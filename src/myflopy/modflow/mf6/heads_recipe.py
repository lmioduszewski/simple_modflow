"""The ``mf.Heads`` recipe: a finished run's heads, resolved at declaration and at build.

``mf.ic(strt=mf.Heads(old_run))`` starts a model from the heads another run ENDED with, on
whatever grid the build supplies. Two steps, each in its own function:

* :func:`stamp_heads`, at DECLARATION (``mf.ic`` calls it, as it resolves every recipe's
  path there, ledger 293): finds the head file the old model's output control names, the
  binary grid file MF6 wrote beside its discretization, each stress period's length and
  the CRS, hashes both files, and checks ``per`` against the file -- so a stopped run, a
  period whose end was not saved or a step MF6 did not converge on is refused at the call.
* :func:`resolve_heads`, at BUILD: reads the head file and the grid file -- and, where the
  run's ``mfsim.nam`` still sits beside them, its period lengths and ``mfsim.lst``, so
  ``per`` and the convergence check follow the run as it then is -- never loading the old
  simulation, and carries the field onto the new grid through the engine
  ``model.hds.remap`` uses, checking the layer map against the built model's layer count
  and leaving cells the new model removes out of what is counted and warned about.

Ledger 596.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np

from myflopy._defaults import DEFAULT_CRS
from myflopy._logging import get_logger
from myflopy.grid_identity import source_fingerprint
from myflopy.modflow.mf6.grid.remap import _field_rows, remap_field
from myflopy.modflow.mf6.grid.voronoi import grid_view_from_modelgrid
from myflopy.modflow.mf6.headsplus import _run_directory, read_saved_field, saved_record
from myflopy.modflow.mf6.output_files import (
    _declared_field_output,
    _declared_models,
    binary_grid_file,
    discretization_file,
    recorded_crs,
)
from myflopy.modflow.mf6.period_end import tdis_period_records
from myflopy.sources import Heads

logger = get_logger(__name__)

#: The name ``mf.ic`` reaches this through, for refusals and the warning.
VERB = "mf.Heads"


def _simulation_directory(folder: Path) -> Path:
    """The directory holding ``mfsim.nam``: ``folder`` or, for a model of a coupled run, its parent."""

    for candidate in (folder, folder.parent):
        if (candidate / "mfsim.nam").is_file():
            return candidate
    raise FileNotFoundError(
        f"{VERB}: {folder} is not a finished MF6 run -- there is no mfsim.nam in it or above "
        "it. Give the run: a Run, a model read from one (run.model(), "
        "mf.load_mf6_run(ws)), or the run's directory."
    )


def _period_lengths(simulation_dir: Path) -> tuple[float, ...]:
    """Each stress period's length from the run's TDIS file, loaded alone (no simulation).

    FloPy reads the TDIS package by itself in milliseconds (measured 2-40 ms on Ten
    Trails runs) where loading the simulation took 1-3 s; the file is named in
    ``mfsim.nam``'s TIMING block.
    """

    import flopy

    namefile = (simulation_dir / "mfsim.nam").read_text(encoding="utf-8", errors="ignore")
    tdis_file = None
    for line in namefile.splitlines():
        fields = line.split("#", 1)[0].split()
        if len(fields) >= 2 and fields[0].upper() == "TDIS6":
            tdis_file = fields[1]
            break
    if tdis_file is None:
        raise ValueError(f"{VERB}: {simulation_dir / 'mfsim.nam'} names no TDIS6 file.")
    simulation = flopy.mf6.MFSimulation(sim_ws=str(simulation_dir), verbosity_level=0)
    tdis = flopy.mf6.ModflowTdis(simulation, filename=tdis_file, pname="tdis")
    tdis.load(strict=False)
    records = tdis_period_records(tdis)
    if not records:
        raise ValueError(f"{VERB}: {simulation_dir / tdis_file} declares no stress periods.")
    return tuple(float(record[0]) for record in records)


def stamp_heads(recipe: Heads) -> Heads:
    """Resolve the run a :class:`~myflopy.sources.Heads` names into files, at DECLARATION.

    Parameters
    ----------
    recipe : Heads
        As ``mf.Heads(...)`` builds it: ``path`` the run's directory or a model's folder.
        One already resolved (``recipe.resolved``, as a reloaded project holds it) is
        returned as it is.

    Returns
    -------
    Heads
        The recipe with ``path`` the head file (absolute), ``grid_path`` the binary grid
        file, ``model``, ``perlen``, ``crs`` and both content hashes filled in.

    Raises
    ------
    FileNotFoundError
        No ``mfsim.nam`` at the run; the head file or the binary grid file is missing.
    ValueError
        The run has no GWF model of that name, several with none named, or a DISU; its
        output control declares no head file; ``per`` is outside the periods, its end was
        not saved, or MF6 did not converge on it.
    """

    if recipe.resolved:
        return recipe
    folder = Path(recipe.path)
    if not folder.exists():
        raise FileNotFoundError(f"{VERB}: no run at {folder}.")
    folder = folder.resolve()
    if folder.is_file():
        raise ValueError(
            f"{VERB}({folder}): give the run -- its directory, a Run or a model -- not a "
            "file; the recipe finds the head file and the grid file itself."
        )
    simulation_dir = _simulation_directory(folder)
    models = _declared_models(simulation_dir)
    flow = [name for kind, name in models if kind.startswith("gwf")]
    model = recipe.model
    if model is None:
        if len(flow) != 1:
            raise ValueError(
                f"{VERB}: the run at {simulation_dir} has "
                + (f"{len(flow)} GWF models ({', '.join(flow)})" if flow else "no GWF model")
                + "; name the one whose heads to read: mf.Heads(run, model=...)."
            )
        model = flow[0]
    elif model not in flow:
        kinds = dict((name, kind) for kind, name in models)
        raise ValueError(
            f"{VERB}(model={model!r}): "
            + (f"that model is {kinds[model].upper()}, and a head is a GWF model's"
               if model in kinds else f"the run has no model of that name ({', '.join(kinds)})")
            + "."
        )
    head_file = _declared_field_output(simulation_dir, model, "head", caller=VERB)
    if not head_file.is_file():
        raise FileNotFoundError(
            f"{VERB}: model {model!r}'s head file {head_file} does not exist -- the run has "
            "not been executed, or saved no heads."
        )
    kind, _dis_file = discretization_file(simulation_dir, model)
    if kind == "disu6":
        raise ValueError(
            f"{VERB}: model {model!r} is DISU, and MF6 writes a DISU's nodes as one layer, so "
            "its heads have no plan-view layer to carry. Carry a DIS or DISV model's heads."
        )
    grid_file, named = binary_grid_file(simulation_dir, model)
    if named == "nogrb":
        raise ValueError(
            f"{VERB}: model {model!r}'s discretization says NOGRB, so MF6 wrote no binary "
            "grid file, and the recipe reads the old mesh from one. Rerun without NOGRB."
        )
    if grid_file is None or not grid_file.is_file():
        raise FileNotFoundError(
            f"{VERB}: model {model!r}'s binary grid file "
            + (f"{grid_file} " if grid_file is not None else "")
            + "does not exist"
            + (" (named by its GRB6 FILEOUT)" if named == "fileout" else "")
            + " -- MF6 writes it when the model runs. The recipe reads the old mesh from it."
        )
    perlen = _period_lengths(simulation_dir)
    import flopy

    reader = flopy.utils.HeadFile(str(head_file), text="head")
    try:
        saved_record(
            reader, per=recipe.per, perlen=perlen, label="head", run_dir=simulation_dir
        )
        if recipe.layers is not None:
            # An old layer the run does not have is knowable here; refused with the
            # engine's own words so the two cannot drift apart.
            _field_rows(recipe.layers, int(reader.nlay), verb=VERB)
    finally:
        reader.close()
    crs = recipe.crs or recorded_crs(simulation_dir, model) or DEFAULT_CRS
    return replace(
        recipe,
        path=head_file,
        model=model,
        grid_path=grid_file,
        perlen=perlen,
        crs=crs,
        source_sha=source_fingerprint(head_file),
        grid_sha=source_fingerprint(grid_file),
    )


def _old_grid(recipe: Heads):
    """The old mesh, read from the run's binary grid file."""

    from flopy.mf6.utils import MfGrdFile

    modelgrid = MfGrdFile(str(recipe.grid_path), verbose=False).modelgrid
    return grid_view_from_modelgrid(modelgrid, crs=recipe.crs, name=recipe.model or "old")


def _current_period_lengths(recipe: Heads) -> tuple[float, ...]:
    """The run's period lengths as they NOW are, where its ``mfsim.nam`` is still beside it.

    The build reads the head file as it now is, so ``per`` (and ``-1`` especially) is
    counted on the same run's TDIS: a run re-run with more periods after the declaration
    would otherwise have its old last period read as its last (the review of ledger 596).
    A bundled copy has no ``mfsim.nam`` beside it, and the lengths stored at the
    declaration -- the ones that copy was written with -- are used.
    """

    folder = Path(recipe.path).parent
    for candidate in (folder, folder.parent):
        if (candidate / "mfsim.nam").is_file():
            return _period_lengths(candidate)
    return tuple(recipe.perlen)


def _layer_rows(recipe: Heads, old_nlay: int, nlay: int) -> list[int]:
    """The OLD layer each layer of the new model reads, checked against ``nlay``."""

    layers = recipe.layers
    if layers is None:
        if old_nlay != nlay:
            example = [min(row, old_nlay - 1) for row in range(nlay)]
            raise ValueError(
                f"{VERB}: the old run has {old_nlay} layer(s) and this model has {nlay}, so "
                "layer for layer does not say which old layer fills each new one. Name one "
                f"old layer per new layer: mf.Heads(..., layers={example}) for instance."
            )
        return list(range(nlay))
    if isinstance(layers, int):
        if nlay != 1:
            raise ValueError(
                f"{VERB}(layers={layers}): this model has {nlay} layers, and one old layer "
                f"fills a ONE-layer model. To fill every layer from it: layers={[layers] * nlay}."
            )
        return [layers]
    if len(layers) != nlay:
        raise ValueError(
            f"{VERB}(layers={list(layers)}): one old layer per layer of this model, which has "
            f"{nlay}; {len(layers)} given."
        )
    return list(layers)


def resolve_heads(recipe: Heads, grid, *, nlay: int, domain=None) -> np.ndarray:
    """The ``(nlay, ncpl)`` starting heads a :class:`~myflopy.sources.Heads` gives, at BUILD.

    Parameters
    ----------
    recipe : Heads
        A recipe :func:`stamp_heads` resolved.
    grid : VoronoiGridPlus
        The new model's built grid.
    nlay : int
        The new model's layer count; the layer map is checked against it.
    domain : array_like, optional
        The new model's idomain, ``(nlay, ncpl)`` or ``(ncpl,)``: a cell at ``<= 0`` is
        never counted, warned about or refused as uncovered. Any other shape is refused,
        naming both. Default ``None``: every cell counts.

    Returns
    -------
    numpy.ndarray
        ``(nlay, grid.ncpl)`` heads, finite in every cell.

    Raises
    ------
    ValueError
        The recipe was not resolved (it reached a build without ``mf.ic``); ``domain`` does
        not fit the model; the layer map does not fit ``nlay``; or as ``model.hds.remap``
        refuses.
    """

    if not recipe.resolved:
        raise ValueError(
            f"{VERB} reached a build unresolved: declare it through "
            "mf.ic(strt=mf.Heads(...)), which finds the run's files."
        )
    active = None
    if domain is not None:
        cells = np.asarray(domain)
        if cells.shape not in ((grid.ncpl,), (int(nlay), grid.ncpl)):
            raise ValueError(
                f"{VERB}: the model's domain has shape {cells.shape}, and this model's built "
                f"grid is ({int(nlay)}, {grid.ncpl}); the context's domain= does not fit the "
                "discretization it is declared beside."
            )
        active = np.broadcast_to(cells.reshape(-1, grid.ncpl), (int(nlay), grid.ncpl)) > 0
    import flopy

    old = _old_grid(recipe)
    reader = flopy.utils.HeadFile(str(recipe.path), text="head")
    try:
        raw, _key, when = read_saved_field(
            reader,
            per=recipe.per,
            perlen=_current_period_lengths(recipe),
            ncpl=int(old.ncpl),
            label="head",
            run_dir=_run_directory(recipe.path),
        )
    finally:
        reader.close()
    rows = _layer_rows(recipe, raw.shape[0], int(nlay))
    field, _report = remap_field(
        raw,
        old,
        grid,
        layers=rows,
        method=recipe.method,
        uncovered=recipe.uncovered,
        water_table=True,
        active=active,
        label="head",
        verb=VERB,
        when=when,
    )
    return field

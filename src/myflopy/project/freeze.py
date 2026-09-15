"""Freeze a resolved simulation so a run can be rebuilt without its sources.

A run's manifest stores a *recipe* -- ``SimulationSpec.to_dict()`` -- and a recipe
cannot hold an array: ``to_dict`` refuses ``numpy.ndarray`` outright. So a model
carrying an inline ``k`` or ``strt`` array writes ``spec_recipe=None``, and
``reopen_run(name).spec`` comes back ``None``: **the run cannot be rebuilt at all**.
Everything expressed as ``mf.ref`` into the package library survives, because the
arrays live in the library; everything inline does not.

This module is the snapshot that closes that. At the end of ``Run.build`` the
resolved simulation is written to ``<run>/.myflopy/frozen.npz``: the structure as
JSON, every array beside it, and the grid's own vertices, cells, CRS and layer
elevations. A reopened run thaws that instead of re-resolving, and reads no source.

**One file, not a directory.** ``os.replace`` of a single file is atomic, so an
interrupted or concurrent freeze can never leave a half-written snapshot, nor delete
a good one it then fails to replace. A directory needed two operations and a shared
staging name, and was neither.

**Written at BUILD, not at ``prepare_run``**, which saves its manifest before
anything is resolved -- and DISCARDED when a run directory is reused, because a
snapshot that outlives the simulation it describes is worse than none: it would
rebuild the previous model in silence.

Measured on a real 9,400-cell, 72-period model: 2.78 MB, against 22.2 MB of MF6
input files and a 1,023 MB run. Binary beats MF6's ~15 characters per number, and a
repetitive array compresses to almost nothing.

The format is deliberately NOT the spec format. ``to_dict``/``from_dict`` carry the
structure, and this wraps them: arrays are lifted out to a token first and put back
after, so the spec serializer's own contract is untouched and nothing here has to
know what a package is.
"""

from __future__ import annotations

import json
import os
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np

from myflopy._logging import get_logger
from myflopy.specs import PackageSpec, SimulationSpec, _json_value

logger = get_logger(__name__)

#: Where a run keeps its snapshot, relative to the run workspace.
FROZEN_FILE = Path(".myflopy") / "frozen.npz"

#: Marks a value that was lifted into the array store.
_TOKEN = "$npz"

#: Marks a caller's own dict that happens to look like a token.
_LITERAL = "$literal"

#: The npz member holding the JSON structure.
_STRUCTURE = "__structure__"

#: The snapshot format's own version, so a reader can refuse what it cannot read
#: rather than mis-reading it. Bump when the on-disk shape changes incompatibly.
FORMAT = 2


class FrozenSnapshotError(ValueError):
    """A snapshot cannot be written, or cannot be trusted to describe this run."""


def _versions() -> dict[str, str]:
    """The toolchain that produced a snapshot, for reading a surprise later."""

    import flopy

    from myflopy import __version__ as myflopy_version

    return {
        "myflopy": str(myflopy_version),
        "flopy": str(flopy.__version__),
        "numpy": str(np.__version__),
    }


def _refuse_object_array(value: np.ndarray, where: str) -> None:
    """Refuse an array numpy can only store by pickling it.

    ``savez`` is permissive and ``load`` is strict: an object-dtype array pickles
    without complaint on the way out, and reading a pickle can execute code, so the
    way back in refuses it. The first version of this module wrote such snapshots
    happily and could never read one -- and FloPy's own list-BC recarray is exactly
    that shape, ``[('cellid', 'O'), ('head', '<f8'), ('boundname', 'O')]``, because a
    cellid is a tuple and a boundname is a string.

    Refusing here puts the failure where it can be acted on: at the build that would
    have written it, naming the package and the option, rather than on a rebuild
    months later as a bare numpy message naming nothing.
    """

    if not value.dtype.hasobject:
        return
    names = getattr(value.dtype, "names", None)
    detail = (
        f" Its object field(s): {[n for n in names if value.dtype[n].hasobject]}."
        if names
        else ""
    )
    raise FrozenSnapshotError(
        f"{where} holds an object-dtype array, which cannot go into a frozen "
        "snapshot: numpy can only store one by pickling, and a snapshot is read back "
        f"with pickling disabled.{detail} Give it a list of records, or a recarray "
        "whose fields are all numeric or fixed-width strings."
    )


def _lift(value: Any, store: dict[str, np.ndarray], where: str) -> Any:
    """Replace every array in ``value`` with a token, collecting them in ``store``.

    Recurses through the containers a package's options actually use -- dicts keyed
    by stress period, lists of records, tuples of cellids -- and leaves everything
    else alone.
    """

    if isinstance(value, np.ndarray):
        _refuse_object_array(value, where)
        key = f"a{len(store)}"
        store[key] = value
        return {_TOKEN: key}
    if isinstance(value, dict):
        if len(value) == 1 and (_TOKEN in value or _LITERAL in value):
            # A caller's own dict that happens to look like a token would be eaten by
            # `_drop` and come back as somebody else's array. Escape it instead.
            return {_LITERAL: {str(k): _lift(v, store, where) for k, v in value.items()}}
        return {
            str(key): _lift(item, store, f"{where}[{key!r}]")
            for key, item in value.items()
        }
    if isinstance(value, (list, tuple)):
        lifted = [_lift(item, store, where) for item in value]
        return type(value)(lifted) if isinstance(value, tuple) else lifted
    return value


def _drop(value: Any, arrays: dict[str, np.ndarray]) -> Any:
    """Put the lifted arrays back where :func:`_lift` took them from."""

    if isinstance(value, dict):
        if len(value) == 1 and _TOKEN in value:
            key = value[_TOKEN]
            if key not in arrays:
                raise FrozenSnapshotError(
                    f"The frozen snapshot references array {key!r}, which it does not "
                    "contain. The snapshot is incomplete and cannot be thawed."
                )
            return arrays[key]
        if len(value) == 1 and _LITERAL in value:
            return {key: _drop(item, arrays) for key, item in value[_LITERAL].items()}
        return {key: _drop(item, arrays) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        dropped = [_drop(item, arrays) for item in value]
        return type(value)(dropped) if isinstance(value, tuple) else dropped
    return value


def _lift_package(package: Any, store: dict[str, np.ndarray], model: str) -> Any:
    """A copy of ``package`` whose options and inputs hold tokens, not arrays.

    A simulation's package tuple holds more than ``PackageSpec``: an ``mf.ref`` is a
    ``PackageRef``, which carries a key and no options at all. Those pass through
    untouched -- their arrays live in the package library, and the run's own copies
    under ``_library/`` are what make them reproducible (ledger 230).
    """

    if not isinstance(package, PackageSpec):
        return package
    where = f"{model}.{package.name}"
    return replace(
        package,
        options=_lift(package.options, store, f"{where} option"),
        inputs=_lift(package.inputs, store, f"{where} input"),
    )


def _drop_package(package: Any, arrays: dict[str, np.ndarray]) -> Any:
    """The inverse of :func:`_lift_package`."""

    if not isinstance(package, PackageSpec):
        return package
    return replace(
        package,
        options=_drop(package.options, arrays),
        inputs=_drop(package.inputs, arrays),
    )


def _grid_block(grid, store: dict[str, np.ndarray], tag: str) -> dict | None:
    """The geometry a thawed model needs, or ``None`` when there is none to keep."""

    if grid is None:
        return None
    try:
        props = grid.get_disv_gridprops()
        verts = np.asarray(grid.verts, dtype=float)
        iverts = [list(map(int, cell)) for cell in grid.iverts]
        # The grid's OWN cell centres, read from the gridprops the DISV was built
        # from -- NOT `grid.centroids`, which are the polygon centroids. For a
        # Voronoi cell those are different points: rebuilding from the centroids put
        # a thawed grid 24.18 length units off its own DISV, which the build-time
        # consistency guard (ledger 249) refused.
        centres = np.asarray(
            [[float(row[1]), float(row[2])] for row in props["cell2d"]], dtype=float
        )
    except Exception:  # noqa: BLE001 - see below; the tuple cannot be closed
        # Deliberately broad, and it only ever costs the snapshot its grid. `grid` is
        # whatever object a caller put on a context -- a user's builder result, a
        # structured grid, a mesh rebuilt from a run folder -- so reading its geometry
        # reaches flopy, numpy and any property it computes lazily. A snapshot without
        # a grid still rebuilds every package, which is the point.
        logger.debug("could not snapshot the grid for %s", tag, exc_info=True)
        return None
    store[f"{tag}.verts"] = verts
    store[f"{tag}.centres"] = centres
    return {
        "verts": f"{tag}.verts",
        "centres": f"{tag}.centres",
        "iverts": iverts,
        "crs": str(getattr(grid, "crs", "") or "") or None,
    }


def _column_label(name: Any) -> int | str:
    """Return ``name`` as JSON carries it, keeping an integer label an integer.

    The elevation table's contract is INTEGER columns -- ``specs._elevation_table_from_layers``
    builds ``{0: top, 1..nlay: bottoms}`` and ``GeoPackageSource.surface_value`` looks
    them up by integer -- so ``str(name)`` thawed ``[0, 1, 2]`` back as ``['0', '1',
    '2']`` and a reopened run died on its first elevation read, which is the one thing
    this block exists to prevent (ledger 265). JSON has no integer KEY but it has an
    integer list ITEM, and this is a list.

    Anything that is not an int or a str after conversion falls back to ``str`` -- the
    old behaviour, and no worse than it. A `MultiIndex` label converts to a tagged
    tuple, which is a dict, and a dict is not hashable as a column name; and a label
    the converter refuses outright would otherwise raise inside the caller's broad
    ``except`` and drop the surfaces entirely. Neither is an improvement on a label
    that is merely spelled differently.
    """

    try:
        converted = _json_value(name)
    except (ValueError, TypeError):
        # `ValueError` is the converter's own refusal; `TypeError` is what a label
        # whose `.item()` will not unwrap raises out of it. Either way the answer is
        # the label spelled as a string, which is what this line did before.
        return str(name)
    return converted if isinstance(converted, (int, str)) and not isinstance(converted, bool) else str(name)


def _surfaces_block(surfaces, store: dict[str, np.ndarray], tag: str) -> dict | None:
    """The per-cell layer elevations, as columns of numbers.

    ``context.surfaces`` is the table ``specs._layer_elevations`` reads -- the single
    function SFR, LAK and ``CellSurfaceOffset`` all resolve elevations through
    (ledger 209) -- so a snapshot that drops it rebuilds today's packages correctly
    and then refuses tomorrow's edits. Only the NUMERIC columns are stored: the
    geometry is the grid's, and is put back in first position on thaw because
    ``gdf_topbtm`` is geometry-first by contract (ledger 196).
    """

    if surfaces is None:
        return None
    try:
        numeric = surfaces.drop(columns="geometry", errors="ignore")
        values = np.asarray(numeric.to_numpy(), dtype=float)
        columns = [_column_label(name) for name in numeric.columns]
    except Exception:  # noqa: BLE001 - see below; the tuple cannot be closed
        # Same reasoning as the grid. `surfaces` is documented loosely ("gdf_topbtm
        # or equivalent") and may be a frame, a LayerBuildResult or a bare array, so
        # reading it reaches pandas, geopandas and numpy. Losing it costs the
        # snapshot its elevations and nothing else.
        logger.debug("could not snapshot the surfaces for %s", tag, exc_info=True)
        return None
    store[f"{tag}.surfaces"] = values
    return {"values": f"{tag}.surfaces", "columns": columns}


def _thaw_grid(block: dict | None, arrays: dict[str, np.ndarray]):
    """Rebuild a grid from :func:`_grid_block`, or ``None`` if it was not kept."""

    if not block:
        return None
    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus

    missing = [key for key in (block["verts"], block["centres"]) if key not in arrays]
    if missing:
        raise FrozenSnapshotError(
            f"The frozen snapshot references grid array(s) {missing}, which it does "
            "not contain. The snapshot is incomplete and cannot be thawed."
        )
    return VoronoiGridPlus(
        verts=arrays[block["verts"]],
        iverts=[list(cell) for cell in block["iverts"]],
        xcyc=arrays[block["centres"]],
        crs=block["crs"],
    )


def _thaw_surfaces(block: dict | None, arrays: dict[str, np.ndarray], grid):
    """Rebuild the elevation table, geometry first, and attach it to ``grid``."""

    if not block or grid is None:
        return None
    import geopandas as gpd

    key = block["values"]
    if key not in arrays:
        raise FrozenSnapshotError(
            f"The frozen snapshot references surfaces array {key!r}, which it does "
            "not contain. The snapshot is incomplete and cannot be thawed."
        )
    values = arrays[key]
    columns = block["columns"]
    frame = gpd.GeoDataFrame(
        {name: values[:, index] for index, name in enumerate(columns)},
        geometry=list(grid.gdf_vorPolys.geometry),
        crs=getattr(grid, "crs", None),
    )
    frame = frame[["geometry", *columns]]
    grid.gdf_topbtm = frame
    return frame


def _open(path: Path) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Read a snapshot into its structure and its arrays."""

    with np.load(path, allow_pickle=False) as handle:
        arrays = {key: handle[key] for key in handle.files if key != _STRUCTURE}
        payload = json.loads(str(handle[_STRUCTURE]))
    return payload, arrays


def freeze_simulation(
    spec: SimulationSpec,
    workspace: Path | str,
    *,
    run_name: str | None = None,
    contexts: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Write ``spec`` to ``<workspace>/.myflopy/frozen.npz`` and return a summary.

    ``run_name`` is stamped in, so a thaw can refuse a snapshot belonging to another
    run rather than rebuilding somebody else's model.

    ``contexts`` maps model name to the context the build RESOLVED, and is what makes
    a deferred grid survive. A declared ``ModelSpec`` carries the grid RECIPE, not a
    grid: ``ModelSpec.build`` resolves it into a local context and hands that to the
    built model, leaving the spec's own ``context.grid`` as ``None``. Freezing the
    declared spec therefore captured no grid at all for exactly the models whose grid
    is expensive to rebuild -- and left the recipe in place, so a reopened run went
    back to the boundary and refinement files it was supposed to no longer need.
    """

    workspace = Path(workspace)
    target = workspace / FROZEN_FILE
    staging = target.with_name(target.name + ".writing.npz")

    store: dict[str, np.ndarray] = {}
    lifted = replace(
        spec,
        packages=tuple(_lift_package(p, store, spec.name) for p in spec.packages),
        models=tuple(
            replace(
                model,
                packages=tuple(
                    _lift_package(p, store, model.name) for p in model.packages
                ),
            )
            for model in spec.models
        ),
    )
    grids, domains, surfaces = {}, {}, {}
    resolved = contexts or {}
    for model in spec.models:
        context = resolved.get(model.name, model.context)
        grids[model.name] = _grid_block(
            getattr(context, "grid", None), store, f"grid.{model.name}"
        )
        domain = getattr(context, "domain", None)
        if domain is not None:
            # `_refuse_object_array` is reached from `_lift` for every package option
            # and `_grid_block`/`_surfaces_block` force `dtype=float`, which left this
            # as the one write into the store with neither guard. `ModelContext.domain`
            # is typed `Any` and validated nowhere, so anything reaches it: a
            # `LayerBuildResult` passed instead of its `.idomain`, or a list holding a
            # `None`, becomes an object array. Measured, that froze "successfully" --
            # 3 arrays, 3,885 bytes -- and made the WHOLE snapshot unreadable, every
            # package in it, with `frozen_summary` reporting `readable: False` and no
            # cause (ledger 265).
            values = np.asarray(domain)
            _refuse_object_array(values, f"the model domain for {model.name!r}")
            store[f"domain.{model.name}"] = values
            domains[model.name] = f"domain.{model.name}"
        surfaces[model.name] = _surfaces_block(
            getattr(context, "surfaces", None), store, f"surf.{model.name}"
        )

    payload = {
        "format": FORMAT,
        "run": run_name,
        "simulation": spec.name,
        "versions": _versions(),
        "spec": lifted.to_dict(),
        "grids": grids,
        "domains": domains,
        "surfaces": surfaces,
    }

    target.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(staging, **{_STRUCTURE: np.array(json.dumps(payload))}, **store)
    os.replace(staging, target)     # one atomic step, whoever else is writing

    size = target.stat().st_size
    logger.info(
        "froze %s: %d array(s), %.2f MB at %s",
        spec.name, len(store), size / 1048576, target,
    )
    return {"path": str(target), "arrays": len(store), "format": FORMAT, "bytes": size}


def discard_snapshot(workspace: Path | str) -> bool:
    """Delete a run's snapshot; ``True`` if there was one.

    Called when a run directory is REUSED. A snapshot that outlives the simulation it
    describes is worse than no snapshot: ``prepare_run(overwrite=True)`` keeps the
    directory, so without this a re-prepared run would thaw and rebuild the PREVIOUS
    model in silence -- and a run staged with ``build=False`` would carry an answer it
    was never asked for.
    """

    path = Path(workspace) / FROZEN_FILE
    if not path.exists():
        return False
    path.unlink()
    logger.debug("discarded the frozen snapshot at %s", path)
    return True


def frozen_summary(workspace: Path | str) -> dict[str, Any] | None:
    """What a run's snapshot says about itself, or ``None`` when it has none."""

    path = Path(workspace) / FROZEN_FILE
    if not path.exists():
        return None
    try:
        payload, _arrays = _open(path)
    except Exception:  # noqa: BLE001 - see below; the tuple cannot be closed
        # Deliberately broad: this is the "what have I got?" verb, reached from
        # `Run.frozen` and from summaries, and a corrupt or truncated file reaches
        # numpy, zipfile and json. Reporting it as unreadable IS the answer; raising
        # from an inspection accessor is not.
        logger.debug("could not read the frozen snapshot at %s", path, exc_info=True)
        return {"path": str(path), "readable": False}
    return {
        "path": str(path),
        "readable": True,
        "format": payload.get("format"),
        "run": payload.get("run"),
        "simulation": payload.get("simulation"),
        "versions": payload.get("versions", {}),
        "models": [model["name"] for model in payload.get("spec", {}).get("models", [])],
        "bytes": path.stat().st_size,
    }


def thaw_simulation(
    workspace: Path | str, *, run_name: str | None = None
) -> SimulationSpec | None:
    """Rebuild the resolved simulation from a run's snapshot, or ``None`` if absent.

    Reads no source file: every array came out of the snapshot and every cell number
    was already resolved when it was taken.
    """

    path = Path(workspace) / FROZEN_FILE
    if not path.exists():
        return None
    payload, arrays = _open(path)

    found = payload.get("format")
    if found != FORMAT:
        raise FrozenSnapshotError(
            f"The frozen snapshot at {path} is format {found!r}; this myflopy reads "
            f"format {FORMAT}. Prepare the run again to write a current one."
        )
    stamped = payload.get("run")
    if run_name is not None and stamped is not None and stamped != run_name:
        raise FrozenSnapshotError(
            f"The frozen snapshot at {path} was written for run {stamped!r}, not "
            f"{run_name!r}. It describes a different simulation and would rebuild that "
            "one. Prepare this run again to write its own."
        )

    spec = SimulationSpec.from_dict(payload["spec"])
    grids = payload.get("grids", {})
    domains = payload.get("domains", {})
    surfaces = payload.get("surfaces", {})
    models = []
    for model in spec.models:
        grid = _thaw_grid(grids.get(model.name), arrays)
        table = _thaw_surfaces(surfaces.get(model.name), arrays, grid)
        domain_key = domains.get(model.name)
        context = replace(
            model.context,
            grid=grid if grid is not None else model.context.grid,
            domain=arrays[domain_key] if domain_key in arrays else model.context.domain,
            surfaces=table if table is not None else model.context.surfaces,
        )
        models.append(
            replace(
                model,
                context=context,
                # The snapshot's grid is authoritative, so drop the RECIPE that would
                # rebuild one. `ModelSpec.build` prefers `self.grid` and re-resolves it
                # unconditionally -- reading the boundary and refinement files again --
                # so leaving it in place meant a deferred-grid run still died on a
                # missing source with the snapshot sitting right there. A reopened run
                # rebuilds what it RAN; re-resolving is what a new prepare_run is for.
                grid=None if grid is not None else model.grid,
                packages=tuple(_drop_package(p, arrays) for p in model.packages),
            )
        )
    return replace(
        spec,
        models=tuple(models),
        packages=tuple(_drop_package(p, arrays) for p in spec.packages),
    )


__all__ = [
    "FORMAT",
    "FROZEN_FILE",
    "FrozenSnapshotError",
    "discard_snapshot",
    "freeze_simulation",
    "frozen_summary",
    "thaw_simulation",
]

"""Native pyEMU parameterization for the modern ``myflopy`` PEST facade.

This module is the engine behind :meth:`PestProject.parameterize`. Where the
legacy path (``parameters.py``) hand-rolled template files and a custom forward
run, this path compiles a small, readable spec straight to
``pyemu.utils.PstFrom.add_parameters`` and lets pyEMU's own
``apply_list_and_array_pars`` drive the forward run.

The user writes::

    cal.parameterize("k", style="constant", bounds=(0.2, 5), physical=(1e-3, 100))
    cal.parameterize("recharge", bounds=(0.5, 1.5))
    cal.parameterize("ghb.cond", bounds=(0.1, 10))

and each call expands into the correct ``add_parameters`` invocation against the
external array/list files produced by ``sim.set_all_data_external()``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from myflopy.modflow.mf6.pest.model_lookup import resolve_transport_model_name
from myflopy.modflow.mf6.pest.zones import resolve_zone_array

# --- recipe registry -------------------------------------------------------
#
# Each recipe knows how to find the external MODFLOW input file(s) a target maps
# to and how pyEMU should read them. Recipes are written for DISV list files,
# which index on ``[layer, cell]`` (index_cols=[0, 1]); ``use_col`` is the
# zero-based column holding the value pyEMU multiplies/offsets, on DISV. A DIS
# cellid is three columns, ``layer row col``, so on a structured grid every
# column after it moves one to the right -- :func:`list_columns` is the ONE rule
# that says where they are on the model's own grid (ledger 420).


@dataclass(frozen=True)
class _Recipe:
    """How one calibration target maps onto an external MODFLOW input file."""

    canonical: str
    family: str  # "array" or "list"
    file_stub: str  # formatted with model name; may contain a "*" glob
    use_col: int | None = None  # list family only
    default_style: str = "constant"
    additive: bool = False  # default to additive offsets (e.g. drain elevation)
    # Which model in the simulation owns the file. A PestProject hangs off the
    # FLOW model, so "transport" targets resolve the GWT sibling by name --
    # everything downstream (pyEMU's own apply_list_and_array_pars) is purely
    # filename-driven and needs no other change.
    model: str = "flow"
    # Array family only: where the griddata lives on the FloPy model, so it can
    # be re-stored one array per layer before externalization (see
    # `relayer_array_target`).
    package: str | None = None
    variable: str | None = None
    # List family only: which columns hold the CELL IDENTITY. Most MF6 list
    # packages write ``layer cell <values...>``, so (0, 1) is the default -- but
    # UZF numbers its own cells first and writes ``iuzno layer cell ...``. pyEMU
    # geolocates a list parameter from these columns, and getting them wrong is
    # silent: every parameter lands on one cell's coordinates and only the
    # geostatistics are wrong (see the `uzf.vks` entry below).
    index_cols: tuple[int, ...] = (0, 1)


# DISV list files are usually written as ``layer cell <values...>``; see
# `_Recipe.index_cols` for the exceptions.
_LIST_INDEX_COLS = (0, 1)

#: How many columns an MF6 cellid takes in a list file, by FloPy's grid type:
#: ``layer row col`` (DIS), ``layer cell`` (DISV), ``node`` (DISU).
_CELLID_WIDTH = {"structured": 3, "vertex": 2, "unstructured": 1}


def list_columns(recipe: _Recipe, grid_type: str) -> tuple[list[int], int]:
    """Where a list recipe's cell identity and value sit on a grid of ``grid_type``.

    A recipe is written for a DISV file, whose cellid is two columns (``layer
    cell``). On DIS the cellid is three (``layer row col``), so the cell identity
    spans one more column and the value is one further right. Before this the DISV
    columns were used on every grid, and on a DIS model a CHD multiplier was
    written onto the COLUMN NUMBER of each record (``usecol:2`` of ``1 1 1 1.0E+02``)
    -- the build succeeded (ledger 420).

    Parameters
    ----------
    recipe : _Recipe
        A list-family recipe; its ``index_cols`` (two adjacent columns, the DISV
        cellid) and ``use_col`` are read.
    grid_type : str
        FloPy's ``modelgrid.grid_type``: ``"structured"`` (DIS), ``"vertex"``
        (DISV) or ``"unstructured"`` (DISU). Any other value raises ``KeyError``.

    Returns
    -------
    tuple of (list of int, int)
        The zero-based ``index_cols`` -- the cellid's columns, pyEMU reading its
        LAST two as the cell's ``(i, j)`` -- and ``use_col``, both on that grid.
    """

    first, last = recipe.index_cols[0], recipe.index_cols[-1]
    if last - first != len(_LIST_INDEX_COLS) - 1:
        raise ValueError(
            f"Recipe {recipe.canonical!r} names index columns {recipe.index_cols}, "
            "which are not one DISV cellid (two adjacent columns)."
        )
    shift = _CELLID_WIDTH[grid_type] - len(_LIST_INDEX_COLS)
    width = _CELLID_WIDTH[grid_type]
    return list(range(first, first + width)), int(recipe.use_col) + shift

_RECIPES: dict[str, _Recipe] = {
    "k": _Recipe("k", "array", "{model}.npf_k.txt", package="npf", variable="k"),
    "k33": _Recipe("k33", "array", "{model}.npf_k33.txt", package="npf", variable="k33"),
    "recharge": _Recipe("recharge", "list", "{model}.rch_stress_period_data_*.txt", use_col=2),
    "chd": _Recipe("chd", "list", "{model}.chd_stress_period_data_*.txt", use_col=2),
    "ghb.cond": _Recipe("ghb.cond", "list", "{model}.ghb_stress_period_data_*.txt", use_col=3),
    "ghb.bhead": _Recipe("ghb.bhead", "list", "{model}.ghb_stress_period_data_*.txt", use_col=2),
    "drn.cond": _Recipe("drn.cond", "list", "{model}.drn_stress_period_data_*.txt", use_col=3),
    "drn.elev": _Recipe(
        "drn.elev", "list", "{model}.drn_stress_period_data_*.txt", use_col=2, additive=True
    ),
    "wel": _Recipe("wel", "list", "{model}.wel_stress_period_data_*.txt", use_col=2),
    # UZF saturated vertical K, from PACKAGEDATA. Two things differ from every
    # recipe above and both are load-bearing:
    #
    # * UZF numbers its own cells, so the row is ``iuzno layer icell2d landflag
    #   ivertcon surfdep vks ...`` -- the cell identity is at columns (1, 2), not
    #   (0, 1). With the default, pyEMU reads the LAYER as the cell number and
    #   every UZF parameter is placed at cell 0's coordinates: the .pst builds,
    #   the forward run is correct, and only the geostatistics are garbage,
    #   surfacing much later as `error inverting cov` in the prior draw.
    # * `boundname` is declared LAST in the MF6 dfn, so enabling it appends a
    #   column and leaves `vks` at 6 (measured both ways, 2026-07-30).
    #
    # Only vks: it is the UZF quantity heads actually respond to (`vks x3` moves
    # heads 1.026 ft on the canonical testing profile, against 0.011 ft for
    # `finf x2`), and thts/eps are 0.995-collinear with each other -- shipping
    # them as a set would build the same ill-posed problem as K/porosity
    # (ledger 112).
    "uzf.vks": _Recipe(
        "uzf.vks", "list", "{model}.uzf_packagedata.txt", use_col=6,
        index_cols=(1, 2),
    ),
    # Transport. MST writes porosity as ONE whole-grid array (no LAYERED
    # keyword), so unlike npf_k there is no per-layer file to select from --
    # see the `layers=` guard in `add_native_parameter`.
    "porosity": _Recipe(
        "porosity", "array", "{model}.mst_porosity.txt", model="transport",
        package="mst", variable="porosity",
    ),
}

# Friendly aliases -> canonical key.
_ALIASES: dict[str, str] = {
    "npf.k": "k",
    "npf_k": "k",
    "hk": "k",
    "kh": "k",
    "npf.k33": "k33",
    "kv": "k33",
    "rch": "recharge",
    "rcha": "recharge",
    "ghb": "ghb.cond",
    "drn": "drn.cond",
    "wel.q": "wel",
    "well": "wel",
    "pumping": "wel",
    "mst.porosity": "porosity",
    "mst_porosity": "porosity",
    "n": "porosity",
    "uzf": "uzf.vks",
    "vks": "uzf.vks",
    "uzf_vks": "uzf.vks",
}

# pyEMU ``par_type`` values that need a cell spatial reference (Phase 2 work for
# Voronoi array parameters). These are accepted for list packages, where no
# spatial reference is required.
_SPATIAL_STYLES = {"grid", "pilotpoints"}


def resolve_target(target: str) -> _Recipe:
    """Return the recipe for a friendly target name, raising a helpful error."""

    key = str(target).strip().lower()
    key = _ALIASES.get(key, key)
    if key not in _RECIPES:
        known = ", ".join(sorted(set(_RECIPES) | set(_ALIASES)))
        raise KeyError(f"Unknown calibration target {target!r}. Known targets: {known}.")
    return _RECIPES[key]


@dataclass
class NativeParameterSpec:
    """A declarative, pyEMU-bound parameter request.

    Parameters
    ----------
    target
        Friendly target name, e.g. ``"k"``, ``"recharge"``, ``"ghb.cond"``.
    style
        Spatial style: ``"constant"`` (one multiplier), ``"zone"`` (one per
        ``zones`` value), ``"grid"`` (one per entry/cell), or ``"pilotpoints"``.
        ``None`` uses the recipe default. ``"grid"``/``"pilotpoints"`` on array
        targets need a cell spatial reference (Phase 2).
    bounds
        ``(lower, upper)`` for the adjustable multiplier (or additive offset).
    physical
        ``(ult_lbound, ult_ubound)`` -- hard limits on the *final* model value
        after all multipliers are applied. Strongly recommended for multipliers.
    transform
        ``"log"`` (default) or ``"none"``. Forced to ``"none"`` for additive.
    additive
        Apply as an additive offset rather than a multiplier. Defaults to the
        recipe (e.g. drain elevation is additive).
    zones
        Zone array for ``style="zone"`` (and to mask inactive cells).
    correlation
        Variogram range (model length units) for ``grid`` geostatistics.
    anisotropy
        Anisotropy ratio of the ``grid`` variogram -- how many times farther the
        property is correlated along the major axis than across it (``1.0`` =
        isotropic). E.g. ``5`` for K in a buried channel correlated 5x farther
        down-valley than across.
    bearing
        Azimuth (degrees) of the anisotropy major axis. Ignored when
        ``anisotropy`` is ``1.0``.
    nugget
        Nugget (unresolved short-scale variance) of the ``grid`` variogram.
    temporal
        Temporal correlation range (days) for time-varying list packages
        (Phase 2 -- recorded but not yet wired).
    name
        Parameter group / name base. Defaults to a slug of ``target``.
    capture
        Also record the *resolved* model-input field (after multipliers are
        applied) as zero-weight observations, so each realization's per-cell
        property field is carried in the ensembles. Enables spatial parameter
        maps via ``IesResults.plot_field`` / ``IesResults.field``.
    """

    target: str
    style: str | None = None
    bounds: tuple[float, float] = (0.5, 2.0)
    physical: tuple[float, float] | None = None
    transform: str = "log"
    additive: bool | None = None
    zones: Any = None
    layers: tuple[int, ...] | None = None
    pp_space: int | None = None
    pp_points: Any = None
    correlation: float | None = None
    anisotropy: float = 1.0
    bearing: float = 0.0
    nugget: float = 0.0
    temporal: float | None = None
    name: str | None = None
    capture: bool = False
    extra: dict = field(default_factory=dict)
    resolved_files: list[str] = field(default_factory=list)

    def __post_init__(self):
        """Resolve the target recipe and fill defaults for style, additive, and name."""

        self.recipe = resolve_target(self.target)
        if self.style is None:
            self.style = self.recipe.default_style
        self.style = str(self.style).strip().lower()
        if self.additive is None:
            self.additive = self.recipe.additive
        if self.name is None:
            self.name = self.recipe.canonical.replace(".", "")
        if self.style == "pilotpoints" and self.recipe.family == "array" and not (
            self.recipe.package and self.recipe.variable
        ):
            # Pilot points multiply a BASE ARRAY read off the model, so the
            # recipe has to say which one. (Until 2026-07-30 that array was a
            # hardcoded `gwf.npf.k` for every target, which is how `k33` came to
            # be overwritten with horizontal K -- 30x wrong, forward run exit 0.
            # This guard replaces a narrower one that only refused TRANSPORT
            # targets and therefore missed k33 entirely.)
            raise NotImplementedError(
                f"style='pilotpoints' needs to know which model array "
                f"{self.recipe.canonical!r} scales, and its recipe declares no "
                "package/variable. Add them to the `_RECIPES` entry, or use "
                "style='grid' (one geostatistically correlated multiplier per "
                "cell), 'zone', or 'constant'."
            )

    @property
    def resolved_transform(self) -> str:
        """Transform actually used (additive parameters cannot be log)."""

        return "none" if self.additive else self.transform


def _flatten_array_file(path: Path, row_length: int = 1) -> None:
    """Rewrite an MF6 external array file as ``row_length`` values per line.

    MF6 wraps large arrays at a fixed number of values per line, leaving a
    short final line. pyEMU reads array files with ``numpy.loadtxt``, which
    rejects that ragged layout (and which also emits a ``%F`` format warning).
    Rewriting the values in rows of one length (free format, scientific
    notation) is read identically by MF6 and cleanly by pyEMU, whatever the cell
    count.

    The length is the grid's: one value per line on DISV, where pyEMU reads row
    ``i`` as cell ``i``, and one ROW of the grid per line on DIS (``ncol``), where
    pyEMU reads the file as ``(nrow, ncol)`` and looks each value's centre up at
    ``[i, j]`` -- a one-column file there made every array parameter on a DIS
    model fail in ``write_array_tpl`` with ``IndexError: index 5 is out of bounds
    for axis 0 with size 5`` (ledger 420). See :func:`array_row_length`.
    """

    values = np.array(Path(path).read_text().split(), dtype=float)
    np.savetxt(path, values.reshape(-1, int(row_length)), fmt="%.10E")


def array_row_length(modelgrid) -> int:
    """How many values one line of an array file pyEMU reads holds, on ``modelgrid``.

    ``ncol`` on a structured (DIS) grid -- pyEMU indexes its spatial reference at
    ``[row, col]`` -- else ``1``: on DISV and DISU it indexes a cell by its line.

    Parameters
    ----------
    modelgrid : flopy.discretization.Grid
        The staged model's FloPy grid; ``grid_type`` and, on DIS, ``ncol`` are read.

    Returns
    -------
    int
    """

    return int(modelgrid.ncol) if modelgrid.grid_type == "structured" else 1


def model_name_for(project, recipe: _Recipe) -> str:
    """Name of the model whose external input files this recipe reads.

    The project hangs off the FLOW model, so a transport target has to find its
    GWT sibling. Nothing else in the pipeline cares: pyEMU's
    ``apply_list_and_array_pars`` matches on FILENAME and never parses a model
    name or package type, so the forward run needs no transport-specific code.
    """

    if recipe.model == "transport":
        return resolve_transport_model_name(project)
    return project.model.name


def staged_model(project, model_name: str):
    """One model of the simulation the build STAGED from the run's files.

    ``PestProject.build`` loads the simulation from its staging copy of the run
    (``project._staged_simulation``), and everything the build reads off a FloPy
    model -- a recipe's package array (:func:`relayer_array_target`, the
    pilot-point base array), the cell centres pilot points are placed on, the
    default spatial reference -- is read THERE: the numbers the template's model
    holds, as written to its files. They were read off ``project.model``, the
    caller's object, until ledger 356, and the relayering changed that object's
    storage.

    Parameters
    ----------
    project : PestProject
        The calibration being built; its ``_staged_simulation`` is read.
    model_name : str
        The model's name in that simulation (the flow model's, or a transport
        sibling's from :func:`model_name_for`).

    Returns
    -------
    flopy.mf6.MFModel
        The staged model.

    Raises
    ------
    RuntimeError
        Outside ``build()``, when there is no staged simulation to read.
    KeyError
        When the run's files hold no model of that name, naming the ones they hold.
    """

    simulation = getattr(project, "_staged_simulation", None)
    if simulation is None:
        raise RuntimeError(
            f"Model {model_name!r} is read from the copy of the run PestProject.build() "
            "stages, and none is staged: call this during build()."
        )
    model = simulation.get_model(model_name)
    if model is None:
        raise KeyError(
            f"The run's files hold no model named {model_name!r} (they hold "
            f"{list(simulation.model_names)}); the PEST build reads the model from them."
        )
    return model


def flopy_model_for(project, recipe: _Recipe):
    """The FloPy model whose packages this recipe reads: the build's STAGED copy's.

    See :func:`staged_model`; a transport target's model is the GWT sibling.

    Parameters
    ----------
    project : PestProject
        The calibration being built.
    recipe : _Recipe
        The target's recipe; its ``model`` says flow or transport.

    Returns
    -------
    flopy.mf6.MFModel
        The staged model, as :func:`staged_model` returns it.

    Raises
    ------
    RuntimeError
        Outside ``build()``.
    KeyError
        When the run's files hold no model of that name.
    """

    return staged_model(project, model_name_for(project, recipe))


def relayer_array_target(project, recipe: _Recipe) -> bool:
    """Re-store an array target one array per layer. Returns whether it applied.

    MF6 griddata supplied as a scalar (``ModflowGwtmst(gwt, porosity=0.25)``)
    externalizes to ONE file holding ``nlay * ncpl`` values, with no LAYERED
    keyword. pyEMU cannot parameterize that on a DISV grid: ``write_array_tpl``
    looks up ``get_xy([i, j])`` for **every** array row against a spatial
    reference with ``ncpl`` entries, so the file raises ``IndexError: index 441
    is out of bounds`` from inside ``add_parameters`` -- an error naming neither
    the target nor the cause.

    ``make_layered()`` + a per-layer ``set_data`` writes the same numbers in the
    shape ``npf_k`` already has (``<model>.mst_porosity_layer1.txt`` ...), which
    MF6 reads identically.

    It acts on the simulation the build STAGED (:func:`flopy_model_for`), never on
    the caller's model, whose storage it used to change (ledger 110, struck by
    ledger 356).

    Skipped when the array is **already** stored one per layer (which ``npf_k``
    starts out as, and which a run whose arrays an earlier build externalized
    holds) and when the model has a single layer -- there a whole-grid file
    already has exactly ``ncpl`` values, so pyEMU is happy and splitting it would
    be busywork. ``store_internal()`` first because FloPy refuses to make
    EXTERNAL data layered, which is the state a run whose arrays are in external
    files is loaded in.
    """

    if recipe.family != "array" or not recipe.package or not recipe.variable:
        return False
    model = flopy_model_for(project, recipe)
    package = model.get_package(recipe.package)
    array = getattr(package, recipe.variable, None) if package is not None else None
    if array is None or not array.supports_layered():
        return False
    # FloPy exposes no public "is this stored per layer" flag; falling back to
    # False when the private accessor moves means attempting the relayer, which
    # is the safe direction.
    storage = getattr(array, "_get_storage_obj", lambda: None)()
    if getattr(storage, "layered", False):
        return False
    nlay = int(model.modelgrid.nlay)
    if nlay <= 1:
        return False
    values = np.asarray(array.get_data(), dtype=float)
    if values.ndim < 2 or values.shape[0] != nlay:
        values = values.reshape(nlay, -1)
    array.store_internal()
    array.make_layered()
    array.set_data([values[layer] for layer in range(nlay)])
    return True


def _select_layer_files(files: list[str], layers, recipe: _Recipe) -> list[str]:
    """Narrow per-layer array files to the requested zero-based layers.

    Refuses rather than falling through when nothing matches. MF6 externalizes
    some griddata as ONE whole-grid array (MST porosity) and some one file per
    layer (NPF K); on the former there is no per-layer file to pick, and the
    previous ``if selected:`` fallthrough silently parameterized every layer
    while the caller believed ``layers=`` had restricted it.
    """

    wanted = {int(layer) for layer in layers}
    selected = [
        name
        for name in files
        if (match := re.search(r"_layer(\d+)\.txt$", name))
        and (int(match.group(1)) - 1) in wanted
    ]
    if selected:
        return selected
    layered = [name for name in files if re.search(r"_layer(\d+)\.txt$", name)]
    if not layered:
        raise ValueError(
            f"layers={sorted(wanted)} cannot be applied to target "
            f"{recipe.canonical!r}: MODFLOW writes it as a single whole-grid "
            f"array ({files[0]}), not one file per layer, so there is no "
            "per-layer file to select. Drop `layers=` (the parameter covers "
            "the whole grid), or use `zones=` to restrict it spatially."
        )
    available = sorted(
        int(re.search(r"_layer(\d+)\.txt$", name).group(1)) - 1 for name in layered
    )
    raise ValueError(
        f"layers={sorted(wanted)} matched no external input file for target "
        f"{recipe.canonical!r}. Available layers: {available}."
    )


def _resolve_files(template_workspace: Path, model_name: str, recipe: _Recipe) -> list[str]:
    """Return external input filenames (relative to the template) for a recipe."""

    stub = recipe.file_stub.format(model=model_name)
    workspace = Path(template_workspace)
    if "*" in stub:
        matches = sorted(p.name for p in workspace.glob(stub))
    else:
        # Multi-layer DISV/DIS arrays are externalized one file per layer, e.g.
        # ``<model>.npf_k_layer1.txt``. Match those precisely so that resolving
        # ``k`` never picks up ``k33`` (``npf_k_layer*`` requires ``_layer``
        # immediately after ``npf_k``).
        #
        # Per-layer FIRST: relayering an array that was already external leaves
        # the old whole-grid file behind, unreferenced by the package but still
        # on disk. Matching the exact name first would resolve to that stale
        # file -- pointing every parameter at values MODFLOW no longer reads.
        base = stub[:-4] if stub.endswith(".txt") else stub
        matches = sorted(p.name for p in workspace.glob(f"{base}_layer*.txt"))
        if not matches and (workspace / stub).exists():
            matches = [stub]
    if not matches:
        raise FileNotFoundError(
            f"No external input file found for target {recipe.canonical!r} "
            f"(looked for {stub!r} or per-layer '{stub[:-4]}_layer*.txt' in "
            f"{template_workspace}). Did the model build with this package, and "
            "was set_all_data_external() applied?"
        )
    return matches


def add_native_parameter(project, spec: NativeParameterSpec):
    """Compile one :class:`NativeParameterSpec` into ``pf.add_parameters``.

    Returns the parameter dataframe pyEMU produced, and records it on the
    project for later inspection via :meth:`PestProject.settings`.
    """

    recipe = spec.recipe
    files = _resolve_files(
        project.template_workspace, model_name_for(project, recipe), recipe
    )
    if spec.layers is not None and recipe.family == "array":
        files = _select_layer_files(files, spec.layers, recipe)
    spec.resolved_files = list(files)

    grid = flopy_model_for(project, recipe).modelgrid
    if recipe.family == "array":
        # Normalize MF6's wrapped array layout so pyEMU can read any cell count,
        # in the shape it reads on this grid (one grid row per line on DIS).
        for filename in files:
            _flatten_array_file(
                Path(project.template_workspace) / filename, array_row_length(grid)
            )

    if recipe.family == "array" and spec.style == "pilotpoints":
        raise NotImplementedError(
            f"style='pilotpoints' on array target {recipe.canonical!r} is not "
            "wired for Voronoi grids yet. Use style='grid' -- one geostatistically "
            "correlated multiplier per cell, drawn against the model's spatial "
            "reference -- or style='constant'/'zone'."
        )

    kwargs: dict[str, Any] = {
        "par_type": spec.style,
        "par_name_base": spec.name,
        "pargp": spec.name,
        "lower_bound": float(spec.bounds[0]),
        "upper_bound": float(spec.bounds[1]),
        "transform": spec.resolved_transform,
    }
    if spec.physical is not None:
        kwargs["ult_lbound"] = float(spec.physical[0])
        kwargs["ult_ubound"] = float(spec.physical[1])
    if spec.additive:
        kwargs["par_style"] = "a"
    if spec.zones is not None:
        # Normalized per family: pyEMU wants per-cell for array targets and
        # (layer, cell) for list targets, and rejects the other with an error
        # that names neither the target nor the shape. See `pest/zones.py`.
        kwargs["zone_array"] = resolve_zone_array(
            spec.zones, family=recipe.family, model=project.model
        )
    if recipe.family == "list":
        index_cols, use_col = list_columns(recipe, grid.grid_type)
        kwargs["index_cols"] = index_cols
        kwargs["use_cols"] = [use_col]
    if spec.correlation is not None and spec.style in _SPATIAL_STYLES:
        kwargs["geostruct"] = project._geostruct_for(spec)
    kwargs.update(spec.extra)

    # List/constant: pass the file list so one parameter scales every file
    # together (recharge across stress periods, or one constant K multiplier
    # across all layers). Grid array parameters instead get an independent
    # parameter set *per layer file*, so each layer's K field varies on its own.
    if recipe.family == "array" and spec.style == "grid" and len(files) > 1:
        frames = []
        for index, filename in enumerate(files, start=1):
            layer_kwargs = dict(kwargs)
            layer_kwargs["par_name_base"] = f"{spec.name}l{index}"
            layer_kwargs["pargp"] = f"{spec.name}l{index}"
            frames.append(project.pf.add_parameters(filename, **layer_kwargs))
        frame = pd.concat(frames, ignore_index=True)
    else:
        filenames = files if len(files) > 1 else files[0]
        frame = project.pf.add_parameters(filenames, **kwargs)
    project._native_parameter_frames[spec.name] = frame
    return frame

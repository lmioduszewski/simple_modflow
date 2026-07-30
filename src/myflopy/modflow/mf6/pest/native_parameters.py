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
# to and how pyEMU should read them. Voronoi/DISV list files index on
# ``[layer, cell]`` (index_cols=[0, 1]); structured grids would use
# ``[layer, row, col]``. ``use_col`` is the zero-based column holding the value
# pyEMU multiplies/offsets.


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


def _flatten_array_file(path: Path) -> None:
    """Rewrite an MF6 external array file as one value per line.

    MF6 wraps large arrays at a fixed number of values per line, leaving a
    short final line. pyEMU reads array files with ``numpy.loadtxt``, which
    rejects that ragged layout (and which also emits a ``%F`` format warning).
    Rewriting the values in a single column (free format, scientific notation)
    is read identically by MF6 and cleanly by pyEMU, regardless of cell count.
    """

    values = np.array(Path(path).read_text().split(), dtype=float)
    np.savetxt(path, values.reshape(-1, 1), fmt="%.10E")


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


def flopy_model_for(project, recipe: _Recipe):
    """The FloPy model object whose packages this recipe reads."""

    if recipe.model == "transport":
        return project.model.sim.get_model(resolve_transport_model_name(project))
    return project.model.gwf


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

    Skipped when the array is **already** stored one per layer (which every
    build after the first sees, and which ``npf_k`` starts out as) and when the
    model has a single layer -- there a whole-grid file already has exactly
    ``ncpl`` values, so pyEMU is happy and splitting it would be busywork.
    ``store_internal()`` first because FloPy refuses to make EXTERNAL data
    layered, which is the state a model loaded from a previous run is in.
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

    if recipe.family == "array":
        # Normalize MF6's wrapped array layout so pyEMU can read any cell count.
        for filename in files:
            _flatten_array_file(Path(project.template_workspace) / filename)

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
        kwargs["index_cols"] = list(recipe.index_cols)
        kwargs["use_cols"] = [recipe.use_col]
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

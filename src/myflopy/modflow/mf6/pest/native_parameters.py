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

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import re

import numpy as np
import pandas as pd


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


# DISV list files are written as ``layer cell <values...>`` -> index_cols=[0, 1].
_LIST_INDEX_COLS = [0, 1]

_RECIPES: dict[str, _Recipe] = {
    "k": _Recipe("k", "array", "{model}.npf_k.txt"),
    "k33": _Recipe("k33", "array", "{model}.npf_k33.txt"),
    "recharge": _Recipe("recharge", "list", "{model}.rch_stress_period_data_*.txt", use_col=2),
    "chd": _Recipe("chd", "list", "{model}.chd_stress_period_data_*.txt", use_col=2),
    "ghb.cond": _Recipe("ghb.cond", "list", "{model}.ghb_stress_period_data_*.txt", use_col=3),
    "ghb.bhead": _Recipe("ghb.bhead", "list", "{model}.ghb_stress_period_data_*.txt", use_col=2),
    "drn.cond": _Recipe("drn.cond", "list", "{model}.drn_stress_period_data_*.txt", use_col=3),
    "drn.elev": _Recipe(
        "drn.elev", "list", "{model}.drn_stress_period_data_*.txt", use_col=2, additive=True
    ),
    "wel": _Recipe("wel", "list", "{model}.wel_stress_period_data_*.txt", use_col=2),
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
        Variogram range for ``grid``/``pilotpoints`` geostatistics (Phase 2).
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
    temporal: float | None = None
    name: str | None = None
    capture: bool = False
    extra: dict = field(default_factory=dict)
    resolved_files: list[str] = field(default_factory=list)

    def __post_init__(self):
        self.recipe = resolve_target(self.target)
        if self.style is None:
            self.style = self.recipe.default_style
        self.style = str(self.style).strip().lower()
        if self.additive is None:
            self.additive = self.recipe.additive
        if self.name is None:
            self.name = self.recipe.canonical.replace(".", "")

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


def _resolve_files(template_workspace: Path, model_name: str, recipe: _Recipe) -> list[str]:
    """Return external input filenames (relative to the template) for a recipe."""

    stub = recipe.file_stub.format(model=model_name)
    workspace = Path(template_workspace)
    if "*" in stub:
        matches = sorted(p.name for p in workspace.glob(stub))
    else:
        candidate = workspace / stub
        if candidate.exists():
            matches = [candidate.name]
        else:
            # Multi-layer DISV/DIS arrays are externalized one file per layer,
            # e.g. ``<model>.npf_k_layer1.txt``. Match those precisely so that
            # resolving ``k`` never picks up ``k33`` (``npf_k_layer*`` requires
            # ``_layer`` immediately after ``npf_k``).
            base = stub[:-4] if stub.endswith(".txt") else stub
            matches = sorted(p.name for p in workspace.glob(f"{base}_layer*.txt"))
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
    files = _resolve_files(project.template_workspace, project.model.name, recipe)
    if spec.layers is not None and recipe.family == "array":
        wanted = {int(layer) for layer in spec.layers}
        selected = [
            name
            for name in files
            if (match := re.search(r"_layer(\d+)\.txt$", name))
            and (int(match.group(1)) - 1) in wanted
        ]
        if selected:
            files = selected
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
        kwargs["zone_array"] = spec.zones
    if recipe.family == "list":
        kwargs["index_cols"] = list(_LIST_INDEX_COLS)
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

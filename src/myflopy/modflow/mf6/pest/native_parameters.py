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
    """

    target: str
    style: str | None = None
    bounds: tuple[float, float] = (0.5, 2.0)
    physical: tuple[float, float] | None = None
    transform: str = "log"
    additive: bool | None = None
    zones: Any = None
    correlation: float | None = None
    temporal: float | None = None
    name: str | None = None
    extra: dict = field(default_factory=dict)

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


def _resolve_files(template_workspace: Path, model_name: str, recipe: _Recipe) -> list[str]:
    """Return external input filenames (relative to the template) for a recipe."""

    stub = recipe.file_stub.format(model=model_name)
    if "*" in stub:
        matches = sorted(p.name for p in Path(template_workspace).glob(stub))
    else:
        candidate = Path(template_workspace) / stub
        matches = [candidate.name] if candidate.exists() else []
    if not matches:
        raise FileNotFoundError(
            f"No external input file found for target {recipe.canonical!r} "
            f"(looked for {stub!r} in {template_workspace}). Did the model build "
            "with this package, and was set_all_data_external() applied?"
        )
    return matches


def add_native_parameter(project, spec: NativeParameterSpec):
    """Compile one :class:`NativeParameterSpec` into ``pf.add_parameters``.

    Returns the parameter dataframe pyEMU produced, and records it on the
    project for later inspection via :meth:`PestProject.settings`.
    """

    recipe = spec.recipe
    files = _resolve_files(project.template_workspace, project.model.name, recipe)

    if recipe.family == "array" and spec.style in _SPATIAL_STYLES:
        raise NotImplementedError(
            f"style={spec.style!r} on array target {recipe.canonical!r} needs a "
            "Voronoi cell spatial reference, which is Phase 2. Use style='constant' "
            "or style='zone' for now, or the legacy KPilotPointParameter path."
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

    # pyEMU accepts a single filename or a list; pass the list for multi-period
    # list packages so one parameter scales every stress period uniformly.
    filenames = files if len(files) > 1 else files[0]
    frame = project.pf.add_parameters(filenames, **kwargs)
    project._native_parameter_frames[spec.name] = frame
    return frame

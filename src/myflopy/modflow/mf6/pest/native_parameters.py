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

import numpy as np
import pandas as pd

from myflopy.modflow.mf6.package_registry import (
    SeveralPackagesError,
    get_package_explorer_spec,
    get_package_input_field_spec,
    model_packages_of_type,
    package_name_of,
    package_named,
    packages_of_type,
)
from myflopy.modflow.mf6.pest.model_lookup import resolve_transport_model_name
from myflopy.modflow.mf6.pest.zones import resolve_zone_array

# --- recipe registry -------------------------------------------------------
#
# Each recipe names the MF6 package TYPE and the FloPy dataset a calibration target
# scales, and pyEMU reads the files FloPy's OWN records name for that dataset on the
# build's staged model (:func:`package_target_files`) -- never a file name guessed
# from the type, which scaled one of two GHBs in silence and found nothing for a GHB
# not named ``ghb`` (ledger 517). ``use_col`` is the zero-based column
# holding the value pyEMU multiplies/offsets, on DISV. A DIS cellid is three columns,
# ``layer row col``, so on a structured grid every column after it moves one to the
# right -- :func:`list_columns` is the ONE rule that says where they are on the
# model's own grid (ledger 420).


@dataclass(frozen=True)
class _Recipe:
    """How one calibration target maps onto a FloPy dataset of an MF6 package TYPE.

    ``package`` is the package TYPE (``"npf"``, ``"riv"``, ...), matched by
    :func:`~myflopy.modflow.mf6.package_registry.packages_of_type` (RCHA counts as
    RCH), never a package NAME; ``variable`` is the FloPy dataset whose external
    files are scaled (``"k"``, ``"stress_period_data"``, ``"packagedata"``); and
    ``field`` is the record field a list target scales. A ``stress_period_data``
    recipe's ``use_col`` and whether it is an ELEVATION come from the package
    registry's ``record_fields`` and ``FieldSpec.kind`` (:attr:`elevation`), so
    neither is restated here.
    """

    canonical: str
    family: str  # "array" or "list"
    package: str  # MF6 package TYPE whose dataset is scaled
    variable: str  # FloPy dataset on that package
    field: str | None = None  # list family: the record field scaled
    use_col: int | None = None  # list family; derived from record_fields when None
    default_style: str = "constant"
    # Which model in the simulation owns the file. A PestProject hangs off the
    # FLOW model, so "transport" targets resolve the GWT sibling by name --
    # everything downstream (pyEMU's own apply_list_and_array_pars) is purely
    # filename-driven and needs no other change.
    model: str = "flow"
    # List family only: which columns hold the CELL IDENTITY. Most MF6 list
    # packages write ``layer cell <values...>``, so (0, 1) is the default -- but
    # UZF numbers its own cells first and writes ``iuzno layer cell ...``. pyEMU
    # geolocates a list parameter from these columns, and getting them wrong is
    # silent: every parameter lands on one cell's coordinates and only the
    # geostatistics are wrong (see the `uzf.vks` entry below).
    index_cols: tuple[int, ...] = (0, 1)

    def __post_init__(self):
        """Derive a ``stress_period_data`` recipe's value column from the registry."""

        if self.family == "list" and self.use_col is None:
            object.__setattr__(
                self, "use_col", self.index_cols[-1] + 1 + _record_offset(self.package, self.field)
            )

    @property
    def elevation(self) -> bool:
        """Whether the scaled field is an ELEVATION (``FieldSpec.kind``, ledger 297).

        An elevation is calibrated by an ADDITIVE offset in model units, never by a
        multiplier, whose prior would scale with the datum (ledger 518).
        """

        if self.field is None:
            return False
        spec = get_package_input_field_spec(self.package, self.field)
        return spec is not None and spec.kind == "elevation"


def _record_offset(package_type: str, field_name: str | None) -> int:
    """Where ``field_name`` sits among the package type's record fields after the cellid."""

    spec = get_package_explorer_spec(package_type)
    fields = tuple(spec.record_fields) if spec is not None else ()
    if field_name not in fields:
        raise ValueError(
            f"Package type {package_type!r} has no record field {field_name!r} in the "
            f"package registry (it lists {list(fields)}); give the recipe a use_col."
        )
    return fields.index(field_name)


# DISV list files are usually written as ``layer cell <values...>``; see
# `_Recipe.index_cols` for the exceptions.
_LIST_INDEX_COLS = (0, 1)

#: How many columns an MF6 cellid takes in a list file, by FloPy's grid type:
#: ``layer row col`` (DIS), ``layer cell`` (DISV), ``node`` (DISU).
_CELLID_WIDTH = {"structured": 3, "vertex": 2, "unstructured": 1}


def list_columns(recipe: _Recipe, grid_type: str, field: str | None = None) -> tuple[list[int], int]:
    """Where a list recipe's cell identity and a value sit on a grid of ``grid_type``.

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
    field : str or None, optional
        Another record field of the recipe's package type whose column is wanted --
        ``"rbot"`` beside ``riv.stage``, which moves the pair (ledger
        519). Default ``None``: the recipe's own value column.

    Returns
    -------
    tuple of (list of int, int)
        The zero-based ``index_cols`` -- the cellid's columns, pyEMU reading its
        LAST two as the cell's ``(i, j)`` -- and the value's column, both on that
        grid.
    """

    first, last = recipe.index_cols[0], recipe.index_cols[-1]
    if last - first != len(_LIST_INDEX_COLS) - 1:
        raise ValueError(
            f"Recipe {recipe.canonical!r} names index columns {recipe.index_cols}, "
            "which are not one DISV cellid (two adjacent columns)."
        )
    shift = _CELLID_WIDTH[grid_type] - len(_LIST_INDEX_COLS)
    width = _CELLID_WIDTH[grid_type]
    column = int(recipe.use_col) if field is None else last + 1 + _record_offset(recipe.package, field)
    return list(range(first, first + width)), column + shift


_RECIPES: dict[str, _Recipe] = {
    "k": _Recipe("k", "array", "npf", "k"),
    "k33": _Recipe("k33", "array", "npf", "k33"),
    "recharge": _Recipe("recharge", "list", "rch", "stress_period_data", "recharge"),
    "chd": _Recipe("chd", "list", "chd", "stress_period_data", "head"),
    "ghb.cond": _Recipe("ghb.cond", "list", "ghb", "stress_period_data", "cond"),
    "ghb.bhead": _Recipe("ghb.bhead", "list", "ghb", "stress_period_data", "bhead"),
    "drn.cond": _Recipe("drn.cond", "list", "drn", "stress_period_data", "cond"),
    "drn.elev": _Recipe("drn.elev", "list", "drn", "stress_period_data", "elev"),
    "riv.cond": _Recipe("riv.cond", "list", "riv", "stress_period_data", "cond"),
    # An additive offset on STAGE that moves RBOT with it (ledgers 290/291's pair
    # rule): pyEMU offsets the stage column, and `forward_run.carry_river_bottom`,
    # run after pyEMU's apply, sets rbot = rbot_org + (stage - stage_org), so the
    # declared depth survives every realization (ledger 519).
    "riv.stage": _Recipe("riv.stage", "list", "riv", "stress_period_data", "stage"),
    "wel": _Recipe("wel", "list", "wel", "stress_period_data", "q"),
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
        "uzf.vks", "list", "uzf", "packagedata", "vks", use_col=6, index_cols=(1, 2),
    ),
    # Transport. MST writes porosity as ONE whole-grid array (no LAYERED
    # keyword), so unlike npf_k there is no per-layer file to select from --
    # see `relayer_array_target` and `select_layers`.
    "porosity": _Recipe("porosity", "array", "mst", "porosity", model="transport"),
}

#: How many of a pilot-point set's nearest points a cell is interpolated from, by
#: default (ledger 526).
DEFAULT_PP_NEAREST = 8

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
    "riv": "riv.cond",
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

#: The four spatial styles :meth:`PestProject.parameterize` builds. Anything else is
#: REFUSED at the call -- pyEMU's own pilot-point spellings (``"pp"``,
#: ``"pilot_points"``, ...) reached pyEMU's structured-grid pilot points and died on
#: a DISV grid with ``'VertexGrid' object has no attribute 'delr'`` (ledger
#: 520).
STYLES = ("constant", "zone", "grid", "pilotpoints")

#: pyEMU's other spellings of pilot points, named by the refusal.
_PYEMU_PILOT_POINT_SPELLINGS = frozenset({
    "pp", "pilot_points", "pilot_point", "pilotpoint", "pilot-point", "pilot-points",
})

# pyEMU ``par_type`` values that need a cell spatial reference (Phase 2 work for
# Voronoi array parameters). These are accepted for list packages, where no
# spatial reference is required.
_SPATIAL_STYLES = {"grid", "pilotpoints"}

#: Default ``bounds`` of a MULTIPLIER: a factor of 2 either way of the input value.
MULTIPLIER_BOUNDS = (0.5, 2.0)

#: Default ``bounds`` of an ADDITIVE offset, in model length units: 2 either way of
#: the input elevation (ledger 518). Symmetric about 0, the value an
#: offset starts at.
OFFSET_BOUNDS = (-2.0, 2.0)

#: The most parameters one CORRELATED ``style="grid"`` block (one layer file of an
#: array target, or one list target) may hold: its prior is drawn from a DENSE
#: ``n x n`` covariance, which pyEMU builds, scales and decomposes in memory, at every
#: ``run_ies()``/``prior()``. Measured (``docs/reproducers/correlated_grid_prior_size.py``,
#: 100 realizations): peak 0.99 GB and 3.7 s at 5,000, 1.92 GB and 7.7 s at 8,000,
#: growing as ``n**2`` -- near 12 GB for one 20,323-cell Ten Trails layer. Pilot points
#: carry a correlated field at a fraction of the parameters (ledger 522).
CORRELATED_GRID_LIMIT = 5000


def resolve_target(target: str) -> _Recipe:
    """Return the recipe for a friendly target name, raising a helpful error."""

    key = str(target).strip().lower()
    key = _ALIASES.get(key, key)
    if key not in _RECIPES:
        known = ", ".join(sorted(set(_RECIPES) | set(_ALIASES)))
        raise KeyError(f"Unknown calibration target {target!r}. Known targets: {known}.")
    return _RECIPES[key]


def _package_slug(package) -> str | None:
    """The ``package=`` argument as a name suffix: ``"crisp"``, ``"riv_crisp"``, or ``None``."""

    if package is None:
        return None
    names = [package] if isinstance(package, str) else list(package)
    return "_".join(str(name).strip().lower() for name in names)


def checked_package_argument(package) -> str | tuple[str, ...] | None:
    """``package=`` as :class:`NativeParameterSpec` holds it: ``None``, a name, or a tuple of names.

    Parameters
    ----------
    package : str or sequence of str or None
        One package NAME (FloPy's ``package_name``, any case), several sharing one
        parameter set, or ``None``.

    Returns
    -------
    str or tuple of str or None
        The names stripped and lower-cased; a sequence as a tuple.

    Raises
    ------
    TypeError
        When it is not a non-empty string, a sequence of them, or ``None``.
    ValueError
        When a sequence is empty or names one package twice.
    """

    if package is None:
        return None
    if isinstance(package, str):
        if not package.strip():
            raise TypeError("package= takes a package's name, not an empty string.")
        return package.strip().lower()
    try:
        names = list(package)
    except TypeError:
        raise TypeError(
            f"package= takes a package's name or a sequence of names, not {package!r}."
        ) from None
    if not names or not all(isinstance(name, str) and name.strip() for name in names):
        raise TypeError(
            f"package= takes a package's name or a non-empty sequence of names, not {package!r}."
        )
    cleaned = tuple(name.strip().lower() for name in names)
    repeated = sorted({name for name in cleaned if cleaned.count(name) > 1})
    if repeated:
        raise ValueError(f"package={package!r} names {repeated[0]!r} twice.")
    return cleaned


def checked_layers(layers) -> tuple[int, ...]:
    """``layers=`` as :class:`NativeParameterSpec` holds it: a tuple of Python ``int``.

    A numpy integer -- ``layers=np.arange(3)`` -- is carried as ``int``: the build's
    metadata is JSON, which refuses ``numpy.int64``, so it crashed the build after
    every parameter was registered.

    Parameters
    ----------
    layers : int or iterable of int
        Zero-based layers; one layer may be given bare. A whole-number float
        (``2.0``) is that layer.

    Returns
    -------
    tuple of int
        The layers, in the order given.

    Raises
    ------
    ValueError
        When a value is not a whole number at or above 0, or the iterable is empty.
    """

    values = [layers] if np.ndim(layers) == 0 else list(np.asarray(layers, dtype=object).ravel())
    checked = []
    for value in values:
        try:
            number = float(value)
        except (TypeError, ValueError):
            number = float("nan")
        if isinstance(value, (bool, np.bool_, str, bytes)) or not np.isfinite(number) or number < 0 \
                or number != int(number):
            raise ValueError(
                f"layers={layers!r}: each layer is a zero-based whole number, and {value!r} "
                "is not."
            )
        checked.append(int(number))
    if not checked:
        raise ValueError(f"layers={layers!r} names no layer; pass None for every layer.")
    return tuple(checked)


@dataclass
class NativeParameterSpec:
    """A declarative, pyEMU-bound parameter request.

    Built by :meth:`PestProject.parameterize`, whose docstring is the reference for
    every field; constructing one directly runs the same checks, bar the two that
    need the project -- a ``name`` already declared (ledger 523) and a ``package``
    the model's packages refuse.

    Parameters
    ----------
    target : str
        Friendly target name, e.g. ``"k"``, ``"recharge"``, ``"ghb.cond"``,
        ``"riv.stage"``; an unknown one raises ``KeyError`` (:func:`resolve_target`).
    style : str or None, optional
        ``"constant"`` (one value), ``"zone"`` (one per ``zones`` label), ``"grid"``
        (one per list entry / cell) or ``"pilotpoints"``, any case. ``None`` (default)
        uses the recipe default, ``"constant"``. Anything else raises ``ValueError``
        naming the four -- pyEMU's own pilot-point spellings (``"pp"``,
        ``"pilot_points"``, ...) name ``"pilotpoints"`` -- and so does
        ``"pilotpoints"`` on a LIST target, which no engine here builds (ledger
        520).
    bounds : tuple of (float, float) or None, optional
        ``(lower, upper)`` of the adjustable multiplier (factors) or additive offset
        (model units). ``None`` (default) is :data:`MULTIPLIER_BOUNDS`
        ``(0.5, 2.0)`` for a multiplier and :data:`OFFSET_BOUNDS` ``(-2.0, 2.0)`` for
        an offset. ``ValueError`` when ``lower >= upper``, when a bound is not finite,
        when the bounds exclude the value the parameter starts at -- ``1.0`` for a
        multiplier (each pilot point too), ``0.0`` for an offset -- or when a
        log-transformed multiplier's lower bound is not above 0 (ledger
        518).
    physical : tuple of (float, float) or None, optional
        ``(ult_lbound, ult_ubound)`` -- hard limits on the *final* model value after
        all multipliers are applied, pilot points included: their step clamps a
        covered cell to the tightest limits of its sets and of the pyEMU multipliers
        on its file (ledger 528). Default ``None``: unclamped.
    transform : str, optional
        ``"log"`` (default) or ``"none"`` (pyEMU also takes ``"fixed"``). Forced to
        ``"none"`` for an additive offset (:attr:`resolved_transform`).
    additive : bool or None, optional
        Apply as an additive offset instead of a multiplier. ``None`` (default)
        follows the scaled field's ``FieldSpec.kind``: an ELEVATION (``chd``,
        ``ghb.bhead``, ``drn.elev``, ``riv.stage``) is an offset, everything else a
        multiplier. ``False`` on an elevation raises ``ValueError`` -- a multiplier on
        an elevation scales with the datum; ``True`` on array pilot points raises
        ``ValueError`` -- they multiply.
    zones : ZoneSpec or array_like of int or None, optional
        Zone labels for ``style="zone"`` (and to mask cells); with array-target
        ``style="pilotpoints"``, one pilot-point set per label ``> 0`` on each
        selected layer (ledger 527); see :meth:`PestProject.parameterize`.
        Default ``None``.
    layers : int or iterable of int or None, optional
        Zero-based layers an ARRAY target is restricted to, held as a tuple of
        ``int`` (:func:`checked_layers`, which refuses anything else). Default
        ``None``: every layer. Given for a LIST target it raises ``ValueError``, naming ``zones=``.
    package : str or tuple of str or None, optional
        Which package(s) of the target's TYPE are scaled, by NAME
        (:func:`checked_package_argument`). ``None`` (default): the model's only
        package of the type; several raise ``SeveralPackagesError`` at
        :meth:`PestProject.parameterize` or :meth:`PestProject.build`. A tuple shares
        one parameter set across the packages it names.
    pp_spacing : float or None, optional
        Pilot-point net spacing, in MODEL LENGTH UNITS, for array-target
        ``style="pilotpoints"`` (ledger 527); see :meth:`PestProject.parameterize`.
        Default ``None``. Array-target pilot points need exactly one of it and
        ``pp_points``, and a spacing that is not finite and above zero raises
        ``ValueError``; any other style refuses it.
    pp_points : geopandas.GeoDataFrame or array_like or None, optional
        Explicit pilot-point locations for array-target ``style="pilotpoints"``;
        see :meth:`PestProject.parameterize`. Default ``None``. Refused beside
        ``pp_spacing`` and on any other style.
    pp_nearest : int, optional
        How many of a set's nearest pilot points each cell is interpolated from
        (ledger 526): a whole number, 1 or more. Default
        :data:`DEFAULT_PP_NEAREST` (``8``). Any other value is refused on every
        style but array-target ``style="pilotpoints"``.
    correlation : float or None, optional
        Variogram range (model length units) for ``style="grid"`` or
        ``"pilotpoints"``. Default ``None``. Given with ``"constant"``/``"zone"`` it
        raises ``ValueError``: neither has anything to correlate.
    anisotropy : float, optional
        Anisotropy ratio of the variogram (``1.0``, the default, is isotropic). Not
        ``1.0`` without ``correlation`` raises ``ValueError``.
    bearing : float, optional
        Azimuth (degrees east of north) of the anisotropy major axis. Default
        ``0.0``. Not ``0.0`` without ``correlation``, or with ``anisotropy == 1.0``,
        raises ``ValueError``: it cannot act.
    nugget : float, optional
        Nugget of the variogram, relative to its contribution of ``1.0``; the prior
        draw scales the two to a sill of 1, so a nugget ``n`` divides every
        correlation by ``1 + n`` and leaves the bound-derived variance alone (ledger
        489). Default ``0.0``. Not ``0.0`` without ``correlation`` raises
        ``ValueError``.
    temporal : float or None, optional
        Reserved; nothing reads it. Default ``None``; any other value raises
        ``ValueError``.
    name : str or None, optional
        Parameter group / name base. Default ``None``: the canonical target with
        dots removed, plus ``_<package>`` when ``package`` is given
        (``"rivcond_crisp"``). Pilot-point sets are named ``<name>l<layer>`` /
        ``<name>l<layer>z<zone>``, lower-cased.
    capture : bool, optional
        Also record the resolved per-cell field as zero-weight observations, for
        ``IesResults.plot_field`` / ``IesResults.field``. Default ``False``.
    extra : dict, optional
        Keyword arguments passed last to pyEMU's ``PstFrom.add_parameters``. Default
        empty. Not empty for array-target pilot points, which never call it, raises
        ``ValueError``.
    resolved_files : list of str, optional
        Filled by :meth:`PestProject.build`: the template-relative files the
        parameter scales. Default empty.
    pilot_point_sets : list of dict, optional
        Filled by :meth:`PestProject.build` for array-target pilot points: one
        summary per pilot-point set (``name``, ``layer``, ``zone``, ``points``,
        ``cells``). Default empty.

    Attributes
    ----------
    recipe : _Recipe
        The target's recipe.
    resolved_targets : list of TargetFile
        Filled by :meth:`PestProject.build`, one per file: its package, layer and
        period, from FloPy's own records.
    resolved_zone_labels : list of int or None
        Filled by :meth:`PestProject.build` when ``zones`` is given: the distinct
        labels of the zone array it resolved (pyEMU's ``zone_array``; for array
        pilot points, the array whose labels ``> 0`` became sets).
    """

    target: str
    style: str | None = None
    bounds: tuple[float, float] | None = None
    physical: tuple[float, float] | None = None
    transform: str = "log"
    additive: bool | None = None
    zones: Any = None
    layers: tuple[int, ...] | None = None
    package: str | tuple[str, ...] | None = None
    pp_spacing: float | None = None
    pp_points: Any = None
    pp_nearest: int = DEFAULT_PP_NEAREST
    correlation: float | None = None
    anisotropy: float = 1.0
    bearing: float = 0.0
    nugget: float = 0.0
    temporal: float | None = None
    name: str | None = None
    capture: bool = False
    extra: dict = field(default_factory=dict)
    resolved_files: list[str] = field(default_factory=list)
    pilot_point_sets: list[dict] = field(default_factory=list)

    def __post_init__(self):
        """Resolve the recipe, fill the defaults, and refuse an argument that cannot act."""

        self.recipe = resolve_target(self.target)
        self.resolved_targets: list[TargetFile] = []
        self.resolved_zone_labels: list[int] | None = None
        recipe = self.recipe
        self.style = _checked_style(self.style if self.style is not None else recipe.default_style)
        self.package = checked_package_argument(self.package)
        if self.layers is not None:
            self.layers = checked_layers(self.layers)
        if self.correlation is not None:
            self.correlation = float(self.correlation)
        if self.additive is None:
            self.additive = recipe.elevation
        elif not self.additive and recipe.elevation:
            raise ValueError(
                f"additive=False on target {recipe.canonical!r}: its {recipe.field} is an "
                "ELEVATION, and a multiplier on an elevation scales with the datum (x1.1 "
                "moves a 400 ft head 40 ft and a 40 ft one 4). It is an additive offset in "
                f"model units: drop additive=False (bounds default to {OFFSET_BOUNDS}, "
                "ledger 518)."
            )
        self.additive = bool(self.additive)
        if self.style == "pilotpoints" and recipe.family == "list":
            raise ValueError(
                f"style='pilotpoints' on list target {recipe.canonical!r}: myflopy's pilot "
                "points interpolate onto an ARRAY, and pyEMU's own reject a list file "
                "(\"unrecognized 'typ' ... not 'pilotpoints'\") on every grid, so it cannot "
                "act. Use style='grid' (one parameter per record; correlation= correlates "
                "them), 'zone' or 'constant' (ledger 520)."
            )
        pilot_array = self.style == "pilotpoints" and recipe.family == "array"
        if pilot_array and self.additive:
            raise ValueError(
                f"additive=True on array-target pilot points ({recipe.canonical!r}): pilot "
                "points MULTIPLY the array in the forward run, so an offset cannot act. "
                "Drop additive=, or use style='constant'/'zone'/'grid'."
            )
        if pilot_array and self.extra:
            raise ValueError(
                f"{sorted(self.extra)} given to array-target pilot points "
                f"({recipe.canonical!r}): extra keywords go to pyEMU's "
                "PstFrom.add_parameters, which myflopy's pilot points never call, so they "
                "cannot act. Drop them."
            )
        if self.layers is not None and recipe.family == "list":
            raise ValueError(
                f"layers= on list target {recipe.canonical!r}: it scales every record of "
                "its package's files, whatever layer each names, so layers= cannot act. "
                "Restrict it with zones= -- a (nlay, ncpl) label array whose labels <= 0 "
                "mark records left alone -- or pick the package with package=."
            )
        if self.temporal is not None:
            raise ValueError(
                f"temporal={self.temporal!r}: nothing reads a temporal correlation yet, so "
                "it cannot act. Drop it."
            )
        self._check_pilot_point_arguments()
        self._refuse_an_idle_variogram()
        self.bounds = self._checked_bounds()
        if self.name is None:
            slug = _package_slug(self.package)
            base = recipe.canonical.replace(".", "")
            self.name = f"{base}_{slug}" if slug else base

    def _refuse_an_idle_variogram(self) -> None:
        """Refuse a variogram argument the spec's style or other arguments leave idle."""

        if self.correlation is not None:
            if self.style not in _SPATIAL_STYLES:
                raise ValueError(
                    f"correlation={self.correlation!r} with style={self.style!r}: a "
                    f"{self.style} parameter has nothing to correlate, so it cannot act. "
                    "Drop it, or use style='grid' or 'pilotpoints'."
                )
            if float(self.anisotropy) == 1.0 and float(self.bearing) != 0.0:
                raise ValueError(
                    f"bearing={self.bearing!r} with anisotropy=1.0: an isotropic variogram "
                    "has no axis to orient, so it cannot act. Give anisotropy=, or drop "
                    "bearing=."
                )
            return
        idle = {
            name: value
            for name, value, default in (
                ("anisotropy", self.anisotropy, 1.0),
                ("bearing", self.bearing, 0.0),
                ("nugget", self.nugget, 0.0),
            )
            if float(value) != default
        }
        if idle:
            raise ValueError(
                f"{', '.join(f'{key}={value!r}' for key, value in idle.items())} without "
                "correlation=: they shape a variogram, and none is built without one, so "
                "they cannot act. Give correlation= (style='grid' or 'pilotpoints'), or "
                "drop them."
            )

    def _checked_bounds(self) -> tuple[float, float]:
        """``bounds`` filled from the parameter's kind, refusing ones that exclude its start."""

        default = OFFSET_BOUNDS if self.additive else MULTIPLIER_BOUNDS
        given = default if self.bounds is None else self.bounds
        try:
            lower, upper = (float(value) for value in given)
        except (TypeError, ValueError):
            raise ValueError(
                f"bounds={self.bounds!r}: give (lower, upper), two numbers."
            ) from None
        what = "offset" if self.additive else "multiplier"
        if not (np.isfinite(lower) and np.isfinite(upper)) or lower >= upper:
            raise ValueError(
                f"bounds={self.bounds!r} for {self.recipe.canonical!r}: give two finite "
                "numbers, lower < upper."
            )
        start = 0.0 if self.additive else 1.0
        if not lower <= start <= upper:
            unit = "model units" if self.additive else "a factor"
            raise ValueError(
                f"bounds={self.bounds!r} for {self.recipe.canonical!r} exclude {start:g}, the "
                f"value every {what} starts at ({unit}; pilot points start each point at "
                "1.0), so PEST would begin outside its own bounds -- PESTPP-IES warns and "
                "holds the draw to them, so the input value is never tried. Give bounds "
                f"around {start:g}, e.g. {default}; to move the input value itself, change "
                "the model (ledger 518)."
            )
        if self.resolved_transform == "log" and lower <= 0.0:
            raise ValueError(
                f"bounds={self.bounds!r} for {self.recipe.canonical!r}: a log-transformed "
                f"multiplier's lower bound must be above 0, and {lower:g} is not. Give one "
                "above 0, or transform='none'."
            )
        return (lower, upper)

    def _check_pilot_point_arguments(self) -> None:
        """Refuse pilot-point placement arguments that cannot act, naming the spelling.

        ``pp_spacing`` replaced ``pp_space`` (median cell widths) with no alias
        (ledger 527): a net on a strongly refined mesh measured in its
        median cell laid 4,501 points per layer on Ten Trails.
        """

        pilot = self.style == "pilotpoints" and self.recipe.family == "array"
        if not pilot:
            given = [
                name for name, value in (("pp_spacing", self.pp_spacing),
                                         ("pp_points", self.pp_points))
                if value is not None
            ]
            if self.pp_nearest != DEFAULT_PP_NEAREST:
                given.append("pp_nearest")
            if given:
                raise ValueError(
                    f"{', '.join(given)} place and interpolate pilot points, and "
                    f"{self.recipe.canonical!r} with style={self.style!r} has none; they act "
                    "only with style='pilotpoints' on an array target (k, k33, porosity)."
                )
            return
        if (self.pp_spacing is None) == (self.pp_points is None):
            raise ValueError(
                f"style='pilotpoints' on {self.recipe.canonical!r} places its points by "
                "pp_spacing= (a regular net, its spacing in MODEL LENGTH UNITS -- e.g. "
                "pp_spacing=1500.0 for a net 1,500 ft apart on a model in feet) or by "
                "pp_points= (explicit locations): give exactly one."
            )
        if self.pp_spacing is not None:
            if isinstance(self.pp_spacing, bool) or not np.isfinite(float(self.pp_spacing)) \
                    or float(self.pp_spacing) <= 0:
                raise ValueError(
                    f"pp_spacing={self.pp_spacing!r} must be a finite length above zero, in "
                    "model length units."
                )
            self.pp_spacing = float(self.pp_spacing)
        if isinstance(self.pp_nearest, bool) or int(self.pp_nearest) != self.pp_nearest \
                or int(self.pp_nearest) < 1:
            raise ValueError(
                f"pp_nearest={self.pp_nearest!r} must be a whole number of pilot points, 1 "
                "or more."
            )
        self.pp_nearest = int(self.pp_nearest)

    @property
    def resolved_transform(self) -> str:
        """Transform actually used (additive parameters cannot be log)."""

        return "none" if self.additive else self.transform


def _checked_style(style) -> str:
    """``style`` lower-cased, refusing anything but the four :data:`STYLES`."""

    value = str(style).strip().lower()
    if value in STYLES:
        return value
    hint = (
        " -- it is pyEMU's spelling of pilot points, whose own engine cannot place them "
        "on a DISV grid; myflopy's is style='pilotpoints'"
        if value in _PYEMU_PILOT_POINT_SPELLINGS else ""
    )
    raise ValueError(
        f"style={style!r} is not one of {list(STYLES)}{hint} (ledger 520)."
    )


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


def pyemu_spatial_reference(modelgrid):
    """The spatial reference pyEMU locates a cell with, on ``modelgrid``.

    The FloPy grid itself on DIS and DISV, which pyEMU indexes at ``[row, col]``
    and ``[cell]``. On DISU a dict ``{node: (x, y)}`` of the zero-based nodes'
    centres, pyEMU's spelling for an unstructured grid: handed the
    ``UnstructuredGrid`` it indexed the one-dimensional centre array at ``[i, j]``
    and every array parameter died in ``write_array_tpl`` with ``IndexError: too
    many indices for array`` (ledger 452). An array file holds one node per
    line there (:func:`array_row_length`), and a list file one ``node`` column
    (:func:`list_columns`), and pyEMU reads both through the dict.

    Parameters
    ----------
    modelgrid : flopy.discretization.Grid
        The staged model's FloPy grid; ``grid_type`` and, on DISU,
        ``xcellcenters``/``ycellcenters`` are read.

    Returns
    -------
    flopy.discretization.Grid or dict
        ``modelgrid`` unchanged, or the node-centre dict on DISU.
    """

    if getattr(modelgrid, "grid_type", None) != "unstructured":
        return modelgrid
    x = np.asarray(modelgrid.xcellcenters, dtype=float).reshape(-1)
    y = np.asarray(modelgrid.ycellcenters, dtype=float).reshape(-1)
    return {node: (float(x[node]), float(y[node])) for node in range(x.size)}


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
    model -- a recipe's package array (:func:`relayer_array_target`), the cell
    centres pilot points are placed on (their base is the template's own array
    file, ledger 525), the
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



# --- which files a target scales: FloPy's own records, by package TYPE ----------


@dataclass(frozen=True)
class TargetFile:
    """One external file a parameter scales, as FloPy's records of the staged model name it.

    Attributes
    ----------
    filename : str
        The file, relative to the simulation's workspace (the template's).
    package : str
        The FloPy ``package_name`` of the package that reads it.
    layer : int or None
        The zero-based layer of a per-layer array file; ``None`` for a whole-grid
        array or a list file.
    period : int or None
        The zero-based stress period of a transient dataset's block; ``None`` for
        data that is not per period (griddata, ``packagedata``).
    """

    filename: str
    package: str
    layer: int | None = None
    period: int | None = None


def select_packages(model, recipe: _Recipe, package=None) -> list:
    """The packages of ``recipe``'s TYPE a parameter scales, chosen by the ledger-352 rule.

    Found by TYPE (:func:`~myflopy.modflow.mf6.package_registry.packages_of_type`,
    RCHA counting as RCH), never by FloPy's name-first ``get_package``; a package is
    never picked in silence.

    Parameters
    ----------
    model : flopy.mf6.MFModel or SimulationBase
        A FloPy model (the build's staged one) or a myflopy model
        (:func:`~myflopy.modflow.mf6.package_registry.model_packages_of_type`).
    recipe : _Recipe
        The target's recipe; its ``package`` is the TYPE.
    package : str or tuple of str or None, optional
        ``None`` (default): the model's only package of the type. A name: that
        package. A tuple of names: each of them, in that order.

    Returns
    -------
    list
        The FloPy packages, in the order named (model order for ``None``).

    Raises
    ------
    ValueError
        When the model holds no package of the type, or a name names none of them
        (listing the ones it holds).
    SeveralPackagesError
        When ``package`` is ``None`` and the model holds several of the type, naming
        each and the spelling that picks one or shares one set.
    """

    if callable(getattr(model, "_packages_of_type", None)):
        everyone = model_packages_of_type(model, recipe.package)
    else:
        everyone = packages_of_type(model, recipe.package)
    kind = recipe.package.upper()
    model_name = getattr(model, "name", "?")
    if not everyone:
        raise ValueError(
            f"Target {recipe.canonical!r} scales the {kind} package's {recipe.variable}, and "
            f"model {model_name!r} in the run's files has no {kind} package. Parameterize a "
            "target the model has."
        )
    if package is None:
        if len(everyone) == 1:
            return list(everyone)
        names = [package_name_of(found) for found in everyone]
        raise SeveralPackagesError(
            f"Model {model_name!r} carries {len(everyone)} {kind} packages "
            f"({', '.join(repr(name) for name in names)}), and target {recipe.canonical!r} "
            "scales one package's records unless told otherwise -- a package is never "
            "picked in silence (ledger 352). Name it: "
            + " or ".join(
                f"cal.parameterize({recipe.canonical!r}, package={name!r})" for name in names
            )
            + f"; or share ONE parameter set: package={tuple(names)!r}."
        )
    names = (package,) if isinstance(package, str) else tuple(package)
    return [package_named(model, recipe.package, name, everyone=everyone) for name in names]


def _storage_files(storage) -> list[tuple[int | None, str | None]]:
    """``(layer, file name)`` of one FloPy data storage: one per layer of a layered array."""

    layer_storage = getattr(storage, "layer_storage", None)
    if layer_storage is not None:
        elements = list(layer_storage.elements())
        if getattr(storage, "layered", False):
            return [(layer, element.fname) for layer, element in enumerate(elements)]
        return [(None, elements[0].fname if elements else None)]
    return [(None, getattr(storage, "fname", None))]


def _dataset_files(data) -> list[tuple[int | None, int | None, str | None]]:
    """``(period, layer, file name)`` for every block of a FloPy dataset.

    Read off the data's STORAGE, never ``get_record()``, which loads the data and, for
    a list of a model kept in a sub-directory, opened it at the model's path where
    MF6 (and FloPy's writer) put it at the simulation's.
    """

    keys = getattr(data, "get_active_key_list", None)
    if callable(keys):
        return [
            (int(key), layer, filename)
            for key, storage in keys()
            for layer, filename in _storage_files(storage)
        ]
    return [(None, layer, filename) for layer, filename in _storage_files(data._get_storage_obj())]


def package_target_files(model, recipe: _Recipe, packages, workspace) -> list[TargetFile]:
    """The external files FloPy's records say ``packages`` read for ``recipe``'s dataset.

    Read off the packages' own data storage -- each period's of a transient dataset
    (``get_active_key_list()``), each layer's of a layered array, one otherwise --
    after ``set_all_data_external()``, which names every file in memory, relative to
    the simulation's directory as MF6 reads it. Never a file name built from the type: FloPy names an
    external file after the package's FILE (``m.crisp_stress_period_data_1.txt``), so
    a glob on the type scaled one of two GHBs and found nothing for a GHB named
    ``edge`` (ledger 517); and a stale whole-grid file a relayering leaves
    behind is never named.

    Parameters
    ----------
    model : flopy.mf6.MFModel
        The model the packages belong to (named by the refusals).
    recipe : _Recipe
        The target's recipe; ``variable`` names the dataset.
    packages : list
        The FloPy packages (:func:`select_packages`).
    workspace : str or pathlib.Path
        The simulation's workspace, where each file must exist.

    Returns
    -------
    list of TargetFile
        Package by package as given, each in period then layer order.

    Raises
    ------
    ValueError
        When a package holds no such dataset (an array-form RCH asked for its list
        ``stress_period_data``), or its data is not in external files.
    FileNotFoundError
        When a file a record names is not in ``workspace``.
    """

    workspace = Path(workspace)
    found: list[TargetFile] = []
    for package in packages:
        name = package_name_of(package)
        data = getattr(package, recipe.variable, None)
        if data is None or not hasattr(data, "_get_storage_obj"):
            raise ValueError(
                f"Target {recipe.canonical!r} scales the {recipe.package.upper()} "
                f"package's {recipe.variable}, and package {name!r} of model "
                f"{getattr(model, 'name', '?')!r} ({type(package).__name__}) has none"
                + (
                    " -- it is the ARRAY form (READASARRAYS), which this target does not "
                    "scale"
                    if str(getattr(package, "package_type", "")).lower() in {"rcha", "evta"}
                    else ""
                )
                + "."
            )
        for period, layer, filename in _dataset_files(data):
            if not filename:
                raise ValueError(
                    f"Package {name!r}'s {recipe.variable} (target {recipe.canonical!r}) "
                    f"is not in an external file{'' if period is None else f' in period {period}'}"
                    "; the build moves every array and list into one "
                    "(set_all_data_external()) before it reads them."
                )
            if not (workspace / filename).is_file():
                raise FileNotFoundError(
                    f"Package {name!r} reads {filename!r} for target {recipe.canonical!r}, "
                    f"and {workspace} holds no such file."
                )
            found.append(TargetFile(str(filename), name, layer, period))
    if not found:
        raise ValueError(
            f"Target {recipe.canonical!r}: package(s) "
            f"{[package_name_of(package) for package in packages]} hold no "
            f"{recipe.variable} records to scale."
        )
    return found


def select_layers(targets: list[TargetFile], layers, recipe: _Recipe) -> list[TargetFile]:
    """Narrow per-layer array files to the requested zero-based layers.

    Refuses rather than falling through when nothing matches. MF6 externalizes
    some griddata as ONE whole-grid array (MST porosity) and some one file per
    layer (NPF K); on the former there is no per-layer file to pick, and the
    previous ``if selected:`` fallthrough silently parameterized every layer
    while the caller believed ``layers=`` had restricted it.

    Parameters
    ----------
    targets : list of TargetFile
        An array target's files (:func:`package_target_files`).
    layers : iterable of int
        The zero-based layers wanted.
    recipe : _Recipe
        The target's recipe, named by the refusals.

    Returns
    -------
    list of TargetFile
        The files of those layers, in the order given.

    Raises
    ------
    ValueError
        When the target is one whole-grid array, or no requested layer has a file
        (naming the layers that do).
    """

    wanted = {int(layer) for layer in layers}
    selected = [target for target in targets if target.layer is not None and target.layer in wanted]
    if selected:
        return selected
    layered = [target for target in targets if target.layer is not None]
    if not layered:
        raise ValueError(
            f"layers={sorted(wanted)} cannot be applied to target "
            f"{recipe.canonical!r}: MODFLOW writes it as a single whole-grid "
            f"array ({targets[0].filename}), not one file per layer, so there is no "
            "per-layer file to select. Drop `layers=` (the parameter covers "
            "the whole grid), or use `zones=` to restrict it spatially."
        )
    raise ValueError(
        f"layers={sorted(wanted)} matched no external input file for target "
        f"{recipe.canonical!r}. Available layers: {sorted({t.layer for t in layered})}."
    )


def target_files(project, spec: NativeParameterSpec) -> list[TargetFile]:
    """The files ``spec`` scales in the build's staged model, recorded on the spec.

    The one rule every recipe goes through (ledger 517): the packages of
    the target's TYPE (:func:`select_packages`, honouring ``package=``), their files
    from FloPy's own records (:func:`package_target_files`), narrowed to ``layers=``
    (:func:`select_layers`). Sets ``spec.resolved_targets`` and
    ``spec.resolved_files``.

    Parameters
    ----------
    project : PestProject
        The calibration being built; its staged simulation is read.
    spec : NativeParameterSpec
        The parameter request.

    Returns
    -------
    list of TargetFile
        The files, as ``spec.resolved_targets`` holds them.
    """

    recipe = spec.recipe
    model = flopy_model_for(project, recipe)
    packages = select_packages(model, recipe, spec.package)
    workspace = Path(model.simulation_data.mfpath.get_sim_path())
    targets = package_target_files(model, recipe, packages, workspace)
    if spec.layers is not None:
        targets = select_layers(targets, spec.layers, recipe)
    spec.resolved_targets = list(targets)
    spec.resolved_files = [target.filename for target in targets]
    return targets



#: What MF6 uses for an NPF array the package leaves out, said when a target names one.
_UNDECLARED_MEANS = {
    ("npf", "k33"): "MF6 then sets K33 equal to K in every cell",
    ("npf", "k22"): "MF6 then sets K22 equal to K in every cell",
}


def refuse_an_undeclared_array(project, recipe: _Recipe, package=None) -> None:
    """Refuse an array target whose package, in the run's files, holds no such array.

    ``cal.parameterize("k33")`` on a model whose NPF declares no ``k33`` died inside
    :func:`relayer_array_target` with ``cannot reshape array of size 1 into shape
    (2,newaxis)`` on two layers, and with "No external input file found" on one --
    neither naming the target nor the cause: there is no ``k33`` file to multiply,
    because MF6 takes K33 from K. It reads the build's STAGED model
    (:func:`flopy_model_for`), so it answers about the run's files.

    Parameters
    ----------
    project : PestProject
        The calibration being built; its staged simulation is read.
    recipe : _Recipe
        The target's recipe. Anything but an ``"array"`` recipe is let through unread.
    package : str or tuple of str or None, optional
        The spec's ``package=`` (:func:`select_packages`). Default ``None``: the
        model's only package of the type.

    Raises
    ------
    ValueError
        When the model has no such package, or the package does not declare the
        array, naming the target, the package and what MF6 uses instead.
    SeveralPackagesError
        When the model holds several packages of the type and ``package`` is ``None``.
    """

    if recipe.family != "array":
        return
    model = flopy_model_for(project, recipe)
    # By TYPE, never FloPy's name-first `get_package` (ledger 352); with no package
    # of the type this raises naming the target.
    for found in select_packages(model, recipe, package):
        array = getattr(found, recipe.variable, None)
        if array is None or not array.has_data():
            break
    else:
        return
    means = _UNDECLARED_MEANS.get((recipe.package, recipe.variable))
    raise ValueError(
        f"Target {recipe.canonical!r} multiplies the {recipe.package.upper()} package's "
        f"{recipe.variable.upper()}, and the {recipe.package.upper()} of model "
        f"{model.name!r} in the run's files does not declare it"
        + (f" ({means})" if means else "")
        + ", so there is no array to multiply. Declare "
        f"`{recipe.variable}=` on the package and run the model again, or parameterize "
        "a target the package declares."
    )


def relayer_array_target(project, recipe: _Recipe, package=None) -> bool:
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
    files is loaded in. The values are read WITH their ``FACTOR``, so the template's
    model holds the numbers the run's did (a whole-grid ``FACTOR 2.0`` was dropped,
    and the template's K was half the run's).

    Parameters
    ----------
    project : PestProject
        The calibration being built; its staged simulation is changed.
    recipe : _Recipe
        The target's recipe; anything but an ``"array"`` recipe is left alone.
    package : str or tuple of str or None, optional
        The spec's ``package=`` (:func:`select_packages`). Default ``None``.

    Returns
    -------
    bool
        Whether any package's array was re-stored.
    """

    if recipe.family != "array":
        return False
    model = flopy_model_for(project, recipe)
    applied = False
    for found in select_packages(model, recipe, package):
        applied = _relayer_one(model, getattr(found, recipe.variable, None)) or applied
    return applied


def _relayer_one(model, array) -> bool:
    """Re-store one package's griddata one array per layer; see :func:`relayer_array_target`."""

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
    # WITH the multiplier: a whole-grid array read `INTERNAL FACTOR 2.0` holds its raw
    # values, and re-storing those dropped the factor -- the template's K was half the
    # run's. The re-stored layers carry FACTOR 1.0.
    values = np.asarray(array.get_data(apply_mult=True), dtype=float)
    if values.ndim < 2 or values.shape[0] != nlay:
        values = values.reshape(nlay, -1)
    array.store_internal()
    array.make_layered()
    array.set_data([values[layer] for layer in range(nlay)])
    return True



def _refuse_a_singular_grid_prior(spec: NativeParameterSpec, modelgrid) -> None:
    """Refuse a correlated ``grid`` parameter on DISU nodes that share a centre.

    On DISU one array file holds every node, so a layered DISU puts one parameter
    per node in ONE correlated set, and the nodes stacked in a column share their
    centre. Their correlation is 1, so the prior covariance is singular and
    PESTPP-IES's draw (pyEMU's ``geostatistical_draws``) died with ``LinAlgError:
    Matrix is not positive definite`` (ledger 452). A ``nugget`` now reaches the
    draw (ledger 489) and makes it drawable, but it would correlate the nodes
    of one column at ``1 / (1 + nugget)`` -- a vertical structure the map-plane
    variogram never described (ledger 119's one covariance per layer) -- so the set
    is refused whatever the nugget. DIS and DISV write one file per layer, each its
    own set, whose centres are distinct.
    """

    if spec.correlation is None or getattr(modelgrid, "grid_type", None) != "unstructured":
        return
    centres = np.column_stack([
        np.asarray(modelgrid.xcellcenters, dtype=float).reshape(-1),
        np.asarray(modelgrid.ycellcenters, dtype=float).reshape(-1),
    ])
    _unique, group, counts = np.unique(centres, axis=0, return_inverse=True, return_counts=True)
    group = np.asarray(group).reshape(-1)
    stacked = counts[group] > 1
    if not stacked.any():
        return
    node = int(np.flatnonzero(stacked)[0])
    others = [int(other) for other in np.flatnonzero(group == group[node]) if other != node]
    raise ValueError(
        f"style='grid' with correlation= on target {spec.recipe.canonical!r} draws one "
        f"correlated multiplier per node of this DISU model, and {int(stacked.sum())} of its "
        f"{len(centres)} nodes share a centre with another -- node {node} at "
        f"({centres[node, 0]:g}, {centres[node, 1]:g}) with node(s) {others[:5]}, as a "
        "layered DISU's stacked nodes do -- so their correlation is 1 and the prior "
        "cannot be drawn. Use style='zone' (a zone per layer of nodes), "
        "style='pilotpoints', style='constant', or style='grid' without correlation=."
    )


def _refuse_a_large_correlated_grid(spec: NativeParameterSpec, frame: pd.DataFrame) -> None:
    """Refuse a correlated ``style="grid"`` block above :data:`CORRELATED_GRID_LIMIT` parameters.

    One correlated block's prior is drawn from a DENSE covariance of all its
    parameters (pyEMU's ``geostatistical_draws``, at every :meth:`PestProject.run_ies`
    and :meth:`PestProject.prior`): 3.3 GB of float64 for one 20,323-cell layer
    before any copy, and the draw raised ``MemoryError`` under a 12 GB address-space
    limit -- with nothing said at build (ledger 522).
    """

    if spec.style != "grid" or spec.correlation is None or frame is None or frame.empty:
        return
    sizes = frame.groupby("pargp").size()
    largest = sizes.idxmax()
    count = int(sizes.max())
    if count <= CORRELATED_GRID_LIMIT:
        return
    gigabytes = count * count * 8 / 1e9
    raise ValueError(
        f"style='grid' with correlation= on target {spec.recipe.canonical!r} puts {count:,} "
        f"parameters in one correlated block ({largest!r}), above the {CORRELATED_GRID_LIMIT:,} "
        f"this build allows: its prior is drawn from a dense {count:,} x {count:,} covariance, "
        f"{gigabytes:.1f} GB before pyEMU's copies, at every run_ies()/prior(). Use "
        "style='pilotpoints' (a correlated field from a few hundred points), restrict the "
        "block with zones=, or drop correlation= (ledger 522)."
    )


def _refuse_a_river_stage_below_its_bottom(project, spec: NativeParameterSpec, grid_type: str) -> None:
    """Refuse ``riv.stage`` on records whose STAGE is below their RBOT.

    MF6 refuses such a record at any icelltype, and the offset moves the pair, so
    every realization would carry the same refusal (ledger 519).
    """

    _index, stage_col = list_columns(spec.recipe, grid_type)
    _index, rbot_col = list_columns(spec.recipe, grid_type, field="rbot")
    for target in spec.resolved_targets:
        rows = _read_list_file(Path(project.template_workspace) / target.filename)
        stage = pd.to_numeric(rows[stage_col])
        rbot = pd.to_numeric(rows[rbot_col])
        below = np.flatnonzero((stage < rbot).to_numpy())
        if below.size:
            row = int(below[0])
            raise ValueError(
                f"riv.stage on package {target.package!r}: {below.size} record(s) of "
                f"{target.filename} hold a stage below their river bottom -- row {row + 1}: "
                f"stage {stage.iloc[row]:g} < rbot {rbot.iloc[row]:g}. MF6 refuses such a "
                "record at any icelltype, and the offset keeps each record's depth, so every "
                "realization would fail. Correct the model's river first."
            )


def _read_list_file(path: Path) -> pd.DataFrame:
    """An MF6 external list file as text columns (a boundname stays text)."""

    return pd.read_csv(path, sep=r"\s+", header=None, dtype=str, keep_default_na=False,
                       comment="#")


def _carry_river_bottom(project, spec: NativeParameterSpec, grid_type: str) -> None:
    """Make every forward run move each scaled river record's RBOT with its STAGE.

    Writes, beside each file, the stage and rbot the template starts from, and adds
    ``forward_run.carry_river_bottom`` as a pre-model command, which pyEMU runs AFTER
    its own ``apply_list_and_array_pars`` (that one is first in the list). Once per
    file, however many ``riv.stage`` parameters scale it.
    """

    _index, stage_col = list_columns(spec.recipe, grid_type)
    _index, rbot_col = list_columns(spec.recipe, grid_type, field="rbot")
    # Emptied by each `PestProject.build`: the set belongs to ONE template's forward run.
    carried = project.__dict__.setdefault("_carried_river_files", set())
    helper = str(Path(__file__).with_name("forward_run.py"))
    for target in spec.resolved_targets:
        if target.filename in carried:
            continue
        rows = _read_list_file(Path(project.template_workspace) / target.filename)
        pair_csv = f"{target.filename}.riv_pair.csv"
        pd.DataFrame({
            "stage": pd.to_numeric(rows[stage_col]).to_numpy(),
            "rbot": pd.to_numeric(rows[rbot_col]).to_numpy(),
        }).to_csv(Path(project.template_workspace) / pair_csv, index=False)
        project.pf.add_py_function(
            helper,
            f"carry_river_bottom(model_file='{target.filename}', pair_csv='{pair_csv}', "
            f"stage_col={stage_col}, rbot_col={rbot_col})",
            is_pre_cmd=True,
        )
        carried.add(target.filename)


def add_native_parameter(project, spec: NativeParameterSpec):
    """Compile one :class:`NativeParameterSpec` into ``pf.add_parameters``.

    Returns the parameter dataframe pyEMU produced, and records it on the
    project for later inspection via :meth:`PestProject.settings` and the build's
    metadata.
    """

    recipe = spec.recipe
    targets = target_files(project, spec)
    files = [target.filename for target in targets]

    grid = flopy_model_for(project, recipe).modelgrid
    if recipe.family == "array":
        # Normalize MF6's wrapped array layout so pyEMU can read any cell count,
        # in the shape it reads on this grid (one grid row per line on DIS).
        for filename in files:
            _flatten_array_file(
                Path(project.template_workspace) / filename, array_row_length(grid)
            )

    if recipe.family == "array" and spec.style == "grid":
        _refuse_a_singular_grid_prior(spec, grid)
    if recipe.canonical == "riv.stage":
        _refuse_a_river_stage_below_its_bottom(project, spec, grid.grid_type)

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
        spec.resolved_zone_labels = sorted(
            int(label) for label in np.unique(np.asarray(kwargs["zone_array"]))
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
    # parameter set *per layer file*, so each layer's K field varies on its own,
    # named by its ZERO-BASED layer as pilot points are (ledger 521).
    if recipe.family == "array" and spec.style == "grid" and len(files) > 1:
        frames = []
        for target in targets:
            suffix = target.layer if target.layer is not None else files.index(target.filename)
            layer_kwargs = dict(kwargs)
            layer_kwargs["par_name_base"] = f"{spec.name}l{suffix}"
            layer_kwargs["pargp"] = f"{spec.name}l{suffix}"
            frames.append(project.pf.add_parameters(target.filename, **layer_kwargs))
        frame = pd.concat(frames, ignore_index=True)
    else:
        filenames = files if len(files) > 1 else files[0]
        frame = project.pf.add_parameters(filenames, **kwargs)
    _refuse_a_large_correlated_grid(spec, frame)
    if recipe.canonical == "riv.stage":
        _carry_river_bottom(project, spec, grid.grid_type)
    project._native_parameter_frames[spec.name] = frame
    return frame

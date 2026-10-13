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

import inspect
import re
import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from myflopy._logging import get_logger
from myflopy.modflow.mf6.boundary_support import bottom_clearance
from myflopy.modflow.mf6.hydraulics import declared_time_unit
from myflopy.modflow.mf6.package_registry import (
    SeveralPackagesError,
    feature_boundnames,
    get_package_explorer_spec,
    get_package_input_field_spec,
    model_packages_of_type,
    package_name_of,
    package_named,
    packages_of_type,
)
from myflopy.modflow.mf6.period_calendar import (
    PeriodCalendar,
    parse_start,
    steady_periods,
    tdis_start_date_time,
)
from myflopy.modflow.mf6.period_end import period_lengths
from myflopy.modflow.mf6.pest.model_lookup import resolve_transport_model_name
from myflopy.modflow.mf6.pest.naming import pest_group_name
from myflopy.modflow.mf6.pest.zones import resolve_zone_array

logger = get_logger(__name__)

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
    # The dataset the package's ARRAY form (READASARRAYS) holds the same quantity in,
    # for a target whose package comes in both forms: ``recharge`` is a list RCH's
    # record field and an array RCH's per-period array. ``None``: one form only. The
    # form is the staged package's (:func:`recipe_for_packages`), never guessed.
    array_dataset: str | None = None
    # Array family only: the dataset is per STRESS PERIOD (an array RCH's
    # ``recharge``), one file per period block, rather than griddata.
    periodic: bool = False
    # A FEATURE target (LAK, SFR; ledger 558): the word for one feature
    # ("lake", "reach") -- what `features=` names, by a ZERO-based number or a
    # boundname -- else None. ``key_cols`` are the columns BEFORE the cellid that name a
    # row uniquely in its file: LAK's (lakeno, iconn), SFR's (ifno). pyEMU keys a row by
    # them, never by its cell, because a lake connects to one cell twice and two
    # reaches can share a cell: a key that is unique per row is one record per row.
    feature: str | None = None
    key_cols: tuple[int, ...] = ()
    # The styles the target takes; () is every style but "feature".
    styles: tuple[str, ...] = ()

    def __post_init__(self):
        """Derive a ``stress_period_data`` recipe's value column from the registry."""

        if self.family == "list" and self.use_col is None:
            object.__setattr__(
                self, "use_col", self.index_cols[-1] + 1 + _record_offset(self.package, self.field)
            )

    @property
    def has_periods(self) -> bool:
        """Whether the target's files are per STRESS PERIOD, so ``periods=`` can group them.

        A list package's ``stress_period_data`` and an array RCH's ``recharge`` are;
        griddata (``k``, ``sto.sy``) and ``packagedata`` (``uzf.vks``) are not.
        """

        return self.variable == "stress_period_data" or self.periodic \
            or self.array_dataset is not None

    @property
    def written_fields(self) -> tuple[str, ...]:
        """The record fields a forward run WRITES for this list target: its field, and
        ``rbot`` beside ``riv.stage``, which the forward run carries with the stage
        (ledger 519). Empty for an array target. What a write-back substitutes
        (ledger 560).
        """

        if self.family != "list" or self.field is None:
            return ()
        return (self.field, "rbot") if self.canonical == "riv.stage" else (self.field,)

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
    # Both forms of RCH (ledger 547): a list RCH's ``recharge`` record
    # field, or an array RCH's (READASARRAYS) per-period ``recharge`` array -- one
    # file per period BLOCK, which FloPy numbers by the block's ORDINAL, never by its
    # period, so a file's period is read from FloPy's records. ``irch`` is never
    # touched. Which form is the staged package's (:func:`recipe_for_packages`).
    "recharge": _Recipe("recharge", "list", "rch", "stress_period_data", "recharge",
                        array_dataset="recharge"),
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
    # Lake-bed leakance, per LAKE (ledger 558). CONNECTIONDATA rows are
    # ``lakeno iconn <cellid> claktype bedleak belev telev connlen connwidth``, so
    # bedleak is one column after the cellid. A row is keyed by (lakeno, iconn): a
    # lake connects to many cells, and to one cell more than once. The default style is
    # "feature", one multiplier per lake: a lake's bed is one sediment, and a
    # leakance the stage data can inform lake by lake. A ``NONE`` bedleak (MF6: the
    # conductance is the aquifer's alone) on a selected lake is refused.
    "lak.bedleak": _Recipe(
        "lak.bedleak", "list", "lak", "connectiondata", "bedleak", use_col=5,
        index_cols=(2, 3), default_style="feature", feature="lake", key_cols=(0, 1),
        styles=("feature", "constant"),
    ),
    # Streambed K, per REACH, from PACKAGEDATA: ``ifno <cellid> rlen rwid rgrd rtp rbth
    # rhk man ncon ustrf ndv [aux] [boundname]``, rhk five columns after the cellid. A
    # reach with cellid ``NONE`` is one column SHORTER, and pyEMU, which reads columns by
    # position, multiplied that reach's Manning's n (0.03 became 0.06, measured), so such
    # a reach is refused. "constant" by default; "grid" is one per reach, keyed
    # (ifno, cellid) so pyEMU locates each reach at its cell; "feature" one per named
    # stream (BOUNDNAME); "zone" one per zones= label, keyed by the cell.
    "sfr.rhk": _Recipe(
        "sfr.rhk", "list", "sfr", "packagedata", "rhk", use_col=8, index_cols=(1, 2),
        feature="reach", key_cols=(0,), styles=("constant", "grid", "feature", "zone"),
    ),
    # Storage (ledger 548). Griddata, as K is: LAYERED writes one file per
    # layer, a scalar one whole-grid file `relayer_array_target` splits. MF6 reads
    # neither in a steady-state period, and SY only in a convertible cell
    # (ICONVERT != 0), so `refuse_storage_that_cannot_act` refuses a target that
    # could not move a head. Under STO's STORAGECOEFFICIENT option `ss` holds a
    # storage coefficient; the multiplier scales whichever it is.
    "sto.ss": _Recipe("sto.ss", "array", "sto", "ss"),
    "sto.sy": _Recipe("sto.sy", "array", "sto", "sy"),
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
    "bedleak": "lak.bedleak",
    "rhk": "sfr.rhk",
    "ss": "sto.ss",
    "specific_storage": "sto.ss",
    "sy": "sto.sy",
    "specific_yield": "sto.sy",
}

#: The styles :meth:`PestProject.parameterize` builds. Anything else is REFUSED at the
#: call -- pyEMU's own pilot-point spellings (``"pp"``, ``"pilot_points"``, ...)
#: reached pyEMU's structured-grid pilot points and died on a DISV grid with
#: ``'VertexGrid' object has no attribute 'delr'`` (ledger 520). ``"feature"`` -- one
#: parameter per lake or named stream -- acts on a FEATURE target alone (ledger
#: 558).
STYLES = ("constant", "zone", "grid", "pilotpoints", "feature")

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


#: FloPy's package types for the ARRAY form (READASARRAYS) of a package type that has
#: two forms; each is its list type's (``package_registry.packages_of_type``).
_ARRAY_FORM_TYPES = frozenset({"rcha", "evta"})


def array_form(recipe: _Recipe) -> _Recipe:
    """The recipe of a two-form target's ARRAY form: its ``array_dataset``, per period.

    Parameters
    ----------
    recipe : _Recipe
        A recipe whose package comes in both forms (``array_dataset`` set).

    Returns
    -------
    _Recipe
        The same target as an ``"array"`` recipe on ``array_dataset``, ``periodic``.

    Raises
    ------
    ValueError
        When the recipe's package has one form only.
    """

    if recipe.array_dataset is None:
        raise ValueError(f"Target {recipe.canonical!r} has no array form.")
    return replace(recipe, family="array", variable=recipe.array_dataset, field=None,
                   use_col=None, array_dataset=None, periodic=True)


def is_array_form(package) -> bool:
    """Whether a FloPy package is the ARRAY form (READASARRAYS) of its type: RCHA or EVTA.

    Parameters
    ----------
    package : flopy.mf6.MFPackage
        A FloPy package; its ``package_type`` is read.

    Returns
    -------
    bool
    """

    return str(getattr(package, "package_type", "")).lower() in _ARRAY_FORM_TYPES


def recipe_for_packages(recipe: _Recipe, packages) -> _Recipe:
    """The recipe in the FORM the packages it scales are written in (ledger 547).

    ``recharge`` scales a list RCH's ``recharge`` records and an array RCH's
    per-period ``recharge`` arrays; which one is the package's, as the staged model
    holds it -- never guessed. A recipe of one form is returned unchanged.

    Parameters
    ----------
    recipe : _Recipe
        The DECLARED recipe (:func:`resolve_target`).
    packages : list
        The FloPy packages it scales (:func:`select_packages`).

    Returns
    -------
    _Recipe
        ``recipe`` for a list package, :func:`array_form` of it for an array one.

    Raises
    ------
    ValueError
        When the packages named to share one parameter set mix the two forms.
    """

    if recipe.array_dataset is None:
        return recipe
    forms = {is_array_form(package) for package in packages}
    if len(forms) > 1:
        names = [f"{package_name_of(package)} ({'array' if is_array_form(package) else 'list'})"
                 for package in packages]
        raise ValueError(
            f"Target {recipe.canonical!r} on packages {names}: one parameter set cannot scale "
            "a list package's records and an array package's arrays together. Parameterize "
            "each with its own call: package=<name>."
        )
    return array_form(recipe) if forms == {True} else recipe


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


def checked_features(features, recipe: _Recipe) -> tuple[int | str, ...] | None:
    """``features=`` as :class:`NativeParameterSpec` holds it: a tuple of ids and boundnames.

    A feature is named as :class:`~myflopy.LakeStageTargets` names a lake: by MF6's
    number for it, ZERO-based, or by its boundname (ledger 558). Text
    that reads as a whole number is that number; other text is a boundname, matched
    later against the package (:func:`resolve_features`), any case.

    Parameters
    ----------
    features : int or str or sequence of int or str or None
        The lakes (``lak.bedleak``) or reaches/streams (``sfr.rhk``) to parameterize;
        one may be given bare. ``None``: every feature the style covers.
    recipe : _Recipe
        The target's recipe; its ``feature`` word ("lake", "reach") is read, and a
        target with none refuses any ``features``.

    Returns
    -------
    tuple of (int or str) or None
        Each number as ``int`` and each boundname stripped, in the order given.

    Raises
    ------
    ValueError
        When the target has no features, the sequence is empty, or an entry is blank,
        a fraction, negative or a bool.
    TypeError
        When an entry is neither a number nor text.
    """

    if features is None:
        return None
    if recipe.feature is None:
        featured = sorted(name for name, entry in _RECIPES.items() if entry.feature)
        raise ValueError(
            f"features={features!r} on target {recipe.canonical!r}: it names lakes or "
            f"reaches, and {recipe.canonical!r} scales a package that has none, so it "
            f"cannot act. Targets with features: {featured}."
        )
    items = [features] if isinstance(features, (str, int, float, np.integer, np.floating)) \
        else list(features)
    if not items:
        raise ValueError(f"features={features!r} names no {recipe.feature}; pass None for all.")
    checked: list[int | str] = []
    for item in items:
        def refuse(why: str, item=item) -> ValueError:
            return ValueError(
                f"features={features!r}: {item!r} {why}. A {recipe.feature} is MF6's "
                f"{recipe.feature} number, ZERO-based, or its boundname."
            )

        if isinstance(item, (bool, np.bool_)):
            raise refuse("is not a number or a boundname")
        if isinstance(item, str):
            text = item.strip()
            if not text:
                raise refuse("is blank")
            try:
                number = float(text)
            except ValueError:
                checked.append(text)
                continue
        elif isinstance(item, (int, float, np.integer, np.floating)):
            number = float(item)
        else:
            raise TypeError(
                f"features={features!r}: {item!r} is neither a {recipe.feature} number nor "
                "a boundname."
            )
        if not np.isfinite(number) or number != np.floor(number):
            raise refuse("is not a whole number")
        if number < 0:
            raise refuse("is negative: the first is 0")
        checked.append(int(number))
    return tuple(checked)


#: The one word ``periods=`` takes: a group per calendar month (ledger 549).
PERIODS_BY_MONTH = "month"

#: A ``periods=`` group label: it ends a PEST parameter name (``<name>_<label>``), so
#: letters, digits and ``_`` only, lower-cased.
_PERIOD_LABEL = re.compile(r"[a-z0-9_]+")


def _checked_period_numbers(label: str, members) -> tuple[int, ...]:
    """One ``periods=`` group's periods as sorted zero-based ``int``, refusing anything else."""

    values = [members] if np.ndim(members) == 0 else list(np.asarray(members, dtype=object).ravel())
    checked = []
    for value in values:
        try:
            number = float(value)
        except (TypeError, ValueError):
            number = float("nan")
        if isinstance(value, (bool, np.bool_, str, bytes)) or not np.isfinite(number) \
                or number < 0 or number != int(number):
            raise ValueError(
                f"periods= group {label!r}: each period is a zero-based whole number, and "
                f"{value!r} is not."
            )
        checked.append(int(number))
    if not checked:
        raise ValueError(f"periods= group {label!r} names no period; drop it.")
    if len(set(checked)) != len(checked):
        twice = sorted({value for value in checked if checked.count(value) > 1})[0]
        raise ValueError(f"periods= group {label!r} names period {twice} twice.")
    return tuple(sorted(checked))


def checked_periods(periods, recipe: _Recipe):
    """``periods=`` as :class:`NativeParameterSpec` holds it (ledger 549).

    Parameters
    ----------
    periods : str or mapping of str to sequence of int or None
        ``None``; ``"month"`` (any case); or ``{label: [zero-based periods]}``.
    recipe : _Recipe
        The target's DECLARED recipe; ``periods`` is refused unless its files are per
        stress period (:attr:`_Recipe.has_periods`).

    Returns
    -------
    str or dict of str to tuple of int or None
        ``None``, :data:`PERIODS_BY_MONTH`, or ``{label: periods}`` with each label
        stripped and lower-cased and its periods sorted, in the order given.

    Raises
    ------
    ValueError
        When the target's files are not per period; a string is not ``"month"``; a
        label is not letters, digits and ``_`` (or two labels are one once
        lower-cased); a period is not a whole number at or above 0; a group is
        empty or names a period twice; two groups share a period; or the mapping
        is empty.
    TypeError
        When it is neither a string, a mapping nor ``None``.
    """

    if periods is None:
        return None
    if not recipe.has_periods:
        raise ValueError(
            f"periods= on target {recipe.canonical!r}: its {recipe.variable} is not per "
            "stress period -- one value per cell for the whole run -- so there is nothing "
            "to group, and it cannot act. Drop periods=; it groups the per-period files of "
            "recharge and the list targets (chd, ghb, drn, riv, wel)."
        )
    if isinstance(periods, str):
        if periods.strip().lower() == PERIODS_BY_MONTH:
            return PERIODS_BY_MONTH
        raise ValueError(
            f"periods={periods!r}: the one word periods= takes is 'month' (a parameter per "
            "calendar month of the transient periods, m01..m12); otherwise give the groups, "
            "periods={'wet': [1, 2, 3], 'dry': [4, 5, 6]} (zero-based stress periods)."
        )
    if not isinstance(periods, Mapping):
        raise TypeError(
            f"periods= takes 'month' or {{label: [zero-based periods]}}, not "
            f"{type(periods).__name__}."
        )
    groups: dict[str, tuple[int, ...]] = {}
    owner: dict[int, str] = {}
    for label, members in periods.items():
        key = label.strip().lower() if isinstance(label, str) else None
        if key is None or not _PERIOD_LABEL.fullmatch(key):
            raise ValueError(
                f"periods= label {label!r}: a label ends the group's PEST parameter names "
                "(<name>_<label>), so it is text of letters, digits and _ only (e.g. "
                "'wet', 'q1')."
            )
        if key in groups:
            raise ValueError(f"periods= labels {label!r} and another are both {key!r}; rename one.")
        numbers = _checked_period_numbers(label, members)
        for number in numbers:
            if number in owner:
                raise ValueError(
                    f"periods= names period {number} in groups {owner[number]!r} and {key!r}; "
                    "a period's value is scaled by one group's parameter."
                )
            owner[number] = key
        groups[key] = numbers
    if not groups:
        raise ValueError("periods={} names no group; pass None for one parameter set over every period.")
    return groups


@dataclass(frozen=True)
class PeriodGroup:
    """One parameter set of a target's per-period files, as ``periods=`` groups them.

    Attributes
    ----------
    label : str or None
        The group's label (``"m03"``, ``"wet"``); ``None`` when ``periods`` is not
        given -- one set over every file.
    periods : tuple of int
        The zero-based stress periods whose values the group's files hold: each
        file's own period and every later period that repeats its block.
    targets : tuple of TargetFile
        The files the group's parameters scale.
    """

    label: str | None
    periods: tuple[int, ...]
    targets: tuple[TargetFile, ...]

    @property
    def name_suffix(self) -> str:
        """``"_<label>"``, appended to the parameter name; ``""`` with no label."""

        return "" if self.label is None else f"_{self.label}"


def month_labels(calendar: PeriodCalendar) -> list[str | None]:
    """Each stress period's ``periods="month"`` group: ``m01``..``m12``, ``None`` when steady.

    A TRANSIENT period is labelled by the calendar month of its MIDPOINT
    (:attr:`PeriodCalendar.months`, ledger 529): Ten Trails starts 2017-09-30, so
    each period's START date names the month before. A steady-state period is in no
    group: it is the long-term average, and only a parameter declared without
    ``periods`` reaches it (the user's rule, Q3).

    Parameters
    ----------
    calendar : PeriodCalendar
        The staged run's calendar (:func:`staged_calendar`); a calendar with no
        steady flags reads every period as transient.

    Returns
    -------
    list of str or None
        One entry per period, zero-based.

    Raises
    ------
    ValueError
        When a transient period is longer than 31 days: it spans more than one
        calendar month, so its month is no month's.
    """

    steady = calendar.steady or (False,) * calendar.nper
    lengths = (calendar.ends - calendar.starts) / pd.Timedelta(days=1)
    labels: list[str | None] = []
    for per in range(calendar.nper):
        if steady[per]:
            labels.append(None)
            continue
        if lengths[per] > 31.0 + 1e-9:
            raise ValueError(
                f"periods='month': transient period {per} runs {lengths[per]:g} days "
                f"({calendar.starts[per].date()} to {calendar.ends[per].date()}), more than "
                "one calendar month, so it is no month's. Give the groups yourself, "
                "periods={'wet': [...], 'dry': [...]}, or split the run into monthly periods "
                "(mf.AverageYear)."
            )
        labels.append(f"m{int(calendar.months[per]):02d}")
    return labels


def file_spans(targets: list[TargetFile], nper: int) -> dict[TargetFile, range]:
    """The stress periods each per-period file holds the values of, by MF6's block rule.

    A file serves its own block's period and every later period of its package up
    to the package's next block (or the last period): MF6 repeats a block until the
    next one. The next block may be EMPTY (``stress_period_data={2: []}``), which has
    no file of its own and still ends the block before it, so it comes from FloPy's
    records (:attr:`TargetFile.next_block`) as well as from the other files. Files
    are numbered by block ORDINAL (an array RCH's ``m.rcha_recharge_2.txt`` is the
    period-3 block's when only periods 1, 3 and 4 have blocks), so this reads each
    file's ``period`` from FloPy's records (:attr:`TargetFile.period`), never from
    its name.

    Parameters
    ----------
    targets : list of TargetFile
        Per-period files, each with its zero-based ``period``.
    nper : int
        The number of stress periods.

    Returns
    -------
    dict of TargetFile to range
        Each file's periods, in the order ``targets`` gives them.
    """

    starts: dict[str, list[int]] = {}
    for target in targets:
        starts.setdefault(target.package, []).append(int(target.period))
    for periods in starts.values():
        periods.sort()
    spans = {}
    for target in targets:
        ends = [period for period in starts[target.package] if period > target.period][:1]
        if target.next_block is not None:
            ends.append(int(target.next_block))
        spans[target] = range(int(target.period), min(ends) if ends else int(nper))
    return spans


def _describe_group(label: str | None, by_month: bool) -> str:
    if label is not None:
        return f"group {label!r}"
    return "a steady-state period, in no seasonal group" if by_month else "in no group"


def group_period_files(spec, targets: list[TargetFile], labels: list[str | None]) -> list[PeriodGroup]:
    """Split a target's period files into ``periods=`` groups (ledger 549).

    MF6 repeats a period's block until the next one, so a file serves its own
    period and every later period of its package without a block of its own. A
    file is in the group of every period it serves -- and when those are not one
    group (a later period in another group, or in none, a steady-state period
    under ``"month"`` included), one multiplier could not be one group's without
    being the other's, so it is REFUSED, naming the fix. Files serving periods in
    no group are left out: a ``"month"`` group never scales the steady period.

    Parameters
    ----------
    spec : NativeParameterSpec
        The parameter (named by the refusals); its ``periods`` decides the mode.
    targets : list of TargetFile
        Every file of the target, each with its ``period`` (:func:`target_files`).
    labels : list of str or None
        Each stress period's group, zero-based, ``len == nper``
        (:func:`month_labels`, or a ``periods=`` mapping spread over the periods).

    Returns
    -------
    list of PeriodGroup
        One per label, in month order for ``"month"`` and the mapping's order
        otherwise; only groups holding a file.

    Raises
    ------
    ValueError
        When a file serves periods of two groups (or of a group and none); when a
        group's periods are served by no file of the target's packages -- before
        their first block, or under an empty block -- so nothing scales them; or
        when no group holds a file.
    """

    by_month = spec.periods == PERIODS_BY_MONTH
    by_package = sorted({target.package for target in targets})
    grouped: dict[str, list[TargetFile]] = {}
    served: dict[str, set[int]] = {}
    for target, span in file_spans(targets, len(labels)).items():
        package = target.package
        own = labels[target.period]
        for period in span:
            if labels[period] != own:
                raise ValueError(
                    f"periods= on {spec.recipe.canonical!r} (name {spec.name!r}): package "
                    f"{package!r} declares no block for period {period} "
                    f"({_describe_group(labels[period], by_month)}), so MF6 repeats period "
                    f"{target.period}'s ({_describe_group(own, by_month)}) -- one file, "
                    f"{target.filename}, serves both, and its multiplier cannot be one's "
                    "without being the other's. Declare a block for every period -- a value "
                    "per period, mf.rch.array(recharge={period: ...}) or "
                    "mf.chd(stress_period_data={period: ...}); mf.AverageYear(...).by_month(...) "
                    "gives each period its own."
                )
        if own is None:
            continue
        grouped.setdefault(own, []).append(target)
        served.setdefault(own, set()).update(span)
    order = sorted(set(label for label in labels if label is not None)) if by_month \
        else list(spec.periods)
    if not by_month:
        empty = [label for label in order if label not in grouped]
        if empty:
            label = empty[0]
            raise ValueError(
                f"periods= group {label!r} of {spec.recipe.canonical!r} (name {spec.name!r}) "
                f"holds periods {list(spec.periods[label])}, which no block of package(s) "
                f"{sorted(by_package)} holding records serves -- they come before its first "
                "block, or under an empty one: MF6 holds no value of them to scale, so its "
                "parameter cannot act. Drop the group, or declare blocks for its periods."
            )
    groups = [
        PeriodGroup(label, tuple(sorted(served[label])), tuple(grouped[label]))
        for label in order if label in grouped
    ]
    if not groups:
        raise ValueError(
            f"periods={spec.periods!r} on {spec.recipe.canonical!r} (name {spec.name!r}) "
            "groups none of its files -- every period it reaches is in no group -- so it "
            "cannot act."
        )
    return groups


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
        (one per list entry / cell), ``"pilotpoints"`` or ``"feature"`` (one per lake
        or named stream, a FEATURE target's alone), any case. ``None`` (default) uses
        the recipe default: ``"feature"`` for ``lak.bedleak``, ``"constant"`` for every
        other target. Anything else raises ``ValueError`` naming the five --
        pyEMU's own pilot-point spellings (``"pp"``, ``"pilot_points"``, ...) name
        ``"pilotpoints"`` -- and so does a style the target does not take:
        ``"pilotpoints"`` on a LIST target, which no engine here builds (ledger
        520; on ``recharge``, whose form is its package's, that refusal comes from
        :meth:`resolve_form` once the package is known), ``"feature"`` on a target
        with no features, and ``"grid"``/``"zone"``/``"pilotpoints"`` on
        ``lak.bedleak`` (ledger 558).
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
        Default ``None``. With ``style="constant"`` on a LIST target, which pyEMU
        scales with one number whatever the labels, it raises ``ValueError`` naming
        ``style="zone"`` -- here, or for ``recharge`` from :meth:`resolve_form` once
        its package's form is known (ledger 569).
    layers : int or iterable of int or None, optional
        Zero-based layers an ARRAY target is restricted to, held as a tuple of
        ``int`` (:func:`checked_layers`, which refuses anything else). Default
        ``None``: every layer. Given for a LIST target it raises ``ValueError``
        naming ``package=`` -- nothing selects a list record by layer, and ``zones=``
        leaves none out (with ``style="zone"`` every label, ``0`` too, is a parameter;
        ledger 559); for ``recharge`` in either form (MF6 applies it to one cell
        per column) it raises naming ``zones=``.
    package : str or tuple of str or None, optional
        Which package(s) of the target's TYPE are scaled, by NAME
        (:func:`checked_package_argument`). ``None`` (default): the model's only
        package of the type; several raise ``SeveralPackagesError`` at
        :meth:`PestProject.parameterize` or :meth:`PestProject.build`. A tuple shares
        one parameter set across the packages it names; refused with ``features`` or
        ``style="feature"``, which number features within ONE package, and at
        :meth:`PestProject.build` when the packages' records differ in width -- one
        declares BOUNDNAMES or AUXILIARY the other does not (ledger 559).
    features : int or str or sequence of int or str or None, optional
        The lakes (``lak.bedleak``) or reaches and named streams (``sfr.rhk``) the
        parameter covers, each MF6's number for it, ZERO-based, or its boundname (any
        case), held as a tuple (:func:`checked_features`). Default ``None``: every
        lake, every reach -- and, for ``style="feature"`` on ``sfr.rhk``, every
        boundname. With ``"feature"`` each entry is ONE parameter; with
        ``"constant"`` they share one; with ``"grid"`` each reach among them is one.
        Refused on a target with no features, with ``style="zone"``, and -- at
        :meth:`PestProject.parameterize` when the model lists its packages, else at
        :meth:`PestProject.build` -- naming no feature of the package, or a feature
        twice (ledger 558).
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
        raises ``ValueError``: neither has anything to correlate. On a LIST
        ``"grid"`` whose parameters pyEMU places at one point in one covariance
        group -- a column's records in two layers on DISV/DISU -- or then splits one
        per group -- every reach of an ``sfr.rhk`` set on DISV, one record per layer
        of a lone DIS column -- :meth:`PestProject.build` raises ``ValueError``
        (ledger 570).
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
    periods : str or mapping of str to sequence of int or None, optional
        Groups the target's per-period files, one parameter set per group (ledger
        549); see :meth:`PestProject.parameterize`. ``None`` (default):
        one set over every period. ``"month"``: one group per calendar month of the
        TRANSIENT periods, by each period's midpoint, labelled ``m01``..``m12``.
        A mapping ``{label: [zero-based periods]}``: explicit groups, each label
        letters, digits and ``_``. Held as ``"month"``, ``None`` or a dict of
        lower-cased labels to sorted tuples of ``int`` (:func:`checked_periods`,
        which refuses anything else, and refuses it on a target whose files are not
        per period).
    name : str or None, optional
        Parameter group / name base. Default ``None``: the canonical target with
        dots removed, plus ``_<package>`` when ``package`` is given
        (``"rivcond_crisp"``). Pilot-point sets are named ``<name>l<layer>`` /
        ``<name>l<layer>z<zone>``, lower-cased.
    capture : bool, optional
        Also record the resolved per-cell field as zero-weight observations, for
        ``IesResults.plot_field`` / ``IesResults.field``. Default ``False``. ``True``
        on a FEATURE target (``lak.bedleak``, ``sfr.rhk``) raises ``ValueError``: its
        rows are lake connections and reaches, not one value per cell (deferred, ledger
        558).
    extra : dict, optional
        Keyword arguments passed to pyEMU's ``PstFrom.add_parameters``. Default
        empty. Not empty for array-target pilot points, which never call it, raises
        ``ValueError``. A key the facade sets itself (``_FACADE_KEYWORDS``:
        ``par_type``, ``par_name_base``, ``pargp``, the bounds, ``transform``,
        ``ult_lbound``/``ult_ubound``, ``par_style``, ``zone_array``, ``geostruct``,
        ``index_cols``, ``use_cols``, ``use_rows``) raises ``TypeError`` naming the
        argument that sets it (ledger 571), and so does a key
        ``add_parameters`` does not take.
    resolved_files : list of str, optional
        Filled by :meth:`PestProject.build`: the template-relative files the
        parameter scales. Default empty.
    pilot_point_sets : list of dict, optional
        Filled by :meth:`PestProject.build` for array-target pilot points: one
        summary per pilot-point set (``name``, ``layer`` -- ``None`` for a set on
        period files --, ``zone``, ``period_group``, ``points``, ``cells``). Default
        empty.

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
    declared_recipe : _Recipe
        The recipe the target names, before its package's FORM is known: for
        ``recharge`` the list form, whose array form :attr:`recipe` becomes when
        the staged RCH is READASARRAYS (:func:`recipe_for_packages`).
    period_groups : list of PeriodGroup
        Filled by :meth:`PestProject.build`: the parameter sets ``periods`` made,
        each with its label, the periods its files serve and the files
        (:func:`group_period_files`); one unlabelled group without ``periods``.
    feature_sets : list of dict
        Filled by :meth:`PestProject.build` for a FEATURE target: one record per
        parameter set the features made -- its PEST group ``pargp``, ``label``
        (``None`` but for ``"feature"``), ``package``, the zero-based ``features``
        (lakes or reaches) it covers and the number of file ``rows`` it scales.
    period_group_pargps : list of list of str
        Filled by :meth:`PestProject.build`, one list per :attr:`period_groups`
        entry in its order: the PEST groups the build gave that group's parameters
        -- ``<name>_<label>``, or one per layer file of a ``"grid"`` array
        (``<name>_<label>l<layer>``) and one per pilot-point set
        (``<name>_<label>[l<layer>][z<zone>]``). ``pargp`` in ``extra`` is refused
        (ledger 571), so these are the groups.
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
    features: tuple[int | str, ...] | None = None
    pp_spacing: float | None = None
    pp_points: Any = None
    pp_nearest: int = DEFAULT_PP_NEAREST
    correlation: float | None = None
    anisotropy: float = 1.0
    bearing: float = 0.0
    nugget: float = 0.0
    periods: str | Mapping[str, Sequence[int]] | None = None
    name: str | None = None
    capture: bool = False
    extra: dict = field(default_factory=dict)
    resolved_files: list[str] = field(default_factory=list)
    pilot_point_sets: list[dict] = field(default_factory=list)

    def __post_init__(self):
        """Resolve the recipe, fill the defaults, and refuse an argument that cannot act."""

        self.declared_recipe = resolve_target(self.target)
        self.resolved_targets: list[TargetFile] = []
        self.resolved_zone_labels: list[int] | None = None
        self.period_groups: list[PeriodGroup] = []
        self.period_group_pargps: list[list[str]] = []
        self.feature_sets: list[dict] = []
        declared = self.declared_recipe
        self.style = _checked_style(self.style if self.style is not None else declared.default_style)
        _refuse_a_style_the_target_does_not_take(declared, self.style)
        # A target of two forms (`recharge`) is held in the form its style can only
        # have -- pilot points exist on arrays alone -- until the package's form is
        # known (`resolve_form`), which refuses pilot points on a list package.
        self.recipe = (
            array_form(declared)
            if declared.array_dataset is not None and self.style == "pilotpoints" else declared
        )
        recipe = self.recipe
        self.package = checked_package_argument(self.package)
        self.features = checked_features(self.features, declared)
        if declared.feature is not None:
            self._refuse_feature_arguments_that_cannot_act()
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
        if self.layers is not None and declared.array_dataset is not None:
            raise ValueError(
                f"layers= on target {declared.canonical!r}: MF6 applies it to ONE cell of "
                "each column (IRCH, or the top active cell), in either form, so it has no "
                "layers to restrict. Zone it by column with zones= -- one label per "
                "cell of a layer: on an array RCH a label <= 0 is left alone, while on a "
                "list RCH, with style='zone', every label -- 0 and negatives too -- is a "
                "parameter of its own (ledger 559)."
            )
        if self.layers is not None and recipe.family == "list":
            raise ValueError(
                f"layers= on list target {recipe.canonical!r}: it scales every record of "
                "its package's files, whatever layer each names, so layers= cannot act. "
                "Nor can zones= leave records out: with style='zone' every label -- 0 "
                "and negatives too -- is a parameter of its own, and on DIS/DISV one "
                "label covers a whole column, every layer alike ('grid' only names its "
                "records by it, 'constant' ignores it). Pick the package with package= "
                "(ledger 559)."
            )
        if recipe.family == "list" and declared.array_dataset is None:
            # A target of two forms (`recharge`) is asked once its form is known.
            self._refuse_zones_a_list_constant_ignores(recipe)
        self.periods = checked_periods(self.periods, declared)
        _refuse_extra_keywords_the_facade_sets(self.extra)
        _refuse_extra_keywords_pyemu_does_not_take(self.extra)
        self._check_pilot_point_arguments()
        self._refuse_an_idle_variogram()
        self.bounds = self._checked_bounds()
        if self.name is None:
            slug = _package_slug(self.package)
            base = recipe.canonical.replace(".", "")
            self.name = f"{base}_{slug}" if slug else base

    def _refuse_feature_arguments_that_cannot_act(self) -> None:
        """Refuse an argument a FEATURE target (LAK, SFR) leaves idle (ledger 558)."""

        recipe = self.declared_recipe
        plural = "lakes" if recipe.feature == "lake" else "reaches"
        if self.zones is not None and self.style != "zone":
            raise ValueError(
                f"zones= on {recipe.canonical!r} with style={self.style!r}: zones= labels "
                f"CELLS, its rows are {plural}, and only style='zone' reads a row's cell "
                "label -- pyEMU ignores zones= on a 'constant' or 'feature' parameter -- "
                f"so it cannot act. Choose the {plural} with features="
                + (", or use style='zone'." if "zone" in recipe.styles else ".")
            )
        if self.features is not None and self.style == "zone":
            raise ValueError(
                f"features= with style='zone' on {recipe.canonical!r}: a zone parameter "
                "covers every reach in its zones= cells, so features= cannot narrow it. "
                "Use zones= alone, or features= with style='constant'/'grid'/'feature'."
            )
        if isinstance(self.package, tuple) and (self.features is not None
                                                 or self.style == "feature"):
            raise ValueError(
                f"package={self.package!r} with "
                + ("features=" if self.features is not None else "style='feature'")
                + f" on {recipe.canonical!r}: {plural} are numbered and named within ONE "
                "package, so a set shared across packages cannot name them. Parameterize "
                "each package with its own call: package=<name>."
            )
        if self.capture:
            raise ValueError(
                f"capture=True on {recipe.canonical!r}: a capture records one value per "
                f"CELL, and its rows are {plural} -- a lake connects to one cell more than "
                "once, two reaches can share a cell -- so it is deferred (ledger "
                "558). Drop capture=; the parameter values are in the "
                "ensembles."
            )

    def _refuse_zones_a_list_constant_ignores(self, recipe: _Recipe) -> None:
        """Refuse ``zones=`` on a LIST ``style="constant"`` parameter (ledger 569).

        pyEMU reads ``zone_array`` only for a ``"zone"`` or ``"grid"`` list
        parameter, so a constant one scaled every record with one number and the zones
        did nothing (measured, ledger 559). A FEATURE target refuses it in
        :meth:`_refuse_feature_arguments_that_cannot_act`.
        """

        if self.zones is None or self.style != "constant" or recipe.feature is not None:
            return
        raise ValueError(
            f"zones= with style='constant' on list target {recipe.canonical!r}: pyEMU reads "
            "zones= only for a 'zone' or 'grid' list parameter, so one "
            + ("offset" if self.additive else "multiplier")
            + " would scale every record of the package and zones= cannot act. Use "
            "style='zone' -- one parameter per label, and on a list target every label, 0 "
            "and negatives too, is one (ledger 559) -- or drop zones= (ledger "
            "569)."
        )

    def _refuse_an_idle_variogram(self) -> None:
        """Refuse a variogram argument the spec's style or other arguments leave idle."""

        if self.correlation is not None:
            if self.style == "feature":
                raise ValueError(
                    f"correlation={self.correlation!r} with style='feature': a feature "
                    "parameter is one multiplier per lake or named stream, keyed by its "
                    "number, not located at a cell, so a variogram has nothing to act "
                    "on. Drop it, or use style='grid' (sfr.rhk: one per reach, at its "
                    "cell)."
                )
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
                    "only with style='pilotpoints' on an array target (k, k33, sto.ss, "
                    "sto.sy, porosity, recharge on an array RCH)."
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

    def resolve_form(self, packages) -> _Recipe:
        """Set :attr:`recipe` to the FORM of the packages it scales, and return it.

        :func:`recipe_for_packages` on :attr:`declared_recipe`; called by
        :meth:`PestProject.parameterize` when the model can list its packages and by
        :meth:`PestProject.build` on the staged model.

        Parameters
        ----------
        packages : list
            The FloPy packages the spec scales (:func:`select_packages`).

        Returns
        -------
        _Recipe
            The recipe now held.

        Raises
        ------
        ValueError
            When the packages mix the two forms, or ``style="pilotpoints"`` was given
            for a target whose package is the LIST form (pilot points interpolate onto
            an array, ledger 520), or ``zones=`` with ``style="constant"`` for one
            (pyEMU ignores it there, ledger 569).
        """

        recipe = recipe_for_packages(self.declared_recipe, packages)
        if recipe.family == "list":
            self._refuse_zones_a_list_constant_ignores(recipe)
        if self.style == "pilotpoints" and recipe.family == "list":
            raise ValueError(
                f"style='pilotpoints' on {recipe.canonical!r}: package(s) "
                f"{[package_name_of(package) for package in packages]} hold it as LIST "
                "records, and myflopy's pilot points interpolate onto an ARRAY. Write the "
                "package in its array form (mf.rch.array, READASARRAYS) to put pilot points "
                "on it, or use style='grid' (one parameter per record; correlation= "
                "correlates them), 'zone' or 'constant' (ledger 520)."
            )
        self.recipe = recipe
        return recipe

    @property
    def resolved_transform(self) -> str:
        """Transform actually used (additive parameters cannot be log)."""

        return "none" if self.additive else self.transform


#: The ``PstFrom.add_parameters`` keywords the facade sets itself, each with the
#: :meth:`PestProject.parameterize` argument (or what else) that sets it. An
#: ``**extra`` key among them is refused (ledger 571).
_FACADE_KEYWORDS = {
    "par_type": "style=",
    "par_name_base": "name=",
    "pargp": "name=",
    "lower_bound": "bounds=",
    "upper_bound": "bounds=",
    "transform": "transform=",
    "ult_lbound": "physical=",
    "ult_ubound": "physical=",
    "par_style": "additive=",
    "zone_array": "zones=",
    "geostruct": "correlation= (with anisotropy=, bearing=, nugget=)",
    "index_cols": "the target: its recipe names the cellid columns on the model's grid",
    "use_cols": "the target: one per column, e.g. 'riv.cond' or 'riv.stage'",
    "use_rows": "features= on lak.bedleak/sfr.rhk, package= on any list target",
}


def _refuse_extra_keywords_the_facade_sets(extra: dict) -> None:
    """Refuse an ``**extra`` key the facade itself hands ``PstFrom.add_parameters``.

    They were applied last and overrode what the facade derived, so an ``index_cols``
    or ``use_rows`` there failed only after the whole staging -- or, under
    ``use_rows``, scaled rows the template check then refused -- and a ``pargp`` or
    ``par_type`` moved what the build's metadata records out from under it (ledger
    571). Each is named with the argument that sets it.
    """

    colliding = [key for key in extra if key in _FACADE_KEYWORDS]
    if colliding:
        raise TypeError(
            "parameterize() sets "
            + ", ".join(f"{key}= through {_FACADE_KEYWORDS[key]}" for key in colliding)
            + " itself, so an extra keyword of that name would override it; give it there "
            "(ledger 571)."
        )


def _refuse_extra_keywords_pyemu_does_not_take(extra: dict) -> None:
    """Refuse an ``**extra`` key that is no parameter of pyEMU's ``PstFrom.add_parameters``.

    ``**extra`` is passed there at :meth:`PestProject.build`, so a misspelt key -- or
    ``temporal=``, which ``periods=`` replaced (ledger 549) -- reached it
    only after the run was staged and died with pyEMU's ``unexpected keyword
    argument``. Asked at the call, of pyEMU's own signature; without pyEMU nothing is
    asked (the build needs it anyway).
    """

    if not extra:
        return
    try:
        import pyemu
    except ImportError:
        return
    accepted = set(inspect.signature(pyemu.utils.PstFrom.add_parameters).parameters)
    unknown = sorted(set(extra) - (accepted - {"self", "filenames"}))
    if unknown:
        raise TypeError(
            f"parameterize() got {', '.join(f'{key}=' for key in unknown)}, a keyword neither "
            "it nor pyEMU's PstFrom.add_parameters (where extra keywords go) takes. Seasonal "
            "or other period groups are periods= ('month' or {label: [periods]})."
        )


def _refuse_a_style_the_target_does_not_take(recipe: _Recipe, style: str) -> None:
    """Refuse a style the recipe's target cannot build (ledger 558).

    ``"feature"`` needs features (lakes, reaches) to make one parameter each; a
    lake-bed leakance is keyed by lake and connection, not by a cell, so pyEMU's
    ``"grid"``/``"zone"`` -- which locate a record at its cell -- and pilot points do
    not act on it.
    """

    allowed = recipe.styles or tuple(name for name in STYLES if name != "feature")
    if style in allowed:
        return
    if style == "feature":
        featured = sorted(name for name, entry in _RECIPES.items() if entry.feature)
        raise ValueError(
            f"style='feature' on target {recipe.canonical!r}: it makes one parameter per "
            f"lake or named stream, and {recipe.canonical!r} scales a package that has "
            f"none, so it cannot act. Targets with features: {featured}."
        )
    raise ValueError(
        f"style={style!r} on target {recipe.canonical!r}: it takes "
        f"{', '.join(repr(name) for name in allowed)}. A lake-bed leakance is keyed by "
        "lake and connection, not by a cell, so a style that locates a record at its "
        "cell cannot act on it: 'feature' is one multiplier per lake (the default), "
        "'constant' one for every lake (ledger 558)."
        if recipe.feature == "lake" else
        f"style={style!r} on target {recipe.canonical!r}: it takes "
        f"{', '.join(repr(name) for name in allowed)} (ledger 558)."
    )


def _checked_style(style) -> str:
    """``style`` lower-cased, refusing anything but the :data:`STYLES`."""

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
    next_block : int or None
        The zero-based period of the package's NEXT block, an EMPTY one included
        (``stress_period_data={2: []}``, which FloPy writes as ``BEGIN period 3`` /
        ``END period 3`` with no file): MF6 repeats this block until then, so the file
        serves ``period`` to ``next_block - 1`` (:func:`file_spans`). ``None`` for the
        package's last block, for data that is not per period, and for a file built by
        hand, whose span :func:`file_spans` then reads off the other files.
    """

    filename: str
    package: str
    layer: int | None = None
    period: int | None = None
    next_block: int | None = None


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


def _block_is_empty(storage) -> bool:
    """Whether a period block's storage holds no data: an EMPTY block, not a missing file.

    FloPy keeps ``stress_period_data={2: []}`` as a block of its own (``BEGIN period
    3`` / ``END period 3``, from which MF6 has no boundary of the package) and
    externalizes nothing for it, so it names no file and that is right.
    """

    has_data = getattr(storage, "has_data", None)
    return callable(has_data) and not has_data()


def _dataset_files(data) -> list[tuple[int | None, int | None, str | None, bool]]:
    """``(period, layer, file name, empty)`` for every block of a FloPy dataset.

    Read off the data's STORAGE, never ``get_record()``, which loads the data and, for
    a list of a model kept in a sub-directory, opened it at the model's path where
    MF6 (and FloPy's writer) put it at the simulation's. An empty period block
    (:func:`_block_is_empty`) is one entry with no file and ``empty`` True.
    """

    keys = getattr(data, "get_active_key_list", None)
    if callable(keys):
        found = []
        for key, storage in keys():
            files = _storage_files(storage)
            if not any(filename for _layer, filename in files) and _block_is_empty(storage):
                found.append((int(key), None, None, True))
                continue
            found.extend((int(key), layer, filename, False) for layer, filename in files)
        return found
    return [(None, layer, filename, False)
            for layer, filename in _storage_files(data._get_storage_obj())]


def package_target_files(model, recipe: _Recipe, packages, workspace) -> list[TargetFile]:
    """The external files FloPy's records say ``packages`` read for ``recipe``'s dataset.

    Read off the packages' own data storage -- each period's of a transient dataset
    (``get_active_key_list()``), each layer's of a layered array, one otherwise --
    after ``set_all_data_external()``, which names every file in memory, relative to
    the simulation's directory as MF6 reads it. Never a file name built from the type: FloPy names an
    external file after the package's FILE (``m.crisp_stress_period_data_1.txt``), so
    a glob on the type scaled one of two GHBs and found nothing for a GHB named
    ``edge`` (ledger 517); and a stale whole-grid file a relayering leaves
    behind is never named. An EMPTY period block (``stress_period_data={2: []}``)
    holds no file and nothing to scale: it is left out, and it ends the block before
    it (:attr:`TargetFile.next_block`), since MF6 has no boundary of the package from
    it on.

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
        Package by package as given, each in period then layer order, every period
        file carrying its package's next block (an empty one included) as
        ``next_block``.

    Raises
    ------
    ValueError
        When a package holds no such dataset, or a block that holds data is not in an
        external file. (The dataset is the package's FORM's -- an array RCH's
        ``recharge`` arrays, a list RCH's records -- since :func:`target_files` hands
        over the recipe :meth:`NativeParameterSpec.resolve_form` chose, ledger 547.)
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
                f"{getattr(model, 'name', '?')!r} ({type(package).__name__}) has none."
            )
        entries = _dataset_files(data)
        blocks = sorted({period for period, _layer, _filename, _empty in entries
                         if period is not None})
        for period, layer, filename, empty in entries:
            if empty:
                continue
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
            later = [block for block in blocks if period is not None and block > period]
            found.append(TargetFile(str(filename), name, layer, period,
                                    later[0] if later else None))
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


class _StagedFlowModel:
    """A staged FloPy model as the calendar's readers read a myflopy model.

    :func:`~myflopy.modflow.mf6.period_calendar.steady_periods`,
    :func:`~myflopy.modflow.mf6.period_end.period_lengths` and
    :func:`~myflopy.modflow.mf6.hydraulics.declared_time_unit` read ``model.sim``
    and ``model._packages_of_type``; a FloPy model has ``simulation`` and finds a
    package by TYPE through :func:`packages_of_type`. The ONE reading of a period's
    dates and steady flag stays theirs (ledger 529).
    """

    def __init__(self, gwf):
        self.gwf = gwf
        self.sim = gwf.simulation
        self.name = gwf.name
        self.model_type = str(getattr(gwf, "model_type", "gwf6") or "gwf6")

    def _packages_of_type(self, package_type: str) -> list:
        return packages_of_type(self.gwf, package_type)


def staged_nper(model) -> int:
    """The number of stress periods of the staged FloPy ``model``'s simulation.

    Parameters
    ----------
    model : flopy.mf6.MFModel
        A model of the build's staged simulation.

    Returns
    -------
    int

    Raises
    ------
    ValueError
        When its simulation's TDIS has no period block.
    """

    lengths = period_lengths(_StagedFlowModel(model))
    if not lengths:
        raise ValueError(f"The run's TDIS (model {model.name!r}) holds no period block.")
    return len(lengths)


def staged_calendar(project, model) -> PeriodCalendar:
    """The :class:`PeriodCalendar` of the staged run, for ``periods="month"``.

    The TDIS lengths, ``TIME_UNITS`` and the flow model's STO steady flags of the
    run's FILES (ledger 356), read by the calendar's own readers
    (:class:`_StagedFlowModel`); the start is the TDIS ``START_DATE_TIME``, or,
    when the TDIS declares none, the calibration's ``start_datetime``.

    Parameters
    ----------
    project : PestProject
        The calibration being built; its ``start_datetime`` is read.
    model : flopy.mf6.MFModel
        The staged FLOW model.

    Returns
    -------
    PeriodCalendar

    Raises
    ------
    ValueError
        When the TDIS's ``START_DATE_TIME`` and the calibration's ``start_datetime``
        are both declared and differ (a month would depend on which is read); when
        neither is; or when the TDIS declares no ``TIME_UNITS`` MF6 reads.
    """

    flow = _StagedFlowModel(model)
    declared = getattr(project, "start_datetime", None)
    declared = None if declared is None else parse_start(
        declared, what="PestProject(start_datetime=)")
    held = tdis_start_date_time(flow)
    start = None if held is None else parse_start(held, what="The TDIS START_DATE_TIME")
    if start is not None and declared is not None and start != declared:
        raise ValueError(
            f"periods='month' dates the periods, and the run's TDIS starts them on {held!r} "
            f"while the calibration's start_datetime is {project.start_datetime!r}: a "
            "period's month would depend on which is read. Make them one -- model.pest() "
            "takes start_datetime from the TDIS when it is not given."
        )
    start = start if start is not None else declared
    if start is None:
        raise ValueError(
            "periods='month' dates the periods, and neither the run's TDIS "
            "(START_DATE_TIME) nor the calibration (start_datetime=) gives a start."
        )
    units = declared_time_unit(flow)
    if units is None:
        raise ValueError(
            "periods='month' dates the periods, and the run's TDIS declares no TIME_UNITS "
            "MF6 reads, so its lengths cannot become dates; declare one, "
            "mf.tdis(time_units='days', ...), or give the groups, periods={label: [...]}."
        )
    perlen = period_lengths(flow)
    return PeriodCalendar(start, perlen, units, steady_periods(flow, len(perlen)))


def period_labels(project, spec: NativeParameterSpec, model) -> list[str | None]:
    """Each stress period's ``periods=`` group, zero-based (``None``: in no group).

    Parameters
    ----------
    project : PestProject
        The calibration being built.
    spec : NativeParameterSpec
        A spec whose ``periods`` is ``"month"`` or a mapping.
    model : flopy.mf6.MFModel
        The staged model the target's packages belong to.

    Returns
    -------
    list of str or None

    Raises
    ------
    ValueError
        When a mapping names a period the run does not have; or what
        :func:`staged_calendar` and :func:`month_labels` refuse.
    """

    if spec.periods == PERIODS_BY_MONTH:
        return month_labels(staged_calendar(project, model))
    nper = staged_nper(model)
    labels: list[str | None] = [None] * nper
    for label, periods in spec.periods.items():
        for period in periods:
            if period >= nper:
                raise ValueError(
                    f"periods= group {label!r} names period {period}, and the run has "
                    f"{nper} (0..{nper - 1})."
                )
            labels[period] = label
    return labels


def period_groups_for(project, spec: NativeParameterSpec, model,
                      targets: list[TargetFile]) -> list[PeriodGroup]:
    """The parameter sets ``spec`` scales its files in: ``periods=`` groups, else one.

    Parameters
    ----------
    project : PestProject
        The calibration being built.
    spec : NativeParameterSpec
        The parameter.
    model : flopy.mf6.MFModel
        The staged model the files belong to.
    targets : list of TargetFile
        The target's files (:func:`package_target_files`).

    Returns
    -------
    list of PeriodGroup
        Without ``periods``: one unlabelled group of every file, its ``periods`` the
        periods they serve (none for data that is not per period).
    """

    if spec.periods is None:
        served: set[int] = set()
        if all(target.period is not None for target in targets):
            for span in file_spans(targets, staged_nper(model)).values():
                served.update(span)
        return [PeriodGroup(None, tuple(sorted(served)), tuple(targets))]
    return group_period_files(spec, targets, period_labels(project, spec, model))


def refuse_storage_that_cannot_act(project, spec: NativeParameterSpec, model, packages,
                                   targets: list[TargetFile]) -> None:
    """Refuse an STO target that could move no head (ledger 548).

    MF6 reads storage only in a TRANSIENT period, and SY only in a convertible cell
    (``ICONVERT != 0``): a parameter on storage of an all-steady run, or on SY of
    cells that are all confined, is pure prior noise to PESTPP-IES. Refused per
    PARAMETER SET -- the whole selection for ``"constant"``, each zone label for
    ``"zone"``, each layer file for ``"grid"``/``"pilotpoints"`` (each a set of its
    own) -- with ``zones=`` labels ``<= 0`` left out of the selection.

    Parameters
    ----------
    project : PestProject
        The calibration being built.
    spec : NativeParameterSpec
        An ``sto.ss``/``sto.sy`` parameter.
    model : flopy.mf6.MFModel
        The staged flow model; ``zones=`` is resolved on its grid.
    packages : list
        Its STO package (:func:`select_packages`).
    targets : list of TargetFile
        The files selected (after ``layers=``).

    Raises
    ------
    ValueError
        As above, naming the layer or zone and the spelling that leaves it out.
    """

    canonical = spec.recipe.canonical
    nper = staged_nper(model)
    steady = steady_periods(_StagedFlowModel(model), nper)
    if steady is not None and all(steady):
        raise ValueError(
            f"Target {canonical!r} scales the STO's {spec.recipe.variable}, and every one of "
            f"the run's {nper} stress period(s) is steady state, where MF6 reads no storage, "
            "so it cannot act. Declare a transient period (mf.sto(transient={1: True}, ...)) "
            "and run the model again, or drop the target."
        )
    if spec.recipe.variable != "sy":
        return
    nlay = int(model.modelgrid.nlay)
    iconvert = None
    for package in packages:
        held = getattr(package, "iconvert", None)
        iconvert = None if held is None else held.array
    iconvert = np.zeros(1, dtype=int) if iconvert is None else np.asarray(iconvert)
    per_layer = np.broadcast_to(iconvert.reshape(-1), (int(model.modelgrid.nnodes),)) \
        if iconvert.size == 1 else iconvert.reshape(-1)
    per_layer = np.asarray(per_layer).reshape(nlay, -1) if per_layer.size % nlay == 0 \
        else per_layer.reshape(1, -1)
    labels = None
    if spec.zones is not None:
        labels = np.asarray(
            resolve_zone_array(spec.zones, family="array", model=model)).reshape(-1)
    covered = []
    for target in targets:
        cells = per_layer.reshape(-1) if target.layer is None else per_layer[target.layer]
        selected = np.ones(cells.size, dtype=bool)
        zone = None
        if labels is not None and labels.size == cells.size:
            selected = labels > 0
            zone = labels
        covered.append((target, cells != 0, selected, zone))

    def refuse(where: str, advice: str) -> None:
        raise ValueError(
            f"Target 'sto.sy' (name {spec.name!r}): every cell of {where} has ICONVERT 0 "
            "(confined), and MF6 reads SY only in a convertible cell, so its parameter "
            f"cannot act. {advice}"
        )

    if spec.style in ("grid", "pilotpoints"):
        for target, convertible, selected, _zone in covered:
            if not (convertible & selected).any():
                refuse(f"layer {target.layer if target.layer is not None else 0} "
                       f"({target.filename}) selected",
                       "Leave it out with layers=, or make the layer convertible.")
    elif spec.style == "zone" and labels is not None:
        for label in sorted({int(value) for value in labels if value > 0}):
            if not any((convertible & (zone == label)).any()
                       for _target, convertible, _selected, zone in covered if zone is not None):
                refuse(f"zone {label} on the layers selected",
                       "Give it label 0 to leave it alone, or select a convertible layer.")
    elif not any((convertible & selected).any() for _t, convertible, selected, _z in covered):
        refuse("the layers and zones selected",
               "Select a convertible layer with layers=, or parameterize sto.ss instead.")


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

    model = flopy_model_for(project, spec.declared_recipe)
    packages = select_packages(model, spec.declared_recipe, spec.package)
    # The package's FORM decides the dataset: an array RCH's per-period `recharge`
    # arrays, a list RCH's records (ledger 547).
    recipe = spec.resolve_form(packages)
    workspace = Path(model.simulation_data.mfpath.get_sim_path())
    targets = package_target_files(model, recipe, packages, workspace)
    if spec.layers is not None:
        targets = select_layers(targets, spec.layers, recipe)
    if recipe.package == "sto":
        refuse_storage_that_cannot_act(project, spec, model, packages, targets)
    spec.period_groups = period_groups_for(project, spec, model, targets)
    spec.period_group_pargps = [[] for _group in spec.period_groups]
    _refuse_one_set_over_records_of_unequal_width(spec, recipe, packages, workspace)
    targets = [target for group in spec.period_groups for target in group.targets]
    if spec.capture and sum(target.period is not None for target in targets) > 1:
        raise NotImplementedError(
            f"capture=True on {recipe.canonical!r} (name {spec.name!r}): it scales "
            f"{len(targets)} per-period files, and a capture records ONE field under one "
            "prefix, so their observation names would collide. Capturing per-period files "
            "is deferred (ledger 547); drop capture=."
        )
    spec.resolved_targets = list(targets)
    spec.resolved_files = [target.filename for target in targets]
    return targets


def _list_record_width(path: Path) -> int | None:
    """The columns of an MF6 external list file's first record, as pyEMU splits it.

    ``None`` for a file holding no record.
    """

    try:
        first = pd.read_csv(path, sep=r"\s+", header=None, dtype=str, keep_default_na=False,
                            comment="#", nrows=1)
    except pd.errors.EmptyDataError:
        return None
    return int(first.shape[1])


def _record_options(package) -> str:
    """What a list package declares that widens its records: BOUNDNAMES, AUXILIARY."""

    said = []
    boundnames = getattr(package, "boundnames", None)
    if boundnames is not None and boundnames.get_data():
        said.append("BOUNDNAMES")
    auxiliary = getattr(package, "auxiliary", None)
    names = auxiliary.get_data() if auxiliary is not None else None
    if names is not None and len(names):
        said.append("AUXILIARY " + " ".join(str(name) for name in list(names[0])[1:]))
    return ", ".join(said) if said else "neither BOUNDNAMES nor AUXILIARY"


def _refuse_one_set_over_records_of_unequal_width(spec: NativeParameterSpec, recipe: _Recipe,
                                                  packages, workspace: Path) -> None:
    """Refuse ONE parameter set over list packages whose records are not one width.

    ``package=("riv", "crisp")`` shares one set, and every ``add_parameters`` call
    hands pyEMU the files of all the packages it names, which pyEMU reads as ONE
    table: when one package declares BOUNDNAMES (or AUXILIARY variables) the other
    does not, its records carry more columns and pyEMU stops with ``shape mismatch
    for array style`` after the whole staging (ledger 559). A ``style="feature"``
    set is one call per package (``package=`` there is one name, ledger 558), so it
    is let through.

    Parameters
    ----------
    spec : NativeParameterSpec
        The parameter request; its ``period_groups`` (one ``add_parameters`` call
        each) and ``name`` are read.
    recipe : _Recipe
        The target's recipe in its package's form. Anything but a ``"list"`` recipe
        is let through.
    packages : list
        The FloPy packages :func:`select_packages` chose; their BOUNDNAMES and
        AUXILIARY options are named in the refusal.
    workspace : Path
        The staged simulation's workspace, which the file names are relative to.

    Raises
    ------
    ValueError
        When the files one call would hand pyEMU hold records of different widths,
        naming each package, its width and its options, and the spelling that
        builds: one ``parameterize`` call per package.
    """

    if recipe.family != "list" or (recipe.feature is not None and spec.style == "feature"):
        return
    options = {package.package_name: package for package in packages}
    for group in spec.period_groups:
        widths: dict[str, int] = {}
        for target in group.targets:
            if target.package in widths:
                continue
            width = _list_record_width(Path(workspace) / target.filename)
            if width is not None:
                widths[target.package] = width
        if len(set(widths.values())) <= 1:
            continue
        said = "; ".join(
            f"{name!r} {width} columns ({_record_options(options[name])})"
            if name in options else f"{name!r} {width} columns"
            for name, width in widths.items()
        )
        calls = " and ".join(f"parameterize({recipe.canonical!r}, package={name!r})"
                             for name in widths)
        raise ValueError(
            f"{recipe.canonical!r} shared over packages {tuple(widths)} ({spec.name!r}): "
            f"their records are not one width -- {said}. One parameter set hands pyEMU "
            "the files of every package it names as ONE table, and pyEMU stops on that "
            "('shape mismatch for array style'). Declare one call per package -- "
            f"{calls} -- or declare BOUNDNAMES and AUXILIARY alike on every package "
            "(ledger 559)."
        )



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
        The target's recipe. Anything but an ``"array"`` recipe is let through unread,
        and so is a per-period one (an array RCH's ``recharge``), whose files
        :func:`package_target_files` finds or refuses.
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

    if recipe.family != "array" or recipe.periodic:
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
        The target's recipe; anything but an ``"array"`` griddata recipe is left
        alone -- a per-period array (an array RCH's ``recharge``) is one value per
        cell of a layer already.
    package : str or tuple of str or None, optional
        The spec's ``package=`` (:func:`select_packages`). Default ``None``.

    Returns
    -------
    bool
        Whether any package's array was re-stored.
    """

    if recipe.family != "array" or recipe.periodic:
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


def _cell_bottoms(model) -> np.ndarray:
    """Every cell's bottom elevation on ``model``'s discretization, indexed by a cellid.

    Read from the discretization package itself -- DIS ``botm`` ``(nlay, nrow, ncol)``,
    DISV ``botm`` ``(nlay, ncpl)``, DISU ``bot`` ``(nodes,)`` -- so a ZERO-based FloPy
    cellid indexes it as written: ``bottoms[(layer, cell)]``, ``bottoms[(node,)]``.

    Parameters
    ----------
    model : flopy.mf6.MFModel
        The build's staged flow model (:func:`staged_model`).

    Returns
    -------
    numpy.ndarray
        The bottoms, float.

    Raises
    ------
    ValueError
        When the model holds no DIS, DISV or DISU package.
    """

    for package in model.packagelist:
        kind = str(getattr(package, "package_type", "")).lower()
        if kind in ("dis", "disv"):
            return np.asarray(package.botm.array, dtype=float)
        if kind == "disu":
            return np.asarray(package.bot.array, dtype=float).reshape(-1)
    raise ValueError(f"Model {model.name!r} holds no DIS, DISV or DISU package to read cell bottoms from.")


def _refuse_a_river_bottom_lowered_through_its_cell(project, spec: NativeParameterSpec, model,
                                                    grid_type: str) -> None:
    """Refuse a ``riv.stage`` lower bound that carries a record's RBOT below its cell's bottom.

    The offset moves each record's ``rbot`` with its stage (ledger 519), so a lower
    bound below a record's clearance -- ``rbot`` minus its cell's bottom -- lowers that
    ``rbot`` through the cell: at ``icelltype != 0`` MF6 refuses the realization by
    name, at ``0`` it over-drains under Normal termination (ledger 290). The floor is
    the cell bottom plus :func:`~myflopy.modflow.mf6.boundary_support.bottom_clearance`,
    since FloPy's writing of the two numbers can round an exact equality below
    (ledger 290). Offsets on one file ADD in the forward run, so the lower bounds of
    every ``riv.stage`` parameter already compiled onto it count too. ``physical=``
    is not credited. A record whose ``rbot`` the run already holds at or below its
    floor is not this check's: no offset puts it there (ledger 568).
    """

    index_cols, _stage_col = list_columns(spec.recipe, grid_type)
    _index, rbot_col = list_columns(spec.recipe, grid_type, field="rbot")
    bottoms = _cell_bottoms(model)
    # Emptied by each `PestProject.build`: the offsets one template's forward run adds.
    floors = project.__dict__.setdefault("_river_offset_lower_bounds", {})
    lower = float(spec.bounds[0])
    offending = 0
    worst = None
    safest = -np.inf
    for target in spec.resolved_targets:
        rows = _read_list_file(Path(project.template_workspace) / target.filename)
        cellids = rows[index_cols].astype(float).astype(int).to_numpy() - 1
        bottom = bottoms[tuple(cellids.T)]
        rbot = pd.to_numeric(rows[rbot_col]).to_numpy(dtype=float)
        floor = bottom + bottom_clearance(bottom)
        earlier = float(floors.get(target.filename, 0.0))
        checked = rbot > floor
        lowest = rbot + earlier + lower
        below = checked & (lowest < floor)
        if checked.any():
            safest = max(safest, float(np.max((floor - rbot - earlier)[checked])))
        offending += int(below.sum())
        if below.any():
            depth = (lowest - floor)[below]
            row = int(np.flatnonzero(below)[np.argmin(depth)])
            if worst is None or depth.min() < worst[0]:
                worst = (float(depth.min()), target, row, tuple(int(value) for value in cellids[row]),
                         float(rbot[row]), float(bottom[row]))
    if offending:
        _depth, target, row, cellid, rbot, bottom = worst
        others = float(floors.get(target.filename, 0.0))
        # Rounded TOWARD zero, so the bound named is safe as printed.
        safe = np.ceil(safest * 1000.0) / 1000.0
        stacked = (f" -- with the {others:g} the riv.stage parameters already on its "
                   "file add" if others else "")
        raise ValueError(
            f"riv.stage {spec.name!r}: its lower bound {lower:g} carries the river bottom of "
            f"{offending} record(s) below their cell's bottom{stacked} -- worst, row "
            f"{row + 1} of {target.filename} (package {target.package!r}), cellid {cellid} "
            f"(zero-based): rbot {rbot:g} over a cell bottom of {bottom:g}. The offset moves "
            "each rbot with its stage, and an rbot below its cell is refused by MF6 at "
            "icelltype != 0 and over-drains at 0 (ledger 290). "
            + (f"The most negative lower bound that keeps every rbot inside its cell is "
               f"{safe:g}: bounds=({safe:g}, {float(spec.bounds[1]):g})."
               if safe < 0 else
               "No negative lower bound keeps every rbot inside its cell: bounds=(0.0, "
               f"{float(spec.bounds[1]):g}) only raises the river.")
            + " (ledger 568)"
        )
    for target in spec.resolved_targets:
        floors[target.filename] = float(floors.get(target.filename, 0.0)) + lower


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


# --- feature targets: one parameter per lake or named stream ---------------------


@dataclass(frozen=True)
class FeatureSet:
    """The features one parameter set of a FEATURE target covers (ledger 558).

    Attributes
    ----------
    label : str or None
        The set's name suffix -- a lake's or stream's boundname, or its number, slugged
        by :func:`~myflopy.modflow.mf6.pest.naming.pest_group_name` -- for
        ``style="feature"``; ``None`` for one set over the selection.
    features : tuple of int
        The ZERO-based lakes or reaches it covers.
    """

    label: str | None
    features: tuple[int, ...]


def _feature_numbers(package) -> list[int]:
    """The ZERO-based feature numbers a LAK or SFR package's ``packagedata`` declares."""

    data = package.packagedata.get_data()
    if data is None:
        return []
    table = data if isinstance(data, pd.DataFrame) else pd.DataFrame(data)
    return [int(value) for value in table.iloc[:, 0]]


def resolve_features(spec: NativeParameterSpec, package) -> list[FeatureSet]:
    """The parameter sets ``spec``'s features make on ``package``, checked against it.

    Read off the package's ``packagedata`` alone -- its features' ZERO-based numbers
    and boundnames (:func:`~myflopy.modflow.mf6.package_registry.feature_boundnames`,
    the one boundname rule observations use too) -- so :meth:`PestProject.parameterize`
    asks it of the caller's model and :meth:`PestProject.build` of the staged one.

    Parameters
    ----------
    spec : NativeParameterSpec
        A FEATURE target's spec; its ``features``, ``style`` and ``name`` are read.
    package : flopy.mf6.ModflowGwflak or flopy.mf6.ModflowGwfsfr
        The ONE package the spec scales.

    Returns
    -------
    list of FeatureSet
        For ``style="feature"``: one per entry of ``features`` -- by default one per
        lake (LAK) or per boundname (SFR, in order of first appearance) -- labelled by
        its boundname, else its number. Otherwise one unlabelled set over the entries
        of ``features``, or over every feature.

    Raises
    ------
    ValueError
        When an entry names no feature of the package (listing those it has), two
        entries name one feature, two labels make one PEST name, or ``style="feature"``
        on an SFR with no boundnames is given no ``features`` (naming
        ``style='grid'``).
    """

    recipe = spec.declared_recipe
    word = recipe.feature
    plural = "lakes" if word == "lake" else "reaches"
    name = package_name_of(package)
    held = _feature_numbers(package)
    data = package.packagedata.get_data()
    named = feature_boundnames(data)
    # Each feature's boundname as the package spells it, for its label and messages;
    # matching is `feature_boundnames`' (any case).
    spelled: dict[int, str] = {}
    if named:
        table = data if isinstance(data, pd.DataFrame) else pd.DataFrame(data)
        for number, bound in zip(table.iloc[:, 0], table["boundname"], strict=True):
            text = "" if bound is None else str(bound).strip()
            if text and text.lower() in named:
                spelled.setdefault(int(number), text)
    own_name = {number: spelled[number] for number in spelled}
    where = f"{recipe.package.upper()} package {name!r}"

    def listing() -> str:
        span = f"{plural} 0 .. {max(held)}" if held else f"no {plural}"
        names = ", ".join(repr(spelled[numbers[0]]) for numbers in named.values()) or "none"
        return f"it has {span} (numbered from 0); boundnames: {names}"

    if spec.features is None:
        if spec.style != "feature":
            return [FeatureSet(None, tuple(held))]
        if word == "lake":
            entries = list(held)
        elif named:
            entries = [spelled[numbers[0]] for numbers in named.values()]
        else:
            raise ValueError(
                f"style='feature' on {recipe.canonical!r} makes one parameter per NAMED "
                f"stream, and {where} declares no BOUNDNAMES, so it names none. Use "
                "style='grid' (one parameter per reach), or name the reaches in the "
                "package (boundnames) or with features=."
            )
    else:
        entries = list(spec.features)

    resolved: list[tuple[str, tuple[int, ...]]] = []
    wrong: list[str] = []
    for entry in entries:
        if isinstance(entry, str):
            numbers = named.get(entry.lower(), [])
            if not numbers:
                wrong.append(f"boundname {entry!r}")
                continue
            label = entry
        else:
            if entry not in held:
                hint = (f" (one past the last: MF6's files number {plural} from 1)"
                        if held and entry == max(held) + 1 else "")
                wrong.append(f"{word} {entry}{hint}")
                continue
            numbers = [entry]
            label = own_name.get(entry, str(entry)) if word == "lake" else str(entry)
        resolved.append((label, tuple(numbers)))
    if wrong:
        raise ValueError(
            f"features= on {recipe.canonical!r} names {', '.join(wrong)}, which "
            f"{where} does not have -- {listing()}. Give each {word} its ZERO-based "
            "number or its boundname (any case)."
        )
    seen: dict[int, str] = {}
    for label, numbers in resolved:
        for number in numbers:
            if number in seen:
                raise ValueError(
                    f"features= on {recipe.canonical!r} names {word} {number} twice -- "
                    f"through {seen[number]!r} and {label!r} -- and a row takes one "
                    "multiplier from one call. Name it once."
                )
            seen[number] = label
    if spec.style != "feature":
        return [FeatureSet(None, tuple(number for _label, numbers in resolved
                                       for number in numbers))]
    slugs: dict[str, str] = {}
    sets = []
    for label, numbers in resolved:
        slug = pest_group_name(label)
        if not re.search(r"[a-z0-9]", slug):
            raise ValueError(
                f"style='feature' on {recipe.canonical!r}: {label!r} holds no letter or "
                "digit, so it makes no PEST name. Name the feature by its number."
            )
        if slug in slugs:
            raise ValueError(
                f"style='feature' on {recipe.canonical!r}: {slugs[slug]!r} and {label!r} "
                f"both make the PEST name {spec.name}_{slug!s}, so their parameters "
                "would be one. Rename one, or name them by number in features=."
            )
        slugs[slug] = label
        sets.append(FeatureSet(slug, numbers))
    return sets


def feature_columns(recipe: _Recipe, grid_type: str, style: str) -> tuple[list[int], int]:
    """pyEMU's ``index_cols`` and the value column of a FEATURE target, by style.

    ``"zone"`` keys a row by its CELL, which pyEMU looks its zone up at; ``"grid"`` by
    its feature key AND its cell, so each reach is one parameter and pyEMU locates it
    at the cell's centre (its last two index columns, or one on DISU); every other
    style by the feature key alone (LAK ``lakeno iconn``, SFR ``ifno``), unique in the
    file (ledger 558).

    Parameters
    ----------
    recipe : _Recipe
        A FEATURE target's recipe.
    grid_type : str
        FloPy's ``modelgrid.grid_type``, as :func:`list_columns` reads it.
    style : str
        The spec's style.

    Returns
    -------
    tuple of (list of int, int)
        The zero-based index columns and the value column, on that grid.
    """

    cellid, use_col = list_columns(recipe, grid_type)
    if style == "zone":
        return cellid, use_col
    if style == "grid":
        return list(recipe.key_cols) + cellid, use_col
    return list(recipe.key_cols), use_col


def _is_none(text) -> bool:
    """Whether a list file's field is MF6's ``NONE``."""

    return str(text).strip().lower() == "none"


def _feature_rows(project, spec: NativeParameterSpec, target: TargetFile, grid_type: str,
                  index_cols: list[int], value_col: int) -> pd.DataFrame:
    """One row per record of a FEATURE target's file: its feature, pyEMU key and value.

    Read from the file pyEMU reads (text, 1-based as written), refusing what pyEMU
    would mis-read: an SFR reach whose cellid is ``NONE`` is one column short, and
    pyEMU, reading by position, multiplied its Manning's n (measured).
    """

    recipe = spec.recipe
    rows = _read_list_file(Path(project.template_workspace) / target.filename)
    if recipe.feature == "reach":
        first_cell = list_columns(recipe, grid_type)[0][0]
        unplaced = [int(row[0]) - 1 for _index, row in rows.iterrows()
                    if _is_none(row[first_cell])]
        if unplaced:
            raise ValueError(
                f"{recipe.canonical!r} on SFR package {target.package!r}: reach(es) "
                f"{unplaced[:10]} (zero-based) have cellid NONE, so their rows in "
                f"{target.filename} are one column shorter than the rest, and pyEMU reads "
                "columns by position -- it multiplied such a reach's Manning's n in place "
                "of its rhk (measured). Give each reach its cell, or take the unconnected "
                "reaches out of the package (ledger 558)."
            )
    cellid = list_columns(recipe, grid_type)[0]
    keys = [tuple(int(float(value)) for value in row[index_cols]) for _index, row in rows.iterrows()]
    return pd.DataFrame({
        "feature": [int(float(value)) - 1 for value in rows[recipe.key_cols[0]]],
        "connection": ([int(float(value)) - 1 for value in rows[recipe.key_cols[1]]]
                       if len(recipe.key_cols) > 1 else None),
        "key": keys,
        # The cell's place in the map plane: the cellid without its layer.
        "column": [tuple(str(value) for value in row[cellid[1:] or cellid])
                   for _index, row in rows.iterrows()],
        "value": rows[value_col].to_numpy(),
    })


def _record_points(modelgrid, keys: Sequence[tuple[int, ...]]) -> np.ndarray:
    """Where pyEMU places each list record of a ``"grid"`` parameter: its cell's centre.

    ``keys`` are the records' pyEMU index columns as the file writes them, ONE-based,
    their cellid last -- pyEMU reads a record's ``(i, j)`` from its last two (DIS: row,
    column), its cell from its last one (DISV) or its node (DISU, through
    :func:`pyemu_spatial_reference`'s dict). Returns an ``(n, 2)`` array of x, y.
    """

    x = np.asarray(modelgrid.xcellcenters, dtype=float)
    y = np.asarray(modelgrid.ycellcenters, dtype=float)
    if modelgrid.grid_type == "structured":
        rows = np.array([key[-2] - 1 for key in keys], dtype=int)
        cols = np.array([key[-1] - 1 for key in keys], dtype=int)
        return np.c_[x[rows, cols], y[rows, cols]]
    cells = np.array([key[-1] - 1 for key in keys], dtype=int)
    return np.c_[x.reshape(-1)[cells], y.reshape(-1)[cells]]


def _refuse_records_pyemu_cannot_correlate(spec: NativeParameterSpec, modelgrid,
                                          keys: Sequence[tuple[int, ...]], labels: Sequence,
                                          *, plural: str, advice: str) -> None:
    """Refuse ``correlation=`` on a LIST ``"grid"`` set whose records pyEMU cannot correlate.

    pyEMU locates each record at its cell's centre (:func:`_record_points`) and draws
    the prior from one covariance per covariance group. Where two records sit at one
    point it splits the group by the index column before the cell's -- a DIS
    cellid's layer, an SFR reach's number -- and only when it HAS one: a DISV or
    DISU cellid has none. Mirrored here, two outcomes are refused at build (ledger
    570; ledger 558 refused the second for ``sfr.rhk`` alone):

    * two parameters at one point in ONE group -- a column's records in two layers
      on DISV/DISU, or two reaches in one DIS cell and layer (records sharing a whole
      cellid share one parameter, pyEMU's rule, so they are one point, not two):
      their correlation is 1, the covariance is singular, and pyEMU's draw at :meth:`PestProject.run_ies`/:meth:`PestProject.prior`
      died with ``LinAlgError: Matrix is not positive definite`` (measured on DISV and
      DISU). A ``nugget`` makes it drawable -- measured from ``1e-10`` up -- but
      correlates the records of one column at ``1 / (1 + nugget)``, a vertical
      structure the map-plane variogram never described, so this is refused WHATEVER
      the nugget, as :func:`_refuse_a_singular_grid_prior` refuses a DISU array
      (ledger 452; review of ledger 570);
    * a split that leaves every record a group of its own -- every reach of a
      ``sfr.rhk`` set on DISV: the prior came out uncorrelated, off-diagonal within
      0.03 of 0 against 0.66-0.83 without the shared cell (measured, ledger 558).

    A DIS split into groups of several records each is let through: pyEMU correlates
    the records of each LAYER and not across layers (measured: 0.44-0.74 within a
    layer, -0.05 to 0.00 across; ``docs/reproducers/correlated_list_grid_shared_points.py``).

    Parameters
    ----------
    spec : NativeParameterSpec
        The parameter; nothing is asked unless its style is ``"grid"`` with
        ``correlation=``.
    modelgrid : flopy.discretization.Grid
        The staged model's grid.
    keys : sequence of tuple of int
        One per PARAMETER, its pyEMU index columns as written in the file (one-based),
        the cellid last.
    labels : sequence
        How each record is named in the refusal, in ``keys``' order.
    plural : str
        The records' noun, plural (``"records"``, ``"reaches"``).
    advice : str
        The spellings the refusal names.

    Raises
    ------
    ValueError
        In either case above.
    """

    if spec.style != "grid" or spec.correlation is None or len(keys) < 2:
        return
    points = _record_points(modelgrid, keys)
    frame = pd.DataFrame({"x": points[:, 0], "y": points[:, 1], "label": list(labels)})
    if not frame.duplicated(["x", "y"]).any():
        return
    split = len(keys[0]) > 2
    frame["group"] = [key[-3] for key in keys] if split else 0
    canonical = spec.recipe.canonical
    stacked = frame[frame.duplicated(["group", "x", "y"], keep=False)]
    if not stacked.empty:
        first = stacked.iloc[0]
        together = list(stacked.loc[(stacked["x"] == first["x"]) & (stacked["y"] == first["y"])
                                    & (stacked["group"] == first["group"]), "label"])
        why = (
            "so their correlation is 1, the prior covariance is singular, and pyEMU's "
            "prior draw at run_ies()/prior() fails with 'Matrix is not positive definite' "
            "(measured)"
            if float(spec.nugget) == 0.0 else
            f"so nugget={float(spec.nugget):g} correlates them at "
            f"{1.0 / (1.0 + float(spec.nugget)):.3g} -- drawable, but a structure along the "
            "column that the map-plane variogram never described, refused whatever the "
            "nugget as a DISU array's stacked nodes are (ledger 452)"
        )
        raise ValueError(
            f"correlation= on {canonical!r} with style='grid': {plural} {together} "
            f"(zero-based) lie at one point -- one model column"
            + (", one layer" if split else "")
            + f" -- in one covariance group ({len(stacked)} selected {plural} share a point "
            f"with another), {why}. {advice} (ledger 570)."
        )
    if split and frame["group"].is_unique:
        shared = frame[frame.duplicated(["x", "y"], keep=False)]
        first = shared.iloc[0]
        together = list(shared.loc[(shared["x"] == first["x"]) & (shared["y"] == first["y"]),
                                   "label"])
        raise ValueError(
            f"correlation= on {canonical!r} with style='grid': {plural} {together} "
            f"(zero-based) lie in one model column ({len(shared)} selected {plural} share a "
            "column with another), so pyEMU places them at one point and then draws EVERY "
            "one's prior uncorrelated -- correlation= would not act (measured). "
            f"{advice} (ledger 558, ledger 570)."
        )


def _refuse_an_uncorrelatable_reach_grid(spec: NativeParameterSpec, grid_type: str,
                                         chosen: pd.DataFrame, modelgrid) -> None:
    """Refuse ``correlation=`` on ``sfr.rhk`` ``"grid"`` where pyEMU cannot correlate reaches.

    On DISU pyEMU cannot locate a reach at all: a reach's key is two columns, and
    pyEMU looks a node up by one. Elsewhere :func:`_refuse_records_pyemu_cannot_correlate`
    asks what it asks of every list target: two selected reaches in one model COLUMN
    sit at one point, and pyEMU then puts every reach on DISV in a covariance group of
    its own (``covgp`` ``<name>_cov<ifno>``), drawing the prior uncorrelated (measured;
    ledger 558).
    """

    if spec.style != "grid" or spec.correlation is None:
        return
    if grid_type == "unstructured":
        raise ValueError(
            f"correlation= on {spec.recipe.canonical!r} with style='grid' on a DISU model: "
            "pyEMU locates a DISU record by its node alone, and a reach is keyed by its "
            "number and node, so it cannot be placed. Drop correlation=, or use "
            "style='feature'/'zone'/'constant'."
        )
    _refuse_records_pyemu_cannot_correlate(
        spec, modelgrid, list(chosen["key"]), [int(reach) for reach in chosen["feature"]],
        plural="reaches",
        advice=("Leave all but one of each such set out with features=, drop correlation=, "
                "or use style='feature' (one per named stream) or 'zone'"),
    )


def _refuse_rows_pyemu_did_not_select(project, frame: pd.DataFrame, expected: set,
                                      spec: NativeParameterSpec) -> None:
    """Refuse a template that scales rows other than ``expected`` (ZERO-based keys).

    pyEMU sets up parameters for EVERY row when no ``use_rows`` matches (measured), so
    a key mismatch would scale every lake in silence; the template it wrote says which
    rows it took.
    """

    for tpl in sorted({str(name) for name in frame["tpl_filename"]}):
        path = Path(tpl) if Path(tpl).is_absolute() else Path(project.template_workspace) / tpl
        written = pd.read_csv(path, skiprows=1, usecols=["sidx"])["sidx"]
        took = {tuple(int(part) for part in str(value).strip("()").split(",") if part.strip())
                for value in written}
        if took != expected:
            raise RuntimeError(
                f"pyEMU's template {path.name} for {spec.recipe.canonical!r} (name "
                f"{spec.name!r}) scales {len(took)} row(s) where {len(expected)} were "
                "selected; refusing a parameter that would scale other rows."
            )


def _add_feature_parameters(project, spec: NativeParameterSpec, targets: list[TargetFile],
                            grid_type: str, kwargs: dict) -> pd.DataFrame:
    """Compile a FEATURE target (``lak.bedleak``, ``sfr.rhk``) into ``add_parameters`` calls.

    One call per :class:`FeatureSet` -- one per lake or named stream for
    ``style="feature"``, one over the selection otherwise -- each with ``use_rows``
    computed from the staged package, never left to pyEMU's all-rows fallback, and
    the template checked to scale exactly those rows (ledger 558).
    """

    recipe = spec.recipe
    model = flopy_model_for(project, spec.declared_recipe)
    index_cols, value_col = feature_columns(recipe, grid_type, spec.style)
    kwargs = dict(kwargs)
    kwargs["index_cols"] = index_cols
    kwargs["use_cols"] = [value_col]
    tables: dict[str, pd.DataFrame] = {}
    sets_by_package: dict[str, list[FeatureSet]] = {}
    for target in targets:
        tables[target.filename] = _feature_rows(project, spec, target, grid_type,
                                                index_cols, value_col)
        if target.package not in sets_by_package:
            package = package_named(model, recipe.package, target.package,
                                    everyone=packages_of_type(model, recipe.package))
            sets_by_package[target.package] = resolve_features(spec, package)
    plural = "lakes" if recipe.feature == "lake" else "reaches"
    calls: list[tuple[str, FeatureSet, list[TargetFile]]] = []
    for target in targets:
        for feature_set in sets_by_package[target.package]:
            calls.append((target.package, feature_set, [target]))
    if spec.style != "feature" and len(targets) > 1:
        # One set shared across the packages a tuple names: ONE call over every file.
        calls = [(", ".join(target.package for target in targets),
                  sets_by_package[targets[0].package][0], list(targets))]
    if spec.style == "feature" and spec.features is None and recipe.feature == "reach":
        for target in targets:
            named = {number for item in sets_by_package[target.package]
                     for number in item.features}
            unnamed = sorted(set(tables[target.filename]["feature"]) - named)
            if unnamed:
                message = (
                    f"style='feature' on {recipe.canonical!r} ({spec.name!r}): {len(unnamed)} "
                    f"reach(es) of SFR package {target.package!r} carry no boundname -- "
                    f"{unnamed[:10]} (zero-based) -- so no named stream covers them and "
                    "they keep their rhk. Name them in features=, or stack a "
                    "style='constant' parameter over every reach."
                )
                logger.warning(message)
                warnings.warn(message, stacklevel=2)
    frames = []
    pargps = spec.period_group_pargps[0] if spec.period_group_pargps else []
    spec.feature_sets = []
    for package, feature_set, files in calls:
        use_rows: set = set()
        count = 0
        for target in files:
            table = tables[target.filename]
            covered = set(feature_set.features) if len(files) == 1 else set(table["feature"])
            chosen = table[table["feature"].isin(covered)]
            if recipe.feature == "lake":
                blank = chosen[[_is_none(value) for value in chosen["value"]]]
                if not blank.empty:
                    first = blank.iloc[0]
                    raise ValueError(
                        f"{recipe.canonical!r} on LAK package {target.package!r}: lake "
                        f"{int(first['feature'])} connection {int(first['connection'])} "
                        f"(both zero-based) has bedleak NONE -- {len(blank)} connection(s) of "
                        "the lakes selected do -- which MF6 reads as no lakebed (the "
                        "conductance is the aquifer's alone), so there is nothing to "
                        "multiply. Give those connections a bedleak, or leave their lakes "
                        "out with features= (ledger 558)."
                    )
            _refuse_an_uncorrelatable_reach_grid(spec, grid_type, chosen, model.modelgrid)
            if chosen.empty:
                raise ValueError(
                    f"{recipe.canonical!r}: the {plural} {list(feature_set.features)[:10]} of "
                    f"package {target.package!r} hold no rows in {target.filename}, so "
                    "the parameter would scale nothing."
                )
            use_rows.update(chosen["key"])
            count += len(chosen)
        call_kwargs = dict(kwargs)
        base = spec.name if feature_set.label is None else f"{spec.name}_{feature_set.label}"
        call_kwargs["par_name_base"] = base
        call_kwargs["pargp"] = base
        call_kwargs["use_rows"] = sorted(use_rows)
        call_kwargs.update(spec.extra)
        filenames = [target.filename for target in files]
        frame = project.pf.add_parameters(
            filenames if len(filenames) > 1 else filenames[0], **call_kwargs)
        # pyEMU holds a key zero-based: one less than the file's when the project's
        # files count from one (`PestProject(zero_based=False)`, the default).
        shift = 0 if getattr(project.pf, "zero_based", False) else 1
        expected = {tuple(value - shift for value in key) for key in use_rows}
        _refuse_rows_pyemu_did_not_select(project, frame, expected, spec)
        frames.append(frame)
        _record_pargps(pargps, frame)
        spec.feature_sets.append({
            "pargp": base,
            "label": feature_set.label,
            "package": package,
            "features": [int(number) for number in feature_set.features],
            "rows": int(count),
        })
    return pd.concat(frames, ignore_index=True)


def _refuse_list_records_pyemu_cannot_correlate(project, spec: NativeParameterSpec, modelgrid,
                                                targets: list[TargetFile],
                                                index_cols: list[int]) -> None:
    """:func:`_refuse_records_pyemu_cannot_correlate` on one ``add_parameters`` call's files.

    pyEMU makes one ``"grid"`` parameter per distinct index key across the call's
    files, so the keys are those, in the order the files hold them.
    """

    if spec.style != "grid" or spec.correlation is None:
        return
    keys: dict[tuple[int, ...], None] = {}
    for target in targets:
        rows = _read_list_file(Path(project.template_workspace) / target.filename)
        for values in rows[index_cols].itertuples(index=False):
            keys.setdefault(tuple(int(float(value)) for value in values), None)
    ordered = list(keys)
    _refuse_records_pyemu_cannot_correlate(
        spec, modelgrid, ordered, [tuple(value - 1 for value in key) for key in ordered],
        plural="records",
        advice=("Drop correlation= (one uncorrelated parameter per record), or use "
                "style='zone' or 'constant'"),
    )


def _record_pargps(pargps: list[str], frame) -> None:
    """Append the PEST groups of one ``add_parameters`` frame to a period group's list."""

    for group in frame["pargp"]:
        if str(group) not in pargps:
            pargps.append(str(group))


def add_native_parameter(project, spec: NativeParameterSpec):
    """Compile one :class:`NativeParameterSpec` into ``pf.add_parameters``.

    Returns the parameter dataframe pyEMU produced, and records it on the
    project for later inspection via :meth:`PestProject.settings` and the build's
    metadata.
    """

    targets = target_files(project, spec)
    recipe = spec.recipe
    files = [target.filename for target in targets]

    staged = flopy_model_for(project, recipe)
    grid = staged.modelgrid
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
        _refuse_a_river_bottom_lowered_through_its_cell(project, spec, staged, grid.grid_type)

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
        # On the STAGED model's grid, like everything else the build reads (ledger 356).
        kwargs["zone_array"] = resolve_zone_array(spec.zones, family=recipe.family, model=staged)
        spec.resolved_zone_labels = sorted(
            int(label) for label in np.unique(np.asarray(kwargs["zone_array"]))
        )
    if recipe.family == "list":
        index_cols, use_col = list_columns(recipe, grid.grid_type)
        kwargs["index_cols"] = index_cols
        kwargs["use_cols"] = [use_col]
    if spec.correlation is not None and spec.style in _SPATIAL_STYLES:
        kwargs["geostruct"] = project._geostruct_for(spec)
    if recipe.feature is not None:
        # `feature` is myflopy's style: each lake or named stream is one pyEMU constant.
        kwargs["par_type"] = "constant" if spec.style == "feature" else spec.style
        frame = _add_feature_parameters(project, spec, targets, grid.grid_type, kwargs)
        _refuse_a_large_correlated_grid(spec, frame)
        project._native_parameter_frames[spec.name] = frame
        return frame
    kwargs.update(spec.extra)

    # One `add_parameters` call per `periods=` group, with that group's files, named
    # `<name>_<label>` and its own PEST group (ledger 549); one call
    # over every file without `periods`. Within a call the file list shares one
    # parameter set (recharge across its period files, a constant K multiplier
    # across all layers) -- except grid ARRAY parameters on LAYER files, which get an
    # independent set *per layer file*, so each layer's K field varies on its own,
    # named by its ZERO-BASED layer as pilot points are (ledger 521). Period files
    # carry no layer and share one set: pyEMU's own rule, one per-cell set across the
    # files of a call (ledger 547).
    frames = []
    for group, pargps in zip(spec.period_groups, spec.period_group_pargps, strict=True):
        name = f"{spec.name}{group.name_suffix}"
        group_files = [target.filename for target in group.targets]
        layers = sorted({target.layer for target in group.targets if target.layer is not None})
        if recipe.family == "array" and spec.style == "grid" and len(layers) > 1:
            for target in group.targets:
                layer_kwargs = dict(kwargs)
                layer_kwargs["par_name_base"] = f"{name}l{target.layer}"
                layer_kwargs["pargp"] = f"{name}l{target.layer}"
                frames.append(project.pf.add_parameters(target.filename, **layer_kwargs))
                _record_pargps(pargps, frames[-1])
            continue
        if recipe.family == "list":
            _refuse_list_records_pyemu_cannot_correlate(project, spec, staged.modelgrid,
                                                        group.targets, kwargs["index_cols"])
        group_kwargs = dict(kwargs)
        group_kwargs["par_name_base"] = name
        group_kwargs["pargp"] = name
        group_kwargs.update(spec.extra)
        frames.append(project.pf.add_parameters(
            group_files if len(group_files) > 1 else group_files[0], **group_kwargs))
        _record_pargps(pargps, frames[-1])
    frame = pd.concat(frames, ignore_index=True)
    _refuse_a_large_correlated_grid(spec, frame)
    if recipe.canonical == "riv.stage":
        _carry_river_bottom(project, spec, grid.grid_type)
    project._native_parameter_frames[spec.name] = frame
    return frame

"""Shared normalization/series helpers for observation targets."""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd

from myflopy._logging import get_logger
from myflopy.modflow.mf6.grid.geometry import where_they_are
from myflopy.modflow.mf6.package_registry import agg_for_kind
from myflopy.modflow.mf6.period_end import (
    at_period_end,
    period_end_rows,
    period_end_steps,
    period_ends,
    period_lengths,
)

logger = get_logger(__name__)

#: How far apart two copies of ONE value may sit and still be the same number,
#: relative to the value (floored at 1). The rows a feature-valued table repeats
#: for each connected cell are bit-identical in every table measured, so this only
#: absorbs a rounding a future reader might introduce (ledger 331).
_COPY_RTOL = 1e-9


def _copy_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Return a shallow copy with a clean integer index."""

    return frame.copy().reset_index(drop=True)


def _normalize_identifier_series(values: pd.Series) -> pd.Series:
    """Normalize identifiers for robust joins and observation naming.

    Parameters
    ----------
    values
        Series of location or column identifiers.

    Returns
    -------
    pandas.Series
        Lower-case, stripped string identifiers.
    """

    return values.astype(str).str.strip().str.lower()


def _shown(value) -> str:
    """A value for a message as the user wrote it -- ``1.5``, not ``np.float64(1.5)``."""

    if isinstance(value, np.generic):
        value = value.item()
    return repr(value)


def _is_blank(value) -> bool:
    """True for an entry that means "not measured": missing, or empty or blank text."""

    if isinstance(value, str):
        return not value.strip()
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        # ``pd.isna`` of a list or an array answers elementwise, which ``bool``
        # refuses; such an entry is not a blank, and the numeric conversion that
        # follows refuses it by name.
        return False


def _coerce_target_values(
    long: pd.DataFrame,
    column: str,
    *,
    owner: str,
    source: pd.Series | str,
) -> pd.Series:
    """A target table's measured values as floats: a BLANK is unmeasured, anything else a number.

    ``pd.to_numeric(errors="coerce")`` read both the same way -- a blank cell and a
    typo like ``"12,5"`` or ``"dry"`` each became ``NaN`` -- so a value that failed
    to parse was quietly unmeasured (ledger 340). The two are told apart BEFORE the
    conversion: an entry that is missing, or text that is empty or blank, is
    ``NaN``; one that parses as a FINITE number is that number; anything else is
    refused -- an infinite value included, which reached the control file as an
    ``obsval`` pyEMU refuses as a NaN, writing ``pst.obs_data.nans.csv`` into the
    caller's working directory (ledger 341).

    Parameters
    ----------
    long : pandas.DataFrame
        The normalized long table: ``column`` plus, for the message, ``name`` and
        ``time`` where it has them.
    column : str
        The column of ``long`` holding the raw values.
    owner : str
        The target class, named in the error (``"HeadTargets"``).
    source : pandas.Series or str
        Where each value came from in the caller's own table -- the value column
        of a LONG table, or the target's own column of a WIDE one -- as a Series
        aligned with ``long``, or one name for every row.

    Returns
    -------
    pandas.Series
        ``column`` as float, ``NaN`` exactly where the entry was blank.

    Raises
    ------
    ValueError
        When any entry is neither blank nor a finite number, naming up to five by
        target, time, the entry and its column.
    """

    raw = long[column]
    blank = raw.map(_is_blank).astype(bool)
    numbers = pd.to_numeric(raw.where(~blank), errors="coerce").astype(float)
    failed = pd.Series(False, index=long.index)
    for index in (~blank & numbers.isna())[lambda held: held].index:
        # What pandas would not read, Python may (``"1_000"``); text that says
        # "nan" is a number that is missing, which is a blank written out.
        try:
            numbers.at[index] = float(str(raw.at[index]).strip())
        except ValueError:
            failed.at[index] = True
    failed |= np.isinf(numbers.to_numpy())
    if failed.any():
        sources = source if isinstance(source, pd.Series) else pd.Series(source, index=long.index)
        offenders = [
            f"{_shown(long.at[index, 'name'])} at time {_shown(long.at[index, 'time'])} "
            f"holds {_shown(raw.at[index])} (column {_shown(sources.at[index])})"
            if {"name", "time"} <= set(long.columns)
            else f"{_shown(raw.at[index])} (column {_shown(sources.at[index])})"
            for index in failed[failed].index[:5]
        ]
        more = int(failed.sum()) - len(offenders)
        raise ValueError(
            f"{owner} values must be finite numbers, or blank where a target was "
            f"not measured, and {int(failed.sum())} are neither: "
            + "; ".join(offenders)
            + (f"; and {more} more" if more else "")
            + ". A blank is read as unmeasured; an entry that is not a number is "
            "not read as one, because a typo would then silently drop the "
            "measurement from PEST."
        )
    return numbers


def _refuse_flows_leaving_the_frame(
    long: pd.DataFrame,
    column: str,
    *,
    owner: str,
    quantity: str,
    leaving: str,
    measured: str,
) -> None:
    """Refuse a flow target of the sign MF6 can never write for it (ledger 350).

    A reach's streamflow (MF6's ``DOWNSTREAM-FLOW``) and a drain's discharge (MF6's
    ``DRN``) are both water LEAVING something -- the reach, the aquifer -- and MF6
    writes that NEGATIVE, from the reach's or the aquifer's point of view; the
    targets compare in that sign on both routes, the user's rule. So a value above
    zero can never be matched and is a sign error: a positive gauged ``Q`` scored
    a residual of ``-2Q`` at a perfect fit (measured, phi 327,106 on three exact
    SFR observations), a misfit only a smaller simulated flow lowers, so it pulls a
    calibration toward drying the reach. Zero is accepted --
    a reach or drain measured dry -- and a blank (``NaN``) is unmeasured.

    Parameters
    ----------
    long : pandas.DataFrame
        The normalized long target table: ``column`` plus, for the message,
        ``name`` and ``time``.
    column : str
        The column of ``long`` holding the observed flows.
    owner : str
        The target class, named in the error (``"SfrFlowTargets"``).
    quantity : str
        What the simulated flow is, for the message (``"MF6's DOWNSTREAM-FLOW"``).
    leaving : str
        What the water leaves, for the message (``"the reach"``).
    measured : str
        What the user measured, for the spelling (``"a gauged streamflow"``).

    Raises
    ------
    ValueError
        When any value is above zero, naming up to five by target, time and value,
        the frame, and the spelling (``-Q``).
    """

    values = pd.to_numeric(long[column], errors="coerce")
    positive = values > 0.0
    if not positive.any():
        return
    offenders = [
        f"{_shown(long.at[index, 'name'])} at time {_shown(long.at[index, 'time'])} "
        f"holds {_shown(long.at[index, column])}"
        if {"name", "time"} <= set(long.columns)
        else _shown(long.at[index, column])
        for index in positive[positive].index[:5]
    ]
    more = int(positive.sum()) - len(offenders)
    raise ValueError(
        f"{owner} values are {quantity}, in MF6's own sign: MF6 writes water leaving "
        f"{leaving} as negative, so the flow is never above 0, and "
        f"{int(positive.sum())} value(s) are: "
        + "; ".join(offenders)
        + (f"; and {more} more" if more else "")
        + f". Enter {measured} Q as -Q (0, a dry reading, is accepted)."
    )


def _observation_label(value: str) -> str:
    """Return a compact, MF6-friendly observation label fragment."""

    normalized = "".join(
        character if str(character).isalnum() else "_"
        for character in str(value).strip()
    )
    normalized = normalized.strip("_")
    return normalized or "obs"


def _clean_observation_output_values(values: pd.Series, *, no_value: float = 0.0) -> pd.Series:
    """Normalize MF6 observation CSV values into ordinary numeric series.

    MF6 writes an observation it cannot take as a huge number -- ``DNODATA``,
    3e30, as its STAGE of an inactive lake. For a FLOW that means no
    contribution, so the default ``no_value=0.0`` stands; an ELEVATION such as a
    stage has no such reading, and a stage caller passes ``NaN`` (ledger 329: an
    inactive lake read as a stage of 0.0 ft; ledger 344: the SFR stage reader
    does too, although MF6 6.8 writes an INACTIVE reach's stage as its bed top,
    measured, not 3e30).
    """

    numeric = pd.to_numeric(values, errors="coerce")
    numeric = numeric.mask(numeric.abs() > 1.0e20, no_value)
    return numeric.astype(float)


#: How many characters of an observation name MF6 keeps (its ``LENOBSNAME``): a
#: longer name is cut, without a word, and the CSV header carries the cut --
#: measured, MF6 6.8: two 42-character names sharing 40 characters wrote two
#: columns under one header, each holding its own record's values (ledger 344).
_MF6_OBSNAME_LENGTH = 40

#: How far apart two of one package's observation files may put one time step's
#: time and still be on MF6's one clock, relative to the time. MF6 writes a TEXT
#: file's ``time`` column at 13 significant digits (``3.747252747253``,
#: ``0.1388889000000E-01`` -- measured, MF6 6.8, whatever the OBS ``DIGITS``), the
#: values beside it with every digit unless ``DIGITS`` says fewer, and a BINARY
#: file's time in full, so the same step read 3.333333333333 in
#: one and 3.3333333333333335 in the other: a half unit in the 13th digit is at most
#: 5e-13 of the time, and this is twice it. Exact equality refused a zone declared
#: over one file of each (ledger 347).
_OBSERVATION_TIME_RTOL = 1e-12


@dataclass(frozen=True)
class _ObservedQuantity:
    """What one target family reads on the observation route, and from which package type.

    Attributes
    ----------
    owner : str
        The target class, named in every message (``"LakeStageTargets"``).
    package_type : str
        The FloPy package type whose OBS declares the observations (``"lak"``).
    obstype : str
        MF6's observation type for the family's quantity, lower-cased
        (``"stage"``, ``"downstream-flow"``, ``"drn"``).
    feature : str
        What a target's id names, for messages: ``"lake"``, ``"reach"``, ``"cell"``.
    file_route : str
        What the family reads when no target is declared, for messages.
    file_route_reads_every_package : bool
        Whether that file route serves a model carrying several packages of the
        type: the DRN budget route sums every drain in a zone's cells whichever
        package holds it; the LAK and SFR files belong to ONE package.
    """

    owner: str
    package_type: str
    obstype: str
    feature: str
    file_route: str
    file_route_reads_every_package: bool = False

    @property
    def cellid(self) -> bool:
        """True when a target's id is a ``(layer, cell)`` pair rather than a feature number."""

        return self.feature == "cell"


@dataclass(frozen=True)
class _DeclaredObservation:
    """One continuous record a package's OBS declares: where MF6 writes it and what it measures.

    Attributes
    ----------
    package : str
        The owning package's name (FloPy's ``package_name``).
    csv : str
        The output file as declared, relative to the model's workspace.
    column : int
        The record's position among that file's records, from 0: MF6 writes a
        file's records as its columns, in declaration order, after ``time``.
    width : int
        How many records the file declares -- its value columns.
    obsname : str
        The observation's name as declared.
    obstype : str
        Its MF6 observation type, lower-cased.
    ident : object
        Its id as FloPy holds it: a number from one or a boundname (LAK, SFR), a
        zero-based cellid tuple (DRN).
    features : tuple
        What the id names in the package: zero-based lake or reach numbers, or
        ``(layer, cell)`` pairs. Empty when it names nothing (a boundname no
        feature carries, or a DRN boundname, which names drains, not a cell).
    """

    package: str
    csv: str
    column: int
    width: int
    obsname: str
    obstype: str
    ident: Any
    features: tuple

    def record(self) -> str:
        """The record as its OBS file reads, and where: `` `drying STAGE 1` in a.csv (package 'lak')``."""

        ident = self.ident
        if isinstance(ident, (tuple, list, np.ndarray)):
            ident = "(" + ", ".join(str(int(value)) for value in ident) + ")"
        return f"`{self.obsname} {self.obstype.upper()} {ident}` in {self.csv} (package {self.package!r})"


def _packages_of_type(model, package_type: str) -> list:
    """Every FloPy package of ``package_type`` in ``model``; ``[]`` when it has none.

    Asked by TYPE: FloPy's ``get_package(name)`` answers a package NAMED like the
    type before the type itself, so ``model.package("lak")`` on a model holding
    LAK packages named ``lak`` and ``lak_2`` returns the first alone. The lookup
    through ``model.package`` comes first because a reloaded run loads its
    packages there; a model that has none of the type -- or a test double with no
    ``package`` -- raises ``AttributeError``, and has nothing declared.
    """

    try:
        first = model.package(package_type)
    except AttributeError:
        return []
    getter = getattr(getattr(model, "gwf", None), "get_package", None)
    found = getter(package_type, type_only=True) if getter is not None else first
    if found is None:
        return []
    return list(found) if isinstance(found, (list, tuple)) else [found]


def _attached_obs_packages(package) -> list:
    """The OBS packages attached to one FloPy package (FloPy's child-package list)."""

    container = getattr(package, "obs", None)
    if container is None:
        return []
    attached = []
    while True:
        try:
            attached.append(container[len(attached)])
        except (ValueError, IndexError):
            # FloPy's ``MFChildPackages`` raises ValueError past its last package
            # (a list, IndexError); that is how its length is read.
            return attached


def _declared_features(package, ident, *, cellid: bool) -> tuple:
    """What a declared id names: feature numbers from 0, or ``(layer, cell)`` pairs."""

    if cellid:
        if isinstance(ident, (tuple, list, np.ndarray)):
            try:
                return (tuple(int(value) for value in ident),)
            except (TypeError, ValueError):
                return ()
        return ()
    text = str(ident).strip()
    try:
        return (int(text) - 1,)  # MF6 numbers lakes and reaches from one
    except ValueError:
        pass
    packagedata = getattr(package, "packagedata", None)
    data = packagedata.get_data() if packagedata is not None else None
    if data is None or "boundname" not in (data.dtype.names or ()):
        return ()
    key = text.lower()
    return tuple(
        int(row[0])
        for row in data
        if row["boundname"] is not None and str(row["boundname"]).strip().lower() == key
    )


def _declared_observations(packages, *, cellid: bool) -> list[_DeclaredObservation]:
    """Every continuous record the OBS packages of ``packages`` declare, resolved."""

    declared: list[_DeclaredObservation] = []
    for package in packages:
        name = str(getattr(package, "package_name", "") or "")
        for obs in _attached_obs_packages(package):
            continuous = obs.continuous.get_data() or {}
            for csv, records in continuous.items():
                rows = list(records)
                for column, row in enumerate(rows):
                    declared.append(
                        _DeclaredObservation(
                            package=name,
                            csv=str(csv),
                            column=column,
                            width=len(rows),
                            obsname=str(row[0]).strip(),
                            obstype=str(row[1]).strip().lower(),
                            ident=row[2],
                            features=_declared_features(package, row[2], cellid=cellid),
                        )
                    )
    return declared


def _listing(items, *, limit: int = 5) -> str:
    """Up to ``limit`` items joined by ``"; "``, and a count of the rest."""

    items = list(items)
    text = "; ".join(str(item) for item in items[:limit])
    if len(items) > limit:
        text += f"; and {len(items) - limit} more"
    return text


def _named_package(model, quantity: _ObservedQuantity, package: str):
    """The FloPy package of the family's type NAMED ``package`` (any case), or a refusal naming the rest."""

    everyone = _packages_of_type(model, quantity.package_type)
    for candidate in everyone:
        if str(getattr(candidate, "package_name", "")).strip().lower() == str(package).strip().lower():
            return candidate
    held = ", ".join(repr(str(getattr(p, "package_name", ""))) for p in everyone) or "none"
    raise ValueError(
        f"{quantity.owner}: package={package!r} names no {quantity.package_type.upper()} "
        f"package of this model; its {quantity.package_type.upper()} packages: {held}."
    )


def _several_packages(quantity: _ObservedQuantity, packages, reason: str, fix: str) -> ValueError:
    """The refusal for a set that cannot tell which of several packages of its type it reads."""

    names = ", ".join(repr(str(getattr(package, "package_name", ""))) for package in packages)
    return ValueError(
        f"{quantity.owner}: this model carries {len(packages)} {quantity.package_type.upper()} "
        f"packages ({names}); {reason}. {fix}"
    )


def _declared_route(
    model,
    quantity: _ObservedQuantity,
    wanted: dict[str, tuple[str, str, Any]],
    *,
    package: str | None = None,
) -> dict[str, _DeclaredObservation] | None:
    """Which declared observation each target reads, or ``None`` for the file route (ledger 344).

    The observation route used to be chosen, and read, by the SPELLING of a
    target's name: the first CSV anywhere in the workspace, in sorted order,
    whose header held every label. So a LAK target read a GWF HEAD column, a LAK
    VOLUME, an SFR ``DOWNSTREAM-FLOW`` or another lake's stage, and a stale CSV
    from a previous run or a backup copy froze PEST's simulated value in every
    realization. It reads the package's DECLARED observations now:

    - a target takes the observation route only when its label is declared, in
      the package's own OBS, as the family's quantity AND for the target's own
      lake, reach or cell; declared as another quantity, for another feature, or
      more than once as the quantity, it RAISES naming the target, its id and
      the records -- a record of another quantity under the same label is
      another observation, and is left alone when the quantity's own is there;
    - the choice is for the whole SET: every target declared, or none -- a set
      of some of each RAISES naming both groups, because the two routes have not
      always read the same numbers -- MF6 time against period until ledger 346,
      the LAK stage file's record clock until ledger 348, an SFR flow's sign
      until ledger 350, a dry lake's stage until ledger 351, and still a dry
      reach's (-1e30 in the SFR stage file) -- and a silent flip moved every
      target of the set;
    - a set with no target declared takes the family's file route, as before.

    Parameters
    ----------
    model : SimulationBase
        The run model: ``model.package`` finds the packages (a model without it,
        or without a package of the type, has nothing declared).
    quantity : _ObservedQuantity
        The family's package type, observation type and messages.
    wanted : dict of str to tuple of (str, str, object)
        Key to ``(target name, label, identity)``: the label MF6 would write (as
        ``_observation_label`` makes it) and the target's own zero-based lake or
        reach number, or ``(layer, cell)``. A DRN zone has one key per cell.
    package : str or None, optional
        The owning package's NAME, any case: only its observations count. Default
        ``None``: the model's one package of the type.

    Returns
    -------
    dict of str to _DeclaredObservation or None
        The declared observation per key, or ``None`` when no target is declared.

    Raises
    ------
    ValueError
        As listed above; when ``package`` names no package of the type in the
        model; when the model carries several packages of the type, ``package``
        is not given and a label is declared in them; or when the file route
        would have to read one of several packages (LAK, SFR) or one DRN package
        apart from the others.
    """

    everyone = _packages_of_type(model, quantity.package_type)
    packages = everyone if package is None else [_named_package(model, quantity, package)]
    declared = _declared_observations(packages, cellid=quantity.cellid)
    by_label: dict[str, list[_DeclaredObservation]] = {}
    for observation in declared:
        by_label.setdefault(observation.obsname.upper(), []).append(observation)
    hits = {key: by_label.get(str(label).upper(), []) for key, (_name, label, _identity) in wanted.items()}
    if len(packages) > 1 and any(hits.values()):
        raise _several_packages(
            quantity,
            packages,
            "the targets' labels are declared in their observations, and "
            + (
                f"a {quantity.package_type.upper()} observation measures ONE package's drain in "
                "its cell"
                if quantity.cellid
                else f"a {quantity.feature} number counts within ONE package"
            ),
            f"Name the package the targets belong to: {quantity.owner}(..., package=<its name>).",
        )

    chosen: dict[str, _DeclaredObservation] = {}
    missing: list[str] = []
    wrong: list[str] = []
    for key, (name, label, identity) in wanted.items():
        found = hits[key]
        own = [observation for observation in found if observation.obstype == quantity.obstype]
        where = f"target {name!r} ({quantity.feature} {_shown(identity)}, label {str(label).upper()})"
        if not found:
            missing.append(key)
        elif not own:
            wrong.append(
                f"{where} is declared as another quantity, not {quantity.obstype.upper()}: "
                + ", ".join(observation.record() for observation in found)
            )
        elif len(own) > 1:
            wrong.append(
                f"{where} is declared {len(own)} times as {quantity.obstype.upper()}: "
                + ", ".join(observation.record() for observation in own)
            )
        elif own[0].features != (identity,):
            named = ", ".join(_shown(feature) for feature in own[0].features) or "nothing the package has"
            wrong.append(
                f"{where} is declared for {quantity.feature} {named}: {own[0].record()}"
            )
        else:
            chosen[key] = own[0]
    if wrong:
        raise ValueError(
            f"{quantity.owner}: a target's observation label names another observation in "
            f"the {quantity.package_type.upper()} package's own OBS, so the observation route "
            f"would read something else in its place -- {_listing(wrong)}. Correct the "
            "target's id, or rename it (or the observation). "
            + (
                "A (layer, cell) is zero-based, as FloPy holds a cellid."
                if quantity.cellid
                else f"A target's {quantity.feature} number is zero-based; MF6's OBS numbers "
                "them from one."
            )
        )
    if chosen and missing:
        targets = list(dict.fromkeys(name for name, _label, _identity in wanted.values()))
        short = {wanted[key][0] for key in missing}
        declared_names = [name for name in targets if name not in short]
        raise ValueError(
            f"{quantity.owner}: the observation route reads a whole target set or none of it, "
            f"and {len(declared_names)} of {len(targets)} targets are declared in the "
            f"{quantity.package_type.upper()} package's observations "
            f"({_listing(repr(name) for name in declared_names)}) while "
            f"{_listing(repr(name) for name in targets if name in short)} "
            f"{'is' if len(short) == 1 else 'are'} not (labels "
            f"{_listing(str(wanted[key][1]).upper() for key in missing)}). The two routes do "
            "not always read the same numbers (a dry SFR reach is still -1e30 in its stage "
            f"file), so the set would be read from {quantity.file_route} in "
            "silence. Declare the missing observations (to_flopy_obs / attach_flopy_obs), or "
            "give each route its own target set."
        )
    if chosen:
        return chosen
    if len(everyone) > 1 and package is None and not quantity.file_route_reads_every_package:
        raise _several_packages(
            quantity,
            everyone,
            f"none of the targets is declared in their observations, and {quantity.file_route}, "
            f"which they would be read from instead, belongs to ONE {quantity.package_type.upper()} "
            "package",
            "Declare the targets' observations on the package they belong to, and name it: "
            f"{quantity.owner}(..., package=<its name>).",
        )
    if len(everyone) > 1 and package is not None:
        raise _several_packages(
            quantity,
            everyone,
            f"package={package!r} declares none of the targets' labels, and "
            + (
                f"{quantity.file_route} sums every {quantity.package_type.upper()} package's "
                "drains in a zone -- it cannot read one package's alone"
                if quantity.file_route_reads_every_package
                else f"{quantity.file_route} cannot be chosen by package"
            ),
            f"Declare the targets' observations on {package!r}"
            + (", or drop package= to read every package's drains."
               if quantity.file_route_reads_every_package else "."),
        )
    return None


def _observation_output_table(path: Path, *, owner: str) -> tuple[list[str], pd.DataFrame]:
    """One MF6 continuous-observation output file: its header and its rows, by position.

    MF6 writes a continuous block as TEXT (``time,NAME1,NAME2`` then one line per
    time step) or, when the block is declared ``FILEOUT <file> BINARY``, as a
    stream: a 100-byte header opening ``cont double`` (or ``cont single``) with
    ``LENOBSNAME`` in bytes 12-15, the int32 record count, each name blank-padded
    to ``LENOBSNAME``, then one ``time`` plus one value per record per time step in
    that precision (``Obs.f90``'s ``build_headers`` and ``ObsUtility.f90``). The two
    are told apart by the file's first bytes -- ``cont`` against ``time`` -- because
    that is what MF6 wrote, whatever the declaration now says. A binary file was
    read as text, and died on ``UnicodeDecodeError`` (ledger 345).

    Parameters
    ----------
    path : pathlib.Path
        The output file, as the OBS declares it.
    owner : str
        The target class, for messages.

    Returns
    -------
    tuple of (list of str, pandas.DataFrame)
        The header -- ``"time"`` and each record's name as MF6 wrote it -- and the
        rows, columns numbered from 0 (``time``) in the file's order.

    Raises
    ------
    ValueError
        When a binary file's header cannot be read, or its values do not fill a
        whole number of rows (MF6 stopped mid-write, or the file is not MF6's).
    """

    with path.open("rb") as handle:
        opening = handle.read(4)
    if opening != b"cont":
        with path.open() as handle:
            header = [text.strip() for text in handle.readline().rstrip("\r\n").split(",")]
        try:
            # ``round_trip``: MF6 writes each value with every digit it holds (``G0``),
            # and pandas' default fast parser read 4 of the 27 DRN values of ledger
            # 346's toy one unit in the last place off -- the observation route and the
            # budget file then disagreed about a number MF6 wrote once.
            frame = pd.read_csv(path, header=None, skiprows=1, float_precision="round_trip")
        except pd.errors.EmptyDataError:
            frame = pd.DataFrame(columns=range(len(header)))
        return header, frame

    raw = path.read_bytes()
    try:
        control = raw[:100].decode("ascii")
        precision = control[5:11].strip().lower()
        length = int(control[11:15])
        count = int(np.frombuffer(raw, dtype="<i4", count=1, offset=100)[0])
    except (UnicodeDecodeError, ValueError, IndexError) as error:
        raise ValueError(
            f"{owner}: {path} opens like MF6's binary observation output ('cont') but its "
            f"header cannot be read ({error})."
        ) from error
    if precision not in {"single", "double"} or length <= 0 or count < 0:
        raise ValueError(
            f"{owner}: {path} opens like MF6's binary observation output but its header "
            f"reads precision {precision!r}, name length {length}, {count} records."
        )
    start = 104 + count * length
    names = [
        raw[104 + index * length:104 + (index + 1) * length].decode("ascii", errors="replace").strip()
        for index in range(count)
    ]
    itemsize = 8 if precision == "double" else 4
    body = raw[start:]
    if len(body) % ((count + 1) * itemsize):
        raise ValueError(
            f"{owner}: {path} holds {len(body)} bytes of {precision}-precision values after "
            f"its header, not a whole number of rows of {count + 1} (time and {count} "
            "records): MF6 stopped writing it part way, or it is not MF6's file."
        )
    values = np.frombuffer(body, dtype="<f8" if itemsize == 8 else "<f4")
    frame = pd.DataFrame(values.astype(float).reshape(-1, count + 1))
    return ["time", *names], frame


def _observation_workspace(model) -> Path:
    """The directory a model's declared observation files are written in."""

    for attribute in ("workspace", "model_output_folder_path"):
        value = getattr(model, attribute, None)
        if value is not None:
            return Path(value)
    raise ValueError("The model has no workspace to read its declared observation files from.")


def _read_declared_observations(
    model,
    chosen: dict[str, _DeclaredObservation],
    *,
    owner: str,
    no_value: float = 0.0,
) -> pd.DataFrame:
    """Each key's declared observation at each stress period's END: ``time`` (the period) plus one column per key.

    Only the files the package declares are opened, never a workspace scan: a
    CSV left by a previous run in a reused directory (``prepare_run(overwrite=True)``
    keeps MF6's old outputs), a backup copy or another package's file is never
    read (ledger 344). A record is found by its POSITION -- MF6 writes a file's
    records as its columns, in declaration order -- and its header is checked
    against the name MF6 keeps (the first 40 characters, upper-cased), so two
    names MF6 cuts to one header are still two columns. A file MF6 wrote as
    BINARY (``FILEOUT <file> BINARY``) is read as such, by the same rules
    (:func:`_observation_output_table`, ledger 345).

    **One row per stress period, the row that ENDS it, keyed by the zero-based
    period** (ledger 346, the user's rule): MF6 writes a continuous observation at
    every time step, labelled by simulated time, and the route handed every step
    on -- so a target keyed by period matched nothing, and under ATS the number of
    rows, and so what PEST's line-positional instruction file read, moved with the
    parameters (12, 13 or 14 rows at K = 10, 1, 2; a K = 2 realization read t = 23
    as the end of the last period). Which row ends each period is
    :func:`~myflopy.modflow.mf6.period_end.period_end_rows`, checked against the
    model's TDIS by the rule every route reads, so a file short of a period's end
    -- a run that stopped, a truncated copy -- is refused naming the period.

    Parameters
    ----------
    model : SimulationBase
        The run model; its workspace holds the files, and its TDIS
        (``model.sim.tdis``) says how long each stress period is.
    chosen : dict of str to _DeclaredObservation
        What :func:`_declared_route` chose, in the order the columns come back.
    owner : str
        The target class, for messages.
    no_value : float, optional
        What MF6's no-data value (3e30) becomes (see
        :func:`_clean_observation_output_values`). Default ``0.0``, for a flow.

    Returns
    -------
    pandas.DataFrame
        ``time`` -- the zero-based stress period, ``0 .. nper - 1``, every period
        once -- and one float column per key of ``chosen``, each the value MF6
        observed at the end of that period.

    Raises
    ------
    FileNotFoundError
        When a declared file is not in the workspace.
    ValueError
        When a file's columns do not line up with its declared records (a count
        or a header that differs), two files hold different times -- compared
        within the 13 significant digits MF6 writes a text file's time at
        (:data:`_OBSERVATION_TIME_RTOL`), so a text file and a binary one of the
        same run agree -- or a binary
        file cannot be read whole (:func:`_observation_output_table`); when the
        model has no TDIS to tell the periods by; or when the rows do not end
        every period the TDIS declares
        (:func:`~myflopy.modflow.mf6.period_end.period_end_rows`).
    """

    workspace = _observation_workspace(model)
    files: dict[str, list[tuple[str, _DeclaredObservation]]] = {}
    for key, observation in chosen.items():
        files.setdefault(observation.csv, []).append((key, observation))
    columns: dict[str, pd.Series] = {}
    time: pd.Series | None = None
    time_file = None
    for csv, entries in files.items():
        path = workspace / csv
        if not path.is_file():
            raise FileNotFoundError(
                f"{owner}: {entries[0][1].record()} declares its output in {csv}, which is "
                f"not in {workspace}. MF6 writes it when the model runs; has this workspace "
                "run since the observation was declared?"
            )
        header, frame = _observation_output_table(path, owner=owner)
        width = entries[0][1].width
        if len(header) != width + 1 or frame.shape[1] != width + 1:
            raise ValueError(
                f"{owner}: {path} holds {len(header) - 1} observation column(s) where its OBS "
                f"declares {width} records ({entries[0][1].package!r}): the file and the "
                "declaration do not line up."
            )
        for key, observation in entries:
            written = header[observation.column + 1].upper()
            if written != observation.obsname.upper()[:_MF6_OBSNAME_LENGTH]:
                raise ValueError(
                    f"{owner}: column {observation.column + 1} of {path} is headed {written!r} "
                    f"where {observation.record()} is declared there: the file and the "
                    "declaration do not line up."
                )
            columns[key] = _clean_observation_output_values(
                frame[observation.column + 1], no_value=no_value
            )
        this_time = pd.to_numeric(frame[0], errors="coerce").astype(float).reset_index(drop=True)
        if time is None:
            time, time_file = this_time, csv
        elif len(time) != len(this_time) or not np.allclose(
            time.to_numpy(), this_time.to_numpy(), rtol=_OBSERVATION_TIME_RTOL, atol=0.0
        ):
            # Within the 13 digits MF6 writes a text file's time at, not bit for bit:
            # a binary file beside it holds the same step's time in full (ledger 347).
            raise ValueError(
                f"{owner}: the declared observation files {time_file} and {csv} hold different "
                "times, so their rows cannot be read side by side."
            )
    result = pd.DataFrame({"time": time if time is not None else pd.Series(dtype=float)})
    for key in chosen:
        result[key] = columns[key].reset_index(drop=True).to_numpy()

    perlen = period_lengths(model)
    packages = sorted({observation.package for observation in chosen.values()})
    described = (
        f"the observation file{'s' if len(files) > 1 else ''} {', '.join(files)} "
        f"(package{'s' if len(packages) > 1 else ''} {', '.join(repr(p) for p in packages)})"
    )
    if perlen is None:
        raise ValueError(
            f"{owner}: cannot tell which row of {described} ends each stress period -- "
            "MF6 labels its rows by simulated time alone, and the model carries no TDIS "
            "saying how long each period is. The observation route reads one row per "
            "stress period, at its end (ledger 346)."
        )
    rows = period_end_rows(result["time"], perlen, file=described)
    periods = list(rows)
    reduced = result.iloc[[rows[per] for per in periods]].reset_index(drop=True)
    reduced["time"] = periods
    return reduced


def _period_key_hints(keys, perlen) -> list[str]:
    """For each key that is a stress period's END in simulated time, the period key to use.

    The observation route was keyed by MF6's simulated time until ledger 346, so a
    target set written for it -- ``10, 20, 30`` for three 10-day periods -- names
    period ends, not periods. This says which period each such key ends; a key
    that is no period's end gets no hint.

    Parameters
    ----------
    keys : iterable of float
        The keys that name no stress period of the model.
    perlen : sequence of float or None
        Each period's length (:func:`~myflopy.modflow.mf6.period_end.period_lengths`);
        ``None`` or empty gives no hints.

    Returns
    -------
    list of str
        One ``"10 is the end of period 0 in simulated time -- key it 0"`` per key
        that is a period's end, in the order given.
    """

    if not perlen:
        return []
    ends = np.cumsum(np.asarray(list(perlen), dtype=float))
    hints = []
    for key in keys:
        number = float(key)
        # A hint, not a match: a key typed to six or seven figures still names its end.
        near = np.flatnonzero(np.isclose(ends, number, rtol=1e-6, atol=1e-9))
        if near.size:
            period = int(near[0])
            hints.append(
                f"{_shown(key)} is the end of period {period} in simulated time -- key it {period}"
            )
    return hints


def _refuse_keys_off_the_periods(
    targets: pd.DataFrame,
    periods,
    *,
    model,
    caller: str,
) -> None:
    """Refuse a target keyed by a period the observation route does not have (ledger 346).

    The observation route is one row per stress period now, keyed ``0 .. nper - 1``,
    where it was one row per time step keyed by MF6's simulated time -- the key its
    docs taught (``10, 20, 30``). Such a target used to match; now it would match
    nothing, and ``compare()`` would answer ``NaN`` for it while ``stats()``
    counted one row fewer. So it is refused, naming the period each simulated-time
    key ends. A time key that is ALSO a period -- any key below ``nper``: ``1`` on
    daily periods, and ``2, 4, 6`` on twelve 2-day periods, the ENDS of periods 0,
    1, 2 -- cannot be told from one and reads that period in silence (ledger 347);
    a set is refused only when one of its keys names no period.

    Parameters
    ----------
    targets : pandas.DataFrame
        The target rows, with ``name``, ``time`` (the key as given) and ``per``
        (read by :func:`_numeric_periods`).
    periods : iterable of int
        The periods the route's table holds.
    model : SimulationBase
        The run model, whose TDIS gives the hints.
    caller : str
        What is reading, named in the error (``"DrnFlowTargets.compare(...)"``).

    Raises
    ------
    ValueError
        When any target's ``per`` is not one of ``periods``, naming up to five by
        target and key, and the period each simulated-time key ends.
    """

    held = sorted({int(per) for per in periods})
    off = targets.loc[~targets["per"].isin(held)]
    if off.empty:
        return
    listing = [
        f"{_shown(row.name)} at time {_shown(row.time)}"
        for row in off.loc[:, ["name", "time"]].drop_duplicates().itertuples(index=False)
    ]
    hints = _period_key_hints(pd.unique(off["time"]), period_lengths(model))
    span = f"0 .. {held[-1]}" if held else "none"
    raise ValueError(
        f"{caller}: the observation route is keyed by zero-based stress period, one row per "
        f"period at its end (ledger 346), and this model's periods are {span}; "
        f"{len(listing)} target row(s) name none of them: {_listing(listing)}. "
        + (f"{_listing(hints)}. " if hints else "")
        + "The route was keyed by MF6's simulated time until ledger 346; key each target by "
        "its stress period, 0, 1, 2, ..."
    )


def _refuse_observation_kind(kind, quantity: _ObservedQuantity) -> str:
    """``kind`` upper-cased, when it is the family's quantity; anything else raises (ledger 344)."""

    written = str(kind).strip().upper()
    if written != quantity.obstype.upper():
        raise ValueError(
            f"{quantity.owner} writes {quantity.obstype.upper()} observations, the quantity it "
            f"reads back; kind={kind!r} would be declared under the targets' names and refused "
            "when read (a label declared as another quantity is not the target's). Declare "
            "another observation type with FloPy directly, under a name of its own."
        )
    return written


def _refuse_shared_labels(names, quantity: _ObservedQuantity) -> None:
    """Refuse two target names MF6 would read as one observation label (ledger 344).

    MF6 upper-cases an observation's name, so ``"Pond A"`` and ``"pond-a"`` --
    written ``Pond_A`` and ``pond_a`` -- are one label: both were accepted, MF6
    wrote the header ``time,POND_A,POND_A``, and both targets read the first lake.
    """

    owners: dict[str, list[str]] = {}
    for name in names:
        owners.setdefault(_observation_label(name).upper(), []).append(str(name))
    clashes = {label: named for label, named in owners.items() if len(named) > 1}
    if clashes:
        raise ValueError(
            f"{quantity.owner}: two target names would be written as one MF6 observation label "
            "(letters and digits only, any case): "
            + _listing(f"{label} from {', '.join(repr(name) for name in named)}"
                       for label, named in clashes.items())
            + ". Rename one of them."
        )


def _ensure_unique_observation_names(
    frame: pd.DataFrame,
    *,
    base_column: str = "name",
    layer_column: str = "layer",
    cell_column: str = "cell",
) -> pd.Series:
    """Build deterministic unique observation names for FloPy OBS records."""

    counts = frame[base_column].value_counts(dropna=False)
    names: list[str] = []
    used: dict[str, int] = {}
    for row in frame.itertuples(index=False):
        base_name = _observation_label(getattr(row, base_column))
        if counts.get(getattr(row, base_column), 0) > 1:
            base_name = f"{base_name}_lay{int(getattr(row, layer_column))}"
        candidate = base_name
        if candidate in used:
            candidate = f"{candidate}_c{int(getattr(row, cell_column))}"
        if candidate in used:
            used[candidate] += 1
            candidate = f"{candidate}_{used[candidate]}"
        else:
            used[candidate] = 1
        names.append(candidate)
    return pd.Series(names, index=frame.index, name="obsname")


def _load_locations(
    locations: str | Path | gpd.GeoDataFrame | pd.DataFrame | pd.Series | dict | list | tuple,
    *,
    name_column: str,
    layer_column: str | None,
    group_column: str | None,
    weight_column: str | None,
) -> pd.DataFrame:
    """Load and validate a head-target locations layer."""

    if isinstance(locations, gpd.GeoDataFrame):
        frame = locations.copy()
    elif isinstance(locations, pd.Series):
        index_name = locations.index.name or name_column
        frame = locations.rename("cell").rename_axis(index_name).reset_index()
        if index_name != name_column:
            frame = frame.rename(columns={index_name: name_column})
    elif isinstance(locations, dict):
        lower_keys = {str(key).strip().lower() for key in locations.keys()}
        if {"name", "cell"} & lower_keys or {"name", "node"} & lower_keys:
            frame = pd.DataFrame(locations)
        else:
            frame = pd.DataFrame({name_column: list(locations.keys()), "cell": list(locations.values())})
    elif isinstance(locations, (list, tuple)):
        frame = pd.DataFrame(locations)
    elif isinstance(locations, pd.DataFrame):
        frame = locations.copy()
    else:
        frame = gpd.read_file(locations)

    required = [name_column]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(
            "Head target locations are missing required columns: "
            + ", ".join(missing)
        )

    keep_columns = [name_column]
    cell_column = None
    for candidate in ("cell", "node"):
        if candidate in frame.columns:
            cell_column = candidate
            keep_columns.append(candidate)
            break
    for column in (layer_column, group_column, weight_column):
        if column is not None and column in frame.columns:
            keep_columns.append(column)
    if isinstance(frame, gpd.GeoDataFrame):
        if frame.geometry.is_empty.any():
            raise ValueError("Head target locations contain empty geometries.")
        keep_columns.append(frame.geometry.name)
    elif cell_column is None:
        raise ValueError(
            "Non-geospatial head target locations must include a 'cell' or 'node' column."
        )
    frame = frame.loc[:, list(dict.fromkeys(keep_columns))].copy()

    frame = frame.rename(columns={name_column: "name"})
    if cell_column is not None:
        frame = frame.rename(columns={cell_column: "cell"})
    if layer_column and layer_column in frame.columns:
        frame = frame.rename(columns={layer_column: "layer"})
    if group_column and group_column in frame.columns:
        frame = frame.rename(columns={group_column: "group"})
    if weight_column and weight_column in frame.columns:
        frame = frame.rename(columns={weight_column: "weight"})

    if "layer" not in frame.columns:
        frame["layer"] = 0
    else:
        frame["layer"] = pd.to_numeric(frame["layer"], errors="coerce").fillna(0).astype(int)

    if "group" not in frame.columns:
        frame["group"] = pd.NA
    if "weight" not in frame.columns:
        frame["weight"] = np.nan

    frame["name"] = frame["name"].astype(str)
    frame["name_key"] = _normalize_identifier_series(frame["name"])
    if "cell" in frame.columns:
        frame["cell"] = pd.to_numeric(frame["cell"], errors="coerce")
        if frame["cell"].isna().any():
            raise ValueError("Head target cell/node values must be numeric.")
        frame["cell"] = frame["cell"].astype(int)
    if isinstance(frame, gpd.GeoDataFrame):
        frame["x"] = frame.geometry.centroid.x
        frame["y"] = frame.geometry.centroid.y
    else:
        frame["x"] = pd.to_numeric(frame.get("x", np.nan), errors="coerce")
        frame["y"] = pd.to_numeric(frame.get("y", np.nan), errors="coerce")
    return frame


def _normalize_values_frame(
    values: str | Path | pd.DataFrame | pd.Series | dict | list | tuple,
    *,
    time_column: str,
    value_column: str,
    times=None,
    owner: str = "HeadTargets",
) -> pd.DataFrame:
    """Normalize long or wide target tables to a consistent long format.

    A blank value is ``NaN`` (unmeasured) and any other value must be a number
    (:func:`_coerce_target_values`, which names ``owner`` when it refuses one).
    """

    if isinstance(values, pd.Series):
        index_name = values.index.name or time_column
        frame = values.rename(value_column).rename_axis(index_name).reset_index()
        if index_name != time_column:
            frame = frame.rename(columns={index_name: time_column})
    elif isinstance(values, dict):
        frame = pd.DataFrame(values)
        if time_column not in frame.columns and times is not None:
            frame.insert(0, time_column, list(times))
    elif isinstance(values, (list, tuple)):
        frame = pd.DataFrame(values)
        if time_column not in frame.columns and times is not None:
            if frame.shape[1] == 1:
                frame.columns = [value_column]
            if frame.shape[0] != len(list(times)):
                raise ValueError(
                    f"Provided times length {len(list(times))} does not match values length {frame.shape[0]}."
                )
            frame.insert(0, time_column, list(times))
    elif isinstance(values, pd.DataFrame):
        frame = values.copy()
    else:
        frame = pd.read_csv(values)

    if time_column not in frame.columns:
        raise ValueError(
            f"Head target values must include the time/per column {time_column!r}."
        )

    lower_columns = {str(column).strip().lower(): column for column in frame.columns}
    name_column = lower_columns.get("name")
    explicit_value_column = lower_columns.get(value_column.strip().lower())

    if name_column is not None and explicit_value_column is not None:
        long = frame.rename(
            columns={
                time_column: "time",
                name_column: "name",
                explicit_value_column: "head_target",
            }
        ).loc[:, ["time", "name", "head_target"]]
        source = explicit_value_column
    elif name_column is None and explicit_value_column is not None and len(frame.columns) == 2:
        long = frame.rename(columns={time_column: "time", explicit_value_column: "head_target"}).loc[
            :, ["time", "head_target"]
        ]
        source = explicit_value_column
    else:
        value_columns = [column for column in frame.columns if column != time_column]
        long = frame.melt(
            id_vars=[time_column],
            value_vars=value_columns,
            var_name="name",
            value_name="head_target",
        ).rename(columns={time_column: "time"})
        source = long["name"]

    if "name" in long.columns:
        long["name"] = long["name"].astype(str)
        long["name_key"] = _normalize_identifier_series(long["name"])
    long["head_target"] = _coerce_target_values(long, "head_target", owner=owner, source=source)
    return _copy_frame(long)


#: The reader behind each dependent-variable table, what to call its file, and the
#: output-control word that schedules it -- for the period-end rule's refusal.
_FIELD_FILES = {
    "all_heads": ("hds", "the head file", "HEAD"),
    "all_conc": ("conc", "the concentration file", "CONCENTRATION"),
}


def _simulated_field_by_period(
    model: Any,
    *,
    table_attribute: str = "all_heads",
    store_column: str = "elev",
) -> pd.DataFrame:
    """One simulated value per period/layer/cell: each period's END (ledger 336).

    The record read for a period is the one that ENDS it, by the one rule every
    target uses (:func:`~myflopy.modflow.mf6.period_end.period_ends`), on the
    dependent-variable file's own clock. A period whose end the file does not
    hold -- output control saved ``FIRST`` or ``STEPS`` short of it -- is refused
    rather than read at an earlier step, which it used to be: the last saved step
    was taken, wherever it fell.

    ``table_attribute``/``store_column`` select the field: heads
    (``all_heads``/``elev``) or concentration (``all_conc``/``conc``). The output
    column is ``sim_head`` for BOTH -- it is a frame-internal name meaning "the
    simulated value", and keeping it shared is what lets the calibration-plot and
    PEST machinery in ``calcs/calibration.py`` serve concentration unchanged
    (ledger 102).

    Parameters
    ----------
    model : SimulationBase
        A run model holding the dependent-variable file.
    table_attribute : {"all_heads", "all_conc"}, default "all_heads"
        The model's ``kstpkper``/``layer``/``cell`` table to read.
    store_column : str, default "elev"
        That table's value column (``"elev"`` for heads, ``"conc"``).

    Returns
    -------
    pandas.DataFrame
        ``per``, ``layer``, ``cell``, ``sim_head``: one row per cell per period
        the file saved, the value at that period's end.

    Raises
    ------
    ValueError
        When a period the file saved lacks its final step, naming the period,
        the steps saved and the ``SAVE HEAD LAST`` setting that fixes it.
    """

    frame = getattr(model, table_attribute).reset_index().copy()
    reader_attribute, file, saves = _FIELD_FILES[table_attribute]
    ends = period_ends(
        model,
        reader=lambda: getattr(model, reader_attribute),
        file=file,
        saves=saves,
        fallback=frame["kstpkper"],
    )
    frame = at_period_end(frame, ends).copy()
    frame["per"] = frame["kstpkper"].map(lambda item: int(item[1]))
    frame = frame.sort_values("per", kind="stable").rename(columns={store_column: "sim_head"})
    return frame.loc[:, ["per", "layer", "cell", "sim_head"]].reset_index(drop=True)


def _simulated_heads_by_period(model: Any) -> pd.DataFrame:
    """Simulated heads per period/layer/cell (the heads instance of the above)."""

    return _simulated_field_by_period(model)


def _default_obs_csv_name(filename: str | Path, default_name: str) -> str:
    """Derive a default MF6 observation CSV name from an OBS control filename."""

    if filename is None:
        return default_name
    stem = Path(filename).stem
    return f"{stem}.csv"


def _coerce_label_list(values, *, count: int, default_prefix: str) -> list[str]:
    """Normalize optional names/groups to a length-matched list."""

    if values is None:
        return [f"{default_prefix}_{idx + 1:03d}" for idx in range(count)]
    if isinstance(values, str):
        if count != 1:
            raise ValueError(f"{default_prefix} label scalar only valid for one target.")
        return [values]
    items = list(values)
    if len(items) != count:
        raise ValueError(f"Expected {count} {default_prefix} labels, received {len(items)}.")
    return [str(item) for item in items]


def _coerce_scalar_or_list(values, *, count: int, default=None):
    """Broadcast scalars or validate sequence lengths."""

    if values is None:
        return [default] * count
    if isinstance(values, str) or not isinstance(values, (list, tuple, np.ndarray, pd.Series)):
        return [values] * count
    items = list(values)
    if len(items) != count:
        raise ValueError(f"Expected {count} values, received {len(items)}.")
    return items


def _build_values_from_cells(
    values,
    *,
    names: list[str],
    times,
    time_column: str,
    value_column: str,
) -> pd.DataFrame:
    """Normalize cell-based target values into the standard long-format input."""

    if isinstance(values, pd.Series):
        if len(names) != 1:
            raise ValueError("Series input for from_cells(...) only supports one target name.")
        frame = values.rename(names[0]).reset_index()
        frame.columns = [time_column, names[0]]
        return frame
    if isinstance(values, dict):
        frame = pd.DataFrame(values)
        if time_column not in frame.columns:
            if times is None:
                raise ValueError("times is required when values dict does not include the time column.")
            frame.insert(0, time_column, list(times))
        return frame
    if isinstance(values, pd.DataFrame):
        frame = values.copy()
        if time_column not in frame.columns:
            if times is None:
                raise ValueError("times is required when values DataFrame does not include the time column.")
            frame.insert(0, time_column, list(times))
        unnamed = [column for column in frame.columns if column != time_column]
        if unnamed == list(range(len(unnamed))):
            rename = {column: names[idx] for idx, column in enumerate(unnamed)}
            frame = frame.rename(columns=rename)
        return frame

    array = np.asarray(values, dtype=object)
    if array.ndim == 1:
        if len(names) != 1:
            raise ValueError("One-dimensional values for from_cells(...) require exactly one target.")
        if times is None:
            times = list(range(len(array)))
        return pd.DataFrame({time_column: list(times), names[0]: array.tolist()})
    if array.ndim != 2:
        raise ValueError("from_cells(...) values must be 1D or 2D.")
    if array.shape[1] != len(names):
        raise ValueError(
            f"from_cells(...) expected {len(names)} value columns, received {array.shape[1]}."
        )
    if times is None:
        times = list(range(array.shape[0]))
    return pd.DataFrame(array.tolist(), columns=names).assign(**{time_column: list(times)})[
        [time_column, *names]
    ]


def _normalize_lake_stage_values(
    values: str | Path | pd.DataFrame | pd.Series | dict | list | tuple | None,
    *,
    time_column: str,
    value_column: str,
    times=None,
    owner: str = "LakeStageTargets",
) -> pd.DataFrame:
    """Normalize lake-stage target values to a long format.

    A blank stage is ``NaN`` (unmeasured) and any other entry must be a number
    (:func:`_coerce_target_values`).
    """

    if values is None:
        return pd.DataFrame(columns=["time", "name", "stage"])
    if isinstance(values, pd.Series):
        index_name = values.index.name or time_column
        frame = values.rename(value_column).rename_axis(index_name).reset_index()
        if index_name != time_column:
            frame = frame.rename(columns={index_name: time_column})
    elif isinstance(values, dict):
        frame = pd.DataFrame(values)
        if time_column not in frame.columns and times is not None:
            frame.insert(0, time_column, list(times))
    elif isinstance(values, (list, tuple)):
        frame = pd.DataFrame(values)
        if time_column not in frame.columns and times is not None:
            if frame.shape[1] == 1:
                frame.columns = [value_column]
            if frame.shape[0] != len(list(times)):
                raise ValueError(
                    f"Provided times length {len(list(times))} does not match values length {frame.shape[0]}."
                )
            frame.insert(0, time_column, list(times))
    elif isinstance(values, pd.DataFrame):
        frame = values.copy()
    else:
        frame = pd.read_csv(values)

    if time_column not in frame.columns:
        raise ValueError(f"Lake-stage target values must include {time_column!r}.")
    lower_columns = {str(column).strip().lower(): column for column in frame.columns}
    name_column = lower_columns.get("name")
    explicit_value_column = lower_columns.get(value_column.strip().lower())
    if name_column is not None and explicit_value_column is not None:
        long = frame.rename(
            columns={
                time_column: "time",
                name_column: "name",
                explicit_value_column: "stage",
            }
        ).loc[:, ["time", "name", "stage"]]
        source = explicit_value_column
    elif name_column is None and explicit_value_column is not None and len(frame.columns) == 2:
        long = frame.rename(columns={time_column: "time", explicit_value_column: "stage"}).loc[
            :, ["time", "stage"]
        ]
        source = explicit_value_column
    else:
        value_columns = [column for column in frame.columns if column != time_column]
        long = frame.melt(
            id_vars=[time_column],
            value_vars=value_columns,
            var_name="name",
            value_name="stage",
        ).rename(columns={time_column: "time"})
        source = long["name"]
    if "name" in long.columns:
        long["name"] = long["name"].astype(str)
    long["stage"] = _coerce_target_values(long, "stage", owner=owner, source=source)
    return _copy_frame(long)


def _normalize_named_integer_locations(
    locations,
    *,
    id_column: str,
    target_name: str,
) -> pd.DataFrame:
    """Normalize simple ``name -> integer id`` location definitions."""

    source = {} if locations is None else locations
    if isinstance(source, pd.Series):
        source = source.to_dict()
    elif isinstance(source, pd.DataFrame):
        required = {"name", id_column}
        if required - set(source.columns):
            missing = sorted(required - set(source.columns))
            raise ValueError(f"{target_name} locations DataFrame is missing columns: {', '.join(missing)}")
        source = source.loc[:, ["name", id_column]].drop_duplicates(subset=["name"]).set_index("name")[id_column].to_dict()
    elif isinstance(source, (list, tuple)):
        source_frame = pd.DataFrame(source)
        required = {"name", id_column}
        if required - set(source_frame.columns):
            missing = sorted(required - set(source_frame.columns))
            raise ValueError(f"{target_name} location records are missing columns: {', '.join(missing)}")
        source = (
            source_frame.loc[:, ["name", id_column]]
            .drop_duplicates(subset=["name"])
            .set_index("name")[id_column]
            .to_dict()
        )

    normalized = {
        str(name): int(value)
        for name, value in dict(source).items()
    }
    return pd.DataFrame([{"name": name, id_column: value} for name, value in normalized.items()])


def _normalize_named_series_values(
    values: str | Path | pd.DataFrame | pd.Series | dict | list | tuple | None,
    *,
    time_column: str,
    value_column: str,
    output_column: str,
    times=None,
    owner: str = "target",
) -> pd.DataFrame:
    """Normalize long or wide named-series targets to one long table.

    A blank value is ``NaN`` (unmeasured) and any other entry must be a number
    (:func:`_coerce_target_values`, which names ``owner`` when it refuses one).
    """

    if values is None:
        return pd.DataFrame(columns=["time", "name", output_column])
    if isinstance(values, pd.Series):
        index_name = values.index.name or time_column
        frame = values.rename(value_column).rename_axis(index_name).reset_index()
        if index_name != time_column:
            frame = frame.rename(columns={index_name: time_column})
    elif isinstance(values, dict):
        frame = pd.DataFrame(values)
        if time_column not in frame.columns and times is not None:
            frame.insert(0, time_column, list(times))
    elif isinstance(values, (list, tuple)):
        frame = pd.DataFrame(values)
        if time_column not in frame.columns and times is not None:
            if frame.shape[1] == 1:
                frame.columns = [value_column]
            if frame.shape[0] != len(list(times)):
                raise ValueError(
                    f"Provided times length {len(list(times))} does not match values length {frame.shape[0]}."
                )
            frame.insert(0, time_column, list(times))
    elif isinstance(values, pd.DataFrame):
        frame = values.copy()
    else:
        frame = pd.read_csv(values)

    if time_column not in frame.columns:
        raise ValueError(f"Target values must include {time_column!r}.")

    lower_columns = {str(column).strip().lower(): column for column in frame.columns}
    name_column = lower_columns.get("name")
    explicit_value_column = lower_columns.get(value_column.strip().lower())

    if name_column is not None and explicit_value_column is not None:
        long = frame.rename(
            columns={
                time_column: "time",
                name_column: "name",
                explicit_value_column: output_column,
            }
        ).loc[:, ["time", "name", output_column]]
        source = explicit_value_column
    elif name_column is None and explicit_value_column is not None and len(frame.columns) == 2:
        long = frame.rename(columns={time_column: "time", explicit_value_column: output_column}).loc[
            :, ["time", output_column]
        ]
        source = explicit_value_column
    else:
        value_columns = [column for column in frame.columns if column != time_column]
        long = frame.melt(
            id_vars=[time_column],
            value_vars=value_columns,
            var_name="name",
            value_name=output_column,
        ).rename(columns={time_column: "time"})
        source = long["name"]
    if "name" in long.columns:
        long["name"] = long["name"].astype(str)
    long[output_column] = _coerce_target_values(long, output_column, owner=owner, source=source)
    return _copy_frame(long)


def _numeric_periods(frame: pd.DataFrame, *, caller: str) -> pd.Series:
    """Zero-based stress periods read from a table's ``per`` (else ``time``) column.

    The ONE reading of a target's key as a period, for every family and route and
    for PEST's build (ledger 340). A period is a whole number: ``1.0`` is period 1,
    and a key with a fraction is REFUSED. It used to be truncated, so a target at
    ``10.9`` read period 10's value, and one at ``3.3333`` period 3's, with nothing
    said -- while the PEST build, matching by the key's spelling, dropped both.

    Parameters
    ----------
    frame : pandas.DataFrame
        A target or result table with a ``per`` or ``time`` column (``per``
        wins when it has both).
    caller : str
        What is reading it, named in the error (``"HeadTargets.compare(...)"``).

    Returns
    -------
    pandas.Series
        The periods as ``int``, aligned with ``frame``.

    Raises
    ------
    ValueError
        When the table has neither column; or when a key is not a finite number,
        or has a fraction -- naming up to five such keys.
    """

    if "per" in frame.columns:
        column = "per"
    elif "time" in frame.columns:
        column = "time"
    else:
        raise ValueError(f"{caller} requires a 'per' or 'time' column.")
    raw = frame[column]
    per = pd.to_numeric(raw, errors="coerce").astype(float)
    finite = np.isfinite(per.to_numpy())
    if not finite.all():
        labels = ", ".join(_shown(value) for value in pd.unique(raw[~finite])[:5])
        raise ValueError(
            f"{caller} reads `{column}` as a zero-based stress period, and "
            f"{labels} {'is' if (~finite).sum() == 1 else 'are'} not a number."
        )
    whole = per.to_numpy() == np.floor(per.to_numpy())
    if not whole.all():
        labels = ", ".join(_shown(value) for value in pd.unique(raw[~whole])[:5])
        raise ValueError(
            f"{caller} reads `{column}` as a zero-based stress period, and {labels} "
            f"{'is' if (~whole).sum() == 1 else 'are'} not a whole number: a period "
            "has no fraction, so there is no period to read for it. Key each target "
            "by its stress period, 0, 1, 2, ... (1.0 is period 1)."
        )
    return per.astype(int)


def _measured_periods(
    values: pd.DataFrame,
    *,
    value_column: str,
    caller: str,
    name_column: str = "name",
) -> pd.Series:
    """Each target row's period, refusing two measurements of one target in one period.

    PEST holds ONE observation per target per period, and the build assigned each
    row in turn, so the LAST of two readings in one period silently won there while
    ``compare()`` and ``stats()`` counted both (ledger 340). Which one a period's
    observation should be -- a mean, the last reading -- is a choice about the data,
    so it is the caller's: two rows HOLDING a value for one (target, period) are
    refused, naming them. A blank row is not a measurement and never collides.

    Parameters
    ----------
    values : pandas.DataFrame
        A target set's normalized values table (``time``, ``name_column`` and
        ``value_column``), one row per target per time as given -- NOT a table
        already joined to the locations, which repeats a row for each location
        row sharing the name (two screens of one well).
    value_column : str
        The measured value's column (``"head_target"``, ``"stage"``, ...).
    caller : str
        What is reading it, named in the errors.
    name_column : str, default "name"
        The column identifying the target; compared stripped and lower-cased.

    Returns
    -------
    pandas.Series
        The zero-based periods as ``int``, aligned with ``values``
        (:func:`_numeric_periods`, whose refusals apply).

    Raises
    ------
    ValueError
        As :func:`_numeric_periods`; or when two rows hold a value for one target
        in one period, naming up to five such keys.
    """

    per = _numeric_periods(values, caller=caller)
    if value_column in values.columns:
        held = values[value_column].notna()
    else:
        held = pd.Series(True, index=values.index)
    keys = pd.DataFrame(
        {
            "key": _normalize_identifier_series(values[name_column]),
            "name": values["name"] if "name" in values.columns else values[name_column],
            "per": per,
        }
    )[held.to_numpy()]
    repeated = keys.duplicated(subset=["key", "per"], keep=False)
    if repeated.any():
        groups = keys[repeated].groupby(["key", "per"], sort=True)
        listing = [
            f"{_shown(group['name'].iloc[0])} at period {period}: {len(group)} values"
            for (_key, period), group in list(groups)[:5]
        ]
        more = groups.ngroups - len(listing)
        raise ValueError(
            f"{caller}: two or more measurements of one target in one stress period -- "
            + "; ".join(listing)
            + (f"; and {more} more" if more else "")
            + ". PEST holds ONE observation per target per period, so it would keep "
            "whichever row came last. Reduce them to one value per period first (a "
            "mean, the last reading): which one is a choice about your data."
        )
    return per


def _aggregate_named_series(
    frame: pd.DataFrame,
    *,
    id_column: str,
    value_column: str,
    kind: str,
    ends: dict[int, tuple[int, int]] | None = None,
) -> pd.DataFrame:
    """Reduce a result table to ONE value per ``per`` and integer id, by what the value IS.

    Rows ADD only when the value is extensive -- ledger 297's rule, asked of
    :func:`~myflopy.modflow.mf6.package_registry.agg_for_kind` rather than restated
    here. A flow's rows are the cells or connections of one feature in one time
    step, and they add: that is what a seepage zone or a reach's outflow IS. A
    stage table repeats a feature's ONE stage on every cell the feature connects
    to, so its rows of one id in one record are copies, and adding them read a
    pond at 526.00 ft as 401,339.97 ft -- 763 connected cells (ledger 331). For a
    non-extensive value the copies are reduced to the one number they hold, which
    is what the registry's ``"mean"`` gives for copies that agree; copies that
    DISAGREE are refused, because averaging them would invent a value the table
    never held and a table like that is not the one this assumes.

    **Time steps never add, whatever the kind** (ledgers 331, 336). A table that
    carries ``kstpkper`` is first reduced to the record that ENDS each period --
    the one rule every target reads,
    :func:`~myflopy.modflow.mf6.period_end.period_end_steps` -- so a period's
    value is its end-of-period value. Adding a period's saved steps multiplied an
    SFR stage, an SFR flow and a DRN zone's flow by the number of steps saved
    (3.0x at ``nstp=3``). Every result table carries it -- the LAK stage table since
    ledger 348, whose ``per`` was the index of the saved record until then (ledger
    301's stage clock); a table without it (a test double's) keeps every row.

    A row whose VALUE is missing is kept, and an id with no value in a period comes
    back ``NaN`` rather than absent or 0.0: a lake with no stage (ledger 329) is a
    row in the stage table, and dropping it deleted that period from the wide
    series -- which ``write_named_series_targets`` hands PEST, whose instruction
    file advances one line per row of the BASE run, so a forward run one row short
    fails at end of file (``docs/reproducers/pest_no_value_fields.py``).
    ``sum(min_count=1)`` keeps an all-``NaN`` group ``NaN``; a plain ``sum``
    makes it 0.0.

    Parameters
    ----------
    frame : pandas.DataFrame
        One row per (record, id, cell): a ``per`` (or ``time``) column holding
        zero-based integers, ``id_column``, ``value_column`` and, optionally,
        ``kstpkper`` tuples on the clock of the file the table was read from.
    id_column : str
        The integer feature id each row belongs to (``"lake"``, ``"reach"``).
    value_column : str
        The column holding the value to reduce (``"stage"``, ``"sim_flow"``).
    kind : str
        What the value is, in :class:`~myflopy.modflow.mf6.package_registry.FieldSpec`'s
        words: ``"extensive"`` (a flux -- rows of one step add), or
        ``"elevation"`` / ``"intensive"`` (a level -- rows are copies of one
        number). Required, so a new caller has to say.
    ends : dict of int to tuple of (int, int), optional
        The record that ends each period, ``{kper: (kstp, kper)}``, as
        ``period_end.period_ends`` reads it off the table's own file
        and checks it against the TDIS. ``None`` (the default) takes each
        period's last ``kstpkper`` in the table, unchecked -- for a table with no
        file behind it. Ignored for a table without ``kstpkper``.

    Returns
    -------
    pandas.DataFrame
        ``per``, ``id_column`` and ``value_column``, one row per (``per``, id),
        sorted by both; ``NaN`` where the id holds no value in that ``per``.

    Raises
    ------
    ValueError
        If ``kind`` is not one of the three; if a ``per`` or ``time`` value is not
        numeric; or, for a non-extensive ``kind``, if the rows of one id in one
        ``per`` disagree by more than ``1e-9`` of the value (at least ``1e-9``) or
        some hold a value and some do not -- naming the id, the ``per`` and the
        spread.
    """

    combine = agg_for_kind(kind)
    data = frame.copy()
    data["per"] = _numeric_periods(data, caller="_aggregate_named_series(...)")
    data[id_column] = pd.to_numeric(data[id_column], errors="coerce")
    data[value_column] = pd.to_numeric(data[value_column], errors="coerce")
    data = data.dropna(subset=[id_column]).copy()
    data[id_column] = data[id_column].astype(int)
    if "kstpkper" in data.columns and not data.empty:
        if ends is None:
            ends = period_end_steps(data["kstpkper"])
        data = at_period_end(data, ends)
    grouped = data.groupby(["per", id_column], dropna=False)[value_column]
    if combine == "sum":
        reduced = grouped.sum(min_count=1)
    else:
        reduced = _one_value_per_group(grouped, id_column=id_column, value_column=value_column)
    return reduced.reset_index().sort_values(["per", id_column]).reset_index(drop=True)


def _stage_at_period_end(
    model: Any,
    frame: pd.DataFrame,
    *,
    reader: Any,
    file: str,
    id_column: str,
) -> pd.DataFrame:
    """One stage per feature per stress period: the one that ENDS it (ledgers 331, 336, 348).

    The stage-file route of a lake or reach target. A LAK or SFR stage table holds
    one row per feature-cell per saved step, labelled from the stage file's own
    clock -- which MF6 writes when it saves HEAD, not BUDGET (ledger 301). The
    period-end record is read off that file and checked against the TDIS
    (``period_end.period_ends``), so a period whose end was not saved is refused
    rather than read at an earlier step; the rows of that record are then reduced
    to the feature's one stage (a stage is not added, ledger 297). Both families
    call this: it was ``SfrStageTargets``' alone until ledger 348, while the LAK
    table, keyed by saved record, kept every record and a lake target keyed by
    period read a step inside an EARLIER period.

    Parameters
    ----------
    model : SimulationBase
        The run model the table was read from.
    frame : pandas.DataFrame
        ``model.lak.results.stage.get()`` or ``model.sfr.results.stage.get()``:
        ``kstpkper``, ``per``, ``id_column`` and ``stage``. A table without
        ``kstpkper`` (a test double's) keeps every row.
    reader : callable
        Takes no arguments and returns the stage file's open FloPy reader
        (``lambda: model.fl.lak.output.stage()``); called only when the model has
        a TDIS.
    file : str
        What the file is, for the refusal (``"the LAK stage file"``).
    id_column : str
        The feature id column, ``"lake"`` or ``"reach"``.

    Returns
    -------
    pandas.DataFrame
        ``per``, ``id_column``, ``stage``: one row per feature per period, ``NaN``
        where it has no stage.

    Raises
    ------
    ValueError
        When a period the stage file saved lacks its final step -- naming the
        period, the steps saved and ``SAVE HEAD LAST`` -- or a feature's rows in one
        record disagree.
    """

    ends = period_ends(
        model,
        reader=reader,
        file=file,
        saves="HEAD",
        fallback=frame["kstpkper"] if "kstpkper" in frame.columns else (),
    )
    return _aggregate_named_series(
        frame, id_column=id_column, value_column="stage", kind="elevation", ends=ends
    )


def _one_value_per_group(grouped, *, id_column: str, value_column: str) -> pd.Series:
    """The one number each group's rows are copies of; refuse a group whose copies differ.

    Taken as the FIRST copy, so it is the table's number to the bit; for copies
    that agree that is the registry's ``"mean"`` (ledger 331). A group that is
    ``NaN`` throughout stays ``NaN``.

    Parameters
    ----------
    grouped : pandas.core.groupby.SeriesGroupBy
        The value column grouped by ``(per, id)``.
    id_column : str
        The id level's name, used to name an offending group (``"lake"``).
    value_column : str
        The value's name, used in the message (``"stage"``).

    Returns
    -------
    pandas.Series
        One value per ``(per, id)``, indexed by both.

    Raises
    ------
    ValueError
        If any group's copies span more than ``_COPY_RTOL`` of the value (floored
        at 1), or some hold a value and some do not; the message names up to
        five such groups by id, ``per`` and spread.
    """

    first = grouped.first()
    low = grouped.min()
    high = grouped.max()
    rows = grouped.size()
    held = grouped.count()
    spread = high - low
    tolerance = _COPY_RTOL * np.maximum(1.0, high.abs().fillna(0.0))
    partial = (held > 0) & (held < rows)
    differ = (spread > tolerance) | partial
    if differ.any():
        offenders = []
        for (per, feature) in differ[differ].index[:5]:
            if partial.loc[(per, feature)]:
                offenders.append(
                    f"{id_column} {feature} at per={per}: {held.loc[(per, feature)]} of "
                    f"{rows.loc[(per, feature)]} rows hold a {value_column}"
                )
            else:
                offenders.append(
                    f"{id_column} {feature} at per={per}: {rows.loc[(per, feature)]} rows "
                    f"span {low.loc[(per, feature)]:.10g}..{high.loc[(per, feature)]:.10g} "
                    f"(spread {spread.loc[(per, feature)]:.6g})"
                )
        more = int(differ.sum()) - len(offenders)
        listing = "; ".join(offenders) + (f"; and {more} more" if more else "")
        raise ValueError(
            f"The rows of one {id_column} in one saved record should be copies of ONE "
            f"{value_column} -- the {id_column}'s own, repeated on each cell it connects "
            f"to -- and these are not: {listing}. A {value_column} is not added "
            "(ledger 297), and averaging values that differ would invent one "
            "(ledger 331)."
        )
    return first


def _default_compare_stats(frame: pd.DataFrame, *, residual_column: str = "residual") -> pd.DataFrame:
    """Return basic residual statistics for a compare-style frame."""

    data = frame.dropna(subset=[residual_column]).copy()
    if data.empty:
        return pd.DataFrame([{"n": 0, "mean_error": np.nan, "mae": np.nan, "rmse": np.nan}])
    residual = data[residual_column].to_numpy(dtype=float)
    return pd.DataFrame(
        [
            {
                "n": int(len(data)),
                "mean_error": float(np.mean(residual)),
                "mae": float(np.mean(np.abs(residual))),
                "rmse": float(np.sqrt(np.mean(residual**2))),
            }
        ]
    )


def _normalize_compare_time(frame: pd.DataFrame) -> pd.DataFrame:
    """Return one canonical ``time`` column after target/result merges."""

    result = frame.copy()
    if "time" not in result.columns:
        candidates = [
            result[column]
            for column in ("time_x", "time_y")
            if column in result.columns
        ]
        if candidates:
            time = candidates[0].copy()
            for candidate in candidates[1:]:
                time = time.combine_first(candidate)
            result["time"] = time
    return result.drop(columns=["time_x", "time_y"], errors="ignore")


def _in_grid_crs(features, model):
    """``features`` in the GRID's CRS -- the same object when nothing converts.

    The grid decides where a cell is, so a feature set read in whatever CRS its
    file carries is converted before it is matched. Returning the SAME object
    when there is nothing to convert is what lets a caller tell the two apart
    without comparing CRSs itself.
    """

    grid_crs = getattr(getattr(model, "vor", None), "crs", None)
    if (
        not isinstance(features, gpd.GeoDataFrame)
        or grid_crs is None
        or features.crs is None
        or features.crs == grid_crs
    ):
        return features
    return features.to_crs(grid_crs)


def _placement_note(features, model, *, label: str = "targets") -> str:
    """Where the features are and where the grid is -- what an empty match is read with.

    ``label`` names them in the sentence (``"targets"``, ``"zones"``). The bounds
    are the SOURCE frame's own, in the CRS the sentence names, because a mismatch
    is two boxes and reprojecting one of them first would hide it. The wording
    itself lives in ``grid.geometry.where_they_are``, so the GeoPackage source,
    the selector and these observations all say it the same way.
    """

    return where_they_are(
        features,
        model.vor.gdf_vorPolys,
        noun=label,
        source_crs=features.crs if isinstance(features, gpd.GeoDataFrame) else None,
        grid_crs=getattr(model.vor, "crs", None),
    )


def _report_missed_features(
    owner,
    missed: list[str],
    *,
    total: int,
    subject: str,
    nothing_left: str,
    where,
    log,
) -> None:
    """Ledger 211 for observation targets: raise on a total miss, warn on a partial one.

    ``subject`` is the singular noun (``"head target"``, ``"DRN zone"``) and
    ``nothing_left`` finishes "so there is nothing to ..."; ``total`` is how many
    features were looked for, counted the same way ``missed`` is. A total miss
    means the features and the model are in different places, so it raises rather
    than handing back an empty set whose comparison is silently about nothing. A
    PARTIAL miss can be deliberate -- one file serving several grids -- so it
    warns, by name, through ``warnings.warn`` as well as ``log`` (ledger 342):
    ``import myflopy`` installs a ``NullHandler``, so a log line alone reaches
    nobody in a default session -- ledger 288's finding, for the ``.gpkg``
    sources. The message names the features, so two different misses are two
    warnings rather than one Python deduplicates.

    ``where`` is CALLED, and only for the raise, since the warning does not use
    it: the placement note scans both extents, and the deduped path below returns
    without ever needing the sentence. ``log`` is required rather than defaulted
    to this module's logger, so a miss is always reported under the name of the
    module that found it.

    ``owner`` carries the once-only set: every consumer resolves the same features
    again, so without it one miss is said once per call. Pass ``None`` where there
    is nothing to remember it on, and every call warns.
    """

    listing = ", ".join(missed[:10])
    if len(missed) > 10:
        listing += f", and {len(missed) - 10} more"
    if len(missed) == total:
        raise ValueError(
            f"No {subject} lands on an active cell, so there is nothing to "
            f"{nothing_left}: {listing}. {where()}"
        )
    if owner is not None:
        warned = getattr(owner, "_warned_misses", None)
        if warned is None:
            warned = set()
            owner._warned_misses = warned
        if listing in warned:
            return
        warned.add(listing)
    message = f"{len(missed)} of {total} {subject}s land on no active cell and are left out: {listing}"
    log.warning(message)
    warnings.warn(message, UserWarning, stacklevel=4)

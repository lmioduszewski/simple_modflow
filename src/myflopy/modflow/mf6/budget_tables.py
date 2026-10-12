"""Table-building helpers for :mod:`myflopy.modflow.mf6.budget`.

The functions in this module keep raw binary budget access, dataframe shaping,
and package-output reshaping separate from the public wrapper classes.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from myflopy.modflow.mf6.grid.geometry import (
    footprint_cells,
    read_footprints,
    refuse_footprint_misses,
)
from myflopy.modflow.mf6.package_explorer_utils import _model_budget_record_frame
from myflopy.modflow.mf6.package_registry import _PACKAGE_EXPLORER_SPECS, one_package

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase


def budget_types(model: SimulationBase) -> list[str]:
    """Return record types present in the model budget file.

    Parameters
    ----------
    model
        Parent model providing the raw MF6 budget reader.
    """

    return [record.astype(str) for record in model._get_budget_reader().get_unique_record_names()]


def raw_budget(model: SimulationBase, gwf_package: str | None = None):
    """Return the raw model budget reader or one package's raw record list.

    Parameters
    ----------
    model
        Parent model providing access to the groundwater budget file.
    gwf_package
        Optional groundwater package filter such as ``"drn"`` or ``"rch"``.
    """

    budget_reader = model._get_budget_reader()
    if gwf_package is None:
        return budget_reader
    return budget_reader.get_data(text=gwf_package)


def _cell_based_budget_packages() -> frozenset[str]:
    """Packages whose budget ``node2`` is zero-based, from the package registry.

    Derived rather than hardcoded: a hardcoded ``{"drn", "ghb", "rch"}`` silently
    left ``chd``/``riv``/``wel``/``evt`` on MF6's raw 1-based node ids for as long
    as those packages have existed (see ``tests/test_budget_node_basing.py``).

    Scope narrowed 2026-07-27. This used to govern the ``node`` column too, on
    the belief -- stated in the registry and in this docstring -- that the
    advanced packages "carry feature ids, not model cells". **That was measured
    false.** In the *model* budget file every record's ``node`` is a 1-based
    model cell, including ``SFR``/``LAK``/``UZF-GWRCH`` (it is ``node2`` that
    holds the feature id there; the feature-first layout belongs to the separate
    *package-output* budget file). Excluding those three therefore handed back
    node ids exactly one cell high. ``node`` is now zero-based once, universally,
    in :func:`~myflopy.modflow.mf6.package_budget._model_budget_record_frame`,
    and this set governs only the package-specific ``node2``.

    That leaves the descriptor field ``zero_base_budget_nodes`` misnamed -- it
    now answers a ``node2`` question. Renaming it moves ``api_snapshot.json`` and
    the descriptor payoff count, so it is deferred (ledger 93).
    """

    return frozenset(
        name
        for name, spec in _PACKAGE_EXPLORER_SPECS.items()
        if spec.zero_base_budget_nodes
    )


def _zero_base_budget_frame(frame: pd.DataFrame, gwf_package: str) -> pd.DataFrame:
    """Normalize a package budget's ``node2`` column to zero-based indexing.

    The ``node`` index arrives already zero-based from
    :func:`~myflopy.modflow.mf6.package_budget._model_budget_record_frame`, which
    is the single place MF6's 1-based cell ids are shifted. Callers must NOT
    re-normalize -- doing so is what made ``group.bud('drn')`` come back a cell
    low (fixed 2026-07-18).

    Notes
    -----
    Only the packages listed by :func:`_cell_based_budget_packages` have a
    ``node2`` that myflopy zero-bases; for SFR/LAK/UZF it is a feature id and is
    left exactly as MF6 wrote it.
    """

    normalized = frame.copy()
    if gwf_package in _cell_based_budget_packages() and "node2" in normalized.columns:
        normalized["node2"] = pd.to_numeric(normalized["node2"], errors="coerce") - 1
    return normalized


def budget_df(model: SimulationBase, gwf_package: str) -> pd.DataFrame:
    """Return one concatenated dataframe for a package budget across all periods.

    Parameters
    ----------
    model
        Parent model supplying the raw budget file and stress-period metadata.
    gwf_package
        Groundwater package whose budget should be tabulated.

    Returns
    -------
    pandas.DataFrame
        Concatenated package budget with zero-based node indexing where
        appropriate.
    """

    if gwf_package is None:
        raise ValueError("gwf_package is required to build a budget dataframe.")

    records = raw_budget(model, gwf_package)
    # Label each record from the SAME selection that produced it, not from the
    # file's unique save times. A model with TWO packages of one type -- two RIV6
    # entries in the name file, say -- writes TWO records per saved time, so a
    # 38-period model held 76 RIV records against 38 unique times and this raised
    # IndexError on record 38. Ledger 301, same shape as package_output_df's.
    kstpkper_all = package_output_kstpkper(
        model._get_budget_reader(), bud_type=gwf_package
    )
    frames = []
    for per, record in enumerate(records):
        # Shared with the modern explorer path so the two cannot disagree about
        # node basing again, and so imeth=1 full-array records (STO-SS,
        # STORAGE-AQUEOUS, ...) build a real per-cell table instead of raising
        # pandas' opaque "Must pass 2-d input".
        frame = _model_budget_record_frame(model, record, budget_text=gwf_package)
        frame["kstpkper"] = [kstpkper_all[per] for _ in range(len(frame))]
        frames.append(frame)

    combined = pd.concat(frames)
    combined = combined.set_index(["node", "kstpkper"])
    return _zero_base_budget_frame(combined, gwf_package.lower())


def coerce_plot_times(times) -> pd.DatetimeIndex | None:
    """Normalize optional plotting times to a ``DatetimeIndex``.

    Parameters
    ----------
    times
        ``None``, an existing ``DatetimeIndex``, or an iterable coercible to
        datetimes.
    """

    if times is None:
        return None
    if isinstance(times, pd.DatetimeIndex):
        return times
    return pd.DatetimeIndex(times)


def budget_obs_df(
    budget,
    *,
    model: SimulationBase,
    shp_gpkg,
    q: str = "q",
    name_field: str = "name",
    multiplier: float = 24 * 60 * 60,
) -> pd.DataFrame:
    """Sum one package's budget flows over observation areas, per saved time.

    Every rule here is a shared one (ledger 476): the areas are read through
    :func:`~myflopy.modflow.mf6.grid.geometry.read_footprints` (reprojected to the
    grid's CRS; a file carrying none is taken to be in it; a multi-part polygon is
    one area), a cell belongs to an area it COVERS by area
    (:func:`~myflopy.modflow.mf6.grid.geometry.footprint_cells`), an area on no cell
    is refused naming both frames, every record on a cell is ADDED (a flow is
    extensive, ledger 297), and the rows are the BUDGET file's own saved times
    (ledger 301), never the head file's.

    Parameters
    ----------
    budget : Budget
        The package accessor (``model.bud("drn")``); its :attr:`Budget.df` holds the
        flows, indexed by ``(node, kstpkper)``.
    model : SimulationBase
        The model whose grid (``model.vor``) the areas are matched against.
    shp_gpkg : str or pathlib.Path
        A polygon file (its first layer) whose rows are the observation areas.
    q : str, default "q"
        The column of :attr:`Budget.df` to sum; one it lacks raises ``KeyError``
        naming the columns it has.
    name_field : str, default "name"
        The attribute column naming each area; a missing column raises
        ``ValueError``.
    multiplier : float, default 86400
        The DIVISOR applied to every sum.

    Returns
    -------
    pandas.DataFrame
        Indexed by the budget file's own ``kstpkper``, in file order: ``total``
        (every record of the package), then one column per area. Flows are negated,
        so water leaving the aquifer reads positive.
    """

    frame = budget.df
    if q not in frame.columns:
        raise KeyError(
            f"q={q!r} is not a column of this budget; it holds {list(frame.columns)}."
        )
    flows = frame[q].astype(float) * -1
    times = list(dict.fromkeys(flows.index.get_level_values("kstpkper")))
    by_time = flows.groupby(level="kstpkper", sort=False)
    table = {"total": by_time.sum().reindex(times) / multiplier}

    polys = model.vor.gdf_vorPolys
    areas = read_footprints(
        Path(shp_gpkg), crs=model.vor.crs, id_field=name_field, noun="observation area"
    )
    area_cells = footprint_cells(polys, areas.geometry)
    refuse_footprint_misses(
        polys,
        areas,
        area_cells,
        noun="observation areas",
        consequence=", so they would read no flow",
        source_crs=areas.attrs.get("source_crs"),
        grid_crs=model.vor.crs,
    )
    nodes = flows.index.get_level_values("node")
    for label, cells in zip(areas.index, area_cells, strict=True):
        inside = flows[nodes.isin(cells)]
        table[label] = (
            inside.groupby(level="kstpkper", sort=False).sum().reindex(times, fill_value=0.0)
            / multiplier
        )
    result = pd.DataFrame(table)
    result.index = pd.Index(times, name="kstpkper", tupleize_cols=False)
    return result


def package_output_types(package) -> list[str]:
    """Return record types present in one package output budget file."""

    return [record.astype(str) for record in package.output.budget().get_unique_record_names()]


def package_output_budget(package):
    """Return the raw FloPy package output budget reader."""

    return package.output.budget()


def package_output_kstpkper(reader, *, bud_type: str) -> list[tuple[int, int]]:
    """Return the zero-based ``(kstp, kper)`` of every ``bud_type`` record in ``reader``.

    Taken from the budget file's own ``recordarray`` (1-based on disk) under
    FloPy's own ``get_indices`` mask -- the same selection, in the same file
    order, that ``get_data(text=bud_type)`` returns. So the two zip exactly.
    """

    indices = reader.get_indices(text=bud_type)
    if indices is None:
        return []
    recordarray = reader.recordarray
    return [
        (int(recordarray["kstp"][idx]) - 1, int(recordarray["kper"][idx]) - 1)
        for idx in indices
    ]


def package_output_df(
    model: SimulationBase,
    package,
    *,
    bud_type: str,
):
    """Return one concatenated dataframe for a package-output budget type.

    Each record is labelled with the ``(kstp, kper)`` the package-output budget
    file recorded for it -- never with ``model.kstpkper``, which is the HEAD
    file's clock.

    Notes
    -----
    An advanced package's budget file and the model's dependent-variable file
    are written on two INDEPENDENT MF6 clocks. One OC ``saverecord`` sets
    ``icbcfl`` from ``SAVE BUDGET`` and ``idvsave`` from ``SAVE HEAD``
    (``gwf.f90``), and LAK/SFR gate their BUDGET output on the first while
    gating their STAGE output on the second. Measured: ``SAVE HEAD LAST`` +
    ``SAVE BUDGET ALL`` at ``nstp=5`` writes 2 head records against 10 LAK
    budget records, and ``SAVE HEAD STEPS 1`` + ``SAVE BUDGET STEPS 3`` writes
    two of each naming DIFFERENT time steps. Pairing the two positionally --
    which this function did until ledger 301 -- crashed on the first and
    silently mislabelled on the second.

    Parameters
    ----------
    model
        Parent model. Not consulted for labelling; see the note above.
    package
        FloPy package object exposing ``output.budget()``.
    bud_type
        Output record type to concatenate across saved times.
    """

    del model  # the package budget file carries its own clock; see Notes
    reader = package.output.budget()
    records = reader.get_data(text=bud_type)
    labels = package_output_kstpkper(reader, bud_type=bud_type)
    if len(records) != len(labels):
        # Unreachable through FloPy's own mask (get_data and get_indices select
        # from one recordarray); a guard against a reader whose record order
        # stops tracking its index.
        raise ValueError(
            f"{bud_type!r}: the package budget file yielded {len(records)} "
            f"records but {len(labels)} (kstp, kper) labels"
        )

    frames = []
    for record, label in zip(records, labels, strict=True):
        frame = pd.DataFrame(record)
        frame["kstpkper"] = [label for _ in range(len(frame))]
        frames.append(frame)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames)


def lake_stage_array(model: SimulationBase):
    """Return LAK stage output reshaped to ``(nper, nlakes)``."""

    lak = one_package(model, "lak")  # the model's only LAK, or refused by name (ledger 352)
    nlakes = lak.nlakes.data
    stages = lak.output.stage()
    return stages.get_alldata().reshape(-1, nlakes)


def sfr_stage_df(model: SimulationBase) -> pd.DataFrame:
    """Return SFR stage output as a reach-by-period dataframe."""

    sfr = one_package(model, "sfr")  # the model's only SFR, or refused by name (ledger 352)
    nreaches = sfr.nreaches.data
    stages = sfr.output.stage().get_alldata()
    data = stages.reshape(model.nper, nreaches).transpose()
    frame = pd.DataFrame(data)
    frame.index.name = "reaches"
    return frame

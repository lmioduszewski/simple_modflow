"""Table-building helpers for :mod:`myflopy.modflow.mf6.budget`.

The functions in this module keep raw binary budget access, dataframe shaping,
and package-output reshaping separate from the public wrapper classes.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from myflopy.modflow.mf6.boundaries import Boundaries
from myflopy.modflow.utils.datatypes.readers import read_shp_gpkg

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase

idxx = pd.IndexSlice


def budget_types(model: "SimulationBase") -> list[str]:
    """Return record types present in the model budget file.

    Parameters
    ----------
    model
        Parent model providing the raw MF6 budget reader.
    """

    return [record.astype(str) for record in model._get_budget_reader().get_unique_record_names()]


def raw_budget(model: "SimulationBase", gwf_package: str | None = None):
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


def _zero_base_budget_frame(frame: pd.DataFrame, gwf_package: str) -> pd.DataFrame:
    """Normalize package budget node columns to zero-based indexing when needed.

    Notes
    -----
    Only groundwater budget packages whose node ids map directly to model cells
    are normalized here. The level names are preserved so downstream code sees
    stable ``node`` / ``kstpkper`` columns after ``reset_index()``.
    """

    normalized = frame.copy()
    if gwf_package in {"drn", "ghb", "rch"}:
        index_names = normalized.index.names
        new_index = pd.MultiIndex.from_tuples([(int(node) - 1, kstpkper) for node, kstpkper in normalized.index])
        new_index = new_index.set_names(index_names)
        normalized.index = new_index
        if "node2" in normalized.columns:
            normalized["node2"] = pd.to_numeric(normalized["node2"], errors="coerce") - 1
    return normalized


def budget_df(model: "SimulationBase", gwf_package: str) -> pd.DataFrame:
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
    kstpkper_all = model._get_budget_kstpkper()
    frames = []
    for per, record in enumerate(records):
        frame = pd.DataFrame(record)
        frame["kstpkper"] = [kstpkper_all[per] for _ in range(len(frame))]
        frames.append(frame)

    combined = pd.concat(frames)
    node_column = str(combined.columns[0])
    combined = combined.set_index([combined.columns[0], "kstpkper"])
    combined.index = combined.index.set_names([node_column, "kstpkper"])
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
    model: "SimulationBase",
    shp_gpkg: Path,
    q: str = "q",
    name_field: str = "name",
    multiplier: float = 24 * 60 * 60,
) -> pd.DataFrame:
    """Aggregate package budget flows by observation polygon and stress period.

    Parameters
    ----------
    budget
        Parent :class:`~myflopy.modflow.mf6.budget.Budget` accessor.
    model
        Parent model providing the grid and stress-period selectors.
    shp_gpkg
        Polygon dataset defining the observation areas.
    q
        Budget value column to aggregate, usually ``"q"``.
    name_field
        Attribute field used to label each observation polygon.
    multiplier
        Unit-conversion divisor applied after summing flows.
    """

    full_idx = range(model.vor.ncpl)
    df = budget.df[~budget.df.index.duplicated(keep="first")]
    per0 = model.kstpkper[0]

    pkg_cells = df.loc[idxx[:, per0], :].index.get_level_values(0).to_list()
    b_obs = read_shp_gpkg(shp_gpkg)
    crs = b_obs.crs.to_epsg()
    b_obs["cells"] = Boundaries(model, model.vor, shp_gpkg=shp_gpkg, crs=crs).intersections

    def drop_non_pkg_cells(cells):
        return [cell for cell in cells if cell in pkg_cells]

    b_obs["cells"] = b_obs["cells"].apply(drop_non_pkg_cells)

    budget_dict: dict[str, dict] = {"total": {}}
    for obs in b_obs.iterrows():
        budget_dict[obs[1].loc[name_field]] = {}
        obs_cells = obs[1].loc["cells"]
        for kstpkper in model.kstpkper:
            try:
                flows = df.loc[idxx[:, kstpkper], :][q].droplevel(1) * -1
                budget_dict["total"][kstpkper] = flows.sum() / multiplier
                flows = flows.reindex(full_idx, fill_value=0)
                budget_dict[obs[1].loc[name_field]][kstpkper] = flows.loc[obs_cells].sum() / multiplier
            except KeyError:
                continue

    return pd.DataFrame.from_dict(budget_dict).reset_index(drop=True)


def package_output_types(package) -> list[str]:
    """Return record types present in one package output budget file."""

    return [record.astype(str) for record in package.output.budget().get_unique_record_names()]


def package_output_budget(package):
    """Return the raw FloPy package output budget reader."""

    return package.output.budget()


def package_output_df(
    model: "SimulationBase",
    package,
    *,
    bud_type: str,
):
    """Return one concatenated dataframe for a package-output budget type.

    Parameters
    ----------
    model
        Parent model supplying stress-period metadata.
    package
        FloPy package object exposing ``output.budget()``.
    bud_type
        Output record type to concatenate across periods.
    """

    records = package.output.budget().get_data(text=bud_type)
    if len(records) != len(model.kstpkper):
        raise ValueError("length of periods and budget array list do not match")

    frames = []
    for idx, per in enumerate(model.kstpkper):
        frame = pd.DataFrame(records[idx])
        frame["kstpkper"] = [per for _ in range(len(frame))]
        frames.append(frame)
    return pd.concat(frames)


def lake_stage_array(model: "SimulationBase"):
    """Return LAK stage output reshaped to ``(nper, nlakes)``."""

    nlakes = model.lak.nlakes.data
    stages = model.lak.output.stage()
    return stages.get_alldata().reshape(-1, nlakes)


def sfr_stage_df(model: "SimulationBase") -> pd.DataFrame:
    """Return SFR stage output as a reach-by-period dataframe."""

    nreaches = model.sfr.nreaches.data
    stages = model.sfr.output.stage().get_alldata()
    data = stages.reshape(model.nper, nreaches).transpose()
    frame = pd.DataFrame(data)
    frame.index.name = "reaches"
    return frame

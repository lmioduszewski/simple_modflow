"""Observation-oriented helpers for :mod:`myflopy.modflow.mf6.headsplus`.

These functions keep observation lookups and table shaping separate from the
core binary-head reader so ``HeadsPlus`` can focus on reading and exposing
model-aware head data.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from myflopy.modflow.utils.validators import valid_list_of_cell_ints

if TYPE_CHECKING:
    from myflopy.modflow.mf6.headsplus import HeadsPlus

idxx = pd.IndexSlice


def sort_dict_by_keys(dict_to_sort: dict | None = None) -> dict:
    """Return ``dict_to_sort`` ordered by sorted key values."""

    if dict_to_sort is None:
        return {}
    sorted_keys = sorted(dict_to_sort.keys())
    return {key: dict_to_sort[key] for key in sorted_keys}


def get_obs_cells(
    heads: "HeadsPlus",
    locs: Path | None,
    crs: str | None = None,
    loc_name_field: str = "ExploName",
) -> dict:
    """Return the single model cell containing each observation point.

    Parameters
    ----------
    heads
        Parent :class:`HeadsPlus` instance providing ``vor`` and default CRS
        information.
    locs
        Path to a point shapefile or GeoPackage containing observation
        locations. If omitted, ``heads.obs_path`` is used.
    crs
        Optional CRS override for the input file. Defaults to ``heads.crs``.
    loc_name_field
        Attribute field containing observation names.
    """

    locs = heads.obs_path if locs is None else locs
    crs = heads.crs if crs is None else crs
    if locs is None:
        raise ValueError("No observation path found. Provide locs or set heads.obs_path first.")

    obs_dict = heads.vor.get_vor_cells_as_dict(
        locs=locs,
        crs=crs,
        predicate="contains",
        loc_name_field=loc_name_field,
    )
    obs_dict = sort_dict_by_keys(obs_dict)
    obs_dict = {key: value for key, value in obs_dict.items() if len(value) > 0}

    normalized: dict[str, int] = {}
    for obs_name, cell_ids in obs_dict.items():
        assert len(cell_ids) == 1, f"more than one cell found for {obs_name}. Fix to make it one cell"
        normalized[obs_name] = cell_ids[0][0]
    return normalized


def get_obs_heads(
    heads: "HeadsPlus",
    locs: Path | list[int] | None = None,
    crs: str | None = None,
    loc_name_field: str = "ExploName",
    long_format: bool = False,
) -> pd.DataFrame:
    """Return observation heads for explicit cells or point locations.

    Parameters
    ----------
    heads
        Parent :class:`HeadsPlus` instance supplying the heads dataframe and
        observation metadata cache.
    locs
        Either a path to point locations or a list of zero-based cell ids.
    crs
        Optional CRS override for point inputs.
    loc_name_field
        Attribute field containing point names when ``locs`` is spatial data.
    long_format
        If ``True``, return a long/tidy dataframe indexed by location, layer,
        and stress period.
    """

    crs = heads.crs if crs is None else crs

    if isinstance(locs, Path):
        heads.obs_path = locs

    if isinstance(locs, list):
        if valid_list_of_cell_ints(heads.model, locs):
            obs_cells = locs
            heads._obs = locs
        else:
            raise ValueError("locs must be a valid list of zero-based cell ids.")
    elif not heads.obs and heads.obs_path is not None:
        resolved = get_obs_cells(heads, heads.obs_path, crs=crs, loc_name_field=loc_name_field)
        heads._obs.update(resolved)
        obs_cells = list(heads.obs.values())
    else:
        raise ValueError("No observation path found or observation cells provided.")

    all_heads = heads.all_heads.copy()
    obs_heads = all_heads.loc[idxx[:, :, obs_cells], :]
    obs_reset_idx = obs_heads.reset_index()
    obs_heads = obs_reset_idx.pivot(
        index=["layer", "kstpkper"],
        columns="cell",
        values="elev",
    )

    new_cols = []
    for col in obs_heads.columns:
        if isinstance(heads.obs, dict):
            obs_name = next(key for key, value in heads.obs.items() if value == col)
        else:
            obs_name = col
        new_cols.append(obs_name)
    obs_heads.columns = new_cols
    obs_heads = obs_heads[sorted(obs_heads.columns)]

    if long_format:
        obs_heads = (
            obs_heads.melt(ignore_index=False, var_name="locs", value_name="elev")
            .reset_index()
            .set_index(["locs", "layer", "kstpkper"])
        )

    return obs_heads

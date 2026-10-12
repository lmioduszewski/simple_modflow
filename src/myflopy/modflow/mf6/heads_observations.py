"""Observation-oriented helpers for :mod:`myflopy.modflow.mf6.headsplus`.

These functions keep observation lookups and table shaping separate from the
core binary-head reader so ``HeadsPlus`` can focus on reading and exposing
model-aware head data.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd
import pyproj

from myflopy.modflow.utils.datatypes.readers import read_vector_layer
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
    heads: HeadsPlus,
    locs: str | os.PathLike | None,
    crs: str | None = None,
    loc_name_field: str = "ExploName",
) -> dict:
    """Return the single model cell holding each observation point.

    Parameters
    ----------
    heads
        Parent :class:`HeadsPlus` instance providing ``vor``.
    locs
        Path to a point file (its first layer). ``None`` reads ``heads.obs_path``.
    crs
        The CRS a file carrying NONE is in (ledger 476); ``None`` takes such a
        file to be in the grid's. A file carrying its own CRS is reprojected to the
        grid's, and a ``crs`` naming another one raises ``ValueError``.
    loc_name_field
        Attribute field containing observation names.
    """

    locs = heads.obs_path if locs is None else locs
    if locs is None:
        raise ValueError("No observation path found. Provide locs or set heads.obs_path first.")
    frame = read_vector_layer(Path(locs), None)
    if crs is not None:
        named = pyproj.CRS.from_user_input(crs)
        if frame.crs is None:
            frame = frame.set_crs(named)
        elif frame.crs != named:
            raise ValueError(
                f"crs={crs!r} names the frame of a points file that carries none, and "
                f"{Path(locs).name} carries its own ({frame.crs.to_string()}), which it "
                "is reprojected from. Drop crs=, or fix the file's CRS."
            )
    # The grid's own selector: the shared cell rule (a point is the one cell holding
    # it, ledgers 210-212), the frame put in the grid's CRS, a miss logged by name.
    cells = heads.vor.get_vor_cells_as_series(frame, name_field=loc_name_field).to_dict()
    several = {name: list(found) for name, found in cells.items() if len(found) > 1}
    if several:
        listed = ", ".join(f"{name!r} ({len(found)} cells)" for name, found in several.items())
        raise ValueError(
            f"An observation location is ONE cell, and {listed} reach several -- a "
            "line or polygon, not a point. Give each location as a point."
        )
    return {name: int(found[0]) for name, found in sort_dict_by_keys(cells).items() if found}


def get_obs_heads(
    heads: HeadsPlus,
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
        A path (``str`` or path-like) to point locations, resolved now, or a list
        of zero-based cell ids; ``None`` reuses the resolved locations, else
        ``heads.obs_path``.
    crs
        The CRS of a points file that carries none (see :func:`get_obs_cells`).
    loc_name_field
        Attribute field containing point names when ``locs`` is spatial data.
    long_format
        If ``True``, return a long/tidy dataframe indexed by location, layer,
        and stress period.

    Returns
    -------
    pandas.DataFrame
        Every saved time step, every layer. A cell MF6 wrote as dry (-1e30) or
        inactive (+1e30) reads ``NaN`` (ledger 495), as in
        :meth:`HeadsPlus.long`.
    """

    if isinstance(locs, (str, os.PathLike)):
        # A path given HERE is read now, whatever the reader resolved before: a
        # str was not recognised as one and the constructor's file was read in
        # its place, and a second path raised (ledger 476).
        heads.obs_path = Path(locs)
        heads._obs = get_obs_cells(heads, heads.obs_path, crs=crs, loc_name_field=loc_name_field)
    elif isinstance(locs, list):
        if not valid_list_of_cell_ints(heads.model, locs):
            raise ValueError("locs must be a valid list of zero-based cell ids.")
        heads._obs = {cell: cell for cell in locs}
    elif locs is not None:
        raise TypeError(
            f"locs must be a path to a points file or a list of cell ids, not "
            f"{type(locs).__name__}."
        )
    elif not heads.obs:
        if heads.obs_path is None:
            raise ValueError("No observation path found or observation cells provided.")
        heads._obs = get_obs_cells(heads, heads.obs_path, crs=crs, loc_name_field=loc_name_field)
    named_cells = dict(heads.obs)

    obs_cells = sorted(set(named_cells.values()))
    # `long()`, not `all_heads`: the masked series (`headsplus._is_sentinel`, ledger
    # 315), so a dry (-1e30) or inactive (+1e30) cell reads NaN (ledger 495).
    obs_heads = heads.long().loc[idxx[:, :, obs_cells]].reset_index().pivot(
        index=["layer", "kstpkper"],
        columns="cell",
        values=heads.store_column,
    )
    obs_heads = pd.DataFrame(
        {name: obs_heads[cell] for name, cell in sort_dict_by_keys(named_cells).items()},
        index=obs_heads.index,
    )

    if long_format:
        obs_heads = (
            obs_heads.melt(ignore_index=False, var_name="locs", value_name="elev")
            .reset_index()
            .set_index(["locs", "layer", "kstpkper"])
        )

    return obs_heads

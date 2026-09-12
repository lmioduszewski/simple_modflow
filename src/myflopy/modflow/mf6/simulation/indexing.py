from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from myflopy._logging import get_logger
from myflopy.modflow.mf6.grid.geometry import (
    miss_reason,
    selected_cells,
    where_they_are,
)
from myflopy.modflow.utils.datatypes.readers import (
    assign_voronoi_cells_to_layers,
    read_shp_gpkg,
)

logger = get_logger(__name__)


def build_idomain(vor, idomain_path: Path):
    """
    Build the MF6 idomain array from a layer-tagged polygon dataset.

    Each feature's ``layer`` tag says which layers it switches OFF, and the cells
    it reaches come from the shared coverage rule: a polygon makes inactive the
    cells it COVERS by area, not every cell it happens to touch. The whole file is
    resolved in ONE selector call, so a single row landing off the grid is a
    partial miss and not a whole-input one.

    Ledger 211: a file whose features land on NO cell raises. It used to come back
    as a domain with every cell active in every layer -- the model then runs to
    normal termination with the idomain doing nothing -- and the cause is almost
    always a location or a CRS, so the message carries both extents and both CRSs.
    A file that lands only PARTLY keeps what landed; the selector names the rest
    in a warning.
    """
    idomain_gdf = read_shp_gpkg(idomain_path)
    idomain_gdf, col_names = assign_voronoi_cells_to_layers(
        idomain_gdf, vor, return_col_names=True
    )

    idomain_by_layer_dict = {}
    for col in col_names:
        # Every non-null value is the LIST of cells that feature reaches, so
        # `selected_cells` -- the one copy of that flattening -- is the whole job.
        idomain_by_layer_dict[int(col)] = selected_cells(idomain_gdf[col].dropna())

    # A column that is empty is a layer no feature named; a file where EVERY
    # column is empty landed nowhere at all (ledger 211), which used to answer
    # "every cell active in every layer" and let the model run regardless.
    if col_names and not any(idomain_by_layer_dict.values()):
        polys = vor.gdf_vorPolys
        reasons = [
            f"{index!r} ({miss_reason(polys, geometry, target='grid cell')})"
            for index, geometry in idomain_gdf.geometry.items()
        ]
        listing = ", ".join(reasons[:10])
        if len(reasons) > 10:
            listing += f", and {len(reasons) - 10} more"
        raise ValueError(
            f"No feature in '{getattr(idomain_path, 'name', idomain_path)}' lands on a "
            "grid cell, so every cell would stay active and the idomain would do "
            f"nothing. {listing}. "
            + where_they_are(
                idomain_gdf,
                polys,
                source_crs=idomain_gdf.crs,
                no_crs_subject="the file has NO CRS, so its coordinates",
            )
        )

    idomain = []
    for layer in range(vor.nlay):
        inactive_cells = idomain_by_layer_dict.get(layer, set())
        layer_idomain = [0 if idx in inactive_cells else 1 for idx in range(vor.ncpl)]
        idomain.append(layer_idomain)

    return idomain, idomain_gdf


def coerce_per_dates(per_dates):
    """Coerce a list of period dates to a pandas ``DatetimeIndex`` (asserting an existing index)."""

    if per_dates is not None:
        if isinstance(per_dates, list):
            try:
                per_dates = pd.to_datetime(per_dates)
            except ValueError:
                logger.warning(
                    'could not parse the given period dates as datetimes; '
                    'leaving them as provided'
                )
        else:
            assert isinstance(per_dates, pd.DatetimeIndex), 'Cannot recognize valid dates in provide period dates'
    return per_dates


def build_model_times(gwf):
    """Map each ``(timestep, period)`` to its cumulative simulation time (``totim``)."""

    times = gwf.modeltime.totim
    steps = gwf.modeltime.kper_kstp
    return {(step[1], step[0]): float(time) for step, time in zip(steps, times, strict=False)}


def build_ncpl_arr(modelgrid) -> np.ndarray:
    """The per-layer cell count as an ``(nlay,)`` array (broadcasting a scalar ``ncpl``)."""

    return (
        np.full(modelgrid.nlay, modelgrid.ncpl, dtype=int)
        if isinstance(modelgrid.ncpl, (int, np.integer))
        else np.asarray(modelgrid.ncpl, dtype=int)
    )


def build_offsets(ncpl_arr: np.ndarray) -> np.ndarray:
    """Cumulative node offsets per layer (a leading 0), for flat node <-> layer indexing."""

    return np.concatenate(([0], np.cumsum(ncpl_arr)))


def build_node_to_lni(ncpl_arr: np.ndarray, offsets: np.ndarray) -> dict[int, tuple[int, int]]:
    """Map each flat node number to its ``(layer, in-layer index)`` pair."""

    out = {}
    for layer, ncpl in enumerate(ncpl_arr):
        start = offsets[layer]
        for idx_in_layer in range(int(ncpl)):
            out[start + idx_in_layer] = (layer, idx_in_layer)
    return out

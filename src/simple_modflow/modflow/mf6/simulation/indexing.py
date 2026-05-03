from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from simple_modflow.modflow.utils.datatypes.readers import (
    assign_voronoi_cells_to_layers,
    read_shp_gpkg,
)


def build_idomain(vor, idomain_path: Path):
    """
    Build the MF6 idomain array from a layer-tagged polygon dataset.
    """
    idomain_gdf = read_shp_gpkg(idomain_path)
    idomain_gdf, col_names = assign_voronoi_cells_to_layers(
        idomain_gdf, vor, return_col_names=True
    )

    idomain_by_layer_dict = {}
    for col in col_names:
        values = idomain_gdf[col].dropna()
        cell_lists = values[values.apply(lambda value: isinstance(value, list))]
        all_cells = sorted(set(cell for sublist in cell_lists for cell in sublist))
        idomain_by_layer_dict[int(col)] = all_cells

    idomain = []
    for layer in range(vor.nlay):
        inactive_cells = idomain_by_layer_dict.get(layer, set())
        layer_idomain = [0 if idx in inactive_cells else 1 for idx in range(vor.ncpl)]
        idomain.append(layer_idomain)

    return idomain, idomain_gdf


def coerce_per_dates(per_dates):
    if per_dates is not None:
        if isinstance(per_dates, list):
            try:
                per_dates = pd.to_datetime(per_dates)
            except ValueError:
                print("Can't convert provided period dates to Pandas DateTime")
        else:
            assert isinstance(per_dates, pd.DatetimeIndex), 'Cannot recognize valid dates in provide period dates'
    return per_dates


def build_model_times(gwf):
    times = gwf.modeltime.totim
    steps = gwf.modeltime.kper_kstp
    return {(step[1], step[0]): float(time) for step, time in zip(steps, times)}


def build_ncpl_arr(modelgrid) -> np.ndarray:
    return (
        np.full(modelgrid.nlay, modelgrid.ncpl, dtype=int)
        if isinstance(modelgrid.ncpl, (int, np.integer))
        else np.asarray(modelgrid.ncpl, dtype=int)
    )


def build_offsets(ncpl_arr: np.ndarray) -> np.ndarray:
    return np.concatenate(([0], np.cumsum(ncpl_arr)))


def build_node_to_lni(ncpl_arr: np.ndarray, offsets: np.ndarray) -> dict[int, tuple[int, int]]:
    out = {}
    for layer, ncpl in enumerate(ncpl_arr):
        start = offsets[layer]
        for idx_in_layer in range(int(ncpl)):
            out[start + idx_in_layer] = (layer, idx_in_layer)
    return out

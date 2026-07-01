"""Regression tests for grid top/bottom overlap adjustment.

`get_vor_cells_as_series` returns one row per geometry, each holding a *list* of
intersected cell ids. `adjust_top_btm_overlaps` must flatten those lists into
individual cells; iterating them directly passed a whole list (or a non-int id)
into `find_adjacent_cells`, which slices `iac[:cell_id]` and raised
"slice indices must be integers".
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import box

import myflopy as mf
from myflopy.layers import Array
from myflopy.modflow.mf6.canonical_example import (
    CanonicalModelConfig,
    rectangular_voronoi,
)
from myflopy.modflow.mf6.grid.geometry import find_adjacent_cells


@pytest.fixture(scope="module")
def vor3():
    """A small 3-layer Voronoi grid with top/botm populated."""
    vor = rectangular_voronoi(CanonicalModelConfig(nrow=8, ncol=8, nlay=3, nper=1))
    ncpl = int(vor.ncpl)
    (
        mf.LayerStack(vor, top=Array(np.full(ncpl, 100.0)), length_units="feet")
        .add("a", thickness=20.0)
        .add("b", thickness=20.0)
        .add("c", thickness=20.0)
        .build(attach=True)
    )
    return vor


def test_find_adjacent_cells_tolerates_nonint_ids(vor3):
    """numpy-int and float-typed cell ids from spatial joins must not break slicing."""
    from_int = find_adjacent_cells(vor3, 5)
    from_np = find_adjacent_cells(vor3, np.int64(5))
    from_float = find_adjacent_cells(vor3, 5.0)
    assert from_int == from_np == from_float
    assert all(isinstance(cell, int) for cell in from_int)


def test_adjust_top_btm_overlaps_flattens_cell_lists(vor3, tmp_path):
    """A polygon selecting several cells is flattened, not passed as a list."""
    xmin, ymin, xmax, ymax = vor3.gdf_vorPolys.total_bounds
    cx, cy = (xmin + xmax) / 2, (ymin + ymax) / 2
    w, h = (xmax - xmin) * 0.25, (ymax - ymin) * 0.25
    gpkg = tmp_path / "selection.gpkg"
    gpd.GeoDataFrame(
        {"geometry": [box(cx - w, cy - h, cx + w, cy + h)]}, crs=vor3.crs
    ).to_file(gpkg)

    # Before the fix this raised "slice indices must be integers".
    out = vor3.adjust_top_btm_overlaps(shp=gpkg, buffer=1, min_sep=1)
    assert out.shape[0] == int(vor3.ncpl)
    assert 0 in out.columns and 1 in out.columns

    # The adjustment must actually change layer-1 bottoms, not silently no-op
    # (the chained-inplace .update() would drop the change under copy-on-write).
    baseline = vor3.reconcile_surfaces(min_sep=1)
    assert not out[1].equals(baseline[1])

"""Regression test: GHB.from_vector tolerates missing optional attribute fields.

`from_polygons` used to read every field (name/elev/height/cond/layer/min_elev)
unconditionally, so a GHB shapefile that only carried the fields it uses (e.g.
an elevation-only boundary with name/elev/cond/layer) raised ``KeyError: 'height'``.
The optional fields (height_over_btm, min_elev, elevation) are now read defensively.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
from shapely.geometry import box

import myflopy as mf
from myflopy.layers import Array
from myflopy.modflow.mf6.canonical_example import (
    CanonicalModelConfig,
    rectangular_voronoi,
)
from myflopy.modflow.mf6.ghb import GHB


def test_ghb_from_vector_tolerates_missing_optional_fields(tmp_path):
    vor = rectangular_voronoi(CanonicalModelConfig(nrow=8, ncol=8, nlay=3, nper=1))
    ncpl = int(vor.ncpl)
    (
        mf.LayerStack(vor, top=Array(np.full(ncpl, 100.0)), length_units="feet")
        .add("a", thickness=20.0)
        .add("b", thickness=20.0)
        .add("c", thickness=20.0)
        .build(attach=True)
    )

    # A GHB strip along the left edge with ONLY name/elev/cond/layer -- no
    # 'height' or 'min_elev' columns (the case that used to KeyError).
    xmin, ymin, xmax, ymax = vor.gdf_vorPolys.total_bounds
    strip = box(xmin, ymin, xmin + (xmax - xmin) * 0.15, ymax)
    shp = tmp_path / "ghb.gpkg"
    gpd.GeoDataFrame(
        {"name": ["west"], "elev": [95.0], "cond": [100.0], "layer": [1],
         "geometry": [strip]},
        crs=vor.crs,
    ).to_file(shp)

    ghb = GHB(vor=vor, shp_gpkg=shp, uid="name")
    spd = ghb.from_vector()  # must not raise KeyError on the missing fields

    assert 0 in spd
    # any produced records use the elevation attribute as the boundary head
    assert all(rec[1] == 95.0 for rec in spd[0])


def test_ghb_from_vector_selects_all_cells_in_polygon_not_just_edges(tmp_path):
    """A GHB *zone* polygon must apply to every cell it covers, not only the
    perimeter cells. `from_polygons` used to hard-code edges_only=True, which
    silently dropped the interior cells of a wide polygon -- e.g. it produced 24
    of the 46 discharge cells the legacy `get_ghb_from_shp` selected, letting the
    modeled water table mound tens of feet where those boundaries were missing.
    """
    vor = rectangular_voronoi(CanonicalModelConfig(nrow=8, ncol=8, nlay=3, nper=1))
    ncpl = int(vor.ncpl)
    (
        mf.LayerStack(vor, top=Array(np.full(ncpl, 100.0)), length_units="feet")
        .add("a", thickness=20.0)
        .add("b", thickness=20.0)
        .add("c", thickness=20.0)
        .build(attach=True)
    )

    # A central block wide enough to contain interior (non-perimeter) cells.
    xmin, ymin, xmax, ymax = vor.gdf_vorPolys.total_bounds
    bx0, by0 = xmin + (xmax - xmin) * 0.30, ymin + (ymax - ymin) * 0.30
    bx1, by1 = xmin + (xmax - xmin) * 0.70, ymin + (ymax - ymin) * 0.70
    shp = tmp_path / "ghb_zone.gpkg"
    gpd.GeoDataFrame(
        {"name": ["zone"], "elev": [95.0], "cond": [100.0], "layer": [1],
         "geometry": [box(bx0, by0, bx1, by1)]},
        crs=vor.crs,
    ).to_file(shp)

    all_cells = GHB(vor=vor, shp_gpkg=shp, uid="name").from_vector()          # default
    edge_cells = GHB(vor=vor, shp_gpkg=shp, uid="name").from_polygons(edges_only=True)

    all_ids = {rec[0] for rec in all_cells[0]}
    edge_ids = {rec[0] for rec in edge_cells[0]}

    # The default (all-cells) is a strict superset of the perimeter-only selection.
    assert edge_ids < all_ids
    assert len(all_ids) > len(edge_ids)  # interior cells are included by default

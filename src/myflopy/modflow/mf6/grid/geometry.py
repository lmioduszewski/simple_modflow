from __future__ import annotations

import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely as shp
from shapely.geometry import LineString, Polygon
from shapely.prepared import prep

from myflopy.modflow.utils.datatypes.readers import read_shp_gpkg


def get_voronoi_polygons(vor) -> list[Polygon]:
    """
    Build Shapely polygons for each Voronoi cell from vertex arrays.
    """
    polygons = []
    for region in vor.iverts:
        vertices = [vor.verts[j] for j in region]
        polygons.append(Polygon(vertices))
    return polygons


def generate_grid_coordinates(vor, grid_spacing: int) -> tuple[list[int], list[int]]:
    """
    Generate evenly spaced coordinates spanning the grid extent.
    """
    if hasattr(vor, 'x_vor') and hasattr(vor, 'y_vor'):
        xs = np.asarray(vor.x_vor)
        ys = np.asarray(vor.y_vor)
    else:
        verts = np.asarray(vor.verts)
        xs = verts[:, 0]
        ys = verts[:, 1]

    x_min = xs.min()
    x_max = xs.max()
    y_min = ys.min()
    y_max = ys.max()
    x_coords = []
    y_coords = []
    for x in range(int(x_min), int(x_max) + grid_spacing, grid_spacing):
        for y in range(int(y_min), int(y_max) + grid_spacing, grid_spacing):
            x_coords.append(x)
            y_coords.append(y)
    return x_coords, y_coords


def get_centroids(vor) -> tuple[list[float], list[float]]:
    """
    Return centroid coordinates for each Voronoi polygon.
    """
    if vor._gdf_vorPolys is None:
        vor.get_gdf_vorPolys()

    centroids = vor._gdf_vorPolys.geometry.centroid
    return centroids.x.tolist(), centroids.y.tolist()


def find_adjacent_polygons(gdf: gpd.GeoDataFrame) -> list[list[int]]:
    """
    Determine which polygons in the grid share an edge.
    """
    adjacent = [[] for _ in range(len(gdf))]
    gdf_sindex = gdf.sindex

    for i, poly in gdf.geometry.items():
        candidates = list(gdf_sindex.intersection(poly.bounds))

        for j in candidates:
            if i == j:
                continue

            shared_edge = poly.intersection(gdf.iloc[j].geometry)
            if isinstance(shared_edge, LineString):
                adjacent[i].append(j)

    return adjacent


def find_adjacent_cells(vor, cell_id: int) -> list[int]:
    """
    Return neighboring cell ids from the DISU connectivity vectors.
    """
    start_index = int(np.sum(vor.iac[:cell_id]))
    num_connections = int(vor.iac[cell_id])
    connections = vor.ja[start_index:start_index + num_connections]
    return [int(cell) for cell in connections if int(cell) != cell_id]


def calculate_distance(gdf: gpd.GeoDataFrame, poly_idx1: int, poly_idx2: int) -> float:
    """
    Calculate centroid-to-shared-face distance for two adjacent polygons.
    """
    poly1 = gdf.iloc[poly_idx1].geometry
    poly2 = gdf.iloc[poly_idx2].geometry
    shared_edge = poly1.intersection(poly2)
    centroid = poly1.centroid
    return centroid.distance(shared_edge)


def shared_face_length(poly1, poly2) -> float:
    """
    Calculate the length of the shared face between adjacent polygons.
    """
    intersection = poly1.intersection(poly2)
    if not isinstance(intersection, LineString):
        return 0
    return intersection.length


def get_gdf_vor_polys(vor, crs=None) -> gpd.GeoDataFrame:
    """
    Build and cache the GeoDataFrame of Voronoi polygons.
    """
    crs = vor.crs if crs is None else crs
    vertices_by_cells = []
    xvertices_by_cells = []
    yvertices_by_cells = []
    polygons = []

    for cell in vor.iverts:
        thiscell = []
        thiscellx = []
        thiscelly = []
        for vidx in cell:
            thisvert = vor.verts[vidx].tolist()
            thiscell.append(thisvert)
            thiscellx.append(thisvert[0])
            thiscelly.append(thisvert[1])
        xvertices_by_cells.append(thiscellx)
        yvertices_by_cells.append(thiscelly)
        vertices_by_cells.append(thiscell)

    for cell in vertices_by_cells:
        polygons.append(shp.Polygon(cell))

    gdf = gpd.GeoDataFrame(geometry=polygons, crs=crs)
    vor.gdf_vorPolys = gdf
    vor.x_coords_by_node = xvertices_by_cells
    vor.y_coords_by_node = yvertices_by_cells
    return gdf


def get_gdf_latlon(vor):
    """
    Return the Voronoi polygons in latitude/longitude coordinates.
    """
    if vor._gdf_latlon is None:
        vor._gdf_latlon = vor.gdf_vorPolys.to_crs(vor.crs_latlon)
    return vor._gdf_latlon


def get_latlon(vor):
    """
    Return GeoJSON-like lat/lon geometry for the Voronoi grid.
    """
    if vor._latlon is None:
        vor._latlon = json.loads(vor.gdf_latlon["geometry"].to_json())
    return vor._latlon


def get_grid_centroid(vor):
    """
    Return the centroid of the Voronoi domain in lat/lon coordinates.
    """
    return vor.gdf_latlon.union_all().centroid


def get_overlapping_area(vor, shp_gpkg=None, cell_list=None):
    """
    Return the area of cells overlapping a geometry or explicit cell ids.
    """
    if shp_gpkg is not None:
        cells = vor.get_vor_cells_as_series(shp_gpkg).to_list()
    elif cell_list is not None:
        cells = cell_list
    else:
        raise ValueError('You must provide a shp_gpkg or cell_list')
    return vor.gdf_vorPolys.loc[cells].union_all().area


def get_vor_idx_from_geometry(
    vor,
    shp_to_query: shp = None,
    gdf_to_query: gpd.GeoDataFrame = None,
    crs: str = "EPSG:2927",
    crs_latlon: str = "EPSG:4326",
    name_col: str = None,
    predicate: str = "intersects",
) -> dict:
    """
    Determine intersecting Voronoi cells for each geometry in a query layer.
    """
    if shp_to_query:
        gdf_query = read_shp_gpkg(shp_to_query).to_crs(crs)
    elif gdf_to_query is not None:
        gdf_query = gdf_to_query.to_crs(crs)
    else:
        raise ValueError('Provide shp_to_query or gdf_to_query')

    vor_idx_dict = {}
    latlon_vor = vor.gdf_vorPolys["geometry"].to_crs(crs_latlon)
    latlon_query = gdf_query.to_crs(crs_latlon)

    for idx in latlon_query.index:
        cell_intersections = (
            gpd.GeoSeries(latlon_query.iloc[idx]["geometry"])
            .sindex.query(latlon_vor, predicate=predicate)[0]
            .tolist()
        )
        key = gdf_query.iloc[idx][name_col] if name_col is not None else idx
        vor_idx_dict[key] = cell_intersections

    return vor_idx_dict


def get_vor_idx_from_geometry_idx(
    vor,
    gdf_to_query: gpd.GeoDataFrame = None,
    idx: int = 0,
    predicate: str = "intersects",
) -> list[int]:
    """
    Determine intersecting Voronoi cells for a single geometry in a GeoDataFrame.
    """
    return (
        vor.gdf_vorPolys["geometry"]
        .sindex.query(gdf_to_query['geometry'].iloc[idx], predicate=predicate, sort=True)
        .tolist()
    )


def set_k_vor(vor, k_dict: dict = None, k_default=100) -> list:
    """
    Set hydraulic conductivity values on the Voronoi polygon GeoDataFrame.
    """
    vor.gdf_vorPolys["Kh"] = k_default
    for key in k_dict.keys():
        vor.gdf_vorPolys.loc[vor.gdf_vorPolys.index.isin(k_dict[key]), ["Kh"]] = key
    return vor.gdf_vorPolys["Kh"].to_list()


def generate_grid_around_point(center_point: shp.Point, spacing: float, size: int, crs: str) -> gpd.GeoSeries:
    """
    Generate a GeoSeries of grid points around a center point.
    """
    minx = center_point.x - size / 2 * spacing
    miny = center_point.y - size / 2 * spacing
    maxx = center_point.x + size / 2 * spacing
    maxy = center_point.y + size / 2 * spacing

    x_coords = list(range(int(minx), int(maxx) + 1, spacing))
    y_coords = list(range(int(miny), int(maxy) + 1, spacing))
    points = [shp.Point(x, y) for x in x_coords for y in y_coords]
    return gpd.GeoSeries(points).set_crs(crs)


def generate_grid_polygons(center_point, spacing, size, gap):
    """
    Generate square polygons in a grid pattern around a center point.
    """
    minx = center_point.x - size / 2 * spacing
    miny = center_point.y - size / 2 * spacing
    maxx = center_point.x + size / 2 * spacing
    maxy = center_point.y + size / 2 * spacing

    x_coords = list(range(int(minx), int(maxx) + 1, spacing))
    y_coords = list(range(int(miny), int(maxy) + 1, spacing))

    polygons = []
    for x in x_coords[:-1]:
        for y in y_coords[:-1]:
            polygons.append(shp.box(x + gap / 2, y + gap / 2, x + spacing - gap / 2, y + spacing - gap / 2))

    return shp.MultiPolygon(polygons)


def voronoi_refine_by_point(point: shp.Point, spacing: int, tri) -> gpd.GeoDataFrame:
    """
    Add a small refinement stencil around a point to a Triangle grid.
    """
    polypoints = [
        (point.x, point.y),
        (point.x + spacing, point.y),
        (point.x + spacing, point.y - spacing),
        (point.x, point.y - spacing),
    ]
    poly_main = shp.Polygon(polypoints)
    transform_dist = spacing * 2

    poly_e = shp.transform(poly_main, lambda x: x + [transform_dist, 0])
    poly_w = shp.transform(poly_main, lambda x: x - [transform_dist, 0])
    poly_n = shp.transform(poly_main, lambda x: x + [0, transform_dist])
    poly_s = shp.transform(poly_main, lambda x: x - [0, transform_dist])
    poly_ne = shp.transform(poly_main, lambda x: x + [transform_dist, transform_dist])
    poly_nw = shp.transform(poly_main, lambda x: x + [-transform_dist, transform_dist])
    poly_se = shp.transform(poly_main, lambda x: x + [transform_dist, -transform_dist])
    poly_sw = shp.transform(poly_main, lambda x: x + [-transform_dist, -transform_dist])
    all_polys = [poly_main, poly_e, poly_w, poly_n, poly_s, poly_ne, poly_nw, poly_se, poly_sw]

    gdf_all_polys = gpd.GeoDataFrame(geometry=all_polys)
    for poly_idx in range(len(gdf_all_polys)):
        tri.add_polygon(gdf_all_polys.loc[poly_idx, 'geometry'])

    return gdf_all_polys


def reconcile_surfaces(vor, df: pd.DataFrame = None, min_sep=0.1, trigger_sep=1, which='bottom'):
    """
    Adjust stacked surface elevations to preserve layer ordering.
    """
    df = vor.gdf_topbtm.copy() if df is None else df
    if isinstance(df, gpd.GeoDataFrame):
        df = df.drop(columns='geometry').map(lambda x: pd.to_numeric(x, errors='coerce'))
    else:
        df = df.map(lambda x: pd.to_numeric(x, errors='coerce'))

    labels = list(df.columns)
    df = df.loc[:, labels]
    for i, label in enumerate(labels):
        if i == 0:
            continue
        diffs = df.diff(axis=1)
        diff_list = list(diffs[diffs[label] >= -trigger_sep].index)
        if which == 'bottom':
            df.iloc[diff_list, i] = df.iloc[diff_list, (i - 1)] - min_sep
        elif which == 'top':
            df.iloc[diff_list, (i - 1)] = df.iloc[diff_list, i] + min_sep
        else:
            raise ValueError(f'which arg {which} is not valid. Must be "top" or "bottom"')

    return df


def adjust_cells_by_id(
    vor,
    cell_ids: list,
    adjustment: int | float,
    df: pd.DataFrame = None,
    layer: int = 0,
    reconcile: bool = True,
):
    """
    Adjust elevations for selected cells and optionally reconcile surfaces.
    """
    df = vor.gdf_topbtm.copy() if df is None else df.copy()
    df.loc[cell_ids, layer] = df.loc[cell_ids, layer] + adjustment
    return vor.reconcile_surfaces(df) if reconcile else df


def adjust_top_btm_overlaps(
    vor,
    elev_df=None,
    shp_path: Path = None,
    buffer=1,
    layer_bottom_name=1,
    min_sep=None,
) -> pd.DataFrame:
    """
    Lower selected bottoms so they overlap adjacent cell tops consistently.
    """
    elev_df = vor.reconcile_surfaces(min_sep=min_sep) if elev_df is None else elev_df
    cells_to_adjust = vor.get_vor_cells_as_series(shp_path)
    new_bottoms = {}

    for cell_id in cells_to_adjust:
        adjacent_cells = vor.find_adjacent_cells(cell_id)
        ja_cell_tops = elev_df[0].loc[adjacent_cells]
        ja_min = ja_cell_tops.min()
        new_bottoms[cell_id] = ja_min - buffer

    new_bottoms = pd.Series(new_bottoms, name=layer_bottom_name)
    elev_df[layer_bottom_name].update(new_bottoms)
    return vor.reconcile_surfaces(df=elev_df)

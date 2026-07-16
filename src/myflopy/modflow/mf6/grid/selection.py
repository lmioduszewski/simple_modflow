from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely as shp
from shapely.geometry import Polygon

from myflopy.modflow.utils.datatypes.readers import read_shp_gpkg


def get_vor_cells_as_series(
    gdf_vor_polys,
    overlapping_geometry: shp.Polygon | shp.Point | gpd.GeoSeries | Path = None,
    predicate: str = 'intersects',
    return_dict: bool = False,
    name_field: str = 'ExploName',
) -> pd.Series | dict:
    """
    Identify Voronoi cells matching a spatial predicate against provided geometries.
    """
    names = None

    if isinstance(overlapping_geometry, (shp.Polygon, shp.Point, shp.LineString, shp.MultiPolygon)):
        geometries = gpd.GeoSeries([overlapping_geometry])
        names = geometries.name
    elif isinstance(overlapping_geometry, (gpd.GeoDataFrame, gpd.GeoSeries)):
        geometries = overlapping_geometry.geometry
        if name_field in overlapping_geometry.columns:
            names = overlapping_geometry[name_field]
    elif isinstance(overlapping_geometry, Path):
        geoms = gpd.read_file(overlapping_geometry)
        geometries = geoms.geometry
        if name_field in geoms.columns:
            names = geoms[name_field]
    else:
        return ValueError('Unsupported data type for overlapping_geometry.')

    if names is None:
        names = list(range(len(geometries)))
    assert len(names) == len(geometries), 'length of names must match length of geometries'

    if geometries.crs is None:
        geometries.crs = gdf_vor_polys.crs
    elif geometries.crs != gdf_vor_polys.crs:
        print('Warning: Provided geometries are not in the same CRS as the Voronoi cells. Reprojecting geometries to match Voronoi cells.')
        geometries = geometries.to_crs(gdf_vor_polys.crs)

    joined = gpd.sjoin(gdf_vor_polys, geometries.to_frame('geometry'), predicate=predicate, how='inner')
    group_join = joined.groupby('index_right').apply(lambda frame: frame.index.to_list(), include_groups=False)
    group_join.name = 'cells'
    if not isinstance(names, list):
        names_idx = names.loc[group_join.index]
        group_join.index = names_idx

    if return_dict:
        return group_join.to_dict()
    return group_join


def get_vor_cells_as_dict(
    vor,
    locs: Path,
    crs: str = None,
    predicate: str = 'intersects',
    loc_name_field: str = None,
    return_gdf: bool = False,
) -> dict:
    """
    Provide a shapefile (locs) and get back the voronoi cells that contain them.
    """
    crs = vor.crs if crs is None else crs
    gdf_locs = gpd.read_file(locs).to_crs(crs)

    loc_vor_cell_dict = {}
    for idx in gdf_locs.index:
        location_name = idx if loc_name_field is None else gdf_locs.iloc[idx][loc_name_field]

        vor_cells = get_vor_cells_as_series(vor.gdf_vorPolys, gdf_locs.geometry[idx], predicate)
        if not isinstance(vor_cells, pd.Series):
            print(f'Warning: {location_name} does not intersect any Voronoi cells. Skipping it.')
            continue
        loc_vor_cell_dict[location_name] = vor_cells.tolist()

    if return_gdf:
        return loc_vor_cell_dict, gdf_locs
    return loc_vor_cell_dict


def get_model_boundary_polygons(gdf_vor_polys) -> dict:
    """
    Return a dict of polygons that form the model domain boundary.
    """
    print("Deprecated. Only works if the convex hull is equivalent to the actual grid boundary")
    print("Use get_grid_edge")
    grid = shp.MultiPolygon(gdf_vor_polys.geometry.to_list())
    polygons = gdf_vor_polys.geometry.to_list()
    convex_hull_boundary = grid.convex_hull.boundary

    boundary_polygons = []
    boundary_polygons_idx = []
    for idx, polygon in enumerate(polygons):
        if convex_hull_boundary.dwithin(polygon, 0.01):
            boundary_polygons.append(polygon)
            boundary_polygons_idx.append(idx)
    return dict(zip(boundary_polygons_idx, boundary_polygons, strict=False))


def get_grid_edge_cells(vor, idomain: list = None, idomain_path: Path = None, include_interiors: bool = True) -> list:
    """
    Get edge cells for the voronoi grid, optionally excluding inactive cells.
    """
    df = vor.gdf_vorPolys
    idomain_path = vor.idomain_path if idomain_path is None else idomain_path

    adj_lengths = np.fromiter((len(vor.adjacent_cells_idx[idx]) for idx in df.index), dtype=int)
    face_counts = df.geometry.exterior.apply(lambda geom: len(geom.coords) - 1).to_numpy()
    edge_mask = face_counts > adj_lengths
    edges = df.index[edge_mask].tolist()

    if idomain_path or idomain:
        if idomain_path:
            idomain_geom = read_shp_gpkg(idomain_path).union_all()
            icells = get_vor_cells_as_series(vor.gdf_vorPolys, idomain_geom)[0]
        else:
            icells = idomain

        icell_set = set(icells)
        remaining_cells = df.loc[~df.index.isin(icell_set)]
        union_geom = remaining_cells.geometry.union_all()

        if isinstance(union_geom, Polygon):
            edge_geoms = [union_geom.exterior]
            interiors = union_geom.interiors if include_interiors else []
        elif isinstance(union_geom, shp.MultiPolygon):
            edge_geoms = [poly.exterior for poly in union_geom.geoms]
            interiors = [ring for poly in union_geom.geoms for ring in poly.interiors] if include_interiors else []
        else:
            edge_geoms, interiors = [], []

        edge_cells = [pd.Series(get_vor_cells_as_series(vor.gdf_vorPolys, geom)[0]) for geom in edge_geoms]
        if include_interiors and interiors:
            edge_cells += [pd.Series(get_vor_cells_as_series(vor.gdf_vorPolys, ring)[0]) for ring in interiors]

        flat = pd.concat(edge_cells)
        flat.index = flat.values

        print('getting grid edge cells excluding idomain cells')
        return flat.loc[~flat.index.isin(icell_set)].index.tolist()

    print('geting all edge cells')
    return edges

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus as Vor

from pathlib import Path

import geopandas as gpd
import pandas as pd
import shapely as shp


def read_gpkg(gpkg_path: Path) -> gpd.GeoDataFrame:
    """
    Reads a GPKG file, iterates through all features, and returns a geopandas.GeoDataFrame with all found geometries
    Useful when you have multiple layers and geometries in a file.

    :param gpkg_path: gpkg file path
    :return: geopandas.GeoDataFrame with all found geometries
    """
    num_features = 0
    layer_num = 0
    layers = []
    gpkg = gpd.read_file(gpkg_path)
    crs = gpkg.crs
    try:
        for idx, row in gpkg.iterrows():
            if row.geometry is None:
                print(f'skipping row {idx} because the geometry is None')
                continue
            if isinstance(row.geometry, shp.Polygon | shp.Point | shp.MultiLineString | shp.LineString):
                layers.append(row)
                num_features += 1
            elif isinstance(row.geometry, shp.MultiPolygon):
                for geom in row.geometry.geoms:
                    new_row = row.copy()
                    new_row.geometry = geom
                    layers.append(new_row)
                    num_features += 1
            else:
                raise TypeError(f'Unexpected geometry type: {type(row.geometry)}')
            layer_num += 1
    # if layer number is not valid, end the while loop by setting layer to False
    except:
        layer = False
        if layer_num == 0:
            raise ValueError('Could not read gpkg file')

    print(f'Imported {num_features} features from {gpkg_path}')
    gdf = gpd.GeoDataFrame.from_records(data=layers)
    gdf.set_geometry('geometry', inplace=True)
    gdf.crs = crs

    return gdf


def read_shp_gpkg(files: list | Path) -> gpd.GeoDataFrame:
    """
    return a GeoDataFrame with all found geometries from a provided shapefile, geopackage, or list thereof.
    :param files: shapefile, geopackage, or list thereof
    :return: GeoDataFrame with all found geometries
    """
    if isinstance(files, Path):
        files = [files]
    if isinstance(files, list):
        assert all(isinstance(path, Path) for path in files), 'all items in list must be Path objects'
    else:
        raise TypeError('files must be Path objects')

    gdfs = []
    for file in files:
        assert isinstance(file, Path)
        if file.suffix == '.shp':
            d = gpd.read_file(file)
        elif file.suffix == '.gpkg':
            d = read_gpkg(file)
        else:
            raise ValueError(f'File type {file.suffix} not supported')
        gdfs.append(d)

    gdf_join = pd.concat(gdfs)
    gdf_join = gdf_join.reset_index(drop=True)
    gdf_join = gpd.GeoDataFrame(gdf_join, geometry='geometry')

    return gdf_join


def _parse_string_for_ints(s: str) -> list | None:
    """
    Parses a given string to extract integers. The function attempts to convert
    each space-separated value in the provided string into an integer. If any value
    cannot be converted to an integer or if the input string results in an invalid
    conversion, the function will return None and will give an error message.

    :param s: The input string to parse and convert into integers.
    :type s: str
    :return: A list of integers if conversion is successful, otherwise None.
    :rtype: list | None
    """
    try:
        ints = [int(x) for x in s.strip().split() if x.strip() not in ["", ","]]
        assert all(isinstance(i, int) for i in ints), \
            ("couldn't convert all values to integers \n "
             f"got: {ints}")
        return ints

    except ValueError(f'cannot convert all values to integers. \n'
                      f'double check: {s}'):
        return None


def parse_geodataframe_for_layer_ints(
        gdf: gpd.GeoDataFrame,
        col_name: str = 'layer'
) -> gpd.GeoDataFrame:
    """
    Processes a GeoDataFrame to parse integer layers from a given column and creates new columns
    for each unique integer layer. The function first parses the specified column for integers and
    identifies all unique integers across the entire column. For each unique integer, a new column
    is added to the GeoDataFrame where its value represents the integer if it exists in the parsed
    data for that row; otherwise, the value will be None.

    :param gdf: A GeoDataFrame that contains the column to parse.
    :param col_name: The name of the column in the GeoDataFrame to parse for integer values. Defaults to 'layer'.
    :return: A copy of the GeoDataFrame with additional columns, where each column corresponds
        to a unique integer parsed from the target column. The values in these columns are the
        integer if it exists for the given row, or None if it does not.
    """
    parsed = gdf[col_name].apply(_parse_string_for_ints)
    unique_layers = sorted(set(i for sublist in parsed.dropna() for i in sublist))

    for layer in unique_layers:
        gdf[str(layer)] = parsed.apply(lambda x: layer if x and layer in x else None)

    return gdf


def get_layer_col_names(gdf: gpd.GeoDataFrame) -> list:
    """
    Extract column names from a GeoDataFrame representing layer names as integers.

    This function scans through all column names of the given GeoDataFrame and
    selects the ones that are entirely numeric. These are interpreted as
    representing layer names. The function also validates that all these numeric
    layer names can be successfully converted to integers.

    :param gdf: Input GeoDataFrame from which to extract numeric layer column names
    :type gdf: gpd.GeoDataFrame
    :return: List of numeric column names from the GeoDataFrame, converted as strings
    :rtype: list
    """
    names = [col for col in gdf.columns if col.isdigit()]
    assert all(isinstance(int(col), int) for col in names), \
        f'all layer names must be integers. got: {names}'
    return names


def assign_voronoi_cells_to_layers(
        gdf: gpd.GeoDataFrame,
        vor: Vor,
        layer_col: str = 'layer',
        return_col_names: bool = False,
) -> gpd.GeoDataFrame | tuple[gpd.GeoDataFrame, list]:
    """
    Uses parse_geodataframe_for_layer_ints() to create layer columns,
    and replaces the values in those columns with Voronoi cell indices
    from the geometry in each row if the layer applies to that row.

    :param return_col_names: if True will also return a list of column names as second return value
    :param gdf: Input GeoDataFrame with geometries and a string column of layer IDs.
    :param vor: voronoi grid object.
    :param layer_col: Name of the column containing the space-separated layer strings.
    :return: Modified GeoDataFrame with columns for each unique layer, each filled
             with Voronoi cell indices or None for the geometry in each row.
    """
    # Parse the layers and create one column per layer
    gdf = parse_geodataframe_for_layer_ints(gdf, col_name=layer_col)

    # Identify which columns are layer columns (integers as strings)
    layer_columns = get_layer_col_names(gdf)

    # Replace each value with voronoi cell indices where appropriate
    for layer in layer_columns:
        gdf[layer] = gdf.apply(
            lambda row: list(vor.get_vor_cells_as_series(row.geometry)) if pd.notna(row[layer]) else None,
            axis=1
        )
    if return_col_names:
        return gdf, layer_columns
    else:
        return gdf


def explode_dxf_polygons(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Take a GeoDataFrame of CAD-imported geometries (often GeometryCollections /
    invalid multipolygons) and return a cleaned GeoDataFrame of single-part polygons.
    """

    # 0) Drop Z/M values (DXF often has 3D coords)
    gdf = gdf.copy()
    gdf.geometry = gdf.geometry.map(shp.force_2d)

    # 1) Fix invalid geometries (splits weird rings into valid Polygons/MultiPolygons)
    gdf.geometry = gdf.geometry.map(shp.make_valid)

    # 2) Extract polygonal parts from GeometryCollections or MultiPolygons
    parts, rows = [], []
    for row in gdf.itertuples(index=False):
        geom = row.geometry
        for part in shp.get_parts(geom):
            if part.geom_type in ("Polygon", "MultiPolygon"):
                parts.append(part)
                rows.append(row)

    poly_gdf = gpd.GeoDataFrame(rows, geometry=parts, crs=gdf.crs)

    # 3) Explode multipolygons into single polygons
    poly_gdf = poly_gdf.explode(index_parts=False, ignore_index=True)

    return poly_gdf

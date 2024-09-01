import pandas as pd
from pathlib import Path
import shapely as shp
import geopandas as gpd


def read_gpkg(filepath: Path, crs="EPSG:2927"):
    """
    Reads a GPKG file, iterates through all features, and returns a geopandas.GeoDataFrame with all found geometries

    :param filepath: gpkg file path
    :param crs: coordinate reference system, defaults to EPSG:2927
    :return: geopandas.GeoDataFrame with all found geometries
    """
    layer = True
    layer_num = 0
    layers = []
    num_features = 0

    while layer:
        try:
            f = gpd.read_file(filepath, layer=layer_num)
            for g in f.geometry:
                if isinstance(g, shp.Polygon):
                    layers.append(g)
                    num_features += 1
                elif isinstance(g, shp.MultiPolygon):
                    for geom in g.geoms:
                        layers.append(geom)
                        num_features += 1
                else:
                    raise TypeError(f'Unexpected geometry type: {type(g)}')
            layer_num += 1
        except:
            layer = False
            if layer_num == 0:
                raise ValueError('Could not read gpkg file')

    print(f'Imported {num_features} features from {filepath}')
    gdf = gpd.GeoDataFrame(geometry=layers, crs=crs)

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

    return gdf_join

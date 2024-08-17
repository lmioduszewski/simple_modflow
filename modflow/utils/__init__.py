from simple_modflow.modflow.utils.surfaces import InterpolatedSurface
import geopandas as gpd
from pathlib import Path

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
            geometry = gpd.read_file(filepath, layer=layer_num)
            for feature in geometry.geometry:
                layers.append(feature)
                num_features += 1
            layer_num += 1
        except:
            layer = False
            if layer_num == 0:
                raise ValueError('Could not read gpkg file')
    print(f'Imported {num_features} features from {filepath}')
    gdf = gpd.GeoDataFrame(geometry=layers, crs=crs)
    return gdf
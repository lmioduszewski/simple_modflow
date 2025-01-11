from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
    from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor

import rasterio
import shapely as shp
from pathlib import Path
import geopandas as gpd
import pickle
import pandas as pd


class RasterData:

    def __init__(self, raster_path: Path, point: shp.Point = None, vor: Vor = None):

        self.raster_path = raster_path
        self.point = point
        self.vor = vor

    @property
    def raster_elevs(self):
        with rasterio.open(self.raster_path) as src:
            elevations = src.read(1)
        return elevations

    def sample_raster(self):
        """
        Returns the raster values for the given iterable of Shapely points.
        """

        if self.point:
            print(f'sampling raster from single point {self.point}')
            points = [self.point]
        elif self.vor:
            print(f'sampling raster at all centroids in voronoi grid')
            points = self.vor.gdf_vorPolys.centroid.to_list()
        else:
            raise ValueError('No valid sample points provided')

        with rasterio.open(self.raster_path) as src:
            # Convert Shapely points to (x, y) coordinate pairs
            coords = [(point.x, point.y) for point in points]

            # rasterio.sample returns an iterator of arrays (one per band)
            # For a single-band raster, each returned array has 1 value: e.g. array([value])
            # Extract that single value with val[0].
            values = [val[0] for val in src.sample(coords)]

        return pd.Series(values)


if __name__ == '__main__':

    prism_raster = Path(r"C:\Users\lukem\mf6\Cumberland general\PRISM_ppt_30yr_normal_800mM4_annual_bil\PRISM_ppt_30yr_normal_annual_inches_EPSG_2926.tif")
    weather_stn = Path(r"C:\Users\lukem\mf6\Cumberland general\landsburg.gpkg")
    vor_path_v2c = Path(r"C:\Users\lukem\mf6\Cumberland general\cumberland_v2c.vor")

    with open(vor_path_v2c, 'rb') as file:
        vor: Vor = pickle.load(file)
    point = gpd.read_file(weather_stn).geometry.iloc[0]

    rs = RasterData(prism_raster, vor=vor)
    print(rs.sample_raster())


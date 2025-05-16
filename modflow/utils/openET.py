from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
    from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor

import requests
import glob
import rasterio as rio
import pandas as pd
import geopandas as gpd
from shapely import Polygon
from pathlib import Path


class OpenETRequest:
    """Class to request ET data from OpenET"""

    def __init__(self, date_range=None, api_key=None, geometry=None, interval="monthly", units="in"):
        self.stack_url = "https://openet-api.org/raster/export/stack"
        self.download_url = "https://openet-api.org/account/storage"
        self.date_range = date_range
        self.api_key = api_key
        self.geometry = geometry
        self.interval = interval
        self.units = units

        self.header = {"Authorization": api_key}
        self.date_range = date_range

    def request_args(self):
        # endpoint arguments
        args = {
            "date_range": self.date_range,
            "interval": self.interval,
            "geometry": self.geometry,
            "model": "Ensemble",
            "variable": "ET",
            "reference_et": "gridMET",
            "units": self.units,
            "encrypt": False
        }
        return args

    def request(self):
        """Requests a stack of OpenET raster data based on class inputs"""
        req = requests.post(
            headers=self.header,
            json=self.request_args(),
            url=self.stack_url
        )
        return req

    def dl(self):
        """Requests download links from OpenET API based on the request. Run self.request first"""
        dl_req = requests.get(
            headers=self.header,
            url=self.download_url
        )
        return dl_req

    def dl_links(self):
        """get download links"""
        dl = self.dl()
        return dl.json()


class OpenETtoVor:
    """Class to assign OpenET data to Voronoi Grid cells"""

    def __init__(self, vor: Vor, raster_dir: Path):
        self.vor = vor
        self.raster_dir = raster_dir
        self._vor_et = None

    @property
    def raster_paths(self):
        gtif = glob.glob(pathname='*.tif', root_dir=self.raster_dir)
        gtiff = glob.glob(pathname='*.tiff', root_dir=self.raster_dir)
        g = gtif + gtiff
        tif_path_list = [self.raster_dir / f for f in g]
        return tif_path_list

    @property
    def vor_et(self):
        """returns a dictionary where the keys are the OpenET raster paths, and values are
        a Pandas Series where the index is the voronoi cell index and the data are area-weighted
        ET values for each voronoi cell"""
        if self._vor_et is None:

            et_per_vorcell = {}

            for i, path in enumerate(self.raster_paths):

                print(f'working on {path}, number {i + 1} of {len(self.raster_paths)}', end='\r')

                with rio.open(path) as src:

                    src: rio.DatasetReader
                    crs = src.crs
                    band = src.read(1)
                    band = band.astype('float32')

                    # get polygon equivalents of each raster pixel
                    poged = [(p, v) for p, v in rio.features.shapes(band, transform=src.transform)]
                    poged = [[Polygon(poged[i][0]['coordinates'][0]), poged[i][1]] for i in range(len(poged))]

                    # make DataFrames of polygonized raster and voronoi polygons
                    df = pd.DataFrame(poged)
                    df.columns = ['geometry', 'ET']
                    gdf = gpd.GeoDataFrame(df.iloc[:, 1], geometry=df.iloc[:, 0], crs=crs)
                    gdf = gdf.to_crs(crs=self.vor.crs)
                    vor_polys = self.vor.gdf_vorPolys.copy()
                    vor_polys['idx'] = vor_polys.index

                    # get intersections of polygonized raster and voronoi polygons
                    intersections = gpd.overlay(gdf, vor_polys, how='intersection')
                    intersections["int_area"] = intersections.area
                    # Weighted "mass" = raster_val * intersection area
                    intersections["weighted_val"] = intersections['ET'] * intersections["int_area"]

                    # Group by the Voronoi polygon IDX
                    grouped = intersections.groupby("idx").agg(
                        total_area=("int_area", "sum"),
                        sum_weighted_val=("weighted_val", "sum")
                    )
                    grouped["mean_ET"] = grouped["sum_weighted_val"] / grouped["total_area"]
                    et = grouped.reset_index().set_index('idx')['mean_ET']
                    et_per_vorcell[path] = et

            self._vor_et = et_per_vorcell

        return self._vor_et


if __name__ == "__main__":

    import pickle
    raster_dir = Path(r"C:\Users\lukem\mf6\Cumberland general\ET")
    model_path = Path(r"C:\Users\lukem\mf6\cumb_v7b\cumb_v7b.model")
    with open(model_path, 'rb') as file:
        model = pickle.load(file)
    vor = model.vor
    et = OpenETtoVor(vor, raster_dir)

    et_path = Path().home() / 'mf6' / 'et_per_cell_pits.et'
    with open(et_path, 'wb') as file:
        pickle.dump(et.vor_et, file)

    """api_key = '8n9LsycWsdg6EQ2RGgD8O3mBKQBRQNN1GBGXMK2JaMpUbMwBfCfohscxqevS'
    date_range = ["2024-01-01", "2024-02-01"]
    geometry = [-121.973571242, 47.326251526, -121.897909136, 47.327471047,
                -121.897390248, 47.259707461, -121.970796765, 47.257356792]

    r = OpenETRequest(date_range=date_range, api_key=api_key, geometry=geometry)
    # r.request()
    print(r.dl_links())"""
    pass

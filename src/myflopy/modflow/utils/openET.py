from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus as Vor

import glob
import re
from pathlib import Path

import geopandas as gpd
import pandas as pd
import rasterio as rio
import rasterio.features as features
import requests
from shapely import Polygon


def extract_dates_from_paths(paths) -> list[Path]:
    """
    Extracts dates in 'YYYY_MM_DD' format from a list of file paths.

    This function iterates over a list of file paths, searches each path for
    a date string in the format 'YYYY_MM_DD', and appends any matches to a
    list. Dates are identified using a regular expression pattern.

    :param paths: A list of file paths where the function will search for
                  date strings.
    :type paths: list[Path]

    :return: A list of date strings in the 'YYYY_MM_DD' format extracted
             from the input paths.
    :rtype: list[str]
    """
    pattern = re.compile(r"\d{4}_\d{2}_\d{2}")
    dates = []

    for path in paths:
        path_str = str(path)
        match = pattern.search(path_str)
        if match:
            dates.append(match.group(0))
    return dates


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
                    poged = [(p, v) for p, v in features.shapes(band, transform=src.transform)]
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

    def avg_months(self, months: list = None):
        """
        Calculate monthly averages from a set of evapotranspiration data.

        This function processes evapotranspiration data retrieved from the
        OpenETtoVor instance and groups them by month. It averages the data
        for each specified month and fills any missing or zero values in
        columns with the column-wise mean. Finally, it outputs the monthly
        averages in a dictionary, ordered from September to October.

        :param months: List of integers representing the months for each raster file
            for which to calculate the monthly averages. Each integer should correspond
            to the month for a specific file, in order of the raster files in the directory.
        :return: Dictionary with keys as months and values as Pandas Series
            containing the average values for each grid cell.
        :rtype: dict
        """
        vor_et = OpenETtoVor(self.vor, self.raster_dir).vor_et
        # get month integers from file names
        if months is None:
            dates = extract_dates_from_paths(self.raster_paths)
            months = [int(m[5:7]) for m in dates]
        keys = list(vor_et.keys())
        et_dict = {}
        for month in pd.Series(months).unique().tolist():
            et_dict[month] = []
        for i, month in enumerate(months):
            ets = vor_et[keys[i]]
            ets[ets < 0] = 0
            et_dict[month].append(ets)
        et_concat = {}
        for month in months:
            et_concat[month]: pd.DataFrame = pd.concat(et_dict[month], axis=1)

        mon_avgs = {}
        for month in pd.Series(months).unique().tolist():

            mo = et_concat[month].reindex(list(range(self.vor.ncpl)), fill_value=0)
            mo.columns = list(range(len(mo.columns)))
            for column in mo.columns:
                m_col = mo.loc[:, column]
                col_mean = m_col[m_col != 0].mean()
                m_col[m_col == 0] = col_mean
            mon_avg = mo.mean(axis=1)
            mon_avgs[month] = mon_avg
        key_order = [10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        mon_avgs = {key: mon_avgs[key] for key in key_order if key in mon_avgs}

        return mon_avgs


if __name__ == "__main__":
    import pickle

    raster_dir = Path(r"C:\Users\lukem\mf6\Cumberland general\ET\alt")
    vor_v12 = Path(r"C:\Users\lukem\mf6\Cumberland general\cumb_v14c_existing_algomesh_v2.1.vor")
    with open(vor_v12, 'rb') as file:
        vor = pickle.load(file)
    et = OpenETtoVor(vor, raster_dir).avg_months()

    et_path = Path().home() / 'mf6' / 'et_new_v13.et'
    with open(et_path, 'wb') as file:
        pickle.dump(et, file)

    """api_key = '8n9LsycWsdg6EQ2RGgD8O3mBKQBRQNN1GBGXMK2JaMpUbMwBfCfohscxqevS'
    date_range = ["2024-01-01", "2024-12-01"]
    geometry = [-121.96498100745063, 47.263702752801436, -121.89545157703465, 47.263702752801436,
                -121.89545157703465, 47.32354206806991, -121.96498100745063, 47.32354206806991, 
                -121.96498100745063, 47.263702752801436]

    r = OpenETRequest(date_range=date_range, api_key=api_key, geometry=geometry)
    # r.request()
    print(r.dl_links())"""


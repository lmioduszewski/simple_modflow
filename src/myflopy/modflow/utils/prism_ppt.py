from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus as Vor

from pathlib import Path

from shapely import Point

from myflopy import read_shp_gpkg
from myflopy.modflow.utils.raster import RasterData


class PrismPrecipScaling:
    """class to determine the precipitation scaling for a voronoi grid based on a
    weather station and the PRISM 30-year normal"""

    def __init__(
            self,
            vor: Vor,
            prism_path: Path,
            weather_station: Path,
    ):
        """
        Use to determine the precipitation scaling for a voronoi grid based on PRISM 30-year normal and
        a weather station location
        :param vor: Voronoi grid
        :param prism_path: path to PRISM 30-year normal tif
        :param weather_station: path to weather station point shapefile or geopackage
        """
        self.vor = vor
        self.prism_path = prism_path
        self._weather_station = None

        self.weather_station = weather_station

    @property
    def weather_station(self):
        """The weather-station location point used as the PRISM scaling reference."""

        return self._weather_station

    @weather_station.setter
    def weather_station(self, weather_station):
        """Set the station from a shapefile/GeoPackage ``Path`` (first point; CRS must match the grid)."""

        assert isinstance(weather_station, Path), "weather_station must be a Path object"
        weather_station = read_shp_gpkg(weather_station)
        if len(weather_station) > 1:
            print(f'More than one point provided for weather station. Assuming 1st weather station: '
                  f'{weather_station.geometry.iloc[0]} is what you want')
        assert weather_station.crs == self.vor.crs, \
            f'weather_station crs {weather_station.crs} does not match voronoi grid {self.vor.crs}'
        station = weather_station.geometry.iloc[0]
        assert isinstance(station, Point), "station must be a Point object"
        self._weather_station = station

    @property
    def scaling(self):
        """returns a Pandas Series where the values are the scaling factors to use for precipitation for each
        voronoi cell, listed by voronoi cell number index"""

        weather_station_ppt = RasterData(self.prism_path, point=self.weather_station).sample_raster()
        voronoi_ppt = RasterData(self.prism_path, vor=self.vor).sample_raster()
        scaling = voronoi_ppt / weather_station_ppt.iloc[0]
        return scaling


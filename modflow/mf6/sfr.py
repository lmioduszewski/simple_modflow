import flopy
import pandas as pd
import geopandas as gpd
from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor
from pathlib import Path
import shapely as shp

class SFR:

    def __init__(
            self,
            vor: Vor = None,
            stream_path: Path = None
    ):
        self.vor = vor
        self.stream_path = stream_path
        self.stream_gdf = gpd.read_file(stream_path)
        self.stream_geom = self.stream_gdf.geometry
        self.stream_points = {}
        self.vor_intersect = {}

        for i, geom in enumerate(self.stream_geom):
            x, y = geom.coords.xy
            self.stream_points[i] = [shp.Point(point) for point in list(zip(x,y))]
            self.vor_intersect[i] = vor.get_vor_cells_as_series(geom)
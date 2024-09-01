from simple_modflow.modflow.utils.surfaces import InterpolatedSurface
from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor
import figs as f
import pandas as pd
import geopandas as gpd


class SurfaceData(InterpolatedSurface):

    def __init__(self):
        super().__init__()



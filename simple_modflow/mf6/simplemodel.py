from pathlib import Path
import flopy
from voronoiplus import VoronoiGridPlus as Vor
from voronoiplus import TriangleGrid as Triangle
from boundaries import Boundaries
import numpy as np
from shapely import Polygon
from headsplus import HeadsPlus as hp
import shapely as shp
import mfsimbase as mf


class SimpleModel(mf.SimulationBase):
    """Class for a simple modflow model for quick results and analysis."""

    def __init__(
            self,
            vor: Vor,
            nper: int = 1,
            initial_sat_thickness: float = 1,
            k: int = 100,
            top=None,
            bottom=None,
            rch_dict: dict = None,
            boundary_conductance: int = 1000
    ):
        super().__init__()
        self.vor = vor
        self.k = k
        self.top = top
        self.bottom = bottom
        self.rch_dict = rch_dict
        self.initial_sat_thickness = initial_sat_thickness
        self.boundary_cells = list(self.vor.get_model_boundary_polygons().keys())
        self.tdis = mf.TemporalDiscretization(model=self, nper=nper)
        self.disu = mf.DisuGrid(vor=self.vor, model=self, top=self.top, bottom=self.bottom)
        self.drain_stress_period_data = Boundaries(vor=self.vor).get_drn_stress_period_data(
            cells=self.boundary_cells,
            bottom_addition=0,
            conductance=boundary_conductance,
        )
        self.ic = mf.InitialConditions(model=self, vor=self.vor, initial_sat_thickness=self.initial_sat_thickness)
        self.k = mf.KFlow(model=self, k_list=self.k)
        self.oc = mf.OutputControl(model=self)
        self.drn = mf.Drains(model=self, stress_period_data=self.drain_stress_period_data)
        if rch_dict:
            self.rch = mf.Recharge(model=self, vor=self.vor, rch_dict=self.rch_dict)


if __name__ == "__main__":
    """Example simple model"""
    tri = Triangle(
        model_ws=Path.cwd(),
        angle=30
    )
    tri.add_circle(radius=100, center_coords=(50, 200), )
    tri.add_circle(radius=5, center_coords=(50, 200), )
    tri.add_circle(radius=5, center_coords=(50, 150), )
    tri.add_regions(((50, 150), (0, 200)),
                    maximum_areas=(0.5, 10))
    tri.model_ws = Path.cwd().joinpath('sample_model_output')
    tri.build()
    # get bottom elevs
    strike = 30  # Strike given in degrees from north
    dip = 20  # Dip given in degrees from horizontal
    known_point = (50, 200, 0)  # Known point (x, y, elevation)
    pixel_size = 1  # Pixel size
    bottom_raster_path = Path.cwd().joinpath('sample_model_output', 'bottom_raster.tif')
    top_raster_path = Path.cwd().joinpath('sample_model_output', 'top_raster.tif')
    vor = Vor(tri)
    bottom_elevs = vor.get_raster_from_strike_dip(strike, dip, known_point, pixel_size, bottom_raster_path)
    top_elevs = vor.get_raster_from_strike_dip(0, 0, (0,0,50), 1, top_raster_path)
    vor = Vor(tri, rasters=[bottom_raster_path, top_raster_path])
    nper = 31
    center_cells = [32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56,
                    57, 58, 59, 60, 61, 62, 63, 98, 123, 124, 163, 164, 165, 193, 194, 195, 196, 210, 211, 212, 213,
                    234, 235, 236, 314, 337, 374, 375, 388, 392]
    rch_trans = [np.random.random() + 4 for per in range(nper)]
    rch_dict = {}
    for per in range(nper):
        cell_list = []
        for cell in range(vor.ncpl):
            if cell in center_cells:
                cell_list.append([cell, rch_trans[per]])
            else:
                cell_list.append([cell, 0.01])

        rch_dict[per] = cell_list
    model = SimpleModel(
        vor,
        #bottom=bottom_elevs['elev'].to_list(),
        #top=50,
        nper=nper,
        rch_dict=rch_dict
    )
    model.run_simulation()
    hds = hp(hds_path=model.model_output_folder_path.joinpath('mf6_model.hds'), vor=vor)
    hds.plot_choropleth((19, 0), zoom=19, plot_mounding=True).show()
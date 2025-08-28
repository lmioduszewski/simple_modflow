from pathlib import Path
import flopy
from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor
from simple_modflow.modflow.mf6.voronoiplus import TriangleGrid as Triangle
from simple_modflow.modflow.mf6.boundaries import Boundaries
import numpy as np
from shapely import Polygon
from simple_modflow.modflow.mf6.headsplus import HeadsPlus as hp
import shapely as shp
from simple_modflow.modflow.mf6 import mfsimbase as mf
import pickle
from simple_modflow.modflow.mf6.drn import DRN


class SimpleModel(mf.SimulationBase):
    """Class for a simple modflow model for quick results and analysis."""

    def __init__(
            self,
            vor: Vor,
            name: str = 'simplemodel',
            nper: int = 1,
            initial_sat_thickness: float = 1,
            iheads=None,
            k: int | float | list = 100,
            k33_vert: int | float | list = None,
            top=None,
            bottom=None,
            rch_dict: dict = None,
            boundary_conductance: int = 1000,
            nlay=1,
            disv_disu='disv'
    ):
        super().__init__(nper=nper, name=name, vor=vor)
        self.vor = vor
        self.nper = nper
        self.name = name
        self.k = k
        self.iheads = iheads
        self.nlay = nlay
        self.top = top
        self.bottom = bottom
        self.rch_dict = rch_dict
        self.initial_sat_thickness = initial_sat_thickness
        self.boundary_cells = self.vor.get_grid_edge()
        self.tdis = mf.TemporalDiscretization(model=self, per_len=30)
        if disv_disu == 'disu':
            self.disu = mf.DisuGrid(vor=self.vor, model=self, top=self.top, bottom=self.bottom)
        if disv_disu == 'disv':
            self.disv = mf.DisvGrid(vor=self.vor, model=self, nlay=nlay, top=self.top, bottom=self.bottom)
        self.drain_stress_period_data = []
        for layer in range(self.nlay):
            bottom_addition = 20 if layer == 2 else 0.1
            self.drain_stress_period_data += DRN(model=self, vor=self.vor).get_drn_stress_period_data(
                cells=self.boundary_cells,
                bottom_addition=bottom_addition,
                conductance=boundary_conductance,
                disMf=disv_disu,
                layer=layer
            )
        self.ic = mf.InitialConditions(model=self, vor=self.vor, initial_sat_thickness=self.initial_sat_thickness,
                                       nlay=self.nlay, strt=self.iheads)
        self.k = mf.KFlow(model=self, k=self.k, k33_vert=k33_vert)
        self.oc = mf.OutputControl(model=self)
        self.drn = mf.Drains(model=self, stress_period_data=self.drain_stress_period_data)
        self.sto = mf.Storage(model=self, specific_yield=0.2, specific_storage=0.0001,
                              sto_transient={1: True}, sto_steady={0: True})
        if rch_dict:
            self.rch = mf.Recharge(model=self, vor=self.vor, rch_dict=self.rch_dict)


if __name__ == "__main__":
    """Example simple model"""
    tri = Triangle(
        model_ws=Path.cwd(),
        angle=30
    )
    tri.add_circle(radius=10_000, center_coords=(0, 0))
    # tri.add_circle(radius=1000, center_coords=(0, 0))
    """ssb = tri.add_circle(2, (0, 0))
    ssb2 = tri.add_circle(2, (50, 0))
    ssb3 = tri.add_circle(2, (0, 50))
    ssb4 = tri.add_circle(2, (50, 50))"""
    outer_circle = tri.add_circle(50, (0, 0))

    # xs = np.arange(-150, 250, 100)
    """xs = np.arange(-50, 50, 50)
    xx, yy = np.meshgrid(xs, xs)
    coords = list(zip(xx.ravel(), yy.ravel()))
    uics = []
    for coord in coords:
        uics.append(tri.add_circle(10, coord))"""


    #ssb2 = tri.add_rectangle(10, 300, origin=(150, -150), max_area=50)
    #ssb3 = tri.add_rectangle(10, 300, origin=(0, -150), max_area=50)
    tri.add_region(point=(-400, -400), maximum_area=10_000)
    #tri.add_region(point=(-2000, -2000), maximum_area=8000)
    tri.add_region(point=(-10, 0), maximum_area=1000)
    tri.model_ws = Path.cwd().joinpath('sample_model_output')
    tri.build()

    # get bottom elevs
    pixel_size = 1000  # Pixel size
    l1_botm = Path.cwd().joinpath('sample_model_output', 'l1_botm.tif')
    l2_botm = Path.cwd().joinpath('sample_model_output', 'l2_botm.tif')
    l3_botm = Path.cwd().joinpath('sample_model_output', 'l3_botm.tif')
    top_raster_path = Path.cwd().joinpath('sample_model_output', 'top_raster.tif')
    vor = Vor(tri)
    vor.get_disu_connectivity()
    top_elevs = vor.get_raster_from_strike_dip(0, 0, (0, 0, 660), pixel_size, top_raster_path)
    l1_bottom_elevs = vor.get_raster_from_strike_dip(0, 0, (0, 0, 640), pixel_size, l1_botm)
    l2_bottom_elevs = vor.get_raster_from_strike_dip(0, 0, (0, 0, 550), pixel_size, l2_botm)
    l3_bottom_elevs = vor.get_raster_from_strike_dip(0, 0, (0, 0, 420), pixel_size, l3_botm)

    #with open(Path(r"C:\Users\lukem\mf6\simplemodel\iheads.hds"), 'rb') as file:
        #iheads = pickle.load(file)  # initial heads import

    vor = Vor(tri, rasters=[top_raster_path, l1_botm, l2_botm, l3_botm])
    surface_elevs = vor.reconcile_surfaces()
    vor.gdf_topbtm.loc[:, 0:] = surface_elevs
    nper = 1

    rch_dict = {per: [] for per in range(nper)}
    c = vor.get_vor_cells_as_series(outer_circle).to_list()
    ssb_cells = c

    facility_area = vor.get_overlapping_area(cell_list=ssb_cells)
    rch_trans = [(2_000 * 1) / facility_area for i in range(nper)]  # 100 gpm for one UIC

    for per in range(nper):
        cell_list = []
        for cell in range(vor.ncpl):
            if cell in ssb_cells:
                cell_list.append([(2, cell), rch_trans[per]])
            else:
                # cell_list.append([(0, cell), 0.001])
                pass

        rch_dict[per] += cell_list

    botm1 = surface_elevs.loc[:, 1].values
    botm2 = surface_elevs.loc[:, 2].values
    botm3 = surface_elevs.loc[:, 3].values
    botms = [botm1, botm2, botm3]

    model = SimpleModel(
        name='ssb_flat',
        vor=vor,
        k=[1000, 40, 40],
        bottom=botms,
        top=vor.gdf_topbtm.loc[:, 0].to_list(),
        nper=nper,
        rch_dict=rch_dict,
        nlay=3,
        boundary_conductance=5000,
    )

    model_file_path = model.model_output_folder_path / f'{model.name}.model'
    with open(model_file_path, 'wb') as file:
        pickle.dump(model, file)
    # print(rch_dict)
    model.run_simulation()
    vor.show()
    # model.choro(kstpkper=(9, 0)).plot()
    # model.xsect(cells=[2247, 2245], spacing=1).show()
    # model.surf.hds(per=0).plot()

import flopy.utils.binaryfile
import shapely as shp
from pathlib import Path
import geopandas as gpd
import pandas as pd
from simple_modflow.modflow.mf6.headsplus import HeadsPlus as hp
import simple_modflow.modflow.mf6.mfsimbase as mf
import grid
import simple_modflow.modflow.mf6.mfsimbase as mf
import grid
from simple_modflow.modflow.mf6.headsplus import HeadsPlus as hp
from simple_modflow.modflow.mf6.boundaries import Boundaries
from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor
import pickle
from simple_modflow.modflow.mf6.boundaries import Boundaries

idxx = pd.IndexSlice
# Conversion factors
inches_to_feet = 1 / 12

class RechargeFromShp(Boundaries):

    def __init__(
            self,
            model: mf.SimulationBase = None,
            vor: Vor = None,
            shp: Path = None,
            uid: str = None,
            crs: int = 2927,
            rch_cols = None

    ):
        """
        class to set up recharge for a modflow 6 model
        :param model: model to which this boundary applies
        :param vor: voronoi grid to which this boundary apples
        :param shp: path to shapefile that holds the polygons for the boundary
        :param uid: the field name in the shapefile attribute table that holds the unique ids, one for each polygon
        :param crs: coordinate reference system for boundary, should be integer EPSG code.
        :param rch_cols: field names corresponding to the recharge data in the shapefile attribute table
        :param bound_type: arbitary identifier for this boundary type
        """
        super().__init__(model, vor, shp, uid, crs)
        self.bound_type = 'rch'



    def get_rch(
            self,
            cell_ids: dict = None,
            recharges: dict = None,
            grid_type:str = 'disv',
            background_rch: int | float = None,
            nper: int = 1
    ) -> dict:
        """
        get a recharge dictionary to pass to flopy in setting of a recharge package. Assumes recharge only applied to
        top layer
        :param cell_ids: dictionary where each key is an arbitrary name given each recharge area and the values
        are a list of cell ids in that area where recharge will be applied. Cell id is the cell2d number.
        :param nper: number of stress periods for model
        :param recharges: dictionary where each key is an arbitary name for each recharge area. Must match the keys
        in the cell_ids dict. The dictionary values are each a list of recharge. Length of the list must equal to the
        number of stress periods.
        :param grid_type: string identifying grid type - 'disv' or 'disu'
        :return: recharge dictionary of stress period data to pass to flopy
        """
        rch_dict = {}
        nper = nper if self.nper is None else self.nper
        assert nper == len(list(recharges.values())[0]), 'Number of periods and length of recharge values must match'
        for per in range(nper):
            cell_list = []
            all_rch_cells = []
            for name, cell_nums in cell_ids.items():
                all_rch_cells += cell_nums
                recharge = recharges[name][per]
                for cell in cell_nums:
                    cell_id = cell if grid_type == 'disu' else (0, cell)
                    cell_list.append([cell_id, recharge])
            if background_rch is not None:
                for cell in range(self.vor.ncpl):
                    if cell not in all_rch_cells:
                        cell_id = cell if grid_type == 'disu' else (0, cell)
                        cell_list.append([cell_id, background_rch])
            rch_dict[per] = cell_list
        return rch_dict


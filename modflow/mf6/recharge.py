import flopy.utils.binaryfile
import shapely as shp
from pathlib import Path
import geopandas as gpd
import pandas as pd
from simple_modflow.modflow.mf6.headsplus import HeadsPlus as hp
import simple_modflow.modflow.mf6.mfsimbase as mf
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
            rch_fields: list | slice = None,
            rch_fields_to_pers: list = None

    ):
        """
        class to set up recharge for a modflow 6 model
        :param model: model to which this boundary applies
        :param vor: voronoi grid to which this boundary apples
        :param shp: path to shapefile that holds the polygons for the boundary
        :param uid: the field name in the shapefile attribute table that holds the unique ids, one for each polygon
        :param crs: coordinate reference system for boundary, should be integer EPSG code.
        :param rch_fields: field names corresponding to the recharge data in the shapefile attribute table
        :param rch_fields_to_pers: list of indices of length nper that correspond to the fields in rch_fields. Defines which field should be used for each stress period.
        :param bound_type: arbitary identifier for this boundary type
        """
        super().__init__(model, vor, shp, uid, crs)
        self.bound_type = 'rch'
        self.fields = rch_fields
        self.rch_fields_to_pers = rch_fields_to_pers
        self._cell_ids = None
        self._rch_fields = None
        self.uid = uid
        self._recharges = None

    @property
    def cell_ids(self):
        """gets cell ids for each recharge polygon in the shapefile as a dict. Keys are
        unique ids (uids)"""
        if self._cell_ids is None:
            cell_ids = self.intersections_no_duplicates
            cell_ids = cell_ids.set_index(self.uid)['no_dup'].to_dict()
            self._cell_ids = cell_ids
        return self._cell_ids

    @property
    def rch_fields(self):
        """gets a DataFrame of just the recharge data fields from the shapefile.
        Used to build the recharge dict for input into a flopy modflow model"""
        if self._rch_fields is None:
            rch_fields = self.gdf.loc[:, self.fields]
            self._rch_fields = pd.concat([self.gdf[self.uid], rch_fields], axis=1).set_index(self.uid)
        return self._rch_fields

    @property
    def recharges(self):
        if self._recharges is None:
            recharges = {}
            nper = self.nper
            uids = self.gdf[self.uid].to_list()
            rch_fields_dict = self.rch_fields.to_dict()
            rch_fields_cols = self.rch_fields.columns
            for uid in uids:
                recharges[uid] = []
            for per in range(nper):
                for uid in uids:
                    field_for_per = rch_fields_cols[self.rch_fields_to_pers[per]]
                    recharges[uid].append(rch_fields_dict[field_for_per][uid])
            self._recharges = recharges
        return self._recharges

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


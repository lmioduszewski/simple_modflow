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
inches_to_feet = 1 / 12

class GeneralHeadBoundary(Boundaries):

    def __init__(
            self,
            model: mf.SimulationBase = None,
            vor: Vor = None,
            shp: Path = None,
            line_shp: Path = None,
            uid: str = None,
            crs: int = 2927,
            rch_fields: list | slice = None,
            rch_fields_to_pers: list = None,
            xlsx_rch: Path = None,
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
        :param xlsx_rch: path to an excel file which contains the recharge data for each uid polygon. optional, otherwise data will be taken froim the shapefile attribute table. if excel is provided, it will be prioritized over the shapefile
        """
        super().__init__(model, vor, shp, uid, crs)
        self.bound_type = 'rch'
        self.fields = rch_fields
        self.rch_fields_to_pers = rch_fields_to_pers
        self.xlsx_rch = xlsx_rch
        self._cell_ids = None
        self._rch_fields = None
        self.uid = uid
        self._fields_to_pers = None
        self._recharges = None
        self.line_shp = line_shp
        self._line_ghb = None
        self.line_fields = {'cond': ['cond_strt', 'cond_end'], 'elev': ['elev_strt', 'elev_end'], 'layer': 'layer'}

    @property
    def cell_ids(self):
        """gets cell ids for each recharge polygon in the shapefile as a dict. Keys are
        unique ids (uids)"""
        if self._cell_ids is None:
            cell_ids = self.intersections_no_duplicates
            cell_ids = cell_ids['no_dup'].to_dict()
            self._cell_ids = cell_ids
        return self._cell_ids

    @property
    def rch_fields(self):
        """gets a DataFrame of just the recharge data fields from the shapefile or from excel file if proivded.
        Used to build the recharge dict for input into a flopy modflow model"""
        if self._rch_fields is None:
            if self.xlsx_rch:
                """use excel if it exists, otherwise get from shapefile"""
                rch_fields = pd.read_excel(self.xlsx_rch).set_index(self.uid)
                assert len(rch_fields) == len(self.gdf), 'number of rows in excel file and number of shapefile polys must be the equal'
            else:
                rch_fields = self.gdf.loc[:, self.fields]
            if self.rch_in_vol:  # if recharge is in volumes, divide by the voronoi area of each polygon
                rch_fields = rch_fields.apply(lambda row: row / self.gdf.area)
            self._rch_fields = rch_fields
        return self._rch_fields

    @property
    def fields_to_pers(self):
        """creates a list of indices and values that correspond to the columns/fields of the rch_fields DataFrame.
        In the case that the length of the fields is not long enough, -1 or -2 is added with correspond to either,
        use the last index given for the remaining stress periods or set the remaining stress periods to a background
        recharge, defined by setting the background recharge class attribute and setting apply_background_rch to True."""
        if self._fields_to_pers is None:
            fields_to_pers = self.rch_fields_to_pers.copy()
            fields_to_pers = [] if fields_to_pers is None else fields_to_pers
            for per in range(self.nper):
                if len(fields_to_pers) <= per:
                    if self.apply_background_rch:
                        fields_to_pers.append(-2)
                    else:
                        fields_to_pers.append(-1)
            self._fields_to_pers = fields_to_pers
        return self._fields_to_pers

    @property
    def recharges(self):
        """gets a dict of recharge values for each recharge area (each uid) for each stress period. Pass to
        get_rch() to generate a recharge dict to pass to the flopy recharge class."""
        nper = self.nper
        uids = self.gdf.index.to_list()
        scaled_rch_fields = self.rch_fields.mul(self.shp_to_vor_poly_scale, axis=0)
        rch_fields_dict = scaled_rch_fields.to_dict()
        rch_fields_cols = self.rch_fields.columns

        if self._recharges is None:
            recharges = {}
            for uid in uids:
                recharges[uid] = []
            for per in range(nper):
                if self.fields_to_pers[per] == -1:
                    for uid in uids:
                        #  if -1 then apply the last field indicated in the provided rch_fields_to_pers agrument
                        field_for_per = rch_fields_cols[self.rch_fields_to_pers[-1]]
                        recharges[uid].append(rch_fields_dict[field_for_per][uid])
                elif self.fields_to_pers[per] == -2:
                    for uid in uids:
                        recharges[uid].append(self.background_rch)
                else:
                    for uid in uids:
                        field_for_per = rch_fields_cols[self.fields_to_pers[per]]
                        recharges[uid].append(rch_fields_dict[field_for_per][uid])
            self._recharges = recharges
        return self._recharges

    @property
    def line_ghb(self):
        if self._line_ghb is None:
            assert isinstance(self.line_shp, Path), 'No line shapefile provided. Set the line_shp attribute'
            line_ghb = Boundaries(self.model, self.vor, self.line_shp, self.uid, self.crs, bound_type='ghb')
            self._line_ghb = line_ghb
        return self._line_ghb

    def add_line_ghb(self, existing_ghb_dict: dict=None):
        """adds a line based ghb boundary to an existing ghb dict"""
        elev_strt = self.line_ghb.gdf[self.line_fields['elev'][0]]
        elev_end = self.line_ghb.gdf[self.line_fields['elev'][1]]
        elev_delta = elev_end - elev_strt
        cond_strt = self.line_ghb.gdf[self.line_fields['cond'][0]]
        cond_end = self.line_ghb.gdf[self.line_fields['cond'][1]]
        cond_delta = cond_end - cond_strt
        ghb_line_lists = []
        for idx, delta in enumerate(elev_delta):
            line_len = self.line_ghb.gdf.geometry[idx].length
            sorted_cells = self.line_ghb.sorted_cells_along_line(idx)
            line_layer = self.line_ghb.gdf['layer'].iloc[idx] - 1
            for cell in sorted_cells.index:
                len_ratio = sorted_cells.loc[cell, 'distance_along_line'] / line_len
                elev_delta_ratio = elev_delta.iloc[idx] * len_ratio
                this_elev = elev_strt.iloc[idx] + elev_delta_ratio
                cond_delta_ratio = cond_delta.iloc[idx] * len_ratio
                this_cond = cond_strt.iloc[idx] + cond_delta_ratio
                this_cell = (line_layer, cell)
                this_ghb_list = [this_cell, this_elev, this_cond]
                ghb_line_lists.append(this_ghb_list)
        # TODO this only works assumiung constant ghb boundary. Adds the ghb data to stress period 1. Need to make more generic
        to_add = {0: ghb_line_lists}
        updated_ghb_dict = self.add_to_ghb_dict(existing_ghb_dict, to_add)
        return updated_ghb_dict

    def add_to_ghb_dict(self, existing_ghb_dict: dict=None, dict_to_add:dict =None):
        for per, ghb_list_of_lists in dict_to_add.items():
            assert isinstance(ghb_list_of_lists, list)
            updated_ghb_dict = existing_ghb_dict.copy()
            existing_list_of_lists: list = existing_ghb_dict[per]
            existing_cells = [existing_list[0] for existing_list in existing_list_of_lists]
            for ghb_list in ghb_list_of_lists:
                assert isinstance(ghb_list, list)
                assert ghb_list[0] not in existing_cells, (f'{ghb_list[0]} already in ghb dict. Remove duplicate cells.')
                updated_ghb_dict[per] = updated_ghb_dict[per] + [ghb_list]
        return updated_ghb_dict

    def get_rch(
            self,
            cell_ids: dict = None,
            recharges: dict = None,
            grid_type:str = 'disv',
            background_rch: int | float = None
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
        cell_ids = self.cell_ids if cell_ids is None else cell_ids
        recharges = self.recharges if recharges is None else recharges
        nper = self.nper
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


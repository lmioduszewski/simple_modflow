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


class GHB(Boundaries):

    def __init__(
            self,
            model: mf.SimulationBase = None,
            vor: Vor = None,
            shp_gpkg: Path = None,
            line_shp: Path = None,
            uid: str = 'name',
            crs: int = 2927,
            fields: list | slice = None,
            fields_to_pers: list = None,
            xlsx: Path = None,
            idomain: list[int] | pd.Series = None,
            idomain_path: Path = None,
            verbose: bool = False,
    ):
        """
        class to set up general head boundaries for a modflow 6 model
        :param model: model to which this boundary applies
        :param vor: voronoi grid to which this boundary apples
        :param shp_gpkg: path to shapefile that holds the polygons for the boundary
        :param uid: the field name in the shapefile attribute table that holds the unique ids, one for each polygon
        :param crs: coordinate reference system for boundary, should be integer EPSG code.
        :param fields: field names corresponding to the data in the shapefile attribute table or Excel file
        :param fields_to_pers: list of indices of length nper that correspond to the fields in fields.
        Defines which field should be used for each stress period.
        :param xlsx: path to an Excel file which contains the recharge data for each uid polygon.
        optional, otherwise data will be taken from the shapefile attribute table. if excel is provided,
        it will be prioritized over the shapefile
        """
        super().__init__(model, vor, shp_gpkg, uid, crs, idomain=idomain, idomain_path=idomain_path)
        self.bound_type = 'rch'
        self._fields = fields
        # self.fields_to_pers = fields_to_pers
        self.xlsx = xlsx
        self._cell_ids = None
        self._fields = None
        self.uid = uid
        self._fields_to_pers = None
        self._recharges = None
        self.line_shp = line_shp
        self._line_ghb = None
        self.line_fields = {'cond': ['cond_strt', 'cond_end'], 'elev': ['elev_strt', 'elev_end'], 'layer': 'layer'}
        self.verbose = verbose

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
    def fields(self):
        """gets a DataFrame of just the data fields from the shapefile or from excel file if provided.
        Used to build the dict for input into a flopy modflow model"""
        if self._fields is None:
            if self.xlsx:
                """use excel if it exists, otherwise get from shapefile"""
                fields = pd.read_excel(self.xlsx).set_index(self.uid)
                assert len(fields) == len(
                    self.gdf), 'number of rows in excel file and number of shapefile polys must be the equal'
            else:
                fields = self.gdf.loc[:, self._fields]
            self._fields = fields
        return self._fields

    @property
    def fields_to_pers(self):
        """creates a list of indices and values that correspond to the columns/fields of the fields DataFrame.
        In the case that the length of the fields is not long enough, -1 is added which tells the class to fill the
        remaining stress periods with the last field value given"""
        if self._fields_to_pers is None:
            fields_to_pers = self.fields_to_pers.copy()
            fields_to_pers: list = [] if fields_to_pers is None else fields_to_pers
            for per in range(self.nper):
                if len(fields_to_pers) <= per:
                    fields_to_pers.append(-1)
            self._fields_to_pers = fields_to_pers
        return self._fields_to_pers

    @property
    def data(self):
        """gets a dict of values for each boundary area (each uid) for each stress period. Pass to
        self.get() to generate a dict to pass to the flopy modflow boundary class."""
        nper = self.nper
        uids = self.gdf.index.to_list()
        scaled_rch_fields = self.fields.mul(self.shp_to_vor_poly_scale, axis=0)
        rch_fields_dict = scaled_rch_fields.to_dict()
        rch_fields_cols = self.fields.columns

        if self._recharges is None:
            recharges = {}
            for uid in uids:
                recharges[uid] = []
            for per in range(nper):
                if self.fields_to_pers[per] == -1:
                    for uid in uids:
                        #  if -1 then apply the last field indicated in the provided fields_to_pers argument
                        field_for_per = rch_fields_cols[self.fields_to_pers[-1]]
                        recharges[uid].append(rch_fields_dict[field_for_per][uid])
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

    def add_line_ghb(self, existing_ghb_dict: dict = None):
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
        updated_ghb_dict = self.add_to_dict(existing_ghb_dict, to_add)
        return updated_ghb_dict

    def get_from_poly(
            self,
            grid_type: str = 'disv',
            fields: dict = None,
            elev_reference: dict = None,
            reference_offset: float | int = 0
    ) -> dict:
        """
        Get a data dict for a flopy model from a polygon file (shapefile or geopackage)
        :param grid_type: default is 'disv'
        :param fields: a dict of custom field names in the geometry file, including 'name', 'height_over_btm',
        'conductance', 'layer', elevation', and 'min_elev'
        :return: a dict of drn data
        """
        nper = self.nper if self.nper is not None else 1
        if fields is None:
            fields = {
                'name': 'name',
                'elevation': 'elev',
                'height_over_btm': 'height',
                'conductance': 'cond',
                'layer': 'layer',
                'min_elev': 'min_elev',
            }
        # get bottoms of model layers from voronoi grid
        if self.vor is not None:
            lyr_botms = self.vor.gdf_topbtm.drop('geometry', axis='columns').iloc[:, 1:]
        else:
            lyr_botms = None
        # check to see if the geodataframe index has already been set to the correct 'name' field
        if self.gdf.index.name != fields['name']:
            gdf = self.gdf.set_index(fields['name'])
        else:
            gdf = self.gdf
        cells = self.edge_intersections.to_dict()
        boundary_dict = {}

        for per in range(nper):  # FOR EACH STRESS PERIOD
            cell_list = []

            for name, cell_nums in cells.items():  # FOR EACH GHB BOUDARY
                boundary_height = gdf.loc[name, fields['height_over_btm']]
                conductance = gdf.loc[name, fields['conductance']]
                layer = gdf.loc[name, fields['layer']]
                min_elev = gdf.loc[name, fields['min_elev']]
                elev = gdf.loc[name, fields['elevation']]
                min_elev = min_elev if min_elev is not None else 0
                # adjust layer number for zero-based indexing
                layer_idx = layer - 1

                for cell in cell_nums:  # FOR EACH CELL IN THIS GHB BOUNDARY
                    if self.inactive_cells is not None and cell in self.inactive_cells:
                        continue  # skip this if this cell is inactive
                    if elev is not None:
                        if name in elev_reference.keys():
                            try:
                                e = elev_reference[name][per]
                            except:
                                boundary_elev = elev + reference_offset
                                e = None
                            if isinstance(e, float | int):
                                if self.verbose:
                                    print(f'using {name} as an elevation reference for cell {cell}')
                                boundary_elev = e + reference_offset
                            else:
                                boundary_elev = elev + reference_offset
                        else:
                            boundary_elev = elev
                    elif boundary_height is not None:
                        boundary_elev = lyr_botms.iloc[cell, layer_idx] + boundary_height
                    else:
                        boundary_elev = lyr_botms.iloc[cell, layer_idx]
                    if boundary_elev < min_elev:
                        boundary_elev = min_elev  # adjusts drn elev to minimum allowed if specified
                    cell_id = cell if grid_type == 'disu' else (layer_idx, cell)
                    cell_list.append([cell_id, boundary_elev, conductance])

            boundary_dict[per] = cell_list

        return boundary_dict

    def add_to_dict(self, existing_dict: dict = None, dict_to_add: dict = None):
        assert self._verify_boundary_dict_structure(dict_to_add), 'dict to add not valid or missing'
        if existing_dict is not None:
            # 'period_list' is a list of list. Each list starts with a cellid and then the boundary data for that cellid
            for per, period_list in dict_to_add.items():
                # dict where cellids are the keys and boundary data are the values
                cell_dict = {cell_list[0]: cell_list[1:] for cell_list in period_list}
                ex_cell_dict = {ex_cell_list[0]: ex_cell_list[1:] for ex_cell_list in existing_dict[per]}
                # check if existing boundary dict already contains the cells to add/update
                for cell_id, vals in cell_dict.items():
                    if cell_id in ex_cell_dict.keys():
                        # if so, warns about overwrite, then overwrites
                        print(f'{cell_id} already in dict. Overwriting to {vals}.')
                    ex_cell_dict[cell_id] = vals
                # convert dict back to a nested list for this period in the loop
                existing_dict[per] = [[k] + v for k, v in ex_cell_dict.items()]

        else:
            existing_dict = dict_to_add
        # existing_dict has now been updated...so return
        return existing_dict

    def _verify_boundary_dict_structure(self, dict_to_check: dict = None):
        if dict_to_check is None:
            return dict_to_check
        for per, period_list in dict_to_check.items():
            assert per in range(self.model.nper), f'stress period {per} not in model'
            cell_dict = {cell_list[0]: cell_list[1:] for cell_list in period_list}
            for k in cell_dict.keys():
                assert k[0] in range(self.model.modelgrid.nlay), f'layer {k[0]} not in model grid'
                assert k[1] in range(self.model.modelgrid.ncpl), f'cell {k} not a valid model cell id'
            # TODO check cell_dict values based on type of boundary
        return True


if __name__ == '__main__':
    model_path_v4b_lowK = Path(r"C:\Users\lukem\mf6\cumb_v4b\cumb_v4b.model")
    with open(model_path_v4b_lowK, 'rb') as file:
        model = pickle.load(file)
    ghb_path = Path(r"C:\Users\lukem\mf6\Cumberland general\Boundaries\ghb.gpkg")
    ghb = GHB(shp_gpkg=ghb_path, model=model, crs=2926)
    print(ghb.get_from_poly()[0])


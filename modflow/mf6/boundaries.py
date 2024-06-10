import pandas as pd
from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor
from pathlib import Path
import numpy as np
from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
import geopandas as gpd

idxx = pd.IndexSlice
# Conversion factors
inches_to_feet = 1 / 12

def remove_duplicates(lst: list, seen: set = None):
    """Removes duplicates from a list"""
    seen = set() if seen is None else seen
    new_lst = []
    for num in lst:
        if num not in seen:
            new_lst.append(num)
            seen.add(num)
    return new_lst


class Boundaries:

    def __init__(
            self,
            model: SimulationBase = None,
            vor: Vor = None,
            shp: Path = None,
            uid: str = None,
            crs: int = 2927,
            bound_type: str = None

    ):
        """
        Base class for boundary conditions. Shouldn't need to instantiate. Instead use the boundary condition
        classes that inherit from this.
        :param model: model to which this boundary applies
        :param vor: voronoi grid to which this boundary apples
        :param shp: path to shapefile that holds the polygons for the boundary
        :param uid: the field name in the shapefile attribute table that holds the unique ids, one for each polygon
        :param crs: coordinate reference system for boundary, should be integer EPSG code.
        :param bound_type: arbitary identifier for this boundary type
        """

        self.model = model
        self.vor = vor
        self.bound_type = bound_type
        if shp is not None:
            self.gdf = gpd.read_file(shp)
            self.gdf.to_crs(inplace=True, epsg=crs)
        self.nper = model.nper
        self.uid = uid
        self._intersections = None
        self._intersections_no_duplicates = None
        self._vor_bound_polys = None
        self._rch_scale = None

    @property
    def intersections(self):
        """gets a DataFrame with unique ids (uid) for each shapefile polygon and the associated
        intersecting voronoi grid cells"""
        if self._intersections is None:
            vor_polys = self.vor.gdf_vorPolys
            df_intersect = self.gdf.geometry.apply(
                lambda geom: vor_polys[vor_polys.intersects(geom)].index.tolist())
            df_intersect = pd.concat([self.gdf[self.uid], df_intersect], axis=1)
            df_intersect.columns = [self.uid, 'intersect']
            self._intersections = df_intersect
        return self._intersections

    @property
    def intersections_no_duplicates(self):
        """gets a DataFrame of intersecting cells with duplicate cells removed"""
        if self._intersections_no_duplicates is None:
            seen = set()
            no_dups = self.intersections.copy()
            no_dups['len'] = no_dups['intersect'].apply(lambda x: len(x))
            no_dups['no_dup'] = no_dups.sort_values(by='len').loc[:, 'intersect'].apply(
                lambda x: remove_duplicates(x, seen))
            no_dups.drop(['len', 'intersect'], inplace=True, axis='columns')
            self._intersections_no_duplicates = no_dups
        return self._intersections_no_duplicates

    @property
    def vor_bound_polys(self):
        """gets the intersecting voronoi polygons equivalent to the shapefile polygons"""
        if self._vor_bound_polys is None:
            vor_polys = self.intersections_no_duplicates.copy()
            vor_polys['geometry'] = vor_polys['no_dup'].apply(lambda x: self.vor.gdf_vorPolys.loc[x].unary_union)
            self._vor_bound_polys = gpd.GeoDataFrame(vor_polys, geometry='geometry').drop(columns='no_dup')
        return self._vor_bound_polys

    @property
    def shp_to_vor_poly_scale(self):
        """gets a DataFrame giving the scaling between the areas of the shapefile vs. voronoi polys"""
        if self._rch_scale is None:
            rch_scale = self.gdf.area / self.vor_bound_polys.area
            rch_scale.index = self.gdf[self.uid]
            self._rch_scale = rch_scale
        return self._rch_scale

    def get_drn_stress_period_data(
        self, 
        cells: list,
        bottom_addition: float = 0,
        conductance: float = 100,
        disMf: str = 'disu',
        bottoms: dict = None,
        ) -> list:
        """Returns a list of lists. Each nested list corresponds to the DRN package
        boundary data for a particular voronoi cell in the grid, which includes cell
        ID, elevation of drain, and conductance. Can be passed to the flopy DRN package.

        Args:
            cells (list): list of cell IDs in this drain
            bottoms (dict): dict of bottoms, keys are the cell indices, values are the bottom elevations
            bottom_addition (float): height above the bottom of cell for the drain. This is added to the bottom of cell elevation derived from the Voronoi grid object.
            conductance (float): conductance for this drain cell
            disMf (str, optional): Either 'disu' or 'disv' works, and refers to the MODFLOW6 discretization package being used. Defaults to 'disu'.

        Returns:
            list: List of lists that contain the data for this drain and can be passed to the flopy DRN package
        """
        
        if self.vor is None:
            return print("No voronoi grid defined")
        drn_values = []
        for cell in cells:
            if bottoms:
                thisdrn = [cell, (bottoms[cell] + bottom_addition), conductance]
            elif self.vor.gdf_topbtm is not None:
                try:
                    thisdrn = [cell, (self.vor.gdf_topbtm.loc[cell, "bottom"] + bottom_addition), conductance]
                except:
                    print("can't get bottom elevations for drains. Assuming bottom elev is zero")
                    thisdrn = [cell, bottom_addition, conductance]
            else:
                thisdrn = [cell, bottom_addition, conductance]
            if disMf == "disv":
                thisdrn = [0] + thisdrn  # add layer num for disv grid
            drn_values.append(thisdrn)
        return drn_values


    def get_drn_from_shp(
            self,
            shapefile_path: Path,
            grid_type: str = 'disv',
            nper: int = 1,
            fields: dict = None
    ):
        nper = nper if self.nper is None else self.nper
        if fields is None:
            fields = {
                'name': 'name',
                'height_over_btm': 'height',
                'conductance': 'cond',
                'layer': 'layer'
            }
        # get drain cells based on shapefile
        drn_cells, gdf_drn = self.vor.get_vor_cells_as_dict(
            locs=shapefile_path,
            predicate='intersects',
            loc_name_field=fields['name'],
            return_gdf=True
        )
        # get bottoms of model layers from voronoi grid
        lyr_botms = self.vor.gdf_topbtm.drop('geometry', axis='columns').iloc[:, 1:]
        gdf_drn = gdf_drn.set_index(fields['name'])
        drn_dict = {}
        for per in range(nper):
            cell_list = []
            for name, cell_nums in drn_cells.items():
                boundary_height = gdf_drn.loc[name, fields['height_over_btm']]
                conductance = gdf_drn.loc[name, fields['conductance']]
                layer = gdf_drn.loc[name, fields['layer']]
                # adjust layer number for zero-based indexing
                layer_idx = layer - 1
                for cell in cell_nums:
                    boundary_elev = lyr_botms.iloc[cell, layer_idx] + boundary_height
                    cell_id = cell if grid_type == 'disu' else (layer_idx, cell)
                    cell_list.append([cell_id, boundary_elev, conductance])
            drn_dict[per] = cell_list
        return drn_dict

    def get_rch_dict(
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
        print(nper)
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

    def get_ghb_from_shp(
            self,
            shapefile_path: Path,
            grid_type: str = 'disv',
            nper: int = 1,
            fields: dict = None
    ):
        """
        Returns a dictionary to use as input into the ghb flopy constructor. keys of the dict are stress periods.
        :param shapefile_path: path to shapefile with ghb information
        :param grid_type: disv or disu
        :param nper: number of stress periods in the model
        :param fields: field names in the shapefile corresponding to name, elevation, conductance, and layer of the ghb
        :return: dict where keys are stress periods and values are the ghb data for flopy
        """
        if fields is None:
            fields = {
                'name': 'name',
                'elevation': 'elev',
                'conductance': 'cond',
                'layer': 'layer'
            }
        ghb_cells, gdf_ghb = self.vor.get_vor_cells_as_dict(
            locs=shapefile_path,
            predicate='intersects',
            loc_name_field=fields['name'],
            return_gdf=True
        )
        gdf_ghb = gdf_ghb.set_index(fields['name'])
        ghb_dict = {}
        for per in range(nper):
            cell_list = []
            for name, cell_nums in ghb_cells.items():
                boundary_head = gdf_ghb.loc[name, fields['elevation']]
                conductance = gdf_ghb.loc[name, fields['conductance']]
                layer = gdf_ghb.loc[name, fields['layer']]
                layer_idx = layer - 1
                for cell in cell_nums:
                    cell_id = cell if grid_type == 'disu' else (layer_idx, cell)
                    cell_list.append([cell_id, boundary_head, conductance])
            ghb_dict[per] = cell_list
        return ghb_dict

    def get_chd_from_shp(
            self,
            shapefile_path: Path,
            grid_type: str = 'disv',
            nper: int = 1,
            fields: dict = None
    ):
        """
        Returns a dictionary to use as input into the chd flopy constructor. keys of the dict are stress periods.
        :param shapefile_path: path to shapefile with chd information
        :param grid_type: disv or disu
        :param nper: number of stress periods in the model
        :param fields: field names in the shapefile corresponding to name, elevation, and layer of the chd
        :return: dict where keys are stress periods and values are the chd data for flopy
        """
        if fields is None:
            fields = {
                'name': 'name',
                'elevation': 'elev',
                'layer': 'layer'
            }
        chd_cells, gdf_chd = self.vor.get_vor_cells_as_dict(
            locs=shapefile_path,
            predicate='intersects',
            loc_name_field=fields['name'],
            return_gdf=True
        )
        gdf_chd = gdf_chd.set_index(fields['name'])
        chd_dict = {}
        for per in range(nper):
            cell_list = []
            for name, cell_nums in chd_cells.items():
                boundary_head = gdf_chd.loc[name, fields['elevation']]
                layer = gdf_chd.loc[name, fields['layer']]
                layer_idx = layer - 1
                for cell in cell_nums:
                    cell_id = cell if grid_type == 'disu' else (layer_idx, cell)
                    cell_list.append([cell_id, boundary_head])
            chd_dict[per] = cell_list
        return chd_dict

    def get_k_from_shp(
            self,
            shapefile_path: Path,
            grid_type: str = 'disv',
            fields: dict = None,
            nlay: int = 1,
            return_array: bool = True
    ):
        if fields is None:
            fields = {
                'name': 'name',
                'k': 'k',
                'layer': 'layer'
            }
        k_cells, gdf_k = self.vor.get_vor_cells_as_dict(
            locs=shapefile_path,
            predicate='intersects',
            loc_name_field=fields['name'],
            return_gdf=True
        )
        gdf_k = gdf_k.set_index(fields['name'])

        k_midx = pd.MultiIndex.from_product(
            iterables=[list(range(nlay)), list(range(self.vor.ncpl))],
            names=['layer', 'cell'])
        k_df = pd.DataFrame(index=k_midx, columns=['k'])
        k_lists = []
        for name, cell_nums in k_cells.items():
            k = gdf_k.loc[name, fields['k']]
            layer = gdf_k.loc[name, fields['layer']]
            layer_idx = layer - 1
            k_df.loc[idxx[layer_idx, cell_nums], 'k'] = k
        if return_array:
            for layer in range(nlay):
                k_lists.append(k_df.loc[layer].squeeze().tolist())
            k_array = np.array(k_lists)
            return k_array
        else:
            return k_df
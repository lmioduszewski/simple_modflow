from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
    from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor

import pandas as pd
# from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor
from pathlib import Path
import numpy as np
# from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
import geopandas as gpd
import shapely as shp

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
            shp_gpkg: Path = None,
            uid: str = None,
            crs: int = None,
            bound_type: str = None,
            idomain: list[int] | pd.Series = None

    ):
        """
        Base class for boundary conditions. Shouldn't need to instantiate. Instead use the boundary condition
        classes that inherit from this.
        :param model: model to which this boundary applies
        :param vor: voronoi grid to which this boundary apples
        :param shp_gpkg: path to shapefile that holds the polygons for the boundary
        :param uid: the field name in the shapefile attribute table that holds the unique ids, one for each polygon. required
        :param crs: coordinate reference system for boundary, should be integer EPSG code.
        :param bound_type: arbitary identifier for this boundary type
        :param idomain: list of integers (1 or 0), one for each cell in the grid. If 0, that cell index is inactive
        """

        self.model = model
        try:
            self.vor = self.model.vor if vor is None else vor
            self.crs = self.vor.crs if crs is None else crs
            self.nper = self.model.nper
        except:
            self.vor = None
            self.nper = None
        self.bound_type = bound_type
        self.uid = uid
        self.crs = crs
        if shp_gpkg is not None:
            self.gdf = gpd.read_file(shp_gpkg)
            if self.uid is not None:
                self.gdf = self.gdf.set_index(self.uid)
            self.gdf.to_crs(inplace=True, epsg=crs)
        self._intersections = None
        self._edge_intersections = None
        self._intersections_no_duplicates = None
        self._vor_bound_polys = None
        self._rch_scale = None
        self._sorted_cells_along_line = None
        self._inactive_cells = None
        self.idomain = idomain

    @property
    def inactive_cells(self):
        if self._inactive_cells is None:
            if self.idomain is None:
                return None
            elif self.idomain is not None and self.vor is not None:
                assert len(self.idomain) == self.vor.ncpl, 'idomain length must be equal to num cells in vor grid'
                if isinstance(self.idomain,list):
                    self.idomain = pd.Series(self.idomain)
                assert isinstance(self.idomain,pd.Series), 'error: could not make idomain a pd.Series'
                inactive_cells = self.idomain[self.idomain == 0].index.tolist()
                self._inactive_cells = inactive_cells
            else:
                print('no active grid object provided. Cannot determine inactive cells')
        return self._inactive_cells

    @property
    def intersections(self):
        """gets a DataFrame with unique ids (uid) for each shapefile polygon and the associated
        intersecting voronoi grid cells"""
        if self._intersections is None:
            vor_polys = self.vor.gdf_vorPolys
            df_intersect = self.gdf.geometry.apply(
                lambda geom: vor_polys[vor_polys.intersects(geom)].index.tolist())
            df_intersect.name = 'intersect'
            self._intersections = df_intersect
        return self._intersections

    @property
    def edge_intersections(self):
        """gets a DataFrame with unique ids (uid) for each shapefile polygon and the associated
        intersecting voronoi grid cells, but then filters for only those on a grid edge"""
        if self._edge_intersections is None:
            edge_cells = self.vor.get_grid_edge()
            intersections = self.intersections
            filtered = intersections.apply(lambda x: [cell for cell in x if cell in edge_cells])
            self._edge_intersections = filtered
        return self._edge_intersections

    @property
    def intersections_no_duplicates(self):
        """gets a DataFrame of intersecting cells with duplicate cells removed"""
        if self._intersections_no_duplicates is None:
            seen = set()
            no_dups = self.intersections.copy()
            lens = no_dups.apply(lambda x: len(x)).sort_values()
            lens.name = 'len'
            no_dups = pd.concat([lens, no_dups], axis=1)
            no_dups['no_dup'] = no_dups.loc[:, 'intersect'].apply(
                lambda x: remove_duplicates(x, seen))
            no_dups.drop(['len', 'intersect'], inplace=True, axis='columns')
            self._intersections_no_duplicates = no_dups
        return self._intersections_no_duplicates

    @property
    def vor_bound_polys(self):
        """gets the intersecting voronoi polygons equivalent to the shapefile polygons"""
        if self._vor_bound_polys is None:
            vor_polys = self.intersections_no_duplicates.copy()
            vor_polys['geometry'] = vor_polys['no_dup'].apply(lambda x: self.vor.gdf_vorPolys.loc[x].union_all())
            self._vor_bound_polys = gpd.GeoDataFrame(vor_polys, geometry='geometry').drop(columns='no_dup')
        return self._vor_bound_polys

    @property
    def shp_to_vor_poly_scale(self):
        """gets a DataFrame giving the scaling between the areas of the shapefile vs. voronoi polys"""
        if self._rch_scale is None:
            rch_scale = self.gdf.area / self.vor_bound_polys.area
            self._rch_scale = rch_scale
        return self._rch_scale

    def sorted_cells_along_line(self, idx=0):
        """
        method to get a list of voronoi model cells along the length of a line in order from
        the beginning to end of the linestring.
        :param idx: index of the line, if there was more than one provided in the shapefile
        :return: dataframe of intersecting cells, sorted with distance along line and centroids
        """

        line_gdf = self.gdf
        line = line_gdf.geometry[idx]
        assert isinstance(line, shp.geometry.linestring.LineString), 'geometry not a LineString!'
        intersecting_cells = self.vor.gdf_vorPolys.loc[self.intersections.to_list()[idx], :]
        intersecting_cells['centroid'] = intersecting_cells.centroid

        # Calculate the distance along the line for each cell centroid
        intersecting_cells['distance_along_line'] = intersecting_cells.centroid.apply(
            lambda point: line.project(point)
        )
        # Sort the cells based on the distance along the line
        sorted_cells = intersecting_cells.sort_values('distance_along_line')
        sorted_cells.index.name = 'cell'

        return sorted_cells

    def get_drn_stress_period_data(
            self, *args, **kwargs
    ):

        print('Deprecated. Use method of same name from the DRN class')

        return NotImplementedError

    def get_drn_from_shp(
            self, *args, **kwargs
    ):

        print('Deprecated. Use method of same name from the DRN class')

        return NotImplementedError

    def get_rch_dict(
            self,
            zone_cell_id_dict: dict = None,
            zone_rch_dict: dict = None,
            grid_type: str = 'disv',
            background_rch: int | float = None,
            nper: int = 1,
            shift: int = 0
    ) -> dict:
        """
        get a recharge dictionary to pass to flopy in setting of a recharge package. Assumes recharge only applied to
        top layer
        :param zone_cell_id_dict: dictionary where each key is an arbitrary name given each recharge area and the values
        are a list of cell ids in that area where recharge will be applied. Cell id is the cell2d number.
        :param nper: number of stress periods for model
        :param zone_rch_dict: dictionary where each key is an arbitary name for each recharge area. Must match the keys
        in the rch_zone_dict dict. The dictionary values are each a list of recharge. Length of the list must equal to the
        number of stress periods.
        :param shift: number of stress periods to shift each recharge, to delay it if needed
        :param grid_type: string identifying grid type - 'disv' or 'disu'
        :return: recharge dictionary of stress period data to pass to flopy
        """
        rch_dict = {}
        nper = nper if self.nper is None else self.nper
        assert nper == len(list(zone_rch_dict.values())[0]), "nper and length of rch dict must match"
        for per in range(nper):
            if per + shift >= nper:
                continue
            cell_list = []
            all_rch_cells = []
            for name, cell_nums in zone_cell_id_dict.items():
                if isinstance(cell_nums, int):
                    cell_nums = [cell_nums]
                cell_nums = [cell for cell in cell_nums if cell not in self.inactive_cells]
                all_rch_cells += cell_nums
                recharge = zone_rch_dict[name][per]
                for cell in cell_nums:
                    cell_id = cell if grid_type == 'disu' else (0, cell)
                    cell_list.append([cell_id, recharge])
            if background_rch is not None:
                back_cells = [cell for cell in list(range(self.vor.ncpl)) if cell not in self.inactive_cells]
                for cell in back_cells:
                    if cell not in all_rch_cells:
                        cell_id = cell if grid_type == 'disu' else (0, cell)
                        cell_list.append([cell_id, background_rch])

            if per == 0 and shift > 0:
                for s in range(shift):
                    rch_dict[s] = cell_list  # duplicate the first stress period to meet specified shift
            rch_dict[per + shift] = cell_list
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
        crs = self.vor.crs if self.vor is not None else None
        k_cells, gdf_k = self.vor.get_vor_cells_as_dict(
            locs=shapefile_path,
            predicate='intersects',
            loc_name_field=fields['name'],
            return_gdf=True,
            crs=crs
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

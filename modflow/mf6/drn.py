import flopy.utils.binaryfile
import shapely as shp
from pathlib import Path
import geopandas as gpd
import pandas as pd
import simple_modflow.modflow.mf6.mfsimbase as mf
from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor
import pickle
from simple_modflow.modflow.mf6.boundaries import Boundaries

idxx = pd.IndexSlice
inches_to_feet = 1 / 12


class DRN(Boundaries):

    def __init__(
            self,
            model: mf.SimulationBase = None,
            vor: Vor = None,
            shp_gpkg: Path = None,
            uid: str = None,
            crs: int = 2927,
            grid_type: str = 'disv',
    ):
        super().__init__(model, vor, shp_gpkg, uid, crs)
        self.bound_type = 'drn'
        self.grid_type = grid_type.lower()




    def get_drn_stress_period_data(
            self,
            cells: list,
            bottom_addition: float | int = 0,
            conductance: float | int | list = 100,
            disMf: str = 'disv',
            bottoms: dict = None,
            layer: int = None,
    ) -> list:
        """Returns a list of lists. Each nested list corresponds to the DRN package
        boundary data for a particular voronoi cell in the grid, which includes cell
        ID, elevation of drain, and conductance. Can be passed to the flopy DRN package.

        Args:
            cells (list): list of cell IDs in this drain
            bottoms (dict): dict of bottoms, keys are the cell indices, values are the bottom elevations
            bottom_addition (float): height above the bottom of cell for the drain. This is added to the bottom of cell elevation derived from the Voronoi grid object.
            conductance (float | int | list): conductance for the cell of this DRN, can also provide a list
            of conductance values that is of equal length as the list of cells. Function will create a
            dict of cell:conductance key:value pairs.
            disMf (str, optional): Either 'disu' or 'disv' works, and refers to the MODFLOW6 discretization package being used. Defaults to 'disu'.
            layer: layer to apply drains

        Returns:
            list: List of lists that contain the data for this drain and can be passed to the flopy DRN package
        """

        if self.vor is None:
            return print("No voronoi grid defined")
        drn_values = []

        # Set up conductance dict
        if isinstance(conductance, int | float):
            conductance = {cell: conductance for cell in cells}
        elif isinstance(conductance, list):
            assert len(conductance) == len(cells), 'conductance list length must equal number of cells in DRN'
            conductance = dict(zip(cells, conductance))

        # generate list of lists for all DRN boundary cells
        for cell in cells:
            cell_id = cell if disMf == 'disu' else (layer, cell)
            if bottoms:
                thisdrn = [cell_id, (bottoms[cell] + bottom_addition), conductance[cell]]
            elif self.vor.gdf_topbtm is not None:
                try:
                    #  need a better way
                    thisdrn = [cell_id, (self.vor.gdf_topbtm.loc[cell, layer + 1] + bottom_addition), conductance [cell]]
                except:
                    print("can't get bottom elevations for drains. Assuming bottom elev is zero")
                    thisdrn = [cell_id, bottom_addition, conductance[cell]]
            else:
                thisdrn = [cell_id, bottom_addition, conductance[cell]]
            """if disMf == "disv":
                thisdrn = [0] + thisdrn  # add layer num for disv grid"""
            drn_values.append(thisdrn)
        return drn_values

    def get_drn_from_poly(
            self,
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
        """# get drain cells based on shapefile
        drn_cells, gdf_drn = self.vor.get_vor_cells_as_dict(
            locs=shapefile_path,
            predicate='intersects',
            loc_name_field=fields['name'],
            return_gdf=True
        )"""
        # get bottoms of model layers from voronoi grid
        lyr_botms = self.vor.gdf_topbtm.drop('geometry', axis='columns').iloc[:, 1:]
        # check to see if the geodataframe index has already been set to the correct 'name' field
        if self.gdf.index.name != fields['name']:
            gdf_drn = self.gdf.set_index(fields['name'])
        else:
            gdf_drn = self.gdf
        drn_cells = self.edge_intersections.to_dict()
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



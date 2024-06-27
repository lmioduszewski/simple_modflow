import numpy as np
import pandas as pd
from .mf2Dplots import WaterLevelPlot
from .voronoiplus import VoronoiGridPlus as vgp
import flopy.utils.binaryfile as bf
from pathlib import Path
import plotly.graph_objs as go
import geopandas as gpd
from . import mf2Dplots
import figs
from scipy.interpolate import RBFInterpolator

idxx = pd.IndexSlice  # for easy index slicing in a MultiIndex DataFrame
crs_latlon = "EPSG:4326"

class HeadsPlus(bf.HeadFile):

    def __init__(
            self,
            hds_path: Path,
            vor: vgp = None,
            obs_path: Path = None
    ):
        """
        Class to do stuff with a MODFLOW heads file. This class subclasses the flopy.utils.binaryfile.HeadFile class.
        :param hds_path: path to the heads file
        :param vor: voronoi grid representing model grid for the heads file, optional
        """
        super().__init__(filename=hds_path)

        self.hds = bf.HeadFile(filename=hds_path)
        self.kstpkper = self.get_kstpkper()
        self.vor = vor

        self.modebar = go.layout.Modebar(
            add=[
                'togglespikelines',
                'hovercompare',
                'togglehover',
                'drawline',
                'drawopenpath',
                'drawclosedpath',
                'drawcircle',
                'drawrect',
                'eraseshape'
            ]
        )
        self.fig_layout = go.Layout(  # default fig layout
            dragmode='pan',
            modebar=self.modebar
        )
        self.fig_config = {  # default fig config
            'scrollZoom': True
        }
        self.vor = vor
        self.obs_heads_df = None
        self._all_heads = None
        self.nper = pd.DataFrame(self.hds.get_kstpkper()).iloc[:, 1].max() + 1
        self.numstp = pd.DataFrame(self.hds.get_kstpkper()).iloc[:, 0].max() + 1
        self.vor_list = self.vor.gdf_vorPolys.geometry.to_list()
        self.cell_list = [i for i in range(len(self.vor_list))]
        self.area_list = [cell.area for cell in self.vor_list]
        self.x_list = [cell.centroid.xy[0][0] for cell in self.vor_list]
        self.y_list = [cell.centroid.xy[1][0] for cell in self.vor_list]
        self._obs = {}
        self._obs_heads = None
        self.obs_path = obs_path


    @property
    def all_heads(self):
        if self._all_heads is None:
            self._all_heads = self.get_all_heads()
        return self._all_heads

    @property
    def obs(self):
        return self._obs

    @property
    def obs_path(self):
        return self._obs_path

    @obs_path.setter
    def obs_path(self, val):
        self._obs_path = val

    @property
    def obs_heads(self):
        if self._obs_heads is None:
            self._obs_heads = self.get_obs_heads()
        return self._obs_heads

    def get_all_heads(self):
        """Method to get all heads for this model and store in
            a dataframe"""

        vor_cell_list = list(self.vor.gdf_vorPolys.index)

        """set generic MultiIndex for all stress periods and all cells"""
        hds_mdx = pd.MultiIndex.from_product(
            iterables=[
                self.kstpkper,
                list(range(self.nlay)),
                vor_cell_list
            ],
            names=['kstpkper', 'layer', 'cell']
        )
        """Set up a MultiIndex DataFrame to hold the heads
            for all cells and stress periods"""
        df_heads = pd.DataFrame(
            index=hds_mdx,
            columns=['elev']
        )
        """get data for each stress period"""
        for kstpkper in self.kstpkper:
            spHds = pd.DataFrame(self.get_data(kstpkper=kstpkper).squeeze().transpose())
            """copy and paste this stress period data to the MultiIndex DataFrame"""
            for layer in range(self.nlay):
                df_heads.loc[idxx[kstpkper, layer], :] = spHds.iloc[:, layer].values

        return df_heads

    @staticmethod
    def sort_dict_by_keys(
            dict_to_sort: dict = None
    ):
        """Sorts the given dict by its keys and returns the sorted dict"""
        sorted_keys = sorted(dict_to_sort.keys())
        sorted_dict = {i: dict_to_sort[i] for i in sorted_keys}
        return sorted_dict

    def get_obs_cells(self, locs: Path, crs: str = "EPSG:2927", loc_name_field='ExploName'):
        """
        method to get cells that contain certain observation locations. Locations
        should be points in a shapefile
        :param locs: Path for shapefile with locations of obs as points
        :param crs: crs of shapefile. Defaults to EPSG_2927 (South WA)
        :param loc_name_field: field name in the shapefile attribute table containing observation names.
        :return: dict where observation names are keys and the lists of cells containing them are the values.
        """
        locs = self.obs_path if locs is None else locs
        if locs is None:
            return print('no obs path found')
        obs_dict = self.vor.get_vor_cells_as_dict(
            locs=locs,
            crs=crs,
            predicate='contains',
            loc_name_field=loc_name_field)
        obs_dict = self.sort_dict_by_keys(obs_dict)
        for obs, cell_ids in obs_dict.items():
            assert len(cell_ids) == 1, f'more than one cell found for {obs}. Fix to make it one cell'
            obs_dict[obs] = cell_ids[0]
        return obs_dict

    def get_obs_heads(self, locs: Path = None, crs: str = "EPSG:2927", loc_name_field='ExploName'):

        if locs is not None:
            new_obs = self.get_obs_cells(locs, crs, loc_name_field)
            self._obs.update(new_obs)
        if not self.obs:
            d = self.get_obs_cells(self.obs_path, crs, loc_name_field)
            self._obs.update(d)
        obs_cell_dict = self.obs
        assert obs_cell_dict, 'no observations found'
        obs_cells = list(obs_cell_dict.values())
        all_heads = self.all_heads.copy()
        obs_heads = all_heads.loc[idxx[:, :, obs_cells], :]
        obs_reset_idx = obs_heads.reset_index()
        for ob, cell in obs_cell_dict.items():
            obs_reset_idx.loc[obs_reset_idx['cell'] == cell, 'cell'] = ob
        obs_heads = obs_reset_idx.pivot(
            index=['layer', 'kstpkper'],
            columns='cell',
            values='elev'
        )
        obs_heads.columns.name = 'obs'

        return obs_heads


    def plot_heads(
            self,
            locs: Path,
            crs: str = "EPSG:2927",
            normalized_head_per_obs: dict = None,
            to_date_range: tuple = None,
            loc_name_field='ExploName'
    ):
        """Plots heads for specified locations in the model.
            Locations should be specified as the model cell/node to 
            plot

            Args:
                locs (Path, optional): Path to shapefile of points to plots heads. Defaults to
                 None.
                normalized_head_per_obs (Dict, optional): Dict where the keys are the loc names
                 and the values are the elevations to normalize to
                to_date_range (tuple, optional): tuple of starting and ending dates. Assumes
                daily data.
            Returns:
                Returns the fig object and plots it
            """

        if locs is not None:

            heads = self.get_all_heads().to_dict()['elev']
            obs_dict = self.vor.get_vor_cells_as_dict(locs=locs,
                                                      crs=crs,
                                                      predicate='contains',
                                                      loc_name_field=loc_name_field
                                                      )
            obs_dict = self.sort_dict_by_keys(obs_dict)

            obs_heads_dict = {}  # define the empty dict

            for obs in obs_dict:  # add all heads for obs to the new dict
                obs_heads_dict[obs] = []
                for stp_per in self.kstpkper:
                    for layer in range(self.nlay):
                        obs_heads_dict[obs] += [heads[(stp_per, layer, obs_dict[obs][0])]]
            fig = go.Figure(
                layout=self.fig_layout,
            )
            for obs in obs_dict:
                if normalized_head_per_obs is None:
                    y_norm = obs_heads_dict[obs]
                else:
                    this_normalized_obs_head = normalized_head_per_obs[obs]
                    y_norm = np.array(obs_heads_dict[obs]) - np.array(
                        [this_normalized_obs_head] * len(obs_heads_dict[obs]))
                    y_norm = list(y_norm)
                if to_date_range is None:
                    x = list(range(len(self.kstpkper)))
                else:
                    x = pd.date_range(to_date_range[0], to_date_range[1])
                fig.add_scattergl(
                    y=y_norm,
                    x=x,
                    name=obs
                )
            df = pd.DataFrame(obs_heads_dict)
            try:  # save the sorted obs heads dict to self
                df = df[sorted(df.columns)]
                self.obs_heads_df = df
            except:
                self.obs_heads_df = df

            fig.show(renderer='browser', config=self.fig_config)
            return fig

        elif locs is None:
            return None

    def plot_choropleth(
            self,
            stp_per_to_plot: tuple = (0, 0),
            plot_mounding: bool = False,
            zmin=None,
            zmax=None,
            zoom=18,
            custom_hover: dict = None,
            bottom=None,
            bottom_array=None,
            all_layers=False,
            layer=1,
            obs: Path = None,
            obs_name: str = 'ExploName'
    ):
        """Plot heads for a specified time step and stress period on a
            choropleth map. Heads may show saturated thickness (mounding) or
            show elevation head.

            :param stp_per_to_plot (tuple, optional): Tuple defining time step and stress period to plot. Defaults to (0,0).
            :param plot_mounding (boolean, optional): Boolean to determine whether to plot elevation head or mounding (sat thickness)
            :param zoom: define zoom level of plot. Default is 18.
            """

        stp_to_plot = stp_per_to_plot[0]
        per_to_plot = stp_per_to_plot[1]
        kstpkper_key = f"sp{per_to_plot}ts{stp_to_plot}"
        choro_dict = {}
        choro_heads = {}
        vor = self.vor

        """If plot_mounding == True then plot head over cell bottom"""
        if bottom:
            bottom_elev = bottom
        elif plot_mounding:
            if bottom_array is not None:
                bottom_elev = bottom_array
            else:
                bottom_elev = vor.gdf_topbtm.loc[:, layer].to_numpy()
        else:
            bottom_elev = 0

        if all_layers is True:
            pass
        choro_heads[kstpkper_key] = self.all_heads.loc[idxx[(stp_to_plot, per_to_plot), layer], :]
        choro_dict[kstpkper_key] = choro_heads[kstpkper_key]['elev'] - bottom_elev

        if zmax is None:
            zmax = choro_dict[kstpkper_key].max()
        if zmin is None:
            zmin = choro_dict[kstpkper_key].min()

        fig_mbox = mf2Dplots.ChoroplethPlot(vor=vor, zoom=zoom)

        # create lists for hover data and get hover template
        head_list = choro_heads[kstpkper_key]['elev'].to_list()
        hover_dict = {
            'Cell No.': self.cell_list,
            'Area': self.area_list,
            'x': self.x_list,
            'y': self.y_list,
        }
        for lyr in range(self.nlay):
            lyr_heads = self.all_heads.loc[idxx[(stp_to_plot, per_to_plot), lyr], 'elev'].to_list()
            hover_dict[f'Layer {lyr+1} Heads'] = lyr_heads
        if plot_mounding:
            mounding_list = choro_dict[kstpkper_key].to_list()
            hover_dict['Mounding'] = mounding_list
        if custom_hover:
            for name, data in custom_hover.items():
                hover_dict[str(name)] = data

        custom_data, hover_template = figs.create_hover(hover_dict)

        fig_mbox.add_choroplethmapbox(
            geojson=vor.latslons,
            featureidkey="id",
            locations=vor.gdf_latlon.index.to_list(),
            z=choro_dict[kstpkper_key],
            hovertemplate=hover_template,
            customdata=custom_data,
            colorscale="earth",
            zmax=zmax,
            zmin=zmin,
        )
        if obs:
            obs = gpd.read_file(obs).to_crs(crs_latlon)
            fig_mbox.add_scattermapbox(
                lat=obs.geometry.y,
                lon=obs.geometry.x,
                text=obs[obs_name],
                hoverinfo='text',
                marker_color='red'
            )

        return fig_mbox

    def interpolate_and_save_raster(
            self, x, y, z, output_tif='interpd.tif', grid_res=100, save=True):
        """
        Interpolates irregular x, y points with associated z-values to a gridded surface
        and writes the output to a GeoTIFF raster with the specified CRS.

        Parameters:
        x (array-like): Array of x coordinates.
        y (array-like): Array of y coordinates.
        z (array-like): Array of z values.
        output_tif (str): Path to the output GeoTIFF file.
        crs (str): Coordinate reference system (CRS) in EPSG format (e.g., 'EPSG:2927').
        grid_res (int): Resolution of the grid (number of grid points along each axis).

        Returns:
        None
        """
        # Define the extent of the grid
        xmin, ymin = np.min(x), np.min(y)
        xmax, ymax = np.max(x), np.max(y)
        crs = self.vor.crs

        # Create grid coordinates
        grid_x, grid_y = np.meshgrid(
            np.linspace(xmin, xmax, grid_res),
            np.linspace(ymin, ymax, grid_res)
        )

        # Interpolate data using RBFInterpolator
        interpolator = RBFInterpolator(np.column_stack((x, y)), z, neighbors=10)
        grid_z = interpolator(np.column_stack((grid_x.ravel(), grid_y.ravel()))).reshape(grid_x.shape)

        # Define the transform and CRS
        xres = (xmax - xmin) / (grid_res - 1)
        yres = (ymax - ymin) / (grid_res - 1)
        transform = from_bounds(xmin, ymin, xmax, ymax, grid_res, grid_res)

        # Flip the grid_z array vertically
        grid_z = np.flipud(grid_z)

        # Save the grid as a GeoTIFF
        if save:
            with rasterio.open(output_tif, 'w', driver='GTiff',
                               height=grid_res, width=grid_res,
                               count=1, dtype=grid_z.dtype,
                               crs=crs, transform=transform) as dst:
                dst.write(grid_z, 1)

        return grid_x, grid_y, grid_z

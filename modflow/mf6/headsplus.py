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

idxx = pd.IndexSlice  # for easy index slicing in a MultiIndex DataFrame
crs_latlon = "EPSG:4326"


class HeadsPlus(bf.HeadFile):

    def __init__(
            self,
            hds_path: Path = None,
            model=None,
            vor: vgp = None,
            obs_path: Path = None
    ):
        """
        Class to do stuff with a MODFLOW heads file. This class subclasses the flopy.utils.binaryfile.HeadFile class.
        :param hds_path: path to the heads file
        :param vor: voronoi grid representing model grid for the heads file, optional
        """

        from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
        if model is None:
            self.model = None
        elif isinstance(model, SimulationBase):
            self.model = model
        else:
            raise ValueError("model must be an instance of SimulationBase")
        if hds_path is None:
            if self.model is None:
                raise ValueError("Must provide heads file or model")
            self.hds_path = self.model.model_output_folder_path / f'{self.model.name}.hds'
        else:
            self.hds_path = hds_path

        super().__init__(filename=self.hds_path)

        self.hds = bf.HeadFile(filename=self.hds_path)
        self.kstpkper = self.get_kstpkper()
        self.vor = self.model.vor if vor is None else vor
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
                # df_heads.loc[idxx[kstpkper, layer], :] = spHds.iloc[:, layer].values
                assert self.vor.ncpl == spHds.shape[0], (
                    'Are you using the wrong voronoi grid??? \n'
                    f'The provided vor grid has {self.vor.ncpl} cells, but there are {spHds.shape[0]} heads '
                    f'in the model'
                )
                df_heads.loc[idxx[kstpkper, layer, :]] = spHds.iloc[:, layer].values.reshape(-1, 1)
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
        #  pivot so obs cells are now columns
        obs_heads = obs_reset_idx.pivot(
            index=['layer', 'kstpkper'],
            columns='cell',
            values='elev'
        )
        new_cols = []
        # change cell nums to obs names based on obs dict
        for col in obs_heads.columns:
            obs_name = next(k for k, v in self.obs.items() if v == col)
            new_cols.append(obs_name)
        obs_heads.columns = new_cols
        obs_heads.sort_index(axis=1, inplace=True)

        return obs_heads

    def plot_heads(
            self,
            locs: Path | int | list,
            crs: str = "EPSG:2927",
            layer: int = 0,
            loc_name_field='ExploName'
    ):
        """Plots heads for specified locations in the model.
            Locations should be specified as the model cell/node to 
            plot

            Args:
                locs (Path, optional): Path to shapefile of points to plots heads. If int (or list of ints) is provided, it
                corresponds to a cell idx Defaults to None.
            Returns:
                Returns the fig object and plots it
            """
        if locs is not None:
            fig = figs.Fig()
            heads = self.all_heads
            if isinstance(locs, Path):
                obs_dict = self.vor.get_vor_cells_as_dict(
                    locs=locs,
                    crs=crs,
                    predicate='contains',
                    loc_name_field=loc_name_field
                )
                obs_df = pd.DataFrame.from_dict(obs_dict).transpose()
                obs_locs = obs_df.index
            elif isinstance(locs, int):
                obs_locs = [locs]
            elif isinstance(locs, list):
                obs_locs = locs
            for obs_loc in obs_locs:
                if isinstance(locs, Path):
                    obs_heads = heads.loc[idxx[:, layer, obs_df.loc[obs_loc]], 'elev']
                elif isinstance(locs, int | list):
                    obs_heads = heads.loc[idxx[:, layer, obs_loc], 'elev']
                fig.add_scattergl(
                    x=list(range(len(self.kstpkper))),
                    y=obs_heads,
                    name=obs_loc
                )
            return fig.show()

        return None

    def plot_choropleth(self, *args, **kwargs):
        """Plot heads for a specified time step and stress period on a
            choropleth map. Heads may show saturated thickness (mounding) or
            show elevation head.

            :param obs_name:
            :param obs:
            :param layer:
            :param zmax:
            :param zmin:
            :param bottom:
            :param custom_hover:
            :param all_layers:
            :param bottom_array:
            :param stp_per_to_plot: (tuple, optional) Tuple defining time step and stress period to plot. Defaults to (0,0).
            :param plot_mounding: (boolean, optional) Boolean to determine whether to plot elevation head or mounding (sat thickness)
            :param zoom: define zoom level of plot. Default is 18.
            """
        fig = self.choropleth(*args, **kwargs)
        fig.show()

    def choropleth(
            self,
            kstpkper: tuple = (0, 0),
            plot_mounding: bool = False,
            zmin=None,
            zmax=None,
            zoom=13,
            custom_hover: dict = None,
            bottom=None,
            bottom_array=None,
            all_layers: bool = False,
            layer: int = 1,
            obs: Path = None,
            obs_name: str = 'ExploName'
    ):
        """Instantiate choropleth figure for a specified time step and stress period on a
            choropleth map. Heads may show saturated thickness (mounding) or
            show elevation head.

            :param obs_name: name field in the obs file, defaults to 'ExploName'
            :param obs: file containing locations to plot as scatter points on the map
            :param layer: the layer to plot
            :param zmax: defines the maximum z on the colorscale
            :param zmin: defines the minimum z on the colorscale
            :param zoom: define zoom level of plot. Default is 13.
            :param bottom:
            :param custom_hover: a dict where the keys are the names/labels what will be displayed on hover
                                 and the values for each key is a list of data for all model cells, corresponding
                                 to that label.
            :param all_layers:
            :param bottom_array:
            :param kstpkper: (tuple, optional) Tuple defining time step and stress period to plot. Defaults to (0,0).
            :param plot_mounding: (boolean, optional) Boolean to determine whether to plot elevation head
                                  or mounding (sat thickness)
            :param zoom: define zoom level of plot. Default is 18.
            """

        stp_to_plot = kstpkper[0]
        per_to_plot = kstpkper[1]
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
            hover_dict[f'Layer {lyr + 1} Heads'] = lyr_heads
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

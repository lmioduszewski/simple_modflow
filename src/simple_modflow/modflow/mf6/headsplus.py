from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
    from .voronoiplus import VoronoiGridPlus as Vor
    from figs import Fig

import pandas as pd
import flopy.utils.binaryfile as bf
from pathlib import Path
import geopandas as gpd
from . import mf2Dplots
import figs
from simple_modflow.modflow.utils.datatypes.datalists import convert_nested_to_int
from simple_modflow.modflow.utils.validators import valid_list_of_cell_ints

idxx = pd.IndexSlice  # for easy index slicing in a MultiIndex DataFrame
crs_latlon = "EPSG:4326"


def multimodel_plot_heads(models: list[SimulationBase], locs: int | list[int] | Path, **kwargs):
    figs = []
    for model in models:
        fig = model.hds.plot_heads(locs=locs, plot_fig=False, return_fig=True, **kwargs)
        for trc in fig.data:
            trc.update(name=f'{model.name}-{trc["name"]}')
        figs.append(fig)

    multifig: Fig = figs[0]

    for i, fig in enumerate(figs[1:]):
        for trace in fig.data:
            multifig.add_trace(trace)

    return multifig.show()


class HeadsPlus(bf.HeadFile):

    def __init__(
            self,
            hds_path: Path = None,
            model=None,
            vor: Vor = None,
            obs_path: Path = None
    ):
        """
        Class to do stuff with a MODFLOW heads file. This class subclasses the flopy.utils.binaryfile.HeadFile class.
        :param hds_path: path to the heads file
        :param vor: voronoi grid representing model grid for the heads file, optional
        """
        from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
        if hds_path is None:
            if model is None:
                raise ValueError("Must provide heads file or model")
            assert isinstance(model, SimulationBase), 'no valid model provided'
            self.hds_path = model.model_output_folder_path / f'{model.name}.hds'
        else:
            self.hds_path = hds_path

        super().__init__(filename=self.hds_path)

        if model is None:
            self.model = None
        elif isinstance(model, SimulationBase):
            self.model = model
        else:
            raise ValueError("model must be an instance of SimulationBase")

        self.hds = bf.HeadFile(filename=self.hds_path)
        self.kstpkper = convert_nested_to_int(self.get_kstpkper())
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
        self.crs = self.vor.crs

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

    def get_obs_cells(self, locs: Path, crs: str = None, loc_name_field='ExploName'):
        """
        method to get cells that contain certain observation locations. Locations
        should be points in a shapefile
        :param locs: Path for shapefile with locations of obs as points
        :param crs: crs of shapefile. Gets from Voronoi grid object if no crs provided
        :param loc_name_field: field name in the shapefile attribute table containing observation names.
        :return: dict where observation names are keys and the lists of cells containing them are the values.
        """
        locs = self.obs_path if locs is None else locs
        crs = self.crs if crs is None else crs
        if locs is None:
            return print('no obs path found')
        obs_dict = self.vor.get_vor_cells_as_dict(
            locs=locs,
            crs=crs,
            predicate='contains',
            loc_name_field=loc_name_field)
        obs_dict = self.sort_dict_by_keys(obs_dict)
        # remove dict entries where the loc was not contained in a cell (outside the grid)
        obs_dict = {key: value for key, value in obs_dict.items() if len(value) > 0}
        for obs, cell_ids in obs_dict.items():
            assert len(cell_ids) == 1, f'more than one cell found for {obs}. Fix to make it one cell'
            obs_dict[obs] = cell_ids[0]
        return obs_dict

    def get_obs_heads(
            self, locs: Path | list[int] = None,
            crs: str = None,
            loc_name_field='ExploName',
            long_format=False
    ):
        """
        Retrieves and processes observation heads data for specific observation locations
        and formats it based on user requirements. This function supports both wide and
        long formats for returned data and allows users to specify locations via a file path
        or a list of cell integers.

        :param locs: Path to observation locations file or a list of integer cell indices.
        :param crs: Coordinate reference system (CRS) to use for spatial data transformations.
        :param loc_name_field: Field name in input data representing observation locations.
        :param long_format: If True, returns data in long format; otherwise, wide format.
        :return: A pandas DataFrame containing formatted observation head data.
        :rtype: pandas.DataFrame
        :raises ValueError: If no observation path is provided or no observation cells are given.
        """
        crs = self.crs if crs is None else crs

        if isinstance(locs, Path):
            self.obs_path = locs

        if isinstance(locs, list):
            if valid_list_of_cell_ints(self.model, locs):
                obs_cells = locs
                self._obs = locs

        elif not self.obs and self.obs_path is not None:
            d = self.get_obs_cells(self.obs_path, crs, loc_name_field)
            self._obs.update(d)
            obs_cells = list(self.obs.values())
        else:
            return ValueError('no obs path found or obs cells provided')

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
            if isinstance(self.obs, dict):
                obs_name = next(k for k, v in self.obs.items() if v == col)
            else:
                obs_name = col
            new_cols.append(obs_name)
        obs_heads.columns = new_cols
        obs_heads.sort_index(axis=1, inplace=True)
        obs_heads = obs_heads[sorted(obs_heads.columns)]

        if long_format:
            obs_heads = obs_heads.melt(ignore_index=False, var_name='locs', value_name='elev'
                                       ).reset_index().set_index(['locs', 'layer', 'kstpkper'])

        return obs_heads

    def plot_heads(
            self,
            locs: Path | int | list,
            crs: str = None,
            layer: int = 0,
            loc_name_field='ExploName',
            plot_fig: bool = True,
            return_fig: bool = False,
            show_dates: bool = False,
            show_times: bool = False,
            start_period: int = 0,
    ):
        """
        Plots head values at specified observation locations over the stress periods or
        dates of a model simulation. The function supports plotting for individual
        locations or multiple locations specified by input parameters. Results can be
        visualized directly or returned for further use.

        :param show_times: if True, will use model times as x-value for each plotted kstpkper
        :param start_period: stress period index to start plotting from. Defaults to 0.
        :param locs: Path to a file with observation locations, an integer representing
            a single location index, or a list of location indices.
        :param crs: Coordinate Reference System (CRS) as a string. Defaults to the
            object's CRS if None.
        :param layer: The specific layer of the model for which the heads should be
            plotted. Defaults to 0.
        :param loc_name_field: The field name in the input locs file that specifies
            location names. Used when locs is provided as a file path.
        :param plot_fig: Boolean flag to indicate whether the generated plot should be
            displayed. Defaults to True.
        :param return_fig: Boolean flag to indicate whether the generated plot object
            should be returned. Defaults to False.
        :param show_dates: Boolean flag to indicate whether to use actual model period
            dates on the x-axis of the plot. If False, stress period indices are used.
            Defaults to False.
        :return: Returns the generated plot object if return_fig is True. Otherwise,
            returns None.
        """
        crs = self.crs if crs is None else crs
        if locs is not None:
            fig = figs.Fig()
            heads = self.all_heads

            # drops stress periods less than start_period if it's greater than 0
            if start_period > 0:
                per_tuples = heads.index.get_level_values('kstpkper')
                mask = [j >= start_period for (i, j) in per_tuples]
                heads = heads[mask]

            if isinstance(locs, Path):
                obs_dict = self.vor.get_vor_cells_as_dict(
                    locs=locs,
                    crs=crs,
                    predicate='contains',
                    loc_name_field=loc_name_field
                )
                # remove dict entries where the loc was not contained in a cell (outside the grid)
                obs_dict = {key: value for key, value in obs_dict.items() if len(value) > 0}
                obs_df = pd.DataFrame.from_dict(obs_dict).transpose()
                obs_locs = obs_df.index
            elif isinstance(locs, int):
                obs_locs = [locs]
            elif isinstance(locs, list):
                obs_locs = locs
            else:
                raise ValueError('locs must be a file path, integer, or list of integers')
            for obs_loc in obs_locs:
                if isinstance(locs, Path):
                    obs_heads = heads.loc[idxx[:, layer, obs_df.loc[obs_loc]], 'elev']
                elif isinstance(locs, int | list):
                    obs_heads = heads.loc[idxx[:, layer, obs_loc], 'elev']
                # plot actual period dates on x-axis if provided with model
                if show_times:
                    periods = heads.index.get_level_values('kstpkper').unique().to_list()
                    xs = [round(self.model.times[per]) for per in periods]
                elif show_dates:
                    xs = self.model.per_dates if self.model.per_dates is not None else list(range(len(self.kstpkper)))
                else:
                    xs = list(range(len(self.kstpkper)))

                fig.add_scattergl(
                    x=xs,
                    y=obs_heads,
                    name=obs_loc
                )

            if plot_fig:
                fig.show()
            if return_fig:
                return fig

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

        fig_mbox.add_choroplethmap(
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
            fig_mbox.add_scattermap(
                lat=obs.geometry.y,
                lon=obs.geometry.x,
                text=obs[obs_name],
                hoverinfo='text',
                marker_color='red'
            )

        return fig_mbox

from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
    from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor

from pandas import IndexSlice as idxx
from figs import Fig, create_hover
import dash
from dash import dcc
from dash import html, Input, Output
import dash_bootstrap_components as dbc
import numpy as np
import shapely as shp
from pathlib import Path
from simple_modflow.modflow.utils.datatypes.readers import read_shp_gpkg
import geopandas as gpd
import json
import pandas as pd
import plotly.graph_objs as go
from simple_modflow.modflow.utils.animations import Animation


class Choro:

    def __init__(
            self,
            model: SimulationBase = None,
            vor: Vor = None,
            kstpkper: tuple = None,
            per: int = None,
            layer: int = 0,
            choro_type: str = 'hds',
            custom_hover: dict = None,
            custom_zs: list = None,
            zmin: float | int = None,
            zmax: float | int = None,
            zoom: int = 13,
            show_layer_elevs: bool = True,
            show_mounding: bool = False,
            hover_heads: bool = True,
            hover_ks: bool = False,
            locs: Path = None,
            rch_scale: float = None,
            bgs: bool = False

    ):
        """
        Class defining the basic choropleth plots generated from a modflow model.
        :param model:
        :param vor:
        :param kstpkper:
        :param per: can just provide stress period. appropriate kstpkper tuple will be determined, will throw an
        error if more than one valid kstpkper in the model output exists with the provided per
        :param layer: what layer to plot, a zero-index. 0 equals layer 1.
        :param choro_type:
        :param custom_hover:
        :param custom_zs:
        :param zmin:
        :param zmax:
        :param zoom:
        :param show_layer_elevs: Default is True. To show elevations of all layers on hover
        :param show_mounding: if True, colorscale will be mounding over given Layer
        :param hover_heads:
        :param hover_ks:
        :param locs: specify the path of a shapefile or geopackage with location points to show on choropleth
        :param rch_scale: value to scale z-values in choropleth. For example to convert units in recharge to another L/T
        :param bgs: if True, will take precedence, and will plot water levels in given Layer relative to top of model
        """

        self._vor = None
        self.model = model
        self.vor = model.vor if model is not None else vor
        if model is not None:
            self.kstpkper = self.model.hds.kstpkper[0] if kstpkper is None else kstpkper
        else:
            self.kstpkper = None
        self._per = None
        self.per = per
        self._layer = layer
        self.choro_type = choro_type
        self._custom_hover = custom_hover
        self._custom_zs = custom_zs
        self._zmin = zmin
        self._zmax = zmax
        self.zoom = zoom
        self.show_layer_elevs = show_layer_elevs
        self.show_mounding = show_mounding
        self.hover_heads = hover_heads
        self.hover_ks = hover_ks
        self.locs = locs
        self.rch_scale = rch_scale
        self.bgs = bgs
        self._show_mounding_above_ground = False

        self.fig = Fig()
        self.vor_list = self.vor.gdf_vorPolys.geometry.to_list()
        self.cell_list = [i for i in range(len(self.vor_list))]
        self.area_list = [cell.area for cell in self.vor_list]
        self.x_list = [cell.centroid.xy[0][0] for cell in self.vor_list]
        self.y_list = [cell.centroid.xy[1][0] for cell in self.vor_list]

        self._hover_dict = self.hover_dict_default.copy()
        self._all_heads = None
        self._all_ks = None
        self._hover_dict = self.hover_dict_default

    @property
    def per(self):
        return self._per

    @per.setter
    def per(self, per):
        if per is not None:
            kstpkper = self.model.kstpkper
            per_idx = [i for i, period in enumerate(list(zip(*kstpkper))[1]) if period == per]
            assert len(per_idx) == 1, f'more than one kstpkper with stress period - {per}. Provide unique kstpkper.'
            per_tpl = kstpkper[per_idx[0]]
            self.kstpkper = per_tpl
        self._per = per

    @property
    def all_heads(self):
        if self._all_heads is None:
            self._all_heads = self.model.hds.all_heads
        return self._all_heads

    @property
    def all_ks(self):
        if self._all_ks is None:
            self._all_ks = self.model.gwf.npf.k.data
        return self._all_ks

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, model):
        # assert isinstance(model, SimulationBase), 'model must be an instance of SimulationBase'
        self._model = model

    @property
    def hover_dict_default(self):
        """The hover_dict default"""
        hover_dict_default = {
            'Cell No.': self.cell_list,
            'Area': self.area_list,
            'x': self.x_list,
            'y': self.y_list,
        }
        return hover_dict_default

    @property
    def vor(self):
        return self._vor

    @vor.setter
    def vor(self, vor):
        # assert isinstance(vor, Vor), 'vor must be an instance of VoronoiGridPlus'
        self._vor = vor

    @property
    def kstpkper(self):
        return self._kstpkper

    @kstpkper.setter
    def kstpkper(self, kstpkper):
        if kstpkper is not None:
            assert isinstance(kstpkper, tuple), 'kstpkper must be an instance of tuple'
            assert kstpkper in self.model.hds.kstpkper, f'kstpkper {kstpkper} invalid, not listed in hds file'
        self._kstpkper = kstpkper

    @property
    def nlay(self):
        return self.model.gwf.modelgrid.nlay

    @property
    def hover_dict(self):

        if self.model is not None:
            if self.choro_type == 'hds' or self.hover_heads is True:
                for lyr in range(self.nlay):
                    lyr_heads = self.all_heads.loc[idxx[self.kstpkper, lyr], 'elev'].to_list()
                    self._hover_dict[f'Layer {lyr + 1} Heads'] = lyr_heads

            if self.choro_type == 'ks' or self.hover_ks is True:
                for lyr in range(self.nlay):
                    lyr_ks = self.all_ks[lyr].tolist()
                    self._hover_dict[f'Layer {lyr + 1} Kh'] = lyr_ks

            if self.choro_type == 'rch':
                self._hover_dict['Recharge'] = self.output_rch_zs

        if self.show_layer_elevs:
            layer_nums = self.vor.gdf_topbtm.columns[2:].to_list()
            if self.model is not None:
                botms = self.model.gwf.modelgrid.botm
                top = self.model.gwf.modelgrid.top
            else:
                botms = self.vor.gdf_topbtm.iloc[:, 2:].to_numpy().reshape(-1, 1).transpose()
                top = self.vor.gdf_topbtm.iloc[:, 1].to_numpy().reshape(-1, 1).transpose()[
                    0]  # TODO why do i have to add [0]
            layer_nums = list(range(len(botms)))
            self._hover_dict.update(
                {f'Top of Model': np.round(top, 2)})
            self._hover_dict.update(
                {
                    f'Layer {lyr + 1} Bottom': np.round(botm, 2) for lyr, botm in enumerate(botms)
                }
            )
        if self.show_mounding:
            if self.layer == -1:
                self._show_mounding_above_ground = True
                self.layer = 0
            if self._show_mounding_above_ground is True:
                layer_bottom = self.model.gwf.modelgrid.top.transpose()
            else:
                layer_bottom = self.model.gwf.modelgrid.botm[self.layer].transpose()
            z_hd = self.all_heads.loc[idxx[self.kstpkper, self.layer], 'elev'].reset_index(drop=True)
            zs = z_hd - layer_bottom
            # remove negative mounding values
            zs = zs.mask(zs < 0, 0)
            self._hover_dict.update(
                {
                    f'Layer {self.layer + 1} Mounding': zs
                }
            )
        else:
            self._hover_dict.update(
                {
                    f'zs': self.zs
                }
            )
        if self._custom_hover:
            for key in self._custom_hover.keys():
                self._hover_dict.update(
                    {
                        f'{key}': self._custom_hover[key]
                    }
                )
        if self.bgs:
            z_hd = self.all_heads.loc[idxx[self.kstpkper, self.layer], 'elev'].reset_index(drop=True)
            model_top = self.model.gwf.modelgrid.top.transpose()
            zs = z_hd - model_top
            self._hover_dict.update(
                {
                    f'Layer {self.layer + 1} Below Ground': zs
                }
            )

        return self._hover_dict

    @property
    def layer(self):
        return self._layer

    @layer.setter
    def layer(self, layer):
        assert int(layer) in list(range(self.model.gwf.modelgrid.nlay)), \
            'layer must be an integer, from 0 up to 1 less than the number of model layers'
        self._layer = layer

    @property
    def custom_zs(self):
        return self._custom_zs

    @custom_zs.setter
    def custom_zs(self, custom_zs):
        if not isinstance(self.custom_zs, list):
            raise ValueError('custom_zs must be an instance of list')
        assert len(custom_zs) == self.vor.ncpl, 'customs zs must be provided for every cell'
        self._custom_zs = custom_zs

    @property
    def zs(self):
        """Defines the z values of the choropleth plot, which will be represented by a varying colorscale"""
        if self.custom_zs is not None:
            return self.custom_zs

        if self.choro_type == 'hds' and self.model is not None:
            if self.show_mounding is True:
                # if layer is specified as -1, show mounding above ground
                if self.layer == -1:
                    self._show_mounding_above_ground = True
                    self.layer = 0
                if self._show_mounding_above_ground is True:
                    layer_bottom = self.model.gwf.modelgrid.top.transpose()
                else:
                    layer_bottom = self.model.gwf.modelgrid.botm[self.layer].transpose()
                z_hd = self.all_heads.loc[idxx[self.kstpkper, self.layer], 'elev'].reset_index(drop=True)
                z_hd.loc[z_hd == 1e+30] = np.nan  # make modflow empty elevations NaN
                zs = z_hd - layer_bottom
                # remove negative mounding values
                zs = zs.mask(zs < 0, 0)
            elif self.bgs is True:
                z_hd = self.all_heads.loc[idxx[self.kstpkper, self.layer], 'elev'].reset_index(drop=True)
                model_top = self.model.gwf.modelgrid.top.transpose()
                zs = z_hd - model_top
            else:
                zs = self.all_heads.loc[idxx[self.kstpkper, self.layer], 'elev'].reset_index(drop=True)
                zs.loc[zs == 1e+30] = np.nan  # make modflow empty elevations NaN

        elif self.choro_type == 'ks' and self.model is not None:
            zs = self.all_ks[self.layer].tolist()

        elif self.choro_type == 'output_rch' and self.model is not None:
            zs = self.output_rch_zs

        elif self.choro_type == 'input_rch' and self.model is not None:
            zs = self.input_rch_zs

        else:
            zs = self.vor.gdf_vorPolys.index.to_list()

        return zs

    @property
    def output_rch_zs(self):
        """gets model output recharge values for the specified kstpkper"""
        df = self.model.bud('rch').df
        df = df[~df.index.duplicated()].reset_index()  # remove duplicates...TODO why are there duplicates sometimes
        areas = pd.Series(self.model.vor.area_list, name='area')
        areas.index.name = 'node'
        m = pd.merge(df, areas, on='node')
        full_idx = range(self.model.vor.ncpl)
        m['rch_per_period'] = (m['q'] / m['area'])  # calc L/T recharge
        m = m[m.loc[:, 'kstpkper'] == self.kstpkper]
        m.loc[:, 'node'] = m['node'] - 1  # adjust for zero-based index
        m = m.set_index('node')
        m = m.reindex(full_idx, fill_value=0)
        zs = m['rch_per_period'].to_list()
        if self.rch_scale is not None:
            zs = (np.array(zs) * self.rch_scale).tolist()
        return zs

    @property
    def input_rch_zs(self):
        """gets model input recharge values for the specified kstpkper"""
        d = pd.DataFrame(self.model.gwf.rch.stress_period_data.data[37])
        d['cell'] = d['cellid'].apply(lambda x: x[1])
        d = d.set_index('cell')
        d = d[~d.index.duplicated()]  # remove duplicates if they exist
        full_idx = range(self.model.vor.ncpl)
        d = d.reindex(full_idx, fill_value=0)  # fill in missing cell indices
        zs = d['recharge'].to_list()
        if self.rch_scale is not None:
            zs = (np.array(zs) * self.rch_scale).tolist()
        return zs

    @property
    def colorscale(self):
        if self.choro_type == 'hds':
            return 'earth'
        elif self.choro_type == 'ks':
            return 'earth'
        elif self.choro_type == 'rch':
            return 'earth'
        else:
            return 'earth'

    @property
    def locs(self):
        return self._locs

    @locs.setter
    def locs(self, locs):
        if locs is not None:
            assert isinstance(locs, Path), 'locs must be a Path object'
            try:
                gdf = gpd.read_file(locs)
                gdf = gdf.to_crs(epsg=4326)  # convert to lat/lon
                locs = gdf
            except ValueError:
                print(f'Unable to read {locs}')
        self._locs = locs

    def update_layout(self):

        # Set up default choropleth map styles
        if self.vor:
            map_center = {"lat": self.vor.grid_centroid.y, "lon": self.vor.grid_centroid.x}
        else:
            map_center = None
        self.fig.update_layout(
            margin={"r": 0, "t": 20, "l": 0, "b": 0},
            map_style="carto-voyager",
            map_zoom=self.zoom,
            map_center=map_center
        )

    def add_choropleth(self):
        """creates a choropleth map based on the provided params and adds to the fig"""
        custom_data, hover_template = create_hover(self.hover_dict)
        self.update_layout()
        self.fig.add_trace(self.get_choropleth())

    def get_choropleth(self):

        custom_data, hover_template = create_hover(self.hover_dict)
        choropleth = go.Choroplethmap(
            geojson=self.vor.latlon,
            featureidkey="id",
            locations=self.vor.gdf_latlon.index.to_list(),
            z=self.zs,
            hovertemplate=hover_template,
            customdata=custom_data,
            colorscale=self.colorscale,
            zmax=self._zmax,
            zmin=self._zmin,
        )
        return choropleth

    def add_locs(self, name_field='ExploName'):
        if self.locs is not None:
            geoms = self.locs.geometry
            for idx, row in self.locs.iterrows():
                geom = row.geometry
                if isinstance(geom, shp.Polygon):
                    coords = geom.exterior.xy
                    mode = 'lines'
                elif isinstance(geom, shp.Point):
                    coords = geom.xy
                    mode = 'markers'
                self.fig.add_scattermap(
                    mode=mode,
                    lat=coords[1].tolist(),
                    lon=coords[0].tolist(),
                    name=row[name_field],
                    marker_color='black',
                    showlegend=False,
                )

    @property
    def choropleth(self):
        self.add_choropleth()
        if self.locs is not None:
            self.add_locs()
        return self.fig

    @property
    def ani(self):
        """get animation frames for a choropleth plot"""

        frames = []
        zmin = 1_000_000
        zmax = 0

        for i, per in enumerate(self.model.kstpkper):
            if i > 5:
                continue
            print(f'reading kstpkper {per}', end='\r')
            self.kstpkper = per
            choropleth = self.get_choropleth()
            frame = go.Frame(data=choropleth, name=str(per), baseframe=str(self.model.kstpkper[0]))
            frames.append(frame)
            frame_zmin = round(choropleth.z.min())
            frame_zmax = round(choropleth.z.max())
            zmin = frame_zmin if frame_zmin < zmin else zmin
            zmax = frame_zmax if frame_zmax > zmax else zmax

        # make zmin and zmax the same for all frames
        for frame in frames:
            frame.data[0]['zmin'] = zmin
            frame.data[0]['zmax'] = zmax


        self.fig = Fig(
            data=frames[0].data,
            frames=frames,
            layout=go.Layout(
                updatemenus=[
                    dict(
                        type="buttons",
                        buttons=[dict(label="Play", method="animate", args=[None])],
                    ),
                ],
                # sliders=sliders,
            ),
        )
        # self.fig.update_layout(updatemenus=Animation(self.model).updatemenus)
        # self.fig.update_layout(sliders=Animation(self.model).sliders)

        return self.fig

    def plot(self):
        fig = self.choropleth
        self.fig.show()

    def dash_selector(self):

        self.add_choropleth()

        app = dash.Dash()
        app.layout = html.Div(
            [
                dbc.Row(
                    dbc.Col(
                        [
                            dcc.Graph(
                                figure=self.fig,
                                className="flex-grow-1",
                                style={"height": "95vh"},
                                id="fig",
                            )
                        ],
                        class_name="h-100 d-flex flex-column",
                        style={"height": "95vh"},
                    ),
                    style={"height": "95vh"},
                ),
                dcc.Store(id="selected"),
                html.Div(
                    id='cell_print'
                )
            ],
            style={"height": "95vh"},
        )

        @app.callback(
            Output(component_id="selected", component_property="data"),
            Output(component_id='cell_print', component_property="children"),
            Input(component_id="fig", component_property="selectedData"),
            prevent_initial_callbacks=True,
        )
        def on_select(selectedData):
            if not selectedData:
                return None, None
            selected_cells = []
            for cell in selectedData["points"]:
                selected_cells.append(int(cell["location"]))

            return selected_cells, str(selected_cells)

        app.run(debug=True, port=8050, jupyter_mode='external', use_reloader=False)

from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.grid.voronoi import VoronoiGridPlus as Vor
    from simple_modflow.modflow.mf6.simulation.base import SimulationBase

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
import rasterio
from rasterio.warp import transform_bounds
from PIL import Image
import base64, mimetypes


class Choro:

    def __init__(
            self,
            model: SimulationBase = None,
            vor: Vor = None,
            kstpkper: tuple = None,
            per: int = None,
            layer: int = 0,
            type: str = 'hds',
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
            bgs: bool = False,
            hillshade_path: Path = None,
            colorscale: str = None,
            logscale: bool = False,
            **kwargs

    ):
        """
        Class defining the basic choropleth plots generated from a modflow model.
        :param model: a modflow simulation object. SimulationBase is the base class for all modflow simulations.
        :param vor: a VoronoiGridPlus object.
        :param kstpkper: tuple of stress period and time step to plot
        :param per: can just provide stress period. appropriate kstpkper tuple will be determined, will throw an
        error if more than one valid kstpkper in the model output exists with the provided per
        :param layer: what layer to plot, a zero-index. 0 equals layer 1.
        :param type: type of choropleth to plot - options are 'hds', 'ks', 'input_rch', 'output_rch'
        :param custom_hover: custom dictionary of hover labels to use for choropleth.
        Must be same length as no. of cells in model.
        :param custom_zs: custom list of z values to use for choropleth. Must be same length as no. of cells in model.
        :param zmin: minimum z value to use for choropleth colorscale.
        :param zmax: maximum z value to use for choropleth colorscale.
        :param zoom: zoom level for choropleth map. Default is 13.
        :param show_layer_elevs: Default is True. To show elevations of all layers on hover
        :param show_mounding: if True, colorscale will be mounding over given Layer
        :param hover_heads: boolean to show heads on hover.
        :param hover_ks: boolean to show Kh on hover.
        :param locs: specify the path of a shapefile or geopackage with location points to show on choropleth
        :param rch_scale: value to scale z-values in choropleth. For example to convert units in recharge to another L/T
        :param bgs: if True, will take precedence, and will plot water levels in given Layer relative to top of model
        :param hillshade_path: path to hillshade raster.
        :param colorscale: colorscale to use for choropleth. Default is 'earth'.
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
        self.type = type
        self._custom_hover = custom_hover
        self._custom_zs = custom_zs
        self._zmin = zmin
        self._locs = None
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
        self._colorscale = None
        self.logscale = logscale
        self.kwargs = kwargs

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

        self.hillshade_path = hillshade_path
        self.colorscale = colorscale

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
            if self.type == 'hds' or self.hover_heads is True:
                for lyr in range(self.nlay):
                    lyr_heads = self.all_heads.loc[idxx[self.kstpkper, lyr], 'elev'].to_list()
                    self._hover_dict[f'Layer {lyr + 1} Heads'] = lyr_heads

            if self.type == 'ks' or self.hover_ks is True:
                for lyr in range(self.nlay):
                    lyr_ks = self.all_ks[lyr].tolist()
                    self._hover_dict[f'Layer {lyr + 1} Kh'] = lyr_ks

            if self.type == 'rch':
                self._hover_dict['Recharge'] = self.output_rch_zs

        if self.show_layer_elevs:
            layer_nums = self.vor.gdf_topbtm.columns[2:].to_list()
            if self.model is not None:
                botms = self.model.gwf.modelgrid.botm
                top = self.model.gwf.modelgrid.top
            else:
                botms = self.vor.gdf_topbtm.iloc[:, 2:].to_numpy().reshape(-1, self.vor.nlay).transpose()
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

        if self.type == 'hds' and self.model is not None:
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

        elif self.type == 'ks' and self.model is not None:
            zs = self.all_ks[self.layer].tolist()

        elif self.type == 'output_rch' and self.model is not None:
            zs = self.output_rch_zs

        elif self.type == 'input_rch' and self.model is not None:
            zs = self.input_rch_zs

        else:
            zs = self.vor.gdf_vorPolys.index.to_list()

        if self.logscale:
            zs = np.log10(zs)

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
        """
        Gets the current colorscale used by an instance. If no colorscale
        has been explicitly set, it defaults to 'earth'.

        :return: The current colorscale.
        :rtype: str
        """
        if self._colorscale is None:
            self._colorscale = 'earth'
        return self._colorscale

    @colorscale.setter
    def colorscale(self, colorscale):
        valid_colorscales = [
            'Blackbody', 'Bluered', 'Blues', 'Cividis', 'Earth', 'Electric',
            'Greens', 'Greys', 'Hot', 'Jet', 'Picnic', 'Portland', 'Rainbow',
            'RdBu', 'Reds', 'Viridis', 'YlGnBu', 'YlOrRd']
        if isinstance(colorscale, str):
            if colorscale.lower() in [scale.lower() for scale in valid_colorscales]:
                self._colorscale = colorscale
            else:
                print(f'colorscale {colorscale} not recognized, using default: "earth".\n'
                      f'colorscale options are: {valid_colorscales}')
                self._colorscale = 'earth'

    @property
    def locs(self):
        return self._locs

    @locs.setter
    def locs(self, locs):
        if locs is not None:
            if isinstance(locs, Path):
                try:
                    gdf = gpd.read_file(locs)
                    gdf = gdf.to_crs(epsg=4326)  # convert to lat/lon
                    locations = gdf
                except ValueError:
                    print(f'Unable to read {locs}')
            elif isinstance(locs, gpd.GeoDataFrame):
                locs = locs.to_crs(epsg=4326)
                locations = locs
            else:
                print('could not determine locs; locs should be path or geodataframe')
                locations = None
            self._locs = locations

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
            **self.kwargs,
        )
        return choropleth

    def add_locs(self, name_field='ExploName'):
        if self.locs is not None:
            for idx, row in self.locs.iterrows():
                try:
                    name = row[name_field]
                except:
                    name = idx
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
                    name=name,
                    marker_color='black',
                    showlegend=False,
                )

    @property
    def choropleth(self):
        self.add_choropleth()
        if self.locs is not None:
            self.add_locs()
        if self.hillshade_path is not None:
            self.add_hillshade(self.hillshade_path)
        return self.fig

    def add_hillshade(
            self,
            tif_path: Path = None,
    ):
        hillshade = ChoroplethHillshadeBackground(tif_path, self.model)
        if hillshade.png_path.exists():
            print('hillshade png already exists, using existing file')
        else:
            hillshade.save_geotiff_as_png()
        hillshade.update_fig_layout(self.fig)

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
        """if self.hillshade_path is not None:
            self.add_hillshade(self.hillshade_path)"""
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


class ChoroplethHillshadeBackground:

    def __init__(
            self,
            tif_path: Path = None,
            model: SimulationBase = None,
            png_path: Path = None,
    ):

        self.tif_path = tif_path
        self.model = model
        if png_path is None:
            if model is not None:
                self.png_path = model.model_output_folder_path / 'hillshade.png'
            else:
                self.png_path = Path('hillshade.png')
        self._bounds = None

    @property
    def bounds(self):
        """
        Retrieves the geographical bounds of the raster file in EPSG:4326
        coordinate reference system. It reads the bounds from
        the source file and transforms them into the desired CRS.

        :return: A tuple containing the bounds in the format
                 (left, bottom, right, top).
        :rtype: tuple[float, float, float, float]
        """
        if self._bounds is None:
            with rasterio.open(self.tif_path) as src:
                left, bottom, right, top = transform_bounds(
                    src.crs, "EPSG:4326",
                    src.bounds.left, src.bounds.bottom,
                    src.bounds.right, src.bounds.top
                )
            self._bounds = (left, bottom, right, top)
        return self._bounds

    def save_geotiff_as_png(self):
        """
           Read a GeoTIFF and save as PNG suitable for Plotly image layer.
           - Handles 1, 2, 3, 4 bands.
           Returns: ((left, bottom, right, top) in EPSG:4326)
           """
        with rasterio.open(self.tif_path) as src:
            data = src.read()  # (bands, rows, cols)

        bands = data.shape[0]

        # Normalize to 0–255 uint8 if needed
        def norm255(a):
            a = a.astype(np.float32)
            amin, amax = np.nanmin(a), np.nanmax(a)
            if not np.isfinite(amin) or not np.isfinite(amax) or amax == amin:
                # fall back to zeros if degenerate
                return np.zeros_like(a, dtype=np.uint8)
            return np.clip(255 * (a - amin) / (amax - amin), 0, 255).astype(np.uint8)

        if bands == 1:
            # (H, W) grayscale
            arr = norm255(data[0])
            Image.fromarray(arr, mode="L").save(self.png_path)
        elif bands == 2:
            # Treat band2 as alpha if that makes sense; otherwise just ignore it
            rgb = np.stack([norm255(data[0])] * 3, axis=-1)  # fake RGB from band1
            alpha = norm255(data[1])
            rgba = np.dstack([rgb, alpha])
            Image.fromarray(rgba, mode="RGBA").save(self.png_path)
        elif bands == 3:
            # RGB
            rgb = np.moveaxis(data[:3], 0, 2)
            if rgb.dtype != np.uint8:
                rgb = np.dstack([norm255(data[0]), norm255(data[1]), norm255(data[2])])
            Image.fromarray(rgb, mode="RGB").save(self.png_path)
        else:
            # 4+ bands -> RGBA using first 4
            r, g, b, a = (data[0], data[1], data[2], data[3])
            rgba = np.dstack([*(norm255(x) if x.dtype != np.uint8 else x for x in (r, g, b, a))])
            Image.fromarray(rgba, mode="RGBA").save(self.png_path)

        return

    def file_to_data_uri(self):
        """
        Converts the file specified by the instance attribute `png_path` to a Data URI format.

        This method reads the content of a file, encodes its bytes in base64, deduces
        its MIME type, and constructs a Data URI string in the format:
        `data:[MIME-type];base64,[base64-string]`. By default, if the MIME type is not
        detected, it assumes `image/png`.

        :return: A string representing the file as a Data URI in base64-encoded format.
        :rtype: str
        """
        mime = mimetypes.guess_type(self.png_path)[0] or "image/png"
        with open(self.png_path, "rb") as f:
            b64 = base64.b64encode(f.read()).decode("ascii")
        return f"data:{mime};base64,{b64}"

    def update_fig_layout(self, fig: Fig, zoom: int = 12):
        """
        Updates the layout of a choropleth map figure to show a png
        hillshade background.

        :param fig: The choropleth map figure to update, which should be an instance
            of Choro.
        :param zoom: The zoom level of the map. Defaults to 12.
        :type zoom: int
        :return: None
        """
        left, bottom, right, top = self.bounds
        coords = [[left, top], [right, top], [right, bottom], [left, bottom]]

        fig.update_layout(
            map=dict(
                style="white-bg",  # minimalist canvas
                layers=[
                    dict(
                        sourcetype="image",
                        source=self.file_to_data_uri(),  # local path or URL
                        coordinates=coords,
                        below="traces",  # rendered beneath all data traces
                        opacity=1.0,
                        visible=True,
                        name='hillshade'
                    )
                ],
                center={"lat": (top + bottom) / 2, "lon": (left + right) / 2},
                zoom=zoom
            ))
        fig.update_traces(marker=dict(opacity=0.6), selector=dict(type="choroplethmap"))

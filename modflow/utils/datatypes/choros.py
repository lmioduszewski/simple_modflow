from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
    from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor

from pandas import IndexSlice as idxx
from figs import Fig, create_hover


class Choro:

    def __init__(
            self,
            model: SimulationBase = None,
            vor: Vor = None,
            kstpkper: tuple = None,
            layer: int = 0,
            choro_type: str = 'hds',
            custom_hover: dict = None,
            zmin: float | int = None,
            zmax: float | int = None,
            zoom: int = 13,
            show_layer_elevs: bool = False,
            show_mounding: bool = False,
            hover_heads: bool = True,
            hover_ks: bool = False,
            locs=None

    ):
        """Class defining the basic choropleth plots generated from a modflow model."""

        self.model = model
        self.vor = model.vor if model is not None else vor
        self.kstpkper = self.model.hds.kstpkper[0] if kstpkper is None else kstpkper
        self._layer = layer
        self.choro_type = choro_type
        self._custom_hover = custom_hover
        self._zmin = zmin
        self._zmax = zmax
        self.zoom = zoom
        self.show_layer_elevs = show_layer_elevs
        self.show_mounding = show_mounding
        self.hover_heads = hover_heads
        self.hover_ks = hover_ks
        self.locs = locs

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
        assert isinstance(kstpkper, tuple), 'kstpkper must be an instance of tuple'
        assert kstpkper in self.model.hds.kstpkper, f'kstpkper {kstpkper} invalid, not listed in hds file'
        self._kstpkper = kstpkper

    @property
    def nlay(self):
        return self.model.gwf.modelgrid.nlay

    @property
    def hover_dict(self):

        if self.choro_type == 'hds' or self.hover_heads is True:
            for lyr in range(self.nlay):
                lyr_heads = self.all_heads.loc[idxx[self.kstpkper, lyr], 'elev'].to_list()
                self._hover_dict[f'Layer {lyr + 1} Heads'] = lyr_heads

        if self.choro_type == 'ks' or self.hover_ks is True:
            for lyr in range(self.nlay):
                lyr_ks = self.all_ks[lyr].tolist()
                self._hover_dict[f'Layer {lyr + 1} Kh'] = lyr_ks

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
    def zs(self):
        if self.choro_type == 'hds':
            zs = self.all_heads.loc[idxx[self.kstpkper, self.layer], 'elev'].to_list()
        elif self.choro_type == 'ks':
            zs = self.all_ks[self.layer].tolist()
        else:
            zs = None
        return zs

    @property
    def colorscale(self):
        if self.choro_type == 'hds':
            return 'earth'
        elif self.choro_type == 'ks':
            return 'earth'

    def update_layout(self):

        # Set up default choropleth map styles
        if self.vor:
            mapbox_center = {"lat": self.vor.grid_centroid.y, "lon": self.vor.grid_centroid.x}
        else:
            mapbox_center = None
        self.fig.update_layout(
            margin={"r": 0, "t": 20, "l": 0, "b": 0},
            mapbox_style="carto-positron",
            mapbox_zoom=self.zoom,
            mapbox_center=mapbox_center
        )

    def add_choropleth(self):
        """creates a choropleth map based on the provided params and adds to the fig"""
        custom_data, hover_template = create_hover(self.hover_dict)
        self.update_layout()
        self.fig.add_choroplethmapbox(
            geojson=self.vor.latslons,
            featureidkey="id",
            locations=self.vor.gdf_latlon.index.to_list(),
            z=self.zs,
            hovertemplate=hover_template,
            customdata=custom_data,
            colorscale=self.colorscale,
            zmax=self._zmax,
            zmin=self._zmin,
        )

    def plot(self):
        self.add_choropleth()
        self.fig.show()

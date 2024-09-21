from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
    from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor

from pandas import IndexSlice as idxx
from figs import Fig, create_hover
from simple_modflow.modflow.utils.datatypes.surface_data import ModelSurface
from shapely.geometry import LineString


class XSection:

    def __init__(
            self,
            model: SimulationBase = None,
            kstpkper: tuple = None,
            layer: int = None,
            cell: int = None,
            xy: str = 'x',
            points: list | LineString = None,
            spacing: int = 10,
            num_points: int = 100,
    ):
        self._model = model
        self._vor = None
        self._kstpkper = kstpkper
        self._layer = layer
        self._cell = cell
        self._xy = xy
        self._points = points
        self.spacing = spacing
        self._num_points = num_points

    @property
    def model(self):
        return self._model

    @property
    def vor(self):
        if self._vor is None:
            self._vor = self.model.vor
        return self._vor

    @property
    def surface(self):
        surface = ModelSurface(model=self.model)
        surface = surface.hds(layer=self.layer, kstpkper=self.kstpkper)
        return surface

    @property
    def kstpkper(self):
        return self._kstpkper

    @kstpkper.setter
    def kstpkper(self, kstpkper):
        assert kstpkper in self.model.hds.kstpkper, f'{kstpkper} is not a valid kstpkper'
        self._kstpkper = kstpkper

    @property
    def layer(self):
        return self._layer

    @layer.setter
    def layer(self, layer):
        assert layer in list(range(self.model.gwf.modelgrid.nlay)), f'{layer} is not a valid layer'
        self._layer = layer

    @property
    def cell(self):
        return self._cell

    @cell.setter
    def cell(self, cell):
        assert cell in self.model.vor.gdf_vorPolys.index.to_list(), f'{cell} is not a valid cell'
        self._cell = cell

    @property
    def xy(self):
        if self._xy is None:
            self._xy = 'x'
        return self._xy

    @xy.setter
    def xy(self, xy):
        assert xy in ['x', 'y'], f'{xy} is not x or y'
        self._xy = xy

    @property
    def x_min_max(self):
        xmin, ymin, xmax, ymax = self.vor.get_domain().bounds
        return xmin, xmax

    @property
    def y_min_max(self):
        xmin, ymin, xmax, ymax = self.vor.get_domain().bounds
        return ymin, ymax

    @property
    def num_points(self):
        return self._num_points

    @num_points.setter
    def num_points(self, num_points):
        assert isinstance(num_points, int), f'{num_points} is not an integer'
        self._num_points = num_points

    @property
    def points(self):
        """Get points for cross-section"""

        if self.cell:
            if self.xy == 'x':
                xmin, xmax = self.x_min_max
                y = self.vor.centroids_y[self.cell]
                xsect = LineString(((xmin, y), (xmax, y)))

            if self.xy == 'y':
                ymin, ymax = self.y_min_max
                x = self.vor.centroids_x[self.cell]
                xsect = LineString(((ymin, x), (ymax, x)))

        length = xsect.length
        spacing = length / (self.num_points - 1)

        # Generate equally spaced points along the LineString
        points = [xsect.interpolate(spacing * i) for i in range(self.num_points)]
        self._points = points

        return self._points

    @property
    def xs(self):
        return

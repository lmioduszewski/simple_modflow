from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
    from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor

from pandas import IndexSlice as idxx
from figs import Fig, create_hover
from simple_modflow.modflow.utils.datatypes.surface_data import ModelSurface
from shapely.geometry import LineString
from simple_modflow.modflow.utils.surfaces import InterpolatedSurface
import figs as f
import plotly.graph_objs as go
import numpy as np


class XSection:

    def __init__(
            self,
            model: SimulationBase = None,
            kstpkper: tuple = None,
            layer: int = 0,
            cells: int | list[int] = None,
            x_or_y: str = 'x',
            spacing: int = 10,
            num_points: int = 100,
            extrapolate_beyond_section_ends: bool = False,
            surf_type: str = 'hds'
    ):
        """
        Use to plot a cross-section of heads through a model. Can be used to create an animation
        of head changes for all stress periods. The cross-section line can be defined by providing
        one cell (the 'cells' parameter) or as two ends by providing two cells to the 'cells'
        parameter. If just one cell is given, 'x_or_y' parameter defines whether the cross-section
        is vertical (along 'y' axis) or horizontal (along 'x' axis).

        Examples:

            Show an animated cross-section of all stress periods:

                XSection(model, cells=[1653, 651, 1241]).ani.show() ...OR...
                XSection(model, cells=69, layer=2, x_or_y='y').ani.show()

            Show just a cross-section of one stress period, no animation:

                XSection(model, cells=[1653, 651, 1241], kstpkper=(9, 50)).show()

        :param model: model (SimulationBase object) instance
        :param kstpkper: defaults to the first model stress period if not provided
        :param layer: defaults to 0
        :param cells: defines cross-section location. Can provide any number of cells
        :param x_or_y: only used if one cell is given, defines whether
        the cross-section is vertical (along 'y' axis) or horizontal (along 'x' axis).
        :param spacing: x distance between points on the plot
        :param num_points: number of points in the cross-section plot
        :param extrapolate_beyond_section_ends:
        :param surf_type: can be hds (default) or lyr (for model layers)

        """
        self._model = model
        self._vor = None
        self._kstpkper = kstpkper
        self._layer = layer
        self._cells = cells
        self._x_or_y = x_or_y
        self.spacing = spacing
        self._num_points = num_points
        self._points = None
        self._extrapolate_beyond_section_ends = extrapolate_beyond_section_ends
        self._xsect_linestring = None
        self._xs_as_length = None
        self._x_min_max = None
        self._y_min_max = None
        self.surf_type = surf_type

    @property
    def model(self):
        return self._model

    @property
    def vor(self):
        if self._vor is None:
            self._vor = self.model.vor
        return self._vor

    """@property
    def surface(self):
        surface = ModelSurface(model=self.model)
        surface = surface.hds(layer=self.layer, kstpkper=self.kstpkper)
        return surface"""

    @property
    def kstpkper(self):
        if self._kstpkper is None:
            self._kstpkper = self.model.kstpkper[0]
        return self._kstpkper

    @kstpkper.setter
    def kstpkper(self, kstpkper):
        assert kstpkper in self.model.kstpkper, f'{kstpkper} is not a valid kstpkper'
        self._kstpkper = kstpkper

    @property
    def layer(self):
        return self._layer

    @layer.setter
    def layer(self, layer):
        assert layer in list(range(self.model.gwf.modelgrid.nlay)), f'{layer} is not a valid layer'
        self._layer = layer

    @property
    def cells(self):
        return self._cells

    @cells.setter
    def cells(self, cells):
        if isinstance(cells, int):
            cells = [cells]
        assert all(isinstance(cell, int) for cell in cells), 'cells must be one or more integers'
        assert cells in self.model.vor.gdf_vorPolys.index.to_list(), f'{cells} is not a valid cell'
        self._cells = cells

    @property
    def x_or_y(self):
        if self._x_or_y is None:
            self._x_or_y = 'x'
        return self._x_or_y

    @x_or_y.setter
    def x_or_y(self, xy):
        assert xy in ['x', 'y'], f'{xy} is not x or y'
        self._x_or_y = xy

    @property
    def x_min_max(self):
        if self._x_min_max is None:
            xmin, ymin, xmax, ymax = self.vor.get_domain().bounds
            self._x_min_max = xmin, xmax
        return self._x_min_max

    @property
    def y_min_max(self):
        if self._y_min_max is None:
            xmin, ymin, xmax, ymax = self.vor.get_domain().bounds
            self._y_min_max = ymin, ymax
        return self._y_min_max

    @property
    def extrapolate_beyond_section_ends(self):
        return self._extrapolate_beyond_section_ends

    @extrapolate_beyond_section_ends.setter
    def extrapolate_beyond_section_ends(self, val):
        assert isinstance(val, bool), f'{val} is not a bool'
        self._extrapolate_beyond_section_ends = val

    @property
    def num_points(self):
        return self._num_points

    @num_points.setter
    def num_points(self, num_points):
        assert isinstance(num_points, int), f'{num_points} is not an integer'
        self._num_points = num_points

    @property
    def xsect_linestring(self):

        if self._xsect_linestring is None:

            # if 'cells' is just one point, use it define a vertical or horizontal section line
            if len(self.cells) == 1:
                cell = self.cells[0]
                if self.x_or_y == 'x':
                    xmin, xmax = self.x_min_max
                    y = self.vor.centroids_y[cell]
                    linestring = LineString(((xmin, y), (xmax, y)))

                elif self.x_or_y == 'y':
                    ymin, ymax = self.y_min_max
                    x = self.vor.centroids_x[cell]
                    linestring = LineString(((x, ymin), (x, ymax)))

                else:
                    raise ValueError(f'x_or_y must be either x or y')

            # if 'cells' > 1, then use the cells to define a cross-section line
            elif len(self.cells) > 1:
                cell_centroids = [[self.vor.centroids_x[cell], self.vor.centroids_y[cell]] for cell in self.cells]
                linestring = LineString(cell_centroids)
                if self.extrapolate_beyond_section_ends:
                    raise NotImplementedError

            else:
                raise ValueError(f'{self.cells} has no cells given to define the cross-section')

            self._xsect_linestring = linestring

        return self._xsect_linestring

    @property
    def points(self):
        """Get points for cross-section from self.xsect_linestring based on number of points"""

        if self._points is None:

            linestring = self.xsect_linestring
            length = linestring.length
            spacing = length / (self.num_points - 1)

            # Generate equally spaced points along the LineString
            points = [linestring.interpolate(spacing * i) for i in range(self.num_points)]
            self._points = points

        return self._points

    @property
    def memfile(self):
        """Returns a rasterio memfile of an interpolated surface
        at a particular stress period and layer for the given model"""
        interp = InterpolatedSurface(
            model=self.model,
            layer=self.layer,
            kstpkper=self.kstpkper,
            surf_type=self.surf_type
        )
        memfile = interp.memfile

        return memfile

    @property
    def xsect(self):
        """Opens a rasterio memfile to get points
        elevations at those points to draw a cross-section"""
        with self.memfile.open() as dataset:
            # Use the sample method to extract the elevation along the profile line
            points = [(point.x, point.y) for point in self.points]
            elevations = list(dataset.sample(points))
            elevations = [e[0] for e in elevations]

        # drop the x-section ends [1:-1] then return
        return points[1:-1], elevations[1:-1]

    @property
    def xs_as_length(self):
        """gets xs for cross-section as length along self.xsect_linestring"""

        if self._xs_as_length is None:

            length = self.xsect_linestring.length
            x_start = 0
            xs = np.linspace(x_start, length, self.num_points)
            xs = xs[1:-1]  # drop cross-section ends
            self._xs_as_length = xs

        return self._xs_as_length

    @property
    def fig(self):
        """returns figure of the cross-section"""

        fig = f.Fig()
        points, elevations = self.xsect

        fig.add_scatter(
            x=self.xs_as_length,
            y=elevations,
        )

        return fig

    @property
    def ani(self):
        """get animation frames for cross-section"""

        frames = []
        y_max = 0
        y_min = 1_000_000

        for per in self.model.kstpkper:

            print(f'reading kstpkper {per}', end='\r')
            self.kstpkper = per
            points, elevations = self.xsect

            if np.max(elevations) > y_max:
                y_max = np.max(elevations)
            if np.min(elevations) < y_min:
                y_min = np.min(elevations)

            # define frame for this stress period and append to the frames list
            frame = go.Frame(data=go.Scatter(
                x=self.xs_as_length,
                y=elevations,
                name=f'{per}'

            ), name=f'{per}')
            frames.append(frame)

        y_max = y_max + ((y_max - y_min) * 0.05)  # add a buffer of 5% of the total y-span to y max

        for frame in frames:
            frame.update(dict1={
                'layout': {
                    'yaxis': {
                        'range': [y_min, y_max]
                    }
                }
            })

        # define figure and update layout to include buttons and slider
        fig = f.Fig(
            data=self.fig.data,
            frames=frames
        )
        fig.update_layout(
            yaxis={
                'range': [y_min, y_max]
            },
            updatemenus=[{
                'type': 'buttons',
                'buttons': [
                    {'args': [None, {'frame': {'duration': 125, 'redraw': False},
                                     'transition': {'duration': 0, 'easing': 'quad-in'},
                                     'fromcurrent': True,
                                     'mode': 'afterall'}],
                     'label': 'Play',
                     'method': 'animate'},
                    {'args': [[None], {'mode': 'immediate',
                                       'frame': {'duration': 0, 'redraw': False}
                                       }],
                     'label': 'Pause',
                     'method': 'animate'}
                ]
            }])

        fig.update_layout(
            sliders=[{
                'active': 0,
                'yanchor': 'top',
                'xanchor': 'left',
                'currentvalue': {
                    'font': {'size': 20},
                    'prefix': 'Frame:',
                    'visible': True,
                    'xanchor': 'right'
                },
                'transition': {'duration': 0,
                               'easing': 'linear'},
                'pad': {'b': 10, 't': 50},
                'len': 0.9,
                'x': 0.1,
                'y': 0,
                'steps': [{
                    'args': [[f'{per}'],
                             {'frame': {'duration': 100, 'redraw': False},
                              'mode': 'immediate',
                              'transition': {
                                  'duration': 0,
                                  'easing': 'linear'
                              }}],
                    'label': f'{per}',
                    'method': 'animate'} for per in self.model.kstpkper]}
            ]
        )

        return fig

    def show(self):
        """shows the figure"""
        self.fig.show()

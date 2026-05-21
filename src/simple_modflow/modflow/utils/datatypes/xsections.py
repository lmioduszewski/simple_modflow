from __future__ import annotations
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.grid.voronoi import VoronoiGridPlus as Vor
    from simple_modflow.modflow.mf6.simulation.base import SimulationBase

from pandas import IndexSlice as idxx
from figs import Fig, create_hover
from simple_modflow.modflow.utils.datatypes.surface_data import ModelSurface
from shapely.geometry import LineString
from simple_modflow.modflow.utils.surfaces import InterpolatedSurface
import figs as f
import plotly.graph_objs as go
import numpy as np
from shapely import line_locate_point
from simple_modflow.modflow.utils.animations import Animation
import shapely as shp
from flopy.mf6 import MFSimulation


class XSection:

    def __init__(
            self,
            model: SimulationBase = None,
            per: int = None,
            kstpkper: tuple = None,
            layer: int | list[int] = 0,
            cells: int | list[int] = None,
            x_or_y: str = 'x',
            spacing: int = 10,
            num_points: int = 100,
            extrapolate_beyond_section_ends: bool = False,
            surf_type: str = 'hds',
            interpolate: bool = False,
            use_rbf: bool = True,
            interpolator: str = None,
            section_name: str = None,
            clip: shp.Polygon = None,
            show_model_top: bool = True,
            show_model_btm: bool = False,
            **kwargs
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
        :param section_name: name of section to show on figure (optional), will default to model name
        :param per: stress period number, O-based index; will take presedence over kstpkper if provided
        :param kstpkper: defaults to the first model stress period if not provided
        :param layer: defaults to 0
        :param cells: defines cross-section location. Can provide any number of cells
        :param x_or_y: only used if one cell is given, defines whether
        the cross-section is vertical (along 'y' axis) or horizontal (along 'x' axis).
        :param spacing: x distance between points on the plot
        :param num_points: number of points in the cross-section plot
        :param extrapolate_beyond_section_ends: not implemented
        :param surf_type: can be hds (default) or lyr (for model layers)
        :param interpolate: whether to interpolate the cross-section
        :param use_rbf: Defaults to True, rbf is an interpolation method, use this if having issues
        :param clip: shapely Polygon object to clip the cross-section to
        :param interpolator: define interpolation method. See InterpolatedSurface class for options.
        :param kwargs: additional keyword arguments to pass to InterpolatedSurface class

        """
        self._model = model
        self._vor = None
        self._kstpkper = self.model.kstpkper[per] if per is not None else kstpkper
        self.interpolator = interpolator
        self._layer = None
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
        self.interpolate = True if self.interpolator is not None else interpolate
        self.use_rbf = use_rbf
        self._all_heads = None
        self.section_name = model.name if section_name is None else section_name
        self._clip = clip
        self._kwargs = kwargs
        self.show_model_top = show_model_top
        self.show_model_btm = show_model_btm
        self._overlapping_cells = None
        self._xs = None
        self._model_top = None

        self.layer = layer

    @property
    def model(self):
        return self._model

    @property
    def all_heads(self):
        if self._all_heads is None:
            self._all_heads = self.model.hds.all_heads
        return self._all_heads

    @property
    def model_top(self):
        if self._model_top is None:
            sim = MFSimulation.load(
                sim_name=self.model.sim.name_file.filename[:-4],
                sim_ws=self.model.model_output_folder_path,
                load_only=['disv']
            )
            top = pd.Series(sim.gwf[0].disv.top.data)
            self._model_top = top
        return self._model_top

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
        if isinstance(layer, int):
            layer = [layer]
        assert isinstance(layer, list), f'{layer} is not a integer or list'
        assert all(lyr in list(range(self.model.gwf.modelgrid.nlay)) for lyr in layer), \
            f'one of {layer} is not a valid layer'
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

    @property
    def overlapping_cells(self):

        if self._overlapping_cells is None:
            ov = self.vor.get_vor_cells_as_series(self.xsect_linestring)[0]
            self._overlapping_cells = ov
        return self._overlapping_cells

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
        lyr = self.layer[0]
        interp = InterpolatedSurface(
            model=self.model,
            layer=lyr,
            kstpkper=self.kstpkper,
            surf_type=self.surf_type,
            use_rbf=self.use_rbf,
            clip=self._clip,
            interpolator=self.interpolator,
            **self._kwargs
        )
        memfile = interp.memfile

        return memfile

    @property
    def xs(self):
        """
        Retrieves or calculates the xs property. This is derived based on the
        intersection of the cross-section linestring and the centroidal points
        of Voronoi polygons for the model's geometry. It sorts and filters the
        values by overlapping cells.

        :return: A Pandas Series containing xs values, which are calculated
            by locating points along the cross-section line string with respect
            to the Voronoi polygon centroids. These are filtered and sorted for
            overlapping cells.
        :rtype: pandas.Series
        """
        if self._xs is None:
            vor = self.model.vor
            linestring = self.xsect_linestring
            xs = line_locate_point(
                linestring, vor.gdf_vorPolys.centroid).loc[self.overlapping_cells].sort_values()
            self._xs = xs
        return self._xs

    @property
    def xsect(self):
        """
        Provides the cross-section data for a given profile line.

        This property calculates and returns the cross-section coordinates and
        their corresponding elevations based on whether interpolation is enabled
        or not. When interpolation is enabled, it interpolates elevations along
        the profile line. Otherwise, it retrieves the elevations of overlapping
        cells in the dataset.

        :return: A tuple containing two lists:
            1. The list of coordinates (x, y) along the cross-section line
            2. Corresponding elevation values for those coordinates
        :rtype: tuple[list[tuple[float, float]], list[float]] or None
        """

        if self.interpolate is True:
            # TODO make work for multiple layers
            with self.memfile.open() as dataset:
                # Use the sample method to extract the elevation along the profile line
                points = [(point.x, point.y) for point in self.points]
                elevations = list(dataset.sample(points))
                elevations = [e[0] for e in elevations]

            # drop the x-section ends [1:-1] then return
            return points[1:-1], elevations[1:-1]

        elif self.interpolate is False:
            # if no interpolation, just get head elevations of each overlapping cell
            xs = self.xs
            ys = [self.all_heads.loc[idxx[self.kstpkper, lyr, xs.index.to_list()]] for lyr in self.layer]
            ys_layers = []
            # get head elevations for each layer for each overlapping cell
            for y_lyr in ys:
                ylist = [y[0] for y in y_lyr.values]
                ys_layers.append(ylist)
            # return distance along xsection line for each cell centroid and head of each cell
            return xs.values, ys_layers
        else:
            return None

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

    def to_frame(
            self,
            include_model_top: bool | None = None,
            include_model_btm: bool | None = None,
    ) -> pd.DataFrame:
        """
        Return the cross-section as a long-form DataFrame for external plotting.
        """
        include_model_top = self.show_model_top if include_model_top is None else include_model_top
        include_model_btm = self.show_model_btm if include_model_btm is None else include_model_btm

        rows = []
        points, elevations = self.xsect

        if self.interpolate is True:
            for distance, elevation in zip(self.xs_as_length, elevations):
                rows.append(
                    {
                        "distance": float(distance),
                        "elevation": float(elevation),
                        "series": self.section_name,
                        "kind": "profile",
                    }
                )
        else:
            for i, lyr in enumerate(self.layer):
                series_name = f"Lyr {lyr} hds - {self.section_name}"
                for distance, elevation in zip(points, elevations[i]):
                    rows.append(
                        {
                            "distance": float(distance),
                            "elevation": float(elevation),
                            "series": series_name,
                            "kind": "profile",
                            "layer": lyr,
                        }
                    )

        if include_model_top:
            xs = self.xs
            ys = self.model_top.loc[self.xs.index.to_list()].to_list()
            for distance, elevation in zip(xs, ys):
                rows.append(
                    {
                        "distance": float(distance),
                        "elevation": float(elevation),
                        "series": "model top",
                        "kind": "model_top",
                    }
                )

        if include_model_btm:
            xs = self.xs
            btm_layers = pd.DataFrame(self.model.gwf.modelgrid.botm.T)
            for lyr in btm_layers.columns:
                if lyr not in self.layer:
                    continue
                ys = btm_layers.loc[xs.index.to_list(), lyr].to_list()
                for distance, elevation in zip(xs, ys):
                    rows.append(
                        {
                            "distance": float(distance),
                            "elevation": float(elevation),
                            "series": f"Lyr {lyr} Btm",
                            "kind": "model_bottom",
                            "layer": lyr,
                        }
                    )

        return pd.DataFrame(rows)

    def plot_mpl(self, **kwargs):
        """
        Plot the cross-section with the figs matplotlib cross-section helper.
        """
        from figs.mpl import plot_cross_section

        data = self.to_frame()
        kwargs.setdefault("title", self.section_name)
        kwargs.setdefault("ylabel", "Elevation (ft)")
        return plot_cross_section(
            data=data,
            x="distance",
            y="elevation",
            series_col="series",
            **kwargs,
        )

    @property
    def fig(self):
        """returns figure of the cross-section"""

        fig = f.Fig()
        points, elevations = self.xsect

        if self.interpolate is True:
            fig.add_scatter(x=self.xs_as_length, y=elevations, name=self.section_name)

        elif self.interpolate is False:
            for i, lyr in enumerate(self.layer):
                fig.add_scatter(x=points, y=elevations[i], name=f'Lyr {lyr} hds - {self.section_name}')

        if self.show_model_top:
            xs = self.xs
            # ys = self.vor.gdf_topbtm.loc[xs.index.to_list(), 0].to_list()
            ys = self.model_top.loc[self.xs.index.to_list()].to_list()
            fig.add_scatter(x=xs, y=ys, mode='lines', name='model top')
        if self.show_model_btm:
            xs = self.xs
            btm_layers = pd.DataFrame(self.model.gwf.modelgrid.botm.T)
            for lyr in btm_layers.columns:
                if lyr in self.layer:
                    ys = btm_layers.loc[xs.index.to_list(), lyr].to_list()
                    fig.add_scatter(
                        x=xs, y=ys, mode='lines',
                        name=f'Lyr {lyr} Btm',
                        line=dict(color='black', width=1, dash='dash')
                    )

        return fig

    @property
    def ani(self):
        """get animation frames for cross-section"""

        frames = []
        y_max = 0
        y_min = 1_000_000

        print(f'reading {self.model.name} data...', end='\n')
        if self.interpolator:
            print(f'using {self.interpolator} interpolation method')

        for per in self.model.kstpkper:

            print(f'reading kstpkper {per}', end='\r')
            try:
                # print(f'reading kstpkper {per}', end='\r')
                self.kstpkper = per
                points, elevations = self.xsect
            except:
                continue

            # if minimum and max y-values for this period are greater than the previous max and min,
            # then update the max and min values for the animation
            per_y_max = np.array([float(np.max(elev)) for elev in elevations]).max()
            y_max = per_y_max if per_y_max > y_max else y_max
            per_y_min = np.array([float(np.min(elev)) for elev in elevations]).min()
            y_min = per_y_min if per_y_min < y_min else y_min

            # define frame for this stress period and append to the frames list
            if self.interpolate is True:
                # TODO make work for multiple layers
                frame = go.Frame(data=[
                    go.Scatter(
                        x=self.xs_as_length,
                        y=elevations,
                        name=f'{per}')
                ],
                    name=f'{per}')
            elif self.interpolate is False:
                frame = go.Frame(data=[], name=f'{per}')
                for i, lyr in enumerate(self.layer):
                    tr = go.Scatter(
                        x=points, y=elevations[i],
                        name=f'Lyr {lyr} hds - {self.section_name}'
                    )
                    frame.data += (tr,)  # added comma so tr is treated as a tuple

            if self.show_model_top:
                ys = self.model_top.loc[self.xs.index.to_list()].to_list()
                xs = self.xs
                # ys = self.vor.gdf_topbtm.loc[xs.index.to_list(), 0].to_list()
                model_top = go.Scatter(x=xs, y=ys, mode='lines', name='model top')
                frame.data += (model_top,)
                y_max = np.max(ys)  # y_max is top of model if show_model_top is True

            if self.show_model_btm:
                xs = self.xs
                btm_layers = pd.DataFrame(self.model.gwf.modelgrid.botm.T)
                for lyr in btm_layers.columns:
                    ys = btm_layers.loc[xs.index.to_list(), lyr].to_list()
                    y_min = np.min(ys) if np.min(ys) < y_min else y_min
                    model_btm = go.Scatter(
                        x=xs, y=ys, mode='lines',
                        name=f'Lyr {lyr} Btm',
                        line=dict(color='black', width=1, dash='dash')
                    )
                    frame.data += (model_btm,)

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
            updatemenus=Animation(self.model).updatemenus)
        fig.update_layout(
            sliders=Animation(self.model).sliders,
            xaxis=dict(uirevision="lock"),
            yaxis=dict(uirevision="lock"),
        )
        return fig

    def show(self):
        """shows the figure"""
        self.fig.show()


class MultiModelXSection:

    def __init__(
            self,
            models: list[SimulationBase],
            section_names: list[str] = None,
            cells: int | list[int] = None,
            per: int = None,
            kstpkper: tuple = None,
            layer: int = 0,
            x_or_y: str = 'x',
            spacing: int = 10,
            num_points: int = 100,
            extrapolate_beyond_section_ends: bool = False,
            surf_type: str = 'hds',
            interpolate: bool = False,
            interpolator: str = None,
            use_rbf: bool = True,
            clip: shp.Polygon = None,
            **kwargs
    ):
        """
        set up a list of XSection objects, one for each model provided in 'models' arg. For each model, the
        remaining args will be applied in creating a list of XSection objects. Can then show the animation figure
        using .show() method.
        :param models: list of SimulationBase objects for which to create XSection objects. Should have the same
        model grid or there will be errors or will return erroneous results.
        :param section_names: list of section names to show on figure (optional), will default to model names
        :param cells: defines cross-section location. Can provide any number of cells
        :param per: stress period number, O-based index; will take presedence over kstpkper if provided
        :param kstpkper: defaults to the first model stress period if not provided
        :param layer: defaults to 0
        :param x_or_y: only used if one cell is given, defines whether
        the cross-section is vertical (along 'y' axis) or horizontal (along 'x' axis).
        :param spacing: x distance between points on the plot
        :param num_points: number of points in the cross-section plot
        :param extrapolate_beyond_section_ends: not yet implemented
        :param surf_type: can be hds (default) or lyr (for model layers)
        :param interpolate: if True, interpolate between cells along cross section line
        :param use_rbf: Defaults to True, rbf is an interpolation method, use this if having issues
        :param interpolator: define interpolation method. See InterpolatedSurface class for options.
        :param kwargs: additional keyword arguments to pass to XSection class

        Example usage: MultiModelXSection(models=[model7a, model7b, model7c], cells=[23444, 15525, 16264]).show()
        """

        self.xsect_class_objs = []
        if section_names is not None:
            assert (isinstance(section_names, list)), 'section_names must be a list'
            assert (len(section_names) == len(models)), 'section_names must have same length as models'
            section_names = [str(name) for name in section_names]
        else:
            section_names = [model.name for model in models]
        for i, model in enumerate(models):
            self.xsect_class_objs.append(
                XSection(
                    model=model,
                    section_name=section_names[i],
                    cells=cells,
                    per=per,
                    kstpkper=kstpkper,
                    layer=layer,
                    x_or_y=x_or_y,
                    spacing=spacing,
                    num_points=num_points,
                    extrapolate_beyond_section_ends=extrapolate_beyond_section_ends,
                    surf_type=surf_type,
                    interpolate=interpolate,
                    interpolator=interpolator,
                    use_rbf=use_rbf,
                    clip=clip,
                    **kwargs
                )
            )

        self._anis = None

    @property
    def anis(self):
        """get a list of animation figures for each XSection object"""

        if self._anis is None:
            anis = [xsect_obj.ani for xsect_obj in self.xsect_class_objs]
        self._anis = anis

        return self._anis

    @property
    def fig(self):
        """create and show the animation figure, including all models provided to MultiModelXSection class"""

        anis = self.anis
        animation_fig = anis[0]

        if len(anis) == 1:
            return animation_fig
        else:
            print(f'\n{len(anis)} models for xsection animation')
            for ani in anis[1:]:
                animation_fig.add_traces(ani.data)
                for i, frame in enumerate(animation_fig.frames):
                    frame.data += tuple(ani.frames[i].data)

            return animation_fig

    def show(self):
        self.fig.show()


if __name__ == '__main__':
    import pickle
    from pathlib import Path

    model_path_v7b_et = Path(r"C:\Users\lukem\mf6\cum7bET\cum7bET.model")
    with open(model_path_v7b_et, 'rb') as file:
        model7b: SimulationBase = pickle.load(file)

    fig = XSection(
        model=model7b,
        cells=[16264, 15525, 23444],
        section_name='dev',
        interpolator='cloughTocher2D',
        resolution=300
    ).ani.show()

from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase

import numpy as np
import plotly.graph_objs as go
import shapely as shp
from flopy.mf6 import MFSimulation
from pandas import IndexSlice as idxx
from shapely import line_locate_point
from shapely.geometry import LineString

from myflopy import viz as f
from myflopy.modflow.mf6.grid.interpolated_surface import InterpolatedSurface
from myflopy.modflow.utils.animations import Animation
from myflopy.viz import Fig


def _normalize_section_line(line) -> LineString:
    """Coerce a section-line spec to a shapely ``LineString``.

    Accepts a ``LineString``, a sequence of ``(x, y)`` pairs, or a flopy-style
    ``{"line": [...]}`` dict, so every ``xs(line=...)`` entry point shares one
    normalization.
    """

    if isinstance(line, LineString):
        return line
    if isinstance(line, dict) and "line" in line:
        return _normalize_section_line(line["line"])
    try:
        coords = [(float(x), float(y)) for x, y in line]
    except (TypeError, ValueError) as err:
        raise TypeError(
            "line must be a LineString, a sequence of (x, y) pairs, or a "
            f"{{'line': [...]}} dict; got {line!r}."
        ) from err
    if len(coords) < 2:
        raise ValueError("line must contain at least two (x, y) points.")
    return LineString(coords)


def combined_section_frame(sections: dict[str, XSection]) -> pd.DataFrame:
    """Concatenate section profile tables into one long overlay frame.

    Each section contributes its head-profile series (labeled by its
    ``section_name``) at its *current* ``kstpkper``; the model top is included
    once, from the first section. This is the shared data step behind
    :func:`render_xsections` and the ``kind="xs"`` composers.
    """

    if not sections:
        raise ValueError("at least one section is required.")
    frames = [
        section.to_frame(include_model_top=(index == 0), include_model_btm=False)
        for index, section in enumerate(sections.values())
    ]
    return pd.concat(frames, ignore_index=True)


def render_xsections(
        sections: dict[str, XSection],
        *,
        backend: str = "plotly",
        title: str | None = None,
):
    """Render one or more :class:`XSection` objects as a single overlay figure.

    The shared renderer behind the ``xs`` grammar verb on model, group, and
    diff heads surfaces: each section contributes its head-profile series
    (labeled by its ``section_name``); the model top is drawn once from the
    first section. ``backend="plotly"`` returns a ``viz.Fig``; ``backend="mpl"``
    a matplotlib figure.
    """

    data = combined_section_frame(sections)

    normalized = str(backend).lower()
    if normalized in ("mpl", "matplotlib", "static"):
        from myflopy.viz import plot_cross_section

        result = plot_cross_section(
            data=data,
            x="distance",
            y="elevation",
            series_col="series",
            title=title or "Cross section",
            ylabel="Elevation",
        )
        # figs returns (fig, ax); the grammar contract is a bare figure
        return result[0] if isinstance(result, tuple) else result
    if normalized in ("plotly", "interactive"):
        fig = Fig()
        for series_name, sub in data.groupby("series", sort=False):
            fig.add_scatter(
                x=sub["distance"].to_numpy(),
                y=sub["elevation"].to_numpy(),
                mode="lines",
                name=str(series_name),
            )
        fig.update_layout(
            title=title or "Cross section",
            xaxis_title="Distance",
            yaxis_title="Elevation",
        )
        return fig
    raise ValueError(f"backend must be 'plotly' or 'mpl', got {backend!r}.")


class XSection:

    def __init__(
            self,
            model: SimulationBase = None,
            per: int = None,
            kstpkper: tuple = None,
            layer: int | list[int] = 0,
            cells: int | list[int] = None,
            line=None,
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
            animation_kstpkpers=None,
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
        :param line: defines the cross-section location directly as a shapely
        LineString, a sequence of (x, y) pairs, or a flopy-style {'line': [...]}
        dict. Takes precedence over 'cells' when provided.
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
        self._xsect_linestring = None if line is None else _normalize_section_line(line)
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
        self.animation_kstpkpers = list(model.kstpkper if animation_kstpkpers is None else animation_kstpkpers)
        self._overlapping_cells = None
        self._xs = None
        self._model_top = None

        self.layer = layer

    @property
    def model(self):
        """The MODFLOW simulation this cross-section is taken from."""

        return self._model

    @property
    def all_heads(self):
        """The model's full dependent-variable table (cached on first access).

        Named ``all_heads`` for history, but it is whatever field the model
        actually has: heads on GWF, concentration on GWT, temperature on GWE.
        Reading ``model.hds`` here made ``model.conc.xs()`` and
        ``model.temp.xs()`` raise "is a GWT model; '.hds' is only available on
        GWF models" -- so the ``xs`` verb was broken on the transport readers
        from the day §6.1/6.2 shipped them, while the docs advertised the full
        grammar (fixed 2026-07-27, ledger 99).

        The column NAME differs per kind (``elev``/``conc``/``temp``) but never
        matters: the one consumer reads values positionally.
        """

        if self._all_heads is None:
            # `_field_reader` is the ungated kind-neutral reader; fall back to
            # `.hds` for duck-typed stand-ins that are not a SimulationBase.
            reader = getattr(self.model, "_field_reader", None)
            self._all_heads = (
                reader.all_values if reader is not None else self.model.hds.all_heads
            )
        return self._all_heads

    @property
    def model_top(self):
        """Model-top elevation per cell, read from the DISV grid (cached).

        Loaded from the written simulation's ``disv`` package the first time it is
        needed, so it reflects the model as run.
        """

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
        """The model's Voronoi grid (cached)."""

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
        """The ``(timestep, stress_period)`` sampled (defaults to the model's first)."""

        if self._kstpkper is None:
            self._kstpkper = self.model.kstpkper[0]
        return self._kstpkper

    @kstpkper.setter
    def kstpkper(self, kstpkper):
        """Set the sampled timestep, asserting it exists in the model output."""

        assert kstpkper in self.model.kstpkper, f'{kstpkper} is not a valid kstpkper'
        self._kstpkper = kstpkper

    @property
    def layer(self):
        """The layer(s) drawn in the section, stored as a list of zero-based indices."""

        return self._layer

    @layer.setter
    def layer(self, layer):
        """Set the section layer(s); a bare int is wrapped to a list and each is validated."""

        if isinstance(layer, int):
            layer = [layer]
        assert isinstance(layer, list), f'{layer} is not a integer or list'
        assert all(lyr in list(range(self.model.gwf.modelgrid.nlay)) for lyr in layer), \
            f'one of {layer} is not a valid layer'
        self._layer = layer

    @property
    def cells(self):
        """The cell(s) that define the section line (a single cell or an ordered path)."""

        return self._cells

    @cells.setter
    def cells(self, cells):
        """Set the defining cells; a bare int is wrapped to a list and validated against the grid."""

        if isinstance(cells, int):
            cells = [cells]
        assert all(isinstance(cell, int) for cell in cells), 'cells must be one or more integers'
        assert cells in self.model.vor.gdf_vorPolys.index.to_list(), f'{cells} is not a valid cell'
        self._cells = cells

    @property
    def x_or_y(self):
        """Orientation of a single-cell section line: ``'x'`` (E-W) or ``'y'`` (N-S); default ``'x'``."""

        if self._x_or_y is None:
            self._x_or_y = 'x'
        return self._x_or_y

    @x_or_y.setter
    def x_or_y(self, xy):
        """Set the single-cell section orientation (must be ``'x'`` or ``'y'``)."""

        assert xy in ['x', 'y'], f'{xy} is not x or y'
        self._x_or_y = xy

    @property
    def x_min_max(self):
        """``(xmin, xmax)`` of the model domain, used to span an E-W section line (cached)."""

        if self._x_min_max is None:
            xmin, ymin, xmax, ymax = self.vor.get_domain().bounds
            self._x_min_max = xmin, xmax
        return self._x_min_max

    @property
    def y_min_max(self):
        """``(ymin, ymax)`` of the model domain, used to span a N-S section line (cached)."""

        if self._y_min_max is None:
            xmin, ymin, xmax, ymax = self.vor.get_domain().bounds
            self._y_min_max = ymin, ymax
        return self._y_min_max

    @property
    def extrapolate_beyond_section_ends(self):
        """Whether the section is extended past its end cells (not yet implemented)."""

        return self._extrapolate_beyond_section_ends

    @extrapolate_beyond_section_ends.setter
    def extrapolate_beyond_section_ends(self, val):
        """Set the beyond-ends extrapolation flag (must be a bool)."""

        assert isinstance(val, bool), f'{val} is not a bool'
        self._extrapolate_beyond_section_ends = val

    @property
    def num_points(self):
        """Number of equally spaced sample points taken along the section line."""

        return self._num_points

    @property
    def overlapping_cells(self):
        """Grid cells the section line passes through, ordered along the line (cached)."""

        if self._overlapping_cells is None:
            ov = self.vor.get_vor_cells_as_series(self.xsect_linestring)[0]
            self._overlapping_cells = ov
        return self._overlapping_cells

    @num_points.setter
    def num_points(self, num_points):
        """Set the number of sample points along the section (must be an integer)."""

        assert isinstance(num_points, int), f'{num_points} is not an integer'
        self._num_points = num_points

    @property
    def xsect_linestring(self):
        """The section's ``LineString`` in model coordinates, built from ``cells`` (cached).

        A single cell yields a full-width horizontal or vertical line through the
        cell centroid (per ``x_or_y``); multiple cells yield a polyline through
        their centroids in order.
        """

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
                    raise ValueError('x_or_y must be either x or y')

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
            for distance, elevation in zip(self.xs_as_length, elevations, strict=False):
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
                for distance, elevation in zip(points, elevations[i], strict=False):
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
            for distance, elevation in zip(xs, ys, strict=False):
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
                for distance, elevation in zip(xs, ys, strict=False):
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
        from myflopy.viz import plot_cross_section

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

        periods = self.animation_kstpkpers
        for per in periods:

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

        # define figure and update layout to include buttons and slider
        fig = f.Fig(
            data=self.fig.data,
            frames=frames
        )
        fig.update_layout(
            yaxis={
                'range': [y_min, y_max]
            },
            updatemenus=Animation(self.model, periods=periods).updatemenus)
        fig.update_layout(
            sliders=Animation(self.model, periods=periods).sliders,
            uirevision="lock",
            xaxis=dict(uirevision="lock"),
            yaxis=dict(uirevision="lock"),
        )
        return fig

    def show(self):
        """shows the figure"""
        self.fig.show()



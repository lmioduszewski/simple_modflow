from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase

import pandas as pd

from myflopy._logging import get_logger
from myflopy.modflow.mf6.package_plotting import (
    _apply_backend,
    refuse_noun_parameters,
    resolve_noun_hover,
)
from myflopy.modflow.utils.datatypes.choros import Choro

logger = get_logger(__name__)


class Inputs:
    def __init__(self, model: SimulationBase):
        """Bind the legacy input-explorer namespace to ``model``."""

        self.model = model

    @property
    def rch(self):
        """Recharge input helper for the bound model."""

        return RchInput(self.model)

    @property
    def uzf(self):
        """UZF input helper for the bound model."""

        return UzfInput(self.model)


class UzfInput:
    def __init__(self, model: SimulationBase):
        """Bind the UZF input helper to ``model``."""

        self.model = model

    def finf(self, per: int | str = 0, multiplier: float = 12 * 30):
        """
        Calculate the modified infiltration rate for UZF (Unsaturated Zone Flow) package
        based on the given period and multiplier. The output is adjusted to include all
        model cells even if they have zero infiltration.

        :param per: The stress period index for which the infiltration rates need
            to be calculated. Defaults to 0.
        :type per: int
        :param multiplier: A scaling factor applied to the infiltration rate.
            Defaults to 12 * 30 - converts ft/day to in/month.
        :type multiplier: float
        :return: A pandas Series containing the adjusted infiltration rates for all
            model cells, indexed by cell ID.
        :rtype: pd.Series
        """
        uzf = self.model.gwf.uzf
        ncpl = self.model.modelgrid.ncpl
        packagedata = pd.DataFrame(uzf.packagedata.get_data()).set_index('ifno')

        if isinstance(per, str):
            if per.lower() == 'all':
                dfs = []
                for k, v in uzf.perioddata.data.items():
                    df = pd.DataFrame(v).loc[:, ['ifno', 'finf']]
                    df.columns = ['ifno', k]
                    df = df.set_index('ifno')
                    dfs.append(df)

                perioddata = pd.concat(dfs, axis=1).reset_index()
            else:
                raise ValueError(f'Invalid period value: {per}. "all" is only accepted for this argument.')

        elif isinstance(per, int):
            perioddata = pd.DataFrame(uzf.perioddata.data[per]).loc[:, ['ifno', 'finf']]

        # create 'cell' column - converts ifno to cellid
        perioddata['cell'] = perioddata.ifno.apply(lambda x: packagedata.loc[x].cellid[1])
        perioddata = perioddata.drop(columns='ifno')

        # create dataframe with just cellid and uzf infiltration data
        # then reindex so all model cells are included even if zero finf
        finf = perioddata.set_index('cell').reindex(list(range(ncpl)), fill_value=0)
        logger.info('scaling UZF infiltration rates by %s', multiplier)
        finf = finf * multiplier

        return finf

    def map(
            self,
            per: int = 0,
            multiplier: float = 12 * 30,
            *,
            # -- colour ----------------------------------------------------------
            zmin: float | None = None,
            zmax: float | None = None,
            colorscale: str | list | tuple | None = None,
            logscale: bool = False,
            # -- contours --------------------------------------------------------
            contours: bool | str = False,
            contour_values=None,
            contour_levels: int | float | list = 10,
            contour_color: str = "black",
            contour_width: float = 1.5,
            contour_name: str | None = None,
            contour_clip: bool = True,
            contour_resolution: int = 150,
            contour_method: str = "linear",
            # -- highlighting ----------------------------------------------------
            select=None,
            select_style: str = "outline",
            select_color: str | None = None,
            # -- overlays and framing --------------------------------------------
            locs=None,
            hillshade_path=None,
            fit_bounds: bool = True,
            bounds_padding: float = 0.05,
            # -- hover -----------------------------------------------------------
            hover=None,
            # -- renderer --------------------------------------------------------
            backend: str = "plotly",
            **trace_kwargs,
    ):
        """UZF infiltration (``finf``) for one period, as a choropleth.

        A RECORD noun, and a legacy one: ``model.inputs.uzf`` predates the
        ``model.packages.<pkg>.inputs.<field>`` grammar and survives because UZF
        is indexed by ``ifno`` rather than by cellid, so its period table has to
        be joined back to the package data before it can be coloured. That join
        is :meth:`finf`, and it is why this signature carries ``multiplier`` on
        top of the shared drawing parameters -- the free verb never converts
        units, because it is handed the values directly. Cells with no UZF object
        are filled with zero rather than left blank, so the picture is of the
        whole grid.

        Every parameter is named rather than swept into ``**kwargs``: PyCharm and
        Pylance read the ``def`` line and never run the module, so a parameter
        that arrives through a tail is one no editor can ever offer (plan 8.8).

        ``per`` and ``multiplier`` stay positional-or-keyword, because that is how
        this method has always been called (``uzf.map(3, 1.0)``); everything added
        below is keyword-only.

        The per-layer arguments the free verb has -- ``kstpkper``,
        ``per_timestep``, ``bgs``, ``hover_layers``, ``hover_surfaces``,
        ``show_layer_elevs``, ``show_mounding`` -- are deliberately absent. Each
        was measured to leave a record noun's figure byte-identical, and
        ``show_mounding`` is worse than inert: it injects a head-derived row into
        an elevation tooltip on every cell. ``values``/``type`` and the legacy
        hover trio raise instead of being forwarded: infiltration is what this
        noun draws, and overriding the numbers under a tooltip that still reports
        ``finf`` produces a figure whose cells and whose labels disagree.

        Parameters
        ----------
        per : int, default 0
            Stress period to read, zero-based -- both the UZF period whose
            ``finf`` is drawn and the period the underlying map reads for its
            tooltip. Unlike :meth:`finf`, ``"all"`` is not accepted here: a
            choropleth colours one number per cell, so pick a period (or build
            several maps and hand them to ``model.plot.animate``).
        multiplier : float, default 360
            Scale every infiltration rate before drawing. The default ``12 * 30``
            converts MODFLOW's ft/day into inches/month, which is the unit
            recharge is usually reviewed in; pass ``1.0`` to see the raw model
            units, or ``365.25`` for ft/year. Applied to every cell, including
            the zero-filled ones, so it never changes which cells are blank.
        zmin, zmax : float, optional
            Fixed color-scale limits. Set both to hold the scale steady across
            frames or panels; ``zmin >= zmax`` raises rather than rendering one flat
            color. With neither, the range comes from the data.
        colorscale : str or list of (float, str), optional
            A Plotly colorscale name, or explicit stops. **Pass stops for a diverging
            scale** -- names round-trip through a plotly-to-matplotlib table that
            maps ``'rdbu'`` to the REVERSED colormap, so a named diverging scale
            renders mirrored between backends (ledger 69/70).
        logscale : bool, default False
            Color on a log scale. Non-positive values are masked.
        contours : bool or str, default False
            Overlay contour lines. ``True`` contours the mapped values; a string
            names a different field to contour instead.
        contour_values : sequence of float, optional
            Contour a supplied array rather than the mapped one.
        contour_levels : int or float or list of float, default 10
            A count of levels, a fixed interval, or explicit level values.
        contour_color : str, default 'black'
            Line colour for the contour trace.
        contour_width : float, default 1.5
            Line width for the contour trace.
        contour_name : str, optional
            Legend name for the contour trace.
        contour_clip : bool, default True
            Clip contours to the active domain instead of the full grid extent.
        contour_resolution : int, default 150
            Grid size used to build the contours, per axis. Higher is smoother and
            slower. **Only acts under ``contour_method="cubic"``**, which is the one
            that interpolates onto a square grid before contouring; the default
            linear method triangulates the cell centres directly, so there is no
            grid for this to size and passing it changes nothing (measured: 338
            contour points at 40, 150 and 400 alike).
        contour_method : {'linear', 'cubic'}, default 'linear'
            How the scattered cell values become contours. ``'linear'``
            triangulates the cell centres and contours the triangulation --
            fast, exact at the centres, and faceted. ``'cubic'`` interpolates onto
            a ``contour_resolution``-square grid first (Clough-Tocher) and contours
            that -- smoother, slower, and able to overshoot between cells.
            ``'tri'``/``'tricontour'`` and ``'clough'``/``'clough_tocher'``/
            ``'cloughtocher'`` are accepted as aliases. Anything else raises naming
            both; ``'nearest'`` in particular was documented here for a while and
            has never been implemented.
        select : sequence of int, ndarray, str, Path or geometry, optional
            Cells to highlight. Cell indices, a boolean mask of length ``ncpl``, the
            name of a registered model region (``"all_streams"``), a path to a
            vector file, or a shapely/GeoPandas geometry to intersect. Highlights
            nothing (with a warning) when the selection is empty.
        select_style : {'outline', 'dim', 'both'}, default 'outline'
            How the highlight is drawn. ``'outline'`` traces the dissolved boundary
            of the selection and leaves every cell at full opacity. ``'dim'`` fades
            the *unselected* cells to 20% instead, which suits a bare grid where the
            selection is the subject, but on a field it costs a measured 6.9x of
            readable contrast everywhere you did not select -- and one box/lasso
            gesture in the browser overwrites it. ``'both'`` draws each.
        select_color : str, optional
            Highlight colour, defaulting to ``viz.PALETTE.highlight``. Worth setting
            when the default red collides with a red-blue diverging colorscale.
        locs : Path or GeoDataFrame, optional
            Point locations to mark -- wells, observations, samples. A path is read
            as a vector file.
        hillshade_path : Path, optional
            A hillshade GeoTIFF to draw beneath the cells for topographic context.
        fit_bounds : bool, default True
            Fit the initial view to the grid extent rather than using ``zoom``.
        bounds_padding : float, default 0.05
            Fractional padding around the fitted bounds.
        hover : HoverSpec, optional
            Replace the sectioned hover outright. See
            :mod:`myflopy.modflow.utils.datatypes.hover`.
        backend : {'plotly', 'mpl'}, default 'plotly'
            Which renderer draws the map. ``'plotly'`` returns the interactive
            ``Choro`` picture -- pan, zoom, hover, a basemap. ``'mpl'`` returns a
            static :class:`matplotlib.figure.Figure` instead, for a report, a
            multi-panel figure of your own, or anywhere a live figure is not wanted.
            Accepts ``'interactive'`` and ``'matplotlib'``/``'static'`` as aliases;
            anything else raises rather than being ignored.

            The switch changes the RENDERER, never the subject: both backends draw
            this same map. Two differences are worth knowing before you rely on
            one. The Matplotlib branch draws in **model coordinates** with no
            basemap, so ``bgs`` and ``zoom`` have nothing to act on there; and a
            NAMED diverging colorscale renders mirrored between the two, because the
            name round-trips through a plotly-to-matplotlib table that maps
            ``'rdbu'`` to the reversed colormap -- pass explicit stops when the
            direction carries meaning (ledger 69/70).

            ``backend='mpl'`` and ``.plot_mpl()`` on the returned picture are the
            same renderer reached two ways. Prefer the parameter: it is the spelling
            the whole grammar shares, so it also works on the nouns
            (``model.hds.map(backend='mpl')``) and on the composers
            (``mosaic``/``animate``), where there is no intermediate picture to call
            a method on.

        Returns
        -------
        Choro or matplotlib.figure.Figure
            A :class:`~myflopy.viz.Picture` under the default backend --
            ``.show()``, ``.save(path)``, ``.html(path)`` -- or a bare
            Matplotlib figure with ``backend="mpl"``.

        See Also
        --------
        finf : the per-cell infiltration series behind the picture.

        Examples
        --------
        >>> model.inputs.uzf.map()                       # in/month, period 0
        >>> model.inputs.uzf.map(per=2)
        >>> model.inputs.uzf.map(3, 1.0)                 # period 3, raw ft/day
        >>> model.inputs.uzf.map(logscale=True, colorscale="Blues")
        >>> model.inputs.uzf.map(contours=True, contour_levels=6, zmax=4)
        >>> model.inputs.uzf.map(backend="mpl").savefig("finf.png")
        """

        refuse_noun_parameters("inputs.uzf", "finf", trace_kwargs)
        hover = resolve_noun_hover(hover, trace_kwargs)
        finf = self.finf(per, multiplier).finf
        choro = self.model.plot.map(
            per=per,
            custom_zs=finf.to_list(),
            zmin=zmin,
            zmax=zmax,
            colorscale=colorscale,
            logscale=logscale,
            contours=contours,
            contour_values=contour_values,
            contour_levels=contour_levels,
            contour_color=contour_color,
            contour_width=contour_width,
            contour_name=contour_name,
            contour_clip=contour_clip,
            contour_resolution=contour_resolution,
            contour_method=contour_method,
            select=select,
            select_style=select_style,
            select_color=select_color,
            locs=locs,
            hillshade_path=hillshade_path,
            fit_bounds=fit_bounds,
            bounds_padding=bounds_padding,
            hover=hover,
            **trace_kwargs,
        )
        return _apply_backend(choro, backend)


class RchInput:

    def __init__(
            self,
            model: SimulationBase,
    ):
        """Bind the RCH input helper to ``model`` (default period 0)."""

        self.model = model
        self._df = None
        self.per = 0  # default stress period

    def df(self, per: int|str = None):
        """The RCH stress-period-data table for ``per`` indexed by cell number."""

        per = self.per if per is None else per
        df = pd.DataFrame(self.model.gwf.rch.stress_period_data.data[per])
        df['cellid'] = df['cellid'].apply(lambda x: x[1])  # extract cell numbers
        df = df.set_index('cellid')
        return df

    def map(
            self,
            per: tuple = None,
            multiplier: float | int = 1
    ):
        """Plot the recharge choropleth for one period (filling zero-recharge cells)."""

        model = self.model
        vor = model.vor

        df = pd.DataFrame(model.gwf.rch.stress_period_data.data[per])
        df['cellid'] = df['cellid'].apply(lambda x: x[1])  # extract cell numbers
        df = df.set_index('cellid')
        c_hov = df.reindex(list(range(vor.ncpl)), fill_value=0)  # insert missing cell nums
        c_hov = c_hov * multiplier
        c_hov = c_hov.to_dict(orient='list')
        return Choro(vor=vor, custom_hover=c_hov, custom_zs=c_hov[list(c_hov.keys())[0]])


class DrnInput:

    def __init__(
            self,
            model: SimulationBase,
    ):
        """Bind the DRN input helper to ``model``."""

        self.model = model

    def map(
            self,
            per: tuple = None,
    ):
        """Highlight the DRN cells on the grid for one stress period."""

        model = self.model
        vor = model.vor

        drn_cells = pd.DataFrame(model.gwf.drn.stress_period_data.data[per])['cellid'].apply(lambda x: x[1]).to_list()
        return vor.plot.map(select=drn_cells)

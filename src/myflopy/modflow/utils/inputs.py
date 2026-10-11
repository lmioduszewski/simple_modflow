from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase

import pandas as pd

from myflopy._logging import get_logger
from myflopy.modflow.mf6.package_explorer_utils import cell_numbering
from myflopy.modflow.mf6.package_plotting import (
    _apply_backend,
    refuse_noun_parameters,
    resolve_noun_hover,
)
from myflopy.modflow.mf6.package_registry import (
    SeveralPackagesError,
    model_packages_of_type,
    one_package,
    package_name_of,
)
from myflopy.modflow.utils.datatypes.choros import Choro

logger = get_logger(__name__)


def _legacy_list_package(model: SimulationBase, package_type: str):
    """The ONE list package (RCH, DRN) a legacy helper here reads, by TYPE (ledger 352).

    These read FloPy's ``model.gwf.<type>``, which answers a package NAMED like the
    type first -- beside a second package of the type it read that one alone, in
    silence. Several are refused, naming the view that reads every package.
    """

    found = model_packages_of_type(model, package_type)
    if len(found) == 1:
        return found[0]
    kind = package_type.upper()
    if not found:
        raise AttributeError(f"Model {getattr(model, 'name', '?')!r} has no {kind} package.")
    names = ", ".join(repr(package_name_of(package)) for package in found)
    raise SeveralPackagesError(
        f"Model {getattr(model, 'name', '?')!r} carries {len(found)} {kind} packages ({names}), "
        f"and this legacy helper reads ONE. model.{package_type}.inputs reads every one, each "
        "row tagged with its package's name (pname)."
    )


def _cell_numbers(model: SimulationBase, cellids: pd.Series):
    """Each FloPy cellid's in-layer cell number, by the one cellid rule (ledger 419).

    Parameters
    ----------
    model : SimulationBase
        The model whose grid numbers the cells.
    cellids : pandas.Series
        FloPy cellids: ``(layer, row, col)`` on DIS, ``(layer, cell)`` on DISV,
        ``(node,)`` on DISU.

    Returns
    -------
    numpy.ndarray
        The cell number the model's grid, heads table and cell polygons use:
        ``row * ncol + col`` on DIS, the cell on DISV, the node on DISU. These
        helpers read ``cellid[1]`` before -- a DIS cellid's ROW.
    """

    return cell_numbering(model).layer_cells(cellids.tolist())[1]


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
        uzf = one_package(self.model, "uzf")  # several: refused, naming model.uzf(package=...)
        # The one cellid rule (ledger 419): a DIS cellid's cell is
        # row * ncol + col, and its second item -- read here before -- is the ROW.
        numbering = cell_numbering(self.model)
        ncpl = numbering.ncpl
        packagedata = pd.DataFrame(uzf.packagedata.get_data()).set_index('ifno')
        cell_of_ifno = pd.Series(
            numbering.layer_cells(packagedata['cellid'].tolist())[1], index=packagedata.index
        )

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
        perioddata['cell'] = perioddata.ifno.map(cell_of_ifno)
        perioddata = perioddata.drop(columns='ifno')

        # create dataframe with just cellid and uzf infiltration data
        # then reindex so all model cells are included even if zero finf
        finf = perioddata.set_index('cell').reindex(list(range(ncpl)), fill_value=0)
        # A run LOADED from its files holds the rates as text (FloPy keeps the
        # column object-typed for time-series names), which no float multiplies.
        finf = finf.apply(pd.to_numeric)
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
            features=None,
            hillshade_path=None,
            fit_bounds: bool = True,
            bounds_padding: float = 0.05,
            extent=None,
            # -- hover -----------------------------------------------------------
            hover=None,
            show_layer_elevs: bool | None = None,
            hover_surfaces: bool | None = None,
            # -- renderer --------------------------------------------------------
            backend: str = "plotly",
            **trace_kwargs,
    ):
        """UZF infiltration (``finf``) for one period, as a choropleth.

        A RECORD noun, and a legacy one: ``model.inputs.uzf`` predates the
        ``model.<pkg>.inputs.<field>`` grammar and survives because UZF
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

            .. generated from myflopy.plot.map -- edit there, then run scripts/derive_docstrings.py
        zmin, zmax : float, optional
            Fixed color-scale limits. Set both to hold the scale steady across
            frames or panels; ``zmin >= zmax`` raises rather than rendering one flat
            color. With neither, the range comes from the data.
        colorscale : str or list of (float, str), optional
            Any Plotly colorscale name -- ``'earth'`` (the default), ``'Viridis'``,
            ``'Turbo'``, ``'RdBu'``, ...; append ``_r`` to reverse one -- or explicit
            ``[[position, color], ...]`` stops. The static backend builds its
            colormap from the same stops, so a name renders the same way on both.
            An unknown name warns and falls back to ``'earth'``.
        logscale : bool, default False
            Color on a log10 scale, for a field that spans orders of magnitude
            (hydraulic conductivity, recharge). The colorbar is labelled in real
            units at whole decades on both backends and the hover keeps the real
            values. ``zmin``/``zmax`` (and a ``zmid``) are then in DATA units --
            ``zmin=1e-4, zmax=10`` -- and a non-positive limit raises. Non-positive
            values draw as gaps.
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
            vector file, or a shapely/GeoPandas geometry, by coverage. Highlights
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
        features : path, GeoDataFrame, geometry, Feature, or list/dict of those, optional
            Custom geometry drawn over the cells -- wells, pond footprints, a
            boundary. A path to any vector file GeoPandas reads (name a GeoPackage
            layer with ``Feature(path, layer=...)``), a ``GeoDataFrame`` or
            ``GeoSeries``, a shapely geometry, or a :class:`myflopy.plot.Feature`
            carrying its own ``color``/``width``/``size``/``fill``/``label``/
            ``hover``. A list draws several layers; a dict names them. Points
            become markers, lines polylines, polygons outlines (or fills). Each
            layer gets a legend entry, a stable house color and a hover listing its
            attribute columns. A source without a CRS is taken to be in the grid's.
            Drawn by both backends and kept by ``mosaic`` and ``animate``.
        hillshade_path : Path, optional
            A hillshade GeoTIFF to draw beneath the cells for topographic context.
        fit_bounds : bool, default True
            Fit the initial view to the grid extent rather than using ``zoom``. The
            fit is for the figure's own size, or plotly's default 700x450 when it
            names none -- a larger display only adds margin around the grid -- and
            ``.save()`` to an image refits for the export's size.
        bounds_padding : float, default 0.05
            Fractional padding around the fitted bounds.
        extent : tuple, points, geometry, GeoDataFrame, area or str, optional
            The window the map OPENS on, on both backends -- the one way to frame a
            static (matplotlib) map, which cannot be zoomed. Four numbers
            ``(xmin, ymin, xmax, ymax)`` in the grid's CRS (shapely's order, NOT
            FloPy's ``(xmin, xmax, ymin, ymax)``); a run of ``(x, y)`` points (their
            bounding box); any shapely geometry -- a polygon, a transect line -- by its
            bounds; a GeoDataFrame, reprojected from its own CRS; or the name of an
            area of interest (``model.aoi``). It sets the VIEW only: no cell is
            cropped and the colour scale still spans the whole field, so pass
            ``zmin``/``zmax`` to tighten it for a close-up. ``bounds_padding`` pads it
            as it pads the grid. A box off the grid raises, printing both boxes and
            both CRSs. In a ``mosaic``, one extent frames every panel and a LIST of
            them frames each panel on its own.
        hover : HoverSpec, optional
            Replace the sectioned hover outright. See
            :mod:`myflopy.modflow.utils.datatypes.hover`.
        hover_surfaces : bool, optional
            Add the model-top / layer-bottom table to the sectioned hover.
        show_layer_elevs : bool, optional
            Add model-top and per-layer-bottom rows to the hover. Defaults to
            whether the grid actually carries layer elevations (``vor.gdf_topbtm``),
            because forcing it on a grid without them raises.
        backend : {'plotly', 'mpl'}, default 'plotly'
            Which renderer draws the map. ``'plotly'`` returns the interactive
            ``Choro`` picture -- pan, zoom, hover, a basemap. ``'mpl'`` returns a
            static :class:`matplotlib.figure.Figure` instead, for a report, a
            multi-panel figure of your own, or anywhere a live figure is not wanted.
            Accepts ``'interactive'`` and ``'matplotlib'``/``'static'`` as aliases;
            anything else raises rather than being ignored.

            The switch changes the RENDERER, never the subject: both backends draw
            this same map, and a ``colorscale`` -- named or as stops -- renders the
            same ramp on each (the static colormap is built from Plotly's own
            stops). One difference is worth knowing before you rely on it: the
            Matplotlib branch draws in **model coordinates** with no basemap, so
            ``bgs`` and ``zoom`` have nothing to act on there, while ``features``,
            the ``select`` outline and a log colorbar all carry across. Pass explicit
            stops for a diverging scale when the direction carries meaning: the
            stops state which end is red, which a bare name cannot (ledger 253).

            ``backend='mpl'`` and ``.plot_mpl()`` on the returned picture are the
            same renderer reached two ways. Prefer the parameter: it is the spelling
            the whole grammar shares, so it also works on the nouns
            (``model.hds.map(backend='mpl')``) and on the composers
            (``mosaic``/``animate``), where there is no intermediate picture to call
            a method on.

            .. end generated

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
            features=features,
            hillshade_path=hillshade_path,
            fit_bounds=fit_bounds,
            bounds_padding=bounds_padding,
            extent=extent,
            hover=hover,
            show_layer_elevs=show_layer_elevs,
            hover_surfaces=hover_surfaces,
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
        df = pd.DataFrame(_legacy_list_package(self.model, "rch").stress_period_data.data[per])
        df['cellid'] = _cell_numbers(self.model, df['cellid'])
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

        df = pd.DataFrame(_legacy_list_package(model, "rch").stress_period_data.data[per])
        df['cellid'] = _cell_numbers(model, df['cellid'])
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

        drn_cells = _cell_numbers(
            model, pd.DataFrame(_legacy_list_package(model, "drn").stress_period_data.data[per])['cellid']
        ).tolist()
        return vor.plot.map(select=drn_cells)

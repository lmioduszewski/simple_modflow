"""Input explorer classes behind model.packages."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.modflow.mf6.package_explorer_utils import (
    _filter_normalized_table,
)
from myflopy.modflow.mf6.package_plotting import (
    FieldMappable,
    SpatialView,
    _apply_backend,
    _as_layer_cell_property,
    build_cell_input_map_payload,
    refuse_noun_parameters,
    resolve_noun_hover,
)
from myflopy.modflow.mf6.package_registry import (
    FieldSpec,
    get_default_package_colorscale,
    get_default_package_value_column,
    get_package_explorer_spec,
    get_package_input_field_names,
    get_package_input_field_spec,
)
from myflopy.modflow.mf6.package_tables import (
    build_cell_package_input_table,
    build_uzf_field_input_table,
    build_uzf_field_input_wide_table,
    summarize_input_table,
)
from myflopy.modflow.utils.datatypes.hover import cell_input_hover


class CellPackageInputsExplorer(FieldMappable):
    """Normalized input explorer for one cell-based MF6 stress-period package.

    Each registry-backed field is a first-class node (``ghb.inputs.cond``) with
    the full unified grammar; the namespace verbs take ``field=`` as sugar --
    ``inputs.map(field="cond")`` is exactly ``inputs.cond.map()``, and
    ``inputs.map()`` draws the package's default field.
    """

    def __init__(self, model: SimulationBase, package_name: str):
        """Bind a cell-based input explorer to ``model`` for one package (name lowercased)."""

        self.model = model
        self.package_name = str(package_name).lower()

    @property
    def _default_field(self) -> str:
        """The registry-declared preferred input field (drives bare ``inputs.map()``)."""

        return get_default_package_value_column(self.package_name)

    def _field_names(self) -> list[str]:
        """The registry-declared input field names for this package."""

        return get_package_input_field_names(self.package_name)

    def __getattr__(self, field_name: str) -> CellPackageInputFieldExplorer:
        """Return a field-specific explorer for registry-backed input fields."""

        field_spec = get_package_input_field_spec(self.package_name, field_name)
        if field_spec is None:
            raise AttributeError(
                f"{type(self).__name__!s} has no input field {field_name!r}"
            )
        return CellPackageInputFieldExplorer(self, field_spec)

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return a normalized input table for this package.

        Parameters
        ----------
        per
            Optional zero-based stress period.
        layer
            Optional zero-based layer or layers to keep.
        cells
            Optional zero-based cell ids to keep.
        """

        frame = build_cell_package_input_table(
            self.model,
            self.package_name,
            per=per,
            layer=layer,
            cells=cells,
        )
        return _filter_normalized_table(frame, per=per, layer=layer, cells=cells)

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of available inputs for this package."""

        frame = self.get()
        value_columns = [
            column
            for column in frame.columns
            if column not in {"model", "package", "per", "layer", "cell"}
        ]
        return summarize_input_table(
            frame, label=f"{self.package_name}.inputs", value_columns=value_columns
        )

    @property
    def default(self) -> CellPackageInputFieldExplorer:
        """Return the registry-defined preferred input field."""

        return getattr(self, get_default_package_value_column(self.package_name))

    # map/plot/xs/mosaic/animate come from FieldMappable: they dispatch to the
    # (default) field node, so ``inputs.map()`` == ``inputs.<default>.map()``.



def _sibling_input_fields(package_name: str, field_name: str) -> tuple[str, ...]:
    """The package's OTHER input fields, for an input map's hover.

    Reading a GHB map means asking "what head, against what conductance", and a
    hover carrying only the coloured field answers half of it. The payload
    already has them -- `build_cell_input_map_payload` aggregates every field in
    the selection -- so this only tells the hover spec to render them.

    Names the payload does not carry are dropped by `Fields.build`, so this is
    safe on the explorers whose payload holds one array.
    """

    return tuple(
        name for name in get_package_input_field_names(package_name)
        if name != field_name
    )


def _input_context_fields(hover: Mapping[str, Sequence[Any]]) -> tuple[str, ...]:
    """Cell context worth showing for THIS map: always the layer, records if they merged.

    Which layer you are looking at is always worth saying -- the commonest reason
    a boundary map comes back empty is that the package has no records in the
    default layer 0. The record count only earns its line when some cell actually
    aggregated more than one record, which is the only thing that explains a
    summed value; on a model with one record per cell it would read "records 1"
    on every cell forever.
    """

    counts = hover.get("Record Count") or ()
    merged = any(int(count) > 1 for count in counts)
    return ("Layer", "Record Count") if merged else ("Layer",)


class CellPackageInputFieldExplorer(SpatialView):
    """Field-specific view over a cell package's normalized input table."""

    def __init__(self, inputs: CellPackageInputsExplorer, field_spec: FieldSpec):
        """Pin a package inputs explorer to one field described by ``field_spec``."""

        self.inputs = inputs
        self.field_spec = field_spec
        self.field_name = field_spec.name

    @property
    def model(self):
        """Return the underlying model."""

        return self.inputs.model

    @property
    def package_name(self) -> str:
        """Return the underlying package name."""

        return self.inputs.package_name

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return rows for this input field, preserving package metadata."""

        frame = self.inputs.get(per=per, layer=layer, cells=cells)
        required_columns = ["model", "package", "per", "layer", "cell", self.field_name]
        if frame.empty:
            return pd.DataFrame(columns=required_columns)
        missing = [column for column in required_columns if column not in frame.columns]
        if missing:
            raise KeyError(
                f"Input field {self.field_name!r} is missing required columns: {missing}"
            )
        return frame.loc[:, required_columns].copy()

    def summary(self) -> pd.DataFrame:
        """Return a compact summary for this input field."""

        return summarize_input_table(
            self.get(),
            label=f"{self.package_name}.inputs.{self.field_name}",
            value_columns=[self.field_name],
        )

    def map(
        self,
        *,
        # -- which records ---------------------------------------------------
        per: int = 0,
        layer: int = 0,
        # -- table to cells --------------------------------------------------
        multiplier: float = 1.0,
        fill_value: float | None = None,
        agg: str | None = None,
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
        """This input field as a choropleth: one colour per cell.

        A RECORD noun. The package's period table is reduced to one value per
        cell (``agg``), cells the package does not touch take ``fill_value``, and
        the result is drawn like any other per-cell array. That reduction is why
        this signature carries ``multiplier``/``fill_value``/``agg`` on top of
        the shared drawing parameters -- the free verb never does it, because it
        is handed the values directly.

        Every parameter is named rather than swept into ``**kwargs``: PyCharm and
        Pylance read the ``def`` line and never run the module, so a parameter
        that arrives through a tail is one no editor can ever offer (plan 8.8).

        The per-layer arguments the free verb has -- ``kstpkper``,
        ``per_timestep``, ``bgs``, ``hover_layers``, ``hover_surfaces``,
        ``show_layer_elevs``, ``show_mounding`` -- are deliberately absent. Each
        was measured to leave a record noun's figure byte-identical, and
        ``show_mounding`` is worse than inert: it injects a head-derived row into
        an elevation tooltip on every cell.

        Parameters
        ----------
        per : int, default 0
            Stress period to read, zero-based. A package whose records do not
            change with time has them all under period 0.
        layer : int, default 0
            Zero-based layer. Records in other layers are not drawn; pass the
            layer the boundary actually sits in, which
            ``.get()`` will show you in its ``layer`` column.
        multiplier : float, default 1.0
            Scale every value before drawing -- unit conversion, or a sign flip
            to draw an outflow as positive. Applied before ``agg``.
        fill_value : float, optional
            Value given to cells this package has no record for. Defaults to the
            field's own declared fill (usually ``0.0``); pass ``float("nan")`` to
            leave untouched cells uncoloured instead of colouring them zero.
        agg : {'sum', 'mean', 'min', 'max', 'first', 'last'}, optional
            How several records landing in ONE cell are reduced to one number.
            Defaults to the field's declared aggregation -- ``sum`` for an
            extensive quantity such as a conductance or a flow, ``mean`` for an
            intensive one such as an elevation or a head. Worth setting
            explicitly when a vector source put many features in one cell.
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
        get : the records behind the picture, as a DataFrame.
        summary : the same records reduced to one row.

        Examples
        --------
        >>> model.packages.drn.inputs.elev.map()
        >>> model.packages.ghb.inputs.cond.map(per=3, agg="sum")
        >>> model.packages.drn.inputs.cond.map(logscale=True, colorscale="Viridis")
        >>> model.packages.ghb.inputs.bhead.map(fill_value=float("nan"))
        >>> model.packages.drn.inputs.elev.map(select="all_streams")
        >>> model.packages.ghb.inputs.cond.map(backend="mpl").savefig("cond.png")
        """

        refuse_noun_parameters(
            f"{self.package_name}.inputs.{self.field_name}",
            self.field_name,
            trace_kwargs,
        )
        hover = resolve_noun_hover(hover, trace_kwargs)
        selected = self.inputs.get(per=per, layer=layer)
        values, cell_hover = build_cell_input_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            value_column=self.field_name,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=self.field_spec.fill_value if fill_value is None else fill_value,
            agg=self.field_spec.agg if agg is None else agg,
        )
        choro = self.model.plot.map(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=cell_hover,
            hover_heads=False,
            hover_ks=False,
            # The noun's OWN spec is the picture's BASE hover, so it goes in
            # `hover_spec=`; `hover=` is the call-site override slot and stays the
            # caller's. Putting the default in the override slot renders the same
            # but leaves `choro.hover_spec` None, which disagrees with `model.hds`
            # and with what this noun did before the 8.8 conversion.
            hover=hover,
            hover_spec=cell_input_hover(
                self.field_name,
                extra_fields=_sibling_input_fields(self.package_name, self.field_name),
                context_fields=_input_context_fields(cell_hover),
            ),
            zmin=zmin,
            zmax=zmax,
            colorscale=(
                colorscale
                or self.field_spec.colorscale
                or get_default_package_colorscale(self.package_name)
            ),
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
            **trace_kwargs,
        )
        return _apply_backend(choro, backend)


class UzfFieldInputsExplorer(SpatialView):
    """Normalized explorer for one UZF perioddata field."""

    def __init__(self, model: SimulationBase, field_name: str):
        """Bind an explorer to one UZF perioddata field (e.g. ``finf``, ``pet``)."""

        self.model = model
        self.field_name = str(field_name)

    #: This explorer is UZF-only by construction, but the hover helper asks every
    #: explorer the same question, so it answers it the same way.
    package_name = "uzf"

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return the normalized UZF field table for the selected rows."""

        frame = build_uzf_field_input_table(
            self.model,
            self.field_name,
            per=per,
            layer=layer,
            cells=cells,
        )
        return _filter_normalized_table(frame, per=per, layer=layer, cells=cells)

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of this UZF field's available input data."""

        frame = self.get()
        return summarize_input_table(
            frame,
            label=f"uzf.inputs.{self.field_name}",
            value_columns=[self.field_name],
        )

    def wide(
        self,
        *,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return a wide DataFrame with one row per UZF record and one column per period."""

        return build_uzf_field_input_wide_table(
            self.model,
            self.field_name,
            layer=layer,
            cells=cells,
        )

    def map(
        self,
        *,
        # -- which records ---------------------------------------------------
        per: int = 0,
        layer: int = 0,
        # -- table to cells --------------------------------------------------
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
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
        """This UZF perioddata field as a choropleth: one colour per cell.

        A RECORD noun. UZF perioddata is one row per UZF cell per stress period
        -- ``finf``, ``pet``, ``extdp``, ``extwc``, ``ha``, ``hroot``,
        ``rootact`` -- so one period is selected (``per``), the rows landing in a
        single grid cell are reduced to one number (``agg``), and cells with no
        UZF object at all take ``fill_value``. That reduction is why this
        signature carries ``multiplier``/``fill_value``/``agg`` on top of the
        shared drawing parameters: the free verb never does it, because it is
        handed the values directly.

        Worth knowing about UZF specifically: a vertical column of UZF objects
        stacked under one cell all share that cell, so ``agg`` is doing real work
        on a multi-layer unsaturated column even when ``layer`` is pinned --
        ``sum`` is right for a flux such as ``finf``, ``mean`` for a property
        such as ``extdp``.

        Every parameter is named rather than swept into ``**kwargs``: PyCharm and
        Pylance read the ``def`` line and never run the module, so a parameter
        that arrives through a tail is one no editor can ever offer (plan 8.8).

        The per-layer arguments the free verb has -- ``kstpkper``,
        ``per_timestep``, ``bgs``, ``hover_layers``, ``hover_surfaces``,
        ``show_layer_elevs``, ``show_mounding`` -- are deliberately absent. Each
        was measured to leave a record noun's figure byte-identical, and
        ``show_mounding`` is worse than inert: it injects a head-derived row into
        an infiltration tooltip on every cell.

        Parameters
        ----------
        per : int, default 0
            Stress period to read, zero-based. UZF perioddata is genuinely
            time-varying -- an infiltration series is the usual reason the
            package exists -- so this is the selector you will reach for most.
            ``.wide()`` shows one column per period if you are not sure which one
            you want.
        layer : int, default 0
            Zero-based layer. UZF objects below this layer are not drawn; layer 0
            holds the land-surface objects that receive ``finf``, which is why it
            is the default. ``.get()`` shows the ``layer`` column if the package
            was built with a deeper unsaturated column.
        multiplier : float, default 1.0
            Scale every value before drawing -- unit conversion is the common
            case, e.g. ``multiplier=365.25`` to read a ft/day infiltration rate
            as ft/year. Applied before ``agg``.
        fill_value : float, default 0.0
            Value given to cells with no UZF object. Zero reads naturally for a
            flux ("no infiltration here"); pass ``float("nan")`` to leave the
            non-UZF part of the grid uncoloured instead, which is the honest
            choice for a property such as ``extdp`` where zero is a real and
            different statement.
        agg : {'sum', 'mean', 'min', 'max', 'first', 'last'}, default 'sum'
            How several UZF records landing in ONE grid cell are reduced to one
            number. ``sum`` is right for an extensive quantity (a rate applied
            over the cell); switch to ``mean`` for an intensive one such as
            ``extdp``, ``extwc`` or ``rootact``, where summing a stacked column
            produces a number with no physical meaning.
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
        get : the records behind the picture, as a DataFrame.
        wide : one row per UZF record, one column per period.
        summary : the same records reduced to one row.

        Examples
        --------
        >>> model.packages.uzf.inputs.finf.map()
        >>> model.packages.uzf.inputs.finf.map(per=5, multiplier=365.25)
        >>> model.packages.uzf.inputs.pet.map(per=3, agg="mean")
        >>> model.packages.uzf.inputs.extdp.map(agg="mean", fill_value=float("nan"))
        >>> model.packages.uzf.inputs.finf.map(logscale=True, colorscale="Blues")
        >>> model.packages.uzf.inputs.pet.map(backend="mpl").savefig("pet.png")
        """

        refuse_noun_parameters(
            f"uzf.inputs.{self.field_name}",
            self.field_name,
            trace_kwargs,
        )
        hover = resolve_noun_hover(hover, trace_kwargs)
        selected = self.get(per=per, layer=layer)
        values, cell_hover = build_cell_input_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            value_column=self.field_name,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        choro = self.model.plot.map(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=cell_hover,
            hover_heads=False,
            hover_ks=False,
            # The noun's OWN spec is the picture's BASE hover, so it goes in
            # `hover_spec=`; `hover=` is the call-site override slot and stays the
            # caller's. Putting the default in the override slot renders the same
            # but leaves `choro.hover_spec` None, which disagrees with `model.hds`
            # and with what this noun did before the 8.8 conversion.
            hover=hover,
            hover_spec=cell_input_hover(
                self.field_name,
                extra_fields=_sibling_input_fields(self.package_name, self.field_name),
                context_fields=_input_context_fields(cell_hover),
            ),
            zmin=zmin,
            zmax=zmax,
            colorscale=(
                colorscale
                or get_default_package_colorscale(f"uzf_{self.field_name}")
                or "earth"
            ),
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
            **trace_kwargs,
        )
        return _apply_backend(choro, backend)


class UzfInputsNamespace(FieldMappable):
    """Namespace for normalized UZF input explorers.

    Each perioddata field (``finf``, ``pet``, ...) is a first-class node with
    the full unified grammar; the namespace verbs take ``field=`` as sugar
    (default field: ``finf``).
    """

    _default_field = "finf"

    def __init__(self, model: SimulationBase):
        """Bind the UZF inputs namespace to ``model``."""

        self.model = model

    def _field_names(self) -> list[str]:
        """The registry-declared UZF perioddata field names."""

        return get_package_input_field_names("uzf")

    @property
    def fields(self) -> pd.DataFrame:
        """Return registry metadata for supported UZF perioddata fields."""

        spec = get_package_explorer_spec("uzf")
        if spec is None:
            return pd.DataFrame(
                columns=["field", "label", "colorscale", "fill_value", "agg"]
            )
        rows = [
            {
                "field": field_spec.name,
                "label": field_spec.label,
                "colorscale": field_spec.colorscale,
                "fill_value": field_spec.fill_value,
                "agg": field_spec.agg,
            }
            for field_spec in spec.inputs.values()
        ]
        return pd.DataFrame(rows)

    def summary(self) -> pd.DataFrame:
        """Return one compact summary row per supported UZF input field."""

        frames = []
        for field_name in self.fields["field"].tolist():
            frames.append(self._field(field_name).summary())
        if not frames:
            return pd.DataFrame(
                columns=[
                    "label",
                    "records",
                    "periods",
                    "layers",
                    "cells",
                    "value_columns",
                ]
            )
        summary = pd.concat(frames, ignore_index=True)
        summary["field_name"] = summary["value_columns"].apply(
            lambda values: values[0] if values else None
        )
        field_metadata = self.fields.rename(columns={"field": "field_name"})
        return summary.merge(field_metadata, on="field_name", how="left")

    def _field(self, field_name: str) -> UzfFieldInputsExplorer:
        """Return one registry-backed UZF perioddata field explorer."""

        field_spec = get_package_input_field_spec("uzf", field_name)
        if field_spec is None:
            raise AttributeError(
                f"{type(self).__name__!s} has no UZF input field {field_name!r}"
            )
        return UzfFieldInputsExplorer(self.model, field_spec.name)

    def __getattr__(self, field_name: str) -> UzfFieldInputsExplorer:
        """Return a registry-backed UZF perioddata field explorer."""

        return self._field(field_name)

    @property
    def default(self) -> UzfFieldInputsExplorer:
        """Return the preferred UZF input field."""

        return self.finf

    # map/plot/xs/mosaic/animate come from FieldMappable (field= sugar).

    @property
    def finf(self) -> UzfFieldInputsExplorer:
        """Return the preferred infiltration-rate explorer."""

        return self._field("finf")

    @property
    def pet(self) -> UzfFieldInputsExplorer:
        """Return the potential evapotranspiration explorer."""

        return self._field("pet")

    @property
    def extdp(self) -> UzfFieldInputsExplorer:
        """Return the ET extinction-depth explorer."""

        return self._field("extdp")

    @property
    def extwc(self) -> UzfFieldInputsExplorer:
        """Return the ET extinction-water-content explorer."""

        return self._field("extwc")

    @property
    def ha(self) -> UzfFieldInputsExplorer:
        """Return the surface-depression-storage-depth explorer."""

        return self._field("ha")

    @property
    def hroot(self) -> UzfFieldInputsExplorer:
        """Return the root-zone-thickness explorer."""

        return self._field("hroot")

    @property
    def rootact(self) -> UzfFieldInputsExplorer:
        """Return the root-activity explorer."""

        return self._field("rootact")


class StaticArrayFieldExplorer(SpatialView):
    """Explorer for static layer/cell arrays such as IC, NPF, and STO fields."""

    def __init__(
        self,
        model: SimulationBase,
        package_name: str,
        field_name: str,
        *,
        label: str | None = None,
        colorscale: str = "Viridis",
    ):
        """Bind a static layer/cell array field (IC/NPF/STO) with a label and colorscale."""

        self.model = model
        self.package_name = str(package_name).lower()
        self.field_name = str(field_name)
        self.label = label or f"{self.package_name}.{self.field_name}"
        self.colorscale = colorscale

    def _array(self) -> np.ndarray:
        """Read this field from its package and broadcast it to a ``(nlay, ncpl)`` array."""

        package = self.model.package(self.package_name)
        data = getattr(package, self.field_name)
        values = getattr(data, "array", None)
        if values is None:
            values = getattr(data, "data", data)
        nlay = int(
            getattr(self.model.gwf.modelgrid, "nlay", getattr(self.model, "nlay", 1))
        )
        ncpl = int(self.model.vor.ncpl)
        return _as_layer_cell_property(values, nlay=nlay, ncpl=ncpl, label=self.label)

    def get(
        self,
        *,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return this static array field as a normalized layer/cell table."""

        arr = self._array()
        rows = []
        for layer_index in range(arr.shape[0]):
            for cell in range(arr.shape[1]):
                rows.append(
                    {
                        "model": self.model.name,
                        "package": self.package_name,
                        "field": self.field_name,
                        "layer": layer_index,
                        "cell": cell,
                        self.field_name: arr[layer_index, cell],
                    }
                )
        frame = pd.DataFrame(rows)
        return _filter_normalized_table(frame, layer=layer, cells=cells)

    def summary(self) -> pd.DataFrame:
        """Return a compact summary for this array field."""

        return summarize_input_table(
            self.get(),
            label=f"{self.package_name}.arrays.{self.field_name}",
            value_columns=[self.field_name],
        )

    def long(
        self,
        *,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.Series:
        """Return this field as a series indexed by ``layer/cell``."""

        frame = self.get(layer=layer, cells=cells)
        if frame.empty:
            empty_index = pd.MultiIndex.from_arrays([[], []], names=["layer", "cell"])
            return pd.Series([], index=empty_index, dtype=float, name=self.field_name)
        series = frame.set_index(["layer", "cell"])[self.field_name].sort_index()
        series.name = self.field_name
        return series

    def stack(self, **kwargs) -> pd.Series:
        """Alias for :meth:`long`."""

        return self.long(**kwargs)

    def wide(
        self,
        *,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Pivot this field to one row per cell and one column per layer."""

        frame = self.get(layer=layer, cells=cells)
        if frame.empty:
            return pd.DataFrame(columns=["cell"])
        wide = frame.pivot_table(
            index="cell", columns="layer", values=self.field_name, aggfunc="first"
        )
        wide.columns = [f"layer_{int(column)}" for column in wide.columns]
        return wide.reset_index()

    def map(
        self,
        *,
        # -- which slice of the array ----------------------------------------
        per: int = 0,
        layer: int = 0,
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
        """One layer of this static array as a choropleth: one colour per cell.

        A STATIC ARRAY noun -- ``npf.k``, ``npf.k33``, ``ic.strt``, ``sto.sy``.
        The package already stores exactly one value per cell per layer, so
        unlike a record noun there is nothing to reduce: no ``multiplier``, no
        ``fill_value``, no ``agg``. Pick a ``layer`` and the array is drawn as it
        stands.

        **There is no time axis here, and ``per`` therefore does nothing.** A
        static array is written once in the model's input file and never varies
        by stress period; ``per`` is accepted only so that a caller sweeping the
        same keyword across a row of mixed nouns -- an input record map beside a
        ``k`` map -- does not have to special-case this one. It is documented
        below and then discarded.

        Every parameter is named rather than swept into ``**kwargs``: PyCharm and
        Pylance read the ``def`` line and never run the module, so a parameter
        that arrives through a tail is one no editor can ever offer (plan 8.8).

        The per-layer arguments the free verb has -- ``kstpkper``,
        ``per_timestep``, ``bgs``, ``hover_layers``, ``hover_surfaces``,
        ``show_layer_elevs``, ``show_mounding`` -- are deliberately absent. They
        describe a simulated per-layer FIELD read from an output file; this
        picture is model input, drawn one layer at a time, and each was measured
        to leave it byte-identical. ``show_mounding`` is worse than inert: it
        injects a head-derived row into the tooltip of every cell, on a picture
        that is not a head.

        Parameters
        ----------
        per : int, default 0
            Accepted and ignored. A static array has no stress-period dimension,
            so every value of ``per`` draws the same picture. Present for
            uniformity with the record nouns (``ghb.inputs.cond.map(per=3)``),
            where it selects; if you want a quantity that varies with time, you
            want a ``.inputs`` or ``.results`` noun instead of an array field.
        layer : int, default 0
            Zero-based model layer to draw. This is the ONLY selector that acts
            here: a static array is ``(nlay, ncpl)``, and a constant supplied for
            the whole model is broadcast across layers, so layer 0 and layer 3
            look identical for a uniform ``k`` and quite different for one built
            from a per-unit list. Use ``.wide()`` to see every layer side by side
            before choosing.
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
        get : this array as a normalized layer/cell table.
        wide : the same values pivoted to one column per layer.
        summary : the array reduced to one row (count, range, layers).

        Examples
        --------
        >>> model.packages.npf.k.map()
        >>> model.packages.npf.k33.map(layer=2, logscale=True)
        >>> model.packages.ic.strt.map(contours=True, contour_levels=20)
        >>> model.packages.sto.sy.map(zmin=0.0, zmax=0.3, colorscale="Blues")
        >>> model.packages.npf.k.map(select="all_streams", select_color="red")
        >>> model.packages.sto.ss.map(backend="mpl").savefig("ss.png")
        """

        refuse_noun_parameters(
            f"{self.package_name}.{self.field_name}",
            self.field_name,
            trace_kwargs,
        )
        hover = resolve_noun_hover(hover, trace_kwargs)
        del per
        selected = self.get(layer=layer)
        values, cell_hover = build_cell_input_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            value_column=self.field_name,
            per=None,
            layer=layer,
            agg="first",
        )
        choro = self.model.plot.map(
            per=0,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=cell_hover,
            hover_heads=False,
            hover_ks=False,
            # The noun's OWN spec is the picture's BASE hover, so it goes in
            # `hover_spec=`; `hover=` is the call-site override slot and stays the
            # caller's. Putting the default in the override slot renders the same
            # but leaves `choro.hover_spec` None, which disagrees with `model.hds`
            # and with what this noun did before the 8.8 conversion.
            hover=hover,
            hover_spec=cell_input_hover(
                self.field_name,
                extra_fields=_sibling_input_fields(self.package_name, self.field_name),
                context_fields=_input_context_fields(cell_hover),
            ),
            zmin=zmin,
            zmax=zmax,
            colorscale=colorscale or self.colorscale,
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
            **trace_kwargs,
        )
        return _apply_backend(choro, backend)


__all__ = [
    "CellPackageInputsExplorer",
    "CellPackageInputFieldExplorer",
    "UzfFieldInputsExplorer",
    "UzfInputsNamespace",
    "StaticArrayFieldExplorer",
]

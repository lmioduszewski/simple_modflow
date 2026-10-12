"""Input explorer classes behind model.packages."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.modflow.mf6.hydraulics import model_transmissivity
from myflopy.modflow.mf6.package_explorer_utils import (
    _filter_normalized_table,
    views_nlay,
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

        One row per record of every package of this type (each tagged ``pname``,
        its FloPy name in upper case): the record's own fields (``cond``,
        ``bhead``, ``boundname``, auxiliaries, ...) plus ``layer``, ``cell``,
        ``package``, ``per`` and ``model``. An array-form RCH/EVT gives one row
        per column per period. The three selectors combine with AND.

        Parameters
        ----------
        per : int, optional
            Zero-based stress period: the records MF6 RUNS in it. A block replaces
            every record before it and a period with no block of its own reuses
            the last block before it, so this reads the last block at or before
            ``per`` -- a GHB written once, in period 0, gives its records at every
            ``per`` -- and none under an EMPTY block, each row stamped ``per``. An
            array-form RCH/EVT takes each array from the last block that set it.
            Reading never changes the package (it used to ADD an empty block
            ``per``, which MF6 then ran as "no records from here on"; ledger
            449). A non-integer raises ``TypeError``, a period outside
            ``0 .. nper - 1`` ``IndexError``. Default ``None``: every block as
            WRITTEN, each row's ``per`` the block it came from (an empty block
            gives no row).
        layer : int or iterable of int, optional
            Zero-based layer(s) to keep. Default ``None``: every layer.
        cells : int or iterable of int, optional
            Zero-based cell id(s) WITHIN a layer (``0`` to ``ncpl - 1``) to keep;
            a bare int works as well as a list. Default ``None``: every cell.
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
            if column not in {"model", "package", "pname", "per", "layer", "cell"}
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
        """Return rows for this input field, preserving package metadata.

        Columns ``model``, ``package``, ``pname`` (where the table carries it),
        ``per``, ``layer``, ``cell`` and this field, one row per record; the
        selectors are those of the package's ``inputs.get`` and combine with AND.

        Parameters
        ----------
        per : int, optional
            Zero-based stress period: the records MF6 runs in it -- the last block
            at or before ``per``, none under an empty block -- as the package's
            ``inputs.get`` reads it (ledger 449). Default ``None``: every
            block as written.
        layer : int or iterable of int, optional
            Zero-based layer(s) to keep. Default ``None``: every layer.
        cells : int or iterable of int, optional
            Zero-based cell id(s) within a layer to keep. Default ``None``: every
            cell.
        """

        frame = self.inputs.get(per=per, layer=layer, cells=cells)
        required_columns = ["model", "package", "per", "layer", "cell", self.field_name]
        if frame.empty:
            return pd.DataFrame(columns=required_columns)
        missing = [column for column in required_columns if column not in frame.columns]
        if missing:
            raise KeyError(
                f"Input field {self.field_name!r} is missing required columns: {missing}"
            )
        # `pname` says which package a row came from, where a type has several
        # (ledger 229); a table built without it simply lacks the column.
        keep = required_columns[:2] + [c for c in ("pname",) if c in frame.columns] + required_columns[2:]
        return frame.loc[:, keep].copy()

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
        ``per_timestep``, ``bgs``, ``hover_layers``, ``show_mounding`` -- are
        deliberately not named. Each was measured to leave a record noun's figure
        byte-identical, and ``show_mounding`` is worse than inert: it injects a
        head-derived row into an elevation tooltip on every cell. Passed anyway,
        each raises ``TypeError`` (ledger 465).
        ``show_layer_elevs``/``hover_surfaces`` are named, and act: they read the
        GRID's layer frame (ledger 258).

        Parameters
        ----------
        per : int, default 0
            Stress period to read, zero-based: the records MF6 runs in it -- the
            last block at or before ``per``, so a package written once, in period
            0, draws the same at every ``per`` -- as the package's ``inputs.get``
            reads it (ledger 449). Under an empty block every cell is
            ``fill_value``.
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
        **trace_kwargs
            The open tail, forwarded to ``model.plot.map``. ``values=``,
            ``type=`` and the legacy ``custom_hover``/``hover_heads``/``hover_ks``
            raise ``TypeError``, and ``hover_spec=`` is read as a legacy spelling
            of ``hover`` (``hover`` wins when both are given). A key naming a
            :func:`myflopy.plot.map` parameter this noun does not name reaches it
            as that parameter, unrefused -- ``zoom``, ``hover_fields``,
            ``animation_kstpkpers``, ``rch_scale`` (which acts on a recharge
            map only). The per-layer ``kstpkper``, ``per_timestep``, ``bgs``,
            ``hover_layers`` and ``show_mounding`` raise ``TypeError`` (ledger
            465; ``show_mounding=True`` had kept every colour and added a
            head-derived row to every tooltip, ledger 441). Anything else
            rides to the ``go.Choroplethmap`` trace (``zmid``, ``colorbar``), where
            Plotly refuses an unknown name with ``ValueError`` when the figure is
            built.

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
        get : the records behind the picture, as a DataFrame.
        summary : the same records reduced to one row.

        Examples
        --------
        >>> model.drn.inputs.elev.map()
        >>> model.ghb.inputs.cond.map(per=3, agg="sum")
        >>> model.drn.inputs.cond.map(logscale=True, colorscale="Viridis")
        >>> model.ghb.inputs.bhead.map(fill_value=float("nan"))
        >>> model.drn.inputs.elev.map(select="all_streams")
        >>> model.ghb.inputs.cond.map(backend="mpl").savefig("cond.png")
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
            show_layer_elevs=show_layer_elevs,
            hover_surfaces=hover_surfaces,
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
            features=features,
            hillshade_path=hillshade_path,
            fit_bounds=fit_bounds,
            bounds_padding=bounds_padding,
            extent=extent,
            **trace_kwargs,
        )
        return _apply_backend(choro, backend)


class UzfFieldInputsExplorer(SpatialView):
    """Normalized explorer for one UZF perioddata field."""

    def __init__(self, model: SimulationBase, field_name: str, *, package: str | None = None):
        """Bind an explorer to one UZF perioddata field (e.g. ``finf``, ``pet``).

        ``package`` is the UZF package's NAME (``None``: the model's only one;
        several raise ``SeveralPackagesError``, ledger 352).
        """

        self.model = model
        self.field_name = str(field_name)
        self.package = package

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
        """Return the normalized UZF field table for the selected rows.

        Columns ``model``, ``package`` (``"uzf"``), ``per``, ``ifno``, ``layer``,
        ``cell`` and this field: one row per UZF object listed in a period block,
        placed at its own cell through ``packagedata``. The selectors combine
        with AND.

        Parameters
        ----------
        per : int, optional
            Zero-based stress period: what UZF HOLDS in it. A UZF PERIOD block sets
            only the objects it lists and leaves the rest as an earlier block set
            them (MF6's ``uzf_rp``, measured through MF6), so each object's row
            comes from the last block at or before ``per`` that lists it, stamped
            ``per``; an empty block changes nothing. Reading never changes the
            package (it used to ADD an empty block ``per``; ledger 449). A
            non-integer raises ``TypeError``, a period outside ``0 .. nper - 1``
            ``IndexError``. Default ``None``: every block as WRITTEN, one row per
            object a block lists.
        layer : int or iterable of int, optional
            Zero-based layer(s) of the UZF object's cell to keep. Default
            ``None``: every layer.
        cells : int or iterable of int, optional
            Zero-based cell id(s) within a layer to keep. Default ``None``: every
            cell.
        """

        frame = build_uzf_field_input_table(
            self.model,
            self.field_name,
            per=per,
            layer=layer,
            cells=cells,
            package=self.package,
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
        """Return a wide DataFrame with one row per UZF record and one column per period.

        Rows are every UZF object of the package (``ifno``, ``layer``, ``cell``)
        that the selectors keep; columns ``per_<n>`` are one per PERIOD block the
        package writes (a period MF6 fills from an earlier block has no column),
        ``NaN`` where a block does not list the object. A selection keeping no
        object returns an empty frame with the three key columns.

        Parameters
        ----------
        layer : int or iterable of int, optional
            Zero-based layer(s) of the object's cell to keep. Default ``None``:
            every layer.
        cells : int or iterable of int, optional
            Zero-based cell id(s) within a layer to keep. Default ``None``: every
            cell.
        """

        return build_uzf_field_input_wide_table(
            self.model,
            self.field_name,
            layer=layer,
            cells=cells,
            package=self.package,
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
        ``per_timestep``, ``bgs``, ``hover_layers``, ``show_mounding`` -- are
        deliberately not named. Each was measured to leave a record noun's figure
        byte-identical, and ``show_mounding`` is worse than inert: it injects a
        head-derived row into an infiltration tooltip on every cell. Passed anyway,
        each raises ``TypeError`` (ledger 465).
        ``show_layer_elevs``/``hover_surfaces`` are named, and act: they read the
        GRID's layer frame (ledger 258).

        Parameters
        ----------
        per : int, default 0
            Stress period to read, zero-based: what UZF holds in it, each object
            from the last block at or before ``per`` that lists it (:meth:`get`;
            ledger 449). UZF perioddata is genuinely time-varying -- an
            infiltration series is the usual reason the package exists -- so this
            is the selector you will reach for most. ``.wide()`` shows one column
            per WRITTEN block if you are not sure which one you want.
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
        **trace_kwargs
            The open tail, forwarded to ``model.plot.map``. ``values=``,
            ``type=`` and the legacy ``custom_hover``/``hover_heads``/``hover_ks``
            raise ``TypeError``, and ``hover_spec=`` is read as a legacy spelling
            of ``hover`` (``hover`` wins when both are given). A key naming a
            :func:`myflopy.plot.map` parameter this noun does not name reaches it
            as that parameter, unrefused -- ``zoom``, ``hover_fields``,
            ``animation_kstpkpers``, ``rch_scale`` (which acts on a recharge
            map only). The per-layer ``kstpkper``, ``per_timestep``, ``bgs``,
            ``hover_layers`` and ``show_mounding`` raise ``TypeError`` (ledger
            465; ``show_mounding=True`` had kept every colour and added a
            head-derived row to every tooltip, ledger 441). Anything else
            rides to the ``go.Choroplethmap`` trace (``zmid``, ``colorbar``), where
            Plotly refuses an unknown name with ``ValueError`` when the figure is
            built.

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
        get : the records behind the picture, as a DataFrame.
        wide : one row per UZF record, one column per period.
        summary : the same records reduced to one row.

        Examples
        --------
        >>> model.uzf.inputs.finf.map()
        >>> model.uzf.inputs.finf.map(per=5, multiplier=365.25)
        >>> model.uzf.inputs.pet.map(per=3, agg="mean")
        >>> model.uzf.inputs.extdp.map(agg="mean", fill_value=float("nan"))
        >>> model.uzf.inputs.finf.map(logscale=True, colorscale="Blues")
        >>> model.uzf.inputs.pet.map(backend="mpl").savefig("pet.png")
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
            show_layer_elevs=show_layer_elevs,
            hover_surfaces=hover_surfaces,
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
            features=features,
            hillshade_path=hillshade_path,
            fit_bounds=fit_bounds,
            bounds_padding=bounds_padding,
            extent=extent,
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

    def __init__(self, model: SimulationBase, *, package: str | None = None):
        """Bind the UZF inputs namespace to ``model``'s UZF package ``package``.

        ``package`` is the package's NAME (``None``: the model's only UZF).
        """

        self.model = model
        self.package = package

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
        return UzfFieldInputsExplorer(self.model, field_spec.name, package=self.package)

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

    #: The free verb's ``type`` for a DERIVED field this noun draws through
    #: ``model.plot.map`` (transmissivity); ``None`` for a stored array.
    _map_type: str | None = None

    def _map_per(self, per):
        """The ``per`` a derived field's map reads (``None``: it has no time)."""

        return None

    def __init__(
        self,
        model: SimulationBase,
        package_name: str,
        field_name: str,
        *,
        label: str | None = None,
        colorscale: str = "Viridis",
        siblings: Sequence[str] = (),
    ):
        """Bind a static layer/cell array field (IC/NPF/STO) with a label and colorscale.

        ``siblings`` are the package's other DECLARED array fields (``k22``/``k33``
        beside ``k``, ``sy`` beside ``ss``). They are read at the SAME layer and
        listed under the mapped value in the hover -- reading a K map means asking
        "and what is K33 here", and a hover carrying only the coloured field
        answers half of it. The same rule the list packages already follow through
        :func:`_sibling_input_fields`, and declared by the namespace that owns the
        field table rather than guessed here.
        """

        self.model = model
        self.package_name = str(package_name).lower()
        self.field_name = str(field_name)
        self.label = label or f"{self.package_name}.{self.field_name}"
        self.colorscale = colorscale
        self.siblings = tuple(str(name) for name in siblings)

    #: What reading a SIBLING array may raise. It is a convenience on someone
    #: else's map, so an unreadable one is dropped rather than taking the figure
    #: down: a missing attribute, a shape `_as_layer_cell_property` refuses, or a
    #: package object that does not answer `.array`/`.data`.
    _SIBLING_UNAVAILABLE = (AttributeError, IndexError, KeyError, TypeError, ValueError)

    def _sibling_arrays(self, layer: int) -> dict[str, list]:
        """``{field: per-cell values}`` for this package's other arrays at ``layer``.

        Only fields with VALUES at this layer. The presence test is the finite
        data, not the attribute: FloPy answers ``npf.k22`` with a live object on
        every model, and its ``.array`` is ``None`` on the many that never set
        it -- which reached the hover as a row of blanks under a real field name.
        A row that is empty on every cell is the same defect as naming a
        parameter that cannot act.
        """

        package = self.model.package(self.package_name)
        ncpl = int(self.model.vor.ncpl)
        arrays: dict[str, list] = {}
        for name in self.siblings:
            if getattr(package, name, None) is None:
                continue
            try:
                values = type(self)(self.model, self.package_name, name)._array()
            except self._SIBLING_UNAVAILABLE:
                continue
            if int(layer) >= values.shape[0]:
                continue
            at_layer = np.asarray(values[int(layer)], dtype=float)[:ncpl]
            if not np.isfinite(at_layer).any():
                continue
            arrays[name] = at_layer.tolist()
        return arrays

    def _array(self) -> np.ndarray:
        """Read this field from its package and broadcast it to a ``(nlay, ncpl)`` array."""

        package = self.model.package(self.package_name)
        data = getattr(package, self.field_name)
        values = getattr(data, "array", None)
        if values is None:
            values = getattr(data, "data", data)
        # The views' layers, not FloPy's: on a DISU FloPy infers two from IHC and
        # the node array was broadcast onto both, inventing a layer 1 (ledger 448).
        nlay = views_nlay(self.model.gwf)
        ncpl = int(self.model.vor.ncpl)
        return _as_layer_cell_property(values, nlay=nlay, ncpl=ncpl, label=self.label)

    def get(
        self,
        *,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return this static array field as a normalized layer/cell table.

        Columns ``model``, ``package``, ``field``, ``layer``, ``cell`` and this
        field, one row per layer and cell of the ``(nlay, ncpl)`` array -- a
        constant or a per-layer list declared for the model is broadcast to every
        cell first. Values are in the package's own units. The selectors combine
        with AND.

        Parameters
        ----------
        layer : int or iterable of int, optional
            Zero-based layer(s) to keep; one the model lacks selects nothing.
            Default ``None``: every layer.
        cells : int or iterable of int, optional
            Zero-based cell id(s) within a layer to keep. Default ``None``: every
            cell.
        """

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
        """Return this field as a series indexed by ``layer/cell``.

        Sorted by layer then cell, named after the field; a selection keeping
        nothing returns an empty float series with the same two index levels.

        Parameters
        ----------
        layer : int or iterable of int, optional
            Zero-based layer(s) to keep, as in :meth:`get`. Default ``None``:
            every layer.
        cells : int or iterable of int, optional
            Zero-based cell id(s) within a layer to keep. Default ``None``: every
            cell.
        """

        frame = self.get(layer=layer, cells=cells)
        if frame.empty:
            empty_index = pd.MultiIndex.from_arrays([[], []], names=["layer", "cell"])
            return pd.Series([], index=empty_index, dtype=float, name=self.field_name)
        series = frame.set_index(["layer", "cell"])[self.field_name].sort_index()
        series.name = self.field_name
        return series

    def stack(self, **kwargs) -> pd.Series:
        """Alias for :meth:`long`.

        Parameters
        ----------
        **kwargs
            Passed to :meth:`long` unchanged: ``layer=`` and ``cells=``, by
            keyword. Anything else -- ``per=`` included -- raises ``TypeError``.
        """

        return self.long(**kwargs)

    def wide(
        self,
        *,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Pivot this field to one row per cell and one column per layer.

        A ``cell`` column, then ``layer_<k>`` for each layer kept. A selection
        keeping nothing returns an empty frame with a ``cell`` column.

        Parameters
        ----------
        layer : int or iterable of int, optional
            Zero-based layer(s) to keep as columns. Default ``None``: every layer.
        cells : int or iterable of int, optional
            Zero-based cell id(s) within a layer to keep as rows. Default
            ``None``: every cell.
        """

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
        """One layer of this static array as a choropleth: one colour per cell.

        A STATIC ARRAY noun -- ``npf.k``, ``npf.k33``, ``ic.strt``, ``sto.sy``.
        The package already stores exactly one value per cell per layer, so
        unlike a record noun there is nothing to reduce: no ``multiplier``, no
        ``fill_value``, no ``agg``. Pick a ``layer`` and the array is drawn as it
        stands.

        **There is no time axis here, and ``per`` therefore does nothing** -- with
        one exception, ``npf.transmissivity``, whose thickness is the SATURATED
        one and so moves with the heads. A static array is written once in the
        model's input file and never varies by stress period; ``per`` is accepted
        only so that a caller sweeping the same keyword across a row of mixed
        nouns -- an input record map beside a ``k`` map -- does not have to
        special-case this one.

        Every parameter is named rather than swept into ``**kwargs``: PyCharm and
        Pylance read the ``def`` line and never run the module, so a parameter
        that arrives through a tail is one no editor can ever offer (plan 8.8).

        The per-layer arguments the free verb has -- ``kstpkper``,
        ``per_timestep``, ``bgs``, ``hover_layers``, ``show_mounding`` -- are
        deliberately not named. They describe a simulated per-layer FIELD read
        from an output file; this picture is model input, drawn one layer at a
        time, and on a stored array (``npf.k``) each was measured to leave it
        byte-identical. On ``npf.transmissivity``, which draws through the free
        map's own transmissivity hover, ``hover_layers`` does reshape the
        tooltip (measured). ``show_mounding`` is worse than inert: it injects a
        head-derived row into the tooltip of every cell, on a picture that is not
        a head. Passed anyway, each raises ``TypeError`` -- bar ``hover_layers`` on
        ``npf.transmissivity``, which acts there and is passed on (ledger
        465).
        ``show_layer_elevs``/``hover_surfaces`` are named, and act: they read the
        GRID's layer frame (ledger 258).

        Parameters
        ----------
        per : int, default 0
            Accepted and ignored -- except by ``npf.transmissivity``, where it
            picks the stress period whose LAST saved heads set the saturated
            thickness (``model.plot.map(type="trans", kstpkper=...)`` reaches any
            other time step). A static array has no stress-period dimension,
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
        **trace_kwargs
            The open tail, forwarded to ``model.plot.map``. ``values=``,
            ``type=`` and the legacy ``custom_hover``/``hover_heads``/``hover_ks``
            raise ``TypeError``, and ``hover_spec=`` is read as a legacy spelling
            of ``hover`` (``hover`` wins when both are given). A key naming a
            :func:`myflopy.plot.map` parameter this noun does not name reaches it
            as that parameter, unrefused -- ``zoom``, ``hover_fields``,
            ``animation_kstpkpers``, ``rch_scale`` (which acts on a recharge
            map only). The per-layer ``kstpkper``, ``per_timestep``, ``bgs``,
            ``hover_layers`` and ``show_mounding`` raise ``TypeError`` -- bar
            ``hover_layers`` on ``npf.transmissivity``, which reshapes its
            tooltip and is passed on (ledger 465; ``show_mounding=True``
            had kept every colour and added a head-derived row to every tooltip,
            ledger 441). Anything else
            rides to the ``go.Choroplethmap`` trace (``zmid``, ``colorbar``), where
            Plotly refuses an unknown name with ``ValueError`` when the figure is
            built.

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
        get : this array as a normalized layer/cell table.
        wide : the same values pivoted to one column per layer.
        summary : the array reduced to one row (count, range, layers).

        Examples
        --------
        >>> model.npf.k.map()
        >>> model.npf.k33.map(layer=2, logscale=True)
        >>> model.ic.strt.map(contours=True, contour_levels=20)
        >>> model.sto.sy.map(zmin=0.0, zmax=0.3, colorscale="Blues")
        >>> model.npf.k.map(select="all_streams", select_color="red")
        >>> model.sto.ss.map(backend="mpl").savefig("ss.png")
        """

        refuse_noun_parameters(
            f"{self.package_name}.{self.field_name}",
            self.field_name,
            trace_kwargs,
            # A derived field (transmissivity) keeps the free verb's per-layer
            # hover, which `hover_layers` reshapes; a plain array has none.
            honours=("hover_layers",) if self._map_type is not None else (),
        )
        hover = resolve_noun_hover(hover, trace_kwargs)
        drawing = dict(
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
            features=features,
            hillshade_path=hillshade_path,
            fit_bounds=fit_bounds,
            bounds_padding=bounds_padding,
            extent=extent,
            show_layer_elevs=show_layer_elevs,
            hover_surfaces=hover_surfaces,
        )
        if self._map_type is not None:
            # A DERIVED field the free verb already draws: this noun is that
            # map's grammar, so it draws THROUGH it -- one picture, one hover.
            choro = self.model.plot.map(
                per=self._map_per(per),
                layer=layer,
                type=self._map_type,
                hover=hover,
                **drawing,
                **trace_kwargs,
            )
            return _apply_backend(choro, backend)
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
        # The package's other arrays at this layer, listed under the mapped value.
        cell_hover = {**cell_hover, **self._sibling_arrays(int(layer))}
        siblings = tuple(name for name in self.siblings if name in cell_hover)
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
                extra_fields=siblings,
                context_fields=_input_context_fields(cell_hover),
            ),
            **drawing,
            **trace_kwargs,
        )
        return _apply_backend(choro, backend)


class TransmissivityMaxExplorer(StaticArrayFieldExplorer):
    """``model.npf.transmissivity_max`` -- K times the FULL layer thickness.

    The most a layer could conduct, so it needs no heads and draws on a model that
    has not been run. Every verb reads :func:`~myflopy.modflow.mf6.hydraulics
    .model_transmissivity`, and ``map`` draws through ``model.plot.map(type=
    "trans_max")``, so this noun, its mosaic and the free map are one picture
    (ledger 308).
    """

    _map_type = "trans_max"
    _field = "transmissivity_max"
    _label = "Max transmissivity (K x full thickness)"

    def __init__(self, model: SimulationBase):
        """Bind this transmissivity field of ``model``'s NPF package."""

        super().__init__(model, "npf", self._field, label=self._label, colorscale="Viridis")

    def _result(self, kstpkper=None):
        """The :class:`~myflopy.modflow.mf6.hydraulics.Transmissivity` at ``kstpkper``."""

        return model_transmissivity(self.model)

    def _array(self) -> np.ndarray:
        """``(nlay, ncpl)`` transmissivity (the last saved time, when saturated)."""

        return self._result().transmissivity

    def _spatial_layers(self) -> list[int]:
        """Every layer the views number (the base asks ``get()``, which would compute them all)."""

        return list(range(views_nlay(self.model.gwf)))

    def _times(self, per):
        """``[(per, kstpkper)]`` this table covers: one untimed entry for the maximum."""

        return [(None, None)]

    def get(
        self,
        *,
        per: int | Iterable[int] | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Transmissivity, K and thickness as a normalized layer/cell table.

        ``per`` selects stress periods for the SATURATED noun -- each at its last
        saved time step, every saved period by default -- and is ignored by the
        maximum, which has no time. Built with numpy rather than one row at a
        time: a saturated table is ``periods x layers x cells`` long.

        Columns ``model``, ``package``, ``field``, (saturated only) ``per`` and
        ``kstp``, ``layer``, ``cell``, the transmissivity (length^2/time in the
        model's units, e.g. ft2/d), ``k`` (length/time) and ``thickness``
        (length).

        Parameters
        ----------
        per : int or iterable of int, optional
            On ``npf.transmissivity``: the zero-based stress period(s), each read
            at the LAST saved head step of the period; a period with no saved
            heads, or an empty iterable, raises ``ValueError``, and an unrun model
            raises naming ``transmissivity_max``. Default ``None``: every period
            with saved heads. On ``npf.transmissivity_max``: accepted and
            IGNORED -- the maximum has no time, so every value returns the same
            table, with no ``per``/``kstp`` columns.
        layer : int or iterable of int, optional
            Zero-based layer(s) to keep. Default ``None``: every layer.
        cells : int or iterable of int, optional
            Zero-based cell id(s) within a layer to keep. Default ``None``: every
            cell.
        """

        frames = []
        for period, key in self._times(per):
            result = self._result(key)
            nlay, ncpl = result.transmissivity.shape
            frame = pd.DataFrame(
                {
                    "model": self.model.name,
                    "package": self.package_name,
                    "field": self.field_name,
                    "layer": np.repeat(np.arange(nlay), ncpl),
                    "cell": np.tile(np.arange(ncpl), nlay),
                    self.field_name: result.transmissivity.reshape(-1),
                    "k": result.k.reshape(-1),
                    "thickness": result.thickness.reshape(-1),
                }
            )
            if key is not None:
                frame.insert(3, "per", int(period))
                frame.insert(4, "kstp", int(key[0]))
            frames.append(frame)
        frame = pd.concat(frames, ignore_index=True)
        return _filter_normalized_table(frame, layer=layer, cells=cells)


class TransmissivityExplorer(TransmissivityMaxExplorer):
    """``model.npf.transmissivity`` -- K times the SATURATED thickness.

    The default meaning of transmissivity here (the user's decision, ledger 308):
    what the layer conducts at a stress period, following MF6's own ICELLTYPE rule
    -- ``min(head, top) - bottom`` for a convertible layer, the full thickness for a
    confined one. It needs heads, so a model that has not been run raises and
    names ``transmissivity_max``. ``per`` selects the period (its last saved time
    step); ``mosaic(by="layer")`` and ``animate()`` over periods work as on any
    noun.
    """

    _map_type = "trans"
    _field = "transmissivity"
    _label = "Transmissivity (K x saturated thickness)"

    def _map_per(self, per):
        """The map reads the period asked for -- this field has a time. ``None``
        is the free map's default time, as on every other noun."""

        return None if per is None else int(per)

    def long(
        self,
        *,
        per: int | Iterable[int] | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.Series:
        """Transmissivity indexed by ``per/layer/cell`` -- the PERIOD is part of the
        key, since the saturated thickness moves with the heads (the inherited
        static form indexed by ``layer/cell`` alone repeated every key once per
        period, ledger 311).

        Parameters
        ----------
        per : int or iterable of int, optional
            Zero-based stress period(s), each at its last saved head step, as in
            :meth:`get`; a period with no saved heads raises ``ValueError``.
            Default ``None``: every period with saved heads.
        layer : int or iterable of int, optional
            Zero-based layer(s) to keep. Default ``None``: every layer.
        cells : int or iterable of int, optional
            Zero-based cell id(s) within a layer to keep. Default ``None``: every
            cell.
        """

        frame = self.get(per=per, layer=layer, cells=cells)
        series = frame.set_index(["per", "layer", "cell"])[self.field_name].sort_index()
        series.name = self.field_name
        return series

    def stack(self, **kwargs) -> pd.Series:
        """Alias for :meth:`long`.

        Parameters
        ----------
        **kwargs
            Passed to :meth:`long` unchanged: ``per=``, ``layer=`` and ``cells=``,
            by keyword. Anything else raises ``TypeError``.
        """

        return self.long(**kwargs)

    def wide(
        self,
        *,
        per: int | Iterable[int] | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """One row per ``(per, cell)``, one column per layer -- never one period
        silently standing for all of them.

        Columns ``per``, ``cell``, then ``layer_<k>`` for each layer kept.

        Parameters
        ----------
        per : int or iterable of int, optional
            Zero-based stress period(s), each at its last saved head step, as in
            :meth:`get`. Default ``None``: every period with saved heads.
        layer : int or iterable of int, optional
            Zero-based layer(s) to keep as columns. Default ``None``: every layer.
        cells : int or iterable of int, optional
            Zero-based cell id(s) within a layer to keep as rows. Default
            ``None``: every cell.
        """

        frame = self.get(per=per, layer=layer, cells=cells)
        wide = frame.pivot_table(
            index=["per", "cell"], columns="layer", values=self.field_name, aggfunc="first"
        )
        wide.columns = [f"layer_{int(column)}" for column in wide.columns]
        return wide.reset_index()

    def _saved_times(self) -> list[tuple[int, int]]:
        """Each saved period's LAST ``(kstp, kper)``, in order; refuses an unrun model."""

        try:
            keys = [tuple(int(v) for v in key) for key in self.model.hds.kstpkper]
        except (OSError, ValueError, KeyError) as error:
            raise ValueError(
                "npf.transmissivity is the SATURATED transmissivity and needs heads; "
                "this model has none yet. Run it, or use "
                "model.npf.transmissivity_max (K x full thickness)."
            ) from error
        last: dict[int, tuple[int, int]] = {}
        for key in keys:
            last[key[1]] = key
        return [last[period] for period in sorted(last)]

    def _times(self, per):
        """``[(per, kstpkper)]`` for the requested periods (all saved, by default)."""

        saved = {key[1]: key for key in self._saved_times()}
        wanted = sorted(saved) if per is None else self._resolve_periods(per)
        missing = [p for p in wanted if p not in saved]
        if missing:
            raise ValueError(
                f"Stress period(s) {missing} have no saved heads. Available: {sorted(saved)}"
            )
        return [(p, saved[p]) for p in wanted]

    def _result(self, kstpkper=None):
        """Saturated transmissivity at ``kstpkper`` (the last saved time by default)."""

        key = kstpkper if kstpkper is not None else self._saved_times()[-1]
        heads = np.vstack([
            self.model.hds.array(layer=layer, kstpkper=key, masked=False)
            for layer in range(views_nlay(self.model.gwf))
        ])
        return model_transmissivity(self.model, heads=heads)

    def _spatial_periods(self) -> list[int]:
        """The stress periods with saved heads."""

        return [key[1] for key in self._saved_times()]


__all__ = [
    "CellPackageInputsExplorer",
    "CellPackageInputFieldExplorer",
    "UzfFieldInputsExplorer",
    "UzfInputsNamespace",
    "StaticArrayFieldExplorer",
    "TransmissivityExplorer",
    "TransmissivityMaxExplorer",
]

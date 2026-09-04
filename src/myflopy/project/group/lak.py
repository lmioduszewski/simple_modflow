"""Grouped LAK budget/stage/connection results + outputs + accessor."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from myflopy.modflow.mf6.package_explorer import (
    LakStageResultsExplorer,
    _exchange_colorscale,
    _normalize_connection_type_filter,
    _symmetric_color_limit,
    build_cell_input_map_payload,
    build_group_input_compare_map_payload,
    build_lak_budget_result_table,
    build_lak_connection_table,
    build_lak_q_map_payload,
    get_default_group_compare_colorscale,
)
from myflopy.modflow.mf6.package_plotting import (
    _apply_backend,
    refuse_noun_parameters,
    resolve_noun_hover,
)
from myflopy.modflow.mf6.package_surface_water import join_lak_stage
from myflopy.modflow.utils.datatypes.hover import cell_input_hover, compare_hover, lak_hover
from myflopy.project.group._shared import (
    _ensure_group_map_compatible,
    _filter_group_input_table,
    _resolve_group_compare_target,
)
from myflopy.project.group.results import GroupCellPackageResults, GroupCellPackageResultsNamespace
from myflopy.project.group.spatial import _GroupSpatialView

if TYPE_CHECKING:
    from myflopy.project.group.core import ModelGroup


class GroupLakBudgetResults(GroupCellPackageResults):
    """Grouped LAK exchange accessor with area-normalized maps."""

    def __init__(self, group: ModelGroup, *, budget_text: str = "GWF", value_name: str = "q_lake"):
        """Bind a grouped LAK exchange accessor (defaults to the ``GWF`` term's ``q_lake``).

        LAK's exchange is feature-referenced, so the emitted column is
        ``q_lake`` (accessor stays ``results.q``).
        """

        super().__init__(group, "lak", budget_text=budget_text, value_name=value_name)

    def get(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
        connection_type: str | Sequence[str] | None = None,
    ) -> pd.DataFrame:
        """Return aligned LAK exchange rows, including ``q_per_area``."""

        frames = []
        for current_model_name, model in self.group.models.items():
            frame = build_lak_budget_result_table(
                model,
                budget_text=self.budget_text,
                value_name=self.value_name,
            )
            frame["model"] = current_model_name
            frames.append(frame)

        if not frames:
            return pd.DataFrame()
        combined = pd.concat(frames, ignore_index=True)
        combined = _filter_group_input_table(
            combined,
            model_name=model_name,
            per=per,
            layer=layer,
            cells=cells,
        )
        connection_types = _normalize_connection_type_filter(connection_type)
        if connection_types is not None and "claktype" in combined.columns:
            combined = combined.loc[combined["claktype"].astype("string").str.upper().isin(connection_types)].copy()
        return combined.reset_index(drop=True)

    def compare(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
        connection_type: str | Sequence[str] | None = None,
    ) -> pd.DataFrame:
        """Compare grouped LAK exchange rows against the reference model."""

        data = self.get(per=per, layer=layer, cells=cells, connection_type=connection_type)
        if data.empty:
            return pd.DataFrame()

        reference = self.group.reference
        q = self.value_name  # feature-referenced exchange column (q_lake)
        value_columns = [q, "q_per_area"] if "q_per_area" in data.columns else [q]
        key_columns = [column for column in data.columns if column not in {"model", *value_columns}]
        ref = data[data["model"] == reference].drop(columns=["model"]).copy()
        rename_map = {column: f"reference_{column}" for column in value_columns}
        ref = ref.rename(columns=rename_map)
        comp = data[data["model"] != reference].merge(ref, on=key_columns, how="inner")
        comp["reference_model"] = reference
        comp[f"{q}_diff"] = comp[q].astype(float) - comp[f"reference_{q}"].astype(float)
        if "q_per_area" in value_columns:
            comp["q_per_area_diff"] = comp["q_per_area"].astype(float) - comp["reference_q_per_area"].astype(float)
        ordered = ["model", "reference_model", *key_columns, q, f"reference_{q}", f"{q}_diff"]
        if "q_per_area" in value_columns:
            ordered.extend(["q_per_area", "reference_q_per_area", "q_per_area_diff"])
        result = comp[ordered]
        return _filter_group_input_table(
            result,
            model_name=model_name,
            per=per,
            layer=layer,
            cells=cells,
        )

    def _spatial_map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        model=None,
        connection_type: str | Sequence[str] | None = None,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Render one model's LAK exchange (area-normalized) as a ``Choro``."""

        del agg
        target_name = self._group_target(model)
        target_model = self.group.models[target_name]
        selected = join_lak_stage(
            target_model,
            self.get(model_name=target_name, per=per, layer=layer, connection_type=connection_type),
            per=per,
        )
        values, hover = build_lak_q_map_payload(
            selected,
            ncpl=target_model.vor.ncpl,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            value_column=self.value_name,
        )
        absmax = _symmetric_color_limit(values)
        kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
        kwargs.setdefault("zmax", absmax if absmax > 0 else None)
        kwargs.setdefault("zmid", 0.0)
        kwargs.setdefault("hover_spec", lak_hover())
        return target_model.plot.map(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            # match the SFR convention: gaining (negative q) blue, losing red
            # LAK exchange is feature-referenced (gaining is POSITIVE), so blue
            # sits at the positive end -- matching the single-model LAK map. Using
            # the raw blue-at-negative scale here inverted the group map (gaining
            # drew red) while the single-model map drew it blue.
            colorscale=colorscale or _exchange_colorscale("feature"),
            **kwargs,
        )

    def compare_map(
        self,
        *,
        model_name: str | None = None,
        per: int = 0,
        layer: int = 0,
        connection_type: str | Sequence[str] | None = None,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Build a grouped LAK diff map using normalized exchange per area."""

        del agg
        target_name = _resolve_group_compare_target(self.group, model_name)
        _ensure_group_map_compatible(self.group, target_name)
        reference_model = self.group.models[self.group.reference]
        data = self.get(per=per, layer=layer, connection_type=connection_type)

        def _normalized_by_cell(frame: pd.DataFrame, current_model_name: str) -> pd.DataFrame:
            """Per-cell exchange per unit connection area (``sum(q)/sum(flow_area)``) for one model."""

            selected = frame[frame["model"] == current_model_name].copy()
            if selected.empty:
                return pd.DataFrame(columns=["cell", "q_per_area"])
            q = self.value_name
            selected[q] = pd.to_numeric(selected[q], errors="coerce")
            selected["flow_area"] = pd.to_numeric(selected["flow_area"], errors="coerce")
            grouped = selected.groupby("cell", as_index=False).agg({q: "sum", "flow_area": "sum"})
            grouped["q_per_area"] = np.where(grouped["flow_area"] > 0.0, grouped[q] / grouped["flow_area"], np.nan)
            return grouped[["cell", "q_per_area"]]

        reference = _normalized_by_cell(data, self.group.reference).rename(
            columns={"q_per_area": "reference_q_per_area"}
        )
        target = _normalized_by_cell(data, target_name)
        comparison = target.merge(reference, on="cell", how="inner")
        comparison["model"] = target_name
        comparison["reference_model"] = self.group.reference
        comparison["per"] = int(per)
        comparison["layer"] = int(layer)
        comparison["q_per_area_diff"] = (
            comparison["q_per_area"].astype(float) - comparison["reference_q_per_area"].astype(float)
        )
        values, hover, absmax = build_group_input_compare_map_payload(
            comparison,
            ncpl=reference_model.vor.ncpl,
            value_column="q_per_area",
            diff_column="q_per_area_diff",
            model_name=target_name,
            reference_model=self.group.reference,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg="sum",
        )
        kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
        kwargs.setdefault("zmax", absmax if absmax > 0 else None)
        kwargs.setdefault("zmid", 0.0)
        kwargs.setdefault(
            "hover_spec",
            compare_hover(
                "q_per_area",
                "q_per_area_diff",
                title="Δ lake exchange vs reference",
                units={
                    "q_per_area": "ft/d",
                    "reference_q_per_area": "ft/d",
                    "q_per_area_diff": "ft/d",
                },
                labels={"q_per_area_diff": "Δ q / area"},
            ),
        )
        return reference_model.plot.map(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or get_default_group_compare_colorscale(),
            **kwargs,
        )


class GroupLakOutputs:
    """Lake-output accessor for :class:`ModelGroup`."""

    def __init__(self, group: ModelGroup):
        """Bind the grouped lake-output accessor to ``group``."""

        self.group = group

    def stage(self) -> pd.DataFrame:
        """Return lake stages for all models as one aligned long-format table."""

        rows: list[pd.DataFrame] = []
        for model_name, model in self.group.models.items():
            stage_data = np.asarray(model.outputs.lak.stage.get(), dtype=float)
            if stage_data.ndim == 1:
                stage_data = stage_data.reshape(-1, 1)
            periods = list(model.kstpkper)
            frame = pd.DataFrame(stage_data, columns=list(range(stage_data.shape[1])))
            frame["kstpkper"] = periods[: len(frame)]
            frame = frame.melt(id_vars="kstpkper", var_name="lake", value_name="stage")
            frame["model"] = model_name
            rows.append(frame)

        if not rows:
            return pd.DataFrame(columns=["model", "kstpkper", "lake", "stage"])
        return pd.concat(rows, ignore_index=True)[["model", "kstpkper", "lake", "stage"]]


class GroupLakStageResults(_GroupSpatialView):
    """Grouped accessor for lake stages and stage comparisons.

    Inherits the :class:`SpatialView` grammar (``map``/``mosaic``/``animate``)
    with the model axis being the group's members; each panel delegates to the
    single-model LAK stage explorer so grouped stage maps match the single-model
    ``lak.results.stage.map`` exactly.
    """

    def __init__(self, group: ModelGroup):
        """Bind the grouped lake stage-results accessor to ``group``."""

        self.group = group

    def _spatial_value_label(self) -> str:
        """The mapped quantity's label -- lake stage."""

        return "stage"

    def _spatial_map(self, *, per: int = 0, layer: int = 0, model=None, **kwargs):
        """Render one model's lake stage broadcast to connected cells."""

        target_model = self.group.models[self._group_target(model)]
        return LakStageResultsExplorer(target_model).map(
            per=int(per), layer=int(layer), backend="plotly", **kwargs
        )

    def get(
        self,
        *,
        model_name: str | None = None,
        lake: int | None = None,
        per: int | None = None,
    ) -> pd.DataFrame:
        """Return aligned lake stages for all models."""

        rows: list[pd.DataFrame] = []
        for current_model_name, model in self.group.models.items():
            stage_data = np.asarray(model.outputs.lak.stage.get(), dtype=float)
            if stage_data.ndim == 1:
                stage_data = stage_data.reshape(-1, 1)
            periods = pd.DataFrame(stage_data, columns=list(range(stage_data.shape[1])))
            periods["per"] = periods.index.astype(int)
            frame = periods.melt(id_vars="per", var_name="lake", value_name="stage")
            frame["lake"] = frame["lake"].astype(int)
            frame["model"] = current_model_name
            rows.append(frame[["model", "per", "lake", "stage"]])

        if not rows:
            return pd.DataFrame(columns=["model", "per", "lake", "stage"])
        combined = pd.concat(rows, ignore_index=True)
        if model_name is not None:
            combined = combined.loc[combined["model"] == str(model_name)].copy()
        if lake is not None:
            combined = combined.loc[combined["lake"] == int(lake)].copy()
        if per is not None:
            combined = combined.loc[combined["per"] == int(per)].copy()
        return combined.reset_index(drop=True)

    def compare(
        self,
        *,
        model_name: str | None = None,
        lake: int | None = None,
        per: int | None = None,
    ) -> pd.DataFrame:
        """Compare lake stages against the reference model."""

        data = self.get(lake=lake, per=per)
        if data.empty:
            return pd.DataFrame(columns=["model", "reference_model", "per", "lake", "stage", "reference_stage", "stage_diff"])

        reference = self.group.reference
        ref = (
            data.loc[data["model"] == reference, ["per", "lake", "stage"]]
            .rename(columns={"stage": "reference_stage"})
            .copy()
        )
        comp = data.loc[data["model"] != reference].merge(ref, on=["per", "lake"], how="inner")
        comp["reference_model"] = reference
        comp["stage_diff"] = comp["stage"].astype(float) - comp["reference_stage"].astype(float)
        if model_name is not None:
            comp = comp.loc[comp["model"] == str(model_name)].copy()
        return comp[["model", "reference_model", "per", "lake", "stage", "reference_stage", "stage_diff"]]

    # NOTE: the series view is the unified grammar's ``plot()`` (SpatialView) --
    # one line per model and lake.


class GroupLakConnections:
    """Grouped accessor for lake-connection geometry."""

    def __init__(self, group: ModelGroup):
        """Bind the grouped lake-connection-geometry accessor to ``group``."""

        self.group = group

    def get(
        self,
        *,
        model_name: str | None = None,
        lake: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """Return aligned LAK connection rows for all models."""

        frames = []
        for current_model_name, model in self.group.models.items():
            frame = build_lak_connection_table(model)
            frame["model"] = current_model_name
            frames.append(frame)

        if not frames:
            return pd.DataFrame()
        combined = pd.concat(frames, ignore_index=True)
        combined = _filter_group_input_table(
            combined,
            model_name=model_name,
            per=None,
            layer=layer,
            cells=cells,
        )
        if lake is not None and "lake" in combined.columns:
            combined = combined.loc[combined["lake"] == int(lake)].copy()
        return combined.reset_index(drop=True)

    def map(
        self,
        *,
        # -- which model, which connections -----------------------------------
        model_name: str | None = None,
        lake: int | None = None,
        layer: int = 0,
        # -- table to cells ---------------------------------------------------
        value_column: str = "connection_area",
        agg: str = "sum",
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        # -- colour -----------------------------------------------------------
        zmin: float | None = None,
        zmax: float | None = None,
        colorscale: str | list | tuple | None = None,
        logscale: bool = False,
        # -- contours ---------------------------------------------------------
        contours: bool | str = False,
        contour_values=None,
        contour_levels: int | float | list = 10,
        contour_color: str = "black",
        contour_width: float = 1.5,
        contour_name: str | None = None,
        contour_clip: bool = True,
        contour_resolution: int = 150,
        contour_method: str = "linear",
        # -- highlighting -----------------------------------------------------
        select=None,
        select_style: str = "outline",
        select_color: str | None = None,
        # -- overlays and framing ---------------------------------------------
        locs=None,
        hillshade_path=None,
        fit_bounds: bool = True,
        bounds_padding: float = 0.05,
        # -- hover ------------------------------------------------------------
        hover=None,
        # -- renderer ---------------------------------------------------------
        backend: str = "plotly",
        **trace_kwargs,
    ):
        """LAK connection geometry for ONE member model, as a choropleth.

        A GROUP noun drawing a single picture. The group's axis is the model
        list, but a choropleth has one grid, so ``model_name`` picks the member
        and everything else describes that member's map. The rows come from
        :meth:`get` -- one per LAK connection -- and are reduced to one value
        per cell by ``agg``, which is why this signature carries
        ``value_column``/``agg``/``multiplier``/``fill_value`` on top of the
        shared drawing parameters. Use ``mosaic`` on a results noun when you
        want every member side by side; connection geometry is static input, so
        the useful comparison is one member at a time against another.

        Every parameter is named rather than swept into ``**kwargs``: PyCharm
        and Pylance read the ``def`` line and never run the module, so a
        parameter that arrives through a tail is one no editor can ever offer
        (plan 8.8). ``backend`` in particular used to reach the renderer only
        through that tail.

        The per-layer arguments the free verb has -- ``kstpkper``,
        ``per_timestep``, ``bgs``, ``hover_layers``, ``hover_surfaces``,
        ``show_layer_elevs``, ``show_mounding`` -- are deliberately absent.
        Connection geometry has no time axis and no per-layer profile behind
        it, so each was measured to leave the figure byte-identical, and
        ``show_mounding`` is worse than inert: it injects a head-derived row
        into a geometry tooltip on every cell.

        Parameters
        ----------
        model_name : str, optional
            Which member of the group to draw. Defaults to the group's
            reference model (``group.reference``). Pass a member name to draw
            that model's connections instead -- the usual reason is to look at
            the same lake on a refined grid. A name that is not in the group
            raises ``KeyError`` rather than silently falling back.
        lake : int, optional
            Draw only one lake's connections, by zero-based lake number
            (``ifno``). Defaults to ``None``, meaning every lake in the
            package, which on a multi-lake model paints them all into one
            colour range. Set it when one lake's geometry is the question, or
            when a large lake's areas are flattening a small one.
        layer : int, default 0
            Zero-based layer whose connections are drawn. A lake connects
            downward through several layers; only the rows in this layer are
            mapped, so pass the layer the connections you care about actually
            sit in -- ``get()``'s ``layer`` column lists them.
        value_column : str, default 'connection_area'
            Which connection column to colour by. ``connection_area`` is the
            physical interface area (plan-view cell area for a vertical
            connection, ``connwidth * (telev - belev)`` for a horizontal one).
            Any other numeric column of :meth:`get` works -- ``belev``,
            ``telev``, ``connlen``, ``connwidth`` -- and a name that is not in
            the table raises ``KeyError`` naming it.
        agg : {'sum', 'mean', 'min', 'max', 'first', 'last'}, default 'sum'
            How several connections landing in ONE cell are reduced to one
            number. ``sum`` is right for ``connection_area``, which is
            extensive -- a cell touched by two lakes offers both areas. Switch
            to ``mean`` for an intensive column such as ``belev`` or
            ``connlen``, where adding two elevations produces a meaningless
            number.
        multiplier : float, default 1.0
            Scale every value before drawing -- unit conversion (``0.3048`` for
            feet to metres on an elevation column), or a sign flip. Applied
            before ``agg``.
        fill_value : float, default 0.0
            Value given to cells with no LAK connection in this layer. The
            default colours them zero, which reads as "no lake here" on an
            area map; pass ``float("nan")`` to leave them uncoloured instead,
            which is what you want on an elevation column where zero is a real
            elevation and would distort the colour range.
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
        get : the connection rows behind the picture, as a DataFrame, aligned
            across every model in the group.

        Examples
        --------
        >>> group.packages.lak.connections.map()
        >>> group.packages.lak.connections.map(model_name="refined")
        >>> group.packages.lak.connections.map(lake=0, layer=1)
        >>> group.packages.lak.connections.map(
        ...     value_column="belev", agg="mean", fill_value=float("nan")
        ... )
        >>> group.packages.lak.connections.map(logscale=True, colorscale="Viridis")
        >>> group.packages.lak.connections.map(backend="mpl").savefig("lak_conn.png")
        """

        refuse_noun_parameters(
            "group.packages.lak.connections",
            value_column,
            trace_kwargs,
        )
        hover = resolve_noun_hover(hover, trace_kwargs)
        if "per" in trace_kwargs:
            # Every SIBLING noun takes `per`, so reaching for it here is the
            # natural mistake -- but connection geometry is static and the
            # forward hardcodes `per=0`, so the name collides in the tail and
            # dies as "ModelPlots.map() got multiple values for keyword
            # argument 'per'", naming a class the caller never typed.
            raise TypeError(
                "group.packages.lak.connections.map() has no per=: LAK "
                "connection geometry does not change with time, so there is "
                "one picture for every stress period. Use "
                "`group.packages.lak.results.<noun>` for a quantity that does."
            )
        target_name = self.group.reference if model_name is None else str(model_name)
        if target_name not in self.group.models:
            raise KeyError(f"Model {target_name!r} is not in the group.")
        target_model = self.group.models[target_name]
        selected = self.get(model_name=target_name, lake=lake, layer=layer)
        if value_column not in selected.columns:
            raise KeyError(f"LAK connection column {value_column!r} was not found.")
        values, hover_table = build_cell_input_map_payload(
            selected,
            ncpl=target_model.vor.ncpl,
            value_column=value_column,
            per=None,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        choro = target_model.plot.map(
            per=0,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover_table,
            hover_heads=False,
            hover_ks=False,
            # The noun's OWN spec is the picture's BASE hover, so it goes in
            # `hover_spec=`; `hover=` is the call-site override slot and stays the
            # caller's. Putting the default in the override slot renders the same
            # but leaves `choro.hover_spec` None, which disagrees with `model.hds`
            # and with what this noun did before the 8.8 conversion.
            hover=hover,
            hover_spec=cell_input_hover(value_column),
            zmin=zmin,
            zmax=zmax,
            colorscale=colorscale or "earth",
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


class GroupLakResultsNamespace(GroupCellPackageResultsNamespace):
    """Namespace for grouped LAK result accessors."""

    def _field_names(self):
        """The mappable grouped LAK result fields: exchange ``q`` and ``stage``."""

        return ["q", "stage"]

    @property
    def stage(self) -> GroupLakStageResults:
        """Return grouped LAK stage helpers."""

        return GroupLakStageResults(self._result_accessor.group)

    @property
    def q(self) -> GroupLakBudgetResults:
        """Return grouped LAK exchange helpers with area-normalized map behavior."""

        return self._result_accessor


class GroupLakPackageAccessor:
    """Namespace for grouped LAK geometry and result helpers."""

    def __init__(self, group: ModelGroup, results_namespace: GroupLakResultsNamespace):
        """Bind the grouped LAK package accessor (connections + results) to ``group``."""

        self.group = group
        self._results_namespace = results_namespace

    @property
    def connections(self) -> GroupLakConnections:
        """Return grouped LAK connection-geometry helpers."""

        return GroupLakConnections(self.group)

    @property
    def results(self) -> GroupLakResultsNamespace:
        """Return grouped LAK result helpers."""

        return self._results_namespace




from myflopy.modflow.mf6.package_plotting import NOUN_MAP_PARAMS  # noqa: E402
from myflopy.plot import inherit_map_docs  # noqa: E402


# --- the shared map reference, spliced onto this module's noun ----------------
#
# `myflopy.plot` (layer 7) splices every noun BELOW it at import. This module is
# layer 9, so reaching down for the reference is this module's job: doing it from
# `plot` pulls in `myflopy.project` (layer 13) while `simulation.base` is still
# half-built, and the import fails outright. So the reference is fetched HERE
# instead, at the bottom of the module -- module-level (this file is layer 9 and
# both targets are below it, so the ratchet stays where it was) but placed after
# the class, which must exist before it can be spliced.
def _inherit_group_lak_docs() -> None:
    """Give ``GroupLakConnections.map`` the free verb's entries for what it takes."""

    # `keep` is exactly the shared drawing set. Deliberately NOT `| {"per",
    # "layer"}`:
    #   `per`   -- this noun has no `per` parameter (connection geometry is
    #              static, so the forward hardcodes `per=0`). Inheriting its
    #              entry documented a knob that does not exist, and
    #              `map(per=3)` reached the tail and died as
    #              `ModelPlots.map() got multiple values for keyword argument
    #              'per'` -- a TypeError naming a class the caller never typed,
    #              the exact defect `refuse_noun_parameters` was written against.
    #   `layer` -- documented locally above, with what it means for a LAK
    #              connection. Inheriting it too put two `layer` entries in one
    #              Parameters block.
    inherit_map_docs(GroupLakConnections.map, keep=set(NOUN_MAP_PARAMS))


_inherit_group_lak_docs()

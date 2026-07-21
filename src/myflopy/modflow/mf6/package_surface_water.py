"""LAK, SFR, and combined surface-water package explorers."""

from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from myflopy import viz as figs
from myflopy._deprecation import deprecated_instance_getattr
from myflopy.viz import mpl_axes

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.modflow.mf6.package_budget import (
    build_lak_budget_result_table,
    build_lak_budget_term_table,
    build_lak_stage_change_table,
    build_lak_stage_result_table,
    build_sfr_budget_result_table,
    build_sfr_budget_term_table,
    build_sfr_long_profile_table,
    build_sfr_stage_result_table,
    build_surface_water_exchange_cell_table,
)
from myflopy.modflow.mf6.package_explorer_utils import (
    _default_show_layer_elevs,
    _filter_normalized_table,
    _normalize_connection_type_filter,
)
from myflopy.modflow.mf6.package_plotting import (
    FieldMappable,
    SpatialView,
    _apply_backend,
    _blue_white_red_diverging_colorscale,
    _exchange_colorscale,
    _symmetric_color_limit,
    build_cell_input_map_payload,
    build_lak_q_map_payload,
    build_sfr_q_map_payload,
    build_surface_water_q_map_payload,
)
from myflopy.modflow.mf6.package_registry import (
    get_default_budget_term,
    get_package_result_spec,
)
from myflopy.modflow.mf6.package_results import (
    CellBudgetResultsExplorer,
    PackageBudgetTermExplorer,
    StageResultsExplorer,
)
from myflopy.modflow.mf6.package_tables import (
    build_lak_connection_table,
    build_lak_input_table,
    build_sfr_input_table,
    summarize_input_table,
)
from myflopy.modflow.utils.datatypes.hover import (
    cell_input_hover,
    lak_hover,
    sfr_hover,
    surface_water_hover,
)


def _join_feature_stage(frame, stage_table, *, per=None):
    """Join per-cell feature stage onto an exchange frame for the hover.

    ``stage_table`` is a normalized stage result table (one row per feature-cell
    per period). Cells touched by more than one feature take the mean stage,
    matching :class:`StageResultsExplorer`'s series aggregate.
    """

    if frame.empty or "stage" in frame.columns or stage_table.empty:
        return frame
    selected = stage_table
    if per is not None and "per" in selected.columns:
        selected = selected[selected["per"] == int(per)]
    if selected.empty or "stage" not in selected.columns:
        return frame
    mapping = (
        pd.to_numeric(selected["stage"], errors="coerce")
        .groupby(selected["cell"])
        .mean()
    )
    joined = frame.copy()
    joined["stage"] = joined["cell"].map(mapping)
    return joined


def join_sfr_stage(model, frame, *, per=None):
    """Add SFR reach stage to an exchange frame; a no-op if stage is unavailable."""

    try:
        stage_table = build_sfr_stage_result_table(model)
    except Exception:
        return frame
    return _join_feature_stage(frame, stage_table, per=per)


def join_lak_stage(model, frame, *, per=None):
    """Add lake stage to an exchange frame; a no-op if stage is unavailable."""

    try:
        stage_table = build_lak_stage_result_table(model)
    except Exception:
        return frame
    return _join_feature_stage(frame, stage_table, per=per)


class SfrBudgetResultsExplorer(CellBudgetResultsExplorer):
    """SFR-specific result explorer with reach-profile helpers."""

    def __init__(
        self,
        model: SimulationBase,
        *,
        budget_text: str = "SFR",
        value_name: str = "q_gwf",
    ):
        """Bind an SFR budget-result explorer (defaults to the ``SFR`` term's ``q_gwf``).

        ``value_name`` is the emitted column (``q_gwf`` -- the SFR exchange is
        aquifer-referenced); the accessor stays ``results.q`` via ``result_name``.
        """

        super().__init__(model, "sfr", budget_text, value_name, result_name="q")

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return the normalized SFR budget-result table for selected rows."""

        frame = build_sfr_budget_result_table(
            self.model,
            budget_text=self.budget_text,
            value_name=self.value_name,
        )
        return _filter_normalized_table(frame, per=per, layer=layer, cells=cells)

    def profile(self, *, per: int = 0) -> pd.DataFrame:
        """Return one reach-ordered profile table for the selected period."""

        frame = self.get(per=per)
        if frame.empty:
            return frame
        return frame.sort_values(["reach", "cell"]).reset_index(drop=True)

    def map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        backend: str = "plotly",
        **kwargs,
    ):
        """Build an SFR exchange choropleth normalized by total reach length.

        Notes
        -----
        The mapped value is ``sum(q) / sum(rlen)`` within each cell. This makes
        SFR exchange maps less sensitive to cells that only appear larger
        because they contain longer stream reaches. The ``agg`` argument is
        accepted for API compatibility but is not used because the
        normalization is computed explicitly from total exchange and total
        length per cell.

        MF6 reports the SFR ``GWF`` budget term as flow from the stream reach
        to the groundwater cell. Positive values therefore indicate losing
        reaches, while negative values indicate gaining reaches. The default
        diverging colorscale is defined explicitly so gaining reaches plot blue
        and losing reaches plot red.
        """

        del agg
        selected = join_sfr_stage(self.model, self.get(per=per, layer=layer), per=per)
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(self.model))
        values, hover = build_sfr_q_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
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
        kwargs.setdefault("hover_spec", sfr_hover())
        choro = self.model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            # SFR's cell record is flow FROM reach TO cell (the "gwf" frame), so
            # gaining is negative. The orientation is derived from the registry,
            # not assumed, so it cannot drift from the data.
            colorscale=colorscale or _exchange_colorscale(_exchange_frame("sfr")),
            **kwargs,
        )
        return _apply_backend(choro, backend)

    def plot_profile(
        self,
        *,
        per: int = 0,
        x: str = "distance",
        plot_fig: bool = False,
        return_fig: bool = True,
    ):
        """Plot one SFR exchange profile by reach or cumulative stream distance.

        Parameters
        ----------
        per
            Zero-based stress period to plot.
        x
            Either ``"distance"`` for cumulative stream distance or
            ``"reach"`` for raw reach number.
        plot_fig
            If ``True``, call ``show()`` on the created figure.
        return_fig
            If ``True``, return the created figure.
        """

        frame = self.profile(per=per)
        fig = figs.Fig()
        if frame.empty:
            if plot_fig:
                fig.show()
            if return_fig:
                return fig
            return None

        x_column = "distance_mid" if str(x).lower() == "distance" else "reach"
        if x_column not in frame.columns:
            raise KeyError(f"Profile x-axis column {x_column!r} was not found.")
        fig.add_scattergl(
            x=frame[x_column],
            y=frame[self.value_name],
            mode="lines+markers",
            name=f"SFR {self.value_name} per {per}",
            customdata=np.column_stack([frame["reach"], frame["cell"]]),
            hovertemplate=(
                "reach=%{customdata[0]}<br>"
                "cell=%{customdata[1]}<br>"
                f"{self.value_name}=%{{y}}<extra></extra>"
            ),
        )
        fig.update_layout(
            xaxis_title="Stream Distance" if x_column == "distance_mid" else "Reach",
            yaxis_title=self.value_name,
            title=f"SFR {self.value_name} profile (per={per})",
        )
        if plot_fig:
            fig.show()
        if return_fig:
            return fig
        return None


class LakBudgetResultsExplorer(CellBudgetResultsExplorer):
    """LAK-specific result explorer with area-normalized exchange maps."""

    def __init__(
        self,
        model: SimulationBase,
        *,
        budget_text: str = "GWF",
        value_name: str = "q_lake",
    ):
        """Bind a LAK budget-result explorer (defaults to the ``GWF`` exchange ``q_lake``).

        ``value_name`` is the emitted column (``q_lake`` -- LAK's exchange is
        feature-referenced); the accessor stays ``results.q`` via ``result_name``.
        """

        super().__init__(model, "lak", budget_text, value_name, result_name="q")

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
        connection_type: str | Iterable[str] | None = None,
    ) -> pd.DataFrame:
        """Return the normalized LAK budget-result table for selected rows."""

        frame = build_lak_budget_result_table(
            self.model,
            budget_text=self.budget_text,
            value_name=self.value_name,
        )
        frame = _filter_normalized_table(frame, per=per, layer=layer, cells=cells)
        connection_types = _normalize_connection_type_filter(connection_type)
        if connection_types is not None and "claktype" in frame.columns:
            frame = frame.loc[
                frame["claktype"].astype("string").str.upper().isin(connection_types)
            ].copy()
        return frame.reset_index(drop=True)

    def map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        connection_type: str | Iterable[str] | None = None,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        backend: str = "plotly",
        **kwargs,
    ):
        """Build a LAK exchange choropleth normalized by flow-surface area.

        Notes
        -----
        The mapped value is ``sum(q) / sum(flow_area)`` within each cell. This
        yields a signed exchange intensity in length-per-time units instead of
        raw volumetric exchange, which would otherwise scale with lake
        connection area.
        """

        del agg
        selected = join_lak_stage(
            self.model,
            self.get(per=per, layer=layer, connection_type=connection_type),
            per=per,
        )
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(self.model))
        values, hover = build_lak_q_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
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
        choro = self.model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            # NOT the SFR frame: LAK's GWF record comes from the LAK package
            # budget, written from the LAKE's perspective (the "feature" frame),
            # so gaining is POSITIVE. This map once shipped inverted -- losing
            # lakes drew blue -- because a frame literal was copied from SFR;
            # deriving it from the registry removes that whole failure mode.
            colorscale=colorscale or _exchange_colorscale(_exchange_frame("lak")),
            **kwargs,
        )
        return _apply_backend(choro, backend)

    def budget_summary(
        self,
        *,
        per: int | None = None,
        connection_type: str | Iterable[str] | None = None,
    ) -> pd.DataFrame:
        """Summarize lake-groundwater exchange by period, lake, and connection type.

        Parameters
        ----------
        per
            Optional zero-based stress period filter. When omitted, all periods
            are summarized.

        Returns
        -------
        pandas.DataFrame
            Summary table with signed volumetric exchange ``q`` and
            area-normalized exchange ``q_per_area`` grouped by period, lake,
            and connection type.
        """

        frame = self.get(per=per, connection_type=connection_type)
        q = self.value_name  # feature-referenced exchange column (q_lake)
        if frame.empty:
            return pd.DataFrame(
                columns=[
                    "per",
                    "lake",
                    "claktype",
                    "record_count",
                    q,
                    "flow_area",
                    "q_per_area",
                ]
            )
        summary = (
            frame.groupby(["per", "lake", "claktype"], dropna=False, as_index=False)
            .agg(
                record_count=("cell", "size"),
                **{q: (q, "sum")},
                flow_area=("flow_area", "sum"),
            )
            .sort_values(["per", "lake", "claktype"])
            .reset_index(drop=True)
        )
        summary["q_per_area"] = np.where(
            pd.to_numeric(summary["flow_area"], errors="coerce") > 0.0,
            pd.to_numeric(summary[q], errors="coerce")
            / pd.to_numeric(summary["flow_area"], errors="coerce"),
            np.nan,
        )
        return summary

    def plot_budget(
        self,
        *,
        per: int = 0,
        connection_type: str | Iterable[str] | None = None,
        value: str = "q_per_area",
        ax=None,
        return_fig: bool = True,
    ):
        """Plot a compact LAK budget summary by connection type for one period.

        Parameters
        ----------
        per
            Zero-based stress period to plot.
        value
            Summary field to visualize. Supported values are ``"q"``,
            ``"flow_area"``, and ``"q_per_area"``.
        ax
            Optional Matplotlib axes object to draw onto.
        return_fig
            If ``True``, return the created figure.
        """

        summary = self.budget_summary(per=per, connection_type=connection_type)
        # Accept the accessor-facing "q" as an alias for the frame-named column.
        if value == "q":
            value = self.value_name
        if value not in {self.value_name, "flow_area", "q_per_area"}:
            raise ValueError(
                f"value must be one of: 'q', {self.value_name!r}, 'flow_area', 'q_per_area'"
            )
        if ax is None:
            fig, ax = mpl_axes(figsize=(8, 4))
        else:
            fig = ax.figure
        if summary.empty:
            ax.set_title(f"LAK {value} summary (per={per})")
            ax.set_xlabel("Connection Type")
            ax.set_ylabel(value)
            if return_fig:
                return fig
            return None

        summary = summary.copy()
        labels = summary.apply(
            lambda row: f"Lake {int(row['lake'])}\n{row['claktype']}",
            axis=1,
        )
        colors = []
        if value == "flow_area":
            colors = ["#4c78a8" for _ in range(len(summary))]
        else:
            for current in pd.to_numeric(summary[value], errors="coerce").fillna(0.0):
                colors.append("#1f77b4" if current >= 0.0 else "#d62728")
        ax.bar(labels, summary[value].astype(float), color=colors)
        ax.axhline(0.0, color="black", linewidth=1.0, alpha=0.6)
        ax.set_title(f"LAK {value} summary (per={per})")
        ax.set_xlabel("Lake / Connection Type")
        ax.set_ylabel(value)
        ax.tick_params(axis="x", rotation=0)
        fig.tight_layout()
        if return_fig:
            return fig
        return None


class LakStageResultsExplorer(StageResultsExplorer):
    """LAK stage explorer -- the unified grammar's ``plot()`` draws one line
    per lake by stress period (replaced the old ``plot_timeseries``)."""

    def __init__(self, model: SimulationBase):
        """Bind a LAK stage-result explorer using the lake stage table builder."""

        super().__init__(model, "lak", build_lak_stage_result_table)


class LakStageChangeExplorer:
    """Explorer for lake-stage changes between stress periods."""

    def __init__(self, model: SimulationBase):
        """Bind the lake stage-change explorer to ``model``."""

        self.model = model

    def get(
        self,
        *,
        lake: int | None = None,
        per0: int | None = None,
        per1: int | None = None,
    ) -> pd.DataFrame:
        """Return stage-change rows for one lake and/or one period transition."""

        frame = build_lak_stage_change_table(self.model)
        if lake is not None:
            frame = frame.loc[frame["lake"] == int(lake)].copy()
        if per0 is not None:
            frame = frame.loc[frame["per0"] == int(per0)].copy()
        if per1 is not None:
            frame = frame.loc[frame["per1"] == int(per1)].copy()
        return frame.reset_index(drop=True)

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of available lake-stage transitions."""

        frame = self.get()
        if frame.empty:
            return pd.DataFrame(
                [
                    {
                        "label": "lak.results.stage_change",
                        "records": 0,
                        "lakes": 0,
                        "transitions": 0,
                    }
                ]
            )
        return pd.DataFrame(
            [
                {
                    "label": "lak.results.stage_change",
                    "records": int(len(frame)),
                    "lakes": int(frame["lake"].nunique()),
                    "transitions": int(
                        frame[["per0", "per1"]].drop_duplicates().shape[0]
                    ),
                }
            ]
        )

    def plot(
        self,
        *,
        lake: int | None = None,
        backend: str = "plotly",
        title: str | None = None,
    ):
        """Series panel: stage change by stress-period transition, per lake.

        The grammar's ``plot`` verb for this node (its x-axis is the period
        *transition*, not the period, so it does not use the generic series
        engine). ``backend="plotly"`` returns a ``viz.Fig``; ``backend="mpl"``
        a matplotlib figure.
        """

        from myflopy.viz import Fig

        frame = self.get(lake=lake)
        heading = title or "LAK stage change by transition"
        normalized = str(backend).lower()
        if normalized in ("plotly", "interactive"):
            fig = Fig()
            for lake_id, group in frame.groupby("lake", dropna=False):
                labels = [
                    f"{int(start)}->{int(end)}"
                    for start, end in zip(group["per0"], group["per1"], strict=False)
                ]
                fig.add_scatter(
                    x=labels,
                    y=group["stage_change"].astype(float).to_numpy(),
                    mode="lines+markers",
                    name=f"Lake {int(lake_id)}",
                )
            fig.update_layout(
                title=heading,
                xaxis_title="Stress-Period Transition",
                yaxis_title="Stage Change",
            )
            return fig
        if normalized not in ("mpl", "matplotlib", "static"):
            raise ValueError(f"backend must be 'plotly' or 'mpl', got {backend!r}.")
        fig, ax = mpl_axes(figsize=(8, 4))
        for lake_id, group in frame.groupby("lake", dropna=False):
            labels = [
                f"{int(start)}->{int(end)}"
                for start, end in zip(group["per0"], group["per1"], strict=False)
            ]
            ax.plot(
                labels,
                group["stage_change"].astype(float).to_numpy(),
                marker="o",
                linewidth=2.0,
                label=f"Lake {int(lake_id)}",
            )
        ax.axhline(0.0, color="black", linewidth=1.0, alpha=0.6)
        ax.set_title(heading)
        ax.set_xlabel("Stress-Period Transition")
        ax.set_ylabel("Stage Change")
        if not frame.empty:
            ax.legend()
        fig.tight_layout()
        return fig


class LakConnectionsExplorer:
    """Explorer for LAK connection geometry and exchange interface area."""

    def __init__(self, model: SimulationBase):
        """Bind the LAK connection-geometry explorer to ``model``."""

        self.model = model

    def get(
        self,
        *,
        lake: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return the normalized LAK connection table."""

        frame = build_lak_connection_table(self.model)
        frame = _filter_normalized_table(frame, per=None, layer=layer, cells=cells)
        if lake is not None and "lake" in frame.columns:
            frame = frame.loc[frame["lake"] == int(lake)].copy()
        return frame.reset_index(drop=True)

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of LAK connection geometry."""

        frame = self.get()
        if frame.empty:
            return pd.DataFrame(
                [
                    {
                        "label": "lak.connections",
                        "records": 0,
                        "lakes": 0,
                        "layers": 0,
                        "cells": 0,
                        "total_connection_area": 0.0,
                    }
                ]
            )
        return pd.DataFrame(
            [
                {
                    "label": "lak.connections",
                    "records": int(len(frame)),
                    "lakes": int(frame["lake"].nunique()),
                    "layers": int(frame["layer"].nunique()),
                    "cells": int(frame["cell"].nunique()),
                    "total_connection_area": float(
                        pd.to_numeric(frame["connection_area"], errors="coerce").sum()
                    ),
                }
            ]
        )

    def map(
        self,
        *,
        lake: int | None = None,
        layer: int = 0,
        value_column: str = "connection_area",
        agg: str = "sum",
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        colorscale: str | None = None,
        backend: str = "plotly",
        **kwargs,
    ):
        """Build a choropleth of LAK connection geometry by cell.

        Parameters
        ----------
        lake
            Optional zero-based lake id filter.
        layer
            Zero-based model layer to render.
        value_column
            Connection field to map. Common choices are ``"connection_area"``
            and ``"connwidth"``.
        agg
            Aggregation passed through to the generic cell-input map builder.
        multiplier
            Optional scalar multiplier applied to the mapped values.
        fill_value
            Fill value for cells without lake connections.
        colorscale
            Optional choropleth colorscale override.
        """

        selected = self.get(lake=lake, layer=layer)
        if value_column not in selected.columns:
            raise KeyError(f"LAK connection column {value_column!r} was not found.")
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(self.model))
        values, hover = build_cell_input_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            value_column=value_column,
            per=None,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        kwargs.setdefault("hover_spec", cell_input_hover(value_column))
        choro = self.model.cor(
            per=0,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or "earth",
            **kwargs,
        )
        return _apply_backend(choro, backend)


class SfrStageResultsExplorer(StageResultsExplorer):
    """SFR stage explorer with reach-profile helpers."""

    def __init__(self, model: SimulationBase):
        """Bind an SFR stage-result explorer using the reach stage table builder."""

        super().__init__(model, "sfr", build_sfr_stage_result_table)

    def profile(self, *, per: int = 0) -> pd.DataFrame:
        """Return one reach-ordered stage profile for the selected period."""

        frame = self.get(per=per)
        if frame.empty:
            return frame
        return frame.sort_values(["reach", "cell"]).reset_index(drop=True)

    def plot_profile(
        self,
        *,
        per: int = 0,
        x: str = "distance",
        plot_fig: bool = False,
        return_fig: bool = True,
    ):
        """Plot one SFR stage profile by reach or cumulative stream distance."""

        frame = self.profile(per=per)
        fig = figs.Fig()
        if frame.empty:
            if plot_fig:
                fig.show()
            if return_fig:
                return fig
            return None

        x_column = "distance_mid" if str(x).lower() == "distance" else "reach"
        if x_column not in frame.columns:
            raise KeyError(f"Profile x-axis column {x_column!r} was not found.")
        fig.add_scattergl(
            x=frame[x_column],
            y=frame["stage"],
            mode="lines+markers",
            name=f"SFR stage per {per}",
            customdata=np.column_stack([frame["reach"], frame["cell"]]),
            hovertemplate=(
                "reach=%{customdata[0]}<br>"
                "cell=%{customdata[1]}<br>"
                "stage=%{y}<extra></extra>"
            ),
        )
        fig.update_layout(
            xaxis_title="Stream Distance" if x_column == "distance_mid" else "Reach",
            yaxis_title="stage",
            title=f"SFR stage profile (per={per})",
        )
        if plot_fig:
            fig.show()
        if return_fig:
            return fig
        return None


class LakResultsNamespace(FieldMappable):
    """Namespace for LAK result explorers.

    Mappable fields: ``q`` (lake-groundwater exchange, the default) and
    ``stage`` -- use ``results.map(field="stage")`` or ``results.stage.map()``.
    ``stage_change`` stays a first-class accessor (``results.stage_change.plot()``)
    but is a per-transition series, not a spatial field, so it is not in
    ``field_names()``.
    """

    _default_field = "q"

    def _field_names(self):
        """The mappable LAK result fields: exchange ``q`` and ``stage``."""

        return ["q", "stage"]

    def __init__(self, model: SimulationBase):
        """Bind the LAK results namespace to ``model``."""

        self.model = model

    @property
    def stage(self) -> LakStageResultsExplorer:
        """Return the lake stage explorer mapped to connected cells."""

        return LakStageResultsExplorer(self.model)

    @property
    def stage_change(self) -> LakStageChangeExplorer:
        """Return the lake-stage change explorer."""

        return LakStageChangeExplorer(self.model)

    @property
    def q(self) -> LakBudgetResultsExplorer:
        """Return the lake-groundwater exchange result explorer.

        Exchange maps are normalized to lake connection area, so
        ``map(...)`` renders ``sum(q) / sum(flow_area)`` by cell.
        """

        budget_text, value_name = get_default_budget_term("lak") or ("GWF", "q_lake")
        return LakBudgetResultsExplorer(
            self.model, budget_text=budget_text, value_name=value_name
        )


class LakBudgetNamespace:
    """Namespace for all MF6-defined LAK package-output budget terms."""

    def __init__(self, model: SimulationBase):
        """Bind the LAK budget-term namespace to ``model``."""

        self.model = model

    @property
    def types(self) -> list[str]:
        """Return the LAK package-output budget term names available for the model."""

        return [
            str(value).strip().upper() for value in self.model.outputs.lak.bud.types
        ]

    def get(
        self,
        *,
        term: str | Iterable[str] | None = None,
        per: int | Iterable[int] | None = None,
        lakes: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return a canonical long dataframe of LAK package-output budget terms.

        Parameters
        ----------
        term
            Optional LAK budget term filter such as ``"GWF"`` or
            ``["GWF", "STORAGE"]``.
        per
            Optional zero-based stress period filter.
        lakes
            Optional iterable of zero-based lake ids to keep.
        """

        frame = build_lak_budget_term_table(self.model, term=term)
        if per is not None:
            if isinstance(per, Iterable) and not isinstance(per, (str, bytes)):
                periods = {int(value) for value in per}
                frame = frame.loc[frame["per"].isin(periods)].copy()
            else:
                frame = frame.loc[frame["per"] == int(per)].copy()
        if lakes is not None and "lake" in frame.columns:
            if isinstance(lakes, Iterable) and not isinstance(lakes, (str, bytes)):
                lake_ids = {int(value) for value in lakes}
            else:
                lake_ids = {int(lakes)}
            frame = frame.loc[
                pd.to_numeric(frame["lake"], errors="coerce").isin(lake_ids)
            ].copy()
        return frame.reset_index(drop=True)

    def summary(
        self,
        *,
        term: str | Iterable[str] | None = None,
        per: int | Iterable[int] | None = None,
        lakes: Iterable[int] | None = None,
        by: list[str] | tuple[str, ...] | None = None,
    ) -> pd.DataFrame:
        """Summarize LAK budget terms by selected grouping columns."""

        frame = self.get(term=term, per=per, lakes=lakes)
        group_columns = list(by) if by is not None else ["per", "lake", "term"]
        if frame.empty:
            columns = [*group_columns, "record_count", "q"]
            return pd.DataFrame(columns=columns)

        agg_map: dict[str, tuple[str, str]] = {
            "record_count": ("q", "size"),
            "q": ("q", "sum"),
        }
        for column in ("FLOW-AREA", "flow_area", "VOLUME"):
            if column in frame.columns:
                agg_map[column] = (column, "sum")
        summary = (
            frame.groupby(group_columns, dropna=False, as_index=False)
            .agg(**agg_map)
            .sort_values(group_columns)
            .reset_index(drop=True)
        )
        if "flow_area" in summary.columns:
            summary["q_per_area"] = np.where(
                pd.to_numeric(summary["flow_area"], errors="coerce") > 0.0,
                pd.to_numeric(summary["q"], errors="coerce")
                / pd.to_numeric(summary["flow_area"], errors="coerce"),
                np.nan,
            )
        return summary

    def wide(
        self,
        *,
        term: str | Iterable[str] | None = None,
        per: int | Iterable[int] | None = None,
        lakes: Iterable[int] | None = None,
        index: list[str] | tuple[str, ...] = ("per", "lake"),
        values: str = "q",
    ) -> pd.DataFrame:
        """Pivot the long LAK budget table to one wide table by term."""

        frame = self.get(term=term, per=per, lakes=lakes)
        if frame.empty:
            return pd.DataFrame()
        if values not in frame.columns:
            raise KeyError(f"LAK budget column {values!r} was not found.")
        wide = frame.pivot_table(
            index=list(index),
            columns="term",
            values=values,
            aggfunc="sum",
        ).sort_index()
        if isinstance(wide.columns, pd.Index):
            wide.columns.name = None
        return wide.reset_index()

    @property
    def gwf(self) -> PackageBudgetTermExplorer:
        """Lake-groundwater exchange term helper."""

        return PackageBudgetTermExplorer(self, term="GWF", label="lak.budget.gwf")

    @property
    def storage(self) -> PackageBudgetTermExplorer:
        """Lake storage term helper."""

        return PackageBudgetTermExplorer(
            self, term="STORAGE", label="lak.budget.storage"
        )

    @property
    def runoff(self) -> PackageBudgetTermExplorer:
        """Lake runoff term helper."""

        return PackageBudgetTermExplorer(self, term="RUNOFF", label="lak.budget.runoff")

    @property
    def rainfall(self) -> PackageBudgetTermExplorer:
        """Lake rainfall term helper."""

        return PackageBudgetTermExplorer(
            self, term="RAINFALL", label="lak.budget.rainfall"
        )

    @property
    def evaporation(self) -> PackageBudgetTermExplorer:
        """Lake evaporation term helper."""

        return PackageBudgetTermExplorer(
            self, term="EVAPORATION", label="lak.budget.evaporation"
        )

    @property
    def withdrawal(self) -> PackageBudgetTermExplorer:
        """Lake withdrawal term helper."""

        return PackageBudgetTermExplorer(
            self, term="WITHDRAWAL", label="lak.budget.withdrawal"
        )

    @property
    def constant(self) -> PackageBudgetTermExplorer:
        """Lake constant-stage balancing flow term helper."""

        return PackageBudgetTermExplorer(
            self, term="CONSTANT", label="lak.budget.constant"
        )

    @property
    def ext_inflow(self) -> PackageBudgetTermExplorer:
        """External inflow term helper."""

        return PackageBudgetTermExplorer(
            self, term="EXT-INFLOW", label="lak.budget.ext_inflow"
        )

    @property
    def ext_outflow(self) -> PackageBudgetTermExplorer:
        """External outflow term helper."""

        return PackageBudgetTermExplorer(
            self, term="EXT-OUTFLOW", label="lak.budget.ext_outflow"
        )

    @property
    def from_mvr(self) -> PackageBudgetTermExplorer:
        """Mover inflow term helper."""

        return PackageBudgetTermExplorer(
            self, term="FROM-MVR", label="lak.budget.from_mvr"
        )

    @property
    def to_mvr(self) -> PackageBudgetTermExplorer:
        """Mover outflow term helper."""

        return PackageBudgetTermExplorer(self, term="TO-MVR", label="lak.budget.to_mvr")

    @property
    def flow_ja_face(self) -> PackageBudgetTermExplorer:
        """Lake-to-lake outlet/routing connection term helper."""

        return PackageBudgetTermExplorer(
            self, term="FLOW-JA-FACE", label="lak.budget.flow_ja_face"
        )

    @property
    def auxiliary(self) -> PackageBudgetTermExplorer:
        """Auxiliary term helper."""

        return PackageBudgetTermExplorer(
            self, term="AUXILIARY", label="lak.budget.auxiliary"
        )

    @property
    def mvr(self) -> PackageBudgetTermExplorer:
        """Combined mover-related LAK budget term helper."""

        return PackageBudgetTermExplorer(
            self, term=["FROM-MVR", "TO-MVR"], label="lak.budget.mvr"
        )

    @property
    def lake_fluxes(self) -> PackageBudgetTermExplorer:
        """Combined lake-level flux term helper excluding connection-level GWF rows."""

        return PackageBudgetTermExplorer(
            self,
            term=[
                "EXT-INFLOW",
                "RUNOFF",
                "RAINFALL",
                "EVAPORATION",
                "WITHDRAWAL",
                "STORAGE",
                "CONSTANT",
                "EXT-OUTFLOW",
                "FROM-MVR",
                "TO-MVR",
            ],
            label="lak.budget.lake_fluxes",
        )


class SfrBudgetNamespace:
    """Namespace for all MF6-defined SFR package-output budget terms."""

    def __init__(self, model: SimulationBase):
        """Bind the SFR budget-term namespace to ``model``."""

        self.model = model

    @property
    def types(self) -> list[str]:
        """Return the SFR package-output budget term names available for the model."""

        return [
            str(value).strip().upper() for value in self.model.outputs.sfr.bud.types
        ]

    def get(
        self,
        *,
        term: str | Iterable[str] | None = None,
        per: int | Iterable[int] | None = None,
        reaches: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return a canonical long dataframe of SFR package-output budget terms."""

        frame = build_sfr_budget_term_table(self.model, term=term)
        if per is not None:
            if isinstance(per, Iterable) and not isinstance(per, (str, bytes)):
                periods = {int(value) for value in per}
                frame = frame.loc[frame["per"].isin(periods)].copy()
            else:
                frame = frame.loc[frame["per"] == int(per)].copy()
        if reaches is not None and "reach" in frame.columns:
            if isinstance(reaches, Iterable) and not isinstance(reaches, (str, bytes)):
                reach_ids = {int(value) for value in reaches}
            else:
                reach_ids = {int(reaches)}
            frame = frame.loc[
                pd.to_numeric(frame["reach"], errors="coerce").isin(reach_ids)
            ].copy()
        return frame.reset_index(drop=True)

    def summary(
        self,
        *,
        term: str | Iterable[str] | None = None,
        per: int | Iterable[int] | None = None,
        reaches: Iterable[int] | None = None,
        by: list[str] | tuple[str, ...] | None = None,
    ) -> pd.DataFrame:
        """Summarize SFR budget terms by selected grouping columns."""

        frame = self.get(term=term, per=per, reaches=reaches)
        group_columns = list(by) if by is not None else ["per", "reach", "term"]
        if frame.empty:
            columns = [*group_columns, "record_count", "q"]
            return pd.DataFrame(columns=columns)

        agg_map: dict[str, tuple[str, str]] = {
            "record_count": ("q", "size"),
            "q": ("q", "sum"),
        }
        if "rlen" in frame.columns and "reach" in group_columns:
            agg_map["rlen"] = ("rlen", "first")
        for column in ("FLOW-AREA", "VOLUME"):
            if column in frame.columns:
                agg_map[column] = (column, "sum")
        summary = (
            frame.groupby(group_columns, dropna=False, as_index=False)
            .agg(**agg_map)
            .sort_values(group_columns)
            .reset_index(drop=True)
        )
        if "rlen" in summary.columns:
            summary["q_per_length"] = np.where(
                pd.to_numeric(summary["rlen"], errors="coerce") > 0.0,
                pd.to_numeric(summary["q"], errors="coerce")
                / pd.to_numeric(summary["rlen"], errors="coerce"),
                np.nan,
            )
        return summary

    def wide(
        self,
        *,
        term: str | Iterable[str] | None = None,
        per: int | Iterable[int] | None = None,
        reaches: Iterable[int] | None = None,
        index: list[str] | tuple[str, ...] = ("per", "reach"),
        values: str = "q",
    ) -> pd.DataFrame:
        """Pivot the long SFR budget table to one wide table by term."""

        frame = self.get(term=term, per=per, reaches=reaches)
        if frame.empty:
            return pd.DataFrame()
        if values not in frame.columns:
            raise KeyError(f"SFR budget column {values!r} was not found.")
        wide = frame.pivot_table(
            index=list(index),
            columns="term",
            values=values,
            aggfunc="sum",
        ).sort_index()
        if isinstance(wide.columns, pd.Index):
            wide.columns.name = None
        return wide.reset_index()

    @property
    def gwf(self) -> PackageBudgetTermExplorer:
        """Stream-groundwater exchange term helper."""

        return PackageBudgetTermExplorer(self, term="GWF", label="sfr.budget.gwf")

    @property
    def flow_ja_face(self) -> PackageBudgetTermExplorer:
        """Reach-to-reach routing connection term helper."""

        return PackageBudgetTermExplorer(
            self, term="FLOW-JA-FACE", label="sfr.budget.flow_ja_face"
        )

    @property
    def ext_inflow(self) -> PackageBudgetTermExplorer:
        """External inflow term helper."""

        return PackageBudgetTermExplorer(
            self, term="EXT-INFLOW", label="sfr.budget.ext_inflow"
        )

    @property
    def runoff(self) -> PackageBudgetTermExplorer:
        """Runoff term helper."""

        return PackageBudgetTermExplorer(self, term="RUNOFF", label="sfr.budget.runoff")

    @property
    def rain(self) -> PackageBudgetTermExplorer:
        """Rainfall term helper."""

        return PackageBudgetTermExplorer(self, term="RAIN", label="sfr.budget.rain")

    @property
    def evaporation(self) -> PackageBudgetTermExplorer:
        """Evaporation term helper."""

        return PackageBudgetTermExplorer(
            self, term="EVAPORATION", label="sfr.budget.evaporation"
        )

    @property
    def ext_outflow(self) -> PackageBudgetTermExplorer:
        """External outflow term helper."""

        return PackageBudgetTermExplorer(
            self, term="EXT-OUTFLOW", label="sfr.budget.ext_outflow"
        )

    @property
    def storage(self) -> PackageBudgetTermExplorer:
        """Storage term helper."""

        return PackageBudgetTermExplorer(
            self, term="STORAGE", label="sfr.budget.storage"
        )

    @property
    def from_mvr(self) -> PackageBudgetTermExplorer:
        """Mover inflow term helper."""

        return PackageBudgetTermExplorer(
            self, term="FROM-MVR", label="sfr.budget.from_mvr"
        )

    @property
    def to_mvr(self) -> PackageBudgetTermExplorer:
        """Mover outflow term helper."""

        return PackageBudgetTermExplorer(self, term="TO-MVR", label="sfr.budget.to_mvr")

    @property
    def auxiliary(self) -> PackageBudgetTermExplorer:
        """Auxiliary term helper."""

        return PackageBudgetTermExplorer(
            self, term="AUXILIARY", label="sfr.budget.auxiliary"
        )

    @property
    def mvr(self) -> PackageBudgetTermExplorer:
        """Combined mover-related SFR budget term helper."""

        return PackageBudgetTermExplorer(
            self, term=["FROM-MVR", "TO-MVR"], label="sfr.budget.mvr"
        )

    @property
    def stream_fluxes(self) -> PackageBudgetTermExplorer:
        """Combined reach-level flux term helper excluding GWF and routing rows."""

        return PackageBudgetTermExplorer(
            self,
            term=[
                "EXT-INFLOW",
                "RUNOFF",
                "RAIN",
                "EVAPORATION",
                "EXT-OUTFLOW",
                "STORAGE",
                "FROM-MVR",
                "TO-MVR",
            ],
            label="sfr.budget.stream_fluxes",
        )


def _exchange_frame(package: str) -> str:
    """The declared reference frame for a package's exchange ``q``.

    Reads ``ResultSpec.reference_frame`` from the registry so a map's colour
    orientation is derived from the same field that documents the sign, never a
    literal at the call site that could fall out of step with the data.
    """

    spec = get_package_result_spec(package, "q")
    return spec.reference_frame if spec is not None else "gwf"


def _signed_exchange_colors() -> tuple[str, str]:
    """Return the ``(gaining, losing)`` bar colors for signed exchange plots.

    Both are read off ``_blue_white_red_diverging_colorscale`` rather than
    written out, so discrete bars and continuous maps cannot drift apart. The
    colours are frame-independent (gaining is always blue, losing always red);
    which *sign* of ``q`` is gaining is a per-frame fact the caller applies. The
    signed SFR profile bars are in the "gwf" frame, where negative ``q`` is the
    reach gaining (BLUE) and positive ``q`` is the reach losing (RED).
    """

    scale = _blue_white_red_diverging_colorscale()
    return str(scale[0][1]), str(scale[-1][1])


class SfrProfileView:
    """The SFR long profile: one merged reach-ordered table and its figure.

    Follows the house view shape (``docs/view_layer_conventions.md``) -- a noun
    reached from the results namespace, carrying ``get`` for the frame, ``plot``
    for the figure, and ``summary`` for the compact digest. Calling the view
    rebinds the stress period, so these are the same figure::

        model.packages.sfr.results.profile.plot(per=3)
        model.packages.sfr.results.profile(per=3).plot()

    The table merges reach geometry, streambed elevations, simulated stage, and
    stream-groundwater exchange, so it feeds custom analysis as readily as it
    feeds ``plot``.
    """

    def __init__(self, model: SimulationBase, *, per: int = 0):
        """Bind the profile view to ``model`` at stress period ``per``."""

        self.model = model
        self.per = int(per)

    def __repr__(self) -> str:
        """Show the bound period, since it is the view's only state."""

        return f"{type(self).__name__}(per={self.per})"

    def __call__(self, *, per: int) -> SfrProfileView:
        """Return an equivalent view bound to stress period ``per``."""

        return type(self)(self.model, per=per)

    def _period(self, per: int | None) -> int:
        """Resolve an explicit ``per`` against the period bound to this view."""

        return self.per if per is None else int(per)

    def get(self, *, per: int | None = None) -> pd.DataFrame:
        """Return the merged reach-ordered profile table."""

        return build_sfr_long_profile_table(self.model, per=self._period(per))

    def summary(self, *, per: int | None = None) -> pd.DataFrame:
        """Return a compact digest of the profile's stage and exchange fields."""

        frame = self.get(per=per)
        return summarize_input_table(
            frame,
            label="sfr.results.profile",
            value_columns=[
                column
                for column in ("stage", "streambed_top", "q_gwf", "q_per_length")
                if column in frame.columns
            ],
        )

    def plot(
        self,
        *,
        per: int | None = None,
        x: str = "distance",
        include_stage: bool = True,
        include_streambed: bool = True,
        include_exchange: bool = True,
        signed_exchange: bool = True,
        plot_fig: bool = False,
        return_fig: bool = True,
    ):
        """Plot the long profile with the common hydrologic overlays.

        Parameters
        ----------
        per
            Zero-based stress period; defaults to the period bound to the view.
        x
            Either ``"distance"`` for cumulative stream distance or ``"reach"``
            for raw reach number.
        include_stage
            Whether to show simulated stream stage.
        include_streambed
            Whether to show streambed top and bottom elevations.
        include_exchange
            Whether to show stream-groundwater exchange on a secondary axis.
        signed_exchange
            When ``True`` (the default) draw exchange as per-reach bars colored
            by sign -- blue where the reach gains (NEGATIVE q_gwf: SFR's cell
            record is aquifer-referenced, MF6's raw sign), red where it loses,
            matching the SFR map colorscale. When ``False`` draw one unsigned line.
        plot_fig
            If ``True``, call ``show()`` on the created figure.
        return_fig
            If ``True``, return the created figure.
        """

        period = self._period(per)
        frame = self.get(per=period)
        fig = figs.Fig()
        if frame.empty:
            if plot_fig:
                fig.show()
            if return_fig:
                return fig
            return None

        x_column = "distance_mid" if str(x).lower() == "distance" else "reach"
        if x_column not in frame.columns:
            raise KeyError(f"Long-profile x-axis column {x_column!r} was not found.")

        customdata = np.column_stack([frame["reach"], frame["cell"]])
        if include_streambed and "streambed_top" in frame.columns:
            fig.add_scattergl(
                x=frame[x_column],
                y=frame["streambed_top"],
                mode="lines",
                name="Streambed Top",
                line={"color": "#8c564b", "width": 2},
                customdata=customdata,
                hovertemplate=(
                    "reach=%{customdata[0]}<br>"
                    "cell=%{customdata[1]}<br>"
                    "streambed_top=%{y}<extra></extra>"
                ),
            )
        if include_streambed and "streambed_bottom" in frame.columns:
            fig.add_scattergl(
                x=frame[x_column],
                y=frame["streambed_bottom"],
                mode="lines",
                name="Streambed Bottom",
                line={"color": "#c49c94", "width": 2, "dash": "dash"},
                customdata=customdata,
                hovertemplate=(
                    "reach=%{customdata[0]}<br>"
                    "cell=%{customdata[1]}<br>"
                    "streambed_bottom=%{y}<extra></extra>"
                ),
            )
        if include_stage and "stage" in frame.columns:
            fig.add_scattergl(
                x=frame[x_column],
                y=frame["stage"],
                mode="lines+markers",
                name="Stage",
                line={"color": "#1f77b4", "width": 3},
                marker={"size": 7},
                customdata=customdata,
                hovertemplate=(
                    "reach=%{customdata[0]}<br>"
                    "cell=%{customdata[1]}<br>"
                    "stage=%{y}<extra></extra>"
                ),
            )
        if include_exchange and "q_gwf" in frame.columns:
            self._add_exchange_trace(
                fig,
                frame=frame,
                x_column=x_column,
                customdata=customdata,
                signed=signed_exchange,
            )
            fig.update_layout(
                yaxis2={
                    "title": "Exchange q",
                    "overlaying": "y",
                    "side": "right",
                    "showgrid": False,
                    "zeroline": True,
                }
            )

        fig.update_layout(
            xaxis_title="Stream Distance" if x_column == "distance_mid" else "Reach",
            yaxis_title="Elevation / Stage",
            title=f"SFR long profile (per={period})",
        )
        if plot_fig:
            fig.show()
        if return_fig:
            return fig
        return None

    @staticmethod
    def _add_exchange_trace(fig, *, frame, x_column, customdata, signed: bool) -> None:
        """Add the exchange trace, as signed bars or as one unsigned line."""

        hovertemplate = (
            "reach=%{customdata[0]}<br>"
            "cell=%{customdata[1]}<br>"
            "q=%{y}<extra></extra>"
        )
        if not signed:
            fig.add_scattergl(
                x=frame[x_column],
                y=frame["q_gwf"],
                mode="lines+markers",
                name="Exchange q",
                line={"color": "#d62728", "width": 2},
                marker={"size": 6, "symbol": "diamond"},
                yaxis="y2",
                customdata=customdata,
                hovertemplate=hovertemplate,
            )
            return

        gaining_color, losing_color = _signed_exchange_colors()
        q = pd.to_numeric(frame["q_gwf"], errors="coerce").to_numpy(float)
        positions = pd.to_numeric(frame[x_column], errors="coerce").to_numpy(float)
        # One bar per reach, sized just under the reach spacing so adjacent
        # bars read as separate reaches. A single reach has no spacing to
        # measure, so fall back to plotly's own default width.
        spacing = np.diff(np.sort(positions[np.isfinite(positions)]))
        spacing = spacing[spacing > 0.0]
        width = float(np.median(spacing)) * 0.9 if spacing.size else None
        fig.add_bar(
            x=frame[x_column],
            y=frame["q_gwf"],
            name="Exchange q (blue gains, red loses)",
            marker={
                # SFR's q keeps MF6's raw sign (the "gwf" frame): NEGATIVE q is
                # the reach gaining -> blue; positive q is losing -> red.
                "color": np.where(q < 0.0, gaining_color, losing_color).tolist(),
                "line": {"width": 0},
            },
            opacity=0.55,
            width=width,
            yaxis="y2",
            customdata=customdata,
            hovertemplate=hovertemplate,
        )


class SfrResultsNamespace(FieldMappable):
    """Namespace for SFR result explorers.

    Fields: ``q`` (stream-groundwater exchange, the default) and ``stage``.
    Use ``results.map(field="stage")`` or ``results.stage.map()``.
    """

    _default_field = "q"

    def _field_names(self):
        """The mappable SFR result fields: exchange ``q`` and ``stage``."""

        return ["q", "stage"]

    def __init__(self, model: SimulationBase):
        """Bind the SFR results namespace to ``model``."""

        self.model = model

    @property
    def stage(self) -> StageResultsExplorer:
        """Return the stream stage explorer mapped to reach cells."""

        return SfrStageResultsExplorer(self.model)

    @property
    def q(self) -> SfrBudgetResultsExplorer:
        """Return the stream-groundwater exchange result explorer."""

        budget_text, value_name = get_default_budget_term("sfr") or ("SFR", "q_gwf")
        return SfrBudgetResultsExplorer(
            self.model, budget_text=budget_text, value_name=value_name
        )

    @property
    def profile(self) -> SfrProfileView:
        """Return the long-profile view: ``.get()`` for the table, ``.plot()``.

        Named ``profile`` to match ``results.stage.profile`` one level down --
        inside an SFR namespace a profile can only be longitudinal, so "long"
        carried no information. This namespace-level view merges every field;
        the field-level ones cover a single field.
        """

        return SfrProfileView(self.model)

    # -- backing methods for the pre-view spellings ---------------------------
    # ``profile`` replaced two older names when the derived tables became view
    # objects. These bodies preserve the OLD return values exactly (a DataFrame
    # and an unsigned-line figure); the mapping below is the only place the old
    # spellings appear, and they resolve solely through __getattr__ so they stay
    # out of dir()/completion (D12).
    #
    # Deliberately named ``_legacy_*`` rather than echoing the old spelling: a
    # private member called ``_long_profile_frame`` would still surface the
    # retired name in IDE completion, which is exactly what D12 exists to stop.
    def _legacy_profile_frame(self, *, per: int = 0) -> pd.DataFrame:
        """Back the retired frame spelling; returns what it always returned."""

        return self.profile.get(per=per)

    def _legacy_profile_plot(self, **kwargs):
        """Back the retired plot spelling.

        Pins ``signed_exchange=False`` so it keeps drawing the single unsigned
        line it always drew; signed bars are the new ``profile.plot()`` default.
        """

        kwargs.setdefault("signed_exchange", False)
        return self.profile.plot(**kwargs)

    __getattr__ = deprecated_instance_getattr(
        {
            "long_profile": (
                "_legacy_profile_frame",
                "model.packages.sfr.results.profile.get",
                "0.1.0",
            ),
            "plot_long_profile": (
                "_legacy_profile_plot",
                "model.packages.sfr.results.profile.plot",
                "0.1.0",
            ),
        },
        "myflopy.modflow.mf6.package_surface_water.SfrResultsNamespace",
    )


class SurfaceWaterExchangeResultsExplorer(SpatialView):
    """Combined SFR/LAK exchange explorer with one shared physical sign scale."""

    def __init__(self, model: SimulationBase):
        """Bind the combined SFR/LAK exchange explorer to ``model``."""

        self.model = model

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        include: str | Iterable[str] | None = None,
        lak_connection_type: str | Iterable[str] | None = None,
    ) -> pd.DataFrame:
        """Return combined SFR/LAK exchange rows in a unified L/T convention."""

        return build_surface_water_exchange_cell_table(
            self.model,
            per=per,
            layer=layer,
            include=include,
            lak_connection_type=lak_connection_type,
        )

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of combined surface-water exchange rows."""

        frame = self.get()
        return summarize_input_table(
            frame,
            label="surface_water.results.q",
            value_columns=["exchange_intensity"],
        )

    def map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        include: str | Iterable[str] | None = None,
        lak_connection_type: str | Iterable[str] | None = None,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        colorscale: str | None = None,
        backend: str = "plotly",
        **kwargs,
    ):
        """Build one combined SFR/LAK exchange map with a shared L/T scale.

        Notes
        -----
        The mapped value uses a unified physical sign convention across SFR and
        LAK:

        - positive = groundwater gaining into the surface-water feature
        - negative = surface-water losing to groundwater
        """

        selected = self.get(
            per=per,
            layer=layer,
            include=include,
            lak_connection_type=lak_connection_type,
        )
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(self.model))
        values, hover = build_surface_water_q_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
        )
        absmax = _symmetric_color_limit(values)
        kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
        kwargs.setdefault("zmax", absmax if absmax > 0 else None)
        kwargs.setdefault("zmid", 0.0)
        kwargs.setdefault("hover_spec", surface_water_hover())
        choro = self.model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            # This map draws ``exchange_intensity``, myflopy's OWN unified field,
            # which is normalized so POSITIVE = the feature gains (see
            # build_surface_water_exchange_cell_table). That matches the
            # "feature" orientation -- blue at the positive end -- regardless of
            # each source package's own raw frame.
            colorscale=colorscale or _exchange_colorscale("feature"),
            **kwargs,
        )
        return _apply_backend(choro, backend)


class SurfaceWaterResultsNamespace(FieldMappable):
    """Namespace for combined surface-water result explorers (field: ``q``)."""

    _default_field = "q"

    def _field_names(self):
        """The single mappable combined surface-water field: exchange ``q``."""

        return ["q"]

    def __init__(self, model: SimulationBase):
        """Bind the combined surface-water results namespace to ``model``."""

        self.model = model

    @property
    def q(self) -> SurfaceWaterExchangeResultsExplorer:
        """Return one shared SFR/LAK exchange explorer."""

        return SurfaceWaterExchangeResultsExplorer(self.model)


class SurfaceWaterInputFieldExplorer(SpatialView):
    """One cell-mapped LAK or SFR input field."""

    def __init__(self, inputs: SurfaceWaterInputsNamespace, field_name: str):
        """Pin a LAK/SFR inputs namespace to one mappable numeric input field."""

        self.inputs = inputs
        self.model = inputs.model
        self.package_name = inputs.package_name
        self.field_name = str(field_name)

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return rows containing this mapped input field."""

        frame = self.inputs.get(per=per, layer=layer, cells=cells)
        metadata = [
            column
            for column in ("model", "package", "per", "lake", "reach", "layer", "cell")
            if column in frame.columns
        ]
        return frame.loc[:, [*metadata, self.field_name]].copy()

    def summary(self) -> pd.DataFrame:
        """Return a compact one-field summary of this LAK/SFR input field."""

        return summarize_input_table(
            self.get(),
            label=f"{self.package_name}.inputs.{self.field_name}",
            value_columns=[self.field_name],
        )

    def map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str | None = None,
        colorscale: str | None = None,
        backend: str = "plotly",
        **kwargs,
    ):
        """Build a Plotly choropleth for this LAK or SFR input field."""

        selected = self.get(per=per, layer=layer)
        if agg is None:
            agg = (
                "sum"
                if self.field_name in {"connection_area", "rlen", "inflow", "runoff"}
                else "first"
            )
        values, hover = build_cell_input_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            value_column=self.field_name,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(self.model))
        kwargs.setdefault("hover_spec", cell_input_hover(self.field_name))
        choro = self.model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or "earth",
            **kwargs,
        )
        return _apply_backend(choro, backend)


class SurfaceWaterInputsNamespace(FieldMappable):
    """Consistent input exploration namespace for LAK and SFR.

    Every numeric input column is a first-class field node (``sfr.inputs.rhk``,
    ``lak.inputs.connection_area``) with the unified grammar; the namespace
    verbs take ``field=`` as sugar over them.
    """

    def __init__(self, model: SimulationBase, package_name: str):
        """Bind a LAK or SFR inputs namespace to ``model`` (``package_name`` lowercased)."""

        self.model = model
        self.package_name = str(package_name).lower()

    @property
    def _default_field(self) -> str:
        """Preferred input field: LAK ``connection_area``, SFR ``rhk``."""

        return "connection_area" if self.package_name == "lak" else "rhk"

    def _field_names(self) -> list[str]:
        """The mappable numeric input field names discovered for this package."""

        return self.fields["field"].tolist()

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return normalized package inputs mapped to groundwater cells."""

        builder = (
            build_lak_input_table
            if self.package_name == "lak"
            else build_sfr_input_table
        )
        return _filter_normalized_table(
            builder(self.model), per=per, layer=layer, cells=cells
        )

    @property
    def fields(self) -> pd.DataFrame:
        """Return numeric input fields that can be mapped."""

        frame = self.get()
        excluded = {
            "model",
            "package",
            "per",
            "lake",
            "reach",
            "layer",
            "cell",
            "ifno",
            "iconn",
            "ncon",
            "ndv",
            "nlakeconn",
        }
        fields = [
            column
            for column in frame.columns
            if column not in excluded and pd.api.types.is_numeric_dtype(frame[column])
        ]
        return pd.DataFrame({"field": fields})

    def summary(self) -> pd.DataFrame:
        """Return one stacked summary row per mappable input field for this package."""

        frames = [
            getattr(self, field).summary() for field in self.fields["field"].tolist()
        ]
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    def __getattr__(self, field_name: str) -> SurfaceWaterInputFieldExplorer:
        """Resolve ``inputs.<field>`` to a field-pinned explorer (else ``AttributeError``)."""

        if field_name not in set(self.fields["field"].tolist()):
            raise AttributeError(
                f"{type(self).__name__!s} has no input field {field_name!r}"
            )
        return SurfaceWaterInputFieldExplorer(self, field_name)

    @property
    def default(self) -> SurfaceWaterInputFieldExplorer:
        """Return the preferred package input field."""

        return getattr(self, self._default_field)

    # map/plot/xs/mosaic/animate come from FieldMappable (field= sugar).


class LakPackageExplorer:
    """Top-level LAK package explorer namespace."""

    def __init__(self, model: SimulationBase):
        """Bind the top-level LAK explorer (inputs/connections/budget/results) to ``model``."""

        self.model = model

    @property
    def inputs(self) -> SurfaceWaterInputsNamespace:
        """Return mapped LAK starting-stage, connection, and period inputs."""

        return SurfaceWaterInputsNamespace(self.model, "lak")

    @property
    def connections(self) -> LakConnectionsExplorer:
        """Return LAK connection-geometry exploration helpers."""

        return LakConnectionsExplorer(self.model)

    @property
    def budget(self) -> LakBudgetNamespace:
        """Return LAK package-output budget helpers for all MF6-defined terms."""

        return LakBudgetNamespace(self.model)

    @property
    def results(self) -> LakResultsNamespace:
        """Return the LAK result exploration namespace."""

        return LakResultsNamespace(self.model)


class SfrPackageExplorer:
    """Top-level SFR package explorer namespace."""

    def __init__(self, model: SimulationBase):
        """Bind the top-level SFR explorer (inputs/budget/results) to ``model``."""

        self.model = model

    @property
    def inputs(self) -> SurfaceWaterInputsNamespace:
        """Return mapped SFR reach-hydraulic and period inputs."""

        return SurfaceWaterInputsNamespace(self.model, "sfr")

    @property
    def budget(self) -> SfrBudgetNamespace:
        """Return SFR package-output budget helpers for all MF6-defined terms."""

        return SfrBudgetNamespace(self.model)

    @property
    def results(self) -> SfrResultsNamespace:
        """Return the SFR result exploration namespace."""

        return SfrResultsNamespace(self.model)


class SurfaceWaterPackageExplorer:
    """Top-level combined surface-water explorer namespace."""

    def __init__(self, model: SimulationBase):
        """Bind the top-level combined SFR/LAK explorer to ``model``."""

        self.model = model

    @property
    def results(self) -> SurfaceWaterResultsNamespace:
        """Return combined SFR/LAK result explorers."""

        return SurfaceWaterResultsNamespace(self.model)


__all__ = [
    "SfrBudgetResultsExplorer",
    "LakBudgetResultsExplorer",
    "LakStageResultsExplorer",
    "LakStageChangeExplorer",
    "LakConnectionsExplorer",
    "SfrStageResultsExplorer",
    "SfrProfileView",
    "LakResultsNamespace",
    "LakBudgetNamespace",
    "SfrBudgetNamespace",
    "SfrResultsNamespace",
    "SurfaceWaterExchangeResultsExplorer",
    "SurfaceWaterResultsNamespace",
    "SurfaceWaterInputFieldExplorer",
    "SurfaceWaterInputsNamespace",
    "LakPackageExplorer",
    "SfrPackageExplorer",
    "SurfaceWaterPackageExplorer",
]

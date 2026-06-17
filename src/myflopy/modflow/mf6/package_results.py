"""Generic result explorer classes behind model.packages."""

from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.modflow.mf6.package_registry import (
    get_default_package_colorscale,
    get_package_explorer_spec,
    get_package_result_spec,
)
from myflopy.modflow.mf6.package_explorer_utils import (
    _default_show_layer_elevs,
    _filter_normalized_table,
    _normalize_iterable_filter,
    _normalize_term_filter,
)
from myflopy.modflow.mf6.package_tables import (
    summarize_input_table,
)
from myflopy.modflow.mf6.package_budget import (
    build_budget_result_table,
)
from myflopy.modflow.mf6.package_plotting import (
    MappedFieldVisualizationMixin,
    _symmetric_color_limit,
    build_cell_input_map_payload,
)


class CellBudgetResultsExplorer(MappedFieldVisualizationMixin):
    """Normalized explorer for one cell-based package result term."""

    def __init__(
        self,
        model: "SimulationBase",
        package_name: str,
        budget_text: str,
        value_name: str,
    ):
        self.model = model
        self.package_name = str(package_name).lower()
        self.budget_text = str(budget_text)
        self.value_name = str(value_name)

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return the normalized result table for this budget term."""

        frame = build_budget_result_table(
            self.model,
            budget_text=self.budget_text,
            package_name=self.package_name,
            value_name=self.value_name,
        )
        return _filter_normalized_table(frame, per=per, layer=layer, cells=cells)

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of the available result rows."""

        frame = self.get()
        return summarize_input_table(
            frame,
            label=f"{self.package_name}.results.{self.value_name}",
            value_columns=[self.value_name],
        )

    def wide(
        self,
        *,
        per: int | Iterable[int] | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
        index: list[str] | tuple[str, ...] = ("layer", "cell"),
        values: str | None = None,
        agg: str = "sum",
    ) -> pd.DataFrame:
        """Pivot this result term to one column per stress period."""

        value_column = self.value_name if values is None else str(values)
        frame = self.get(layer=layer, cells=cells)
        per_values = _normalize_iterable_filter(per)
        if per_values is not None and "per" in frame.columns:
            frame = frame[frame["per"].isin(per_values)]
        if frame.empty:
            return pd.DataFrame(columns=[*index])
        if value_column not in frame.columns:
            raise KeyError(f"Result value column {value_column!r} was not found.")
        missing_index = [column for column in index if column not in frame.columns]
        if missing_index:
            raise KeyError(f"Wide result index columns were not found: {missing_index}")

        wide = frame.pivot_table(
            index=list(index),
            columns="per",
            values=value_column,
            aggfunc=agg,
        )
        wide.columns = [f"per_{int(column)}" for column in wide.columns]
        return wide.reset_index()

    def long(
        self,
        *,
        per: int | Iterable[int] | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
        values: str | None = None,
        agg: str = "sum",
    ) -> pd.Series:
        """Return a long result series indexed by ``kstpkper/layer/cell``."""

        value_column = self.value_name if values is None else str(values)
        frame = self.get(layer=layer, cells=cells)
        per_values = _normalize_iterable_filter(per)
        if per_values is not None and "per" in frame.columns:
            frame = frame[frame["per"].isin(per_values)]
        index_columns = ["kstpkper", "layer", "cell"]
        if frame.empty:
            empty_index = pd.MultiIndex.from_arrays(
                [[] for _ in index_columns],
                names=index_columns,
            )
            return pd.Series([], index=empty_index, dtype=float, name=value_column)
        if value_column not in frame.columns:
            raise KeyError(f"Result value column {value_column!r} was not found.")
        missing_index = [
            column for column in index_columns if column not in frame.columns
        ]
        if missing_index:
            raise KeyError(f"Long result index columns were not found: {missing_index}")

        series = (
            frame.groupby(index_columns, dropna=False)[value_column]
            .agg(agg)
            .sort_index()
        )
        series.name = value_column
        return series

    def stack(self, **kwargs) -> pd.Series:
        """Alias for :meth:`long`."""

        return self.long(**kwargs)

    def plot_timeseries(
        self,
        *,
        cells: int | Iterable[int] | None = None,
        layer: int | Iterable[int] | None = None,
        agg: str = "sum",
        ax=None,
        return_fig: bool = True,
    ):
        """Plot this cell result by stress period for selected cells."""

        selected_cells = [int(cells)] if isinstance(cells, (int, np.integer)) else cells
        frame = self.get(layer=layer, cells=selected_cells)
        result_spec = get_package_result_spec(self.package_name, self.value_name)
        display_label = (
            result_spec.label
            if result_spec is not None and result_spec.label is not None
            else f"{self.package_name.upper()} {self.value_name}"
        )
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 4))
        else:
            fig = ax.figure
        if frame.empty:
            ax.set_title(f"{display_label} by stress period")
            ax.set_xlabel("Stress Period")
            ax.set_ylabel(self.value_name)
            if return_fig:
                return fig
            return None

        grouped_keys = [
            column for column in ("layer", "cell") if column in frame.columns
        ]
        for key, group in frame.groupby(grouped_keys, dropna=False):
            if not isinstance(key, tuple):
                key = (key,)
            key_map = dict(zip(grouped_keys, key, strict=False))
            series = (
                group.groupby("per", as_index=False)[self.value_name]
                .agg(agg)
                .sort_values("per")
            )
            layer_label = (
                f"L{int(key_map['layer'])} "
                if "layer" in key_map and pd.notna(key_map["layer"])
                else ""
            )
            cell_label = (
                f"C{int(key_map['cell'])}"
                if "cell" in key_map and pd.notna(key_map["cell"])
                else "All cells"
            )
            ax.plot(
                series["per"].astype(int).to_numpy(),
                series[self.value_name].astype(float).to_numpy(),
                marker="o",
                linewidth=2.0,
                label=f"{layer_label}{cell_label}",
            )

        ax.set_title(f"{display_label} by stress period")
        ax.set_xlabel("Stress Period")
        ax.set_ylabel(self.value_name)
        ax.legend()
        fig.tight_layout()
        if return_fig:
            return fig
        return None

    def map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Build a choropleth for this cell-based result field."""

        selected = self.get(per=per, layer=layer)
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(self.model))
        values, hover = build_cell_input_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            value_column=self.value_name,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        if self.value_name == "q":
            absmax = _symmetric_color_limit(values)
            kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
            kwargs.setdefault("zmax", absmax if absmax > 0 else None)
            kwargs.setdefault("zmid", 0.0)
        result_spec = get_package_result_spec(self.package_name, self.value_name)
        return self.model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=(
                colorscale
                or (result_spec.colorscale if result_spec is not None else None)
                or ("RdBu" if self.value_name == "q" else None)
                or get_default_package_colorscale(self.package_name)
                or "Viridis"
            ),
            **kwargs,
        )


class StageResultsExplorer(MappedFieldVisualizationMixin):
    """Normalized explorer for cell-mapped stage results such as LAK and SFR."""

    def __init__(self, model: "SimulationBase", package_name: str, builder):
        self.model = model
        self.package_name = str(package_name).lower()
        self._builder = builder

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return the normalized stage table for the selected rows."""

        frame = self._builder(self.model)
        return _filter_normalized_table(frame, per=per, layer=layer, cells=cells)

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of available stage results."""

        frame = self.get()
        return summarize_input_table(
            frame,
            label=f"{self.package_name}.results.stage",
            value_columns=["stage"],
        )

    def map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "first",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Build a stage choropleth mapped to cells."""

        selected = self.get(per=per, layer=layer)
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(self.model))
        values, hover = build_cell_input_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            value_column="stage",
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        return self.model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or "Blues",
            **kwargs,
        )


class CellPackageResultsNamespace:
    """Namespace for cell-based package results represented by one budget term."""

    def __init__(self, model: "SimulationBase", package_name: str):
        self.model = model
        self.package_name = str(package_name).lower()

    @property
    def fields(self) -> pd.DataFrame:
        """Return registry metadata for supported package result fields."""

        spec = get_package_explorer_spec(self.package_name)
        if spec is None:
            return pd.DataFrame(
                columns=[
                    "field",
                    "budget_text",
                    "value_name",
                    "label",
                    "colorscale",
                    "diverging",
                ]
            )
        rows = [
            {
                "field": result_spec.name,
                "budget_text": result_spec.budget_text,
                "value_name": result_spec.value_name,
                "label": result_spec.label,
                "colorscale": result_spec.colorscale,
                "diverging": result_spec.diverging,
            }
            for result_spec in spec.results.values()
        ]
        return pd.DataFrame(rows)

    def summary(self) -> pd.DataFrame:
        """Return one compact summary row per supported package result field."""

        frames = []
        for field_name in self.fields["field"].tolist():
            frames.append(getattr(self, field_name).summary())
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

    def __getattr__(self, result_name: str) -> CellBudgetResultsExplorer:
        """Return a registry-backed named result explorer."""

        result_spec = get_package_result_spec(self.package_name, result_name)
        if result_spec is None:
            raise AttributeError(
                f"{type(self).__name__!s} has no result {result_name!r}"
            )
        return CellBudgetResultsExplorer(
            self.model,
            self.package_name,
            result_spec.budget_text,
            result_spec.value_name,
        )

    @property
    def q(self) -> CellBudgetResultsExplorer:
        """Return the primary package-exchange result explorer."""

        result_spec = get_package_result_spec(self.package_name, "q")
        if result_spec is None:
            return CellBudgetResultsExplorer(
                self.model, self.package_name, self.package_name.upper(), "q"
            )
        return CellBudgetResultsExplorer(
            self.model,
            self.package_name,
            result_spec.budget_text,
            result_spec.value_name,
        )


class UzfResultsNamespace:
    """Namespace for UZF result explorers."""

    def __init__(self, model: "SimulationBase"):
        self.model = model

    @property
    def fields(self) -> pd.DataFrame:
        """Return registry metadata for supported UZF result fields."""

        spec = get_package_explorer_spec("uzf")
        if spec is None:
            return pd.DataFrame(
                columns=[
                    "field",
                    "budget_text",
                    "value_name",
                    "label",
                    "colorscale",
                    "diverging",
                ]
            )
        rows = [
            {
                "field": result_spec.name,
                "budget_text": result_spec.budget_text,
                "value_name": result_spec.value_name,
                "label": result_spec.label,
                "colorscale": result_spec.colorscale,
                "diverging": result_spec.diverging,
            }
            for result_spec in spec.results.values()
        ]
        return pd.DataFrame(rows)

    def summary(self) -> pd.DataFrame:
        """Return one compact summary row per supported UZF result field."""

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

    def _field(self, field_name: str) -> CellBudgetResultsExplorer:
        """Return one registry-backed UZF result explorer."""

        result_spec = get_package_result_spec("uzf", field_name)
        if result_spec is None:
            raise AttributeError(
                f"{type(self).__name__!s} has no UZF result field {field_name!r}"
            )
        return CellBudgetResultsExplorer(
            self.model, "uzf", result_spec.budget_text, result_spec.value_name
        )

    def __getattr__(self, field_name: str) -> CellBudgetResultsExplorer:
        """Return a registry-backed UZF result explorer."""

        return self._field(field_name)

    @property
    def gwrch(self) -> CellBudgetResultsExplorer:
        """Return groundwater recharge from the UZF package."""

        return self._field("gwrch")

    @property
    def sat(self) -> CellBudgetResultsExplorer:
        """Return normalized unsaturated-zone saturation results."""

        return self._field("sat")


class PackageBudgetTermExplorer:
    """Filtered helper for one package budget term or a small term family."""

    def __init__(
        self,
        namespace,
        *,
        term: str | Iterable[str],
        label: str,
    ):
        self._namespace = namespace
        self.term = term
        self.label = label

    @property
    def types(self) -> list[str]:
        """Return the normalized MF6 LAK term names covered by this helper."""

        return _normalize_term_filter(self.term) or []

    def get(
        self,
        *,
        per: int | Iterable[int] | None = None,
        **filters,
    ) -> pd.DataFrame:
        """Return the filtered package budget-term dataframe."""

        return self._namespace.get(term=self.term, per=per, **filters)

    def summary(
        self,
        *,
        per: int | Iterable[int] | None = None,
        by: list[str] | tuple[str, ...] | None = None,
        **filters,
    ) -> pd.DataFrame:
        """Summarize the filtered package budget terms."""

        return self._namespace.summary(term=self.term, per=per, by=by, **filters)

    def wide(
        self,
        *,
        per: int | Iterable[int] | None = None,
        index: list[str] | tuple[str, ...] = ("per", "lake"),
        values: str = "q",
        **filters,
    ) -> pd.DataFrame:
        """Pivot the filtered package budget terms to a wide dataframe."""

        return self._namespace.wide(
            term=self.term, per=per, index=index, values=values, **filters
        )


__all__ = [
    "CellBudgetResultsExplorer",
    "StageResultsExplorer",
    "CellPackageResultsNamespace",
    "UzfResultsNamespace",
    "PackageBudgetTermExplorer",
]

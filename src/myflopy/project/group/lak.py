"""Grouped LAK budget/stage/connection results + outputs + accessor."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from myflopy.modflow.mf6.package_explorer import (
    LakStageResultsExplorer,
    _blue_white_red_diverging_colorscale,
    _normalize_connection_type_filter,
    _symmetric_color_limit,
    build_cell_input_map_payload,
    build_group_input_compare_map_payload,
    build_lak_budget_result_table,
    build_lak_connection_table,
    build_lak_q_map_payload,
    get_default_group_compare_colorscale,
)
from myflopy.modflow.mf6.package_surface_water import join_lak_stage
from myflopy.modflow.utils.datatypes.hover import cell_input_hover, compare_hover, lak_hover
from myflopy.project.group._shared import (
    _default_show_layer_elevs,
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

    def __init__(self, group: ModelGroup, *, budget_text: str = "GWF", value_name: str = "q"):
        """Bind a grouped LAK exchange accessor (defaults to the ``GWF`` term's ``q``)."""

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
        value_columns = ["q", "q_per_area"] if "q_per_area" in data.columns else ["q"]
        key_columns = [column for column in data.columns if column not in {"model", *value_columns}]
        ref = data[data["model"] == reference].drop(columns=["model"]).copy()
        rename_map = {column: f"reference_{column}" for column in value_columns}
        ref = ref.rename(columns=rename_map)
        comp = data[data["model"] != reference].merge(ref, on=key_columns, how="inner")
        comp["reference_model"] = reference
        comp["q_diff"] = comp["q"].astype(float) - comp["reference_q"].astype(float)
        if "q_per_area" in value_columns:
            comp["q_per_area_diff"] = comp["q_per_area"].astype(float) - comp["reference_q_per_area"].astype(float)
        ordered = ["model", "reference_model", *key_columns, "q", "reference_q", "q_diff"]
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
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(target_model))
        values, hover = build_lak_q_map_payload(
            selected,
            ncpl=target_model.vor.ncpl,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
        )
        absmax = _symmetric_color_limit(values)
        kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
        kwargs.setdefault("zmax", absmax if absmax > 0 else None)
        kwargs.setdefault("zmid", 0.0)
        kwargs.setdefault("hover_spec", lak_hover())
        return target_model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            # match the SFR convention: gaining (negative q) blue, losing red
            colorscale=colorscale or _blue_white_red_diverging_colorscale(),
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
            selected["q"] = pd.to_numeric(selected["q"], errors="coerce")
            selected["flow_area"] = pd.to_numeric(selected["flow_area"], errors="coerce")
            grouped = selected.groupby("cell", as_index=False).agg({"q": "sum", "flow_area": "sum"})
            grouped["q_per_area"] = np.where(grouped["flow_area"] > 0.0, grouped["q"] / grouped["flow_area"], np.nan)
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
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(reference_model))
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
        return reference_model.cor(
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
        model_name: str | None = None,
        lake: int | None = None,
        layer: int = 0,
        value_column: str = "connection_area",
        agg: str = "sum",
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        colorscale: str | None = None,
        **kwargs,
    ):
        """Build a grouped LAK connection-geometry map for one selected model."""

        target_name = self.group.reference if model_name is None else str(model_name)
        if target_name not in self.group.models:
            raise KeyError(f"Model {target_name!r} is not in the group.")
        target_model = self.group.models[target_name]
        selected = self.get(model_name=target_name, lake=lake, layer=layer)
        if value_column not in selected.columns:
            raise KeyError(f"LAK connection column {value_column!r} was not found.")
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(target_model))
        values, hover = build_cell_input_map_payload(
            selected,
            ncpl=target_model.vor.ncpl,
            value_column=value_column,
            per=None,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        kwargs.setdefault("hover_spec", cell_input_hover(value_column))
        return target_model.cor(
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



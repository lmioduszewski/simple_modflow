"""Grouped SFR budget/stage results + namespace."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from myflopy.modflow.mf6.package_explorer import (
    SfrStageResultsExplorer,
    _exchange_colorscale,
    _symmetric_color_limit,
    build_group_input_compare_map_payload,
    build_sfr_budget_result_table,
    build_sfr_q_map_payload,
    get_default_group_compare_colorscale,
)
from myflopy.modflow.mf6.package_surface_water import join_sfr_stage
from myflopy.modflow.utils.datatypes.hover import compare_hover, sfr_hover
from myflopy.project.group._shared import (
    _default_show_layer_elevs,
    _ensure_group_map_compatible,
    _filter_group_input_table,
    _reduce_to_period_end,
    _resolve_group_compare_target,
    _stable_compare_keys,
)
from myflopy.project.group.results import GroupCellPackageResults, GroupCellPackageResultsNamespace
from myflopy.project.group.spatial import _GroupSpatialView

if TYPE_CHECKING:
    from myflopy.project.group.core import ModelGroup


class GroupSfrBudgetResults(GroupCellPackageResults):
    """Grouped SFR exchange accessor with reach-length-normalized maps."""

    def __init__(self, group: ModelGroup, *, budget_text: str = "SFR", value_name: str = "q_gwf"):
        """Bind a grouped SFR exchange accessor (defaults to the ``SFR`` term's ``q_gwf``).

        The SFR exchange is aquifer-referenced, so the emitted column is
        ``q_gwf`` (accessor stays ``results.q``).
        """

        super().__init__(group, "sfr", budget_text=budget_text, value_name=value_name)

    def get(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """Return aligned SFR exchange rows, including ``q_per_length``."""

        frames = []
        for current_model_name, model in self.group.models.items():
            frame = build_sfr_budget_result_table(
                model,
                budget_text=self.budget_text,
                value_name=self.value_name,
            )
            frame["model"] = current_model_name
            frames.append(frame)

        if not frames:
            return pd.DataFrame()
        combined = pd.concat(frames, ignore_index=True)
        return _filter_group_input_table(
            combined,
            model_name=model_name,
            per=per,
            layer=layer,
            cells=cells,
        )

    def compare(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """Compare grouped SFR exchange rows against the reference model."""

        data = self.get(per=per, layer=layer, cells=cells)
        if data.empty:
            return pd.DataFrame()
        data = _reduce_to_period_end(data)

        reference = self.group.reference
        q = self.value_name  # aquifer-referenced exchange column (q_gwf)
        value_columns = [q, "q_per_length"] if "q_per_length" in data.columns else [q]
        key_columns = _stable_compare_keys(data, value_columns)
        ref = data[data["model"] == reference][key_columns + value_columns].copy()
        rename_map = {column: f"reference_{column}" for column in value_columns}
        ref = ref.rename(columns=rename_map)
        comp = data[data["model"] != reference].merge(ref, on=key_columns, how="inner")
        comp["reference_model"] = reference
        comp[f"{q}_diff"] = comp[q].astype(float) - comp[f"reference_{q}"].astype(float)
        if "q_per_length" in value_columns:
            comp["q_per_length_diff"] = (
                comp["q_per_length"].astype(float) - comp["reference_q_per_length"].astype(float)
            )
        ordered = [
            "model",
            "reference_model",
            *key_columns,
            q,
            f"reference_{q}",
            f"{q}_diff",
        ]
        if "q_per_length" in value_columns:
            ordered.extend(["q_per_length", "reference_q_per_length", "q_per_length_diff"])
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
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Render one model's SFR exchange (reach-length-normalized) as a ``Choro``."""

        del agg
        target_name = self._group_target(model)
        target_model = self.group.models[target_name]
        selected = join_sfr_stage(
            target_model, self.get(model_name=target_name, per=per, layer=layer), per=per
        )
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(target_model))
        values, hover = build_sfr_q_map_payload(
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
        kwargs.setdefault("hover_spec", sfr_hover())
        return target_model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            # SFR exchange is aquifer-referenced (gaining is NEGATIVE), so blue
            # sits at the negative end -- matching the single-model SFR map.
            colorscale=colorscale or _exchange_colorscale("gwf"),
            **kwargs,
        )

    def compare_map(
        self,
        *,
        model_name: str | None = None,
        per: int = 0,
        layer: int = 0,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Build a grouped SFR diff map using normalized exchange per unit length."""

        del agg
        target_name = _resolve_group_compare_target(self.group, model_name)
        _ensure_group_map_compatible(self.group, target_name)
        reference_model = self.group.models[self.group.reference]
        data = self.get(per=per, layer=layer)

        def _normalized_by_cell(frame: pd.DataFrame, current_model_name: str) -> pd.DataFrame:
            """Per-cell exchange per unit reach length (``sum(q)/sum(rlen)``) for one model."""

            selected = frame[frame["model"] == current_model_name].copy()
            if selected.empty:
                return pd.DataFrame(columns=["cell", "q_per_length"])
            q = self.value_name
            selected[q] = pd.to_numeric(selected[q], errors="coerce")
            selected["rlen"] = pd.to_numeric(selected["rlen"], errors="coerce")
            grouped = selected.groupby("cell", as_index=False).agg({q: "sum", "rlen": "sum"})
            grouped["q_per_length"] = np.where(grouped["rlen"] > 0.0, grouped[q] / grouped["rlen"], np.nan)
            return grouped[["cell", "q_per_length"]]

        reference = _normalized_by_cell(data, self.group.reference).rename(
            columns={"q_per_length": "reference_q_per_length"}
        )
        target = _normalized_by_cell(data, target_name)
        comparison = target.merge(reference, on="cell", how="inner")
        comparison["model"] = target_name
        comparison["reference_model"] = self.group.reference
        comparison["per"] = int(per)
        comparison["layer"] = int(layer)
        comparison["q_per_length_diff"] = (
            comparison["q_per_length"].astype(float) - comparison["reference_q_per_length"].astype(float)
        )
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(reference_model))
        values, hover, absmax = build_group_input_compare_map_payload(
            comparison,
            ncpl=reference_model.vor.ncpl,
            value_column="q_per_length",
            diff_column="q_per_length_diff",
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
                "q_per_length",
                "q_per_length_diff",
                title="Δ stream exchange vs reference",
                units={
                    "q_per_length": "ft²/d",
                    "reference_q_per_length": "ft²/d",
                    "q_per_length_diff": "ft²/d",
                },
                labels={"q_per_length_diff": "Δ q / length"},
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

class GroupSfrStageResults(_GroupSpatialView):
    """Grouped accessor for SFR reach stages and stage comparisons.

    Mirrors :class:`GroupLakStageResults` for streams (keyed by ``reach``),
    reading ``model.outputs.sfr.stage``. Inherits the :class:`SpatialView`
    grammar; each panel delegates to the single-model SFR stage explorer.
    """

    def __init__(self, group: ModelGroup):
        """Bind the grouped SFR reach stage-results accessor to ``group``."""

        self.group = group

    def _spatial_value_label(self) -> str:
        """The mapped quantity's label -- reach stage."""

        return "stage"

    def _spatial_map(self, *, per: int = 0, layer: int = 0, model=None, **kwargs):
        """Render one model's reach stage broadcast to reach cells."""

        target_model = self.group.models[self._group_target(model)]
        return SfrStageResultsExplorer(target_model).map(
            per=int(per), layer=int(layer), backend="plotly", **kwargs
        )

    def get(
        self,
        *,
        model_name: str | None = None,
        reach: int | None = None,
        per: int | None = None,
    ) -> pd.DataFrame:
        """Return aligned SFR reach stages for all models."""

        rows: list[pd.DataFrame] = []
        for current_model_name, model in self.group.models.items():
            stage_data = np.asarray(model.outputs.sfr.stage.get(), dtype=float)
            if stage_data.ndim == 1:
                stage_data = stage_data.reshape(-1, 1)
            periods = pd.DataFrame(stage_data, columns=list(range(stage_data.shape[1])))
            periods["per"] = periods.index.astype(int)
            frame = periods.melt(id_vars="per", var_name="reach", value_name="stage")
            frame["reach"] = frame["reach"].astype(int)
            frame["model"] = current_model_name
            rows.append(frame[["model", "per", "reach", "stage"]])

        if not rows:
            return pd.DataFrame(columns=["model", "per", "reach", "stage"])
        combined = pd.concat(rows, ignore_index=True)
        if model_name is not None:
            combined = combined.loc[combined["model"] == str(model_name)].copy()
        if reach is not None:
            combined = combined.loc[combined["reach"] == int(reach)].copy()
        if per is not None:
            combined = combined.loc[combined["per"] == int(per)].copy()
        return combined.reset_index(drop=True)

    def compare(
        self,
        *,
        model_name: str | None = None,
        reach: int | None = None,
        per: int | None = None,
    ) -> pd.DataFrame:
        """Compare SFR reach stages against the reference model."""

        columns = [
            "model", "reference_model", "per", "reach",
            "stage", "reference_stage", "stage_diff",
        ]
        data = self.get(reach=reach, per=per)
        if data.empty:
            return pd.DataFrame(columns=columns)

        reference = self.group.reference
        ref = (
            data.loc[data["model"] == reference, ["per", "reach", "stage"]]
            .rename(columns={"stage": "reference_stage"})
            .copy()
        )
        comp = data.loc[data["model"] != reference].merge(ref, on=["per", "reach"], how="inner")
        comp["reference_model"] = reference
        comp["stage_diff"] = comp["stage"].astype(float) - comp["reference_stage"].astype(float)
        if model_name is not None:
            comp = comp.loc[comp["model"] == str(model_name)].copy()
        return comp[columns]


class GroupSfrResultsNamespace(GroupCellPackageResultsNamespace):
    """Namespace for grouped SFR result accessors."""

    def _field_names(self):
        """The mappable grouped SFR result fields: exchange ``q`` and ``stage``."""

        return ["q", "stage"]

    @property
    def q(self) -> GroupSfrBudgetResults:
        """Return grouped SFR exchange helpers with normalized map behavior."""

        return self._result_accessor

    @property
    def stage(self) -> GroupSfrStageResults:
        """Return grouped SFR reach-stage helpers."""

        return GroupSfrStageResults(self._result_accessor.group)



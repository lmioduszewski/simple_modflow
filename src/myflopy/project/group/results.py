"""Grouped per-cell package results (base class + namespace)."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

import pandas as pd

from myflopy.modflow.mf6.package_explorer import (
    FieldMappable,
    _symmetric_color_limit,
    build_budget_result_table,
    build_cell_input_map_payload,
    build_group_input_compare_map_payload,
    get_default_group_compare_colorscale,
    get_default_package_colorscale,
)
from myflopy.modflow.utils.datatypes.hover import compare_hover, result_hover
from myflopy.project.group._shared import (
    _default_show_layer_elevs,
    _ensure_group_map_compatible,
    _filter_group_input_table,
    _reduce_to_period_end,
    _resolve_group_compare_target,
    _stable_compare_keys,
)
from myflopy.project.group.spatial import _GroupSpatialView

if TYPE_CHECKING:
    from myflopy.project.group.core import ModelGroup


class GroupCellPackageResults(_GroupSpatialView):
    """Grouped accessor for cell-based package result tables.

    Inherits the unified ``map`` / ``mosaic`` / ``animate`` grammar (``mosaic()``
    defaults to one panel per model).
    """

    def __init__(self, group: ModelGroup, package_name: str, *, budget_text: str, value_name: str):
        """Bind a grouped cell-budget result accessor: one ``budget_text`` term + value column."""

        self.group = group
        self.package_name = str(package_name).lower()
        self.budget_text = str(budget_text)
        self.value_name = str(value_name)

    def get(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """Return aligned package-result rows for all models in the group."""

        frames = []
        for current_model_name, model in self.group.models.items():
            frame = build_budget_result_table(
                model,
                budget_text=self.budget_text,
                package_name=self.package_name,
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
        """Compare grouped package results against the reference model."""

        data = self.get(per=per, layer=layer, cells=cells)
        if data.empty:
            return pd.DataFrame()
        data = _reduce_to_period_end(data)

        reference = self.group.reference
        value_columns = [self.value_name]
        key_columns = _stable_compare_keys(data, value_columns)
        ref = data[data["model"] == reference][key_columns + value_columns].copy()
        rename_map = {column: f"reference_{column}" for column in value_columns}
        ref = ref.rename(columns=rename_map)
        comp = data[data["model"] != reference].merge(ref, on=key_columns, how="inner")
        comp["reference_model"] = reference
        comp[f"{self.value_name}_diff"] = (
            comp[self.value_name].astype(float) - comp[f"reference_{self.value_name}"].astype(float)
        )
        ordered = [
            "model",
            "reference_model",
            *key_columns,
            self.value_name,
            f"reference_{self.value_name}",
            f"{self.value_name}_diff",
        ]
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
        """Render one model's package result field as a ``Choro``."""

        target_name = self._group_target(model)
        target_model = self.group.models[target_name]
        selected = self.get(model_name=target_name, per=per, layer=layer)
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(target_model))
        values, hover = build_cell_input_map_payload(
            selected,
            ncpl=target_model.vor.ncpl,
            value_column=self.value_name,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        # The signed exchange columns name their reference frame (q_gwf / q_lake);
        # "q" is only the legacy raw name. All three take the diverging RdBu scale
        # with a symmetric, zero-centred range.
        is_exchange = self.value_name in {"q", "q_gwf", "q_lake"}
        if is_exchange:
            absmax = _symmetric_color_limit(values)
            kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
            kwargs.setdefault("zmax", absmax if absmax > 0 else None)
            kwargs.setdefault("zmid", 0.0)
        kwargs.setdefault(
            "hover_spec",
            result_hover(
                self.value_name,
                title=f"{self.package_name.upper()} {self.value_name}",
                units={self.value_name: "ft³/d"} if is_exchange else None,
            ),
        )
        return target_model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=(
                colorscale
                or ("RdBu" if is_exchange else None)
                or get_default_package_colorscale(self.package_name)
                or "earth"
            ),
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
        """Build a diff choropleth for one grouped package result field."""

        target_name = _resolve_group_compare_target(self.group, model_name)
        _ensure_group_map_compatible(self.group, target_name)
        reference_model = self.group.models[self.group.reference]
        selected = self.compare(model_name=target_name, per=per, layer=layer)
        diff_column = f"{self.value_name}_diff"
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(reference_model))
        values, hover, absmax = build_group_input_compare_map_payload(
            selected,
            ncpl=reference_model.vor.ncpl,
            value_column=self.value_name,
            diff_column=diff_column,
            model_name=target_name,
            reference_model=self.group.reference,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
        kwargs.setdefault("zmax", absmax if absmax > 0 else None)
        kwargs.setdefault("zmid", 0.0)
        kwargs.setdefault(
            "hover_spec",
            compare_hover(
                self.value_name,
                diff_column,
                title=f"Δ {self.package_name.upper()} {self.value_name} vs reference",
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

    # NOTE: the series view is the unified grammar's ``plot()`` (SpatialView) --
    # one line per model (and per cell with ``cells=[...]``).


class GroupCellPackageResultsNamespace(FieldMappable):
    """Namespace for grouped cell-based package result accessors.

    ``results.map()`` maps the exchange field ``q`` for one model;
    ``results.mosaic(by="model")`` and ``results.animate(...)`` inherit the
    unified grammar. (SFR/LAK add a ``stage`` field via their subclasses.)
    """

    _default_field = "q"

    def _field_names(self):
        """The single mappable grouped result field: exchange ``q``."""

        return ["q"]

    def __init__(self, result_accessor: GroupCellPackageResults):
        """Wrap a grouped cell-budget result accessor as a ``field=``-aware namespace."""

        self._result_accessor = result_accessor

    @property
    def q(self) -> GroupCellPackageResults:
        """Return the primary grouped package-exchange result accessor."""

        return self._result_accessor



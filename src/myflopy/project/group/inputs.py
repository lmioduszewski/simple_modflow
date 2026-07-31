"""Grouped package-input exploration (`GroupPackageInputs` + field views)."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

import pandas as pd

from myflopy._logging import get_logger
from myflopy.modflow.mf6.package_explorer import (
    LeafFieldSugar,
    build_cell_input_map_payload,
    build_cell_package_input_table,
    build_group_input_compare_map_payload,
    get_default_group_compare_colorscale,
    get_default_package_colorscale,
    get_default_package_value_column,
    get_package_input_field_names,
)
from myflopy.modflow.mf6.package_tables import PACKAGE_TABLE_UNAVAILABLE
from myflopy.modflow.utils.datatypes.hover import cell_input_hover, compare_hover
from myflopy.project.group._shared import (
    _default_show_layer_elevs,
    _ensure_group_map_compatible,
    _filter_group_input_table,
    _resolve_group_compare_target,
)
from myflopy.project.group.spatial import _GroupSpatialView

if TYPE_CHECKING:
    from myflopy.project.group.core import ModelGroup

logger = get_logger(__name__)


class GroupPackageInputField(_GroupSpatialView):
    """One input field of a grouped package, pinned for the unified grammar.

    ``group.packages.ghb.inputs.cond`` -- same verbs as the parent accessor
    but every panel/series draws this field.
    """

    def __init__(self, parent: GroupPackageInputs, field_name: str):
        """Pin a grouped package-inputs accessor to one field for the unified grammar."""

        self._parent = parent
        self.group = parent.group
        self.package_name = parent.package_name
        self.field_name = str(field_name).lower()

    def get(self, **kwargs) -> pd.DataFrame:
        """Return the aligned input rows narrowed to this field."""

        frame = self._parent.get(**kwargs)
        keep = [
            column
            for column in ("model", "package", "per", "layer", "cell", self.field_name)
            if column in getattr(frame, "columns", [])
        ]
        return frame[keep] if keep else frame

    def summary(self, **kwargs) -> pd.DataFrame:
        """Return the parent's per-model coverage summary."""

        return self._parent.summary(**kwargs)

    def _spatial_map(self, **kwargs):
        """Draw one model's map for this pinned field (delegates to the parent accessor)."""

        kwargs.setdefault("value_column", self.field_name)
        return self._parent._spatial_map(**kwargs)

    def _spatial_value_label(self) -> str:
        """The mapped quantity's label -- this pinned field name."""

        return self.field_name

    def _series_value_column(self, frame) -> str:
        """The column ``plot()`` draws -- this pinned field name."""

        return self.field_name


class GroupPackageInputs(LeafFieldSugar, _GroupSpatialView):
    """Grouped input accessor for simple cell-based MF6 stress-period packages.

    Inherits the unified grammar (``map``/``plot``/``mosaic``/``animate``;
    ``mosaic()`` defaults to one panel per model). Registry-backed fields are
    first-class nodes (``inputs.cond``) and every verb takes ``field=`` as
    sugar over them; without ``field=`` the package default field is drawn.
    """

    def __init__(self, group: ModelGroup, package_name: str):
        """Bind a grouped BC-package inputs accessor to ``group`` for one package."""

        self.group = group
        self.package_name = str(package_name).lower()

    def _field_names(self) -> list[str]:
        """The registry-declared input field names for this package."""

        return get_package_input_field_names(self.package_name)

    def _field_node(self, name: str) -> GroupPackageInputField:
        """Return this accessor pinned to field ``name``."""

        return GroupPackageInputField(self, name)

    def get(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """Return aligned stress-period input data for all models."""

        frames = []
        for current_model_name, model in self.group.models.items():
            # A group member that does not carry this package is skipped, NOT an
            # error. It used to raise, and because `compare()` (and so
            # `ModelDiff._value_cells_changed`) builds the whole group before
            # filtering to one model, ONE member without the package discarded
            # the comparison for every other member -- reporting models with
            # genuinely different values as identical to the reference.
            try:
                frame = build_cell_package_input_table(model, self.package_name)
            except PACKAGE_TABLE_UNAVAILABLE:
                logger.debug(
                    "%s carries no readable %s; leaving it out of the group table",
                    current_model_name, self.package_name, exc_info=True,
                )
                continue
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

    def summary(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """Return a per-model summary of this package's input coverage.

        Mirrors the single-model ``model.packages.<pkg>.inputs.summary()`` verb:
        one row per model with the number of applied cells and stress periods.
        """

        data = self.get(model_name=model_name, per=per, layer=layer)
        columns = ["model", "package", "cells", "periods"]
        if data.empty:
            return pd.DataFrame(columns=columns)
        rows = []
        for current_model_name, sub in data.groupby("model"):
            rows.append(
                {
                    "model": current_model_name,
                    "package": self.package_name,
                    "cells": int(sub["cell"].nunique()) if "cell" in sub.columns else 0,
                    "periods": int(sub["per"].nunique()) if "per" in sub.columns else 0,
                }
            )
        return pd.DataFrame(rows, columns=columns)

    def compare(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """Compare input values for all models against the reference model."""

        data = self.get(per=per, layer=layer, cells=cells)
        if data.empty:
            return pd.DataFrame()

        reference = self.group.reference
        value_columns = [
            column
            for column in data.columns
            if column not in {"model", "per", "layer", "cell"}
            and pd.api.types.is_numeric_dtype(data[column])
        ]
        key_columns = [column for column in data.columns if column not in {"model", *value_columns}]
        ref = data[data["model"] == reference].drop(columns=["model"]).copy()
        rename_map = {column: f"reference_{column}" for column in value_columns}
        ref = ref.rename(columns=rename_map)
        comp = data[data["model"] != reference].merge(ref, on=key_columns, how="inner")
        comp["reference_model"] = reference
        for column in value_columns:
            comp[f"{column}_diff"] = comp[column].astype(float) - comp[f"reference_{column}"].astype(float)
        ordered = ["model", "reference_model", *key_columns]
        for column in value_columns:
            ordered.extend([column, f"reference_{column}", f"{column}_diff"])
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
        value_column: str | None = None,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Render one model's package-input field as a raw-value ``Choro``."""

        target_name = self._group_target(model)
        target_model = self.group.models[target_name]
        selected = self.get(model_name=target_name, per=per, layer=layer)
        fallback = get_default_package_value_column(self.package_name)
        source = selected if not selected.empty else self.get(model_name=target_name)
        numeric_columns = [
            column
            for column in source.columns
            if column not in {"model", "package", "per", "layer", "cell"}
            and pd.api.types.is_numeric_dtype(source[column])
        ]
        chosen_value_column = (
            value_column
            or (fallback if fallback in source.columns else None)
            or (numeric_columns[0] if numeric_columns else None)
        )
        if chosen_value_column is None:
            raise ValueError(f"Could not infer a numeric value column for package {self.package_name!r}.")
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(target_model))
        values, hover = build_cell_input_map_payload(
            selected,
            ncpl=target_model.vor.ncpl,
            value_column=chosen_value_column,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        kwargs.setdefault("hover_spec", cell_input_hover(chosen_value_column))
        return target_model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or get_default_package_colorscale(self.package_name),
            **kwargs,
        )

    def compare_map(
        self,
        *,
        model_name: str | None = None,
        per: int = 0,
        layer: int = 0,
        value_column: str | None = None,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Build a diff choropleth against the group's reference model.

        Parameters
        ----------
        model_name
            Non-reference model to compare against the reference. When the
            group contains exactly one non-reference model, it is inferred.
        per, layer
            Zero-based stress period and layer to map.
        value_column
            Numeric input field whose difference should be mapped.
        multiplier, fill_value, agg
            Passed through to the shared diff-map payload builder.
        colorscale
            Optional diverging colorscale override.
        kwargs
            Forwarded to ``reference_model.cor(...)``.
        """

        target_name = _resolve_group_compare_target(self.group, model_name)
        _ensure_group_map_compatible(self.group, target_name)
        reference_model = self.group.models[self.group.reference]
        selected = self.compare(model_name=target_name, per=per, layer=layer)
        fallback = get_default_package_value_column(self.package_name)
        source = selected if not selected.empty else self.compare(model_name=target_name)
        numeric_columns = [
            column
            for column in source.columns
            if column not in {"per", "layer", "cell", "model", "reference_model", "package"}
            and pd.api.types.is_numeric_dtype(source[column])
            and not column.startswith("reference_")
            and not column.endswith("_diff")
        ]
        chosen_value_column = (
            value_column
            or (fallback if fallback in source.columns else None)
            or (numeric_columns[0] if numeric_columns else None)
        )
        if chosen_value_column is None:
            raise ValueError(f"Could not infer a numeric value column for package {self.package_name!r}.")
        diff_column = f"{chosen_value_column}_diff"
        if diff_column not in source.columns:
            raise KeyError(f"Difference column {diff_column!r} was not found.")
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(reference_model))
        values, hover, absmax = build_group_input_compare_map_payload(
            selected,
            ncpl=reference_model.vor.ncpl,
            value_column=chosen_value_column,
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
        kwargs.setdefault("hover_spec", compare_hover(chosen_value_column, diff_column))
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



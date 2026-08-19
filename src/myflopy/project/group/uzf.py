"""Grouped UZF inputs/fields/results + accessor."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

import pandas as pd

from myflopy.modflow.mf6.package_explorer import (
    FieldMappable,
    build_cell_input_map_payload,
    build_group_input_compare_map_payload,
    build_uzf_field_input_table,
    get_default_budget_term,
    get_default_group_compare_colorscale,
    get_default_package_colorscale,
    get_package_input_field_spec,
)
from myflopy.modflow.utils.datatypes.hover import cell_input_hover, compare_hover
from myflopy.project.group._shared import (
    _ensure_group_map_compatible,
    _filter_group_input_table,
    _resolve_group_compare_target,
)
from myflopy.project.group.results import GroupCellPackageResults
from myflopy.project.group.spatial import _GroupSpatialView

if TYPE_CHECKING:
    from myflopy.project.group.core import ModelGroup


class GroupUzfFieldAccessor(_GroupSpatialView):
    """Grouped accessor for one UZF perioddata field such as ``finf``."""

    def __init__(self, group: ModelGroup, field_name: str):
        """Bind a grouped UZF perioddata-field accessor to ``group`` for one field."""

        self.group = group
        self.field_name = str(field_name)

    def get(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """Return one aligned UZF perioddata field for all models."""

        frames = []
        for current_model_name, model in self.group.models.items():
            frame = build_uzf_field_input_table(model, self.field_name)
            frame["model"] = current_model_name
            frames.append(frame)

        if not frames:
            return pd.DataFrame(columns=["model", "per", "ifno", "layer", "cell", self.field_name])

        combined = pd.concat(frames, ignore_index=True)
        combined["ifno"] = combined["ifno"].astype(int)
        combined["layer"] = combined["layer"].astype(int)
        combined["cell"] = combined["cell"].astype(int)
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
        """Compare the grouped UZF field against the reference model."""

        data = self.get(per=per, layer=layer, cells=cells)
        if data.empty:
            return pd.DataFrame()

        reference = self.group.reference
        ref = (
            data[data["model"] == reference]
            .rename(columns={self.field_name: f"reference_{self.field_name}"})
            .drop(columns=["model"])
        )
        comp = data[data["model"] != reference].merge(
            ref,
            on=["per", "ifno", "layer", "cell"],
            how="inner",
        )
        comp["reference_model"] = reference
        comp[f"{self.field_name}_diff"] = (
            comp[self.field_name].astype(float) - comp[f"reference_{self.field_name}"].astype(float)
        )
        result = comp[
            [
                "model",
                "reference_model",
                "per",
                "ifno",
                "layer",
                "cell",
                self.field_name,
                f"reference_{self.field_name}",
                f"{self.field_name}_diff",
            ]
        ]
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
        """Render one model's UZF field as a raw-value ``Choro``."""

        target_name = self._group_target(model)
        target_model = self.group.models[target_name]
        selected = self.get(model_name=target_name, per=per, layer=layer)
        values, hover = build_cell_input_map_payload(
            selected,
            ncpl=target_model.vor.ncpl,
            value_column=self.field_name,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        kwargs.setdefault("hover_spec", cell_input_hover(self.field_name))
        return target_model.plot.map(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or get_default_package_colorscale(f"uzf_{self.field_name}") or "earth",
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
        """Build a diff choropleth for one grouped UZF field."""

        target_name = _resolve_group_compare_target(self.group, model_name)
        _ensure_group_map_compatible(self.group, target_name)
        reference_model = self.group.models[self.group.reference]
        selected = self.compare(model_name=target_name, per=per, layer=layer)
        diff_column = f"{self.field_name}_diff"
        values, hover, absmax = build_group_input_compare_map_payload(
            selected,
            ncpl=reference_model.vor.ncpl,
            value_column=self.field_name,
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
        kwargs.setdefault("hover_spec", compare_hover(self.field_name, diff_column))
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


class GroupUzfInputs:
    """Namespace for grouped UZF input fields."""

    def __init__(self, group: ModelGroup):
        """Bind the grouped UZF inputs namespace to ``group``."""

        self.group = group

    def _field(self, field_name: str) -> GroupUzfFieldAccessor:
        """Return one registry-backed grouped UZF field accessor."""

        field_spec = get_package_input_field_spec("uzf", field_name)
        if field_spec is None:
            raise AttributeError(f"{type(self).__name__!s} has no UZF input field {field_name!r}")
        return GroupUzfFieldAccessor(self.group, field_spec.name)

    def __getattr__(self, field_name: str) -> GroupUzfFieldAccessor:
        """Return a registry-backed grouped UZF field accessor."""

        return self._field(field_name)

    @property
    def finf(self) -> GroupUzfFieldAccessor:
        """Return grouped infiltration accessors."""

        return self._field("finf")

    @property
    def pet(self) -> GroupUzfFieldAccessor:
        """Return grouped potential evapotranspiration accessors."""

        return self._field("pet")

    @property
    def extdp(self) -> GroupUzfFieldAccessor:
        """Return grouped ET extinction-depth accessors."""

        return self._field("extdp")

    @property
    def extwc(self) -> GroupUzfFieldAccessor:
        """Return grouped ET extinction-water-content accessors."""

        return self._field("extwc")

    @property
    def ha(self) -> GroupUzfFieldAccessor:
        """Return grouped surface-depression-storage-depth accessors."""

        return self._field("ha")

    @property
    def hroot(self) -> GroupUzfFieldAccessor:
        """Return grouped root-zone-thickness accessors."""

        return self._field("hroot")

    @property
    def rootact(self) -> GroupUzfFieldAccessor:
        """Return grouped root-activity accessors."""

        return self._field("rootact")


class GroupUzfPackageAccessor:
    """Namespace for grouped UZF exploration helpers."""

    def __init__(self, accessor: GroupUzfInputs):
        """Wrap a grouped UZF inputs ``accessor`` and expose ``.inputs`` / ``.results``."""

        self.inputs = accessor

    @property
    def results(self):
        """Return grouped UZF result helpers."""

        budget_text, value_name = get_default_budget_term("uzf_gwrch") or ("UZF-GWRCH", "gwrch")
        return GroupUzfResultsNamespace(
            GroupCellPackageResults(
                self.inputs.group,
                "uzf",
                budget_text=budget_text,
                value_name=value_name,
            )
        )


class GroupUzfResultsNamespace(FieldMappable):
    """Namespace for grouped UZF result accessors.

    Fields: ``gwrch`` (groundwater recharge, the default) and ``sat``.
    """

    _default_field = "gwrch"

    def _field_names(self):
        """The mappable grouped UZF result fields: ``gwrch`` and ``sat``."""

        return ["gwrch", "sat"]

    def __init__(self, gwrch_accessor: GroupCellPackageResults, sat_accessor: GroupCellPackageResults | None = None):
        """Bind the grouped UZF results namespace; ``sat`` is built lazily if not given."""

        self._gwrch = gwrch_accessor
        self._sat = sat_accessor

    @property
    def gwrch(self) -> GroupCellPackageResults:
        """Return grouped groundwater-recharge results from UZF."""

        return self._gwrch

    @property
    def sat(self) -> GroupCellPackageResults:
        """Return grouped UZF saturation results."""

        if self._sat is None:
            budget_text, value_name = get_default_budget_term("uzf_sat") or ("DATA-SAT", "sat")
            self._sat = GroupCellPackageResults(
                self._gwrch.group,
                "uzf",
                budget_text=budget_text,
                value_name=value_name,
            )
        return self._sat



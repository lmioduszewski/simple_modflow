"""Input explorer classes behind model.packages."""

from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.modflow.mf6.package_registry import (
    FieldSpec,
    get_default_package_colorscale,
    get_default_package_value_column,
    get_package_explorer_spec,
    get_package_input_field_names,
    get_package_input_field_spec,
)
from myflopy.modflow.mf6.package_explorer_utils import (
    _default_show_layer_elevs,
    _filter_normalized_table,
)
from myflopy.modflow.mf6.package_tables import (
    build_cell_package_input_table,
    build_uzf_field_input_table,
    build_uzf_field_input_wide_table,
    summarize_input_table,
)
from myflopy.modflow.utils.datatypes.hover import cell_input_hover
from myflopy.modflow.mf6.package_plotting import (
    FieldMappable,
    SpatialView,
    _apply_backend,
    _as_layer_cell_property,
    build_cell_input_map_payload,
)


class CellPackageInputsExplorer(FieldMappable):
    """Normalized input explorer for one cell-based MF6 stress-period package.

    Each registry-backed field is a first-class node (``ghb.inputs.cond``) with
    the full unified grammar; the namespace verbs take ``field=`` as sugar --
    ``inputs.map(field="cond")`` is exactly ``inputs.cond.map()``, and
    ``inputs.map()`` draws the package's default field.
    """

    def __init__(self, model: "SimulationBase", package_name: str):
        self.model = model
        self.package_name = str(package_name).lower()

    @property
    def _default_field(self) -> str:
        return get_default_package_value_column(self.package_name)

    def _field_names(self) -> list[str]:
        return get_package_input_field_names(self.package_name)

    def __getattr__(self, field_name: str) -> "CellPackageInputFieldExplorer":
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
    def default(self) -> "CellPackageInputFieldExplorer":
        """Return the registry-defined preferred input field."""

        return getattr(self, get_default_package_value_column(self.package_name))

    # map/plot/xs/mosaic/animate come from FieldMappable: they dispatch to the
    # (default) field node, so ``inputs.map()`` == ``inputs.<default>.map()``.


class CellPackageInputFieldExplorer(SpatialView):
    """Field-specific view over a cell package's normalized input table."""

    def __init__(self, inputs: CellPackageInputsExplorer, field_spec: FieldSpec):
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
        per: int = 0,
        layer: int = 0,
        multiplier: float = 1.0,
        fill_value: float | None = None,
        agg: str | None = None,
        colorscale: str | None = None,
        backend: str = "plotly",
        **kwargs,
    ):
        """Build a choropleth for this specific input field."""

        selected = self.inputs.get(per=per, layer=layer)
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(self.model))
        values, hover = build_cell_input_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            value_column=self.field_name,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=self.field_spec.fill_value if fill_value is None else fill_value,
            agg=self.field_spec.agg if agg is None else agg,
        )
        kwargs.setdefault("hover_spec", cell_input_hover(self.field_name))
        choro = self.model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=(
                colorscale
                or self.field_spec.colorscale
                or get_default_package_colorscale(self.package_name)
            ),
            **kwargs,
        )
        return _apply_backend(choro, backend)


class UzfFieldInputsExplorer(SpatialView):
    """Normalized explorer for one UZF perioddata field."""

    def __init__(self, model: "SimulationBase", field_name: str):
        self.model = model
        self.field_name = str(field_name)

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
        per: int = 0,
        layer: int = 0,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        backend: str = "plotly",
        **kwargs,
    ):
        """Build a choropleth for the selected UZF field."""

        selected = self.get(per=per, layer=layer)
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(self.model))
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
        kwargs.setdefault("hover_spec", cell_input_hover(self.field_name))
        choro = self.model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale
            or get_default_package_colorscale(f"uzf_{self.field_name}")
            or "earth",
            **kwargs,
        )
        return _apply_backend(choro, backend)


class UzfInputsNamespace(FieldMappable):
    """Namespace for normalized UZF input explorers.

    Each perioddata field (``finf``, ``pet``, ...) is a first-class node with
    the full unified grammar; the namespace verbs take ``field=`` as sugar
    (default field: ``finf``).
    """

    _default_field = "finf"

    def __init__(self, model: "SimulationBase"):
        self.model = model

    def _field_names(self) -> list[str]:
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
        model: "SimulationBase",
        package_name: str,
        field_name: str,
        *,
        label: str | None = None,
        colorscale: str = "Viridis",
    ):
        self.model = model
        self.package_name = str(package_name).lower()
        self.field_name = str(field_name)
        self.label = label or f"{self.package_name}.{self.field_name}"
        self.colorscale = colorscale

    def _array(self) -> np.ndarray:
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
        per: int = 0,
        layer: int = 0,
        colorscale: str | None = None,
        backend: str = "plotly",
        **kwargs,
    ):
        """Build a choropleth for one layer of this static array field."""

        del per
        selected = self.get(layer=layer)
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(self.model))
        values, hover = build_cell_input_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            value_column=self.field_name,
            per=None,
            layer=layer,
            agg="first",
        )
        kwargs.setdefault("hover_spec", cell_input_hover(self.field_name))
        choro = self.model.cor(
            per=0,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or self.colorscale,
            **kwargs,
        )
        return _apply_backend(choro, backend)


__all__ = [
    "CellPackageInputsExplorer",
    "CellPackageInputFieldExplorer",
    "UzfFieldInputsExplorer",
    "UzfInputsNamespace",
    "StaticArrayFieldExplorer",
]

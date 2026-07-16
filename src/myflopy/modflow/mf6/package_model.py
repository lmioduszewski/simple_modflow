"""Top-level model package explorer wiring."""

from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.modflow.mf6.package_inputs import (
    CellPackageInputsExplorer,
    StaticArrayFieldExplorer,
    UzfInputsNamespace,
)
from myflopy.modflow.mf6.package_registry import (
    get_package_explorer_spec,
)
from myflopy.modflow.mf6.package_results import (
    CellPackageResultsNamespace,
    UzfResultsNamespace,
)
from myflopy.modflow.mf6.package_surface_water import (
    LakPackageExplorer,
    SfrPackageExplorer,
    SurfaceWaterPackageExplorer,
)


class PackageExplorer:
    """Namespace for one package's preferred exploration helpers."""

    def __init__(self, model: SimulationBase, package_name: str):
        """Bind a generic package explorer (``.inputs`` / ``.results``) to ``model``."""

        self.model = model
        self.package_name = str(package_name).lower()

    @property
    def inputs(self) -> CellPackageInputsExplorer:
        """Return normalized input helpers for this package."""

        return CellPackageInputsExplorer(self.model, self.package_name)

    @property
    def results(self) -> CellPackageResultsNamespace:
        """Return normalized result helpers for this package."""

        return CellPackageResultsNamespace(self.model, self.package_name)


class StaticArrayPackageExplorer:
    """Namespace for static array fields in one package."""

    def __init__(
        self,
        model: SimulationBase,
        package_name: str,
        fields: Mapping[str, dict[str, str]],
    ):
        """Bind a static-array explorer to ``model`` for a package's declared array ``fields``."""

        self.model = model
        self.package_name = str(package_name).lower()
        self._fields = dict(fields)

    def _available_field_items(self) -> list[tuple[str, dict[str, str]]]:
        """The declared ``(field_name, metadata)`` pairs that actually exist on the package."""

        package = self.model.package(self.package_name)
        available = []
        for field_name, metadata in self._fields.items():
            data = getattr(package, field_name, None)
            if data is None:
                continue
            available.append((field_name, metadata))
        return available

    @property
    def fields(self) -> pd.DataFrame:
        """Return supported static array fields available on this package."""

        return pd.DataFrame(
            [
                {
                    "field": field_name,
                    "label": metadata.get("label"),
                    "colorscale": metadata.get("colorscale", "Viridis"),
                }
                for field_name, metadata in self._available_field_items()
            ]
        ).reindex(columns=["field", "label", "colorscale"])

    def summary(self) -> pd.DataFrame:
        """Return one compact summary row per supported array field."""

        frames = [
            self._field(field_name).summary()
            for field_name in self.fields["field"].tolist()
        ]
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
        return summary.merge(
            self.fields.rename(columns={"field": "field_name"}),
            on="field_name",
            how="left",
        )

    def _field(self, field_name: str) -> StaticArrayFieldExplorer:
        """A field-pinned explorer for one static array (raises if unknown or absent on the package)."""

        metadata = self._fields.get(str(field_name))
        if metadata is None:
            raise AttributeError(
                f"{type(self).__name__!s} has no array field {field_name!r}"
            )
        package = self.model.package(self.package_name)
        if getattr(package, str(field_name), None) is None:
            raise AttributeError(
                f"Package {self.package_name!r} has no available array field {field_name!r}"
            )
        return StaticArrayFieldExplorer(
            self.model,
            self.package_name,
            str(field_name),
            label=metadata.get("label"),
            colorscale=metadata.get("colorscale", "Viridis"),
        )

    def __getattr__(self, field_name: str) -> StaticArrayFieldExplorer:
        """Return a supported static array field explorer."""

        return self._field(field_name)


class UzfPackageExplorer:
    """Top-level UZF package explorer namespace."""

    def __init__(self, model: SimulationBase):
        """Bind the top-level UZF explorer (``.inputs`` / ``.results``) to ``model``."""

        self.model = model

    @property
    def inputs(self) -> UzfInputsNamespace:
        """Return the UZF input exploration namespace."""

        return UzfInputsNamespace(self.model)

    @property
    def results(self) -> UzfResultsNamespace:
        """Return the UZF result exploration namespace."""

        return UzfResultsNamespace(self.model)


class ModelPackages:
    """Preferred package exploration namespace for one model/run.

    Examples
    --------
    ``model.packages.rch.inputs.get()``
        Normalized recharge input table.
    ``model.packages.rch.inputs.map(per=0)``
        Recharge choropleth using the shared choropleth styling.
    ``model.packages.uzf.inputs.finf.summary()``
        Compact summary of UZF infiltration inputs.
    """

    def __init__(self, model: SimulationBase):
        """Bind the preferred ``model.packages`` exploration namespace to ``model``."""

        self.model = model

    def __getattr__(self, package_name: str) -> PackageExplorer:
        """Return a registry-backed generic package explorer."""

        spec = get_package_explorer_spec(package_name)
        if spec is None or spec.kind != "cell_stress":
            raise AttributeError(
                f"{type(self).__name__!s} has no package {package_name!r}"
            )
        return PackageExplorer(self.model, spec.name)

    @property
    def rch(self) -> PackageExplorer:
        """Recharge package exploration helpers."""

        return PackageExplorer(self.model, "rch")

    @property
    def chd(self) -> PackageExplorer:
        """Constant-head package exploration helpers."""

        return PackageExplorer(self.model, "chd")

    @property
    def drn(self) -> PackageExplorer:
        """Drain package exploration helpers."""

        return PackageExplorer(self.model, "drn")

    @property
    def ghb(self) -> PackageExplorer:
        """General-head boundary package exploration helpers."""

        return PackageExplorer(self.model, "ghb")

    @property
    def wel(self) -> PackageExplorer:
        """Well package exploration helpers."""

        return PackageExplorer(self.model, "wel")

    @property
    def uzf(self) -> UzfPackageExplorer:
        """UZF package exploration helpers."""

        return UzfPackageExplorer(self.model)

    @property
    def ic(self) -> StaticArrayPackageExplorer:
        """Initial conditions array exploration helpers."""

        return StaticArrayPackageExplorer(
            self.model,
            "ic",
            {
                "strt": {"label": "Starting head", "colorscale": "Viridis"},
            },
        )

    @property
    def npf(self) -> StaticArrayPackageExplorer:
        """Node property flow array exploration helpers."""

        return StaticArrayPackageExplorer(
            self.model,
            "npf",
            {
                "k": {
                    "label": "Horizontal hydraulic conductivity",
                    "colorscale": "Viridis",
                },
                "k22": {
                    "label": "Horizontal hydraulic conductivity K22",
                    "colorscale": "Viridis",
                },
                "k33": {
                    "label": "Vertical hydraulic conductivity",
                    "colorscale": "Viridis",
                },
            },
        )

    @property
    def sto(self) -> StaticArrayPackageExplorer:
        """Storage package array exploration helpers."""

        return StaticArrayPackageExplorer(
            self.model,
            "sto",
            {
                "ss": {"label": "Specific storage", "colorscale": "Viridis"},
                "sy": {"label": "Specific yield", "colorscale": "Viridis"},
            },
        )

    @property
    def lak(self) -> LakPackageExplorer:
        """LAK package exploration helpers."""

        return LakPackageExplorer(self.model)

    @property
    def sfr(self) -> SfrPackageExplorer:
        """SFR package exploration helpers."""

        return SfrPackageExplorer(self.model)

    @property
    def surface_water(self) -> SurfaceWaterPackageExplorer:
        """Combined SFR/LAK exploration helpers."""

        return SurfaceWaterPackageExplorer(self.model)


__all__ = [
    "PackageExplorer",
    "StaticArrayPackageExplorer",
    "UzfPackageExplorer",
    "ModelPackages",
]

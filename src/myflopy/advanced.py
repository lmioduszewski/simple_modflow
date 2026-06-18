"""Readable package-spec factories for advanced MODFLOW 6 GWF packages."""

from __future__ import annotations

from collections.abc import Iterable
import functools
from typing import Any

import flopy

from myflopy.specs import PackageSpec


def _named_options(name: str, options: dict[str, Any]) -> dict[str, Any]:
    """Return package options with consistent FloPy package naming."""

    return {"pname": name, "filename": f"{{model_name}}.{name}", **options}


def _build_named(constructor, model, **options):
    """Build a package after resolving the model-dependent filename."""

    options = {
        key: value.format(model_name=model.name) if isinstance(value, str) else value
        for key, value in options.items()
    }
    return constructor(model, **options)


def _factory(constructor):
    # functools.partial over a module-level function (not a lambda) so the
    # resulting PackageSpec can be pickled and reused across variants/sessions.
    return functools.partial(_build_named, constructor)


def wel_spec(
    stress_period_data,
    *,
    name: str = "wel",
    auxiliary=None,
    boundnames: bool = True,
    **options,
) -> PackageSpec:
    """Return a WEL package specification."""

    values = _named_options(
        name,
        {
            "stress_period_data": stress_period_data,
            "auxiliary": auxiliary,
            "boundnames": boundnames,
            "save_flows": True,
            **options,
        },
    )
    return PackageSpec(name, _factory(flopy.mf6.ModflowGwfwel), values)


def chd_spec(
    stress_period_data,
    *,
    name: str = "chd",
    boundnames: bool = False,
    **options,
) -> PackageSpec:
    """Return a CHD package specification."""

    values = _named_options(
        name,
        {
            "stress_period_data": stress_period_data,
            "boundnames": boundnames,
            "save_flows": True,
            **options,
        },
    )
    return PackageSpec(name, _factory(flopy.mf6.ModflowGwfchd), values)


def drn_spec(
    stress_period_data,
    *,
    name: str = "drn",
    boundnames: bool = False,
    **options,
) -> PackageSpec:
    """Return a DRN package specification."""

    values = _named_options(
        name,
        {
            "stress_period_data": stress_period_data,
            "boundnames": boundnames,
            "save_flows": True,
            **options,
        },
    )
    return PackageSpec(name, _factory(flopy.mf6.ModflowGwfdrn), values)


def ghb_spec(
    stress_period_data,
    *,
    name: str = "ghb",
    auxiliary=None,
    boundnames: bool = False,
    **options,
) -> PackageSpec:
    """Return a GHB package specification."""

    values = _named_options(
        name,
        {
            "stress_period_data": stress_period_data,
            "auxiliary": auxiliary,
            "boundnames": boundnames,
            "save_flows": True,
            **options,
        },
    )
    return PackageSpec(name, _factory(flopy.mf6.ModflowGwfghb), values)


def rch_spec(
    stress_period_data,
    *,
    name: str = "rch",
    boundnames: bool = False,
    **options,
) -> PackageSpec:
    """Return a list-based RCH package specification."""

    values = _named_options(
        name,
        {
            "stress_period_data": stress_period_data,
            "maxbound": max((len(records) for records in stress_period_data.values()), default=0),
            "boundnames": boundnames,
            "save_flows": True,
            **options,
        },
    )
    return PackageSpec(name, _factory(flopy.mf6.ModflowGwfrch), values)


def uzf_spec(
    packagedata,
    perioddata,
    *,
    nuzfcells: int | None = None,
    name: str = "uzf",
    mover: bool = False,
    simulate_et: bool = False,
    **options,
) -> PackageSpec:
    """Return a UZF package specification from prepared UZF data."""

    count = len(packagedata) if nuzfcells is None else nuzfcells
    values = _named_options(
        name,
        {
            "nuzfcells": count,
            "packagedata": packagedata,
            "perioddata": perioddata,
            "mover": mover,
            "simulate_et": simulate_et,
            "save_flows": True,
            "budget_filerecord": "{model_name}_budget.uzf",
            "budgetcsv_filerecord": f"{{model_name}}_{name}_budget.csv",
            "package_convergence_filerecord": f"{{model_name}}_{name}_package_convergence.csv",
            **options,
        },
    )
    return PackageSpec(name, _factory(flopy.mf6.ModflowGwfuzf), values)


def lak_spec(
    packagedata,
    connectiondata,
    perioddata,
    *,
    nlakes: int | None = None,
    name: str = "lak",
    noutlets: int = 0,
    ntables: int = 0,
    mover: bool = False,
    **options,
) -> PackageSpec:
    """Return a LAK package specification from prepared lake data."""

    count = len(packagedata) if nlakes is None else nlakes
    values = _named_options(
        name,
        {
            "nlakes": count,
            "noutlets": noutlets,
            "ntables": ntables,
            "packagedata": packagedata,
            "connectiondata": connectiondata,
            "perioddata": perioddata,
            "mover": mover,
            "save_flows": True,
            "stage_filerecord": f"{{model_name}}_{name}_stage.lak",
            "budget_filerecord": f"{{model_name}}_{name}_budget.lak",
            "budgetcsv_filerecord": "{model_name}_lake_budget.csv",
            "package_convergence_filerecord": "{model_name}_lake_convergence.csv",
            **options,
        },
    )
    return PackageSpec(name, _factory(flopy.mf6.ModflowGwflak), values)


def sfr_spec(
    packagedata,
    connectiondata,
    perioddata,
    *,
    nreaches: int | None = None,
    name: str = "sfr",
    mover: bool = False,
    diversions=None,
    **options,
) -> PackageSpec:
    """Return an SFR package specification from prepared stream data."""

    count = len(packagedata) if nreaches is None else nreaches
    values = _named_options(
        name,
        {
            "nreaches": count,
            "packagedata": packagedata,
            "connectiondata": connectiondata,
            "perioddata": perioddata,
            "diversions": diversions,
            "mover": mover,
            "save_flows": True,
            "stage_filerecord": f"{{model_name}}_{name}_stage.sfr",
            "budget_filerecord": f"{{model_name}}_{name}_budget.sfr",
            **options,
        },
    )
    return PackageSpec(name, _factory(flopy.mf6.ModflowGwfsfr), values)


def _package_names(packages: Iterable[Any]) -> tuple[str, ...]:
    """Return normalized package names from MVR package records."""

    names = []
    for record in packages:
        value = record[0] if isinstance(record, (list, tuple)) else record
        names.append(str(value))
    return tuple(names)


def mvr_spec(
    packages,
    perioddata,
    *,
    name: str = "mvr",
    maxmvr: int | None = None,
    maxpackages: int | None = None,
    **options,
) -> PackageSpec:
    """Return an MVR spec that declares and validates package dependencies."""

    required = _package_names(packages)
    if len(set(required)) != len(required):
        raise ValueError("MVR packages must contain unique package names.")
    required_set = set(required)
    for period, records in perioddata.items():
        for record in records:
            source, destination = str(record[0]), str(record[2])
            missing = sorted({source, destination} - required_set)
            if missing:
                raise ValueError(
                    f"MVR period {period} references undeclared packages: {', '.join(missing)}"
                )

    count = max((len(records) for records in perioddata.values()), default=0)
    resolved_maxmvr = count if maxmvr is None else maxmvr
    resolved_maxpackages = len(required) if maxpackages is None else maxpackages
    if resolved_maxmvr < count:
        raise ValueError(f"maxmvr must be at least {count}.")
    if resolved_maxpackages < len(required):
        raise ValueError(f"maxpackages must be at least {len(required)}.")
    values = _named_options(
        name,
        {
            "maxmvr": resolved_maxmvr,
            "maxpackages": resolved_maxpackages,
            "packages": packages,
            "perioddata": perioddata,
            "budget_filerecord": f"{{model_name}}.{name}.bud",
            "budgetcsv_filerecord": f"{{model_name}}.{name}.csv",
            **options,
        },
    )
    return PackageSpec(
        name,
        _factory(flopy.mf6.ModflowGwfmvr),
        values,
        requires=required,
    )


__all__ = [
    "chd_spec",
    "drn_spec",
    "ghb_spec",
    "lak_spec",
    "mvr_spec",
    "rch_spec",
    "sfr_spec",
    "uzf_spec",
    "wel_spec",
]

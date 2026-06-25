"""Low-level ``*_spec`` factories: the raw-FloPy escape hatch for package specs.

Each function here wraps a FloPy GWF package constructor in a serializable
:class:`~myflopy.specs.PackageSpec` -- the data form of the package-first API.
They take native FloPy inputs verbatim (``stress_period_data`` dicts, ``recarray``
package/connection data, period data) and add only myflopy's conventions: a
consistent ``pname``/``filename``, ``save_flows`` on by default, and the standard
output file records for the advanced packages.

You normally reach these through the package-first facade's ``.flopy(...)`` form,
e.g. ``mf.ghb.flopy(stress_period_data=...)`` calls :func:`ghb_spec`, and
``mf.uzf.flopy(...)`` calls :func:`uzf_spec`. Use them directly only when you want
the raw FloPy data path and not the high-level GIS-aware builders (``mf.uzf(...)``,
``mf.sfr(...)``, ``mf.lak(...)``) that resolve cells from a grid for you. The
returned spec is built at run time by ``PackageSpec.build(model)``.
"""

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
    """Return a WEL (well) package spec from raw FloPy stress-period data.

    The data form behind ``mf.wel.flopy(...)``. Wraps
    :class:`flopy.mf6.ModflowGwfwel` with ``save_flows`` enabled.

    Parameters
    ----------
    stress_period_data
        FloPy stress-period mapping ``{per: [(cellid, q, [aux...], [bname]), ...]}``.
    name
        Package name used for ``pname`` and the output filename.
    auxiliary
        Optional auxiliary variable name(s).
    boundnames
        Enable named boundaries (default ``True``).
    **options
        Extra keyword options passed straight to the FloPy constructor.
    """

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
    """Return a CHD (constant-head) package spec from raw FloPy stress-period data.

    The data form behind ``mf.chd.flopy(...)``. Wraps
    :class:`flopy.mf6.ModflowGwfchd` with ``save_flows`` enabled.

    Parameters
    ----------
    stress_period_data
        FloPy mapping ``{per: [(cellid, head, [bname]), ...]}``.
    name
        Package name used for ``pname`` and the output filename.
    boundnames
        Enable named boundaries (default ``False``).
    **options
        Extra keyword options passed straight to the FloPy constructor.
    """

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
    """Return a DRN (drain) package spec from raw FloPy stress-period data.

    The data form behind ``mf.drn.flopy(...)``. Wraps
    :class:`flopy.mf6.ModflowGwfdrn` with ``save_flows`` enabled.

    Parameters
    ----------
    stress_period_data
        FloPy mapping ``{per: [(cellid, elev, cond, [bname]), ...]}``.
    name
        Package name used for ``pname`` and the output filename.
    boundnames
        Enable named boundaries (default ``False``).
    **options
        Extra keyword options passed straight to the FloPy constructor.
    """

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
    """Return a GHB (general-head boundary) package spec from raw FloPy data.

    The data form behind ``mf.ghb.flopy(...)``. Wraps
    :class:`flopy.mf6.ModflowGwfghb` with ``save_flows`` enabled.

    Parameters
    ----------
    stress_period_data
        FloPy mapping ``{per: [(cellid, bhead, cond, [aux...], [bname]), ...]}``.
    name
        Package name used for ``pname`` and the output filename.
    auxiliary
        Optional auxiliary variable name(s).
    boundnames
        Enable named boundaries (default ``False``).
    **options
        Extra keyword options passed straight to the FloPy constructor.
    """

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
    """Return a list-based RCH (recharge) package spec from raw FloPy data.

    The data form behind ``mf.rch.flopy(...)``. Builds the list (cell-by-cell)
    form of :class:`flopy.mf6.ModflowGwfrch`; ``maxbound`` is inferred from the
    largest period's record count. (For the array form, or recharge derived from
    GIS/PRISM rasters, use the ``mf.rch(...)`` / ``RCHBuilder`` high-level path.)

    Parameters
    ----------
    stress_period_data
        FloPy mapping ``{per: [(cellid, recharge, [bname]), ...]}``.
    name
        Package name used for ``pname`` and the output filename.
    boundnames
        Enable named boundaries (default ``False``).
    **options
        Extra keyword options passed straight to the FloPy constructor.
    """

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
    """Return a UZF package spec from already-prepared UZF cell data.

    The raw-data form behind ``mf.uzf.flopy(...)``. Expects the FloPy UZF
    ``packagedata`` (one record per UZF cell, including vertical connectivity)
    and ``perioddata`` already assembled -- unlike ``mf.uzf(...)`` /
    ``UZFBuilder``, which derive cell records from a grid and infiltration inputs
    for you. Sets ``nuzfcells`` (defaulting to ``len(packagedata)``), enables
    ``save_flows``, and wires the standard UZF budget/convergence output records.

    Parameters
    ----------
    packagedata
        FloPy UZF packagedata records (one per UZF cell).
    perioddata
        FloPy UZF stress-period data (infiltration, PET, extinction depth, ...).
    nuzfcells
        Cell count; defaults to ``len(packagedata)``.
    name
        Package name used for ``pname`` and the output filenames.
    mover
        Enable MVR participation (rejected infiltration as a provider).
    simulate_et
        Enable evapotranspiration simulation.
    **options
        Extra keyword options passed straight to the FloPy constructor.
    """

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
    """Return a LAK package spec from already-prepared lake data.

    The raw-data form behind ``mf.lak.flopy(...)``. Expects the FloPy LAK
    ``packagedata`` (one record per lake), ``connectiondata`` (lake-to-cell
    connections), and ``perioddata`` already assembled -- unlike ``mf.lak(...)`` /
    ``LAKBuilder``, which resolve lake-cell connections from a grid and lake
    geometry for you. Sets ``nlakes`` (defaulting to ``len(packagedata)``),
    enables ``save_flows``, and wires the standard stage/budget/convergence
    output records.

    Parameters
    ----------
    packagedata
        FloPy LAK packagedata records (one per lake).
    connectiondata
        FloPy LAK connection records (lake-to-cell).
    perioddata
        FloPy LAK stress-period data (status, stage, withdrawals, ...).
    nlakes
        Lake count; defaults to ``len(packagedata)``.
    name
        Package name used for ``pname`` and the output filenames.
    noutlets, ntables
        Number of lake outlets and lake-bathymetry tables.
    mover
        Enable MVR participation (lake outflow as a provider/receiver).
    **options
        Extra keyword options passed straight to the FloPy constructor.
    """

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
    """Return an SFR package spec from already-prepared stream-reach data.

    The raw-data form behind ``mf.sfr.flopy(...)``. Expects the FloPy SFR
    ``packagedata`` (one record per reach), ``connectiondata`` (reach-to-reach
    topology), and ``perioddata`` already assembled -- unlike ``mf.sfr(...)`` /
    ``SFRBuilder``, which build reaches and connectivity from a stream centerline
    on a grid for you. Sets ``nreaches`` (defaulting to ``len(packagedata)``),
    enables ``save_flows``, and wires the stage/budget output records.

    Parameters
    ----------
    packagedata
        FloPy SFR packagedata records (one per reach).
    connectiondata
        FloPy SFR reach-to-reach connection records.
    perioddata
        FloPy SFR stress-period data (inflow, runoff, status, ...).
    nreaches
        Reach count; defaults to ``len(packagedata)``.
    name
        Package name used for ``pname`` and the output filenames.
    mover
        Enable MVR participation (reach outflow as a provider/receiver).
    diversions
        Optional SFR diversion records.
    **options
        Extra keyword options passed straight to the FloPy constructor.
    """

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
    """Return an MVR (water mover) spec, validating its package dependencies.

    The raw-data form behind ``mf.mvr.flopy(...)``. Takes the FloPy MVR
    ``packages`` declarations and ``perioddata`` mover records directly -- unlike
    ``mf.mvr(moves=...)`` / ``MVRBuilder``, which assemble these from semantic
    :class:`~myflopy.modflow.mf6.mvr.Move` objects. Validates that package names
    are unique, that every mover record references a declared package, and that
    ``maxmvr``/``maxpackages`` are large enough (both inferred when omitted). The
    returned spec carries ``requires=`` so the run ordering check ensures the
    moved packages are built before the mover.

    Parameters
    ----------
    packages
        MVR package declarations, e.g. ``[["sfr"], ["lak"]]``.
    perioddata
        Mover records per period: ``{per: [(pname1, id1, pname2, id2, mvrtype,
        value), ...]}``.
    name
        Package name used for ``pname`` and the output filenames.
    maxmvr, maxpackages
        Upper bounds; inferred from ``perioddata``/``packages`` when ``None``.
    **options
        Extra keyword options passed straight to the FloPy constructor.

    Raises
    ------
    ValueError
        On duplicate package names, references to undeclared packages, or
        ``maxmvr``/``maxpackages`` smaller than required.
    """

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

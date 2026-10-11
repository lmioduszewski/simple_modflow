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

import functools
from collections.abc import Iterable
from typing import Any

import flopy

from myflopy.grid_identity import CELLWISE_OPTIONS, INDEXED_PACKAGES, cellwise_options, pin_metadata
from myflopy.modflow.mf6.package_registry import (
    _PACKAGE_EXPLORER_SPECS,
    advanced_output_filerecords,
)
from myflopy.specs import ModelContext, PackageSpec


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
    """A picklable builder that constructs ``constructor``'s package (for use in a ``PackageSpec``)."""

    # functools.partial over a module-level function (not a lambda) so the
    # resulting PackageSpec can be pickled and reused across variants/sessions.
    return functools.partial(_build_named, constructor)



def _names_cells(options, ncpl: int) -> tuple[str, ...]:
    """WHICH resolved options are per-COLUMN arrays, i.e. name cells.

    The ncpl-EXACT spelling of :func:`~myflopy.grid_identity.cellwise_options`, kept as
    a name because `package_api`'s lazy areal pin imports it. The rule itself moved to
    `grid_identity` (import layer 0) when `package_move` needed the same one: npf/ic/sto
    carry per-cell arrays too and nothing was measuring them (ledger 295).
    """

    return cellwise_options(options, ncpl)


def _cellwise_options(options) -> tuple[str, ...]:
    """WHICH options are per-column when there is no ``ncpl`` to measure them against.

    The conservative substitute, valid HERE and not in general: under ``READASARRAYS``
    "one value per column" is the format's own definition, so any array-valued option is
    per-column by construction. It cannot know the length is wrong, only that the value
    addresses cells at all -- which is why `package_move` uses the exact form instead.
    """

    return cellwise_options(options)


def _array_bc_spec(
    constructor, *, name: str, options, context: ModelContext | None = None
) -> PackageSpec:
    """Build one array-form (``READASARRAYS``) boundary spec.

    Sibling of :func:`_list_bc_spec`, for the packages MODFLOW 6 offers in both
    shapes. The list form addresses cells one record at a time; the array form
    supplies one value per column per period and lets MODFLOW 6 resolve which
    layer receives it. Same physics, different input.

    ``readasarrays`` is FloPy's default on ``ModflowGwfrcha``/``ModflowGwfevta``
    and is deliberately not passed -- the class IS the choice.

    Goes through ``_factory``/``_named_options`` like everything else here, so
    ``pname``/``filename`` templating and picklability come for free.

    **The pin is CONDITIONAL**, on whether the resolved payload names cells rather
    than on whether a ``context=`` was passed. `mf.rch.array(recharge=0.001,
    irch=0)` resolves to two scalars and is genuinely grid-free -- pinning it would
    refuse a correct move. `irch="top_active"`, the default, resolves to an
    `(ncpl,)` array derived from the domain, and a same-`ncpl` mesh swap then
    applies the wrong per-cell layer indices with nothing said (ledger 267; ledger
    264 claimed this form was already pinned and it never was).
    """

    values = _named_options(name, dict(options))
    spec = PackageSpec(name, _factory(constructor), values)
    grid = getattr(context, "grid", None) if context is not None else None
    ncpl = getattr(grid, "ncpl", None)
    cellwise = (
        _names_cells(values, int(ncpl)) if ncpl is not None else _cellwise_options(values)
    )
    if not cellwise:
        return spec
    stamped = pin_metadata(context, certifies=cellwise) if ncpl is not None else {}
    # A pin when there was a mesh to pin TO, and the weaker marker otherwise. Both say
    # the same thing -- these options address cells -- but only the pin can say WHICH
    # mesh, so an unpinnable one classifies `unverifiable` rather than `unchanged`
    # (ledger 282). `pin_metadata` returns {} for a grid it cannot fingerprint, and
    # that falls through to the marker rather than losing it.
    return spec.with_metadata(**(stamped or {CELLWISE_OPTIONS: sorted(cellwise)}))


def _list_bc_spec(
    package: str,
    stress_period_data,
    *,
    name: str,
    boundnames: bool,
    **options,
) -> PackageSpec:
    """Build one cell-based list-BC spec from its registry descriptor.

    The seven ``*_spec`` factories below are thin wrappers over this. They keep
    their own names, signatures and docstrings -- those are the public,
    IDE-facing surface and each documents a genuinely different boundary -- but
    the body is here, once.

    That is the point: the bodies had drifted before. ``evt_spec`` and
    ``rch_spec`` each grew a ``maxbound`` that FloPy already computes from
    ``stress_period_data`` at write time, which crashed on any period whose
    record count differed. With one body, "we never set maxbound" is a property
    of the shared implementation rather than something seven functions have to
    independently keep remembering.

    ``options`` is splatted last, so a caller can still override anything --
    including ``save_flows``.
    """

    return PackageSpec(
        name,
        _factory(getattr(flopy.mf6, _PACKAGE_EXPLORER_SPECS[package].flopy_class)),
        _named_options(
            name,
            {
                "stress_period_data": stress_period_data,
                "boundnames": boundnames,
                "save_flows": True,
                **options,
            },
        ),
    )


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
    stress_period_data : dict or list
        FloPy stress-period data, handed to FloPy unchanged: a
        ``{period: [record, ...]}`` mapping with zero-based period keys, or a bare
        list of records, which FloPy files under period 0. Each record is
        ``(cellid, q[, aux...][, boundname])``: ``cellid`` is zero-based --
        ``(layer, cell)`` on DISV, ``(layer, row, col)`` on DIS -- and ``q`` is the
        well rate in model length^3/time, negative for extraction. A period mapped
        to an empty list is written as an empty ``PERIOD`` block. ``maxbound`` is
        never set here; FloPy computes it from these records at write time.
    name : str, default "wel"
        Package name: FloPy's ``pname``, the input file ``<model>.<name>`` (the model
        name is filled in at build), and the returned spec's name.
    auxiliary : str or list of str, optional
        Auxiliary variable name(s), written to MF6's ``AUXILIARY`` option. Each
        record then carries one value per name, after ``q``. ``None`` (the default)
        declares none.
    boundnames : bool, default True
        Write MF6's ``BOUNDNAMES`` option, declaring that records may end with a
        boundary name. A record with no name is still accepted (FloPy stores an
        empty one).
    **options
        Any other keyword of :class:`flopy.mf6.ModflowGwfwel` (``print_input``,
        ``auto_flow_reduce``, ``timeseries``, ``observations``, ...), passed to it
        at build. The names are checked against that class's signature when the
        spec is created, and an unknown one raises ``TypeError`` naming a close
        valid one if any. They are applied last, so they override ``save_flows=True`` and
        the ``pname``/``filename`` set from ``name``. A string value may contain
        ``{model_name}``, which is replaced by the model's name at build.
    """

    return _list_bc_spec(
        "wel",
        stress_period_data,
        name=name,
        boundnames=boundnames,
        auxiliary=auxiliary,
        **options,
    )


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
    stress_period_data : dict or list
        FloPy stress-period data, handed to FloPy unchanged: a
        ``{period: [record, ...]}`` mapping with zero-based period keys, or a bare
        list of records, which FloPy files under period 0. Each record is
        ``(cellid, head[, aux...][, boundname])``: ``cellid`` is zero-based --
        ``(layer, cell)`` on DISV, ``(layer, row, col)`` on DIS -- and ``head`` is
        the specified head, in model length units. A period mapped to an empty
        list is written as an empty ``PERIOD`` block. ``maxbound`` is never set
        here; FloPy computes it from these records at write time.
    name : str, default "chd"
        Package name: FloPy's ``pname``, the input file ``<model>.<name>`` (the model
        name is filled in at build), and the returned spec's name.
    boundnames : bool, default False
        Write MF6's ``BOUNDNAMES`` option, declaring that records may end with a
        boundary name. Leave it ``False`` when the records carry none.
    **options
        Any other keyword of :class:`flopy.mf6.ModflowGwfchd` (``auxiliary``,
        ``print_input``, ``timeseries``, ``observations``, ...), passed to it at
        build. The names are checked against that class's signature when the spec
        is created, and an unknown one raises ``TypeError`` naming a close valid
        one if any. They are applied last, so they override ``save_flows=True`` and the
        ``pname``/``filename`` set from ``name``. A string value may contain
        ``{model_name}``, which is replaced by the model's name at build.
    """

    return _list_bc_spec(
        "chd",
        stress_period_data,
        name=name,
        boundnames=boundnames,
        **options,
    )


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
    stress_period_data : dict or list
        FloPy stress-period data, handed to FloPy unchanged: a
        ``{period: [record, ...]}`` mapping with zero-based period keys, or a bare
        list of records, which FloPy files under period 0. Each record is
        ``(cellid, elev, cond[, aux...][, boundname])``: ``cellid`` is zero-based
        -- ``(layer, cell)`` on DISV, ``(layer, row, col)`` on DIS -- ``elev`` is
        the drain elevation (length) and ``cond`` its conductance
        (length^2/time). A period mapped to an empty list is written as an empty
        ``PERIOD`` block. ``maxbound`` is never set here; FloPy computes it from
        these records at write time.
    name : str, default "drn"
        Package name: FloPy's ``pname``, the input file ``<model>.<name>`` (the model
        name is filled in at build), and the returned spec's name.
    boundnames : bool, default False
        Write MF6's ``BOUNDNAMES`` option, declaring that records may end with a
        boundary name. Leave it ``False`` when the records carry none.
    **options
        Any other keyword of :class:`flopy.mf6.ModflowGwfdrn` (``auxiliary``,
        ``mover``, ``print_input``, ``timeseries``, ``observations``, ...), passed
        to it at build. The names are checked against that class's signature when
        the spec is created, and an unknown one raises ``TypeError`` naming a
        close valid one if any. They are applied last, so they override
        ``save_flows=True`` and the ``pname``/``filename`` set from ``name``. A
        string value may contain ``{model_name}``, which is replaced by the model's
        name at build.
    """

    return _list_bc_spec(
        "drn",
        stress_period_data,
        name=name,
        boundnames=boundnames,
        **options,
    )


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
    stress_period_data : dict or list
        FloPy stress-period data, handed to FloPy unchanged: a
        ``{period: [record, ...]}`` mapping with zero-based period keys, or a bare
        list of records, which FloPy files under period 0. Each record is
        ``(cellid, bhead, cond[, aux...][, boundname])``: ``cellid`` is zero-based
        -- ``(layer, cell)`` on DISV, ``(layer, row, col)`` on DIS -- ``bhead`` is
        the boundary head (length) and ``cond`` the conductance (length^2/time). A
        period mapped to an empty list is written as an empty ``PERIOD`` block.
        ``maxbound`` is never set here; FloPy computes it from these records at
        write time.
    name : str, default "ghb"
        Package name: FloPy's ``pname``, the input file ``<model>.<name>`` (the model
        name is filled in at build), and the returned spec's name.
    auxiliary : str or list of str, optional
        Auxiliary variable name(s), written to MF6's ``AUXILIARY`` option. Each
        record then carries one value per name, after ``cond``. ``None`` (the
        default) declares none.
    boundnames : bool, default False
        Write MF6's ``BOUNDNAMES`` option, declaring that records may end with a
        boundary name. Leave it ``False`` when the records carry none.
    **options
        Any other keyword of :class:`flopy.mf6.ModflowGwfghb` (``mover``,
        ``print_input``, ``timeseries``, ``observations``, ...), passed to it at
        build. The names are checked against that class's signature when the spec
        is created, and an unknown one raises ``TypeError`` naming a close valid
        one if any. They are applied last, so they override ``save_flows=True`` and the
        ``pname``/``filename`` set from ``name``. A string value may contain
        ``{model_name}``, which is replaced by the model's name at build.
    """

    return _list_bc_spec(
        "ghb",
        stress_period_data,
        name=name,
        boundnames=boundnames,
        auxiliary=auxiliary,
        **options,
    )


def riv_spec(
    stress_period_data,
    *,
    name: str = "riv",
    auxiliary=None,
    boundnames: bool = False,
    **options,
) -> PackageSpec:
    """Return a RIV (river) package spec from raw FloPy stress-period data.

    The data form behind ``mf.riv.flopy(...)``. Wraps
    :class:`flopy.mf6.ModflowGwfriv` with ``save_flows`` enabled.

    Parameters
    ----------
    stress_period_data : dict or list
        FloPy stress-period data, handed to FloPy unchanged: a
        ``{period: [record, ...]}`` mapping with zero-based period keys, or a bare
        list of records, which FloPy files under period 0. Each record is
        ``(cellid, stage, cond, rbot[, aux...][, boundname])``: ``cellid`` is
        zero-based -- ``(layer, cell)`` on DISV, ``(layer, row, col)`` on DIS --
        ``stage`` is the river stage (length), ``cond`` the riverbed conductance
        (length^2/time) and ``rbot`` the river-bottom elevation (length). A period
        mapped to an empty list is written as an empty ``PERIOD`` block.
        ``maxbound`` is never set here; FloPy computes it from these records at
        write time.
    name : str, default "riv"
        Package name: FloPy's ``pname``, the input file ``<model>.<name>`` (the model
        name is filled in at build), and the returned spec's name.
    auxiliary : str or list of str, optional
        Auxiliary variable name(s), written to MF6's ``AUXILIARY`` option. Each
        record then carries one value per name, after ``rbot``. ``None`` (the
        default) declares none.
    boundnames : bool, default False
        Write MF6's ``BOUNDNAMES`` option, declaring that records may end with a
        boundary name. Leave it ``False`` when the records carry none.
    **options
        Any other keyword of :class:`flopy.mf6.ModflowGwfriv` (``mover``,
        ``print_input``, ``timeseries``, ``observations``, ...), passed to it at
        build. The names are checked against that class's signature when the spec
        is created, and an unknown one raises ``TypeError`` naming a close valid
        one if any. They are applied last, so they override ``save_flows=True`` and the
        ``pname``/``filename`` set from ``name``. A string value may contain
        ``{model_name}``, which is replaced by the model's name at build.
    """

    return _list_bc_spec(
        "riv",
        stress_period_data,
        name=name,
        boundnames=boundnames,
        auxiliary=auxiliary,
        **options,
    )


def rch_spec(
    stress_period_data,
    *,
    name: str = "rch",
    boundnames: bool = False,
    **options,
) -> PackageSpec:
    """Return a list-based RCH (recharge) package spec from raw FloPy data.

    The data form behind ``mf.rch.flopy(...)``. Builds the list (cell-by-cell)
    form of :class:`flopy.mf6.ModflowGwfrch`. ``maxbound`` is left to FloPy,
    which computes it from ``stress_period_data`` at write time -- so every
    native input shape (a dict of record lists, a bare list, periods set to
    ``None``) passes through verbatim. For the ARRAY form use
    ``mf.rch.array(...)``; for recharge derived from GIS/PRISM rasters use the
    ``mf.rch(...)`` / ``RCHBuilder`` high-level path, which also produces the
    list form.

    Parameters
    ----------
    stress_period_data : dict or list
        FloPy stress-period data, handed to FloPy unchanged: a
        ``{period: [record, ...]}`` mapping with zero-based period keys, or a bare
        list of records, which FloPy files under period 0. Each record is
        ``(cellid, recharge[, aux...][, boundname])``: ``cellid`` is zero-based --
        ``(layer, cell)`` on DISV, ``(layer, row, col)`` on DIS -- and ``recharge``
        is a flux in model length/time. A period mapped to an empty list is written
        as an empty ``PERIOD`` block.
    name : str, default "rch"
        Package name: FloPy's ``pname``, the input file ``<model>.<name>`` (the model
        name is filled in at build), and the returned spec's name.
    boundnames : bool, default False
        Write MF6's ``BOUNDNAMES`` option, declaring that records may end with a
        boundary name. Leave it ``False`` when the records carry none.
    **options
        Any other keyword of :class:`flopy.mf6.ModflowGwfrch` (``auxiliary``,
        ``fixed_cell``, ``print_input``, ``timeseries``, ``observations``, ...),
        passed to it at build. The names are checked against that class's
        signature when the spec is created, and an unknown one raises ``TypeError``
        naming a close valid one if any. They are applied last, so they override
        ``save_flows=True`` and the ``pname``/``filename`` set from ``name``. A
        string value may contain ``{model_name}``, which is replaced by the model's
        name at build.
    """

    return _list_bc_spec(
        "rch",
        stress_period_data,
        name=name,
        boundnames=boundnames,
        **options,
    )


def evt_spec(
    stress_period_data,
    *,
    name: str = "evt",
    nseg: int = 1,
    auxiliary=None,
    boundnames: bool = False,
    **options,
) -> PackageSpec:
    """Return a list-based EVT (evapotranspiration) package spec from raw FloPy data.

    The data form behind ``mf.evt.flopy(...)``. Builds the list (cell-by-cell)
    form of :class:`flopy.mf6.ModflowGwfevt`. ``maxbound`` is left to FloPy,
    which computes it from ``stress_period_data`` at write time -- so every
    native input shape (a dict of record lists, a bare list, periods set to
    ``None``) passes through verbatim. With the default ``nseg=1`` each record is
    ``(cellid, surface, rate, depth)``; for segmented ET pass ``nseg`` and
    include the ``pxdp``/``petm`` (and, with ``surf_rate_specified``, ``petm0``)
    record values FloPy expects.

    Parameters
    ----------
    stress_period_data : dict or list
        FloPy stress-period data, handed to FloPy unchanged: a
        ``{period: [record, ...]}`` mapping with zero-based period keys, or a bare
        list of records, which FloPy files under period 0. Each record is
        ``(cellid, surface, rate, depth[, pxdp..., petm...][, petm0][, aux...]
        [, boundname])``: ``cellid`` is zero-based -- ``(layer, cell)`` on DISV,
        ``(layer, row, col)`` on DIS -- ``surface`` is the ET surface elevation
        (length), ``rate`` the maximum ET flux (length/time) and ``depth`` the
        extinction depth below the surface (length). A period mapped to an empty
        list is written as an empty ``PERIOD`` block.
    name : str, default "evt"
        Package name: FloPy's ``pname``, the input file ``<model>.<name>`` (the model
        name is filled in at build), and the returned spec's name.
    nseg : int, default 1
        MF6's ``NSEG``, the number of ET segments, passed to FloPy. ``1`` is the
        plain linear decay from ``rate`` at the surface to zero at ``depth``, and
        its records carry no ``pxdp``/``petm`` values; with ``nseg > 1`` every
        record must carry ``nseg - 1`` of each.
    auxiliary : str or list of str, optional
        Auxiliary variable name(s), written to MF6's ``AUXILIARY`` option. Each
        record then carries one value per name, before any boundary name. ``None``
        (the default) declares none.
    boundnames : bool, default False
        Write MF6's ``BOUNDNAMES`` option, declaring that records may end with a
        boundary name. Leave it ``False`` when the records carry none.
    **options
        Any other keyword of :class:`flopy.mf6.ModflowGwfevt` (``fixed_cell``,
        ``surf_rate_specified``, ``print_input``, ``timeseries``, ``observations``,
        ...), passed to it at build. The names are checked against that class's
        signature when the spec is created, and an unknown one raises ``TypeError``
        naming a close valid one if any. They are applied last, so they override
        ``save_flows=True`` and the ``pname``/``filename`` set from ``name``. A
        string value may contain ``{model_name}``, which is replaced by the model's
        name at build.
    """

    return _list_bc_spec(
        "evt",
        stress_period_data,
        name=name,
        boundnames=boundnames,
        nseg=nseg,
        auxiliary=auxiliary,
        **options,
    )


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
    packagedata : list
        FloPy UZF ``packagedata`` records, one per UZF object and handed to FloPy
        unchanged: ``(ifno, cellid, landflag, ivertcon, surfdep, vks, thtr, thts,
        thti, eps[, boundname])``, with zero-based ``ifno``, ``cellid`` and
        ``ivertcon`` (FloPy adds one on write, so an ``ivertcon`` of ``-1`` writes
        MF6's ``0``, no underlying UZF object). Its length is the default
        ``nuzfcells``.
    perioddata : dict
        FloPy UZF ``perioddata``, handed to FloPy unchanged:
        ``{period: [(ifno, finf, pet, extdp, extwc, ha, hroot, rootact[, aux...]),
        ...]}`` with zero-based period keys -- infiltration and potential ET in
        length/time, extinction depth in length.
    nuzfcells : int, optional
        MF6's ``NUZFCELLS``. ``None`` (the default) uses ``len(packagedata)``; it is
        not checked against ``packagedata`` when given.
    name : str, default "uzf"
        Package name: FloPy's ``pname``, the input file ``<model>.<name>``, and the
        returned spec's name. It also names the output files written here:
        ``<model>_<name>_budget.uzf``, ``<model>_<name>_budget.csv`` and
        ``<model>_<name>_package_convergence.csv``.
    mover : bool, default False
        Write MF6's ``MOVER`` option, which the package needs before an MVR
        package (``mf.mvr``) can move water to or from it.
    simulate_et : bool, default False
        Write MF6's ``SIMULATE_ET`` option, which turns on ET in the unsaturated
        zone -- and in the GWF cell as well when ``linear_gwet`` or ``square_gwet``
        is also passed through ``**options``.
    **options
        Any other keyword of :class:`flopy.mf6.ModflowGwfuzf` (``simulate_gwseep``,
        ``linear_gwet``, ``ntrailwaves``, ``nwavesets``, ``wc_filerecord``,
        ``observations``, ...), passed to it at build. The names are checked
        against that class's signature when the spec is created, and an unknown one
        raises ``TypeError`` naming a close valid one if any. They are applied last,
        so they override ``save_flows=True``, the ``pname``/``filename`` set from
        ``name`` and the three output file records. A string value may contain
        ``{model_name}``, which is replaced by the model's name at build.
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
            **advanced_output_filerecords("uzf", name),
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
    packagedata : list
        FloPy LAK ``packagedata`` records, one per lake and handed to FloPy
        unchanged: ``(ifno, strt, nlakeconn[, aux...][, boundname])`` -- a
        zero-based lake number, the starting stage (length) and the number of
        cells connected to that lake. Its length is the default ``nlakes``.
    connectiondata : list
        FloPy LAK ``connectiondata`` records, one per lake-cell connection and
        handed to FloPy unchanged: ``(ifno, iconn, cellid, claktype, bedleak,
        belev, telev, connlen, connwidth)``, with zero-based ``ifno``, ``iconn`` and
        ``cellid``.
    perioddata : dict
        FloPy LAK ``perioddata``, handed to FloPy unchanged: ``{period: [(number,
        laksetting...), ...]}`` with zero-based period keys, where each record
        starts with a lake (or outlet) number followed by a setting keyword and its
        value(s) -- ``"status"``, ``"stage"``, ``"rainfall"``, ``"evaporation"``,
        ``"withdrawal"``, ...
    nlakes : int, optional
        MF6's ``NLAKES``. ``None`` (the default) uses ``len(packagedata)``; it is not
        checked against ``packagedata`` when given.
    name : str, default "lak"
        Package name: FloPy's ``pname``, the input file ``<model>.<name>``, and the
        returned spec's name. It also names the output files written here:
        ``<model>_<name>_stage.lak``, ``<model>_<name>_budget.lak``,
        ``<model>_<name>_budget.csv`` and ``<model>_<name>_package_convergence.csv``.
    noutlets : int, default 0
        MF6's ``NOUTLETS``, the number of lake outlets. With a nonzero value, pass
        the outlet records through ``outlets=`` in ``**options``.
    ntables : int, default 0
        MF6's ``NTABLES``, the number of lakes given a stage-volume-area table. With
        a nonzero value, pass the table references through ``tables=`` in
        ``**options``.
    mover : bool, default False
        Write MF6's ``MOVER`` option, which the package needs before an MVR
        package (``mf.mvr``) can move water to or from it.
    **options
        Any other keyword of :class:`flopy.mf6.ModflowGwflak` (``outlets``,
        ``tables``, ``surfdep``, ``length_conversion``, ``time_conversion``,
        ``observations``, ...), passed to it at build. The names are checked
        against that class's signature when the spec is created, and an unknown one
        raises ``TypeError`` naming a close valid one if any. They are applied last,
        so they override ``save_flows=True``, the ``pname``/``filename`` set from
        ``name`` and the four output file records. A string value may contain
        ``{model_name}``, which is replaced by the model's name at build.
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
            **advanced_output_filerecords("lak", name),
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
    packagedata : list
        FloPy SFR ``packagedata`` records, one per reach and handed to FloPy
        unchanged: ``(ifno, cellid, rlen, rwid, rgrd, rtp, rbth, rhk, man, ncon,
        ustrf, ndv[, aux...][, boundname])``, with a zero-based reach number and
        ``cellid``. Its length is the default ``nreaches``.
    connectiondata : list
        FloPy SFR ``connectiondata`` records, handed to FloPy unchanged:
        ``(ifno, ic...)`` -- a reach followed by the reaches it connects to, as
        zero-based reach numbers, negative for a downstream connection (FloPy adds
        one to the magnitude on write, so ``(0, -1)`` writes ``1 -2``).
    perioddata : dict
        FloPy SFR ``perioddata``, handed to FloPy unchanged: ``{period: [(ifno,
        sfrsetting...), ...]}`` with zero-based period keys, where each record
        starts with a reach number followed by a setting keyword and its value(s)
        -- ``"status"``, ``"inflow"``, ``"rainfall"``, ``"evaporation"``,
        ``"runoff"``, ...
    nreaches : int, optional
        MF6's ``NREACHES``. ``None`` (the default) uses ``len(packagedata)``; it is
        not checked against ``packagedata`` when given.
    name : str, default "sfr"
        Package name: FloPy's ``pname``, the input file ``<model>.<name>``, and the
        returned spec's name. It also names the output files written here:
        ``<model>_<name>_stage.sfr`` and ``<model>_<name>_budget.sfr``.
    mover : bool, default False
        Write MF6's ``MOVER`` option, which the package needs before an MVR
        package (``mf.mvr``) can move water to or from it.
    diversions : list, optional
        FloPy SFR ``diversions`` records, handed to FloPy unchanged:
        ``(ifno, idv, iconr, cprior)``. ``None`` (the default) writes none; a reach
        with diversions must also declare their count as ``ndv`` in
        ``packagedata``.
    **options
        Any other keyword of :class:`flopy.mf6.ModflowGwfsfr` (``crosssections``,
        ``initialstages``, ``length_conversion``, ``time_conversion``,
        ``observations``, ...), passed to it at build. The names are checked
        against that class's signature when the spec is created, and an unknown one
        raises ``TypeError`` naming a close valid one if any. They are applied last,
        so they override ``save_flows=True``, the ``pname``/``filename`` set from
        ``name`` and the two output file records. A string value may contain
        ``{model_name}``, which is replaced by the model's name at build.
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
            **advanced_output_filerecords("sfr", name),
            **options,
        },
    )
    return PackageSpec(name, _factory(flopy.mf6.ModflowGwfsfr), values)


def _package_names(packages: Iterable[Any]) -> tuple[str, ...]:
    """Return the package names from MVR package records, refusing a bare string.

    A record is a list or tuple whose first item is the package name. A bare string was
    taken as its own name here while ``str(record[0])`` elsewhere took its first LETTER --
    ``["sfr", "lak"]`` required ``("sfr", "lak")`` and stamped ``["l", "s"]`` as the
    packages the mover addresses by index, so a grid move missed them (ledger 303) -- and
    FloPy reads the same flat list as ONE two-field record and refuses it at build. One
    spelling, named in the refusal (ledger 381).
    """

    if isinstance(packages, (str, bytes)) or any(
        not isinstance(record, (list, tuple)) for record in packages
    ):
        flat = [packages] if isinstance(packages, (str, bytes)) else list(packages)
        spelling = [
            list(record) if isinstance(record, (list, tuple)) else [str(record)]
            for record in flat
        ]
        raise ValueError(
            "MVR packages are one single-item list or tuple per package, not bare names: "
            f"packages={spelling!r}. FloPy reads a flat list of names as ONE record."
        )
    return tuple(str(record[0]) for record in packages)


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
    packages : list
        MVR package declarations, one single-item list or tuple per package that
        provides or receives water: ``[["sfr"], ["lak"]]``. The item is the
        package name, and it is what the returned spec ``requires=`` -- so it must
        match a ``PackageSpec`` name declared earlier in the same model -- and what
        it records as the packages it addresses by index. Names must be unique. A
        bare name (``["sfr", "lak"]``) raises ``ValueError`` naming the record
        spelling: FloPy reads a flat list as ONE record.
    perioddata : dict
        Mover records per period, ``{period: [(pname1, id1, pname2, id2, mvrtype,
        value), ...]}`` with zero-based period keys (a list is not accepted).
        ``pname1`` provides and ``pname2`` receives; both must appear in
        ``packages``. ``id1``/``id2`` are zero-based feature numbers (reach, lake,
        UZF object; FloPy adds one on write), ``mvrtype`` is ``"FACTOR"``,
        ``"EXCESS"``, ``"THRESHOLD"`` or ``"UPTO"``, and ``value`` is its factor or
        rate. A period mapped to an empty list is written as an empty ``PERIOD``
        block.
    name : str, default "mvr"
        Package name: FloPy's ``pname``, the input file ``<model>.<name>``, and the
        returned spec's name. It also names the output files written here:
        ``<model>.<name>.bud`` and ``<model>.<name>.csv``.
    maxmvr : int, optional
        MF6's ``MAXMVR``. ``None`` (the default) uses the largest number of records
        in any one period of ``perioddata``; a smaller value raises ``ValueError``.
    maxpackages : int, optional
        MF6's ``MAXPACKAGES``. ``None`` (the default) uses ``len(packages)``; a
        smaller value raises ``ValueError``.
    **options
        Any other keyword of :class:`flopy.mf6.ModflowGwfmvr` (``print_input``,
        ``print_flows``, ...), passed to it at build. The names are checked against
        that class's signature when the spec is created, and an unknown one raises
        ``TypeError`` naming a close valid one if any. They are applied last, so they
        override the ``pname``/``filename`` set from ``name`` and the two output
        file records. A string value may contain ``{model_name}``, which is
        replaced by the model's name at build. ``modelnames=True`` raises
        ``ValueError``: this builds a MODEL's mover, and MF6 refuses ``MODELNAMES``
        there in either record layout ("MODELNAMES cannot be specified unless the
        mover package is for an exchange", ``docs/reproducers/mvr_modelnames.py``);
        ``modelnames=False`` is MF6's default and is passed on.

    Raises
    ------
    ValueError
        On ``modelnames=True``, a bare package name rather than a record,
        duplicate package names, references to undeclared packages, or
        ``maxmvr``/``maxpackages`` smaller than required.
    """

    if options.get("modelnames"):
        # Read as the single-model layout below, a (mname, pname) record's MODEL name
        # became the package `requires=` and INDEXED_PACKAGES named, and MF6 refuses
        # the option on a model's mover whatever the records say (ledger 413).
        raise ValueError(
            f"modelnames={options['modelnames']!r}: MF6 takes MODELNAMES only on the mover of a "
            "GWF-GWF exchange, and this builds a model's own mover, which MF6 refuses with it "
            "('MODELNAMES cannot be specified unless the mover package is for an exchange'). "
            "Leave modelnames out and name the packages alone: packages=[['sfr'], ['lak']], "
            "records (pname1, id1, pname2, id2, mvrtype, value)."
        )
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
        # Which siblings this package addresses by INDEX. Stamped HERE rather than in
        # `MVRBuilder.build` because `mf.mvr.flopy(...)` reaches this function directly
        # and would otherwise carry no marker at all. Not `requires=`, which already
        # means "build this first" -- reading that as an index reference would refuse a
        # package ordered after a discretization whenever it rebuilt (ledger 303). From
        # `required`, the names `_package_names` read, so the two can never disagree.
        metadata={INDEXED_PACKAGES: sorted(set(required))},
    )


__all__ = [
    "chd_spec",
    "drn_spec",
    "evt_spec",
    "ghb_spec",
    "lak_spec",
    "mvr_spec",
    "rch_spec",
    "riv_spec",
    "sfr_spec",
    "uzf_spec",
    "wel_spec",
]

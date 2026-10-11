"""Package-first helpers for building readable MODFLOW 6 package specs.

The objects in this module are thin facades over the existing builders and
``PackageSpec`` factories. They keep model scripts close to MODFLOW vocabulary:

``mf.drn(...)`` for direct DRN stress-period data, ``mf.drn.gpkg(...)`` for GIS
boundaries, and ``mf.uzf(...)`` for the higher-level UZF builder.
"""

from __future__ import annotations

import functools
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import replace as _replace
from pathlib import Path
from typing import Any, ClassVar

import flopy
import geopandas as gpd
import numpy as np
from shapely.geometry import Point, Polygon

from myflopy._logging import get_logger
from myflopy.advanced import (
    _array_bc_spec,
    _names_cells,
    chd_spec,
    drn_spec,
    evt_spec,
    ghb_spec,
    lak_spec,
    mvr_spec,
    rch_spec,
    riv_spec,
    sfr_spec,
    uzf_spec,
    wel_spec,
)
from myflopy.builders import (
    build_dis,
    build_disu,
    build_disv,
    build_ems,
    build_gwf_gwe_exchange,
    build_gwf_gwf_exchange,
    build_gwf_gwt_exchange,
    build_gwf_prt_exchange,
    build_ic,
    build_ims,
    build_oc,
)
from myflopy.geopackage import (
    CellSurfaceOffset,
    GeoPackageSource,
    RowValue,
    Spread,
    SurfaceResolver,
    array_from_polygons,
    refuse_spread_on_elevation,
)
from myflopy.grid_identity import (
    CELLWISE_OPTIONS,
    GENERATED_SOLVER,
    GRID_PIN,
    GRID_PIN_OPTIONS,
    PIN_VOIDED,
    SOURCE_SHA,
    disv_pin,
    pin_metadata,
    source_fingerprint,
)
from myflopy.modflow.mf6.evapotranspiration import EVTBuilder
from myflopy.modflow.mf6.grid.geometry import cell_shares, top_active_layer
from myflopy.modflow.mf6.grid.surfaces import check_gaps
from myflopy.modflow.mf6.lakes import (
    LAKBuilder,
    LakeConnection,
    LakeOutlet,
    LakeTable,
    LakeTableBuilder,
    _flopy_lak_refusals,
)
from myflopy.modflow.mf6.mvr import Move, MoverConnection, MVRBuilder
from myflopy.modflow.mf6.package_registry import _PACKAGE_EXPLORER_SPECS
from myflopy.modflow.mf6.recharge import RCHBuilder
from myflopy.modflow.mf6.sfr import SFRBuilder, StreamConnection, StreamDiversion
from myflopy.modflow.mf6.uzf import UZFBuilder
from myflopy.sources import BoundarySource, Polygons
from myflopy.specs import (
    _PUBLISHED_BY_DISV,
    _PUBLISHED_BY_KEY,
    ExchangeSpec,
    ModelContext,
    ModelSpec,
    PackageEntry,
    PackageSpec,
    PostBuildHook,
    SimulationSpec,
    _layer_elevations,
    stamp_source_hashes,
)

PathLike = Path | str

logger = get_logger(__name__)


def _model_options(
    *,
    model_nam_file: str | None = None,
    version: str = "mf6",
    exe_name: str = "mf6",
    model_rel_path: str | Path = ".",
    list: str | None = None,
    print_input: bool | None = None,
    print_flows: bool | None = None,
    save_flows: bool | None = None,
    newtonoptions: str | Sequence[str] | None = None,
    nc_mesh2d_filerecord: Any = None,
    nc_structured_filerecord: Any = None,
    nc_filerecord: Any = None,
    extra: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Collect common MODFLOW 6 model constructor options."""

    options = {
        "version": version,
        "exe_name": exe_name,
        "model_rel_path": model_rel_path,
        **({} if extra is None else extra),
    }
    for key, value in {
        "model_nam_file": model_nam_file,
        "list": list,
        "print_input": print_input,
        "print_flows": print_flows,
        "save_flows": save_flows,
        "nc_mesh2d_filerecord": nc_mesh2d_filerecord,
        "nc_structured_filerecord": nc_structured_filerecord,
        "nc_filerecord": nc_filerecord,
    }.items():
        if value is not None:
            options[key] = value
    return options


def gwf(
    name: str,
    *,
    packages: Iterable[PackageSpec] = (),
    context: ModelContext | None = None,
    hooks: Iterable[PostBuildHook] = (),
    grid: Any = None,
    model_nam_file: str | None = None,
    version: str = "mf6",
    exe_name: str = "mf6",
    model_rel_path: str | Path = ".",
    list: str | None = None,
    print_input: bool | None = None,
    print_flows: bool | None = None,
    save_flows: bool | None = None,
    newtonoptions: str | Sequence[str] | None = None,
    nc_mesh2d_filerecord: Any = None,
    nc_structured_filerecord: Any = None,
    nc_filerecord: Any = None,
    **kwargs: Any,
) -> ModelSpec:
    """Return a typed GWF (groundwater-flow) model spec wrapping ``flopy.mf6.ModflowGwf``.

    The primary model constructor of the package-first API. Discretization and
    physics (``mf.disv``, ``mf.npf``, ``mf.ic``, ``mf.sto``, ``mf.oc``) and boundary
    conditions (``mf.chd``, ``mf.ghb``, ``mf.rch``, ``mf.sfr``, ``mf.lak``, ...) are
    passed as ``PackageSpec`` objects through ``packages`` and built after the model
    is created. Geometry rides on the model via ``context`` (not the project).

    Parameters
    ----------
    name : str
        Model name -- the FloPy model name, and the key that IMS solvers
        (``mf.ims(models=(name,))``) and inter-model exchanges reference.
    packages : Iterable[PackageSpec or PackageRef or str], optional
        Package specs to attach, e.g. ``[mf.disv(...), mf.npf(...), mf.ic(...),
        mf.oc(...)]``; an ``mf.ref(key)`` or a bare key string (turned into a
        ``PackageRef``) names a package in the project's library. Built in this
        order after the model. Names must be unique (``ValueError`` otherwise),
        and when every entry is a concrete ``PackageSpec`` each one's ``requires``
        must be present and declared before it. Default ``()``: no packages.
    context : ModelContext, optional
        Geometry/dates for this model -- ``grid``, ``domain`` (idomain),
        ``surfaces``, ``dates``. Reachable afterward as ``model.myflopy_context``;
        GIS-aware helpers (``mf.rch``, ``mf.sfr``, ``mf.lak``, ``mf.X.gpkg``) read it.
        Default ``None`` stores an empty ``ModelContext()``.
    hooks : Iterable[PostBuildHook], optional
        Post-build callbacks, each called as ``callback(model, packages, context)``
        after every package is built; results land in ``BuiltModel.hook_results``
        under the hook's name. Hook names must be unique. Default ``()``.
    grid : GridSpec or GridRef or grid object, optional
        A deferred ``mf.GridSpec.voronoi(...)``, an ``mf.grid_ref(key)`` into the
        project's grid library, or an already-built grid (wrapped with
        ``GridSpec.from_object``). Resolved at build time, before any package
        builds, and it REPLACES ``context.grid`` on the built model. Default
        ``None``: the grid, if any, is ``context.grid``.
    model_nam_file : str, optional
        Name-file name for the model, forwarded to ``ModflowGwf``; omitted when
        ``None``, so FloPy names it after the model.
    version : str, optional
        FloPy's model ``version`` string. Default ``"mf6"``.
    exe_name : str, optional
        Executable name recorded on the FloPy model object. Default ``"mf6"``.
    model_rel_path : str or Path, optional
        The model's folder relative to the simulation workspace. Default ``"."``
        (the workspace itself).
    list : str, optional
        Listing-file name; omitted when ``None``, so MF6 uses the name file's
        base name with ``.lst``.
    print_input : bool, optional
        MF6 ``PRINT_INPUT``: echo package input to the listing file. Omitted
        when ``None`` (MF6's default, off).
    print_flows : bool, optional
        MF6 ``PRINT_FLOWS``: print package flow rates to the listing file.
        Omitted when ``None`` (off).
    save_flows : bool, optional
        MF6 ``SAVE_FLOWS``: write package flow terms to the budget file named
        by output control. Omitted when ``None`` (off).
    newtonoptions : str or Sequence[str], optional
        MF6 Newton-Raphson formulation, e.g. ``"NEWTON"`` or
        ``"NEWTON UNDER_RELAXATION"`` -- strongly recommended for water-table /
        wetting-and-drying models. Omitted when ``None``. A NEWTON model gets a
        BICGSTAB solver when ``mf.simulation`` generates the IMS.
    nc_mesh2d_filerecord : Any, optional
        FloPy's NetCDF layered-mesh output record, passed through when not
        ``None``.
    nc_structured_filerecord : Any, optional
        FloPy's NetCDF structured-grid output record, passed through when not
        ``None``.
    nc_filerecord : Any, optional
        FloPy's NetCDF file record, passed through when not ``None``.
    **kwargs
        Any additional ``flopy.mf6.ModflowGwf`` constructor option, stored on the
        spec and passed to ``ModflowGwf`` at build time. Not checked at the call:
        a name FloPy does not accept fails when the model is built.

    Returns
    -------
    ModelSpec
        A declarative model spec; add it to a :class:`SimulationSpec`
        (``mf.simulation(flow)`` or ``mf.SimulationSpec(...)``).

    Examples
    --------
    >>> ctx = mf.ModelContext(grid=vor, domain=idomain)
    >>> flow = mf.gwf("valley", context=ctx, newtonoptions="NEWTON",
    ...               packages=[mf.disv(...), mf.npf(k=10.0, icelltype=1),
    ...                         mf.ic(strt=100.0), mf.sto(steady_state={0: True}),
    ...                         mf.oc(head_filerecord="valley.hds",
    ...                               saverecord=[("HEAD", "ALL")])])
    >>> sim = mf.simulation(flow)
    """

    options = _model_options(
        model_nam_file=model_nam_file,
        version=version,
        exe_name=exe_name,
        model_rel_path=model_rel_path,
        list=list,
        print_input=print_input,
        print_flows=print_flows,
        save_flows=save_flows,
        nc_mesh2d_filerecord=nc_mesh2d_filerecord,
        nc_structured_filerecord=nc_structured_filerecord,
        nc_filerecord=nc_filerecord,
        extra=kwargs,
    )
    if newtonoptions is not None:
        options["newtonoptions"] = newtonoptions

    return ModelSpec(
        name,
        "gwf",
        packages=tuple(packages),
        context=ModelContext() if context is None else context,
        hooks=tuple(hooks),
        grid=grid,
        options=options,
    )


def gwt(
    name: str,
    *,
    packages: Iterable[PackageSpec] = (),
    context: ModelContext | None = None,
    hooks: Iterable[PostBuildHook] = (),
    grid: Any = None,
    model_nam_file: str | None = None,
    version: str = "mf6",
    exe_name: str = "mf6",
    model_rel_path: str | Path = ".",
    list: str | None = None,
    print_input: bool | None = None,
    print_flows: bool | None = None,
    save_flows: bool | None = True,
    dependent_variable_scaling: bool | None = None,
    nc_mesh2d_filerecord: Any = None,
    nc_structured_filerecord: Any = None,
    nc_filerecord: Any = None,
    **kwargs: Any,
) -> ModelSpec:
    """Return a typed GWT (groundwater solute-transport) model spec.

    The solute-transport counterpart of :func:`gwf`: it wraps
    ``flopy.mf6.ModflowGwt`` and takes the transport packages (``mf.disv`` or a
    shared grid, plus the MST/ADV/DSP/SSM/IC/OC transport packages) through
    ``packages``. Couple it to a flow model with a GWF-GWT exchange
    (``mf.exchange(flow, transport)``) in the same :class:`SimulationSpec`.

    Parameters
    ----------
    name : str
        Model name (the FloPy model name; referenced by IMS and the GWF-GWT exchange).
    packages : Iterable[PackageSpec or PackageRef or str], optional
        Transport package specs to attach (discretization + MST/ADV/DSP/SSM/IC/OC),
        built in this order; an ``mf.ref(key)`` or bare key string names a library
        package. Names must be unique. Default ``()``.
    context : ModelContext, optional
        Geometry/dates, typically the same grid as the paired flow model. Default
        ``None`` stores an empty ``ModelContext()``.
    hooks : Iterable[PostBuildHook], optional
        Post-build callbacks, each called as ``callback(model, packages, context)``
        after every package is built; results land in ``BuiltModel.hook_results``.
        Default ``()``.
    grid : GridSpec or GridRef or grid object, optional
        A deferred ``GridSpec``, an ``mf.grid_ref(key)``, or a built grid (wrapped
        with ``GridSpec.from_object``), resolved at build time and replacing
        ``context.grid``. It must be the paired flow model's mesh: two DISV models
        joined by an exchange on different meshes are refused at build. Default
        ``None``.
    model_nam_file : str, optional
        Name-file name, forwarded to ``ModflowGwt``; omitted when ``None``.
    version : str, optional
        FloPy's model ``version`` string. Default ``"mf6"``.
    exe_name : str, optional
        Executable name recorded on the FloPy model object. Default ``"mf6"``.
    model_rel_path : str or Path, optional
        The model's folder relative to the simulation workspace. Default ``"."``.
    list : str, optional
        Listing-file name; omitted when ``None`` (MF6 names it after the name file).
    print_input : bool, optional
        MF6 ``PRINT_INPUT``: echo package input to the listing file. Omitted when
        ``None`` (off).
    print_flows : bool, optional
        MF6 ``PRINT_FLOWS``: print package flow rates to the listing file. Omitted
        when ``None`` (off).
    save_flows : bool, optional
        MF6 ``SAVE_FLOWS``: write package flow terms to the budget file. Default
        ``True`` here (unlike :func:`gwf`); ``None`` omits it.
    dependent_variable_scaling : bool, optional
        MF6 ``DEPENDENT_VARIABLE_SCALING``: scale the solution and right-hand side
        to avoid very large concentration values. Omitted when ``None``.
    nc_mesh2d_filerecord : Any, optional
        FloPy's NetCDF layered-mesh output record, passed through when not ``None``.
    nc_structured_filerecord : Any, optional
        FloPy's NetCDF structured-grid output record, passed through when not
        ``None``.
    nc_filerecord : Any, optional
        FloPy's NetCDF file record, passed through when not ``None``.
    **kwargs
        Any additional ``flopy.mf6.ModflowGwt`` constructor option, passed to it at
        build time; an unknown name fails there, not at the call.

    Returns
    -------
    ModelSpec

    Examples
    --------
    >>> transport = mf.gwt("transport", packages=[...])
    >>> mf.simulation(flow, transport, exchanges=[mf.exchange(flow, transport)])
    """

    options = _model_options(
        model_nam_file=model_nam_file,
        version=version,
        exe_name=exe_name,
        model_rel_path=model_rel_path,
        list=list,
        print_input=print_input,
        print_flows=print_flows,
        save_flows=save_flows,
        nc_mesh2d_filerecord=nc_mesh2d_filerecord,
        nc_structured_filerecord=nc_structured_filerecord,
        nc_filerecord=nc_filerecord,
        extra=kwargs,
    )
    if dependent_variable_scaling is not None:
        options["dependent_variable_scaling"] = dependent_variable_scaling

    return ModelSpec(
        name,
        "gwt",
        packages=tuple(packages),
        context=ModelContext() if context is None else context,
        hooks=tuple(hooks),
        grid=grid,
        options=options,
    )


def gwe(
    name: str,
    *,
    packages: Iterable[PackageSpec] = (),
    context: ModelContext | None = None,
    hooks: Iterable[PostBuildHook] = (),
    grid: Any = None,
    model_nam_file: str | None = None,
    version: str = "mf6",
    exe_name: str = "mf6",
    model_rel_path: str | Path = ".",
    list: str | None = None,
    print_input: bool | None = None,
    print_flows: bool | None = None,
    save_flows: bool | None = True,
    dependent_variable_scaling: bool | None = None,
    nc_mesh2d_filerecord: Any = None,
    nc_structured_filerecord: Any = None,
    nc_filerecord: Any = None,
    **kwargs: Any,
) -> ModelSpec:
    """Return a typed GWE (groundwater energy / heat-transport) model spec.

    The heat-transport counterpart of :func:`gwf`: it wraps
    ``flopy.mf6.ModflowGwe`` and takes the energy packages (EST/ADV/CND/SSM/IC/OC)
    through ``packages``. Couple it to a flow model with a GWF-GWE exchange
    (``mf.exchange(flow, energy)``) in the same :class:`SimulationSpec`. Use it to
    model aquifer thermal transport (ATES, heat pumps, thermal plumes).

    Parameters
    ----------
    name : str
        Model name (the FloPy model name; referenced by IMS and the GWF-GWE exchange).
    packages : Iterable[PackageSpec or PackageRef or str], optional
        Energy-transport package specs to attach (discretization +
        EST/ADV/CND/SSM/IC/OC), built in this order; an ``mf.ref(key)`` or bare key
        string names a library package. Names must be unique. Default ``()``.
    context : ModelContext, optional
        Geometry/dates, typically the same grid as the paired flow model. Default
        ``None`` stores an empty ``ModelContext()``.
    hooks : Iterable[PostBuildHook], optional
        Post-build callbacks, each called as ``callback(model, packages, context)``
        after every package is built; results land in ``BuiltModel.hook_results``.
        Default ``()``.
    grid : GridSpec or GridRef or grid object, optional
        A deferred ``GridSpec``, an ``mf.grid_ref(key)``, or a built grid (wrapped
        with ``GridSpec.from_object``), resolved at build time and replacing
        ``context.grid``. It must be the paired flow model's mesh: two DISV models
        joined by an exchange on different meshes are refused at build. Default
        ``None``.
    model_nam_file : str, optional
        Name-file name, forwarded to ``ModflowGwe``; omitted when ``None``.
    version : str, optional
        FloPy's model ``version`` string. Default ``"mf6"``.
    exe_name : str, optional
        Executable name recorded on the FloPy model object. Default ``"mf6"``.
    model_rel_path : str or Path, optional
        The model's folder relative to the simulation workspace. Default ``"."``.
    list : str, optional
        Listing-file name; omitted when ``None`` (MF6 names it after the name file).
    print_input : bool, optional
        MF6 ``PRINT_INPUT``: echo package input to the listing file. Omitted when
        ``None`` (off).
    print_flows : bool, optional
        MF6 ``PRINT_FLOWS``: print package flow rates to the listing file. Omitted
        when ``None`` (off).
    save_flows : bool, optional
        MF6 ``SAVE_FLOWS``: write package flow terms to the budget file. Default
        ``True`` here (unlike :func:`gwf`); ``None`` omits it.
    dependent_variable_scaling : bool, optional
        MF6 ``DEPENDENT_VARIABLE_SCALING``: scale the solution and right-hand side
        to avoid very large temperature values. Omitted when ``None``.
    nc_mesh2d_filerecord : Any, optional
        FloPy's NetCDF layered-mesh output record, passed through when not ``None``.
    nc_structured_filerecord : Any, optional
        FloPy's NetCDF structured-grid output record, passed through when not
        ``None``.
    nc_filerecord : Any, optional
        FloPy's NetCDF file record, passed through when not ``None``.
    **kwargs
        Any additional ``flopy.mf6.ModflowGwe`` constructor option, passed to it at
        build time; an unknown name fails there, not at the call.

    Returns
    -------
    ModelSpec

    Examples
    --------
    >>> energy = mf.gwe("heat", context=ctx, packages=[...])
    >>> mf.simulation(flow, energy, exchanges=[mf.exchange(flow, energy)])
    """

    options = _model_options(
        model_nam_file=model_nam_file,
        version=version,
        exe_name=exe_name,
        model_rel_path=model_rel_path,
        list=list,
        print_input=print_input,
        print_flows=print_flows,
        save_flows=save_flows,
        nc_mesh2d_filerecord=nc_mesh2d_filerecord,
        nc_structured_filerecord=nc_structured_filerecord,
        nc_filerecord=nc_filerecord,
        extra=kwargs,
    )
    if dependent_variable_scaling is not None:
        options["dependent_variable_scaling"] = dependent_variable_scaling

    return ModelSpec(
        name,
        "gwe",
        packages=tuple(packages),
        context=ModelContext() if context is None else context,
        hooks=tuple(hooks),
        grid=grid,
        options=options,
    )


def prt(
    name: str,
    *,
    packages: Iterable[PackageSpec] = (),
    context: ModelContext | None = None,
    hooks: Iterable[PostBuildHook] = (),
    grid: Any = None,
    model_nam_file: str | None = None,
    version: str = "mf6",
    exe_name: str = "mf6",
    model_rel_path: str | Path = ".",
    list: str | None = None,
    print_input: bool | None = None,
    print_flows: bool | None = None,
    save_flows: bool | None = None,
    **kwargs: Any,
) -> ModelSpec:
    """Return a typed PRT (particle tracking) model spec.

    The particle-tracking counterpart of :func:`gwf`: it wraps
    ``flopy.mf6.ModflowPrt`` and takes the PRT packages through ``packages`` --
    ``mf.disv``/``mf.dis`` (they dispatch on the PRT model kind), :func:`mip`,
    :func:`prp`, and ``mf.oc``. Couple it to a flow model with a GWF-PRT exchange
    (``mf.exchange(flow, particles)``) in the same :class:`SimulationSpec` -- flows
    pass through the exchange, so no FMI package is needed -- and give the PRT
    model an :func:`ems` solver (PRT is explicit; IMS is rejected). Use it for
    advective pathlines, capture zones, and travel-time analysis. ``context``
    carries the grid/geometry exactly as for a flow model; the higher-level
    :class:`~myflopy.modflow.mf6.prt.PRTProject` wraps the *post-hoc* pattern
    (tracking over an already-run flow model in a separate simulation, where FMI
    and grid copying ARE needed) -- prefer it for that workflow.

    Parameters
    ----------
    name : str
        Model name (also the FloPy model name), referenced by the EMS solver and the
        GWF-PRT exchange.
    packages : Iterable[PackageSpec or PackageRef or str], optional
        PRT package specs to attach (``mf.disv(...)``, ``mf.mip(...)``,
        ``mf.prp(...)``, ``mf.oc(...)``), built in this order; an ``mf.ref(key)``
        or bare key string names a library package. Names must be unique.
        Default ``()``.
    context : ModelContext, optional
        Geometry/dates :class:`ModelContext` shared with the paired flow model.
        Default ``None`` stores an empty ``ModelContext()``.
    hooks : Iterable[PostBuildHook], optional
        Post-build callbacks, each called as ``callback(model, packages, context)``
        after every package is built; results land in ``BuiltModel.hook_results``.
        Default ``()``.
    grid : GridSpec or GridRef or grid object, optional
        A deferred grid (e.g. ``mf.GridSpec.voronoi(...)``), an
        ``mf.grid_ref(key)``, or a built grid (wrapped with
        ``GridSpec.from_object``), resolved at build time and replacing
        ``context.grid``. Default ``None``.
    model_nam_file : str, optional
        Name-file name, forwarded to ``ModflowPrt``; omitted when ``None``.
    version : str, optional
        FloPy's model ``version`` string. Default ``"mf6"``.
    exe_name : str, optional
        Executable name recorded on the FloPy model object. Default ``"mf6"``.
    model_rel_path : str or Path, optional
        The model's folder relative to the simulation workspace. Default ``"."``.
    list : str, optional
        Listing-file name; omitted when ``None`` (MF6 names it after the name file).
    print_input : bool, optional
        MF6 ``PRINT_INPUT``: echo package input to the listing file. Omitted when
        ``None`` (off).
    print_flows : bool, optional
        MF6 ``PRINT_FLOWS``: print package flow rates to the listing file. Omitted
        when ``None`` (off).
    save_flows : bool, optional
        MF6 ``SAVE_FLOWS``: write package flow terms to the budget file. Omitted
        when ``None`` (off).
    **kwargs
        Additional ``flopy.mf6.ModflowPrt`` constructor options, passed to it at
        build time; an unknown name fails there, not at the call.
    """

    return ModelSpec(
        name,
        "prt",
        packages=tuple(packages),
        context=ModelContext() if context is None else context,
        hooks=tuple(hooks),
        grid=grid,
        options=_model_options(
            model_nam_file=model_nam_file,
            version=version,
            exe_name=exe_name,
            model_rel_path=model_rel_path,
            list=list,
            print_input=print_input,
            print_flows=print_flows,
            save_flows=save_flows,
            extra=kwargs,
        ),
    )


_DISCRETIZATION_TYPES = ("dis", "disv", "disu")


def _built_discretization(model):
    """The discretization already built on ``model``, or ``None``.

    Found by ``package_type`` rather than by name. ``model.get_package("dis")`` does a
    FloPy PREFIX match and hands back a ``ModflowGwfdisv`` -- measured -- and a name test
    would miss `mf.disv(stack, name="mesh")` besides. `mf.hfb` carried
    `requires=("disv",)` until ledger 378, wrong for the same reason: it names the
    myflopy SPEC name, so it false-refused a DIS or DISU model, a renamed discretization
    and a `mf.ref` whose key differs. It orders through this function now.
    """

    for package in getattr(model, "packagelist", None) or ():
        if getattr(package, "package_type", None) in _DISCRETIZATION_TYPES:
            return package
    return None


#: How many integers a cellid holds on each discretization -- what FloPy expects when it
#: WRITES a cellid, and so what a cell pair has to carry to reach MF6.
_CELLID_LENGTH = {"dis": 3, "disv": 2, "disu": 1}


def _build_hfb(package_class, model, *, _spec, **options):
    """Build an HFB after the discretization, found by TYPE (ledger 378).

    The ordering rule ``requires=("disv",)`` used to carry, without its defect: it named
    the SPEC ``disv``, so a barrier was refused on a DIS or DISU model, beside a DISV
    declared under another name, and behind an ``mf.ref`` whose resolved name differed.
    The order itself is load-bearing -- measured: FloPy builds an HFB declared ahead of
    its DIS and fails only when WRITING it ("Expected a cellid containing 0 integer(s)
    for grid type DiscretizationType.UNDEFINED"), naming neither the package nor the fix.

    Also refuses a cell pair whose cellids are the wrong SHAPE for the discretization --
    a plan pair ``(layer, cell)`` on a DIS model, say -- which FloPy likewise accepts here
    and reports only at write. Records in a form this does not read (a recarray handed
    to ``mf.hfb.flopy``) are left to FloPy.
    """

    discretization = _built_discretization(model)
    if discretization is None:
        raise ValueError(
            f"Package {_spec.name!r} is a flow barrier between cells, so the model's "
            "discretization has to be declared BEFORE it -- `mf.dis`, `mf.disv` or "
            "`mf.disu`, under any name or through `mf.ref`. Move it ahead of the barrier "
            "in `packages=[...]`."
        )
    expected = _CELLID_LENGTH[discretization.package_type]
    wrong: list[Any] = []
    period_data = options.get("stress_period_data")
    if isinstance(period_data, Mapping):
        for records in period_data.values():
            if not isinstance(records, (list, tuple)):
                continue
            for record in records:
                if not isinstance(record, (list, tuple)) or len(record) < 2:
                    continue
                for cellid in record[:2]:
                    if isinstance(cellid, (list, tuple)) and len(cellid) != expected:
                        wrong.append(tuple(cellid))
    if wrong:
        spelling = {
            "dis": "(layer, row, column)",
            "disv": "(layer, cell)",
            "disu": "(node,)",
        }[discretization.package_type]
        raise ValueError(
            f"Package {_spec.name!r}: the model is discretized by "
            f"{discretization.package_type.upper()}, whose cellid is {spelling}, but "
            f"{len(wrong)} barrier cellid(s) have another shape (first: {wrong[0]}). "
            "`mf.hfb(pairs=[(a, b)])` with plain integers writes DISV `(layer, cell)` "
            f"cellids; on this model pass full cellids, `pairs=[({spelling}, {spelling})]`."
        )
    return package_class(model, **options)


#: ONE partial, not one per declaration: ``specs._wants_spec`` and the option checks
#: cache on the builder OBJECT, and a partial hashes by identity.
_HFB_BUILDER = functools.partial(_build_hfb, flopy.mf6.ModflowGwfhfb)


def _lazy_resolution_context(model, spec_name: str, *, needs_grid: bool = True):
    """The three preconditions a lazy boundary condition needs, in order.

    Each replaces a failure that names nothing. Returns a ``ModelContext`` local to this
    package's resolution -- deliberately NOT published onto ``model.myflopy_context``,
    because publishing would change what every later package reads, including the eager
    ones, and this commit's contract is that the eager path does not move.

    ``needs_grid=False`` is for the AREAL builders, which read the idomain and a cell
    count and nothing else -- measured with a spy grid: `mf.rch(context=...)` touches
    `get_disv_gridprops` (for the pin) and `ncpl`, no geometry, no CRS, no surfaces.
    Both are on the built discretization, so demanding a grid object would make the lazy
    form strictly NARROWER than the eager one, which accepts a bare
    ``ModelContext(domain=idomain)`` and a model declared with no context at all
    (ledger 284). The ORDERING and ncpl-agreement rules stay here either way, which is
    the reason this is a parameter rather than a second function.
    """

    context = getattr(model, "myflopy_context", None)
    grid = getattr(context, "grid", None)
    if grid is None and needs_grid:
        raise ValueError(
            f"Package {spec_name!r} was declared without `context=`, so its grid has to "
            "come from the model at build time -- and this model's context carries none. "
            "Either declare the model with `mf.gwf(..., context=ModelContext(grid=...))`, "
            "give it a deferred grid with `.with_grid(mf.GridSpec.voronoi(...))`, or pass "
            "`context=` to the helper to resolve it eagerly."
        )

    discretization = _built_discretization(model)
    if discretization is None:
        raise ValueError(
            f"Package {spec_name!r} resolves its cells against the model's grid at BUILD "
            "time, so the discretization has to be declared BEFORE it. Move `mf.disv(...)` "
            "ahead of it in `packages=[...]`."
        )

    ncpl = getattr(getattr(discretization, "ncpl", None), "get_data", lambda: None)()
    if ncpl is not None and grid is not None and int(ncpl) != int(grid.ncpl):
        raise ValueError(
            f"Package {spec_name!r} resolves against this model's grid "
            f"({grid.ncpl} cells) at build time, but the discretization "
            f"already built on the model describes {int(ncpl)} cells -- they are not the "
            "same mesh. Declare the discretization from the same grid; `mf.disv(stack)` "
            "resolves against whatever grid the build supplies."
        )

    domain = getattr(context, "domain", None)
    if domain is None:
        # Read the idomain off the PACKAGE, never off `model.modelgrid.idomain`: flopy
        # SYNTHESIZES an all-ones FLOAT array when the discretization declares none,
        # where the package's own field says None -- which is what "no idomain" means in
        # MF6. Load-bearing, measured: a resolver reading only the context wrote
        # boundaries into cells the discretization makes inactive (drn 5 records against
        # a correct 4), and MF6 then refuses the run -- except under `edges_only=True`,
        # where it runs normally on the wrong boundary.
        idomain = getattr(discretization, "idomain", None)
        domain = idomain.get_data() if idomain is not None else None

    if domain is None and grid is None and ncpl is not None:
        # Nothing left to count cells with: no context, and a discretization declared
        # without `idomain=`. The areal builders fall back to `context.grid.ncpl` for
        # their cell list, so both routes are None here and they die five frames down
        # with a message naming neither the package nor the fix. An absent IDOMAIN means
        # every cell ACTIVE in MF6, so saying that outright is not a guess -- it is the
        # same model, written down (ledger 284).
        nlay = getattr(getattr(discretization, "nlay", None), "get_data", lambda: None)()
        domain = np.ones((int(nlay or 1), int(ncpl)), dtype=int)

    surfaces = getattr(context, "surfaces", None)
    if surfaces is None and getattr(grid, "gdf_topbtm", None) is None:
        # The order is PARITY, not preference: `_layer_elevations` consults the grid's
        # own table before anything else for the eager form, and lazy == eager is this
        # commit's contract. Falling back to the built discretization is what makes
        # `CellSurfaceOffset` -- the documented elevation idiom for drn/riv/ghb -- work
        # behind an eager or hand-written DISV, where it raises today.
        top = getattr(discretization, "top", None)
        botm = getattr(discretization, "botm", None)
        if top is not None and botm is not None:
            surfaces = {"top": top.array, "botm": botm.array}

    return _replace(context, domain=domain, surfaces=surfaces)


def _build_lazy_list_bc(package_class, model, *, _spec, stress_period_data, **options):
    """Resolve a :class:`~myflopy.sources.BoundarySource` against the model's grid.

    The whole of a lazy list boundary condition. It is bound to its FloPy class by
    :func:`_lazy_factory`, so the declaration-time option check sees the real package's
    signature rather than a union of all seven (ledger 279).

    Its NAME comes from ``_spec``, never from the recipe (ledger 278): a name frozen in
    the recipe survives a rename, and two renamed forks then collapse into one MF6
    package.
    """

    recipe = stress_period_data
    context = _lazy_resolution_context(model, _spec.name)
    source = GeoPackageSource(
        recipe.path,
        context,
        recipe.nper,
        layer=recipe.layer,
        name_field=recipe.name_field,
        layer_field=recipe.layer_field,
        period_field=recipe.period_field,
        layer_base=recipe.layer_base,
        period_base=recipe.period_base,
        gaps=recipe.gaps,
        snap_within=recipe.snap_within,
        min_above_bottom=recipe.min_above_bottom,
    )
    kwargs = dict(recipe.values)
    kwargs["name"] = _spec.name
    kwargs["boundnames"] = recipe.boundnames
    if _PACKAGE_EXPLORER_SPECS[recipe.package].capabilities.edges_only:
        kwargs["edges_only"] = recipe.edges_only
    eager = getattr(source, recipe.package)(**kwargs, **options)
    # `_bc_from_features` stamps GRID_PIN afresh but never clears a void, so a recipe
    # forked from an edited package would carry a pin AND its own contradiction, and
    # `with_grid` would complain about it twice (ledger 275).
    if PIN_VOIDED in eager.metadata:
        eager = _replace(
            eager,
            metadata={k: v for k, v in eager.metadata.items() if k != PIN_VOIDED},
        )
    # Record what was actually BUILT, so the run freezes records rather than this
    # recipe. Without it a reopened run goes back to the GeoPackage -- the dependence
    # ledger 251's snapshot exists to remove -- and dies if the file has moved.
    # `ModelSpec.build` collects this into `BuiltModel.resolved` and aliases the key
    # for a `mf.ref` (ledger 280).
    eager = _replace(eager, name=_spec.name)
    resolved = getattr(model, "_myflopy_resolved", None)
    if resolved is None:
        resolved = {}
        model._myflopy_resolved = resolved
    resolved[_spec.name] = eager
    return eager.build(model)


#: Keyed by the FLOPY class, not by the package name. The lazy builder is bound to its
#: FloPy class through a ``functools.partial`` so the option check sees the real
#: package's signature, and looking the builder up from that same object is what keeps
#: the partial's ``keywords`` EMPTY -- see :func:`_lazy_areal_factory` for why a
#: keyword there moved the digest (ledger 284).
_AREAL_BUILDERS = {
    flopy.mf6.ModflowGwfrch: RCHBuilder,
    flopy.mf6.ModflowGwfevt: EVTBuilder,
}


#: What each areal declaration key means when the caller said nothing. A recipe records
#: only what DIFFERS from these -- the `_explicit_options` discipline of ledger 258 --
#: so adding a defaulted argument later cannot move a stored digest.
_AREAL_DEFAULTS = {
    "cells": "top_active",
    "layer": 0,
    "name_by_cell": None,
    "boundnames": False,
}


def _explicit_areal(**declared: Any) -> dict[str, Any]:
    """The areal declaration, minus everything left at its default."""

    return {
        key: value
        for key, value in declared.items()
        if key not in _AREAL_DEFAULTS or value != _AREAL_DEFAULTS[key]
    }


def _refuse_caller_typed_cells(
    recipe: dict[str, Any], *, verb: str,
    value_keys: tuple[str, ...] = ("recharge", "rate", "depth", "surface"),
    default_cells: str = "top_active",
) -> None:
    """Refuse the parts of an areal declaration that name cells the caller chose.

    A lazy recipe is a promise that the package can be rebuilt on another mesh, and
    four of these cannot keep it:

    * an explicit ``cells=`` list -- the numbers are the caller's, and **you cannot pin
      cell numbers the caller typed**. Measured: a lazy `cells=[1, 2, 3]` built the SAME
      three records at ncpl 28 and at ncpl 92, with a `rebuilt` verdict and no warning.
    * ``name_by_cell=`` -- keyed by cell id, and on another mesh it does not error, it
      MISSES and falls back to a generated ``<name>_<layer>_<cell>``. Its tuple keys do
      not survive a JSON round trip either: ``_json_value`` stringifies them and
      ``_spec_value`` revives only int-like keys, so ``(0, 3)`` comes back ``'(0, 3)'``.
    * a ``{cellid: value}`` mapping, for the same two reasons.
    * a per-cell value -- a list, tuple, ndarray, pandas Series or any other array-like
      of one or more dimensions, bare or inside a period -- whose length is the
      selected-cell count on the mesh it was written for.

    All four are still available EAGERLY -- `context=` resolves them against a known
    grid and pins the result -- so this refuses a combination, not a capability, and
    says which one to reach for.

    ``value_keys`` is which of the recipe's entries hold VALUES rather than settings, and
    it is a parameter because the areal helpers and ``mf.uzf`` carry different ones (four
    against thirteen). The RULE is identical, which is why there is one function
    (ledger 304). Refusing at the call rather than at build is the rule
    `_lazy_bc_spec` already follows for a missing file and `nper=0`: neither needs a
    grid, so neither belongs at build time.
    """

    def caller_typed(value) -> bool:
        # A PERIOD key is an integer, a numpy one included (`isinstance(np.int64(0), int)`
        # is False, ledger 265); anything else keys by cell.
        if isinstance(value, Mapping):
            return not all(
                isinstance(key, (int, np.integer)) and not isinstance(key, bool)
                for key in value
            )
        if isinstance(value, (str, bytes)):
            return False
        # ANY array-like is per-cell, not only a `Sequence`: an ndarray and a pandas
        # Series are not `collections.abc.Sequence`s, so the lazy form accepted the same
        # per-cell values it refused as a list -- read `rebuilt` across `with_grid` while
        # riding the recipe by POSITION, and made the spec unsaveable, since `to_dict`
        # refuses an ndarray (ledger 379). A 0-d array is a scalar.
        return isinstance(value, Sequence) or np.ndim(value) > 0

    offenders: list[str] = []
    if not isinstance(recipe.get("cells", default_cells), str):
        offenders.append("cells=[...] names cells by number")
    if recipe.get("name_by_cell") is not None:
        offenders.append("name_by_cell= is keyed by cell id")
    for key in value_keys:
        value = recipe.get(key)
        if value is None:
            continue
        if caller_typed(value):
            offenders.append(f"{key}= is per-cell rather than a scalar or {{period: scalar}}")
        elif isinstance(value, Mapping) and any(caller_typed(item) for item in value.values()):
            offenders.append(f"{key}= is per-cell within a period")
    if not offenders:
        return
    raise TypeError(
        f"{verb} was declared without `context=`, so it resolves against whatever grid "
        "the build supplies -- and " + "; ".join(offenders) + ". Nothing can carry those "
        "to another mesh. Pass `context=` to resolve and pin them against a known grid, "
        "or give a scalar / {period: scalar} value here."
    )


#: The declaration keys an areal recipe carries. Everything else a lazy areal spec is
#: handed is a FloPy option -- the split has to be decidable, because ``PackageSpec.build``
#: splats ``inputs`` and ``options`` into one call and the builder sees them merged.
_AREAL_RECIPE_KEYS = frozenset(
    {"nper", "cells", "layer", "name_by_cell", "boundnames", "recharge", "rate", "depth", "surface"}
)


def _build_lazy_areal(package_class, model, *, _spec, **kwargs):
    """Resolve an areal (RCH/EVT) declaration against the model's grid at build time.

    The list form's counterpart to :func:`_build_lazy_list_bc`, and the same shape: the
    declaration rides in ``PackageSpec.inputs``, the NAME comes from ``_spec`` rather
    than from the recipe (ledger 278), and what was actually built is recorded so the
    run freezes RECORDS rather than this recipe (ledger 251/280).

    ``needs_grid=False``: these builders read the idomain and a cell count, both of
    which the built discretization carries.

    ``package_class`` is the FLOPY class, bound by :func:`_lazy_areal_factory` exactly
    as :func:`_build_lazy_list_bc` binds its own, and it chooses the areal builder.
    Binding ``RCHBuilder`` instead narrowed the accepted option set to that builder's
    eight parameters, and the lazy branch then REFUSED ``print_input``, ``print_flows``
    and ``auxiliary`` -- real MF6 options the eager branch accepts, a ``lazy != eager``
    divergence in the narrowing direction. ``_accepted_options`` recurses into a
    partial's bound ARGS and never reads ``mf6_targets`` off the partial itself, so
    setting that attribute on the partial was dead code (ledger 284).
    """

    recipe = {key: value for key, value in kwargs.items() if key in _AREAL_RECIPE_KEYS}
    options = {key: value for key, value in kwargs.items() if key not in _AREAL_RECIPE_KEYS}
    context = _lazy_resolution_context(model, _spec.name, needs_grid=False)
    eager = _AREAL_BUILDERS[package_class](
        context=context, name=_spec.name, options=options, **recipe
    ).build()
    if GRID_PIN not in (eager.metadata or {}):
        # `_ArealBuilder.build` pins through `pin_metadata`, which fingerprints
        # `context.grid` -- and a model declared with no context has none. The built
        # DISV is the same mesh by another route, and `disv_pin` measured EQUAL to
        # `grid_pin` on it, so the resolved package still records which mesh its
        # records belong to instead of silently going unpinned (ledger 282 is what
        # unpinned costs).
        disv = _built_discretization(model)
        pin = disv_pin(disv) if disv is not None else None
        if pin is not None:
            eager = eager.with_metadata(
                **{GRID_PIN: pin, GRID_PIN_OPTIONS: ["stress_period_data"]}
            )
    resolved = getattr(model, "_myflopy_resolved", None)
    if resolved is None:
        resolved = {}
        model._myflopy_resolved = resolved
    resolved[_spec.name] = eager
    return eager.build(model)


#: The declaration keys a UZF recipe carries. Everything else is a FloPy option, and the
#: split has to be decidable because ``PackageSpec.build`` splats ``inputs`` and
#: ``options`` into one call.
_UZF_RECIPE_KEYS = frozenset(
    {
        "nper", "cells", "boundnames", "mover",
        "vks", "thtr", "thts", "thti", "eps", "surfdep",
        "finf", "pet", "extdp", "extwc", "ha", "hroot", "rootact",
        "simulate_et", "linear_gwet", "square_gwet", "simulate_gwseep",
        "unsat_etwc", "unsat_etae", "ntrailwaves", "nwavesets",
    }
)

#: The UZF recipe entries that hold a VALUE rather than a setting -- what
#: :func:`_refuse_caller_typed_cells` has to inspect for per-cell identity.
_UZF_VALUE_KEYS = (
    "vks", "thtr", "thts", "thti", "eps", "surfdep",
    "finf", "pet", "extdp", "extwc", "ha", "hroot", "rootact",
)


def _build_lazy_uzf(model, *, _spec, **kwargs):
    """Resolve a UZF declaration against the model's DOMAIN at build time.

    ``needs_grid=False``, and for UZF that is not a concession -- it is the whole truth.
    Measured with an attribute spy: ``UZFBuilder.build()`` touches exactly ONE grid
    attribute, ``get_disv_gridprops``, and only from inside ``pin_metadata``; its cell set,
    its ``iuzno`` numbering and its ``ivertcon`` all come from the IDOMAIN alone
    (``uzf_cells`` walks each column, takes the first active layer as land surface and
    descends while active). It builds byte-identical ``packagedata`` against a bare
    ``ModelContext(domain=...)``. So it needs strictly LESS than ``mf.rch`` does.
    """

    recipe = {key: value for key, value in kwargs.items() if key in _UZF_RECIPE_KEYS}
    options = {key: value for key, value in kwargs.items() if key not in _UZF_RECIPE_KEYS}
    context = _lazy_resolution_context(model, _spec.name, needs_grid=False)
    eager = UZFBuilder(
        context=context, name=_spec.name, options=options, **recipe
    ).build()
    if GRID_PIN not in (eager.metadata or {}):
        # `UZFBuilder.build` pins through `pin_metadata`, which fingerprints
        # `context.grid` -- and a model declared with no context has none, in which case
        # it stamps the weaker cell-addressing marker instead (ledger 302). The built DISV
        # is the same mesh by another route and `disv_pin` measured EQUAL to `grid_pin` on
        # it, so the RESOLVED package records which mesh its records belong to rather than
        # going out under a marker that means "could not tell".
        disv = _built_discretization(model)
        pin = disv_pin(disv) if disv is not None else None
        if pin is not None:
            eager = _replace(
                eager,
                metadata={
                    key: value
                    for key, value in (eager.metadata or {}).items()
                    if key != CELLWISE_OPTIONS
                },
            ).with_metadata(**{GRID_PIN: pin, GRID_PIN_OPTIONS: ["packagedata"]})
    # Record what was actually BUILT, so the run freezes records rather than this recipe.
    # Without it a reopened run goes back to the recipe and re-resolves against whatever
    # it finds -- the dependence ledger 251's snapshot exists to remove.
    resolved = getattr(model, "_myflopy_resolved", None)
    if resolved is None:
        resolved = {}
        model._myflopy_resolved = resolved
    resolved[_spec.name] = eager
    return eager.build(model)


def _lazy_uzf_spec(*, name: str, recipe: dict[str, Any], options: dict[str, Any]):
    """A lazy UZF spec: the declaration in ``inputs``, resolved at build time.

    FLAT, like the areal recipe and for the same reason -- a recipe entry sharing a name
    with a FloPy option must COLLIDE so ``packages.edit`` is refused by name rather than
    reaching the builder twice (ledger 284). No source spec and no ``source_sha``: unlike
    every other lazy form in this library, a UZF declaration names no FILE, so there is
    nothing for the bundler, the project anchor or the content digest to get wrong.
    """

    verb = "mf.uzf"
    _refuse_caller_typed_cells(
        recipe, verb=verb, value_keys=_UZF_VALUE_KEYS, default_cells="all_active"
    )
    if int(recipe["nper"]) < 1:
        raise ValueError("nper must be at least 1.")
    if recipe.get("mover"):
        # No `uzf_connection` exists -- there is none in the tree and no
        # `MoverConnection("uzf", ...)` is constructed anywhere -- so a UZF mover endpoint
        # can only be a hand-typed `iuzno`. That number is an enumeration over the ACTIVE
        # CELLS, so it means something different on every mesh and is in range on both:
        # measured, the column at (480, 150) is iuzno 18-20 on one mesh and 33-35 on
        # another. Refused rather than moved wrongly (ledger 304); still available with
        # `context=`, which pins it.
        raise TypeError(
            "mf.uzf(mover=True) needs `context=`. A UZF mover endpoint is a hand-typed "
            "`iuzno`, which is an enumeration over the ACTIVE CELLS -- it names a "
            "different column on every mesh and is in range on all of them, so a lazy "
            "declaration cannot carry it. Pass `context=` to resolve and pin the package, "
            "or drop `mover=True`."
        )
    return PackageSpec(
        name,
        _build_lazy_uzf,
        dict(options),
        requires=(),
        inputs=dict(recipe),
    )


@functools.cache
def _lazy_areal_factory(package: str):
    """A per-package lazy areal builder, bound to its builder class.

    Memoized for the reason :func:`_lazy_factory` is: ``_wants_spec`` and
    ``_required_keyword_only`` cache on the builder OBJECT and a ``functools.partial``
    hashes by identity, so a fresh one per declaration would leak a cache entry a call.
    """

    flopy_class = getattr(flopy.mf6, _PACKAGE_EXPLORER_SPECS[package].flopy_class)
    # ONE bound arg and NO keywords. Carrying the package name as a keyword instead
    # measured a MOVED DIGEST across `to_dict`/`from_dict`: `package_digest` pickles
    # `name` and the builder reference into one payload, a declared spec has the same
    # INTERNED string in both places, and pickle emits a memo back-reference for the
    # second -- 168 bytes against 172 for the revived spec, whose two strings are
    # distinct objects. Equal values, different bytes, and `runs.drift()` would report
    # every reloaded package as EDITED (ledger 284).
    return functools.partial(_build_lazy_areal, flopy_class)


def _lazy_areal_spec(package: str, *, name: str, recipe: dict[str, Any], options: dict[str, Any]):
    """A lazy areal spec: the declaration in ``inputs``, resolved at build time."""

    verb = f"mf.{package}"
    _refuse_caller_typed_cells(recipe, verb=verb)
    if int(recipe["nper"]) < 1:
        raise ValueError("nper must be at least 1.")
    return PackageSpec(
        name,
        _lazy_areal_factory(package),
        dict(options),
        requires=(),
        # FLAT, not nested under one key: a recipe entry that shares a name with a FloPy
        # option -- `boundnames` is the live one for both packages -- must COLLIDE, so
        # that `packages.edit(key, boundnames=True)` is refused by name instead of
        # reaching the builder twice and dying inside a `functools.partial`. Nesting
        # hides every key from `_check_input_collisions` (ledger 284).
        # Hashed on the way in, for the same reason `_lazy_bc_spec` hashes its
        # GeoPackage: an `mf.evt(surface=<Surface>)` recipe reads that file at BUILD
        # time, so digesting its PATH churned on a bundle and was blind to an edit.
        inputs=stamp_source_hashes(dict(recipe)),
    )


#: Per-cell VALUE arrays a :class:`~myflopy.sources.Polygons` map may fill, per helper.
#: An allowlist rather than "anything that is a Polygons", because a helper's ``**options``
#: tail also carries things that are not per-cell arrays at all -- ``steady_state``,
#: ``transient``, filerecords -- and silently painting one of those would produce an
#: array where MF6 expects a period dict.
_POLYGON_ARRAYS: dict[str, frozenset[str]] = {
    "npf": frozenset({"k", "k22", "k33"}),
    "ic": frozenset({"strt"}),
    "sto": frozenset({"ss", "sy"}),
}

#: Array parameters that are integer FLAGS, and why a polygon map is refused on them.
#: Every rule in :func:`~myflopy.modflow.mf6.grid.geometry.polygon_cell_values` returns a
#: FLOAT, and two of the three can return one that is in no polygon -- there is no cell
#: type 0.4. ``majority`` and ``centroid`` would in fact be meaningful here, but shipping
#: that means shipping a cast this commit has not measured, so it is refused by name
#: instead of half-supported (ledger 299).
_POLYGON_FLAGS: dict[str, frozenset[str]] = {
    "npf": frozenset({"icelltype"}),
    "sto": frozenset({"iconvert"}),
}


def _split_polygon_recipes(helper: str, values: dict) -> tuple[dict, dict]:
    """Split ``values`` into plain options and :class:`Polygons` recipes.

    A recipe is MOVED out of ``options`` into ``inputs`` rather than copied: the two are
    splatted into one builder call, so a key in both is a ``TypeError`` five frames down
    naming FloPy rather than the declaration.
    """

    fillable = _POLYGON_ARRAYS.get(helper, frozenset())
    flags = _POLYGON_FLAGS.get(helper, frozenset())
    options: dict = {}
    recipes: dict = {}
    for key, value in values.items():
        if not isinstance(value, Polygons):
            options[key] = value
            continue
        if key in flags:
            raise TypeError(
                f"mf.{helper}({key}=mf.Polygons(...)) is refused: {key!r} is an integer "
                "flag, not a measured value, and every polygon rule returns a float -- "
                "there is no cell type 0.4. Paint the flag with a plain array or a "
                "scalar; a polygon map is for a value that varies by mapped unit."
            )
        if key not in fillable:
            raise TypeError(
                f"mf.{helper}({key}=mf.Polygons(...)) is refused: {key!r} is not a "
                "per-cell array this helper fills from a polygon map. It fills "
                + ", ".join(sorted(fillable))
                + "."
            )
        recipes[key] = value
    return options, recipes


def _refuse_unpaintable_array(helper: str, parameter: str, value) -> None:
    """Refuse, at the DECLARATION, a value an array parameter accepts but cannot use.

    ``mf.npf(k=<Surface>)`` was accepted here and died at BUILD inside FloPy with an
    ``MFDataException`` that names neither ``Surface`` nor what to do -- measured, and
    the same for a ``LayerStack``. A slot that says yes to something it cannot use is
    the defect, not the message; this closes the yes.

    Deliberately NARROW. It names the two types measured to fail rather than allowing a
    list of accepted ones, because the array slots genuinely take a scalar, a sequence,
    an ndarray, an ``OPEN/CLOSE`` path string and FloPy's own layered dicts, and a rule
    written the other way round would refuse something that works today.
    """

    token = getattr(type(value), "_SPEC_TOKEN", None)
    if token not in ("$surface", "$layer_stack"):
        return
    noun = "Surface" if token == "$surface" else "LayerStack"
    raise TypeError(
        f"mf.{helper}({parameter}=<{noun}>) is not supported: a {noun} describes an "
        "ELEVATION field and MF6 is handed it unresolved, which fails inside FloPy with "
        "an error naming neither. Pass a scalar, a per-cell array, or "
        f"`mf.Polygons(path, column=...)` to paint {parameter} from a polygon map."
    )


def _array_package_spec(helper: str, name: str, builder, values: dict) -> PackageSpec:
    """A :class:`PackageSpec` for an array package, lazy only if a recipe is present.

    The eager path is byte-identical to what it was -- same builder object, same
    ``options``, no ``inputs`` -- so declaring ``mf.npf(k=10.0)`` cannot move a stored
    digest. Only a declaration that actually carries a ``Polygons`` becomes lazy.
    """

    for parameter, value in values.items():
        _refuse_unpaintable_array(helper, parameter, value)
    options, recipes = _split_polygon_recipes(helper, values)
    if not recipes:
        return PackageSpec(name, builder, options)
    return PackageSpec(
        name,
        _polygon_array_factory(builder),
        options,
        requires=(),
        inputs={key: _stamped_polygons(recipe) for key, recipe in recipes.items()},
    )


def _stamped_polygons(recipe: Polygons) -> Polygons:
    """Resolve the recipe's path and stamp its content hash, both at DECLARATION.

    Resolved HERE because the declaration would otherwise check the path against the
    cwd while the build anchors it to the project root -- two bases for one string, and
    with a same-named file under the root a model BUILDS CLEAN from the wrong one
    (ledger 293). Hashed here for the same reason ``BoundarySource`` is: the digest must
    read the file's CONTENT, not its address.
    """

    path = Path(recipe.path)
    if not path.exists():
        raise FileNotFoundError(f"Polygon source not found: {path}")
    path = path.resolve()
    return _replace(recipe, path=path, source_sha=source_fingerprint(path))


def _build_polygon_arrays(builder, model, *, _spec, **values):
    """Resolve every :class:`Polygons` in ``values`` against the model's grid, then build.

    Bound to its real builder by :func:`_polygon_array_factory` so the declaration-time
    option check still sees that package's own signature -- ``_accepted_options``
    unwraps a ``functools.partial`` and inspects what is inside it, so binding a
    generic wrapper WITHOUT the target would silently turn the typo check off
    (ledger 284's trap, one helper over).
    """

    context = _lazy_resolution_context(model, _spec.name)
    discretization = _built_discretization(model)
    nlay = getattr(getattr(discretization, "nlay", None), "get_data", lambda: None)()
    resolved = {
        key: (
            array_from_polygons(
                value, context.grid, nlay=int(nlay or 1), parameter=f"{_spec.name}.{key}"
            )
            if isinstance(value, Polygons)
            else value
        )
        for key, value in values.items()
    }
    package = builder(model, **resolved)
    # Record an EAGER spec holding the resolved arrays, so the run freezes numbers
    # rather than this recipe. It must go in `_myflopy_resolved` -- the attribute
    # `ModelSpec.build` harvests into `BuiltModel.resolved` and `Run._freeze` passes to
    # `freeze_simulation(resolved=)`. A separate attribute of its own looked equivalent
    # and was dead: the snapshot then held the RECIPE, so a reopened run went back to
    # the polygon file and died once it was gone -- ledger 270's defect exactly, one
    # package over, which is how it was found here.
    store = getattr(model, "_myflopy_resolved", None)
    if store is None:
        store = {}
        model._myflopy_resolved = store
    store[_spec.name] = _replace(_spec, builder=builder, options=dict(resolved), inputs={})
    return package


@functools.cache
def _polygon_array_factory(builder):
    """The array-resolving builder, bound to the real one.

    Memoized for the reason :func:`_lazy_factory` is: ``specs._required_keyword_only``
    and ``specs._wants_spec`` cache on the builder OBJECT, and a ``functools.partial``
    hashes by identity, so a fresh partial per declaration would leak a cache entry per
    call.
    """

    return functools.partial(_build_polygon_arrays, builder)


@functools.cache
def _lazy_factory(package: str):
    """A per-package lazy builder, bound to its FloPy class.

    Memoized because ``specs._required_keyword_only`` and ``specs._wants_spec`` are both
    cached on the builder OBJECT and a ``functools.partial`` hashes by identity -- a
    fresh partial per declaration would add a permanent cache entry per call.
    """

    flopy_class = getattr(flopy.mf6, _PACKAGE_EXPLORER_SPECS[package].flopy_class)
    return functools.partial(_build_lazy_list_bc, flopy_class)


def _lazy_bc_spec(
    package: str,
    path: PathLike,
    *,
    name: str,
    nper: int,
    layer: str | None,
    name_field: str | None,
    layer_field: str | None,
    period_field: str | None,
    layer_base: int,
    period_base: int,
    gaps: str,
    snap_within: float | None,
    # DEFAULTED, unlike `gaps`/`snap_within`: only chd/ghb/drn/riv have an elevation
    # field, so only those four pass it. Required here it broke the other three on the
    # lazy path outright -- `_lazy_bc_spec() missing 1 required keyword-only argument`.
    min_above_bottom: float | None = None,
    values: dict[str, Any],
    # `wel`/`rch`/`evt` have no `edges_only` parameter at all -- the package registry
    # says so (`capabilities.edges_only`) and their signatures agree -- so their
    # wrappers pass nothing and the recipe records the default.
    edges_only: bool = False,
    boundnames: bool,
    options: dict[str, Any],
) -> PackageSpec:
    """A lazy list-BC spec: the recipe in ``inputs``, resolved at build time.

    One helper for all seven so the registry key and the bound FloPy class can never
    disagree -- they are derived from the same ``package`` string.

    The recipe rides under ``stress_period_data``, the option it supplies. That is
    ledger 272's case (b) -- a source shadowing its own option -- and the only shape
    ``_check_input_collisions`` permits there.

    **The declaration-time checks are repeated here on purpose.** They live in
    ``GeoPackageSource.__post_init__``, which the lazy branch never reaches, so without
    them `mf.chd.gpkg("nope.gpkg", nper=0)` is ACCEPTED and raises at BUILD -- measured.
    Neither needs a grid, so neither belongs at build time; the eager form refuses both
    at the call and the lazy form must match it.
    """

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"GeoPackage not found: {path}")
    # RESOLVED, so the recipe names the file this call just validated. A relative path
    # was stored as typed and checked against the CURRENT WORKING DIRECTORY here, while
    # `_anchored_packages` re-based it against the PROJECT ROOT at build -- two bases
    # for one string. Measured: a declaration that passed then failed at build, and
    # where a same-named file happened to sit under the project root, the model BUILT
    # CLEAN from the wrong one (heads 999.0 against the declared 9.0, no warning).
    # Resolving here leaves the anchor only the paths the BUNDLER wrote, which are
    # project-relative by construction. Free of digest consequences because the content
    # hash below is what the digest reads.
    path = path.resolve()
    if nper < 1:
        raise ValueError("nper must be at least 1.")
    # Neither needs a grid, so neither belongs at build time -- the same rule the two
    # checks above already follow. `snap_within` is NOT resolved here: a lazy
    # declaration has no mesh yet, and a number resolved now would carry the OLD grid's
    # cell width across `with_grid`, which is the stale-cellid defect in another costume
    # (ledger 288).
    refuse_spread_on_elevation(package, values)
    check_gaps(gaps, verb=f"mf.{package}.gpkg")
    if snap_within is not None and not (float(snap_within) > 0):
        raise ValueError(
            f"mf.{package}.gpkg: snap_within is a distance in the grid's own coordinate "
            f"units and must be positive, not {snap_within!r}."
        )

    # ONE read, used twice on purpose. It is the DIGEST's stand-in for the path
    # (`BoundarySource.digest_payload`) and the metadata `runs.drift()`'s source column
    # compares against -- and those are different questions: the recipe's copy says what
    # the file said when this recipe was written, while the column asks what it says
    # now. `None` for an unreadable file, exactly as the eager stamp does; neither is a
    # reason to refuse the declaration.
    sha = source_fingerprint(path)
    recipe = BoundarySource(
        path=path,
        package=package,
        nper=nper,
        source_sha=sha,
        values=dict(values),
        layer=layer,
        name_field=name_field,
        layer_field=layer_field,
        period_field=period_field,
        layer_base=layer_base,
        period_base=period_base,
        gaps=gaps,
        snap_within=snap_within,
        min_above_bottom=min_above_bottom,
        edges_only=edges_only,
        boundnames=boundnames,
    )
    return PackageSpec(
        name,
        _lazy_factory(package),
        dict(options),
        requires=(),
        inputs={"stress_period_data": recipe},
        # The same two keys the EAGER form stamps in `GeoPackageSource._metadata`, and
        # they are what `RunLibrary._library_source_memo` keys off. Without them
        # `runs.drift()`'s source column was INVERTED (ledger 293): the eager package,
        # whose records are frozen and whose file is never reopened, reported
        # `source=missing` when its GeoPackage moved, while the lazy one -- which
        # CANNOT BUILD without that file -- reported `—`, meaning "names no source".
        # Excluded from `package_digest` like the rest of `metadata`, so stamping them
        # moves nothing; the recipe's own `source_sha` is what the digest reads.
        metadata={"source_type": "geopackage", "source_path": str(path), SOURCE_SHA: sha},
    )


def _source(
    path: PathLike,
    *,
    context: ModelContext,
    nper: int,
    layer: str | None = None,
    name_field: str | None = "name",
    layer_field: str | None = "layer",
    period_field: str | None = None,
    layer_base: int = 1,
    period_base: int = 0,
    gaps: str = "warn",
    snap_within: float | None = None,
    min_above_bottom: float | None = None,
) -> GeoPackageSource:
    """Return a configured GeoPackage source for package-first helpers."""

    return GeoPackageSource(
        path,
        context,
        nper,
        layer=layer,
        name_field=name_field,
        layer_field=layer_field,
        period_field=period_field,
        layer_base=layer_base,
        period_base=period_base,
        gaps=gaps,
        snap_within=snap_within,
        min_above_bottom=min_above_bottom,
    )


def tdis(
    *,
    nper: int = 1,
    perioddata: Any = ((1.0, 1, 1.0),),
    time_units: str | None = None,
    start_date_time: str | None = None,
    ats: bool | Iterable[int] | Mapping[int, Mapping[str, Any]] | None = None,
    ats_dt0: float | None = None,
    ats_dtmin: float | None = None,
    ats_dtmax: float | None = None,
    ats_dtadj: float = 2.0,
    ats_dtfailadj: float = 5.0,
    ats_perioddata: Any = None,
    filename: str | None = None,
    pname: str | None = None,
    name: str = "tdis",
    **kwargs: Any,
) -> PackageSpec:
    """Return the simulation-wide TDIS (time-discretization) package spec.

    MODFLOW 6 has one TDIS package per simulation. GWF, GWT, GWE, and PRT
    models inside the same ``SimulationSpec`` share this timing. To use
    different timing, build a separate simulation/run.

    Parameters
    ----------
    nper : int, default 1
        Number of stress periods (must match ``len(perioddata)``).
    perioddata : sequence of tuple, default ((1.0, 1, 1.0),)
        One ``(perlen, nstp, tsmult)`` tuple per period -- period length, number of
        time steps, and time-step multiplier.
    time_units : str, optional
        Time unit label, e.g. ``"days"`` / ``"seconds"``.
    start_date_time : str, optional
        ISO start datetime for the simulation.
    ats : bool, iterable of int, or mapping, optional
        Turn on adaptive time stepping. MODFLOW 6 then chooses its own step
        within each ATS period, shrinking it where the solver struggles and
        growing it back where it does not -- and, crucially, **retrying a failed
        step at a smaller step instead of ending the run**.

        ``True`` applies it to every period; an iterable selects periods by
        **zero-based** index (``ats=[7, 9]``); a mapping additionally overrides
        the settings for a period (``ats={7: {"dtmin": 1e-4}}``, keys
        ``dt0``/``dtmin``/``dtmax``/``dtadj``/``dtfailadj``).
    ats_dt0 : float, optional
        Initial step for ATS periods. ``None`` (default) uses the first sub-step
        the equivalent fixed-step period would have taken.
    ats_dtmin : float, optional
        Smallest step allowed. ``None`` -> ``perlen * 1e-5``.
    ats_dtmax : float, optional
        Largest step allowed. ``None`` -> ``perlen``, i.e. a whole period in one
        step when nothing is straining the solver.
    ats_dtadj : float, default 2.0
        Factor the step grows or shrinks by according to solver effort. Must be
        ``0``, ``1``, or greater than 1; ``0``/``1`` disable growth.
    ats_dtfailadj : float, default 5.0
        Divisor used to retry a step that FAILED to converge. Must be ``0`` or
        greater than 1. **``0`` means a failed step ends the run** -- which is
        MODFLOW's behaviour without ATS, and rarely what you want if you enabled
        ATS in the first place.
    ats_perioddata : optional
        Raw FloPy ATS records, for full control -- the escape hatch, like
        ``.flopy(...)`` on the boundary helpers. Mutually exclusive with ``ats``.
    filename, pname : str, optional
        FloPy file name / package name overrides.
    name : str, default "tdis"
        Spec name.
    **kwargs
        Extra ``flopy.mf6.ModflowTdis`` options.

    Returns
    -------
    PackageSpec

    Examples
    --------
    >>> mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])                    # one steady period
    >>> mf.tdis(nper=12, perioddata=[(30.0, 3, 1.1)] * 12, time_units="days")
    >>> mf.tdis(nper=72, perioddata=monthly, time_units="days",
    ...         ats=True, ats_dtfailadj=5.0)      # retry a failed step at dt/5
    """

    from myflopy.modflow.mf6.simulation.discretization import (
        _build_ats_records,
        _resolve_ats_periods,
    )

    if ats is not None and ats is not False and ats_perioddata is not None:
        raise ValueError(
            "tdis: pass either ats= (built for you) or ats_perioddata= (raw FloPy "
            "records), not both."
        )
    if ats is not None and ats is not False:
        periods = _resolve_ats_periods(ats, nper)
        if periods:
            ats_perioddata = _build_ats_records(
                list(perioddata), periods,
                ats_dt0, ats_dtmin, ats_dtmax, ats_dtadj, ats_dtfailadj,
            )

    values = {"nper": nper, "perioddata": perioddata, **kwargs}
    for key, value in {
        "time_units": time_units,
        "start_date_time": start_date_time,
        "ats_perioddata": ats_perioddata,
        "filename": filename,
        "pname": pname,
    }.items():
        if value is not None:
            values[key] = value
    return PackageSpec(name, flopy.mf6.ModflowTdis, values)


# Captured under a private name so `simulation()` can reach the default timing
# builder from inside a module that also exports `tdis` publicly.
_steady_tdis = tdis


def ims(
    *,
    models: Iterable[str] | None = None,
    name: str = "ims",
    print_option: str | None = None,
    complexity: str | None = None,
    csv_output_filerecord: Any = None,
    csv_outer_output_filerecord: Any = None,
    csv_inner_output_filerecord: Any = None,
    no_ptcrecord: Any = None,
    ats_outer_maximum_fraction: float | None = None,
    outer_hclose: float | None = None,
    outer_dvclose: float | None = None,
    outer_rclosebnd: float | None = None,
    outer_maximum: int | None = None,
    under_relaxation: str | None = None,
    under_relaxation_gamma: float | None = None,
    under_relaxation_theta: float | None = None,
    under_relaxation_kappa: float | None = None,
    under_relaxation_momentum: float | None = None,
    backtracking_number: int | None = None,
    backtracking_tolerance: float | None = None,
    backtracking_reduction_factor: float | None = None,
    backtracking_residual_limit: float | None = None,
    inner_maximum: int | None = None,
    inner_hclose: float | None = None,
    inner_dvclose: float | None = None,
    rcloserecord: Any = None,
    linear_acceleration: str | None = None,
    relaxation_factor: float | None = None,
    preconditioner_levels: int | None = None,
    preconditioner_drop_tolerance: float | None = None,
    number_orthogonalizations: int | None = None,
    scaling_method: str | None = None,
    reordering_method: str | None = None,
    filename: str | None = None,
    pname: str | None = None,
    **kwargs: Any,
) -> PackageSpec:
    """Return an IMS package spec and register it to one or more models.

    ``models`` controls which GWF/GWT/GWE models use this solver. PRT models are
    *explicit* in MF6 and cannot sit in an IMS solution group ("Explicit models
    require EMS6") -- give them :func:`ems` instead. A single IMS can be shared:

    ``mf.ims(models=("flow", "transport"))``

    or separate IMS packages can be declared for different models:

    ``mf.ims(name="flow_solver", models=("flow",))``

    ``mf.ims(name="transport_solver", models=("transport",))``

    ``complexity`` (``"SIMPLE"``/``"MODERATE"``/``"COMPLEX"``) presets the solver
    tolerances; bump it (and ``outer_maximum``/``inner_maximum``) for stiff models
    with wetting/drying or many advanced packages (SFR/LAK/UZF). Pair with a GWF
    ``newtonoptions=...`` for robust water-table convergence.

    Parameters
    ----------
    models : Iterable[str], optional
        Names of the models this solver solves; the IMS is registered to each at
        build time. Default ``None`` leaves it unset so :func:`simulation` can
        assign every model no other solver claims; an IMS whose ``models`` is
        still unset at build raises ``ValueError`` there.
    name : str, optional
        Spec name, unique among the simulation's packages, and the default
        ``pname``. Default ``"ims"`` -- two solvers in one simulation need
        distinct names.
    print_option : str, optional
        MF6 ``PRINT_OPTION``: ``"NONE"`` / ``"SUMMARY"`` / ``"ALL"`` solver
        convergence printing. MF6's default is ``"NONE"``.
    complexity : str, optional
        Tolerance preset: ``"SIMPLE"`` / ``"MODERATE"`` / ``"COMPLEX"``; the other
        solver arguments override individual values. MF6 uses ``"SIMPLE"`` when
        it is omitted.
    csv_output_filerecord : Any, optional
        FloPy record naming a CSV file for solver convergence information.
    csv_outer_output_filerecord : Any, optional
        FloPy record naming a CSV file for outer-iteration convergence.
    csv_inner_output_filerecord : Any, optional
        FloPy record naming a CSV file for inner-iteration convergence.
    no_ptcrecord : Any, optional
        FloPy record for MF6 ``NO_PTC`` (``"FIRST"``/``"ALL"``), disabling
        pseudo-transient continuation in steady-state Newton periods.
    ats_outer_maximum_fraction : float, optional
        Fraction (0-0.5) of ``outer_maximum`` that adaptive time stepping uses to
        grow or shrink the step; MF6's default is one third.
    outer_hclose : float, optional
        Outer (nonlinear) head-change criterion, model length units; MF6's
        deprecated spelling of ``outer_dvclose``.
    outer_dvclose : float, optional
        Outer (nonlinear) dependent-variable change criterion, in the dependent
        variable's units (length for head).
    outer_rclosebnd : float, optional
        MF6's deprecated package-residual tolerance; no effect since MF6 6.1.1.
    outer_maximum : int, optional
        Maximum outer (nonlinear) iterations per time step -- raise for stiff
        models.
    under_relaxation : str, optional
        Nonlinear under-relaxation scheme: ``"NONE"`` / ``"SIMPLE"`` /
        ``"COOLEY"`` / ``"DBD"``. MF6 applies none when omitted.
    under_relaxation_gamma : float, optional
        Relaxation factor for ``"SIMPLE"``, or the memory term for
        ``"COOLEY"``/``"DBD"``.
    under_relaxation_theta : float, optional
        Delta-bar-delta learning-rate reduction factor (0-1); ``"DBD"`` only.
    under_relaxation_kappa : float, optional
        Delta-bar-delta learning-rate increment (0-1); ``"DBD"`` only.
    under_relaxation_momentum : float, optional
        Delta-bar-delta momentum fraction (0-1); ``"DBD"`` only.
    backtracking_number : int, optional
        Maximum backtracking iterations for residual reduction; ``0`` disables
        backtracking.
    backtracking_tolerance : float, optional
        Allowed residual increase factor before backtracking (>= 1).
    backtracking_reduction_factor : float, optional
        Step-size reduction per backtracking iteration (0-1).
    backtracking_residual_limit : float, optional
        Residual below which no further backtracking is done.
    inner_maximum : int, optional
        Maximum inner (linear) iterations per outer iteration.
    inner_hclose : float, optional
        Inner (linear) head-change criterion, model length units; MF6's deprecated
        spelling of ``inner_dvclose``.
    inner_dvclose : float, optional
        Inner (linear) dependent-variable change criterion.
    rcloserecord : Any, optional
        FloPy record ``(inner_rclose, rclose_option)``: the linear flow-residual
        tolerance (volume per time, model units) and optionally ``"STRICT"`` /
        ``"L2NORM_RCLOSE"`` / ``"RELATIVE_RCLOSE"``.
    linear_acceleration : str, optional
        Linear solver, ``"CG"`` or ``"BICGSTAB"``. An asymmetric matrix
        (transport, Newton, an implicit lake) needs ``"BICGSTAB"``; MF6 stops
        on CG.
    relaxation_factor : float, optional
        Relaxation factor (0-1) for the MILU(0)/MILUT preconditioners; MF6's
        default is 0.
    preconditioner_levels : int, optional
        ILUT/MILUT level of fill; above 0 selects ILUT. MF6's default is 0.
    preconditioner_drop_tolerance : float, optional
        ILUT/MILUT drop tolerance. MF6's default is 0.
    number_orthogonalizations : int, optional
        Interval at which the residual is recomputed explicitly. MF6's default
        is 0 (never).
    scaling_method : str, optional
        Matrix scaling: ``"NONE"`` / ``"DIAGONAL"`` / ``"L2NORM"``.
    reordering_method : str, optional
        Matrix reordering: ``"NONE"`` / ``"RCM"`` / ``"MD"``.
    filename : str, optional
        File the IMS is written to. Default ``None`` writes ``"<pname>.ims"``,
        which keeps two solvers in one simulation in separate files.
    pname : str, optional
        FloPy package name. Default ``None`` uses ``name``.
    **kwargs
        Any other ``flopy.mf6.ModflowIms`` option, passed to it at build time. A
        name ``ModflowIms`` does not accept raises ``TypeError`` at the call.

    Returns
    -------
    PackageSpec

    Notes
    -----
    Every solver argument mirrors ``flopy.mf6.ModflowIms`` one-to-one and is
    only forwarded when not ``None``, so MF6's own default (or the
    ``complexity`` preset) applies to anything left out.

    Examples
    --------
    >>> mf.ims(models=("flow",), complexity="MODERATE",
    ...        outer_maximum=100, inner_maximum=100, linear_acceleration="BICGSTAB")
    """

    # `models=None` is left OUT rather than stored as None: `simulation()` fills
    # it with whatever the other solvers did not claim, and an absent key is what
    # `_claim_models` and `_ordered_solvers` both test for.
    values: dict[str, Any] = {"pname": name if pname is None else pname, **kwargs}
    if models is not None:
        values["models"] = tuple(models)
    for key, value in {
        "print_option": print_option,
        "complexity": complexity,
        "csv_output_filerecord": csv_output_filerecord,
        "csv_outer_output_filerecord": csv_outer_output_filerecord,
        "csv_inner_output_filerecord": csv_inner_output_filerecord,
        "no_ptcrecord": no_ptcrecord,
        "ats_outer_maximum_fraction": ats_outer_maximum_fraction,
        "outer_hclose": outer_hclose,
        "outer_dvclose": outer_dvclose,
        "outer_rclosebnd": outer_rclosebnd,
        "outer_maximum": outer_maximum,
        "under_relaxation": under_relaxation,
        "under_relaxation_gamma": under_relaxation_gamma,
        "under_relaxation_theta": under_relaxation_theta,
        "under_relaxation_kappa": under_relaxation_kappa,
        "under_relaxation_momentum": under_relaxation_momentum,
        "backtracking_number": backtracking_number,
        "backtracking_tolerance": backtracking_tolerance,
        "backtracking_reduction_factor": backtracking_reduction_factor,
        "backtracking_residual_limit": backtracking_residual_limit,
        "inner_maximum": inner_maximum,
        "inner_hclose": inner_hclose,
        "inner_dvclose": inner_dvclose,
        "rcloserecord": rcloserecord,
        "linear_acceleration": linear_acceleration,
        "relaxation_factor": relaxation_factor,
        "preconditioner_levels": preconditioner_levels,
        "preconditioner_drop_tolerance": preconditioner_drop_tolerance,
        "number_orthogonalizations": number_orthogonalizations,
        "scaling_method": scaling_method,
        "reordering_method": reordering_method,
        "filename": filename,
    }.items():
        if value is not None:
            values[key] = value
    return PackageSpec(name, build_ims, values)


def ems(
    *,
    models: Iterable[str] | None = None,
    name: str = "ems",
    filename: str | None = None,
    pname: str | None = None,
    **kwargs: Any,
) -> PackageSpec:
    """Return an EMS (explicit solver) package spec registered to PRT models.

    The explicit counterpart of :func:`ims`: PRT models are *explicit* in MF6 --
    putting one in an IMS solution group aborts the run with "Explicit models
    require EMS6", and an EMS that is merely constructed (not registered) is
    silently dropped from ``mfsim.nam``. This spec registers properly via the
    :func:`~myflopy.builders.build_ems` engine::

        mf.ems(models=("particles",))

    EMS has no solver settings (it just marks the solution group explicit);
    ``filename``/``pname`` are the only FloPy options.

    Parameters
    ----------
    models : Iterable[str], optional
        Names of the (PRT) models this EMS solves; it is registered to each at
        build time. Default ``None`` leaves it unset for :func:`simulation` to
        fill; an EMS whose ``models`` is still unset at build raises
        ``ValueError``.
    name : str, optional
        Spec name, unique among the simulation's packages, and the default
        ``pname``. Default ``"ems"``.
    filename : str, optional
        File the EMS is written to. Default ``None`` writes ``"<pname>.ems"``.
    pname : str, optional
        FloPy package name. Default ``None`` uses ``name``.
    **kwargs
        Any other ``flopy.mf6.ModflowEms`` option, passed to it at build time. A
        name ``ModflowEms`` does not accept raises ``TypeError`` at the call.

    Returns
    -------
    PackageSpec
    """

    values: dict[str, Any] = {"pname": name if pname is None else pname, **kwargs}
    if models is not None:
        values["models"] = tuple(models)
    if filename is not None:
        values["filename"] = filename
    return PackageSpec(name, build_ems, values)


def _asymmetric_matrix(model: ModelSpec) -> bool:
    """Whether this model's coefficient matrix is asymmetric, so CG cannot solve it.

    Not a tuning preference -- MF6 refuses outright: *"MODEL x PRODUCES AN
    ASYMMETRIC COEFFICIENT MATRIX, BUT THE CONJUGATE GRADIENT METHOD WAS
    SELECTED. USE BICGSTAB INSTEAD."*

    Keyed on the MATRIX, not the model kind, and that distinction is the whole
    point. Transport and energy are asymmetric because they advect. A GWF model
    is asymmetric when it is formulated with NEWTON -- so ``mf.simulation(flow)``
    on an unconfined model failed exactly the way ``mf.simulation(flow,
    transport)`` did, and a rule written against ``model_type`` alone would have
    fixed one and left the other. Both measured: under ``complexity="SIMPLE"``,
    GWT and GWF-with-Newton both fail and both pass with ``BICGSTAB``.

    A GWF model holding a lake declared ``mf.lak(..., implicit=True)`` is asymmetric
    too: MF6's LAK ``IMPLICIT`` adds each stage to the matrix, and MF6 stops on CG with
    it (measured; ledger 326). Only a lake declared inline is seen here -- an
    ``mf.ref`` to a library lake is not resolved yet. The solver this module generates
    is therefore marked ``GENERATED_SOLVER``, and ``SimulationSpec.build_flopy`` gives it
    BICGSTAB once the refs resolve to an implicit lake.
    """

    if model.model_type in ("gwt", "gwe"):
        return True
    if model.model_type != "gwf":
        return False
    # MF6's LAK IMPLICIT makes the matrix asymmetric too, and MF6 stops on CG with it
    # (measured; ledger 326) -- so a model holding an implicit lake gets BICGSTAB.
    implicit_lake = any(
        isinstance(package, PackageSpec) and (package.options or {}).get("implicit") is True
        for package in model.packages
    )
    return implicit_lake or bool(model.options.get("newtonoptions"))


def _partition_simulation_packages(
    packages: Iterable[Any],
) -> tuple[list[PackageSpec], list[PackageSpec], list[Any]]:
    """Split a ``packages=`` list into ``(timing, solvers, everything else)``.

    Identified by BUILDER, the same signal :func:`~myflopy.specs._package_scope`
    uses, so it does not depend on what a package happens to be named.

    A ``PackageRef`` lands in "everything else" and does not suppress a default:
    an unresolved reference carries no builder, so a ref standing for a TDIS is
    indistinguishable from one standing for anything else, and guessing would
    either drop the caller's timing or write two.
    """

    timing: list[PackageSpec] = []
    solvers: list[PackageSpec] = []
    others: list[Any] = []
    for entry in packages:
        if not isinstance(entry, PackageSpec):
            others.append(entry)
        elif entry.builder is flopy.mf6.ModflowTdis:
            timing.append(entry)
        elif entry.builder in (build_ims, build_ems):
            solvers.append(entry)
        else:
            others.append(entry)
    return timing, solvers, others


def _claim_models(
    solvers: list[PackageSpec], models: tuple[ModelSpec, ...]
) -> list[PackageSpec]:
    """Give every model a solver, without overriding one the caller supplied.

    A solver declaring ``models=`` claims exactly those. A solver declaring none
    claims everything left over -- which is what lets a one-model simulation say
    ``mf.ims(complexity="COMPLEX")`` and name nothing. Models still unclaimed get
    a synthesized solver each.

    PRT is never folded into an unnamed IMS: MF6 rejects an explicit model under
    IMS6, so a PRT model always gets its own EMS even when an open IMS is present.
    """

    claimed = {
        name for spec in solvers for name in (spec.options.get("models") or ())
    }
    open_solvers = [spec for spec in solvers if not spec.options.get("models")]
    if len(open_solvers) > 1:
        names = ", ".join(spec.name for spec in open_solvers)
        raise ValueError(
            f"simulation() was given {len(open_solvers)} solvers with no models=: "
            f"{names}. At most one can claim the models the others left -- give "
            "the rest an explicit models=(...)."
        )

    unclaimed = [model for model in models if model.name not in claimed]
    resolved: list[PackageSpec] = []
    for spec in solvers:
        if spec.options.get("models"):
            resolved.append(spec)
            continue
        # The open solver takes every unclaimed model except PRT, which needs an
        # EMS of its own whatever else is going on.
        taken = tuple(m.name for m in unclaimed if m.model_type != "prt")
        if not taken:
            raise ValueError(
                f"Solver '{spec.name}' declares no models= and there are none "
                "left for it to claim."
            )
        resolved.append(spec.with_options(models=taken))
        unclaimed = [m for m in unclaimed if m.name not in taken]

    for model in unclaimed:
        if model.model_type == "prt":
            resolved.append(ems(models=[model.name], name=f"{model.name}_ems"))
            continue
        # NOT `options` -- `simulation()` has its own `**options`.
        solver_options: dict[str, Any] = {"complexity": "SIMPLE"}
        if _asymmetric_matrix(model):
            solver_options["linear_acceleration"] = "BICGSTAB"
        # Marked as myflopy's own choice, so the build may finish it once an `mf.ref`
        # resolves to an implicit lake this function cannot see (ledger 326).
        resolved.append(
            ims(models=[model.name], name=f"{model.name}_ims", **solver_options)
            .with_metadata(**{GENERATED_SOLVER: True})
        )
    return resolved


def _ordered_solvers(
    solvers: list[PackageSpec], models: tuple[ModelSpec, ...]
) -> list[PackageSpec]:
    """Order solutions the way MF6 requires them in ``mfsim.nam``.

    A GWF solution must be listed before the solution of any model it feeds:
    *"The IMS specified for GWF must be listed in mfsim.nam before the IMS for
    GWT"*, and the GWF-PRT pairing says the same about the EMS. Nothing enforced
    it, so a hand-written ``solver=[transport_ims, flow_ims]`` wrote a name file
    MF6 rejects, for a reason that reads like a physics problem.

    Model DECLARATION order is the order the caller already expressed, so sort by
    it rather than trusting the order a solver list happened to use. The sort is
    stable, and a solver naming no known model keeps its place at the end.
    """

    position = {model.name: index for index, model in enumerate(models)}
    unknown = len(position)

    def rank(spec: PackageSpec) -> int:
        named = spec.options.get("models") or ()
        return min((position.get(name, unknown) for name in named), default=unknown)

    return sorted(solvers, key=rank)


def simulation(
    *models: ModelSpec,
    name: str = "sim",
    packages: Iterable[PackageEntry] = (),
    exchanges: Iterable[Any] = (),
    **options: Any,
) -> SimulationSpec:
    """Assemble a :class:`SimulationSpec`, defaulting the timing and solvers.

    The common case -- one or more models that should just run -- needs no
    boilerplate::

        mf.simulation(flow)              # steady, one IMS solving `flow`
        mf.simulation(flow, transport)   # one IMS each, coupled-ready

    Everything simulation-level goes through ONE ``packages=`` list, the way
    model packages go through ``mf.gwf(packages=[...])``::

        mf.simulation(flow, packages=[mf.tdis(nper=12, perioddata=spd)])
        mf.simulation(flow, transport, packages=[tdis, flow_ims, transport_ims])

    Whatever you supply, myflopy fills in the rest: a steady ``mf.tdis()`` when
    the list has no timing, and a solver for every model no supplied solver
    claims -- an ``IMS`` for GWF/GWT/GWE and an ``EMS`` for PRT, which MF6
    requires (it rejects an explicit model under IMS6).

    A solver with no ``models=`` claims every model the others left, so a
    one-model simulation needs no names at all::

        mf.simulation(flow, packages=[mf.ims(complexity="COMPLEX")])

    Solvers are ordered by MODEL DECLARATION order regardless of the order you
    list them: MF6 requires a GWF solution to appear in ``mfsim.nam`` before the
    solution of any model it feeds.

    Parameters
    ----------
    *models : ModelSpec
        One or more model specs (``mf.gwf(...)`` etc.). At least one is required.
    name : str, default "sim"
        Simulation name.
    packages : Iterable[PackageEntry], optional
        Simulation-level packages -- timing (``mf.tdis``) and solvers
        (``mf.ims`` / ``mf.ems``). A model package raises; it belongs on the
        model. A ``mf.ref`` is passed through untouched and does NOT suppress a
        default, because an unresolved reference cannot be identified -- supply
        the timing or solver it stands for explicitly alongside it.
    exchanges : Iterable, optional
        Inter-model exchange specs -- see :func:`exchange`.
    **options
        Additional :class:`SimulationSpec` fields (``workspace``, ``executable``,
        ``run_name``).

    Returns
    -------
    SimulationSpec

    Raises
    ------
    ValueError
        If no model specs are given, if ``packages`` carries more than one
        timing package or more than one unnamed solver, or if a model package is
        handed to the simulation tier.

    Examples
    --------
    >>> mf.simulation(flow)
    >>> mf.simulation(flow, packages=[mf.ims(complexity="COMPLEX")])
    >>> mf.simulation(flow, transport, packages=[tdis, flow_ims, transport_ims])
    """

    if not models:
        raise ValueError("simulation() requires at least one model spec.")

    timing, solvers, others = _partition_simulation_packages(packages)

    if len(timing) > 1:
        names = ", ".join(spec.name for spec in timing)
        raise ValueError(
            f"simulation() was given {len(timing)} timing packages: {names}. "
            "A simulation has exactly one TDIS."
        )
    if not timing:
        timing = [_steady_tdis()]

    solvers = _claim_models(solvers, models)
    return SimulationSpec(
        name,
        models=tuple(models),
        packages=(*timing, *_ordered_solvers(solvers, models), *others),
        exchanges=tuple(exchanges),
        **options,
    )


class _ExchangeHelper:
    """Package-first inter-model exchanges -- the coupling a multi-model run needs.

    ``mf.simulation(flow, transport)`` wires no exchange for you: the models are
    built and nothing moves between them, MF6 reports normal termination, and the
    transport model sees no flow field. The exchange is the piece that does it.

    ``mf.exchange(flow, transport)`` reads the two models' kinds and picks the MF6
    exchange class. The explicit spellings are there when you want to be plain, or
    need the options only GWF-GWF has::

        mf.exchange(flow, transport)                      # -> ModflowGwfgwt
        mf.exchange.gwfgwt(flow, transport)               # the same, said outright
        mf.exchange.gwfgwf(coarse, fine, nexg=12, exchangedata=rows)

    Three of the four carry no information beyond the model pair, which is why
    inference is safe for them. GWF-GWF is the exception -- it takes 18 options,
    and MF6 rejects one written without ``nexg``/``exchangedata`` ("Required
    block \"DIMENSIONS\" not found").

    Order does not matter: whichever model is the GWF is put first, because that
    is the order MF6 wants and getting it wrong is not a choice anyone makes on
    purpose.
    """

    _BUILDERS: ClassVar[dict[tuple[str, str], Any]] = {}

    def __call__(
        self,
        source: ModelSpec,
        target: ModelSpec,
        *,
        name: str | None = None,
        **options: Any,
    ) -> ExchangeSpec:
        """Return the exchange coupling ``source`` and ``target``, chosen by kind.

        Both must be :class:`ModelSpec` objects -- the kinds are what select the
        MF6 class, so a pair of names cannot be resolved. Use an explicit
        spelling (:meth:`gwfgwt` and friends) if all you have are names.

        Parameters
        ----------
        source : ModelSpec
            One model of the pair. Whichever of the two is the GWF model is
            written first (``exgmnamea``), so the order does not matter. A plain
            name string raises ``TypeError``.
        target : ModelSpec
            The other model. The pair must be GWF with GWT, GWE, PRT or GWF;
            any other pairing raises ``ValueError``.
        name : str, optional
            The exchange's name in ``SimulationSpec.exchanges``, which must be
            unique there. Default ``None`` uses ``"<first>-<second>"`` from the
            model names, GWF first.
        **options
            Passed to the chosen FloPy exchange class (``ModflowGwfgwt``,
            ``ModflowGwfgwe``, ``ModflowGwfprt`` or ``ModflowGwfgwf``) at build
            time; not checked at the call. A GWF-GWF pair needs ``nexg`` and
            ``exchangedata`` here (MF6 rejects the file without them), and this
            form does not check for them -- :meth:`gwfgwf` does.

        Returns
        -------
        ExchangeSpec
            Pass it to ``mf.simulation(..., exchanges=[...])``.
        """

        pair = tuple(_exchange_kind(model, "mf.exchange") for model in (source, target))
        models = (source, target)
        if pair[1] == "gwf" and pair[0] != "gwf":
            pair, models = pair[::-1], models[::-1]
        builder = self._BUILDERS.get(pair)
        if builder is None:
            supported = ", ".join(f"{a}-{b}" for a, b in sorted(self._BUILDERS))
            raise ValueError(
                f"No MF6 exchange couples a {pair[0]} model to a {pair[1]} model. "
                f"Supported: {supported}."
            )
        return self._spec(builder, models, name=name, **options)

    def gwfgwt(self, flow, transport, *, name=None, **options) -> ExchangeSpec:
        """Couple a flow model to a solute-transport model (``ModflowGwfgwt``).

        Parameters
        ----------
        flow : ModelSpec or str
            The GWF model, or its name. Written first (``exgmnamea``); the order
            is taken as given, not inferred.
        transport : ModelSpec or str
            The GWT model, or its name (``exgmnameb``).
        name : str, optional
            The exchange's name in ``SimulationSpec.exchanges``, unique there.
            Default ``None`` uses ``"<flow>-<transport>"``.
        **options
            Passed to ``flopy.mf6.ModflowGwfgwt`` at build time (e.g.
            ``filename``, ``pname``); not checked at the call.

        Returns
        -------
        ExchangeSpec
        """

        return self._spec(build_gwf_gwt_exchange, (flow, transport), name=name, **options)

    def gwfgwe(self, flow, energy, *, name=None, **options) -> ExchangeSpec:
        """Couple a flow model to an energy-transport model (``ModflowGwfgwe``).

        Parameters
        ----------
        flow : ModelSpec or str
            The GWF model, or its name. Written first (``exgmnamea``); the order
            is taken as given, not inferred.
        energy : ModelSpec or str
            The GWE model, or its name (``exgmnameb``).
        name : str, optional
            The exchange's name in ``SimulationSpec.exchanges``, unique there.
            Default ``None`` uses ``"<flow>-<energy>"``.
        **options
            Passed to ``flopy.mf6.ModflowGwfgwe`` at build time (e.g.
            ``filename``, ``pname``); not checked at the call.

        Returns
        -------
        ExchangeSpec
        """

        return self._spec(build_gwf_gwe_exchange, (flow, energy), name=name, **options)

    def gwfprt(self, flow, particles, *, name=None, **options) -> ExchangeSpec:
        """Couple a flow model to a particle-tracking model (``ModflowGwfprt``).

        Parameters
        ----------
        flow : ModelSpec or str
            The GWF model, or its name. Written first (``exgmnamea``); the order
            is taken as given, not inferred.
        particles : ModelSpec or str
            The PRT model, or its name (``exgmnameb``).
        name : str, optional
            The exchange's name in ``SimulationSpec.exchanges``, unique there.
            Default ``None`` uses ``"<flow>-<particles>"``.
        **options
            Passed to ``flopy.mf6.ModflowGwfprt`` at build time (e.g.
            ``filename``, ``pname``); not checked at the call.

        Returns
        -------
        ExchangeSpec
        """

        return self._spec(build_gwf_prt_exchange, (flow, particles), name=name, **options)

    def gwfgwf(self, first, second, *, name=None, **options) -> ExchangeSpec:
        """Couple two flow models across a shared interface (``ModflowGwfgwf``).

        Unlike the other three this carries real content, and MF6 refuses an
        exchange written without it: ``nexg`` and ``exchangedata`` describe the
        connected cell pairs, and a file holding only an empty options block
        fails with *Required block "DIMENSIONS" not found*.

        Parameters
        ----------
        first : ModelSpec or str
            One GWF model, or its name (``exgmnamea``); ``exchangedata``'s first
            cellid in each row is in this model.
        second : ModelSpec or str
            The other GWF model, or its name (``exgmnameb``).
        name : str, optional
            The exchange's name in ``SimulationSpec.exchanges``, unique there.
            Default ``None`` uses ``"<first>-<second>"``.
        **options
            Passed to ``flopy.mf6.ModflowGwfgwf`` at build time. ``nexg`` and
            ``exchangedata`` are REQUIRED -- either missing raises ``ValueError``
            at the call; the rest (``auxiliary``, ``xt3d``, ``newton``,
            ``perioddata``, ...) are optional and not checked at the call.

        Returns
        -------
        ExchangeSpec
        """

        missing = [key for key in ("nexg", "exchangedata") if key not in options]
        if missing:
            raise ValueError(
                "A GWF-GWF exchange needs "
                f"{' and '.join(missing)}: MF6 reads the connected cell pairs from "
                'its DIMENSIONS block and rejects a file without them ("Required '
                'block \"DIMENSIONS\" not found").'
            )
        return self._spec(build_gwf_gwf_exchange, (first, second), name=name, **options)

    @staticmethod
    def _spec(builder, models, *, name=None, **options) -> ExchangeSpec:
        """Name the pair and wrap it, so callers never touch ExchangeSpec."""

        names = tuple(_exchange_name(model) for model in models)
        return ExchangeSpec(
            name or "-".join(names), builder, models=names, options=options
        )


_ExchangeHelper._BUILDERS = {
    ("gwf", "gwt"): build_gwf_gwt_exchange,
    ("gwf", "gwe"): build_gwf_gwe_exchange,
    ("gwf", "prt"): build_gwf_prt_exchange,
    ("gwf", "gwf"): build_gwf_gwf_exchange,
}

exchange = _ExchangeHelper()


def _exchange_name(model: ModelSpec | str) -> str:
    """The model NAME an exchange records, from a spec or a plain name."""

    return model if isinstance(model, str) else model.name


def _exchange_kind(model: ModelSpec | str, verb: str) -> str:
    """The model KIND, refusing a bare name -- inference needs the kind."""

    if isinstance(model, str):
        raise TypeError(
            f"{verb}() needs the model SPECS, not their names: the kinds are what "
            f"choose the MF6 exchange class. Got {model!r}. Pass the objects, or "
            f"name the exchange outright with {verb}.gwfgwt(...) and friends."
        )
    return str(getattr(model.model_type, "value", model.model_type))


def _build_lazy_disv(model, *, layer_stack, build_kwargs, disv_name="disv", **options):
    """Resolve a declared ``LayerStack`` against the model's grid, then build the DISV.

    The whole of lazy discretization. ``ModelSpec.build`` resolves a deferred
    ``GridSpec`` and publishes it onto ``model.myflopy_context`` BEFORE any package
    builds, so by the time this runs the grid exists -- which is why laziness here
    needed no new machinery, only this function in place of a raise.

    It PUBLISHES what it resolved: the stack's ``idomain`` becomes the context's
    domain and the elevations it built -- a geometry-first table made from the RESULT,
    carrying its excavation records in ``attrs["excavations"]`` -- become its surfaces,
    so the surface-aware builders (``mf.sfr`` reach tops, ``mf.lak`` layering) and the
    build-time guard all see the mesh this DISV actually describes. Not
    ``grid.gdf_topbtm``: with ``attach=False`` that is whatever was published before.
    Published by REBINDING, because ``ModelContext`` is frozen with slots and an
    in-place write would raise.

    A DECLARED domain that disagrees with what the stack resolved is refused rather
    than overwritten: the caller described two different models and only one of them
    can be built (ledger 270). Declared SURFACES are not compared -- they are replaced
    by what was built. The published table is marked as the DISV's own
    (``attrs[specs._PUBLISHED_BY_KEY]``), so ``specs._layer_elevations`` does not warn that
    it disagrees with a grid table ``attach=False`` left alone.
    """

    context = getattr(model, "myflopy_context", None)
    grid = getattr(context, "grid", None)
    if grid is None:
        raise ValueError(
            "mf.disv() was given a layer stack with no grid, so the grid has to come "
            "from the model at build time -- but this model's context carries none. "
            "Either declare the model with `mf.gwf(..., context=ModelContext(grid=...))`, "
            "give it a deferred grid with `.with_grid(mf.GridSpec.voronoi(...))`, or "
            "bind the stack itself with `stack.for_grid(vor)` before `mf.disv`."
        )

    result = layer_stack.build(grid, **build_kwargs)
    spec = result.to_disv(grid, name=disv_name)
    built = spec.with_options(**options).build(model) if options else spec.build(model)

    declared_domain = getattr(context, "domain", None)
    if declared_domain is not None:
        declared = np.asarray(declared_domain)
        resolved = np.asarray(result.idomain)
        if declared.shape != resolved.shape or not np.array_equal(declared, resolved):
            raise ValueError(
                f"The domain declared on model {getattr(model, 'name', '?')!r} is not the "
                f"one this layer stack resolves: declared {declared.shape}, the stack "
                f"produced {resolved.shape}. A lazy DISV publishes the domain it built, "
                "so declaring a different one describes two models. Drop the declared "
                "`domain=` and let the stack own it, or build the stack yourself and pass "
                "`mf.disv(result)`."
            )
    # The elevations of THIS result, not `grid.gdf_topbtm`: with `attach=False` the
    # grid still holds whatever was published before, and the context published a top of
    # 500.0 beside a DISV top of 100.0. The excavation records ride on the table itself
    # (`attrs`), so the elevations and the records a model reads come from one result
    # (`grid_identity.excavation_cells`; ledger 323).
    table = result._elevation_table(grid)
    table.attrs["excavations"] = dict(getattr(result, "excavations", None) or {})
    # Marked as the DISV's own: with attach=False the grid keeps an older table, and
    # `_layer_elevations` would otherwise warn on every read that the CONTEXT is stale
    # and advise rebuilding it -- which would put the old elevations back.
    table.attrs[_PUBLISHED_BY_KEY] = _PUBLISHED_BY_DISV
    model.myflopy_context = _replace(context, domain=result.idomain, surfaces=table)
    # Record what was actually built, so the run freezes ARRAYS rather than this
    # recipe. Freezing the declared spec stores the `LayerStack`, and a reopened run
    # then goes back to `gis/ground.tif` -- the dependence ledger 251's snapshot
    # exists to remove. `ModelSpec.build` collects this into `BuiltModel.resolved`.
    resolved = getattr(model, "_myflopy_resolved", None)
    if resolved is None:
        resolved = {}
        model._myflopy_resolved = resolved
    resolved[disv_name] = spec.with_options(**options) if options else spec
    return built


# The lazy builder targets exactly what the eager one does, and saying so is what
# restores the DECLARATION-time checks every other package gets. Without it
# `_accepted_options` reads a plain function, returns `None`, and abstains: a typo in
# `mf.disv(stack, ...)` was accepted at the call and surfaced from inside FloPy at
# build. It also makes the builder's SCOPE decidable, so a lazy DISV handed to
# `mf.simulation(packages=[...])` is refused at the call like any other model package.
_build_lazy_disv.mf6_targets = build_disv.mf6_targets


_DISV_BUILD_DEFAULTS = {
    "attach": True,
    "default_min_thickness": 1.0,
    "default_pinch": "floor",
    "reconcile": "bottom",
    "min_sep": 0.1,
    "trigger_sep": None,
    "method": "area",
    "refresh": False,
}


def _disv_from_layers(layers: Any, grid: Any, *, name: str, build_kwargs: dict[str, Any]):
    """Delegate to whichever ``to_disv`` ``layers`` has, forwarding what it accepts.

    Two different objects answer to ``layers``. A ``LayerStack`` has not sampled
    anything yet, so every build kwarg acts. A ``LayerBuildResult`` already HOLDS the
    resolved arrays -- its ``to_disv`` takes only ``name`` -- so a build kwarg there
    cannot act, and silently accepting one is the same defect as a hidden parameter
    (the `surface(backend=)` precedent). It raises instead, naming the ones passed and
    saying where they belong.

    `myflopy.layers` is import layer 6 and this module is layer 5, so the class is not
    imported: the distinguishing feature is which signature the object's own
    ``to_disv`` has.
    """

    import inspect

    accepted = set(inspect.signature(layers.to_disv).parameters)
    non_default = {
        key: value
        for key, value in build_kwargs.items()
        if value != _DISV_BUILD_DEFAULTS[key]
    }
    if "reconcile" not in accepted:
        if non_default:
            raise TypeError(
                f"mf.disv() cannot apply {', '.join(sorted(non_default))} to an "
                "already-built layer result: its arrays were resolved when you called "
                "`stack.build(...)`, so passing them here would describe a different "
                "model than the one you inspected. Pass them to `stack.build(...)` or "
                "`stack.to_disv(...)`, or hand `mf.disv` the STACK instead of its result."
            )
        return layers.to_disv(grid, name=name) if grid is not None else layers.to_disv(name=name)
    if grid is not None:
        return layers.to_disv(grid, name=name, **build_kwargs)
    return layers.to_disv(name=name, **build_kwargs)


def disv(
    layers: Any = None,
    *,
    grid: Any = None,
    nlay: int | None = None,
    ncpl: int | None = None,
    nvert: int | None = None,
    vertices: Any = None,
    cell2d: Any = None,
    top: Any = None,
    botm: Any = None,
    idomain: Any | None = None,
    name: str = "disv",
    attach: bool = True,
    default_min_thickness: float = 1.0,
    default_pinch: str = "floor",
    reconcile: Any = "bottom",
    min_sep: float = 0.1,
    trigger_sep: float | None = None,
    method: str = "area",
    refresh: bool = False,
    **options: Any,
) -> PackageSpec:
    """Vertex discretization (DISV) package: the unstructured grid + layers.

    Every other package in a model comes from ``mf.<name>(...)``; this one used
    to be the exception, reachable only as ``result.to_disv()``. Pass the built
    layer result and it is not::

        layers = stack.build(vor, attach=True)
        mf.disv(layers)                      # the mesh comes with it

    Three forms, told apart by ``layers``:

    - **A built ``LayerBuildResult``** -- EAGER, and exactly ``result.to_disv()``:
      the arrays you inspected with ``qc()`` are the arrays written.
    - **A ``LayerStack``** -- EAGER when it is bound to a grid (``LayerStack(vor,
      ...)``, ``stack.for_grid(vor)``) or ``grid=`` is given, resolved here with
      the build keywords below. LAZY when neither is: the spec holds the stack and
      resolves it at BUILD against whatever grid the model supplies, so one
      declaration builds on any mesh and a layered model can be saved (ledger 270).
    - **``None``** -- hand-written: ``nlay`` through ``botm`` (and optionally
      ``idomain``), typically from ``vor.get_disv_gridprops()``.

    Parameters
    ----------
    layers : LayerStack or LayerBuildResult, optional
        The layering to discretize; see the three forms above. Anything else with no
        ``to_disv`` raises ``AttributeError``. A LAZY spec holds its OWN copy of the
        stack, rebuilt from the stack's payload the way a reload rebuilds it (ledger
        270, amended 2026-09-25): a later ``stack.add(...)``/``replace(...)`` does
        not reach a package already declared, a library copy of it, or its digest --
        declare it again to take an edit. The copy digests as the original. A stack
        holding per-cell arrays (``Surface.from_array``) cannot be written down, so it
        is copied with ``for_grid(None)`` instead. An EAGER spec from a
        ``LayerBuildResult`` holds COPIES of its ``top``/``botm``/``idomain``, so an
        in-place edit of the result after the call (``result.idomain[0, :5] = 0``)
        does not reach it either -- it did, measured, because ``np.asarray`` shares
        memory.
    grid : VoronoiGridPlus, optional
        The mesh to discretize a gridless stack or result on, which makes the form
        EAGER. Wins over the grid the stack or result carries. A STACK resolves on it
        whatever its size; a RESULT's arrays are already per-cell, so a grid with a
        different cell count RAISES ``ValueError`` -- it used to write the new
        ``ncpl`` (92) over a 28-cell ``top``, which FloPy refused only at build.
        Ignored by the hand-written form.
    nlay : int, optional
        Hand-written form: number of layers.
    ncpl : int, optional
        Hand-written form: number of cells per layer.
    nvert : int, optional
        Hand-written form: number of vertices.
    vertices : list of [int, float, float] or ndarray, optional
        Hand-written form: FloPy's vertex list ``[[iv, x, y], ...]``, zero-based,
        coordinates in model length units.
    cell2d : list of [int, float, float, int, int, ...], optional
        Hand-written form: FloPy's cell list ``[[icell2d, xc, yc, ncvert, iv1, ...],
        ...]`` -- ``vor.get_disv_gridprops()["cell2d"]``.
    top : float or array_like of shape (ncpl,), optional
        Hand-written form: model-top elevation, length units.
    botm : array_like of shape (nlay, ncpl), optional
        Hand-written form: layer bottom elevations, length units, top layer first.
    idomain : array_like of int, shape (nlay, ncpl), optional
        Hand-written form: ``1`` active, ``0`` inactive, ``-1`` vertical
        pass-through. Omitted, nothing is written and MF6 treats every cell as
        active. The layered forms derive it from the pinch policy instead.
    name : str, default "disv"
        The SPEC's name -- the key ``with_packages(slot=)``, ``without_packages`` and
        a run's frozen snapshot use. It does not reach FloPy: the built package is
        always named ``disv`` and written to ``<model>.disv`` (measured with
        ``name="dis2"``).
    attach : bool, default True
        Layer-stack forms: publish the sampled elevations onto the grid's
        ``gdf_topbtm`` (measured: set with ``True``, left unset with ``False``), which
        ``mf.sfr`` reach tops, ``mf.lak`` lake-cell layering and the layer-elevation
        hover rows read. Note ``LayerStack.build`` defaults it to ``False``. On the
        LAZY form the MODEL's context gets the elevations the stack just built
        either way, with its excavation records (``attrs["excavations"]``); this
        decides only whether the GRID is published to as well. Until ledger 323
        ``attach=False`` handed the context whatever table the grid already held --
        measured, a top of 500.0 against the DISV's 100.0.
    default_min_thickness : float, default 1.0
        Layer-stack forms: minimum thickness, length units, for units that declare
        none with ``.add(..., min_thickness=)``.
    default_pinch : {"floor", "inactive", "passthrough"}, default "floor"
        Layer-stack forms: thin-layer policy for units that declare none --
        ``"floor"`` holds the cell open at its minimum thickness, ``"inactive"``
        writes idomain 0, ``"passthrough"`` idomain -1 (flow passes vertically and
        the cell can carry no boundary).
    reconcile : {"bottom", "top", True, False, None}, default "bottom"
        Layer-stack forms: what to do where two surfaces cross or come closer than
        ``trigger_sep`` -- ``"bottom"``/``True`` push the lower one down (the model
        top stays where you put it), ``"top"`` raises the upper one (and can move the
        model top), ``False``/``None`` leave them crossed. Full text on
        :meth:`LayerStack.build`.
    min_sep : float, default 0.1
        Layer-stack forms: the separation, length units, reconcile leaves between a
        PINCHING unit's surfaces; a ``"floor"`` unit uses its own minimum thickness.
    trigger_sep : float, optional
        Layer-stack forms: the separation at or below which reconcile acts at all,
        length units. Omitted, ``1.0 / <largest split>``. Never reaches a pinching
        unit at or above its own ``min_thickness``; full text on
        :meth:`LayerStack.build`.
    method : {"area", "centroid"}, default "area"
        Layer-stack forms: how a raster is sampled onto a cell -- area-weighted, or
        the value at its centroid.
    refresh : bool, default False
        Layer-stack forms: rebuild a stale contour surface's GRASS cache before
        sampling; left ``False``, a stale cache is reused with a warning.
    **options
        Any other ``ModflowGwfdisv`` keyword (``length_units``, ``xorigin``,
        ``yorigin``, ``angrot``, ...), applied on every form; one it does not have
        raises here, at the call, lazy form included.

    Returns
    -------
    PackageSpec
        EAGER: ``options`` hold the mesh and the ``top``/``botm``/``idomain``
        arrays, ``inputs`` is empty. LAZY: ``inputs`` holds ``layer_stack`` (the
        package's own copy), ``build_kwargs`` and ``disv_name``, and the arrays
        exist only once a model builds it.

    Raises
    ------
    TypeError
        ``layers`` together with any of ``nlay``/``ncpl``/``nvert``/``vertices``/
        ``cell2d``/``top``/``botm`` (two meshes); the hand-written form missing any
        of those seven; a non-default build keyword on a ``LayerBuildResult``; an
        ``**options`` keyword ``ModflowGwfdisv`` does not accept.
    AttributeError
        ``layers`` that is neither a stack nor a result (it has no ``to_disv``).
    ValueError
        A ``LayerBuildResult`` with a ``grid=`` whose cell count differs from its
        own. At BUILD, for the lazy form: a model whose
        context carries no grid, or a declared ``domain`` that disagrees with the one
        the stack resolves.

    Notes
    -----
    The eight build keywords ``attach`` through ``refresh`` change the arrays a stack
    resolves to -- a stack does not determine its own arrays. On a
    ``LayerBuildResult`` they cannot act, so a non-default one RAISES; the
    hand-written form does NOT look at them (found 2026-09-25, not changed: a
    ``reconcile=False`` there is accepted and does nothing).

    Examples
    --------
    >>> import myflopy as mf
    >>> stack = mf.LayerStack(top=mf.Surface.flat(100.0)).add("sand", thickness=20.0)
    >>> spec = mf.disv(stack)                 # gridless -> LAZY
    >>> sorted(spec.inputs)
    ['build_kwargs', 'disv_name', 'layer_stack']
    >>> _ = stack.add("clay", thickness=30.0)  # the declared package keeps its copy
    >>> [layer.name for layer in spec.inputs["layer_stack"]._layers]
    ['sand']
    >>> one_cell = mf.disv(nlay=1, ncpl=1, nvert=4,
    ...                    vertices=[[0, 0.0, 0.0], [1, 0.0, 10.0],
    ...                              [2, 10.0, 10.0], [3, 10.0, 0.0]],
    ...                    cell2d=[[0, 5.0, 5.0, 4, 0, 1, 2, 3]],
    ...                    top=10.0, botm=[0.0])
    >>> one_cell.options["ncpl"]
    1
    """

    if layers is not None:
        # Delegate rather than re-derive: the result's own `to_disv` is what
        # guarantees the DISV matches the arrays you inspected with `qc()`.
        # Re-deriving here would be a SECOND resolution of the same stack, which
        # measured 8.5 ft of drift on a real model (see LayerStack.to_disv).
        supplied = {
            key: value
            for key, value in {
                "nlay": nlay, "ncpl": ncpl, "nvert": nvert, "vertices": vertices,
                "cell2d": cell2d, "top": top, "botm": botm,
            }.items()
            if value is not None
        }
        if supplied:
            raise TypeError(
                "mf.disv(layers) takes the geometry FROM the layer result; "
                f"also passing {', '.join(sorted(supplied))} would describe a "
                "different mesh. Pass one or the other."
            )
        # Every one of `to_disv`'s build kwargs is forwarded. It forwarded NONE of
        # them, so `mf.disv(stack, default_pinch="inactive")` raised
        # `TypeError: ModflowGwfdisv does not accept 'default_pinch'` and a stack
        # silently used `to_disv`'s defaults. That matters because **a stack does not
        # determine its own arrays**: same stack, same grid, a different `reconcile`
        # measures a different `botm`, and `LayerBuildResult` records no
        # `reconcile`/`min_sep`/`trigger_sep`/`method` at all. Once the DISV resolves
        # lazily on another mesh, forwarding nothing would re-derive with defaults --
        # a different idomain, no warning (ledger 270).
        build_kwargs = {
            "attach": attach,
            "default_min_thickness": default_min_thickness,
            "default_pinch": default_pinch,
            "reconcile": reconcile,
            "min_sep": min_sep,
            "trigger_sep": trigger_sep,
            "method": method,
            "refresh": refresh,
        }
        reachable = grid if grid is not None else getattr(layers, "vor", None)
        if reachable is None and hasattr(layers, "to_disv"):
            # A declared, gridless LayerStack: emit a LAZY spec instead of raising.
            # Duck-typed on `to_disv`, NOT on `_SPEC_TOKEN`: the token means "this
            # serializes", which stopped being a synonym for "this is a layer stack"
            # the moment ledger 274 gave one to `Spread` and `CellSurfaceOffset`.
            # Measured then: `mf.disv(mf.Spread("x"))` turned from an immediate
            # `AttributeError` naming the real problem into an ACCEPTED, persistable
            # spec that died at build talking about layer stacks.
            # The stack rides in `inputs`, which `PackageSpec.build` splats into the
            # builder and which `_json_value` and `freeze._lift` already handle -- so
            # no new `recipe` field is needed (ledger 270). `package_digest` did NOT
            # handle it: a two-layer stack digested one way live and another after a
            # reload, so every such package read as EDITED forever. Fixed in 273.
            #
            # The package holds its OWN copy of the stack (ledger 270, amended): it
            # held the caller's object whenever nothing was stamped, so a later
            # `stack.add(...)` rewrote a package already declared. `stamp_source_hashes`
            # rebuilds every stack it can write down; one holding per-cell arrays
            # (`Surface.from_array`) cannot be written down, so it is copied here.
            held = stamp_source_hashes(layers)
            if held is layers and hasattr(layers, "for_grid"):
                held = layers.for_grid(None)
            return PackageSpec(
                name,
                _build_lazy_disv,
                dict(options),
                requires=(),
                inputs={
                    "layer_stack": held,
                    "build_kwargs": build_kwargs,
                    "disv_name": name,
                },
            )
        spec = _disv_from_layers(layers, grid, name=name, build_kwargs=build_kwargs)
        return spec.with_options(**options) if options else spec

    missing = [
        key
        for key, value in {
            "nlay": nlay, "ncpl": ncpl, "nvert": nvert,
            "vertices": vertices, "cell2d": cell2d, "top": top, "botm": botm,
        }.items()
        if value is None
    ]
    if missing:
        raise TypeError(
            f"mf.disv() is missing {', '.join(missing)}. Pass a built layer "
            "result instead -- mf.disv(stack.build(vor)) -- and the mesh comes "
            "with it."
        )

    values = {
        "nlay": nlay,
        "ncpl": ncpl,
        "nvert": nvert,
        "vertices": vertices,
        "cell2d": cell2d,
        "top": top,
        "botm": botm,
        **options,
    }
    if idomain is not None:
        values["idomain"] = idomain
    return PackageSpec(name, build_disv, values)


def dis(
    *,
    nlay: int,
    nrow: int,
    ncol: int,
    delr: Any,
    delc: Any,
    top: Any,
    botm: Any,
    idomain: Any | None = None,
    name: str = "dis",
    **options: Any,
) -> PackageSpec:
    """Structured (row/column) discretization -- DIS.

    myflopy is Voronoi/DISV-first; ``mf.dis`` exists so a plain rectilinear grid can
    be built with the same package-first grammar, and to interoperate with externally
    supplied structured models. A DIS model builds and runs, but the choropleth-map /
    cross-section / animation viz targets DISV/Voronoi meshes and is **not** wired for
    structured grids -- use FloPy's own plotting, or DISV, for those. Like ``mf.disv``
    it dispatches the FloPy class off the model kind (GWF/GWT/GWE/PRT).

    Parameters
    ----------
    nlay : int
        Number of layers.
    nrow : int
        Number of rows.
    ncol : int
        Number of columns.
    delr : float or array_like of shape (ncol,)
        Column widths along a row, model length units; a scalar is one width for
        every column.
    delc : float or array_like of shape (nrow,)
        Row heights along a column, model length units; a scalar is one height for
        every row.
    top : float or array_like of shape (nrow, ncol)
        Model-top elevation, length units; a scalar is a flat top.
    botm : float or array_like of shape (nlay, nrow, ncol)
        Layer bottom elevations, length units, top layer first; a scalar or a
        per-layer list is constant within each layer (FloPy's array input).
    idomain : array_like of int, shape (nlay, nrow, ncol), optional
        ``1`` active, ``0`` inactive, ``-1`` vertical pass-through. Default
        ``None`` writes nothing, and MF6 treats every cell as active.
    name : str, optional
        The spec's name in the model's package list. It is not passed to FloPy.
        Default ``"dis"``.
    **options
        Extra ``flopy.mf6.Modflow<Kind>dis`` options (e.g. ``length_units``,
        ``xorigin``/``yorigin``/``angrot``), passed at build time; one the DIS
        class does not accept raises ``TypeError`` at the call.

    Returns
    -------
    PackageSpec
    """

    values = {
        "nlay": nlay,
        "nrow": nrow,
        "ncol": ncol,
        "delr": delr,
        "delc": delc,
        "top": top,
        "botm": botm,
        **options,
    }
    if idomain is not None:
        values["idomain"] = idomain
    return PackageSpec(name, build_dis, values)


def disu(
    *,
    nodes: int,
    nja: int,
    top: Any,
    bot: Any,
    area: Any,
    iac: Any,
    ja: Any,
    idomain: Any | None = None,
    name: str = "disu",
    **options: Any,
) -> PackageSpec:
    """Fully unstructured discretization -- DISU (connectivity you supply yourself).

    The most general MF6 grid: an explicit node list plus a hand-built connectivity
    table (``iac``/``ja`` and, via ``**options``, ``ihc``/``cl12``/``hwva``). myflopy
    is Voronoi/DISV-first -- prefer ``mf.disv`` with a ``VoronoiGridPlus``, which
    generates that connectivity for you; reach for ``mf.disu`` only when ingesting an
    existing DISU mesh. As with ``mf.disv`` the choropleth / cross-section / animation
    viz targets DISV/Voronoi meshes and is **not** wired for raw DISU. Dispatches the
    FloPy class off the model kind (GWF/GWT/GWE); MF6 has no ``ModflowPrtdisu``.

    Parameters
    ----------
    nodes : int
        Number of cells.
    nja : int
        Length of ``ja``: every cell's connections counted from both sides, plus
        one entry per cell for the cell itself.
    top : float or array_like of shape (nodes,)
        Top elevation of each cell, model length units.
    bot : float or array_like of shape (nodes,)
        Bottom elevation of each cell, length units.
    area : float or array_like of shape (nodes,)
        Plan-view area of each cell, length units squared.
    iac : array_like of int, shape (nodes,)
        Number of connections of each cell PLUS ONE (the cell itself); sums to
        ``nja``.
    ja : array_like of int, shape (nja,)
        For each cell in turn: the cell's own number followed by the cells it
        connects to. Zero-based -- FloPy adds one on write.
    idomain : array_like of int, shape (nodes,), optional
        ``1`` active, ``0`` inactive. Default ``None`` writes nothing, and MF6
        treats every cell as active.
    name : str, optional
        The spec's name in the model's package list. It is not passed to FloPy.
        Default ``"disu"``.
    **options
        Extra ``flopy.mf6.Modflow<Kind>disu`` options, passed at build time. MF6
        also REQUIRES ``ihc``/``cl12``/``hwva`` (connection direction, per-side
        connection length, width or area), and they go here; this helper does not
        check for them. Also ``nvert``/``vertices``/``cell2d``, ``length_units``,
        ``xorigin``/``angrot``, .... One the DISU class does not accept raises
        ``TypeError`` at the call.

    Returns
    -------
    PackageSpec
    """

    values = {
        "nodes": nodes,
        "nja": nja,
        "top": top,
        "bot": bot,
        "area": area,
        "iac": iac,
        "ja": ja,
        **options,
    }
    if idomain is not None:
        values["idomain"] = idomain
    return PackageSpec(name, build_disu, values)


def ic(*, strt: Any, name: str = "ic", **options: Any) -> PackageSpec:
    """Initial-conditions (IC) package: the starting dependent-variable field.

    Dispatches on the model kind: starting head (GWF), concentration (GWT), or
    temperature (GWE).

    Parameters
    ----------
    strt : float or array-like
        Starting value -- a scalar, a per-cell array ``(ncpl,)``, or a
        ``(nlay, ncpl)`` array. For GWF a sensible head (near the water table)
        helps Newton / under-relaxation runs converge.
    name : str, default "ic"
        Package name.
    **options
        Extra ``flopy.mf6.ModflowGw{f,t,e}ic`` options.

    Returns
    -------
    PackageSpec

    Examples
    --------
    >>> mf.ic(strt=9.0)                 # uniform
    >>> mf.ic(strt=starting_heads)      # per-cell / per-layer array
    """

    return _array_package_spec("ic", name, build_ic, {"strt": strt, **options})


def npf(
    *,
    k: Any,
    k33: Any | None = None,
    icelltype: Any | None = None,
    save_flows: bool = True,
    save_specific_discharge: bool = False,
    name: str = "npf",
    **options: Any,
) -> PackageSpec:
    """Node-property-flow (NPF) package: hydraulic conductivity.

    Parameters
    ----------
    k : float or array-like
        Horizontal hydraulic conductivity -- a scalar, per-cell array, or
        ``(nlay, ncpl)`` array.
    k33 : float or array-like, optional
        Vertical (z) hydraulic conductivity; defaults to isotropic (``k``) in MF6.
    icelltype : int or array-like, optional
        Cell type: 0 = confined, non-zero = convertible (water-table). Use a
        convertible type for unconfined / wetting-and-drying layers.
    save_flows : bool, default True
        Save cell-by-cell flows.
    save_specific_discharge : bool, default False
        Save the Darcy velocity -- required for particle tracking / velocity plots.
    name : str, default "npf"
        Package name.
    **options
        Extra ``flopy.mf6.ModflowGwfnpf`` options.

    Returns
    -------
    PackageSpec

    Examples
    --------
    >>> mf.npf(k=10.0)                       # uniform K
    >>> mf.npf(k=k_array, k33=k_array * 0.1, icelltype=1)   # per-cell, convertible
    """

    values = {
        "k": k,
        "save_flows": save_flows,
        "save_specific_discharge": save_specific_discharge,
        **options,
    }
    if k33 is not None:
        values["k33"] = k33
    if icelltype is not None:
        values["icelltype"] = icelltype
    return _array_package_spec("npf", name, flopy.mf6.ModflowGwfnpf, values)


def sto(
    *,
    ss: Any | None = None,
    sy: Any | None = None,
    iconvert: Any | None = None,
    steady_state: Any | None = None,
    transient: Any | None = None,
    save_flows: bool = True,
    name: str = "sto",
    **options: Any,
) -> PackageSpec:
    """Storage (STO) package: storativity + per-period steady/transient flags.

    A purely steady model still needs STO with ``steady_state={0: True}``.

    Parameters
    ----------
    ss : float or array-like, optional
        Specific storage (used in transient periods).
    sy : float or array-like, optional
        Specific yield (used in transient periods for convertible cells).
    iconvert : int or array-like, optional
        Convertible-storage flag per cell (non-zero enables ``sy``).
    steady_state : dict, optional
        ``{period: True}`` marking steady-state stress periods.
    transient : dict, optional
        ``{period: True}`` marking transient stress periods.
    save_flows : bool, default True
        Save cell-by-cell storage flows.
    name : str, default "sto"
        Package name.
    **options
        Extra ``flopy.mf6.ModflowGwfsto`` options.

    Returns
    -------
    PackageSpec

    Examples
    --------
    >>> mf.sto(steady_state={0: True})                          # steady
    >>> mf.sto(ss=1e-5, sy=0.15, transient={0: True})           # transient
    """

    values = {"save_flows": save_flows, **options}
    for key, value in {
        "ss": ss,
        "sy": sy,
        "iconvert": iconvert,
        "steady_state": steady_state,
        "transient": transient,
    }.items():
        if value is not None:
            values[key] = value
    return _array_package_spec("sto", name, flopy.mf6.ModflowGwfsto, values)


def oc(
    *,
    head_filerecord: Any | None = None,
    budget_filerecord: Any | None = None,
    saverecord: Any | None = None,
    printrecord: Any | None = None,
    name: str = "oc",
    **options: Any,
) -> PackageSpec:
    """Output-control (OC) package: what to save/print and when.

    Dispatches on the model kind (GWF/GWT/GWE/PRT). To get the dependent variable
    and a cell budget on disk you must name the output files **and** request them
    in ``saverecord`` -- MF6 errors if you ask to save without an output file.
    For PRT pass ``track_filerecord=``/``trackcsv_filerecord=`` via ``**options``
    (the track CSV is what the results tier reads).

    Parameters
    ----------
    head_filerecord : str, optional
        GWF output head file name (e.g. ``"model.hds"``). For GWT/GWE pass
        ``concentration_filerecord=`` / ``temperature_filerecord=`` via
        ``**options`` instead (and the matching ``saverecord`` text).
    budget_filerecord : str, optional
        Output cell-budget file name (e.g. ``"model.cbc"``).
    saverecord : list, optional
        What/when to save, e.g. ``[("HEAD", "ALL"), ("BUDGET", "LAST")]``
        (use ``"CONCENTRATION"`` / ``"TEMPERATURE"`` for GWT/GWE).
    printrecord : list, optional
        What/when to print to the listing file.
    name : str, default "oc"
        Package name.
    **options
        Extra ``flopy.mf6.ModflowGw{f,t,e}oc`` options (incl. the GWT/GWE
        ``concentration_filerecord`` / ``temperature_filerecord``).

    Returns
    -------
    PackageSpec

    Examples
    --------
    >>> mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
    ...       saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])
    """

    values = dict(options)
    for key, value in {
        "head_filerecord": head_filerecord,
        "budget_filerecord": budget_filerecord,
        "saverecord": saverecord,
        "printrecord": printrecord,
    }.items():
        if value is not None:
            values[key] = value
    return PackageSpec(name, build_oc, values)


# --- GWT (solute transport) packages -------------------------------------------
# Thin factories, same shape as mf.ic/mf.npf: a PackageSpec around the FloPy class.
# The GeoPackage/registry/results-tier treatment of the transport list BCs
# (cnc/ctp input maps + hover, model.conc/temp) is Phase 6, not here (plan 5.3B).


def adv(*, scheme: str | None = None, name: str = "adv", **options: Any) -> PackageSpec:
    """Advection (ADV) package for a GWT model.

    Parameters
    ----------
    scheme : str, optional
        Advection scheme: ``"central"``, ``"upstream"``, ``"tvd"``, or ``"utvd"``
        (MF6 defaults to upstream). TVD reduces numerical dispersion at more cost.
    name : str, default "adv"
        Package name.
    **options
        Extra ``flopy.mf6.ModflowGwtadv`` options.
    """

    values = dict(options)
    if scheme is not None:
        values["scheme"] = scheme
    return PackageSpec(name, flopy.mf6.ModflowGwtadv, values)


def dsp(
    *,
    alh: Any | None = None,
    ath1: Any | None = None,
    diffc: Any | None = None,
    name: str = "dsp",
    **options: Any,
) -> PackageSpec:
    """Dispersion (DSP) package for a GWT model.

    Parameters
    ----------
    alh : float or array-like, optional
        Longitudinal dispersivity in the horizontal direction.
    ath1 : float or array-like, optional
        Transverse dispersivity in the first horizontal direction.
    diffc : float or array-like, optional
        Effective molecular diffusion coefficient.
    name : str, default "dsp"
        Package name.
    **options
        Extra ``flopy.mf6.ModflowGwtdsp`` options (``alv``/``ath2``/``atv``, xt3d flags).
    """

    values = dict(options)
    for key, value in {"alh": alh, "ath1": ath1, "diffc": diffc}.items():
        if value is not None:
            values[key] = value
    return PackageSpec(name, flopy.mf6.ModflowGwtdsp, values)


def mst(*, porosity: Any, name: str = "mst", **options: Any) -> PackageSpec:
    """Mobile storage and transfer (MST) package for a GWT model.

    Parameters
    ----------
    porosity : float or array-like
        Mobile-domain porosity (required).
    name : str, default "mst"
        Package name.
    **options
        Extra ``flopy.mf6.ModflowGwtmst`` options -- ``decay``/``sorption``/
        ``bulk_density``/``distcoef`` for decay and reactive transport.
    """

    return PackageSpec(name, flopy.mf6.ModflowGwtmst, {"porosity": porosity, **options})


def ist(
    *,
    porosity: Any,
    volfrac: Any,
    zetaim: Any,
    name: str = "ist",
    **options: Any,
) -> PackageSpec:
    """Immobile storage and transfer (IST) package for a GWT model (dual-domain).

    Parameters
    ----------
    porosity : float or array-like
        Immobile-domain porosity (required).
    volfrac : float or array-like
        Volume fraction of the immobile domain (required).
    zetaim : float or array-like
        Mass-transfer rate coefficient between mobile and immobile domains (required).
    name : str, default "ist"
        Package name.
    **options
        Extra ``flopy.mf6.ModflowGwtist`` options (decay/sorption arrays, filerecords).
    """

    return PackageSpec(
        name,
        flopy.mf6.ModflowGwtist,
        {"porosity": porosity, "volfrac": volfrac, "zetaim": zetaim, **options},
    )


def ssm(*, sources: Any, name: str = "ssm", **options: Any) -> PackageSpec:
    """Source-sink mixing (SSM) package for a GWT model.

    Injects solute mass through the flow model's source/sink boundaries (RCH, WEL,
    GHB, ...), so it is effectively required whenever the paired GWF model has any.

    Parameters
    ----------
    sources : list
        Rows of ``[pname, srctype, auxname]`` naming each flow package and how its
        concentration is supplied, e.g. ``[["chd", "AUX", "concentration"]]``.
    name : str, default "ssm"
        Package name.
    **options
        Extra ``flopy.mf6.ModflowGwtssm`` options (``fileinput``, ...).
    """

    return PackageSpec(name, flopy.mf6.ModflowGwtssm, {"sources": sources, **options})


def cnc(*, stress_period_data: Any, name: str = "cnc", **options: Any) -> PackageSpec:
    """Constant-concentration (CNC) list-BC package for a GWT model.

    Parameters
    ----------
    stress_period_data : dict
        FloPy mapping ``{period: [[cellid, conc], ...]}``.
    name : str, default "cnc"
        Package name.
    **options
        Extra ``flopy.mf6.ModflowGwtcnc`` options (``boundnames``, ``auxiliary``, ...).
    """

    return PackageSpec(
        name, flopy.mf6.ModflowGwtcnc, {"stress_period_data": stress_period_data, **options}
    )


def src(*, stress_period_data: Any, name: str = "src", **options: Any) -> PackageSpec:
    """Mass-source loading (SRC) list-BC package for a GWT model.

    Parameters
    ----------
    stress_period_data : dict
        FloPy mapping ``{period: [[cellid, smassrate], ...]}`` (solute mass loading rate).
    name : str, default "src"
        Package name.
    **options
        Extra ``flopy.mf6.ModflowGwtsrc`` options.
    """

    return PackageSpec(
        name, flopy.mf6.ModflowGwtsrc, {"stress_period_data": stress_period_data, **options}
    )


# --- GWE (energy transport) packages -------------------------------------------


def est(
    *,
    porosity: Any,
    heat_capacity_solid: Any,
    density_solid: Any,
    name: str = "est",
    **options: Any,
) -> PackageSpec:
    """Energy storage and transfer (EST) package for a GWE model.

    Parameters
    ----------
    porosity : float or array-like
        Aquifer porosity (required).
    heat_capacity_solid : float or array-like
        Heat capacity of the solid (required).
    density_solid : float or array-like
        Density of the solid (required).
    name : str, default "est"
        Package name.
    **options
        Extra ``flopy.mf6.ModflowGweest`` options; the water-side properties
        (``density_water``, ``heat_capacity_water``, ``latent_heat_vaporization``)
        carry sensible FloPy defaults.
    """

    return PackageSpec(
        name,
        flopy.mf6.ModflowGweest,
        {
            "porosity": porosity,
            "heat_capacity_solid": heat_capacity_solid,
            "density_solid": density_solid,
            **options,
        },
    )


def cnd(
    *,
    ktw: Any | None = None,
    kts: Any | None = None,
    alh: Any | None = None,
    ath1: Any | None = None,
    name: str = "cnd",
    **options: Any,
) -> PackageSpec:
    """Conduction and dispersion (CND) package for a GWE model.

    Parameters
    ----------
    ktw : float or array-like, optional
        Thermal conductivity of water.
    kts : float or array-like, optional
        Thermal conductivity of the solid.
    alh : float or array-like, optional
        Longitudinal dispersivity in the horizontal direction.
    ath1 : float or array-like, optional
        Transverse dispersivity in the first horizontal direction.
    name : str, default "cnd"
        Package name.
    **options
        Extra ``flopy.mf6.ModflowGwecnd`` options.
    """

    values = dict(options)
    for key, value in {"ktw": ktw, "kts": kts, "alh": alh, "ath1": ath1}.items():
        if value is not None:
            values[key] = value
    return PackageSpec(name, flopy.mf6.ModflowGwecnd, values)


def ctp(*, stress_period_data: Any, name: str = "ctp", **options: Any) -> PackageSpec:
    """Constant-temperature (CTP) list-BC package for a GWE model.

    Parameters
    ----------
    stress_period_data : dict
        FloPy mapping ``{period: [[cellid, temp], ...]}``.
    name : str, default "ctp"
        Package name.
    **options
        Extra ``flopy.mf6.ModflowGwectp`` options.
    """

    return PackageSpec(
        name, flopy.mf6.ModflowGwectp, {"stress_period_data": stress_period_data, **options}
    )


def esl(*, stress_period_data: Any, name: str = "esl", **options: Any) -> PackageSpec:
    """Energy source loading (ESL) list-BC package for a GWE model.

    Parameters
    ----------
    stress_period_data : dict
        FloPy mapping ``{period: [[cellid, senerrate], ...]}`` (energy source loading rate).
    name : str, default "esl"
        Package name.
    **options
        Extra ``flopy.mf6.ModflowGweesl`` options.
    """

    return PackageSpec(
        name, flopy.mf6.ModflowGweesl, {"stress_period_data": stress_period_data, **options}
    )


def mip(
    *,
    porosity: Any,
    retfactor: Any | None = None,
    izone: Any | None = None,
    name: str = "mip",
    **options: Any,
) -> PackageSpec:
    """Model input (MIP) package for a PRT model: tracking porosity.

    The one physical input particle tracking needs -- advective velocity is the
    cell flow divided by ``porosity``. Attach it to a :func:`prt` model alongside
    the grid (``mf.disv``/``mf.dis`` dispatch on the PRT kind), :func:`prp`, and
    ``mf.oc``.

    Parameters
    ----------
    porosity : float or array-like
        Effective (transport) porosity -- a scalar or per-cell array.
    retfactor : float or array-like, optional
        Retardation factor (>= 1 slows particles; default no retardation).
    izone : int or array-like, optional
        Zone numbers reported in the track output's ``izone`` column.
    name : str, default "mip"
        Package name.
    **options
        Extra ``flopy.mf6.ModflowPrtmip`` options.

    Returns
    -------
    PackageSpec
    """

    values = {"porosity": porosity, **options}
    if retfactor is not None:
        values["retfactor"] = retfactor
    if izone is not None:
        values["izone"] = izone
    return PackageSpec(name, flopy.mf6.ModflowPrtmip, values)


def prp(
    *,
    packagedata: Any,
    nreleasepts: int | None = None,
    perioddata: Any | None = None,
    boundnames: bool | None = None,
    name: str = "prp",
    **options: Any,
) -> PackageSpec:
    """Particle release point (PRP) package for a PRT model.

    Declares where and when particles are released. ``packagedata`` rows are
    ``(irpt, (layer, cell), x, y, z)`` -- zero-based ``irpt`` and cellid -- with
    an optional trailing ``boundname`` string to group release points (the name
    is echoed -- uppercased -- into the track CSV's ``name`` column, which is what
    ``results.capture.map()`` groups by). ``nreleasepts`` defaults to
    ``len(packagedata)``; ``perioddata`` defaults to ``{0: ["FIRST"]}`` (release
    at the first step, matching :class:`~myflopy.modflow.mf6.prt.PRTProject`);
    ``boundnames`` is auto-enabled when rows carry the sixth (name) field.

    Parameters
    ----------
    packagedata : list
        Release-point rows ``[(irpt, (layer, cell), x, y, z[, boundname]), ...]``.
    nreleasepts : int, optional
        Number of release points; defaults to ``len(packagedata)``.
    perioddata : dict, optional
        ``{kper: settings}`` release schedule (``ALL``/``FIRST``/``LAST``/
        ``FREQUENCY <n>``/``STEPS <s...>``/``FRACTION <f>``); default
        ``{0: ["FIRST"]}``.
    boundnames : bool, optional
        Enable named release points; auto-detected from the row width when
        omitted.
    name : str, default "prp"
        Package name.
    **options
        Extra ``flopy.mf6.ModflowPrtprp`` options (``stoptime``, ``drape``,
        ``local_z``, ``extend_tracking``, ``stop_at_weak_sink``, ...).

    Returns
    -------
    PackageSpec
    """

    rows = list(packagedata)
    if boundnames is None:
        boundnames = any(len(tuple(row)) >= 6 for row in rows)
    values = {
        "packagedata": rows,
        "nreleasepts": len(rows) if nreleasepts is None else nreleasepts,
        "perioddata": {0: ["FIRST"]} if perioddata is None else perioddata,
        # explicit pname: flopy otherwise numbers PRP instances ("prp_0"), which
        # breaks get_package("prp") lookups.
        "pname": name,
        **options,
    }
    if boundnames:
        values["boundnames"] = True
    return PackageSpec(name, flopy.mf6.ModflowPrtprp, values)


class _CHDPackage:
    """Package-first CHD (constant/specified-head boundary) helpers.

    A CHD pins listed cells to a specified head -- used for fixed boundaries (a
    lake edge held at stage, a regional head condition). Three entry points:

    - ``mf.chd(stress_period_data=...)`` -- direct MF6 records;
    - ``mf.chd.gpkg(path, context=, nper=)`` -- build records from GeoPackage
      features mapped onto cells (``edges_only=True`` keeps just the active
      domain's edge cells);
    - ``mf.chd.flopy(...)`` -- the FloPy-native form.

    Examples
    --------
    >>> mf.chd(stress_period_data={0: [[(0, 0), 10.0]]})   # (cellid, head)
    >>> mf.chd.gpkg("bcs.gpkg", layer="fixed_head", context=ctx, nper=1)
    """

    def __call__(
        self,
        *,
        stress_period_data: Any,
        name: str = "chd",
        boundnames: bool = False,
        **options: Any,
    ) -> PackageSpec:
        """Return a CHD package spec from direct MF6 stress-period data.

        Parameters
        ----------
        stress_period_data : dict
            FloPy mapping ``{period: [[cellid, head], ...]}`` where ``cellid`` is
            ``(layer, cell)`` on a DISV grid.
        name : str, default "chd"
            Package name (also its slot in the model).
        boundnames : bool, default False
            Enable named boundaries in the package.
        **options
            Extra ``flopy.mf6.ModflowGwfchd`` options.

        Returns
        -------
        PackageSpec

        Examples
        --------
        >>> mf.chd(stress_period_data={0: [[(0, 0), 10.0]]})
        """

        return chd_spec(stress_period_data, name=name, boundnames=boundnames, **options)

    def flopy(
        self,
        *,
        stress_period_data: Any,
        name: str = "chd",
        boundnames: bool = False,
        **options: Any,
    ) -> PackageSpec:
        """Return a direct (list-based) FloPy-style CHD spec from stress-period data.

        Parameters
        ----------
        stress_period_data : dict
            FloPy mapping ``{period: [[cellid, head], ...]}``.
        name : str, default "chd"
            Package name.
        boundnames : bool, default False
            Enable named boundaries.
        **options
            Extra ``flopy.mf6.ModflowGwfchd`` options.

        Returns
        -------
        PackageSpec

        Examples
        --------
        >>> mf.chd.flopy(stress_period_data={0: [[(0, 0), 10.0]]})
        """

        return chd_spec(stress_period_data, name=name, boundnames=boundnames, **options)

    def gpkg(
        self,
        path: PathLike,
        *,
        context: ModelContext | None = None,
        nper: int,
        head: RowValue = "head",
        layer: str | None = None,
        name_field: str | None = "name",
        layer_field: str | None = "layer",
        period_field: str | None = None,
        layer_base: int = 1,
        period_base: int = 0,
        gaps: str = "warn",
        snap_within: float | None = None,
        min_above_bottom: float | None = None,
        name: str = "chd",
        edges_only: bool = False,
        boundnames: bool = True,
        **options: Any,
    ) -> PackageSpec:
        """Return a CHD package spec built from GeoPackage features mapped onto the grid.

        Parameters
        ----------
        path : Path or str
            GeoPackage/shapefile of boundary features (auto-reprojected to the grid CRS).
        context : ModelContext, optional
            The grid/domain the features are mapped onto. Given, the records resolve
            NOW, against ``context.grid`` (``ValueError`` when it carries none) and
            ``context.domain``, and ``with_grid`` re-declares them from the recorded
            recipe. Default ``None`` makes the spec LAZY: it holds the recipe and
            resolves at build against the model's own grid, which the model's
            context must then carry; the model's discretization must be declared
            before this package.
        nper : int
            Number of stress periods the boundary spans.
        head : str, float, list of str, or CellSurfaceOffset, default "head"
            The specified head. A value spec: an attribute column name, a numeric constant, a
            *list* of column names (indexed by stress period), or a
            :class:`~myflopy.geopackage.CellSurfaceOffset`.
            Two CHD records in one cell is an MF6 error, not a sum -- if features
            share a cell after mapping, the run fails at read time.
        layer : str, optional
            GeoPackage layer/table name (defaults to the first/only layer).
        name_field : str, optional, default "name"
            Attribute column used for boundnames.
        layer_field : str, optional, default "layer"
            Attribute column giving each feature's model layer (see ``layer_base``).
        period_field : str, optional
            Attribute column giving a feature's stress period (``None`` = all periods).
        layer_base : int, default 1
            Index base of ``layer_field`` values (1 = one-based GIS convention).
        period_base : int, default 0
            Index base of ``period_field`` values (0 = zero-based).
        gaps : {"warn", "nearest", "raise", "ignore"}, optional
            What a feature that lands on no active cell means. ``"warn"`` (the
            default) leaves it out and says so -- a ``UserWarning`` as well as the
            log, naming each feature -- in the same message that reports features
            clipped at the grid's edge and elevations below their cell bottom.
            ``"nearest"`` moves a missed POINT onto the nearest usable cell when it
            is within ``snap_within`` and reports the move; a line or polygon is
            never moved, and is reported as under ``"warn"``. ``"raise"`` raises
            ``ValueError`` wherever ``"warn"`` would warn. ``"ignore"`` says
            nothing. In every mode a source with NO feature landing raises
            ``ValueError`` naming both extents and CRSs; an unknown word raises
            ``ValueError`` at the call.
        snap_within : float, optional
            ``gaps="nearest"`` only (ignored otherwise): how far a missed point may
            be moved, in the grid's coordinate units, measured to the cell POLYGON
            rather than its centroid. Must be positive (``ValueError`` at the
            call). Default ``None`` is one cell width -- ``sqrt(area)`` of the cell
            it would move to -- resolved when the records resolve, so on the lazy
            form it follows the mesh across ``with_grid``.
        min_above_bottom : float, optional
            Keep every resolved elevation at least this far above ITS OWN cell bottom,
            per cell. ``None`` (the default) changes nothing but still REPORTS a record
            that falls below -- see ledger 290 for why the detection is default-on and the
            correction is not. ``0.0`` means as close to the bottom as the written file
            can express: FloPy writes ``botm`` fixed-point and a boundary value in E
            format, so a value set exactly equal to its cell bottom can round below it and
            MF6's strict comparison trips (measured, 13 of 24 real runs failed), and a CHD
            set exactly there aborts the solver. On ``riv`` it floors ``rbot`` and raises
            ``stage`` by the same delta, so the declared river depth survives.
        name : str, default "chd"
            Package name.
        edges_only : bool, default False
            Keep only cells on the edge of the ACTIVE domain -- bordering the
            outside of the mesh, or an inactive cell -- read from ``context``
            per layer (a domain-edge head condition).
        boundnames : bool, default True
            Enable named boundaries.
        **options
            Extra ``flopy.mf6.ModflowGwfchd`` options.

        Returns
        -------
        PackageSpec

        Examples
        --------
        >>> mf.chd.gpkg("bcs.gpkg", layer="fixed_head", context=ctx, nper=1)
        """

        if context is None:
            return _lazy_bc_spec(
                "chd",
                path,
                name=name,
                nper=nper,
                layer=layer,
                name_field=name_field,
                layer_field=layer_field,
                period_field=period_field,
                layer_base=layer_base,
                period_base=period_base,
                gaps=gaps,
                snap_within=snap_within,
                min_above_bottom=min_above_bottom,
                values={"head": head},
                edges_only=edges_only,
                boundnames=boundnames,
                options=options,
            )

        return _source(
            path,
            context=context,
            nper=nper,
            layer=layer,
            name_field=name_field,
            layer_field=layer_field,
            period_field=period_field,
            layer_base=layer_base,
            period_base=period_base,
            gaps=gaps,
            snap_within=snap_within,
            min_above_bottom=min_above_bottom,
        ).chd(
            head=head,
            name=name,
            edges_only=edges_only,
            boundnames=boundnames,
            **options,
        )


class _GHBPackage:
    """Package-first GHB (general-head boundary) helpers.

    A GHB is a head-dependent flux that connects model cells to a fixed external
    head through a conductance -- the workhorse for regional underflow / far-field
    boundaries. Three entry points:

    - ``mf.ghb(stress_period_data=...)`` -- direct MF6 records;
    - ``mf.ghb.gpkg(path, context=, nper=)`` -- build records from GeoPackage
      features (auto-reprojected and mapped onto cells);
    - ``mf.ghb.flopy(...)`` -- the FloPy-native form.

    Examples
    --------
    >>> mf.ghb(stress_period_data={0: [[(0, 5), 86.0, 50.0]]})        # (cellid, head, cond)
    >>> mf.ghb.gpkg("bcs.gpkg", layer="underflow", context=ctx, nper=1)
    """

    def __call__(
        self,
        *,
        stress_period_data: Any,
        name: str = "ghb",
        auxiliary: Any = None,
        boundnames: bool = False,
        **options: Any,
    ) -> PackageSpec:
        """Return a GHB package spec from direct MF6 stress-period data.

        Parameters
        ----------
        stress_period_data : dict
            FloPy mapping ``{period: [[cellid, bhead, cond], ...]}`` where ``cellid``
            is ``(layer, cell)`` on a DISV grid, ``bhead`` the boundary head, and
            ``cond`` the conductance.
        name : str, default "ghb"
            Package name.
        auxiliary : optional
            Auxiliary variable name(s) forwarded to FloPy.
        boundnames : bool, default False
            Enable named boundaries.
        **options
            Extra ``flopy.mf6.ModflowGwfghb`` options.

        Returns
        -------
        PackageSpec

        Examples
        --------
        >>> mf.ghb(stress_period_data={0: [[(0, 5), 86.0, 50.0]]})
        """

        return ghb_spec(
            stress_period_data,
            name=name,
            auxiliary=auxiliary,
            boundnames=boundnames,
            **options,
        )

    def flopy(
        self,
        *,
        stress_period_data: Any,
        name: str = "ghb",
        auxiliary: Any = None,
        boundnames: bool = False,
        **options: Any,
    ) -> PackageSpec:
        """Return a direct (list-based) FloPy-style GHB spec from stress-period data.

        Parameters
        ----------
        stress_period_data : dict
            FloPy mapping ``{period: [[cellid, bhead, cond], ...]}``.
        name : str, default "ghb"
            Package name.
        auxiliary : optional
            Auxiliary variable name(s) forwarded to FloPy.
        boundnames : bool, default False
            Enable named boundaries.
        **options
            Extra ``flopy.mf6.ModflowGwfghb`` options.

        Returns
        -------
        PackageSpec

        Examples
        --------
        >>> mf.ghb.flopy(stress_period_data={0: [[(0, 5), 86.0, 50.0]]})
        """

        return ghb_spec(
            stress_period_data,
            name=name,
            auxiliary=auxiliary,
            boundnames=boundnames,
            **options,
        )

    def gpkg(
        self,
        path: PathLike,
        *,
        context: ModelContext | None = None,
        nper: int,
        head: RowValue = "head",
        conductance: RowValue = "conductance",
        layer: str | None = None,
        name_field: str | None = "name",
        layer_field: str | None = "layer",
        period_field: str | None = None,
        layer_base: int = 1,
        period_base: int = 0,
        gaps: str = "warn",
        snap_within: float | None = None,
        min_above_bottom: float | None = None,
        name: str = "ghb",
        edges_only: bool = False,
        boundnames: bool = True,
        **options: Any,
    ) -> PackageSpec:
        """Return a GHB package spec built from GeoPackage features mapped onto the grid.

        Parameters
        ----------
        path : Path or str
            GeoPackage/shapefile of boundary features (auto-reprojected to the grid CRS).
        context : ModelContext, optional
            The grid/domain the features are mapped onto. Given, the records resolve
            NOW, against ``context.grid`` (``ValueError`` when it carries none) and
            ``context.domain``, and ``with_grid`` re-declares them from the recorded
            recipe. Default ``None`` makes the spec LAZY: it holds the recipe and
            resolves at build against the model's own grid, which the model's
            context must then carry; the model's discretization must be declared
            before this package.
        nper : int
            Number of stress periods the boundary spans.
        head : str, float, list of str, or CellSurfaceOffset, default "head"
            The external head the cells are connected to. A value spec: an attribute column
            name, a numeric constant, a *list* of column names (indexed by stress period), or a
            :class:`~myflopy.geopackage.CellSurfaceOffset`. MF6 rejects a head below its cell
            bottom.
        conductance : str, float, list of str, or Spread, default "conductance"
            Conductance of the connection, L^2/T. Wrap it in :class:`~myflopy.geopackage.Spread`
            when the feature spans more than one cell -- a bare column copies the SAME value
            into every cell the feature covers, which is right for an elevation and multiplies a
            conductance by the cell count.
        layer : str, optional
            Name of the layer/table inside the file. Defaults to the first, which is
            the only one in a single-layer GeoPackage or a shapefile.
        name_field : str, optional, default "name"
            Attribute column used for boundnames. Set ``boundnames=False`` if the
            file has no such column, or the build raises for the missing field.
        layer_field : str, optional, default "layer"
            Attribute column giving each feature's model layer. Cells inactive in that layer are
            skipped, as are vertical-passthrough cells (``idomain < 0``). ``None`` puts every
            feature in layer 1.
        period_field : str, optional
            Attribute column giving the stress period a feature belongs to. ``None``
            (the default) applies every feature to all ``nper`` periods. Use it when
            the file holds one row per feature per period; for values that change in
            time on a fixed set of features, pass a *list* of column names as the
            value instead -- it is indexed by period.
        layer_base : int, default 1
            Index base of ``layer_field``. ``1`` matches the usual GIS convention of
            numbering layers from one; pass ``0`` for a zero-based column.
        period_base : int, default 0
            Index base of ``period_field``. ``0`` matches MODFLOW's zero-based
            periods as myflopy counts them.
        gaps : {"warn", "nearest", "raise", "ignore"}, optional
            What a feature that lands on no active cell means. ``"warn"`` (the
            default) leaves it out and says so -- a ``UserWarning`` as well as the
            log, naming each feature -- in the same message that reports features
            clipped at the grid's edge and elevations below their cell bottom.
            ``"nearest"`` moves a missed POINT onto the nearest usable cell when it
            is within ``snap_within`` and reports the move; a line or polygon is
            never moved, and is reported as under ``"warn"``. ``"raise"`` raises
            ``ValueError`` wherever ``"warn"`` would warn. ``"ignore"`` says
            nothing. In every mode a source with NO feature landing raises
            ``ValueError`` naming both extents and CRSs; an unknown word raises
            ``ValueError`` at the call.
        snap_within : float, optional
            ``gaps="nearest"`` only (ignored otherwise): how far a missed point may
            be moved, in the grid's coordinate units, measured to the cell POLYGON
            rather than its centroid. Must be positive (``ValueError`` at the
            call). Default ``None`` is one cell width -- ``sqrt(area)`` of the cell
            it would move to -- resolved when the records resolve, so on the lazy
            form it follows the mesh across ``with_grid``.
        min_above_bottom : float, optional
            Keep every resolved elevation at least this far above ITS OWN cell bottom,
            per cell. ``None`` (the default) changes nothing but still REPORTS a record
            that falls below -- see ledger 290 for why the detection is default-on and the
            correction is not. ``0.0`` means as close to the bottom as the written file
            can express: FloPy writes ``botm`` fixed-point and a boundary value in E
            format, so a value set exactly equal to its cell bottom can round below it and
            MF6's strict comparison trips (measured, 13 of 24 real runs failed), and a CHD
            set exactly there aborts the solver. On ``riv`` it floors ``rbot`` and raises
            ``stage`` by the same delta, so the declared river depth survives.
        name : str, default "ghb"
            Package name.
        edges_only : bool, default False
            Keep only cells on the edge of the ACTIVE domain, per layer (a
            far-field underflow condition).
        boundnames : bool, default True
            Enable named boundaries.
        **options
            Extra ``flopy.mf6.ModflowGwfghb`` options.

        Returns
        -------
        PackageSpec

        Examples
        --------
        >>> mf.ghb.gpkg("bcs.gpkg", layer="underflow", context=ctx, nper=1)
        """

        if context is None:
            return _lazy_bc_spec(
                "ghb",
                path,
                name=name,
                nper=nper,
                layer=layer,
                name_field=name_field,
                layer_field=layer_field,
                period_field=period_field,
                layer_base=layer_base,
                period_base=period_base,
                gaps=gaps,
                snap_within=snap_within,
                min_above_bottom=min_above_bottom,
                values={"head": head, "conductance": conductance},
                edges_only=edges_only,
                boundnames=boundnames,
                options=options,
            )

        return _source(
            path,
            context=context,
            nper=nper,
            layer=layer,
            name_field=name_field,
            layer_field=layer_field,
            period_field=period_field,
            layer_base=layer_base,
            period_base=period_base,
            gaps=gaps,
            snap_within=snap_within,
            min_above_bottom=min_above_bottom,
        ).ghb(
            head=head,
            conductance=conductance,
            name=name,
            edges_only=edges_only,
            boundnames=boundnames,
            **options,
        )


class _DRNPackage:
    """Package-first DRN (drain) helpers.

    A drain removes water from a cell only when head rises above the drain
    elevation (one-way, head-dependent) -- used for toe-of-slope springs, tile
    drains, and seepage faces. Three entry points:

    - ``mf.drn(stress_period_data=...)`` -- direct MF6 records;
    - ``mf.drn.gpkg(path, context=, nper=)`` -- build records from GeoPackage
      features mapped onto cells;
    - ``mf.drn.flopy(...)`` -- the FloPy-native form.

    Examples
    --------
    >>> mf.drn(stress_period_data={0: [[(0, 12), 95.0, 30.0]]})   # (cellid, elev, cond)
    >>> mf.drn.gpkg("bcs.gpkg", layer="springs", context=ctx, nper=1)
    """

    def __call__(
        self,
        *,
        stress_period_data: Any,
        name: str = "drn",
        boundnames: bool = False,
        **options: Any,
    ) -> PackageSpec:
        """Return a DRN package spec from direct MF6 stress-period data.

        Parameters
        ----------
        stress_period_data : dict
            FloPy mapping ``{period: [[cellid, elev, cond], ...]}`` -- drain
            elevation and conductance per ``(layer, cell)``.
        name : str, default "drn"
            Package name.
        boundnames : bool, default False
            Enable named boundaries.
        **options
            Extra ``flopy.mf6.ModflowGwfdrn`` options.

        Returns
        -------
        PackageSpec

        Examples
        --------
        >>> mf.drn(stress_period_data={0: [[(0, 12), 95.0, 30.0]]})
        """

        return drn_spec(stress_period_data, name=name, boundnames=boundnames, **options)

    def flopy(
        self,
        *,
        stress_period_data: Any,
        name: str = "drn",
        boundnames: bool = False,
        **options: Any,
    ) -> PackageSpec:
        """Return a direct (list-based) FloPy-style DRN spec from stress-period data.

        Parameters
        ----------
        stress_period_data : dict
            FloPy mapping ``{period: [[cellid, elev, cond], ...]}``.
        name : str, default "drn"
            Package name.
        boundnames : bool, default False
            Enable named boundaries.
        **options
            Extra ``flopy.mf6.ModflowGwfdrn`` options.

        Returns
        -------
        PackageSpec

        Examples
        --------
        >>> mf.drn.flopy(stress_period_data={0: [[(0, 12), 95.0, 30.0]]})
        """

        return drn_spec(stress_period_data, name=name, boundnames=boundnames, **options)

    def gpkg(
        self,
        path: PathLike,
        *,
        context: ModelContext | None = None,
        nper: int,
        elevation: RowValue = "elevation",
        conductance: RowValue = "conductance",
        layer: str | None = None,
        name_field: str | None = "name",
        layer_field: str | None = "layer",
        period_field: str | None = None,
        layer_base: int = 1,
        period_base: int = 0,
        gaps: str = "warn",
        snap_within: float | None = None,
        min_above_bottom: float | None = None,
        name: str = "drn",
        edges_only: bool = False,
        boundnames: bool = True,
        **options: Any,
    ) -> PackageSpec:
        """Return a DRN package spec built from GeoPackage features mapped onto the grid.

        Parameters
        ----------
        path : Path or str
            GeoPackage/shapefile of drain features (auto-reprojected to the grid CRS).
        context : ModelContext, optional
            The grid/domain the features are mapped onto. Given, the records resolve
            NOW, against ``context.grid`` (``ValueError`` when it carries none) and
            ``context.domain``, and ``with_grid`` re-declares them from the recorded
            recipe. Default ``None`` makes the spec LAZY: it holds the recipe and
            resolves at build against the model's own grid, which the model's
            context must then carry; the model's discretization must be declared
            before this package.
        nper : int
            Number of stress periods.
        elevation : str, float, list of str, or CellSurfaceOffset, default "elevation"
            The drain invert. A value spec: an attribute column name, a numeric constant, a
            *list* of column names (indexed by stress period), or a
            :class:`~myflopy.geopackage.CellSurfaceOffset`.
            ``CellSurfaceOffset("cell_bottom", offset=2.0)`` puts it 2 ft above each
            cell's own floor; ``("model_top", offset=-2.0)`` puts it 2 ft below ground.
            **MF6 rejects an invert below its cell bottom only at ``icelltype != 0``**;
            at ``icelltype = 0`` -- what ``mf.npf(k=...)`` writes, naming no icelltype at
            all -- it ACCEPTS one and over-drains under Normal termination (measured 2.17x
            on a bench, 1.09x on a real model). Either way one is reported now, and
            ``min_above_bottom=`` lifts it (ledger 290).
        conductance : str, float, list of str, or Spread, default "conductance"
            Drain conductance, L^2/T. Wrap it in :class:`~myflopy.geopackage.Spread` when the
            feature spans more than one cell -- a bare column copies the SAME value into every
            cell the feature covers, which is right for an elevation and multiplies a
            conductance by the cell count.
        layer : str, optional
            Name of the layer/table inside the file. Defaults to the first, which is
            the only one in a single-layer GeoPackage or a shapefile.
        name_field : str, optional, default "name"
            Attribute column used for boundnames. Set ``boundnames=False`` if the
            file has no such column, or the build raises for the missing field.
        layer_field : str, optional, default "layer"
            Attribute column giving each feature's model layer. Cells inactive in that layer are
            skipped, as are vertical-passthrough cells (``idomain < 0``). ``None`` puts every
            feature in layer 1.
        period_field : str, optional
            Attribute column giving the stress period a feature belongs to. ``None``
            (the default) applies every feature to all ``nper`` periods. Use it when
            the file holds one row per feature per period; for values that change in
            time on a fixed set of features, pass a *list* of column names as the
            value instead -- it is indexed by period.
        layer_base : int, default 1
            Index base of ``layer_field``. ``1`` matches the usual GIS convention of
            numbering layers from one; pass ``0`` for a zero-based column.
        period_base : int, default 0
            Index base of ``period_field``. ``0`` matches MODFLOW's zero-based
            periods as myflopy counts them.
        gaps : {"warn", "nearest", "raise", "ignore"}, optional
            What a feature that lands on no active cell means. ``"warn"`` (the
            default) leaves it out and says so -- a ``UserWarning`` as well as the
            log, naming each feature -- in the same message that reports features
            clipped at the grid's edge and elevations below their cell bottom.
            ``"nearest"`` moves a missed POINT onto the nearest usable cell when it
            is within ``snap_within`` and reports the move; a line or polygon is
            never moved, and is reported as under ``"warn"``. ``"raise"`` raises
            ``ValueError`` wherever ``"warn"`` would warn. ``"ignore"`` says
            nothing. In every mode a source with NO feature landing raises
            ``ValueError`` naming both extents and CRSs; an unknown word raises
            ``ValueError`` at the call.
        snap_within : float, optional
            ``gaps="nearest"`` only (ignored otherwise): how far a missed point may
            be moved, in the grid's coordinate units, measured to the cell POLYGON
            rather than its centroid. Must be positive (``ValueError`` at the
            call). Default ``None`` is one cell width -- ``sqrt(area)`` of the cell
            it would move to -- resolved when the records resolve, so on the lazy
            form it follows the mesh across ``with_grid``.
        min_above_bottom : float, optional
            Keep every resolved elevation at least this far above ITS OWN cell bottom,
            per cell. ``None`` (the default) changes nothing but still REPORTS a record
            that falls below -- see ledger 290 for why the detection is default-on and the
            correction is not. ``0.0`` means as close to the bottom as the written file
            can express: FloPy writes ``botm`` fixed-point and a boundary value in E
            format, so a value set exactly equal to its cell bottom can round below it and
            MF6's strict comparison trips (measured, 13 of 24 real runs failed), and a CHD
            set exactly there aborts the solver. On ``riv`` it floors ``rbot`` and raises
            ``stage`` by the same delta, so the declared river depth survives.
        name : str, default "drn"
            Package name.
        edges_only : bool, default False
            Keep only cells on the edge of the ACTIVE domain, per layer.
        boundnames : bool, default True
            Enable named boundaries.
        **options
            Extra ``flopy.mf6.ModflowGwfdrn`` options.

        Returns
        -------
        PackageSpec

        Examples
        --------
        >>> mf.drn.gpkg("bcs.gpkg", layer="springs", context=ctx, nper=1)
        """

        if context is None:
            return _lazy_bc_spec(
                "drn",
                path,
                name=name,
                nper=nper,
                layer=layer,
                name_field=name_field,
                layer_field=layer_field,
                period_field=period_field,
                layer_base=layer_base,
                period_base=period_base,
                gaps=gaps,
                snap_within=snap_within,
                min_above_bottom=min_above_bottom,
                values={"elevation": elevation, "conductance": conductance},
                edges_only=edges_only,
                boundnames=boundnames,
                options=options,
            )

        return _source(
            path,
            context=context,
            nper=nper,
            layer=layer,
            name_field=name_field,
            layer_field=layer_field,
            period_field=period_field,
            layer_base=layer_base,
            period_base=period_base,
            gaps=gaps,
            snap_within=snap_within,
            min_above_bottom=min_above_bottom,
        ).drn(
            elevation=elevation,
            conductance=conductance,
            name=name,
            edges_only=edges_only,
            boundnames=boundnames,
            **options,
        )


class _RIVPackage:
    """Package-first RIV (river) helpers.

    A river exchanges water with the aquifer through a streambed conductance:
    head-dependent in both directions while the water table is above the river
    bottom (``rbot``), and a fixed stage-to-rbot leakage once it drops below --
    the standard boundary for streams too small to warrant SFR routing. Three
    entry points:

    - ``mf.riv(stress_period_data=...)`` -- direct MF6 records;
    - ``mf.riv.gpkg(path, context=, nper=)`` -- build records from GeoPackage
      features mapped onto cells;
    - ``mf.riv.flopy(...)`` -- the FloPy-native form.

    Examples
    --------
    >>> mf.riv(stress_period_data={0: [[(0, 7), 98.0, 40.0, 96.5]]})   # (cellid, stage, cond, rbot)
    >>> mf.riv.gpkg("bcs.gpkg", layer="river", context=ctx, nper=1)
    """

    def __call__(
        self,
        *,
        stress_period_data: Any,
        name: str = "riv",
        auxiliary: Any = None,
        boundnames: bool = False,
        **options: Any,
    ) -> PackageSpec:
        """Return a RIV package spec from direct MF6 stress-period data.

        Parameters
        ----------
        stress_period_data : dict
            FloPy mapping ``{period: [[cellid, stage, cond, rbot], ...]}`` --
            river stage, riverbed conductance, and river-bottom elevation per
            ``(layer, cell)`` on a DISV grid.
        name : str, default "riv"
            Package name.
        auxiliary : optional
            Auxiliary variable name(s) forwarded to FloPy.
        boundnames : bool, default False
            Enable named boundaries.
        **options
            Extra ``flopy.mf6.ModflowGwfriv`` options.

        Returns
        -------
        PackageSpec

        Examples
        --------
        >>> mf.riv(stress_period_data={0: [[(0, 7), 98.0, 40.0, 96.5]]})
        """

        return riv_spec(
            stress_period_data,
            name=name,
            auxiliary=auxiliary,
            boundnames=boundnames,
            **options,
        )

    def flopy(
        self,
        *,
        stress_period_data: Any,
        name: str = "riv",
        auxiliary: Any = None,
        boundnames: bool = False,
        **options: Any,
    ) -> PackageSpec:
        """Return a direct (list-based) FloPy-style RIV spec from stress-period data.

        Parameters
        ----------
        stress_period_data : dict
            FloPy mapping ``{period: [[cellid, stage, cond, rbot], ...]}``.
        name : str, default "riv"
            Package name.
        auxiliary : optional
            Auxiliary variable name(s) forwarded to FloPy.
        boundnames : bool, default False
            Enable named boundaries.
        **options
            Extra ``flopy.mf6.ModflowGwfriv`` options.

        Returns
        -------
        PackageSpec

        Examples
        --------
        >>> mf.riv.flopy(stress_period_data={0: [[(0, 7), 98.0, 40.0, 96.5]]})
        """

        return riv_spec(
            stress_period_data,
            name=name,
            auxiliary=auxiliary,
            boundnames=boundnames,
            **options,
        )

    def gpkg(
        self,
        path: PathLike,
        *,
        context: ModelContext | None = None,
        nper: int,
        stage: RowValue = "stage",
        conductance: RowValue = "conductance",
        rbot: RowValue = "rbot",
        layer: str | None = None,
        name_field: str | None = "name",
        layer_field: str | None = "layer",
        period_field: str | None = None,
        layer_base: int = 1,
        period_base: int = 0,
        gaps: str = "warn",
        snap_within: float | None = None,
        min_above_bottom: float | None = None,
        name: str = "riv",
        edges_only: bool = False,
        boundnames: bool = True,
        **options: Any,
    ) -> PackageSpec:
        """Return a RIV package spec built from GeoPackage features mapped onto the grid.

        Parameters
        ----------
        path : Path or str
            GeoPackage/shapefile of river features (auto-reprojected to the grid CRS).
        context : ModelContext, optional
            The grid/domain the features are mapped onto. Given, the records resolve
            NOW, against ``context.grid`` (``ValueError`` when it carries none) and
            ``context.domain``, and ``with_grid`` re-declares them from the recorded
            recipe. Default ``None`` makes the spec LAZY: it holds the recipe and
            resolves at build against the model's own grid, which the model's
            context must then carry; the model's discretization must be declared
            before this package.
        nper : int
            Number of stress periods.
        stage : str, float, list of str, or CellSurfaceOffset, default "stage"
            Water-surface elevation in the river. A value spec: an attribute column name, a
            numeric constant, a *list* of column names (indexed by stress period), or a
            :class:`~myflopy.geopackage.CellSurfaceOffset`.
        conductance : str, float, list of str, or Spread, default "conductance"
            Riverbed conductance, L^2/T. Wrap it in :class:`~myflopy.geopackage.Spread` when the
            feature spans more than one cell -- a bare column copies the SAME value into every
            cell the feature covers, which is right for an elevation and multiplies a
            conductance by the cell count.
        rbot : str, float, list of str, or CellSurfaceOffset, default "rbot"
            Elevation of the riverbed bottom, which must sit below ``stage`` and
            above the cell bottom. A value spec: an attribute column name, a numeric constant, a
            *list* of column names (indexed by stress period), or a
            :class:`~myflopy.geopackage.CellSurfaceOffset`.
        layer : str, optional
            Name of the layer/table inside the file. Defaults to the first, which is
            the only one in a single-layer GeoPackage or a shapefile.
        name_field : str, optional, default "name"
            Attribute column used for boundnames. Set ``boundnames=False`` if the
            file has no such column, or the build raises for the missing field.
        layer_field : str, optional, default "layer"
            Attribute column giving each feature's model layer. Cells inactive in that layer are
            skipped, as are vertical-passthrough cells (``idomain < 0``). ``None`` puts every
            feature in layer 1.
        period_field : str, optional
            Attribute column giving the stress period a feature belongs to. ``None``
            (the default) applies every feature to all ``nper`` periods. Use it when
            the file holds one row per feature per period; for values that change in
            time on a fixed set of features, pass a *list* of column names as the
            value instead -- it is indexed by period.
        layer_base : int, default 1
            Index base of ``layer_field``. ``1`` matches the usual GIS convention of
            numbering layers from one; pass ``0`` for a zero-based column.
        period_base : int, default 0
            Index base of ``period_field``. ``0`` matches MODFLOW's zero-based
            periods as myflopy counts them.
        gaps : {"warn", "nearest", "raise", "ignore"}, optional
            What a feature that lands on no active cell means. ``"warn"`` (the
            default) leaves it out and says so -- a ``UserWarning`` as well as the
            log, naming each feature -- in the same message that reports features
            clipped at the grid's edge and elevations below their cell bottom.
            ``"nearest"`` moves a missed POINT onto the nearest usable cell when it
            is within ``snap_within`` and reports the move; a line or polygon is
            never moved, and is reported as under ``"warn"``. ``"raise"`` raises
            ``ValueError`` wherever ``"warn"`` would warn. ``"ignore"`` says
            nothing. In every mode a source with NO feature landing raises
            ``ValueError`` naming both extents and CRSs; an unknown word raises
            ``ValueError`` at the call.
        snap_within : float, optional
            ``gaps="nearest"`` only (ignored otherwise): how far a missed point may
            be moved, in the grid's coordinate units, measured to the cell POLYGON
            rather than its centroid. Must be positive (``ValueError`` at the
            call). Default ``None`` is one cell width -- ``sqrt(area)`` of the cell
            it would move to -- resolved when the records resolve, so on the lazy
            form it follows the mesh across ``with_grid``.
        min_above_bottom : float, optional
            Keep every resolved elevation at least this far above ITS OWN cell bottom,
            per cell. ``None`` (the default) changes nothing but still REPORTS a record
            that falls below -- see ledger 290 for why the detection is default-on and the
            correction is not. ``0.0`` means as close to the bottom as the written file
            can express: FloPy writes ``botm`` fixed-point and a boundary value in E
            format, so a value set exactly equal to its cell bottom can round below it and
            MF6's strict comparison trips (measured, 13 of 24 real runs failed), and a CHD
            set exactly there aborts the solver. On ``riv`` it floors ``rbot`` and raises
            ``stage`` by the same delta, so the declared river depth survives.
        name : str, default "riv"
            Package name.
        edges_only : bool, default False
            Keep only cells on the edge of the ACTIVE domain, per layer.
        boundnames : bool, default True
            Enable named boundaries.
        **options
            Extra ``flopy.mf6.ModflowGwfriv`` options.

        Returns
        -------
        PackageSpec

        Examples
        --------
        >>> mf.riv.gpkg("bcs.gpkg", layer="river", context=ctx, nper=1)
        """

        if context is None:
            return _lazy_bc_spec(
                "riv",
                path,
                name=name,
                nper=nper,
                layer=layer,
                name_field=name_field,
                layer_field=layer_field,
                period_field=period_field,
                layer_base=layer_base,
                period_base=period_base,
                gaps=gaps,
                snap_within=snap_within,
                min_above_bottom=min_above_bottom,
                values={"stage": stage, "conductance": conductance, "rbot": rbot},
                edges_only=edges_only,
                boundnames=boundnames,
                options=options,
            )

        return _source(
            path,
            context=context,
            nper=nper,
            layer=layer,
            name_field=name_field,
            layer_field=layer_field,
            period_field=period_field,
            layer_base=layer_base,
            period_base=period_base,
            gaps=gaps,
            snap_within=snap_within,
            min_above_bottom=min_above_bottom,
        ).riv(
            stage=stage,
            conductance=conductance,
            rbot=rbot,
            name=name,
            edges_only=edges_only,
            boundnames=boundnames,
            **options,
        )


class _WELPackage:
    """Package-first WEL (well) helpers.

    A WEL applies a specified volumetric flux to cells -- negative to pump out
    (extraction), positive to inject. Three entry points:

    - ``mf.wel(stress_period_data=...)`` -- direct MF6 records;
    - ``mf.wel.gpkg(path, context=, nper=)`` -- build records from GeoPackage
      point features mapped onto cells (with per-period rates);
    - ``mf.wel.flopy(...)`` -- the FloPy-native form.

    Examples
    --------
    >>> mf.wel(stress_period_data={0: [[(0, 42), -500.0]]})   # (cellid, rate); negative = pumping
    >>> mf.wel.gpkg("wells.gpkg", context=ctx, nper=12)
    """

    def __call__(
        self,
        *,
        stress_period_data: Any,
        name: str = "wel",
        auxiliary: Any = None,
        boundnames: bool = True,
        **options: Any,
    ) -> PackageSpec:
        """Return a WEL package spec from direct MF6 stress-period data.

        Parameters
        ----------
        stress_period_data : dict
            FloPy mapping ``{period: [[cellid, rate], ...]}`` -- a volumetric flux
            per ``(layer, cell)``; **negative to pump out**, positive to inject.
        name : str, default "wel"
            Package name.
        auxiliary : optional
            Auxiliary variable name(s) forwarded to FloPy.
        boundnames : bool, default True
            Enable named boundaries.
        **options
            Extra ``flopy.mf6.ModflowGwfwel`` options.

        Returns
        -------
        PackageSpec

        Examples
        --------
        >>> mf.wel(stress_period_data={0: [[(0, 42), -500.0]]})   # extraction
        """

        return wel_spec(
            stress_period_data,
            name=name,
            auxiliary=auxiliary,
            boundnames=boundnames,
            **options,
        )

    def flopy(
        self,
        *,
        stress_period_data: Any,
        name: str = "wel",
        auxiliary: Any = None,
        boundnames: bool = True,
        **options: Any,
    ) -> PackageSpec:
        """Return a direct (list-based) FloPy-style WEL spec from stress-period data.

        Parameters
        ----------
        stress_period_data : dict
            FloPy mapping ``{period: [[cellid, rate], ...]}`` -- **negative to
            pump out**, positive to inject.
        name : str, default "wel"
            Package name.
        auxiliary : optional
            Auxiliary variable name(s) forwarded to FloPy.
        boundnames : bool, default True
            Enable named boundaries.
        **options
            Extra ``flopy.mf6.ModflowGwfwel`` options.

        Returns
        -------
        PackageSpec

        Examples
        --------
        >>> mf.wel.flopy(stress_period_data={0: [[(0, 42), -500.0]]})
        """

        return wel_spec(
            stress_period_data,
            name=name,
            auxiliary=auxiliary,
            boundnames=boundnames,
            **options,
        )

    def gpkg(
        self,
        path: PathLike,
        *,
        context: ModelContext | None = None,
        nper: int,
        rate: RowValue = "rate",
        layer: str | None = None,
        name_field: str | None = "name",
        layer_field: str | None = "layer",
        period_field: str | None = None,
        layer_base: int = 1,
        period_base: int = 0,
        gaps: str = "warn",
        snap_within: float | None = None,
        name: str = "wel",
        boundnames: bool = True,
        **options: Any,
    ) -> PackageSpec:
        """Return a WEL package spec built from GeoPackage point features mapped onto the grid.

        Parameters
        ----------
        path : Path or str
            GeoPackage/shapefile of well point features (auto-reprojected to the grid CRS).
        context : ModelContext, optional
            The grid/domain the wells are mapped onto. Given, the records resolve
            NOW, against ``context.grid`` (``ValueError`` when it carries none) and
            ``context.domain``, and ``with_grid`` re-declares them from the recorded
            recipe. Default ``None`` makes the spec LAZY: it holds the recipe and
            resolves at build against the model's own grid, which the model's
            context must then carry; the model's discretization must be declared
            before this package.
        nper : int
            Number of stress periods.
        rate : str, float, list of str, or Spread, default "rate"
            Volumetric rate, L^3/T, **negative for withdrawal**. A value spec: an attribute
            column name, a numeric constant, a *list* of column names (indexed by stress
            period), or a :class:`~myflopy.geopackage.CellSurfaceOffset`. A list of column names
            gives one rate per stress period; ``period_field`` is the alternative when the file
            holds one row per feature per period. Wrap it in :class:`~myflopy.geopackage.Spread`
            when the feature spans more than one cell -- a bare column copies the SAME value
            into every cell the feature covers, which is right for an elevation and multiplies a
            conductance by the cell count.
        layer : str, optional
            Name of the layer/table inside the file. Defaults to the first, which is
            the only one in a single-layer GeoPackage or a shapefile.
        name_field : str, optional, default "name"
            Attribute column used for boundnames. Set ``boundnames=False`` if the
            file has no such column, or the build raises for the missing field.
        layer_field : str, optional, default "layer"
            Attribute column giving each feature's model layer. Cells inactive in that layer are
            skipped, as are vertical-passthrough cells (``idomain < 0``). ``None`` puts every
            feature in layer 1.
        period_field : str, optional
            Attribute column giving the stress period a feature belongs to. ``None``
            (the default) applies every feature to all ``nper`` periods. Use it when
            the file holds one row per feature per period; for values that change in
            time on a fixed set of features, pass a *list* of column names as the
            value instead -- it is indexed by period.
        layer_base : int, default 1
            Index base of ``layer_field``. ``1`` matches the usual GIS convention of
            numbering layers from one; pass ``0`` for a zero-based column.
        period_base : int, default 0
            Index base of ``period_field``. ``0`` matches MODFLOW's zero-based
            periods as myflopy counts them.
        gaps : {"warn", "nearest", "raise", "ignore"}, optional
            What a feature that lands on no active cell means. ``"warn"`` (the
            default) leaves it out and says so -- a ``UserWarning`` as well as the
            log, naming each feature -- in the same message that reports features
            clipped at the grid's edge and elevations below their cell bottom.
            ``"nearest"`` moves a missed POINT onto the nearest usable cell when it
            is within ``snap_within`` and reports the move; a line or polygon is
            never moved, and is reported as under ``"warn"``. ``"raise"`` raises
            ``ValueError`` wherever ``"warn"`` would warn. ``"ignore"`` says
            nothing. In every mode a source with NO feature landing raises
            ``ValueError`` naming both extents and CRSs; an unknown word raises
            ``ValueError`` at the call.
        snap_within : float, optional
            ``gaps="nearest"`` only (ignored otherwise): how far a missed point may
            be moved, in the grid's coordinate units, measured to the cell POLYGON
            rather than its centroid. Must be positive (``ValueError`` at the
            call). Default ``None`` is one cell width -- ``sqrt(area)`` of the cell
            it would move to -- resolved when the records resolve, so on the lazy
            form it follows the mesh across ``with_grid``.
        name : str, default "wel"
            Package name.
        boundnames : bool, default True
            Enable named boundaries.
        **options
            Extra ``flopy.mf6.ModflowGwfwel`` options.

        Returns
        -------
        PackageSpec

        Examples
        --------
        >>> mf.wel.gpkg("wells.gpkg", context=ctx, nper=12)
        """

        if context is None:
            return _lazy_bc_spec(
                "wel",
                path,
                name=name,
                nper=nper,
                layer=layer,
                name_field=name_field,
                layer_field=layer_field,
                period_field=period_field,
                layer_base=layer_base,
                period_base=period_base,
                gaps=gaps,
                snap_within=snap_within,
                values={"rate": rate},
                boundnames=boundnames,
                options=options,
            )

        return _source(
            path,
            context=context,
            nper=nper,
            layer=layer,
            name_field=name_field,
            layer_field=layer_field,
            period_field=period_field,
            layer_base=layer_base,
            period_base=period_base,
            gaps=gaps,
            snap_within=snap_within,
        ).wel(
            rate=rate,
            name=name,
            boundnames=boundnames,
            **options,
        )


class _RCHPackage:
    """Package-first RCH (recharge) helpers.

    Areally-distributed recharge applied to the top active cell of each column.
    Entry points:

    - ``mf.rch(context=, nper=, recharge=)`` -- high-level builder; picks the
      top-active cells from the domain and applies a scalar/per-cell/``{period:}`` rate;
    - ``mf.rch(stress_period_data=...)`` / ``mf.rch.flopy(...)`` -- direct MF6 records;
    - ``mf.rch.gpkg(path, context=, nper=)`` -- from GeoPackage features.

    Examples
    --------
    >>> ctx = mf.ModelContext(grid=vor, domain=idomain)
    >>> mf.rch(context=ctx, nper=1, recharge=6.0e-4)                 # uniform rate
    >>> mf.rch.flopy(stress_period_data={0: [[(0, 3), 6.0e-4]]})     # explicit cells
    """

    def __call__(
        self,
        *,
        stress_period_data: Any = None,
        context: ModelContext | None = None,
        nper: int | None = None,
        recharge: Any = None,
        cells: str | Sequence[int | tuple[int, int]] = "top_active",
        layer: int = 0,
        name_by_cell: Mapping[int | tuple[int, int], str] | None = None,
        boundnames: bool = False,
        name: str = "rch",
        **options: Any,
    ) -> PackageSpec:
        """Return an RCH package spec, either from direct data or the domain-aware builder.

        Two mutually exclusive forms:

        * **Direct** -- pass ``stress_period_data=`` (like ``mf.drn`` / ``mf.wel``)
          for a list-based RCH spec.
        * **Builder** -- pass ``context=``, ``nper=``, and ``recharge=`` to compute the
          recharge cells from the model domain automatically.

        Parameters
        ----------
        stress_period_data : dict, optional
            Direct form: ``{period: [[cellid, recharge], ...]}``.
        context : ModelContext, optional
            Builder form: carries the grid/domain used to select recharge cells.
            Omitted with ``nper`` and a value given, the spec is LAZY and resolves
            at build against the model's grid; it then refuses (``TypeError``,
            naming ``context=``) anything that names cells -- an explicit ``cells``
            list, ``name_by_cell``, a ``{cellid: value}`` mapping and a per-cell
            value of any array-like kind (list, tuple, ndarray, Series), bare or
            inside a period.
        nper : int, optional
            Builder form: number of stress periods.
        recharge : scalar, sequence, or {period: ...}, optional
            Builder form: the recharge rate (L/T) -- a constant, one value per
            selected cell, or a per-period mapping.
        cells : str or sequence, default "top_active"
            Which cells receive recharge: ``"top_active"`` (top active cell per
            column), ``"all_active"``, ``"surface_only"``, or an explicit list of
            ``int`` / ``(layer, cell)`` ids.
        layer : int, default 0
            Layer used when ``cells`` are given as bare cell ints.
        name_by_cell : mapping, optional
            Optional ``{cell: boundname}`` mapping (builder form).
        boundnames : bool, default False
            Enable named boundaries.
        name : str, default "rch"
            Package name.
        **options
            Extra ``flopy.mf6.ModflowGwfrch`` options.

        Returns
        -------
        PackageSpec

        Raises
        ------
        TypeError
            If neither ``stress_period_data`` nor the full builder trio
            (``context``/``nper``/``recharge``) is given.

        Examples
        --------
        >>> ctx = mf.ModelContext(grid=vor, domain=idomain)
        >>> mf.rch(context=ctx, nper=1, recharge=6.0e-4)              # uniform, top-active cells
        >>> mf.rch(stress_period_data={0: [[(0, 3), 6.0e-4]]})        # explicit cells
        """

        if stress_period_data is not None:
            return rch_spec(
                stress_period_data, name=name, boundnames=boundnames, **options
            )
        if nper is None or recharge is None:
            raise TypeError(
                "mf.rch requires either stress_period_data=, or the builder "
                "arguments nper= and recharge= (plus context= to resolve eagerly)."
            )
        if context is None:
            # LAZY: the same declaration, resolved against whatever grid the build
            # supplies. The trigger is exactly the declaration that RAISED before, which
            # is what keeps every stored digest still (ledger 284).
            return _lazy_areal_spec(
                "rch",
                name=name,
                recipe=_explicit_areal(
                    nper=nper, recharge=recharge, cells=cells, layer=layer,
                    name_by_cell=name_by_cell, boundnames=boundnames,
                ),
                options=options,
            )
        return RCHBuilder(
            context=context,
            nper=nper,
            recharge=recharge,
            cells=cells,
            layer=layer,
            name_by_cell=name_by_cell,
            boundnames=boundnames,
            name=name,
            options=options,
        ).build()

    def array(
        self,
        *,
        recharge: Any,
        context: ModelContext | None = None,
        irch: Any = "top_active",
        fixed_cell: bool = False,
        auxiliary: Any = None,
        auxmultname: str | None = None,
        observations: Any = None,
        print_input: bool = False,
        save_flows: bool = True,
        name: str = "rch",
        **options: Any,
    ) -> PackageSpec:
        """Return an array-form (``READASARRAYS``) RCH spec -- one array per period.

        The same physics as :meth:`__call__`, supplied differently: one value per
        column per period instead of one record per cell. Compact and much faster
        when recharge covers the model -- measured on a 72-period, 9405-cell model,
        0.88 s end to end against 5.35 s for the equivalent list, and the FloPy
        constructor alone 413x faster.

        Two modes, decided by what is left to the grid. **Eager**: pass ``context=``
        (or give ``irch`` explicitly) and the layer index is resolved here, at the
        call. **Lazy**: leave ``irch="top_active"`` with no ``context=``, and it is
        resolved at BUILD against the model's own discretization -- so one
        declaration builds on any mesh and saves with its project (ledger 285).

        Parameters
        ----------
        recharge : float, int, array_like, str, or dict of {int: any of those}
            Recharge rate, length/time in the model's units (ft/d on a model in
            feet and days), applied over each column's plan area. Required.

            A number writes MODFLOW's ``CONSTANT``. An array or list of ``ncpl``
            values (a ``(1, ncpl)`` array works too) gives each column its own rate,
            INCLUDING inactive columns; a wrong length is refused by FloPy at build
            (``MFDataException``), not here. A ``str`` is not a number: FloPy writes
            it as the name of an external array file (``OPEN/CLOSE``), which MF6
            must be able to open from the run's workspace. A ``{period: value}``
            dict, zero-based periods, varies it in time; a period left out repeats
            the one before it (measured: ``{0: 1e-3, 2: 0.0}`` over three periods
            applied 360, 360 and 0 ft3/d). A single value applies to every period.
        context : ModelContext, optional
            Where ``irch="top_active"`` is resolved from: only ``context.domain`` is
            read for it, and it must be ``(nlay, ncpl)``. When the context also
            carries a ``grid``, the resolved index is PINNED to that mesh
            (``metadata['grid_pin']``), so moving the model to another mesh refuses
            it rather than applying one mesh's layer per column to another's cells.
            Omit it, with ``irch`` left at ``"top_active"``, for the LAZY mode. With
            an explicit ``irch`` it resolves nothing and only pins.
        irch : {"top_active"}, array_like of int, dict of {int: array_like of int}, or None, default "top_active"
            The ZERO-based layer each column's recharge enters. FloPy adds one when
            it writes MODFLOW's 1-based IRCH, so 1-based values shift every column
            down a layer, and one pointing past the bottom stops MF6 with
            ``Invalid layer number`` (measured: 2 on a two-layer model).

            ``"top_active"`` is each column's highest ACTIVE cell -- ``idomain >
            0``, the one rule every "top active" in myflopy uses
            (:func:`~myflopy.modflow.mf6.grid.geometry.top_active_layer`). A
            vertical pass-through cell (``idomain == -1``) is skipped: MF6 has
            removed it, and recharge addressed to it is dropped without a word.
            This read ``idomain != 0`` until ledger 319, and on a stack pinched with
            ``pinch="passthrough"`` it lost 42.2% of a model's recharge under Normal
            termination. A column with no active cell gets 0; its cell is inactive
            either way. Resolved at the call with ``context=``, at build without.

            An array or list gives the layer per column, a ``{period: array}`` dict
            per period. ``None`` writes no IRCH, and MF6 then applies recharge to
            layer 1 in every column and SILENTLY skips any column whose layer 1 is
            inactive -- measured, a 10-column model with 4 such holes lost 40% of
            its recharge under Normal termination, and ``fixed_cell`` does not
            rescue it. Correct only where layer 1 is active everywhere.

            **Checked at build against the DISV it is built beside** (ledgers
            321-322). A per-column index -- an array, or ``"top_active"`` resolved
            eagerly or at build -- that names a REMOVED cell (``idomain <= 0``) in a
            column that has an active cell is REFUSED: it was resolved against
            another domain -- a ``context=`` from before the stack was excavated, or,
            declared without one, the MODEL's context ``domain=`` -- and MF6 would
            drop that column's recharge without a word. ``None`` or one integer for
            every column stays legal, and the build WARNS naming the columns and the
            rate lost (exact: MF6's budget shows that much less recharge).
        fixed_cell : bool, default False
            Write ``FIXED_CELL``, MF6's option to keep recharge in the IRCH cell
            even when that cell is inactive -- a convertible cell gone DRY included
            -- where by default MF6 passes it to the highest active cell below.
            What it loses is the dry cells' share. Measured
            (``docs/reproducers/top_active_passthrough.py``) on a model where 207.98
            of 360.00 ft3/d falls on 52 dry layer-1 cells: the default applied all
            360.00, ``fixed_cell=True`` 152.02 -- only the columns whose IRCH cell
            was wet. Neither setting reaches a cell ``idomain`` removed: that is why
            ``"top_active"`` skips a pass-through cell rather than addressing it.
        auxiliary : str or list of str, optional
            Names of auxiliary variables, written in OPTIONS. Their values go in
            FloPy's ``aux=`` keyword through ``**options``:
            ``aux={period: [array, ...]}``, one ``ncpl`` array per name, in order.
        auxmultname : str, optional
            The auxiliary variable, one of ``auxiliary``, that MULTIPLIES the
            recharge rate -- measured, an aux array of 0.5 took 360 ft3/d to 180.
        observations : dict, optional
            FloPy's observation mapping, ``{"file.csv": [(obsname, "rch", cellid),
            ...]}`` with a zero-based ``(layer, cell)`` cellid. MF6 writes the
            recharge into that cell, rate times area, to the CSV (measured 4.375
            ft3/d at 1e-3 ft/d on a 4,375 ft2 cell). An observation names a cell:
            boundnames do not exist under ``READASARRAYS``.
        print_input : bool, default False
            Write ``PRINT_INPUT``: MF6 echoes the arrays to the model listing file.
        save_flows : bool, default True
            Write ``SAVE_FLOWS``: MF6 stores cell-by-cell recharge in the budget file
            ``mf.oc``'s budget record names, which the results tier reads.
        name : str, default "rch"
            The package name, which also names its file (``<model>.<name>``). Two
            recharge packages in one model need two names.
        **options
            Any other ``flopy.mf6.ModflowGwfrcha`` keyword, e.g. ``print_flows=True``
            or ``aux=``. One that class does not have raises here, at the call.

        Returns
        -------
        PackageSpec
            EAGER: ``options['irch']`` holds the resolved array and ``inputs`` is
            empty. LAZY: ``inputs == {'irch': 'top_active'}``, which is also what
            its digest holds.

        Raises
        ------
        ValueError
            ``irch="top_active"`` with a ``context=`` that carries no domain, or a
            domain that is not ``(nlay, ncpl)``; an ``irch`` string other than
            ``"top_active"`` -- at the call with ``context=``, at build without it;
            ``boundnames=True``, which MF6 refuses under ``READASARRAYS``.
        TypeError
            A ``mf.Spread`` or ``mf.CellSurfaceOffset`` in any keyword: both describe
            a per-FEATURE resolution and the array form has no features. Also a
            keyword ``ModflowGwfrcha`` does not accept.

        Notes
        -----
        Array form needs a value for EVERY column, including inactive ones, and
        is unavailable on DISU. For a boundary covering a small part of the grid
        the list form is smaller and faster -- 40 cells of 9405 is ~1 kB as a
        list against ~160 kB as an array.

        Examples
        --------
        >>> import numpy as np
        >>> import myflopy as mf
        >>> idomain = np.array([[1, -1, 0],
        ...                     [1,  1, 1]])
        >>> spec = mf.rch.array(recharge=1e-3, context=mf.ModelContext(domain=idomain))
        >>> spec.options["irch"].tolist()     # the -1 column's recharge enters layer 1
        [0, 1, 1]
        >>> mf.rch.array(recharge={0: 1e-3, 1: 0.0}).inputs    # lazy: resolved at build
        {'irch': 'top_active'}
        """

        if options.pop("boundnames", None):
            raise ValueError(
                "mf.rch.array: MODFLOW 6 does not allow boundnames under READASARRAYS. "
                "Use mf.rch(...) for the list form, which does."
            )
        return _array_bc_result(
            flopy.mf6.ModflowGwfrcha,
            name=name,
            index_name="irch",
            verb="mf.rch.array",
            fields={"recharge": recharge},
            context=context,
            index=irch,
            fixed_cell=fixed_cell,
            auxiliary=auxiliary,
            auxmultname=auxmultname,
            observations=observations,
            print_input=print_input,
            save_flows=save_flows,
            **options,
        )

    def flopy(
        self,
        *,
        stress_period_data: Any,
        name: str = "rch",
        boundnames: bool = False,
        **options: Any,
    ) -> PackageSpec:
        """Return a direct (list-based) FloPy-style RCH spec from stress-period data.

        Parameters
        ----------
        stress_period_data : dict
            FloPy mapping ``{period: [[cellid, recharge], ...]}``.
        name : str, default "rch"
            Package name.
        boundnames : bool, default False
            Enable named boundaries.
        **options
            Extra ``flopy.mf6.ModflowGwfrch`` options.

        Returns
        -------
        PackageSpec

        Examples
        --------
        >>> mf.rch.flopy(stress_period_data={0: [[(0, 3), 6.0e-4]]})
        """

        return rch_spec(stress_period_data, name=name, boundnames=boundnames, **options)

    def gpkg(
        self,
        path: PathLike,
        *,
        context: ModelContext | None = None,
        nper: int,
        recharge: RowValue = "recharge",
        layer: str | None = None,
        name_field: str | None = "name",
        layer_field: str | None = "layer",
        period_field: str | None = None,
        layer_base: int = 1,
        period_base: int = 0,
        gaps: str = "warn",
        snap_within: float | None = None,
        name: str = "rch",
        boundnames: bool = True,
        **options: Any,
    ) -> PackageSpec:
        """Return a list-based RCH package spec built from GeoPackage features.

        Parameters
        ----------
        path : Path or str
            GeoPackage/shapefile of recharge-zone features (auto-reprojected to the grid CRS).
        context : ModelContext, optional
            The grid/domain the features are mapped onto. Given, the records resolve
            NOW, against ``context.grid`` (``ValueError`` when it carries none) and
            ``context.domain``, and ``with_grid`` re-declares them from the recorded
            recipe. Default ``None`` makes the spec LAZY: it holds the recipe and
            resolves at build against the model's own grid, which the model's
            context must then carry; the model's discretization must be declared
            before this package.
        nper : int
            Number of stress periods.
        recharge : str, float, list of str, or Spread, default "recharge"
            Recharge rate, L/T. A value spec: an attribute column name, a numeric constant, a
            *list* of column names (indexed by stress period), or a
            :class:`~myflopy.geopackage.CellSurfaceOffset`.
            A rate is intensive, so it is copied to every cell a feature covers --
            which is correct. Wrap it in :class:`~myflopy.geopackage.Spread` only if
            the column actually holds a total volume to divide up.
            For grid-wide recharge prefer :meth:`array`, which is much faster.
        layer : str, optional
            Name of the layer/table inside the file. Defaults to the first, which is
            the only one in a single-layer GeoPackage or a shapefile.
        name_field : str, optional, default "name"
            Attribute column used for boundnames. Set ``boundnames=False`` if the
            file has no such column, or the build raises for the missing field.
        layer_field : str, optional, default "layer"
            Attribute column giving each feature's model layer. Recharge normally belongs in the
            uppermost active layer; ``mf.rch.array`` resolves that per column with ``irch=``
            instead. ``None`` puts every feature in layer 1.
        period_field : str, optional
            Attribute column giving the stress period a feature belongs to. ``None``
            (the default) applies every feature to all ``nper`` periods. Use it when
            the file holds one row per feature per period; for values that change in
            time on a fixed set of features, pass a *list* of column names as the
            value instead -- it is indexed by period.
        layer_base : int, default 1
            Index base of ``layer_field``. ``1`` matches the usual GIS convention of
            numbering layers from one; pass ``0`` for a zero-based column.
        period_base : int, default 0
            Index base of ``period_field``. ``0`` matches MODFLOW's zero-based
            periods as myflopy counts them.
        gaps : {"warn", "nearest", "raise", "ignore"}, optional
            What a feature that lands on no active cell means. ``"warn"`` (the
            default) leaves it out and says so -- a ``UserWarning`` as well as the
            log, naming each feature -- in the same message that reports features
            clipped at the grid's edge and elevations below their cell bottom.
            ``"nearest"`` moves a missed POINT onto the nearest usable cell when it
            is within ``snap_within`` and reports the move; a line or polygon is
            never moved, and is reported as under ``"warn"``. ``"raise"`` raises
            ``ValueError`` wherever ``"warn"`` would warn. ``"ignore"`` says
            nothing. In every mode a source with NO feature landing raises
            ``ValueError`` naming both extents and CRSs; an unknown word raises
            ``ValueError`` at the call.
        snap_within : float, optional
            ``gaps="nearest"`` only (ignored otherwise): how far a missed point may
            be moved, in the grid's coordinate units, measured to the cell POLYGON
            rather than its centroid. Must be positive (``ValueError`` at the
            call). Default ``None`` is one cell width -- ``sqrt(area)`` of the cell
            it would move to -- resolved when the records resolve, so on the lazy
            form it follows the mesh across ``with_grid``.
        name : str, default "rch"
            Package name.
        boundnames : bool, default True
            Enable named boundaries.
        **options
            Extra ``flopy.mf6.ModflowGwfrch`` options.

        Returns
        -------
        PackageSpec

        Examples
        --------
        >>> mf.rch.gpkg("recharge_zones.gpkg", context=ctx, nper=12)
        """

        if context is None:
            return _lazy_bc_spec(
                "rch",
                path,
                name=name,
                nper=nper,
                layer=layer,
                name_field=name_field,
                layer_field=layer_field,
                period_field=period_field,
                layer_base=layer_base,
                period_base=period_base,
                gaps=gaps,
                snap_within=snap_within,
                values={"recharge": recharge},
                boundnames=boundnames,
                options=options,
            )

        return _source(
            path,
            context=context,
            nper=nper,
            layer=layer,
            name_field=name_field,
            layer_field=layer_field,
            period_field=period_field,
            layer_base=layer_base,
            period_base=period_base,
            gaps=gaps,
            snap_within=snap_within,
        ).rch(
            recharge=recharge,
            name=name,
            boundnames=boundnames,
            **options,
        )


class _EVTPackage:
    """Package-first EVT (evapotranspiration) helpers.

    Areally-distributed ET drawn from the water table: full ``rate`` at the ET
    ``surface``, decaying linearly to zero at extinction ``depth`` below it (the
    default single-segment form; segmented curves via ``nseg``). Entry points:

    - ``mf.evt(context=, nper=, rate=, depth=)`` -- high-level builder; picks the
      top-active cells from the domain and resolves the ET ``surface`` (default
      model top / land surface) through the surface engine, like ``mf.rch(context=...)``;
    - ``mf.evt(stress_period_data=...)`` -- direct MF6 records;
    - ``mf.evt.gpkg(path, context=, nper=)`` -- build records from GeoPackage
      features mapped onto cells (areal polygons, like ``mf.rch.gpkg``);
    - ``mf.evt.flopy(...)`` -- the FloPy-native form.

    Examples
    --------
    >>> ctx = mf.ModelContext(grid=vor, domain=idomain, surfaces=vor.gdf_topbtm)
    >>> mf.evt(context=ctx, nper=1, rate=2.0e-3, depth=2.5)              # top-active cells, surface = cell top
    >>> mf.evt(stress_period_data={0: [[(0, 3), 100.0, 2.0e-3, 2.5]]})   # (cellid, surface, rate, depth)
    >>> mf.evt.gpkg("et_zones.gpkg", context=ctx, nper=12)
    """

    def __call__(
        self,
        *,
        stress_period_data: Any = None,
        context: ModelContext | None = None,
        nper: int | None = None,
        rate: Any = None,
        depth: Any = None,
        surface: Any = None,
        cells: str | Sequence[int | tuple[int, int]] = "top_active",
        layer: int = 0,
        name_by_cell: Mapping[int | tuple[int, int], str] | None = None,
        name: str = "evt",
        nseg: int = 1,
        auxiliary: Any = None,
        boundnames: bool = False,
        **options: Any,
    ) -> PackageSpec:
        """Return an EVT package spec, either from direct data or the domain-aware builder.

        Two mutually exclusive forms:

        * **Direct** -- pass ``stress_period_data=`` for a list-based EVT spec.
        * **Builder** -- pass ``context=``, ``nper=``, ``rate=`` and ``depth=`` to
          compute the ET cells from the model domain automatically and resolve the
          ET ``surface`` (default cell top) through the surface engine.

        Parameters
        ----------
        stress_period_data : dict, optional
            Direct form: ``{period: [[cellid, surface, rate, depth], ...]}`` --
            ET surface elevation, maximum ET rate (L/T), and extinction depth
            per ``(layer, cell)`` on a DISV grid (the ``nseg=1`` record form).
        context : ModelContext, optional
            Builder form: carries the grid/domain used to select ET cells and the
            surfaces used to resolve the ET surface elevation.
            Omitted with ``nper`` and a value given, the spec is LAZY and resolves
            at build against the model's grid; it then refuses (``TypeError``,
            naming ``context=``) anything that names cells -- an explicit ``cells``
            list, ``name_by_cell``, a ``{cellid: value}`` mapping and a per-cell
            value of any array-like kind (list, tuple, ndarray, Series), bare or
            inside a period.
        nper : int, optional
            Builder form: number of stress periods.
        rate, depth : scalar, sequence, or {period: ...}, optional
            Builder form: the maximum ET flux rate (L/T) and extinction depth (L)
            -- a constant, one value per selected cell, or a per-period mapping.
        surface : optional
            Builder form: the ET surface elevation. Defaults to the model top /
            land surface (``CellSurfaceOffset("model_top")``); also accepts a
            numeric constant, a per-cell sequence / ``{cellid: value}`` mapping,
            or a numeric ``CellSurfaceOffset``.
        cells : str or sequence, default "top_active"
            Builder form: which cells receive ET (``"top_active"``,
            ``"all_active"``, or explicit ``int`` / ``(layer, cell)`` ids).
        layer : int, default 0
            Builder form: layer used when ``cells`` are bare cell ints.
        name_by_cell : mapping, optional
            Builder form: optional ``{cell: boundname}`` mapping.
        name : str, default "evt"
            Package name.
        nseg : int, default 1
            Direct form: number of ET segments; segmented records add
            ``pxdp``/``petm`` values (see :func:`myflopy.advanced.evt_spec`). The
            builder form is single-segment only and now RAISES on ``nseg > 1``
            rather than silently building one segment (ledger 283).
        auxiliary : optional
            Auxiliary variable name(s) forwarded to FloPy. Reaches the builder form
            too -- it used to be dropped there while ``mf.rch`` kept it (ledger 283).
        boundnames : bool, default False
            Enable named boundaries.
        **options
            Extra ``flopy.mf6.ModflowGwfevt`` options.

        Returns
        -------
        PackageSpec

        Raises
        ------
        TypeError
            If neither ``stress_period_data`` nor the full builder set
            (``context``/``nper``/``rate``/``depth``) is given.

        Examples
        --------
        >>> ctx = mf.ModelContext(grid=vor, domain=idomain, surfaces=vor.gdf_topbtm)
        >>> mf.evt(context=ctx, nper=1, rate=2.0e-3, depth=2.5)             # top-active cells
        >>> mf.evt(stress_period_data={0: [[(0, 3), 100.0, 2.0e-3, 2.5]]})  # explicit records
        """

        if stress_period_data is not None:
            return evt_spec(
                stress_period_data,
                name=name,
                nseg=nseg,
                auxiliary=auxiliary,
                boundnames=boundnames,
                **options,
            )
        if nper is None or rate is None or depth is None:
            raise TypeError(
                "mf.evt requires either stress_period_data=, or the builder "
                "arguments nper=, rate= and depth= (plus context= to resolve eagerly)."
            )
        if context is None:
            if nseg != 1:
                raise ValueError(
                    "mf.evt(...): the domain-aware builder writes single-segment ET "
                    "only. Use mf.evt.flopy(stress_period_data=..., nseg=...)."
                )
            recipe = _explicit_areal(
                nper=nper, rate=rate, depth=depth, cells=cells, layer=layer,
                name_by_cell=name_by_cell, boundnames=boundnames,
            )
            if surface is not None:
                recipe["surface"] = surface
            lazy_options = dict(options)
            if auxiliary is not None:
                lazy_options["auxiliary"] = auxiliary
            return _lazy_areal_spec("evt", name=name, recipe=recipe, options=lazy_options)
        if nseg != 1:
            # `nseg` and `auxiliary` are NAMED parameters, so neither ever fell into
            # `**options`, and the builder branch assembled its kwargs without them:
            # a caller asking for segmented ET got single-segment ET and a caller
            # asking for an auxiliary variable got none, both without a word
            # (ledger 283). `mf.rch` keeps `auxiliary` because it has no parameter of
            # that name -- the defect was asymmetric inside one helper pair.
            raise ValueError(
                "mf.evt(context=...): the domain-aware builder writes single-segment ET "
                "only -- it has no way to place the pxdp/petm breakpoints a segmented "
                "curve needs. Use mf.evt.flopy(stress_period_data=..., nseg=...) to "
                "supply the segmented records yourself."
            )
        builder_options = dict(options)
        if auxiliary is not None:
            builder_options["auxiliary"] = auxiliary
        builder_kwargs: dict[str, Any] = dict(
            context=context,
            nper=nper,
            rate=rate,
            depth=depth,
            cells=cells,
            layer=layer,
            name_by_cell=name_by_cell,
            boundnames=boundnames,
            name=name,
            options=builder_options,
        )
        if surface is not None:
            builder_kwargs["surface"] = surface
        return EVTBuilder(**builder_kwargs).build()

    def array(
        self,
        *,
        rate: Any,
        depth: Any,
        surface: Any = None,
        context: ModelContext | None = None,
        ievt: Any = "top_active",
        fixed_cell: bool = False,
        auxiliary: Any = None,
        auxmultname: str | None = None,
        observations: Any = None,
        print_input: bool = False,
        save_flows: bool = True,
        name: str = "evt",
        **options: Any,
    ) -> PackageSpec:
        """Return an array-form (``READASARRAYS``) EVT spec -- one array per period.

        The same physics as :meth:`__call__`, one value per column per period
        instead of one record per cell, with the same two modes as
        :meth:`_RCHPackage.array`: **eager** with ``context=`` (resolved at the
        call), **lazy** without it, when ``ievt`` is left at ``"top_active"`` or
        ``surface`` is omitted (resolved at BUILD from the model's discretization,
        ledger 285).

        Parameters
        ----------
        rate : float, int, array_like, str, or dict of {int: any of those}
            Maximum ET rate, length/time in the model's units (ft/d on a model in
            feet and days), taken where the head is at or above ``surface`` and
            falling to zero at ``surface - depth``. Required. Shapes as ``recharge``
            on :meth:`_RCHPackage.array`: a number, ``ncpl`` values, an external
            file name (``str``), or a zero-based ``{period: value}`` dict whose
            missing periods repeat the one before (measured: ``{0: 1e-3, 1: 0.0}``
            took 162.35 ft3/d, then 0).
        depth : float, int, array_like, str, or dict of {int: any of those}
            Extinction depth below ``surface``, length (ft on a model in feet).
            Required. Same shapes as ``rate``.
        surface : float, int, array_like, or dict of {int: any of those}, optional
            The ET surface ELEVATION, length, where the maximum rate applies. Same
            numeric shapes as ``rate``.

            Omitted, it is the MODEL TOP -- resolved at the call from a
            ``context=`` carrying ``grid`` and ``surfaces`` (a ``LayerBuildResult``
            or a geometry-first table), or at build without one. A flat top resolves
            to one number. It used to be left to FloPy, whose default is the datum
            ``0.0``: MF6 then draws the full rate everywhere because every head is
            above the ET surface, and ``depth`` does nothing (ledger 283).

            Two inputs are NOT what they look like (found 2026-09-25, not changed;
            ledger 319): a ``str`` is read as "the model top" whatever it says --
            ``surface="dem.tif"`` wrote the same 162.35 ft3/d as omitting it -- and a
            ``mf.Surface`` is accepted here and refused by FloPy at build
            (``MFDataException``). Sample a raster onto the cells yourself and pass
            the array.
        context : ModelContext, optional
            As on :meth:`_RCHPackage.array`: ``context.domain`` resolves
            ``ievt="top_active"``, and a ``grid`` pins the result to that mesh. With
            ``surface`` omitted it must ALSO carry ``grid`` and ``surfaces``, or the
            call raises. Omit it for the LAZY mode.
        ievt : {"top_active"}, array_like of int, dict of {int: array_like of int}, or None, default "top_active"
            The ZERO-based layer each column's ET is taken from, with exactly the
            rules of ``irch`` on :meth:`_RCHPackage.array`: ``"top_active"`` is the
            highest cell with ``idomain > 0``, so a pass-through cell (``-1``) is
            skipped (ledger 319); 1-based values shift every column down a layer;
            ``None`` takes layer 1 everywhere and silently skips any column whose
            layer 1 is inactive. Checked at build as ``irch`` is (ledgers 321-322):
            a per-column index naming a removed cell in a column with an active cell
            is refused, ``None`` or one layer everywhere is warned with the rate
            lost -- for ET an UPPER bound, "up to" ``rate`` times the lost area at
            the maximum ET rate, because MF6 scales ET by the head's depth below
            ``surface`` and so drops anything from zero to that. A TYPED value is
            written as given in both modes -- with ``surface`` left to the build it
            was dropped, and MF6 took the ET from layer 1 (ledger 324) -- and a numpy
            integer is converted to ``int`` for FloPy.
        fixed_cell : bool, default False
            Write ``FIXED_CELL``, MF6's option to keep ET in the IEVT cell instead of
            passing it to the cell below when that cell is inactive.
        auxiliary : str or list of str, optional
            Names of auxiliary variables, written in OPTIONS; their arrays go in
            FloPy's ``aux=`` through ``**options``, ``aux={period: [array, ...]}``.
        auxmultname : str, optional
            The auxiliary variable, one of ``auxiliary``, that MULTIPLIES the maximum
            ET rate -- measured, an aux array of 0.5 took 162.35 ft3/d to 81.82 (the
            heads rose, so not exactly half).
        observations : dict, optional
            FloPy's observation mapping, ``{"file.csv": [(obsname, "evt", cellid),
            ...]}`` with a zero-based ``(layer, cell)`` cellid; MF6 writes that
            cell's ET to the CSV. No boundnames under ``READASARRAYS``.
        print_input : bool, default False
            Write ``PRINT_INPUT``: MF6 echoes the arrays to the model listing file.
        save_flows : bool, default True
            Write ``SAVE_FLOWS``: MF6 stores cell-by-cell ET in the budget file
            ``mf.oc``'s budget record names.
        name : str, default "evt"
            The package name, which also names its file (``<model>.<name>``).
        **options
            Any other ``flopy.mf6.ModflowGwfevta`` keyword, e.g. ``print_flows=True``
            or ``aux=``; one that class does not have raises here. ``nseg=1`` is
            accepted and dropped.

        Returns
        -------
        PackageSpec
            EAGER: ``options`` holds the resolved ``ievt`` and ``surface``, and
            ``inputs`` is empty. LAZY: ``inputs`` holds what was left to the grid --
            ``{'ievt': 'top_active', 'surface': 'model_top'}`` with both omitted.

        Raises
        ------
        ValueError
            ``nseg`` other than 1 (see Notes); ``boundnames=True``; ``surface``
            omitted with a ``context=`` lacking ``grid`` or ``surfaces``;
            ``ievt="top_active"`` with a ``context=`` carrying no ``(nlay, ncpl)``
            domain; an ``ievt`` string other than ``"top_active"`` -- at the call
            with ``context=``, at build without it.
        TypeError
            A ``mf.Spread`` or ``mf.CellSurfaceOffset`` in any keyword, or a keyword
            ``ModflowGwfevta`` does not accept.

        Notes
        -----
        **Segmented ET cannot be expressed in array form.** ``NSEG`` lives in the
        DIMENSIONS block, which MODFLOW 6 never reads under ``READASARRAYS``, so a
        model with ``nseg > 1`` must use :meth:`__call__` or :meth:`flopy`.
        Passing ``nseg`` here raises rather than silently dropping the segment
        shape -- that is a change to the physics, not to the format.

        Examples
        --------
        >>> import numpy as np
        >>> import myflopy as mf
        >>> ctx = mf.ModelContext(domain=np.array([[1, -1],
        ...                                        [1,  1]]))
        >>> spec = mf.evt.array(rate=2e-3, depth=5.0, surface=100.0, context=ctx)
        >>> spec.options["ievt"].tolist()     # the -1 column's ET comes from layer 1
        [0, 1]
        >>> sorted(mf.evt.array(rate=2e-3, depth=5.0).inputs)   # lazy: both at build
        ['ievt', 'surface']
        """

        if options.pop("nseg", None) not in (None, 1):
            raise ValueError(
                "mf.evt.array: segmented ET (nseg > 1) has no array form -- MODFLOW 6 does "
                "not read the DIMENSIONS block under READASARRAYS, so the segments would be "
                "silently dropped. Use mf.evt(...) or mf.evt.flopy(...) for the list form."
            )
        if options.pop("boundnames", None):
            raise ValueError(
                "mf.evt.array: MODFLOW 6 does not allow boundnames under READASARRAYS. "
                "Use mf.evt(...) for the list form, which does."
            )
        return _array_bc_result(
            flopy.mf6.ModflowGwfevta,
            name=name,
            index_name="ievt",
            verb="mf.evt.array",
            # `_MODEL_TOP` rather than the resolved elevation: with no context there is
            # nothing to resolve it against yet, and leaving the argument out entirely
            # is what wrote MODFLOW's datum (ledger 283).
            fields={"rate": rate, "depth": depth, "surface": _MODEL_TOP if surface is None else surface},
            context=context,
            index=ievt,
            fixed_cell=fixed_cell,
            auxiliary=auxiliary,
            auxmultname=auxmultname,
            observations=observations,
            print_input=print_input,
            save_flows=save_flows,
            **options,
        )

    def flopy(
        self,
        *,
        stress_period_data: Any,
        name: str = "evt",
        nseg: int = 1,
        auxiliary: Any = None,
        boundnames: bool = False,
        **options: Any,
    ) -> PackageSpec:
        """Return a direct (list-based) FloPy-style EVT spec from stress-period data.

        Parameters
        ----------
        stress_period_data : dict
            FloPy mapping ``{period: [[cellid, surface, rate, depth], ...]}``
            (plus ``pxdp``/``petm``/``petm0`` values for segmented forms).
        name : str, default "evt"
            Package name.
        nseg : int, default 1
            Number of ET segments.
        auxiliary : optional
            Auxiliary variable name(s) forwarded to FloPy.
        boundnames : bool, default False
            Enable named boundaries.
        **options
            Extra ``flopy.mf6.ModflowGwfevt`` options.

        Returns
        -------
        PackageSpec

        Examples
        --------
        >>> mf.evt.flopy(stress_period_data={0: [[(0, 3), 100.0, 2.0e-3, 2.5]]})
        """

        return evt_spec(
            stress_period_data,
            name=name,
            nseg=nseg,
            auxiliary=auxiliary,
            boundnames=boundnames,
            **options,
        )

    def gpkg(
        self,
        path: PathLike,
        *,
        context: ModelContext | None = None,
        nper: int,
        surface: RowValue = "surface",
        rate: RowValue = "rate",
        depth: RowValue = "depth",
        layer: str | None = None,
        name_field: str | None = "name",
        layer_field: str | None = "layer",
        period_field: str | None = None,
        layer_base: int = 1,
        period_base: int = 0,
        gaps: str = "warn",
        snap_within: float | None = None,
        name: str = "evt",
        boundnames: bool = True,
        **options: Any,
    ) -> PackageSpec:
        """Return an EVT package spec built from GeoPackage features mapped onto the grid.

        Parameters
        ----------
        path : Path or str
            GeoPackage/shapefile of ET-zone features (auto-reprojected to the grid CRS).
        context : ModelContext, optional
            The grid/domain the features are mapped onto. Given, the records resolve
            NOW, against ``context.grid`` (``ValueError`` when it carries none) and
            ``context.domain``, and ``with_grid`` re-declares them from the recorded
            recipe. Default ``None`` makes the spec LAZY: it holds the recipe and
            resolves at build against the model's own grid, which the model's
            context must then carry; the model's discretization must be declared
            before this package.
        nper : int
            Number of stress periods.
        surface : str, float, list of str, or CellSurfaceOffset, default "surface"
            Elevation ET is measured down from, usually land surface.
            A value spec: an attribute column name, a numeric constant, a
            *list* of column names (indexed by stress period), or a
            :class:`~myflopy.geopackage.CellSurfaceOffset`. ``CellSurfaceOffset("model_top")``
            is the common choice.
        rate : str, float, list of str, default "rate"
            Maximum ET rate, L/T, applied when the head is at ``surface``.
            A value spec: an attribute column name, a numeric constant, a
            *list* of column names (indexed by stress period), or a
            :class:`~myflopy.geopackage.CellSurfaceOffset`.
        depth : str or CellSurfaceOffset, default "depth"
            Extinction depth below the ET surface -- an attribute column or a
            constant.
        layer : str, optional
            Name of the layer/table inside the file. Defaults to the first, which is
            the only one in a single-layer GeoPackage or a shapefile.
        name_field : str, optional, default "name"
            Attribute column used for boundnames. Set ``boundnames=False`` if the
            file has no such column, or the build raises for the missing field.
        layer_field : str, optional, default "layer"
            Attribute column giving each feature's model layer. ET normally belongs in the
            uppermost active layer; ``mf.evt.array`` resolves that per column with ``ievt=``
            instead. ``None`` puts every feature in layer 1.
        period_field : str, optional
            Attribute column giving the stress period a feature belongs to. ``None``
            (the default) applies every feature to all ``nper`` periods. Use it when
            the file holds one row per feature per period; for values that change in
            time on a fixed set of features, pass a *list* of column names as the
            value instead -- it is indexed by period.
        layer_base : int, default 1
            Index base of ``layer_field``. ``1`` matches the usual GIS convention of
            numbering layers from one; pass ``0`` for a zero-based column.
        period_base : int, default 0
            Index base of ``period_field``. ``0`` matches MODFLOW's zero-based
            periods as myflopy counts them.
        gaps : {"warn", "nearest", "raise", "ignore"}, optional
            What a feature that lands on no active cell means. ``"warn"`` (the
            default) leaves it out and says so -- a ``UserWarning`` as well as the
            log, naming each feature -- in the same message that reports features
            clipped at the grid's edge and elevations below their cell bottom.
            ``"nearest"`` moves a missed POINT onto the nearest usable cell when it
            is within ``snap_within`` and reports the move; a line or polygon is
            never moved, and is reported as under ``"warn"``. ``"raise"`` raises
            ``ValueError`` wherever ``"warn"`` would warn. ``"ignore"`` says
            nothing. In every mode a source with NO feature landing raises
            ``ValueError`` naming both extents and CRSs; an unknown word raises
            ``ValueError`` at the call.
        snap_within : float, optional
            ``gaps="nearest"`` only (ignored otherwise): how far a missed point may
            be moved, in the grid's coordinate units, measured to the cell POLYGON
            rather than its centroid. Must be positive (``ValueError`` at the
            call). Default ``None`` is one cell width -- ``sqrt(area)`` of the cell
            it would move to -- resolved when the records resolve, so on the lazy
            form it follows the mesh across ``with_grid``.
        name : str, default "evt"
            Package name.
        boundnames : bool, default True
            Enable named boundaries.
        **options
            Extra ``flopy.mf6.ModflowGwfevt`` options. Note ``nseg`` must stay
            ``1`` here: mapped features yield ``(cellid, surface, rate, depth)``
            records with no ``pxdp``/``petm`` values, so segmented ET is
            rejected (build segmented records yourself and pass them to
            ``mf.evt(...)`` / ``mf.evt.flopy(...)``).

        Returns
        -------
        PackageSpec

        Raises
        ------
        ValueError
            If ``nseg`` other than ``1`` is given.

        Examples
        --------
        >>> mf.evt.gpkg("et_zones.gpkg", context=ctx, nper=12)
        """

        if context is None:
            return _lazy_bc_spec(
                "evt",
                path,
                name=name,
                nper=nper,
                layer=layer,
                name_field=name_field,
                layer_field=layer_field,
                period_field=period_field,
                layer_base=layer_base,
                period_base=period_base,
                gaps=gaps,
                snap_within=snap_within,
                values={"surface": surface, "rate": rate, "depth": depth},
                boundnames=boundnames,
                options=options,
            )

        return _source(
            path,
            context=context,
            nper=nper,
            layer=layer,
            name_field=name_field,
            layer_field=layer_field,
            period_field=period_field,
            layer_base=layer_base,
            period_base=period_base,
            gaps=gaps,
            snap_within=snap_within,
        ).evt(
            surface=surface,
            rate=rate,
            depth=depth,
            name=name,
            boundnames=boundnames,
            **options,
        )


class _UZFPackage:
    """Package-first UZF (unsaturated-zone flow) helpers.

    Call ``mf.uzf(...)`` to build a UZF package from domain-aware inputs (it picks
    the land-surface cells and wires up infiltration/ET for you), or
    ``mf.uzf.flopy(...)`` to pass prepared FloPy ``packagedata``/``perioddata``
    directly. UZF adds a vertical unsaturated column above the water table that
    delays and attenuates recharge and can route rejected infiltration / discharge
    to other packages via MVR.
    """

    def __call__(
        self,
        *,
        context: ModelContext | None = None,
        nper: int,
        vks: Any,
        thtr: Any,
        thts: Any,
        thti: Any,
        cells: str | Sequence[tuple[int, int]] = "all_active",
        eps: Any = 4.0,
        surfdep: Any = 0.001,
        finf: Any = 0.0,
        pet: Any = None,
        extdp: Any = None,
        extwc: Any = None,
        ha: Any = None,
        hroot: Any = None,
        rootact: Any = None,
        name: str = "uzf",
        boundnames: bool = False,
        mover: bool = False,
        simulate_et: bool | None = None,
        linear_gwet: bool = False,
        square_gwet: bool = False,
        simulate_gwseep: bool = False,
        unsat_etwc: bool = False,
        unsat_etae: bool = False,
        ntrailwaves: int = 7,
        nwavesets: int = 40,
        **options: Any,
    ) -> PackageSpec:
        """Build a UZF package from domain-aware inputs (delegates to ``UZFBuilder``).

        Cells are selected automatically from the model domain (one UZF column per
        vertically-connected active cell beneath land surface); pass an explicit
        ``cells`` list to override. The soil properties (``surfdep``/``vks``/``thtr``/
        ``thts``/``thti``/``eps``) are set once per UZF cell; the stress inputs
        (``finf``/``pet``/``extdp``/``extwc``/``ha``/``hroot``/``rootact``) are set
        per period on the land-surface cell of each column.

        A per-cell value (the non-scalar forms below) takes one of three shapes: a
        ``{(layer, cell): value}`` mapping that covers every cell it applies to, an
        array shaped like the domain ``(nlay, ncpl)`` read at each cell, or a 1-D
        sequence with one value per cell it applies to (every UZF cell for a soil
        property, every land-surface cell for a stress input), in feature order.
        Any other shape, or a mapping missing a cell, raises ``ValueError``.

        Parameters
        ----------
        context : ModelContext, optional
            Carries the domain (idomain) the UZF cells are picked and ordered from,
            and the grid the records are pinned to. Given, the package resolves
            NOW; automatic ``cells`` then need ``context.domain`` (``ValueError``
            without it). Default ``None`` makes the spec LAZY: it resolves at build
            from the model's domain -- the context's, else the built
            discretization's ``idomain`` -- and needs no grid. With neither, every
            cell is active only when the model's context carries no grid either;
            a context holding a grid but no domain raises ``ValueError`` at build.
            The lazy form refuses (``TypeError``) an explicit ``cells`` list, a
            per-cell value of any array-like kind (list, tuple, ndarray, Series)
            or a cell mapping for any value, and ``mover=True``; each of those
            works with ``context=``.
        nper : int
            Number of stress periods; must be at least 1 (``ValueError``).
        vks : float or per-cell value
            Saturated vertical hydraulic conductivity of the unsaturated zone,
            length/time (model units), per UZF cell. Must be positive.
        thtr : float or per-cell value
            Residual water content (dimensionless), per UZF cell. Must be less
            than ``thts``.
        thts : float or per-cell value
            Saturated water content (dimensionless), per UZF cell; keep it
            consistent with the storage package's specific yield.
        thti : float or per-cell value
            Initial water content (dimensionless), per UZF cell; must lie between
            ``thtr`` and ``thts`` inclusive.
        cells : {"all_active", "surface_only"} or sequence of (int, int), optional
            ``"all_active"`` (default) makes a UZF cell of each column's top active
            cell (``idomain > 0``) and of every active cell directly below it down
            to the first inactive one; ``"surface_only"`` keeps just the top cell.
            An explicit sequence of zero-based ``(layer, cell)`` ids is used as
            given, in that order; they must be unique and, when a domain is known,
            in range and active (``ValueError`` otherwise).
        eps : float or per-cell value, optional
            Brooks-Corey exponent (dimensionless), per UZF cell; must be positive.
            Default ``4.0``.
        surfdep : float or per-cell value, optional
            Surface depression depth, length units, per UZF cell; must not be
            negative. Default ``0.001``.
        finf : float, str, per-cell value or dict, optional
            Applied infiltration rate, length/time, per land-surface cell and
            period. A ``{period: value}`` dict (zero-based int keys) varies it in
            time: a period it omits gets ``0.0``, and a key outside ``nper`` raises
            ``ValueError``. Any other value applies to every period. A string is
            passed to MF6 as a time-series name. A negative number raises
            ``ValueError``. Default ``0.0``.
        pet : float, str, per-cell value or dict, optional
            Potential ET rate, length/time, in the same forms as ``finf``; used
            only when ET is simulated. Must not be negative. Default ``None``
            writes ``0.0``.
        extdp : float, str, per-cell value or dict, optional
            ET extinction depth below land surface, length units, in the same
            forms as ``finf``; used only when ET is simulated. Must not be
            negative. Default ``None`` writes ``0.0``.
        extwc : float, str, per-cell value or dict, optional
            ET extinction water content (dimensionless), in the same forms as
            ``finf``; MF6 uses it only with ``unsat_etwc``. Default ``None``
            writes ``0.0``.
        ha : float, str, per-cell value or dict, optional
            Air-entry potential (head), length units, in the same forms as
            ``finf``; MF6 uses it only with ``unsat_etae``. Default ``None`` writes
            ``0.0``.
        hroot : float, str, per-cell value or dict, optional
            Root potential (head), length units, in the same forms as ``finf``;
            MF6 uses it only with ``unsat_etae``. Default ``None`` writes ``0.0``.
        rootact : float, str, per-cell value or dict, optional
            Root activity function (root length per soil volume), in the same
            forms as ``finf``; MF6 uses it only with ``unsat_etae``. Default
            ``None`` writes ``0.0``.
        name : str, optional
            Package name: FloPy's ``pname``, the input file ``<model>.<name>``, the
            spec's name, and the stem of the UZF budget and convergence output
            files. Default ``"uzf"``.
        boundnames : bool, optional
            Write MF6 ``BOUNDNAMES`` and name every UZF cell
            ``uzf_<layer>_<cell>``. Default ``False``.
        mover : bool, optional
            Write MF6 ``MOVER`` so an ``mf.mvr`` can move rejected infiltration or
            groundwater discharge to another package. Refused on the lazy form
            (a UZF mover endpoint is a cell enumeration that changes with the
            mesh). Default ``False``.
        simulate_et : bool, optional
            MF6 ``SIMULATE_ET``. Default ``None`` turns it on exactly when ``pet``,
            ``extdp`` or ``extwc`` is given; ``True``/``False`` force it.
        linear_gwet : bool, optional
            MF6 ``LINEAR_GWET``: also take ET from groundwater, with MODFLOW-2005's
            linear formulation. Default ``False``.
        square_gwet : bool, optional
            MF6 ``SQUARE_GWET``: also take ET from groundwater at a constant rate
            down to ``extdp``, smoothed to zero there. Default ``False``.
        simulate_gwseep : bool, optional
            MF6 ``SIMULATE_GWSEEP``: discharge groundwater above land surface.
            MF6 no longer recommends it (a DRN is the suggested route). Default
            ``False``.
        unsat_etwc : bool, optional
            MF6 ``UNSAT_ETWC``: unsaturated-zone ET at the PET rate while the water
            content exceeds ``extwc``. Default ``False``.
        unsat_etae : bool, optional
            MF6 ``UNSAT_ETAE``: unsaturated-zone ET from a capillary-pressure
            (Brooks-Corey) formulation using ``ha``/``hroot``/``rootact``.
            Default ``False``.
        ntrailwaves : int, optional
            MF6 ``NTRAILWAVES``, the number of trailing waves; raise it to lower
            the unsaturated-zone mass-balance error. Default ``7``.
        nwavesets : int, optional
            MF6 ``NWAVESETS``, the number of wave sets. Default ``40``.
        **options
            Extra ``flopy.mf6.ModflowGwfuzf`` options (e.g. ``wc_filerecord``,
            ``observations``), passed at build time. A name that class does not
            accept raises ``TypeError`` -- at the call with ``context=``, but only
            at build on the lazy form. Applied last, so they override the
            ``pname``/``filename`` and output records set from ``name``.

        Returns
        -------
        PackageSpec
            Pass it to ``mf.gwf(packages=[...])``.

        Examples
        --------
        >>> ctx = mf.ModelContext(grid=vor, domain=idomain)
        >>> mf.uzf(context=ctx, nper=1, vks=0.25, thtr=0.08, thts=0.34, thti=0.17,
        ...        finf=3.0e-5, pet=1.0e-4, extdp=7.0)
        """

        recipe = {
            "nper": nper, "cells": cells, "boundnames": boundnames, "mover": mover,
            "vks": vks, "thtr": thtr, "thts": thts, "thti": thti, "eps": eps,
            "surfdep": surfdep, "finf": finf, "pet": pet, "extdp": extdp,
            "extwc": extwc, "ha": ha, "hroot": hroot, "rootact": rootact,
            "simulate_et": simulate_et, "linear_gwet": linear_gwet,
            "square_gwet": square_gwet, "simulate_gwseep": simulate_gwseep,
            "unsat_etwc": unsat_etwc, "unsat_etae": unsat_etae,
            "ntrailwaves": ntrailwaves, "nwavesets": nwavesets,
        }
        if context is None:
            # LAZY: the same declaration, resolved against whatever DOMAIN the build
            # supplies. The trigger is exactly the call that raised `TypeError: missing a
            # required argument` before, so no stored digest can move (ledger 284's rule).
            return _lazy_uzf_spec(name=name, recipe=recipe, options=options)
        return UZFBuilder(
            context=context,
            nper=nper,
            vks=vks,
            thtr=thtr,
            thts=thts,
            thti=thti,
            cells=cells,
            eps=eps,
            surfdep=surfdep,
            finf=finf,
            pet=pet,
            extdp=extdp,
            extwc=extwc,
            ha=ha,
            hroot=hroot,
            rootact=rootact,
            name=name,
            boundnames=boundnames,
            mover=mover,
            simulate_et=simulate_et,
            linear_gwet=linear_gwet,
            square_gwet=square_gwet,
            simulate_gwseep=simulate_gwseep,
            unsat_etwc=unsat_etwc,
            unsat_etae=unsat_etae,
            ntrailwaves=ntrailwaves,
            nwavesets=nwavesets,
            options=options,
        ).build()

    def flopy(
        self,
        *,
        packagedata: Any,
        perioddata: Any,
        nuzfcells: int | None = None,
        name: str = "uzf",
        mover: bool = False,
        simulate_et: bool = False,
        **options: Any,
    ) -> PackageSpec:
        """Return a direct FloPy-style UZF spec from prepared package data (escape hatch).

        Use this when you already have MF6 UZF records; otherwise prefer the
        domain-aware ``mf.uzf(...)`` builder.

        Parameters
        ----------
        packagedata : list
            FloPy UZF ``packagedata`` (one row per UZF cell: iuzno, cellid, ...).
        perioddata : dict
            FloPy UZF ``perioddata`` keyed by stress period.
        nuzfcells : int, optional
            Number of UZF cells (inferred from ``packagedata`` when omitted).
        name : str, default "uzf"
            Package name.
        mover : bool, default False
            Enable MVR routing of rejected infiltration / discharge.
        simulate_et : bool, default False
            Enable evapotranspiration simulation.
        **options
            Extra ``flopy.mf6.ModflowGwfuzf`` options.

        Returns
        -------
        PackageSpec
        """

        return uzf_spec(
            packagedata,
            perioddata,
            nuzfcells=nuzfcells,
            name=name,
            mover=mover,
            simulate_et=simulate_et,
            **options,
        )


class _SFRPackage:
    """Package-first SFR (streamflow-routing) helpers.

    Call ``mf.sfr(...)`` to build an SFR network directly from **stream centerline
    geometry** (a geopackage/shapefile path or a GeoDataFrame): it discretizes each
    line into reaches on the grid, resolves the reach-to-reach topology, and writes
    packagedata/connectiondata/perioddata. Use ``mf.sfr.flopy(...)`` to pass prepared
    FloPy records instead. SFR routes streamflow through the model and can exchange
    water with the aquifer (gaining/losing) and with lakes via MVR.
    """

    def __call__(
        self,
        *,
        context: ModelContext,
        nper: int,
        streams: PathLike | Sequence[PathLike] | gpd.GeoDataFrame,
        stream_id: str | None = None,
        from_node: str | None = None,
        to_node: str | None = None,
        connection_mode: str = "automatic",
        connections: tuple[StreamConnection, ...] = (),
        diversions: tuple[StreamDiversion, ...] = (),
        connection_tolerance: float | None = None,
        reverse_streams: tuple[str, ...] = (),
        reach_layer: int | str | Mapping[str, int] = "top_active",
        width: Any = 10.0,
        gradient: Any = 0.001,
        reach_top: Any = None,
        roughness: Any = 0.03,
        streambed_k: Any = 1.0,
        streambed_thickness: Any = 1.0,
        min_reach_length: float = 0.0,
        inflow: Any = None,
        rainfall: Any = None,
        evaporation: Any = None,
        runoff: Any = None,
        status: Any = None,
        name: str = "sfr",
        mover: bool = False,
        length_conversion: float | None = None,
        time_conversion: float | None = None,
        maximum_picard_iterations: int = 1,
        maximum_iterations: int = 1000,
        maximum_depth_change: float = 0.01,
        **options: Any,
    ) -> PackageSpec:
        """Build an SFR package from stream centerlines (delegates to ``SFRBuilder``).

        Each stream line is cut into reaches where it crosses grid cells, the
        stream-to-stream topology is resolved (``connection_mode``, overridden by
        ``connections``), and the packagedata, connectiondata and perioddata are
        written. Everything is resolved HERE, against ``context``'s grid: every
        error named below is raised by this call, and the package is bound to that
        mesh. ``mf.sfr`` has no lazy form, so after ``with_grid`` it reads
        ``stale`` -- declare it again on the new grid. The returned spec carries an
        ``sfr_index`` (each stream's reaches and their cells) that
        :func:`sfr_connection` reads to name a mover endpoint by stream and place.

        **Reach numbering.** Reaches are numbered from 0, stream by stream, each
        stream's reaches running from its first vertex (the headwater) to its last
        (the outlet). Streams are taken in rounds: first every stream no other
        stream is connected into, in stream-table order, then every stream whose
        upstream streams are all placed, and so on -- so a tributary is numbered
        before the stem it joins. A stretch lying ON a face between two cells goes
        to the lower-numbered cell. Every reach number below means this numbering;
        ``SFRBuilder(...).reaches`` shows it.

        **Per-reach values** (``width``, ``gradient``, ``reach_top``,
        ``roughness``, ``streambed_k``, ``streambed_thickness``) take:

        * a number, for every reach;
        * a string naming a column of the stream table, read per stream; any other
          string is written unchanged to every reach as an MF6 time-series name,
          which MF6's SFR takes for ``streambed_k`` (``RHK``) and ``roughness``
          (``MAN``) only -- for the other four it raises ``ValueError`` naming the
          argument, the reach and its stream (ledger 365);
        * a mapping keyed by reach number, stream id or ``(layer, cell)`` cellid,
          looked up in that order per reach. A reach none of them matches has no
          value, and that raises ``ValueError`` naming the argument, the reach and
          its stream;
        * a sequence of one value per reach, in reach order, or one per stream in
          the stream order above -- NOT the file's order. Any other length raises
          ``ValueError``.

        Every per-reach value must end as a finite number (or a time-series name
        where MF6 takes one): a ``NaN`` in a stream-table column raises
        ``ValueError`` the same way.

        **Per-period values** (``inflow``, ``rainfall``, ``evaporation``,
        ``runoff``, ``status``) take ``None`` (no setting); a mapping whose keys are
        ALL integers, read as ``{period: value}`` with zero-based periods, where a
        period left out -- or at or beyond ``nper`` -- gets no setting; or any other
        value, used for every period. So ``{3: 100.0}`` is period 3, not reach 3:
        key a reach inside a period, ``{0: {3: 100.0}}``. One period's value is a
        number or a string (a time-series name), placed as each entry says; a
        mapping ``{stream_id or reach number: value}``, where a stream id means
        EVERY reach of that stream for ``rainfall``, ``evaporation`` and
        ``status``, its FIRST reach for ``inflow`` and a share by reach length for
        ``runoff`` (a stream not in the table raises ``KeyError`` naming the
        streams); or an iterable of explicit rows ``(reach, value, ...)`` written
        as given. A flat list of numbers is therefore read as rows and raises
        ``TypeError``: unlike ``mf.lak``, it is not one value per period. Period
        and reach keys may be Python or numpy integers.

        Parameters
        ----------
        context : ModelContext
            Required; ``None`` raises ``AttributeError``. Carries the grid the
            reaches are cut on (``context.grid``; ``None`` raises ``ValueError``),
            the ``domain`` read by ``reach_layer="top_active"``, the layer
            elevations read by the default ``reach_top`` (``context.surfaces``,
            else ``grid.gdf_topbtm``), and the ``length_units``/``time_units`` the
            default unit conversions come from. The grid's pin is stamped on the
            package.
        nper : int
            Number of stress periods; the perioddata gets one entry per period
            ``0 .. nper - 1``. Less than 1 raises ``ValueError``.
        streams : str, pathlib.Path, sequence of str or Path, or geopandas.GeoDataFrame
            Stream centerlines, one ``LineString`` per stream; any other geometry,
            a ``MultiLineString`` included, raises ``ValueError``. A path or a list
            of paths is read with :func:`geopandas.read_file` and the tables
            concatenated, so the files must share one CRS (geopandas raises
            ``ValueError`` otherwise). The lines are reprojected to the grid's CRS
            when both carry one. A line runs from its first vertex, the headwater,
            to its last, the outlet (see ``reverse_streams``). A stream crossing no
            grid cell -- or left with no reach by ``min_reach_length`` -- raises
            ``ValueError``.
        stream_id : str or None, optional
            Column holding each stream's id. The values are turned into strings and
            must be unique (``ValueError`` otherwise, and when the column is
            missing); they are the ids every stream-keyed input, ``connections``,
            ``diversions`` and :func:`sfr_connection` use. Default ``None``: each
            file's stem when ``streams`` is paths and no stem repeats across the
            rows read (one feature per file, distinct stems), else ``"stream_0"``,
            ``"stream_1"``, ... by row.
        from_node : str or None, optional
            Column naming each stream's upstream node, read only by
            ``connection_mode="nodes"``, which requires it and ``to_node``
            (``ValueError`` when either is ``None`` or not a column). Default
            ``None``.
        to_node : str or None, optional
            Column naming each stream's downstream node, read only by
            ``connection_mode="nodes"``: a stream flows into the stream whose
            ``from_node`` equals its ``to_node``. No match leaves it unresolved (it
            ends the network); more than one raises ``ValueError``. Default
            ``None``.
        connection_mode : str, optional
            How the stream-to-stream topology is inferred. ``"automatic"``
            (default) connects each stream's last vertex to the nearest other
            stream within ``connection_tolerance`` of it; the two nearest equally
            near raise ``ValueError``, and none leaves it unresolved. ``"nodes"``
            reads ``from_node``/``to_node``. ``"explicit"`` infers nothing, so only
            ``connections`` link streams. Any other value raises ``ValueError``. An
            unresolved stream ends the network; the returned spec's
            ``metadata["unresolved"]`` lists them.
        connections : sequence of StreamConnection, optional
            Explicit stream-to-stream links (:class:`StreamConnection`). Each
            replaces the inferred link for its ``source`` in every mode, and of two
            with the same ``source`` the later wins; ``receiver=None`` marks a true
            outlet. A stream id not in the stream table, or links that route water
            in a cycle, raise ``ValueError``, and so does a link leaving its source
            above the stream's last reach (a ``Point`` ``source_location``): that
            reach would flow into two reaches, which MF6 refuses -- a partial split
            is a ``StreamDiversion`` (ledger 366). Default ``()``.
        diversions : sequence of StreamDiversion, optional
            Partial splits (:class:`StreamDiversion`), each written as an MF6
            ``DIVERSIONS`` record on its source reach plus a ``DIVERSION`` setting
            in each period its ``amount`` covers, plus the reach connection from
            the source reach to the receiving reach that MF6 diverts along --
            written once when a ``StreamConnection`` names the same pair -- so the
            two streams need no other link (ledger 364); the receiving reach
            gets ``USTRF`` 0.0. An unknown stream id or ``priority`` raises
            ``ValueError``, and so do the two placements MF6 cannot run (ledger
            366): a source reach with no downstream reach of its own -- the
            default ``source_location``, the stream's last reach, when the stream
            flows into no other; give a ``Point`` above it -- and a receiving reach
            another reach flows into alone. Default ``()``.
        connection_tolerance : float or None, optional
            Largest distance, in model length units, from a stream's last vertex to
            another stream for ``connection_mode="automatic"`` to connect them.
            Default ``None``: 5% of the square root of the median cell area, and
            never below 1e-9. Unused by the other modes.
        reverse_streams : sequence of str, optional
            Ids of streams whose line is reversed as the stream table is built, for
            lines digitized mouth-first; compared as strings. A bare string is
            iterated character by character, so give even one id as a tuple,
            ``("main",)``. Default ``()``.
        reach_layer : int, str or Mapping[str, int], optional
            Zero-based layer each reach sits in. An integer -- a Python or a numpy
            one -- is every reach's layer; a mapping is
            ``{stream_id: layer}``, where a stream missing from it raises
            ``KeyError``; any other string names a stream-table column holding each
            stream's layer (``ValueError`` when missing). Default
            ``"top_active"``: the top layer with ``idomain > 0`` in the reach's
            column of ``context.domain`` -- a column with none raises
            ``ValueError`` -- or layer 0 when the context has no domain or a
            one-dimensional one.
        width : float, str, Mapping, or sequence, optional
            Reach width, MF6 ``RWID``, in model length units; any per-reach form
            above. Default 10.0. A width that is not positive raises
            ``ValueError``.
        gradient : float, str, Mapping, or sequence, optional
            Streambed slope, MF6 ``RGRD`` (dimensionless); any per-reach form
            above. Default 0.001. A negative gradient raises ``ValueError``.
        reach_top : float, str, Mapping, sequence, or None, optional
            Streambed-top elevation, MF6 ``RTP``, in model length units; any
            per-reach form above, ending as a finite number (MF6 takes no time
            series here). Default ``None``: the top of
            each reach's layer in its cell -- from ``context.surfaces``, else
            ``grid.gdf_topbtm`` (``ValueError`` when neither is available) --
            lowered where needed so no reach is higher than the reach above it in
            the same stream.
        roughness : float, str, Mapping, or sequence, optional
            Manning's roughness coefficient, MF6 ``MAN``; any per-reach form above,
            or a time-series name. Default 0.03. A value that is not positive
            raises ``ValueError``.
        streambed_k : float, str, Mapping, or sequence, optional
            Streambed hydraulic conductivity, MF6 ``RHK``, in model length per
            model time; any per-reach form above, or a time-series name. Default
            1.0. A negative value raises ``ValueError``.
        streambed_thickness : float, str, Mapping, or sequence, optional
            Streambed thickness, MF6 ``RBTH``, in model length units; any per-reach
            form above. Default 1.0. A value that is not positive raises
            ``ValueError``.
        min_reach_length : float, optional
            Shortest piece of a stream, in model length units, that becomes a
            reach: a piece in one cell no longer than ``max(0.01,
            min_reach_length)`` is dropped. A stream clipping the corner of a cell
            yields such a sliver, MF6 divides by reach length, and a 0.005 ft reach
            crashed it while reading the package. The stream stays connected
            through its other reaches, because the topology is rebuilt from the
            surviving order. Default 0.0, which drops only pieces of 0.01 or less.
        inflow : float, str, Mapping, sequence, or None, optional
            Volumetric inflow, MF6 ``INFLOW``, in model length cubed per model
            time; any per-period form above. In one period a number or a
            time-series name enters at reach 0 ONLY -- the headwater of the first
            stream in reach order -- not at every stream's headwater;
            ``{stream_id: q}`` enters at that stream's first reach. Default
            ``None``, no inflow.
        rainfall : float, str, Mapping, sequence, or None, optional
            Rainfall rate on the reach surface, MF6 ``RAINFALL``, in model length
            per model time -- a rate per unit of streambed area, so each reach
            takes it unchanged; any per-period form above. In one period a number
            or a time-series name applies to EVERY reach, and ``{stream_id:
            rate}`` to every reach of that stream (ledger 363); a
            reach-number key names one reach and beats its stream's key whichever
            comes first. Default ``None``.
        evaporation : float, str, Mapping, sequence, or None, optional
            Evaporation rate from the reach surface, MF6 ``EVAPORATION``, in model
            length per model time; placed exactly as ``rainfall`` is -- a number on
            every reach, ``{stream_id: rate}`` on every reach of that stream.
            Default ``None``.
        runoff : float, Mapping, sequence, or None, optional
            Overland runoff, MF6 ``RUNOFF``, in model length cubed per model time.
            In one period a number is EACH stream's total and ``{stream_id:
            total}`` one stream's (a stream not in the table raises ``KeyError``);
            either is split over the stream's reaches in proportion to reach
            length, so refining the grid does not multiply it. An integer key in
            that mapping names a reach and takes its amount as given, as do
            explicit rows. A time-series name cannot be split and raises
            ``ValueError`` unless given per reach -- as rows, ``{period: [(reach,
            "series"), ...]}``, or under a reach-number key. Default ``None``.
        status : str, Mapping, sequence, or None, optional
            MF6 reach ``STATUS`` -- ``"ACTIVE"``, ``"INACTIVE"`` or ``"SIMPLE"`` --
            written as given, not checked here; any per-period form above. In one
            period a string applies to every reach, ``{stream_id: status}`` to every
            reach of that stream, placed as ``rainfall`` is. Default ``None``, no
            setting (MF6 keeps a reach active).
        name : str, optional
            Package name, default ``"sfr"``: FloPy's ``pname``, the input file
            ``<model>.<name>``, the output files ``<model>_<name>_stage.sfr`` and
            ``<model>_<name>_budget.sfr``, and the package a mover endpoint from
            :func:`sfr_connection` names.
        mover : bool, optional
            Write MF6's ``MOVER`` option, which an MVR package (``mf.mvr``) needs
            before it can move water to or from these reaches. Default ``False``.
        length_conversion : float or None, optional
            MF6 ``LENGTH_CONVERSION`` for Manning's equation (meters to model
            length units), written as given. Default ``None``: from
            ``context.length_units`` -- ``"feet"`` (the context's default) 3.28081,
            ``"meters"`` 1.0, ``"centimeters"`` 100.0, also as a singular or an
            abbreviation (``"ft"``, ``"m"``, ``"cm"``) in any case; ``None`` or
            ``"unknown"`` 1.0. Any other unit raises ``ValueError``.
        time_conversion : float or None, optional
            MF6 ``TIME_CONVERSION`` (seconds to model time units), written as
            given. Default ``None``: from ``context.time_units`` -- ``"days"`` (the
            context's default) 86400, ``"hours"`` 3600, ``"minutes"`` 60,
            ``"seconds"`` 1, ``"years"`` 31557600, also as a singular or an
            abbreviation (``"d"``, ``"h"``, ...) in any case; ``None`` or
            ``"unknown"`` 1.0. Any other unit raises ``ValueError``.
        maximum_picard_iterations : int, optional
            MF6 SFR ``MAXIMUM_PICARD_ITERATIONS``, written as given. Default 1.
        maximum_iterations : int, optional
            MF6 SFR ``MAXIMUM_ITERATIONS``, the reach-solution iteration limit,
            written as given. Default 1000.
        maximum_depth_change : float, optional
            MF6 SFR ``MAXIMUM_DEPTH_CHANGE``, in model length units, written as
            given. Default 0.01.
        **options
            Any other keyword of :class:`flopy.mf6.ModflowGwfsfr`
            (``crosssections``, ``initialstages``, ``observations``, ...), collected
            into ``SFRBuilder(options=...)`` and passed on to
            :func:`myflopy.advanced.sfr_spec`, which checks the names against
            FloPy's signature (an unknown one raises ``TypeError``, suggesting the
            closest valid one when there is one) and applies them last, so they override
            ``save_flows=True``, the ``pname``/``filename`` set from ``name`` and
            the stage/budget output file records. ``packagedata``,
            ``connectiondata`` or ``perioddata`` here raises ``TypeError``: the
            builder writes those.

        Returns
        -------
        PackageSpec

        Examples
        --------
        >>> ctx = mf.ModelContext(grid=vor, domain=idomain)
        >>> mf.sfr(context=ctx, nper=1, streams="streams.gpkg", stream_id="name",
        ...        connections=(mf.StreamConnection("trib", "main"),),
        ...        inflow={0: {"trib": 1.0e4}}, width=15.0, gradient=0.001, mover=True)
        """

        builder = SFRBuilder(
            context=context,
            nper=nper,
            streams=streams,
            stream_id=stream_id,
            from_node=from_node,
            to_node=to_node,
            connection_mode=connection_mode,
            connections=connections,
            diversions=diversions,
            connection_tolerance=connection_tolerance,
            reverse_streams=reverse_streams,
            reach_layer=reach_layer,
            width=width,
            gradient=gradient,
            reach_top=reach_top,
            roughness=roughness,
            streambed_k=streambed_k,
            streambed_thickness=streambed_thickness,
            min_reach_length=min_reach_length,
            inflow=inflow,
            rainfall=rainfall,
            evaporation=evaporation,
            runoff=runoff,
            status=status,
            name=name,
            mover=mover,
            length_conversion=length_conversion,
            time_conversion=time_conversion,
            maximum_picard_iterations=maximum_picard_iterations,
            maximum_iterations=maximum_iterations,
            maximum_depth_change=maximum_depth_change,
            options=options,
        )
        # Stash a serializable reach index so movers can be wired semantically
        # (by stream outlet/head, or by the reach in the cell holding a
        # coordinate) -- see
        # mf.sfr_connection.
        return builder.build().with_metadata(sfr_index=_sfr_mover_index(builder))

    def flopy(
        self,
        *,
        packagedata: Any,
        connectiondata: Any,
        perioddata: Any,
        nreaches: int | None = None,
        name: str = "sfr",
        mover: bool = False,
        diversions: Any = None,
        **options: Any,
    ) -> PackageSpec:
        """Return a direct FloPy-style SFR spec from prepared package data (escape hatch).

        Use this when you already have MF6 SFR records; otherwise prefer the
        geometry-driven ``mf.sfr(...)`` builder.

        Parameters
        ----------
        packagedata : list
            FloPy SFR ``packagedata`` (one row per reach).
        connectiondata : list
            FloPy SFR ``connectiondata`` (reach-to-reach topology).
        perioddata : dict
            FloPy SFR ``perioddata`` keyed by stress period.
        nreaches : int, optional
            Total reach count (inferred from ``packagedata`` when omitted).
        name : str, default "sfr"
            Package name.
        mover : bool, default False
            Enable MVR exchange with lakes/UZF.
        diversions : list, optional
            FloPy SFR ``diversions`` records.
        **options
            Extra ``flopy.mf6.ModflowGwfsfr`` options.

        Returns
        -------
        PackageSpec
        """

        return sfr_spec(
            packagedata,
            connectiondata,
            perioddata,
            nreaches=nreaches,
            name=name,
            mover=mover,
            diversions=diversions,
            **options,
        )


class _LAKPackage:
    """Package-first LAK (lake) helpers: ``mf.lak(...)`` and ``mf.lak.flopy(...)``.

    Call ``mf.lak(...)`` to build a LAK package from **lake polygon geometry** (a
    geopackage/shapefile path or a GeoDataFrame): it finds the lake cells, builds
    the lake-aquifer connections, and writes packagedata/connectiondata/perioddata.
    Use ``mf.lak.flopy(...)`` for prepared FloPy records. A lake is a head-dependent
    storage that exchanges water with the aquifer through its bed and can connect to
    streams via MVR and to outlets.

    **Every parameter of ``mf.lak(...)`` is documented on the call itself, one
    entry each**: ``help(mf.lak)`` or ``help(mf.lak.__call__)``, or hover
    ``__call__`` in ``mf.lak.__call__`` in the editor. An editor's signature help
    for ``mf.lak(`` shows only this class text. The rule most worth knowing before
    reading it, under ``connection_modes`` there: MF6 IGNORES the bed elevation
    written on a VERTICAL connection and counts from the connected cell's TOP --
    through that connection, and (without a table) in the lake's volume. In
    ``"bathy"`` the bed still picks each VERTICAL connection's layer AND sets the
    bottom of every HORIZONTAL face written on the lake cell, and MF6 exchanges
    water through those faces from the bed up, even while the stage is below the
    cell top and the lake holds no volume -- though in a TRANSIENT period it lets
    no more out of such a lake than arrives (the entry has the measurements). A
    lake the layer stack was EXCAVATED for (``LayerStack.excavate``) is the way
    out: its vertical cell's top IS the bed, and its walls go into the
    neighbours' cells -- the same entry says how, and ``unexcavated`` what is
    said about a lake that was not.

    See Also
    --------
    LakeTableBuilder : Stage / volume / area tables for ``mf.lak(tables=...)``.
    LakeOutlet : Outlets for ``mf.lak(outlets=...)``.
    LakeConnection : Hand-written connections for ``mf.lak(connection_modes=...)``.

    Examples
    --------
    >>> mf.lak(context=ctx, nper=1, lakes="pond.gpkg", lake_id_field="name",
    ...        starting_stage={"pond": 101.0}, lake_bottom=mf.Raster("bathy.tif"))
    """

    def __call__(
        self,
        *,
        context: ModelContext,
        nper: int,
        lakes: PathLike | Sequence[PathLike] | gpd.GeoDataFrame,
        layer: str | None = None,
        lake_id_field: str | None = None,
        starting_stage: Any = None,
        lake_bottom: Any = None,
        lake_top: Any = None,
        only_vertical: Any = False,
        only_layer: Any = None,
        rectangular_interior: Any = False,
        bed_leakance: Any = 1.0,
        connection_modes: str | Mapping[str, str | Sequence[LakeConnection]] = "bathy",
        unexcavated: str = "warn",
        tables: Mapping[str | int, LakeTable | LakeTableBuilder] | None = None,
        outlets: tuple[LakeOutlet, ...] = (),
        stage: Any = None,
        rainfall: Any = None,
        evaporation: Any = None,
        runoff: Any = None,
        withdrawals: Any = None,
        inflow: Any = None,
        status: Any = None,
        name: str = "lak",
        mover: bool = False,
        boundnames: bool = True,
        length_conversion: float | None = None,
        time_conversion: float | None = None,
        maximum_iterations: int = 100,
        maximum_stage_change: float = 1.0e-5,
        implicit: bool = False,
        **options: Any,
    ) -> PackageSpec:
        """Build a LAK package from lake polygons (delegates to ``LAKBuilder``).

        Lake cells are the ones each polygon footprint COVERS, by AREA -- the
        shared rule a ``.gpkg`` boundary and a head target are resolved with: a
        lake drawn along cell faces reaches only the cells inside it, and an
        overlap that is a sliver of a cell (a CRS round trip leaves these) is
        noise and is dropped, unless the whole lake is that thin. In a MULTIPART
        lake that reaches per part, so a detached part narrower than the floor
        goes with the slivers. A lake covering no cell RAISES, naming why and
        giving both extents and both CRSs. Each lake then gets its connections
        (``connection_modes``), its per-lake inputs, its per-period settings, and
        optionally a table and outlets. Everything is resolved HERE, against
        ``context``'s grid, so the package is bound to that mesh (``mf.lak`` has
        no lazy form).

        **Per-lake values.** ``starting_stage``, ``bed_leakance``,
        ``lake_bottom``, ``lake_top``, ``only_vertical``, ``only_layer`` and
        ``rectangular_interior`` take one value per lake. Where an entry below
        says "any per-lake form", that is: one number for every lake; the name
        of a column of ``lakes``; a numeric string (one number for every lake);
        ``{lake_id: value}`` naming EVERY lake (a missing lake raises; a key is
        read through ``str()``, as the ids are made, so ``101`` names lake
        ``"101"`` -- it used to miss it, and ``only_vertical={101: True}`` was
        ignored silently, ledger 318); a ``pandas.Series`` whose index names
        every lake, read by label, a label naming a lake as a key does (an
        integer-indexed Series used to be read by POSITION, so one out of lake
        order gave each lake another's value; two labels naming one lake raise
        ``ValueError``), except that for ``lake_bottom`` and ``lake_top`` an
        integer-indexed Series that ALSO holds every cell of every lake is read
        by cell, its integers being cell numbers too; or any other sequence of
        one value per lake, in lake order (``LAKBuilder.lake_ids``).

        **Per-period values.** ``stage``, ``rainfall``, ``evaporation``,
        ``runoff``, ``withdrawals``, ``inflow`` and ``status`` each become one MF6
        period setting per lake per period, over ``range(nper)``. Periods are
        ZERO-based. The forms:

        * a scalar: every lake, every period;
        * a list or tuple, or a one-dimensional ``numpy`` array: every lake, by
          position (position 0 = the first stress period; ``None`` skips one);
        * ``{period: value}`` (every key an int): every lake, in the periods
          named;
        * ``{period: {lake_id: value}}``: the lakes named, in the periods named;
        * ``{lake_id: value}``: that lake, every period;
        * ``{lake_id: list}`` or ``{lake_id: {period: value}}``: that lake, per
          period;
        * a ``pandas.Series``: the mapping of its index, so an integer index is
          periods and a lake-id index is lakes.

        A period a value does not reach writes nothing for that lake, and MF6
        KEEPS the previous period's setting (measured: rainfall given for the
        first period only was still applied in the second). A mapping whose keys
        are ALL integers is read as periods and any other by lake id, so with
        integer-like lake ids key lakes by the STRING id. What would never be
        written, or written to a guess, raises ``ValueError``: a period outside
        ``range(nper)`` -- ``rainfall={101: 0.5}`` over lake ``"101"`` wrote
        nothing at all, silently -- a list longer than ``nper``, and an
        all-integer mapping whose every key is ALSO a lake id (give periods as
        a list, or as ``{period: {lake_id: value}}``). A ``numpy`` array or
        ``pandas`` Series was written WHOLE into every period, and MF6 died
        reading it (a segmentation fault); both are read as above now, and an
        array of two or more dimensions raises (ledger 318). Each value is
        checked here as MF6 would check it one run later: a ``status`` other
        than ``ACTIVE``, ``INACTIVE`` or ``CONSTANT`` (any case), a number that
        is not finite, and a negative rainfall, evaporation, runoff, inflow or
        withdrawal all raise ``ValueError``; a string elsewhere is written as
        given, which is how MF6 names a time series. MF6 refused the status and
        each negative forcing by name (measured, all five).

        Parameters
        ----------
        context : ModelContext
            Required (``None`` raises ``AttributeError``). The geometry the lakes
            resolve against. Its ``grid`` is required (``ValueError`` when
            ``None``); the lake cells, their neighbours and the connection
            lengths and widths come from it. The layer elevations come from
            ``surfaces`` when given (a table or a
            ``LayerBuildResult``), else from ``grid.gdf_topbtm``, and one of the
            two must exist (``ValueError``). The EXCAVATION RECORDS that decide
            whether a lake is written by the excavated rule (``connection_modes``)
            come from the same place, never from the geometry: a
            ``LayerBuildResult`` given as ``surfaces`` carries its
            ``excavations``; with no ``surfaces`` the grid carries the ones
            ``LayerStack.build(..., attach=True)`` published beside
            ``gdf_topbtm``; a bare table carries none, unless it IS the grid's
            own ``gdf_topbtm``. ``domain`` is the idomain; ``None``
            means every cell is active. ``length_units``/``time_units`` set the
            conversions when ``length_conversion``/``time_conversion`` are
            ``None``.
        nper : int
            Required. The number of stress periods the per-period values are
            expanded over, from zero. Must be at least 1 (``ValueError``). Give
            the simulation's TDIS count; the two are not compared here.
        lakes : str, pathlib.Path, sequence of str or Path, or geopandas.GeoDataFrame
            Required. The lake polygons, ``Polygon`` or ``MultiPolygon`` only
            (``ValueError`` otherwise, and for a polygon with no area). A path or
            a list of paths is read at ``layer`` (each file's first layer by
            default), and the files are stacked in order; files in different
            CRSs raise ``ValueError: Cannot determine common CRS``. The lakes are
            reprojected to the grid's CRS when both declare one; lakes carrying
            no CRS are read as the grid's. Lake order (MF6's lake numbers, from
            zero) is row order.
        layer : str, optional
            Default ``None``. Which layer of each ``lakes`` FILE to read, by NAME
            (a numpy string is converted with ``str``); ``None`` reads the first,
            which is all ``mf.lak`` could read before this parameter existed (a
            multi-layer GeoPackage reads its first layer behind pyogrio's generic
            warning). A name the file does not hold raises ``ValueError`` listing
            the layers present; so does a ``layer`` with a GeoDataFrame, which has
            none to choose. Read by the shared footprint reader
            (``myflopy.modflow.mf6.grid.geometry.read_footprints``). This is the
            VECTOR layer, not a model layer -- that is ``only_layer`` -- and so an
            integer is refused (``TypeError``) rather than read as a position in
            the file: ``layer=0`` reads as a model layer, and as the file's first
            layer it would have been the default, changing nothing.
            ``geopandas.list_layers(path)`` gives the names.
        lake_id_field : str, optional
            Default ``None``. The attribute naming each lake; its values are
            converted with ``str()`` to become the lake ids (``101`` becomes
            ``"101"``, a float ``101.0`` becomes ``"101.0"``). A name that is not
            a column raises ``ValueError``, and so do repeated ids. ``None`` names
            each lake by its FILE stem when ``lakes`` is paths and every stem is
            unique (one feature per file), and otherwise ``"lake_0"``,
            ``"lake_1"``, ... in row order. Every ``{lake_id: ...}`` mapping,
            ``tables`` key and :class:`~myflopy.modflow.mf6.lakes.LakeOutlet`
            uses these ids.
        starting_stage : float, str, dict of str to float, pandas.Series or sequence of float
            Required (``ValueError: starting_stage is required``). Each lake's
            stage at the start of the simulation (length), in any per-lake form.
            Refused below the lowest ``lake_bottom`` over a generated lake's
            cells, and outside a lake's table (``[minimum_stage,
            maximum_stage]``).
        lake_bottom : float, str, pathlib.Path, Surface, pandas.Series or dict, optional
            Default ``None``. The lakebed elevation (length). For an EXCAVATED
            lake (``connection_modes``) ``None`` takes the excavation's own bed --
            the bed is declared once, on ``LayerStack.excavate`` -- and a
            ``{lake_id: ...}`` mapping that leaves the lake out does the same. A
            value GIVEN for an excavated lake must be that bed: it is sampled as
            below and compared cell by cell, before any capping, and a difference
            beyond ``max(1e-6 * |bed|, 1e-6)`` at any cell raises ``ValueError``
            naming the worst (the stack was cut to another bottom). For every
            other lake it is required by the two generated modes (``ValueError``),
            and it is not read for a lake with explicit connections.
            One value per lake, in any per-lake form, is a flat bed -- what
            ``"rectangular"`` lakes take. ``"bathy"`` wants one per CELL:

            * a :class:`~myflopy.surfaces.Surface` -- ``mf.Raster("bathy.tif")``,
              ``mf.Contours(...)``, or algebra on one (``mf.Raster(dem).below(10)``);
            * a raster path, as a ``Path`` or a ``str`` (a ``str`` naming a column
              of ``lakes`` is read as the column);
            * a cell-indexed ``pd.Series`` (its index is cell numbers, not lake ids).

            Bare, it serves every lake; ``{lake_id: ...}`` gives each its own. A
            Surface or raster is sampled onto the grid once, area-weighted, and
            each lake reads its cells from that. Only the lake's own cells need a
            value: the neighbours ``"bathy"`` compares against are land, and sit at
            the model top wherever the source has none -- a raster clipped to the
            lake, or a Series of the lake's cells. A per-cell bed above its cell's
            model top is capped there, with a warning. A lake cell a raster has no
            data for (one the polygon clips at its edge) is set to the top, with a
            warning; a lake cell missing from a Series raises. What the bed DOES
            in MF6 is under ``connection_modes``: it picks each VERTICAL
            connection's layer, through which MF6 then counts from that cell's
            top, and it is the bottom of each HORIZONTAL face ``"bathy"`` writes on
            the lake cell, through which MF6 exchanges water from the bed up (out
            of the lake, in a transient period, no more than arrives).
        lake_top : float, str, pathlib.Path, Surface, pandas.Series or dict, optional
            Default ``None``. The facility top (length), **required for
            ``"rectangular"``** lakes (``ValueError``) and refused at or below
            ``lake_bottom`` there. It is the top of a rectangular lake's
            HORIZONTAL faces (and the ``telev`` written on its VERTICAL ones,
            which MF6 ignores). Takes the same forms as ``lake_bottom``, but is
            read at the lake's cells only and never filled or capped. Ignored by
            ``"bathy"``.
        only_vertical : bool or dict of str to bool, default False
            ``True``, for every lake, or ``{lake_id: True}`` (a lake left out is
            ``False``) drops every HORIZONTAL connection of a generated lake, in
            either mode: a sealed vault whose only exchange is through its floor.
        only_layer : int, dict of str to int, or None, default None
            A zero-based layer index, for every lake or ``{lake_id: layer}``
            (a lake left out is unrestricted), to which a generated lake's
            connections are confined. ``None`` means every layer a face spans. In
            either generated mode a cell whose bed is not in that layer raises
            ``ValueError: ... does not intersect active cell`` (a ``"rectangular"``
            lake used to get NO vertical connection there, silently; ledger 318).
            On an EXCAVATED lake it must name the top active layer under every
            lake cell (``ValueError`` otherwise), and its walls are confined to
            that layer.
        rectangular_interior : bool or dict of str to bool, default False
            ``"rectangular"`` lakes only. ``True`` (or ``{lake_id: True}``) also
            connects each lake cell toward its neighbours INSIDE the lake,
            modeling a permeable-fill basin (e.g. clean gravel) that is
            hydraulically continuous with the aquifer through every face.
            ``False`` = edges only (a lined trench or vault). It is a physical
            choice about the basin, not a solver setting, although all-faces also
            spreads the lake-stage coupling, which steadies a small basin that
            drains empty. Refused (``ValueError``) on an EXCAVATED lake: its
            floor is open water, with no interior face to conduct through -- a
            gravel-filled basin is aquifer, and is not excavated.
        bed_leakance : float, str, dict, pandas.Series or sequence of float, default 1.0
            The lakebed leakance (1/T: bed hydraulic conductivity divided by bed
            thickness) written on every generated connection, in any per-lake
            form. Two more forms reach one connection TYPE: ``{lake_id:
            {"vertical": v, "horizontal": h}}``, and a key ``(lake_id,
            "vertical")`` or ``(lake_id, "horizontal")``, which wins over the
            lake's own entry. The resolved value must be greater than zero
            (``ValueError``). An explicit
            :class:`~myflopy.modflow.mf6.lakes.LakeConnection` that sets its own
            ``bed_leakance`` wins over all of these.
        connection_modes : str or dict of str to (str or sequence of LakeConnection), default "bathy"
            How each lake's connections are made. A string, in any case, serves
            every lake; ``{lake_id: ...}`` gives each its own and must name every
            lake (``ValueError``). A lake cell can carry a VERTICAL connection and
            HORIZONTAL ones together. What each mode WRITES:

            * ``"bathy"`` (alias ``"automatic"``) -- a natural lake over a bed
              that varies by cell. Per lake cell, ONE VERTICAL connection in the
              highest active layer whose top and bottom bracket that cell's bed,
              written with ``belev = telev = bed`` (a bed in no active layer
              raises ``ValueError``). Then, toward each face neighbour whose bed
              is HIGHER -- another lake cell, or land, whose bed is the per-cell
              source's value there or else the model top -- one HORIZONTAL
              connection per active layer the step spans, written on the LAKE
              CELL ITSELF (the connection's cell id is the lake cell, not the
              neighbour), with ``belev = max(bed, layer bottom)``, ``telev =
              min(neighbour's bed, layer top)``, ``CONNLEN`` the distance from
              the LAKE cell's DISV centre to the line of the shared face
              (``grid.cl12``, MF6's own DISV geometry; until ledger 328 it was the
              full centre-to-centre spacing, about twice that) and ``CONNWIDTH``
              the shared face's length. A flat bed has no steps, so a lake of more
              than one cell gets NO horizontal connection, with a warning.
            * ``"rectangular"`` -- a box: trench, vault, lined basin. Needs a flat
              ``lake_bottom`` and a ``lake_top``. Per lake cell, ONE VERTICAL
              connection in the active layer holding ``lake_bottom``, written
              ``belev = lake_bottom``, ``telev = lake_top``; a cell where
              ``lake_bottom`` falls in no active layer raises ``ValueError``, as in
              ``"bathy"`` (it used to get no VERTICAL connection, silently; ledger
              318). Then,
              toward each neighbour OUTSIDE the lake (every neighbour with
              ``rectangular_interior=True``), one HORIZONTAL connection per active
              layer of the lake cell between ``lake_bottom`` and ``lake_top``,
              again on the lake cell itself, with the same ``CONNLEN`` and
              ``CONNWIDTH``.
            * a sequence of :class:`~myflopy.modflow.mf6.lakes.LakeConnection`
              (dict form only) -- written exactly as given; ``lake_bottom``,
              ``lake_top``, ``only_vertical``, ``only_layer`` and
              ``rectangular_interior`` do not touch that lake.

            **An EXCAVATED lake is written by its own rule, in either generated
            mode.** A lake is excavated when an excavation record reaches
            ``mf.lak`` (see ``context``) and the lake's cells are EXACTLY that
            excavation's cells -- ``LayerStack.excavate(within=<these lakes=
            polygons>)``, which covers cells by the same rule. Chosen from the
            record, never from the geometry (``LAKBuilder.lake_rules`` says which
            rule each lake got). Per lake cell: ONE VERTICAL connection to the top
            ACTIVE cell under it, whose top IS the bed, so MF6 counts from the
            bed; then, toward each face neighbour whose model top -- its bed, for
            a lake cell; the ground, for land; capped at ``lake_top`` for
            ``"rectangular"`` -- is above this cell's bed, one HORIZONTAL
            connection per active layer of the NEIGHBOUR the wall spans, written
            INTO the neighbour's cell (the aquifer the water meets), with ``belev
            = max(bed, layer bottom)``, ``telev = min(wall top, layer top)``,
            ``CONNLEN`` the distance from the NEIGHBOUR's DISV centre to the line
            of the shared face and ``CONNWIDTH`` the face's length -- MF6's own
            DISV geometry, the number ``grid.cl12`` holds in the neighbour's row. A
            cell whose declared bed did not reach below the ground (a rim cell)
            connects at the ground. Measured, a pit lake below the water table in
            uniform flow against the same pit as a K = 1e4 inclusion: through-flow
            1.005 of the reference, where ``"bathy"`` on the un-excavated stack
            gives 0.52 (0.39 before ledger 328). A lake covering only
            part of an excavation, or cells of two, raises ``ValueError``
            whatever ``unexcavated`` says. Declare the model with NEWTON: without
            it MF6 carries nothing through a wall whose cell is dry (the build
            warns). In EVERY mode a HORIZONTAL face thinner than 1e-6 is not
            written: MF6 would read it as the whole cell.

            **MF6 IGNORES the BELEV and TELEV written for a VERTICAL connection**
            and uses the connected cell's TOP (and top + ``surfdep``) instead (MF6
            input guide: "Any value can be specified if CLAKTYPE is VERTICAL";
            ``gwf-lak.f90``). So a bed that lies INSIDE a cell, below its top,
            does not lower the VERTICAL connection: through it MF6 exchanges
            water, and without a table stores it, only from that cell's top up.
            Measured: one VERTICAL connection under a cell whose top is 8.0, lake
            held at 12.0, bed written as 3.0, 6.0 and 0.0; MF6 reported
            ``VOLUME`` 400.0 = 100 ft2 x (12 - 8) all three times. For a
            HORIZONTAL connection MF6 DOES read both, requires them inside the
            cell, and treats ``belev == telev`` as the whole cell -- and
            ``"bathy"`` writes the lake cell's own HORIZONTAL faces from the bed
            up, so through THEM MF6 exchanges water from the bed, below the cell
            top, while the lake holds no volume there. What leaves the lake that
            way is still bounded by what it holds: in a TRANSIENT period MF6 caps
            the water leaving a lake at what is available to it (its storage at
            the start of the step plus what arrives), so a lake below its cell
            top takes water in through such a face freely but lets out only what
            arrives. Measured on the same cell, its VERTICAL connection plus a
            HORIZONTAL face from 3.0 (the bed) to 8.0, lake held at 6.0, below
            the top: ``VOLUME`` 0.0 and the VERTICAL connection 0.0 every time;
            the face carried 3.46 ft3/d out of the lake in a steady period, 0.0
            in a transient one, 2.0 there with an ``INFLOW`` of 2.0, and 5.15
            INTO the lake with the aquifer at 7.0. So the bed decides which LAYER
            a VERTICAL connection lands in and where each HORIZONTAL face starts.
            The elevations written never depend on the (transient) stage: MF6
            wets and dries each face itself, every time step.
        unexcavated : {"warn", "raise", "ignore"}, default "warn"
            What to do about a generated lake NO excavation record reached, when
            MF6 will not see its bed. ``"warn"`` issues a ``UserWarning`` (and a
            log line), ``"raise"`` raises ``ValueError`` with the same text, and
            ``"ignore"`` says nothing; any other value raises ``ValueError`` when
            the builder is constructed. Two findings trigger it. VERTICAL
            connections in a cell whose TOP is above the lake's bed: MF6 counts
            from that top, so each exchanges nothing until the stage reaches it
            and, without a table, the storage between the bed and the top is
            missing; the text gives how many, the worst and the median gap, and
            how many of those tops are above the starting stage and the table's
            highest stage (on a real pond, 577 of 612, worst 12.58). And a lake
            that LOOKS excavated -- more than half its cells have a zero-thickness
            layer above an active cell whose top is the bed -- but carries no
            record, so it was written without walls into its neighbours. A lake whose
            bed IS its cells' tops (a USG export's) triggers neither. Checked when the
            lakes are validated: at the ``mf.lak`` call, or at
            ``LAKBuilder.build``/``validate``. It never softens a partial match with an
            excavation, which always raises. The remedy both findings name: cut the pond
            into the layering and re-declare ``mf.lak`` -- ``stack.excavate(bottom=,
            within=)`` on a ``LayerStack``, rebuilt with ``attach=True`` or its
            ``LayerBuildResult`` passed as ``ModelContext(surfaces=...)``; or, with no
            stack (a DISV built from a loaded run's arrays), ``dug =
            result.excavate(bottom=, within=)`` on the ``LayerBuildResult`` the DISV was
            built from -- a NEW result, to build the DISV and the context from
            (``mf.disv(dug)``, ``ModelContext(surfaces=dug, domain=dug.idomain)``).
        tables : dict of str or int to (LakeTable or LakeTableBuilder), optional
            Default ``None``, the same as ``{}``: no tables. Stage / volume /
            surface-area tables, at most one per lake. A value is a
            :class:`~myflopy.modflow.mf6.lakes.LakeTable` or a
            :class:`~myflopy.modflow.mf6.lakes.LakeTableBuilder`, which is built
            here; anything else raises ``TypeError``. A key is read as a lake id is
            made, through ``str()``, so a lake whose ``lake_id_field`` value is the
            integer ``101`` may be keyed ``101`` or ``"101"``; both at once raise
            ``ValueError``. (An int key used to raise ``unknown lake '101'``,
            printing the int exactly as the id it missed; ledger 318.) A key naming
            no lake raises ``ValueError`` listing the lake ids. Without a table,
            MF6 computes a lake's volume
            and surface area from the plan areas and TOPS of the cells its
            VERTICAL connections reach, not from the bed elevations. A table's
            stages must bracket that lake's ``starting_stage`` (``ValueError``
            otherwise). Each table is written to its own file: its ``filename``,
            or else ``"<name>_<lake_id>.lak.tab"``. To build one from a
            bathymetry raster, from ``mf.Contours`` (via ``.resolve_source()``),
            as a constant-area prism or from explicit rows, see
            ``mf.LakeTableBuilder``.
        outlets : tuple or list of LakeOutlet, default ()
            Lake outlets. Anything in it that is not a
            :class:`~myflopy.modflow.mf6.lakes.LakeOutlet` raises ``ValueError``.
            Each outlet names its ``source`` lake by lake id, and its
            ``receiver`` lake by lake id or ``None`` to discharge out of the
            model. An unknown id raises ``ValueError``. Outlets are numbered by
            position, from zero. ``outlet_type`` picks the discharge law:
            ``"SPECIFIED"`` (the default) uses ``rate``; ``"WEIR"`` uses
            ``invert`` and ``width``; ``"MANNING"`` uses those plus ``roughness``
            and ``slope``. A ``rate`` (a scalar for every period, a list with
            position 0 = the first stress period, or ``{period: value}`` keyed by
            ZERO-based stress period, or an array or Series read as
            ``LakeOutlet.rate`` says) is written as a ``RATE`` period setting, and
            refused (``ValueError``) where it would never be written or is not a
            finite number. MF6 applies it only to a SPECIFIED outlet, where a NEGATIVE rate
            removes water from the lake and a POSITIVE rate ADDS water. ``()``:
            no outlets.
        stage : float, str, list, tuple, dict, numpy.ndarray or pandas.Series, optional
            Default ``None`` (nothing written). The lake stage (length), a
            per-period value. MF6 applies it ONLY to a lake whose ``status`` is
            ``"CONSTANT"``, and holds the lake at it. A number must be finite.
        rainfall : float, str, list, tuple, dict, numpy.ndarray or pandas.Series, optional
            Default ``None``. Rainfall rate on the lake surface (length per
            time), a per-period value; MF6 multiplies it by the lake's surface
            area. Must be >= 0 (``ValueError``; MF6 refuses it by name).
        evaporation : float, str, list, tuple, dict, numpy.ndarray or pandas.Series, optional
            Default ``None``. The MAXIMUM evaporation rate from the lake surface
            (length per time), a per-period value; MF6 lowers it when the lake
            cannot supply it. Must be >= 0 (``ValueError``; MF6 refuses it by
            name).
        runoff : float, str, list, tuple, dict, numpy.ndarray or pandas.Series, optional
            Default ``None``. Runoff into the lake (length cubed per time), a
            per-period value. Must be >= 0 (``ValueError``; MF6 refuses it by
            name).
        withdrawals : float, str, list, tuple, dict, numpy.ndarray or pandas.Series, optional
            Default ``None``. The MAXIMUM withdrawal from the lake (length cubed
            per time), a per-period value, written as MF6's ``WITHDRAWAL``. Must
            be >= 0 (``ValueError``; MF6 refuses it by name): it is a positive
            number that removes water.
        inflow : float, str, list, tuple, dict, numpy.ndarray or pandas.Series, optional
            Default ``None``. External inflow to the lake (length cubed per time),
            a per-period value. Must be >= 0 (``ValueError``; MF6 refuses it by
            name).
        status : str, list, tuple, dict, numpy.ndarray or pandas.Series, optional
            Default ``None`` (MF6's default, ``ACTIVE``). Each lake's status, a
            per-period value: ``"ACTIVE"``, ``"INACTIVE"`` or ``"CONSTANT"`` (the
            last holds the lake at ``stage``), in any case, written as given (MF6
            accepts any case, measured with ``"constant"``). Anything else raises
            ``ValueError`` here; MF6 refused it by name one run later (measured:
            ``Unknown LAK lak status keyword: BOGUS``).
        name : str, default "lak"
            The package name. It is the FloPy ``pname``, the input file
            ``<model>.<name>``, the stem of the four output files
            (``<model>_<name>_stage.lak``, ``_budget.lak``, ``_budget.csv``,
            ``_package_convergence.csv``), the prefix of each default table file
            (``<name>_<lake_id>.lak.tab``) and table child (``<name>_<lake_id>_table``),
            and the package a mover endpoint names (``mf.lak_connection``). Give
            two LAK packages in one model different names.
        mover : bool, default False
            Writes MF6's ``MOVER`` option, so the lakes can give water to, and
            take it from, other packages through ``mf.mvr``; address a lake with
            ``mf.lak_connection(spec, lake_id)``.
        boundnames : bool, default True
            Writes each lake's id as its MF6 boundname (and the ``BOUNDNAMES``
            option), so observations and budgets can name lakes by id. ``False``
            leaves both out.
        length_conversion : float, optional
            Default ``None``. MF6's ``LENGTH_CONVERSION``, from metres to model
            length units; MF6 uses it only for WEIR and MANNING outlet flows.
            ``None`` derives it from ``context.length_units``, read ignoring
            case and surrounding spaces: 3.28081 for ``"feet"``/``"foot"``/
            ``"ft"``, 1.0 for ``"meters"``/``"meter"``/``"m"``, 100.0 for
            ``"centimeters"``/``"centimeter"``/``"cm"``, 1.0 for ``None``,
            ``""`` or ``"unknown"``, and ``ValueError`` for any other name.
        time_conversion : float, optional
            Default ``None``. MF6's ``TIME_CONVERSION``, from seconds to model
            time units; also used only for WEIR and MANNING outlets. ``None``
            derives it from ``context.time_units``, read the same way: 86400.0
            for ``"days"``/``"day"``/``"d"``, 3600.0 for ``"hours"``/``"hour"``/
            ``"hr"``/``"h"``, 60.0 for ``"minutes"``/``"minute"``/``"min"``, 1.0
            for ``"seconds"``/``"second"``/``"sec"``/``"s"``, 31557600.0 for
            ``"years"``/``"year"``/``"yr"``, 1.0 for ``None``, ``""`` or
            ``"unknown"``, and ``ValueError`` for any other name.
        maximum_iterations : int, default 100
            MF6's ``MAXIMUM_ITERATIONS``: the most Newton-Raphson iterations
            MF6 spends on each lake's stage per outer iteration. Raise it only
            when a lake's budget error is large.
        maximum_stage_change : float, default 1e-05
            MF6's ``MAXIMUM_STAGE_CHANGE``: the lake-stage closure tolerance
            (length).
        implicit : bool or numpy.bool_, default False
            Writes MF6's LAK ``IMPLICIT`` option: each lake's stage is solved as one
            more unknown of the groundwater-flow matrix, WITH the heads, instead of by
            MF6's default substitution iteration, which solves the stage and the heads
            by turns. That can take many outer iterations, or fail to converge, when a
            lake and the aquifer are strongly coupled -- a bed that leaks freely, an
            excavated pond exchanging through its floor and walls: on a real one the
            default formulation took 81% more outer iterations and 16 failed time steps,
            and ``IMPLICIT`` 11% and none (the excavation design's runs). The converged
            answer is the same: a transient pit filled at 2,000 ft3/d reached the
            default formulation's stages to 1e-8 ft (measured). MF6's requirements, each
            checked here: the solution's linear acceleration must be BICGSTAB -- MF6
            stops on CG, so an IMS resolving to CG (``linear_acceleration="CG"``, or
            ``complexity="SIMPLE"`` or no complexity with no ``linear_acceleration``) is
            refused when the simulation is built, and the IMS ``mf.simulation`` makes
            for you picks BICGSTAB -- for a lake reached through ``mf.ref`` too, once
            the build resolves it; every lake must have a connection (``ValueError``);
            NEWTON is recommended for a lake that may become perched, and the build
            warns without it. ``maximum_iterations`` and ``maximum_stage_change`` then
            do nothing (MF6 checks the stage in IMS), and the package-convergence file
            is not asked for, since MF6 writes none. **Availability:** MF6's develop
            branch (commit ``e2dfc84c``) and the 6.8 line, in no release before; the
            executable is probed once when the run executes, and one without the option
            is WARNED about, naming it and its version -- the run goes ahead, and MF6
            stops with its own message. FloPy 3.10 has no keyword for it, so myflopy
            writes it into the LAK file after FloPy does (``ModflowGwflakImplicit``). It
            survives a run loaded by ``Run.open`` or ``mf.load_mf6_run`` and written
            again, ``ParallelModelWorkflow.split_model`` and a package artifact; a bare
            ``flopy.mf6.MFSimulation.load`` drops it silently, and ``model.config``
            shows it as the lake's ``implicit`` setting. ``True`` rides in the spec's
            options (it is saved and moves the digest); ``False`` writes and stores
            nothing. A ``numpy.bool_`` is accepted and stored as ``True``/``False``;
            anything else raises ``TypeError``.
        **options : Any
            Any other ``flopy.mf6.ModflowGwflak`` keyword, e.g. ``surfdep``,
            ``print_stage``, ``print_flows``, ``observations`` or
            ``auxiliary``. A name FloPy does not accept raises ``TypeError``
            with a did-you-mean. These override what ``mf.lak`` sets by default
            -- ``save_flows=True`` and the four output file records -- and a
            string value may contain ``{model_name}``, filled in at build.
            Rename the package with ``name=``: a ``pname`` here raises
            ``TypeError``. It used to rename only what FloPy writes, while the
            spec's name, the ``lak_index`` that ``mf.lak_connection`` hands
            ``mf.mvr``, and each table child's ``pname`` and file all kept
            ``name`` (measured: ``pname="zz"`` left all three at ``"lak"``), so a
            mover built from ``mf.lak_connection`` named ``"lak"`` where FloPy
            wrote ``"zz"`` (ledger 318). A ``filename`` here names the LAK INPUT
            file only (default ``<model>.<name>``); the outputs, the table files
            and the name a mover reads follow ``name``. What ``mf.lak`` computes
            itself is not for here: ``nlakes`` raises ``TypeError`` naming it (it
            used to replace the lake count silently, ledger 318), and
            ``noutlets``, ``ntables`` and the three data blocks raise
            ``TypeError: ... got multiple values``.

        Returns
        -------
        PackageSpec
            The LAK package spec, with one child ``laktab`` per table. Its
            metadata carries ``lak_index`` (the package name and each lake id's
            MF6 lake number, which ``mf.lak_connection`` reads), the lake ids,
            each lake's connection mode, the rule each lake was written by
            (``lake_rules``: ``"excavated"``, ``"legacy"`` or ``"explicit"``),
            what each connection is (``lake_connections``: ``{lake_id:
            [[iconn, lake_cell, "floor" | "wall" | "rim"], ...]}`` -- an
            excavated lake's wall is written into its NEIGHBOUR's cell, so the
            row alone does not say which lake cell owns it), and the grid pin
            that refuses a build on another mesh. When the model is built, the
            connections are checked against the built DISV (each cell active,
            each face inside its cell, an excavated lake's vertical cell at its
            bed), and an excavated lake WARNS about a model without NEWTON,
            walls into ``icelltype`` 0 cells, recharge or ET on its vertical
            cells with no ``rainfall``/``evaporation``, and a bed so transparent
            that the mesh sets its infiltration rate.

        Raises
        ------
        ValueError
            For every refusal named above: a missing grid or layer elevations,
            ``nper < 1``, a non-polygon or zero-area lake, a lake covering no
            cell, a missing or repeated lake id, a per-lake mapping missing a
            lake or a sequence of the wrong length, two keys naming one lake
            (``101`` and ``"101"``), an unknown connection mode, a bed in no
            active layer (either generated mode), a missing ``lake_top`` or one
            at or below the bottom (``"rectangular"``), a ``bed_leakance`` of
            zero or less, a ``starting_stage`` below the bed or outside its
            table, an unknown lake id in ``tables``, ``outlets`` or a per-period
            mapping, a non-``LakeOutlet`` outlet, and any ``LakeTableBuilder``
            refusal. For a per-period value or an outlet ``rate``: a period
            outside ``range(nper)``, a list longer than ``nper``, an all-integer
            mapping whose every key is also a lake id, an array of two or more
            dimensions, a ``status`` MF6 does not know, a number that is not
            finite, and a negative rainfall, evaporation, runoff, inflow or
            withdrawal. For an EXCAVATED lake: a lake covering only part of an
            excavation or cells of two, a ``lake_bottom`` that is not the
            excavation's bed, ``rectangular_interior=True``, and an
            ``only_layer`` that is not the top active layer. An
            ``unexcavated`` that is not ``"warn"``, ``"raise"`` or
            ``"ignore"``, and -- with ``unexcavated="raise"`` -- a lake no
            excavation reached whose bed MF6 will not see. With
            ``implicit=True``, a lake with no connection; and when the simulation is
            BUILT, a solution whose linear acceleration resolves to CG.
        TypeError
            A ``layer`` that is not a ``str`` (an integer reads as a model layer;
            that is ``only_layer``); an ``implicit`` that is not a bool; an
            ``**options`` name ``ModflowGwflak`` does not accept, or one ``mf.lak``
            already sets (``implicit`` is not one: it is a named parameter),
            ``nlakes`` or ``pname``; a ``tables`` value that is neither a
            ``LakeTable`` nor a ``LakeTableBuilder``.
        pyogrio.errors.DataSourceError
            A ``lakes`` path that does not exist or is not a vector file.

        See Also
        --------
        LakeTableBuilder : Stage / volume / area tables for ``tables=``.
        LakeOutlet : Outlets for ``outlets=``.
        LakeConnection : Hand-written connections for ``connection_modes=``.
        lak_connection : A lake as a mover endpoint.

        Examples
        --------
        A natural lake whose bed comes from a bathymetry raster, open to movers
        (a flat ``lake_bottom`` here would warn: ``"bathy"`` finds no step in a
        flat bed, so a lake of several cells would get no HORIZONTAL connection):

        >>> ctx = mf.ModelContext(grid=vor, domain=idomain, surfaces=vor.gdf_topbtm)
        >>> mf.lak(context=ctx, nper=1, lakes="lakes.gpkg", lake_id_field="name",
        ...        starting_stage={"valley_lake": 101.0},
        ...        lake_bottom={"valley_lake": mf.Raster("bathy.tif")},
        ...        bed_leakance=0.1, status={"valley_lake": ["ACTIVE"]}, mover=True)

        A natural lake whose bed comes from bathymetry contours, clipped to the pond:

        >>> bed = mf.Contours("pond_contours.gpkg", z="Elev", epsg="2927",
        ...                   region_vector="pond.gpkg")
        >>> mf.lak(context=ctx, nper=1, lakes="pond.gpkg", lake_id_field="name",
        ...        starting_stage={"pond": 101.0}, lake_bottom=bed, bed_leakance=0.1)

        An infiltration trench with its own table, receiving 500 ft3/d in the
        second stress period only:

        >>> mf.lak(context=ctx, nper=3, lakes="trench.gpkg", lake_id_field="name",
        ...        connection_modes="rectangular", starting_stage=96.0,
        ...        lake_bottom=96.0, lake_top=102.0, bed_leakance=0.5,
        ...        inflow={1: 500.0, 2: 0.0},
        ...        tables={"trench": mf.LakeTableBuilder(area=400.0, bottom=96.0,
        ...                                              top=102.0, stage_step=0.5)})
        """

        builder = LAKBuilder(
            context=context,
            nper=nper,
            lakes=lakes,
            layer=layer,
            lake_id_field=lake_id_field,
            starting_stage=starting_stage,
            lake_bottom=lake_bottom,
            lake_top=lake_top,
            only_vertical=only_vertical,
            only_layer=only_layer,
            rectangular_interior=rectangular_interior,
            bed_leakance=bed_leakance,
            connection_modes=connection_modes,
            unexcavated=unexcavated,
            tables={} if tables is None else tables,
            outlets=outlets,
            stage=stage,
            rainfall=rainfall,
            evaporation=evaporation,
            runoff=runoff,
            withdrawals=withdrawals,
            inflow=inflow,
            status=status,
            name=name,
            mover=mover,
            boundnames=boundnames,
            length_conversion=length_conversion,
            time_conversion=time_conversion,
            maximum_iterations=maximum_iterations,
            maximum_stage_change=maximum_stage_change,
            implicit=implicit,
            options=options,
        )
        # Stash lake numbers so movers can target a lake by name (mf.lak_connection).
        return builder.build().with_metadata(
            lak_index={"package": builder.name,
                       "lakes": {str(k): int(v) for k, v in builder.lake_numbers.items()}}
        )

    def flopy(
        self,
        *,
        packagedata: Any,
        connectiondata: Any,
        perioddata: Any,
        nlakes: int | None = None,
        name: str = "lak",
        noutlets: int = 0,
        ntables: int = 0,
        mover: bool = False,
        **options: Any,
    ) -> PackageSpec:
        """Return a direct FloPy-style LAK spec from prepared package data (escape hatch).

        Use this when you already have MF6 LAK records; otherwise prefer the
        geometry-driven ``mf.lak(...)`` builder. Nothing here resolves against a
        grid, and the RECORDS go to ``flopy.mf6.ModflowGwflak`` as given. What is
        compared, at this call, is each COUNT with the block it counts, a
        boundname column with ``boundnames``, and the shape of a period value --
        every case one that FloPy fails at build or MF6 refuses by name when it
        reads the file, a run later, which is where each was caught until ledger
        318. Each case below is measured
        (``docs/reproducers/lak_docstring_claims.py``, section 9). A count held
        in a form this does not count (FloPy's external-file ``{"filename":
        ...}`` form, a DataFrame) is left to FloPy and MF6.

        Parameters
        ----------
        packagedata : list of list
            Required. FloPy LAK ``packagedata``, one row per lake: ``[ifno, strt,
            nlakeconn]`` -- the ZERO-based lake number, its starting stage
            (length) and how many rows it has in ``connectiondata`` -- then any
            ``auxiliary`` values, then a boundname when ``boundnames=True`` is in
            ``**options``. A row longer than that without ``boundnames=True``
            raises ``ValueError`` here; FloPy failed the build on it with an
            ``MFDataException``. The auxiliary values are counted by NAME, in any
            shape FloPy writes the same: ``["a", "b"]``, ``[("a", "b")]``,
            ``["auxiliary", "a", "b"]`` and the recarray FloPy's own
            ``lak.auxiliary.get_data()`` returns each declare two. An
            ``auxiliary`` holding anything but strings is left to FloPy, and the
            boundname check is then skipped.
        connectiondata : list of list
            Required. FloPy LAK ``connectiondata``, one row per lake-aquifer
            connection: ``[ifno, iconn, cellid, claktype, bedleak, belev, telev,
            connlen, connwidth]``, with ``iconn`` zero-based within its lake and
            ``claktype`` one of ``"VERTICAL"``, ``"HORIZONTAL"``, ``"EMBEDDEDH"``
            or ``"EMBEDDEDV"``. MF6 ignores ``belev``/``telev`` on a VERTICAL row
            (``mf.lak``'s ``connection_modes`` entry has the rule), and
            :class:`~myflopy.modflow.mf6.lakes.LakeConnection` describes each
            field.
        perioddata : dict of int to list of list
            Required. FloPy LAK ``perioddata``, keyed by ZERO-based stress period:
            ``{period: [[number, keyword, value], ...]}``, where ``number`` is a
            zero-based lake number for ``STATUS``, ``STAGE``, ``RAINFALL``,
            ``EVAPORATION``, ``RUNOFF``, ``INFLOW``, ``WITHDRAWAL`` and
            ``AUXILIARY``, or an outlet number for ``RATE``, ``INVERT``,
            ``WIDTH``, ``SLOPE`` and ``ROUGH`` (``AUXILIARY`` takes the auxiliary
            variable's name before its value). A period left out, or given an
            empty list, keeps the previous period's settings. A record whose value
            is a ``numpy`` array or ``pandas`` Series raises ``ValueError`` here:
            FloPy writes it whole, as text, and MF6 died reading it (a
            segmentation fault, ledger 318). Give one number per record.
        nlakes : int, optional
            Default ``None``: ``len(packagedata)``. An explicit count REPLACES
            that, so it must EQUAL the number of ``packagedata`` rows: one that
            disagrees raises ``ValueError`` here, where MF6 refused it in either
            direction (measured: ``nlakes=3`` for one row, ``NO DATA SPECIFIED
            FOR LAKE 2``; ``nlakes=1`` for two rows, ``lakeno MUST BE > 0 and <=
            1``).
        name : str, default "lak"
            The package name: the FloPy ``pname``, the input file
            ``<model>.<name>`` and the stem of the four output files, as for
            ``mf.lak``. Give two LAK packages in one model different names.
        noutlets : int, default 0
            MF6's ``NOUTLETS``, which must equal the number of ``outlets`` rows
            passed through ``**options``. FloPy RAISES the written count to that
            number when fewer are declared (measured: ``0`` or ``1`` with two
            rows writes 2), so a count below the rows is accepted; FloPy never
            LOWERS it, so a count above the rows, or one with no ``outlets`` at
            all, raises ``ValueError`` here -- MF6 refused each by name (``NO DATA
            SPECIFIED FOR OUTLET 2``; ``REQUIRED OUTLETS BLOCK NOT FOUND``).
        ntables : int, default 0
            MF6's ``NTABLES``, which must equal the number of ``tables`` rows
            passed through ``**options``, under the same rule: raised to the rows
            given (``0`` with one table writes 1), never lowered, so ``ntables=2``
            with one table (MF6: ``TABLE DATA ARE SPECIFIED 1 TIMES BUT NTABLES
            IS SET TO 2``) and ``ntables=1`` with none raise ``ValueError`` here.
        mover : bool, default False
            Writes MF6's ``MOVER`` option, so this package's lakes and outlets can
            take part in ``mf.mvr``. The spec carries no ``lak_index``, so
            ``mf.lak_connection`` refuses it (``ValueError``); address an endpoint
            as ``mf.MoverConnection(name, number)``.
        **options : Any
            Any other ``flopy.mf6.ModflowGwflak`` keyword: ``outlets``,
            ``tables``, ``boundnames``, ``auxiliary``, ``surfdep``,
            ``length_conversion``, ``time_conversion``, ``observations``, ... A
            name FloPy does not accept raises ``TypeError`` with a did-you-mean.
            These override what this sets by default -- ``save_flows=True``, the
            four output file records and the ``pname``/``filename`` that ``name``
            gives. ``outlets`` and ``tables`` must agree with ``noutlets`` and
            ``ntables`` (above), and a boundname column in ``packagedata`` needs
            ``boundnames=True``; each is checked at this call.

        Returns
        -------
        PackageSpec
            The LAK package spec, holding the records exactly as given. Nothing
            is resolved against a grid and nothing is pinned to one.

        Raises
        ------
        ValueError
            Every disagreement listed in one message: an ``nlakes`` other than the
            number of ``packagedata`` rows; a ``noutlets`` or ``ntables`` above the
            ``outlets`` or ``tables`` rows; a boundname column without
            ``boundnames=True``; an array or Series as a ``perioddata`` value.
        TypeError
            A required block (``packagedata``, ``connectiondata`` or
            ``perioddata``) left out, or an ``**options`` name
            ``ModflowGwflak`` does not accept.

        See Also
        --------
        lak : The geometry-driven builder, which writes these records for you.

        Examples
        --------
        >>> mf.lak.flopy(packagedata=[[0, 101.0, 1]],
        ...              connectiondata=[[0, 0, (0, 42), "VERTICAL", 0.1, 0.0, 0.0, 0.0, 0.0]],
        ...              perioddata={0: [[0, "STATUS", "ACTIVE"]]})

        One lake with a SPECIFIED outlet taking 5 volume units per time, so
        ``noutlets`` matches the one ``outlets`` row:

        >>> mf.lak.flopy(packagedata=[[0, 101.0, 1]],
        ...              connectiondata=[[0, 0, (0, 42), "VERTICAL", 0.1, 0.0, 0.0, 0.0, 0.0]],
        ...              perioddata={0: [[0, "RATE", -5.0]]}, noutlets=1,
        ...              outlets=[[0, 0, -1, "SPECIFIED", 0.0, 0.0, 0.0, 0.0]])

        A count above its block is refused here, not by MF6 a run later:

        >>> mf.lak.flopy(packagedata=[[0, 101.0, 1]],
        ...              connectiondata=[[0, 0, (0, 42), "VERTICAL", 0.1, 0.0, 0.0, 0.0, 0.0]],
        ...              perioddata={0: []}, nlakes=3)
        Traceback (most recent call last):
        ...
        ValueError: mf.lak.flopy: the records disagree with themselves, which FloPy or MF6 ...
        """

        refusals = _flopy_lak_refusals(
            packagedata, perioddata, nlakes=nlakes, noutlets=noutlets, ntables=ntables, options=options
        )
        if refusals:
            raise ValueError(
                "mf.lak.flopy: the records disagree with themselves, which FloPy or MF6 "
                "would refuse a run later (ledger 318):\n  - " + "\n  - ".join(refusals)
            )
        return lak_spec(
            packagedata,
            connectiondata,
            perioddata,
            nlakes=nlakes,
            name=name,
            noutlets=noutlets,
            ntables=ntables,
            mover=mover,
            **options,
        )


class _HFBPackage:
    """Package-first HFB (horizontal flow barrier) helpers.

    A barrier is a thin low-permeability feature -- a fault gouge, a slurry
    cutoff wall, a till finger -- that sits ON the face between two cells rather
    than inside either of them. MODFLOW 6 applies it by reducing the aquifer's
    own intercell conductance rather than by adding flux, which is why it has no
    budget entry of its own: its only trace in the output is a correction to
    ``FLOW-JA-FACE``.

    ``hydchr`` is a hydraulic CHARACTERISTIC, units **1/T** -- the barrier's
    hydraulic conductivity divided by its thickness. A 3 ft slurry wall of
    1e-7 ft/d is ``hydchr = 1e-7 / 3``, not ``1e-7``. Zero is impermeable; a
    NEGATIVE value is a conductance multiplier rather than a barrier
    (``-1.0`` is an exact no-op).

    The entry points differ by where the geometry comes from. FloPy already
    accepts cell pairs, and neither FloPy nor anything else tells you a pair is
    wrong until MODFLOW 6 aborts mid-run -- so resolving a *line* to the faces it
    crosses, and checking the pairs against the grid first, is the point:

    - :meth:`__call__` -- explicit cell pairs;
    - :meth:`line` -- the faces a barrier trace crosses;
    - :meth:`gpkg` -- the same, from a GeoPackage line layer;
    - :meth:`enclose` -- a closed, watertight cutoff wall around a region;
    - :meth:`flopy` -- the FloPy-native form, for a schedule that changes.

    Examples
    --------
    >>> mf.hfb(pairs=[(197, 239), (239, 240)], hydchr=1e-4, context=ctx)
    >>> mf.hfb.gpkg("faults.gpkg", context=ctx, hydchr="hydchr")
    >>> mf.hfb.enclose(pond_footprint, context=ctx, hydchr=1e-8)
    """

    def __call__(
        self,
        *,
        pairs: Sequence[tuple[Any, Any]],
        hydchr: float | Sequence[float],
        layers: int | Sequence[int] | None = None,
        context: ModelContext | None = None,
        first_period: int = 0,
        name: str = "hfb",
        print_input: bool = False,
        **options: Any,
    ) -> PackageSpec:
        """Return an HFB package spec from explicit cell pairs.

        Parameters
        ----------
        pairs : sequence of pairs
            Either ``(cell_a, cell_b)`` zero-based plan-cell integer pairs on a
            DISV mesh -- written as ``(layer, cell)`` into every layer named by
            ``layers`` -- or pairs of full cellids, in which case ``layers`` must
            be ``None`` (``ValueError`` otherwise): ``((lay, cell), (lay, cell))``
            on DISV, ``((lay, row, col), (lay, row, col))`` on DIS and
            ``((node,), (node,))`` on DISU. A mixture of the two forms, or of
            cellid lengths, raises ``ValueError``, and so does a cellid other than
            ``(layer, cell)`` when ``context`` carries a grid. The build refuses
            cellids whose shape is not the model's discretization's, naming the
            spelling. An empty sequence gives a package with no barriers. The
            discretization (any type, any name) must be declared before the
            barrier; the build refuses it otherwise.
        hydchr : float or sequence of float
            Hydraulic characteristic, 1/T (the barrier's K over its thickness;
            ``0`` is impermeable, a negative value is a conductance multiplier).
            One value for every pair, or one per pair -- a length that matches
            neither raises ``ValueError``.
        layers : int or sequence of int, optional
            Zero-based layers each plan pair is written into. ``None`` means layer
            0 only; pass a sequence for a barrier that cuts several layers.
        context : ModelContext, optional
            Carries the grid. Supplied, every pair is validated against the grid
            connectivity BEFORE the run -- a self-pair or two cells that are not
            connected raises ``ValueError``, and a repeated face is kept once
            (MF6 would apply it twice) -- and the records are pinned to that
            mesh. Omitted, no validation is possible and a bad pair surfaces only
            when MODFLOW 6 aborts.
        first_period : int, optional
            Zero-based period the barriers first appear in. MODFLOW 6 reuses the
            set in every later period; use :meth:`flopy` for a schedule that
            changes. Default ``0``.
        name : str, optional
            The spec's name in the model's package list. It is not passed to
            FloPy (put ``pname=`` in ``**options`` for that). Default ``"hfb"``.
        print_input : bool, optional
            Echo the barrier list to the listing file (MF6 ``PRINT_INPUT``);
            written only when ``True``. Default ``False``.
        **options
            Extra ``flopy.mf6.ModflowGwfhfb`` options, passed at build time; a name
            that class does not accept raises ``TypeError`` at the call.

        Returns
        -------
        PackageSpec
        """

        records = _hfb_records(
            pairs, hydchr=hydchr, layers=layers, context=context, verb="mf.hfb"
        )
        return _hfb_spec_from_records(
            records, name=name, first_period=first_period, print_input=print_input,
            context=context, **options
        )

    def line(
        self,
        barrier: Any,
        *,
        context: ModelContext,
        hydchr: float | Sequence[float],
        layers: int | Sequence[int] | None = None,
        first_period: int = 0,
        name: str = "hfb",
        print_input: bool = False,
        **options: Any,
    ) -> PackageSpec:
        """Return an HFB spec from a barrier trace, resolved to the connections it cuts.

        Parameters
        ----------
        barrier : LineString, MultiLineString, GeoSeries or GeoDataFrame
            The trace, reprojected to the grid's CRS when it declares one. It cuts
            every connection whose centre-to-centre segment crosses it; a trace
            that cuts none gives an empty package and a logged warning.
        context : ModelContext
            Required -- the grid is what a line is resolved against; a context
            with no grid raises ``ValueError``. The resolved pairs are validated
            and pinned as for :meth:`__call__`.
        hydchr : float or sequence of float
            Hydraulic characteristic, 1/T: one value, or one per cut connection
            (in the order the trace resolves them).
        layers : int or sequence of int, optional
            Zero-based layers the barrier cuts. ``None`` means layer 0 only.
        first_period : int, optional
            Zero-based period the barriers first appear in; they persist after
            it. Default ``0``.
        name : str, optional
            The spec's name in the model's package list (not the FloPy
            ``pname``). Default ``"hfb"``.
        print_input : bool, optional
            Echo the barrier list to the listing file. Default ``False``.
        **options
            Extra ``flopy.mf6.ModflowGwfhfb`` options, as for :meth:`__call__`.

        Returns
        -------
        PackageSpec
        """

        grid = _hfb_grid(context, verb="mf.hfb.line")
        faces = grid.barrier_faces(barrier)
        if not faces:
            logger.warning("mf.hfb.line: the barrier cuts no connection; HFB will be empty")
        return self(
            pairs=faces,
            hydchr=hydchr,
            layers=layers,
            context=context,
            first_period=first_period,
            name=name,
            print_input=print_input,
            **options,
        )

    def gpkg(
        self,
        path: PathLike,
        *,
        context: ModelContext,
        hydchr: RowValue = "hydchr",
        layer: str | None = None,
        layers: int | Sequence[int] | None = None,
        first_period: int = 0,
        name: str = "hfb",
        print_input: bool = False,
        **options: Any,
    ) -> PackageSpec:
        """Return an HFB spec from a GeoPackage layer of barrier traces.

        Parameters
        ----------
        path : Path or str
            The GeoPackage (any file ``geopandas.read_file`` opens) of barrier
            traces, reprojected to the grid's CRS when both declare one. Each
            feature is resolved like :meth:`line`; if no feature cuts a
            connection the package is empty and a warning is logged.
        context : ModelContext
            Carries the grid the traces resolve against; required (``ValueError``
            when it carries no grid). Pairs are validated and pinned as for
            :meth:`__call__`.
        hydchr : str or float, optional
            An attribute column giving each feature's hydraulic characteristic
            (1/T), applied to every connection that feature cuts, or one number for
            all of them. Default ``"hydchr"``.
        layer : str, optional
            GeoPackage table name; the first (only) table by default.
        layers : int or sequence of int, optional
            Zero-based model layers the barriers cut. ``None`` means layer 0 only.
        first_period : int, optional
            Zero-based period the barriers first appear in; they persist after
            it. Default ``0``.
        name : str, optional
            The spec's name in the model's package list (not the FloPy
            ``pname``). Default ``"hfb"``.
        print_input : bool, optional
            Echo the barrier list to the listing file. Default ``False``.
        **options
            Extra ``flopy.mf6.ModflowGwfhfb`` options, as for :meth:`__call__`.

        Returns
        -------
        PackageSpec
        """

        grid = _hfb_grid(context, verb="mf.hfb.gpkg")
        frame = gpd.read_file(path, layer=layer)
        if grid.crs is not None and frame.crs is not None and frame.crs != grid.crs:
            frame = frame.to_crs(grid.crs)

        pairs: list[tuple[int, int]] = []
        values: list[float] = []
        for _, row in frame.iterrows():
            crossed = grid.barrier_faces(row.geometry)
            value = row[hydchr] if isinstance(hydchr, str) else hydchr
            pairs.extend(crossed)
            values.extend([float(value)] * len(crossed))
        if not pairs:
            logger.warning("mf.hfb.gpkg: no feature in %s cuts a connection", path)

        return self(
            pairs=pairs,
            hydchr=values,
            layers=layers,
            context=context,
            first_period=first_period,
            name=name,
            print_input=print_input,
            **options,
        )

    def enclose(
        self,
        polygon: Any,
        *,
        context: ModelContext,
        hydchr: float,
        layers: int | Sequence[int] | None = None,
        first_period: int = 0,
        name: str = "hfb",
        print_input: bool = False,
        **options: Any,
    ) -> PackageSpec:
        """Return an HFB spec sealing the cells inside ``polygon`` from the rest.

        The barrier is the CUT between the enclosed cell set -- cells whose
        CENTROID is inside -- and its neighbours. :meth:`line` on the same ring
        cuts nearly the same connections (82 against 83 on the canonical ring,
        since ledger 216), but only this refuses a cut that leaks. Before 216,
        ``.line`` took the faces the ring crossed: 53 of them, sealing nothing.

        Parameters
        ----------
        polygon : Polygon, MultiPolygon, GeoSeries or GeoDataFrame
            The region to seal, in the grid's coordinates (a GeoSeries/GeoDataFrame
            is unioned, and reprojected to the grid's CRS when both declare one). A
            cell is inside when its CENTROID is.
        context : ModelContext
            Carries the grid the cut is computed on; required. Pairs are pinned to
            it as for :meth:`__call__`.
        hydchr : float
            Hydraulic characteristic of the wall, 1/T -- the wall's K over its
            thickness, one value for every face of the cut.
        layers : int or sequence of int, optional
            Zero-based layers the wall cuts. ``None`` means layer 0 only.
        first_period : int, optional
            Zero-based period the wall first appears in; it persists after it.
            Default ``0``.
        name : str, optional
            The spec's name in the model's package list (not the FloPy
            ``pname``). Default ``"hfb"``.
        print_input : bool, optional
            Echo the barrier list to the listing file. Default ``False``.
        **options
            Extra ``flopy.mf6.ModflowGwfhfb`` options, as for :meth:`__call__`.

        Returns
        -------
        PackageSpec

        Raises
        ------
        ValueError
            If the resulting cut does not in fact disconnect the interior. A
            leaking cutoff wall is always a bug, so this is an error rather than
            a warning -- unlike :meth:`line`, where an open-ENDED barrier is normal.
            Also raised when the cut is empty -- no cell centroid inside
            ``polygon``, or every one -- or when ``context`` carries no grid.
        """

        from myflopy.modflow.mf6.grid.barriers import is_watertight

        grid = _hfb_grid(context, verb="mf.hfb.enclose")
        faces, interior = grid.enclosed_faces(polygon)
        if not faces:
            raise ValueError(
                "mf.hfb.enclose: the polygon encloses no cell, so there is nothing to seal. "
                "A cell counts as inside when its CENTROID is inside the polygon."
            )
        if not is_watertight(grid, faces, interior[0]):
            raise ValueError(
                f"mf.hfb.enclose: the {len(faces)} resulting faces do not disconnect the "
                f"{len(interior)} enclosed cell(s) from the rest of the grid."
            )
        return self(
            pairs=faces,
            hydchr=hydchr,
            layers=layers,
            context=context,
            first_period=first_period,
            name=name,
            print_input=print_input,
            **options,
        )

    def flopy(
        self,
        *,
        stress_period_data: Any,
        name: str = "hfb",
        print_input: bool = False,
        **options: Any,
    ) -> PackageSpec:
        """Return a direct FloPy-style HFB spec.

        ``stress_period_data`` is ``{period: [[cellid1, cellid2, hydchr], ...]}``.
        An ABSENT period reuses the previous set; an EMPTY one (``{4: []}``)
        removes every barrier.

        Nothing is validated here -- this is the escape hatch, and the pairs go
        to MODFLOW 6 as given. The build checks only that the discretization (any
        type, any name) was declared first and that list-of-records cellids have
        its shape -- both of which FloPy otherwise reports only at write.

        Parameters
        ----------
        stress_period_data : dict
            FloPy's ``{period: [[cellid1, cellid2, hydchr], ...]}``, zero-based
            periods and cellids in the discretization's shape -- ``(layer, cell)``
            on DISV, ``(layer, row, col)`` on DIS, ``(node,)`` on DISU -- and
            ``hydchr`` in 1/T. Handed to ``ModflowGwfhfb`` unchanged and not pinned
            to any mesh.
        name : str, optional
            The spec's name in the model's package list (not the FloPy
            ``pname``). Default ``"hfb"``.
        print_input : bool, optional
            Echo the barrier list to the listing file; written only when ``True``.
            Default ``False``.
        **options
            Extra ``flopy.mf6.ModflowGwfhfb`` options, passed at build time; a name
            that class does not accept raises ``TypeError`` at the call.

        Returns
        -------
        PackageSpec
        """

        values = {"stress_period_data": stress_period_data, **options}
        if print_input:
            values["print_input"] = True
        return PackageSpec(name, _HFB_BUILDER, values)


hfb = _HFBPackage()


def _hfb_grid(context: ModelContext | None, *, verb: str):
    """The grid off a context, or a message naming what is missing."""

    grid = getattr(context, "grid", None)
    if grid is None:
        raise ValueError(
            f"{verb} needs the grid to resolve geometry: pass context=ModelContext(grid=vor)."
        )
    return grid


def _hfb_records(
    pairs: Sequence[tuple[Any, Any]],
    *,
    hydchr: float | Sequence[float],
    layers: int | Sequence[int] | None,
    context: ModelContext | None,
    verb: str,
) -> list[list[Any]]:
    """Expand plan pairs across layers, validate them, and attach ``hydchr``."""

    pairs = list(pairs)
    if not pairs:
        return []

    # A pair is QUALIFIED when both ends are full cellids -- `(layer, cell)` on a DISV,
    # `(layer, row, column)` on a DIS, `(node,)` on a DISU -- and PLAN when both are
    # integers. Only the length 2 was recognised before ledger 378, so a DIS pair
    # died in `int(tuple)`. The cellid's shape against the discretization is checked at
    # build, where the discretization is known (`_build_hfb`).
    def is_cellid(value) -> bool:
        return isinstance(value, (tuple, list))

    qualified = all(is_cellid(a) and is_cellid(b) for a, b in pairs)
    if not qualified and any(is_cellid(a) or is_cellid(b) for a, b in pairs):
        raise ValueError(
            f"{verb}: every pair must be two plan-cell integers `(cell_a, cell_b)` or two "
            "full cellids `((layer, cell), (layer, cell))` -- not a mixture."
        )
    if qualified:
        lengths = {len(cellid) for pair in pairs for cellid in pair}
        if len(lengths) != 1:
            raise ValueError(
                f"{verb}: the pairs mix cellids of {sorted(lengths)} integers; one "
                "discretization writes one shape."
            )
        if layers is not None:
            raise ValueError(
                f"{verb}: pairs are already full cellids, so layers= would be ambiguous."
            )
        if getattr(context, "grid", None) is not None and lengths != {2}:
            raise ValueError(
                f"{verb}: the context's grid is a DISV mesh, whose cellid is "
                f"(layer, cell), but the pairs hold {lengths.pop()} integers each. Give "
                "(layer, cell) pairs, or drop context= for a DIS or DISU model."
            )

    if isinstance(hydchr, (int, float)):
        values = [float(hydchr)] * len(pairs)
    else:
        values = [float(v) for v in hydchr]
        if len(values) != len(pairs):
            raise ValueError(
                f"{verb}: hydchr has {len(values)} value(s) for {len(pairs)} pair(s); "
                "pass one value or one per pair."
            )

    if qualified:
        expanded = [(tuple(a), tuple(b)) for a, b in pairs]
        expanded_values = values
    else:
        layer_list = [0] if layers is None else ([layers] if isinstance(layers, int) else list(layers))
        expanded, expanded_values = [], []
        for (cell_a, cell_b), value in zip(pairs, values, strict=True):
            for layer in layer_list:
                expanded.append(((int(layer), int(cell_a)), (int(layer), int(cell_b))))
                expanded_values.append(value)

    grid = getattr(context, "grid", None)
    if grid is None:
        logger.debug(
            "%s: no grid on the context, so barrier pairs are unvalidated; "
            "MODFLOW 6 will abort at run time on any pair that is not a real connection",
            verb,
        )
    else:
        from myflopy.modflow.mf6.grid.barriers import validate_barrier_pairs

        keep = validate_barrier_pairs(grid, expanded)
        if len(keep) != len(expanded):
            kept = {(tuple(a), tuple(b)) for a, b in keep}
            paired = [
                (pair, value)
                for pair, value in zip(expanded, expanded_values, strict=True)
                if (tuple(pair[0]), tuple(pair[1])) in kept
            ]
            seen: set[tuple] = set()
            expanded, expanded_values = [], []
            for pair, value in paired:
                key = tuple(sorted([tuple(pair[0]), tuple(pair[1])]))
                if key in seen:
                    continue
                seen.add(key)
                expanded.append(pair)
                expanded_values.append(value)

    return [[list(a), list(b), value] for (a, b), value in zip(expanded, expanded_values, strict=True)]


def _hfb_spec_from_records(
    records: list[list[Any]],
    *,
    name: str,
    first_period: int,
    print_input: bool,
    context: ModelContext | None = None,
    **options: Any,
) -> PackageSpec:
    """Wrap resolved barrier records in a ``PackageSpec``.

    ``maxhfb`` is deliberately not exposed: FloPy computes it from the data at
    write time, and a hand-set value that disagrees is the same defect class as
    the ``maxbound`` bug the list-BC builder records.

    The records are PINNED to the mesh they were resolved against when a context
    supplies one. A barrier's records are cell PAIRS, so the range check sees
    nothing to check, and every pair resolved on a coarse mesh is still in range
    on a finer one -- measured, a barrier resolved on 28 cells BUILT CLEAN on 92
    (ledger 267). ``.line``, ``.gpkg`` and ``.enclose`` all funnel through the
    single caller above, so pinning there covers every resolved route; a
    hand-written pair list with no context stays unpinned, like any other cellid
    the caller typed.
    """

    values = {"stress_period_data": {first_period: records}, **options}
    if print_input:
        values["print_input"] = True
    spec = PackageSpec(name, _HFB_BUILDER, values)
    return (
        spec.with_metadata(**pin_metadata(context, certifies=("stress_period_data",)))
        if context is not None
        else spec
    )



def _refuse_per_feature_values(values: dict[str, Any], *, verb: str) -> None:
    """Refuse a ``Spread`` or ``CellSurfaceOffset`` on an ARRAY-form BC.

    Both describe how to resolve a value PER FEATURE -- ``Spread`` splits an
    extensive field across the cells one feature covers, ``CellSurfaceOffset`` reads
    a surface under each cell a feature lands on. The array form has no features:
    it takes one array or scalar per period. So neither can act, and a parameter
    that cannot act is the same defect as a hidden one.

    Measured before this: all four shapes were ACCEPTED, the live object reached
    ``options``, and the package died inside FloPy as an ``MFDataException`` naming
    a data element rather than the argument. Worse, ``Project.prepare_run`` writes
    the manifest (``save_manifest``) BEFORE it builds, and ``build=False`` is a
    supported call -- so a run could record a digest for a package holding one, and
    ledger 274's "no stored digest can move" would have had a hole in it. Refusing
    here is what makes that structural rather than merely unobserved.

    **It is handed the MERGED values, and that is the point** (ledger 283). It used to
    see only ``fields`` -- ``{"recharge"}`` for RCH, ``{"rate", "depth", "surface"}``
    for EVT -- so ``irch=``, ``auxmultname=``, ``observations=`` and every other
    keyword sailed past the check whose docstring says no digest can hold one.
    Measured: ``recharge=Spread(...)`` refused, ``irch=Spread(...)`` accepted AND
    digested. Called after the index is resolved so it sees that too.
    """

    offenders = {
        key: type(value).__name__
        for key, value in values.items()
        if isinstance(value, (Spread, CellSurfaceOffset))
    }
    if not offenders:
        return
    named = ", ".join(f"{key}={kind}(...)" for key, kind in sorted(offenders.items()))
    raise TypeError(
        f"{verb}: {named} describes a per-FEATURE resolution, and the array form has "
        "no features -- it takes one array or scalar per period. Use the list form "
        f"({verb.rsplit('.', 1)[0]}.gpkg(...)) where each feature resolves its own "
        "value, or pass a plain array here."
    )


#: What an omitted `mf.evt.array` surface is, before a grid exists to resolve it against.
#: A string is unambiguous here: the argument otherwise holds numbers.
_MODEL_TOP = "model_top"

#: The VALUE fields of each array form, so a lazy rebuild can split them back out of the
#: options it was handed.
_ARRAY_FIELDS = {"irch": ("recharge",), "ievt": ("rate", "depth", "surface")}


def _array_bc_result(
    constructor,
    *,
    name: str,
    index_name: str,
    verb: str,
    fields: dict[str, Any],
    context: ModelContext | None,
    index: Any,
    **rest: Any,
):
    """One array-form spec, EAGER or LAZY -- the single tail of both `.array` helpers.

    **The trigger is narrower than the list form's, and deliberately so.** A lazy `.gpkg`
    and a lazy `mf.rch(...)` can key on `context is None` because every such declaration
    RAISED before, so no stored run holds one. That is not true here:
    `mf.rch.array(recharge=..., irch=0 | None | <array>)` with no context is legal,
    working and digestible TODAY. Relocating its payload into ``inputs`` would change
    `package_digest`, and `project.runs.drift()` would report untouched packages as
    EDITED the first time anyone re-ran their declaring script -- the exact failure
    ledger 250 kept the pin out of the digest to avoid. So the branch is taken only for
    the parts that genuinely need a grid, which are exactly the parts that raised:
    ``irch``/``ievt`` left at ``"top_active"``, and an omitted EVT ``surface``
    (ledger 285).
    """

    deferred = {}
    if context is None:
        if isinstance(index, str):
            deferred[index_name] = index
        if isinstance(fields.get("surface"), str):
            deferred["surface"] = fields["surface"]
    if deferred:
        options = {
            key: value
            for key, value in {**fields, **rest}.items()
            if key not in deferred and value is not None
        }
        if index is not None and index_name not in deferred:
            # A TYPED index rides beside a deferred EVT surface. It was in neither
            # `fields` nor `rest`, so it was dropped: `ievt=1` wrote no IEVT at all and MF6
            # took the ET from layer 1 (ledger 324). `_build_lazy_array` pops it back out.
            options[index_name] = index
        # The lazy branch returns before `_array_bc_values` runs, so the declaration-time
        # checks that do not need a grid have to be repeated here -- the rule
        # `_lazy_bc_spec` already follows for a missing file and `nper=0`. Without it
        # `mf.rch.array(recharge=mf.Spread(...))` was ACCEPTED at the call and refused at
        # BUILD, quietly undoing ledger 283 for every lazy declaration.
        _refuse_per_feature_values({**options, **deferred}, verb=verb)
        return PackageSpec(
            name,
            _lazy_array_factory(index_name),
            options,
            requires=(),
            inputs=deferred,
        )

    if isinstance(fields.get("surface"), str):
        fields = {**fields, "surface": _model_top_surface(context, verb=verb)}
    values = _array_bc_values(
        fields, context=context, index=index, index_name=index_name, verb=verb, **rest
    )
    return _array_bc_spec(constructor, name=name, options=values, context=context)


def _build_lazy_array(constructor, index_name, model, *, _spec, **kwargs):
    """Resolve an array-form RCH/EVT declaration against the model's grid at build time.

    It re-enters ``_array_bc_values``/``_array_bc_spec`` with a resolved context, so the
    result is the eager one by construction rather than by a second implementation that
    has to be kept in step.
    """

    verb = f"mf.{'rch' if index_name == 'irch' else 'evt'}.array"
    context = _lazy_resolution_context(model, _spec.name, needs_grid=False)
    kwargs = dict(kwargs)
    index = kwargs.pop(index_name, None)
    fields = {key: kwargs.pop(key) for key in _ARRAY_FIELDS[index_name] if key in kwargs}
    if isinstance(fields.get("surface"), str):
        fields["surface"] = _model_top_surface(context, verb=verb)
    values = _array_bc_values(
        fields, context=context, index=index, index_name=index_name, verb=verb, **kwargs
    )
    eager = _array_bc_spec(constructor, name=_spec.name, options=values, context=context)
    if GRID_PIN not in (eager.metadata or {}):
        # No grid on the context means `_array_bc_spec` could only stamp the weaker
        # per-cell MARKER (ledger 282). The built DISV is the same mesh by another
        # route, so the resolved package records WHICH mesh rather than only that it
        # names cells; the marker is dropped because a pin supersedes it.
        disv = _built_discretization(model)
        pin = disv_pin(disv) if disv is not None else None
        if pin is not None:
            metadata = {
                key: value
                for key, value in (eager.metadata or {}).items()
                if key != CELLWISE_OPTIONS
            }
            eager = _replace(
                eager,
                metadata={
                    **metadata,
                    GRID_PIN: pin,
                    GRID_PIN_OPTIONS: sorted(_names_cells(values, int(_ncpl_of(disv)))),
                },
            )
    resolved = getattr(model, "_myflopy_resolved", None)
    if resolved is None:
        resolved = {}
        model._myflopy_resolved = resolved
    resolved[_spec.name] = eager
    return eager.build(model)


def _ncpl_of(discretization) -> int:
    """The cell count a built discretization describes."""

    return int(getattr(getattr(discretization, "ncpl", None), "get_data", lambda: 0)() or 0)


@functools.cache
def _lazy_array_factory(index_name: str):
    """A per-package lazy array builder, bound to its FloPy class.

    Bound arg only, no keywords -- see :func:`_lazy_areal_factory` for the digest that
    a keyword moved. Cached for the reason the other two factories are.
    """

    constructor = (
        flopy.mf6.ModflowGwfrcha if index_name == "irch" else flopy.mf6.ModflowGwfevta
    )
    return functools.partial(_build_lazy_array, constructor, index_name)


def _model_top_surface(context: ModelContext | None, *, verb: str):
    """The ET surface for an array-form EVT that did not name one: the model top.

    **Omitting it used to write the DATUM.** ``surface`` was simply not passed, so
    FloPy's own ``ModflowGwfevta`` default of ``0.0`` reached the file as
    ``surface CONSTANT 0.00000000`` -- while the parameter's documentation promised
    "the model top when omitted". MODFLOW 6 applies the full rate wherever the head is
    at or above the surface, so an ET surface pinned at datum makes the extinction
    depth inert and draws maximum ET everywhere, for the life of the model, under
    Normal termination. Measured on a 4x4 two-layer model with top 100 and heads at 95:
    an ET total of -151.13 against 0.0 for the same declaration naming the top
    (ledger 283). It is ledger 163's ``irch`` trap in the sibling package.

    Resolved through :class:`SurfaceResolver`, the same engine the list form's default
    ``CellSurfaceOffset("model_top")`` uses, so the two shapes of one package agree
    about what "the model top" means and the ``gdf_topbtm`` column layout is read in
    one place (ledger 196/209).

    Refusing when it cannot be resolved is deliberate: silently writing a number the
    caller never gave is what this replaces. ``mf.evt.array`` already refuses
    ``ievt='top_active'`` without a context for the same reason.
    """

    grid = getattr(context, "grid", None)
    ncpl = getattr(grid, "ncpl", None)
    surfaces = _layer_elevations(context) if context is not None else None
    if ncpl is None or surfaces is None:
        raise ValueError(
            f"{verb}: surface= is the ET surface ELEVATION, and omitting it writes "
            "MODFLOW's datum (0.0), not the land surface -- which applies the full ET "
            "rate everywhere and makes depth= inert. Pass surface= explicitly, or pass "
            "context=ModelContext(grid=..., surfaces=...) so the model top can be "
            "resolved."
        )
    return SurfaceResolver(context).surface_column("model_top", layer=0, ncpl=int(ncpl))


def _array_bc_values(
    fields: dict[str, Any],
    *,
    context: ModelContext | None,
    index: Any,
    index_name: str,
    verb: str,
    # Defaulted so the LAZY rebuild can restate only what the caller actually set
    # (ledger 285); the public helpers pass all of them explicitly, as before.
    fixed_cell: bool = False,
    auxiliary: Any = None,
    auxmultname: str | None = None,
    observations: Any = None,
    print_input: bool = False,
    save_flows: bool = True,
    **options: Any,
) -> dict[str, Any]:
    """Assemble the constructor arguments shared by the array-form BC helpers.

    Resolving ``irch``/``ievt`` is the whole of it, and it is not cosmetic:
    without one MODFLOW 6 applies the flux to layer 1 unconditionally and drops
    every column whose layer 1 is inactive, without a word.
    """

    values: dict[str, Any] = dict(fields)
    # `readasarrays` is FloPy's default on both array classes -- the CLASS is the
    # choice, so passing it again would only be a second place to get it wrong.
    values.update(
        fixed_cell=fixed_cell,
        print_input=print_input,
        save_flows=save_flows,
        **options,
    )
    for key, value in (
        ("auxiliary", auxiliary),
        ("auxmultname", auxmultname),
        ("observations", observations),
    ):
        if value is not None:
            values[key] = value

    if index is None:
        logger.debug(
            "%s: %s omitted, so MODFLOW 6 applies the flux to layer 1 in every column",
            verb,
            index_name.upper(),
        )
        _refuse_per_feature_values(values, verb=verb)
        return values

    if isinstance(index, str):
        if index != "top_active":
            raise ValueError(
                f"{verb}: {index_name}= accepts 'top_active', an array, a {{period: array}} "
                f"mapping, or None -- not {index!r}."
            )
        domain = getattr(context, "domain", None)
        if domain is None:
            raise ValueError(
                f"{verb}: {index_name}='top_active' needs the model's active-cell array to "
                f"resolve. Pass context=ModelContext(domain=idomain), give {index_name}= "
                f"explicitly, or pass {index_name}=None to apply the flux to layer 1 "
                "everywhere (which silently skips any column whose layer 1 is inactive)."
            )
        shape = np.shape(domain)
        if len(shape) != 2:
            raise ValueError(
                f"{verb}: expected a (nlay, ncpl) domain to resolve {index_name}, got shape "
                f"{shape}."
            )
        # ZERO-based layer index of the highest ACTIVE cell. FloPy converts to
        # MODFLOW's 1-based IRCH on write by adding one -- measured: passing 1
        # here produced `2` in the file and `Invalid layer number: 3` at the
        # holes. A wholly inactive column takes layer 0; its cell is inactive
        # either way.
        #
        # ACTIVE is `idomain > 0`, the shared rule. This read `!= 0`, so a
        # pass-through cell (-1) counted as active and a column whose top cell was
        # one had its flux written to a cell MF6 has removed -- and MF6 skips that
        # record silently: 42.2% of a pass-through stack's recharge, under Normal
        # termination (ledger 319).
        index = top_active_layer(domain)
        index = np.where(index < 0, 0, index)

    # A numpy INTEGER scalar -- one layer for every column, typed as `np.int64(1)` --
    # died inside FloPy ("object of type 'numpy.int64' has no len()", an
    # MFDataException naming neither argument), eager and lazy alike: numpy's integer
    # scalars subclass no `int` (ledger 265's class; ledger 324).
    if isinstance(index, np.integer):
        index = int(index)
    elif isinstance(index, dict):
        index = {
            period: int(value) if isinstance(value, np.integer) else value
            for period, value in index.items()
        }
    values[index_name] = index
    _refuse_per_feature_values(values, verb=verb)
    return values


class _MVRPackage:
    """Package-first MVR (water mover) helpers.

    Call ``mf.mvr(...)`` to declare **moves** between advanced packages by name and
    id -- e.g. route a stream's outflow into a lake, or rejected UZF infiltration
    into a stream. Each move is an ``mf.Move(provider, receiver)`` of two
    ``mf.MoverConnection(package_name, id)`` endpoints. The provider/receiver packages
    must be declared on the model **with** ``mover=True`` and **before** the mover in
    the package list (myflopy validates this). Use ``mf.mvr.flopy(...)`` for raw
    FloPy ``packages``/``perioddata``.
    """

    def __call__(
        self,
        *,
        nper: int,
        moves: Mapping[int, Sequence[Move]] | Sequence[Move],
        name: str = "mvr",
        print_input: bool = False,
        print_flows: bool = False,
        modelnames: bool = False,
        **options: Any,
    ) -> PackageSpec:
        """Build an MVR package from semantic moves (delegates to ``MVRBuilder``).

        Parameters
        ----------
        nper : int
            Number of stress periods the mover spans; must be at least 1
            (``ValueError`` otherwise). Every period ``0 .. nper-1`` gets a
            ``PERIOD`` block.
        moves : Sequence[Move] or Mapping[int, Sequence[Move]]
            Either a flat sequence of :class:`Move` (applied every period) or a
            ``{period: [Move, ...]}`` mapping (zero-based periods) for time-varying
            routing; a period the mapping leaves out moves nothing, and a key
            outside ``nper`` raises ``ValueError``. Each :class:`Move` connects a
            provider :class:`MoverConnection` to a receiver one
            (``MoverConnection(package_name, id)``, a zero-based feature number),
            optionally with a ``value`` (FACTOR fraction by default). No move at
            all raises ``ValueError``. Every package the moves name becomes a
            ``requires`` of the spec, so it must be declared earlier in the model.
        name : str, optional
            Package name: FloPy's ``pname``, the input file ``<model>.<name>``, the
            budget outputs ``<model>.<name>.bud``/``.csv``, and the spec's name.
            Default ``"mvr"``.
        print_input : bool, optional
            MF6 ``PRINT_INPUT``: echo the mover records to the listing file.
            Default ``False``.
        print_flows : bool, optional
            MF6 ``PRINT_FLOWS``: print the moved rates to the listing file.
            Default ``False``.
        modelnames : bool, optional
            MF6 ``MODELNAMES``: records name a model before each package, for a
            mover between models. The records built from :class:`Move` carry
            package names only, so leave it ``False`` here. Default ``False``.
        **options
            Extra ``flopy.mf6.ModflowGwfmvr`` options, passed at build time; a name
            that class does not accept raises ``TypeError`` at the call. Applied
            last, so they override the ``pname``, ``filename`` and budget file
            records set from ``name``.

        Returns
        -------
        PackageSpec

        Examples
        --------
        >>> # send the main stem's outflow (reach 0) into the lake (lake 0)
        >>> mf.mvr(nper=1, moves=(mf.Move(mf.MoverConnection("sfr", 0),
        ...                               mf.MoverConnection("lak", 0)),))
        """

        return MVRBuilder(
            nper=nper,
            moves=moves,
            name=name,
            print_input=print_input,
            print_flows=print_flows,
            modelnames=modelnames,
            options=options,
        ).build()

    def flopy(
        self,
        *,
        packages: Any,
        perioddata: Any,
        name: str = "mvr",
        maxmvr: int | None = None,
        maxpackages: int | None = None,
        **options: Any,
    ) -> PackageSpec:
        """Return a direct FloPy-style MVR spec from package and period data (escape hatch).

        Use this when you already have MF6 mover records; otherwise prefer the
        semantic ``mf.mvr(moves=...)`` form.

        Parameters
        ----------
        packages : list
            FloPy MVR ``packages`` (the provider/receiver package names).
        perioddata : dict
            FloPy MVR ``perioddata`` keyed by stress period.
        name : str, default "mvr"
            Package name.
        maxmvr : int, optional
            Maximum number of movers (auto-sized from ``perioddata`` when omitted).
        maxpackages : int, optional
            Maximum number of packages involved (auto-sized when omitted).
        **options
            Extra ``flopy.mf6.ModflowGwfmvr`` options.

        Returns
        -------
        PackageSpec
        """

        return mvr_spec(
            packages,
            perioddata,
            name=name,
            maxmvr=maxmvr,
            maxpackages=maxpackages,
            **options,
        )


def _sfr_mover_index(builder: SFRBuilder) -> dict[str, Any]:
    """Serializable reach lookups for wiring movers semantically (see ``sfr_connection``).

    ``cells`` is the grid cell each reach occupies, ``cell_rings`` those cells'
    outlines and ``crs`` the frame they are drawn in -- which is what lets
    ``sfr_connection(at=(x, y))`` put a point in its CELL by the shared rule
    instead of comparing distances to reach centroids (ledger 239). Only the cells
    the stream's reaches occupy are carried: on a 97x97 grid (9,409 cells) with a
    97-reach stream the whole index is 9,992 bytes of JSON where the grid's own
    outlines are 897,304 -- deep-copied by every fork and hashed by every
    ``PackageLibrary.digests()`` -- and every cell that would add is one this
    stream has no reach in.

    Keys are strings throughout. ``PackageSpec.to_dict`` writes JSON, where an
    integer key comes back a string, and ledger 179 promotes those back to ``int``
    only for ``options``/``inputs``, never ``metadata`` -- so an int-keyed lookup
    worked until the spec had been saved and reloaded.
    """

    reaches = builder.reaches
    polys = builder.grid.gdf_vorPolys
    cells = {str(int(row.rno)): int(row.cellid[1]) for row in reaches.itertuples()}
    rings = {
        str(cell): [[float(x), float(y)] for x, y in polys.geometry.loc[cell].exterior.coords]
        for cell in sorted(set(cells.values()))
    }
    return {
        "package": builder.name,
        "outlets": {sid: int(builder.outlet_reach(sid)) for sid in builder.stream_ids},
        "heads": {sid: int(builder.stream_reaches[sid][0]) for sid in builder.stream_ids},
        "by_stream": {sid: [int(r) for r in builder.stream_reaches[sid]] for sid in builder.stream_ids},
        "cells": cells,
        "cell_rings": rings,
        "bounds": [float(value) for value in polys.total_bounds],
        "crs": None if polys.crs is None else str(polys.crs),
    }


def _reach_at_point(index: Mapping[str, Any], stream_id: str, at: Any) -> int:
    """The reach of ``stream_id`` whose grid CELL holds the point ``at``.

    An ``(x, y)`` names a place on the grid, so the grid decides which cell it is
    in -- :func:`myflopy.modflow.mf6.grid.geometry.cell_shares`, the same rule that
    puts a head target or a GeoPackage feature in its cell -- and the reach follows
    from the cell. Comparing the point to each reach piece's CENTROID instead chose
    a reach the point was not in at all: on a 4x4 grid of 25 ft cells with reaches
    in cells 4-6, a point 3 ft inside cell 5's own reach (9.5 ft from that reach's
    centroid) picked cell 6's, whose 10 ft piece was centred 8 ft nearer
    (ledger 239).

    A point on the face between two of this stream's reaches lies in BOTH cells.
    Ledger 215 gives the water at a reach boundary to the DOWNSTREAM reach -- it
    enters there at that reach's upstream end, where the upstream reach would carry
    it along a stretch above the confluence -- so the eligible set is narrowed to
    that one cell before ``cell_shares``' general "lowest-numbered cell wins a tie"
    can hand it upstream. ``SFRBuilder._reach_at`` states the same rule and cannot
    be called here: it needs the reach GEOMETRY and the ``sequence`` column, and a
    serialized index carries neither. Position in ``by_stream`` stands in for
    ``sequence`` and is equal to it, because ``SFRBuilder.reaches`` sorts a
    stream's pieces by projection along its line before numbering them.
    """

    if "cell_rings" not in index:
        raise ValueError(
            "This SFR reach index predates locating a point by cell; rebuild the spec "
            "with mf.sfr(...) to use at=(x, y)."
        )
    reaches = index["by_stream"].get(stream_id)
    if not reaches:
        raise KeyError(f"No stream {stream_id!r}; have {sorted(index['by_stream'])}.")
    x, y = (float(value) for value in at)
    point = Point(x, y)
    # ``by_stream`` runs upstream -> downstream, so a later reach overwrites an
    # earlier one: where two reaches of one stream share a cell the downstream one
    # takes it, which is ledger 215 again.
    order = {int(rno): position for position, rno in enumerate(reaches)}
    cell_reach = {int(index["cells"][str(int(rno))]): int(rno) for rno in reaches}

    def usable(candidates: list[int]) -> list[int]:
        """Which of these cells this stream's water may enter -- at a shared face, the downstream one."""

        reaching = [cell for cell in candidates if cell in cell_reach]
        if len(reaching) < 2:
            return reaching
        return [max(reaching, key=lambda cell: order[cell_reach[cell]])]

    polys = gpd.GeoDataFrame(
        geometry=[Polygon(ring) for ring in index["cell_rings"].values()],
        index=[int(cell) for cell in index["cell_rings"]],
    )
    hits = cell_shares(polys, point, usable=usable)
    if not hits:
        minx, miny, maxx, maxy = index["bounds"]
        frame = index.get("crs")
        where = f" in {frame}" if frame else ""
        if not (minx <= x <= maxx and miny <= y <= maxy):
            raise ValueError(
                f"No reach of stream {stream_id!r} at ({x}, {y}): that point is outside "
                f"the grid, which spans ({minx}, {miny}) to ({maxx}, {maxy}){where} -- "
                f"at= is a coordinate in the grid's own CRS."
            )
        listed = sorted(cell_reach)
        shown = ", ".join(str(cell) for cell in listed[:10])
        if len(listed) > 10:
            shown += f", and {len(listed) - 10} more"
        # Name the nearest reach cell: it is usually the reach the caller meant,
        # and it is what the old nearest-centroid rule handed back in silence.
        nearest = min(listed, key=lambda cell: polys.geometry.loc[cell].distance(point))
        centre = polys.geometry.loc[nearest].centroid
        raise ValueError(
            f"No reach of stream {stream_id!r} at ({x}, {y}): the point is on no cell "
            f"carrying one of its reaches, which are in cells {shown}. The nearest is "
            f"cell {nearest}, holding reach {cell_reach[nearest]}, centred at "
            f"({centre.x:g}, {centre.y:g})."
        )
    return cell_reach[hits[0][0]]


def sfr_connection(sfr_spec: PackageSpec, stream_id: str, at: Any = "downstream") -> MoverConnection:
    """Return a mover endpoint for one SFR reach, chosen *semantically*.

    Instead of hard-coding a raw reach number, point at the reach you mean and let
    myflopy resolve the index from the stream geometry built by ``mf.sfr(...)``.

    Parameters
    ----------
    sfr_spec : PackageSpec
        The spec returned by ``mf.sfr(...)``: it carries the ``sfr_index``
        metadata this reads. A spec without one (``mf.sfr.flopy(...)``) raises
        ``ValueError``.
    stream_id : str
        The stream's id (its ``stream_id`` attribute, e.g. ``"main_stem"``);
        converted with ``str()``. An unknown id raises ``KeyError`` listing the
        ones the spec has.
    at : str or tuple of (float, float), optional
        Which reach of that stream:

        - ``"downstream"`` / ``"outlet"`` / ``"end"`` (default
          ``"downstream"``) -- the stream's final reach;
        - ``"upstream"`` / ``"head"`` / ``"start"`` -- the stream's first reach;
        - an ``(x, y)`` coordinate -- the reach of that stream in the grid CELL
          that holds the point, in the grid's own CRS. A point on the boundary
          between two of its reaches goes to the DOWNSTREAM one (ledger 215); a
          point on no cell this stream has a reach in raises ``ValueError``,
          rather than resolving to whichever reach happens to be nearest.

        The words are case-insensitive; any other string raises ``ValueError``.

    Returns
    -------
    MoverConnection
        Use it as the source/receiver of an ``mf.Move``.

    Examples
    --------
    >>> sfr = mf.sfr(context=ctx, nper=1, streams="streams.gpkg", ...)
    >>> lak = mf.lak(context=ctx, nper=1, lakes="lakes.gpkg", lake_id_field="name", ...)
    >>> mf.mvr(nper=1, moves=(
    ...     mf.Move(mf.sfr_connection(sfr, "main_stem"),            # the stream outlet
    ...             mf.lak_connection(lak, "valley_lake")),))
    >>> mf.sfr_connection(sfr, "main_stem", at=(1500.0, 800.0))     # the reach whose cell holds it
    """

    index = sfr_spec.metadata.get("sfr_index")
    if not index:
        raise ValueError("sfr_connection() needs a spec from mf.sfr(...); this spec has no reach index.")
    sid = str(stream_id)
    if isinstance(at, str):
        key = {"downstream": "outlets", "outlet": "outlets", "end": "outlets",
               "upstream": "heads", "head": "heads", "start": "heads"}.get(at.lower())
        if key is None:
            raise ValueError(
                f"Unknown reach location {at!r}; use 'downstream'/'upstream' or an (x, y) point."
            )
        if sid not in index[key]:
            raise KeyError(f"No stream {sid!r}; have {sorted(index['outlets'])}.")
        rno = index[key][sid]
    else:
        rno = _reach_at_point(index, sid, at)
    return MoverConnection(index["package"], int(rno))


def lak_connection(lak_spec: PackageSpec, lake_id: str) -> MoverConnection:
    """Return a mover endpoint for one lake, by name.

    ``mf.Move`` takes a package name and a zero-based lake NUMBER, which is the
    lake's row in ``mf.lak``'s ``lakes``. This looks the number up by lake id, so a
    reordered lake file cannot send the water to the wrong lake.

    Parameters
    ----------
    lak_spec : PackageSpec
        The spec returned by ``mf.lak(...)``: it carries the ``lak_index`` metadata
        this reads. A spec from ``mf.lak.flopy(...)`` has none.
    lake_id : str
        The lake's id (its ``lake_id_field`` value, e.g. ``"valley_lake"``).
        Converted with ``str()``, so an integer ``101`` finds lake ``"101"``.

    Returns
    -------
    MoverConnection
        The package name and the lake's zero-based MF6 lake number. The package
        needs ``mover=True`` for MF6 to move water through it.

    Raises
    ------
    ValueError
        If ``lak_spec`` carries no lake index (it did not come from ``mf.lak``).
    KeyError
        If no lake has that id; the message lists the ids there are.

    See Also
    --------
    lak : Builds the spec, and numbers its lakes in row order.
    Move : Pairs this endpoint with another in ``mf.mvr``.

    Examples
    --------
    >>> mf.Move(mf.sfr_connection(sfr, "main_stem"), mf.lak_connection(lak, "valley_lake"))
    """

    index = lak_spec.metadata.get("lak_index")
    if not index:
        raise ValueError("lak_connection() needs a spec from mf.lak(...); this spec has no lake index.")
    lakes = index["lakes"]
    if str(lake_id) not in lakes:
        raise KeyError(f"No lake {lake_id!r}; have {sorted(lakes)}.")
    return MoverConnection(index["package"], int(lakes[str(lake_id)]))


chd = _CHDPackage()
ghb = _GHBPackage()
drn = _DRNPackage()
riv = _RIVPackage()
wel = _WELPackage()
rch = _RCHPackage()
evt = _EVTPackage()
uzf = _UZFPackage()
sfr = _SFRPackage()
lak = _LAKPackage()
mvr = _MVRPackage()


__all__ = [
    "adv",
    "chd",
    "cnc",
    "cnd",
    "ctp",
    "dis",
    "disu",
    "disv",
    "drn",
    "dsp",
    "ems",
    "esl",
    "est",
    "evt",
    "ghb",
    "gwe",
    "gwf",
    "gwt",
    "hfb",
    "ic",
    "ims",
    "ist",
    "lak",
    "mip",
    "mst",
    "mvr",
    "npf",
    "oc",
    "prp",
    "prt",
    "rch",
    "riv",
    "sfr",
    "simulation",
    "src",
    "ssm",
    "sto",
    "tdis",
    "uzf",
    "wel",
]

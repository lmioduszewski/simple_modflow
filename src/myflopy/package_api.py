"""Package-first helpers for building readable MODFLOW 6 package specs.

The objects in this module are thin facades over the existing builders and
``PackageSpec`` factories. They keep model scripts close to MODFLOW vocabulary:

``mf.drn(...)`` for direct DRN stress-period data, ``mf.drn.gpkg(...)`` for GIS
boundaries, and ``mf.uzf(...)`` for the higher-level UZF builder.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

import flopy
import geopandas as gpd

from myflopy.advanced import (
    chd_spec,
    drn_spec,
    ghb_spec,
    lak_spec,
    mvr_spec,
    rch_spec,
    sfr_spec,
    uzf_spec,
    wel_spec,
)
from myflopy.builders import build_ims
from myflopy.geopackage import GeoPackageSource, RowValue
from myflopy.modflow.mf6.lakes import LAKBuilder, LakeConnection, LakeOutlet, LakeTable, LakeTableBuilder
from myflopy.modflow.mf6.mvr import MVRBuilder, Move, MoverConnection
from myflopy.modflow.mf6.recharge import RCHBuilder
from myflopy.modflow.mf6.sfr import SFRBuilder, StreamConnection, StreamDiversion
from myflopy.modflow.mf6.uzf import UZFBuilder
from myflopy.specs import (
    ModelContext,
    ModelSpec,
    PackageSpec,
    PostBuildHook,
    SimulationSpec,
)


PathLike = Path | str


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
    """Return a typed GWF model spec wrapping ``flopy.mf6.ModflowGwf``.

    The model-level arguments mirror FloPy's GWF constructor. Package specs
    such as ``mf.disv(...)`` and ``mf.npf(...)`` are supplied through
    ``packages`` and are built after the model is created.
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
    save_flows: bool | None = None,
    dependent_variable_scaling: bool | None = None,
    nc_mesh2d_filerecord: Any = None,
    nc_structured_filerecord: Any = None,
    nc_filerecord: Any = None,
    **kwargs: Any,
) -> ModelSpec:
    """Return a typed GWT (groundwater transport) model spec.

    The solute-transport counterpart of :func:`gwf`: it wraps
    ``flopy.mf6.ModflowGwt`` and takes the transport packages (``mf.disv`` or a
    shared grid, plus the MST/ADV/DSP/SSM/IC/OC transport packages) through
    ``packages``. Couple it to a flow model with a GWF-GWT exchange
    (:func:`build_gwf_gwt_exchange`) in the same :class:`SimulationSpec`.

    Examples
    --------
    >>> transport = mf.gwt("transport", packages=[...])
    >>> mf.SimulationSpec("sim", models=(flow, transport),
    ...                   packages=[mf.tdis(...), mf.ims(models=("flow",)),
    ...                             mf.ims(models=("transport",)),
    ...                             mf.build_gwf_gwt_exchange("flow", "transport")])
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
    save_flows: bool | None = None,
    dependent_variable_scaling: bool | None = None,
    nc_mesh2d_filerecord: Any = None,
    nc_structured_filerecord: Any = None,
    nc_filerecord: Any = None,
    **kwargs: Any,
) -> ModelSpec:
    """Return a typed GWE (groundwater energy / heat transport) model spec.

    The heat-transport counterpart of :func:`gwf`: it wraps
    ``flopy.mf6.ModflowGwe`` and takes the energy packages (EST/ADV/CND/SSM/IC/OC)
    through ``packages``. Couple it to a flow model with a GWF-GWE exchange
    (:func:`build_gwf_gwe_exchange`) in the same :class:`SimulationSpec`. Use it to
    model aquifer thermal transport (ATES, heat pumps, thermal plumes).
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
    ``flopy.mf6.ModflowPrt`` and takes the PRT packages (MIP/PRP/OC) through
    ``packages``. Couple it to a flow model with a GWF-PRT exchange
    (:func:`build_gwf_prt_exchange`) in the same :class:`SimulationSpec` so
    particles advect through that model's flow field. Use it for advective
    pathlines, capture zones, and travel-time analysis. ``context`` carries the
    grid/geometry exactly as for a flow model; the higher-level
    :class:`~myflopy.modflow.mf6.prt.PRTProject` wraps this for common setups.

    Parameters
    ----------
    name
        Model name (also the FloPy model name).
    packages
        PRT package specs to attach (e.g. ``mf.disv`` plus MIP/PRP/OC).
    context
        Geometry/dates :class:`ModelContext` shared with the paired flow model.
    grid
        Optional deferred grid (e.g. ``mf.GridSpec.voronoi(...)``).
    **kwargs
        Additional ``ModflowPrt`` options.
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
    )


def tdis(
    *,
    nper: int = 1,
    perioddata: Any = ((1.0, 1, 1.0),),
    time_units: str | None = None,
    start_date_time: str | None = None,
    ats_perioddata: Any = None,
    filename: str | None = None,
    pname: str | None = None,
    name: str = "tdis",
    **kwargs: Any,
) -> PackageSpec:
    """Return the simulation-wide TDIS package spec.

    MODFLOW 6 has one TDIS package per simulation. GWF, GWT, GWE, and PRT
    models inside the same ``SimulationSpec`` share this timing. To use
    different timing, build a separate simulation/run.

    ``perioddata`` is one ``(perlen, nstp, tsmult)`` tuple per period. This wraps
    ``flopy.mf6.ModflowTdis``.

    Examples
    --------
    >>> mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])                    # one steady period
    >>> mf.tdis(nper=12, perioddata=[(30.0, 3, 1.1)] * 12, time_units="days")
    """

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


# Captured under a private name so simulation()'s ``tdis=`` keyword can shadow
# the public ``tdis`` builder without losing access to the default.
_steady_tdis = tdis


def ims(
    *,
    models: Iterable[str],
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

    ``models`` controls which GWF/GWT/GWE/PRT models use this solver. A single
    IMS can be shared:

    ``mf.ims(models=("flow", "transport"))``

    or separate IMS packages can be declared for different models:

    ``mf.ims(name="flow_solver", models=("flow",))``

    ``mf.ims(name="transport_solver", models=("transport",))``

    ``complexity`` (``"SIMPLE"``/``"MODERATE"``/``"COMPLEX"``) presets the solver
    tolerances; bump it (and ``outer_maximum``/``inner_maximum``) for stiff models
    with wetting/drying or many advanced packages (SFR/LAK/UZF). Pair with a GWF
    ``newtonoptions=...`` for robust water-table convergence.

    All named arguments mirror ``flopy.mf6.ModflowIms``. Extra FloPy keyword
    arguments can still be passed through ``**kwargs``.

    Examples
    --------
    >>> mf.ims(models=("flow",), complexity="MODERATE",
    ...        outer_maximum=100, inner_maximum=100, linear_acceleration="BICGSTAB")
    """

    values = {"models": tuple(models), "pname": name if pname is None else pname, **kwargs}
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


def simulation(
    *models: ModelSpec,
    name: str = "sim",
    tdis: PackageSpec | None = None,
    solver: PackageSpec | Iterable[PackageSpec] | None = None,
    complexity: str | None = "SIMPLE",
    exchanges: Iterable[Any] = (),
    packages: Iterable[PackageSpec] = (),
    **options: Any,
) -> SimulationSpec:
    """Assemble a :class:`SimulationSpec` with sensible timing/solver defaults.

    The common case -- one or more models that should just run -- needs no
    timing or solver boilerplate::

        mf.simulation(flow)              # steady, one IMS solving `flow`
        mf.simulation(flow, transport)   # one IMS each (coupled-ready)

    Defaults: a single steady stress period (``mf.tdis()``) and one ``IMS`` per
    model, each solving only that model -- which is what both single-model and
    coupled GWF/GWT/GWE/PRT runs require. Override any of it::

        mf.simulation(flow, tdis=mf.tdis(nper=12, perioddata=spd))
        mf.simulation(flow, transport, solver=[flow_ims, transport_ims])

    Extra simulation-level packages are appended via ``packages=``; couplings
    via ``exchanges=``.
    """

    if not models:
        raise ValueError("simulation() requires at least one model spec.")

    timing = _steady_tdis() if tdis is None else tdis
    if solver is None:
        solvers: list[PackageSpec] = [
            ims(models=[m.name], name=f"{m.name}_ims", complexity=complexity)
            for m in models
        ]
    elif isinstance(solver, PackageSpec):
        solvers = [solver]
    else:
        solvers = list(solver)

    return SimulationSpec(
        name,
        models=list(models),
        packages=[timing, *solvers, *packages],
        exchanges=tuple(exchanges),
        **options,
    )


def disv(
    *,
    nlay: int,
    ncpl: int,
    nvert: int,
    vertices: Any,
    cell2d: Any,
    top: Any,
    botm: Any,
    idomain: Any | None = None,
    name: str = "disv",
    **options: Any,
) -> PackageSpec:
    """Vertex discretization (DISV) package: the unstructured grid + layers.

    DISV defines the model geometry from a vertex/cell2d mesh plus per-layer ``top``
    and ``botm`` (and optional ``idomain``). With a myflopy grid you rarely type the
    mesh by hand -- pull it from the grid and the layer arrays from a ``LayerStack``::

        gp = vor.get_disv_gridprops()        # ncpl, vertices, cell2d
        layers = stack.build()               # top, botm, idomain, nlay
        mf.disv(nlay=layers.nlay, ncpl=gp["ncpl"], nvert=len(gp["vertices"]),
                vertices=gp["vertices"], cell2d=gp["cell2d"],
                top=layers.top, botm=layers.botm, idomain=layers.idomain)

    Parameters
    ----------
    nlay, ncpl, nvert
        Number of layers, cells-per-layer, and vertices.
    vertices, cell2d
        FloPy vertex and cell2d definition lists (``vor.get_disv_gridprops()``).
    top, botm, idomain
        Model-top (ncpl,), layer bottoms (nlay, ncpl), and active-domain array
        (idomain 0 = inactive / pinched out).
    """

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
    return PackageSpec(name, flopy.mf6.ModflowGwfdisv, values)


def ic(*, strt: Any, name: str = "ic", **options: Any) -> PackageSpec:
    """Initial-conditions (IC) package: the starting head field.

    ``strt`` is the starting head -- a scalar, a per-cell array, or a ``(nlay, ncpl)``
    array (one value per cell per layer). For Newton/under-relaxation runs a sensible
    starting head (e.g. near the water table) helps convergence.

    Examples
    --------
    >>> mf.ic(strt=9.0)                 # uniform
    >>> mf.ic(strt=starting_heads)      # per-cell / per-layer array
    """

    return PackageSpec(name, flopy.mf6.ModflowGwfic, {"strt": strt, **options})


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

    ``k`` is horizontal K and ``k33`` vertical K -- each a scalar, per-cell array, or
    ``(nlay, ncpl)`` array. ``icelltype`` controls confined (0) vs convertible
    (non-zero, water-table) behavior. ``save_specific_discharge=True`` is needed for
    particle tracking / Darcy-velocity plots.

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
    return PackageSpec(name, flopy.mf6.ModflowGwfnpf, values)


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

    Provide ``ss`` (specific storage) and ``sy`` (specific yield) for transient
    periods, and mark each period steady or transient. A purely steady model still
    needs STO with ``steady_state={0: True}``.

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
    return PackageSpec(name, flopy.mf6.ModflowGwfsto, values)


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

    To get heads and a cell budget on disk you must name the output files
    (``head_filerecord`` / ``budget_filerecord``) **and** request them in
    ``saverecord`` -- MF6 errors if you ask to save heads without a head file.

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
    return PackageSpec(name, flopy.mf6.ModflowGwfoc, values)


class _CHDPackage:
    """Package-first CHD (constant/specified-head boundary) helpers.

    A CHD pins listed cells to a specified head -- used for fixed boundaries (a
    lake edge held at stage, a regional head condition). Three entry points:

    - ``mf.chd(stress_period_data=...)`` -- direct MF6 records;
    - ``mf.chd.gpkg(path, context=, nper=)`` -- build records from GeoPackage
      features mapped onto cells (``edges_only=True`` keeps just boundary cells);
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
        """Return a CHD package spec from direct stress-period data."""

        return chd_spec(stress_period_data, name=name, boundnames=boundnames, **options)

    def gpkg(
        self,
        path: PathLike,
        *,
        context: ModelContext,
        nper: int,
        head: RowValue = "head",
        layer: str | None = None,
        name_field: str | None = "name",
        layer_field: str | None = "layer",
        period_field: str | None = None,
        layer_base: int = 1,
        period_base: int = 0,
        name: str = "chd",
        edges_only: bool = False,
        boundnames: bool = True,
        **options: Any,
    ) -> PackageSpec:
        """Return a CHD package spec from GeoPackage features."""

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
    - same as ``()`` here (the FloPy-native form).

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

        ``stress_period_data`` is the FloPy mapping ``{period: [[cellid, bhead,
        cond], ...]}`` (cellid is ``(layer, cell)`` on DISV).
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
        context: ModelContext,
        nper: int,
        head: RowValue = "head",
        conductance: RowValue = "conductance",
        layer: str | None = None,
        name_field: str | None = "name",
        layer_field: str | None = "layer",
        period_field: str | None = None,
        layer_base: int = 1,
        period_base: int = 0,
        name: str = "ghb",
        edges_only: bool = False,
        boundnames: bool = True,
        **options: Any,
    ) -> PackageSpec:
        """Return a GHB package spec from GeoPackage features."""

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

        ``stress_period_data`` is ``{period: [[cellid, elev, cond], ...]}``.
        """

        return drn_spec(stress_period_data, name=name, boundnames=boundnames, **options)

    def gpkg(
        self,
        path: PathLike,
        *,
        context: ModelContext,
        nper: int,
        elevation: RowValue = "elevation",
        conductance: RowValue = "conductance",
        layer: str | None = None,
        name_field: str | None = "name",
        layer_field: str | None = "layer",
        period_field: str | None = None,
        layer_base: int = 1,
        period_base: int = 0,
        name: str = "drn",
        edges_only: bool = False,
        boundnames: bool = True,
        **options: Any,
    ) -> PackageSpec:
        """Return a DRN package spec from GeoPackage features."""

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
        ).drn(
            elevation=elevation,
            conductance=conductance,
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
        """Return a WEL package spec from direct stress-period data."""

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
        context: ModelContext,
        nper: int,
        rate: RowValue = "rate",
        layer: str | None = None,
        name_field: str | None = "name",
        layer_field: str | None = "layer",
        period_field: str | None = None,
        layer_base: int = 1,
        period_base: int = 0,
        name: str = "wel",
        boundnames: bool = True,
        **options: Any,
    ) -> PackageSpec:
        """Return a WEL package spec from GeoPackage features."""

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
        """Return an RCH package spec.

        Two forms, mirroring the other stress packages:

        * Direct ``stress_period_data=`` (like ``mf.drn`` / ``mf.wel``) returns a
          list-based RCH spec.
        * The high-level builder form (``context=``, ``nper=``, ``recharge=``)
          computes the cells from the model domain; ``recharge`` is a scalar,
          per-cell sequence, or ``{period: ...}`` mapping.
        """

        if stress_period_data is not None:
            return rch_spec(
                stress_period_data, name=name, boundnames=boundnames, **options
            )
        if context is None or nper is None or recharge is None:
            raise TypeError(
                "mf.rch requires either stress_period_data=, or the builder "
                "arguments context=, nper=, and recharge=."
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

    def flopy(
        self,
        *,
        stress_period_data: Any,
        name: str = "rch",
        boundnames: bool = False,
        **options: Any,
    ) -> PackageSpec:
        """Return a direct FloPy-style RCH spec from stress-period data."""

        return rch_spec(stress_period_data, name=name, boundnames=boundnames, **options)

    def gpkg(
        self,
        path: PathLike,
        *,
        context: ModelContext,
        nper: int,
        recharge: RowValue = "recharge",
        layer: str | None = None,
        name_field: str | None = "name",
        layer_field: str | None = "layer",
        period_field: str | None = None,
        layer_base: int = 1,
        period_base: int = 0,
        name: str = "rch",
        boundnames: bool = True,
        **options: Any,
    ) -> PackageSpec:
        """Return a list-based RCH package spec from GeoPackage features."""

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
        ).rch(
            recharge=recharge,
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
        context: ModelContext,
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
        ``cells`` list to override. Soil-hydraulic inputs (``vks``/``thtr``/``thts``/
        ``thti``/``eps``) and the stress inputs (``finf``/``pet``/``extdp``...) accept
        a scalar (applied everywhere), a per-cell sequence, or a ``{period: ...}``
        mapping for transient values.

        Parameters
        ----------
        context
            :class:`ModelContext` carrying the grid + domain (idomain) used to pick
            and order UZF cells. Required for ``cells="all_active"/"surface_only"``.
        nper
            Number of stress periods.
        vks, thtr, thts, thti, eps
            Saturated K of the unsaturated zone, residual / saturated / initial water
            content, and the Brooks-Corey exponent.
        cells
            ``"all_active"`` (default), ``"surface_only"``, or an explicit list of
            ``(layer, cell)`` ids.
        finf, pet, extdp, extwc
            Infiltration rate, PET, ET extinction depth, and extinction water content
            (scalar / per-cell / ``{period: ...}``).
        mover
            Enable MVR so rejected infiltration / groundwater discharge can be moved
            to another package.

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
        """Return a direct FloPy-style UZF spec from prepared package data."""

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
        """Build an SFR package from stream geometry (delegates to ``SFRBuilder``).

        Reaches are cut where each stream line crosses grid cells; ``connection_mode
        ="automatic"`` infers the topology from geometry. When tributaries meet
        ambiguously, disambiguate with explicit ``connections`` (and ``diversions``).
        Reach properties (``width``/``gradient``/``roughness``/``streambed_k``/
        ``streambed_thickness``) accept a scalar, one value per reach, or one per
        stream; ``reach_top`` defaults from the grid surface when available.

        Parameters
        ----------
        context
            :class:`ModelContext` carrying the grid the reaches are placed on.
        nper
            Number of stress periods.
        streams
            Stream centerlines: a geopackage/shapefile path (or list), or a
            ``GeoDataFrame`` of ``LineString`` rows. ``stream_id`` names the id column.
        connection_mode, connections, diversions
            ``"automatic"`` topology from geometry, plus explicit
            :class:`StreamConnection`/:class:`StreamDiversion` overrides for
            confluences/splits.
        inflow, rainfall, evaporation, runoff, status
            Per-period boundary inputs, keyed by stream id (e.g.
            ``inflow={0: {"main": 1.0e4}}``).
        mover
            Enable MVR so this stream can give/receive water from lakes or UZF.
        length_conversion, time_conversion
            Manning's-equation unit conversions (e.g. ``3.28081`` ft/m, ``86400`` s/day).

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
        # (by stream outlet/head or by nearest reach to a coordinate) -- see
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
        """Return a direct FloPy-style SFR spec from prepared package data."""

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
    """Package-first LAK (lake) helpers.

    Call ``mf.lak(...)`` to build a LAK package from **lake polygon geometry** (a
    geopackage/shapefile path or a GeoDataFrame): it finds the lake cells, builds
    the lake-aquifer connections, and writes packagedata/connectiondata/perioddata.
    Use ``mf.lak.flopy(...)`` for prepared FloPy records. A lake is a head-dependent
    storage that exchanges water with the aquifer through its bed and can connect to
    streams via MVR and to outlets.
    """

    def __call__(
        self,
        *,
        context: ModelContext,
        nper: int,
        lakes: PathLike | Sequence[PathLike] | gpd.GeoDataFrame,
        lake_id_field: str | None = None,
        starting_stage: Any = None,
        lake_bottom: Any = None,
        lake_top: Any = None,
        only_vertical: Any = False,
        only_layer: Any = None,
        bed_leakance: Any = 1.0,
        connection_modes: str | Mapping[str, str | Sequence[LakeConnection]] = "bathy",
        tables: Mapping[str, LakeTable | LakeTableBuilder] | None = None,
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
        **options: Any,
    ) -> PackageSpec:
        """Build a LAK package from lake geometry (delegates to ``LAKBuilder``).

        Lake cells are taken from the polygon footprint(s); ``connection_modes`` picks
        how the bed connections are generated per lake:

        * ``"bathy"`` (default; alias ``"automatic"``) -- natural lakes. Needs a
          *per-cell* ``lake_bottom`` (a raster ``Path`` or a cell-indexed ``Series``).
          One vertical connection per cell, plus a horizontal connection toward each
          neighbor whose lake bottom is *higher* (the exposed step), clipped to every
          layer the exposed face spans.
        * ``"rectangular"`` -- box facilities (infiltration trenches / vaults). Needs a
          flat ``lake_bottom`` and an explicit ``lake_top`` per lake. Vertical
          connections everywhere (flat-bottom infiltration); horizontal connections
          only on the *perimeter* (cells bordering non-lake cells), spanning
          ``lake_bottom -> lake_top`` clipped to each layer. Set ``only_vertical`` for a
          vault (concrete walls -> no horizontals).

        Connection ``telev``/``belev`` are pure geometry and never depend on the
        (transient) starting stage -- MF6 wets/dries each face per timestep.

        Per-lake inputs (``starting_stage``/``lake_bottom``/``lake_top``/
        ``bed_leakance``/``status``/forcings) are keyed by lake id (the
        ``lake_id_field`` attribute, e.g. ``{"valley_lake": 101.0}``).

        Parameters
        ----------
        context
            :class:`ModelContext` with the grid + domain + surfaces (lake-cell layer
            assignment reads cell-bottom elevations).
        nper
            Number of stress periods.
        lakes
            Lake polygons: a geopackage/shapefile path (or list), or a GeoDataFrame.
            ``lake_id_field`` is the attribute naming each lake.
        starting_stage, lake_bottom, bed_leakance
            Initial stage, lake-bottom elevation, and bed leakance per lake.
        lake_top
            Facility top per lake -- **required for ``"rectangular"``** lakes; the
            horizontal-connection ``telev``. Ignored by ``"bathy"``.
        only_vertical
            ``True`` (or a per-lake mapping) to build vault lakes with vertical
            connections only (no horizontal sidewall exchange).
        only_layer
            Restrict a lake's connections to a single model layer (special cases only);
            an int or per-lake mapping. Default ``None`` = all layers a face spans.
        outlets, tables
            Optional lake outlets and stage-volume-area tables.
        mover
            Enable MVR so the lake can give/receive water from streams or UZF.

        Returns
        -------
        PackageSpec

        Examples
        --------
        >>> ctx = mf.ModelContext(grid=vor, domain=idomain, surfaces=vor.gdf_topbtm)
        >>> mf.lak(context=ctx, nper=1, lakes="lakes.gpkg", lake_id_field="name",
        ...        starting_stage={"valley_lake": 101.0}, lake_bottom={"valley_lake": 96.0},
        ...        bed_leakance=0.1, status={"valley_lake": ["ACTIVE"]}, mover=True)
        """

        builder = LAKBuilder(
            context=context,
            nper=nper,
            lakes=lakes,
            lake_id_field=lake_id_field,
            starting_stage=starting_stage,
            lake_bottom=lake_bottom,
            lake_top=lake_top,
            only_vertical=only_vertical,
            only_layer=only_layer,
            bed_leakance=bed_leakance,
            connection_modes=connection_modes,
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
        """Return a direct FloPy-style LAK spec from prepared package data."""

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
        nper
            Number of stress periods.
        moves
            Either a flat sequence of :class:`Move` (applied every period) or a
            ``{period: [Move, ...]}`` mapping for time-varying routing. Each
            :class:`Move` connects a provider :class:`MoverConnection` to a receiver
            one (``MoverConnection(package_name, id)``), optionally with a ``value``
            (FACTOR fraction by default).

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
        """Return a direct FloPy-style MVR spec from package and period data."""

        return mvr_spec(
            packages,
            perioddata,
            name=name,
            maxmvr=maxmvr,
            maxpackages=maxpackages,
            **options,
        )


def _sfr_mover_index(builder: SFRBuilder) -> dict[str, Any]:
    """Serializable reach lookups for wiring movers semantically (see ``sfr_connection``)."""

    reaches = builder.reaches
    return {
        "package": builder.name,
        "outlets": {sid: int(builder.outlet_reach(sid)) for sid in builder.stream_ids},
        "heads": {sid: int(builder.stream_reaches[sid][0]) for sid in builder.stream_ids},
        "by_stream": {sid: [int(r) for r in builder.stream_reaches[sid]] for sid in builder.stream_ids},
        "centroids": {
            int(row.rno): [float(row.geometry.centroid.x), float(row.geometry.centroid.y)]
            for row in reaches.itertuples()
        },
    }


def sfr_connection(sfr_spec: PackageSpec, stream_id: str, at: Any = "downstream") -> MoverConnection:
    """Return a mover endpoint for one SFR reach, chosen *semantically*.

    Instead of hard-coding a raw reach number, point at the reach you mean and let
    myflopy resolve the index from the stream geometry built by ``mf.sfr(...)``.

    Parameters
    ----------
    sfr_spec
        The spec returned by ``mf.sfr(...)`` (carries the reach index).
    stream_id
        The stream's id (its ``stream_id`` attribute, e.g. ``"main_stem"``).
    at
        Which reach of that stream:

        - ``"downstream"`` / ``"outlet"`` (default) -- the stream's final reach;
        - ``"upstream"`` / ``"head"`` -- the stream's first reach;
        - an ``(x, y)`` coordinate -- the reach of that stream nearest the point.

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
    >>> mf.sfr_connection(sfr, "main_stem", at=(1500.0, 800.0))     # nearest reach to a point
    """

    index = sfr_spec.metadata.get("sfr_index")
    if not index:
        raise ValueError("sfr_connection() needs a spec from mf.sfr(...); this spec has no reach index.")
    sid = str(stream_id)
    if isinstance(at, str):
        key = {"downstream": "outlets", "outlet": "outlets", "end": "outlets",
               "upstream": "heads", "head": "heads", "start": "heads"}.get(at.lower())
        if key is None:
            raise ValueError(f"Unknown reach location {at!r}; use 'downstream'/'upstream' or an (x, y) point.")
        if sid not in index[key]:
            raise KeyError(f"No stream {sid!r}; have {sorted(index['outlets'])}.")
        rno = index[key][sid]
    else:
        x, y = at
        centroids = index["centroids"]
        candidates = index["by_stream"].get(sid)
        if not candidates:
            raise KeyError(f"No stream {sid!r}; have {sorted(index['by_stream'])}.")
        rno = min(candidates, key=lambda r: (centroids[r][0] - x) ** 2 + (centroids[r][1] - y) ** 2)
    return MoverConnection(index["package"], int(rno))


def lak_connection(lak_spec: PackageSpec, lake_id: str) -> MoverConnection:
    """Return a mover endpoint for one lake, by name.

    Parameters
    ----------
    lak_spec
        The spec returned by ``mf.lak(...)``.
    lake_id
        The lake's id (its ``lake_id_field`` value, e.g. ``"valley_lake"``).

    Returns
    -------
    MoverConnection

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
wel = _WELPackage()
rch = _RCHPackage()
uzf = _UZFPackage()
sfr = _SFRPackage()
lak = _LAKPackage()
mvr = _MVRPackage()


__all__ = [
    "chd",
    "disv",
    "drn",
    "ghb",
    "gwe",
    "gwf",
    "gwt",
    "ic",
    "ims",
    "lak",
    "mvr",
    "npf",
    "oc",
    "prt",
    "rch",
    "sfr",
    "simulation",
    "sto",
    "tdis",
    "uzf",
    "wel",
]

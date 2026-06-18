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
from myflopy.modflow.mf6.mvr import MVRBuilder, Move
from myflopy.modflow.mf6.recharge import RCHBuilder
from myflopy.modflow.mf6.sfr import SFRBuilder, StreamConnection, StreamDiversion
from myflopy.modflow.mf6.uzf import UZFBuilder
from myflopy.specs import ModelContext, ModelSpec, PackageSpec, PostBuildHook


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
    """Return a typed GWT model spec wrapping ``flopy.mf6.ModflowGwt``."""

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
    """Return a typed GWE model spec wrapping ``flopy.mf6.ModflowGwe``."""

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
    """Return a typed PRT model spec wrapping ``flopy.mf6.ModflowPrt``."""

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

    This wraps ``flopy.mf6.ModflowTdis``.
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

    All named arguments mirror ``flopy.mf6.ModflowIms``. Extra FloPy keyword
    arguments can still be passed through ``**kwargs``.
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
    """Return a DISV package spec wrapping ``flopy.mf6.ModflowGwfdisv``."""

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
    """Return an IC package spec wrapping ``flopy.mf6.ModflowGwfic``."""

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
    """Return an NPF package spec wrapping ``flopy.mf6.ModflowGwfnpf``."""

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
    """Return an STO package spec wrapping ``flopy.mf6.ModflowGwfsto``."""

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
    """Return an OC package spec wrapping ``flopy.mf6.ModflowGwfoc``."""

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
    """Package-first CHD helpers."""

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
    """Package-first GHB helpers."""

    def __call__(
        self,
        *,
        stress_period_data: Any,
        name: str = "ghb",
        auxiliary: Any = None,
        boundnames: bool = False,
        **options: Any,
    ) -> PackageSpec:
        """Return a GHB package spec from direct stress-period data."""

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
    """Package-first DRN helpers."""

    def __call__(
        self,
        *,
        stress_period_data: Any,
        name: str = "drn",
        boundnames: bool = False,
        **options: Any,
    ) -> PackageSpec:
        """Return a DRN package spec from direct stress-period data."""

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
    """Package-first WEL helpers."""

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
    """Package-first RCH helpers."""

    def __call__(
        self,
        *,
        context: ModelContext,
        nper: int,
        recharge: Any,
        cells: str | Sequence[int | tuple[int, int]] = "top_active",
        layer: int = 0,
        name_by_cell: Mapping[int | tuple[int, int], str] | None = None,
        boundnames: bool = False,
        name: str = "rch",
        **options: Any,
    ) -> PackageSpec:
        """Return a high-level RCH package spec from the model domain."""

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
    """Package-first UZF helpers."""

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
        """Return a high-level UZF package spec from domain-aware inputs."""

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
    """Package-first SFR helpers."""

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
        length_conversion: float = 1.0,
        time_conversion: float = 1.0,
        maximum_picard_iterations: int = 1,
        maximum_iterations: int = 1000,
        maximum_depth_change: float = 0.01,
        **options: Any,
    ) -> PackageSpec:
        """Return a high-level SFR package spec from stream geometry."""

        return SFRBuilder(
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
        ).build()

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
    """Package-first LAK helpers."""

    def __call__(
        self,
        *,
        context: ModelContext,
        nper: int,
        lakes: PathLike | Sequence[PathLike] | gpd.GeoDataFrame,
        lake_id_field: str | None = None,
        starting_stage: Any = None,
        lake_bottom: Any = None,
        bed_leakance: Any = 1.0,
        connection_modes: str | Mapping[str, str | Sequence[LakeConnection]] = "automatic",
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
        length_conversion: float = 1.0,
        time_conversion: float = 1.0,
        maximum_iterations: int = 100,
        maximum_stage_change: float = 1.0e-5,
        **options: Any,
    ) -> PackageSpec:
        """Return a high-level LAK package spec from lake geometry."""

        return LAKBuilder(
            context=context,
            nper=nper,
            lakes=lakes,
            lake_id_field=lake_id_field,
            starting_stage=starting_stage,
            lake_bottom=lake_bottom,
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
        ).build()

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
    """Package-first MVR helpers."""

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
        """Return a high-level MVR package spec from semantic moves."""

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
    "sto",
    "tdis",
    "uzf",
    "wel",
]

# `myflopy` Preferred API

`myflopy` uses immutable specifications to describe models and `Project` /
`Run` objects to materialize them.

## Core Flow

```python
import myflopy as mf

project = mf.Project("modeling/my-project")
project.add_simulation(simulation_spec)
run = project.run("baseline", "simulation-name")
```

- `PackageSpec` describes one FloPy package.
- `ModelSpec` describes one GWF, GWT, GWE, or PRT model.
- `SimulationSpec` groups models, solvers, time discretization, and exchanges.
- `ModelContext` carries grids, surfaces, dates, domains, and metadata.
- `PostBuildHook` handles observations, targets, and validation after packages exist.
- `Project` owns related run workspaces.
- `Run` owns one built or completed simulation.

## GeoPackage Inputs

`GeoPackageSource` is the public GIS-to-model API. It does not require a live
model and returns reusable package specifications or arrays.

```python
context = mf.ModelContext(
    grid=vor,
    surfaces={"top": top, "bottom": bottom},
    dates=dates,
    domain=idomain,
)

chd = mf.GeoPackageSource(
    "inputs/constant_heads.gpkg",
    context,
    nper=12,
).chd(head=["head_jan", "head_feb", "head_mar"])

ghb = mf.GeoPackageSource(
    "inputs/general_heads.gpkg",
    context,
    nper=12,
).ghb(head="head", conductance="conductance")

drn = mf.GeoPackageSource("inputs/drains.gpkg", context, nper=12).drn()
wel = mf.GeoPackageSource("inputs/wells.gpkg", context, nper=12).wel()
rch = mf.GeoPackageSource("inputs/recharge.gpkg", context, nper=12).rch()

k = mf.GeoPackageSource(
    "inputs/hydraulic_conductivity.gpkg",
    context,
    nper=12,
).k_array(value="k", nlay=4, defaults=[10.0, 5.0, 1.0, 0.5])
```

External GIS layer values are one-based by default. Stress periods are
zero-based by default. Use `layer_base` and `period_base` when inputs use
different conventions.

A single value field applies to every stress period. A sequence of fields maps
values by stress period, repeating the last field when needed:

```python
source.chd(head=["head_0", "head_1", "head_2"])
```

For long-format GeoPackages containing a period column:

```python
source = mf.GeoPackageSource(
    "inputs/wells.gpkg",
    context,
    nper=12,
    period_field="period",
)
wel = source.wel(rate="rate")
```

## Compose a GWF Model

```python
flow = mf.ModelSpec(
    "flow",
    "gwf",
    context=context,
    packages=(
        disv,
        ic,
        npf,
        sto,
        chd,
        ghb,
        drn,
        wel,
        rch,
        uzf,
        lak,
        sfr,
        mvr,
        oc,
    ),
    hooks=(
        mf.PostBuildHook("observations", attach_observations),
        mf.PostBuildHook("validation", validate_model),
    ),
)
```

MVR should be declared after the packages it references. Prefer `mf.mvr(...)`
with semantic endpoints from the package helpers where possible.

## Advanced Packages

Prefer package-first helpers for advanced packages. Their `.flopy(...)` methods
are the explicit escape hatch when you already have prepared FloPy-style data;
the older `*_spec` factories remain available for compatibility.

```python
lak = mf.lak(
    context=context,
    nper=nper,
    lakes=lake_polygons,
    lake_id_field="lake_id",
    starting_stage="stage",
    mover=True,
)
sfr = mf.sfr(context=context, nper=nper, streams=streams, stream_id="stream_id")
mvr = mf.mvr(
    nper=nper,
    moves=[
        mf.Move(
            source=sfr_builder.connection("main_stem"),
            receiver=lak_builder.connection("north_lake"),
            value=0.25,
        )
    ],
)

lak_from_prepared_data = mf.lak.flopy(
    packagedata=lake_packagedata,
    connectiondata=lake_connections,
    perioddata=lake_perioddata,
    mover=True,
)
```

## Simulation and Scenarios

```python
baseline = mf.SimulationSpec(
    "canonical",
    models=(flow,),
    packages=(tdis, ims),
)

npf = baseline.model("flow").package("npf")
high_k_flow = baseline.model("flow").with_package(
    npf.with_options(k=k * 2.0)
)
high_k = baseline.with_model(high_k_flow)

baseline_run = project.run("baseline", baseline)
high_k_run = project.run("high-k", high_k)
```

GWT, GWE, and PRT models belong beside GWF models within `SimulationSpec` and
are connected with explicit `ExchangeSpec` values.

# Core Specifications

`myflopy` targets FloPy 3.10 and uses its public MF6 constructors directly.
The core specification API organizes those constructors without hiding them.

## Vocabulary

- `PackageSpec` describes one model-level or simulation-level package.
- `ModelSpec` describes one GWF, GWT, GWE, or PRT model and its packages.
- `ModelContext` carries grid, surfaces, dates, domain information, and metadata.
- `PostBuildHook` runs after all packages on a model have been built.
- `ExchangeSpec` connects named models within a simulation.
- `SimulationSpec` owns the coupled models, simulation packages, and exchanges.
- `Run` owns a workspace, a built or loaded simulation, and model access.
- `ModelView` is the preferred myflopy object for a GWF model, with helpers
  such as `cor()`, `hds`, `bud()`, `packages`, and `outputs`.
- `BuiltSimulation` is the lower-level container for raw FloPy objects.

Specs are immutable. Methods such as `with_package()` return a changed copy,
which makes scenario variants explicit without mutating the baseline.

## Package Builder Contract

Package builders hold their complete configuration and expose a zero-argument
`build()` method:

```python
builder = UZFBuilder(
    context=context,
    nper=12,
    cells="top_active",
    vks=soil.values("vks"),
    thtr=0.05,
    thts=0.30,
    thti=0.15,
    finf=recharge.period_values("recharge"),
    pet=climate.period_values("pet"),
)

uzf = builder.build()
```

Calling `build()` never supplies or overrides configuration. Scenario variants
are created from a changed builder configuration, then built independently.
This contract applies to package-data builders such as `UZFBuilder`,
`SFRBuilder`, and `LakeBuilder`.

Specs and runs have different responsibilities. `PackageSpec.build(target)`
needs the target FloPy object. `SimulationSpec.build()` creates a `Run`, using
the workspace stored on the spec or an explicit workspace argument. Use
`SimulationSpec.build_flopy()` only when you want the raw in-memory FloPy
objects without the run/model workflow wrapper.

## Coupled GWF and GWT Example

```python
from myflopy import (
    ExchangeSpec,
    ModelSpec,
    SimulationSpec,
    build_gwf_gwt_exchange,
    ims,
    load_run,
    npf,
    tdis,
)

flow_npf = npf(k=10.0)
flow = ModelSpec("flow", "gwf", packages=(flow_npf,))
transport = ModelSpec("transport", "gwt")

simulation = SimulationSpec(
    "coupled",
    models=(flow, transport),
    packages=(
        tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
        ims(name="flow_solver", models=("flow",)),
        ims(name="transport_solver", models=("transport",)),
    ),
    exchanges=(
        ExchangeSpec(
            "flow_transport",
            build_gwf_gwt_exchange,
            models=("flow", "transport"),
        ),
    ),
)

run = simulation.build("runs/coupled")
flow_model = run.model("flow")
```

For a single-model simulation, the model name is optional:

```python
model = run.model()
model.gwf       # raw flopy.mf6.ModflowGwf
model.cor()     # myflopy map/dashboard helper
model.hds       # heads helper
```

The same shape is used when reopening an existing workspace:

```python
run = load_run("runs/coupled")
model = run.model("flow")
```

`tdis(...)` is simulation-wide. All models inside one `SimulationSpec` share
the same stress-period timing. `ims(...)` is different: it creates a solver
package and registers it to the named model or models. You can share one solver:

```python
shared_solver = ims(models=("flow", "transport"), print_option="SUMMARY")
```

or split solvers by model:

```python
flow_solver = ims(name="flow_solver", models=("flow",), complexity="SIMPLE")
transport_solver = ims(name="transport_solver", models=("transport",), complexity="COMPLEX")
```

Package alternatives are ordinary values:

```python
high_k_flow = flow.with_package(flow_npf.with_options(k=100.0))
high_k_simulation = simulation.with_model(high_k_flow)
```

## Model Context and Hooks

Keep domain information beside the build instructions instead of hiding it in
FloPy package arguments:

```python
context = ModelContext(
    grid=vor,
    surfaces=surfaces,
    dates=dates,
    domain=idomain,
    metadata={"scenario": "baseline"},
)

flow = ModelSpec(
    "flow",
    "gwf",
    context=context,
    hooks=(PostBuildHook("observations", attach_observations),),
)
```

A hook receives the completed FloPy model, its built package dictionary, and
the model context. Hook return values are retained in
`BuiltModel.hook_results`.

Package-first helpers return ordinary `PackageSpec` values:

```python
flow = ModelSpec(
    "flow",
    "gwf",
    context=ModelContext(grid=vor, surfaces=surfaces, dates=dates),
    packages=(
        wel(stress_period_data=well_data),
        ghb(stress_period_data=ghb_data),
        uzf.flopy(packagedata=uzf_packagedata, perioddata=uzf_perioddata),
        lak.flopy(
            packagedata=lake_packagedata,
            connectiondata=lake_connections,
            perioddata=lake_perioddata,
            mover=True,
        ),
        sfr.flopy(
            packagedata=sfr_packagedata,
            connectiondata=sfr_connections,
            perioddata=sfr_perioddata,
            mover=True,
        ),
        mvr(
            nper=nper,
            moves=[
                Move(
                    source=MoverConnection("sfr", outlet),
                    receiver=MoverConnection("lak", 0),
                )
            ],
        ),
    ),
)
```

`MVRBuilder` declares source and destination package dependencies from its
`Move` records. `MoverConnection` is the low-level endpoint, while builders such
as `SFRBuilder` and `LAKBuilder` expose helpers like
`sfr.connection("main_stem")` and `lak.connection("north_lake")` so examples can
stay keyed to stable feature IDs.

Advanced package builders keep their complete configuration and return normal
package specs without controlling package attachment:

```python
uzf = UZFBuilder(
    context=context,
    nper=nper,
    cells="all_active",
    vks=vks,
    thtr=thtr,
    thts=thts,
    thti=thti,
    finf=finf,
).build()
```

The default `cells="all_active"` builds a vertical UZF chain from the
uppermost active cell through each directly connected active underlying layer.
The uppermost cell has `landflag=1`; underlying cells have `landflag=0` and
are linked with `ivertcon`. Period inputs such as `finf` and `pet` are applied
only to land-surface cells. An inactive layer stops the chain.

The resulting values are normal package specs and can be replaced, disabled,
validated, and reused like any other package.

### Stream Networks

`SFRBuilder` resolves a grid-independent stream network before mapping streams
to generated SFR reaches:

```python
sfr = SFRBuilder(
    context=context,
    nper=nper,
    streams="streams.gpkg",
    stream_id="stream_id",
    connection_mode="automatic",
    connection_tolerance=25.0,
    width="width",
    reach_top="bed_elevation",
    roughness="mannings",
    streambed_k="streambed_k",
    streambed_thickness="bed_thickness",
    mover=True,
)

flow = flow.with_package(sfr.build())
```

Automatic mode connects a stream's downstream endpoint to another stream
geometry within the configured tolerance. Explicit connections override
automatic decisions without referring to generated cells or reach numbers:

```python
connection = StreamConnection(
    source="tributary_a",
    receiver="main_stem",
    source_location="downstream",
    receiver_location="nearest",
)
```

`source_location` describes where water leaves the source stream.
`receiver_location` describes where it enters the receiving stream. Both may
also be geographic `Point` values. Set `receiver=None` to declare an outlet.
Node-field topology is available with `connection_mode="nodes"`,
`from_node="from_node"`, and `to_node="to_node"`.

Diversions are always explicit:

```python
diversion = StreamDiversion(
    source="main_stem",
    receiver="canal",
    source_location="downstream",
    receiver_location="nearest",
    priority="FRACTION",
    amount={0: 0.25, 1: 0.40},
)
```

Inspect `sfr.network` before grid mapping and `sfr.reaches` afterward.
Multiple reaches may occupy the same groundwater cell without losing their
individual routing identities.

### Lakes

`LAKBuilder` combines stable lake identities, groundwater connections, period
settings, outlets, and optional lake tables:

```python
lak = LAKBuilder(
    context=context,
    nper=nper,
    lakes="lakes.gpkg",
    lake_id_field="lake_name",
    starting_stage="starting_stage",
    lake_bottom="bottom_elevation",
    bed_leakance="bed_leakance",
    connection_modes={
        "natural_lake": "automatic",
        "infiltration_trench": "rectangular",
        "managed_facility": (
            LakeConnection(
                cellid=(0, 42),
                connection_type="VERTICAL",
                bottom_elevation=95.0,
                top_elevation=95.0,
            ),
        ),
    },
    mover=True,
)

flow = flow.with_package(lak.build())
```

`"automatic"` builds natural-lake bottom and shoreline connections.
`"rectangular"` builds a box-like facility with bottom connections and
vertical sidewall connections across every face of its footprint cells.
A sequence of `LakeConnection` objects supplies explicit connections for one
lake. Every generated connection resolves one MF6 `bedleak` value; mappings
may provide different values by lake or by `(lake_id, connection_type)`.

Lake tables remain separate, reusable values:

```python
table = LakeTableBuilder(
    area=10_000.0,
    bottom=95.0,
    top=105.0,
    stage_step=0.5,
)

lak = lak.with_updates(tables={"natural_lake": table})
```

Tables may instead be derived from bathymetry:

```python
table = LakeTableBuilder(
    dem=Path("lake_bathymetry.tif"),
    footprint=Path("lake_footprint.gpkg"),
    stages=np.arange(95.0, 105.5, 0.5),
)
```

`LakeTableBuilder.build()` returns an immutable `LakeTable`. `LAKBuilder`
attaches configured tables as LAK child packages when its package spec builds.
Semantic `LakeOutlet` objects refer to stable source and receiver lake IDs
rather than generated lake numbers.

## GeoPackage Sources

Use package-first `.gpkg(...)` helpers for common GIS boundaries. They map
features through the grid carried by `ModelContext` and return package specs
directly:

```python
chd_package = chd.gpkg("inputs/chd.gpkg", context=context, nper=12, head=["head_0", "head_1"])
ghb_package = ghb.gpkg("inputs/ghb.gpkg", context=context, nper=12, head="head", conductance="conductance")
drn_package = drn.gpkg("inputs/drn.gpkg", context=context, nper=12, elevation="elevation", conductance="conductance")
wel_package = wel.gpkg("inputs/wel.gpkg", context=context, nper=12, rate="rate")
rch_package = rch.gpkg("inputs/rch.gpkg", context=context, nper=12, recharge="recharge")

source = GeoPackageSource("inputs/k_zones.gpkg", context, nper=12)
k = source.k_array(value="k", nlay=4, defaults=1.0)
```

For recharge that is not coming from GIS, use `rch(...)`. It selects the top
active cell in each active column by default:

```python
flat_rch = rch(
    context=context,
    nper=12,
    recharge=1.0e-4,
)

cell_rch = rch(
    context=context,
    nper=12,
    recharge={
        (0, 10): 1.0e-4,
        (0, 11): 1.2e-4,
    },
)

transient_cell_rch = rch(
    context=context,
    nper=2,
    recharge={
        0: {(0, 10): 1.0e-4},
        1: {(0, 10): 2.0e-4},
    },
)

direct_rch = rch.flopy(stress_period_data={0: [[(0, 10), 1.0e-4]]})
```

Boundary elevations and heads can also be derived from model surfaces:

```python
drn_package = drn.gpkg(
    "inputs/drn.gpkg",
    context=context,
    nper=12,
    elevation=CellSurfaceOffset(reference="cell_top", offset=-2.0),
    conductance="conductance",
)

ghb_package = ghb.gpkg(
    "inputs/ghb.gpkg",
    context=context,
    nper=12,
    head=CellSurfaceOffset(
        reference="cell_bottom",
        offset="height_above_bottom",
        minimum="min_elev",
    ),
    conductance="conductance",
)
```

Use `period_field` for long-format inputs containing one feature row per
stress period. The source filters inactive cells using `ModelContext.domain`.

The same structure supports GWE and PRT models. They belong beside GWF and GWT
inside `SimulationSpec`, rather than under a separate tracking subsystem.

## Projects and Runs

`Project` owns reusable specs and creates named `Run` workspaces:

```python
from myflopy import Project

project = Project("modeling/my_project")
project.add_simulation(simulation)

baseline = project.run("baseline", "coupled")
reopened = project.reopen_run("baseline")
```

A `Run` owns one concrete workspace and provides `build()`, `write()`,
`execute()`, and `open()`. Its `run.json` manifest records lifecycle status and
a readable summary of the simulation spec. Reopening uses FloPy 3.10's native
`MFSimulation.load`; Python builder callables are intentionally not serialized.

Project spec registries are therefore in-memory configuration. Project and run
manifests, MF6 input files, and results are durable and can be reopened later:

```python
project = Project.reopen("modeling/my_project")
baseline = project.reopen_run("baseline")
```

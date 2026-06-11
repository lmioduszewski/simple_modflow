# `simple_modflow` Preferred API

This document describes the preferred modern API for `simple_modflow`.

It is the counterpart to [codebase_structure.md](C:/Users/lukem/Python/Projects/simple_modflow/docs/codebase_structure.md):  
- `codebase_structure.md` tells you where code lives  
- `preferred_api.md` tells you how the package is intended to be used

The package still carries some compatibility surfaces from older workflows. Those
older names are not necessarily wrong, but new code should prefer the patterns
described here.

## Quick Resources

- [simple_modflow_api_pamphlet.pdf](C:/Users/lukem/Python/Projects/simple_modflow/docs/simple_modflow_api_pamphlet.pdf)
  A short visually organized workflow summary, arranged one topic per page.
- [targets_api_quickstart.ipynb](C:/Users/lukem/Python/Projects/simple_modflow/examples/mf6/notebooks/targets_api_quickstart.ipynb)
  A short notebook focused on the flexible target API, model-bound target registry, FloPy observation generation, and PEST touchpoints.
- [mp3du_quickstart.md](C:/Users/lukem/Python/Projects/simple_modflow/docs/mp3du_quickstart.md)
  Focused quickstart for the supported MP3DU particle-tracking API.

## Cheat Sheet

### Import

```python
import simple_modflow as mf
```

### Build a model

```python
model = mf.SimulationBase(
    name="my_model",
    mf_folder_path=workspace,
    vor=vor,
    nper=2,
)
```

### Build grid from geometry

```python
tri = mf.TriangleGrid(model_ws=str(workspace / "mesh"))
tri.set_domain_file(domain_path)
tri.add_line_feature(stream_path, buffer=50, max_area=3000, label="stream")
tri.add_region_file(lake_path, max_area=2500, label="lake")
tri.build_mesh(profile="balanced")

vor = mf.VoronoiGridPlus(tri)
```

### Attach core packages

```python
DisvGrid(...)
TemporalDiscretization(...)
InitialConditions(...)
KFlow(...)
Storage(...)
OutputControl(...)
```

### Build package data from vector inputs

```python
drn_dict = mf.DRNFromVector(model=model, vor=vor, shp_gpkg=drn_path, uid="name").from_vector()
ghb_dict = mf.GHBFromVector(model=model, vor=vor, shp_gpkg=ghb_path, uid="name").from_vector()
chd_dict = mf.CHDFromVector(model=model, vor=vor, shp_gpkg=chd_path, uid="name").from_vector()
rch_dict = mf.RCHFromVector(model=model, vor=vor, shp_gpkg=rch_path, uid="zone").from_vector()
k_array = mf.KFromVector(model=model, vor=vor, shp_gpkg=k_path, uid="name").from_vector(defaults=[10.0])
```

### Attach package wrappers

```python
CHD(model=model, stress_period_data=chd_dict)
mf.modflow.mf6.simulation.packages.GHB(model=model, stress_period_data=ghb_dict)
mf.modflow.mf6.simulation.packages.Drains(model=model, stress_period_data=drn_dict)
Recharge(model=model, vor=vor, rch_dict=rch_dict)
```

### UZF / LAK / SFR / MVR

```python
uzf = UZFPackageData(..., add_uzf=True)
lake_connections = LakeConnectionData(...)
connectiondata = lake_connections.connection_data
lake_packagedata = LakePackageData(...)
lake_perioddata = LakePeriodData(...)
LAKPackage(..., packagedata=lake_packagedata.packagedata, connectiondata=connectiondata, perioddata=lake_perioddata.perioddata, mover=True)

sfr = SFR(..., mover=True, add_sfr=True)
mvr_perioddata = {0: [["sfr", sfr.stream_reaches[0][-1] - 1, "lak", 0, "FACTOR", 0.25]]}
model.validate_surface_water(
    nlakes=1,
    lak_packagedata=lake_packagedata.packagedata,
    lak_connectiondata=connectiondata,
    lak_perioddata=lake_perioddata.perioddata,
    sfr=sfr,
    maxmvr=1,
    maxpackages=2,
    mvr_packages=[["sfr"], ["lak"]],
    mvr_perioddata=mvr_perioddata,
    raise_on_error=True,
)
MVR(model=model, maxmvr=1, maxpackages=2, packages=[["sfr"], ["lak"]], perioddata=mvr_perioddata)
```

### Run the model

```python
success, buff = model.run_simulation()
```

### Calibration / PEST first slice

```python
targets = mf.HeadTargets(
    locations=head_points_gpkg,
    values=head_targets_csv,
    name_column="name",
    layer_column="layer",
    time_column="per",
)

pest = mf.PestProject(
    model=model,
    name="cumberland_pest",
    workspace=model.workspace / "pest",
    start_datetime="2024-01-01",
)

pest.add_parameter(
    mf.KPilotPointParameter(
        name="hk",
        source=mf.VectorParameterSource(path=hk_gpkg, value_column="k", zone_column="unit"),
        bounds=(0.25, 4.0),
        bounds_mode="multiplier",
        transform="log",
        pp_spacing=1500.0,
        geostruct=mf.ExpGeoStruct(range=3000.0, transform="log"),
    )
)

pest.add_parameter(
    mf.DrainElevationParameter(
        name="drn_elev",
        source=mf.VectorParameterSource(
            path=drn_gpkg,
            value_column="elev",
            feature_id_column="name",
            group_column="group",
            layer_column="layer",
        ),
        bounds=(-10.0, 10.0),
        bounds_mode="absolute",
    )
)

pest.add_parameter(
    mf.DrainConductanceParameter(
        name="drn_cond",
        source=mf.VectorParameterSource(
            path=drn_gpkg,
            value_column="cond",
            feature_id_column="name",
            group_column="group",
            layer_column="layer",
        ),
        bounds=(0.25, 4.0),
        bounds_mode="multiplier",
        transform="log",
    )
)

pest.add_observation(mf.HeadTargetObservationSpec(targets=targets))
pst = pest.build_pst("first_slice.pst")
```

### Outputs

```python
model.hds
model.hds.long()
model.hds.wide()
model.hds.map(contours=True)
model.hds.map(contours="top")
model.hds.map(contours="bottom", contour_levels=10)
model.hds.map(contours=True, contour_method="linear")
model.hds.map(contours=True, contour_method="cubic", contour_resolution=150)
model.hds.map(contours=True, contour_levels=[100, 105, 110])
model.bud("drn")
model.packages.rch.inputs.get()
model.packages.rch.inputs.map(per=0)
model.packages.wel.inputs.q.get()
model.packages.wel.inputs.q.map(per=0)
model.packages.ic.strt.get()
model.packages.ic.strt.long()
model.packages.ic.strt.map(layer=0)
model.packages.npf.k.get()
model.packages.npf.k.wide()
model.packages.npf.k.map(layer=0)
model.packages.npf.k.map(layer=0, contours=True)
model.packages.npf.k.map(layer=0, contours=True, contour_method="cubic")
model.packages.sto.ss.get()
model.packages.sto.sy.get()
model.packages.sto.ss.map(layer=0)
model.packages.uzf.inputs.finf.get()
model.packages.uzf.inputs.finf.map(per=0)
model.packages.uzf.inputs.fields
model.packages.uzf.inputs.summary()
model.packages.uzf.inputs.pet.get()
model.packages.uzf.inputs.rootact.map(per=0)
model.packages.chd.results.q.get()
model.packages.chd.results.q.wide()
model.packages.chd.results.q.long()
model.packages.chd.results.q.stack()
model.packages.chd.results.q.plot_timeseries(cells=[1, 3])
model.packages.chd.results.q.map(per=0)
model.packages.uzf.results.gwrch.get()
model.packages.uzf.results.gwrch.wide()
model.packages.uzf.results.gwrch.long()
model.packages.uzf.results.gwrch.plot_timeseries(cells=[0, 1])
model.packages.uzf.results.gwrch.map(per=0)
model.packages.uzf.results.fields
model.packages.uzf.results.summary()
model.packages.uzf.results.sat.get()
model.packages.uzf.results.sat.plot_timeseries(cells=[0, 1])
model.packages.uzf.results.sat.map(per=0)
model.packages.lak.connections.get()
model.packages.lak.connections.map()
model.packages.lak.budget.types
model.packages.lak.budget.get()
model.packages.lak.budget.get(term="GWF", per=0)
model.packages.lak.budget.summary()
model.packages.lak.budget.wide(index=["per", "lake"])
model.packages.lak.budget.gwf.get(per=0)
model.packages.lak.budget.storage.get()
model.packages.lak.budget.runoff.get()
model.packages.lak.budget.rainfall.get()
model.packages.lak.budget.evaporation.get()
model.packages.lak.budget.withdrawal.get()
model.packages.lak.budget.constant.get()
model.packages.lak.budget.ext_inflow.get()
model.packages.lak.budget.ext_outflow.get()
model.packages.lak.budget.from_mvr.get()
model.packages.lak.budget.to_mvr.get()
model.packages.lak.budget.flow_ja_face.get()
model.packages.lak.budget.mvr.get()
model.packages.lak.budget.lake_fluxes.get()
model.packages.lak.results.stage.get()
model.packages.lak.results.stage.map(per=0)
model.packages.lak.results.stage.plot_timeseries()
model.packages.lak.results.stage_change.get()
model.packages.lak.results.stage_change.plot_timeseries()
model.packages.lak.results.q.get()
model.packages.lak.results.q.map(per=0)
model.packages.lak.results.q.map(per=0, connection_type="VERTICAL")
model.packages.lak.results.q.map(per=0, connection_type="HORIZONTAL")
model.packages.lak.results.q.budget_summary(per=0)
model.packages.lak.results.q.plot_budget(per=0)
model.packages.surface_water.results.q.get()
model.packages.surface_water.results.q.map(per=0)
model.packages.sfr.budget.types
model.packages.sfr.budget.get()
model.packages.sfr.budget.get(term="GWF", per=0)
model.packages.sfr.budget.summary()
model.packages.sfr.budget.wide(index=["per", "reach"])
model.packages.sfr.budget.gwf.get(per=0)
model.packages.sfr.budget.flow_ja_face.get()
model.packages.sfr.budget.ext_inflow.get()
model.packages.sfr.budget.runoff.get()
model.packages.sfr.budget.rain.get()
model.packages.sfr.budget.evaporation.get()
model.packages.sfr.budget.ext_outflow.get()
model.packages.sfr.budget.storage.get()
model.packages.sfr.budget.from_mvr.get()
model.packages.sfr.budget.to_mvr.get()
model.packages.sfr.budget.mvr.get()
model.packages.sfr.budget.stream_fluxes.get()
model.packages.sfr.results.stage.get()
model.packages.sfr.results.stage.map(per=0)
model.packages.sfr.results.stage.profile(per=0)
model.packages.sfr.results.stage.plot_profile(per=0)
model.packages.sfr.results.q.get()
model.packages.sfr.results.q.map(per=0)
model.packages.sfr.results.q.profile(per=0)
model.packages.sfr.results.q.plot_profile(per=0)
model.packages.sfr.results.long_profile(per=0)
model.packages.sfr.results.plot_long_profile(per=0)
model.outputs.lak.stage.get()
model.outputs.sfr.stage.get()
```

### Project-managed runs

```python
catalog = mf.ProjectCatalog(project_root)

model = mf.SimulationBase(
    vor=vor,
    nper=1,
    project_catalog=catalog,
    run_id="source_run",
    model_spec="base_model",
)
```

### Reopen existing runs lazily

```python
archive = mf.explore_runs(root)
run = archive.open("run_id")

run.summary()
run.vor
run.hds
run.bud("drn")
run.packages.rch.inputs.get()
run.packages.rch.inputs.map(per=0)
run.packages.uzf.results.gwrch.get()
run.packages.uzf.results.gwrch.wide()
run.packages.uzf.results.gwrch.plot_timeseries(cells=[0, 1])
run.packages.uzf.results.gwrch.map(per=0)
run.packages.uzf.results.sat.get()
run.packages.sfr.results.q.profile(per=0)
run.packages.sfr.results.q.plot_profile(per=0)
run.packages.sfr.results.long_profile(per=0)
run.packages.sfr.results.plot_long_profile(per=0)
run.outputs.lak.stage.get()
run.uzf
run.load_all()
```

### Compare multiple models

```python
group = mf.ModelGroup(
    {
        "predev": pre_path,
        "postdev": post_path,
    },
    reference="predev",
    shared_grid=True,
)

group.hds.get()
group.hds.compare()
group.bud("rch").get()
group.bud("rch").compare()
group.packages.rch.inputs.get()
group.packages.rch.inputs.compare()
group.packages.rch.inputs.map(model_name="postdev", per=0)
group.packages.rch.inputs.compare_map(model_name="postdev", per=0)
group.packages.uzf.inputs.finf.get()
group.packages.uzf.inputs.finf.compare()
group.packages.uzf.inputs.finf.map(model_name="postdev", per=0)
group.packages.uzf.inputs.finf.compare_map(model_name="postdev", per=0)
group.packages.chd.results.q.get()
group.packages.chd.results.q.compare()
group.packages.chd.results.q.map(model_name="postdev", per=0)
group.packages.chd.results.q.compare_map(model_name="postdev", per=0)
group.packages.chd.results.q.subplot_map(per=0)
group.packages.chd.results.q.plot_timeseries(cells=[1, 3])
group.packages.uzf.results.gwrch.get()
group.packages.uzf.results.gwrch.compare()
group.packages.uzf.results.gwrch.subplot_map(per=0)
group.packages.uzf.results.gwrch.plot_timeseries(cells=[0, 1])
group.packages.uzf.results.sat.get()
group.packages.lak.connections.get()
group.packages.lak.connections.map(model_name="postdev")
group.packages.lak.results.stage.get()
group.packages.lak.results.stage.compare()
group.packages.lak.results.stage.plot_timeseries()
group.packages.lak.results.q.get()
group.packages.lak.results.q.compare()
group.packages.lak.results.q.compare_map(model_name="postdev", per=0)
group.packages.lak.results.q.compare_map(model_name="postdev", per=0, connection_type="HORIZONTAL")
group.packages.lak.results.q.subplot_map(per=0)
group.packages.sfr.results.q.get()
group.packages.sfr.results.q.compare()
group.packages.sfr.results.q.compare_map(model_name="postdev", per=0)
group.packages.sfr.results.q.subplot_map(per=0)
group.packages.surface_water.results.q.subplot_map(per=0)
group.rch.get()
group.rch.compare()
group.outputs.lak.stage()
```

## Design Goals

The preferred API is built around a few principles:

- Keep FloPy as the low-level engine.
- Use `simple_modflow` for higher-level workflow APIs.
- Make model creation, run management, comparison, and inspection feel natural.
- Keep indexing exposed by the public API zero-based.
- Prefer lazy loading when reopening existing MF6 runs from disk.
- Keep package-building helpers organized by role:
  - model object
  - vector/material builders
  - package wrappers
  - outputs
  - project/run/artifact helpers

## Public Entry Points

The top-level import should usually be:

```python
import simple_modflow as mf
```

The most important top-level objects are:

- `mf.SimulationBase`
- `mf.TriangleGrid`
- `mf.VoronoiGridPlus`
- `mf.ProjectCatalog`
- `mf.ModelGroup`
- `mf.explore_runs(...)`

Preferred vector/material builders:

- `mf.DRNFromVector`
- `mf.GHBFromVector`
- `mf.CHDFromVector`
- `mf.RCHFromVector`
- `mf.KFromVector`

Surface-water validation:

- `mf.validate_surface_water_configuration(...)`

## Core Mental Model

The preferred workflow revolves around four layers:

1. `SimulationBase`
   The main live model object for building and running MF6 models.

2. Grid builders
   `TriangleGrid` and `VoronoiGridPlus` are the preferred path for building and
   working with refined Voronoi-style grids.

3. Project/run helpers
   `ProjectCatalog`, package artifacts, and related-run helpers are used when
   you want provenance, reuse, and scenario management.

4. Lazy loaded runs
   `explore_runs(...)` and `LoadedMf6Run` let you reopen MF6 directories
   without paying the full cost up front.

## Building a New Model

The preferred way to build a new model is:

```python
model = mf.SimulationBase(
    name="my_model",
    mf_folder_path=workspace,
    vor=vor,
    nper=2,
)
```

Then attach core model pieces explicitly:

```python
from simple_modflow.modflow.mf6.simulation.discretization import DisvGrid, TemporalDiscretization
from simple_modflow.modflow.mf6.simulation.packages import (
    InitialConditions,
    KFlow,
    Storage,
    OutputControl,
)

DisvGrid(...)
TemporalDiscretization(...)
InitialConditions(...)
KFlow(...)
Storage(...)
OutputControl(...)
```

Why this is preferred:

- `SimulationBase` is the central object used for both live builds and most
  model-facing exploration code.
- core package wrappers keep your defaults and naming conventions consistent.
- the model object stays compatible with project/catalog workflows.

## Grid Workflow

Preferred grid path:

1. Build or refine geometry with `TriangleGrid`
2. Convert to `VoronoiGridPlus`
3. Attach that grid to `SimulationBase`

Typical flow:

```python
tri = mf.TriangleGrid(...)
tri.set_domain_...
tri.add_region_...
tri.add_line_feature(...)
tri.build_mesh(profile="balanced")

vor = mf.VoronoiGridPlus(tri)
model = mf.SimulationBase(vor=vor, ...)
```

Use `TriangleGrid` when you want:

- file-driven domain setup
- code-defined domain setup
- refinement zones
- line-feature refinement
- mesh cleanup and optimization

## Preferred Vector and Material Builders

When package or property data comes from shapefiles or geopackages, prefer the
dedicated builder classes.

### Boundary Builders

- `DRNFromVector`
- `GHBFromVector`
- `CHDFromVector`
- `RCHFromVector`

Typical pattern:

```python
builder = mf.DRNFromVector(model=model, vor=vor, shp_gpkg=path, uid="name")
drn_dict = builder.from_vector(...)
```

Then attach the actual MF6 package with the package wrapper:

```python
from simple_modflow.modflow.mf6.simulation.packages import Drains

Drains(model=model, stress_period_data=drn_dict)
```

### Material Property Builders

For vector-driven hydraulic-property assignment, prefer:

- `KFromVector`

Typical pattern:

```python
k_builder = mf.KFromVector(model=model, vor=vor, shp_gpkg=k_path, uid="name")
k_array = k_builder.from_vector(defaults=[18.0])
```

## Package Wrappers

The thin wrappers in `simulation/packages.py` are still preferred for common
package attachment because they encode package defaults and integrate with the
artifact workflow.

Common wrappers:

- `InitialConditions`
- `KFlow`
- `Storage`
- `OutputControl`
- `Recharge`
- `Drains`
- `GHB`
- `CHD`
- `LAK`
- `UZF`
- `MVR`

These wrappers are not meant to replace FloPy conceptually. They are meant to:

- keep your MF6 package defaults consistent
- keep filenames and naming conventions predictable
- optionally capture reusable package artifacts

## Surface-Water Workflow

For coupled lake/stream/mover models, the preferred sequence is:

1. Build `LAK` inputs
2. Build and validate `SFR`
3. Validate the combined `LAK/SFR/MVR` configuration
4. Attach `MVR`
5. Run the model

Typical shape:

```python
lake_connections = LakeConnectionData(...)
connectiondata = lake_connections.connection_data
lake_packagedata = LakePackageData(...)
lake_perioddata = LakePeriodData(...)

LAKPackage(...)

sfr = SFR(..., add_sfr=True, mover=True)

report = model.validate_surface_water(
    nlakes=1,
    lak_packagedata=lake_packagedata.packagedata,
    lak_connectiondata=connectiondata,
    lak_perioddata=lake_perioddata.perioddata,
    sfr=sfr,
    maxmvr=1,
    maxpackages=2,
    mvr_packages=[["sfr"], ["lak"]],
    mvr_perioddata=mvr_perioddata,
    raise_on_error=True,
)

MVR(...)
```

### Validation

Preferred validation entry points:

- `model.validate_surface_water(...)`
- `mf.validate_surface_water_configuration(...)`

The validation layer is intended to catch common issues early:

- non-zero-based or out-of-range ids
- invalid SFR reach lengths/widths/gradients
- mismatched SFR connection counts
- inconsistent LAK packagedata/connectiondata
- MVR references to packages or ids that do not exist

Use the report directly in notebooks:

```python
report.summary_frame()
report.to_frame()
```

## Project, Runs, and Artifacts

Use `ProjectCatalog` when you want run tracking and reusable package pieces.

Preferred pattern:

1. Create the model with `SimulationBase`
2. Pass run metadata directly into `SimulationBase`
3. Optionally create package artifacts during package construction
4. Reuse those artifacts in related runs

Typical example:

```python
catalog = mf.ProjectCatalog(project_root)

model = mf.SimulationBase(
    vor=vor,
    nper=1,
    project_catalog=catalog,
    run_id="source_run",
    model_spec="base_model",
)
```

Then define packages normally, optionally with `artifact_id=...`.

Related runs can reference stored package artifacts through:

- `package_versions`
- `apply_registered_package_artifacts()`
- automatic artifact application before `run_simulation()` when configured

Use this layer when you want:

- scenario management
- provenance
- reusable package definitions
- related-run comparisons

## Reopening Existing MF6 Runs

Preferred entry point:

```python
archive = mf.explore_runs(root)
run = archive.open("run_id")
```

The returned run object is designed to feel like a `SimulationBase`-style model
object, but lazily.

Preferred usage:

```python
run.summary()
run.vor
run.hds
run.bud("drn")
run.outputs.lak.stage.get()
run.uzf
run.load_all()
```

Important behavior:

- `summary()` should stay cheap
- `vor` loads the grid lazily
- `hds` and budget accessors load outputs lazily
- package properties like `run.uzf` load packages lazily
- `load_all()` forces the full FloPy load

## Outputs vs Package Access

Keep this distinction in mind:

- `run.uzf`, `run.rch`, `run.sfr`
  Package access

- `run.packages.rch.inputs`, `run.packages.chd.inputs`, `run.packages.uzf.inputs.finf`
  Preferred normalized input exploration access

- `run.packages.chd.results.q`, `run.packages.uzf.results.gwrch`,
  `run.packages.uzf.results.sat`, `run.packages.lak.results.stage`,
  `run.packages.sfr.results.stage`, `run.packages.sfr.results.q`
  Preferred normalized result exploration access

- `run.outputs.uzf`, `run.outputs.lak`, `run.outputs.sfr`
  Output access

Examples:

```python
run.packages.rch.inputs.get()
run.packages.rch.inputs.map(per=0)
run.packages.wel.inputs.q.get()
run.packages.wel.inputs.q.map(per=0)
run.packages.uzf.inputs.finf.get()
run.packages.uzf.inputs.finf.map(per=0)
run.packages.uzf.inputs.pet.get()
run.packages.uzf.inputs.rootact.map(per=0)
run.packages.chd.results.q.get()
run.packages.chd.results.q.wide()
run.packages.chd.results.q.long()
run.packages.chd.results.q.map(per=0)
run.packages.uzf.results.gwrch.get()
run.packages.uzf.results.gwrch.stack()
run.packages.uzf.results.gwrch.map(per=0)
run.packages.uzf.results.sat.get()
run.packages.uzf.results.sat.map(per=0)
run.packages.lak.results.stage.get()
run.packages.sfr.results.stage.get()
run.packages.sfr.results.stage.profile(per=0)
run.packages.sfr.results.stage.plot_profile(per=0)
run.packages.sfr.results.q.get()
run.packages.sfr.results.q.map(per=0)
run.packages.sfr.results.q.profile(per=0)
run.packages.sfr.results.q.plot_profile(per=0)
run.packages.sfr.results.long_profile(per=0)
run.packages.sfr.results.plot_long_profile(per=0)
run.outputs.lak.stage.get()
run.outputs.sfr.stage.get()
run.outputs.uzf.some_output_helper
```

For SFR specifically, the preferred higher-level overview is:

```python
run.packages.sfr.results.long_profile(per=0)
run.packages.sfr.results.plot_long_profile(per=0)
```

That merged profile aligns stream stage, stream-groundwater exchange, and
streambed geometry along cumulative stream distance, which tends to be more
useful in notebooks than looking at stage and exchange separately.

SFR exchange maps are normalized by stream length by default. In other words,
`model.packages.sfr.results.q.map(...)` shows `sum(q) / sum(rlen)` within each
cell rather than raw `q`, which makes choropleths less biased toward cells with
longer reaches.

This split is intentional:

- package access is about inputs/attached package objects
- `run.packages...` is about tidy input/result tables and choropleth-ready maps
- output access is about result files and package-specific outputs

## Model Groups

Use `ModelGroup` when you want multi-model comparisons with a consistent API.

Preferred creation patterns:

```python
group = mf.ModelGroup(
    {
        "predev": pre_model,
        "postdev": post_model,
    },
    reference="predev",
)
```

or directly from run directories:

```python
group = mf.ModelGroup(
    {
        "predev": pre_path,
        "postdev": post_path,
    },
    reference="predev",
    shared_grid=True,
)
```

Preferred grouped access:

- `group.hds.get()`
- `group.hds.compare()`
- `group.bud("drn").get()`
- `group.bud("drn").compare()`
- `group.packages.rch.inputs.get()`
- `group.packages.rch.inputs.compare()`
- `group.packages.rch.inputs.map(model_name="scenario", per=0)`
- `group.packages.rch.inputs.compare_map(model_name="scenario", per=0)`
- `group.packages.wel.inputs.get()`
- `group.packages.wel.inputs.compare()`
- `group.packages.uzf.inputs.finf.get()`
- `group.packages.uzf.inputs.finf.compare()`
- `group.packages.uzf.inputs.pet.get()`
- `group.packages.uzf.inputs.pet.compare()`
- `group.packages.uzf.inputs.finf.map(model_name="scenario", per=0)`
- `group.packages.uzf.inputs.finf.compare_map(model_name="scenario", per=0)`
- `group.packages.chd.results.q.get()`
- `group.packages.chd.results.q.compare()`
- `group.packages.chd.results.q.map(model_name="scenario", per=0)`
- `group.packages.chd.results.q.compare_map(model_name="scenario", per=0)`
- `group.packages.uzf.results.gwrch.get()`
- `group.packages.uzf.results.gwrch.compare()`
- `group.packages.uzf.results.sat.get()`
- `group.packages.sfr.results.q.get()`
- `group.packages.sfr.results.q.compare()`
- `group.packages.sfr.results.q.compare_map(model_name="scenario", per=0)`
- `group.rch.get()`
- `group.rch.compare()`
- `group.uzf.finf.get()`
- `group.uzf.finf.compare()`
- `group.outputs.lak.stage()`

Mental model:

- `get()` returns aligned raw values across all models
- `compare()` returns differences relative to the group reference model
- `map()` renders one model's raw input values on the grid
- `compare_map()` renders one model's difference from the reference model

## Preferred Naming

Preferred builder names:

- `DRNFromVector`
- `GHBFromVector`
- `CHDFromVector`
- `RCHFromVector`
- `KFromVector`

Older compatibility names still exist in places, for example:

- `DRN`
- `GHB`
- `RechargeFromShp`

New code should prefer the `*FromVector` names.

## Compatibility and Legacy Notes

Some older helpers remain for compatibility.

That is expected.

The preferred direction is:

- `SimulationBase` as the main live model object
- vector/material builders for GIS-driven inputs
- package wrappers for model attachment
- `ProjectCatalog` for managed run workflows
- `LoadedMf6Run`/`RunExplorer` for lazy reopening
- `ModelGroup` for multi-model comparison

If you are deciding between an older helper and one of the patterns above,
prefer the modern pattern unless you are specifically maintaining legacy code.

## Best Example Notebooks

The most useful notebooks for the preferred API are:

- [project_artifact_workflow.ipynb](C:/Users/lukem/Python/Projects/simple_modflow/examples/mf6/notebooks/project_artifact_workflow.ipynb)
  Project-managed artifact workflow.
- [optimized_related_run_artifact_workflow.ipynb](C:/Users/lukem/Python/Projects/simple_modflow/examples/mf6/notebooks/optimized_related_run_artifact_workflow.ipynb)
  Related runs on the same optimized mesh.
- [refined_end_to_end_preferred_api_workflow.ipynb](C:/Users/lukem/Python/Projects/simple_modflow/examples/mf6/notebooks/refined_end_to_end_preferred_api_workflow.ipynb)
  Refined single-model workflow with preferred vector builders plus validated
  `UZF/LAK/SFR/MVR`.
- [pest_first_slice_workflow.ipynb](C:/Users/lukem/Python/Projects/simple_modflow/examples/mf6/notebooks/pest_first_slice_workflow.ipynb)
  First reusable PEST/pyEMU workflow with GIS-defined `K`, GIS-defined drains,
  `HeadTargets`, and forward-run parameter application.
- [cumberland_predev_pest_first_slice_workflow.ipynb](C:/Users/lukem/Python/Projects/simple_modflow/examples/mf6/notebooks/cumberland_predev_pest_first_slice_workflow.ipynb)
  Cumberland-specific first-slice calibration example using the current
  pre-development `K` and drain GIS inputs on one steady-state stress period.
- [cumberland_observed_snapshot_pest_first_slice_workflow.ipynb](C:/Users/lukem/Python/Projects/simple_modflow/examples/mf6/notebooks/cumberland_observed_snapshot_pest_first_slice_workflow.ipynb)
  Cumberland-specific first-slice calibration example that uses a real
  observed-head snapshot from `calib_observations.xlsx` instead of synthetic
  pseudo-targets.
- [cumberland_steady_snapshot_pest_first_slice.py](C:/Users/lukem/Python/Projects/simple_modflow/examples/mf6/cumberland_steady_snapshot_pest_first_slice.py)
  Command-line Cumberland steady-state snapshot calibration example using the
  current first-slice `PestProject` API with one steady-state stress period, a
  representative observed-head snapshot, low-weight supplemental wells,
  timestamped artifact workspaces, and optional fast/parallel run modes.
- [cumberland_transient_observed_pest_first_slice.py](C:/Users/lukem/Python/Projects/simple_modflow/examples/mf6/cumberland_transient_observed_pest_first_slice.py)
  Command-line Cumberland transient calibration example using the current
  first-slice `PestProject` API with one initial steady-state period, 12 monthly
  transient periods, real observed heads, low-weight one-time supplemental
  wells, timestamped artifact workspaces, and optional fast/parallel run modes.
- [import_existing_runs_workflow.ipynb](C:/Users/lukem/Python/Projects/simple_modflow/examples/mf6/notebooks/import_existing_runs_workflow.ipynb)
  Lazy run reopening and archive browsing.

## Interactive Results and Particle Tracking

Use `model.visualize` for the canonical interactive result-export surface:

```python
model.visualize.cross_section_slider_html(line, "cross_section.html")
model.visualize.head_map_slider_html("head_map.html", layer=0)
model.visualize.head_layer_mosaic_slider_html(
    "head_mosaic.html",
    dpi=240,
    panel_figsize=(7, 5),
)

## Unified Parallel Model Splitting

Use the same API for a single GWF model or a coupled GWF-GWT/GWE simulation.
`simple_modflow` inspects the simulation and selects FloPy's appropriate
splitting implementation internally.

```python
parallel = model.parallel.split_model(
    workspace="parallel_run",
    nparts=8,
    active_only=True,
)

parallel.plot_partitions()
parallel.summary()

# Validate the split simulation before introducing MPI.
parallel.run_serial()
parallel.compare_heads()

# Requires mpiexec and a parallel-enabled MODFLOW 6 executable.
parallel.run(processors=8)

# Reconstruct partitioned results onto the original model grid.
heads = parallel.results.heads()
```

For the shortest workflow, preparation and execution can be combined:

```python
parallel = model.parallel.run(
    workspace="parallel_run",
    nparts=8,
    processors=8,
    validate_serial=True,
)
```

For a hydrologically informed custom partition, pass an integer mask instead of
`nparts`. Automatic masks require the optional `parallel` dependencies:

```bash
python -m pip install "simple_modflow[parallel]"
```

model.visualize.plotly_cross_section_animation(
    cells=[10, 20, 30],
    output_path="cross_section_plotly.html",
)
model.visualize.plotly_head_map_animation(
    layer=0,
    output_path="head_map_plotly.html",
)
```

Use MF6 PRT as the preferred integrated particle-tracking engine:

```python
release_points = mf.PRTReleasePoints.from_cells(model, [100, 120, 140])
prt = model.particle_tracking.prt(
    workspace=model.workspace.parent / "prt",
    release_points=release_points,
)
result = prt.run()

result.plot_map()
result.export_3d_html("prt_pathlines.html")
```

MP3DU remains available through `model.particle_tracking.mp3du(...)`.

See:

- [interactive_visualization_and_prt.md](C:/Users/lukem/Python/Projects/simple_modflow/docs/interactive_visualization_and_prt.md)
- [princeton_2026_visualization_review.md](C:/Users/lukem/Python/Projects/simple_modflow/docs/princeton_2026_visualization_review.md)

## Where To Look Next

- For file locations: [codebase_structure.md](C:/Users/lukem/Python/Projects/simple_modflow/docs/codebase_structure.md)
- For notebook entry points: [examples/mf6/README.md](C:/Users/lukem/Python/Projects/simple_modflow/examples/mf6/README.md)
`model.packages.lak.results.q.map(...)` uses lake-groundwater exchange normalized by
total lake connection flow area in each cell, so the mapped quantity is a
signed length-per-time exchange intensity instead of raw volumetric `q`.
By default it sums both vertical and horizontal lake connections within each
model cell. Pass `connection_type="VERTICAL"` or `connection_type="HORIZONTAL"`
to isolate one connection family.

`model.packages.surface_water.results.q.map(...)` combines SFR and LAK on one
shared map. It converts both packages to one physical sign convention before
plotting:

- positive = groundwater gaining into the surface-water feature
- negative = surface-water losing to groundwater

This means SFR exchange is internally sign-flipped relative to the raw MF6
`SFR/GWF` budget term so the combined map can share one `L/T` color scale with
LAK.

Grouped LAK exchange comparison maps use the same area-normalized quantity, and
diverging signed maps keep `0` centered in white with symmetric color limits.

For grouped scenario review, the preferred quick-look panels are:

- `group.packages.lak.results.q.subplot_map(...)`
- `group.packages.sfr.results.q.subplot_map(...)`
- `group.packages.surface_water.results.q.subplot_map(...)`

These draw one map per model and keep a shared symmetric `L/T` color range
across all panels.

`model.packages.lak.connections.map(...)` maps connection geometry such as total
connection area by cell, while `model.packages.lak.results.stage.plot_timeseries()`,
`model.packages.lak.results.stage_change.plot_timeseries()`, and
`model.packages.lak.results.q.plot_budget(...)` provide quick lake-centric
views that complement the choropleths. For grouped workflows,
`group.packages.lak.results.stage.plot_timeseries()` is the preferred way to
compare lake stage trajectories across scenarios.

`model.packages.sfr.results.q.map(...)` uses MF6 `GWF` exchange directly, where
positive values mean flow from the stream reach to groundwater. The default map
colors are therefore reversed so gaining reaches plot blue and losing reaches
plot red.

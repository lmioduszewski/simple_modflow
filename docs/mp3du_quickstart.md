# MP3DU Quickstart

The supported particle-tracking API in `myflopy` is:

- `ParticleTrackingInput`
- `prepare_particle_tracking`
- `run_particle_tracking`

Those are the only recommended entry points for new MP3DU workflows.

## Recommended One-Liner

```python
from pathlib import Path
import pickle
import myflopy as mf

with open(Path(r"C:\path\to\model.model"), "rb") as f:
    model = pickle.load(f)

result = mf.run_particle_tracking(
    model=model,
    particles=Path(r"C:\path\to\particle_selection.gpkg"),
    porosity=0.2,
    output_path=Path(r"C:\path\to\mp3du_run"),
)
```

`particles=` accepts either:

- a point shapefile or geopackage with MP3DU particle fields already present
- a polygon shapefile or geopackage with no particle attributes, in which case `myflopy` maps geometry to model cells and generates the required particle file
- a list of zero-based model cell IDs

## Zero-Based Cell IDs

```python
result = mf.run_particle_tracking(
    model=model,
    particles=[9151, 9152, 9747, 10347],
    porosity=0.2,
    output_path=Path(r"C:\path\to\mp3du_run"),
)
```

Generated MP3DU start files preserve the original zero-based IDs in a `Cell0` field and write one-based `P3D_CellID` values for MP3DU itself.

## Prepare First, Run Later

Use `prepare_particle_tracking(...)` when you want the configured tracker plus generated inputs before deciding whether to execute MP3DU.

```python
tracker, result = mf.prepare_particle_tracking(
    model=model,
    particles=Path(r"C:\path\to\particle_selection.gpkg"),
    porosity=0.2,
    output_path=Path(r"C:\path\to\mp3du_run"),
    execute=False,
    convert_output=False,
)
```

Then run explicitly:

```python
final_result = tracker.run(execute=True, convert_output=True)
```

## Full Control

```python
from myflopy.modflow.mp3du import ParticleTrackingInput

tracker = ParticleTrackingInput(
    model=model,
    particle_shp=Path(r"C:\path\to\particle_selection.gpkg"),
    porosities_by_layer=[0.2 for _ in range(model.gwf.modelgrid.nlay)],
    output_path=Path(r"C:\path\to\mp3du_run"),
    direction="FORWARD",
    simulation_end_time=5000.0,
    iface_overrides={"DRN": 7, "RCH": 6, "SFR": 6},
)
result = tracker.run()
```

## Outputs

`run_particle_tracking(...)` and `tracker.run(...)` return a `ParticleTrackingResult` with:

- `json_file`
- `path_file`
- `output_json`
- `start_cell_diagnostics`
- `endpoint_summary`
- `diagnostics_file`

The result still supports dictionary-style access for backward compatibility:

```python
result["endpoint_summary"]
```

## Known Limitations

- Some MP3DU `IFACE` package names are not accepted by the current executable build, even if the MF6 model contains those packages.
- Particles that enter sink/source-controlled flow paths may terminate quickly. `endpoint_summary` is the first place to check why.
- Legacy FloPy helpers for the native MODFLOW `PRT` model still exist in `myflopy.modflow.mp3du.legacy_prt`, but they are not part of the supported MP3DU API and remain constrained by current MODFLOW PRT support for many DISV grids.

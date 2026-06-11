# Interactive Visualization and MF6 PRT

`simple_modflow` supports three complementary interactive result workflows:

1. standalone Matplotlib frame sliders
2. Plotly animations with sliders
3. FloPy VTK / PyVista 3D particle scenes

The Matplotlib and PyVista HTML exports are standalone files. They do not need
JupyterLab or a running Python kernel after export.

## Install 3D Support

```powershell
python -m pip install -e ".[viz3d]"
```

The `viz3d` extra installs PyVista and the Trame packages needed by
`Plotter.export_html()`.

## Model-Bound Visualization API

Use `model.visualize` as the preferred entry point:

```python
# Existing simple_modflow/FloPy Matplotlib cross-section style, through time.
model.visualize.cross_section_slider_html(
    line,
    "cross_section.html",
)

# FloPy PlotMapView head maps through time.
model.visualize.head_map_slider_html(
    "head_map.html",
    layer=0,
)

# Multi-layer head mosaics through time.
model.visualize.head_layer_mosaic_slider_html(
    "head_mosaic.html",
    layers=[0, 1, 2],
    ncols=3,
    dpi=240,
    panel_figsize=(7, 5),
)

# Existing Plotly cross-section animation.
model.visualize.plotly_cross_section_animation(
    cells=[10, 20, 30],
    output_path="cross_section_plotly.html",
)

# Existing Plotly map animation.
model.visualize.plotly_head_map_animation(
    layer=0,
    output_path="head_map_plotly.html",
    frame_stride=2,
    max_frames=12,
    zmin=100,
    zmax=125,
    include_plotlyjs="cdn",
)
```

Layer mosaics default to `200` DPI. Use `dpi=` for the most direct resolution
control. `panel_figsize=(width, height)` controls each layer panel's size in
inches, so final pixel dimensions are approximately
`panel_figsize * rows/columns * dpi`. Higher DPI and larger panels produce
sharper maps but increase render time, memory use, and HTML or PNG size.

Plotly cross-section animations preserve the current interactive x- and y-axis
zoom while playing or selecting frames. Double-click the plot to deliberately
reset the axes to their initial ranges. Plotly HTML exports also inherit the
standard `figs` configuration, including scroll-wheel zoom and a hidden Plotly
modebar logo. Pass `config={...}` to override or extend those defaults.

Plotly maps calculate an initial center and zoom from padded WGS84 model bounds,
but do not constrain navigation; users can zoom out or pan anywhere. Use
`fit_bounds=False, zoom=13` to supply the initial zoom directly, or adjust the
automatically calculated home view with `bounds_padding=0.1`.

Standalone Plotly head-map HTML exports use a trace-only replay controller:
each frame restyles the existing map trace's values and hover data while the
layout and global coloraxis remain unchanged. This avoids the colorbar flashing
and unreliable second-play behavior seen when Plotly's native map-frame engine
redraws choropleth traces. The returned Python figure retains normal Plotly
frames for notebook inspection; use `output_path=...` for the repeatable
standalone animation.

Plotly frame selection is applied before traces are built, and the DISV GeoJSON
is stored only on the base trace rather than repeated in every animation frame.
Each animation frame explicitly updates that base trace with its own result
values and hover table. Standard head-map hover includes time step, stress
period, cell number, coordinates, heads by layer, model top, and layer bottoms.
The shared colorscale and colorbar live on a persistent layout-level
`coloraxis`, while map view remains locked so user pan and zoom are preserved.
For renderer compatibility and reliable repeated playback, each frame carries a
complete choropleth trace, including geometry and hover metadata. Provide
`zmin` and `zmax` to set the fixed global colorscale range explicitly; otherwise
it is calculated across the selected animation frames.
For large standalone exports, `include_plotlyjs="cdn"` avoids embedding the
Plotly JavaScript library while preserving the model geometry and results.

The standalone Matplotlib slider exporter pre-renders figures into PNG frames
and embeds those frames plus a browser-side slider into one HTML file. This is
the standalone equivalent of an `ipywidgets.Image` plus `IntSlider` callback.

For large models, prefer resumable external frames:

```python
model.visualize.head_map_slider_html(
    "head_map.html",
    embed_frames=False,
    frame_stride=2,
    max_frames=12,
    resume=True,
    progress=True,
)
```

External-frame exports use atomic file writes, a frame manifest, and an
exclusive output lock. Interrupted exports can reuse completed frames with
`resume=True`; concurrent processes are prevented from writing the same export.

## Canonical MF6 PRT Workflow

MF6 PRT is the preferred integrated particle-tracking engine. MP3DU remains
available through the same model-bound namespace.

```python
import simple_modflow as mf

release_points = mf.PRTReleasePoints.from_cells(
    model,
    cells=[100, 120, 140],
    layer=0,
    local_z=0.5,
)

prt = model.particle_tracking.prt(
    workspace=model.workspace.parent / "prt",
    release_points=release_points,
    porosity=0.25,
)

result = prt.run()

pathlines = result.pathlines
terminal_points = result.terminal_points

result.plot_map()
result.export_3d_html("prt_pathlines.html", vertical_exaggeration=5)
```

`PRTProject` creates a separate MF6 PRT simulation, copies the GWF DIS or DISV
grid and TDIS timing, references the GWF head and budget outputs through FMI,
and writes a PRT track CSV.

Completed PRT runs can be reopened:

```python
result = mf.open_prt_run(model, "path/to/prt_workspace")
```

The existing MP3DU workflow remains available:

```python
tracker, result = model.particle_tracking.mp3du(
    particles=[100, 120, 140],
    execute=False,
)
```

## Particle Visualization

PRT pathlines, MODPATH pathlines, and other FloPy-compatible pathline tables use
the same viewers:

```python
mf.plot_particle_pathlines(model, pathlines)

scene = mf.build_particle_tracking_scene(
    model,
    pathlines,
    vertical_exaggeration=5,
)
scene.export_html("particle_scene.html")
```

The 3D implementation follows the FloPy pattern:

```python
vtk = flopy.export.vtk.Vtk(model=model.gwf, binary=False)
vtk.add_model(model.gwf)
vtk.add_pathline_points(pathlines)
meshes = vtk.to_pyvista()
```

`simple_modflow` wraps this sequence so pathline results can be reviewed and
exported consistently.

## Testing Contract

The visualization and PRT tests cover:

- self-contained Matplotlib slider HTML without `ipywidgets`
- FloPy-style head maps, layer mosaics, and cross sections
- Plotly cross-section and map animations with all frames and sliders
- real FloPy VTK conversion
- real PyVista/Trame standalone HTML export
- DIS and DISV PRT release-point definitions
- PRT input writing and completed-run reopening
- a real GWF-to-PRT run followed by shared PyVista scene construction

# Princeton 2026 FloPy Notebook Visualization Review

This review covers the notebooks in:

`C:\Users\lukem\Python\Projects\modflow-training-princeton2026\examples\notebooks`

It records the patterns that informed the canonical `myflopy`
visualization and MF6 PRT APIs.

## Slider Patterns

### `flopy-intro-gwt-A.ipynb`

- Uses `ipywidgets.Image` and `IntSlider`.
- A Python callback renders each Matplotlib figure to an in-memory PNG.
- Demonstrates six-panel head-map mosaics through time with
  `flopy.plot.PlotMapView`.
- Demonstrates head cross sections through time with
  `flopy.plot.PlotCrossSection`.

### `flopy-intro-gwt-B.ipynb`

- Uses the same `Image` plus `IntSlider` callback pattern.
- Demonstrates concentration maps and concentration cross sections through
  time.
- Uses consistent color normalization and contour overlays.

### `parallel.ipynb`

- Uses `widgets.SelectionSlider` and `widgets.interact` for SFT concentration
  through time.
- Also demonstrates static head maps, specific-discharge vectors, and
  `PlotCrossSection`.

These patterns require a live Jupyter kernel. The corresponding standalone
`myflopy` path pre-renders the Matplotlib frames and embeds them into a
browser-side HTML slider:

- `model.visualize.head_map_slider_html(...)`
- `model.visualize.head_layer_mosaic_slider_html(...)`
- `model.visualize.cross_section_slider_html(...)`

The cross-section exporter deliberately reuses
`plot_model_cross_section(...)`, preserving the existing `myflopy`
Matplotlib cross-section style.

## Plotly Patterns Already in myflopy

Before this review, `myflopy` already supported:

- `model.plot.section(...).ani` for Plotly cross-section animations
- `model.plot.map(...).ani` for Plotly map animations

Those paths are now exposed consistently through:

- `model.visualize.plotly_cross_section_animation(...)`
- `model.visualize.plotly_head_map_animation(...)`

The Plotly map animation was corrected to include every available model frame
and an enabled slider.

## MF6 PRT Patterns

### `prt_voronoi.ipynb`

- Builds a separate PRT simulation with `ModflowPrt`.
- Copies a DISV grid into `ModflowPrtdisv`.
- Defines release points with `ModflowPrtprp`.
- References GWF heads and budgets through `ModflowPrtfmi`.
- Reads the PRT track CSV with pandas.
- Plots pathlines on a FloPy map with scalar fields sampled along paths.

### `prt_backward.ipynb`

- Compares MODPATH 7 and MF6 PRT backward tracking.
- Reverses GWF head and budget files before PRT backward tracking.
- Reads the PRT track CSV.
- Uses `PlotMapView.plot_pathline`.
- Builds 3D scenes with FloPy `Vtk`, `add_model`,
  `add_pathline_points`, and `to_pyvista`.
- Adds explicit PyVista meshes for wells, rivers, and model layers.

### `prt_watertable.ipynb`

- Builds a structured-grid PRT simulation.
- Plots pathlines colored by travel time.
- Creates 3D model, boundary, and pathline meshes with FloPy VTK and PyVista.

### `xt3d-whirls.ipynb`

- Uses FloPy VTK and PyVista for full 3D particle-path review.
- Builds reusable PyVista scene functions with custom camera perspectives.

These patterns are represented by:

- `PRTReleasePoints`
- `PRTProject`
- `PRTRunResults`
- `open_prt_run(...)`
- `model.particle_tracking.prt(...)`
- `result.plot_map(...)`
- `result.scene(...)`
- `result.export_3d_html(...)`

MP3DU remains available through `model.particle_tracking.mp3du(...)`.

## Animation Notebooks

The following notebooks use `matplotlib.animation.FuncAnimation` rather than
sliders:

- `density-bubble.ipynb`
- `density-henry-hilleke.ipynb`
- `gwe-ates.ipynb`
- `gwe-stallman.ipynb`

Those are useful references for video/GIF output, but they solve a different
interaction problem than a user-controlled standalone HTML slider.

## Canonical Direction

The preferred visualization structure is:

- Matplotlib for report-style static frames and standalone frame-slider HTML
- Plotly for vector-based interactive maps and cross-section animations
- MF6 PRT for integrated particle tracking
- MP3DU as an alternative tracking engine
- FloPy VTK plus PyVista for shared 3D particle/result scenes

The result-viewing layer should remain independent of the particle-tracking
engine whenever the engine can provide a FloPy-compatible pathline table.

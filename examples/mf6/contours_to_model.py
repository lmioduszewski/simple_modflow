"""Rasters + contours -> layer stack -> model: the copy-me starting point.

The path most real projects actually take: you have a **DEM raster** for the
ground surface, **digitized contours** for some geologic contacts, maybe another
raster for bedrock, and you need a multi-layer DISV model out of the mixture.

    raster (.tif) ------------------\\
                                     >--  Surface --LayerStack--> top/botm/idomain
    contours (.gpkg) --GRASS-------/                        |
                                            ModelContext ---+
                                                  |
                          mf.gwf(context=ctx, packages=[...]) --> Project --> run

**Source kinds mix freely.** ``Surface`` is one type whatever produced it, so
``top=`` and ``bottom=`` do not care whether a contact came from a raster, from
contours, or from algebra on another surface. Differing resolutions need no
pre-alignment either -- every surface is area-weighted onto the same Voronoi
cells.

Run it directly and it builds its own raster + contour fixtures, so it works with
no data of your own; the comments mark exactly where to swap yours in.

**GRASS prerequisites.** ``mf.Contours`` shells out to GRASS GIS, a **system**
dependency that ``pip`` does not supply: install it with your package manager
(Linux/macOS) or via OSGeo4W (Windows). Discovery is then automatic -- the
launcher is found on ``PATH`` or in the OSGeo4W/QGIS bundles, and the GRASS
Python bindings (inside the install, at ``<prefix>/etc/python``) are located by
asking that launcher, so no ``PYTHONPATH`` is needed. Override only if GRASS
lives somewhere undiscoverable::

    export GRASS_BIN=/usr/bin/grass     # or ...\\grass84.bat on Windows

If you would rather not depend on GRASS at all, ``mf.Points(xs, ys, zs)``
interpolates the same contour vertices with no external tool -- see
``surfaces_without_grass()`` at the bottom.

Reference: ``docs/model_building_cheatsheet.md``.
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import geopandas as gpd
import numpy as np
from shapely.geometry import LineString, Polygon

# Make the in-repo src/ importable when myflopy is not pip-installed.
_SRC = Path(__file__).resolve().parents[2] / "src"
if _SRC.exists() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import myflopy as mf

CRS = "EPSG:2927"
WIDTH, HEIGHT = 3000.0, 2000.0


# --------------------------------------------------------------------------- #
# Fixtures -- delete this whole section and point at your own files.
# --------------------------------------------------------------------------- #
def write_fixtures(root: Path) -> dict[str, Path]:
    """Write a boundary, a DEM raster, a contour set, and a bedrock raster."""

    gpd.GeoDataFrame(
        {"name": ["domain"]},
        geometry=[Polygon([(0, 0), (WIDTH, 0), (WIDTH, HEIGHT), (0, HEIGHT)])],
        crs=CRS,
    ).to_file(root / "boundary.gpkg", driver="GPKG")

    def raster(surface, name, res=25.0):
        """Write `surface` as a north-up GeoTIFF, standing in for a real DEM."""

        import rasterio
        from rasterio.transform import from_origin

        nx, ny = int(WIDTH / res), int(HEIGHT / res)
        data = np.array(
            [[surface(x, y) for x in np.linspace(0, WIDTH, nx)]
             for y in np.linspace(HEIGHT, 0, ny)],
            dtype="float32",
        )
        path = root / f"{name}.tif"
        with rasterio.open(path, "w", driver="GTiff", height=ny, width=nx,
                           count=1, dtype="float32", crs=CRS,
                           transform=from_origin(0, HEIGHT, res, res)) as dst:
            dst.write(data, 1)
        return path

    def contours(surface, levels, name):
        """Digitize `surface` as labelled LineStrings, the way a map would."""

        rows = []
        for z in levels:
            # Walk west to east and solve surface(x, y) = z for y. `surface`
            # takes REAL coordinates, same as the raster writer above.
            pts = []
            for x in np.linspace(0, WIDTH, 60):
                lo, hi = surface(x, 0.0), surface(x, HEIGHT)
                y = (z - lo) / max(hi - lo, 1e-9) * HEIGHT
                if 0.0 <= y <= HEIGHT:
                    pts.append((x, y))
            if len(pts) > 1:
                rows.append({"Elev": float(z), "geometry": LineString(pts)})
        if not rows:
            raise ValueError(f"{name}: no contour crossed the domain -- check levels")
        gpd.GeoDataFrame(rows, crs=CRS).to_file(root / f"{name}.gpkg", driver="GPKG")
        return root / f"{name}.gpkg"

    return {
        "boundary": root / "boundary.gpkg",
        # a DEM, as you'd get from LiDAR
        "ground_dem": raster(lambda x, y: 200.0 - 0.020 * x + 0.008 * y, "ground_dem"),
        # a contact digitized as labelled contour lines
        "clay_top": contours(lambda x, y: 150.0 - 0.015 * x + 0.005 * y,
                             np.arange(100, 166, 5), "clay_top_contours"),
        # and another raster underneath
        "bedrock": raster(lambda x, y: 60.0 - 0.008 * x + 0.003 * y, "bedrock"),
    }


# --------------------------------------------------------------------------- #
def build(root: Path):
    """Build the project, returning (project, simulation, layers, vor)."""

    src = write_fixtures(root)                      # <- your files instead

    # 1. GRID -- eagerly, because the GIS package helpers resolve cells at
    #    declaration time. `boundary=` needs a ShapeSource, not a path string.
    vor = mf.GridSpec.voronoi(
        name="valley",
        boundary=mf.ShapeSource(src["boundary"], crs=CRS),
        crs=CRS,
        boundary_max_area=40_000.0,
    ).resolve(workspace=root / "_grid")

    # 2. SOURCES -> SURFACES. Three different kinds, one `Surface` type.
    #    Rasters are read directly. `mf.Contours` interpolates through GRASS
    #    ONCE, caches to `<name>.interp.tif` beside the source, and reuses it
    #    forever after -- so only the first build pays for it. `region_vector`
    #    (or `region_raster`) is REQUIRED: it defines the GRASS region.
    ground = mf.Raster(src["ground_dem"])
    clay_top = mf.Contours(src["clay_top"], z="Elev", epsg="2927",
                           resolution=25, region_vector=src["boundary"])
    bedrock = mf.Raster(src["bedrock"])

    # 3. LAYER STACK -- top once, then each layer by bottom OR thickness.
    #    Note the sources are mixed and nothing here has to know: raster top,
    #    contoured contact, a DERIVED contact 20 ft below it, raster base.
    stack = (
        mf.LayerStack(vor, top=ground, length_units="feet")
        .add("upper_sand", bottom=clay_top, min_thickness=2.0, pinch="inactive")
        .add("clay", bottom=clay_top.below(20), pinch="passthrough")
        .add("lower_aquifer", bottom=bedrock, min_thickness=5.0, pinch="inactive")
    )
    print(stack.qc())                    # READ THIS before trusting the build
    layers = stack.build(attach=True)    # attach=True publishes vor.gdf_topbtm

    # 4. CONTEXT -- geometry the package builders read; rides on the MODEL.
    ctx = mf.ModelContext(grid=vor, domain=layers.idomain, surfaces=vor.gdf_topbtm)

    # 5. PACKAGES -- one flat declarative list.
    gp = vor.get_disv_gridprops()
    cx, _ = vor.centroids
    west = int(min(range(vor.ncpl), key=lambda i: cx[i]))
    east = int(max(range(vor.ncpl), key=lambda i: cx[i]))

    flow = mf.gwf(
        "valley",
        context=ctx,
        save_flows=True,
        packages=[
            mf.disv(nlay=layers.nlay, ncpl=gp["ncpl"], nvert=len(gp["vertices"]),
                    vertices=gp["vertices"], cell2d=gp["cell2d"],
                    top=layers.top, botm=layers.botm, idomain=layers.idomain,
                    length_units="FEET"),
            mf.ic(strt=float(np.nanmean(layers.top)) - 10.0),
            # one value per layer: sand / clay aquitard / lower aquifer
            mf.npf(k=[25.0, 0.05, 40.0], k33=[2.5, 0.005, 4.0],
                   icelltype=1, save_flows=True),
            mf.sto(steady_state={0: True}),
            mf.chd(stress_period_data={0: [
                [(0, west), float(layers.top[west]) - 5.0],
                [(0, east), float(layers.top[east]) - 15.0],
            ]}),                          # or mf.chd.gpkg("bcs.gpkg", context=ctx, nper=1)
            mf.rch.flopy(stress_period_data={
                0: [[(0, c), 4.0e-4] for c in range(vor.ncpl)]
            }),
            mf.oc(head_filerecord="valley.hds", budget_filerecord="valley.cbc",
                  saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
        ],
    )

    # 6. SIMULATION + PROJECT
    simulation = mf.SimulationSpec(
        "baseline",
        models=(flow,),
        packages=(
            mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
            mf.ims(models=("valley",), complexity="COMPLEX"),
        ),
    )
    project = mf.Project(root / "project", name="valley")
    project.add_simulation(simulation)
    return project, simulation, layers, vor


def surfaces_without_grass(contours_gpkg: Path) -> mf.Surface:
    """The GRASS-free fallback: interpolate the contour VERTICES as points.

    Coarser than ``r.surf.contour`` -- it knows nothing about the contours being
    nested lines -- but it needs no external tool and is often good enough for a
    first pass or for a surface with dense control.
    """

    gdf = gpd.read_file(contours_gpkg)
    xs, ys, zs = [], [], []
    for z, geom in zip(gdf["Elev"], gdf.geometry, strict=True):
        for x, y in geom.coords:
            xs.append(x)
            ys.append(y)
            zs.append(float(z))
    return mf.Points(xs, ys, zs, method="linear")


def main() -> int:
    root = Path(tempfile.mkdtemp(prefix="contours_to_model_"))
    project, _, layers, _ = build(root)

    run = project.prepare_run("baseline", "baseline")
    success, report = run.execute()
    print(f"\nworkspace: {run.workspace}")
    print(f"converged: {success}")
    if not success:
        print("\n".join(report[-25:]))
        return 1

    model = run.model("valley")
    heads = model.hds.array(layer=0)
    print(f"layer-0 heads: {heads.min():.1f} .. {heads.max():.1f}")

    # Everything drawn is a Picture: .fig / .show() / .save(p) / .html(p).
    # `.html()` is used here because it needs nothing extra; `.save("x.png")`
    # rasterizes through kaleido, which requires Chrome (`plotly_get_chrome`).
    model.plot.map(layer=0, contours=True).html(root / "heads.html")
    layers.plot.section(y=HEIGHT / 2).save(root / "section.png")   # Matplotlib
    print(f"wrote {root / 'heads.html'} and {root / 'section.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Rasters + contours -> layer stack -> model: the copy-me starting point.

The path most real projects actually take: you have a **DEM raster** for the
ground surface, **digitized contours** for some geologic contacts, maybe another
raster for bedrock, and you need a multi-layer DISV model out of the mixture.

    Project                       <- made FIRST; paths resolve against its root
      inputs/*.tif, *.gpkg        <- your reference data
        raster ---------\\
                          >-- Surface -- LayerSurfaces      (no grid needed yet)
        contours --GRASS/                     |
      add_grid("base", GridSpec) -> resolve --+
                                              v
                          to_disv(vor) -> add_package("disv/base")
                                              |
                    mf.gwf(packages=[mf.ref("disv/base"), ...]) -> run

**The Project comes first.** It is the root that relative data paths resolve
against and the registry that ``mf.ref``/``mf.grid_ref`` resolve through, so
declaring against it from the start makes the whole model a portable recipe that
``save()``/``load()`` round-trips.

**The layering does not need the grid.** Surfaces are lazy and ``LayerSurfaces``
is just an ordered list of them, so it is declared before the grid exists.
``Surface`` is one type whatever produced it -- raster, contours, or algebra on
another surface -- and differing resolutions need no pre-alignment, since every
surface is area-weighted onto the same Voronoi cells.

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
    """Write a boundary, a DEM raster, a contour set, and a bedrock raster.

    ``root`` is the project's ``inputs/`` directory, so everything lands where
    a relative ``ShapeSource("inputs/...")`` will find it.
    """

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
    """Build the project, returning (project, simulation, layering, vor)."""

    # ------------------------------------------------------------------ #
    # 1. PROJECT FIRST. It is the root that relative data paths resolve
    #    against and the registry that mf.ref/mf.grid_ref resolve through,
    #    so everything below is declared against it.
    # ------------------------------------------------------------------ #
    project = mf.Project(root / "valley", name="valley")
    project.layout.ensure()                      # specs/ inputs/ runs/

    src = write_fixtures(project.layout.inputs_dir)   # <- your files instead

    # ------------------------------------------------------------------ #
    # 2. SOURCES -> SURFACES. Three kinds, one `Surface` type, NO grid yet:
    #    surfaces are lazy, so nothing is read or interpolated here.
    #    `mf.Contours` runs GRASS ONCE, caches to `<name>.interp.tif` beside
    #    the source, and reuses it forever after. `region_vector` (or
    #    `region_raster`) is REQUIRED: it defines the GRASS region.
    # ------------------------------------------------------------------ #
    ground = mf.Raster(src["ground_dem"])
    clay_top = mf.Contours(src["clay_top"], z="Elev", epsg="2927",
                           resolution=25, region_vector=src["boundary"])
    bedrock = mf.Raster(src["bedrock"])

    # ------------------------------------------------------------------ #
    # 3. LAYERING -- still no grid. `LayerSurfaces` is an ordered list of
    #    surfaces: [0] is the model top, each one after is a layer bottom.
    #    Sources are mixed and nothing here has to know: raster top,
    #    contoured contact, a DERIVED contact 20 ft below it, raster base.
    # ------------------------------------------------------------------ #
    layering = mf.LayerSurfaces(
        [ground, clay_top, clay_top.below(20), bedrock],
        labels=["ground", "upper_sand", "clay", "lower_aquifer"],
    )

    # ------------------------------------------------------------------ #
    # 4. GRID -- the recipe is registered on the project, with a path
    #    RELATIVE to the project root, then resolved into a real grid.
    # ------------------------------------------------------------------ #
    project.add_grid("base", mf.GridSpec.voronoi(
        name="valley",
        boundary=mf.ShapeSource("inputs/boundary.gpkg", crs=CRS),
        crs=CRS,
        boundary_max_area=40_000.0,
    ))
    vor = project.grids["base"].resolve(
        project_root=project.root,                     # <- anchors the path
        workspace=project.root / "grids" / "base",
    )

    # ------------------------------------------------------------------ #
    # 5. BIND the layering to the grid -> a PackageSpec -> register it.
    #    Putting DISV in the library rather than inline is what lets the
    #    project be saved: the simulation then holds only a reference, and
    #    the arrays ride in a pickle sidecar beside the library JSON.
    # ------------------------------------------------------------------ #
    print(layering.thickness_report(vor, minimum_thickness=2.0))
    project.add_package("disv/base", layering.to_disv(
        vor, pinch_out=True, minimum_thickness=2.0,
        length_units="FEET", attach=True))          # attach -> vor.gdf_topbtm
    project.add_package("npf/base", mf.npf(
        k=[25.0, 0.05, 40.0], k33=[2.5, 0.005, 4.0],
        icelltype=1, save_flows=True))

    # ------------------------------------------------------------------ #
    # 6. CONTEXT -- geometry the package builders read; rides on the MODEL.
    # ------------------------------------------------------------------ #
    ctx = mf.ModelContext(grid=vor, surfaces=vor.gdf_topbtm)

    # ------------------------------------------------------------------ #
    # 7. PACKAGES -- one flat declarative list, library entries by name.
    # ------------------------------------------------------------------ #
    top = np.asarray(vor.gdf_topbtm["ground"], dtype=float)
    cx, _ = vor.centroids
    west = int(min(range(vor.ncpl), key=lambda i: cx[i]))
    east = int(max(range(vor.ncpl), key=lambda i: cx[i]))

    flow = mf.gwf(
        "valley",
        context=ctx,
        save_flows=True,
        packages=[
            mf.ref("disv/base"),          # <- the layer stack, by reference
            mf.ic(strt=float(np.nanmean(top)) - 10.0),
            mf.ref("npf/base"),
            mf.sto(steady_state={0: True}),
            mf.chd(stress_period_data={0: [
                [(0, west), float(top[west]) - 5.0],
                [(0, east), float(top[east]) - 25.0],
            ]}),                          # or mf.chd.gpkg("bcs.gpkg", context=ctx, nper=1)
            mf.rch.flopy(stress_period_data={
                0: [[(0, c), 4.0e-4] for c in range(vor.ncpl)]
            }),
            mf.oc(head_filerecord="valley.hds", budget_filerecord="valley.cbc",
                  saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
        ],
    )

    # ------------------------------------------------------------------ #
    # 8. SIMULATION, registered on the project.
    # ------------------------------------------------------------------ #
    simulation = mf.SimulationSpec(
        "baseline",
        models=(flow,),
        packages=(
            mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
            mf.ims(models=("valley",), complexity="COMPLEX"),
        ),
    )
    project.add_simulation(simulation)
    return project, simulation, layering, vor


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
    project, _, layering, vor = build(root)

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
    # `LayerSurfaces` has no `.plot` namespace; the `LayerStack` facade does.
    # `from_modflow` takes a flopy model (anything with `.modelgrid`), so read
    # the geometry straight back off the built model to draw it.
    stack = mf.LayerStack.from_modflow(vor, model.gwf, resample=False)
    stack.build().plot.section(y=HEIGHT / 2).save(root / "section.png")
    print(f"wrote {root / 'heads.html'} and {root / 'section.png'}")

    # The payoff of declaring against the project: it is a portable recipe.
    # The simulation holds only mf.ref/mf.grid_ref, so it stays JSON-clean;
    # the array-laden DISV rides in a pickle sidecar beside its library JSON.
    print(f"validate(): {project.validate() or 'serializable'}")
    project.save()
    reloaded = mf.Project.load(project.root)
    rerun = reloaded.prepare_run("reloaded", "baseline")
    print(f"reloaded from disk and rebuilt: nlay={rerun.model('valley').nlay}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

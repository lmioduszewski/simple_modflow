"""Pattern 4 -- layer surfaces -> a real DISV grid -> a model that solves.

Build a Voronoi mesh once, describe the layering as a stack of `Surface`s, and
hand it to `mf.disv`. Each surface is independent and swappable -- flat, a
raster, GRASS-interpolated contours, or scattered points -- so the *layering*
becomes a reusable, swappable ingredient of a variant, exactly like packages
and grids.

    layers = mf.LayerSurfaces([top, layer1_botm, layer2_botm])
    disv   = layers.to_disv(vor)          # samples + reconciles -> DISV spec

The surfaces here are flat / point-interpolated so the example needs no data
files. In a real model any surface drops in unchanged:

    mf.Surface.raster("ground.tif")
    mf.Surface.from_contours("aq_botm.gpkg", z="elev", region_raster="dom.tif")

Run:
    python 04_layer_surfaces_to_disv.py
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import myflopy as mf

HERE = Path(__file__).resolve().parent


def voronoi_grid() -> mf.VoronoiGridPlus:
    """Build one Voronoi mesh over a 1000 x 600 rectangular domain."""

    tri = mf.TriangleGrid(model_ws=tempfile.mkdtemp(), angle=30)
    tri.set_domain_rectangle(x_dist=1000, y_dist=600, origin=(0, 0), max_area=12000)
    tri.build()
    return mf.VoronoiGridPlus(tri)


def flow_model(name: str, vor: mf.VoronoiGridPlus, layers: mf.LayerSurfaces):
    """A steady GWF model whose DISV comes straight from `layers.to_disv(vor)`."""

    cx, _ = vor.centroids
    west = int(min(range(vor.ncpl), key=lambda i: cx[i]))   # inflow cell
    east = int(max(range(vor.ncpl), key=lambda i: cx[i]))   # outflow cell

    flow = mf.gwf(
        name,
        grid=vor,
        packages=[
            layers.to_disv(vor),                            # <-- the capstone call
            mf.ic(strt=45.0),
            mf.npf(k=10.0, icelltype=1),
            mf.chd(stress_period_data={0: [[(0, west), 50.0], [(0, east), 40.0]]}),
            mf.oc(head_filerecord=f"{name}.hds", saverecord=[("HEAD", "ALL")]),
        ],
    )
    return mf.SimulationSpec(
        name,
        models=[flow],
        packages=[
            mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
            mf.ims(models=[name], complexity="SIMPLE"),
        ],
    )


vor = voronoi_grid()

# Two layerings on the SAME mesh -- swap one surface to make a variant.
flat = mf.LayerSurfaces(
    [mf.Surface.flat(50.0), mf.Surface.flat(20.0), mf.Surface.flat(0.0)],
    labels=["top", "botm1", "botm2"],
)
sloped = mf.LayerSurfaces(
    [
        mf.Surface.flat(50.0),
        # middle surface dips west->east; reconcile keeps it between neighbours
        mf.Surface.from_points([0, 1000, 0, 1000], [0, 0, 600, 600], [30, 10, 30, 10]),
        mf.Surface.flat(0.0),
    ],
    labels=["top", "botm1", "botm2"],
)

project = mf.Project(HERE / "runs" / "surfaces_demo", name="surfaces_demo")
flat_run = project.run("flat_layers", flow_model("flat_layers", vor, flat))
sloped_run = project.run("sloped_layers", flow_model("sloped_layers", vor, sloped))


def layer1_botm(run):
    botm = run.built.built_model("flat_layers" if run is flat_run else "sloped_layers")
    return botm.package("disv").botm.array[0]


for label, run in [("flat", flat_run), ("sloped", sloped_run)]:
    b1 = layer1_botm(run)
    print(
        f"{label:7s} success={run.success}  ncpl={vor.ncpl}  "
        f"layer1 botm: {b1.min():.1f}..{b1.max():.1f}"
    )

print("-> same mesh, swapped layer surfaces, both models solve")

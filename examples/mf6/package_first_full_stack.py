"""Package-first full-stack reference: a multi-layer valley with the works.

This is the copy-me template for building a real MODFLOW 6 model with myflopy's
**package-first API** -- the preferred way to assemble models. One model carries
the full surface-water + unsaturated stack: **GHB, RCH, DRN, UZF, SFR, LAK and
MVR**, on a multi-layer DISV (Voronoi) grid, inside a :class:`myflopy.Project`.

How it fits together (read top to bottom)::

    Project            durable workspace + run lifecycle (holds NO geometry)
      SimulationSpec   one MF6 simulation: tdis + ims solver + the model(s)
        ModelSpec      = mf.gwf(name, context=ctx, packages=[...])
          ModelContext grid + domain + surfaces  (rides on the MODEL)
          packages     disv, npf/ic/sto/oc, ghb/drn/rch, uzf/sfr/lak, mvr

The grid is built **eagerly** and wrapped in a ``ModelContext`` because the GIS
package helpers (``mf.uzf``/``mf.sfr``/``mf.lak``/``mf.X.gpkg``) resolve cells the
moment you call them. Layers are authored with the ``mf.LayerStack`` facade. The
boundary geometry here is borrowed from the canonical valley so the model is
self-contained and converges; the comments mark exactly where you would swap in
your own geopackages, rasters and Excel/CSV tables.

Run it directly to build + execute, or import :func:`build_full_stack_project`.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

# Make the in-repo src/ importable when myflopy is not pip-installed.
_SRC = Path(__file__).resolve().parents[2] / "src"
if _SRC.exists() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

import myflopy as mf
from myflopy.layers import Array
from myflopy.modflow.mf6.canonical_example import (
    CanonicalModelConfig,
    _write_surface_water_inputs,
    rectangular_voronoi,
)


def build_full_stack_project(root: Path | str, *, ncol: int = 22, nrow: int = 16):
    """Assemble the full package-first stack and return a built :class:`mf.Project`.

    Parameters
    ----------
    root
        Project root directory (runs land under ``<root>/runs/<name>``).
    ncol, nrow
        Grid size of the demonstration valley. Small by default so the example
        stays fast; raise for a finer grid.

    Returns
    -------
    (project, simulation)
        The :class:`myflopy.Project` (call ``project.prepare_run(...)``) and the
        :class:`myflopy.SimulationSpec` registered on it.
    """

    root = Path(root)
    config = CanonicalModelConfig(nrow=nrow, ncol=ncol, nlay=2, nper=1)
    nper = config.nper

    # ------------------------------------------------------------------ #
    # 1. PROJECT -- the durable workspace + run/scenario lifecycle.
    # ------------------------------------------------------------------ #
    project = mf.Project(root, name="full_stack_demo")

    # ------------------------------------------------------------------ #
    # 2. GRID -- built EAGERLY (the GIS package helpers resolve cells now).
    #    Swap this for `mf.GridSpec.voronoi(boundary=<gpkg>, refinement=<gpkg>)
    #    .resolve(root / "_grid")` to drive the grid from your geopackages.
    # ------------------------------------------------------------------ #
    vor = rectangular_voronoi(config)
    ncpl = int(vor.ncpl)
    centers = np.asarray(vor.points, dtype=float)
    xn = centers[:, 0] / (config.ncol * config.cell_size)   # 0 up-valley .. 1 mouth
    yn = centers[:, 1] / (config.nrow * config.cell_size)   # 0 .. 1 across valley

    # ------------------------------------------------------------------ #
    # 3. LAYERS -- authored with the LayerStack FACADE (compiles to the
    #    LayerSurfaces engine). Here the top is an array; your real model would
    #    use mf.Raster("ground.tif") / mf.Contours("base.shp") etc.
    # ------------------------------------------------------------------ #
    ground = 150.0 - 52.0 * xn + 8.0 * (2.0 * np.abs(yn - 0.5))   # down-valley slope + gentle U
    stack = (
        mf.LayerStack(vor, top=Array(ground), length_units="feet")
        .add("alluvium", thickness=40.0)     # upper unconfined aquifer
        .add("basin_fill", thickness=60.0)   # lower aquifer
    )
    # build(attach=True) also publishes the elevations onto vor.gdf_topbtm so the
    # grid-aware builders (SFR reach tops, LAK lake-cell layering) can read them.
    layers = stack.build(attach=True)   # -> disv-ready top / botm / idomain

    # ------------------------------------------------------------------ #
    # 4. CONTEXT -- geometry the package builders read. Rides on the MODEL.
    # ------------------------------------------------------------------ #
    ctx = mf.ModelContext(grid=vor, domain=layers.idomain, surfaces=vor.gdf_topbtm)

    gp = vor.get_disv_gridprops()
    regional = 146.0 - 60.0 * xn                       # regional water table
    strt = np.vstack([np.maximum(layers.botm[k] + 1.0, regional) for k in range(config.nlay)])
    k_alluvium = 12.0 + 18.0 * xn + 30.0 * np.exp(-((yn - 0.5) / 0.15) ** 2)  # paleochannel ribbon
    k = np.vstack([k_alluvium, k_alluvium * 0.4])

    # ------------------------------------------------------------------ #
    # 5. BOUNDARY DATA -- here computed/borrowed; swap in your gpkg/Excel.
    #    The lake + stream geometry comes from the canonical valley so the demo
    #    is self-contained: in your model these are mf.lak(lakes="lakes.gpkg")
    #    and mf.sfr(streams="streams.gpkg").
    # ------------------------------------------------------------------ #
    sw = _write_surface_water_inputs(root / "_inputs", config)

    # GHB: regional underflow leaving at the valley mouth (east edge).
    mouth = sorted(range(ncpl), key=lambda c: xn[c])[-max(3, ncpl // 40):]
    ghb_data = {0: [[(config.nlay - 1, c), float(86.0), 50.0] for c in mouth]}

    # DRN: toe-of-slope springs along the valley floor mid-section.
    drn_cells = [c for c in range(ncpl) if 0.30 < xn[c] < 0.45 and abs(yn[c] - 0.5) < 0.12]
    drn_data = {0: [[(0, c), float(regional[c] - 1.0), 30.0] for c in drn_cells]}

    # RCH: mountain-front recharge on the two valley walls.
    rch_rate = np.where(np.abs(yn - 0.5) > 0.30, 6.0e-4, 1.5e-4)
    rch_data = {0: [[(0, c), float(rch_rate[c])] for c in range(ncpl)]}

    # UZF: unsaturated zone + ET across the dry valley floor (away from water).
    floor = [c for c in range(ncpl) if np.abs(yn[c] - 0.5) < 0.30 and xn[c] > 0.15]
    uzf_cells = [(0, c) for c in floor]
    finf = {0: [3.0e-5] * len(uzf_cells)}
    pet = {0: [1.0e-4] * len(uzf_cells)}

    # ------------------------------------------------------------------ #
    # 6. MODEL -- one flat, declarative package list. THE readable payoff.
    # ------------------------------------------------------------------ #
    # The stream + lake are pulled out as handles so the mover can reference them
    # SEMANTICALLY (mf.sfr_connection / mf.lak_connection) instead of raw indices.
    sfr = mf.sfr(context=ctx, nper=nper, streams=sw["streams"],
                 connection_mode="automatic",
                 connections=(mf.StreamConnection("north_trib", "main_stem"),
                              mf.StreamConnection("south_trib", "main_stem")),
                 width=18.0, gradient=0.0012,
                 roughness=0.030, streambed_k=0.05, streambed_thickness=1.5,
                 inflow={0: {"north_trib": 15000.0, "south_trib": 10000.0}},
                 length_conversion=3.28081, time_conversion=86_400.0, mover=True)
    lak = mf.lak(context=ctx, nper=nper, lakes=sw["lakes"], lake_id_field="name",
                 starting_stage={"valley_lake": 101.0}, lake_bottom={"valley_lake": 96.0},
                 bed_leakance=0.11, connection_modes="automatic",
                 status={"valley_lake": ["ACTIVE"]}, mover=True,
                 length_conversion=3.28081, time_conversion=86_400.0)

    flow = mf.gwf(
        "valley",
        context=ctx,
        newtonoptions="UNDER_RELAXATION",
        save_flows=True,
        packages=[
            mf.disv(
                nlay=config.nlay, ncpl=gp["ncpl"], nvert=len(gp["vertices"]),
                vertices=gp["vertices"], cell2d=gp["cell2d"],
                top=layers.top, botm=layers.botm, idomain=layers.idomain,
                length_units="FEET",
            ),
            mf.ic(strt=strt),
            mf.npf(k=k, k33=k * 0.1, save_flows=True),
            mf.sto(steady_state={0: True}),
            mf.ghb(stress_period_data=ghb_data),       # or mf.ghb.gpkg("bcs.gpkg", layer="underflow", context=ctx, nper=nper)
            mf.drn(stress_period_data=drn_data),        # or mf.drn.gpkg(...)
            mf.rch.flopy(stress_period_data=rch_data),  # or mf.rch(context=ctx, nper=nper, recharge=<array>)
            mf.uzf(context=ctx, nper=nper, cells=uzf_cells,
                   vks=0.25, thtr=0.08, thts=0.34, thti=0.17, eps=4.0,
                   finf=finf, pet=pet, extdp=7.0),
            sfr,
            lak,
            # MVR: move part of the main stem's OUTFLOW (its final reach, resolved
            # from geometry) into the lake -- no hard-coded reach numbers. Pass
            # at=(x, y) instead to target the reach nearest a coordinate.
            mf.mvr(nper=nper, moves=(mf.Move(mf.sfr_connection(sfr, "main_stem"),
                                             mf.lak_connection(lak, "valley_lake"),
                                             value=0.5),)),
            mf.oc(head_filerecord="valley.hds", budget_filerecord="valley.cbc",
                  saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
        ],
    )

    # ------------------------------------------------------------------ #
    # 7. SIMULATION + register on the project.
    # ------------------------------------------------------------------ #
    simulation = mf.SimulationSpec(
        "baseline",
        models=(flow,),
        packages=(
            mf.tdis(nper=nper, perioddata=[(1.0, 1, 1.0)] * nper),
            mf.ims(models=("valley",), complexity="COMPLEX",
                   outer_maximum=250, inner_maximum=250,
                   linear_acceleration="BICGSTAB"),
        ),
    )
    project.add_simulation(simulation)
    return project, simulation


def main() -> int:
    import tempfile

    root = Path(tempfile.mkdtemp(prefix="full_stack_demo_"))
    project, _ = build_full_stack_project(root)
    run = project.prepare_run("baseline", "baseline")   # builds in memory
    success, report = run.execute()                     # writes + runs MF6
    print(f"workspace: {run.workspace}")
    print(f"converged: {success}")
    if not success:
        print("\n".join(report[-25:]))
    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main())

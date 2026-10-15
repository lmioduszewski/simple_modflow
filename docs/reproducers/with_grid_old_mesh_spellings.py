"""Which mesh `with_grid` measures a per-cell array against, by how the grid was declared.

Ledger 599 (the old-mesh rule itself landed as ledger 597's `_declared_mesh`). `with_grid` took the OLD mesh -- the cell count a per-cell array is measured
against, and the pin it stamps -- from ``context.grid`` alone, so the spelling
``docs/package_api_reference.md`` teaches, ``mf.gwf("m", grid=<built grid>)``, read an
``(nlay, ncpl)`` ``ic.strt`` ``unchanged`` and BUILT on the old mesh's values. This
script declares one model three ways on a 4 x 7 ``rectangular_voronoi`` mesh (ncpl 28,
two layers, a lazy ``mf.disv``, ``strt`` an ``(2, 28)`` array), moves it to a 7 x 4 mesh
of the SAME cell count -- the case a range check cannot see -- and reports, per spelling,
``move_report`` before and after, whether the build refuses, and ``Run.resolution`` on
the unmoved spec. Then three conditions the entry names:

* an ``mf.ref`` library package holding the same ``strt`` (the saveable form "Starting
  heads from a solved run" teaches), on ``grid=<built>`` and on ``grid=mf.grid_ref``;
* ``GridSpec.resolve()`` with no workspace, run in an empty working directory;
* the same recipe meshed three times (a 3,000 x 2,000 ft box with a buffered diagonal
  refinement, at two sizes): how many distinct grid pins -- the remap meshes a recipe
  once for its overlay and the build meshes it again.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/with_grid_old_mesh_spellings.py

It uses only the public API (plus ``grid_pin``), so it measures whichever tree is on
``PYTHONPATH``; "before" is the parent commit's ``src`` (``git archive <parent> src | tar
-x -C <dir>``, then ``PYTHONPATH=<dir>/src``). Needs no MODFLOW binary.

MEASURED (2026-10-15)                          before (cb8073e5)        now (ledgers 597 + 599)
  grid=<built grid>       ic before / after    unchanged / unchanged    unverifiable / unverifiable
                          build                BUILT                    refused
                          run.resolution ic    unchanged                unverifiable
  grid=GridSpec.from_object                    (same as grid=<built grid>, both columns)
  context=ModelContext()  ic before / after    unverifiable / unverif.  unverifiable / unverifiable
                          build                refused                  refused
                          run.resolution ic    unchanged                unverifiable
  mf.ref library ic       grid=<built>         unchanged, BUILT (2,28)  unverifiable, refused
                          grid=grid_ref        unchanged, BUILT (2,28)  unverifiable, refused
  resolve() no workspace  in an empty cwd      meshed into ./_triangle  ValueError, cwd empty
  distinct pins, 3 meshes ncpl 1178 / 3054     1 / 1                    1 / 1

The ``mf.ref`` rows are classified through ``library=project`` (ledger 597): before it the
spec saw a reference with no options. The pins are one machine and one Triangle build;
nothing checks the two meshes of a remap agree.
"""

from __future__ import annotations

import os
import pathlib
import tempfile
import warnings

import geopandas as gpd
import numpy as np
from shapely.geometry import LineString, box

import myflopy as mf
from myflopy.grid_identity import grid_pin
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, rectangular_voronoi

warnings.simplefilter("ignore")


def _stack():
    stack = mf.LayerStack(top=mf.Surface.flat(100.0))
    stack.add("upper", thickness=20.0)
    stack.add("lower", thickness=30.0)
    return stack


def _sim(model):
    return mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])


def _builds(sim, project, run):
    try:
        built = project.prepare_run(run, sim)
    except ValueError:
        return "refused"
    strt = built.built.simulation.get_model("m").get_package("ic").strt.array
    return f"BUILT {strt.shape}"


def spellings(root):
    a = rectangular_voronoi(CanonicalModelConfig(nrow=4, ncol=7, nlay=2, nper=1))
    b = rectangular_voronoi(CanonicalModelConfig(nrow=7, ncol=4, nlay=2, nper=1))
    assert a.ncpl == b.ncpl == 28
    strt = np.full((2, a.ncpl), 90.0)
    declare = {
        "grid=<built grid>": lambda pk: mf.gwf("m", grid=a, packages=pk),
        "grid=GridSpec.from_object": lambda pk: mf.gwf(
            "m", grid=mf.GridSpec.from_object(a), packages=pk
        ),
        "context=ModelContext(grid=)": lambda pk: mf.gwf(
            "m", context=mf.ModelContext(grid=a), packages=pk
        ),
    }
    project = mf.Project(root / "p", name="p")
    for number, (label, make) in enumerate(declare.items()):
        sim = _sim(make([mf.disv(_stack()), mf.npf(k=1.0), mf.ic(strt=strt),
                         mf.oc(head_filerecord="m.hds")]))
        before = {r["package"]: r["how"] for r in sim.move_report()}["ic"]
        moved = sim.with_grid(b)
        after = {r["package"]: r["how"] for r in moved.move_report()}["ic"]
        resolution = {
            r["package"]: r["how"]
            for r in project.prepare_run(f"plain{number}", sim).resolution
        }["ic"]
        print(f"{label:28s} ic {before} / {after}; build {_builds(moved, project, f'moved{number}')}; "
              f"run.resolution ic {resolution}")

    project.add_package("ic_regrid", mf.ic(strt=strt))
    project.add_grid("mesh", mf.GridSpec.from_object(a))
    for number, (label, grid) in enumerate((("grid=<built>", a), ("grid=grid_ref", mf.grid_ref("mesh")))):
        sim = _sim(mf.gwf("m", grid=grid, packages=[mf.disv(_stack()), mf.npf(k=1.0),
                                                    mf.ref("ic_regrid"), mf.oc(head_filerecord="m.hds")]))
        moved = sim.with_grid(b, library=project)          # ledger 597: refs need it
        how = {r["package"]: r["how"] for r in moved.move_report(library=project)}["ic_regrid"]
        built = _builds(moved, project, f"ref{number}")
        print(f"mf.ref library ic, {label:14s} {how}, {built}")


def no_workspace(root):
    gpd.GeoDataFrame({"n": ["d"]}, geometry=[box(0, 0, 600, 600)], crs=2927).to_file(root / "b.gpkg")
    spec = mf.GridSpec.voronoi(boundary=mf.GeoPackageSourceSpec(path=root / "b.gpkg"), max_area=40_000.0)
    cwd = root / "cwd"
    cwd.mkdir()
    start = os.getcwd()
    os.chdir(cwd)
    try:
        try:
            spec.resolve()
            outcome = "meshed"
        except ValueError:
            outcome = "ValueError"
        print(f"resolve() with no workspace: {outcome}; cwd holds {sorted(p.name for p in cwd.iterdir())}")
    finally:
        os.chdir(start)


def determinism(root):
    gpd.GeoDataFrame({"n": ["d"]}, geometry=[box(0, 0, 3000, 2000)], crs=2927).to_file(root / "box.gpkg")
    gpd.GeoDataFrame(
        {"n": ["r"]}, geometry=[LineString([(100, 100), (2900, 1800)]).buffer(60)], crs=2927
    ).to_file(root / "river.gpkg")
    for max_area in (40_000.0, 4_000.0):
        spec = mf.GridSpec.voronoi(
            boundary=mf.GeoPackageSourceSpec(path=root / "box.gpkg"),
            refinement=mf.GeoPackageSourceSpec(path=root / "river.gpkg"),
            max_area=max_area, refinement_max_area=max_area / 10,
        )
        grids = [spec.resolve(workspace=root / f"m{max_area:.0f}_{i}") for i in range(3)]
        print(f"max_area={max_area:.0f}: ncpl {grids[0].ncpl}, distinct pins over 3 meshes: "
              f"{len({grid_pin(grid) for grid in grids})}")


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as scratch:
        root = pathlib.Path(scratch)
        spellings(root)
        no_workspace(root)
        determinism(root)

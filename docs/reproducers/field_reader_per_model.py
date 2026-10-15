"""What ``model.hds`` costs, and what a re-run into a reused directory does to the head file.

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python \
        docs/reproducers/field_reader_per_model.py [run_workspace]

Needs ``mf6`` on PATH (measured with 6.8.0.dev0). Part 1 runs MF6 four times on a 16-cell
toy in a temporary directory (removed afterwards), about 2 s. Part 2 runs only when a
COMPLETED run's workspace is given -- the figures in ledger 598 are Ten Trails',
``/home/lukem/models/mf6/tentrails/runs/sw_mounding_flow_pitk_9232026_QpogK15_exc_low_Qvr``
(20,323 cells per layer, 38 saved records) -- and READS it through
``mf.load_mf6_run``: it writes nothing there.

Run it against the code before ledger 598 (``PYTHONPATH=<checkout of 430154dd>/src``) for
the "before" column. Measured 2026-10-15, part 2 twice per tree, at load averages of ~12
and ~17-18 (other sessions' work) -- the clocks move with the load:

                                              before 598          after 598
  1. re-run: same inode, new size and mtime      yes                 yes
     readers ``model.hds.map()`` builds           5                   0 (after the first
                                                                      ``model.hds``)
     ResourceWarnings, 5 accesses                 10 (2 per reader)   0
     a reader held across the re-run              reads the new       raises RuntimeError
                                                  bytes, old labels
  2. first ``model.hds``                          0.76-0.99 s         0.006-0.016 s
     each later ``model.hds``                     0.35-1.02 s         0.0001-0.0003 s
     ``model.hds.map()``                          2.73-4.77 s         0.35-2.13 s
     ``nper`` / ``numstp``                        37 / 64             38 / 65

The "before" cost per access IS the reader's five per-cell lists: built cell by cell they
take 0.76-0.82 s on their own, against 8-14 ms for FloPy's index of the file (the
second run's figures).
"""

from __future__ import annotations

import gc
import os
import sys
import tempfile
import time
import warnings
from pathlib import Path

import flopy
import numpy as np

import myflopy as mf
from myflopy.modflow.mf6 import headsplus
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid


def _simulation(vor, nper, west):
    props = vor.get_disv_gridprops()
    edges = vor.gdf_vorPolys.geometry.bounds
    cells = np.flatnonzero(edges.minx <= edges.minx.min() + 1e-6)
    model = mf.gwf("m", grid=vor, packages=[
        mf.disv(nlay=1, ncpl=props["ncpl"], nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=100.0, botm=0.0),
        mf.npf(k=10.0), mf.ic(strt=50.0),
        mf.chd(stress_period_data={0: [[(0, int(c)), west] for c in cells]}),
        mf.oc(head_filerecord="m.hds", saverecord=[("HEAD", "ALL")]),
    ])
    return mf.simulation(model, name="s", packages=[
        mf.tdis(nper=nper, perioddata=[(1.0, 1, 1.0)] * nper), mf.ims(models=["m"])])


def _count_readers(call):
    built = []
    original = headsplus.DependentVariableFile.__init__

    def counting(self, *args, **kwargs):
        built.append(self)
        original(self, *args, **kwargs)

    headsplus.DependentVariableFile.__init__ = counting
    try:
        call()
    finally:
        headsplus.DependentVariableFile.__init__ = original
    return len(built)


def part1(root: Path) -> None:
    vor = irregular_voronoi_grid(nrow=4, ncol=4, cell_size=100.0)
    project = mf.Project(root / "p", name="p")
    run = project.prepare_run("r", _simulation(vor, 2, 50.0), overwrite=True)
    assert run.execute()[0]
    model = run.model()
    path = Path(run.workspace) / "m.hds"

    model.hds          # the first reader, whatever it costs
    print(f"1. readers one model.hds.map() builds: {_count_readers(lambda: model.hds.map())}")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ResourceWarning)
        for _ in range(5):
            model.hds
        gc.collect()
    print(f"   ResourceWarnings over 5 accesses: "
          f"{sum(issubclass(w.category, ResourceWarning) for w in caught)}")

    held = model.hds
    before = os.stat(path)
    rerun = project.prepare_run("r", _simulation(vor, 3, 60.0), overwrite=True)
    assert rerun.execute()[0]
    after = os.stat(path)
    print(f"   re-run: same inode {before.st_ino == after.st_ino}, size {before.st_size} -> "
          f"{after.st_size}, mtime changed {before.st_mtime_ns != after.st_mtime_ns}")
    try:
        values = held.get_data()
        print(f"   held reader: {len(held.kstpkper)} records, reads {values.ravel()[0]:.1f} "
              "(the re-run's west head is 60.0, the original's 50.0)")
    except RuntimeError as error:
        print(f"   held reader: RuntimeError: {str(error)[:70]}...")
    print(f"   model.hds afterwards: {len(model.hds.kstpkper)} records; model.kstpkper "
          f"{model.kstpkper}")


def part2(workspace: Path) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = mf.load_mf6_run(workspace, verbosity_level=0)
        model.vor.gdf_vorPolys
        start = time.perf_counter()
        reader = model.hds
        first = time.perf_counter() - start
        later = []
        for _ in range(3):
            start = time.perf_counter()
            model.hds
            later.append(time.perf_counter() - start)
        maps = []
        for _ in range(2):
            start = time.perf_counter()
            model.hds.map()
            maps.append(time.perf_counter() - start)
    polygons = model.vor.gdf_vorPolys.geometry.to_list()
    start = time.perf_counter()
    _lists = ([cell.area for cell in polygons], [cell.centroid.xy[0][0] for cell in polygons],
              [cell.centroid.xy[1][0] for cell in polygons])
    lists = time.perf_counter() - start
    start = time.perf_counter()
    with flopy.utils.HeadFile(str(reader.output_path)):
        index = time.perf_counter() - start
    print(f"2. {workspace.name}: {model.vor.ncpl} cells per layer, {len(reader.kstpkper)} records")
    print(f"   the reader's per-cell lists, built cell by cell: {lists:.3f} s; FloPy's index "
          f"of the file: {index:.4f} s")
    print(f"   first model.hds {first:.3f} s; later {', '.join(f'{t:.4f}' for t in later)} s")
    print(f"   model.hds.map() {', '.join(f'{t:.2f}' for t in maps)} s; "
          f"nper {reader.nper}, numstp {reader.numstp}")


def main() -> None:
    with tempfile.TemporaryDirectory() as root, warnings.catch_warnings():
        warnings.filterwarnings("ignore", "Run 'r' built but its frozen snapshot", UserWarning)
        part1(Path(root))
    if len(sys.argv) > 1:
        part2(Path(sys.argv[1]))
    print(f"load average: {os.getloadavg()[0]:.1f}")


if __name__ == "__main__":
    main()

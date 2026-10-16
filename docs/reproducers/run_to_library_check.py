"""What ``run.to_library`` checks a run's record against, and why (ledger 613).

``run.to_library(project)`` registers a run's RECORD -- its library copies, its simulation
from the frozen snapshot, its mesh -- and first BUILDS and WRITES that import beside the
run's own files and compares them file by file. This script measures the three things that
chose that check, on a 5x5 one-layer DISV toy whose ``ic.strt`` and TDIS are edited after
the run wrote them:

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python \\
        docs/reproducers/run_to_library_check.py [--run <run directory>]

(a) ``ModelDiff`` between the edited run's files and the record written beside them: it
    compares options and list-BC cells, not array packages;
(b) the record compared IN MEMORY with the files read back: FloPy writes 8 decimals;
(c) ``to_library``'s own check: the rows it reports;
(d) an inline ``mf.disv(stack)`` in the frozen snapshot against its declaration -- why an
    inline package is taken from the manifest's recipe.

``--run`` imports a COPY of a real run (inputs and record only, never outputs) into an
empty scratch project and times it. It writes only under a temporary directory.

MEASURED (2026-10-16, FloPy 3.10.0)
  (a) ModelDiff: its package tier has one row, chd identical, and no ic row; its config
      tier differs on tdis period[0].perlen/nstp only -- the ic.strt edit is invisible
  (b) k = 1/3 in memory against the file's 0.33333333: 25 of 25 cells "differ"
  (c) to_library raises naming exactly m.ic (GRIDDATA) and s.tdis (PERIODDATA)
  (d) snapshot digest != declaration digest for the inline lazy disv; the recipe's matches
  --run on copies of the 23 Ten Trails runs with a run.json (20,323 cells x 5 layers):
      19 hold a snapshot and copies: 18 import "ok" and ss raises on ss.tdis PERIODDATA
      alone (perlen 1 vs 31, nstp 1 vs 10, changed outside myflopy after the write);
      4 hold no snapshot and no recipe and are refused; 3.3-4.8 s each.
"""

from __future__ import annotations

import argparse
import shutil
import tempfile
import time
import warnings
from pathlib import Path

import flopy
import numpy as np

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.run_record import _scan_simulation_inputs
from myflopy.specs import package_digest


def _project(root: Path, *, inline_disv: bool = False):
    grid = irregular_voronoi_grid(nrow=5, ncol=5, cell_size=100.0, seed=1)
    stack = mf.LayerStack(top=mf.Surface.flat(100.0))
    stack.add("a", thickness=10.0)
    project = mf.Project(root, name="toy")
    if not inline_disv:
        project.add_package("dis/base", mf.disv(stack))
    project.add_package("npf/base", mf.npf(k=1.0 / 3.0, icelltype=0))
    project.add_package("ic/base", mf.ic(strt=95.0))
    project.add_package("chd/base", mf.chd(stress_period_data={
        0: [[(0, 0), 95.0], [(0, grid.ncpl - 1), 94.0]]}))
    project.add_package("oc/all", mf.oc(head_filerecord="m.hds", saverecord=[("HEAD", "ALL")]))
    keys = ["npf/base", "ic/base", "chd/base", "oc/all"]
    packages = ([mf.disv(stack)] if inline_disv else [mf.ref("dis/base")]) + [
        mf.ref(key) for key in keys]
    project.add_simulation(mf.simulation(
        mf.gwf("m", context=mf.ModelContext(grid=grid), packages=packages), name="s",
        packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])]))
    return project


def checks(root: Path) -> None:
    project = _project(root / "toy")
    run = project.prepare_run("edited", "s")
    run.write(silent=True)
    record = project.prepare_run("record", "s")
    record.write(silent=True)
    # (b) the record in memory against the files read back.
    loaded = flopy.mf6.MFSimulation.load(sim_ws=str(run.workspace), verbosity_level=0)
    k_file = np.asarray(loaded.get_model("m").npf.k.array).ravel()
    k_memory = np.full(k_file.shape, 1.0 / 3.0)
    print(f"(b) in memory: {(k_file != k_memory).sum()} of {k_file.size} cells 'differ'"
          f" ({k_memory[0]!r} against {k_file[0]!r})")
    # The edits: an array through FloPy, then TDIS outside myflopy.
    gwf = run.simulation.get_model("m")
    gwf.ic.strt.set_data(90.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run.write(silent=True)
    tdis = run.workspace / "s.tdis"
    tdis.write_text(tdis.read_text().replace("1.00000000  1", "2.00000000  3", 1))
    # (a) ModelDiff between the two folders.
    diff = mf.load_mf6_run(run.workspace).diff(record.workspace)
    print(f"(a) ModelDiff summary: {diff.summary().to_dict('records')[:3]}")
    config = diff.config.summary()
    print(f"    config tier: {config.to_dict('records') if len(config) else '(empty)'}")
    settings = diff.config.settings()
    differing = settings[~settings["identical"]] if "identical" in settings else settings
    print(f"    differing settings: {sorted(set(map(str, differing.iloc[:, :3].values.tolist())))}")
    # (c) to_library's own check.
    fresh = mf.Project(root / "fresh", name="fresh")   # one Project per root
    try:
        project.runs["edited"].to_library(fresh)
    except ValueError as error:
        rows = [line.strip() for line in str(error).splitlines() if line.startswith("  ")]
        print(f"(c) to_library refused: {rows}")
    # (d) an inline lazy package in the snapshot against its declaration.
    inline = _project(root / "inline", inline_disv=True)
    declared = inline.simulations["s"].models[0].packages[0]
    run = inline.prepare_run("r", "s")
    run.write(silent=True)
    thawed = mf.load_run(run.workspace, load=False).spec.models[0].packages[0]
    rebuilt = mf.load_run(run.workspace, load=False)._declared_spec().models[0].packages[0]
    print(f"(d) declaration {package_digest(declared)}, snapshot {package_digest(thawed)}, "
          f"_declared_spec {package_digest(rebuilt)}")


def a_real_run(source: Path, root: Path) -> None:
    destination = root / "real" / "runs" / source.name
    destination.mkdir(parents=True)
    for relative in _scan_simulation_inputs(source)[0]:
        (destination / relative).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source / relative, destination / relative)
    shutil.copy2(source / "run.json", destination / "run.json")
    for sub in ("_library", ".myflopy"):
        if (source / sub).is_dir():
            shutil.copytree(source / sub, destination / sub)
    project = mf.Project(root / "real", name="real")
    start = time.perf_counter()
    try:
        report = project.runs[source.name].to_library(project)
        print(f"  {report} in {time.perf_counter() - start:.1f} s")
    except (ValueError, KeyError) as error:
        print(f"  refused in {time.perf_counter() - start:.1f} s: {str(error)[:400]}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run", type=Path, help="a run directory to import a COPY of")
    args = parser.parse_args()
    root = Path(tempfile.mkdtemp(prefix="run_to_library_"))
    try:
        checks(root)
        if args.run is not None:
            print(f"--run {args.run}:")
            a_real_run(args.run, root)
    finally:
        shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    main()

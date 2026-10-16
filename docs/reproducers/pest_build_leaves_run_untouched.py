"""A PEST build leaves the source run and the model object untouched -- and what the staging copy costs.

Ledger 356. ``PestProject.build()`` re-stored the caller's model (each array target
made layered, every array and list moved into an external file) and wrote it back
into the RUN directory, which ``PstFrom`` copied into the template; pyEMU's
``PstFrom.log`` went into the working directory. It now stages a copy of the run
beside the template and does all of that there. This script measures both, through
the public API, on whichever tree is on ``PYTHONPATH``:

    /home/lukem/python/envs/gw/bin/python docs/reproducers/pest_build_leaves_run_untouched.py \\
        [--large <run directory>] [--canonical full|validation|testing|none]

The "before" column is this script run against the parent commit's ``src``
(``git archive a58a8c3 src | tar -x -C <dir>``, then ``PYTHONPATH=<dir>/src``). Needs
``mf6`` (a 6.8 build: LAK ``IMPLICIT``). It writes only under a temporary directory,
which it removes; ``--large`` is COPIED there first and only the copy is built on.

THE SETUP. ``tests/test_pest_build_leaves_run_untouched.py``'s model: a 12 x 12
irregular Voronoi grid of 20 ft cells, two layers (100 - 82 - 40 ft) with a pond
excavated to 82 ft holding an ``IMPLICIT`` lake (inflow 5), NEWTON, scalar K 10 and K33
2 (so both are relayered, ledger 110), CHD 94 ft on the west and east edges of layer 2,
recharge on every seventh cell; one steady and one two-step transient period. The
calibration: ``k`` constant (``capture=True``), ``chd`` constant, ``k33`` pilot points at
three cells, one head target. Built from a scratch working directory.

  A. The ``ModelView`` the run hands back: the run directory hashed before and after,
     the model object's storage (every dataset of every package: layered, and each
     layer internal / constant / external) and simulation path, the working directory,
     ``IMPLICIT`` in the template's LAK file, and MF6 run directly in a copy of the
     template against the source run -- heads, cell budget, lake stage and lake budget
     compared byte for byte, and the template heads' SHA-256 printed so the two trees'
     templates can be compared with each other.
  B. The same model with NO lake inflow, so ``mf.lak`` writes empty period blocks.
  C. A copy of A's run whose arrays are already external (FloPy load + ``keep_implicit``
     + ``set_all_data_external`` + MF6 run), opened with ``mf.load_mf6_run``.
  D. Cost: ``cal.build()`` timed on a fresh copy of the canonical model's run (each
     build gets its own copy, because the old build writes into the run directory),
     ``load_mf6_run`` + ``k`` constant + one head target; the staging split (copy /
     load + externalize) is read from the build's own INFO log line. ``--large`` adds the
     same for a real run directory.
  E. (ledger 357, the review of 356) A's model run under a RELATIVE path from a scratch
     working directory -- every canonical PEST notebook's ``Path('../artifacts/...')`` --
     built into the DEFAULT, therefore relative, template; then built again, reporting
     whether the first build's ``PstFrom.log`` (in the template since 356) was closed
     before the rebuild deleted the template. "During" is d8b297d's ``src`` (ledger 356).

MEASURED (mf6 6.8.0.dev0, pyEMU 1.4.0, FloPy 3.10.0, 2026-10-02)       before          after
  A  run directory: files changed / added / removed                 12 / 22 / 0     0 / 0 / 0
     (7 of the 12 with their arrays moved out; the others re-stamped by FloPy -- only
     when the build writes in a later SECOND than the run did, so 7 to 12: the review
     measured 7 / 22 / 0)
     m.disv lines after the build (502 before it)                     30              502
     model object: storage entries changed                            18              0
     working directory after the build                           PstFrom.log        empty
     template LAK holds IMPLICIT                                      yes             yes
     template's MF6 = source run, byte for byte (5 files)             yes             yes
     template heads SHA-256                                     d6ff442748604d56 on both
  B  build                                                     CRASH MFDataException  built; = source,
                                                               (LAK "perioddata")     byte for byte
     run directory: changed / added / removed                       0 / 21 / 0      0 / 0 / 0
  C  build                                                     CRASH ValueError:      built; = source,
                                                               'k33' ... 'npf' ...    byte for byte
                                                               not there (pilot pts)
     run directory: changed / added / removed                       0 / 0 / 0       0 / 0 / 0
  D  canonical full (100x100x4, nper 6, 115 MB), two builds      31.2 s, 30.4 s   9.4 s, 9.2 s
       of which the copy / loading and externalizing it                              0.03 s / 8.4 s
     Ten Trails sw_mounding_flow_half (411 MB), --large          CRASH (load_mf6_run: 8.3 s, 8.2 s
                                                               no file for 'k')
       of which the copy / loading and externalizing it                              0.09 s / 5.4 s
     (load average 0.2-1.0; every copy on tmpfs)
     review (ledger 357), one build each, load average 0.1-1.3:   31.27 s          10.16 s
  E  (before / during: d8b297d / after: ledger 357)
     build                                                 before: built   during: CRASH
                                                           FileNotFoundError moving the log
                                                           to rel/p/runs/r.pest/s1 -- the
                                                           template was made INSIDE the
                                                           staging copy    after: built
     run directory: changed / added / removed                   7 / 22 / 0     0 / 0 / 0
     working directory after the build                  PstFrom.log, rel      rel
     template's MF6 = source run, byte for byte                    yes            yes
     second build: the first build's log closed first             no             yes

Why each is so, and what was decided: ledgers 356 and 357.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import logging
import os
import shutil
import subprocess
import tempfile
import time
import warnings
from pathlib import Path

import flopy
import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.lak_implicit import has_implicit, keep_implicit

OUTPUTS = ("m.hds", "m.cbc", "m_lak_stage.lak", "m_lak_budget.lak", "m_lak_budget.csv",
           "m.lst", "mfsim.lst")
ANSWERS = ("m.hds", "m.cbc", "m_lak_stage.lak", "m_lak_budget.lak", "m_lak_budget.csv")
PP_CELLS = [3, 70, 140]
MF6 = shutil.which("mf6") or "mf6"


def _simulation(*, inflow):
    grid = irregular_voronoi_grid(nrow=12, ncol=12, cell_size=20.0, refined=False)
    xmin, ymin, xmax, ymax = grid.gdf_vorPolys.total_bounds
    pond = gpd.GeoDataFrame({"name": ["pit"]}, crs=grid.crs, geometry=[
        Point((xmin + xmax) / 2, (ymin + ymax) / 2).buffer(30.0, 32)])
    result = (mf.LayerStack(top=mf.Flat(100.0)).add("upper", bottom=mf.Flat(82.0))
              .add("lower", bottom=mf.Flat(40.0))
              .excavate(bottom=82.0, within=pond, name="pit").build(grid))
    context = mf.ModelContext(grid=grid, domain=result.idomain, surfaces=result)
    xy = np.array([(c[1], c[2]) for c in grid.get_disv_gridprops()["cell2d"]])
    edge = np.flatnonzero((xy[:, 0] < xmin + 20) | (xy[:, 0] > xmax - 20))
    lake = mf.lak(context=context, nper=2, lakes=pond, lake_id_field="name",
                  starting_stage=95.0, bed_leakance=3000.0, implicit=True,
                  unexcavated="raise", inflow=inflow)
    recharge = [[(0, int(c)), 1e-3] for c in range(0, grid.ncpl, 7) if result.idomain[0][c] > 0]
    model = mf.gwf("m", context=context, packages=[
        mf.disv(result), mf.npf(k=10.0, k33=2.0, icelltype=1), mf.ic(strt=95.0),
        mf.chd(stress_period_data={0: [[(1, int(c)), 94.0] for c in edge]}),
        mf.sto(ss=1e-5, sy=0.2, iconvert=1, steady_state={0: True}, transient={1: True}),
        mf.rch(stress_period_data={0: recharge}),
        lake,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ], newtonoptions="NEWTON")
    return mf.simulation(model, name="s", packages=[
        mf.tdis(nper=2, perioddata=[(1.0, 1, 1.0), (10.0, 2, 1.0)])])


def _quietly(call, *args, **kwargs):
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        return call(*args, **kwargs)


def _run(root: Path, *, inflow):
    run = _quietly(mf.Project(root / "p", name="p").run, "r", _simulation(inflow=inflow))
    assert run.success
    return run


def _calibrate(model, workspace: Path | None):
    cal = model.pest("s1", workspace=workspace, start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), capture=True)
    cal.parameterize("chd", style="constant", bounds=(-1.0, 1.0))
    grid = model.gwf.modelgrid
    points = np.c_[np.asarray(grid.xcellcenters)[PP_CELLS], np.asarray(grid.ycellcenters)[PP_CELLS]]
    cal.parameterize("k33", style="pilotpoints", pp_points=points, bounds=(0.5, 2.0))
    cal.observe(mf.HeadTargets(
        locations=pd.DataFrame({"name": ["W"], "cell": [30], "layer": [1]}),
        values=pd.DataFrame({"time": [0, 1], "W": [94.5, 94.6]}),
    ))
    _quietly(cal.build, "s1.pst")
    return cal


def _tree(directory) -> dict[str, str]:
    directory = Path(directory)
    return {str(p.relative_to(directory)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(directory.rglob("*")) if p.is_file()}


def _tree_change(before, after) -> str:
    changed = sum(before[k] != after[k] for k in before.keys() & after.keys())
    return f"{changed} / {len(after.keys() - before.keys())} / {len(before.keys() - after.keys())}"


def _storage(simulation) -> dict:
    found = {}
    for model_name in simulation.model_names:
        for package in simulation.get_model(model_name).packagelist:
            for block in package.blocks.values():
                for name, data in block.datasets.items():
                    held = getattr(data, "_data_storage", None)
                    for key, storage in (held.items() if isinstance(held, dict) else [(None, held)]):
                        if storage is None:
                            continue
                        if hasattr(storage, "layer_storage"):
                            form = (bool(getattr(storage, "layered", False)), tuple(
                                str(layer.data_storage_type).rsplit(".", 1)[-1]
                                for layer in storage.layer_storage.elements()))
                        else:
                            form = (str(storage.data_storage_type).rsplit(".", 1)[-1],)
                        found[(model_name, package.filename, name, key)] = form
    return found


def _rerun(template: Path, destination: Path) -> Path:
    shutil.copytree(template, destination)
    for name in OUTPUTS:
        (destination / name).unlink(missing_ok=True)
    result = subprocess.run([MF6], cwd=destination, capture_output=True, text=True)
    assert "Normal termination" in result.stdout, result.stdout[-1500:]
    return destination


def _same_answers(source: Path, rerun: Path) -> str:
    same = [name for name in ANSWERS if (rerun / name).read_bytes() == (source / name).read_bytes()]
    return "yes" if len(same) == len(ANSWERS) else f"NO ({len(same)} of {len(ANSWERS)} equal)"


def _build_in(working: Path, model, workspace: Path):
    """Build from ``working`` as the working directory; (calibration or None, crash or None)."""

    working.mkdir()
    previous = os.getcwd()
    os.chdir(working)
    try:
        return _calibrate(model, workspace), None
    except Exception as error:  # noqa: BLE001 - a reproducer reports whatever the old tree raised
        return None, f"CRASH {type(error).__name__}: {str(error).splitlines()[0][:100]}"
    finally:
        os.chdir(previous)


def _report(label: str, value) -> None:
    print(f"   {label:<50}{value}")


def _report_template(cal, crash, source: Path, rerun: Path) -> None:
    if crash:
        _report("build", crash)
        return
    template = cal.template_workspace.resolve()
    rerun = _rerun(template, rerun)
    _report("template LAK holds IMPLICIT", has_implicit(template / "m.lak"))
    _report("template's MF6 = source run, byte for byte", _same_answers(source, rerun))
    _report("template heads SHA-256", hashlib.sha256((rerun / "m.hds").read_bytes()).hexdigest()[:16])
    _report("source heads SHA-256", hashlib.sha256((source / "m.hds").read_bytes()).hexdigest()[:16])


def section_a(root: Path):
    print("A. the ModelView a run hands back")
    run = _run(root / "a", inflow=5.0)
    model = run.model()
    workspace = Path(run.workspace).resolve()
    before, storage, path = _tree(workspace), _storage(model.sim), str(model.sim.sim_path)
    cal, crash = _build_in(root / "a_cwd", model, root / "a_pest")
    after_storage = _storage(model.sim)
    changed = sum(storage.get(k) != after_storage.get(k) for k in storage.keys() | after_storage.keys())
    _report("run directory: files changed / added / removed", _tree_change(before, _tree(workspace)))
    _report("m.disv lines after the build", len((workspace / "m.disv").read_text().splitlines()))
    _report("model object: storage entries changed", changed)
    _report("simulation path moved", str(model.sim.sim_path) != path)
    _report("working directory after the build", sorted(p.name for p in (root / "a_cwd").iterdir()) or "empty")
    _report_template(cal, crash, workspace, root / "a_rerun")
    return workspace


def section_b(root: Path):
    print("B. no lake inflow: mf.lak writes empty period blocks")
    run = _run(root / "b", inflow=None)
    workspace = Path(run.workspace).resolve()
    before = _tree(workspace)
    cal, crash = _build_in(root / "b_cwd", run.model(), root / "b_pest")
    _report("run directory: files changed / added / removed", _tree_change(before, _tree(workspace)))
    _report_template(cal, crash, workspace, root / "b_rerun")


def section_c(root: Path, source: Path):
    print("C. a run whose arrays are already external, opened with load_mf6_run")
    external = root / "c_external"
    shutil.copytree(source, external)
    simulation = flopy.mf6.MFSimulation.load(sim_ws=str(external), verbosity_level=0)
    keep_implicit(simulation)
    simulation.set_all_data_external()
    simulation.write_simulation(silent=True)
    for name in OUTPUTS:
        (external / name).unlink(missing_ok=True)
    subprocess.run([MF6], cwd=external, capture_output=True, text=True, check=True)
    model = mf.load_mf6_run(external, verbosity_level=0)
    before = _tree(external)
    cal, crash = _build_in(root / "c_cwd", model, root / "c_pest")
    _report("run directory: files changed / added / removed", _tree_change(before, _tree(external)))
    _report_template(cal, crash, external, root / "c_rerun")


def section_e(root: Path):
    print("E. relative paths, as every canonical notebook gives them; then a second build")
    working = root / "e_cwd"
    working.mkdir()
    previous = os.getcwd()
    os.chdir(working)
    try:
        run = _run(Path("rel"), inflow=5.0)
        model = run.model()
        _report("model workspace, as the run gives it", model.workspace)
        before = _tree(run.workspace)
        try:
            cal, crash = _calibrate(model, None), None
        except Exception as error:  # noqa: BLE001 - a reproducer reports whatever the tree raised
            # An OSError's second path is the destination it could not reach.
            detail = getattr(error, "filename2", None) or str(error).splitlines()[0][:100]
            cal, crash = None, f"CRASH {type(error).__name__}: {detail}"
        _report("run directory: files changed / added / removed", _tree_change(before, _tree(run.workspace)))
        _report("working directory after the build", sorted(p.name for p in working.iterdir()))
        # Resolved while the working directory is still the one the paths are relative to.
        _report_template(cal, crash, Path(run.workspace).resolve(), root / "e_rerun")
        if cal is not None:
            _report("default template (relative)", cal.template_workspace)
            first = getattr(cal.pf.logger, "f", None)
            _quietly(cal.build, "s1.pst")
            _report("second build: the first's log closed", getattr(first, "closed", None))
    finally:
        os.chdir(previous)


class _StagingLog(logging.Handler):
    """Collects the build's own INFO line saying how long the staging took."""

    def __init__(self):
        super().__init__(logging.INFO)
        self.lines: list[str] = []

    def emit(self, record):
        if "Staged" in record.getMessage():
            self.lines.append(record.getMessage())


def _timed_build(model, workspace: Path, *, observe_cell: int | None) -> float:
    """``cal.build()`` of ``k`` constant (and one head target), timed."""

    cal = model.pest("t", workspace=workspace, start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0))
    if observe_cell is not None:
        cal.observe(mf.HeadTargets(
            locations=pd.DataFrame({"name": ["W"], "cell": [observe_cell], "layer": [0]}),
            values=pd.DataFrame({"time": [0], "W": [np.nan]}),
        ))
    started = time.perf_counter()
    _quietly(cal.build, "t.pst")
    return time.perf_counter() - started


def _megabytes(directory: Path) -> float:
    return sum(p.stat().st_size for p in Path(directory).rglob("*") if p.is_file()) / 1e6


def section_d(root: Path, *, canonical: str, large: Path | None, repeat: int):
    print(f"D. build time; load average {os.getloadavg()[0]:.1f}")
    handler = _StagingLog()
    logging.getLogger("myflopy").addHandler(handler)
    logging.getLogger("myflopy").setLevel(logging.INFO)
    try:
        if canonical != "none":
            from myflopy.modflow.mf6.canonical_example import (
                CanonicalModelConfig,
                build_canonical_model,
            )

            config = {"full": CanonicalModelConfig, "validation": CanonicalModelConfig.validation,
                      "testing": CanonicalModelConfig.testing}[canonical]()
            label = f"canonical {canonical} ({config.nrow}x{config.ncol}x{config.nlay}, nper {config.nper})"
            for index in range(repeat):
                # A fresh model each time: the old build writes into the run directory.
                where = root / f"canonical_{index}"
                model = _quietly(build_canonical_model, where / "model", config=config)
                success, report = _quietly(model.run_simulation)
                assert success, report[-10:]
                size = _megabytes(model.workspace)
                elapsed = _timed_build(model, where / "pest" / "t", observe_cell=config.ncpl // 2)
                _report(f"{label} build {index + 1}", f"{elapsed:6.2f} s   (run directory {size:.0f} MB)")
                shutil.rmtree(where, ignore_errors=True)
        if large is not None:
            for index in range(repeat):
                local = root / f"large_{index}"
                started = time.perf_counter()
                shutil.copytree(large, local)
                _report(f"large {large.name}: copying it here (reads the disk)",
                        f"{time.perf_counter() - started:6.2f} s")
                model = mf.load_mf6_run(local, verbosity_level=0)
                try:
                    elapsed = _timed_build(model, root / f"large_{index}.pest" / "t", observe_cell=None)
                except Exception as error:  # noqa: BLE001 - a reproducer reports whatever the old tree raised
                    _report(f"large {large.name} build {index + 1}", f"CRASH {type(error).__name__}: "
                            f"{str(error).splitlines()[0][:90]}")
                else:
                    _report(f"large {large.name} build {index + 1}",
                            f"{elapsed:6.2f} s   (run directory {_megabytes(local):.0f} MB)")
                shutil.rmtree(local, ignore_errors=True)
                shutil.rmtree(root / f"large_{index}.pest", ignore_errors=True)
    finally:
        logging.getLogger("myflopy").removeHandler(handler)
    for line in handler.lines:
        print(f"      staging: {line.split(': ', 1)[-1]}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--large", type=Path, default=None, help="a real run directory to time (copied first)")
    parser.add_argument("--canonical", default="full", choices=["full", "validation", "testing", "none"])
    parser.add_argument("--sections", default="ABCDE")
    parser.add_argument("--repeat", type=int, default=2, help="timed builds per model in D")
    args = parser.parse_args()
    print(f"myflopy: {mf.__file__}")
    with tempfile.TemporaryDirectory(dir=os.environ.get("REPRO_TMP")) as scratch:
        root = Path(scratch)
        source = section_a(root) if "A" in args.sections else None
        if "B" in args.sections:
            section_b(root)
        if "C" in args.sections and source is not None:
            section_c(root, source)
        if "E" in args.sections:
            section_e(root)
        if "D" in args.sections:
            section_d(root, canonical=args.canonical, large=args.large, repeat=args.repeat)


if __name__ == "__main__":
    main()

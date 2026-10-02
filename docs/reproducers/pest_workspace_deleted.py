"""A PEST build deletes its ``workspace=`` only when a build wrote it there (ledger 360).

``model.pest(name, workspace=W)`` makes ``W`` the TEMPLATE directory, and ``build()``
deleted it WHATEVER it held -- ``_clear_stale_template`` ran ``shutil.rmtree(W)`` and
pyEMU's ``PstFrom(remove_existing=True)`` removed it again -- so a ``W`` that was an
existing folder went, with everything in it. This script drives each case through the
public API on whichever tree is on ``PYTHONPATH``:

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/pest_workspace_deleted.py

The "before" column is this script run against the parent commit's ``src``
(``git archive 5376c7c src | tar -x -C <dir>``, then ``PYTHONPATH=<dir>/src``). Needs
``mf6`` and pyEMU. It writes only under a temporary directory, which it removes.

THE SETUP. A 4 x 4 irregular Voronoi grid of 10 ft cells, one layer from 10 to 0 ft,
K 10 ft/d, CHD at 9 ft on the west edge and 8 ft on the east, one steady period, heads
saved; run through ``mf.Project(root / "p", name="p").run("r", ...)``, so the run
directory is ``root/p/runs/r``. The calibration: ``k`` constant and one head target.
Each case runs its own copy of that model under its own ``root``.

  A. ``W`` = the Project root ``root/p``, which holds the run (the reviewer's case).
  B. ``W`` = a data folder beside the project holding two files of the user's.
  C. ``W`` = ``root``, which contains the run, with a build marker FORGED in it
     naming it (the after tree's marker; on the before tree it is just a file).
  D. ``W`` = a template the previous build wrote there: built again.
  E. ``W`` = an empty folder.
  F. The staging directory beside ``W`` (``W.staging``) holding a file of the user's,
     ``W`` itself absent.

MEASURED (mf6 6.8.0.dev0, pyEMU 1.4.0, FloPy 3.10.0, 2026-10-02)      before                 after
  A  W = Project root holding the run    root/p DELETED (0 of 14 files      refused, ValueError, at
                                         left); build then raised           model.pest(); root/p
                                         FileNotFoundError "No MF6          byte-identical
                                         simulation in .../runs/r"
  B  W = a non-empty data folder         the user's 2 files DELETED, the     refused, FileExistsError;
                                         template built there (41 files)     byte-identical
  C  W = root, contains the run,         root DELETED (0 of 15 files, the    refused, ValueError, at
     marker forged                       run included); FileNotFoundError    model.pest(); byte-identical
  D  W = a previous build's template     rebuilt                             rebuilt
  E  W = an empty folder                 built                               built
  F  W.staging = a user folder           DELETED, then used and removed      refused, FileExistsError;
                                                                             byte-identical

THE REVIEW OF LEDGER 360 (ledger 361). "Before" is f4817ce's ``src`` (ledger 360), the
same way. I-L start a parallel run with ``true`` as the PEST++ executable -- it exits at
once, so pyEMU's ``start_workers`` does everything but run PEST -- and the results reader
stubbed.

  G. ``W`` = a template the previous build wrote, emptied with ``rm W/*`` (which skips
     dot-files) and refilled with two files of the user's.
  H. ``W`` = a ``cp -r`` backup of a template built before ledger 360 (no marker),
     kept under its own name, ``cal``, in another folder.
  I. ``run_ies(workers=2, master_dir=<a results folder holding a file of the user's>)``.
  J. ``run_ies(workers=2)`` with a folder of the user's named ``worker_1`` beside the
     default master.
  K. Two calibrations of one model: the second's ``prior(workers=2)``, with the first's
     ``worker_0`` (a copy of its template, as ``start_workers`` makes) beside it.
  L. ``run_ies(workers=2, master_dir=<the run directory>)``.

MEASURED (same versions, 2026-10-02)              before (f4817ce)                    after
  G  emptied with rm *, refilled         the 2 user files DELETED, template  refused, FileExistsError;
                                         rebuilt (dot marker survived)       byte-identical
  H  cp -r backup of a pre-360 template  rebuilt: 3 of its 42 files gone or  refused, FileExistsError;
                                         rewritten, cal.phi.actual.csv one   byte-identical
  I  master_dir = a results folder       the user's file DELETED, the        refused, FileExistsError;
                                         template copied in                  byte-identical
  J  a user's worker_1 beside it         DELETED                             refused, FileExistsError;
                                                                             byte-identical
  K  another calibration's worker_0      DELETED (0 of 42 files left)        refused, FileExistsError;
                                                                             byte-identical
  L  master_dir = the run directory      9 of its 14 files replaced by the   refused, ValueError;
                                         template's, 28 added                byte-identical
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import json
import shutil
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid


def _simulation():
    grid = irregular_voronoi_grid(nrow=4, ncol=4, cell_size=10.0, refined=False)
    props = grid.get_disv_gridprops()
    x = np.array([cell[1] for cell in props["cell2d"]])
    west, east = np.flatnonzero(x < x.min() + 5.0), np.flatnonzero(x > x.max() - 5.0)
    chd = [[(0, int(c)), 9.0] for c in west] + [[(0, int(c)), 8.0] for c in east]
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=10.0, botm=[0.0]),
        mf.npf(k=10.0, icelltype=0), mf.ic(strt=8.5),
        mf.chd(stress_period_data={0: chd}),
        mf.oc(head_filerecord="m.hds", saverecord=[("HEAD", "ALL")]),
    ])
    return mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])


def _quietly(call, *args, **kwargs):
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        return call(*args, **kwargs)


def _run(root: Path):
    run = _quietly(mf.Project(root / "p", name="p").run, "r", _simulation())
    assert run.success
    return run


def _build(model, workspace: Path | None, *, name: str = "cal"):
    cal = model.pest(name, workspace=workspace, start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0))
    cal.observe(mf.HeadTargets(
        locations=pd.DataFrame({"name": ["W"], "cell": [5], "layer": [0]}),
        values=pd.DataFrame({"time": [0], "W": [8.6]}),
    ))
    _quietly(cal.build, f"{name}.pst", noptmax=0)
    return cal


def _tree(directory: Path) -> dict[str, str]:
    directory = Path(directory)
    if not directory.exists():
        return {}
    return {str(p.relative_to(directory)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(directory.rglob("*")) if p.is_file()}


def _state(watched: Path, before: dict[str, str]) -> str:
    """What became of ``watched``, whose files were ``before``: kept means the same bytes."""

    after = _tree(watched)
    if not watched.exists():
        return f"{watched.name} DELETED (0 of {len(before)} files left)"
    if before == after:
        return f"{watched.name} byte-identical ({len(after)} files)"
    kept = sum(after.get(path) == digest for path, digest in before.items())
    return f"{watched.name} CHANGED ({kept} of {len(before)} files kept, {len(after)} now)"


def _attempt(model, workspace: Path, watched: Path, scratch: Path) -> str:
    """Build at ``workspace``; say what happened to it and to ``watched``."""

    before = _tree(watched)
    try:
        cal = _build(model, workspace)
    except Exception as error:  # noqa: BLE001 - a reproducer reports whatever either tree raised
        first = str(error).splitlines()[0].replace(str(scratch), "<tmp>")
        outcome = f"{type(error).__name__}: {first[:150]}"
    else:
        outcome = f"built ({len(cal.pst.par_names)} par, {cal.pst.nobs} obs)"
    return f"{outcome}\n       {_state(watched, before)}"


def _attempt_run(method, watched: Path, scratch: Path, **kwargs) -> str:
    """Start a parallel run with ``true`` as PEST++ (it exits at once, so pyEMU's
    ``start_workers`` does everything but run PEST); say what happened to ``watched``."""

    import myflopy.modflow.mf6.pest.ies as ies

    reader, ies.IesResults = ies.IesResults, (lambda workspace, **_: Path(workspace))
    before = _tree(watched)
    try:
        _quietly(method, exe="true", **kwargs)
    except Exception as error:  # noqa: BLE001 - a reproducer reports whatever either tree raised
        first = str(error).splitlines()[0].replace(str(scratch), "<tmp>")
        outcome = f"{type(error).__name__}: {first[:150]}"
    else:
        outcome = "ran"
    finally:
        ies.IesResults = reader
    return f"{outcome}\n       {_state(watched, before)}"


def main():
    with tempfile.TemporaryDirectory(prefix="pest_workspace_") as scratch:
        scratch = Path(scratch)

        root = scratch / "A"
        model = _run(root).model()
        print("A  W = the Project root holding the run\n    ", _attempt(model, root / "p", root / "p", scratch))

        root = scratch / "B"
        model = _run(root).model()
        data = root / "data"
        data.mkdir()
        (data / "wells.csv").write_text("name,x,y\nW,1,2\n", encoding="utf-8")
        (data / "notes.txt").write_text("field notes\n", encoding="utf-8")
        print("B  W = a non-empty data folder\n    ", _attempt(model, data, data, scratch))

        root = scratch / "C"
        model = _run(root).model()
        (root / "myflopy_pest_build.json").write_text(json.dumps(
            {"written_by": "myflopy PEST build", "role": "template",
             "directory": str(root.resolve())}), encoding="utf-8")
        print("C  W = root (contains the run), marker forged\n    ", _attempt(model, root, root, scratch))

        root = scratch / "D"
        model = _run(root).model()
        template = root / "tmpl"
        _build(model, template)
        print("D  W = a template the previous build wrote\n    ", _attempt(model, template, root / "p", scratch))

        root = scratch / "E"
        model = _run(root).model()
        empty = root / "empty"
        empty.mkdir()
        print("E  W = an empty folder\n    ", _attempt(model, empty, root / "p", scratch))

        root = scratch / "F"
        model = _run(root).model()
        staging = root / "tmpl.staging"
        staging.mkdir()
        (staging / "mine.txt").write_text("the user's\n", encoding="utf-8")
        print("F  W.staging = a user folder\n    ", _attempt(model, root / "tmpl", staging, scratch))

        # -- the review of ledger 360 (ledger 361): G-H the build, I-L the parallel run --

        root = scratch / "G"
        model = _run(root).model()
        template = root / "calib"
        _build(model, template)
        for entry in template.iterdir():
            if not entry.name.startswith("."):   # what `rm calib/*` removes
                shutil.rmtree(entry) if entry.is_dir() else entry.unlink()
        (template / "survey.csv").write_text("name,head\nW,8.6\n", encoding="utf-8")
        (template / "report.txt").write_text("the user's report\n", encoding="utf-8")
        print("G  W = a template emptied with `rm W/*`, refilled with 2 user files\n    ",
              _attempt(model, template, template, scratch))

        root = scratch / "H"
        model = _run(root).model()
        default = _build(model, None).template_workspace
        for marker in default.glob("*myflopy_pest_build.json"):
            marker.unlink()   # a template built before ledger 360
        backup = root / "backups" / "cal"
        shutil.copytree(default, backup)
        (backup / "cal.phi.actual.csv").write_text("iteration,mean\n0,1.0\n", encoding="utf-8")
        print("H  W = a `cp -r` backup of a pre-360 template, under its own name elsewhere\n    ",
              _attempt(model, backup, backup, scratch))

        root = scratch / "I"
        model = _run(root).model()
        cal = _build(model, None)
        results = root / "results"
        results.mkdir()
        (results / "earlier_run.csv").write_text("a,b\n1,2\n", encoding="utf-8")
        print("I  run_ies(workers=2, master_dir=<a results folder of the user's>)\n    ",
              _attempt_run(cal.run_ies, results, scratch, workers=2, master_dir=results))

        root = scratch / "J"
        model = _run(root).model()
        cal = _build(model, None)
        worker = cal.template_workspace.parent / "worker_1"
        worker.mkdir()
        (worker / "mine.txt").write_text("the user's\n", encoding="utf-8")
        print("J  run_ies(workers=2), a user folder `worker_1` beside the default master\n    ",
              _attempt_run(cal.run_ies, worker, scratch, workers=2))

        root = scratch / "K"
        model = _run(root).model()
        first = _build(model, None, name="first")
        worker = first.template_workspace.parent / "worker_0"
        shutil.copytree(first.template_workspace, worker)   # what start_workers makes of it
        second = _build(model, None, name="second")
        print("K  a 2nd calibration's prior(workers=2) beside the 1st's worker_0\n    ",
              _attempt_run(second.prior, worker, scratch, workers=2))

        root = scratch / "L"
        model = _run(root).model()
        cal = _build(model, None)
        run = Path(model.workspace)
        print("L  run_ies(workers=2, master_dir=<the run directory>)\n    ",
              _attempt_run(cal.run_ies, run, scratch, workers=2, master_dir=run))


if __name__ == "__main__":
    main()

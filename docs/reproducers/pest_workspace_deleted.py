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
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import json
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


def _build(model, workspace: Path):
    cal = model.pest("cal", workspace=workspace, start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0))
    cal.observe(mf.HeadTargets(
        locations=pd.DataFrame({"name": ["W"], "cell": [5], "layer": [0]}),
        values=pd.DataFrame({"time": [0], "W": [8.6]}),
    ))
    _quietly(cal.build, "cal.pst", noptmax=0)
    return cal


def _tree(directory: Path) -> dict[str, str]:
    directory = Path(directory)
    if not directory.exists():
        return {}
    return {str(p.relative_to(directory)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(directory.rglob("*")) if p.is_file()}


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
    after = _tree(watched)
    if not watched.exists():
        state = f"{watched.name} DELETED (0 of {len(before)} files left)"
    elif before == after:
        state = f"{watched.name} byte-identical ({len(after)} files)"
    else:
        kept = len(before.keys() & after.keys())
        state = f"{watched.name} CHANGED ({kept} of {len(before)} files kept, {len(after)} now)"
    return f"{outcome}\n       {state}"


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
        (root / ".myflopy_pest_build.json").write_text(json.dumps(
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


if __name__ == "__main__":
    main()

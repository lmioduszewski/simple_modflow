"""A run's record against its files: what MF6 runs, and what each check costs (ledger 611).

A run folder holds the RECORD (``run.spec``, its frozen snapshot, its library copies) and
the MF6 input FILES. This script runs the six scenarios that decided the rule on a 2x2
steady model with no boundary -- its heads equal ``ic.strt``, so the head file says which
``strt`` MF6 ran -- and, given ``--run``, times the checks on a COPY of a real run:

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python \\
        docs/reproducers/run_record_vs_files.py [--run <run directory>]

It writes only under a temporary directory, which it removes; ``--run`` is copied there
first and never touched.

MEASURED (2026-10-16, FloPy 3.10.0, MF6 6.x; "before" is 55947911)     before    after
  S1 completed run, strt edited through FloPy, execute()            ran 1.0   ran 50.0
  S2 then write(): record says 1.0, files 50; marked?               no        run.inputs: edited (flow.ic)
  S4 load_run(load=False).execute() after S2                        ran 1.0, rewrote flow.ic   ran 50.0
  S5 project.runs["baseline"].execute() after S2                    ran 1.0, rewrote flow.ic   ran 50.0
  S6 flow.ic hand-edited to 77; load_run(load=False).execute()      ran 1.0, rewrote flow.ic   ran 77.0
  S6' write() with nothing held over those files                    rebuilt the record         refused
  --run on a copy of Ten Trails sw_mounding_flow (38 periods, 20,323 cells x 5 layers):
    hashing the 17 input files 0.07 s; fingerprint 0.10-0.24 s built / 0.06-0.11 s
    loaded; execute()'s check 0.19-0.40 s; write() +0.2-0.4 s on FloPy's ~2.3 s; full
    check 3.3-4.2 s, 0 rows (three runs on a shared machine).
    The input scan cost 1.3 s per read until it looked for FILEIN/OPEN/CLOSE first.
"""

from __future__ import annotations

import argparse
import json
import shutil
import tempfile
import time
import warnings
from pathlib import Path

import flopy

import myflopy as mf
from myflopy import workspace as workspace_module


def _simulation() -> mf.SimulationSpec:
    return mf.SimulationSpec(
        "baseline",
        models=[mf.gwf("flow", packages=[
            mf.dis(nlay=1, nrow=2, ncol=2, delr=1.0, delc=1.0, top=10.0, botm=0.0),
            mf.ic(strt=1.0), mf.ref("npf/base"),
            mf.oc(head_filerecord="flow.hds", saverecord=[("HEAD", "ALL")]),
        ])],
        packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]), mf.ims(models=["flow"])],
    )


def _fresh(root: Path, tag: str):
    project = mf.Project(root / tag, name="p")
    project.add_package("npf/base", mf.npf(k=1.0))
    run = project.prepare_run("baseline", _simulation())
    assert run.execute()[0]
    return project, run


def _ran(run) -> float:
    heads = flopy.utils.HeadFile(run.workspace / "flow.hds")
    try:
        return float(heads.get_data()[0, 0, 0])
    finally:
        heads.close()


def _strt_on_disk(run) -> str:
    return next(line.strip() for line in (run.workspace / "flow.ic").read_text().splitlines()
                if "CONSTANT" in line)


def scenarios(root: Path) -> None:
    project, run = _fresh(root, "s1")
    run.simulation.get_model("flow").ic.strt.set_data(50.0)
    run.execute()
    print(f"S1  edited in memory, execute(): MF6 ran {_ran(run)}")
    print(f"S2  run.inputs after that write: {run.inputs}")
    reopened = mf.load_run(run.workspace, load=False)
    reopened.execute()
    print(f"S4  load_run(load=False).execute(): ran {_ran(reopened)}, file {_strt_on_disk(run)}")
    entry = mf.Project(root / "s1", name="p").runs["baseline"]
    entry.execute()
    print(f"S5  project.runs['baseline'].execute(): ran {_ran(entry)}")

    _, run6 = _fresh(root, "s6")
    path = run6.workspace / "flow.ic"
    path.write_text(path.read_text().replace("1.00000000", "77.00000000"))
    reopened = mf.load_run(run6.workspace, load=False)
    print(f"S6  hand edit: {reopened.inputs}")
    try:
        reopened.write(silent=True)
        print("S6' write() with nothing held: rebuilt the record over the files")
    except ValueError as error:
        print(f"S6' write() with nothing held: refused ({str(error)[:80]}...)")
    reopened.execute()
    print(f"S6  load_run(load=False).execute(): ran {_ran(reopened)}, file {_strt_on_disk(run6)}")


def costs(source: Path, root: Path) -> None:
    copy = root / "run_copy"
    shutil.copytree(source, copy, symlinks=True)

    def timed(label, function, repeat=1):
        start = time.perf_counter()
        for _ in range(repeat):
            out = function()
        print(f"  {label}: {(time.perf_counter() - start) / repeat:.3f} s")
        return out

    hashes = timed("hash every input file", lambda: workspace_module._input_file_hashes(copy), 3)
    print(f"    ({len(hashes)} input files)")
    loaded = timed("load_run(load=True)", lambda: mf.load_run(copy))
    timed("fingerprint, loaded", lambda: workspace_module._simulation_fingerprint(
        loaded.simulation), 3)
    unfingerprinted = mf.load_run(copy, load=False)
    print(f"    run.inputs before any write here: {unfingerprinted.inputs.verdict}")
    rows = timed("full check (record built + written + compared)",
                 unfingerprinted._record_differences)
    print(f"    {len(rows)} difference row(s): {rows[:3]}")
    built = mf.load_run(copy, load=False)
    timed("run.build()", built.build)
    timed("fingerprint, built", lambda: workspace_module._simulation_fingerprint(
        built.simulation), 3)
    timed("run.write()", lambda: built.write(silent=True))
    print(f"    run.inputs after it: {built.inputs.verdict}")
    timed("execute()'s own check (fingerprint + hashes)", lambda: (
        workspace_module._simulation_fingerprint(built.simulation),
        workspace_module._input_file_hashes(copy)))
    manifest = json.loads((copy / "run.json").read_text())
    print(f"    recorded {len(manifest['metadata']['inputs']['sha'])} input hashes")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run", type=Path, default=None,
                        help="a run directory to time the checks on (a copy is used)")
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="record-vs-files-") as scratch, \
            warnings.catch_warnings():
        warnings.simplefilter("ignore")
        root = Path(scratch)
        scenarios(root)
        if args.run is not None:
            print(f"costs on a copy of {args.run}:")
            costs(args.run, root)


if __name__ == "__main__":
    main()

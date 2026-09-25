"""What PEST++ does with a forward-run CSV that is not shaped like its base run.

Ledger 329. ``LakeStageTargets.simulated_series`` says ``NaN`` where a lake has no
stage, and that series is what ``write_named_series_targets`` writes for PEST. pyEMU
builds the instruction file from the BASE run's CSV (``pst_utils.csv_to_ins_file``):
one ``l1`` line advance per row, and each value found by counting commas. So the
forward run must match the base run row for row and column for column, and every
field must be a number PEST++ can read. This runs PEST++ on each way that can fail.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/pest_no_value_fields.py

THE SETUP. A one-parameter "model" (a shell script) that writes ``out.csv``. The base
file, from which the instruction file is built, is always::

    time,foo
    0,100.0
    1,101.0
    2,102.0

and each case's forward run writes something else. ``pestpp-glm`` with
``NOPTMAX 0`` runs it once. Needs ``pestpp-glm`` on PATH.

MEASURED (pestpp-glm, pyEMU 1.4.0, 2026-09-25):
  empty field   ``1,``                      run FAILS: "error tokenizing output line"
  nan           ``1,nan``                   run FAILS: "casting 'nan' to double yielded
                                            denormal value"
  1e30          ``1,1e+30``                 runs; residual -1e30, phi 1e60
  one row short (row 1 dropped)             run FAILS: "EOF encountered when executing
                                            line advance instruction"
  extra column  ``time,bar,foo`` (bar=7)    runs, and reads ``bar`` as ``foo``: every
                                            modelled value 7.0, no error
"""

from __future__ import annotations

import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile

import pyemu

BASE = "time,foo\n0,100.0\n1,101.0\n2,102.0\n"
CASES = {
    "empty field": "time,foo\n0,100.0\n1,\n2,102.0\n",
    "nan": "time,foo\n0,100.0\n1,nan\n2,102.0\n",
    "1e30": "time,foo\n0,100.0\n1,1e+30\n2,102.0\n",
    "one row short": "time,foo\n0,100.0\n2,102.0\n",
    "extra column": "time,bar,foo\n0,7.0,100.0\n1,7.0,101.0\n2,7.0,102.0\n",
}


def _case(root: pathlib.Path, label: str, forward: str) -> None:
    ws = root / label.replace(" ", "_")
    ws.mkdir()
    (ws / "out.csv").write_text(BASE)
    (ws / "par.dat").write_text("1.0\n")
    (ws / "par.dat.tpl").write_text("ptf ~\n~ p1 ~\n")
    (ws / "forward.csv").write_text(forward)
    # PEST++ deletes the model's output files before each run; the "model" restores one.
    (ws / "run.sh").write_text("#!/bin/sh\ncp forward.csv out.csv\n")
    os.chmod(ws / "run.sh", 0o755)
    pyemu.pst_utils.csv_to_ins_file(str(ws / "out.csv"), ins_filename=str(ws / "out.csv.ins"),
                                    prefix="s")
    pst = pyemu.Pst.from_io_files([str(ws / "par.dat.tpl")], [str(ws / "par.dat")],
                                  [str(ws / "out.csv.ins")], [str(ws / "out.csv")], pst_path=".")
    pst.model_command = ["./run.sh"]
    pst.control_data.noptmax = 0
    pst.observation_data.loc[:, "weight"] = 1.0
    pst.write(str(ws / "t.pst"))
    done = subprocess.run(["pestpp-glm", "t.pst"], cwd=ws, capture_output=True, text=True,
                          timeout=120)
    print(f"== {label}: {'runs' if done.returncode == 0 else 'run FAILS'}")
    reasons = sorted({match.group(1).strip()
                      for match in re.finditer(r"output file line: \d+ : ([^\n]*)", done.stdout)})
    for reason in reasons:
        print(f"   {reason[:140]}")
    rei = ws / "t.base.rei" if (ws / "t.base.rei").exists() else ws / "t.rei"
    if done.returncode == 0 and rei.exists():
        rows = [line.split() for line in rei.read_text().splitlines()
                if line.strip().startswith("s_usecol")]
        for name, _, measured, modelled, *_ in rows:
            print(f"   {name}: measured {float(measured):g}, modelled {float(modelled):g}")


def main() -> None:
    if shutil.which("pestpp-glm") is None:
        sys.exit("pestpp-glm is not on PATH")
    root = pathlib.Path(tempfile.mkdtemp(prefix="pest_no_value_fields_"))
    for label, forward in CASES.items():
        _case(root, label, forward)


if __name__ == "__main__":
    main()

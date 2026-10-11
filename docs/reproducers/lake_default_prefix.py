"""A lake set's default prefix takes a PESTPP-IES phi-factor tag of its own: measured.

Ledger 372. A target set with no group label is in its set's PREFIX as its PEST
observation group (ledger 358), and PESTPP-IES maps an ``ies_phi_factor_file`` tag to
every group whose name CONTAINS it. A lake set's default prefix was ``stage``, inside an
SFR stage set's ``sfr_stage``, so with both registered under their defaults no tag could
match the lake's group alone (ledger 359 measured the refusal). The default is
``lak_stage`` now. This script registers heads, a lake and a reach-stage set under their
DEFAULT prefixes, tags each measured group by its own name, and runs PESTPP-IES.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/lake_default_prefix.py

It uses only the public API, so it measures whichever tree is on ``PYTHONPATH`` -- the
"before" column below is this script run against the parent commit's ``src``
(``git archive <parent> src | tar -x -C <dir>``, then ``PYTHONPATH=<dir>/src``) -- an
ABSOLUTE path: PESTPP-IES runs the forward run in the template, where a relative one
resolves to nothing and the venv's editable install is imported instead. Needs
``mf6`` and ``pestpp-ies``. It reuses ``pest_observation_groups.py``'s model (imported
from it): a 6 x 6 irregular Voronoi grid, one confined layer, a four-reach stream and a
two-cell lake, three 10-day periods; one parameter, ``k``. Heads W2/W1/W10 unlabelled,
the lake ``pond`` measured at periods 1 and 2, reach 1 ``r1`` at periods 0 and 1.

MEASURED (mf6 6.8.0.dev0, PEST++ 5.2.16, pyEMU 1.4.0, 2026-10-11)
                                     before (parent)               after
  control-file groups                hds, hds_unmeasured, stage,   hds, hds_unmeasured,
                                     stage_unmeasured, sfr_stage,  lak_stage,
                                     sfr_stage_unmeasured          lak_stage_unmeasured,
                                                                   sfr_stage,
                                                                   sfr_stage_unmeasured
  tags = the measured groups         hds, sfr_stage, stage         hds, lak_stage, sfr_stage
  pestpp-ies NOPTMAX=1, 4 reals      exit 1, "group 'SFR_STAGE'    exit 0
                                     mapped to multiple tags:
                                     SFR_STAGE STAGE"
"""

from __future__ import annotations

import pathlib
import shutil
import subprocess
import sys
import tempfile

import pandas as pd

import myflopy as mf

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from pest_observation_groups import IES, _built, _calibration, _heads, _toy  # noqa: E402


def main():
    root = pathlib.Path(tempfile.mkdtemp(prefix="lake_default_prefix_"))
    try:
        model = _toy(root)
        cal = _calibration(model, root / "pest", "d")
        cal.observe(_heads([None, None, None]))
        cal.observe(mf.LakeStageTargets(locations={"pond": 0}, values=pd.DataFrame(
            {"time": [1.0, 2.0], "pond": [99.0, 98.0]})))
        cal.observe(mf.SfrStageTargets(locations={"r1": 1}, values=pd.DataFrame(
            {"time": [0.0, 1.0], "r1": [95.0, 96.0]})))
        pst = _built(cal, "d")
        obs = pst.observation_data
        groups = sorted(set(obs["obgnme"]))
        measured = sorted(set(obs.loc[obs["weight"] > 0, "obgnme"]))
        print(f"control-file groups: {groups}")
        print(f"tags (the measured groups): {measured}")
        work = pathlib.Path(cal.template_workspace)
        share = 1.0 / len(measured)
        (work / "phifac.csv").write_text("".join(f"{tag},{share}\n" for tag in measured))
        pst.control_data.noptmax = 1
        pst.pestpp_options.update({"ies_num_reals": 4, "ies_phi_factor_file": "phifac.csv"})
        pst.write(str(work / "d.pst"), version=2)
        done = subprocess.run([IES, "d.pst"], cwd=work, capture_output=True, text=True,
                              timeout=900)
        rec = (work / "d.rec").read_text() if (work / "d.rec").exists() else ""
        said = sorted({line.strip() for line in (rec + done.stdout).splitlines()
                       if "tags" in line and "rror" in line})
        print(f"pestpp-ies NOPTMAX=1, 4 reals: exit {done.returncode} {said[:1]}")
    finally:
        shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    main()

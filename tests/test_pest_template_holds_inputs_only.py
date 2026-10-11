"""A PEST template holds the simulation's inputs and nothing else of the run (ledger 383).

The build's staging copy (ledgers 356-357) copied the WHOLE run directory, and
``PstFrom`` copied that into the template -- so every template, and every IES worker
pyEMU copies from it, carried the run's MF6 outputs (head, budget and listing files:
3.2 MB of the canonical run's 3.9 MB), myflopy's run metadata (``run.json``, a pickled
model) and any PEST workspace nested in the run directory (a pre-ledger-107 build put
``<run>/pest/`` there). None of it is read: the build reads every base observation
from the RUN directory's output files, and each realization's model writes its own.

The staging copy is now what ``mfsim.nam`` leads MF6 to read -- name-file entries,
``FILEIN`` and ``OPEN/CLOSE`` -- and the end-to-end test proves, with a real
PESTPP-IES prior run (NOPTMAX=-1) of ``test_pest_handed_targets.py``'s toy, that the
forward run in such a template reproduces the base observations the build read from
the run, in every family the toy carries.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6.pest.project import BUILD_MARKER, _simulation_input_files
from tests.test_pest_handed_targets import DRN_CELLS, run_toy

pyemu = pytest.importorskip("pyemu")


def _write(path: Path, text: str = "x\n") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def test_the_input_files_are_what_mfsim_nam_leads_mf6_to_read(tmp_path):
    """Name-file entries, FILEIN (scanned in turn) and OPEN/CLOSE (copied unread), quoted
    names included; never LIST, FILEOUT, a name nothing gives, an absolute name or one
    outside the directory."""

    _write(tmp_path / "mfsim.nam", """# a comment naming FILEIN not_me.txt
BEGIN options
  HPC6 FILEIN sim.hpc
END options
BEGIN timing
  TDIS6 s.tdis
END timing
BEGIN models
  gwf6 'sub dir/m.nam' m
END models
BEGIN exchanges
END exchanges
BEGIN solutiongroup 1
  MXITER 2
  ims6 m.ims m
END solutiongroup
""")
    _write(tmp_path / "sim.hpc")
    _write(tmp_path / "s.tdis", "BEGIN options\n  ATS6 FILEIN s.ats\nEND options\n")
    _write(tmp_path / "s.ats")
    _write(tmp_path / "m.ims")
    _write(tmp_path / "sub dir/m.nam", """BEGIN options
  LIST m.lst
END options
BEGIN packages
  NPF6 m.npf npf
  CHD6 m.chd chd  ! trailing comment
  OBS6 m.obs obs
END packages
""")
    _write(tmp_path / "m.npf", "BEGIN griddata\n  k\n    OPEN/CLOSE  'arrays/k 1.txt'  FACTOR 1.0\nEND griddata\n")
    _write(tmp_path / "arrays/k 1.txt", "FILEIN never_scanned.txt\n")
    _write(tmp_path / "m.chd", ("BEGIN options\n  TS6 filein m.ts\n  OBS6 FILEIN m.chd.obs\n"
                                "  head_filerecord FILEOUT m.hds\nEND options\n"
                                "BEGIN period 1\n  open/close m.chd_1.txt\nEND period\n"
                                "  OPEN/CLOSE /absolute/elsewhere.txt\n  OPEN/CLOSE ../outside.txt\n"))
    _write(tmp_path / "m.ts", "BEGIN timeseries\nEND timeseries\n")
    _write(tmp_path / "m.chd.obs", "BEGIN continuous FILEOUT m.chd.csv\nEND continuous\n")
    _write(tmp_path / "m.chd_1.txt")
    _write(tmp_path / "m.obs", "BEGIN continuous FILEOUT heads.csv\nEND continuous\n")
    (tmp_path / "m.hds").write_bytes(b"\0" * 64)
    for stray in ("m.lst", "m.chd.csv", "heads.csv", "never_scanned.txt", "not_me.txt",
                  "run.json", "pest/old.pst", "2"):
        _write(tmp_path / stray)
    _write(tmp_path.parent / "outside.txt")

    assert [str(path) for path in _simulation_input_files(tmp_path)] == [
        "mfsim.nam", "s.tdis", "sub dir/m.nam", "m.ims", "sim.hpc", "s.ats", "m.npf",
        "m.chd", "m.obs", "arrays/k 1.txt", "m.ts", "m.chd.obs", "m.chd_1.txt",
    ]


def test_a_binary_input_is_copied_and_not_scanned(tmp_path):
    """FMI's flow files and a binary array are inputs too; their bytes are never read
    for names."""

    _write(tmp_path / "mfsim.nam", "BEGIN models\n  gwt6 t.nam t\nEND models\n")
    _write(tmp_path / "t.nam", "BEGIN packages\n  FMI6 t.fmi fmi\nEND packages\n")
    _write(tmp_path / "t.fmi", "BEGIN packagedata\n  GWFHEAD FILEIN flow.hds\nEND packagedata\n")
    (tmp_path / "flow.hds").write_bytes(b"\0FILEIN ghost.txt\0")
    _write(tmp_path / "ghost.txt")

    assert [str(path) for path in _simulation_input_files(tmp_path)] == [
        "mfsim.nam", "t.nam", "t.fmi", "flow.hds",
    ]


# -- end to end: the toy, built and run through PESTPP-IES -------------------------------


def _targets():
    heads = mf.HeadTargets(
        locations=pd.DataFrame({"name": ["W", "V"], "cell": [22, 3], "layer": [0, 0]}),
        values=pd.DataFrame({"time": [0, 1, 2], "W": [95.0, 96.0, 97.0], "V": [99.0, None, 98.0]}))
    flow = mf.SfrFlowTargets(locations={"a": 0}, values=pd.DataFrame({"time": [1], "a": [-7.0]}))
    seep = mf.DrnFlowTargets(locations={"zone": DRN_CELLS},
                             values=pd.DataFrame({"time": [2], "zone": [-3.0]}))
    pond = mf.LakeStageTargets(locations={"pond": 0},
                               values=pd.DataFrame({"time": [1], "pond": [99.0]}))
    return heads, flow, seep, pond


@pytest.fixture(scope="module")
def staged(tmp_path_factory):
    """The toy run, holding what a run directory collects besides its inputs, built."""

    root = tmp_path_factory.mktemp("inputs_only")
    model = run_toy(root)
    run = Path(model.workspace)
    # A PEST workspace from before ledger 107, a build's template, and a file of the user's.
    _write(run / "pest" / "old.pst", "pcf\n")
    (run / "pest" / "old.bin").write_bytes(b"\0" * 100_000)
    _write(run / "tpl" / BUILD_MARKER, "{}\n")
    _write(run / "tpl" / "tpl.pst", "pcf\n")
    _write(run / "notes.txt", "mine\n")
    cal = model.pest("inputs", workspace=root / "pest", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    for targets in _targets():
        cal.observe(targets)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build("inputs.pst", noptmax=0)
    return model, cal


def test_the_template_holds_no_output_metadata_or_pest_workspace_of_the_run(staged):
    model, cal = staged
    run, template = Path(model.workspace), cal.template_workspace
    inputs = {str(path) for path in _simulation_input_files(run)}
    in_run = {str(path.relative_to(run)) for path in run.rglob("*") if path.is_file()}
    in_template = {str(path.relative_to(template)) for path in template.rglob("*") if path.is_file()}

    left_behind = in_run - inputs
    # What a run collects besides its inputs, every kind of it.
    assert {"m.hds", "m.cbc", "m.lst", "mfsim.lst", "m.lak.stg", "m.sfr.stg", "m.disv.grb",
            "run.json", "pest/old.pst", "pest/old.bin", f"tpl/{BUILD_MARKER}", "tpl/tpl.pst",
            "notes.txt"} <= left_behind
    assert not left_behind & in_template, sorted(left_behind & in_template)
    assert not {"pest", "tpl"} & {Path(name).parts[0] for name in in_template if "/" in name}
    # Every input reached it, and the model the template runs reads only what it holds.
    assert inputs <= in_template
    assert {str(path) for path in _simulation_input_files(template)} <= in_template


@pytest.mark.slow
def test_a_pestpp_ies_run_of_that_template_reproduces_the_base_observations(staged):
    """PESTPP-IES's prior (NOPTMAX=-1) runs the model in the template: its BASE
    realization (every multiplier 1) reads, in every family, what the build read from
    the run's own output files."""

    _model, cal = staged
    template = cal.template_workspace
    expected = {}
    for ins in sorted(template.glob("*.csv.ins")):
        table = ins.name[: -len(".ins")]
        read = pyemu.pst_utils.InstructionFile(str(ins)).read_output_file(str(template / table))
        expected.update(read["obsval"].astype(float).to_dict())
    families = {name.split("_otype")[0] for name in expected}
    assert families == {"oname:hds", "oname:sfr_flow", "oname:drn_flow", "oname:lak_stage"}

    results = cal.prior(reals=3)
    base = pd.read_csv(Path(results.workspace) / "inputs.0.obs.csv", index_col=0)
    base = base.loc[[label for label in base.index if str(label) == "base"]].iloc[0]
    names = sorted(expected)
    # PESTPP-IES writes its ensemble files to six significant digits.
    assert np.allclose(base[names].to_numpy(dtype=float),
                       [expected[name] for name in names], rtol=1e-5, atol=0.0)

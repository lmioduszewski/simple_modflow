"""Groups and diffs keep both runs of one model and take each run's recorded
CRS (ledgers 373, 374).

Every model here is a real MF6 run on a 2-layer, 3 x 4 structured (DIS) grid,
written and run by FloPy in ``tmp_path`` and opened with ``load_mf6_run`` --
the route a user's two scenario runs take.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from types import SimpleNamespace

import flopy
import pytest

import myflopy as mf
from myflopy.project.group._shared import _diff_member_keys

NROW, NCOL = 3, 4


def _dis_run(workspace: Path, *, name: str = "m", east_head: float = 1.0) -> Path:
    """Write and run a DIS model whose head falls west (5.0) to east (``east_head``)."""

    workspace.mkdir(parents=True, exist_ok=True)
    sim = flopy.mf6.MFSimulation(sim_name="sim", sim_ws=str(workspace), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=2, nrow=NROW, ncol=NCOL, delr=10.0, delc=10.0, top=10.0, botm=[0.0, -10.0]
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    records = []
    for layer in range(2):
        for row in range(NROW):
            records.append(((layer, row, 0), 5.0))
            records.append(((layer, row, NCOL - 1), east_head))
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=records)
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        budget_filerecord=f"{name}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    sim.write_simulation(silent=True)
    success, report = sim.run_simulation(silent=True)
    assert success, "\n".join(report[-20:])
    return workspace


@pytest.fixture(scope="module")
def two_runs(tmp_path_factory):
    """Two runs of ONE model (both named ``m``) whose east boundary differs."""

    root = tmp_path_factory.mktemp("two_runs")
    return _dis_run(root / "base"), _dis_run(root / "wet", east_head=3.0)


# --- 373: two runs of one model -----------------------------------------


def test_a_sequence_of_two_same_named_models_is_refused_naming_the_mapping(two_runs):
    base, wet = (mf.load_mf6_run(ws) for ws in two_runs)
    assert base.name == wet.name == "m"

    with pytest.raises(ValueError) as excinfo:
        mf.ModelGroup([base, wet])
    message = str(excinfo.value)
    assert "'m'" in message
    assert "ModelGroup({'base': m1, 'wet': m2}" in message

    # The mapping spelling it names keeps both.
    group = mf.ModelGroup({"base": base, "wet": wet})
    assert list(group.models) == ["base", "wet"]


def test_model_diff_of_two_runs_of_one_model_compares_them(two_runs):
    base, wet = (mf.load_mf6_run(ws) for ws in two_runs)

    diff = base.diff(wet)

    # Keyed by the folders that tell the runs apart, the reference first.
    assert diff.reference == "base"
    assert diff.model_names == ["wet"]
    heads = diff.hds.summary()
    assert list(heads["model"]) == ["wet"]
    # The east boundary moved 1.0 -> 3.0: the heads differ by exactly that there.
    assert heads["max_abs_diff"].iloc[0] == pytest.approx(2.0)
    report = diff.report()
    assert "reference: `base`" in report and "`wet`" in report


def test_model_diff_takes_a_workspace_path_of_the_same_model(two_runs):
    base_ws, wet_ws = two_runs
    diff = mf.load_mf6_run(base_ws).diff(str(wet_ws))
    assert (diff.reference, diff.model_names) == ("base", ["wet"])


def test_diff_keys_prefer_names_then_the_distinguishing_folder_then_a_fixed_pair(tmp_path):
    def model(name, workspace):
        return SimpleNamespace(name=name, workspace=Path(workspace))

    # Distinct names: unchanged.
    assert _diff_member_keys([model("a", tmp_path / "x"), model("b", tmp_path / "x")]) == ["a", "b"]
    # A loaded run's own folder.
    assert _diff_member_keys(
        [model("m", tmp_path / "base"), model("m", tmp_path / "wet")]
    ) == ["base", "wet"]
    # A built model lives in <folder>/<model name>: the parent tells them apart.
    assert _diff_member_keys(
        [model("m", tmp_path / "reference" / "m"), model("m", tmp_path / "variant" / "m")]
    ) == ["reference", "variant"]
    # One workspace twice.
    same = tmp_path / "base"
    assert _diff_member_keys([model("m", same), model("m", same)]) == ["reference", "other"]
    assert _diff_member_keys([model("m", same)] * 3) == ["reference", "other_1", "other_2"]


# --- 374: each run's recorded CRS ----------------------------------------


@pytest.fixture
def recorded_run(two_runs, tmp_path):
    workspace = tmp_path / "utm"
    shutil.copytree(two_runs[0], workspace)
    (workspace / "run.json").write_text(json.dumps({"crs": {"m": "EPSG:32610"}}), encoding="utf-8")
    return workspace


def test_a_group_of_workspace_paths_takes_each_runs_recorded_crs(recorded_run):
    assert mf.load_mf6_run(recorded_run)._crs == "EPSG:32610"
    assert mf.ModelGroup({"utm": recorded_run}).models["utm"]._crs == "EPSG:32610"
    assert mf.ModelGroup([str(recorded_run)]).models["m"]._crs == "EPSG:32610"
    # An explicit crs= still wins.
    explicit = mf.ModelGroup({"utm": recorded_run}, crs="EPSG:2927")
    assert explicit.models["utm"]._crs == "EPSG:2927"


def test_model_diff_takes_the_other_runs_recorded_crs(two_runs, recorded_run):
    diff = mf.load_mf6_run(two_runs[0]).diff(recorded_run)
    assert diff.model_names == ["utm"]
    assert diff.group.models["utm"]._crs == "EPSG:32610"
    forced = mf.load_mf6_run(two_runs[0]).diff(recorded_run, crs="EPSG:2927")
    assert forced.group.models["utm"]._crs == "EPSG:2927"

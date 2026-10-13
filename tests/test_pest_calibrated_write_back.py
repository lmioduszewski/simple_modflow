"""A calibrated realization comes back as library packages (ledger 560).

``IesResults.calibrated(realization="base")`` replays that realization's own forward
run -- its parameter values filled into every template at full precision, then
``forward_run.py`` cut before the model command -- and reads back what its MODFLOW 6
read. ``CalibratedInputs.to_library`` substitutes those values into the packages the
RUN built and forks the run's simulation onto them.

The model is a package-first ``Project``: two layers on a 6x6 irregular Voronoi grid,
``npf`` and an array RCH as library packages (``mf.ref``), a GHB declared from a
GeoPackage against the model's context (eager, so pinned), two RIVs (``riv`` and a
boundnamed ``crisp``), a CHD and a
steady + transient period pair. The calibration: pilot-point ``k`` on layer 0 with its
field captured, a constant ``k33``, ``recharge`` in two period groups, ``ghb.cond`` and
one ``riv.stage`` offset shared by both RIVs (which carries ``rbot``, ledger 519) --
against four head
targets, one PESTPP-IES iteration of six realizations on two workers.
"""

from __future__ import annotations

import json
import shutil
import warnings
from pathlib import Path
from types import SimpleNamespace

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import LineString

import myflopy as mf
from myflopy.grid_identity import GRID_PIN, MOVE_RECIPE, SOURCE_SHA, package_move
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.pest.calibrated import CalibratedInputs, replay_script
from myflopy.modflow.mf6.pest.ies import IesResults

pyemu = pytest.importorskip("pyemu")
flopy = pytest.importorskip("flopy")

NLAY = 2
PERIODS = [(1.0, 1, 1.0), (30.0, 2, 1.0)]


def _grid():
    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    return grid, props, xy


def _project(root: Path):
    """The module docstring's project, run; ``(project, run, flopy-backed model)``."""

    grid, props, xy = _grid()
    ncpl = grid.ncpl
    xmin, ymin, xmax, ymax = grid.gdf_vorPolys.total_bounds
    west = [int(c) for c in np.flatnonzero(xy[:, 0] < xmin + 100)]
    middle = [int(c) for c in np.flatnonzero(np.abs(xy[:, 0] - (xmin + xmax) / 2) < 50)]
    quarter = [int(c) for c in np.flatnonzero(np.abs(xy[:, 0] - (xmin + 3 * xmax) / 4) < 40)]
    east_line = root / "east.gpkg"
    gpd.GeoDataFrame(
        {"head": [101.0], "cond": [30.0], "layer": [1], "name": ["east"]},
        geometry=[LineString([(xmax - 20.0, ymin + 1.0), (xmax - 20.0, ymax - 1.0)])],
        crs=grid.crs,
    ).to_file(east_line, layer="east")
    context = mf.ModelContext(grid=grid)
    project = mf.Project(root / "proj", name="p")
    k = np.vstack([np.linspace(5.0, 15.0, ncpl), np.full(ncpl, 2.0)])
    project.add_package("npf/base", mf.npf(k=k, k33=[1.0, 0.2], icelltype=0))
    project.add_package("rch/base", mf.rch.array(recharge={0: 1e-3, 1: 2e-3}, irch=0))
    model = mf.gwf("m", context=context, packages=[
        mf.disv(nlay=NLAY, ncpl=ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=110.0,
                botm=[50.0, 0.0]),
        mf.ref("npf/base"),
        mf.sto(ss=1e-5, sy=0.2, iconvert=0, steady_state={0: True}, transient={1: True}),
        mf.ic(strt=100.0),
        mf.chd(stress_period_data={0: [[(0, c), 100.0] for c in west]}),
        mf.ghb.gpkg(east_line, layer="east", context=context, nper=1, head="head",
                    conductance="cond"),
        mf.riv(stress_period_data={0: [[(0, c), 99.0, 40.0, 97.0, "main"] for c in middle]},
               boundnames=True),
        # A second RIV, as Ten Trails carries: named, with boundnames (an identity field).
        # Both carry boundnames: pyEMU refuses one parameter over list files of
        # different widths ("shape mismatch for array style").
        mf.riv(stress_period_data={0: [[(0, c), 100.5, 20.0, 99.5, "crisp"] for c in quarter]},
               name="crisp", boundnames=True),
        mf.ref("rch/base"),
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=len(PERIODS), perioddata=PERIODS, time_units="days")])
    project.add_simulation(simulation)
    run = project.prepare_run("r", simulation)
    run.execute()
    assert run.success, "\n".join(run.report[-20:])
    return project, run, run.model(), xy


@pytest.fixture(scope="module")
def calibration(tmp_path_factory):
    root = tmp_path_factory.mktemp("writeback")
    project, run, model, xy = _project(root)
    picks = [3, 10, 20, 27]
    targets = mf.HeadTargets(
        locations=gpd.GeoDataFrame(
            {"name": ["W1", "W2", "W3", "W4"], "layer": [0, 0, 1, 0]},
            geometry=gpd.points_from_xy(xy[picks, 0], xy[picks, 1])),
        values=pd.DataFrame({"per": [0, 1], "W1": [100.5, 100.8], "W2": [100.3, 100.6],
                             "W3": [100.2, 100.4], "W4": [100.9, 101.0]}),
        time_column="per",
    )
    cal = model.pest("wb", start_datetime="2024-01-01")
    cal.parameterize("k", style="pilotpoints", layers=[0], pp_spacing=250.0,
                     bounds=(0.2, 5.0), capture=True)
    cal.parameterize("k33", style="constant", bounds=(0.2, 5.0))
    cal.parameterize("recharge", style="constant", periods={"a": [0], "b": [1]},
                     bounds=(0.5, 2.0))
    cal.parameterize("ghb.cond", style="constant", bounds=(0.2, 5.0))
    cal.parameterize("riv.stage", style="constant", bounds=(-1.0, 1.0),
                     package=("riv", "crisp"))
    cal.observe(targets)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build("wb.pst", noptmax=0)
        ies = cal.run_ies(reals=6, iterations=1, workers=2, ies_lambda_mults=1.0,
                          lambda_scale_fac=1.0)
    return SimpleNamespace(project=project, run=run, model=model, cal=cal, ies=ies,
                           root=root)


@pytest.fixture(scope="module")
def base(calibration):
    return calibration.ies.calibrated()


def _base_par(ies: IesResults) -> pd.Series:
    path = ies.workspace / f"{ies.case}.{ies.posterior_iteration}.base.par"
    values = pyemu.pst_utils.read_parfile(str(path))["parval1"].astype(float)
    values.index = values.index.str.lower()
    return values


def _value_of(ies: IesResults, prefix: str) -> float:
    values = _base_par(ies)
    return float(values[[name for name in values.index if name.startswith(prefix)]].iloc[0])


def _copy_master(ies: IesResults, target: Path) -> Path:
    shutil.copytree(ies.workspace, target)
    return target


# -- reconstruction -------------------------------------------------------------------


def _metadata_records_where_every_parameter_lands(calibration):
    metadata = json.loads(
        (calibration.ies.workspace / "myflopy_pest_metadata.json").read_text())
    assert metadata["version"] == 2
    records = {record["name"]: record for record in metadata["parameters"]}
    assert set(records) == {"k", "k33", "recharge", "ghbcond", "rivstage_riv_crisp"}
    assert all(record["model"] == "m" for record in records.values())
    assert records["k"]["variable"] == "k" and records["k"]["fields"] == []
    assert records["recharge"]["variable"] == "recharge"
    assert records["ghbcond"]["fields"] == ["cond"]
    assert records["rivstage_riv_crisp"]["fields"] == ["stage", "rbot"]
    assert records["rivstage_riv_crisp"]["packages"] == ["crisp", "riv"]


def _replayed_k_matches_the_realizations_captured_field(calibration, base):
    """The base row's captured k (6 significant digits in obs.csv) is the replayed field."""

    observations = calibration.ies.posterior._df.loc["base"]
    captured = observations[observations.index.str.contains("kfieldl0")]
    cells = [int(name.split("_i:")[1].split("_")[0]) for name in captured.index]
    k = base.get("npf", "k")
    assert k.shape == (NLAY, len(cells))
    replayed = k[0, cells]
    assert np.allclose(replayed, captured.to_numpy(dtype=float), rtol=1e-5, atol=0)
    # The pilot points moved the field: not the run's own linspace any more.
    assert not np.allclose(replayed, np.linspace(5.0, 15.0, len(cells))[cells])


def _base_reads_the_full_precision_par_file(calibration, base):
    """k33 is org x the base.par multiplier to the file's 9 digits -- not par.csv's 6.

    Reading ``par.csv`` (6 significant digits) or filling with pyEMU's
    ``write_input_files`` (4) moves the multiplier by ~1e-6 or more relative.
    """

    multiplier = _value_of(calibration.ies, "pname:k33")
    k33 = base.get("npf", "k33")
    # pyEMU writes the file at 9 significant digits (%15.8E): half a unit is 5e-9.
    assert np.allclose(k33[0], 1.0 * multiplier, rtol=6e-9, atol=0)
    assert np.allclose(k33[1], 0.2 * multiplier, rtol=6e-9, atol=0)
    assert base.parameter_file.name.endswith(".base.par")


def _period_groups_and_list_targets_come_back_with_their_own_values(calibration, base):
    ies = calibration.ies
    recharge = base.get("rch", "recharge")
    assert sorted(recharge) == [0, 1]
    assert np.allclose(recharge[0], 1e-3 * _value_of(ies, "pname:recharge_a"), rtol=6e-9)
    assert np.allclose(recharge[1], 2e-3 * _value_of(ies, "pname:recharge_b"), rtol=6e-9)
    ghb = base.get("ghb", "stress_period_data")
    assert np.allclose(ghb["cond"], 30.0 * _value_of(ies, "pname:ghbcond"), rtol=6e-9)
    assert np.allclose(ghb["bhead"], 101.0)
    riv = base.get("riv", "stress_period_data")
    offset = _value_of(ies, "pname:rivstage")
    assert np.allclose(riv["stage"], 99.0 + offset, rtol=0, atol=1e-6)
    # The forward run carries rbot with the stage (ledger 519), and so does this.
    assert np.allclose(riv["rbot"], 97.0 + offset, rtol=0, atol=1e-6)
    crisp = base.get("crisp", "stress_period_data")
    assert np.allclose(crisp["stage"], 100.5 + offset, rtol=0, atol=1e-6)
    assert np.allclose(crisp["rbot"], 99.5 + offset, rtol=0, atol=1e-6)
    assert set(crisp["boundname"]) == {"crisp"}
    # ... and the package written back carries both: stage AND rbot moved, cond kept.
    records = base.packages()["riv"].options["stress_period_data"][0]
    assert np.allclose([r[1] for r in records], 99.0 + offset, rtol=0, atol=1e-6)
    assert np.allclose([r[3] for r in records], 97.0 + offset, rtol=0, atol=1e-6)
    assert np.allclose([r[2] for r in records], 40.0)
    crisp_records = base.packages()["crisp"].options["stress_period_data"][0]
    assert [r[4] for r in crisp_records] == ["crisp"] * len(crisp_records)
    assert np.allclose([r[3] for r in crisp_records], 99.5 + offset, rtol=0, atol=1e-6)
    summary = base.summary().set_index("package")
    assert summary.loc["riv", "change"] == "offset"
    assert summary.loc["riv", "min"] == pytest.approx(offset, abs=1e-6)
    assert summary.loc["k33" == summary["variable"], "change"].tolist() == ["ratio"]


def _another_realization_is_read_from_the_ensemble(calibration):
    ies = calibration.ies
    other = ies.calibrated("1")
    assert other.realization == "1"
    frame = ies.par_ensemble(ies.posterior_iteration)._df
    row = frame.loc[[label for label in frame.index if str(label) == "1"][0]]
    multiplier = float(row[[n for n in row.index if n.startswith("pname:k33")]].iloc[0])
    assert np.allclose(other.get("npf", "k33")[0], multiplier, rtol=6e-9)


@pytest.mark.slow
def test_the_canonical_testing_model_comes_back_as_its_base_realization_read_it(tmp_path):
    """The canonical TESTING model (``build_canonical_calibration_demo``, 21 x 21, every
    package) under canonical_06's calibration: pilot-point ``k`` on layer 0, a list-RCH
    ``recharge`` constant. Not a Project run, so it answers ``get``/``summary`` and
    refuses ``packages`` by name."""

    from myflopy.modflow.mf6.canonical_calibration import build_canonical_calibration_demo
    from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        demo = build_canonical_calibration_demo(tmp_path / "model",
                                                config=CanonicalModelConfig.testing())
        cal = demo.model.pest("cw", start_datetime="2024-01-01")
        cal.parameterize("k", style="pilotpoints", pp_spacing=800.0, layers=[0],
                         bounds=(0.2, 5.0), capture=True)
        cal.parameterize("recharge", style="constant", bounds=(0.5, 2.0))
        cal.observe(demo.head_targets)
        cal.build("cw.pst", noptmax=0)
        ies = cal.run_ies(reals=4, iterations=1, workers=2, ies_lambda_mults=1.0,
                          lambda_scale_fac=1.0)
    calibrated = ies.calibrated()
    observations = ies.posterior._df.loc["base"]
    captured = observations[observations.index.str.contains("kfieldl0")]
    cells = [int(name.split("_i:")[1].split("_")[0]) for name in captured.index]
    k = calibrated.get("npf", "k")
    assert np.allclose(k[0, cells], captured.to_numpy(dtype=float), rtol=1e-5, atol=0)
    multiplier = _value_of(ies, "pname:recharge")
    summary = calibrated.summary().set_index("package")
    assert summary.loc["rch", "min"] == pytest.approx(multiplier, rel=1e-6)
    with pytest.raises(FileNotFoundError, match="not run through a Project"):
        calibrated.packages()


# -- the library ------------------------------------------------------------------------


def _write_back_reproduces_the_base_realizations_simulated_heads(calibration, base):
    """The round trip: a model rebuilt from to_library simulates the base row's heads."""

    project, ies = calibration.project, calibration.ies
    keys = base.to_library(project)
    assert keys == {"npf/base": "npf/base/wb", "rch/base": "rch/base/wb",
                    "ghb": "ghb/wb", "riv": "riv/wb", "crisp": "crisp/wb"}
    assert "s/wb" in project.simulations
    # A fork, never an edit: the run's packages are what they were.
    assert project.packages["npf/base"].options["k"][0, 0] == 5.0
    run = project.prepare_run("calibrated", project.simulations["s/wb"])
    run.execute()
    assert run.success
    template = calibration.cal.template_workspace
    map_csv = template / "hds_head_target_map.csv"
    from myflopy.modflow.mf6.pest.forward_run import _write_head_target_csv

    table = run.workspace / "hds_simulated_heads.csv"
    _write_head_target_csv(model_name="m", mapping_csv=str(map_csv), output_csv=str(table),
                           perlen=[p[0] for p in PERIODS], sim_ws=str(run.workspace),
                           head_file="m.hds")
    instruction = pyemu.pst_utils.InstructionFile(str(template / "hds_simulated_heads.csv.ins"))
    simulated = instruction.read_output_file(str(table))["obsval"]
    # The base realization's simulated values at full precision: PESTPP-IES's residual
    # file for the iteration (obs.csv holds 6 significant digits, 5e-4 ft at 100 ft).
    residuals = pyemu.pst_utils.read_resfile(
        str(ies.workspace / f"{ies.case}.{ies.posterior_iteration}.base.rei"))
    expected = residuals.loc[simulated.index, "modelled"].astype(float)
    assert np.allclose(simulated.to_numpy(dtype=float), expected.to_numpy(), atol=1e-6, rtol=0)
    # ... and it is not the uncalibrated run's.
    start = pyemu.pst_utils.read_resfile(str(ies.workspace / f"{ies.case}.0.base.rei"))
    assert not np.allclose(start.loc[simulated.index, "modelled"].to_numpy(dtype=float),
                           expected.to_numpy(), atol=1e-3, rtol=0)


def _provenance_rides_in_metadata_and_drift_follows_the_parameter_file(calibration, base):
    project = calibration.project
    packages = base.packages()
    npf = packages["npf/base"]
    stamp = npf.metadata["calibrated_from"]
    assert stamp["pest"] == "wb" and stamp["realization"] == "base"
    assert stamp["iteration"] == calibration.ies.posterior_iteration
    assert stamp["run"] == "r" and stamp["source_key"] == "npf/base"
    assert stamp["parameters"] == ["k", "k33"]
    assert any("apply_pilotpoints_to_array" in line for line in stamp["replayed"])
    assert npf.metadata["source_path"] == str(base.parameter_file)
    assert npf.metadata[SOURCE_SHA]
    # Metadata is not digested: two write-backs of one realization digest alike.
    from myflopy.specs import package_digest

    again = calibration.ies.calibrated().packages()["npf/base"]
    assert package_digest(again) == package_digest(npf)
    # drift()'s source column reads the parameter file, for a run that used the key.
    keys = base.to_library(project, simulation="s/drift", keys={"npf/base": "npf/drift",
                           "rch/base": "rch/drift", "ghb": "ghb/drift", "riv": "riv/drift",
                           "crisp": "crisp/drift"})
    assert keys["npf/base"] == "npf/drift"
    project.prepare_run("drift", project.simulations["s/drift"], build=False)

    def _source():
        frame = project.refresh_runs().drift()
        return frame.loc[(frame["run"] == "drift") & (frame["key"] == "npf/drift"),
                         "source"].iloc[0]

    assert _source() == "ok"
    saved = base.parameter_file.read_bytes()
    try:
        base.parameter_file.write_bytes(saved + b"\n")
        assert _source() == "edited"
    finally:
        base.parameter_file.write_bytes(saved)


def _values_only_substitution_keeps_the_pin_and_drops_the_move_recipe(calibration, base):
    """A calibrated GHB keeps its pin (cells checked) and reads stale on a move."""

    declared = calibration.run.spec.models[0].package("ghb")
    assert GRID_PIN in declared.metadata and MOVE_RECIPE in declared.metadata
    ghb = base.packages()["ghb"]
    assert ghb.metadata[GRID_PIN] == declared.metadata[GRID_PIN]
    assert MOVE_RECIPE not in ghb.metadata
    assert package_move(ghb)[0] == "stale"
    # The declared GeoPackage is kept beside the parameter file that now decides it.
    assert ghb.metadata["calibrated_from"]["declared_source_path"] == \
        declared.metadata["source_path"]
    records = ghb.options["stress_period_data"][0]
    assert [r[0] for r in records] == [r[0] for r in declared.options["stress_period_data"][0]]


def test_a_list_record_on_another_cell_is_refused():
    """A cellid that differs from the run's record means the template is not this run's."""

    from myflopy.modflow.mf6.pest.calibrated import _substituted_records

    # FloPy's own layout: an object cellid beside float fields.
    read = np.empty(2, dtype=[("cellid", "O"), ("stage", "<f8"), ("cond", "<f8"),
                              ("rbot", "<f8")]).view(np.recarray)
    read["cellid"][0], read["cellid"][1] = (0, 9), (1, 0)
    read["stage"], read["cond"], read["rbot"] = 99.5, 40.0, 97.5
    run = {0: [[(0, 9), 99.0, 40.0, 97.0], [(0, 15), 99.0, 40.0, 97.0]]}
    with pytest.raises(ValueError, match="record 1's cellid is"):
        _substituted_records(run, {0: read}, ("stage", "rbot"), "riv")
    run[0][1][0] = (1, 0)
    substituted = _substituted_records(run, {0: read}, ("stage", "rbot"), "riv")
    assert substituted[0][1] == [(1, 0), 99.5, 40.0, 97.5]


# -- refusals ---------------------------------------------------------------------------


def _an_unknown_realization_is_refused_listing_the_labels(calibration):
    with pytest.raises(KeyError, match="No realization in iteration 1's parameter ensemble 'nope'"):
        calibration.ies.calibrated("nope")
    with pytest.raises(ValueError, match="not one of this run's"):
        calibration.ies.calibrated(iteration=7)


def _a_template_built_before_the_records_is_refused(calibration, tmp_path):
    master = _copy_master(calibration.ies, tmp_path / "old")
    path = master / "myflopy_pest_metadata.json"
    metadata = json.loads(path.read_text())
    metadata["version"] = 1
    for record in metadata["parameters"]:
        record.pop("variable")
    path.write_text(json.dumps(metadata))
    with pytest.raises(ValueError, match=r"Rebuild it with cal.build\(\)"):
        IesResults(master).calibrated()


def _a_tied_parameter_is_refused(calibration, tmp_path):
    master = _copy_master(calibration.ies, tmp_path / "tied")
    pst = pyemu.Pst(str(master / "wb.pst"))
    names = list(pst.parameter_data.index)
    pst.parameter_data.loc[names[0], "partrans"] = "tied"
    pst.parameter_data.loc[names[0], "partied"] = names[1]
    pst.write(str(master / "wb.pst"), version=2)
    with pytest.raises(ValueError, match="TIED"):
        IesResults(master).calibrated()


def _a_replay_that_fails_carries_its_stderr(calibration, tmp_path):
    master = _copy_master(calibration.ies, tmp_path / "broken")
    (master / "mult2model_info.csv").unlink()
    with pytest.raises(RuntimeError, match="Replaying the forward run .* failed"):
        IesResults(master).calibrated()


def _a_run_without_a_frozen_snapshot_is_refused(calibration, tmp_path):
    run_copy = tmp_path / "r"
    shutil.copytree(calibration.run.workspace, run_copy)
    (run_copy / ".myflopy" / "frozen.npz").unlink()
    master = _copy_master(calibration.ies, tmp_path / "master")
    path = master / "myflopy_pest_metadata.json"
    metadata = json.loads(path.read_text())
    metadata["original_workspace"] = str(run_copy)
    path.write_text(json.dumps(metadata))
    with pytest.raises(ValueError, match="has no frozen snapshot"):
        IesResults(master).calibrated().packages()


def _to_library_refuses_a_collision_another_project_and_an_unknown_label(
        calibration, base, tmp_path):
    project = calibration.project
    base.to_library(project, simulation=False, keys={
        "npf/base": "npf/c1", "rch/base": "rch/c1", "ghb": "ghb/c1", "riv": "riv/c1",
        "crisp": "crisp/c1"})
    with pytest.raises(KeyError, match="already has"):
        base.to_library(project, simulation=False, keys={
            "npf/base": "npf/c1", "rch/base": "rch/c1", "ghb": "ghb/c1", "riv": "riv/c1",
        "crisp": "crisp/c1"})
    base.to_library(project, simulation=False, overwrite=True, keys={
        "npf/base": "npf/c1", "rch/base": "rch/c1", "ghb": "ghb/c1", "riv": "riv/c1",
        "crisp": "crisp/c1"})
    with pytest.raises(KeyError, match="did not touch"):
        base.to_library(project, keys={"chd": "chd/x"})
    with pytest.raises(ValueError, match="owns no run"):
        base.to_library(mf.Project(tmp_path / "other", name="o"))


# -- the replay's cut, and the fork rule ------------------------------------------------


def test_the_replay_cuts_main_before_the_model_command(tmp_path):
    script = tmp_path / "forward_run.py"
    script.write_text(
        "import pyemu\n"
        "def main():\n"
        "    prepare('a.csv')\n"
        "    pyemu.os_utils.run(r'mf6')\n"
        "    after()\n"
        "if __name__ == '__main__':\n"
        "    main()\n"
    )
    source, kept = replay_script(script)
    assert kept == ["prepare('a.csv')"]
    assert "os_utils" not in source.split("def main")[1].split("if __name__")[0]
    script.write_text("def main():\n    prepare()\n")
    with pytest.raises(RuntimeError, match="never calls pyemu.os_utils.run"):
        replay_script(script)


def test_calibrated_inputs_is_exported():
    assert CalibratedInputs.__module__ == "myflopy.modflow.mf6.pest.calibrated"


# -- the three tests that share the IES run ---------------------------------------------
# One module fixture runs PESTPP-IES; under xdist's worksteal every worker that takes a
# test of this module runs it again, so the checks are grouped into three tests.


@pytest.mark.slow
def test_the_base_realization_comes_back_as_its_mf6_read_it(calibration, base):
    _metadata_records_where_every_parameter_lands(calibration)
    _replayed_k_matches_the_realizations_captured_field(calibration, base)
    _base_reads_the_full_precision_par_file(calibration, base)
    _period_groups_and_list_targets_come_back_with_their_own_values(calibration, base)
    _another_realization_is_read_from_the_ensemble(calibration)


@pytest.mark.slow
def test_to_library_reproduces_the_base_realization_and_records_where_it_came_from(
        calibration, base):
    _write_back_reproduces_the_base_realizations_simulated_heads(calibration, base)
    _provenance_rides_in_metadata_and_drift_follows_the_parameter_file(calibration, base)
    _values_only_substitution_keeps_the_pin_and_drops_the_move_recipe(calibration, base)


@pytest.mark.slow
def test_every_refusal_names_what_to_do(calibration, base, tmp_path):
    _an_unknown_realization_is_refused_listing_the_labels(calibration)
    _a_template_built_before_the_records_is_refused(calibration, tmp_path)
    _a_tied_parameter_is_refused(calibration, tmp_path)
    _a_replay_that_fails_carries_its_stderr(calibration, tmp_path)
    _a_run_without_a_frozen_snapshot_is_refused(calibration, tmp_path)
    _to_library_refuses_a_collision_another_project_and_an_unknown_label(
        calibration, base, tmp_path)

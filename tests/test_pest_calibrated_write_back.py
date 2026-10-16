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
from myflopy.modflow.mf6.pest.calibrated import CalibratedInputs
from myflopy.modflow.mf6.pest.ies import IesResults
from myflopy.modflow.mf6.pest.replay import replay_script

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
        {"head": [101.0], "cond": [30.0], "layer": [1], "name": ["East_L1"]},
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
        mf.riv(stress_period_data={0: [[(0, c), 100.5, 20.0, 99.5, "Crisp_Creek"] for c in quarter]},
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
        cal.build("wb.pst")
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
    # FloPy reads a boundname lower-cased; the run declared it mixed-case (Ten Trails
    # carries `ghb_L3`), and the write-back keeps the run's spelling.
    assert set(crisp["boundname"]) == {"crisp_creek"}
    # ... and the package written back carries both: stage AND rbot moved, cond kept.
    records = base.packages()["riv"].options["stress_period_data"][0]
    assert np.allclose([r[1] for r in records], 99.0 + offset, rtol=0, atol=1e-6)
    assert np.allclose([r[3] for r in records], 97.0 + offset, rtol=0, atol=1e-6)
    assert np.allclose([r[2] for r in records], 40.0)
    crisp_records = base.packages()["crisp"].options["stress_period_data"][0]
    assert [r[4] for r in crisp_records] == ["Crisp_Creek"] * len(crisp_records)
    assert np.allclose([r[3] for r in crisp_records], 99.5 + offset, rtol=0, atol=1e-6)
    summary = base.summary().set_index("package")
    assert summary.loc["riv", "change"] == "offset"
    assert summary.loc["riv", "min"] == pytest.approx(offset, abs=1e-6)
    assert summary.loc["k33" == summary["variable"], "change"].tolist() == ["ratio"]


def _calibrated_inputs_parameters_are_the_replays(calibration, base):
    """``CalibratedInputs.parameters``: every parameter, in the control file's order, as replayed."""

    pst = calibration.cal.pst
    parameters = base.parameters
    assert list(parameters.index) == [str(name).lower() for name in pst.par_names]
    expected = _base_par(calibration.ies).reindex(parameters.index)
    np.testing.assert_array_equal(parameters.to_numpy(), expected.to_numpy())
    parameters.iloc[0] = -1.0  # a copy: the realization's values do not move
    assert base.parameters.iloc[0] != -1.0


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
        cal.build("cw.pst")
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
    # The rebuilt run's outputs, read by the calibration's own readers (ledger 605).
    check = calibration.cal.run_template(calibration.root / "check", values=base,
                                         outputs_from=run)
    simulated = check["simulated"]
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


def test_a_boundname_differing_only_in_case_is_the_same_record():
    """FloPy reads a boundname lower-cased; the run's mixed-case spelling is kept.

    Ten Trails' GHB carries ``ghb_L3``: compared verbatim, every such package was refused
    as "not this run's" (ledger 562).
    """

    from myflopy.modflow.mf6.pest.calibrated import _substituted_records

    read = np.empty(1, dtype=[("cellid", "O"), ("bhead", "<f8"), ("cond", "<f8"),
                              ("boundname", "O")]).view(np.recarray)
    read["cellid"][0], read["bhead"], read["cond"], read["boundname"][0] = \
        (2, 11902), 555.0, 2.5, "ghb_l3"
    run = {0: [[(2, 11902), 555.0, 1.0, "ghb_L3"]]}
    substituted = _substituted_records(run, {0: read}, ("cond",), "ghb")
    assert substituted[0] == [[(2, 11902), 555.0, 2.5, "ghb_L3"]]
    read["boundname"][0] = "ghb_l4"
    with pytest.raises(ValueError, match="record 0's boundname is 'ghb_l3'"):
        _substituted_records(run, {0: read}, ("cond",), "ghb")


def test_a_number_flopy_reads_as_text_is_a_number_and_a_keyword_is_kept():
    """FloPy reads LAK's ``bedleak`` (it may say ``NONE``) and SFR's ``man`` (it may name a
    time series) as TEXT. A spelled number is a value -- written back, or for an unwritten
    column not taken for identity -- and a ``NONE`` the forward run left on a lake it did
    not select stays the run's (ledger 563). Before it, every SFR write-back was refused
    (``man`` ``0.03`` against ``'0.03'``) and a ``NONE`` raised inside ``float()``.
    """

    from myflopy.modflow.mf6.pest.calibrated import (
        CalibratedInputs,
        _Dataset,
        _records_frame,
        _substituted_records,
    )

    # FloPy's own LAK layout: claktype and bedleak are object columns.
    read = np.empty(3, dtype=[("ifno", "<i8"), ("iconn", "<i8"), ("cellid", "O"),
                              ("claktype", "O"), ("bedleak", "O"),
                              ("belev", "<f8")]).view(np.recarray)
    read["ifno"], read["iconn"], read["belev"] = [0, 1, 1], [0, 0, 1], 0.0
    for index, (cell, bedleak) in enumerate([(6, "1.0"), (13, "0.3"), (14, "none")]):
        read["cellid"][index], read["claktype"][index] = (0, cell), "vertical"
        read["bedleak"][index] = bedleak
    run = [[0, 0, (0, 6), "VERTICAL", 0.1, 0.0], [1, 0, (0, 13), "vertical", 0.3, 0.0],
           [1, 1, (0, 14), "vertical", "NONE", 0.0]]
    assert _substituted_records(run, read, ("bedleak",), "lak") == [
        [0, 0, (0, 6), "VERTICAL", 1.0, 0.0], [1, 0, (0, 13), "vertical", 0.3, 0.0],
        [1, 1, (0, 14), "vertical", "NONE", 0.0]]
    assert _records_frame(read, ("bedleak",))["bedleak"].tolist() == [1.0, 0.3, "none"]
    dataset = _Dataset(model="m", package="lak", variable="connectiondata", family="list",
                       fields=("bedleak",), parameters=["lakbedleak"], additive=False)
    flat = CalibratedInputs._flat(read, dataset)
    assert flat[:2].tolist() == [1.0, 0.3] and np.isnan(flat[2])
    # A keyword where the run declared a number: the files are not the run's.
    run[2][4] = 0.35
    with pytest.raises(ValueError, match=r"record 2's bedleak is 0.35 in the run and "
                                         r"'none' .*not a number"):
        _substituted_records(run, read, ("bedleak",), "lak")

    # FloPy's SFR layout: man is an object column the calibration does not write.
    read = np.empty(1, dtype=[("ifno", "<i8"), ("cellid", "O"), ("rhk", "<f8"),
                              ("man", "O"), ("boundname", "O")]).view(np.recarray)
    read["ifno"], read["rhk"] = 0, 1.5
    read["cellid"][0], read["man"][0], read["boundname"][0] = (0, 10), "0.03", "creek_b"
    run = [[0, (0, 10), 0.5, 0.03, "Creek_B"]]
    assert _substituted_records(run, read, ("rhk",), "sfr") == [
        [0, (0, 10), 1.5, 0.03, "Creek_B"]]
    # A time-series NAME there is identity, compared as text.
    run[0][3] = "tsman"
    with pytest.raises(ValueError, match="record 0's man is 'tsman'"):
        _substituted_records(run, read, ("rhk",), "sfr")


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


def _a_copy_of_the_run(calibration, tmp_path, name, *, fingerprinted=True):
    """A copy of the run in the project's ``runs/`` and a master pointing at it."""

    run_copy = calibration.project.runs_dir / name
    shutil.copytree(calibration.run.workspace, run_copy)
    if not fingerprinted:
        manifest = run_copy / "run.json"
        data = json.loads(manifest.read_text())
        data["metadata"].pop("inputs")
        manifest.write_text(json.dumps(data))
    master = _copy_master(calibration.ies, tmp_path / f"master_{name}")
    path = master / "myflopy_pest_metadata.json"
    metadata = json.loads(path.read_text())
    metadata["original_workspace"] = str(run_copy)
    path.write_text(json.dumps(metadata))
    return run_copy, master


def _edit(path: Path, old: str, new: str) -> None:
    text = path.read_text()
    assert old in text, (path, old)
    path.write_text(text.replace(old, new, 1))


_KEYS = {"npf/base": "npf/d", "rch/base": "rch/d", "ghb": "ghb/d", "riv": "riv/d",
         "crisp": "crisp/d"}


def _files_that_depart_from_the_record_are_refused(calibration, tmp_path):
    """PEST calibrated the run's FILES; the write-back substitutes into its RECORD (ledger
    611). A substituted package whose file departs is refused by packages(); any
    departure refuses the simulation fork, naming simulation=False -- ledger 562's shape."""

    project = calibration.project
    made = []
    try:
        # A file of a package the calibration substitutes (npf/base: k and k33).
        run_copy, master = _a_copy_of_the_run(calibration, tmp_path, "r_npf")
        made.append(run_copy)
        _edit(run_copy / "m.npf", "0.20000000", "0.30000000")
        with pytest.raises(ValueError, match=r"package\(s\) this calibration substitutes "
                                             r"\(m\.npf: changed outside myflopy"):
            IesResults(master).calibrated().packages()

        # A file of a package it does NOT substitute: the packages come back, the fork is
        # refused naming simulation=False, and nothing is registered before it.
        run_copy, master = _a_copy_of_the_run(calibration, tmp_path, "r_tdis")
        made.append(run_copy)
        _edit(run_copy / "s.tdis", "30.00000000", "31.00000000")
        inputs = IesResults(master).calibrated()
        assert set(inputs.packages()) == set(_KEYS)
        with pytest.raises(ValueError, match=r"s\.tdis: changed outside myflopy.*"
                                             r"Pass simulation=False"):
            inputs.to_library(project, keys=_KEYS, simulation="s/d")
        assert not set(_KEYS.values()) & set(project.packages)
        assert "s/d" not in project.simulations
        inputs.to_library(project, keys=_KEYS, simulation=False)
        assert set(_KEYS.values()) <= set(project.packages)
        for key in _KEYS.values():
            del project.packages[key]

        # A run written before runs recorded their inputs: the record is built and
        # written beside its files, and the same narrow rule applies to what that finds.
        run_copy, master = _a_copy_of_the_run(calibration, tmp_path, "r_old",
                                              fingerprinted=False)
        made.append(run_copy)
        assert set(IesResults(master).calibrated().packages()) == set(_KEYS)
        _edit(run_copy / "m.npf", "0.20000000", "0.30000000")
        with pytest.raises(ValueError, match=r"m\.npf: not what the record writes -- "
                                             r"GRIDDATA: 1 of"):
            IesResults(master).calibrated().packages()
    finally:
        for path in made:
            shutil.rmtree(path, ignore_errors=True)
        project.refresh_runs()


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


def _to_library_refuses_a_simulation_the_run_did_not_run(calibration, base):
    """The fork keeps the LIBRARY's other packages, so one edited since the run is refused
    -- before any package is registered (ledger 562)."""

    project = calibration.project
    original = project.simulations["s"]
    moved = original.with_packages(
        mf.chd(stress_period_data={0: [[(0, 0), 90.0]]}), slot="chd", on="m")
    keys = {"npf/base": "npf/moved", "rch/base": "rch/moved", "ghb": "ghb/moved",
            "riv": "riv/moved", "crisp": "crisp/moved"}
    project.simulations["s"] = moved
    try:
        with pytest.raises(ValueError, match="inline package 'chd' changed"):
            base.to_library(project, keys=keys, simulation="s/moved")
        assert not set(keys.values()) & set(project.packages)
        assert "s/moved" not in project.simulations
        # The packages alone are still registered on request.
        base.to_library(project, keys=keys, simulation=False)
        assert set(keys.values()) <= set(project.packages)
    finally:
        project.simulations["s"] = original


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
    _calibrated_inputs_parameters_are_the_replays(calibration, base)
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
    _files_that_depart_from_the_record_are_refused(calibration, tmp_path)
    _to_library_refuses_a_collision_another_project_and_an_unknown_label(
        calibration, base, tmp_path)
    _to_library_refuses_a_simulation_the_run_did_not_run(calibration, base)

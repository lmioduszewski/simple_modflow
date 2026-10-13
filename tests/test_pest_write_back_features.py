"""Lake-bed leakance and streambed K come back through the write-back, beside a screened well.

Integration of ledger 560 (the write-back) with ledger 558 (``lak.bedleak``/``sfr.rhk``,
``style="feature"``) and ledger 556 (screened head targets). The model is
``tests/pest_feature_models.py``'s, declared package-first on a ``Project`` so the
write-back has a run to substitute into: DISV 2 x 4 x 5, LAK ``lak`` (lake 0 ``Keevie``,
three connections, two of them to cell 7; lake 1 ``marjorie``) and SFR ``sfr`` (reaches
0-2 ``creek_a``, 3-4 ``Creek_B``), both raw (``mf.lak.flopy``/``mf.sfr.flopy``) so their
records are the run's own. The calibration: one leakance multiplier per lake, one rhk
multiplier per named stream and a constant K on layer 1, against an unscreened well and
a well screened across both layers -- whose T-weighted head the forward run computes
AFTER the model command, from the realization's K files, so the replay (cut before that
command) never runs it. One PESTPP-IES iteration of four realizations on two workers.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from tests.pest_feature_models import LAKE_CONNECTIONS, NCOL, NLAY, NROW, REACHES

pyemu = pytest.importorskip("pyemu")
flopy = pytest.importorskip("flopy")
from flopy.utils.gridutil import get_disv_kwargs  # noqa: E402

K = (10.0, 1.0)
LAKES = ["Keevie", "marjorie"]
#: Two head targets: one in layer 0 between the lakes and the streams, one screened
#: across both layers (top 100, bottoms 50/0) beside the second stream.
SCREENED = {"name": "S", "cell": 12, "screen_top": 80.0, "screen_bottom": 20.0}
PLAIN = {"name": "W", "cell": 8, "layer": 0}


def _project(root: Path):
    """The module docstring's project, run: ``(project, run)``."""

    grid = get_disv_kwargs(NLAY, NROW, NCOL, np.full(NCOL, 10.0), np.full(NROW, 10.0),
                           100.0, [50.0, 0.0])
    connections = []
    for lake, number, cell, kind, bedleak in LAKE_CONNECTIONS:
        geometry = (0.0, 0.0, 0.0, 0.0) if kind == "vertical" else (90.0, 100.0, 5.0, 10.0)
        connections.append([lake, number, (0, cell), kind, bedleak, *geometry])
    lakes = [[lake, 101.0, sum(row[0] == lake for row in LAKE_CONNECTIONS), LAKES[lake]]
             for lake in range(2)]
    downstream = {0: [-1], 1: [0, -2], 2: [1], 3: [-4], 4: [3]}
    reaches = [
        [number, (0, cell), 10.0, 2.0, 0.001, 95.0 - number * 0.1, 1.0, rhk, man,
         len(downstream[number]), 1.0, 0, bound]
        for number, (cell, bound, rhk, man) in enumerate(REACHES)
    ]
    project = mf.Project(root / "proj", name="p")
    model = mf.gwf("m", save_flows=True, packages=[
        mf.disv(nlay=grid["nlay"], ncpl=grid["ncpl"], nvert=grid["nvert"],
                vertices=grid["vertices"], cell2d=grid["cell2d"], top=100.0,
                botm=[50.0, 0.0]),
        mf.npf(k=list(K), icelltype=0),
        mf.ic(strt=100.0),
        mf.chd(stress_period_data={0: [[(0, 0), 100.0], [(0, 4), 99.0]]}),
        mf.lak.flopy(nlakes=2, packagedata=lakes, connectiondata=connections,
                     perioddata={0: [[0, "status", "active"], [1, "status", "active"]]},
                     boundnames=True),
        mf.sfr.flopy(nreaches=len(reaches), packagedata=reaches,
                     connectiondata=[[number, *links] for number, links in downstream.items()],
                     perioddata={0: [[0, "inflow", 10.0], [3, "inflow", 5.0]]},
                     boundnames=True),
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=2, perioddata=[(1.0, 1, 1.0)] * 2, time_units="days"),
        mf.ims(linear_acceleration="BICGSTAB")])
    project.add_simulation(simulation)
    run = project.prepare_run("r", simulation)
    run.execute()
    assert run.success, "\n".join(run.report[-20:])
    return project, run


def _targets():
    rows = [PLAIN, {**SCREENED, "layer": np.nan}]
    locations = pd.DataFrame(rows)
    locations.loc[0, ["screen_top", "screen_bottom"]] = np.nan
    return mf.HeadTargets(locations=locations, values=pd.DataFrame(
        {"time": [0, 1], "W": [99.6, 99.6], "S": [99.5, 99.5]}))


@pytest.fixture(scope="module")
def calibration(tmp_path_factory):
    root = tmp_path_factory.mktemp("writeback_features")
    project, run = _project(root)
    model = run.model()
    cal = model.pest("feat", start_datetime="2024-01-01")
    cal.parameterize("bedleak", bounds=(0.1, 10.0))
    cal.parameterize("rhk", style="feature", bounds=(0.2, 5.0))
    cal.parameterize("k", style="constant", layers=[1], bounds=(0.2, 5.0))
    cal.observe(_targets())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build("feat.pst", noptmax=0)
        ies = cal.run_ies(reals=4, iterations=1, workers=2, ies_lambda_mults=1.0,
                          lambda_scale_fac=1.0)
    return SimpleNamespace(project=project, run=run, cal=cal, ies=ies, base=ies.calibrated())


def _base_par(ies) -> pd.Series:
    path = ies.workspace / f"{ies.case}.{ies.posterior_iteration}.base.par"
    values = pyemu.pst_utils.read_parfile(str(path))["parval1"].astype(float)
    values.index = values.index.str.lower()
    return values


def _group_value(calibration, group: str) -> float:
    """The base realization's value of the one parameter in PEST group ``group``."""

    data = calibration.cal.pst.parameter_data
    names = [str(name).lower() for name in data.index[data["pargp"] == group]]
    assert len(names) == 1, (group, names)
    return float(_base_par(calibration.ies)[names[0]])


def _per_lake(calibration) -> list[float]:
    return [_group_value(calibration, f"lakbedleak_{name.lower()}") for name in LAKES]


def _per_reach(calibration) -> list[float]:
    streams = {"creek_a": _group_value(calibration, "sfrrhk_creek_a"),
               "creek_b": _group_value(calibration, "sfrrhk_creek_b")}
    return [streams[bound.lower()] for _, bound, _, _ in REACHES]


def _the_metadata_records_the_feature_datasets(calibration):
    metadata = json.loads(
        (calibration.ies.workspace / "myflopy_pest_metadata.json").read_text())
    records = {record["name"]: record for record in metadata["parameters"]}
    assert (records["lakbedleak"]["variable"], records["lakbedleak"]["fields"]) == \
        ("connectiondata", ["bedleak"])
    assert (records["sfrrhk"]["variable"], records["sfrrhk"]["fields"]) == \
        ("packagedata", ["rhk"])
    assert records["lakbedleak"]["style"] == records["sfrrhk"]["style"] == "feature"


def _each_lake_and_stream_comes_back_with_its_own_multiplier(calibration):
    base = calibration.base
    lakes, reaches = _per_lake(calibration), _per_reach(calibration)
    # Each feature drew its own value: one multiplier per lake, one per named stream.
    assert len(set(lakes)) == 2 and len(set(reaches)) == 2
    lak = base.get("lak", "connectiondata")
    expected = [bedleak * lakes[lake] for lake, _, _, _, bedleak in LAKE_CONNECTIONS]
    assert np.allclose(lak["bedleak"], expected, rtol=1e-8, atol=0)
    # Only the written field moves: the connection geometry is the run's.
    assert np.allclose(lak["connlen"], [0.0, 5.0, 5.0, 0.0, 0.0])
    sfr = base.get("sfr", "packagedata")
    assert np.allclose(sfr["rhk"], [rhk * factor for (_, _, rhk, _), factor
                                    in zip(REACHES, reaches, strict=True)], rtol=1e-8, atol=0)
    # FloPy reads `man` as text (it may name a time series); the file's is the run's.
    assert np.allclose(sfr["man"].astype(float), [man for _, _, _, man in REACHES])
    summary = base.summary().set_index("package")
    assert summary.loc["lak", "variable"] == "connectiondata"
    assert summary.loc["lak", "change"] == "ratio"
    assert summary.loc["lak", "min"] == pytest.approx(min(lakes), rel=1e-8)
    assert summary.loc["sfr", "max"] == pytest.approx(max(reaches), rel=1e-8)


def _the_package_written_back_keeps_every_other_field_and_the_runs_spelling(calibration):
    packages = calibration.base.packages()
    lakes, reaches = _per_lake(calibration), _per_reach(calibration)
    connections = packages["lak"].options["connectiondata"]
    for record, (lake, number, cell, kind, bedleak) in zip(
            connections, LAKE_CONNECTIONS, strict=True):
        assert record[:4] == [lake, number, (0, cell), kind]
        assert record[4] == pytest.approx(bedleak * lakes[lake], rel=1e-8)
    for record, (cell, bound, rhk, man), factor in zip(
            packages["sfr"].options["packagedata"], REACHES, reaches, strict=True):
        assert record[1] == (0, cell) and record[8] == man
        # FloPy reads `Creek_B` lower-cased; the run's spelling is what is written back.
        assert record[12] == bound
        assert record[7] == pytest.approx(rhk * factor, rel=1e-8)
    assert packages["lak"].metadata["calibrated_from"]["parameters"] == ["lakbedleak"]


def _the_replay_never_reached_the_screened_writer(calibration):
    replayed = calibration.base.replayed
    assert replayed, "the replay ran no statement"
    assert not any("_write_head_target_csv" in line for line in replayed)
    forward_run = (calibration.cal.template_workspace / "forward_run.py").read_text()
    assert "k_sources=" in forward_run  # the screened writer is in the forward run


def _the_rebuilt_model_simulates_the_base_realizations_heads(calibration, tmp_path):
    project, ies = calibration.project, calibration.ies
    keys = calibration.base.to_library(project)
    assert keys == {"npf": "npf/feat", "lak": "lak/feat", "sfr": "sfr/feat"}
    run = project.prepare_run("calibrated", project.simulations["s/feat"])
    run.execute()
    assert run.success, "\n".join(run.report[-20:])
    # What the rebuilt run's MF6 read: the calibrated leakance, rhk and K.
    sim = flopy.mf6.MFSimulation.load(sim_ws=str(run.workspace), verbosity_level=0,
                                      load_only=["lak", "sfr", "npf"])
    gwf = sim.get_model("m")
    lakes, reaches = _per_lake(calibration), _per_reach(calibration)
    assert np.allclose(gwf.get_package("lak").connectiondata.get_data()["bedleak"].astype(float),
                       [b * lakes[lake] for lake, _, _, _, b in LAKE_CONNECTIONS], rtol=1e-8)
    assert np.allclose(gwf.get_package("sfr").packagedata.get_data()["rhk"],
                       [rhk * f for (_, _, rhk, _), f in zip(REACHES, reaches, strict=True)],
                       rtol=1e-8)
    k = np.asarray(gwf.npf.k.array, dtype=float).reshape(NLAY, -1)
    assert np.allclose(k[1], K[1] * _group_value(calibration, "k"), rtol=1e-8)
    # The screened well is weighted by the REBUILT model's K, read as the forward run reads
    # a realization's K files.
    k_sources = []
    for layer in range(NLAY):
        path = tmp_path / f"k_layer{layer}.txt"
        np.savetxt(path, k[layer])
        k_sources.append({"layer": layer, "file": str(path), "factor": 1.0})
    template = calibration.cal.template_workspace
    from myflopy.modflow.mf6.pest.forward_run import _write_head_target_csv

    table = tmp_path / "hds_simulated_heads.csv"
    _write_head_target_csv(model_name="m", mapping_csv=str(template / "hds_head_target_map.csv"),
                           output_csv=str(table), perlen=[1.0, 1.0], sim_ws=str(run.workspace),
                           head_file="m.hds", k_sources=k_sources)
    instruction = pyemu.pst_utils.InstructionFile(str(template / "hds_simulated_heads.csv.ins"))
    simulated = instruction.read_output_file(str(table))["obsval"]
    residuals = pyemu.pst_utils.read_resfile(
        str(ies.workspace / f"{ies.case}.{ies.posterior_iteration}.base.rei"))
    expected = residuals.loc[simulated.index, "modelled"].astype(float)
    assert sorted(simulated.index) == sorted(
        f"oname:hds_otype:lst_usecol:{well}_per:{per}" for well in "ws" for per in (0, 1))
    assert np.allclose(simulated.to_numpy(dtype=float), expected.to_numpy(), atol=1e-6, rtol=0)
    # S is ONE observation of both layers' heads, each weighted by K x its 30 ft of screen
    # (both layers confined) -- with the calibrated layer-1 K.
    heads = flopy.utils.HeadFile(str(run.workspace / "m.hds"))
    try:
        for per in (0, 1):
            h = np.asarray(heads.get_data(kstpkper=(0, per))).reshape(NLAY, -1)[:, SCREENED["cell"]]
            by_hand = float(np.dot(k[:, SCREENED["cell"]], h) / k[:, SCREENED["cell"]].sum())
            assert simulated[f"oname:hds_otype:lst_usecol:s_per:{per}"] == pytest.approx(
                by_hand, abs=1e-9)
            assert abs(h[0] - h[1]) > 1e-6  # the layers differ, so the weighting matters
    finally:
        heads.close()


@pytest.mark.slow
def test_lake_leakance_and_streambed_k_come_back_beside_a_screened_well(calibration, tmp_path):
    _the_metadata_records_the_feature_datasets(calibration)
    _each_lake_and_stream_comes_back_with_its_own_multiplier(calibration)
    _the_package_written_back_keeps_every_other_field_and_the_runs_spelling(calibration)
    _the_replay_never_reached_the_screened_writer(calibration)
    _the_rebuilt_model_simulates_the_base_realizations_heads(calibration, tmp_path)

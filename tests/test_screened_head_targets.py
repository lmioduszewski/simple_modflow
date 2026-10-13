"""A well screened across layers is ONE observation of its transmissivity-weighted head.

Ledger 556: a location row giving ``screen_top``/``screen_bottom`` reads
the heads of every ACTIVE layer its screen opens into, weighted by K x the
saturated length of the screen inside each layer (MF6's ICELLTYPE rule, the one
``hydraulics.transmissivity`` states), and the PEST forward run weights with the
REALIZATION's K files. Before it, a well in two layers was two full-weight
observations of the one measured series (ledger 234), and screen elevations were
dropped.

Ledger 557: ``compare()``/``simulated_heads()`` and the forward-run
writer read through ONE period-end reader, which reads the head file the OC
DECLARES (``compare()`` opened ``<model>.hds`` through a ``HeadsPlus`` built twice
per call).

The model is FloPy-written and run with MF6: DISV 3 x 5 x 5, top 100, bottoms
60/30/0, K 10/1/5 (LAYERED constants), ICELLTYPE 1/0/0, a CHD in layer 0 and a
pumping well in layer 2, so the three layers' heads differ at every cell.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6.hydraulics import screened_head
from myflopy.modflow.mf6.observations.conc import ConcTargets
from myflopy.modflow.mf6.pest.forward_run import _write_head_target_csv

pyemu = pytest.importorskip("pyemu")
flopy = pytest.importorskip("flopy")
from flopy.utils.gridutil import get_disv_kwargs  # noqa: E402

NLAY, NROW, NCOL = 3, 5, 5
NCPL = NROW * NCOL
TOP, BOTM = 100.0, (60.0, 30.0, 0.0)
K = (10.0, 1.0, 5.0)
WELL_CELL = 12  # the centre; cell centre (250, 250)


def _write_run(workspace: Path, *, head_name: str = "m.hds", nper: int = 2, chd: float = 90.0) -> Path:
    sim = flopy.mf6.MFSimulation("s", sim_ws=str(workspace), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=[(1.0, 1, 1.0)] * nper, time_units="days")
    flopy.mf6.ModflowIms(sim, complexity="SIMPLE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m", save_flows=True)
    flopy.mf6.ModflowGwfdisv(gwf, **get_disv_kwargs(
        NLAY, NROW, NCOL, np.full(NCOL, 100.0), np.full(NROW, 100.0), TOP, list(BOTM)))
    flopy.mf6.ModflowGwfnpf(gwf, k=list(K), icelltype=[1, 0, 0])
    flopy.mf6.ModflowGwfic(gwf, strt=90.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, c), chd] for c in (0, 4, 20, 24)])
    flopy.mf6.ModflowGwfwel(gwf, stress_period_data={
        per: [[(2, WELL_CELL), -400.0 * (per + 1)]] for per in range(nper)})
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord=head_name, budget_filerecord="m.cbc",
                           saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])
    sim.write_simulation(silent=True)
    assert sim.run_simulation(silent=True)[0]
    return workspace


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    return _write_run(tmp_path_factory.mktemp("screens") / "run")


def _heads(workspace: Path, head_name: str = "m.hds") -> dict[int, np.ndarray]:
    """The file's own heads, per period, ``(nlay, ncpl)``."""

    hds = flopy.utils.HeadFile(str(workspace / head_name))
    try:
        return {kper: np.asarray(hds.get_data(kstpkper=(0, kper))).reshape(NLAY, -1)
                for kper in range(len(hds.get_kstpkper()))}
    finally:
        hds.close()


def _by_hand(h: np.ndarray, k, top: float, bottom: float, cell: int = WELL_CELL) -> float:
    """The T-weighted head written out longhand: layer 0 convertible, 1-2 confined."""

    tops = (TOP, *BOTM[:-1])
    weights, values = [], []
    for layer in range(NLAY):
        upper, lower = min(top, tops[layer]), max(bottom, BOTM[layer])
        if upper <= lower:
            continue
        length = upper - lower
        if layer == 0:  # icelltype 1: only the saturated part
            length = min(max(min(h[layer, cell], upper) - lower, 0.0), length)
        weights.append(k[layer] * length)
        values.append(h[layer, cell])
    return float(np.dot(weights, values) / np.sum(weights))


def _screened(top=70.0, bottom=20.0, *, layer=None, name="W1", extra=None):
    """One screened well at the centre cell, measured 85.0 in both periods."""

    row = {"name": name, "cell": WELL_CELL, "screen_top": top, "screen_bottom": bottom}
    if layer is not None:
        row["layer"] = layer
    rows = [row, *(extra or [])]
    names = list(dict.fromkeys(r["name"] for r in rows))
    values = pd.DataFrame({"time": [0, 1], **{n: [85.0, 85.0] for n in names}})
    return mf.HeadTargets(locations=pd.DataFrame(rows), values=values)


# --- (a) the pure computation, by hand ------------------------------------------


def test_equal_open_lengths_weight_by_k():
    assert screened_head([100.0, 90.0], [10.0, 1.0], [80.0, 70.0], [85.0, 75.0],
                         icelltype=[0, 0]) == pytest.approx((10 * 100.0 + 90.0) / 11)


def test_a_convertible_layer_whose_head_is_below_the_screen_weighs_nothing():
    assert screened_head([65.0, 90.0], [10.0, 1.0], [80.0, 70.0], [85.0, 75.0],
                         icelltype=[1, 0]) == pytest.approx(90.0)


def test_a_confined_layer_uses_the_whole_open_length_whatever_the_head():
    assert screened_head([65.0, 90.0], [10.0, 1.0], [80.0, 70.0], [85.0, 75.0],
                         icelltype=[0, 0]) == pytest.approx((10 * 65.0 + 90.0) / 11)


def test_a_partly_saturated_convertible_layer_weighs_its_saturated_part():
    # layer 0 open 80-90, head 85: 5 ft saturated; layer 1 confined, open 10 ft
    assert screened_head([85.0, 60.0], [1.0, 1.0], [80.0, 50.0], [90.0, 60.0],
                         icelltype=[1, 0]) == pytest.approx((5 * 85.0 + 10 * 60.0) / 15)


def test_all_dry_reads_the_deepest_open_layer():
    """Under NEWTON a dry cell holds a finite head below its bottom (12,290 of 90,446
    on Ten Trails): the deepest one is still a head that moves with the parameters."""

    assert screened_head([79.0, 49.0], [10.0, 1.0], [80.0, 50.0], [85.0, 60.0],
                         icelltype=[1, 1]) == pytest.approx(49.0)
    assert screened_head([-1e30, -1e30], [10.0, 1.0], [80.0, 50.0], [85.0, 60.0],
                         icelltype=[1, 1]) == -1e30


def test_a_thickstrt_layer_uses_the_fixed_length():
    assert screened_head([10.0, 90.0], [1.0, 1.0], [80.0, 70.0], [90.0, 75.0],
                         icelltype=[-1, 0], fixed_thickness=[10.0, np.nan]) == pytest.approx(
        (10 * 10.0 + 5 * 90.0) / 15)


# --- (b) one observation, compared on MF6's own heads ---------------------------


def test_a_screened_well_is_one_label_on_one_row_per_open_layer(run):
    model = mf.load_mf6_run(run)
    labels = _screened().observation_labels(model)
    assert labels["obsname"].tolist() == ["W1", "W1", "W1"]
    assert labels["layer"].tolist() == [0, 1, 2]
    assert labels["open_top"].tolist() == [70.0, 60.0, 30.0]
    assert labels["open_bottom"].tolist() == [60.0, 30.0, 20.0]
    assert labels["icelltype"].tolist() == [1, 0, 0]


def test_compare_reads_the_t_weighted_head_mf6_wrote(run):
    model = mf.load_mf6_run(run)
    heads = _heads(run)
    frame = _screened().compare(model)
    assert len(frame) == 2  # ONE observation per period, not one per layer
    for per in (0, 1):
        expected = _by_hand(heads[per], K, 70.0, 20.0)
        assert frame.loc[frame["per"] == per, "sim_head"].item() == pytest.approx(expected, rel=1e-12)
        # ... and it is none of the layers' own heads
        assert not np.isclose(expected, heads[per][:, WELL_CELL]).any()
    simulated = _screened().simulated_heads(model)
    assert list(simulated.columns) == ["per", "W1"]
    assert simulated["W1"].to_numpy() == pytest.approx(frame["sim_head"].to_numpy())


def test_a_screen_above_the_water_table_in_layer_0_weighs_its_saturated_part(run):
    model = mf.load_mf6_run(run)
    head0 = _heads(run)[0][0, WELL_CELL]
    assert 60.0 < head0 < 99.0
    frame = _screened(99.0, 40.0).compare(model)
    assert frame.loc[frame["per"] == 0, "sim_head"].item() == pytest.approx(
        _by_hand(_heads(run)[0], K, 99.0, 40.0), rel=1e-12)


def test_a_given_layer_inside_the_screen_is_accepted(run):
    model = mf.load_mf6_run(run)
    frame = _screened(layer=1).compare(model)
    assert frame["layer"].tolist() == [1, 1]
    assert frame["sim_head"].iloc[0] == pytest.approx(_by_hand(_heads(run)[0], K, 70.0, 20.0))


def test_two_layer_rows_without_screens_stay_two_observations(run):
    """Ledger 234 holds for a set with only layer columns."""

    model = mf.load_mf6_run(run)
    targets = mf.HeadTargets(
        locations=pd.DataFrame({"name": ["W1", "W1"], "cell": [WELL_CELL] * 2, "layer": [0, 2]}),
        values=pd.DataFrame({"time": [0], "W1": [85.0]}),
    )
    assert targets.observation_labels(model)["obsname"].tolist() == ["W1_lay0", "W1_lay2"]
    assert len(targets.compare(model)) == 2


def test_a_screened_and_an_unscreened_well_in_one_set(run):
    model = mf.load_mf6_run(run)
    targets = _screened(extra=[{"name": "P2", "cell": 7, "layer": 1,
                                "screen_top": np.nan, "screen_bottom": np.nan}])
    frame = targets.compare(model).sort_values(["name", "per"])
    heads = _heads(run)
    assert frame.loc[frame["name"] == "P2", "sim_head"].tolist() == pytest.approx(
        [heads[0][1, 7], heads[1][1, 7]])
    assert frame.loc[frame["name"] == "P2", "layer"].tolist() == [1, 1]


def test_from_cells_and_from_records_carry_screens(run):
    model = mf.load_mf6_run(run)
    by_cells = mf.HeadTargets.from_cells(cells=[WELL_CELL], names=["W1"], values=[[85.0], [85.0]],
                                         screen_tops=[70.0], screen_bottoms=[20.0])
    records = mf.HeadTargets.from_records([
        {"name": "W1", "cell": WELL_CELL, "time": t, "head": 85.0,
         "screen_top": 70.0, "screen_bottom": 20.0} for t in (0, 1)])
    expected = _screened().compare(model)["sim_head"].to_numpy()
    assert by_cells.compare(model)["sim_head"].to_numpy() == pytest.approx(expected)
    assert records.compare(model)["sim_head"].to_numpy() == pytest.approx(expected)


def test_a_screened_point_with_no_layer_matches_its_cell(run):
    """A point location: ``match_to_model`` resolves it with no layer to test."""

    import geopandas as gpd

    model = mf.load_mf6_run(run)
    points = gpd.GeoDataFrame(
        {"name": ["W1"], "screen_top": [70.0], "screen_bottom": [20.0]},
        geometry=gpd.points_from_xy([250.0], [250.0]),
    )
    targets = mf.HeadTargets(locations=points, values=pd.DataFrame({"time": [0], "W1": [85.0]}))
    frame = targets.compare(model)
    assert frame["cell"].tolist() == [WELL_CELL]
    assert frame["layer"].isna().all()
    assert frame["sim_head"].item() == pytest.approx(_by_hand(_heads(run)[0], K, 70.0, 20.0))


def test_from_hydrographs_carries_screens(run, tmp_path):
    workspace = _write_run(tmp_path / "dated")
    sim = flopy.mf6.MFSimulation.load(sim_ws=str(workspace), verbosity_level=0)
    sim.tdis.start_date_time.set_data("2020-01-01")
    sim.tdis.write()
    model = mf.load_mf6_run(workspace)
    locations = pd.DataFrame({"name": ["W1"], "cell": [WELL_CELL],
                              "screen_top": [70.0], "screen_bottom": [20.0]})
    readings = pd.DataFrame({"date": ["2020-01-02", "2020-01-03"], "W1": [85.0, 86.0]})
    targets = mf.HeadTargets.from_hydrographs(locations, readings, model=model)
    assert targets.screened
    frame = targets.compare(model)
    assert frame["sim_head"].tolist() == pytest.approx(
        [_by_hand(_heads(workspace)[per], K, 70.0, 20.0) for per in (0, 1)])


# --- the one reader reads the head file the OC declares -------------------------


def test_compare_reads_the_oc_declared_head_file_not_a_stale_model_hds(tmp_path):
    """R1-obs left this open: ``HeadsPlus`` opened ``<model>.hds`` whatever the OC said."""

    workspace = _write_run(tmp_path / "renamed", head_name="renamed.hds")
    # A stale <model>.hds holding other numbers, as a reused run directory keeps.
    other = _write_run(tmp_path / "other", chd=80.0)
    shutil.copyfile(other / "m.hds", workspace / "m.hds")
    stale = _heads(workspace, "m.hds")
    heads = _heads(workspace, "renamed.hds")
    assert not np.isclose(stale[0][1, 7], heads[0][1, 7])
    model = mf.load_mf6_run(workspace)
    targets = mf.HeadTargets(
        locations=pd.DataFrame({"name": ["P"], "cell": [7], "layer": [1]}),
        values=pd.DataFrame({"time": [0, 1], "P": [85.0, 85.0]}),
    )
    assert targets.compare(model)["sim_head"].tolist() == pytest.approx(
        [heads[0][1, 7], heads[1][1, 7]])
    assert targets.simulated_heads(model)["P"].tolist() == pytest.approx(
        [heads[0][1, 7], heads[1][1, 7]])


# --- (d) refusals ---------------------------------------------------------------


@pytest.mark.parametrize("top, bottom, match", [
    (70.0, np.nan, "one end of a screen"),
    (np.nan, 20.0, "one end of a screen"),
    (20.0, 70.0, "top must be above its bottom"),
    (20.0, 20.0, "top must be above its bottom"),
    ("deep", 20.0, "finite elevations"),
])
def test_a_malformed_screen_is_refused(top, bottom, match):
    with pytest.raises(ValueError, match=match):
        _screened(top, bottom)


def test_one_screen_column_without_its_twin_is_refused():
    with pytest.raises(ValueError, match="not 'screen_bottom'"):
        mf.HeadTargets(locations=pd.DataFrame({"name": ["W"], "cell": [1], "screen_top": [5.0]}),
                       values=pd.DataFrame({"time": [0], "W": [1.0]}))


def test_a_screened_name_on_two_rows_is_refused():
    with pytest.raises(ValueError, match="SCREENED target is one row"):
        _screened(extra=[{"name": "W1", "cell": WELL_CELL, "screen_top": 15.0,
                          "screen_bottom": 5.0}])


def test_a_layer_the_screen_does_not_open_into_is_refused(run):
    model = mf.load_mf6_run(run)
    with pytest.raises(ValueError, match=r"gives layer 0 .* opens into layer\(s\) \[1, 2\]"):
        _screened(55.0, 20.0, layer=0).compare(model)


def test_a_screen_above_the_model_is_refused(run):
    model = mf.load_mf6_run(run)
    with pytest.raises(ValueError, match="entirely above the model top"):
        _screened(150.0, 120.0).compare(model)


def test_a_screen_reaching_past_the_model_is_clipped_with_a_warning(run):
    model = mf.load_mf6_run(run)
    with pytest.warns(UserWarning, match="5 below the model bottom"):
        labels = _screened(70.0, -5.0).observation_labels(model)
    assert labels["open_bottom"].tolist() == [60.0, 30.0, 0.0]


def test_conc_targets_refuse_screens():
    with pytest.raises(ValueError, match="screened concentration is not read"):
        ConcTargets(locations=pd.DataFrame({"name": ["C"], "cell": [1], "screen_top": [5.0],
                                            "screen_bottom": [1.0]}),
                    values=pd.DataFrame({"time": [0], "C": [1.0]}))


def test_to_flopy_obs_refuses_a_screened_set(run):
    with pytest.raises(ValueError, match="MF6's OBS has no multi-layer head"):
        _screened().to_flopy_obs(mf.load_mf6_run(run))


def test_the_map_refuses_two_unscreened_rows_under_one_label(run, tmp_path):
    mapping = pd.DataFrame({"obsname": ["W", "W"], "name": ["W", "W"], "layer": [0, 1],
                            "cell": [1, 1]})
    mapping.to_csv(tmp_path / "map.csv", index=False)
    with pytest.raises(ValueError, match="without a screen interval"):
        _write_head_target_csv("m", tmp_path / "map.csv", tmp_path / "out.csv", sim_ws=run)


# --- (c) PEST: the forward run weights with the REALIZATION's K -----------------


def _forward_run(cal, pst, worker: Path, values: dict) -> Path:
    shutil.rmtree(worker, ignore_errors=True)
    shutil.copytree(cal.template_workspace, worker)
    par = pst.parameter_data.copy()
    for key, value in values.items():
        par.loc[par["pargp"] == key, "parval1"] = value
    for tpl, inp in zip(pst.template_files, pst.input_files, strict=True):
        pyemu.pst_utils.write_to_template(par["parval1"], str(worker / tpl), str(worker / inp))
    done = subprocess.run([sys.executable, "forward_run.py"], cwd=worker, capture_output=True,
                          text=True)
    assert done.returncode == 0, done.stdout[-2000:] + done.stderr[-2000:]
    return worker


def test_the_forward_run_weights_with_the_realizations_k(run, tmp_path):
    model = mf.load_mf6_run(run)
    cal = model.pest("scr", workspace=tmp_path / "scr", start_datetime="2024-01-01")
    cal.observe(_screened())
    cal.parameterize("k", layers=[2], bounds=(0.1, 20.0))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pst = cal.build("scr.pst", noptmax=0)
    template = cal.template_workspace
    mapping = pd.read_csv(template / "hds_head_target_map.csv")
    assert mapping["obsname"].tolist() == ["W1"] * 3
    # PEST sees ONE observation per period for the well, not one per layer.
    measured = pst.observation_data.loc[pst.observation_data["weight"] > 0]
    assert len(measured) == 2
    base = pd.read_csv(template / "hds_simulated_heads.csv")
    assert base["W1"].tolist() == pytest.approx(
        [_by_hand(_heads(run)[per], K, 70.0, 20.0) for per in (0, 1)])

    group = next(g for g in pst.parameter_data["pargp"].unique() if g.startswith("k"))
    worker = _forward_run(cal, pst, tmp_path / "worker", {group: 10.0})
    k_file = np.array((worker / "m.npf_k_layer3.txt").read_text().split(), dtype=float)
    assert k_file == pytest.approx(10.0 * K[2])  # pyEMU rewrote the file MF6 read
    written = pd.read_csv(worker / "hds_simulated_heads.csv")
    heads = _heads(worker)
    realized = (K[0], K[1], 10.0 * K[2])
    for per in (0, 1):
        expected = _by_hand(heads[per], realized, 70.0, 20.0)
        at_build_k = _by_hand(heads[per], K, 70.0, 20.0)
        assert written.loc[per, "W1"] == pytest.approx(expected, rel=1e-10)
        assert abs(expected - at_build_k) > 1e-3  # the build's K would read another head


def test_a_screened_set_from_hydrographs_builds_and_its_forward_run_reads_it(tmp_path):
    """Where three batches meet: a DATED screened set (``from_hydrographs``; the clock
    check of ledgers 543 and 546, the period-end dates of 544) is built -- its K files found through
    ``_dataset_files``, which grew a fourth field for an empty period block -- and its
    map, a well named ``007`` on three rows, is read back as text (ledger 341's reader,
    which ledger 550 widened) by the forward run, which writes the base run's heads.
    """

    workspace = _write_run(tmp_path / "dated")
    sim = flopy.mf6.MFSimulation.load(sim_ws=str(workspace), verbosity_level=0)
    sim.tdis.start_date_time.set_data("2020-01-01")
    sim.tdis.write()
    model = mf.load_mf6_run(workspace)
    locations = pd.DataFrame({"name": ["007"], "cell": [WELL_CELL],
                              "screen_top": [70.0], "screen_bottom": [20.0]})
    readings = pd.DataFrame({"date": ["2020-01-02", "2020-01-03"], "007": [85.0, 86.0]})
    targets = mf.HeadTargets.from_hydrographs(locations, readings, model=model)
    cal = model.pest("hy", workspace=tmp_path / "hy", start_datetime="2020-01-01")
    cal.observe(targets)
    cal.parameterize("k", layers=[2], bounds=(0.1, 20.0))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pst = cal.build("hy.pst", noptmax=0)
    expected = [_by_hand(_heads(workspace)[per], K, 70.0, 20.0) for per in (0, 1)]
    measured = pst.observation_data.loc[pst.observation_data["weight"] > 0]
    assert measured["obsval"].tolist() == [85.0, 86.0]
    base = pd.read_csv(cal.template_workspace / "hds_simulated_heads.csv", dtype={"007": float})
    assert list(base.columns) == ["per", "007"]
    assert base["007"].tolist() == pytest.approx(expected)
    worker = _forward_run(cal, pst, tmp_path / "worker", {})
    written = pd.read_csv(worker / "hds_simulated_heads.csv", dtype={"007": float})
    assert list(written.columns) == ["per", "007"]
    assert written["007"].tolist() == pytest.approx(expected, rel=1e-10)

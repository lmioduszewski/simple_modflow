"""``periods=`` gives a target one parameter set per group of stress periods (ledger 549).

The Ten Trails recalibration fits an average-year cycle (``mf.AverageYear``: one
steady period, then calendar months) with one recharge multiplier per calendar month
shared by every cycle. Nothing could say that: one multiplier spanned every period file,
and ``temporal=`` was read by nothing (finding build-8). ``periods="month"`` now labels
each TRANSIENT period by the calendar month of its MIDPOINT (``m01``..``m12``), a
mapping gives explicit groups, and each group is one ``add_parameters`` call named
``<name>_<label>``. A steady-state period is in no seasonal group (the user's answer to
Q3): only a parameter declared without ``periods`` reaches it. A period that repeats
another group's block -- MF6's rule for a period with no block of its own -- is refused,
since one file cannot carry two groups' multipliers.

The models are FloPy-written from an ``mf.AverageYear(2, start="2017-10-01")`` -- the
TDIS starts 2017-09-30, as Ten Trails' does, then 24 monthly periods -- with an array
RCH and a CHD carrying a block for every period, and variants whose blocks repeat.
Every forward run is a real MF6 run in a copy of the template.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
import warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6.period_calendar import PeriodCalendar
from myflopy.modflow.mf6.pest import native_parameters as np_mod

pyemu = pytest.importorskip("pyemu")
flopy = pytest.importorskip("flopy")
from flopy.utils.gridutil import get_disv_kwargs  # noqa: E402

NROW, NCOL = 3, 4
NCPL = NROW * NCOL
YEAR = mf.AverageYear(2, start="2017-10-01", nstp=1)
#: Each calendar month's recharge (January = 1), ft/d.
MONTHLY = {month: 1e-4 * month for month in range(1, 13)}


def _run(workspace: Path, *, rch: dict | None = None, perioddata=None, steady=True,
         start=None, chd_every_period=True) -> Path:
    """One layer of NROW x NCOL cells; ``rch`` maps a period to its block's recharge."""

    perioddata = YEAR.perioddata if perioddata is None else perioddata
    sim = flopy.mf6.MFSimulation("s", sim_ws=str(workspace), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=len(perioddata), perioddata=perioddata, time_units="days",
                          start_date_time=YEAR.start_date_time if start is None else start)
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m", save_flows=True)
    flopy.mf6.ModflowGwfdisv(gwf, **get_disv_kwargs(
        1, NROW, NCOL, np.full(NCOL, 10.0), np.full(NROW, 10.0), 100.0, [0.0]))
    flopy.mf6.ModflowGwfnpf(gwf, k=5.0, icelltype=0)
    flopy.mf6.ModflowGwfic(gwf, strt=100.0)
    flopy.mf6.ModflowGwfsto(gwf, sy=0.1, ss=1e-5, iconvert=0,
                            steady_state={0: True} if steady else None, transient={1: True})
    heads = {per: [[(0, 0), 100.0 - 0.01 * per]] for per in range(len(perioddata))} \
        if chd_every_period else [[(0, 0), 100.0]]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=heads)
    rch = YEAR.by_month(MONTHLY) if rch is None else rch
    flopy.mf6.ModflowGwfrcha(gwf, recharge=rch)
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord="m.hds", budget_filerecord="m.cbc",
                           saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])
    sim.write_simulation(silent=True)
    assert sim.run_simulation(silent=True)[0]
    return workspace


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    root = tmp_path_factory.mktemp("seasonal")
    by_month = YEAR.by_month(MONTHLY)
    return {
        "year": _run(root / "year"),
        # Blocks at periods 0 and 1 only: periods 2.. repeat October's block.
        "repeats": _run(root / "repeats", rch={0: by_month[0], 1: by_month[1]}),
        # Blocks at periods 0, 2 and 3: FloPy numbers them _1, _2, _3 (ordinal).
        "ordinal": _run(root / "ordinal", rch={0: 1e-4, 2: 2e-4, 3: 3e-4},
                        perioddata=[(1.0, 1, 1.0)] * 4, chd_every_period=False),
        # One steady day, then a 365-day transient period.
        "long": _run(root / "long", rch={0: 1e-4, 1: 2e-4},
                     perioddata=[(1.0, 1, 1.0), (365.0, 1, 1.0)], chd_every_period=False),
    }


def _cal(runs, tmp_path, kind="year", name="sea"):
    model = mf.load_mf6_run(runs[kind])
    cal = model.pest(name, workspace=tmp_path / name)
    points = gpd.GeoDataFrame({"name": ["W"], "layer": [0]},
                              geometry=gpd.points_from_xy([15.0], [15.0]))
    cal.observe(mf.HeadTargets(locations=points, values=pd.DataFrame(
        {"time": [1], "W": [99.0]})))
    return cal


def _build(cal, name="sea"):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return cal.build(f"{name}.pst")


def _forward_run(cal, pst, tmp_path, values) -> Path:
    worker = tmp_path / "worker"
    shutil.rmtree(worker, ignore_errors=True)
    shutil.copytree(cal.template_workspace, worker)
    par = pst.parameter_data.copy()
    for key, value in values.items():
        names = par.index[par["pargp"] == key] if key in set(par["pargp"]) else [key]
        par.loc[names, "parval1"] = value
    for tpl, inp in zip(pst.template_files, pst.input_files, strict=True):
        pyemu.pst_utils.write_to_template(par["parval1"], str(worker / tpl), str(worker / inp))
    done = subprocess.run([sys.executable, "forward_run.py"], cwd=worker, capture_output=True,
                          text=True)
    assert done.returncode == 0, done.stdout[-2000:] + done.stderr[-2000:]
    assert "Normal termination" in (worker / "mfsim.lst").read_text()
    return worker


def _recharge_by_period(worker: Path) -> np.ndarray:
    """MF6's RCH budget IN per stress period (one step each), from its cell-budget file."""

    budget = flopy.utils.CellBudgetFile(str(worker / "m.cbc"), precision="double")
    return np.array([
        float(sum(record["q"].clip(min=0).sum()
                  for record in budget.get_data(kstpkper=kstpkper, text="RCH")))
        for kstpkper in budget.get_kstpkper()
    ])


# -- months -------------------------------------------------------------------------------


def test_a_period_s_month_is_its_midpoint_s():
    """A TDIS starting 2017-09-30 with a 31-day period: its START is September, its
    midpoint (and every day but the first) October -- m10, where a start rule says m09."""

    calendar = PeriodCalendar("2017-09-30", [31.0, 30.0], "days", steady=(False, False))
    assert np_mod.month_labels(calendar) == ["m10", "m11"]
    assert calendar.starts[0].month == 9
    hours = PeriodCalendar("2017-01-01", [24.0 * 31, 24.0 * 28], "hours")
    assert np_mod.month_labels(hours) == ["m01", "m02"]


def test_month_labels_leave_the_steady_period_out_and_refuse_a_long_period():
    calendar = PeriodCalendar("2017-09-30", [1.0, 31.0], "days", steady=(True, False))
    assert np_mod.month_labels(calendar) == [None, "m10"]
    year = PeriodCalendar("2017-09-30", [1.0, 365.0], "days", steady=(True, False))
    with pytest.raises(ValueError, match=r"transient period 1 runs 365 days.*no month's"):
        np_mod.month_labels(year)


def test_month_gives_one_parameter_per_month_that_reaches_exactly_its_months(runs, tmp_path):
    """Twelve groups; m03 = 2.0 doubles the two March periods and nothing else, and the
    steady period -- in no seasonal group -- is untouched."""

    cal = _cal(runs, tmp_path)
    spec = cal.parameterize("recharge", name="rchm", periods="month", bounds=(0.5, 4.0))
    pst = _build(cal)
    groups = sorted(set(pst.parameter_data["pargp"]))
    assert groups == [f"rchm_m{month:02d}" for month in range(1, 13)]
    march = YEAR.months_of(0)[3], YEAR.months_of(1)[3]
    by_label = {group.label: group for group in spec.period_groups}
    assert by_label["m03"].periods == march
    assert 0 not in {target.period for target in spec.resolved_targets}

    base = _forward_run(cal, pst, tmp_path / "base", {})
    doubled = _forward_run(cal, pst, tmp_path / "march", {"rchm_m03": 2.0})
    ratio = _recharge_by_period(doubled) / _recharge_by_period(base)
    expected = np.ones(YEAR.nper)
    expected[list(march)] = 2.0
    assert np.allclose(ratio, expected)


def test_a_plain_parameter_reaches_the_steady_period_and_stacks_on_the_months(runs, tmp_path):
    """Q3: a non-seasonal multiplier reaches every period, the steady one included; a
    month's multiplier stacks on it in its own months only."""

    cal = _cal(runs, tmp_path)
    cal.parameterize("recharge", name="rch", bounds=(0.5, 4.0))
    cal.parameterize("recharge", name="rchm", periods="month", bounds=(0.5, 4.0))
    pst = _build(cal)
    base = _forward_run(cal, pst, tmp_path / "base", {})
    worker = _forward_run(cal, pst, tmp_path / "both", {"rch": 3.0, "rchm_m03": 2.0})
    ratio = _recharge_by_period(worker) / _recharge_by_period(base)
    expected = np.full(YEAR.nper, 3.0)
    expected[[YEAR.months_of(0)[3], YEAR.months_of(1)[3]]] = 6.0
    assert np.allclose(ratio, expected)
    assert np.isclose(ratio[0], 3.0)


def test_pilot_points_per_month_are_a_set_per_month(runs, tmp_path):
    cal = _cal(runs, tmp_path)
    spec = cal.parameterize("recharge", style="pilotpoints", name="rpp", periods="month",
                            pp_spacing=20.0)
    pst = _build(cal)
    assert [item["name"] for item in spec.pilot_point_sets] == [
        f"rpp_m{month:02d}" for month in range(1, 13)]
    manifest = json.loads((cal.template_workspace / "pilot_point_sets.json").read_text())
    march = next(item for item in manifest["sets"] if item["name"] == "rpp_m03")
    assert march["periods"] == list(np.sort([YEAR.months_of(0)[3], YEAR.months_of(1)[3]]))
    assert len(march["array_files"]) == 2
    base = _forward_run(cal, pst, tmp_path / "base", {})
    worker = _forward_run(cal, pst, tmp_path / "march", {"rpp_m03": 1.5})
    ratio = _recharge_by_period(worker) / _recharge_by_period(base)
    expected = np.ones(YEAR.nper)
    expected[[YEAR.months_of(0)[3], YEAR.months_of(1)[3]]] = 1.5
    assert np.allclose(ratio, expected)


def test_explicit_groups(runs, tmp_path):
    cal = _cal(runs, tmp_path)
    wet = [YEAR.months_of(cycle)[month] for cycle in (0, 1) for month in (11, 12, 1)]
    dry = [YEAR.months_of(cycle)[month] for cycle in (0, 1) for month in (7, 8)]
    spec = cal.parameterize("recharge", name="rch", periods={"Wet": wet, "dry": dry},
                            bounds=(0.2, 5.0))
    pst = _build(cal)
    assert sorted(set(pst.parameter_data["pargp"])) == ["rch_dry", "rch_wet"]
    assert [(group.label, group.periods) for group in spec.period_groups] == [
        ("wet", tuple(sorted(wet))), ("dry", tuple(sorted(dry)))]
    base = _forward_run(cal, pst, tmp_path / "base", {})
    worker = _forward_run(cal, pst, tmp_path / "groups", {"rch_wet": 2.0, "rch_dry": 0.5})
    ratio = _recharge_by_period(worker) / _recharge_by_period(base)
    expected = np.ones(YEAR.nper)
    expected[wet] = 2.0
    expected[dry] = 0.5
    assert np.allclose(ratio, expected)


def test_a_list_target_is_grouped_by_month_too(runs, tmp_path):
    """The CHD carries a block every period: ``chd`` offsets per month, the steady
    period's head left alone."""

    cal = _cal(runs, tmp_path)
    spec = cal.parameterize("chd", name="chdm", periods="month")
    pst = _build(cal)
    assert len(spec.period_groups) == 12 and len(spec.resolved_targets) == 24
    assert 0 not in {target.period for target in spec.resolved_targets}
    worker = _forward_run(cal, pst, tmp_path, {"chdm_m10": 1.25})
    october = {YEAR.months_of(0)[10], YEAR.months_of(1)[10]}
    steady = np_mod.TargetFile("m.chd_stress_period_data_1.txt", "chd_0", None, 0)
    for target in [steady, *spec.resolved_targets]:
        head = pd.read_csv(worker / target.filename, sep=r"\s+", header=None).iloc[0, 2]
        moved = 1.25 if target.period in october else 0.0
        assert np.isclose(head, 100.0 - 0.01 * target.period + moved), target


def test_the_build_metadata_records_each_group(runs, tmp_path):
    cal = _cal(runs, tmp_path)
    cal.parameterize("recharge", name="rchm", periods="month")
    _build(cal)
    metadata = json.loads((cal.template_workspace / "myflopy_pest_metadata.json").read_text())
    (record,) = metadata["parameters"]
    assert record["periods"] == "month" and record["family"] == "array"
    october = next(group for group in record["period_groups"] if group["label"] == "m10")
    assert october["pargps"] == ["rchm_m10"]
    assert october["periods"] == [YEAR.months_of(0)[10], YEAR.months_of(1)[10]]
    assert len(october["files"]) == 2


def test_each_group_records_the_pest_groups_its_parameters_are_in(runs, tmp_path):
    """Review of R2-arrays: a group's record named ``<name>_<label>`` as its one PEST
    group, which no parameter of a zoned pilot-point set is in (``rpp_m10z1``,
    ``rpp_m10z2``): a write-back reading it found no parameters. It names each group's
    real PEST groups, every one of which the control file holds."""

    cal = _cal(runs, tmp_path)
    zones = np.r_[np.ones(NCPL // 2, int), np.full(NCPL - NCPL // 2, 2)]
    cal.parameterize("recharge", style="pilotpoints", name="rpp", periods="month",
                     pp_spacing=15.0, zones=zones)
    cal.parameterize("recharge", style="zone", name="rz", zones=zones,
                     periods={"autumn": [YEAR.months_of(0)[10], YEAR.months_of(0)[11]]})
    pst = _build(cal)
    metadata = json.loads((cal.template_workspace / "myflopy_pest_metadata.json").read_text())
    records = {record["name"]: record for record in metadata["parameters"]}
    held = set(pst.parameter_data["pargp"])
    october = next(group for group in records["rpp"]["period_groups"] if group["label"] == "m10")
    assert october["pargps"] == ["rpp_m10z1", "rpp_m10z2"]
    (autumn,) = records["rz"]["period_groups"]
    assert autumn["pargps"] == ["rz_autumn"]
    for record in records.values():
        named = [pargp for group in record["period_groups"] for pargp in group["pargps"]]
        assert sorted(named) == sorted(record["pargps"]) and set(named) <= held


# -- files numbered by ordinal, repeated blocks -------------------------------------------


def test_groups_follow_the_period_flopy_records_not_the_file_number(runs, tmp_path):
    """Blocks at periods 0, 2, 3 are files _1, _2, _3: a group of period 2 scales _2,
    and period 1 -- which repeats period 0's block -- belongs with period 0."""

    cal = _cal(runs, tmp_path, kind="ordinal")
    spec = cal.parameterize("recharge", name="rch", periods={"a": [0, 1], "b": [2]},
                            bounds=(0.2, 5.0))
    pst = _build(cal)
    assert [(g.label, g.periods, [t.filename for t in g.targets]) for g in spec.period_groups] == [
        ("a", (0, 1), ["m.rcha_recharge_1.txt"]), ("b", (2,), ["m.rcha_recharge_2.txt"])]
    worker = _forward_run(cal, pst, tmp_path, {"rch_b": 3.0})
    assert np.allclose(np.loadtxt(worker / "m.rcha_recharge_2.txt"), 3.0 * 2e-4)
    assert np.allclose(np.loadtxt(worker / "m.rcha_recharge_3.txt"), 3e-4)


@pytest.mark.parametrize("periods, message", [
    ({"a": [0], "b": [1]}, r"no block for period 1 \(group 'b'\).*repeats period 0's \(group 'a'\)"),
    ({"a": [1]}, r"no block for period 1 \(group 'a'\).*repeats period 0's \(in no group\)"),
])
def test_a_period_repeating_another_group_s_block_is_refused(runs, tmp_path, periods, message):
    cal = _cal(runs, tmp_path, kind="ordinal")
    cal.parameterize("recharge", name="rch", periods=periods)
    with pytest.raises(ValueError, match=message):
        _build(cal)


def test_month_on_blocks_that_repeat_across_months_is_refused_naming_the_fix(runs, tmp_path):
    cal = _cal(runs, tmp_path, kind="repeats")
    cal.parameterize("recharge", name="rchm", periods="month")
    with pytest.raises(ValueError, match=r"(?s)no block for period 2 \(group 'm11'\).*"
                                         r"repeats period 1's \(group 'm10'\).*by_month"):
        _build(cal)


def test_a_steady_block_repeated_into_a_month_is_refused(runs, tmp_path):
    """The ``ordinal`` run's CHD has ONE block, at the steady period 0, which MF6
    repeats into every later period: a group there would move the steady head too."""

    cal = _cal(runs, tmp_path, kind="ordinal")
    cal.parameterize("chd", name="chdm", periods={"a": [1]})
    with pytest.raises(ValueError, match=r"period 1 \(group 'a'\).*period 0's \(in no group\)"):
        _build(cal)


def test_month_on_a_period_longer_than_a_month_is_refused(runs, tmp_path):
    cal = _cal(runs, tmp_path, kind="long")
    cal.parameterize("recharge", name="rchm", periods="month")
    with pytest.raises(ValueError, match=r"transient period 1 runs 365 days"):
        _build(cal)


def test_a_group_naming_a_period_the_run_lacks_is_refused(runs, tmp_path):
    cal = _cal(runs, tmp_path, kind="ordinal")
    cal.parameterize("recharge", name="rch", periods={"a": [2, 9]})
    with pytest.raises(ValueError, match=r"group 'a' names period 9, and the run has 4"):
        _build(cal)


def test_a_start_date_disagreeing_with_the_tdis_is_refused_for_months(runs, tmp_path):
    model = mf.load_mf6_run(runs["year"])
    cal = model.pest("sea", workspace=tmp_path / "sea", start_datetime="2020-01-01")
    cal.parameterize("recharge", name="rchm", periods="month")
    with pytest.raises(ValueError, match=r"TDIS starts them on '2017-09-30'.*'2020-01-01'"):
        _build(cal)


# -- the argument itself --------------------------------------------------------------------


@pytest.mark.parametrize("kwargs, error, message", [
    ({"target": "k", "periods": "month"}, ValueError, r"periods= on target 'k'.*not per stress"),
    ({"target": "sto.sy", "periods": {"a": [1]}}, ValueError, r"periods= on target 'sto.sy'"),
    ({"target": "uzf.vks", "periods": "month"}, ValueError, r"periods= on target 'uzf.vks'"),
    ({"target": "recharge", "periods": "months"}, ValueError, r"the one word periods= takes is 'month'"),
    ({"target": "recharge", "periods": {"wet season": [1]}}, ValueError, r"letters, digits and _"),
    ({"target": "recharge", "periods": {"a": [1], "A": [2]}}, ValueError, r"both 'a'"),
    ({"target": "recharge", "periods": {"a": [1, 2], "b": [2]}}, ValueError,
     r"period 2 in groups 'a' and 'b'"),
    ({"target": "recharge", "periods": {"a": []}}, ValueError, r"names no period"),
    ({"target": "recharge", "periods": {"a": [1.5]}}, ValueError, r"zero-based whole number"),
    ({"target": "recharge", "periods": {"a": [1, 1]}}, ValueError, r"names period 1 twice"),
    ({"target": "recharge", "periods": {}}, ValueError, r"names no group"),
    ({"target": "recharge", "periods": [1, 2]}, TypeError, r"takes 'month' or"),
])
def test_periods_that_cannot_act_are_refused_at_the_call(kwargs, error, message):
    with pytest.raises(error, match=message):
        np_mod.NativeParameterSpec(**kwargs)


def test_periods_is_held_normalized():
    assert np_mod.NativeParameterSpec(target="rch", periods=" Month ").periods == "month"
    spec = np_mod.NativeParameterSpec(target="ghb", periods={"Q1": (3, 1, np.int64(2))})
    assert spec.periods == {"q1": (1, 2, 3)}
    assert all(type(value) is int for value in spec.periods["q1"])


def test_temporal_is_gone(runs, tmp_path):
    """Replaced by periods= with no alias (the user is the sole user, Q10)."""

    cal = _cal(runs, tmp_path)
    with pytest.raises(TypeError, match="temporal"):
        cal.parameterize("recharge", temporal=30.0)


# -- a real PESTPP-IES run ------------------------------------------------------------------


@pytest.mark.slow
def test_pestpp_ies_runs_monthly_recharge_and_recharge_pilot_points(runs, tmp_path):
    """PESTPP-IES runs a template carrying a plain recharge multiplier, twelve monthly
    ones and pilot points on the array RCH -- every realization's forward run applies
    them -- and its posterior holds every one of them."""

    model = mf.load_mf6_run(runs["year"])
    cal = model.pest("ies", workspace=tmp_path / "ies")
    points = gpd.GeoDataFrame({"name": ["W", "E"], "layer": [0, 0]},
                              geometry=gpd.points_from_xy([15.0, 35.0], [15.0, 15.0]))
    periods = list(range(1, YEAR.nper, 2))
    cal.observe(mf.HeadTargets(locations=points, values=pd.DataFrame(
        {"time": periods, "W": [100.05] * len(periods), "E": [100.1] * len(periods)})))
    cal.parameterize("recharge", name="rch", bounds=(0.5, 2.0))
    cal.parameterize("recharge", name="rchm", periods="month", bounds=(0.5, 2.0))
    cal.parameterize("recharge", style="pilotpoints", name="rchpp", pp_spacing=20.0,
                     bounds=(0.5, 2.0))
    _build(cal, "ies")
    ies = cal.run_ies(reals=4, iterations=1, workers=2, noise=False,
                      ies_lambda_mults=1.0, lambda_scale_fac=1.0)
    groups = set(cal.pst.parameter_data["pargp"])
    assert {"rch", "rchpp"} | {f"rchm_m{month:02d}" for month in range(1, 13)} <= groups
    posterior = ies.posterior_parameters._df
    assert set(cal.pst.adj_par_names) <= {str(name).lower() for name in posterior.columns}
    assert ies.posterior_iteration >= 1


def test_settings_show_the_period_groups(runs, tmp_path):
    cal = _cal(runs, tmp_path)
    cal.parameterize("recharge", name="rchm", periods="month")
    cal.parameterize("recharge", style="pilotpoints", name="rpp", pp_spacing=20.0)
    _build(cal)
    text = str(cal.settings())
    assert "periods='month'" in text
    assert "rpp: " in text and "(period files)" in text
    assert cal.settings().parameter_frame()["periods"].tolist() == ["month", None]

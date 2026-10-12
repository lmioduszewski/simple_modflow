"""``mf.AverageYear``: one steady period, then repeats of the twelve calendar months (ledger 530).

The period-to-month map an average-year calibration needs -- TDIS, STO flags, every
per-period forcing, and the last cycle the targets are compared on -- stated once.
What these tests pin:

- the periods: non-leap calendar month lengths from ``start``'s month, a steady
  period of ``steady_length`` before it, and a TDIS start that puts the first month
  on ``start`` so every period ends on a month end;
- ``by_month`` keys on CALENDAR month whatever ``start`` is (an October cycle puts
  January at periods 4 and 16), gives the steady period the month-length-weighted
  annual mean (the user's rule, Q3), and refuses what it cannot place;
- ``from_model`` recovers the cycle from a run's own TDIS and STO, live and read back
  from its files, and names the first period that breaks the pattern;
- it works with ``mf.rch.array``, ``mf.sto`` and ``mf.lak`` per-period rainfall and
  evaporation -- each checked against what MF6 WROTE to its budget;
- ``equilibrium`` reads the last two cycles at period end: on a 3-cycle run the two
  agree to the tolerance, a 2-cycle run started far off is flagged, and its numbers
  are the head file's at the period's last step, read by hand.
"""

from __future__ import annotations

import contextlib
import io
import shutil
import warnings
from pathlib import Path

import flopy
import geopandas as gpd
import matplotlib.figure
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import Point

import myflopy as mf
from myflopy import viz
from myflopy.cycles import NONLEAP_DAYS, CycleEquilibrium
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.mf6.observations.hydrographs import model_clock
from myflopy.modflow.mf6.period_calendar import MonthlyCycle, monthly_cycle

#: Calendar month lengths of a non-leap October-September year.
WATER_YEAR = [31, 30, 31, 31, 28, 31, 30, 31, 30, 31, 31, 30]
#: Monthly recharge (length/day), January = 1: wet winters, dry summers.
RATES = {m: 2e-3 * (1.0 + np.cos(2 * np.pi * (m - 1) / 12)) for m in range(1, 13)}
CELL = 100.0
NCPL = 25


# -- the periods ---------------------------------------------------------------------


def test_the_periods_are_calendar_months_after_a_steady_period():
    year = mf.AverageYear(2, start="2016-10-01", nstp=10, tsmult=1.1)
    assert year.nper == 25 and year.months == (10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9)
    assert [row[0] for row in year.perioddata] == [1.0] + WATER_YEAR * 2
    assert year.perioddata[0] == (1.0, 1, 1.0) and year.perioddata[1] == (31.0, 10, 1.1)
    assert year.start_date_time == "2016-09-30"
    periods = year.periods
    transient = periods[~periods["steady"]]
    assert (transient["end"] == transient["end"].dt.normalize()).all()
    assert (transient["end"].dt.day == 1).all()  # each ends at 00:00 after its month's last day
    assert list(transient["start"].dt.month) == list(transient["month"])
    assert list(periods["cycle"][:2].astype("object")) == [pd.NA, 0] and periods["cycle"].iloc[-1] == 1
    assert year.steady_state == {0: True} and year.transient == {1: True}
    assert list(year.calendar.months[1:]) == list(transient["month"])


def test_without_a_steady_period_every_period_is_transient():
    year = mf.AverageYear(1, start="2020-01-01", steady=False, nstp={m: 1 + (m == 2) for m in range(1, 13)})
    assert year.nper == 12 and year.start_date_time == "2020-01-01"
    assert year.steady_state is None and year.transient == {0: True}
    assert year.nstp == {m: 1 + (m == 2) for m in range(1, 13)}
    assert year.perioddata[1] == (28.0, 2, 1.0)  # February, 28 days in a leap year too


def test_the_tdis_is_this_year_and_options_pass_through():
    year = mf.AverageYear(1, start="2016-10-01", steady_length=2.5, nstp=3)
    spec = year.tdis(ats=True, ats_dtmin=1e-3)
    assert spec.options["nper"] == 13 and spec.options["time_units"] == "days"
    assert spec.options["start_date_time"] == "2016-09-28T12:00:00"
    assert list(spec.options["perioddata"]) == year.perioddata
    assert len(spec.options["ats_perioddata"]) == 13
    for owned in ("nper", "perioddata", "time_units", "start_date_time"):
        with pytest.raises(ValueError, match=rf"sets {owned} \(from"):
            year.tdis(**{owned: None})


# -- by_month ------------------------------------------------------------------------


def test_by_month_keys_on_the_calendar_month_whatever_the_start():
    """An October cycle: January is periods 4 and 16, October 1 and 13."""

    year = mf.AverageYear(2, start="2016-10-01", nstp=1)
    placed = year.by_month({m: float(m) for m in range(1, 13)})
    assert sorted(placed) == list(range(25))
    assert [p for p, v in placed.items() if v == 1.0 and p] == [4, 16]
    assert [p for p, v in placed.items() if v == 10.0] == [1, 13]
    assert placed[0] == pytest.approx(sum(m * NONLEAP_DAYS[m] for m in range(1, 13)) / 365)
    assert year.months_of() == {m: 13 + i for i, m in enumerate(year.months)}
    assert year.months_of(0)[1] == 4 and year.months_of(-2) == year.months_of(0)
    assert year.by_month(pd.Series({m: float(m) for m in range(1, 13)})) == placed


def test_by_month_places_arrays_and_lake_dicts_and_averages_them():
    year = mf.AverageYear(1, start="2016-10-01", nstp=1)
    arrays = {m: np.full(4, float(m)) for m in range(1, 13)}
    placed = year.by_month(arrays)
    assert placed[4] is arrays[1], "the same object in every period of its month"
    mean = sum(m * NONLEAP_DAYS[m] for m in range(1, 13)) / 365
    assert np.allclose(placed[0], mean)
    lakes = year.by_month({m: {"pit": 0.001 * m, "pond": 0.0} for m in range(1, 13)})
    assert lakes[4] == {"pit": 0.001, "pond": 0.0}
    assert lakes[0]["pit"] == pytest.approx(0.001 * mean) and lakes[0]["pond"] == 0.0
    records = {m: [[(0, 1), float(m)]] for m in range(1, 13)}
    assert year.by_month(records, steady=[[(0, 1), 0.5]])[0] == [[(0, 1), 0.5]]


@pytest.mark.parametrize(
    "values, kwargs, error, match",
    [
        ([1.0] * 12, {}, TypeError, "ambiguous between calendar order and cycle order"),
        ({m: 1.0 for m in range(1, 12)}, {}, ValueError, "missing: Dec"),
        ({**{m: 1.0 for m in range(1, 13)}, 13: 1.0}, {}, ValueError, r"\[13\] is not"),
        ({**{m: 1.0 for m in range(1, 13)}, 3: None}, {}, ValueError, "Mar hold None"),
        ({m: [[(0, 1), 1.0]] for m in range(1, 13)}, {}, ValueError, r"steady=<value>"),
        ({m: (1.0 if m < 6 else np.ones(3)) for m in range(1, 13)}, {}, ValueError, "steady=<value>"),
        ({m: {"a": 1.0} if m < 6 else {"b": 1.0} for m in range(1, 13)}, {}, ValueError, "same"),
    ],
)
def test_what_by_month_cannot_place_is_refused(values, kwargs, error, match):
    with pytest.raises(error, match=match):
        mf.AverageYear(1, start="2016-10-01", nstp=1).by_month(values, **kwargs)


def test_a_steady_value_without_a_steady_period_is_refused():
    year = mf.AverageYear(1, start="2016-10-01", steady=False, nstp=1)
    with pytest.raises(ValueError, match="no steady period"):
        year.by_month({m: 1.0 for m in range(1, 13)}, steady=1.0)
    assert sorted(year.by_month({m: 1.0 for m in range(1, 13)})) == list(range(12))


@pytest.mark.parametrize(
    "args, kwargs, error, match",
    [
        ((0,), dict(start="2016-10-01", nstp=1), ValueError, "at least one cycle"),
        ((True,), dict(start="2016-10-01", nstp=1), TypeError, "whole number of years"),
        ((1.0,), dict(start="2016-10-01", nstp=1), TypeError, "whole number of years"),
        ((1,), dict(start="2016-10-15", nstp=1), ValueError, r"start='2016-10-01'"),
        ((1,), dict(start="2016-10-01T06:00", nstp=1), ValueError, r"start='2016-10-01'"),
        ((1,), dict(start="October", nstp=1), ValueError, "with its year"),
        ((1,), dict(start="2016-10-01", nstp=0), ValueError, "must be >= 1"),
        ((1,), dict(start="2016-10-01", nstp={m: 1 for m in range(1, 12)}), ValueError, "missing: Dec"),
        ((1,), dict(start="2016-10-01", nstp={**{m: 1 for m in range(1, 13)}, 0: 1}), ValueError, "months 1-12"),
        ((1,), dict(start="2016-10-01", nstp=1.5), TypeError, "an int or"),
        ((1,), dict(start="2016-10-01", nstp=1, tsmult=0.0), ValueError, "tsmult"),
        ((1,), dict(start="2016-10-01", nstp=1, steady_length=0.0), ValueError, "steady_length"),
        ((1,), dict(start="2016-10-01", nstp=1, steady=False, steady_length=1.0), ValueError,
         "no steady period"),
        ((1,), dict(start="2016-10-01", nstp=1, steady="yes"), TypeError, "True or False"),
    ],
)
def test_what_cannot_be_an_average_year_is_refused(args, kwargs, error, match):
    with pytest.raises(error, match=match):
        mf.AverageYear(*args, **kwargs)


def test_months_of_refuses_a_cycle_the_year_does_not_run():
    year = mf.AverageYear(2, start="2016-10-01", nstp=1)
    with pytest.raises(IndexError, match="runs 2 cycle"):
        year.months_of(2)
    with pytest.raises(TypeError):
        year.months_of(True)


# -- runs ----------------------------------------------------------------------------


def _grid(n=5, d=CELL):
    verts = np.array([[i * d, j * d] for j in range(n + 1) for i in range(n + 1)], dtype=float)
    iverts = [
        [j * (n + 1) + i, (j + 1) * (n + 1) + i, (j + 1) * (n + 1) + i + 1, j * (n + 1) + i + 1]
        for j in range(n)
        for i in range(n)
    ]
    xcyc = np.array([[(i + 0.5) * d, (j + 0.5) * d] for j in range(n) for i in range(n)])
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)


def _run_square(root: Path, year, *, strt=50.0, sy=0.05, tdis=None):
    """A 5 x 5 unconfined aquifer drained by a GHB along its west edge, recharged by month."""

    vor = _grid()
    context = mf.ModelContext(grid=vor, domain=np.ones((1, NCPL), dtype=int))
    model = mf.gwf("flow", context=context, packages=[
        mf.disv(nlay=1, top=100.0, botm=[0.0], **vor.get_disv_gridprops()),
        mf.ic(strt=strt),
        mf.npf(k=5.0, icelltype=1),
        mf.sto(sy=sy, ss=1e-5, iconvert=1, steady_state=year.steady_state, transient=year.transient),
        mf.ghb(stress_period_data={0: [[(0, c), 50.0, 100.0] for c in range(0, NCPL, 5)]}),
        mf.rch.array(recharge=year.by_month({m: np.full(NCPL, r) for m, r in RATES.items()}),
                     irch=0, context=context),
        mf.oc(head_filerecord="flow.hds", budget_filerecord="flow.cbc",
              saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[tdis or year.tdis()])
    shutil.rmtree(root, ignore_errors=True)
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        run = mf.Project(root, name="p").prepare_run("r", simulation)
        success, report = run.execute(silent=True)
    assert success, "\n".join(report[-15:])
    return run


def _pond(grid):
    xmin, ymin, xmax, ymax = grid.gdf_vorPolys.total_bounds
    return gpd.GeoDataFrame({"name": ["pit"]}, crs=grid.crs,
                            geometry=[Point((xmin + xmax) / 2, (ymin + ymax) / 2).buffer(30.0, 32)])


def _run_lake(root: Path, year):
    """An excavated pond with monthly rainfall and evaporation, NEWTON, CHD in the lower layer."""

    grid = irregular_voronoi_grid(nrow=12, ncol=12, cell_size=20.0, refined=False)
    pond = _pond(grid)
    result = (mf.LayerStack(top=mf.Flat(100.0)).add("upper", bottom=mf.Flat(82.0))
              .add("lower", bottom=mf.Flat(40.0)).excavate(bottom=82.0, within=pond, name="pit")
              .build(grid))
    context = mf.ModelContext(grid=grid, domain=result.idomain, surfaces=result)
    xy = np.array([(c[1], c[2]) for c in grid.get_disv_gridprops()["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 20) | (xy[:, 0] > xmax - 20))
    rain = {m: 0.004 * (1.0 + np.cos(2 * np.pi * (m - 1) / 12)) for m in range(1, 13)}
    evap = {m: {"pit": 0.003 * (1.0 - np.cos(2 * np.pi * (m - 1) / 12))} for m in range(1, 13)}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        lake = mf.lak(context=context, nper=year.nper, lakes=pond, lake_id_field="name",
                      starting_stage=95.0, bed_leakance=1.0, unexcavated="raise",
                      rainfall=year.by_month(rain), evaporation=year.by_month(evap))
    model = mf.gwf("m", context=context, newtonoptions="NEWTON", packages=[
        mf.disv(result), mf.npf(k=10.0, k33=10.0, icelltype=1), mf.ic(strt=95.0),
        mf.chd(stress_period_data={0: [[(1, int(c)), 94.0] for c in edge]}),
        mf.sto(ss=1e-5, sy=0.2, iconvert=1, steady_state=year.steady_state, transient=year.transient),
        lake,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        year.tdis(), mf.ims(complexity="COMPLEX", linear_acceleration="BICGSTAB")])
    shutil.rmtree(root, ignore_errors=True)
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        run = mf.Project(root, name="p").prepare_run("r", simulation)
        success, report = run.execute(silent=True)
    assert success, "\n".join(report[-15:])
    return run, rain, evap


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    root = tmp_path_factory.mktemp("average_year")
    three = mf.AverageYear(3, start="2016-10-01", nstp=2)
    unsettled = mf.AverageYear(2, start="2016-10-01", steady=False, nstp=2)
    lake_year = mf.AverageYear(2, start="2016-10-01", nstp=1)
    lake_run, rain, evap = _run_lake(root / "lake", lake_year)
    return {
        "three": (three, _run_square(root / "three", three)),
        "unsettled": (unsettled, _run_square(root / "unsettled", unsettled, strt=99.0, sy=0.3)),
        "lake": (lake_year, lake_run, rain, evap),
    }


@pytest.mark.slow
def test_from_model_recovers_the_cycle_live_and_from_the_files(runs):
    for key in ("three", "unsettled"):
        year, run = runs[key]
        assert mf.AverageYear.from_model(run.model("flow")) == year
        assert mf.AverageYear.from_model(mf.load_mf6_run(run.workspace)) == year
    year, run = runs["three"]
    assert run.model("flow").calendar.steady == year.calendar.steady  # mf.sto took the flags


@pytest.mark.slow
def test_a_runs_cycle_is_read_in_one_place(runs):
    """Ledger 546: ``monthly_cycle`` is the one reading of a run's cycle, below the
    observation layer -- ``from_model`` is it, and ``from_hydrographs`` places its
    average-year targets on its last cycle without importing ``myflopy.cycles``."""

    for key in ("three", "unsettled"):
        year, run = runs[key]
        for model in (run.model("flow"), mf.load_mf6_run(run.workspace)):
            cycle = monthly_cycle(model)
            assert isinstance(cycle, MonthlyCycle) and cycle.nper == year.nper
            assert [cycle.months_of(c) for c in range(year.cycles)] == [
                year.months_of(c) for c in range(year.cycles)
            ]
            assert model_clock(model, mode="average_year", owner="HeadTargets")[1] == year.months_of(-1)


@pytest.mark.slow
def test_from_model_names_the_period_that_breaks_the_pattern(tmp_path):
    year = mf.AverageYear(1, start="2016-10-01", nstp=1)
    edited = list(year.perioddata)
    edited[5] = (29.0, 1, 1.0)  # February
    run = _run_square(tmp_path / "edited", year, tdis=mf.tdis(
        nper=year.nper, perioddata=edited, time_units="days", start_date_time=year.start_date_time))
    with pytest.raises(ValueError, match="Period 5 is 29 days long where an AverageYear's Feb is 28"):
        mf.AverageYear.from_model(run.model("flow"))
    short = mf.tdis(nper=12, perioddata=year.perioddata[:12], time_units="days",
                    start_date_time=year.start_date_time)
    run = _run_square(tmp_path / "short", mf.AverageYear(1, start="2016-10-01", nstp=1), tdis=short)
    with pytest.raises(ValueError, match="11 transient period"):
        mf.AverageYear.from_model(run.model("flow"))


@pytest.mark.slow
def test_monthly_recharge_reaches_mf6_as_rate_times_area(runs):
    """MF6's own RCH budget, per period: the month's rate over the 25 cells' area."""

    year, run = runs["three"]
    reader = run.model("flow")._get_budget_reader()
    month_of = {per: month for cycle in range(year.cycles) for month, per in year.months_of(cycle).items()}
    mean = sum(RATES[m] * NONLEAP_DAYS[m] for m in range(1, 13)) / 365
    for kstpkper in reader.get_kstpkper():
        per = int(kstpkper[1])
        q = sum(float(np.sum(r["q"])) for r in reader.get_data(text="RCH", kstpkper=kstpkper))
        rate = mean if per == 0 else RATES[month_of[per]]
        assert q == pytest.approx(rate * NCPL * CELL * CELL, rel=1e-6)


@pytest.mark.slow
def test_monthly_lake_rainfall_and_evaporation_reach_mf6(runs):
    """MF6's LAK budget: rainfall and evaporation per period in proportion to the month's rate."""

    year, run, rain, evap = runs["lake"]
    model = run.model("m")
    budget = model.fl.lak.output.budget()
    month_of = {per: month for cycle in range(year.cycles) for month, per in year.months_of(cycle).items()}
    ratios = {"RAINFALL": [], "EVAPORATION": []}
    for kstpkper in budget.get_kstpkper():
        per = int(kstpkper[1])
        if per == 0:
            continue
        month = month_of[per]
        for text, rate in (("RAINFALL", rain[month]), ("EVAPORATION", evap[month]["pit"])):
            q = sum(float(np.sum(r["q"])) for r in budget.get_data(text=text, kstpkper=kstpkper))
            if rate > 0:
                ratios[text].append(abs(q) / rate)
            else:
                assert q == 0.0
    # The pit's surface area, the same every month (vertical walls, stage inside the pit).
    for text, values in ratios.items():
        assert np.ptp(values) <= 1e-6 * np.mean(values), text
    assert ratios["RAINFALL"][0] == pytest.approx(ratios["EVAPORATION"][0], rel=1e-6)


@pytest.mark.slow
def test_three_cycles_settle_and_the_numbers_are_the_head_files(runs):
    year, run = runs["three"]
    model = run.model("flow")
    equilibrium = year.equilibrium(model, tolerance=0.01)
    summary = equilibrium.summary()
    assert list(summary["source"]) == ["heads"] and summary["n"].iloc[0] == NCPL
    assert summary["max_abs_change"].iloc[0] < 0.01 and summary["passes"].iloc[0]
    assert equilibrium.converged

    # The table holds the head file's value at each period's LAST step (nstp=2), by hand.
    heads = flopy.utils.HeadFile(Path(run.workspace) / "flow.hds")
    table = equilibrium.get()
    row = table[(table["cell"] == 12) & (table["month"] == 1)].iloc[0]
    previous, last = year.months_of(1)[1], year.months_of(2)[1]
    assert (row["per_previous"], row["per_last"]) == (previous, last)
    assert row["previous"] == heads.get_data(kstpkper=(1, previous)).ravel()[12]
    assert row["last"] == heads.get_data(kstpkper=(1, last)).ravel()[12]
    assert row["change"] == row["last"] - row["previous"]
    assert len(table) == 12 * NCPL

    # At targets: the set's own simulated heads, the same numbers as the cells'.
    targets = mf.HeadTargets(locations=pd.DataFrame({"name": ["w12"], "layer": [0], "cell": [12]}),
                             values=pd.DataFrame({"per": [last], "w12": [60.0]}), time_column="per")
    at = year.equilibrium(model, targets=targets).get()
    assert list(at["source"].unique()) == ["HeadTargets"]
    jan = at[at["month"] == 1].iloc[0]
    assert (jan["name"], jan["previous"], jan["last"]) == ("w12", row["previous"], row["last"])
    assert targets.simulated_heads(model).set_index("per").loc[last, "w12"] == row["last"]


@pytest.mark.slow
def test_a_run_that_has_not_settled_is_flagged(runs):
    """Two cycles from heads 49 ft too high: the second year still drains the first's."""

    year, run = runs["unsettled"]
    equilibrium = year.equilibrium(run.model("flow"), tolerance=0.01)
    summary = equilibrium.summary()
    assert summary["max_abs_change"].iloc[0] > 1.0
    assert not summary["passes"].iloc[0] and not equilibrium.converged
    assert summary["annual_mean_drift"].iloc[0] > 1.0
    with pytest.raises(ValueError, match="No tolerance was given"):
        year.equilibrium(run.model("flow")).converged


def _execute(root: Path, simulation):
    shutil.rmtree(root, ignore_errors=True)
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        run = mf.Project(root, name="p").prepare_run("r", simulation)
        success, report = run.execute(silent=True)
    assert success, "\n".join(report[-15:])
    return run


@pytest.mark.slow
def test_a_cell_dry_in_one_cycle_only_is_unsettled(tmp_path):
    """Ledger 531: a thin top layer wet at the first cycle's period ends drains
    dry (MF6 writes -1e30) by the second's. Its cells have a value in ONE cycle only;
    dropped, the report read every remaining cell and passed a loose tolerance."""

    year = mf.AverageYear(2, start="2016-10-01", steady=False, nstp=2)
    vor = _grid()
    context = mf.ModelContext(grid=vor, domain=np.ones((2, NCPL), dtype=int))
    model = mf.gwf("flow", context=context, packages=[
        mf.disv(nlay=2, top=100.0, botm=[85.0, 0.0], **vor.get_disv_gridprops()),
        mf.ic(strt=99.0), mf.npf(k=5.0, icelltype=1),
        mf.sto(sy=0.3, ss=1e-5, iconvert=1, steady_state=year.steady_state, transient=year.transient),
        mf.ghb(stress_period_data={0: [[(1, c), 50.0, 100.0] for c in range(0, NCPL, 5)]}),
        mf.rch.array(recharge=year.by_month({m: np.full(NCPL, r) for m, r in RATES.items()}),
                     irch=0, context=context),
        mf.oc(head_filerecord="flow.hds", saverecord=[("HEAD", "LAST")]),
    ])
    run = _execute(tmp_path / "dry", mf.simulation(model, name="s", packages=[year.tdis()]))
    heads = flopy.utils.HeadFile(Path(run.workspace) / "flow.hds")
    october = year.months_of(0)[10], year.months_of(1)[10]
    top = [heads.get_data(kstpkper=(1, per))[0].ravel() for per in october]
    one_sided = int(np.sum((np.abs(top[0]) < 1e29) != (np.abs(top[1]) < 1e29)))
    assert one_sided > 0, "the fixture must dry cells between the cycles"

    equilibrium = year.equilibrium(run.model("flow"), tolerance=1e6)
    summary = equilibrium.summary().iloc[0]
    assert summary["max_abs_change"] < 1e6  # every two-sided change is inside the tolerance
    assert summary["n_one_sided"] > 0 and not summary["passes"] and not equilibrium.converged
    table = equilibrium.get()
    rows = table[(table["layer"] == 0) & (table["month"] == 10)]
    assert len(rows) == one_sided and rows["change"].isna().all()
    assert rows["previous"].notna().all() and rows["last"].isna().all()


@pytest.mark.slow
def test_a_structured_model_is_compared_too(tmp_path):
    """Ledger 531: the head file is FloPy's reader of the OC's file, so a DIS
    model is read; `model.hds` (DISV-only) raised AttributeError on it."""

    year = mf.AverageYear(2, start="2016-10-01", nstp=1)
    model = mf.gwf("flow", packages=[
        mf.dis(nlay=1, nrow=3, ncol=4, delr=100.0, delc=100.0, top=100.0, botm=[0.0]),
        mf.ic(strt=60.0), mf.npf(k=5.0, icelltype=1),
        mf.sto(sy=0.1, ss=1e-5, iconvert=1, steady_state=year.steady_state, transient=year.transient),
        mf.chd(stress_period_data={0: [[(0, 0, 0), 50.0]]}),
        mf.rch.array(recharge=year.by_month({m: 1e-3 * m for m in range(1, 13)}), irch=0),
        mf.oc(head_filerecord="flow.hds", saverecord=[("HEAD", "LAST")]),
    ])
    run = _execute(tmp_path / "dis", mf.simulation(model, name="s", packages=[year.tdis()]))
    table = year.equilibrium(run.model("flow")).get()
    assert len(table) == 12 * 12 and set(table["cell"]) == set(range(12))
    heads = flopy.utils.HeadFile(Path(run.workspace) / "flow.hds")
    row = table[(table["cell"] == 7) & (table["month"] == 3)].iloc[0]
    assert row["last"] == heads.get_data(kstpkper=(0, row["per_last"]))[0, 1, 3]


@pytest.mark.slow
def test_conc_targets_are_compared_on_the_transport_model(tmp_path):
    """Ledger 531: a transport model's calendar carries no STO flags, and
    comparing them refused every ConcTargets call; with targets=None it is refused by
    name, since it has no heads."""

    year = mf.AverageYear(2, start="2016-10-01", nstp=1)
    vor = _grid()
    context = mf.ModelContext(grid=vor, domain=np.ones((1, NCPL), dtype=int))
    disv = dict(nlay=1, top=100.0, botm=[0.0], **vor.get_disv_gridprops())
    flow = mf.gwf("flow", context=context, packages=[
        mf.disv(**disv), mf.ic(strt=50.0), mf.npf(k=5.0, icelltype=1),
        mf.sto(sy=0.05, ss=1e-5, iconvert=1, steady_state=year.steady_state, transient=year.transient),
        mf.ghb(stress_period_data={0: [[(0, c), 50.0, 100.0, 0.0] for c in range(0, NCPL, 5)]},
               auxiliary=["CONCENTRATION"], pname="ghb"),
        mf.rch.array(recharge=year.by_month({m: np.full(NCPL, r) for m, r in RATES.items()}),
                     irch=0, context=context),
        mf.oc(head_filerecord="flow.hds", budget_filerecord="flow.cbc",
              saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")]),
    ])
    transport = mf.gwt("tport", context=context, packages=[
        mf.disv(**disv), mf.ic(strt=0.0), mf.adv(scheme="upstream"), mf.mst(porosity=0.3),
        mf.ssm(sources=[("ghb", "AUX", "CONCENTRATION")]),
        mf.cnc(stress_period_data={0: [[(0, 24), 100.0]]}),
        mf.oc(concentration_filerecord="tport.ucn", saverecord=[("CONCENTRATION", "LAST")]),
    ])
    run = _execute(tmp_path / "gwt", mf.simulation(
        flow, transport, name="s", packages=[year.tdis()], exchanges=[mf.exchange(flow, transport)]))
    model = run.model("tport")
    assert model.calendar.steady is None
    conc = mf.ConcTargets(locations=pd.DataFrame({"name": ["c23"], "layer": [0], "cell": [23]}),
                          values=pd.DataFrame({"per": [24], "c23": [1.0]}), time_column="per")
    table = year.equilibrium(model, targets=conc).get()
    simulated = conc.simulated_heads(model).set_index("per")["c23"]
    assert len(table) == 12 and (table["source"] == "ConcTargets").all()
    assert list(table["last"]) == [simulated.loc[per] for per in table["per_last"]]
    with pytest.raises(ValueError, match="compares a FLOW model's heads"):
        year.equilibrium(model)


@pytest.mark.slow
def test_equilibrium_reads_every_lake_and_agrees_with_a_lake_target(runs):
    year, run, _, _ = runs["lake"]
    model = run.model("m")
    table = year.equilibrium(model).get()
    assert set(table["source"]) == {"heads", "lake stage (lak)"}
    lake = table[table["source"] == "lake stage (lak)"]
    assert set(lake["name"]) == {"lak:0"} and len(lake) == 12
    series = mf.LakeStageTargets(locations={"pit": 0}).simulated_series(model).set_index("time")["pit"]
    row = lake[lake["month"] == 7].iloc[0]
    assert row["previous"] == series.loc[row["per_previous"]]
    assert row["last"] == series.loc[row["per_last"]]


@pytest.mark.slow
def test_equilibrium_refusals(runs, tmp_path):
    year, run = runs["three"]
    model = run.model("flow")
    with pytest.raises(ValueError, match="runs one cycle"):
        mf.AverageYear(1, start="2016-10-01", nstp=2).equilibrium(model)
    with pytest.raises(ValueError, match=r"not .*AverageYear\(4.*from_model"):
        mf.AverageYear(4, start="2016-10-01", nstp=2).equilibrium(model)
    with pytest.raises(ValueError, match="it starts 2016-09-30 where the year starts 2017-09-30"):
        mf.AverageYear(3, start="2017-10-01", nstp=2).equilibrium(model)
    with pytest.raises(ValueError, match="tolerance"):
        year.equilibrium(model, tolerance=-1.0)
    with pytest.raises(TypeError, match="neither simulated_heads nor simulated_series"):
        year.equilibrium(model, targets=[object()])


@pytest.mark.slow
def test_the_plot_draws_the_last_two_cycle_pairs(runs):
    year, run = runs["three"]
    equilibrium = year.equilibrium(run.model("flow"), tolerance=0.01)
    fig = equilibrium.plot()
    assert isinstance(fig, viz.Fig)
    names = [trace.name for trace in fig.data]
    assert names == ["heads: cycle 1 vs 0", "heads: cycle 2 vs 1"]
    assert list(fig.data[0].x) == ["Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                   "Jul", "Aug", "Sep"]
    assert max(fig.data[1].y) < max(fig.data[0].y), "the later pair moves less"
    assert isinstance(equilibrium.plot(backend="mpl"), matplotlib.figure.Figure)
    assert isinstance(equilibrium, CycleEquilibrium)

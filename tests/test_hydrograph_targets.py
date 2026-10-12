"""Dated readings binned onto stress periods: ``from_hydrographs`` (ledgers 543, 544).

A target set is keyed by stress period (ledger 340) and a date is refused as a key
(ledger 499). ``HeadTargets.from_hydrographs`` / ``LakeStageTargets.from_hydrographs``
bin dated readings onto a model's own calendar (ledger 529) by one rule
(``observations/hydrographs.py``). What these tests pin:

- history mode, ``rule="nearest_end"``: each period's target is the reading
  nearest its END -- on Ten Trails' TDIS (start 2017-09-30) a reading on 31
  October at noon is period 0's, where containment would give it to period 1; a
  tie goes to the earlier reading; a reading midway between two ends goes to the
  period holding it; ``median``/``mean`` read the readings inside ``(start, end]``;
- measured on MF6: readings that ARE the model's own heads at every time step bin
  back to its period-end heads exactly under ``nearest_end``, and miss them under
  ``median``;
- average-year mode: per calendar month, the median over years of each year's
  monthly median (the user's rule, Q4), on the LAST cycle's periods of an
  ``mf.AverageYear`` model -- and a model that is not one is refused naming it;
- the CLOCK STAMP: a set binned on a 3-cycle run is refused at ``compare()`` on a
  2-cycle one and by the PEST build on another model, naming ``from_hydrographs``;
- every refusal the design lists, by message, and out-of-span readings WARNED by
  name (ledger 288);
- a real PESTPP-IES run: the control file holds the binning report's values, the
  build metadata the period-end dates and the binning, and the review draws the
  ensembles against those dates (ies-12) -- while a template without them keeps
  period numbers;
- lake stages bin by the same rule.
"""

from __future__ import annotations

import contextlib
import io
import json
import shutil
import warnings
from pathlib import Path

import flopy
import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6.observations.hydrographs import (
    REPORT_COLUMNS,
    bin_readings,
    read_readings,
)
from myflopy.modflow.mf6.period_calendar import PeriodCalendar
from myflopy.modflow.mf6.pest.ies import IesResults, _panel_axis
from tests.test_average_year import _grid, _run_lake

#: Ten Trails' TDIS: START_DATE_TIME 2017-09-30, monthly periods from October.
TT_START = "2017-09-30"
TT_LENGTHS = [31, 30, 31, 31, 28, 31, 30, 31, 30, 31, 31, 30]
NCPL = 25


def _calendar() -> PeriodCalendar:
    return PeriodCalendar(TT_START, TT_LENGTHS, "days")


def _long(rows) -> pd.DataFrame:
    """Readings as ``read_readings`` returns them, from ``(name, date, value)`` rows."""

    frame = pd.DataFrame(rows, columns=["name", "date", "value"])
    frame["date"] = pd.DatetimeIndex([pd.Timestamp(value) for value in frame["date"]])
    return frame


def _row(binning, name, per):
    report = binning.report
    return report.loc[(report["name"] == name) & (report["per"] == per)].iloc[0]


# -- history: the rule, on a calendar -----------------------------------------------


def test_nearest_end_puts_a_reading_on_the_last_day_of_a_month_in_its_own_period():
    """Period 0 ends at 00:00 on 31 October: a reading that noon is 12 h after it.

    A containment rule ((start, end]) gives it to period 1; the reading nearest
    period 0's END is it, and that is the value MF6 reports for period 0.
    """

    binning = bin_readings(
        _long([("W1", "2017-10-20", 1.0), ("W1", "2017-10-31 12:00", 2.0), ("W1", "2017-11-20", 3.0)]),
        _calendar(),
    )
    first = _row(binning, "W1", 0)
    assert first["reading"] == pd.Timestamp("2017-10-31 12:00") and first["value"] == 2.0
    assert first["offset_days"] == pytest.approx(0.5)
    assert first["n_readings"] == 2  # 20 Oct and 31 Oct noon are both nearest period 0's end
    assert _row(binning, "W1", 1)["value"] == 3.0
    assert list(binning.values.columns) == ["time", "name", "value"]
    assert binning.values.to_dict("list") == {"time": [0, 1], "name": ["W1", "W1"], "value": [2.0, 3.0]}
    assert list(binning.report.columns) == list(REPORT_COLUMNS)
    assert binning.report["month"].tolist() == [10, 11]


def test_a_tie_goes_to_the_earlier_reading():
    binning = bin_readings(
        _long([("W1", "2017-11-01", 5.0), ("W1", "2017-10-30", 4.0)]), _calendar()
    )
    chosen = _row(binning, "W1", 0)
    assert chosen["value"] == 4.0 and chosen["offset_days"] == pytest.approx(-1.0)


def test_a_reading_midway_between_two_ends_goes_to_the_period_holding_it():
    # Period 1 runs 2017-10-31 to 2017-11-30: its midpoint is 15 November 00:00.
    binning = bin_readings(_long([("W1", "2017-11-15", 7.0)]), _calendar())
    assert binning.values["time"].tolist() == [1]


def test_median_and_mean_reduce_the_readings_inside_the_period():
    rows = [
        ("W1", "2017-10-05", 1.0),
        ("W1", "2017-10-20", 2.0),
        ("W1", "2017-10-31", 6.0),        # exactly period 0's end: inside (start, end]
        ("W1", "2017-10-31 12:00", 50.0),  # past it: period 1
    ]
    median = bin_readings(_long(rows), _calendar(), rule="median")
    mean = bin_readings(_long(rows), _calendar(), rule="mean")
    assert _row(median, "W1", 0)["value"] == 2.0 and _row(median, "W1", 0)["n_readings"] == 3
    assert _row(mean, "W1", 0)["value"] == pytest.approx(3.0)
    assert _row(mean, "W1", 1)["value"] == 50.0
    assert pd.isna(_row(median, "W1", 0)["reading"])


def test_max_offset_leaves_a_far_period_unmeasured():
    rows = [("W1", "2017-10-29", 1.0), ("W1", "2017-11-12", 2.0)]
    binning = bin_readings(_long(rows), _calendar(), max_offset="5D")
    assert binning.values["time"].tolist() == [0]  # 12 Nov is 18 days from either end


def test_readings_outside_the_span_are_left_out_and_warned_by_name():
    rows = [("W1", "2017-09-01", 1.0), ("W1", "2017-10-30", 2.0),
            ("W2", "2016-01-01", 3.0), ("W2", "2030-01-01", 4.0), ("W2", "2017-12-30", 5.0)]
    with pytest.warns(UserWarning, match=r"3 reading\(s\) fall outside the model's span.*W1 \(1\), W2 \(2\)"):
        binning = bin_readings(_long(rows), _calendar())
    assert binning.values["value"].tolist() == [2.0, 5.0]


@pytest.mark.filterwarnings("ignore:.*fall outside the model's span:UserWarning")
def test_a_location_with_no_reading_left_is_named():
    rows = [("W1", "2017-10-30", 2.0), ("W2", "2016-01-01", 3.0)]
    with pytest.warns(UserWarning, match=r"1 of 2 location\(s\) have no reading .*: W2"):
        bin_readings(_long(rows), _calendar())


def test_no_reading_in_the_span_is_refused_printing_both_spans():
    with pytest.raises(ValueError, match=r"no reading falls in the model's span, 2017-09-30T00:00:00 to "
                                         r"2018-09-30T00:00:00 \(the readings run 2001-01-01"):
        bin_readings(_long([("W1", "2001-01-01", 1.0)]), _calendar())


# -- average year: the rule, on a calendar ------------------------------------------


def _year_calendar():
    year = mf.AverageYear(2, start="2016-10-01", nstp=1)
    return year.calendar, year.months_of(-1)


def test_average_year_is_the_median_over_years_of_each_years_monthly_median():
    """A January read daily in one year counts once, as a year read once does (Q4)."""

    calendar, months = _year_calendar()
    rows = [("W1", day, 10.0) for day in pd.date_range("2001-01-01", "2001-01-31")]
    rows += [("W1", "2002-01-15", 1.0), ("W1", "2003-01-15", 2.0)]
    binning = bin_readings(_long(rows), calendar, mode="average_year", months=months)
    row = _row(binning, "W1", months[1])
    assert row["value"] == 2.0, "median of the yearly medians [10, 1, 2]"
    assert row["n_years"] == 3 and row["n_readings"] == 33 and row["month"] == 1
    assert np.median([10.0] * 31 + [1.0, 2.0]) == 10.0, "a pooled median would be the daily year's"
    mean = bin_readings(_long(rows), calendar, mode="average_year", rule="mean", months=months)
    assert _row(mean, "W1", months[1])["value"] == pytest.approx(13.0 / 3)
    assert months[1] == 16, "January of the LAST of two October cycles after the steady period"


def test_min_years_and_between_choose_what_counts():
    calendar, months = _year_calendar()
    rows = [("W1", "2001-01-15", 1.0), ("W1", "2002-01-15", 3.0), ("W1", "2002-02-15", 9.0),
            ("W1", "2010-01-15", 100.0)]
    binning = bin_readings(_long(rows), calendar, mode="average_year", months=months,
                           min_years=2, between=(None, "2005-12-31"))
    assert binning.values.to_dict("list") == {"time": [months[1]], "name": ["W1"], "value": [2.0]}
    february = _row(binning, "W1", months[2])
    assert not february["measured"] and pd.isna(february["value"])


# -- refusals -------------------------------------------------------------------------


@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({"mode": "weekly"}, r"mode='weekly'; it is 'history' or 'average_year'"),
        ({"rule": "last"}, r"rule='last'; mode='history' takes 'nearest_end', 'median', 'mean'"),
        ({"between": ("2001-01-01", None)}, r"between is read in mode='average_year' only"),
        ({"min_years": 2}, r"min_years is read in mode='average_year' only"),
        ({"max_offset": 5}, r"max_offset=5 has no unit"),
        ({"max_offset": "-2D"}, r"must be a positive duration"),
        ({"rule": "median", "max_offset": "2D"}, r"acts under rule='nearest_end' only"),
        ({"mode": "average_year", "rule": "nearest_end"}, r"use rule='median'"),
        ({"mode": "average_year", "max_offset": "2D"}, r"max_offset is read in mode='history'"),
        ({"mode": "average_year", "min_years": 0}, r"min_years=0 must be a whole number from 1"),
        ({"mode": "average_year", "between": ("2005-01-01", "2001-01-01")}, r"ends before it starts"),
        ({"mode": "average_year", "months": None}, r"needs the twelve months' periods"),
    ],
)
def test_an_option_the_mode_cannot_act_on_is_refused(kwargs, match):
    calendar, months = _year_calendar()
    options = {"months": months, **kwargs}
    with pytest.raises(ValueError, match=match):
        bin_readings(_long([("W1", "2017-10-30", 1.0)]), calendar, **options)


@pytest.mark.parametrize(
    "readings, match",
    [
        (pd.DataFrame({"date": ["2017-10-30", "next week"], "W1": [1.0, 2.0]}),
         r"1 row\(s\) of 'date' are not a date.*row 1"),
        (pd.DataFrame({"date": ["2017-10-30T00:00:00+02:00"], "W1": [1.0]}), r"carries a time zone"),
        (pd.DataFrame({"date": ["October"], "W1": [1.0]}), r"with its year"),
        (pd.DataFrame({"date": ["2017-10-30", "2017-10-30"], "W1": [1.0, 2.0]}),
         r"1 reading\(s\) are given twice at one instant with different values: W1 at 2017-10-30T00:00:00 \(1 and 2\)"),
        (pd.DataFrame({"date": ["2017-10-30"], "name": ["W1"], "head": [1.0], "weight": [2.0]}),
         r"carry a 'weight' column"),
        (pd.DataFrame({"when": ["2017-10-30"], "W1": [1.0]}), r"have no 'date' column"),
        (pd.DataFrame({"date": ["2017-10-30"], "W1": ["dry"]}), r"must be finite numbers"),
    ],
)
def test_readings_that_cannot_be_placed_are_refused(readings, match):
    with pytest.raises(ValueError, match=match):
        read_readings(readings, value_column="head", owner="HeadTargets")


def test_equal_duplicate_readings_collapse_and_dates_may_be_the_index():
    frame = pd.DataFrame({"W1": [1.0, 1.0, np.nan]},
                         index=pd.to_datetime(["2017-10-30", "2017-10-30", "2017-11-30"]))
    long = read_readings(frame, value_column="head", owner="HeadTargets")
    assert long.to_dict("list") == {"date": [pd.Timestamp("2017-10-30")], "name": ["W1"], "value": [1.0]}


def test_a_date_key_on_a_period_keyed_set_names_from_hydrographs():
    """Ledger 499's refusal now names the constructor that bins dates."""

    targets = mf.HeadTargets(
        locations=pd.DataFrame({"name": ["W1"], "cell": [0]}),
        values=pd.DataFrame({"time": pd.to_datetime(["2018-01-31"]), "W1": [5.0]}),
    )
    with pytest.raises(ValueError, match=r"HeadTargets.from_hydrographs\(locations, readings, model=...\)"):
        targets.to_wide()


# -- on MF6 ---------------------------------------------------------------------------


def _square(root: Path, *, tdis, recharge, sto, save="LAST"):
    """A 5 x 5 unconfined aquifer drained along its west edge, recharged per period."""

    vor = _grid()
    context = mf.ModelContext(grid=vor, domain=np.ones((1, NCPL), dtype=int))
    model = mf.gwf("flow", context=context, packages=[
        mf.disv(nlay=1, top=100.0, botm=[0.0], **vor.get_disv_gridprops()),
        mf.ic(strt=50.0),
        mf.npf(k=5.0, icelltype=1),
        mf.sto(sy=0.05, ss=1e-5, iconvert=1, **sto),
        mf.ghb(stress_period_data={0: [[(0, c), 50.0, 100.0] for c in range(0, NCPL, 5)]}),
        mf.rch.array(recharge={p: np.full(NCPL, r) for p, r in recharge.items()}, irch=0, context=context),
        mf.oc(head_filerecord="flow.hds", budget_filerecord="flow.cbc",
              saverecord=[("HEAD", save), ("BUDGET", "LAST")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[tdis])
    shutil.rmtree(root, ignore_errors=True)
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        return mf.Project(root, name="p").prepare_run("r", simulation)


def _execute(run):
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        success, report = run.execute(silent=True)
    assert success, "\n".join(report[-15:])
    return run


def _history_run(root: Path, *, save="ALL", nstp=6, start=TT_START, units="days"):
    """Ten Trails' clock: 12 monthly periods from 2017-09-30, a wet winter and dry summer."""

    rates = {p: 2e-3 * (1.0 + np.cos(2 * np.pi * (p + 9) / 12)) for p in range(12)}
    tdis = mf.tdis(nper=12, perioddata=[(float(n), nstp, 1.0) for n in TT_LENGTHS],
                   time_units=units, start_date_time=start)
    return _square(root, tdis=tdis, recharge=rates, sto={"transient": {0: True}}, save=save)


def _year_run(root: Path, year):
    rates = {m: 2e-3 * (1.0 + np.cos(2 * np.pi * (m - 1) / 12)) for m in range(1, 13)}
    by_month = year.by_month(rates)
    return _square(root, tdis=year.tdis(), recharge=by_month,
                   sto={"steady_state": year.steady_state, "transient": year.transient})


LOCATIONS = pd.DataFrame({"name": ["W1", "W2"], "cell": [12, 7], "layer": [0, 0]})


@pytest.fixture(scope="module")
def models(tmp_path_factory):
    root = tmp_path_factory.mktemp("hydrographs")
    three = mf.AverageYear(3, start="2016-10-01", nstp=1)
    two = mf.AverageYear(2, start="2016-10-01", nstp=1)
    return {
        "history": _execute(_history_run(root / "history")).model("flow"),
        "three": (three, _execute(_year_run(root / "three", three)).model("flow")),
        # Never run: the clock refusal happens before anything is read.
        "two": _year_run(root / "two", two).model("flow"),
    }


def _heads_by_step(model) -> pd.DataFrame:
    """Every saved head of cells 12 and 7, dated start + totim: the model's own hydrograph."""

    hds = flopy.utils.HeadFile(str(Path(model.workspace) / "flow.hds"))
    start = model.calendar.start
    rows = []
    for totim in hds.get_times():
        data = np.asarray(hds.get_data(totim=totim)).reshape(-1)
        rows.append({"date": start + pd.to_timedelta(totim, unit="D"), "W1": data[12], "W2": data[7]})
    return pd.DataFrame(rows)


@pytest.mark.slow
def test_the_models_own_hydrograph_bins_back_to_its_period_end_heads(models, tmp_path):
    """Readings that ARE the run's heads at every time step: nearest_end returns the period-end head.

    This is what makes nearest_end the history default: the simulated value is the
    END-of-period head (ledger 336), and a monthly median of a seasonal signal
    differs from it -- measured here on the run's own heads, against the same file.
    """

    model = models["history"]
    readings = _heads_by_step(model)
    path = tmp_path / "hydrographs.csv"
    readings.to_csv(path, index=False)  # the synthetic WIDE CSV
    nearest = mf.HeadTargets.from_hydrographs(LOCATIONS, path, model=model)
    compared = nearest.compare(model)
    assert len(compared) == 24 and compared["residual"].abs().max() < 1e-6
    assert (nearest.binning["offset_days"].abs() < 1e-6).all()
    median = mf.HeadTargets.from_hydrographs(LOCATIONS, path, model=model, rule="median")
    lagged = median.compare(model)["residual"].abs()
    assert lagged.max() > 0.05, "a monthly median of a seasonal signal is not the period-end head"
    # The CSV read back: a period's end head, from the file by hand.
    hds = flopy.utils.HeadFile(str(Path(model.workspace) / "flow.hds"))
    by_hand = np.asarray(hds.get_data(kstpkper=(5, 3))).reshape(-1)[12]
    w1 = compared.loc[(compared["name"] == "W1") & (compared["per"] == 3)].iloc[0]
    assert w1["head_target"] == pytest.approx(by_hand) and w1["sim_head"] == pytest.approx(by_hand)


@pytest.mark.slow
def test_average_year_targets_land_on_the_last_cycle_and_another_cycle_count_is_refused(models, tmp_path):
    year, model = models["three"]
    days = pd.date_range("2000-01-01", "2015-12-31", freq="10D")
    readings = pd.DataFrame({"date": days, "W1": 50.0 + np.cos(2 * np.pi * days.dayofyear / 365.0)})
    path = tmp_path / "sixteen_years.csv"
    readings.to_csv(path, index=False)  # a synthetic WIDE CSV, dates as text
    targets = mf.HeadTargets.from_hydrographs(LOCATIONS, path, model=model, mode="average_year")
    january = readings.loc[days.month == 1].assign(year=days[days.month == 1].year)
    assert targets.binning.set_index("month").loc[1, "value"] == pytest.approx(
        january.groupby("year")["W1"].median().median()
    )
    last = year.months_of(-1)
    assert sorted(targets.binning["per"]) == sorted(last.values())
    assert set(targets.compare(model)["per"]) == set(last.values())
    assert targets.calendar.same_clock(model.calendar)
    with pytest.raises(ValueError, match=r"HeadTargets.compare\(...\): .*this model runs 25 stress "
                                         r"period\(s\) where the targets were binned on 37.*"
                                         r"HeadTargets.from_hydrographs\(..., model=<this model>\)"):
        targets.compare(models["two"])
    with pytest.raises(ValueError, match=r"simulated_heads"):
        targets.simulated_heads(models["two"])


@pytest.mark.slow
def test_a_model_that_is_not_an_average_year_is_refused_naming_it(models):
    with pytest.raises(ValueError, match=r"this model is not an average-year cycle \(Period 0, the first "
                                         r"month, starts at 2017-09-30.*mf.AverageYear\(...\)"):
        mf.HeadTargets.from_hydrographs(
            LOCATIONS, pd.DataFrame({"date": ["2001-01-05"], "W1": [1.0]}),
            model=models["history"], mode="average_year",
        )


@pytest.mark.slow
def test_a_well_not_among_the_locations_is_refused(models):
    with pytest.raises(ValueError, match=r"name 1 well\(s\) not among the locations .*: w9"):
        mf.HeadTargets.from_hydrographs(
            LOCATIONS, pd.DataFrame({"date": ["2017-10-30"], "W1": [1.0], "W9": [2.0]}),
            model=models["history"],
        )


@pytest.mark.slow
@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({"start": None}, r"declares no START_DATE_TIME.*mf.tdis\(start_date_time="),
        ({"units": "unknown"}, r"declares no TIME_UNITS MF6 reads.*mf.tdis\(time_units='days'"),
    ],
)
def test_a_model_whose_periods_have_no_dates_is_refused(tmp_path, kwargs, match):
    model = _history_run(tmp_path / "undated", save="LAST", nstp=1, **kwargs).model("flow")
    with pytest.raises(ValueError, match=r"HeadTargets.from_hydrographs\(...\) dates the stress periods .*" + match):
        mf.HeadTargets.from_hydrographs(
            LOCATIONS, pd.DataFrame({"date": ["2017-10-30"], "W1": [1.0]}), model=model
        )


@pytest.mark.slow
def test_a_pest_build_reads_the_binned_targets_and_the_review_draws_their_dates(models, tmp_path):
    """The control file holds the binning's values; PESTPP-IES runs; the review is on dates (ies-12)."""

    model = models["history"]
    readings = _heads_by_step(model)
    readings[["W1", "W2"]] += 0.25  # off the model, so a residual exists
    targets = mf.HeadTargets.from_hydrographs(LOCATIONS, readings, model=model)
    cal = model.pest("dated", workspace=tmp_path / "pest")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    cal.observe(targets)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build("dated.pst", noptmax=0)
    obs = cal.pst.observation_data
    measured = obs.loc[obs["weight"] > 0]
    assert len(measured) == 24
    snapshot = pd.read_csv(cal.template_workspace / "hds_target_values.csv")
    by_obs = snapshot.set_index("obsnme")
    report = targets.binning.set_index(["name", "per"])
    for obsnme, row in by_obs.iterrows():
        assert obs.loc[obsnme, "obsval"] == pytest.approx(report.loc[(row["name"], row["per"]), "value"])
    metadata = json.loads((cal.template_workspace / "myflopy_pest_metadata.json").read_text())
    entry = metadata["observation_sets"][0]
    ends = [stamp.isoformat() for stamp in model.calendar.ends]
    assert entry["period_end_dates"] == ends and entry["observation_kind"] == "level"
    assert entry["binning"] == {"mode": "history", "rule": "nearest_end", "max_offset": None,
                                "between": None, "min_years": None}

    review = cal.prior(reals=3)
    assert isinstance(review, IesResults)
    panels = review._location_panels(None, 12)
    assert len(panels) == 2
    for panel in panels:
        assert panel["dates"] is not None
        periods = panel["times"].astype(int)
        assert list(pd.DatetimeIndex(panel["dates"])) == list(model.calendar.ends[periods])
    figure = review.plot_prior_vs_obs()
    assert pd.Timestamp(figure.data[0].x[0]) == model.calendar.ends[0]
    assert figure.layout.xaxis2.title.text == "date (period end)"

    # A template built before the dates were recorded reviews on period numbers.
    path = review.workspace / "myflopy_pest_metadata.json"
    stored = json.loads(path.read_text())
    for set_entry in stored["observation_sets"]:
        set_entry.pop("period_end_dates", None)
    path.write_text(json.dumps(stored))
    old = IesResults(review.workspace, case_name=review.case)
    old_panels = old._location_panels(None, 12)
    assert all(panel["dates"] is None for panel in old_panels)
    assert list(_panel_axis(old_panels[0])) == list(range(12))
    assert old.plot_vs_obs().layout.xaxis2.title.text == "stress period"


@pytest.mark.slow
def test_a_build_on_a_model_of_another_clock_is_refused(models, tmp_path):
    year, three = models["three"]
    targets = mf.HeadTargets.from_hydrographs(
        LOCATIONS, pd.DataFrame({"date": ["2003-01-05", "2004-01-05"], "W1": [1.0, 2.0]}),
        model=three, mode="average_year",
    )
    cal = models["history"].pest("other", workspace=tmp_path / "pest")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    cal.observe(targets)
    with pytest.raises(ValueError, match=r"HeadTargets \(PEST build, prefix 'hds'\): .*runs 12 stress "
                                         r"period\(s\) where the targets were binned on 37"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cal.build("other.pst", noptmax=0)


@pytest.mark.slow
def test_lake_stages_bin_by_the_same_rule(models, tmp_path):
    year = mf.AverageYear(2, start="2016-10-01", nstp=1)
    run, _, _ = _run_lake(tmp_path / "lake", year)
    model = run.model("m")
    days = pd.date_range("2001-01-01", "2010-12-31", freq="7D")
    stages = pd.DataFrame({"date": days, "pit": 95.0 + 0.5 * np.sin(2 * np.pi * days.dayofyear / 365.0)})
    lake = mf.LakeStageTargets.from_hydrographs({"pit": 0}, stages, model=model, mode="average_year")
    direct = bin_readings(
        read_readings(stages, value_column="stage", owner="LakeStageTargets"), model.calendar,
        mode="average_year", months=year.months_of(-1),
    )
    pd.testing.assert_frame_equal(lake.binning, direct.report)
    compared = lake.compare(model)
    assert set(compared["per"]) == set(year.months_of(-1).values())
    assert compared["sim_stage"].notna().all()
    _, three = models["three"]  # 37 periods where the lake run has 25: refused before LAK is read
    for read in (lake.compare, lake.simulated_series):
        with pytest.raises(ValueError, match=r"LakeStageTargets.(compare|simulated_series)\(...\): .*runs 37 "
                                             r"stress period\(s\) where the targets were binned on 25.*"
                                             r"LakeStageTargets.from_hydrographs\(..., model=<this model>\)"):
            read(three)
    with pytest.raises(ValueError, match=r"not among the locations \(matched exactly, case counts\): PIT"):
        mf.LakeStageTargets.from_hydrographs({"pit": 0}, stages.rename(columns={"pit": "PIT"}), model=model)
    history = mf.LakeStageTargets.from_hydrographs(
        {"pit": 0}, pd.DataFrame({"date": ["2017-01-31 18:00"], "stage": [95.1]}), model=model,
    )
    # One series (date and stage only) names the one lake; 31 January 18:00 is nearest
    # January's end, which on this cycle is period 4's.
    assert history.binning[["name", "per", "value"]].values.tolist() == [["pit", 4, 95.1]]

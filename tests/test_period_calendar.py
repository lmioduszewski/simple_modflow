"""The one reading of a model's clock: ``PeriodCalendar`` / ``mf.period_calendar`` / ``model.calendar``.

Ledger 529. MF6 keeps no dates: a TDIS holds a START_DATE_TIME string, a
TIME_UNITS word and the period lengths, and every place that needs a period's date
must put those together the same way. What these tests pin:

- a period's calendar month is its MIDPOINT's -- the Ten Trails TDIS starts
  2017-09-30, so its first period's START is in September while the period is
  October (31 days, ending 2017-10-31);
- MF6's own seconds per unit (a year is 365.25 days) and FloPy's 8-decimal write
  absorbed by ``same_clock`` and the millisecond rounding;
- the steady flags follow MF6's STO rule, measured on MF6: a period BEFORE the first
  STO block is TRANSIENT, and a model with no STO is steady throughout;
- a live model and the same run read back from its files give the same calendar;
- ``model.pest()`` takes its start date from the one reader.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.mf6.period_calendar import (
    PeriodCalendar,
    steady_periods,
    tdis_start_date_time,
)

pytestmark = pytest.mark.filterwarnings("ignore::UserWarning")


def _grid(n=3, d=10.0):
    verts = np.array([[i * d, j * d] for j in range(n + 1) for i in range(n + 1)], dtype=float)
    iverts = [
        [j * (n + 1) + i, (j + 1) * (n + 1) + i, (j + 1) * (n + 1) + i + 1, j * (n + 1) + i + 1]
        for j in range(n)
        for i in range(n)
    ]
    xcyc = np.array([[(i + 0.5) * d, (j + 0.5) * d] for j in range(n) for i in range(n)])
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)


def _run(root: Path, *, sto=None, tdis=None, nper=3):
    """A 3 x 3 unconfined model with recharge and one CHD cell, run by MF6."""

    vor = _grid()
    context = mf.ModelContext(grid=vor, domain=np.ones((1, 9), dtype=int))
    packages = [
        mf.disv(nlay=1, top=10.0, botm=[0.0], **vor.get_disv_gridprops()),
        mf.ic(strt=5.0),
        mf.npf(k=1.0, icelltype=1),
        mf.chd(stress_period_data={0: [[(0, 0), 6.0]]}),
        mf.rch.array(recharge={0: np.full(9, 1e-3)}, irch=0, context=context),
        mf.oc(head_filerecord="flow.hds", budget_filerecord="flow.cbc",
              saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")]),
    ]
    if sto is not None:
        packages.append(mf.sto(sy=0.2, ss=1e-5, iconvert=1, **sto))
    model = mf.gwf("flow", context=context, packages=packages)
    tdis = tdis or mf.tdis(nper=nper, perioddata=[(10.0, 1, 1.0)] * nper, time_units="days",
                           start_date_time="2016-09-30")
    simulation = mf.simulation(model, name="s", packages=[tdis])
    shutil.rmtree(root, ignore_errors=True)
    run = mf.Project(root, name="p").prepare_run("r", simulation)
    success, report = run.execute(silent=True)
    assert success, "\n".join(report[-15:])
    return run


def _storage_flow(model) -> list[float]:
    """Sum of |storage| flow per period, from MF6's own budget file."""

    reader = model._get_budget_reader()
    flows = []
    for kstpkper in reader.get_kstpkper():
        records = reader.get_data(text="STO-SY", kstpkper=kstpkper) + reader.get_data(
            text="STO-SS", kstpkper=kstpkper
        )
        # Storage is written as a full array (imeth 1), not as records.
        flows.append(float(sum(np.abs(np.asarray(r, dtype=float)).sum() for r in records)))
    return flows


# -- the dates -------------------------------------------------------------------


def test_a_period_is_the_month_of_its_midpoint_not_its_start():
    """Ten Trails' TDIS: 2017-09-30, then 31, 30, 31 days -- October, November, December."""

    calendar = PeriodCalendar("2017-09-30", [31.0, 30.0, 31.0], "days")
    table = calendar.table()
    assert list(table["start"].dt.month) == [9, 10, 11]  # the start-date rule: wrong
    assert list(table["month"]) == [10, 11, 12]
    assert list(table["end"]) == [pd.Timestamp("2017-10-31"), pd.Timestamp("2017-11-30"),
                                  pd.Timestamp("2017-12-31")]
    assert table["mid"][0] == pd.Timestamp("2017-10-15 12:00")
    assert table["steady"].isna().all()


def test_units_convert_at_mf6s_own_seconds():
    hours = PeriodCalendar("2020-01-01", [12.0, 36.0], "hours")
    assert list(hours.ends) == [pd.Timestamp("2020-01-01 12:00"), pd.Timestamp("2020-01-03")]
    years = PeriodCalendar("2020-01-01", [1.0], "years")  # MF6: 365.25 days
    assert years.ends[0] == pd.Timestamp("2020-01-01") + pd.Timedelta(days=365.25)
    assert PeriodCalendar("2020-01-01", [1.0], "d").time_units == "days"


def test_flopys_written_lengths_are_the_same_clock():
    """1/24 d is written 0.04166667: the date rounds to the hour, and the clocks agree."""

    exact = PeriodCalendar("2020-01-01", [1 / 24] * 48, "days")
    written = PeriodCalendar("2020-01-01", [0.04166667] * 48, "days")
    assert written.ends[0] == pd.Timestamp("2020-01-01 01:00")
    assert exact.same_clock(written) and written.same_clock(exact)
    assert not exact.same_clock(PeriodCalendar("2020-01-01", [1 / 24 + 1e-6] * 48, "days"))
    assert not exact.same_clock(PeriodCalendar("2020-01-02", [1 / 24] * 48, "days"))
    assert not exact.same_clock(PeriodCalendar("2020-01-01", [1.0] * 48, "hours"))
    with pytest.raises(TypeError, match="two PeriodCalendars"):
        exact.same_clock("2020-01-01")


@pytest.mark.parametrize(
    "kwargs, match",
    [
        (dict(start="October", perlen=[1.0]), "is not one"),
        (dict(start=None, perlen=[1.0]), "takes a date"),
        (dict(start=20200101, perlen=[1.0]), "takes a date"),
        (dict(start="2020-01-01T00:00:00+02:00", perlen=[1.0]), "time zone"),
        (dict(start="2020-01-01", perlen=[]), "holds no stress period"),
        (dict(start="2020-01-01", perlen=[1.0, -1.0]), "period 1 has length -1.0"),
        (dict(start="2020-01-01", perlen=[float("nan")]), "finite number"),
        (dict(start="2020-01-01", perlen=[1.0], time_units=None), "names no unit"),
        (dict(start="2020-01-01", perlen=[1.0], time_units="unknown"), "names no unit"),
        (dict(start="2020-01-01", perlen=[1.0, 2.0], steady=[True]), "1 flags for 2 periods"),
    ],
)
def test_what_cannot_be_a_calendar_is_refused(kwargs, match):
    with pytest.raises(ValueError, match=match):
        PeriodCalendar(**kwargs)


# -- read off a model --------------------------------------------------------------


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    root = tmp_path_factory.mktemp("calendar")
    return {
        "none": _run(root / "none"),
        "late": _run(root / "late", sto=dict(transient={1: True})),
        "first": _run(root / "first", sto=dict(transient={0: True})),
        "steady_first": _run(root / "steady_first", sto=dict(steady_state={0: True}, transient={1: True})),
    }


def test_a_period_before_the_first_sto_block_is_transient(runs):
    """MF6 6.8, measured: STO with TRANSIENT from period 2 alone runs period 1 transient.

    Its storage flows equal those of a STO declaring TRANSIENT from period 1, to the
    bit, where a steady first period carries none. The calendar reads it so.
    """

    late, first, steady = (_storage_flow(runs[k].model("flow")) for k in ("late", "first", "steady_first"))
    assert late == first and late[0] > 0.1
    assert steady[0] == 0.0
    assert runs["late"].model("flow").calendar.steady == (False, False, False)
    assert runs["first"].model("flow").calendar.steady == (False, False, False)
    assert runs["steady_first"].model("flow").calendar.steady == (True, False, False)
    assert runs["none"].model("flow").calendar.steady == (True, True, True)  # no STO: steady


def test_a_live_model_and_its_files_give_the_same_calendar(runs):
    run = runs["steady_first"]
    live = run.model("flow").calendar
    loaded = mf.load_mf6_run(run.workspace).calendar
    assert live == loaded and live.same_clock(loaded)
    assert mf.period_calendar(run.model("flow")) == live
    assert live.start == pd.Timestamp("2016-09-30") and live.time_units == "days"
    assert list(live.ends) == [pd.Timestamp("2016-10-10"), pd.Timestamp("2016-10-20"),
                               pd.Timestamp("2016-10-30")]


def test_reading_the_flags_adds_no_sto_block(runs):
    """Ledger 449: asking FloPy for a period's block creates it. The reader does not ask."""

    model = runs["late"].model("flow")
    sto = model.fl.sto
    before = sorted(sto.transient.get_active_key_dict())
    assert steady_periods(model, 3) == (False, False, False)
    assert sorted(sto.transient.get_active_key_dict()) == before
    assert sorted(sto.steady_state.get_active_key_dict()) == []


def test_a_transport_model_has_no_flags_of_its_own():
    class Gwt:
        model_type = "gwt6"

    assert steady_periods(Gwt(), 4) is None


def test_a_tdis_without_a_start_or_a_unit_is_refused_naming_the_fix(tmp_path):
    no_start = _run(tmp_path / "a", tdis=mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)], time_units="days"))
    with pytest.raises(ValueError, match=r"no START_DATE_TIME.*mf\.tdis\(start_date_time="):
        no_start.model("flow").calendar
    assert tdis_start_date_time(no_start.model("flow")) is None
    no_unit = _run(tmp_path / "b", tdis=mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)],
                                                start_date_time="2020-01-01"))
    with pytest.raises(ValueError, match=r"no TIME_UNITS.*mf\.tdis\(time_units='days'"):
        mf.period_calendar(no_unit.model("flow"))
    assert tdis_start_date_time(no_unit.model("flow")) == "2020-01-01"


def test_pest_takes_its_start_from_the_one_reader(runs):
    model = runs["none"].model("flow")
    assert not hasattr(model, "_tdis_start_datetime")
    assert str(model.pest("calib").start_datetime) == "2016-09-30"

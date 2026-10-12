"""Observation-side hygiene before new target families land (ledgers 499 .. 505).

Every test runs on the real MF6 toy of ``tests/test_pest_handed_targets.py`` -- a
6 x 6 Voronoi grid, one layer, three 10-day periods of three steps, a four-reach
stream, a three-cell DRN zone, a two-cell lake named ``pond``, no MF6 OBS package --
so every family reads its FILE route, keyed by zero-based stress period.

What each test pins, in the order the defects were found:

* a date, or a duration, in a target's key was read as a period numbered in
  nanoseconds (ledger 499);
* ``compare()``/``stats()`` on every FILE route answered ``NaN`` for a key past the
  run's periods and dropped it from ``stats()`` (ledger 500);
"""

from __future__ import annotations

import pandas as pd
import pytest

import myflopy as mf
from tests.test_pest_handed_targets import _build, _head_targets, run_toy

pytest.importorskip("pyemu")


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    """The toy model of ``test_pest_handed_targets``, run."""

    return run_toy(tmp_path_factory.mktemp("target_hygiene"))


# -- ledger 499: a date is not a period ---------------------------------

_DATED = r"reads `time` as a zero-based stress period NUMBER, and {} (is|are) a date or a duration"


def test_a_date_key_is_refused_by_compare_stats_to_wide_and_the_build(toy, tmp_path):
    """pandas reads a datetime as nanoseconds since 1970 -- finite and whole -- so the
    target was keyed by period 1706659200000000000: compare() NaN, stats() n = 0."""

    dated = pd.DataFrame({"time": pd.to_datetime(["2024-01-31", "2024-02-20"]), "W": [95.0, 96.0]})
    heads = _head_targets([("W", 22)], dated)
    for call in (lambda: heads.compare(toy), lambda: heads.stats(toy), heads.to_wide):
        with pytest.raises(ValueError, match=_DATED.format("2024-01-31, 2024-02-20")):
            call()
    with pytest.raises(ValueError, match=_DATED.format("2024-01-31, 2024-02-20")):
        _build(toy, tmp_path / "pest", observe=(heads,))

    # A Series indexed by date reaches the same table as Timestamp OBJECTS.
    series = pd.Series([99.0], index=pd.DatetimeIndex(["2024-01-31"], name="time"), name="pond")
    pond = mf.LakeStageTargets(locations={"pond": 0}, values=series.reset_index().astype(object))
    with pytest.raises(ValueError, match=_DATED.format("2024-01-31")):
        pond.compare(toy)
    # A duration too, which pandas reads as nanoseconds the same way.
    flow = mf.SfrFlowTargets(locations={"a": 0},
                             values=pd.DataFrame({"time": pd.to_timedelta([10], unit="D"), "a": [-7.0]}))
    with pytest.raises(ValueError, match=_DATED.format("10 days")):
        flow.compare(toy)


def test_a_negative_key_is_refused(toy):
    heads = _head_targets([("W", 22)], pd.DataFrame({"time": [0, -1], "W": [95.0, 96.0]}))
    with pytest.raises(ValueError, match=r"and -1 is negative: the first stress period is 0"):
        heads.stats(toy)


# -- ledger 500: every file route refuses a period the run lacks --------

_OFF = (
    r"{caller}: a target is keyed by zero-based stress period and read at its end, and "
    r"the periods {route} holds are 0 \.\. 2; 1 target row\(s\) name none of them: "
    r"'{name}' at time 7"
)


def test_every_file_route_refuses_a_key_past_the_runs_periods(toy):
    """Keys 0 and 7 on a three-period run: compare() gave period 7 a NaN and stats()
    counted one row, with no message -- while the PEST build refused the same key."""

    def values(name, value):
        return pd.DataFrame({"time": [0, 7], name: [value, value]})

    cases = (
        (_head_targets([("W", 22)], values("W", 95.0)), "HeadTargets", "the head file", "W"),
        (mf.LakeStageTargets(locations={"pond": 0}, values=values("pond", 99.0)),
         "LakeStageTargets", "the LAK stage file", "pond"),
        (mf.SfrStageTargets(locations={"g": 1}, values=values("g", 95.0)),
         "SfrStageTargets", "the SFR stage file", "g"),
        (mf.SfrFlowTargets(locations={"a": 0}, values=values("a", -7.0)),
         "SfrFlowTargets", "the SFR budget file", "a"),
    )
    for targets, owner, route, name in cases:
        refusal = _OFF.format(caller=rf"{owner}\.compare\(\.\.\.\)", route=route, name=name)
        with pytest.raises(ValueError, match=refusal):
            targets.compare(toy)
        with pytest.raises(ValueError, match=refusal):
            targets.stats(toy)
    # The periods the run holds still compare.
    assert _head_targets([("W", 22)], pd.DataFrame({"time": [0, 2], "W": [95.0, 96.0]})).stats(toy)["n"].item() == 2

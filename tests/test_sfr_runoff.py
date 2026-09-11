"""SFR runoff is a TOTAL per stream, spread along it by reach length.

MF6's RUNOFF is a volume per time on each reach, and a scalar was written to EVERY
reach -- so a stream's runoff scaled with its reach count, and refining the grid
under a stream multiplied it. Decided 2026-09-11: a number is each stream's total,
``{stream: total}`` gives each named stream its own, and either is split by reach
length. RAINFALL and EVAPORATION stay per-reach rates.
"""

from __future__ import annotations

import pytest

from tests.test_sfr_builder import _builder


def _runoff(builder, period=0):
    return {int(r[0]): float(r[2]) for r in builder.perioddata[period] if r[1] == "RUNOFF"}


def _per_stream(builder, amounts):
    return {s: sum(amounts[r] for r in reaches) for s, reaches in builder.stream_reaches.items()}


def test_a_number_is_each_streams_total():
    builder = _builder(runoff=10.0)
    assert _per_stream(builder, _runoff(builder)) == pytest.approx({"tributary": 10.0, "main": 10.0})


def test_it_is_spread_by_reach_length():
    builder = _builder(runoff=10.0)
    amounts = _runoff(builder)
    for reaches in builder.stream_reaches.values():
        lengths = builder.reaches.loc[reaches, "rlen"].astype(float)
        assert [amounts[r] for r in reaches] == pytest.approx(list(10.0 * lengths / lengths.sum()))


def test_a_mapping_gives_each_named_stream_its_own_total():
    builder = _builder(runoff={"main": 50.0, "tributary": 10.0})
    assert _per_stream(builder, _runoff(builder)) == pytest.approx({"tributary": 10.0, "main": 50.0})


def test_explicit_reach_rows_are_left_alone():
    assert _runoff(_builder(runoff={0: [(0, 3.0)]})) == {0: 3.0}


def test_a_time_series_runoff_is_refused():
    with pytest.raises(ValueError, match="cannot be spread"):
        _builder(runoff="runoff_series").perioddata


def test_an_unknown_stream_is_named():
    with pytest.raises(KeyError, match="the streams are"):
        _builder(runoff={"nowhere": 1.0}).perioddata


def test_rainfall_is_still_a_rate_on_every_reach():
    builder = _builder(rainfall=0.001)
    rows = [r for r in builder.perioddata[0] if r[1] == "RAINFALL"]
    assert len(rows) == len(builder.reaches) and {r[2] for r in rows} == {0.001}

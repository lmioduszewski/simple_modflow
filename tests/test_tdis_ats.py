"""``mf.tdis(ats=...)``: adaptive time stepping on the package-first path.

The resolution and record-building engine already existed and is covered by
``test_ats_discretization.py``; what is pinned here is that the canonical
``mf.tdis`` reaches it, and that the ``MAXATS`` dimension is sized -- FloPy leaves
it at 1 however many records it is given, so before this the records were written
and MODFLOW 6 read only the first.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

import myflopy as mf

PERIODS = [[1.0, 10, 1.1]] * 4


def _records(spec):
    """The ATS records a tdis spec carries, or None."""

    return spec.options.get("ats_perioddata")


def test_no_ats_by_default():
    """A plain tdis writes no ATS records at all."""

    assert _records(mf.tdis(nper=4, perioddata=PERIODS)) is None
    assert _records(mf.tdis(nper=4, perioddata=PERIODS, ats=False)) is None


def test_ats_true_covers_every_period():
    """``ats=True`` builds one record per stress period."""

    records = _records(mf.tdis(nper=4, perioddata=PERIODS, ats=True))
    assert len(records) == 4
    assert [r[0] for r in records] == [0, 1, 2, 3], "iperats stays zero-based"


def test_ats_selects_periods_by_zero_based_index():
    """An iterable picks periods, counted from zero as everywhere in myflopy."""

    records = _records(mf.tdis(nper=5, perioddata=[[1.0, 10, 1.1]] * 5, ats=[1, 3]))
    assert [r[0] for r in records] == [1, 3]


def test_ats_mapping_overrides_one_period():
    """A mapping sets per-period values without disturbing the others."""

    records = _records(
        mf.tdis(nper=4, perioddata=PERIODS, ats={2: {"dtmin": 1.5e-4, "dtmax": 0.4}})
    )
    assert len(records) == 1
    iperats, _dt0, dtmin, dtmax, _dtadj, _dtfailadj = records[0]
    assert iperats == 2
    assert dtmin == pytest.approx(1.5e-4)
    assert dtmax == pytest.approx(0.4)


def test_ats_defaults_are_derived_from_the_period():
    """dt0/dtmin/dtmax come from the period's own perlen/nstp/tsmult when unset."""

    records = _records(mf.tdis(nper=4, perioddata=PERIODS, ats=True))
    _iperats, dt0, dtmin, dtmax, dtadj, dtfailadj = records[0]
    assert 0 < dt0 < 1.0, "dt0 is the first sub-step of the fixed-step period"
    assert dtmin == pytest.approx(1.0 * 1e-5)
    assert dtmax == pytest.approx(1.0), "a whole period when nothing strains the solver"
    assert dtadj == pytest.approx(2.0)
    assert dtfailadj == pytest.approx(5.0)


def test_ats_step_controls_are_forwarded():
    """The growth and failure-retry factors reach the records."""

    records = _records(
        mf.tdis(nper=4, perioddata=PERIODS, ats=True, ats_dtadj=3.0, ats_dtfailadj=10.0)
    )
    assert all(r[4] == pytest.approx(3.0) for r in records)
    assert all(r[5] == pytest.approx(10.0) for r in records)


def test_ats_and_ats_perioddata_together_are_refused():
    """Two ways to say the same thing would silently pick one; refuse instead."""

    with pytest.raises(ValueError, match="not both"):
        mf.tdis(nper=4, perioddata=PERIODS, ats=True,
                ats_perioddata=[[0, 0.1, 1e-5, 1.0, 2.0, 5.0]])


def test_raw_ats_perioddata_still_works():
    """The escape hatch is untouched -- records pass through verbatim."""

    raw = [[0, 0.25, 1e-5, 1.0, 2.0, 5.0]]
    assert _records(mf.tdis(nper=4, perioddata=PERIODS, ats_perioddata=raw)) == raw


def test_an_out_of_range_period_is_refused():
    """Zero-based indexing makes an off-by-one easy; catch it at construction."""

    with pytest.raises(ValueError, match="out of range"):
        mf.tdis(nper=4, perioddata=PERIODS, ats=[4])


def _written_ats(tmp: Path, tdis):
    """Build a minimal simulation with ``tdis`` and return its written .ats text."""

    model = mf.gwf("m", packages=[
        mf.dis(nlay=1, nrow=1, ncol=2, delr=1.0, delc=1.0, top=10.0, botm=[0.0]),
        mf.ic(strt=5.0),
        mf.npf(k=1.0),
        mf.sto(ss=1e-5, sy=0.1, iconvert=1, transient={0: True}),
        mf.oc(head_filerecord="m.hds", saverecord=[("HEAD", "LAST")]),
    ])
    sim = mf.simulation(model, name="s", tdis=tdis,
                        solver=mf.ims(models=["m"], complexity="simple"))
    built = sim.build(tmp)
    flopy_sim = built.sim if hasattr(built, "sim") else built
    for attr in ("simulation", "_simulation", "mfsim"):
        flopy_sim = getattr(flopy_sim, attr, flopy_sim)
    flopy_sim.write_simulation(silent=True)
    files = list(tmp.glob("*.ats"))
    return files[0].read_text() if files else ""


def test_maxats_is_sized_to_the_record_count():
    """FloPy leaves MAXATS at 1, so MODFLOW 6 would read only the first record."""

    with tempfile.TemporaryDirectory() as directory:
        text = _written_ats(Path(directory), mf.tdis(nper=4, perioddata=PERIODS, ats=True))
    assert text, "no .ats file was written"
    line = next(row for row in text.splitlines() if "MAXATS" in row.upper())
    assert line.split()[-1] == "4", f"MAXATS not sized to the records: {line!r}"


def test_maxats_is_sized_for_raw_records_too():
    """The escape hatch had the same defect; it is fixed at the build, not in tdis."""

    raw = [[0, 0.25, 1e-5, 1.0, 2.0, 5.0], [2, 0.25, 1e-5, 1.0, 2.0, 5.0]]
    with tempfile.TemporaryDirectory() as directory:
        text = _written_ats(
            Path(directory), mf.tdis(nper=4, perioddata=PERIODS, ats_perioddata=raw)
        )
    line = next(row for row in text.splitlines() if "MAXATS" in row.upper())
    assert line.split()[-1] == "2"

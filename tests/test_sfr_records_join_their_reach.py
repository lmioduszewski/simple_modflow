"""An SFR exchange record is joined to its OWN reach, never to its cell (ledger 470).

MODFLOW 6 writes each reach's stream-aquifer exchange to the GWF budget as an ``SFR``
record whose ``node`` is the cell and whose ``node2`` is the reach. The view joined
the records to the reach table on the CELL, so a cell of ``k`` reaches carried
``k * k`` rows, each reach tagged with every flow in its cell: on the testing canonical
model cell 195 holds reaches 22, 37 and 38, and the profile gave each of them the
cell's -493.56 where MF6's own SFR budget gives -20.91, -471.96 and -0.69. These pin
the view against MODFLOW 6's SFR package budget file (its ``GWF`` term: ``node`` the
reach, ``q`` reach-referenced, so the opposite sign).
"""

from __future__ import annotations

from pathlib import Path

import flopy
import numpy as np
import pytest

pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def sfr_budget_gwf(canonical_run):
    """MF6's own per-reach exchange: ``{(kstp, kper): {reach: q_gwf}}``, aquifer-referenced."""

    name = canonical_run.fl.sfr.budget_filerecord.get_data()[0][0]
    path = Path(canonical_run.workspace) / name
    budget = flopy.utils.CellBudgetFile(str(path))
    out = {}
    for kstpkper in budget.get_kstpkper():
        record = budget.get_data(text="GWF", kstpkper=kstpkper)[0]
        # The SFR file's GWF term is REACH-referenced (positive into the reach); the
        # view keeps the GWF budget's aquifer-referenced sign, the opposite one.
        out[kstpkper] = {int(n) - 1: -float(q) for n, q in zip(record["node"], record["q"], strict=True)}
    return out


def _multi_reach_cells(frame):
    per_cell = frame.groupby(["layer", "cell"])["reach"].nunique()
    return per_cell[per_cell > 1]


def test_every_record_is_one_row_with_its_own_reach_flow(canonical_run, sfr_budget_gwf):
    """One row per record, and each reach reads MF6's own exchange for it."""

    frame = canonical_run.sfr.results.q.get()
    assert not _multi_reach_cells(frame).empty, "the fixture must hold a cell of several reaches"
    for kstpkper, flows in sfr_budget_gwf.items():
        at = frame[frame["kstpkper"] == kstpkper]
        assert len(at) == len(flows), kstpkper
        assert at["reach"].is_unique
        got = dict(zip(at["reach"].astype(int), at["q_gwf"].astype(float), strict=True))
        np.testing.assert_allclose(
            [got[reach] for reach in sorted(flows)], [flows[reach] for reach in sorted(flows)],
            rtol=1e-6, atol=1e-6,
        )


def test_the_profile_gives_each_reach_its_own_flow(canonical_run, sfr_budget_gwf):
    """``sfr.results.profile`` reads each reach's flow, not its cell's total."""

    last = max(sfr_budget_gwf)
    profile = canonical_run.sfr.results.profile.get(per=last[1]).set_index("reach")
    flows = sfr_budget_gwf[last]
    np.testing.assert_allclose(
        profile.loc[sorted(flows), "q_gwf"].to_numpy(float),
        [flows[reach] for reach in sorted(flows)], rtol=1e-6, atol=1e-6,
    )


def test_surface_water_sums_each_reach_once(canonical_run, sfr_budget_gwf):
    """A multi-reach cell's surface-water ``q``/``rlen`` are its reaches' sums, once each."""

    frame = canonical_run.sfr.results.q.get()
    last = max(sfr_budget_gwf)
    layer, cell = _multi_reach_cells(frame).index[0]
    reaches = frame.loc[(frame["layer"] == layer) & (frame["cell"] == cell), "reach"].unique()
    flows = sfr_budget_gwf[last]
    table = canonical_run.surface_water.results.q.get(per=last[1])
    row = table[(table["source"] == "sfr") & (table["layer"] == layer) & (table["cell"] == cell)]
    assert len(row) == 1
    assert row["q"].iloc[0] == pytest.approx(sum(flows[int(r)] for r in reaches), rel=1e-6)
    rlen = canonical_run.sfr.results.profile.get(per=last[1]).set_index("reach")["rlen"]
    assert row["rlen"].iloc[0] == pytest.approx(float(rlen.loc[[int(r) for r in reaches]].sum()))

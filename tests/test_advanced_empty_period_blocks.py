"""``mf.sfr``/``mf.lak`` write no EMPTY period block (ledger 468).

The builders gave every period an entry, ``[]`` where it sets nothing, and FloPy wrote
each as an empty block. MF6 keeps an advanced package's settings until a block changes
them, so the empty block changed nothing -- measured byte-identical below -- while
FloPy 3.10 reads the first period's records for its dtype, so a package whose period 0
set nothing made ``gwf.plot()`` and the package's ``.plot()`` raise ``AttributeError``.
Cases and the run come from ``docs/reproducers/advanced_empty_period_blocks.py``.
"""

from __future__ import annotations

import filecmp
import importlib.util
from pathlib import Path

import pytest

_PATH = Path(__file__).resolve().parents[1] / "docs" / "reproducers" / "advanced_empty_period_blocks.py"
_SPEC = importlib.util.spec_from_file_location("advanced_empty_period_blocks", _PATH)
repro = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(repro)


def test_a_period_that_sets_nothing_has_no_entry():
    sfr = repro._sfr_model(False, {1: {"s0": 50.0}})[0][-1]
    assert sfr.options["perioddata"] == {1: [[0, "INFLOW", 50.0]]}
    assert repro._sfr_model(False, None)[0][-1].options["perioddata"] == {}
    lak = repro._lak_model(False, None)[0][-1]
    assert list(lak.options["perioddata"]) == [1]


@pytest.mark.parametrize(
    ("factory", "inflow"),
    [(repro._sfr_model, {1: {"s0": 50.0}}), (repro._sfr_model, None), (repro._lak_model, None)],
    ids=["sfr-period-1", "sfr-none", "lak-period-1"],
)
def test_the_run_is_the_same_and_flopy_can_plot_it(tmp_path, factory, inflow):
    without, plots, blocks, kind = repro.run(factory, False, inflow, tmp_path)
    assert set(plots.values()) == {"ok"}, plots
    with_empty, plots_empty, blocks_empty, _ = repro.run(factory, True, inflow, tmp_path)
    assert blocks_empty == 3 and blocks < 3
    assert "AttributeError" in plots_empty["gwf.plot()"], "the defect the omission avoids"
    for name in ("m.cbc", f"m.{kind}.cbc", "m.hds", f"m.{kind}.stg"):
        assert filecmp.cmp(without / name, with_empty / name, shallow=False), name

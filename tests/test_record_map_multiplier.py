"""A record map's ``multiplier`` scales each record before ``agg`` and never the fill (ledger 473).

It acted AFTER ``agg`` and on the ``fill_value`` cells, so ``fill_value=1,
multiplier=2`` drew 2.0 where there is no record -- and 1.0 when the selection held no
record at all -- and under ``multiplier=-1`` ``agg="max"`` drew minus the largest
record rather than the largest negated one. Every record map shares the payload
builders, so a cell-results noun, the SFR/LAK exchange maps and a group comparison
are each checked once here; the ORDER needs a cell of two records and is checked in
``tests/test_input_fields_combine_by_kind.py``.
"""

from __future__ import annotations

import numpy as np
import pytest

pytestmark = pytest.mark.slow


def _z(picture, ncpl):
    for trace in picture.fig.data:
        values = getattr(trace, "z", None)
        if values is not None and len(values) == ncpl:
            return np.asarray(values, dtype=float)
    raise AssertionError("no per-cell trace")


def test_a_cell_with_no_record_keeps_the_fill(canonical_run):
    ncpl = canonical_run.vor.ncpl
    q = canonical_run.drn.results.q
    table = q.get(per=0, layer=0)
    holds = set(table["cell"].astype(int))
    drawn = _z(q.map(multiplier=2.0, fill_value=1.0), ncpl)
    empty = [cell for cell in range(ncpl) if cell not in holds]
    assert set(drawn[empty].tolist()) == {1.0}
    plain = _z(q.map(), ncpl)
    held = sorted(holds)
    assert drawn[held] == pytest.approx(2.0 * plain[held])


@pytest.mark.parametrize("package", ["sfr", "lak", "surface_water"])
def test_the_exchange_maps_keep_the_fill(canonical_run, package):
    ncpl = canonical_run.vor.ncpl
    view = getattr(canonical_run, package).results.q
    plain = _z(view.map(fill_value=1.0), ncpl)
    scaled = _z(view.map(multiplier=3.0, fill_value=1.0), ncpl)
    filled = plain == 1.0
    assert filled.any()
    assert set(scaled[filled].tolist()) == {1.0}
    assert scaled[~filled] == pytest.approx(3.0 * plain[~filled], nan_ok=True)


def test_a_group_comparison_keeps_the_fill(canonical_run):
    import myflopy as mf
    from myflopy.project.group.core import ModelGroup

    group = ModelGroup({"a": canonical_run, "b": mf.load_mf6_run(canonical_run.workspace)},
                       reference="a")
    ncpl = canonical_run.vor.ncpl
    drawn = _z(group.packages.drn.results.q.compare_map(model_name="b", multiplier=2.0,
                                                        fill_value=1.0), ncpl)
    holds = set(canonical_run.drn.results.q.get(per=0, layer=0)["cell"].astype(int))
    empty = [cell for cell in range(ncpl) if cell not in holds]
    assert set(drawn[empty].tolist()) == {1.0}

"""A budget record carrying its value as an auxiliary reads ONE value column (ledger 471).

NPF's ``DATA-SAT`` writes each cell's saturation as the auxiliary ``sat`` beside a
placeholder flow of 0.0. The reader renamed that placeholder ``sat`` too, so the table
held TWO ``sat`` columns and ``model/group/diff ... uzf.results.sat`` ``map``,
``compare`` and ``compare_map`` raised. The value is the auxiliary; the placeholder goes.
"""

from __future__ import annotations

import numpy as np
import pytest

pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def pair(canonical_run):
    import myflopy as mf
    from myflopy.project.group.core import ModelGroup

    reloaded = mf.load_mf6_run(canonical_run.workspace)
    return ModelGroup({"a": canonical_run, "b": reloaded}, reference="a")


def test_sat_is_one_column_holding_mf6s_saturation(canonical_run):
    frame = canonical_run.uzf.results.sat.get()
    assert list(frame.columns).count("sat") == 1
    record = canonical_run._get_budget_reader().get_data(text="DATA-SAT")[-1]
    last = frame[frame["kstpkper"] == frame["kstpkper"].iloc[-1]].sort_values("node")
    np.testing.assert_allclose(last["sat"].to_numpy(float), record["sat"])
    assert record["sat"].max() > 0.0  # the placeholder flow is 0.0 everywhere


def test_sat_maps_on_one_model(canonical_run):
    picture = canonical_run.uzf.results.sat.map()
    drawn = np.asarray(picture.fig.data[0].z, dtype=float)
    assert np.nanmax(drawn) > 0.0


def test_sat_compares_and_maps_on_a_group_and_a_diff(pair):
    view = pair.packages.uzf.results.sat
    compared = view.compare()
    assert not compared.empty
    view.map().fig
    view.compare_map().fig
    diff = pair.diff().packages.uzf.results.sat
    assert not diff.get().empty
    diff.map().fig

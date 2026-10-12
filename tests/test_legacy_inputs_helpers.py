"""``model.inputs.rch.map`` and ``model.inputs.uzf.finf`` (ledger 476).

``rch.map(multiplier=2)`` doubled each boundname ("zonezone") and ``2.5`` raised;
``uzf.finf`` raised ``ValueError`` on a stacked unsaturated column, whose lower
object MF6 applies no infiltration to.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, build_canonical_model


@pytest.fixture
def unrun(tmp_path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return build_canonical_model(tmp_path / "canonical", config=CanonicalModelConfig.testing())


def _with_boundnames(model):
    rch = model.package("rch")
    rch.boundnames.set_data(True)
    rch.stress_period_data.set_data({
        period: [(cell, rate, "zone") for cell, rate in zip(
            records["cellid"], records["recharge"], strict=True)]
        for period, records in rch.stress_period_data.get_data().items()
    })


@pytest.mark.parametrize("multiplier", [2, 2.5])
def test_rch_map_scales_the_numbers_and_never_a_boundname(unrun, multiplier):
    plain = unrun.inputs.rch.map(0)._custom_hover["recharge"]
    _with_boundnames(unrun)
    hover = unrun.inputs.rch.map(0, multiplier=multiplier)._custom_hover
    named = [name for name in hover["boundname"] if name]
    assert named and set(named) == {"zone"}
    assert np.allclose(hover["recharge"], np.asarray(plain) * multiplier)


def test_uzf_finf_reads_a_stacked_columns_land_surface_object(unrun):
    uzf = unrun.package("uzf")
    before = unrun.inputs.uzf.finf(0)
    packagedata = uzf.packagedata.get_data().copy()
    layer, top_cell = packagedata["cellid"][0]
    moved_from = int(packagedata["cellid"][1][1])
    packagedata["cellid"][1] = (layer + 1, top_cell)  # object 1 now under object 0
    packagedata["landflag"][1] = 0
    uzf.packagedata.set_data(packagedata)
    after = unrun.inputs.uzf.finf(0)
    assert after.loc[top_cell, "finf"] == before.loc[top_cell, "finf"]
    assert before.loc[moved_from, "finf"] > 0
    assert after.loc[moved_from, "finf"] == 0
    packagedata["landflag"][1] = 1  # two land-surface objects in one column
    uzf.packagedata.set_data(packagedata)
    with pytest.raises(ValueError, match="more than one land-surface UZF object"):
        unrun.inputs.uzf.finf(0)

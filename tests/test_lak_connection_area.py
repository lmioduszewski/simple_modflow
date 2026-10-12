"""A LAK connection's ``connection_area`` reads the connection as MODFLOW 6 does (ledger 474).

``build_lak_connection_table`` computed ``connwidth * (telev - belev)`` for every
connection but a vertical one -- embedded ones included, whose ``connwidth``,
``belev`` and ``telev`` MF6 never reads (their exchange area comes from the lake
TABLE) -- and read a HORIZONTAL connection declared with ``belev == telev`` as zero
area, where MF6 resets the two to the cell's bottom and top (the rule
``package_tables.lake_bottoms`` reads, ledger 351). Pinned on a loaded copy of the
canonical run whose LAK connections are rewritten in memory: nothing is run, the
table reads the package alone.
"""

from __future__ import annotations

import shutil

import numpy as np
import pytest

pytestmark = pytest.mark.slow


@pytest.fixture
def rewired(canonical_run, tmp_path):
    import myflopy as mf

    workspace = tmp_path / "run"
    shutil.copytree(canonical_run.workspace, workspace)
    model = mf.load_mf6_run(workspace)
    data = model.fl.lak.connectiondata.get_data().copy()
    kinds = np.char.lower(data["claktype"].astype(str))
    horizontal = int(np.flatnonzero(kinds == "horizontal")[0])
    embedded = int(np.flatnonzero(kinds == "horizontal")[1])
    data["telev"][horizontal] = data["belev"][horizontal]
    data["claktype"][embedded] = "embeddedh"
    model.fl.lak.connectiondata.set_data(data)
    return model, data, horizontal, embedded


def test_connection_area_follows_mf6s_reading(rewired):
    from myflopy.modflow.mf6.package_tables import (
        _cell_top_and_bottom,
        build_lak_connection_table,
    )

    model, data, horizontal, embedded = rewired
    table = build_lak_connection_table(model)
    top, bottom = _cell_top_and_bottom(model)(tuple(data["cellid"][horizontal]))
    assert top > bottom
    assert table["connection_area"].iloc[horizontal] == pytest.approx(
        float(data["connwidth"][horizontal]) * (top - bottom)
    )
    vertical = table["claktype"].eq("VERTICAL").to_numpy()
    areas = np.asarray(model.vor.area_list, dtype=float).reshape(-1)
    assert table.loc[vertical, "connection_area"].to_numpy() == pytest.approx(
        areas[table.loc[vertical, "cell"].to_numpy(int)]
    )
    others = ~vertical & ~np.isin(np.arange(len(table)), [horizontal, embedded])
    expected = data["connwidth"][others] * np.clip(data["telev"][others] - data["belev"][others], 0, None)
    assert table.loc[others, "connection_area"].to_numpy() == pytest.approx(expected)


def test_an_embedded_connection_has_no_area(rewired):
    from myflopy.modflow.mf6.package_tables import build_lak_connection_table

    model, _, _, embedded = rewired
    table = build_lak_connection_table(model)
    assert table["claktype"].iloc[embedded] == "EMBEDDEDH"
    assert np.isnan(table["connection_area"].iloc[embedded])

"""Every package of a type is in its table, each row naming its package.

Two packages of one type -- two CHDs, a list RCH beside an array one -- gave a table
of 0 rows, silently: FloPy's ``get_package`` returned a LIST, which the table
builder read as an array package with no fields. With explicit names it returned
the first package only, and an array RCH (FloPy's ``rcha``) was never asked for:
measured, a third of the applied recharge left out. Decided 2026-09-11: the table
holds them all, with a ``pname`` column (ledger 229).
"""

from __future__ import annotations

import numpy as np
import pytest

import myflopy as mf
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.mf6.package_tables import build_cell_package_input_table


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    """A list RCH named for its type beside an array one -- myflopy's usual case."""

    n, d = 3, 10.0
    verts = np.array([[i * d, j * d] for j in range(n + 1) for i in range(n + 1)], dtype=float)
    iverts = [
        [j * (n + 1) + i, (j + 1) * (n + 1) + i, (j + 1) * (n + 1) + i + 1, j * (n + 1) + i + 1]
        for j in range(n)
        for i in range(n)
    ]
    xcyc = np.array([[(i + 0.5) * d, (j + 0.5) * d] for j in range(n) for i in range(n)])
    vor = VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)
    context = mf.ModelContext(grid=vor, domain=np.ones((1, n * n), dtype=int))
    model = mf.gwf(
        "flow",
        context=context,
        packages=[
            mf.disv(nlay=1, top=10.0, botm=[0.0], **vor.get_disv_gridprops()),
            mf.ic(strt=5.0),
            mf.npf(k=1.0),
            mf.chd(stress_period_data={0: [[(0, 0), 6.0]]}, name="chd"),
            mf.chd(stress_period_data={0: [[(0, 8), 4.0]]}, name="chd_east"),
            mf.rch(stress_period_data={0: [[(0, 4), 1e-4]]}, name="rch"),
            mf.rch.array(recharge={0: np.full(9, 5e-5)}, name="rch_array", context=context),
            mf.oc(),
        ],
    )
    run = mf.Project(tmp_path_factory.mktemp("same"), name="same").prepare_run(
        "r", mf.simulation(model, name="same")
    )
    success, report = run.execute(silent=True)
    assert success, report
    return run


def test_both_chds_are_in_the_head_table(run):
    table = run.model("flow").packages.chd.inputs.head.get()
    assert sorted(zip(table["pname"], table["cell"], table["head"])) == [
        ("CHD", 0, 6.0), ("CHD_EAST", 8, 4.0),
    ]


def test_every_recharge_package_is_in_the_table_array_form_included(run):
    table = run.model("flow").packages.rch.inputs.recharge.get()
    by_package = table.groupby("pname")["recharge"].agg(["count", "sum"])
    assert by_package.loc["RCH", "count"] == 1
    assert by_package.loc["RCH_ARRAY", "count"] == 9
    assert table["recharge"].sum() == pytest.approx(1e-4 + 9 * 5e-5)


def test_a_reopened_run_gives_the_same_table(run):
    table = mf.load_mf6_run(run.workspace).packages.rch.inputs.recharge.get()
    assert set(table["pname"]) == {"RCH", "RCH_ARRAY"}
    assert table["recharge"].sum() == pytest.approx(1e-4 + 9 * 5e-5)


def test_one_package_can_still_be_asked_for_by_its_own_name(run):
    table = build_cell_package_input_table(run.model("flow"), "rch_array")
    assert set(table["pname"]) == {"RCH_ARRAY"}
    assert len(table) == 9


def test_a_map_adds_the_packages_that_share_a_cell(run):
    """MF6 applies both recharge packages to cell 4, and the map sums per cell."""

    z = run.model("flow").packages.rch.inputs.recharge.map(per=0).fig.data[0].z
    assert z[4] == pytest.approx(1e-4 + 5e-5)

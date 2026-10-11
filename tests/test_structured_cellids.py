"""A model's cellids name their cells by ONE rule in every table, diff and region
(ledger 419).

``CellNumbering`` reads a FloPy cellid as the ``(layer, cell)`` the model's views
key by: DIS ``(layer, row, col)`` is ``(layer, row * ncol + col)``, DISV ``(layer,
cell)`` is itself, DISU ``(node,)`` is ``(0, node)`` -- the numbering the model's
own heads table and cell polygons use. ``split_cellid_columns`` read every tuple as
``(layer, cell)``, so a DIS CHD table named its ROWS as its cells and two identical
DIS runs diffed as changed.

Every model here is a real MF6 run, written and run by FloPy in a temporary
directory and opened with ``load_mf6_run``.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import flopy
import numpy as np
import pandas as pd
import pytest
from flopy.utils.gridutil import get_disu_kwargs

import myflopy as mf

# The module, not its new names: docs/reproducers/structured_cellids.py imports this
# file's run builders against the tree before the fix too.
from myflopy.modflow.mf6 import package_explorer_utils as cellids
from myflopy.modflow.mf6.simulation.regions import RegionRegistry, normalize_region_inputs

NLAY, NROW, NCOL = 2, 3, 4
#: Two nearly sealed barriers across the column-1 | column-2 face: row 1 of layer 0
#: (cells 5 | 6) and row 2 of layer 1 (cells 9 | 10). Read as (layer, col) they named
#: cells 1 | 2 of each layer -- row 0's face, which no barrier crosses.
BARRIERS = [((0, 1, 1), (0, 1, 2)), ((1, 2, 1), (1, 2, 2))]


def _write_and_run(sim) -> None:
    sim.write_simulation(silent=True)
    success, report = sim.run_simulation(silent=True)
    assert success, "\n".join(report[-20:])


def _dis_run(
    workspace: Path,
    *,
    changed: dict | None = None,
    uzf: bool = False,
    hfb: bool = False,
    rch: bool = False,
) -> Path:
    """A 2 x 3 x 4 DIS model: CHD 5.0 on the west column, 1.0 on the east one.

    ``changed`` maps a ``(layer, row, col)`` to the head that one CHD record holds
    instead. ``uzf=True`` adds a UZF object on every top-layer cell between the CHD
    columns. ``hfb=True`` adds the barriers of :data:`BARRIERS`. ``rch=True`` adds
    a list RCH of 0.01 on every top-layer cell.
    """

    workspace.mkdir(parents=True, exist_ok=True)
    sim = flopy.mf6.MFSimulation(sim_name="sim", sim_ws=str(workspace), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, complexity="COMPLEX")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m", save_flows=True, newtonoptions="NEWTON")
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=NLAY, nrow=NROW, ncol=NCOL, delr=10.0, delc=10.0, top=10.0,
        botm=[0.0, -10.0],
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0, icelltype=[1, 0])
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    heads = {}
    for layer in range(NLAY):
        for row in range(NROW):
            heads[(layer, row, 0)] = 5.0
            heads[(layer, row, NCOL - 1)] = 1.0
    heads.update(changed or {})
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=list(heads.items()))
    if hfb:
        flopy.mf6.ModflowGwfhfb(
            gwf, stress_period_data={0: [(*pair, 1e-6) for pair in BARRIERS]}
        )
    if rch:
        flopy.mf6.ModflowGwfrch(
            gwf,
            stress_period_data=[((0, row, col), 0.01) for row in range(NROW) for col in range(NCOL)],
        )
    if uzf:
        cells = [(0, row, col) for row in range(NROW) for col in range(1, NCOL - 1)]
        packagedata = [
            [index, cellid, 1, -1, 0.5, 0.1, 0.05, 0.3, 0.2, 3.5]
            for index, cellid in enumerate(cells)
        ]
        flopy.mf6.ModflowGwfuzf(
            gwf,
            simulate_et=False,
            nuzfcells=len(cells),
            ntrailwaves=7,
            nwavesets=40,
            packagedata=packagedata,
            perioddata={0: [[index, 0.001, 0, 0, 0, 0, 0, 0] for index in range(len(cells))]},
        )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord="m.hds",
        budget_filerecord="m.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    _write_and_run(sim)
    return workspace


def _disu_run(workspace: Path) -> Path:
    """The same 2 x 3 x 4 layout written as a DISU: CHD on node ``(node,)`` cellids."""

    workspace.mkdir(parents=True, exist_ok=True)
    kwargs = get_disu_kwargs(
        NLAY, NROW, NCOL, np.full(NCOL, 10.0), np.full(NROW, 10.0), 10.0, [0.0, -10.0],
        return_vertices=True,
    )
    sim = flopy.mf6.MFSimulation(sim_name="sim", sim_ws=str(workspace), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m", save_flows=True)
    flopy.mf6.ModflowGwfdisu(gwf, **kwargs)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    records = []
    for layer in range(NLAY):
        for row in range(NROW):
            first = layer * NROW * NCOL + row * NCOL
            records.append(((first,), 5.0))
            records.append(((first + NCOL - 1,), 1.0))
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=records)
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord="m.hds",
        budget_filerecord="m.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    _write_and_run(sim)
    return workspace


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    root = tmp_path_factory.mktemp("structured_cellids")
    return SimpleNamespace(
        base=_dis_run(root / "base"),
        twin=_dis_run(root / "twin"),
        one_cell=_dis_run(root / "one_cell", changed={(1, 2, 3): 2.0}),
        disu=_disu_run(root / "disu"),
    )


# --- package input tables ---------------------------------------------------


def test_a_dis_package_table_names_the_grids_own_cells(runs):
    model = mf.load_mf6_run(runs.base)
    table = model.chd.inputs.get()

    west = [row * NCOL for row in range(NROW)]
    east = [row * NCOL + NCOL - 1 for row in range(NROW)]
    for layer in range(NLAY):
        rows = table[table["layer"] == layer]
        assert sorted(rows.loc[rows["head"] == 5.0, "cell"]) == west
        assert sorted(rows.loc[rows["head"] == 1.0, "cell"]) == east

    # The cells the table names hold the heads the CHD fixes, in the heads table.
    heads = model.hds.get(per=0)
    joined = table.merge(heads, on=["layer", "cell"], suffixes=("_chd", ""))
    assert len(joined) == len(table) == 2 * NLAY * NROW
    assert joined["head"].to_numpy() == pytest.approx(joined["head_chd"].to_numpy(), abs=1e-6)


def test_a_self_diff_of_two_identical_dis_runs_is_empty(runs):
    base, twin = mf.load_mf6_run(runs.base), mf.load_mf6_run(runs.twin)
    summary = mf.ModelGroup({"base": base, "twin": twin}).diff().summary()

    row = summary.loc[summary["package"] == "chd"].iloc[0]
    assert row["value_cells_changed"] == 0
    assert bool(row["identical"])


def test_a_one_cell_change_diffs_to_that_one_cell(runs):
    base, other = mf.load_mf6_run(runs.base), mf.load_mf6_run(runs.one_cell)
    diff = mf.ModelGroup({"base": base, "other": other}).diff()

    row = diff.summary().loc[lambda frame: frame["package"] == "chd"].iloc[0]
    assert row["value_cells_changed"] == 1
    values = diff.packages.chd.inputs.values()
    assert len(values) == 2 * NLAY * NROW  # one row per CHD cell: no cartesian join
    changed = values[values["head_diff"] != 0.0]
    assert changed[["layer", "cell", "head", "reference_head"]].values.tolist() == [
        [1, 2 * NCOL + 3, 2.0, 1.0]
    ]


def test_a_cellid_off_the_structured_grid_is_refused_not_aliased():
    numbering = cellids.CellNumbering("structured", nlay=NLAY, ncpl=NROW * NCOL, ncol=NCOL)
    frame = pd.DataFrame({"cellid": [(0, 0, 0), (0, 0, NCOL)], "head": [1.0, 2.0]})
    # Column NCOL of row 0 would be cell NCOL -- row 1's first cell -- if not refused.
    with pytest.raises(ValueError, match=r"off the grid.*\(layer, row, col\)"):
        cellids.split_cellid_columns(frame, numbering)
    with pytest.raises(ValueError, match=r"2 item\(s\).*\(layer, row, col\)"):
        cellids.split_cellid_columns(pd.DataFrame({"cellid": [(0, 5)]}), numbering)
    assert numbering.layer_cell((1, 2, 3)) == (1, 2 * NCOL + 3)


def test_a_reach_connected_to_no_cell_splits_to_no_cell():
    # MF6 writes it NONE; FloPy reads it back as the string 'none'. int('none')
    # raised before, and a one-item tuple of it would be refused as a short cellid.
    numbering = cellids.CellNumbering("structured", nlay=NLAY, ncpl=NROW * NCOL, ncol=NCOL)
    frame = pd.DataFrame({"cellid": [(0, 1, 2), "none", None]})
    split = cellids.split_cellid_columns(frame, numbering)
    assert split["layer"].tolist()[0] == 0 and split["cell"].tolist()[0] == NCOL + 2
    assert split["cell"].isna().tolist() == [False, True, True]
    assert split["layer"].isna().tolist() == [False, True, True]


# --- DISU: (node,) is cell `node` of layer 0, as its heads table numbers it --


def test_a_disu_package_table_and_budget_number_a_node_as_its_heads_table_does(runs):
    model = mf.load_mf6_run(runs.disu)
    assert cellids.cell_numbering(model) == cellids.CellNumbering(
        "unstructured", nlay=1, ncpl=NLAY * NROW * NCOL
    )

    heads = model.hds.get(per=0)
    assert set(heads["layer"]) == {0}
    assert sorted(heads["cell"]) == list(range(NLAY * NROW * NCOL))

    table = model.chd.inputs.get()
    nodes = sorted(
        layer * NROW * NCOL + row * NCOL + col
        for layer in range(NLAY) for row in range(NROW) for col in (0, NCOL - 1)
    )
    assert set(table["layer"]) == {0}
    assert sorted(table["cell"]) == nodes
    joined = table.merge(heads, on=["layer", "cell"], suffixes=("_chd", ""))
    assert joined["head"].to_numpy() == pytest.approx(joined["head_chd"].to_numpy(), abs=1e-6)

    # The cell-budget tables read MF6's node through the same rule: FloPy infers two
    # layers from IHC, and node 12 read (1, 0) here beside (0, 12) in the heads table.
    budget = model.chd.results.q.get()
    assert set(budget["layer"]) == {0}
    assert sorted(budget["cell"]) == nodes
    assert model.node_to_lni[12] == (0, 12)


def test_a_disu_region_reads_a_node_and_refuses_a_layer_cell_pair(runs):
    model = mf.load_mf6_run(runs.disu)
    region = model.add_region("r", cellids=[(13,), 14])
    assert region.cells == [13, 14]
    heads = model.region_heads("r")
    assert sorted(heads["cell"]) == [13, 14]

    with pytest.raises(ValueError, match=r"\(node,\)"):
        model.add_region("pair", cellids=[(1, 1)])
    with pytest.raises(ValueError, match="off the grid"):
        model.add_region("past", cellids=[(NLAY * NROW * NCOL,)])
    with pytest.raises(ValueError, match="one layer, 0"):
        model.add_region("layer", cellids=[13], layer=1)
    assert not {"pair", "past", "layer"} & set(model.regions.list_regions())


# --- regions bound every layer and cell when the grid is known ---------------


def test_a_dis_region_refuses_a_layer_or_cell_off_the_grid(runs):
    model = mf.load_mf6_run(runs.base)
    with pytest.raises(ValueError, match=r"off the grid.*layers 0\.\.1"):
        model.add_region("r", cellids=[(5, 2, 3)])
    with pytest.raises(ValueError, match=r"Cell\(s\) \[12\] are off the grid"):
        model.add_region("r", cellids=[NROW * NCOL])
    with pytest.raises(ValueError, match=r"Layer\(s\) \[5\] are off the grid"):
        model.add_region("r", cellids=[1], layer=5)
    assert "r" not in model.regions
    assert model.add_region("r", cellids=[(1, 2, 3)]).layer == 1


def test_a_vertex_grid_region_refuses_a_cell_off_the_grid():
    numbering = cellids.CellNumbering("vertex", nlay=2, ncpl=9)
    with pytest.raises(ValueError, match="off the grid"):
        normalize_region_inputs([(0, 9)], numbering=numbering)
    with pytest.raises(ValueError, match="off the grid"):
        normalize_region_inputs([9], numbering=numbering)
    with pytest.raises(ValueError, match="off the grid"):
        normalize_region_inputs([(2, 0)], numbering=numbering)
    assert normalize_region_inputs([(1, 8)], numbering=numbering)[1:] == ([8], 1)

    grid = SimpleNamespace(grid_type="vertex", nlay=2, ncpl=9)
    regions = RegionRegistry(SimpleNamespace(modelgrid=grid))
    with pytest.raises(ValueError, match="off the grid"):
        regions.add_from_cells("r", [(0, 9)])
    assert "r" not in regions


# --- UZF objects on a DIS grid -------------------------------------------------


def test_a_dis_uzf_maps_each_object_to_the_grids_own_cell(tmp_path):
    model = mf.load_mf6_run(_dis_run(tmp_path / "uzf", uzf=True))
    mapping = model.outputs.uzf.ifno_to_cellid
    expected = [row * NCOL + col for row in range(NROW) for col in range(1, NCOL - 1)]
    assert mapping["layer"].tolist() == [0] * len(expected)
    assert mapping["cellid"].tolist() == expected


# --- a diff needs two models ---------------------------------------------------


def test_a_one_model_group_diff_is_refused_naming_a_second_model(runs):
    model = mf.load_mf6_run(runs.base)
    with pytest.raises(ValueError, match=r"holds only 'base'.*second model"):
        mf.ModelGroup({"base": model}).diff()


# --- review (ledger 419): every other reader of a cellid reads the rule --


def _face_flow(workspace: Path, node_a: int, node_b: int) -> float:
    """MF6's own FLOW-JA-FACE from ``node_a`` to ``node_b``, through its IA/JA."""

    grid = flopy.mf6.utils.MfGrdFile(str(workspace / "m.dis.grb"))
    flows = np.ravel(
        flopy.utils.CellBudgetFile(str(workspace / "m.cbc")).get_data(text="FLOW-JA-FACE")[0]
    )
    for position in range(grid.ia[node_a], grid.ia[node_a + 1]):
        if grid.ja[position] == node_b:
            return float(flows[position])
    raise AssertionError(f"nodes {node_a} and {node_b} are not connected")


def test_a_dis_barrier_names_its_own_face_and_reads_its_own_flow(tmp_path):
    workspace = _dis_run(tmp_path / "hfb", hfb=True)
    model = mf.load_mf6_run(workspace)

    table = model.hfb.inputs.get()
    assert table[["layer", "cell1", "cell2"]].values.tolist() == [[0, 5, 6], [1, 9, 10]]
    assert table["length"].tolist() == [10.0, 10.0]

    # The flow is MF6's across the barrier's own face. Read on row 0's open face it
    # was 5.56 and 14.16 ft3/d through barriers passing 8.5e-5 and 2.5e-4.
    q = model.hfb.results.q.get()
    ncpl = NROW * NCOL
    for (layer, cell1, cell2), flow in zip(
        table[["layer", "cell1", "cell2"]].values.tolist(), q["q_cell1_to_cell2"], strict=True
    ):
        truth = _face_flow(workspace, layer * ncpl + cell1, layer * ncpl + cell2)
        assert flow == pytest.approx(truth, rel=1e-6)
        assert abs(truth) < 1e-3


def test_the_legacy_input_helpers_name_the_grids_own_cells(tmp_path):
    model = mf.load_mf6_run(_dis_run(tmp_path / "legacy", uzf=True, rch=True))

    finf = model.inputs.uzf.finf(per=0, multiplier=1.0)
    expected = [row * NCOL + col for row in range(NROW) for col in range(1, NCOL - 1)]
    assert len(finf) == NROW * NCOL
    assert finf.index[finf.iloc[:, 0] != 0].tolist() == expected
    assert finf.iloc[expected, 0].tolist() == pytest.approx([0.001] * len(expected))

    # cellid[1] made every row's four records one duplicated cell (the row).
    assert sorted(model.inputs.rch.df(per=0).index) == list(range(NROW * NCOL))


def test_a_dis_mover_keeps_its_own_cells_together_when_split(tmp_path):
    from myflopy.modflow.mf6.parallel import _mover_column_groups

    sim = flopy.mf6.MFSimulation(sim_ws=str(tmp_path))
    flopy.mf6.ModflowTdis(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
    flopy.mf6.ModflowGwfdis(gwf, nlay=NLAY, nrow=NROW, ncol=NCOL, top=10.0, botm=[0.0, -10.0])
    reaches = [
        (index, cellid, 10.0, 1.0, 0.001, 9.0, 1.0, 1.0, 0.03, 0, 1.0, 0)
        for index, cellid in enumerate([(0, 2, 3), (0, 1, 2)])
    ]
    flopy.mf6.ModflowGwfsfr(
        gwf, pname="sfr", nreaches=2, packagedata=reaches, connectiondata=[(0,), (1,)]
    )
    flopy.mf6.ModflowGwfmvr(
        gwf, maxmvr=1, maxpackages=1, packages=["sfr"],
        perioddata={0: [("sfr", 0, "sfr", 1, "FACTOR", 1.0)]},
    )
    # The flat index the partition mask uses; the cellids' last items (3, 2) are
    # COLUMNS, which kept two cells of row 0 together instead.
    assert [sorted(int(cell) for cell in group) for group in _mover_column_groups(gwf)] == [
        [1 * NCOL + 2, 2 * NCOL + 3]
    ]

"""A period read of an input view reads the block MODFLOW 6 has IN FORCE, and changes nothing.

Ledger 449. Two defects, one cause -- the views asked FloPy for a period's block
with ``get_data(key=per)``:

* FloPy CREATES the storage it is asked about, and writes it as an EMPTY ``PERIOD``
  block, which MF6 reads as "no records from here on". A GHB written once, in period
  0, read at ``per=3`` and rewritten ran with 0 records in periods 3-5 against 84,
  under Normal termination: a view changed the model.
* A period with no block of its own read NO records, where MF6 runs the last block
  before it in every later period.

Every claim about what MF6 runs is checked against MF6 itself: the tables are compared
with the records MF6 writes to the budget file, period by period.
"""

from __future__ import annotations

import shutil
import warnings

import flopy
import numpy as np
import pytest

import myflopy as mf
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, build_canonical_model
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.mf6.package_explorer_utils import (
    feature_records_in_force,
    period_block_in_force,
    written_period_blocks,
)
from myflopy.modflow.utils.inputs import DrnInput
from myflopy.project.model_group import ModelGroup

# --- a 3 x 3 model whose packages switch blocks ---------------------------------

NPER = 5
#: GHB: two cells in period 0, an EMPTY block in period 2, one cell from period 3.
GHB_BLOCKS = {
    0: [[(0, 2), 5.0, 1.0], [(0, 5), 5.0, 1.0]],
    2: [],
    3: [[(0, 6), 5.5, 1.0]],
}
#: The GHB cells MF6 runs in each period, by MF6's own reuse rule.
GHB_IN_FORCE = {0: [2, 5], 1: [2, 5], 2: [], 3: [6], 4: [6]}
#: HFB: a wall between columns 0 and 1, removed by an empty block, one face back.
HFB_BLOCKS = {
    0: [[(0, 0), (0, 1), 1e-9], [(0, 3), (0, 4), 1e-9], [(0, 6), (0, 7), 1e-9]],
    2: [],
    3: [[(0, 3), (0, 4), 1e-9]],
}
HFB_IN_FORCE = {0: 3, 1: 3, 2: 0, 3: 1, 4: 1}


def _grid(n=3, d=10.0):
    verts = np.array([[i * d, j * d] for j in range(n + 1) for i in range(n + 1)], dtype=float)
    iverts = [
        [j * (n + 1) + i, (j + 1) * (n + 1) + i, (j + 1) * (n + 1) + i + 1, j * (n + 1) + i + 1]
        for j in range(n)
        for i in range(n)
    ]
    xcyc = np.array([[(i + 0.5) * d, (j + 0.5) * d] for j in range(n) for i in range(n)])
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    """The 3 x 3 model, run by MF6: GHB/HFB blocks switch, RCH is array form."""

    vor = _grid()
    context = mf.ModelContext(grid=vor, domain=np.ones((1, 9), dtype=int))
    model = mf.gwf(
        "flow",
        context=context,
        packages=[
            mf.disv(nlay=1, top=10.0, botm=[0.0], **vor.get_disv_gridprops()),
            mf.ic(strt=5.0),
            mf.npf(k=1.0),
            mf.chd(stress_period_data={0: [[(0, 0), 6.0], [(0, 8), 4.0]]}),
            mf.ghb(stress_period_data=GHB_BLOCKS),
            mf.hfb.flopy(stress_period_data=HFB_BLOCKS),
            mf.rch.array(recharge={0: np.full(9, 5e-5), 2: np.full(9, 1e-4)}, irch=0, context=context),
            mf.oc(
                head_filerecord="flow.hds", budget_filerecord="flow.cbc",
                saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
            ),
        ],
    )
    simulation = mf.simulation(
        model, name="blocks", packages=[mf.tdis(nper=NPER, perioddata=[(1.0, 1, 1.0)] * NPER)]
    )
    run = mf.Project(tmp_path_factory.mktemp("blocks"), name="blocks").prepare_run("r", simulation)
    success, report = run.execute(silent=True)
    assert success, "\n".join(report[-15:])
    return run


def _budget_cells(model, text):
    """``{period: sorted zero-based cells}`` of the records MF6 wrote for ``text``."""

    reader = model._get_budget_reader()
    found = {}
    for kstp, kper in reader.get_kstpkper():
        records = reader.get_data(text=text, kstpkper=(kstp, kper))
        found[int(kper)] = sorted(int(node) - 1 for record in records for node in record["node"])
    return found


def test_the_ghb_table_at_each_period_is_what_mf6_ran(run):
    """`get(per=N)` holds exactly the records MF6 wrote to the budget for period N."""

    model = run.model("flow")
    ran = _budget_cells(model, "GHB")
    assert ran == GHB_IN_FORCE  # MF6's own rule, measured: reuse, and an empty block is none
    read = {per: sorted(model.ghb.inputs.get(per=per)["cell"].tolist()) for per in range(NPER)}
    assert read == ran
    assert set(model.ghb.inputs.get(per=4)["per"]) == {4}


def test_a_period_with_no_block_reads_the_block_before_it(run):
    """Period 1 has no GHB block: MF6 runs period 0's, and so does the table."""

    model = run.model("flow")
    period_1 = model.ghb.inputs.get(per=1).drop(columns="per")
    period_0 = model.ghb.inputs.get(per=0).drop(columns="per")
    assert len(period_1) == 2
    assert period_1.equals(period_0)
    assert model.ghb.inputs.cond.map(per=1).fig.data[0].z[2] == pytest.approx(1.0)


def test_no_per_lists_the_blocks_as_written(run):
    model = run.model("flow")
    table = model.ghb.inputs.get()
    assert sorted(table["per"].unique().tolist()) == [0, 3]  # the empty block 2 holds no row
    assert len(table) == 3


def test_an_array_package_reads_each_array_from_the_last_block_that_set_it(run):
    model = run.model("flow")
    budget = model._get_budget_reader()
    for per, rate in ((1, 5e-5), (3, 1e-4)):
        table = model.rch.inputs.get(per=per)
        assert len(table) == 9
        assert set(table["per"]) == {per}
        assert table["recharge"].tolist() == pytest.approx([rate] * 9)
        ran = budget.get_data(text="RCH", kstpkper=(0, per))[0]
        # 100 ft2 cells; MF6 applies none on the two CHD cells, 0 and 8
        applied = {int(n) - 1: float(q) for n, q in zip(ran["node"], ran["q"], strict=True)}
        assert [applied[cell] for cell in range(1, 8)] == pytest.approx([rate * 100.0] * 7)


def test_the_barriers_listed_at_each_step_are_those_mf6_ran(run):
    """HFB blocks replace one another too: the listed faces are the sealed ones."""

    model = run.model("flow")
    flows = model.hfb.results.q.get()
    listed = flows.groupby("per").size().reindex(range(NPER), fill_value=0).to_dict()
    assert listed == HFB_IN_FORCE
    assert flows["q_cell1_to_cell2"].abs().max() < 1e-6  # every listed face is sealed
    # ...and a face the table drops IS open to MF6: the period-0 wall, in period 2.
    ia_ja_reader = model.hfb.results._connectivity()
    ia, ja = ia_ja_reader
    record = model._get_budget_reader().get_data(text="FLOW-JA-FACE", kstpkper=(0, 2))[0]
    position = next(k for k in range(ia[0], ia[1]) if ja[k] == 1)
    assert abs(float(np.ravel(record)[position])) > 1e-3
    assert len(model.hfb.get(per=1)) == 3
    assert len(model.hfb.segments(per=2)) == 0
    assert len(model.hfb.segments(per=4)) == 1


def test_a_disagreeing_per_and_kstpkper_are_refused_naming_both(run):
    model = run.model("flow")
    with pytest.raises(ValueError, match=r"kstpkper=\(0, 3\).*per=2"):
        model.hfb.results.q.map(per=2, kstpkper=(0, 3))


def test_the_results_map_of_a_period_with_no_barrier_draws_the_grid(run):
    """An empty selection used to raise ``KeyError: 'layer'`` from a column-less frame."""

    model = run.model("flow")
    assert model.hfb.results.q.get(per=2).columns.tolist()[-1] == "q_cell1_to_cell2"
    model.hfb.results.q.map(per=2).fig  # noqa: B018 - builds without raising


def test_a_period_outside_the_simulation_is_refused(run):
    model = run.model("flow")
    with pytest.raises(IndexError, match="nper=5"):
        model.ghb.inputs.get(per=NPER)
    with pytest.raises(IndexError, match="nper=5"):
        model.ghb.inputs.get(per=-1)
    with pytest.raises(TypeError, match="stress period"):
        model.ghb.inputs.get(per="3")


def test_the_one_reader_answers_without_adding_a_block(run):
    model = run.model("flow")
    records = model.package("ghb").stress_period_data
    before = (sorted(records.get_active_key_dict()), dict(records.empty_keys))
    assert period_block_in_force(records, 1)[0] == 0
    assert period_block_in_force(records, 2) == (2, None)
    assert sorted(written_period_blocks(records)) == [0, 2, 3]
    assert (sorted(records.get_active_key_dict()), dict(records.empty_keys)) == before


def test_storage_flopy_does_not_write_is_no_block(run, tmp_path):
    """Empty storage is a block only when FloPy WRITES it (ledger 449).

    Any ``get_data(key=N)`` -- a user's own, FloPy's -- leaves empty storage at ``N``.
    With no block header there FloPy writes nothing and MF6 runs the block before it;
    under a header FloPy still holds (it outlives ``remove_transient_key``) it writes
    an EMPTY block and MF6 runs none. The view must read each as MF6 runs it.
    """

    workspace = tmp_path / "copy"
    shutil.copytree(run.workspace, workspace)
    model = mf.load_mf6_run(workspace, model_name="flow")
    records = model.package("ghb").stress_period_data
    records.get_data(key=1)  # no header for period 1: FloPy will write nothing there
    records.remove_transient_key(3)
    records.get_data(key=3)  # period 3's header outlives the removal: an EMPTY block
    read = {per: sorted(model.ghb.inputs.get(per=per)["cell"].tolist()) for per in range(NPER)}
    model.sim.write_simulation(silent=True)
    success, report = model.sim.run_simulation(silent=True)
    assert success, report[-10:]
    ran = _budget_cells(model, "GHB")
    assert ran == {0: [2, 5], 1: [2, 5], 2: [], 3: [], 4: []}
    assert read == ran


def test_the_group_reads_each_member_at_the_block_in_force(run):
    """`group.packages.ghb.inputs.get(per=1)` read no rows: the table held blocks only."""

    model = run.model("flow")
    group = ModelGroup({"a": model, "b": run.model("flow")}, reference="a")
    inputs = group.packages.ghb.inputs
    assert sorted(inputs.get(per=1)["cell"].tolist()) == [2, 2, 5, 5]
    assert inputs.get(per=2).empty
    assert set(inputs.get(per=4)["per"]) == {4}
    assert len(group.packages.ghb.inputs.compare(per=1)) == 2
    group.packages.ghb.inputs.cond.map(per=1).fig  # noqa: B018
    # A mistyped period is refused, not read as "no member carries GHB".
    with pytest.raises(TypeError, match="stress period"):
        inputs.get(per="1")


# --- the canonical model: every view at every period, then MF6 ---------------------


def _canonical(path):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=r"mf\.lak: lake")
        return build_canonical_model(path, config=CanonicalModelConfig.testing())


def _block_bookkeeping(model):
    """Every transient dataset's stored and empty period keys, package by package."""

    state = {}
    for package in model.gwf.packagelist:
        for name in ("stress_period_data", "perioddata"):
            dataset = getattr(package, name, None)
            if dataset is not None and hasattr(dataset, "get_active_key_dict"):
                state[(package.package_name, name)] = (
                    sorted(dataset.get_active_key_dict()),
                    sorted(k for k, v in (dataset.empty_keys or {}).items() if v),
                )
    return state


@pytest.mark.slow
def test_reading_every_input_view_at_every_period_changes_no_model(tmp_path):
    """Read every list and UZF input view at every period; MF6 then runs the same model."""

    model = _canonical(tmp_path / "read")
    for period in range(1, model.nper):  # MF6 runs GHB's period-0 block in all six
        model.package("ghb").stress_period_data.remove_transient_key(period)
    before = _block_bookkeeping(model)

    for pkg in ("chd", "ghb", "rch", "evt", "wel", "drn", "riv"):
        view = getattr(model, pkg).inputs
        for period in range(model.nper):
            assert len(view.get(per=period)), (pkg, period)
    for field in ("finf", "pet", "extdp", "extwc"):
        for period in range(model.nper):
            assert len(getattr(model.uzf.inputs, field).get(per=period)), (field, period)
    for period in range(model.nper):
        model.inputs.uzf.finf(per=period)
    model.ghb.inputs.cond.map(per=3).fig  # noqa: B018
    model.uzf.inputs.finf.map(per=3).fig  # noqa: B018

    assert _block_bookkeeping(model) == before
    model.sim.write_simulation(silent=True)
    success, report = model.sim.run_simulation(silent=True)
    assert success, report[-10:]
    ran = _budget_cells(model, "GHB")
    assert [len(ran[period]) for period in range(model.nper)] == [84] * model.nper
    read = model.ghb.inputs.get(per=5)
    ncpl = model.vor.ncpl
    assert sorted((read["layer"] * ncpl + read["cell"]).tolist()) == ran[5]


@pytest.mark.slow
def test_a_uzf_block_sets_only_the_objects_it_lists(tmp_path):
    """UZF's rule is per OBJECT: a block listing ten leaves the rest as they were.

    Measured through MF6: every land-surface object the period-2 block does not
    list still infiltrates in period 3, and the ten it sets to zero do not.
    """

    model = _canonical(tmp_path / "uzf")
    perioddata = model.package("uzf").perioddata
    first = perioddata.get_data(key=0).copy()
    for period in range(1, model.nper):
        perioddata.remove_transient_key(period)
    listed = first[:10].copy()
    listed["finf"] = 0.0
    perioddata.set_data({0: first, 2: listed})
    model.package("uzf").budget_filerecord.set_data([("uzf_test.cbc",)])

    table = model.uzf.inputs.finf.get(per=3).set_index("ifno")["finf"]
    assert len(table) == len(first)
    assert (table.loc[:9] == 0.0).all()
    assert table.loc[10:].tolist() == pytest.approx(first["finf"][10:].tolist())
    frame = feature_records_in_force(perioddata, 3, id_field="ifno")
    assert len(frame) == len(first)
    legacy = model.inputs.uzf.finf(per=np.int64(3), multiplier=1)["finf"]
    by_cell = model.uzf.inputs.finf.get(per=3).set_index("cell")["finf"]
    assert (legacy.loc[by_cell.index] == by_cell).all()  # the legacy helper reads the same rule

    model.sim.write_simulation(silent=True)
    success, report = model.sim.run_simulation(silent=True)
    assert success, report[-10:]
    budget = flopy.utils.CellBudgetFile(str(model.sim.sim_path / "uzf_test.cbc"))
    def infiltrating(period):
        records = budget.get_data(text="INFILTRATION", kstpkper=(0, period))[0]
        return {int(r["node"]) - 1: float(r["q"]) for r in records}

    before, after = infiltrating(0), infiltrating(3)
    assert all(after[i] == 0.0 for i in range(10))
    # the land-surface objects the block left alone infiltrate as they did
    unlisted = [i for i in range(10, len(first)) if before[i] > 0.0]
    assert len(unlisted) > 100
    assert all(after[i] > 0.0 for i in unlisted)


@pytest.mark.slow
def test_the_legacy_record_helpers_read_the_block_in_force(tmp_path):
    """`model.inputs.rch.df/map` and `DrnInput.map` read ``.data[per]``: ``KeyError``."""

    model = _canonical(tmp_path / "legacy")
    for name in ("rch", "drn"):
        for period in range(1, model.nper):
            model.package(name).stress_period_data.remove_transient_key(period)
    before = _block_bookkeeping(model)
    assert model.inputs.rch.df(per=3).equals(model.inputs.rch.df(per=0))
    model.inputs.rch.map().fig  # noqa: B018 - per=None reads self.per, as df() does
    model.inputs.rch.map(per=4).fig  # noqa: B018
    DrnInput(model).map(per=5).fig  # noqa: B018 - not bound on model.inputs; built directly
    assert _block_bookkeeping(model) == before

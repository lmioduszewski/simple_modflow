"""The HEAD file's clock and a budget file's clock are two different clocks.

One MF6 ``OC`` ``saverecord`` sets two INDEPENDENT flags -- ``idvsave`` from
``SAVE HEAD`` and ``icbcfl`` from ``SAVE BUDGET`` (``gwf.f90``). LAK and SFR gate
their STAGE file on the first and their BUDGET file on the second, so an ordinary
OC such as ``SAVE HEAD LAST`` + ``SAVE BUDGET ALL`` writes a different number of
records to each. Until ledger 301 ``budget_tables.package_output_df`` labelled
package-budget records by zipping them positionally against ``model.kstpkper``
-- the HEAD clock -- which crashed when the counts differed and silently
mislabelled when they happened to match.

Every ``saverecord`` elsewhere in this suite is symmetric (both ``ALL`` or both
``LAST``), which is exactly why the whole class survived. These fixtures are
deliberately asymmetric, and each test asserts the mismatch ACTUALLY HOLDS before
asserting the behaviour -- otherwise a later fixture change turns them into tests
of nothing (ledger 96's lesson).
"""

from __future__ import annotations

import flopy
import pandas as pd
import pytest
from flopy.utils import CellBudgetFile

import myflopy as mf
from myflopy.modflow.mf6.budget_tables import package_output_df
from myflopy.modflow.mf6.output_files import _declared_output_path


def _build_and_run(workspace, name, saverecord, *, nstp, perlen, budget_suffix="cbb"):
    """A 5x5 one-lake model whose OC saves heads and budgets on different clocks."""

    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=str(workspace), exe_name="mf6")
    flopy.mf6.ModflowTdis(
        sim, nper=2, time_units="days", perioddata=[(perlen, nstp, 1.0)] * 2
    )
    flopy.mf6.ModflowIms(sim, complexity="MODERATE", linear_acceleration="BICGSTAB")
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, save_flows=True, newtonoptions="NEWTON"
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=5, ncol=5, delr=100.0, delc=100.0, top=10.0, botm=[-10.0]
    )
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=10.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-5, sy=0.2, transient={0: True})
    flopy.mf6.ModflowGwfic(gwf, strt=8.0)
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=(
            [[(0, row, 0), 9.0] for row in range(5)]
            + [[(0, row, 4), 7.0] for row in range(5)]
        ),
    )
    flopy.mf6.ModflowGwflak(
        gwf,
        boundnames=True,
        save_flows=True,
        budget_filerecord=f"{name}_lak.bud",
        stage_filerecord=f"{name}_lak.stg",
        nlakes=1,
        noutlets=0,
        ntables=0,
        packagedata=[[0, 8.5, 1, "lake1"]],
        connectiondata=[[0, 0, (0, 2, 2), "VERTICAL", 1.0, 0.0, 0.0, 0.0, 0.0]],
        perioddata={0: [[0, "RAINFALL", 0.001]]},
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        # NOT ``<name>.cbc``: the modeller names this file, and a reader that
        # assumes the default extension cannot open it at all (ledger 301).
        budget_filerecord=f"{name}.{budget_suffix}",
        saverecord=saverecord,
    )
    sim.write_simulation(silent=True)
    success, buff = sim.run_simulation(silent=True)
    assert success, "\n".join(buff[-30:])
    return workspace


@pytest.fixture(scope="session")
def budget_more_often_than_heads(tmp_path_factory):
    """``SAVE HEAD LAST`` + ``SAVE BUDGET ALL``: 2 head records, 10 LAK budget records."""

    workspace = tmp_path_factory.mktemp("clock_budget_more")
    _build_and_run(
        workspace, "a", [("HEAD", "LAST"), ("BUDGET", "ALL")], nstp=5, perlen=10.0
    )
    return workspace, "a"


@pytest.fixture(scope="session")
def same_count_different_steps(tmp_path_factory):
    """``STEPS 1`` vs ``STEPS 3``: two records each, naming DIFFERENT time steps."""

    workspace = tmp_path_factory.mktemp("clock_same_count")
    _build_and_run(
        workspace,
        "b",
        {0: [("HEAD", "STEPS", 1), ("BUDGET", "STEPS", 3)]},
        nstp=3,
        perlen=9.0,
    )
    return workspace, "b"


@pytest.fixture(scope="session")
def two_packages_of_one_type(tmp_path_factory):
    """TWO ``RIV6`` entries -- so the budget file holds 2 RIV records per saved time."""

    workspace = tmp_path_factory.mktemp("clock_two_riv")
    name = "c"
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=str(workspace), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=2, time_units="days", perioddata=[(10.0, 1, 1.0)] * 2)
    flopy.mf6.ModflowIms(sim, complexity="MODERATE", linear_acceleration="BICGSTAB")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=5, ncol=5, delr=100.0, delc=100.0, top=10.0, botm=[-10.0]
    )
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=10.0)
    flopy.mf6.ModflowGwfic(gwf, strt=8.0)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, row, 0), 9.0] for row in range(5)]
    )
    # Two RIV packages, as a model with a main stem and a named tributary has.
    flopy.mf6.ModflowGwfriv(
        gwf,
        pname="riv",
        filename=f"{name}.riv",
        stress_period_data=[[(0, 1, col), 8.0, 100.0, 7.0] for col in range(1, 4)],
    )
    flopy.mf6.ModflowGwfriv(
        gwf,
        pname="trib",
        filename=f"{name}.trib",
        stress_period_data=[[(0, 3, col), 8.0, 100.0, 7.0] for col in range(1, 5)],
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        budget_filerecord=f"{name}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    sim.write_simulation(silent=True)
    success, buff = sim.run_simulation(silent=True)
    assert success, "\n".join(buff[-30:])
    return workspace, name


@pytest.mark.slow
def test_two_packages_of_one_type_do_not_overrun_the_budget_clock(
    two_packages_of_one_type,
):
    """``model.bud("riv").df`` on a model with two RIV packages.

    The file's UNIQUE save times are one per record only when the model has one
    package of that type. With two, ``budget_df`` zipped 2N records against N
    times and raised ``IndexError`` on record N -- measured on a real 38-period
    model carrying a main stem and a named tributary (76 records, 38 times).
    """

    workspace, name = two_packages_of_one_type
    model = mf.load_mf6_run(str(workspace))
    reader = model._get_budget_reader()
    riv_records = len(reader.get_data(text="RIV"))
    save_times = len(reader.get_kstpkper())
    assert riv_records == 2 * save_times, "fixture no longer has two RIV packages"

    frame = model.bud("riv").df
    assert len(set(frame.index.get_level_values(1))) == save_times
    # Both packages' cells are present: 3 main-stem + 4 tributary per saved time.
    assert sorted(frame.groupby(level=1).size().unique()) == [7]


def _lak_budget_clock(workspace, name):
    """The (kstp, kper) the LAK package budget file itself recorded, zero-based."""

    reader = CellBudgetFile(str(workspace / f"{name}_lak.bud"), precision="double")
    return sorted({tuple(int(value) for value in pair) for pair in reader.get_kstpkper()})


@pytest.mark.slow
def test_a_lak_budget_is_labelled_by_its_own_clock_not_the_head_file(
    same_count_different_steps,
):
    """Equal record counts, different time steps -- the silent half of ledger 301.

    Reverting ``package_output_df`` to zip against ``model.kstpkper`` makes this
    fail: no length check can see it, and nothing warns.
    """

    workspace, name = same_count_different_steps
    model = mf.load_mf6_run(str(workspace))
    head_clock = sorted(tuple(int(v) for v in pair) for pair in model.kstpkper)
    lak_clock = _lak_budget_clock(workspace, name)

    assert len(head_clock) == len(lak_clock), "fixture no longer has equal counts"
    assert head_clock != lak_clock, "fixture no longer separates the two clocks"

    frame = model.lak.budget.get()
    labelled = sorted({tuple(value) for value in frame.loc[frame.term == "GWF", "kstpkper"]})
    assert labelled == lak_clock
    assert set(frame["per"]) == {kper for _, kper in lak_clock}


@pytest.mark.slow
def test_a_lak_budget_saved_more_often_than_heads_still_reads(
    budget_more_often_than_heads,
):
    """More budget records than head records -- the crashing half of ledger 301.

    This is the user-reported ``AttributeError: 'list' object has no attribute
    'empty'``: the length check raised, ``LakBudget.get`` warned and returned the
    raw record list, and ``build_lak_budget_result_table`` called ``.empty`` on it.
    """

    workspace, name = budget_more_often_than_heads
    model = mf.load_mf6_run(str(workspace))
    head_clock = [tuple(int(v) for v in pair) for pair in model.kstpkper]
    lak_clock = _lak_budget_clock(workspace, name)
    assert len(lak_clock) > len(head_clock), "fixture no longer separates the two clocks"

    frame = model.lak.budget.get()
    assert isinstance(frame, pd.DataFrame)
    gwf = frame.loc[frame.term == "GWF"]
    assert sorted({tuple(value) for value in gwf["kstpkper"]}) == lak_clock

    # The quieter half: the sibling terms used to be dropped one by one by an
    # ``isinstance`` guard that turned the crash into an empty table.
    assert {"GWF", "STORAGE", "RAINFALL"} <= set(frame["term"])
    for term in ("STORAGE", "RAINFALL"):
        assert (frame["term"] == term).sum() == len(lak_clock), term


@pytest.mark.slow
def test_a_package_budget_get_returns_a_dataframe_or_raises_never_a_list(
    budget_more_often_than_heads,
):
    """``return_df=True`` is a contract. The raw list is ``return_df=False``."""

    workspace, _ = budget_more_often_than_heads
    model = mf.load_mf6_run(str(workspace))

    assert isinstance(model.outputs.lak.bud.get("GWF"), pd.DataFrame)
    assert isinstance(model.outputs.lak.bud.get("GWF", return_df=False), list)
    assert isinstance(package_output_df(model, model.fl.lak, bud_type="GWF"), pd.DataFrame)


@pytest.mark.slow
def test_a_failure_inside_the_table_builder_propagates_instead_of_a_raw_list(
    budget_more_often_than_heads, monkeypatch
):
    """The fallback is gone, and nothing can make ``get`` change its return type.

    With the clock fixed nothing raises here in practice, so this forces the
    failure. It certifies the CONTRACT rather than the arithmetic: ``get`` used
    to catch ``ValueError``, warn, and hand back a raw record list -- which
    ``build_lak_budget_result_table`` then called ``.empty`` on. A caller cannot
    guard against a silently different type, and ``return_df=False`` is the
    supported way to ask for the list.
    """

    workspace, _ = budget_more_often_than_heads
    model = mf.load_mf6_run(str(workspace))

    def _explode(*args, **kwargs):
        raise ValueError("boom")

    monkeypatch.setattr(
        "myflopy.modflow.mf6.budget._package_output_df", _explode, raising=True
    )
    with pytest.raises(ValueError, match="boom"):
        model.outputs.lak.bud.get("GWF")
    # ...and still reachable deliberately.
    assert isinstance(model.outputs.lak.bud.get("GWF", return_df=False), list)


@pytest.mark.slow
def test_a_model_budget_file_the_oc_renamed_is_still_found(
    budget_more_often_than_heads,
):
    """``BUDGET FILEOUT <name>.cbb`` -- the reader must not assume ``.cbc``.

    Every cell-budget surface was dead on a real 38-period model for this one
    reason: ``model.bud(pkg)``, ``model.budget.<term>``,
    ``packages.<pkg>.results.q`` and HFB flow all raised ``FileNotFoundError``
    for a ``.cbc`` MF6 never wrote (ledger 301).
    """

    workspace, name = budget_more_often_than_heads
    assert not (workspace / f"{name}.cbc").exists()
    assert (workspace / f"{name}.cbb").exists()

    assert _declared_output_path(workspace, name, "budget").name == f"{name}.cbb"

    model = mf.load_mf6_run(str(workspace))
    assert "CHD" in [str(record).strip() for record in model.budget.types]
    assert not model.chd.results.q.get().empty


def test_a_workspace_whose_oc_names_no_budget_falls_back_to_the_default(tmp_path):
    """No OC file, or no ``BUDGET FILEOUT`` in it -- keep the historical default."""

    assert _declared_output_path(tmp_path, "m", "budget") == tmp_path / "m.cbc"

    (tmp_path / "m.nam").write_text(
        "BEGIN packages\n  OC6  m.oc  oc\nEND packages\n", encoding="utf-8"
    )
    (tmp_path / "m.oc").write_text(
        "BEGIN options\n  HEAD  FILEOUT  m.hds\nEND options\n", encoding="utf-8"
    )
    assert _declared_output_path(tmp_path, "m", "budget") == tmp_path / "m.cbc"
    assert _declared_output_path(tmp_path, "m", "head").name == "m.hds"


@pytest.mark.slow
def test_the_head_clock_and_the_budget_clock_do_not_share_a_cache_slot(
    budget_more_often_than_heads,
):
    """Both orders, on one object. They used to share ``_kstpkper``.

    Whichever accessor was touched first answered for both for the life of the
    model, so the same call gave different answers depending on what a notebook
    had happened to read earlier.
    """

    workspace, name = budget_more_often_than_heads
    expected_head = len([p for p in CellBudgetFile(
        str(workspace / f"{name}.cbb"), precision="double").get_kstpkper()])

    head_first = mf.load_mf6_run(str(workspace))
    head_n = len(head_first.kstpkper)
    budget_after = len(head_first._get_budget_kstpkper())

    budget_first = mf.load_mf6_run(str(workspace))
    budget_n = len(budget_first._get_budget_kstpkper())
    head_after = len(budget_first.kstpkper)

    assert head_n != budget_n, "fixture no longer separates the two clocks"
    assert (head_n, budget_after) == (head_after, budget_n)
    assert budget_n == expected_head

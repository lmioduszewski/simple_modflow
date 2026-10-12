"""The group and diff tiers compare what MODFLOW 6 wrote (ledger 484 and its siblings).

Real MF6 runs, small enough to build per module:

``toy``      a 1 x 5 DIS row, one convertible layer, CHD at 5.0 on the west cell and
             a well on the east cell: at -0.1 every cell is wet; at -50 MF6 writes
             four cells DRY (-1e30) -- read with FloPy, never through the code under
             test, before it is used.
``named``    a 3 x 3 DIS model whose DRN is NAMED ``seeps`` (conductance 10 against
             30) and whose GHB is named ``lake_edge``: a package found by its NAME
             prefix was never compared, and the GHB counted as a LAK.
``movers``   a 3 x 3 DIS model whose two drains (``drn_a``, ``drn_b``) and a well hand
             their water to the mover, all into a one-reach SFR; the second run has
             ``drn_a``'s conductance 30 for 10. MF6 writes ``DRN-TO-MVR``/``WEL-TO-MVR``
             in the groundwater budget file and SFR's ``FROM-MVR`` in its own; the
             diff read ``FROM-MVR``/``TO-MVR`` from the former and was always empty
             (ledger 486).
``steps``    the canonical model on its testing grid (six 30-day periods, a 72-
             connection lake, an SFR-to-lake mover) at one time step per period and
             at two: a period is its END on both, and the LAK join keyed on the
             exact ``kstpkper`` and on MF6's wetted ``FLOW-AREA`` (ledger 485).
"""

from __future__ import annotations

import contextlib
import dataclasses
import io
import warnings

import flopy
import numpy as np
import pytest

import myflopy as mf
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig
from tests.conftest import _build_canonical

# -- the runs ------------------------------------------------------------------------


def _write_run(sim, workspace):
    sim.write_simulation()
    ok, _ = sim.run_simulation(silent=True)
    assert ok, f"{workspace} did not run"
    return mf.load_mf6_run(workspace)


def _toy(workspace, pumping):
    """The 1 x 5 row; ``pumping`` at the east cell decides which cells go dry."""

    sim = flopy.mf6.MFSimulation(sim_name="t", sim_ws=str(workspace), exe_name="mf6",
                                 verbosity_level=0)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, complexity="SIMPLE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="t", save_flows=True)
    flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=1, ncol=5, delr=10.0, delc=10.0,
                            top=10.0, botm=0.0)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0, icelltype=1)
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 0), 5.0]])
    flopy.mf6.ModflowGwfwel(gwf, stress_period_data=[[(0, 0, 4), pumping]])
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord="t.cbc", head_filerecord="t.hds",
                           saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])
    return _write_run(sim, workspace)


def _two_drains(workspace, second_cond):
    """A 3 x 3 model of two DRN packages on ONE cell (``d1`` cond 10, ``d2`` ``second_cond``)."""

    sim = flopy.mf6.MFSimulation(sim_name="t", sim_ws=str(workspace), exe_name="mf6",
                                 verbosity_level=0)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, complexity="SIMPLE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="t", save_flows=True)
    flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=3, ncol=3, delr=10.0, delc=10.0,
                            top=10.0, botm=0.0)
    flopy.mf6.ModflowGwfnpf(gwf, k=10.0)
    flopy.mf6.ModflowGwfic(gwf, strt=9.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, r, 0), 9.0] for r in range(3)])
    flopy.mf6.ModflowGwfdrn(gwf, pname="d1", filename="t1.drn",
                            stress_period_data=[[(0, 1, 2), 5.0, 10.0]])
    flopy.mf6.ModflowGwfdrn(gwf, pname="d2", filename="t2.drn",
                            stress_period_data=[[(0, 1, 2), 6.0, second_cond]])
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord="t.cbc", head_filerecord="t.hds",
                           saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])
    return _write_run(sim, workspace)


def _named(workspace, cond):
    """The 3 x 3 model of a DRN named ``seeps`` and a GHB named ``lake_edge``."""

    sim = flopy.mf6.MFSimulation(sim_name="n", sim_ws=str(workspace), exe_name="mf6",
                                 verbosity_level=0)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, complexity="SIMPLE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="n", save_flows=True)
    flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=3, ncol=3, delr=10.0, delc=10.0,
                            top=10.0, botm=0.0)
    flopy.mf6.ModflowGwfnpf(gwf, k=10.0)
    flopy.mf6.ModflowGwfic(gwf, strt=9.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, r, 0), 9.0] for r in range(3)])
    flopy.mf6.ModflowGwfdrn(gwf, pname="seeps", filename="n.drn",
                            stress_period_data=[[(0, 1, 2), 5.0, cond]])
    flopy.mf6.ModflowGwfghb(gwf, pname="lake_edge", filename="n.ghb",
                            stress_period_data=[[(0, 0, 2), 8.0, 1.0]])
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord="n.cbc", head_filerecord="n.hds",
                           saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])
    return _write_run(sim, workspace)


@pytest.fixture(scope="module")
def toys(tmp_path_factory):
    """``wet`` (every cell wet), ``dry`` (four cells dry) and ``dry_too`` (the same)."""

    root = tmp_path_factory.mktemp("diff_sentinels")
    runs = {"wet": _toy(root / "wet", -0.1), "dry": _toy(root / "dry", -50.0),
            "dry_too": _toy(root / "dry_too", -50.0)}
    for name, wet_cells in (("wet", 5), ("dry", 1), ("dry_too", 1)):
        heads = flopy.utils.HeadFile(str(root / name / "t.hds")).get_data().ravel()
        assert int((np.abs(heads) < 1e29).sum()) == wet_cells, (name, heads)
        assert np.all(heads[np.abs(heads) >= 1e29] == -1e30)
    return runs


@pytest.fixture(scope="module")
def named(tmp_path_factory):
    root = tmp_path_factory.mktemp("diff_named")
    return {"a": _named(root / "a", 10.0), "b": _named(root / "b", 30.0)}


def _movers(workspace, drn_a_cond, *, sfr_budget=True):
    """The 3 x 3 model of two drains and a well moving water into one SFR reach."""

    sim = flopy.mf6.MFSimulation(sim_name="mv", sim_ws=str(workspace), exe_name="mf6",
                                 verbosity_level=0)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, complexity="SIMPLE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="mv", save_flows=True)
    flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=3, ncol=3, delr=10.0, delc=10.0,
                            top=10.0, botm=0.0)
    flopy.mf6.ModflowGwfnpf(gwf, k=10.0)
    flopy.mf6.ModflowGwfic(gwf, strt=9.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, r, 0), 9.0] for r in range(3)])
    flopy.mf6.ModflowGwfdrn(gwf, mover=True, pname="drn_a",
                            stress_period_data=[[(0, 1, 2), 5.0, drn_a_cond]])
    flopy.mf6.ModflowGwfdrn(gwf, mover=True, pname="drn_b", filename="mv_b.drn",
                            stress_period_data=[[(0, 0, 2), 5.0, 10.0]])
    flopy.mf6.ModflowGwfwel(gwf, mover=True, pname="wel",
                            stress_period_data=[[(0, 2, 2), -1.0]])
    flopy.mf6.ModflowGwfsfr(
        gwf, mover=True, pname="sfr", nreaches=1,
        budget_filerecord="mv.sfr.bud" if sfr_budget else None,
        packagedata=[[0, (0, 1, 1), 10.0, 1.0, 0.001, 8.0, 1.0, 0.0, 0.03, 0, 1.0, 0]],
        connectiondata=[[0]], perioddata=[[0, "INFLOW", 0.0]],
    )
    flopy.mf6.ModflowGwfmvr(
        gwf, maxmvr=3, maxpackages=4, packages=[["drn_a"], ["drn_b"], ["wel"], ["sfr"]],
        perioddata={0: [["drn_a", 0, "sfr", 0, "FACTOR", 1.0],
                        ["drn_b", 0, "sfr", 0, "FACTOR", 1.0],
                        ["wel", 0, "sfr", 0, "FACTOR", 1.0]]},
    )
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord="mv.cbc", head_filerecord="mv.hds",
                           saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])
    return _write_run(sim, workspace)


def _moved(path, text):
    """``{package name: q}`` of the records named ``text`` in a budget file, by FloPy alone."""

    budget = flopy.utils.CellBudgetFile(str(path), precision="auto")
    moved = {}
    for index, header in enumerate(budget.recordarray):
        if header["text"].decode().strip() == text:
            record = budget.get_record(index)
            moved[header["paknam2"].decode().strip()] = float(record["q"].sum())
    return moved


@pytest.fixture(scope="module")
def movers(tmp_path_factory):
    root = tmp_path_factory.mktemp("diff_movers")
    runs = {"a": _movers(root / "a", 10.0), "b": _movers(root / "b", 30.0),
            "unwritten": _movers(root / "unwritten", 30.0, sfr_budget=False)}
    for name in ("a", "b", "unwritten"):
        held = {raw.decode().strip() for raw in flopy.utils.CellBudgetFile(
            str(root / name / "mv.cbc"), precision="auto").recordarray["text"]}
        assert {"DRN-TO-MVR", "WEL-TO-MVR"} <= held
        assert not {"FROM-MVR", "TO-MVR"} & held
    runs["moved"] = {name: _moved(root / name / "mv.cbc", "DRN-TO-MVR") for name in ("a", "b")}
    runs["received"] = {name: _moved(root / name / "mv.sfr.bud", "FROM-MVR")
                        for name in ("a", "b")}
    return runs


@pytest.fixture(scope="module")
def steps(tmp_path_factory):
    """The canonical model run at one and at two time steps per period, reloaded."""

    root = tmp_path_factory.mktemp("diff_steps")
    runs = {}
    for name, nstp in (("one", 1), ("two", 2)):
        config = dataclasses.replace(CanonicalModelConfig.testing(), steps_per_period=nstp)
        with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = _build_canonical(root / name, config)
            ok, _ = model.run_simulation()
        assert ok, f"{name} did not run"
        runs[name] = mf.load_mf6_run(model.workspace)
    saved = {name: {tuple(k) for k in run.kstpkper} for name, run in runs.items()}
    assert saved["one"] == {(0, p) for p in range(6)}
    assert saved["two"] == {(k, p) for p in range(6) for k in (0, 1)}
    return runs


# -- (a) MF6's no-value numbers are a gap, not a difference --------------------------


def test_a_dry_cell_in_one_member_is_a_gap_not_1e30(toys):
    diff = mf.ModelGroup({"wet": toys["wet"], "dry": toys["dry"]}, reference="wet").diff()
    table = diff.hds.get()
    assert len(table) == 5
    assert np.isfinite(table["diff"]).sum() == 1
    assert table["elev"].isna().sum() == 4
    assert float(np.nanmax(np.abs(table["diff"]))) < 1.0

    summary = diff.hds.summary().iloc[0]
    assert summary["n"] == 1
    assert summary["n_one_sided"] == 4
    assert summary["max_abs_diff"] < 1.0
    assert not bool(summary["within_tolerance"])
    assert "4 cell(s) with no value in one model only" in diff.report(results=True)


def test_a_cell_dry_in_both_members_agrees(toys):
    diff = mf.ModelGroup({"dry": toys["dry"], "dry_too": toys["dry_too"]},
                         reference="dry").diff()
    summary = diff.hds.summary().iloc[0]
    assert (summary["n"], summary["n_one_sided"]) == (1, 0)
    assert bool(summary["within_tolerance"])


def test_a_group_compare_and_its_map_leave_a_dry_cell_blank(toys):
    group = mf.ModelGroup({"wet": toys["wet"], "dry": toys["dry"]}, reference="wet")
    assert group.hds.get(per=0)["elev"].isna().sum() == 4
    compared = group.hds.compare(per=0)
    assert compared["diff"].isna().sum() == 4
    picture = group.diff().hds.map("dry", per=0, layer=0)
    drawn = np.asarray(picture.fig.data[0].z, dtype=float)
    assert np.isnan(drawn).sum() == 4
    assert np.nanmax(np.abs(drawn)) < 1.0


def test_a_long_and_wide_read_a_dry_cell_as_nan(toys):
    hds = toys["dry"].hds
    assert hds.all_values["elev"].min() == -1e30     # the file's own numbers
    assert hds.long().isna().sum() == 4
    assert np.nanmin(hds.long().to_numpy()) > -1.0
    wide = hds.wide()
    assert wide.filter(like="kstpkper_").isna().to_numpy().sum() == 4


# -- (d) a package is found by its TYPE, whatever it is named ------------------------


def test_a_package_named_otherwise_is_compared(named):
    diff = mf.ModelGroup(named, reference="a").diff()
    assert diff.package_names == ["chd", "drn", "ghb"]
    assert diff.connection_package_names == []
    assert "## `b` -- differs from reference" in diff.report()
    drn = diff.summary().set_index("package").loc["drn"]
    assert int(drn["value_cells_changed"]) == 1
    assert not bool(drn["identical"])


def test_diff_package_answers_package_types_only(named):
    diff = mf.ModelGroup(named, reference="a").diff()
    assert diff.package("DRN").package_name == "drn"
    for name in ("_diff", "__init__", "nope"):
        with pytest.raises(AttributeError, match="does not diff package"):
            diff.package(name)


# -- (b) a period is its end, and a connection is paired by what it IS ---------------


def test_lak_exchange_pairs_every_connection_at_each_period_end(steps):
    group = mf.ModelGroup({"one": steps["one"], "two": steps["two"]}, reference="one")
    compared = group.packages.lak.results.q.compare()
    assert compared.groupby("per").size().to_dict() == {p: 72 for p in range(6)}
    assert set(compared["kstpkper"]) == {(1, p) for p in range(6)}
    diff_rows = group.diff().packages.lak.results.q.get(per=1)
    assert set(diff_rows["kstpkper"]) == {(1, 1)}


def test_a_connection_whose_wetted_area_moved_still_pairs(steps, monkeypatch):
    """MF6's ``FLOW-AREA`` moves with the stage; it is a value, never a join key."""

    import myflopy.project.group.lak as group_lak

    real = group_lak.build_lak_budget_result_table

    def wetter(model, **kwargs):
        table = real(model, **kwargs)
        if model is steps["two"]:
            for column in ("flow_area", "FLOW-AREA"):
                table[column] = table[column] * 1.5
        return table

    monkeypatch.setattr(group_lak, "build_lak_budget_result_table", wetter)
    group = mf.ModelGroup({"one": steps["one"], "two": steps["two"]}, reference="one")
    compared = group.packages.lak.results.q.compare(per=3)
    assert len(compared) == 72
    wet = compared["reference_flow_area"] > 0
    assert np.allclose(compared.loc[wet, "flow_area_diff"],
                       0.5 * compared.loc[wet, "reference_flow_area"])


def test_group_budget_compares_each_period_end(steps):
    group = mf.ModelGroup({"one": steps["one"], "two": steps["two"]}, reference="one")
    compared = group.bud("drn").compare()
    assert set(compared["kstpkper"]) == {(1, p) for p in range(6)}
    assert compared.groupby("per").size().nunique() == 1
    assert "kstpkper" not in compared.columns[3:]


def test_listing_budget_diff_has_no_bookkeeping_terms(steps):
    diff = mf.ModelGroup({"one": steps["one"], "two": steps["two"]}, reference="one").diff()
    terms = set(diff.bud.summary()["term"])
    assert not terms & {"time_step", "stress_period", "tslen"}
    assert {"TOTAL_IN", "PERCENT_DISCREPANCY"} <= terms
    rows = diff.bud.get()
    assert set(rows["per"]) == set(range(6))
    assert set(rows["kstpkper"]) == {(1, p) for p in range(6)}


# -- (c) the mover diff reads what MF6 writes ------------------------------------------


def test_mover_diff_reads_each_providers_term_from_the_groundwater_budget(movers):
    moved = movers["moved"]
    assert moved["a"]["DRN_A"] < 0 and moved["b"]["DRN_A"] < moved["a"]["DRN_A"]
    mvr = mf.ModelGroup({"a": movers["a"], "b": movers["b"]}, reference="a").diff()
    drn = mvr.packages.mvr.results.get(package="drn").set_index("pname")
    assert set(drn.index) == {"DRN_A", "DRN_B"}
    for pname in ("DRN_A", "DRN_B"):
        assert drn.loc[pname, "reference_q"] == pytest.approx(moved["a"][pname])
        assert drn.loc[pname, "q"] == pytest.approx(moved["b"][pname])
    wel = mvr.packages.mvr.results.get(package="wel")
    assert len(wel) == 1 and wel["q"].iloc[0] < 0


def test_mover_diff_reads_a_receivers_term_from_its_own_budget_file(movers):
    received = movers["received"]
    mvr = mf.ModelGroup({"a": movers["a"], "b": movers["b"]}, reference="a").diff()
    sfr = mvr.packages.mvr.results.get(package="sfr", direction="from_mvr")
    assert len(sfr) == 1
    assert sfr["reference_q"].iloc[0] == pytest.approx(received["a"]["SFR"])
    assert sfr["q"].iloc[0] == pytest.approx(received["b"]["SFR"])
    assert sfr["q_diff"].iloc[0] > 50.0     # drn_a's conductance tripled

    summary = mvr.packages.mvr.results.summary().set_index(["package", "direction"])
    assert set(summary.index) == {("drn", "to_mvr"), ("wel", "to_mvr"),
                                  ("sfr", "from_mvr"), ("sfr", "to_mvr")}
    assert not bool(summary.loc[("drn", "to_mvr"), "within_tolerance"])
    assert summary.loc[("sfr", "from_mvr"), "argmax_feature"] == 0


def test_mover_diff_refuses_a_term_mf6_does_not_write(movers):
    mvr = mf.ModelGroup({"a": movers["a"], "b": movers["b"]},
                        reference="a").diff().packages.mvr.results
    with pytest.raises(ValueError, match="DRN-TO-MVR.*direction='to_mvr'"):
        mvr.get(package="drn", direction="from_mvr")
    with pytest.raises(ValueError, match="UZF's own budget file"):
        mvr.get(package="uzf")
    with pytest.raises(ValueError, match="'to_mvr'"):
        mvr.get(package="sfr", direction="TO-MVR")
    with pytest.raises(ValueError, match="not a cell"):
        mvr.get(package="sfr", layer=0)


def test_the_canonical_stream_feeds_its_lake_through_the_mover(steps):
    mvr = mf.ModelGroup({"one": steps["one"], "two": steps["two"]},
                        reference="one").diff().packages.mvr.results
    lake = mvr.get(package="lak", direction="from_mvr")
    assert len(lake) == 6 and (lake["reference_q"] > 0).all()
    stream = mvr.get(package="sfr")
    assert (stream["reference_q"] < 0).any()


def test_mover_summary_says_it_leaves_out_a_uzf_mover(movers, monkeypatch):
    import myflopy.project.model_results_diff as results_diff

    real = results_diff._declares_mover
    monkeypatch.setattr(results_diff, "_declares_mover",
                        lambda model, kind: kind == "uzf" or real(model, kind))
    mvr = mf.ModelGroup({"a": movers["a"], "b": movers["b"]},
                        reference="a").diff().packages.mvr.results
    with pytest.warns(UserWarning, match="leaves out UZF's mover terms"):
        mvr.summary()


# -- (e) a stage diff answers the verbs its namespace lists ---------------------------


@pytest.mark.parametrize("package, entity", [("lak", "lake"), ("sfr", "reach")])
def test_a_stage_diff_is_mapped_on_the_reference_features_cells(steps, package, entity):
    from myflopy.modflow.mf6.package_budget import (
        build_lak_stage_result_table,
        build_sfr_stage_result_table,
    )

    diff = mf.ModelGroup({"one": steps["one"], "two": steps["two"]}, reference="one").diff()
    results = getattr(diff.packages, package).results
    assert "stage" in results.field_names()
    stage = results.stage
    compared = stage.get(per=3).set_index(entity)["stage_diff"]
    table = {"lak": build_lak_stage_result_table,
             "sfr": build_sfr_stage_result_table}[package](steps["one"])
    cells = table.loc[table["layer"] == 0, [entity, "cell"]].drop_duplicates()
    expected = cells.assign(diff=cells[entity].map(compared)).groupby("cell")["diff"].mean()

    drawn = np.asarray(results.map(field="stage", per=3, layer=0).fig.data[0].z, dtype=float)
    assert np.allclose(drawn[expected.index.to_numpy()], expected.to_numpy(), equal_nan=True)
    others = np.setdiff1d(np.arange(drawn.size), expected.index.to_numpy())
    assert np.all(drawn[others] == 0.0)
    assert np.abs(expected).max() > 0

    assert stage.plot() is not None
    assert stage.mosaic() is not None
    assert stage.map(per=3, backend="mpl") is not None


def test_a_receiver_with_no_budget_file_is_refused_not_read_as_nothing(movers):
    group = mf.ModelGroup({"a": movers["a"], "c": movers["unwritten"]}, reference="a")
    mvr = group.diff().packages.mvr.results
    with pytest.raises(ValueError, match="'c' \\(sfr\\).*budget_filerecord="):
        mvr.get(package="sfr", direction="from_mvr")
    with pytest.warns(UserWarning, match="leaves out SFR's mover terms"):
        summary = mvr.summary()
    assert set(summary["package"]) == {"drn", "wel"}


def test_group_budget_pairs_two_packages_on_one_cell_in_file_order(tmp_path):
    a, b = _two_drains(tmp_path / "a", 10.0), _two_drains(tmp_path / "b", 40.0)
    by_package = {
        name: _moved(tmp_path / name / "t.cbc", "DRN") for name in ("a", "b")
    }
    compared = mf.ModelGroup({"a": a, "b": b}, reference="a").bud("drn").compare()
    assert len(compared) == 2
    for row, pname in zip(compared.itertuples(), ("D1", "D2"), strict=True):
        assert row.reference_q == pytest.approx(by_package["a"][pname])
        assert row.q == pytest.approx(by_package["b"][pname])

"""The group and diff tiers compare what MODFLOW 6 wrote (ledger 484 and its siblings).

Real MF6 runs, small enough to build per module:

``toy``      a 1 x 5 DIS row, one convertible layer, CHD at 5.0 on the west cell and
             a well on the east cell: at -0.1 every cell is wet; at -50 MF6 writes
             four cells DRY (-1e30) -- read with FloPy, never through the code under
             test, before it is used.
``named``    a 3 x 3 DIS model whose DRN is NAMED ``seeps`` (conductance 10 against
             30) and whose GHB is named ``lake_edge``: a package found by its NAME
             prefix was never compared, and the GHB counted as a LAK.
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

"""A USG export's boundary layers carry EVERY period's cells, and a switch-off shows.

``usg/export._bc_frames`` took its cells from the FIRST defined period alone, so a
cell that only a later period gives a record never reached ``boundaries.gpkg``, and
a period switched off with ``ITMP = 0`` left no trace there; ``boundary_frame()``
dropped that period's (empty) rows too, so its rows could not tell a switch-off from
a reuse (ledger 430).

The model here is ``test_usg_export_adjacency``'s one-layer square grid with a DRN of
four periods:

* period 0 defines drains on cells 0 and 1;
* period 1 defines cells 0, 1 and 2 -- cell 2 is NEW, with another elevation;
* period 2 is ``ITMP = 0``: the drains are switched off;
* period 3 reuses period 2 (``ITMP < 0``), so they stay off.
"""

from __future__ import annotations

import numpy as np
import pytest

gpd = pytest.importorskip("geopandas")

from myflopy.modflow.usg.export import (  # noqa: E402
    ExportManifest,
    _bc_frames,
    _export_boundaries,
)
from myflopy.modflow.usg.packages import PeriodRecords  # noqa: E402
from tests.test_usg_export_adjacency import CRS, square_grid, usg_model  # noqa: E402

NPER = 4


def _drain_model():
    grid = square_grid()
    model = usg_model(grid, [0, 1])
    model.disu.nper = NPER
    model.disu.perlen = np.ones(NPER)
    model.disu.nstp = np.ones(NPER, int)
    model.disu.tsmult = np.ones(NPER)
    model.disu.steady = np.zeros(NPER, bool)
    model.boundaries = {"DRN": PeriodRecords(
        ftype="DRN",
        periods={
            0: np.array([[1, 7.5, 100.0], [2, 7.5, 100.0]]),
            1: np.array([[1, 7.5, 100.0], [2, 7.5, 100.0], [3, 6.0, 50.0]]),
            2: np.empty((0, 3)),
        },
        reused=(3,),
    )}
    return model


def test_a_cell_only_a_later_period_names_is_exported():
    points, lines = _bc_frames(_drain_model(), "DRN", CRS)
    assert sorted(points["old_cell"]) == [0, 1, 2]
    assert sorted(lines["old_cell"]) == [0, 1, 2]
    new = points.set_index("old_cell").loc[2]
    assert new["elevation"] == 6.0 and new["conductance"] == 50.0
    assert new["first_period"] == 1


def test_each_feature_says_which_periods_it_acts_in():
    """Periods 2 and 3 are off -- the switch-off and its reuse -- for every feature."""

    points, _ = _bc_frames(_drain_model(), "DRN", CRS)
    periods = points.set_index("old_cell")["periods"].to_dict()
    assert periods == {0: "0-1", 1: "0-1", 2: "1"}


def test_the_cycle_columns_are_model_periods_and_a_switch_off_is_nan():
    points, _ = _bc_frames(_drain_model(), "DRN", CRS)
    table = points.set_index("old_cell")
    assert [c for c in table.columns if c.startswith("elevation_m")] == [
        "elevation_m01", "elevation_m02", "elevation_m03", "elevation_m04"
    ]
    assert np.isnan(table.loc[2, "elevation_m01"]) and table.loc[2, "elevation_m02"] == 6.0
    assert table[["elevation_m03", "elevation_m04"]].isna().all().all()


def test_the_manifest_names_the_switched_off_periods(tmp_path):
    manifest = ExportManifest(directory=tmp_path, crs=CRS, nper=NPER)
    _export_boundaries(_drain_model(), tmp_path, CRS, manifest)
    assert manifest.layers["boundaries.gpkg:drn"]["switched_off_periods"] == [2, 3]
    assert manifest.layers["boundaries.gpkg:drn_lines"]["switched_off_periods"] == [2, 3]
    assert any("2 of 4 period(s) (zero-based: 2-3)" in note for note in manifest.notes)
    written = gpd.read_file(tmp_path / "boundaries.gpkg", layer="drn")
    assert sorted(written["old_cell"]) == [0, 1, 2]


def test_boundary_frame_tells_a_switch_off_from_a_reuse():
    frame = _drain_model().boundary_frame("DRN")
    assert sorted(frame["period"].unique()) == [0, 1]
    assert frame.attrs["defined_periods"] == (0, 1, 2)
    assert frame.attrs["empty_periods"] == (2,)
    assert frame.attrs["reused_periods"] == (3,)


def test_a_package_that_never_moves_writes_no_cycle_columns(tmp_path):
    """One defined period reused through the rest: nothing to cycle, nothing off."""

    model = _drain_model()
    model.boundaries["DRN"] = PeriodRecords(
        ftype="DRN", periods={0: np.array([[1, 7.5, 100.0], [2, 7.5, 100.0]])},
        reused=(1, 2, 3),
    )
    points, _ = _bc_frames(model, "DRN", CRS)
    assert not [c for c in points.columns if "_m0" in c]
    assert set(points["periods"]) == {"0-3"}
    manifest = ExportManifest(directory=tmp_path, crs=CRS, nper=NPER)
    _export_boundaries(model, tmp_path, CRS, manifest)
    assert manifest.layers["boundaries.gpkg:drn"]["switched_off_periods"] == []
    assert manifest.layers["boundaries.gpkg:drn"]["features_in_some_periods"] == 0
    assert not manifest.notes


def test_the_manifest_counts_features_in_force_in_only_some_periods(tmp_path):
    """Cell 2 acts in period 1 alone of the package's 0-1: a layer read with nper= would
    apply it in period 0 too, and only its own row said otherwise (Review of U6, ledger
    430). The two drains on cells 0 and 1 act in every period the package does."""

    manifest = ExportManifest(directory=tmp_path, crs=CRS, nper=NPER)
    _export_boundaries(_drain_model(), tmp_path, CRS, manifest)
    assert manifest.layers["boundaries.gpkg:drn"]["features_in_some_periods"] == 1
    assert manifest.layers["boundaries.gpkg:drn_lines"]["features_in_some_periods"] == 1
    assert any("1 of 3 feature(s) are in force in only some of the periods the package "
               "is (0-1)" in note for note in manifest.notes)


def _two_records_on_one_cell(ftype, first, second):
    """Records ``first`` and ``second`` on cell 0 (one-based node 1), one more on cell 1."""

    grid = square_grid()
    model = usg_model(grid, [0, 1])
    model.boundaries = {ftype: PeriodRecords(
        ftype=ftype, periods={0: np.array([[1, *first], [1, *second], [2, *first]])},
    )}
    return model


def _two_drains_on_one_cell(second_elevation, second_conductance=50.0):
    return _two_records_on_one_cell("DRN", (7.5, 100.0), (second_elevation, second_conductance))


def test_two_records_on_one_cell_are_one_feature_combined_by_kind():
    """Two drains on one cell in one period -- MF6 applies both -- stopped the export
    on ``float()`` of a Series (ledger 467). Conductance adds; the elevation is the
    CONDUCTANCE-weighted mean (ledger 493): 7.5 x 100 + 7.0 x 50 over 150 is
    7.333, where the plain mean 7.25 drained 112.5 for MF6's 100.0. A warning names the
    cell, because the weighted drain is MF6's only while the head is not between."""

    with pytest.warns(UserWarning, match=r"DRN: 1 cell\(s\) .*elev values differ .*\[\(1, 0\)\].*conductance-weighted"):
        points, lines = _bc_frames(_two_drains_on_one_cell(7.0), "DRN", CRS)
    table = points.set_index("old_cell")
    assert sorted(table.index) == [0, 1]
    assert table.loc[0, "conductance"] == 150.0
    assert table.loc[0, "elevation"] == pytest.approx((7.5 * 100.0 + 7.0 * 50.0) / 150.0)
    assert table.loc[1, "conductance"] == 100.0 and table.loc[1, "elevation"] == 7.5
    assert lines.set_index("old_cell").loc[0, "conductance"] == 150.0
    assert points.attrs == {"records_combined": 1, "cells_inexact": 1, "cells_unweighted": 0}


def test_a_ghb_s_weighted_head_is_exact_so_nothing_warns(recwarn):
    """GHB is linear in the head, so the conductance-weighted head is MF6's boundary at
    every head: the cell is combined and not counted as inexact."""

    points, _ = _bc_frames(_two_records_on_one_cell("GHB", (9.0, 100.0), (6.0, 50.0)), "GHB", CRS)
    assert not [w for w in recwarn if "several records" in str(w.message)]
    row = points.set_index("old_cell").loc[0]
    assert row["conductance"] == 150.0 and row["head"] == pytest.approx(8.0)
    assert points.attrs == {"records_combined": 1, "cells_inexact": 0, "cells_unweighted": 0}


def test_a_river_s_stage_and_rbot_are_both_weighted_by_conductance():
    # Only the differing rbot warns: copies differing in stage alone are exact at every
    # head (ledger 494).
    with pytest.warns(UserWarning, match=r"RIV: 1 cell\(s\) .*whose rbot values differ"):
        points, _ = _bc_frames(
            _two_records_on_one_cell("RIV", (9.0, 100.0, 7.5), (8.0, 50.0, 7.0)), "RIV", CRS
        )
    row = points.set_index("old_cell").loc[0]
    assert row["conductance"] == 150.0
    assert row["stage"] == pytest.approx((9.0 * 100.0 + 8.0 * 50.0) / 150.0)
    assert row["rbot"] == pytest.approx((7.5 * 100.0 + 7.0 * 50.0) / 150.0)


def test_river_copies_differing_in_stage_alone_are_exact_so_nothing_warns(recwarn):
    """``sum C_i (s_i - max(h, r))`` is ``C (s_w - max(h, r))`` when the copies share
    their rbot, so the weighted record is MF6's at every head -- measured below at a head
    between the two stages and one below the bed (Review of S1, ledger 494)."""

    points, _ = _bc_frames(_two_records_on_one_cell("RIV", (9.0, 100.0, 7.0), (8.0, 50.0, 7.0)), "RIV", CRS)
    assert not [w for w in recwarn if "several records" in str(w.message)]
    assert points.attrs == {"records_combined": 1, "cells_inexact": 0, "cells_unweighted": 0}


def test_a_value_no_copy_disagrees_with_is_written_exactly():
    """The weighting touches only the cells whose copies DIFFER: a cell holding one
    record, and copies that agree, keep their value to the bit. Weighting every cell
    wrote 7.1 x 0.7 / 0.7 = 7.1000000000000005 for a lone drain's 7.1 (Review of S1,
    ledger 494)."""

    grid = square_grid()
    model = usg_model(grid, [0, 1])
    model.boundaries = {"DRN": PeriodRecords(ftype="DRN", periods={0: np.array([
        [1, 7.5, 100.0], [1, 7.0, 50.0],   # cell 0: copies that differ, weighted
        [2, 7.1, 0.7],                     # cell 1: one record
        [3, 7.1, 0.3], [3, 7.1, 1.1],      # cell 2: copies that agree
    ])})}
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        points, _ = _bc_frames(model, "DRN", CRS)
    table = points.set_index("old_cell")
    assert table.loc[1, "elevation"] == 7.1
    assert table.loc[2, "elevation"] == 7.1
    assert table.loc[0, "elevation"] == pytest.approx(22.0 / 3.0)


def test_copies_whose_conductances_total_zero_take_the_plain_mean_and_say_so():
    model = _two_records_on_one_cell("DRN", (7.5, 0.0), (7.0, 0.0))
    with pytest.warns(UserWarning, match=r"conductances total zero.*plain MEAN.*\[\(1, 0\)\]"):
        points, _ = _bc_frames(model, "DRN", CRS)
    row = points.set_index("old_cell").loc[0]
    assert row["conductance"] == 0.0 and row["elevation"] == 7.25
    assert points.attrs == {"records_combined": 1, "cells_inexact": 0, "cells_unweighted": 1}


def test_copies_that_agree_combine_without_a_warning(tmp_path, recwarn):
    manifest = ExportManifest(directory=tmp_path, crs=CRS, nper=1)
    _export_boundaries(_two_drains_on_one_cell(7.5), tmp_path, CRS, manifest)
    assert not [w for w in recwarn if "several records" in str(w.message)]
    assert manifest.layers["boundaries.gpkg:drn"]["records_combined"] == 1
    assert any("1 record(s) shared a cell" in note and "conductance-weighted" in note
               for note in manifest.notes)
    written = gpd.read_file(tmp_path / "boundaries.gpkg", layer="drn").set_index("old_cell")
    assert written.loc[0, "conductance"] == 150.0 and written.loc[0, "elevation"] == 7.5


# --- The exported one-record cell is MF6's two-record cell (ledger 493) --- #

def _mf6_flow(workspace, package, records, hold):
    """MF6's flow through ``package``'s ``records`` on one cell whose head a GHB of
    huge conductance holds at ``hold``."""

    import flopy

    sim = flopy.mf6.MFSimulation(sim_name="s", sim_ws=str(workspace), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, complexity="SIMPLE", inner_dvclose=1e-10, outer_dvclose=1e-10)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m", save_flows=True)
    flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=1, ncol=1, delr=10.0, delc=10.0, top=20.0, botm=0.0)
    flopy.mf6.ModflowGwfic(gwf, strt=hold)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfghb(gwf, stress_period_data={0: [[(0, 0, 0), hold, 1e9]]}, pname="hold")
    constructor = getattr(flopy.mf6, f"ModflowGwf{package}")
    constructor(gwf, stress_period_data={0: [[(0, 0, 0), *r] for r in records]}, pname="bc",
                filename=f"m.bc.{package}")
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord="m.cbc", saverecord=[("BUDGET", "ALL")])
    sim.write_simulation(silent=True)
    ok, _ = sim.run_simulation(silent=True)
    assert ok
    budget = flopy.utils.CellBudgetFile(str(workspace / "m.cbc"))
    records_of = [r for r in budget.get_data(text=package.upper())]
    # The holding GHB is the first GHB record; the boundary under test is the last.
    return float(records_of[-1]["q"].sum())


@pytest.mark.skipif(__import__("shutil").which("mf6") is None, reason="needs mf6 on PATH")
@pytest.mark.parametrize(
    "ftype, first, second, hold, column_order",
    [
        ("DRN", (7.5, 100.0), (7.0, 50.0), 8.0, ("elevation", "conductance")),
        ("GHB", (9.0, 100.0), (6.0, 50.0), 7.0, ("head", "conductance")),
        ("RIV", (9.0, 100.0, 7.5), (8.0, 50.0, 7.0), 8.5, ("stage", "conductance", "rbot")),
        ("RIV", (9.0, 100.0, 7.5), (8.0, 50.0, 7.0), 6.0, ("stage", "conductance", "rbot")),
        # One rbot, two stages: exact between the stages and below the bed alike
        # (ledger 494).
        ("RIV", (9.0, 100.0, 7.0), (8.0, 50.0, 7.0), 8.5, ("stage", "conductance", "rbot")),
        ("RIV", (9.0, 100.0, 7.0), (8.0, 50.0, 7.0), 6.0, ("stage", "conductance", "rbot")),
    ],
)
def test_the_exported_record_moves_what_mf6_s_two_records_move(tmp_path, ftype, first, second, hold, column_order):
    """Measured through MF6 at heads where both copies act: the one exported record's
    flow is the two records' flow. The plain mean it replaced drained 112.5 for 100.0
    (DRN) and leaked 0.0 for 25.0 (RIV at head 8.5)."""

    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        points, _ = _bc_frames(_two_records_on_one_cell(ftype, first, second), ftype, CRS)
    row = points.set_index("old_cell").loc[0]
    exported = [tuple(float(row[c]) for c in column_order)]
    package = ftype.lower()
    two = _mf6_flow(tmp_path / "two", package, [first, second], hold)
    one = _mf6_flow(tmp_path / "one", package, exported, hold)
    assert abs(two) > 1.0
    assert one == pytest.approx(two, rel=1e-6)

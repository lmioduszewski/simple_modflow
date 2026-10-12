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


def _two_drains_on_one_cell(second_elevation):
    grid = square_grid()
    model = usg_model(grid, [0, 1])
    model.boundaries = {"DRN": PeriodRecords(
        ftype="DRN",
        periods={0: np.array([[1, 7.5, 100.0], [1, second_elevation, 50.0], [2, 7.5, 100.0]])},
    )}
    return model


def test_two_records_on_one_cell_are_one_feature_combined_by_kind():
    """Two drains on one cell in one period -- MF6 applies both -- stopped the export
    on ``float()`` of a Series (ledger 467). Conductance adds; the elevation is
    the mean, and a warning names the cell whose copies differ."""

    with pytest.warns(UserWarning, match=r"DRN: 1 cell\(s\) .*elev values differ .*\[\(1, 0\)\]"):
        points, lines = _bc_frames(_two_drains_on_one_cell(7.0), "DRN", CRS)
    table = points.set_index("old_cell")
    assert sorted(table.index) == [0, 1]
    assert table.loc[0, "conductance"] == 150.0 and table.loc[0, "elevation"] == 7.25
    assert table.loc[1, "conductance"] == 100.0 and table.loc[1, "elevation"] == 7.5
    assert lines.set_index("old_cell").loc[0, "conductance"] == 150.0
    assert points.attrs == {"records_combined": 1, "cells_disagreeing": 1}


def test_copies_that_agree_combine_without_a_warning(tmp_path, recwarn):
    manifest = ExportManifest(directory=tmp_path, crs=CRS, nper=1)
    _export_boundaries(_two_drains_on_one_cell(7.5), tmp_path, CRS, manifest)
    assert not [w for w in recwarn if "several records" in str(w.message)]
    assert manifest.layers["boundaries.gpkg:drn"]["records_combined"] == 1
    assert any("1 record(s) shared a cell" in note for note in manifest.notes)
    written = gpd.read_file(tmp_path / "boundaries.gpkg", layer="drn").set_index("old_cell")
    assert written.loc[0, "conductance"] == 150.0 and written.loc[0, "elevation"] == 7.5

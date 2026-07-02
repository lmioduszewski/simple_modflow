"""Tests for Adaptive Time Stepping (ATS) support in ``TemporalDiscretization``.

Three layers:
  * pure-helper unit tests (period resolution, off-by-one, auto-derivation);
  * write-and-inspect integration tests that assert the emitted MF6 ``.tdis`` +
    ``.tdis.ats`` files are correct (right ``MAXATS``, right *1-based* periods);
  * a slow end-to-end test that runs MF6 and proves ATS actually subdivides the
    targeted period.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pytest

from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.modflow.mf6.simulation.discretization import (
    DisvGrid,
    TemporalDiscretization,
    _build_ats_records,
    _first_timestep_length,
    _resolve_ats_periods,
)


# --------------------------------------------------------------------------- #
# Fixtures / helpers                                                          #
# --------------------------------------------------------------------------- #
def _grid() -> VoronoiGridPlus:
    """A tiny 6-cell, single-layer Voronoi grid."""
    vertices = np.array(
        [
            [0, 0], [1, 0], [2, 0], [3, 0],
            [0, 1], [1, 1], [2, 1], [3, 1],
            [0, 2], [1, 2], [2, 2], [3, 2],
        ],
        dtype=float,
    )
    # clockwise winding: MF6 DISV requires cell2d vertices in clockwise order
    # (counter-clockwise yields "CELL2D area less than zero").
    iverts = [
        [4, 5, 1, 0], [5, 6, 2, 1], [6, 7, 3, 2],
        [8, 9, 5, 4], [9, 10, 6, 5], [10, 11, 7, 6],
    ]
    centers = np.array(
        [[0.5, 0.5], [1.5, 0.5], [2.5, 0.5], [0.5, 1.5], [1.5, 1.5], [2.5, 1.5]]
    )
    grid = VoronoiGridPlus(verts=vertices, iverts=iverts, xcyc=centers)
    grid.gdf_topbtm = gpd.GeoDataFrame(
        {0: [10.0] * 6, 1: [0.0] * 6},
        geometry=grid.gdf_vorPolys.geometry,
        crs=grid.crs,
    )
    return grid


def _writable_model(tmp_path, nper, period_data, **ats_kwargs):
    """Build a minimal, writable single-layer model with a TDIS (optionally ATS)."""
    grid = _grid()
    model = SimulationBase(
        name="atsmod", mf_folder_path=tmp_path, nper=nper, vor=grid
    )
    DisvGrid(vor=grid, model=model, top=[10.0] * 6, bottom=[0.0] * 6, nlay=1)
    import flopy

    flopy.mf6.ModflowGwfic(model.gwf, strt=8.0)
    flopy.mf6.ModflowGwfnpf(model.gwf, icelltype=1, k=1.0)
    td = TemporalDiscretization(
        model=model, period_data=period_data, **ats_kwargs
    )
    return model, td


def _ats_file_records(sim_ws, name="atsmod"):
    """Parse the written ``.tdis.ats`` file -> (maxats, [iperats ints])."""
    path = sim_ws / f"{name}.tdis.ats"
    text = path.read_text()
    maxats = None
    periods = []
    in_pd = False
    for line in text.splitlines():
        s = line.strip()
        if s.upper().startswith("MAXATS"):
            maxats = int(s.split()[1])
        elif s.upper().startswith("BEGIN PERIODDATA"):
            in_pd = True
        elif s.upper().startswith("END PERIODDATA"):
            in_pd = False
        elif in_pd and s and not s.startswith("#"):
            periods.append(int(float(s.split()[0])))
    return maxats, periods


# --------------------------------------------------------------------------- #
# Unit tests: pure helpers                                                    #
# --------------------------------------------------------------------------- #
def test_first_timestep_length_uniform_when_no_multiplier():
    assert _first_timestep_length(10.0, 5, 1.0) == pytest.approx(2.0)
    assert _first_timestep_length(1.0, 1, 1.1) == pytest.approx(1.0)


def test_first_timestep_length_matches_geometric_series():
    perlen, nstp, tsmult = 1.0, 10, 1.1
    dt0 = _first_timestep_length(perlen, nstp, tsmult)
    # the geometric series of steps must sum back to perlen
    total = sum(dt0 * tsmult ** i for i in range(nstp))
    assert total == pytest.approx(perlen)
    assert dt0 == pytest.approx(0.06274, abs=1e-4)


def test_resolve_ats_periods_true_selects_all():
    assert _resolve_ats_periods(True, 3) == {1: {}, 2: {}, 3: {}}


def test_resolve_ats_periods_none_and_false_select_nothing():
    assert _resolve_ats_periods(None, 5) == {}
    assert _resolve_ats_periods(False, 5) == {}


def test_resolve_ats_periods_iterable_and_mapping():
    assert _resolve_ats_periods([2, 4], 5) == {2: {}, 4: {}}
    resolved = _resolve_ats_periods({8: {"dtmin": 1e-4}}, 10)
    assert resolved == {8: {"dtmin": 1e-4}}


@pytest.mark.parametrize("bad", [0, 6, -1])
def test_resolve_ats_periods_out_of_range_raises(bad):
    with pytest.raises(ValueError, match="out of range"):
        _resolve_ats_periods([bad], 5)


def test_build_ats_records_offbyone_and_human_mapping():
    period_data = [[1.0, 10, 1.1]] * 4
    flopy_recs, human_recs = _build_ats_records(
        period_data, {2: {}, 4: {}},
        dt0=None, dtmin=None, dtmax=None, dtadj=2.0, dtfailadj=5.0,
    )
    # FloPy records carry the 0-based index (period-1); human carry 1-based.
    assert [r[0] for r in flopy_recs] == [1, 3]
    assert [r[0] for r in human_recs] == [2, 4]


def test_build_ats_records_autoderives_defaults():
    period_data = [[2.0, 4, 1.0]]  # perlen=2, nstp=4, tsmult=1 -> first step 0.5
    flopy_recs, _ = _build_ats_records(
        period_data, {1: {}},
        dt0=None, dtmin=None, dtmax=None, dtadj=2.0, dtfailadj=5.0,
    )
    _, dt0, dtmin, dtmax, dtadj, dtfailadj = flopy_recs[0]
    assert dt0 == pytest.approx(0.5)          # perlen/nstp
    assert dtmin == pytest.approx(2.0 * 1e-5)  # perlen * 1e-5
    assert dtmax == pytest.approx(2.0)         # perlen
    assert dtadj == 2.0 and dtfailadj == 5.0


def test_build_ats_records_scalar_and_per_period_overrides():
    period_data = [[1.0, 10, 1.1]] * 3
    flopy_recs, _ = _build_ats_records(
        period_data, {1: {}, 2: {"dtmin": 9e-3, "dtadj": 3.0}},
        dt0=0.02, dtmin=1e-3, dtmax=0.5, dtadj=2.0, dtfailadj=4.0,
    )
    # period 1: scalar defaults
    assert flopy_recs[0][1:] == [0.02, 1e-3, 0.5, 2.0, 4.0]
    # period 2: per-period overrides win over scalars
    assert flopy_recs[1][2] == 9e-3   # dtmin overridden
    assert flopy_recs[1][4] == 3.0    # dtadj overridden
    assert flopy_recs[1][1] == 0.02   # dt0 still from scalar


# --------------------------------------------------------------------------- #
# Integration tests: write + inspect the emitted MF6 files                    #
# --------------------------------------------------------------------------- #
def test_no_ats_leaves_tdis_plain(tmp_path):
    model, td = _writable_model(tmp_path, 3, [[1.0, 5, 1.1]] * 3)
    model.sim.write_simulation(silent=True)
    assert td.ats is None
    assert td.ats_perioddata is None
    tdis_text = (model.model_output_folder_path / "atsmod.tdis").read_text()
    assert "ATS6" not in tdis_text.upper()
    assert not (model.model_output_folder_path / "atsmod.tdis.ats").exists()


def test_ats_true_writes_every_period_with_correct_maxats(tmp_path):
    model, td = _writable_model(tmp_path, 4, [[1.0, 10, 1.1]] * 4, ats=True)
    model.sim.write_simulation(silent=True)
    tdis_text = (model.model_output_folder_path / "atsmod.tdis").read_text()
    assert "ATS6" in tdis_text.upper()
    maxats, periods = _ats_file_records(model.model_output_folder_path)
    assert maxats == 4
    assert periods == [1, 2, 3, 4]
    # human-facing record carries 1-based periods
    assert [r[0] for r in td.ats_perioddata] == [1, 2, 3, 4]


def test_ats_list_writes_only_selected_1based_periods(tmp_path):
    # The off-by-one, end to end: ask for MF6 periods 2 and 4.
    model, td = _writable_model(
        tmp_path, 5, [[1.0, 10, 1.1]] * 5, ats=[2, 4]
    )
    model.sim.write_simulation(silent=True)
    maxats, periods = _ats_file_records(model.model_output_folder_path)
    assert maxats == 2
    assert periods == [2, 4]
    # and the model exposes the resolved records
    assert model.ats_perioddata == td.ats_perioddata
    assert [r[0] for r in td.ats_perioddata] == [2, 4]


def test_ats_mapping_overrides_reach_the_file(tmp_path):
    model, td = _writable_model(
        tmp_path, 3, [[1.0, 10, 1.1]] * 3,
        ats={3: {"dtmin": 1.5e-4, "dtmax": 0.4}},
    )
    model.sim.write_simulation(silent=True)
    maxats, periods = _ats_file_records(model.model_output_folder_path)
    assert maxats == 1
    assert periods == [3]
    # dtmin/dtmax overrides captured in the human record
    _, _dt0, dtmin, dtmax, _adj, _fail = td.ats_perioddata[0]
    assert dtmin == pytest.approx(1.5e-4)
    assert dtmax == pytest.approx(0.4)


# --------------------------------------------------------------------------- #
# Slow: run MF6 and prove ATS subdivides the targeted period                  #
# --------------------------------------------------------------------------- #
def test_ats_subdivides_targeted_period_when_run(tmp_path):
    """MF6 honours ATS: the ATS period is forced to >=4 steps, others take 1."""
    import flopy
    import re

    grid = _grid()
    model = SimulationBase(
        name="atsrun", mf_folder_path=tmp_path, nper=3, vor=grid
    )
    DisvGrid(vor=grid, model=model, top=[10.0] * 6, bottom=[0.0] * 6, nlay=1)
    gwf = model.gwf
    # TDIS must exist before packages that key data by stress period (STO/WEL).
    # All periods 1 step; ATS on period 2 with dtmax=0.3 forces >=4 substeps.
    TemporalDiscretization(
        model=model,
        period_data=[[1.0, 1, 1.0]] * 3,
        ats={2: {"dt0": 0.3, "dtmin": 1e-4, "dtmax": 0.3, "dtadj": 1.0}},
    )
    flopy.mf6.ModflowGwfic(gwf, strt=8.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=1.0, save_flows=True)
    flopy.mf6.ModflowGwfsto(
        gwf, ss=1e-4, sy=0.1, iconvert=1, transient={0: True}
    )
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: [[(0, 0), 8.0]]})
    # a strong injection only during the (ATS) middle period
    flopy.mf6.ModflowGwfwel(
        gwf, stress_period_data={1: [[(0, 5), 20.0]], 2: []}
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord="atsrun.hds",
        saverecord=[("HEAD", "ALL")],
    )
    model.sim.write_simulation(silent=True)
    success, buff = model.sim.run_simulation(silent=True)
    assert success, "\n".join(buff[-25:])

    lst = (model.model_output_folder_path / "atsrun.lst").read_text()
    steps = re.findall(
        r"END OF TIME STEP\s+(\d+),\s+STRESS PERIOD\s+(\d+)", lst
    )
    per_counts: dict[int, int] = {}
    for step, per in steps:
        per_counts[int(per)] = max(per_counts.get(int(per), 0), int(step))
    assert per_counts.get(1) == 1, per_counts   # non-ATS period: 1 step
    assert per_counts.get(3) == 1, per_counts   # non-ATS period: 1 step
    assert per_counts.get(2, 0) >= 4, per_counts  # ATS period subdivided

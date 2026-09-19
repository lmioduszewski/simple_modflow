"""`min_above_bottom=` keeps an elevation inside its cell, and the detection is default-on.

**What MF6 does with a below-bottom elevation depends on `icelltype`, and both branches
are live in this library** (ledger 290). At `icelltype = 0` -- what `mf.npf(k=...)` writes,
because it names no icelltype at all -- all four packages ACCEPT one, terminate normally
and over-drain. At `icelltype != 0` -- what `CanonicalModelConfig` and 17 taught sites use
-- MF6 REFUSES by name, and the refusal is a per-period data read, so an offender in
period 3 dies after MF6 has solved periods 1 and 2.

Nothing reported it before: no warning, no metadata key, no validator. The silence is what
is constant across models; the magnitude is not.
"""

from __future__ import annotations

import json
import warnings

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import LineString, Point, box

import myflopy as mf
from myflopy.specs import package_digest


@pytest.fixture(scope="module")
def sloping():
    """A 2-layer context whose layer-1 bottom SLOPES, so one flat elevation cannot suit it.

    A flat base would make a scalar floor indistinguishable from a per-cell one, which is
    the whole distinction under test: `CellSurfaceOffset.minimum` is per FEATURE and
    cannot track a bottom that moves across the cells one feature covers.
    """

    import tempfile
    from pathlib import Path

    workspace = Path(tempfile.mkdtemp())
    boundary = workspace / "b.gpkg"
    gpd.GeoDataFrame(
        {"name": ["d"]}, geometry=[box(0.0, 0.0, 600.0, 600.0)], crs=2927
    ).to_file(boundary)
    grid = mf.GridSpec.voronoi(
        boundary=mf.GeoPackageSourceSpec(path=boundary), max_area=30_000.0
    ).resolve(workspace=workspace / "g")

    ncpl = int(grid.ncpl)
    xs = np.asarray(grid.gdf_vorPolys.geometry.centroid.x, dtype=float)
    top = np.full(ncpl, 100.0)
    botm = np.stack([60.0 - 0.05 * xs, np.zeros(ncpl)])
    context = mf.ModelContext(
        grid=grid, domain=np.ones((2, ncpl), dtype=int),
        surfaces={"top": top, "botm": botm},
    )
    drains = workspace / "d.gpkg"
    gpd.GeoDataFrame(
        {"name": ["drn"], "layer": [1], "elev": [45.0], "conductance": [100.0]},
        geometry=[LineString([(20.0, 300.0), (580.0, 300.0)])], crs=2927,
    ).to_file(drains)
    return context, drains, botm[0]


def _elevations(spec):
    records = spec.options["stress_period_data"][0]
    return np.array([r[0][1] for r in records]), np.array([float(r[1]) for r in records])


# --- detection, which is default-on ------------------------------------------


def test_a_below_bottom_elevation_is_reported_without_being_asked(sloping):
    """The half that reaches a model nobody thought to floor.

    The keyword alone only helps someone who already knows they have the problem; at
    `icelltype = 0` nothing would ever tell them.
    """

    context, drains, bottoms = sloping

    with pytest.warns(UserWarning, match="BELOW their own cell bottom") as caught:
        spec = mf.drn.gpkg(
            drains, context=context, nper=1, elevation="elev", conductance="conductance"
        )

    cells, elevations = _elevations(spec)
    assert (elevations < bottoms[cells]).sum() > 0, "the fixture must actually offend"
    assert spec.metadata["coverage"]["below_bottom"] == int((elevations < bottoms[cells]).sum())
    assert spec.metadata["coverage"]["worst_below_bottom"] < 0
    assert "min_above_bottom" in str(caught[0].message)


def test_the_detection_says_nothing_when_every_elevation_is_inside_its_cell(sloping):
    """The refusing direction: a clean declaration must stay silent, or it is noise."""

    context, _, _ = sloping
    path = _one_line(context, elev=None)

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        mf.drn.gpkg(
            path, context=context, nper=1, conductance="conductance",
            elevation=mf.CellSurfaceOffset("cell_bottom", offset=5.0),
        )


def _one_line(context, *, elev):
    import tempfile
    from pathlib import Path

    path = Path(tempfile.mkdtemp()) / "one.gpkg"
    columns = {"name": ["drn"], "layer": [1], "conductance": [100.0]}
    if elev is not None:
        columns["elev"] = [elev]
    gpd.GeoDataFrame(
        columns, geometry=[LineString([(20.0, 300.0), (580.0, 300.0)])], crs=2927
    ).to_file(path)
    return path


# --- the floor ---------------------------------------------------------------


def test_the_floor_is_PER_CELL_where_minimum_is_per_feature(sloping):
    """This is the capability that did not exist.

    `CellSurfaceOffset("cell_bottom", offset=F, minimum=<declared>)` is exact sugar for
    `min_above_bottom=F` when the declared elevation is a constant or a single column --
    measured, 0 of 7 mismatches. What it cannot do is floor a value that came from
    somewhere else: `minimum=` is read once per feature, so one number serves every cell
    a feature covers, and `minimum=<CellSurfaceOffset>` raises `TypeError`.
    """

    context, drains, bottoms = sloping

    with pytest.warns(UserWarning):
        floored = mf.drn.gpkg(
            drains, context=context, nper=1, elevation="elev",
            conductance="conductance", min_above_bottom=1.0,
        )

    cells, elevations = _elevations(floored)
    assert (elevations < bottoms[cells]).sum() == 0
    assert np.all(elevations >= bottoms[cells] + 1.0 - 1e-9)
    # Per cell: the lifted values differ from one another, which one scalar cannot do.
    assert len(set(np.round(elevations, 6))) > 1
    assert floored.metadata["coverage"]["floored"] > 0


def test_a_zero_clearance_still_lands_STRICTLY_above_the_bottom(sloping):
    """`min_above_bottom=0.0` cannot mean "exactly at the bottom".

    FloPy writes `botm` fixed-point and a list-BC value in E format, so a value set
    exactly equal to its cell bottom can round BELOW it in the written file and MF6's
    strict comparison trips -- measured, 13 of 24 real runs FAILED at bottoms drawn from
    [40, 600], and which way the discarded digits fall is a coin flip. A CHD set exactly
    at its bottom is worse: it aborts the SOLVER, at round magnitudes too. So the
    clearance is nudged relative to the magnitude.
    """

    context, drains, bottoms = sloping

    with pytest.warns(UserWarning):
        spec = mf.drn.gpkg(
            drains, context=context, nper=1, elevation="elev",
            conductance="conductance", min_above_bottom=0.0,
        )

    cells, elevations = _elevations(spec)
    lifted = elevations > 45.0
    assert lifted.any()
    assert np.all(elevations[lifted] > bottoms[cells][lifted]), "strictly above, never equal"


def test_a_negative_clearance_is_refused(sloping):
    context, drains, _ = sloping

    with pytest.raises(ValueError, match="CLEARANCE"):
        mf.drn.gpkg(
            drains, context=context, nper=1, elevation="elev",
            conductance="conductance", min_above_bottom=-1.0,
        )


# --- RIV is a PAIR ----------------------------------------------------------


def test_riv_lifts_the_pair_so_the_declared_river_depth_survives(sloping):
    """Where the STAGE is below the floor too, it is carried up with the bed.

    MF6 checks `rbot` and enforces `stage >= rbot` separately, so raising `rbot` alone
    turns a model that ran into a hard failure, and clamping both collapses the depth.
    Measured on one record: as declared q=-868.57; rbot alone REFUSED at both
    icelltypes; both clamped q=-380.00 with the river depth collapsed to 0.000 from a
    declared 1.0; pair-lifted q=-108.57.

    The bottom SLOPES here, so one declaration meets the floor differently in each cell
    -- which is what separates this from its sibling below.
    """

    context, _, bottoms = sloping
    import tempfile
    from pathlib import Path

    path = Path(tempfile.mkdtemp()) / "r.gpkg"
    gpd.GeoDataFrame(
        {"name": ["riv"], "layer": [1], "stage": [40.0], "rbot": [35.0],
         "conductance": [100.0]},
        geometry=[LineString([(20.0, 300.0), (580.0, 300.0)])], crs=2927,
    ).to_file(path)

    with pytest.warns(UserWarning):
        spec = mf.riv.gpkg(
            path, context=context, nper=1, stage="stage", rbot="rbot",
            conductance="conductance", min_above_bottom=1.0,
        )

    records = spec.options["stress_period_data"][0]
    cells = np.array([r[0][1] for r in records])
    stage = np.array([float(r[1]) for r in records])
    rbot = np.array([float(r[3]) for r in records])

    assert (rbot < bottoms[cells]).sum() == 0
    assert np.all(stage >= rbot)
    moved = stage > 40.0
    assert moved.any(), "the fixture must contain a stage that had to move"
    assert np.allclose((stage - rbot)[moved], 5.0), "its declared 5 ft depth must survive"
    assert np.allclose(stage[~moved], 40.0), "and a stage that did not have to move stays"


def test_a_river_whose_stage_is_already_inside_its_cell_keeps_it(sloping):
    """The half ledger 290 generalised past: a bed below the cell, a legal stage.

    That is what an ordinary river looks like, and the stage is the water surface it
    imposes -- lifting it by the bed's delta moves the answer. Measured on a one-cell
    bench (bottom 48, top 60), declared stage 55 / rbot 40: lifting the pair wrote
    stage 63, three feet ABOVE land surface, and turned q = 0.00 into +436.37; with a
    deeper bed (stage 50 / rbot 20) it reversed the sign, -272.73 to +1254.55. Lifting
    the bed alone reproduced the declared answer exactly and ran at both icelltypes.
    """

    context, _, bottoms = sloping
    import tempfile
    from pathlib import Path

    path = Path(tempfile.mkdtemp()) / "ordinary.gpkg"
    gpd.GeoDataFrame(
        {"name": ["riv"], "layer": [1], "stage": [75.0], "rbot": [28.0],
         "conductance": [100.0]},
        geometry=[LineString([(20.0, 300.0), (580.0, 300.0)])], crs=2927,
    ).to_file(path)

    with pytest.warns(UserWarning):
        spec = mf.riv.gpkg(
            path, context=context, nper=1, stage="stage", rbot="rbot",
            conductance="conductance", min_above_bottom=1.0,
        )

    records = spec.options["stress_period_data"][0]
    cells = np.array([r[0][1] for r in records])
    stage = np.array([float(r[1]) for r in records])
    rbot = np.array([float(r[3]) for r in records])

    assert (rbot < bottoms[cells]).sum() == 0, "every bed still comes inside its cell"
    assert np.allclose(stage, 75.0), "and not one stage moved"


def test_a_river_stage_below_its_own_bottom_is_refused(sloping):
    """MF6 refuses `stage < rbot` at BOTH icelltypes, and `mf.riv.gpkg`'s own docstring
    already promised the invariant while nothing enforced it."""

    context, _, _ = sloping
    import tempfile
    from pathlib import Path

    path = Path(tempfile.mkdtemp()) / "bad.gpkg"
    gpd.GeoDataFrame(
        {"name": ["riv"], "layer": [1], "stage": [30.0], "rbot": [40.0],
         "conductance": [100.0]},
        geometry=[LineString([(20.0, 300.0), (580.0, 300.0)])], crs=2927,
    ).to_file(path)

    with pytest.raises(ValueError, match="river stage .* BELOW its river bottom"):
        mf.riv.gpkg(
            path, context=context, nper=1, stage="stage", rbot="rbot",
            conductance="conductance",
        )


# --- where it is and is not offered -----------------------------------------


@pytest.mark.parametrize("helper_name", ["wel", "rch", "evt"])
def test_it_is_not_offered_where_no_field_is_an_elevation(sloping, helper_name):
    """`wel.rate`, `rch.recharge` and `evt`'s rate/depth are not elevations. `evt.surface`
    IS one, but the user excluded evt from the floor deliberately."""

    _, drains, _ = sloping

    # Declared LAZILY so the option check runs at spec construction rather than after the
    # file's required columns are read -- the point is that the keyword is not offered,
    # not that this file lacks a `rate` column.
    with pytest.raises(TypeError, match="does not accept 'min_above_bottom'"):
        getattr(mf, helper_name).gpkg(drains, nper=1, min_above_bottom=1.0)

    # And the three must still WORK without it: making it required on the shared lazy
    # builder broke them outright, which no test in this file would otherwise have seen.
    assert getattr(mf, helper_name).gpkg(drains, nper=1).inputs


def test_the_floor_is_applied_AFTER_a_nearest_snap(sloping):
    """A snapped point lands in a DIFFERENT cell with a different bottom, so a floor
    applied before the snap would use the wrong one (ledger 288 added the snap)."""

    context, _, bottoms = sloping
    import tempfile
    from pathlib import Path

    grid = context.grid
    xmax = float(grid.gdf_vorPolys.total_bounds[2])
    path = Path(tempfile.mkdtemp()) / "near.gpkg"
    gpd.GeoDataFrame(
        {"name": ["just_out"], "layer": [1], "elev": [-500.0], "conductance": [100.0]},
        geometry=[Point(xmax + 3.0, 300.0)], crs=2927,
    ).to_file(path)

    with pytest.warns(UserWarning):
        spec = mf.drn.gpkg(
            path, context=context, nper=1, elevation="elev", conductance="conductance",
            gaps="nearest", min_above_bottom=1.0,
        )

    cells, elevations = _elevations(spec)
    assert len(cells) == 1
    # Floored against the bottom of the cell it SNAPPED INTO.
    assert elevations[0] == pytest.approx(bottoms[cells][0] + 1.0)


# --- it survives being written down -----------------------------------------


def test_the_recipe_carries_the_clearance_and_round_trips(sloping):
    _, drains, _ = sloping
    spec = mf.drn.gpkg(drains, nper=1, elevation="elev", min_above_bottom=2.5)

    assert spec.inputs["stress_period_data"].min_above_bottom == 2.5
    revived = type(spec).from_dict(json.loads(json.dumps(spec.to_dict())))
    assert revived.inputs["stress_period_data"].min_above_bottom == 2.5
    assert package_digest(revived) == package_digest(spec)


def test_a_defaulted_clearance_moves_no_digest(sloping):
    """It rides in a recipe whose `to_dict` omits defaults, so adding the field cannot
    report a stored package as drifted."""

    _, drains, _ = sloping
    without = mf.drn.gpkg(drains, nper=1, elevation="elev")

    assert "min_above_bottom" not in without.to_dict()["inputs"]["stress_period_data"]

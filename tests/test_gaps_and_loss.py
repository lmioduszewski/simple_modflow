"""`gaps=` says what a feature that reached no cell means, and the loss is accounted for.

A feature that landed nowhere was reported through `logger.warning` -- and `import
myflopy` installs a `NullHandler` on the `myflopy` logger while root has none, so
`Logger.callHandlers` counts a handler as FOUND and never reaches `logging.lastResort`.
Measured: zero bytes on stderr in a bare session and in a real ipykernel. The motivating
case is two wells with one off the grid: **90% of the declared pumping gone, no python
warning** (ledger 288).

`nearest` moves POINTS only, which is measured rather than cautious -- see
`test_nearest_does_not_move_a_feature_with_extent`. An areal feature is CLIPPED to the
part inside the grid, which `cell_shares` already did exactly; what was missing is that
the clipped-away part was never reported.
"""

from __future__ import annotations

import json
import warnings

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import LineString, Point, box

import myflopy as mf
from myflopy.geopackage import Spread
from myflopy.specs import package_digest
from tests.test_gpkg_on_cell_faces import _cells, _context

# The 4x4 grid of `_context()` spans (0, 0)-(100, 100) in 25 ft cells.
JUST_OUTSIDE = Point(103.0, 62.5)       # 3 ft east of the grid, beside cell 11
FAR = Point(5_000.0, 5_000.0)


def _wells(tmp_path, *features, names, rate=-500.0):
    path = tmp_path / "w.gpkg"
    gpd.GeoDataFrame(
        {"name": list(names), "layer": [1] * len(features), "q": [rate] * len(features)},
        geometry=list(features),
        crs="EPSG:2927",
    ).to_file(path, layer="bc")
    return path


def _total(spec):
    return sum(float(record[1]) for record in spec.options["stress_period_data"][0])


# --- the four words ----------------------------------------------------------


def test_the_default_warns_and_says_what_was_not_written(tmp_path):
    """The motivating case: one of two wells off the grid, 90% of the pumping gone."""

    path = _wells(tmp_path, Point(62.5, 62.5), FAR, names=["on_grid", "off_grid"])

    with pytest.warns(UserWarning, match="1 of 2 features") as caught:
        spec = mf.wel.gpkg(path, layer="bc", context=_context(), nper=1, rate="q")

    assert _total(spec) == pytest.approx(-500.0)
    assert "rate=q -500" in str(caught[0].message)


def test_raise_refuses_the_same_thing_the_default_warns_about(tmp_path):
    path = _wells(tmp_path, Point(62.5, 62.5), FAR, names=["on_grid", "off_grid"])

    with pytest.raises(ValueError, match="1 of 2 features"):
        mf.wel.gpkg(path, layer="bc", context=_context(), nper=1, rate="q", gaps="raise")


def test_ignore_is_silent_but_still_records_what_happened(tmp_path):
    """Silent is not the same as unrecorded: a warning is best-effort -- it dedupes, and
    on a lazy declaration it fires inside `prepare_run` where no stacklevel points at
    the caller -- so the numbers are stamped into metadata as well."""

    path = _wells(tmp_path, Point(62.5, 62.5), FAR, names=["on_grid", "off_grid"])

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        spec = mf.wel.gpkg(path, layer="bc", context=_context(), nper=1, rate="q", gaps="ignore")

    assert spec.metadata["coverage"]["missed"] == 1
    assert spec.metadata["coverage"]["lost"] == {"rate=q": pytest.approx(-500.0)}


@pytest.mark.parametrize("gaps", ["warn", "ignore", "nearest", "raise"])
def test_every_mode_still_refuses_a_source_that_lands_nowhere(tmp_path, gaps):
    """The all-missed refusal is ORTHOGONAL to the mode.

    The raster sibling's `"ignore"` returns BEFORE its all-missing raise, because
    `Surface.sample` has a documented next stage that fills the gaps. A boundary package
    has no such caller: MF6 rejects an empty list package naming the file, or runs
    without the boundary under `maxbound=`. Reversing ledger 211's central decision in
    one word is what this refuses to do.
    """

    path = _wells(tmp_path, FAR, names=["gone"])

    with pytest.raises(ValueError, match="lands on an active cell"):
        mf.wel.gpkg(path, layer="bc", context=_context(), nper=1, rate="q", gaps=gaps)


def test_an_unknown_mode_is_refused_against_the_one_vocabulary(tmp_path):
    path = _wells(tmp_path, Point(62.5, 62.5), names=["on_grid"])

    with pytest.raises(ValueError, match="unknown gaps mode"):
        mf.wel.gpkg(path, layer="bc", context=_context(), nper=1, rate="q", gaps="nearset")


# --- nearest -----------------------------------------------------------------


def test_nearest_moves_a_point_that_only_just_missed(tmp_path):
    path = _wells(tmp_path, Point(62.5, 62.5), JUST_OUTSIDE, names=["on_grid", "just_outside"])

    with pytest.warns(UserWarning, match="MOVED"):
        spec = mf.wel.gpkg(
            path, layer="bc", context=_context(), nper=1, rate="q", gaps="nearest"
        )

    assert len(spec.options["stress_period_data"][0]) == 2
    assert spec.metadata["coverage"]["snapped"] == 1
    assert spec.metadata["coverage"]["snap_max"] == pytest.approx(3.0)


def test_snap_within_refuses_a_point_beyond_the_tolerance(tmp_path):
    """The tolerance is the discriminator, and it is the ONLY one that works.

    Snap-vector bearing coherence cannot gate this: legitimate one-sided digitizing slop
    measures R = 1.000, identical to a uniform 700 ft datum shift. Magnitude separates
    them -- one cell width accepts 12 of 12 legitimate slop and refuses 77 of 78
    systematic misses (ledger 288).
    """

    path = _wells(tmp_path, Point(62.5, 62.5), JUST_OUTSIDE, names=["on_grid", "just_outside"])

    with pytest.warns(UserWarning, match="land on no active cell"):
        spec = mf.wel.gpkg(
            path, layer="bc", context=_context(), nper=1, rate="q",
            gaps="nearest", snap_within=1.0,
        )

    assert len(spec.options["stress_period_data"][0]) == 1


@pytest.mark.parametrize(
    "feature",
    [
        # PARALLEL is the case the guard actually decides: translating this one back
        # onto the grid DOES resolve -- measured, 3 cells summing to exactly 1.0000 --
        # so without the points-only rule it would be snapped, and whether a line is
        # recovered would depend on which way it happens to point. A line pointing AWAY
        # resolves to zero cells either way, which is why testing only that one passes
        # with the whole rule deleted.
        pytest.param(LineString([(103.0, 30.0), (103.0, 70.0)]), id="line-PARALLEL-to-edge"),
        pytest.param(LineString([(103.0, 62.5), (203.0, 62.5)]), id="line-pointing-away"),
        pytest.param(box(103.0, 40.0, 150.0, 80.0), id="polygon"),
    ],
)
def test_nearest_does_not_move_a_feature_with_extent(tmp_path, feature):
    """Points only, and the measurement is why.

    Translating by the minimal vector recovers a POINT exactly. It recovers a LINE only
    when the line happens to lie PARALLEL to the boundary it missed -- the same line
    POINTING at the grid lands with its tip touching and resolves to ZERO cells. A
    polygon translated by exactly its gap only touches; with 0.01 ft of overshoot the
    shares sum to 0.0001, so a `Spread` conductance would write 0.01% of what was
    declared. Orientation luck is not a rule (ledger 288).

    An areal feature that PARTLY overlaps is a different thing entirely and is not a
    miss at all -- it is clipped, and reported as clipped.
    """

    path = tmp_path / "f.gpkg"
    gpd.GeoDataFrame(
        {"name": ["inside", "outside"], "layer": [1, 1], "q": [-500.0, -500.0]},
        geometry=[Point(62.5, 62.5), feature], crs="EPSG:2927",
    ).to_file(path, layer="bc")

    with pytest.warns(UserWarning, match="land on no active cell"):
        spec = mf.wel.gpkg(
            path, layer="bc", context=_context(), nper=1, rate="q", gaps="nearest"
        )

    assert spec.metadata["coverage"]["snapped"] == 0
    assert spec.metadata["coverage"]["missed"] == 1


# --- clipping, which is what an areal feature gets instead --------------------


def test_a_polygon_partly_off_the_grid_is_clipped_and_the_loss_is_named(tmp_path):
    """`cell_shares` already clipped exactly right -- `sum(share)` equals the true
    in-grid fraction. What was missing is that the discarded part was never reported:
    a half-off polygon under `Spread(mode="clip")` wrote 497.1 of a declared 1000.0
    with no warning at any level."""

    path = tmp_path / "z.gpkg"
    gpd.GeoDataFrame(
        {"name": ["half_off"], "layer": [1], "conductance": [1000.0], "elevation": [5.0]},
        geometry=[box(75.0, 25.0, 125.0, 75.0)], crs="EPSG:2927",
    ).to_file(path, layer="bc")

    with pytest.warns(UserWarning, match="CLIPPED") as caught:
        spec = mf.drn.gpkg(
            path, layer="bc", context=_context(), nper=1,
            conductance=Spread("conductance"), elevation="elevation",
        )

    assert spec.metadata["coverage"]["clipped"] == 1
    assert spec.metadata["coverage"]["lost"]["conductance=conductance"] == pytest.approx(500.0)
    assert "50.0% inside" in str(caught[0].message)


def test_an_intensive_field_is_clipped_without_a_meaningless_total(tmp_path):
    """A conductance NOT wrapped in `Spread` is copied whole to every covered cell, by
    the documented rule -- so declared-minus-written would report a clean build as a
    -300% "loss". Ledger 169 measured `corr(cond, cell area) ~ 0` over 205 real records:
    a conductance is not a leakance. The clip is reported; no total is."""

    path = tmp_path / "z.gpkg"
    gpd.GeoDataFrame(
        {"name": ["half_off"], "layer": [1], "conductance": [1000.0], "elevation": [5.0]},
        geometry=[box(75.0, 25.0, 125.0, 75.0)], crs="EPSG:2927",
    ).to_file(path, layer="bc")

    with pytest.warns(UserWarning, match="CLIPPED"):
        spec = mf.drn.gpkg(path, layer="bc", context=_context(), nper=1)

    assert spec.metadata["coverage"]["clipped"] == 1
    assert "lost" not in spec.metadata["coverage"]


def test_a_source_that_lands_entirely_inside_says_nothing(tmp_path):
    """The refusing direction: a clean build must stay silent, or the warning is noise."""

    path = _wells(tmp_path, Point(62.5, 62.5), names=["on_grid"])

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        spec = mf.wel.gpkg(path, layer="bc", context=_context(), nper=1, rate="q")

    assert _cells(spec) == [10]


# --- it survives being written down ------------------------------------------


def test_the_recipe_carries_the_mode_and_round_trips(tmp_path):
    """A lazy declaration has no grid, so `snap_within` stays None in the recipe and is
    resolved at build: a number resolved at declaration would carry the OLD mesh's cell
    width across `with_grid`."""

    path = _wells(tmp_path, Point(62.5, 62.5), names=["on_grid"])
    spec = mf.wel.gpkg(path, layer="bc", nper=1, rate="q", gaps="nearest")

    recipe = spec.inputs["stress_period_data"]
    assert recipe.gaps == "nearest"
    assert recipe.snap_within is None

    revived = type(spec).from_dict(json.loads(json.dumps(spec.to_dict())))
    assert revived.inputs["stress_period_data"].gaps == "nearest"
    assert package_digest(revived) == package_digest(spec)


def test_a_numpy_snap_within_can_still_be_written_down(tmp_path):
    """`snap_within=np.sqrt(vor.gdf_vorPolys.area.median())` is the obvious spelling.

    `BoundarySource.to_dict` compared each field against its default WITHOUT
    normalizing, so a numpy-typed value went into the payload and it then failed both
    `json.dumps` and `yaml.safe_dump` -- a recipe that could not be saved, frozen or
    round-tripped (ledger 288).
    """

    path = _wells(tmp_path, Point(62.5, 62.5), names=["on_grid"])
    spec = mf.wel.gpkg(
        path, layer="bc", nper=1, rate="q", gaps="nearest", snap_within=np.float64(4.0)
    )

    payload = spec.to_dict()
    json.dumps(payload)
    assert type(payload["inputs"]["stress_period_data"]["snap_within"]) is float


def test_the_coverage_record_does_not_move_a_digest(tmp_path):
    """It rides in `metadata`, which `package_digest` excludes -- so recording what
    happened cannot report a stored run as drifted on the day it lands."""

    path = _wells(tmp_path, Point(62.5, 62.5), FAR, names=["on_grid", "off_grid"])
    context = _context()

    with pytest.warns(UserWarning):
        warned = mf.wel.gpkg(path, layer="bc", context=context, nper=1, rate="q")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        quiet = mf.wel.gpkg(path, layer="bc", context=context, nper=1, rate="q", gaps="ignore")

    assert warned.metadata["coverage"] == quiet.metadata["coverage"]
    assert package_digest(warned) == package_digest(quiet)

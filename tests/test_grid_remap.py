"""A per-cell array is CARRIED onto a new mesh by area overlay (ledger 300).

B3's C6. Before this, a model whose arrays came from nowhere reproducible -- a
PEST-calibrated K field, an array read out of an imported model -- could not move at all:
`with_grid` stamped it `unverifiable` and the build refused it, because an array names its
cells by POSITION and every array is still the right length on a mesh of the same size.

Two things this file is careful about, both measured:

* **A remap is the WORSE route where a source exists.** Re-declaring a polygon map on the
  new mesh (ledger 299) held to 7.1% worst-case flux error against this engine's 12.5% at
  best. So `remap=` is opt-in, the default is unchanged, and the warning says so.
* **K has no safe default**, because the right average depends on the FLOW DIRECTION and
  not on the field. The engine therefore reports how far the candidate rules are apart
  rather than pretending the choice is settled.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

import myflopy as mf
from myflopy.grid_identity import CELLWISE_OPTIONS, MOVED_FROM, REMAPPED_FROM, package_move
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, rectangular_voronoi
from myflopy.modflow.mf6.grid.geometry import (
    REMAP_RULES,
    cell_overlay,
    refuse_a_different_frame,
)
from myflopy.modflow.mf6.grid.remap import RemapRefused, remap_arrays
from myflopy.specs import ModelContext, package_digest


def _grid(n, cell=100.0):
    return rectangular_voronoi(
        CanonicalModelConfig(nrow=n, ncol=n, nlay=1, nper=1, cell_size=cell)
    )


@pytest.fixture(scope="module")
def meshes():
    """A fine mesh and a strictly coarser one over the SAME footprint."""

    return _grid(12), _grid(6, cell=200.0)


@pytest.fixture(scope="module")
def lens(meshes):
    """A strong K contrast: a clay lens inside sand. The case the rules disagree on."""

    fine, _ = meshes
    values = np.full(fine.ncpl, 10.0)
    centres = fine.gdf_vorPolys.geometry.centroid
    inside = (
        (centres.x > 400) & (centres.x < 800) & (centres.y > 400) & (centres.y < 800)
    ).to_numpy()
    values[inside] = 1e-3
    assert inside.sum() > 4, "the fixture must actually contain a lens"
    return values


def _stack():
    stack = mf.LayerStack(top=mf.Surface.flat(100.0))
    stack.add("sand", thickness=50.0)
    return stack


def _sim(grid, *, k, strt=40.0):
    return mf.simulation(
        mf.gwf("m", context=ModelContext(grid=grid),
               packages=[mf.disv(_stack()), mf.npf(k=k), mf.ic(strt=strt)]),
        name="s", packages=[mf.tdis(nper=1)],
    )


def _verdict(sim, package):
    return next(row["how"] for row in sim.move_report() if row["package"] == package)


# ------------------------------------------------------------------------- the point


def test_an_array_package_moves_and_builds_only_with_remap(meshes, lens, tmp_path):
    """`unverifiable` -> `remapped`, and the build stops refusing it."""

    fine, coarse = meshes
    sim = _sim(fine, k=lens, strt=np.full(fine.ncpl, 40.0))

    assert _verdict(sim.with_grid(coarse), "npf") == "unverifiable"

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        moved = sim.with_grid(coarse, remap=True)
    assert _verdict(moved, "npf") == "remapped"
    assert _verdict(moved, "ic") == "remapped"

    built = moved.build(tmp_path / "remapped")
    written = built.simulation.get_model("m").get_package("npf").k.get_data()
    assert written.shape == (1, coarse.ncpl)
    assert np.all(np.isfinite(written))
    # Bracketed by the source, because every rule here is an average of values present --
    # with a float tolerance, because the geometric mean of a uniform region is
    # `exp(sum(a*log(v))/sum(a))` and came back 10.000000000000007 for a source 10.0.
    assert written.min() >= lens.min() * (1 - 1e-9)
    assert written.max() <= lens.max() * (1 + 1e-9)


def test_remap_is_off_by_default(meshes, lens):
    """The existing behaviour must not move: a remap is measurably the inferior route."""

    fine, coarse = meshes
    moved = _sim(fine, k=lens).with_grid(coarse)
    assert _verdict(moved, "npf") == "unverifiable"
    assert REMAPPED_FROM not in (moved.model("m").packages[1].metadata or {})


def test_a_constant_field_survives_exactly(meshes):
    """Nothing about averaging may disturb a field that is already uniform."""

    fine, coarse = meshes
    overlay = cell_overlay(fine.gdf_vorPolys, coarse.gdf_vorPolys)
    for rule in REMAP_RULES:
        painted, reached = overlay.reduce(np.full(fine.ncpl, 7.0), rule=rule)
        assert reached.all()
        assert np.allclose(painted, 7.0, rtol=0, atol=1e-12), rule


def test_a_partially_covered_cell_uses_the_COVERED_area(meshes):
    """The denominator that makes a partly covered cell right rather than diluted.

    With the new cell's OWN area as the denominator a missing weight acts as a value of
    zero, and a constant 7.0 field came back between 0.0066 and 6.9878.
    """

    fine, _ = meshes
    grown = _grid(8, cell=200.0)          # 1,600 wide against the source's 1,200
    overlay = cell_overlay(fine.gdf_vorPolys, grown.gdf_vorPolys)
    painted, reached = overlay.reduce(np.full(fine.ncpl, 7.0), rule="mean")
    assert reached.sum() < grown.ncpl, "this fixture must have uncovered cells to be a test"
    assert np.allclose(painted[reached], 7.0, rtol=0, atol=1e-12)


# -------------------------------------------------------------------------- the rules


def test_the_rules_are_ordered_and_distinct(meshes, lens):
    """harmonic <= geometric <= mean, per cell, and far apart on a real contrast.

    The ordering is a mathematical fact, so it is the cheapest available proof that the
    three branches are three formulas rather than one answering for all of them.
    """

    fine, coarse = meshes
    overlay = cell_overlay(fine.gdf_vorPolys, coarse.gdf_vorPolys)
    mean, _ = overlay.reduce(lens, rule="mean")
    geometric, _ = overlay.reduce(lens, rule="geometric")
    harmonic, _ = overlay.reduce(lens, rule="harmonic")

    assert np.all(harmonic <= geometric + 1e-12)
    assert np.all(geometric <= mean + 1e-12)
    assert np.nanmax(mean / harmonic) > 100.0, "the fixture must make the choice matter"

    # STRICT somewhere, on the cells that actually hold both materials. The `<=` above is
    # satisfied by equality, so on its own it cannot tell one branch answering for two
    # from three formulas -- a mutation making `harmonic` return the geometric mean
    # passed every other test in this file.
    # A RELATIVE threshold, not `!= 0`: on a cell holding one material the three agree to
    # float noise (the harmonic mean of identical 10.0s came back 9.999999999999998), and
    # selecting those made the strict comparison below fail on correct code.
    mixed = mean > harmonic * 1.01
    assert mixed.any(), "no cell holds both materials -- the fixture is not a contrast"
    assert np.all(harmonic[mixed] < geometric[mixed])
    assert np.all(geometric[mixed] < mean[mixed])


def test_the_spread_report_finds_the_cells_where_the_choice_matters(meshes, lens):
    """What makes a conductivity default honest: the engine says where it is a real choice."""

    fine, coarse = meshes
    overlay = cell_overlay(fine.gdf_vorPolys, coarse.gdf_vorPolys)
    ratio, worst = overlay.rule_spread(lens)
    assert worst is not None
    assert np.nanmax(ratio) > 100.0
    # And on a uniform field it must report no spread at all, or the warning cries wolf.
    flat_ratio, _ = overlay.rule_spread(np.full(fine.ncpl, 3.0))
    assert np.allclose(flat_ratio, 1.0, rtol=1e-9)


@pytest.mark.parametrize("rule", ["geometric", "harmonic"])
def test_a_non_positive_value_is_refused_by_the_log_and_reciprocal_rules(meshes, rule):
    """A head may legitimately be negative, so these two must refuse rather than NaN."""

    fine, coarse = meshes
    overlay = cell_overlay(fine.gdf_vorPolys, coarse.gdf_vorPolys)
    values = np.full(fine.ncpl, 5.0)
    values[0] = -1.0
    with pytest.raises(ValueError, match="positive"):
        overlay.reduce(values, rule=rule)


def test_reduce_refuses_a_two_dimensional_array(meshes):
    """Per-layer masking belongs to the caller, so the shape is refused by name."""

    fine, coarse = meshes
    overlay = cell_overlay(fine.gdf_vorPolys, coarse.gdf_vorPolys)
    with pytest.raises(ValueError, match="one value per OLD cell"):
        overlay.reduce(np.ones((2, fine.ncpl)))


# -------------------------------------------------------------------- the frame guard


def test_the_frame_guard_refuses_a_different_crs(meshes):
    """`disjoint` cannot see this: the overlay finds pieces and reports full coverage."""

    fine, coarse = meshes
    relabelled = coarse.gdf_vorPolys.copy().set_crs("EPSG:2926", allow_override=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pytest.raises(ValueError, match="different coordinate systems"):
            refuse_a_different_frame(fine.gdf_vorPolys, relabelled, verb="remap")


def test_the_frame_guard_refuses_a_disjoint_pair(meshes):
    """A mesh in the same frame but another PLACE pins identically under grid_pin."""

    fine, coarse = meshes
    moved = coarse.gdf_vorPolys.copy()
    moved["geometry"] = moved.geometry.translate(1e6, 1e6)
    with pytest.raises(ValueError, match="share NO area"):
        refuse_a_different_frame(fine.gdf_vorPolys, moved, verb="remap")


def test_the_frame_guard_allows_the_same_frame(meshes):
    fine, coarse = meshes
    refuse_a_different_frame(fine.gdf_vorPolys, coarse.gdf_vorPolys, verb="remap")


def test_a_same_crs_unit_rescale_is_NOT_caught_and_that_is_recorded(meshes):
    """The guard's blind spot, pinned so it cannot be mistaken for coverage.

    A mesh rescaled feet -> metres while still declaring EPSG:2927 sits inside the old
    footprint, every cell is covered, and the remap is self-consistent -- it simply
    describes a tenth of the domain. Nothing geometric tells that from a legitimate
    SUBDOMAIN refinement, so a test that fired on it would fire on correct work. Ledger
    211 records the same blind spot for a `.gpkg` source. This test exists so the limit
    is a decision and not a surprise.
    """

    fine, coarse = meshes
    metres = coarse.gdf_vorPolys.copy()
    metres["geometry"] = metres.geometry.scale(1 / 3.28084, 1 / 3.28084, origin=(0, 0))
    refuse_a_different_frame(fine.gdf_vorPolys, metres, verb="remap")


# ---------------------------------------------------------------------- what it refuses


@pytest.mark.parametrize("flag", ["icelltype", "iconvert", "irch", "ievt"])
def test_an_integer_flag_is_refused_by_name(meshes, flag):
    """There is no cell type 0.4, and irch/ievt are zero-based LAYER INDICES."""

    fine, coarse = meshes
    with pytest.raises(RemapRefused, match="integer flag"):
        remap_arrays(
            {flag: np.ones(fine.ncpl)},
            old_polys=fine.gdf_vorPolys, new_polys=coarse.gdf_vorPolys,
            cellwise=(flag,), package="npf",
        )


def test_an_extensive_array_is_refused_by_name(meshes):
    """A rate's denominator is the OLD cell's area; ledger 300 defers it deliberately."""

    fine, coarse = meshes
    with pytest.raises(RemapRefused, match="no measured rule"):
        remap_arrays(
            {"recharge": np.full(fine.ncpl, 1e-4)},
            old_polys=fine.gdf_vorPolys, new_polys=coarse.gdf_vorPolys,
            cellwise=("recharge",), package="rcha",
        )


def test_a_refusal_is_not_swallowed_into_a_silent_no_op(meshes):
    """`remap=True` on a package it cannot serve must RAISE, not move on quietly.

    Catching it and leaving the package stamped `unverifiable` would tell the caller the
    move was fine -- ledger 295's failure again, one verb over.
    """

    fine, coarse = meshes
    sim = mf.simulation(
        mf.gwf("m", context=ModelContext(grid=fine),
               packages=[mf.disv(_stack()),
                         mf.npf(k=10.0, icelltype=np.ones(fine.ncpl, dtype=int))]),
        name="s", packages=[mf.tdis(nper=1)],
    )
    with pytest.raises(RemapRefused, match="integer flag"):
        sim.with_grid(coarse, remap=True)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [({"rule": "nearest"}, "is not a rule"), ({"uncovered": "warn"}, "is not a mode")],
)
def test_an_unknown_keyword_value_is_refused(meshes, kwargs, match):
    fine, coarse = meshes
    with pytest.raises(ValueError, match=match):
        remap_arrays(
            {"k": np.full(fine.ncpl, 1.0)},
            old_polys=fine.gdf_vorPolys, new_polys=coarse.gdf_vorPolys,
            cellwise=("k",), package="npf", **kwargs,
        )


@pytest.mark.parametrize("tier", ["simulation", "model"])
@pytest.mark.parametrize(
    ("kwargs", "match"),
    [({"rule": "geometrc"}, "is not a rule"), ({"uncovered": "neerest"}, "is not a mode")],
)
def test_a_typo_is_refused_at_the_call_with_nothing_to_remap(meshes, tier, kwargs, match):
    """The check lived in `remap_arrays`, reached only for a package holding an ARRAY.

    `ic(strt=40.0)` and `npf(k=1.0)` carry none, so `rule="geometrc"` was accepted.
    """

    fine, coarse = meshes
    sim = _sim(fine, k=1.0)
    target = sim if tier == "simulation" else sim.model("m")
    with pytest.raises(ValueError, match=match):
        target.with_grid(coarse, remap=True, **kwargs)


@pytest.mark.parametrize("tier", ["simulation", "model"])
@pytest.mark.parametrize(
    "kwargs", [{"rule": "harmonic"}, {"uncovered": "raise"}, {"rule": "geometric"}]
)
def test_rule_or_uncovered_without_remap_is_refused_naming_remap(
    meshes, lens, tier, kwargs
):
    """Without `remap=True` both were IGNORED: a `rule=` read as a rule applied."""

    fine, coarse = meshes
    sim = _sim(fine, k=lens)
    target = sim if tier == "simulation" else sim.model("m")
    with pytest.raises(ValueError, match=r"only acts? with remap=True"):
        target.with_grid(coarse, **kwargs)
    # A simulation with no gridded model moves nothing, and still refuses.
    if tier == "simulation":
        bare = mf.simulation(mf.gwf("m", packages=[mf.npf(k=1.0)]), name="s",
                             packages=[mf.tdis(nper=1)])
        with pytest.raises(ValueError, match="remap=True"):
            bare.with_grid(coarse, **kwargs)


def test_a_deferred_grid_model_refuses_remap_by_name(meshes):
    """There are no values on a mesh yet, so there is nothing to carry.

    At the MODEL tier this raises. At the simulation tier the same model is skipped
    entirely instead, because `on=None` moves every model that HAS a grid (ledger 268) --
    so there is no remap to refuse. Both are pinned: the skip is existing behaviour and
    would otherwise read as this refusal failing to fire.
    """

    _, coarse = meshes
    model = mf.gwf("m", packages=[mf.npf(k=10.0)])
    with pytest.raises(ValueError, match="no grid to carry them from"):
        model.with_grid(coarse, remap=True)

    sim = mf.simulation(model, name="s", packages=[mf.tdis(nper=1)])
    assert sim.with_grid(coarse, remap=True).model("m").packages[0].options["k"] == 10.0


# --------------------------------------------------------------- uncovered new cells


def test_an_uncovered_cell_is_filled_from_the_nearest_and_reported_in_detail(meshes, lens):
    """The user's decided default, and the warning has to carry three numbers."""

    fine, _ = meshes
    grown = _grid(8, cell=200.0)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        moved, report = remap_arrays(
            {"k": lens}, old_polys=fine.gdf_vorPolys, new_polys=grown.gdf_vorPolys,
            cellwise=("k",), package="npf",
        )
    assert report["filled"] > 0
    assert report["filled_area"] > 0
    assert report["distance"] > 0
    assert np.all(np.isfinite(moved["k"]))

    said = "\n".join(str(entry.message) for entry in caught)
    assert "npf" in said, "the warning must name the package -- warnings dedupe"
    assert "filled from the nearest one" in said
    assert f"{report['filled']} new cells" in said


def test_uncovered_raise_refuses_instead(meshes, lens):
    fine, _ = meshes
    grown = _grid(8, cell=200.0)
    with pytest.raises(RemapRefused, match="no old cell under them"):
        remap_arrays(
            {"k": lens}, old_polys=fine.gdf_vorPolys, new_polys=grown.gdf_vorPolys,
            cellwise=("k",), package="npf", uncovered="raise",
        )


def test_the_warning_reaches_warnings_not_only_the_logger(meshes, lens, caplog):
    """`import myflopy` installs a NullHandler, so a log line alone reaches NOBODY."""

    fine, coarse = meshes
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        remap_arrays(
            {"k": lens}, old_polys=fine.gdf_vorPolys, new_polys=coarse.gdf_vorPolys,
            cellwise=("k",), package="npf",
        )
    assert any(issubclass(entry.category, UserWarning) for entry in caught)


def test_the_warning_points_at_re_declaration_as_the_better_route(meshes, lens):
    """Measured 7.1% against 12.5%; a reader with a source should be told."""

    fine, coarse = meshes
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        remap_arrays(
            {"k": lens}, old_polys=fine.gdf_vorPolys, new_polys=coarse.gdf_vorPolys,
            cellwise=("k",), package="npf",
        )
    said = "\n".join(str(entry.message) for entry in caught)
    assert "re-declaring" in said and "INTERPOLATED" in said


# ------------------------------------------------------------------ idomain and layers


def test_an_inactive_old_cell_does_not_poison_an_active_new_one(meshes):
    """Mask the WEIGHTS, never pass a NaN value through.

    Measured on a pinching wedge: handing inactive cells through as NaN left ten reached
    new cells NaN, three of them ACTIVE on the new mesh -- and `reached` said nothing.
    """

    fine, coarse = meshes
    values = np.full(fine.ncpl, 10.0)
    values[: fine.ncpl // 2] = np.nan          # as an inactive cell would arrive
    active = np.ones(fine.ncpl, dtype=bool)
    active[: fine.ncpl // 2] = False

    overlay = cell_overlay(fine.gdf_vorPolys, coarse.gdf_vorPolys)
    painted, reached = overlay.reduce(values, rule="mean", active=active)
    assert reached.any()
    assert not np.isnan(painted[reached]).any()


def test_each_layer_is_remapped_separately(meshes):
    fine, coarse = meshes
    source = np.vstack([np.full(fine.ncpl, 2.0), np.full(fine.ncpl, 8.0)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        moved, _ = remap_arrays(
            {"k": source}, old_polys=fine.gdf_vorPolys, new_polys=coarse.gdf_vorPolys,
            cellwise=("k",), package="npf",
        )
    assert moved["k"].shape == (2, coarse.ncpl)
    assert np.allclose(moved["k"][0], 2.0)
    assert np.allclose(moved["k"][1], 8.0)


# -------------------------------------------------------------- stamps and neighbours


def test_a_remapped_package_clears_MOVED_FROM_and_stamps_REMAPPED_FROM(meshes, lens):
    """`moved_package_complaints` fires on MOVED_FROM alone -- leaving it refuses the fix."""

    fine, coarse = meshes
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        moved = _sim(fine, k=lens).with_grid(coarse, remap=True)
    npf = next(p for p in moved.model("m").packages if p.name == "npf")
    assert REMAPPED_FROM in npf.metadata
    assert MOVED_FROM not in npf.metadata
    assert npf.metadata[REMAPPED_FROM]["rule"] == "geometric"


def test_a_pinned_package_is_left_alone(meshes, tmp_path):
    """It is `stale` and has a better route: re-declaration is exact (ledger 296)."""

    import geopandas as gpd
    from shapely.geometry import box

    fine, coarse = meshes
    path = tmp_path / "bc.gpkg"
    xmin, ymin, xmax, ymax = fine.gdf_vorPolys.total_bounds
    gpd.GeoDataFrame(
        {"name": ["west"], "layer": [1], "head": [9.0]},
        geometry=[box(xmin, ymin, xmin + (xmax - xmin) * 0.2, ymax)], crs=fine.crs,
    ).to_file(path, layer="chd", driver="GPKG")

    context = ModelContext(grid=fine)
    sim = mf.simulation(
        mf.gwf("m", context=context,
               packages=[mf.disv(_stack()), mf.npf(k=10.0),
                         mf.chd.gpkg(path, context=context, nper=1, layer="chd")]),
        name="s", packages=[mf.tdis(nper=1)],
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        moved = sim.with_grid(coarse, remap=True)
    chd = next(p for p in moved.model("m").packages if p.name == "chd")
    assert REMAPPED_FROM not in (chd.metadata or {})


def test_the_stamp_is_invisible_to_the_digest(meshes, lens):
    """Metadata is excluded from `package_digest` by contract -- assert it, do not trust it."""

    fine, coarse = meshes
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        moved = _sim(fine, k=lens).with_grid(coarse, remap=True)
    npf = next(p for p in moved.model("m").packages if p.name == "npf")
    from dataclasses import replace

    assert package_digest(npf) == package_digest(replace(npf, metadata={}))


def test_package_move_reports_remapped_before_the_cellwise_branch(meshes, lens):
    """The stamp has to win, or the verdict contradicts the build that follows."""

    fine, _ = meshes
    stamped = mf.npf(k=lens).with_metadata(
        **{REMAPPED_FROM: {"was": "abc", "filled": 3}, CELLWISE_OPTIONS: ["k"]}
    )
    how, why = package_move(stamped, ncpl=fine.ncpl)
    assert how == "remapped"
    assert "3 cells filled" in why


def test_a_remapped_model_CAN_be_saved_by_the_documented_route(meshes, lens, tmp_path):
    """Two reasons it cannot be saved in one step, and the route that gets there anyway.

    Ledger 300 first recorded this as a flat "cannot be saved", which overstated it. The
    two obstacles are real and both predate C6:

    * the remapped values are an ``ndarray``, which ``to_dict`` refuses (ledger 193) --
      and the model was ALREADY unsaveable before the move, which is precisely why it
      needed remapping rather than re-declaring: a package WITH a source reads ``rebuilt``
      and never reaches the engine;
    * the remap needs real cell polygons so it refuses a ``mf.grid_ref`` target, leaving
      a built grid object in ``grid=`` -- itself unsaveable (ledger 257).

    Both have documented answers. This pins the combination, because nothing else here
    would fail if the second ``with_grid`` stopped preserving the remap -- which it does
    only because ``REMAPPED_FROM`` is tested BEFORE the cellwise branch in
    ``package_move``.
    """

    fine, coarse = meshes
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        moved = _sim(fine, k=lens).with_grid(coarse, remap=True)
        # Re-point at the library entry. The verdict and the values both have to survive.
        repointed = moved.with_grid(mf.grid_ref("mesh"))

    assert _verdict(repointed, "npf") == "remapped"
    npf = next(p for p in repointed.model("m").packages if p.name == "npf")
    assert isinstance(npf.options["k"], np.ndarray)

    project = mf.Project(tmp_path / "saved", name="p")
    project.add_grid("mesh", coarse)
    project.add_package("npf/remapped", npf)
    project.add_simulation(
        repointed.with_model(
            repointed.model("m").with_packages(mf.ref("npf/remapped"), slot="npf")
        )
    )
    report = project.save()
    assert report.simulations == 1, report
    assert report.skipped == ()

    reloaded = mf.Project.load(tmp_path / "saved")
    run = reloaded.prepare_run("r", reloaded.simulation("s"))
    run.build()
    written = run.built.simulation.get_model("m").get_package("npf").k.get_data()
    assert np.all(np.isfinite(written))
    assert written.shape == (1, coarse.ncpl)


def test_the_remap_refuses_a_grid_ref_target_by_name(meshes, lens):
    """It needs polygons, and a GridRef is resolved at build. Refusing beats skipping.

    Silently skipping would leave the package stamped `unverifiable` while the caller had
    asked for a remap -- ledger 295's failure again.
    """

    fine, _ = meshes
    with pytest.raises(TypeError, match="not GridRef"):
        _sim(fine, k=lens).with_grid(mf.grid_ref("mesh"), remap=True)

"""A SOLVED field is carried onto another mesh: ``model.hds.remap(grid)``.

Starting heads for a model on a new grid, from the heads an old run ENDED with --
``mf.ic(strt=old.hds.remap(new_vor))`` -- and the same for a concentration or a
temperature. ``with_grid(remap=True)`` carries a model's INPUT ``ic.strt``; nothing
carried a run's OUTPUT until this, and the last time it was done by hand one cell of
~18,000 came back NaN and was patched by hand (ledger 248).

What this file pins, each against the way it went wrong or could:

* **The default is linear, and it is exact for a plane** -- including the band of new
  cells outside the hull of the old centres, which on a jittered mesh is much of the
  boundary (38 of 144 cells), where the area mean it first fell back to was off by 2.79 --
  held to the range the old layer holds.
* **It never interpolates across a gap**: a triangle bridging an inactive strip is not
  used, and the cells either side keep their own side's value.
* **A dry head takes the water table below it; an inactive cell is a gap**, and a gap a
  new cell sits over only is filled from the nearest old cell WITH a warning.
* **The time is the END of the period, checked for that period alone**, so a run that
  stopped later still reads, and a missing end names the fix.
* **Nothing comes back NaN.**
"""

from __future__ import annotations

import inspect
import warnings

import numpy as np
import pytest

import myflopy as mf
from myflopy._numbers import is_no_value
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, rectangular_voronoi
from myflopy.modflow.mf6.grid.remap import (
    DEFAULT_FIELD_METHOD,
    DEFAULT_UNCOVERED,
    RemapRefused,
    no_value_to_nan,
    remap_field,
)
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.mf6.headsplus import DependentVariableFile, HeadsPlus
from myflopy.modflow.mf6.period_end import period_end_record


def _grid(n, cell=100.0):
    return rectangular_voronoi(
        CanonicalModelConfig(nrow=n, ncol=n, nlay=1, nper=1, cell_size=cell)
    )


@pytest.fixture(scope="module")
def coarse():
    """A jittered 6x6 mesh over a 1,200 ft square: its boundary cells are NOT centred on
    the boundary, so the hull of its centres leaves a rim of the square uncovered."""

    return _grid(6, cell=200.0)


@pytest.fixture(scope="module")
def fine():
    return _grid(12)


@pytest.fixture(scope="module")
def grown():
    """A mesh past the coarse one's footprint: 1,600 ft square, so cells lie OUTSIDE it."""

    return _grid(8, cell=200.0)


def _plane(grid, a=50.0, b=0.01, c=-0.02):
    x, y = np.asarray(grid.points, dtype=float).T
    return a + b * x + c * y


def _shifted(grid, dx, *, crs=None):
    """The same mesh moved by ``dx`` (and optionally relabelled), as a new built grid."""

    verts = np.asarray(grid.verts, dtype=float).copy()
    verts[:, 0] += dx
    points = np.asarray(grid.points, dtype=float).copy()
    points[:, 0] += dx
    return VoronoiGridPlus(
        verts=verts, iverts=grid.iverts, xcyc=points,
        crs=grid.crs if crs is None else crs, name="shifted",
    )


# -- the engine: no MF6 ----------------------------------------------------------------


@pytest.mark.parametrize("pair", ["refine", "coarsen"])
def test_a_planar_field_comes_across_exactly_under_linear(coarse, fine, pair):
    """A uniform gradient is reproduced at every new centre, held to the old range.

    The rim is the point: ``coarse`` is a jittered mesh, so many of the fine mesh's
    boundary cells sit outside the hull of the coarse centres (38 of 144). They take the
    nearest old cell's value moved along that cell's own least-squares gradient -- exact
    for a plane -- clamped to the range the old layer holds, so the plane comes back
    exactly, clipped at the extreme corners. The area mean they first took was off by
    2.79 there.
    """

    old, new = (coarse, fine) if pair == "refine" else (fine, coarse)
    field = np.vstack([_plane(old), _plane(old) - 1.0])
    out, report = remap_field(field, old, new)
    for row, shift in ((0, 0.0), (1, -1.0)):
        expected = np.clip(_plane(new) + shift, field[row].min(), field[row].max())
        np.testing.assert_allclose(out[row], expected, atol=1e-9)
    if pair == "refine":
        assert report["rim"] > 0, "the rim this test exists for was not exercised"


def test_the_extension_never_leaves_the_old_range(coarse, fine):
    """A plane fitted across a sharp front would go negative; the layer range stops it.

    Measured without the clamp: -0.055 from a field whose least value was 4.7e-9.
    """

    front = np.exp(-np.asarray(coarse.points, dtype=float)[:, 0] / 60.0)
    out, report = remap_field(front[None, :], coarse, fine)
    assert report["rim"] > 0, "no cell took the extension"
    assert out.min() >= front.min() and out.max() <= front.max()


def test_the_area_mean_is_not_exact_for_a_plane_which_is_why_linear_is_the_default(
    coarse, fine
):
    """``method="area"`` paints every new cell inside an old one with that cell's value."""

    out, _ = remap_field(_plane(coarse)[None, :], coarse, fine, method="area")
    assert np.abs(out[0] - _plane(fine)).max() > 1.0
    assert DEFAULT_FIELD_METHOD == "linear"


@pytest.mark.parametrize("method", ["linear", "area"])
def test_the_same_mesh_returns_the_same_field(coarse, method):
    field = np.vstack([_plane(coarse), _plane(coarse) + 3.0])
    out, report = remap_field(field, coarse, coarse, method=method)
    np.testing.assert_allclose(out, field, atol=1e-9)
    assert report["outside"] == 0


@pytest.mark.parametrize("method", ["linear", "area"])
def test_a_constant_field_stays_constant(coarse, grown, method):
    """Even over cells outside the old footprint, which take the nearest value."""

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        out, _ = remap_field(np.full((2, coarse.ncpl), 7.0), coarse, grown, method=method)
    np.testing.assert_allclose(out, 7.0)


def test_linear_never_interpolates_across_an_inactive_strip(fine):
    """A triangle joining cells either side of an inactive strip is not used.

    West of the strip the field is 10, east 20. Interpolating through a bridging triangle
    would put values between the two on cells that sit wholly on one side.
    """

    x = np.asarray(fine.points, dtype=float)[:, 0]
    field = np.where(x < 600.0, 10.0, 20.0)
    strip = (x > 400.0) & (x < 800.0)
    field[strip] = 1.0e30                                    # INACTIVE
    new = _grid(6, cell=200.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)         # cells over the strip alone
        out, report = remap_field(field[None, :], fine, new)
    new_x = np.asarray(new.points, dtype=float)[:, 0]
    polys = new.gdf_vorPolys.geometry.bounds
    west = (polys.maxx.to_numpy() <= 400.0)
    east = (polys.minx.to_numpy() >= 800.0)
    assert west.any() and east.any()
    np.testing.assert_allclose(out[0][west], 10.0)
    np.testing.assert_allclose(out[0][east], 20.0)
    assert np.isin(np.round(out[0], 9), [10.0, 20.0]).all(), (
        "a value between the two sides means a triangle across the strip was used: "
        f"{sorted(set(np.round(out[0], 3)) - {10.0, 20.0})} at x={new_x[~np.isin(np.round(out[0], 9), [10.0, 20.0])]}"
    )
    # The cells over the strip ALONE are the guessed ones, counted -- independently here.
    valid_area = fine.gdf_vorPolys.geometry[~strip].union_all()
    over_strip_only = new.gdf_vorPolys.geometry.intersection(valid_area).area.to_numpy() < 1e-6
    assert over_strip_only.any(), "no new cell sits over the strip alone; the count is untested"
    assert report["no_value_under"].get(0, 0) == int(over_strip_only.sum())
    assert report["outside"] == 0


# -- no-value cells -----------------------------------------------------------------------


def test_a_dry_head_takes_the_water_table_and_an_inactive_one_is_a_gap():
    raw = np.array([
        [-1e30, -1e30, -1e30, 5.0, 1e30],     # layer 0
        [  4.0,  1e30, -1e30, 4.0, 1e30],     # layer 1
        [  3.0,   2.0,  1e30, 3.0, 1e30],     # layer 2
    ])
    field, filled = no_value_to_nan(raw, water_table=True)
    # cell 0: dry over a wet layer 1 -> 4.0; cell 1: dry over INACTIVE over 2.0 -> 2.0;
    # cell 2: dry over dry over inactive -> no water anywhere below -> a gap.
    np.testing.assert_array_equal(field[0, :4], [4.0, 2.0, np.nan, 5.0])
    assert np.isnan(field[:, 4]).all()                      # inactive stays a gap
    assert np.isnan(field[1, 2])
    assert filled == 2
    concentration, filled = no_value_to_nan(raw, water_table=False)
    assert np.isnan(concentration[0, :3]).all() and filled == 0


def test_cells_over_no_value_only_are_filled_from_the_nearest_with_one_warning(
    coarse, fine, caplog
):
    field = np.vstack([_plane(coarse), _plane(coarse)])
    x = np.asarray(coarse.points, dtype=float)[:, 0]
    field[0, x < 300.0] = 1.0e30
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        out, report = remap_field(field, coarse, fine, label="head", when="the end of period 4")
    remap_warnings = [w for w in caught if "remap" in str(w.message)]
    assert len(remap_warnings) == 1
    text = str(remap_warnings[0].message)
    assert "the end of period 4" in text and "old layer 0" in text
    assert "inactive, or dry with no water below" in text
    assert report["no_value_under"][0] > 0 and report["outside"] == 0
    assert np.isfinite(out).all()
    assert any("remap" in record.getMessage() for record in caplog.records)


def test_cells_outside_the_old_model_warn_with_the_count_area_and_distance(coarse, grown):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        out, report = remap_field(_plane(coarse)[None, :], coarse, grown, layers=0)
    assert out.shape == (grown.ncpl,)
    assert report["outside"] > 0 and report["outside_area"] > 0 and report["distance"] > 0
    (warning,) = [w for w in caught if "remap" in str(w.message)]
    assert f"{report['outside']:,} new cells" in str(warning.message)
    assert "uncovered='raise'" in str(warning.message)
    assert np.isfinite(out).all()


def test_uncovered_raise_refuses_naming_the_counts(coarse, grown):
    with pytest.raises(RemapRefused, match=r"lie outside the old model's cells.*uncovered='nearest'"):
        remap_field(_plane(coarse)[None, :], coarse, grown, uncovered="raise")


def test_a_layer_with_no_value_at_all_is_refused_naming_layers(coarse, fine):
    field = np.vstack([np.full(coarse.ncpl, 1.0e30), _plane(coarse)])
    with pytest.raises(RemapRefused, match=r"old layer 0 holds no head at all.*layers=\[\.\.\.\]"):
        remap_field(field, coarse, fine)
    out, _ = remap_field(field, coarse, fine, layers=[1, 1])          # mapped round it
    np.testing.assert_allclose(out[0], out[1])


# -- the arguments ------------------------------------------------------------------------


def test_layers_is_the_new_to_old_layer_map(coarse, fine):
    field = np.vstack([_plane(coarse), _plane(coarse) + 10.0, _plane(coarse) + 20.0])
    out, report = remap_field(field, coarse, fine, layers=[0, 0, 2, -1])
    assert out.shape == (4, fine.ncpl)
    np.testing.assert_allclose(out[0], out[1])
    np.testing.assert_allclose(out[2], out[3])
    np.testing.assert_allclose(out[2] - out[0], 20.0)
    assert report["rows"] == [0, 0, 2, 2]
    one, _ = remap_field(field, coarse, fine, layers=1)
    assert one.shape == (fine.ncpl,)


@pytest.mark.parametrize(
    "layers, error, message",
    [
        (3, ValueError, "outside the old model's 3 layers"),
        (-4, ValueError, "outside the old model's 3 layers"),
        (True, TypeError, "layer index"),
        ([0, True], TypeError, "layer index"),
        (1.0, TypeError, "layer index"),
        ("01", TypeError, "layer index"),
        ([], ValueError, "empty"),
    ],
)
def test_a_layer_that_is_not_one_is_refused(coarse, fine, layers, error, message):
    field = np.vstack([_plane(coarse)] * 3)
    with pytest.raises(error, match=message):
        remap_field(field, coarse, fine, layers=layers)


@pytest.mark.parametrize(
    "bad, says", [(dict(method="nearest"), "is not a method"), (dict(uncovered="warn"), "is not a mode")]
)
def test_a_choice_that_is_not_one_is_refused_even_with_nothing_to_fill(coarse, bad, says):
    with pytest.raises(ValueError, match=says):
        remap_field(_plane(coarse)[None, :], coarse, coarse, **bad)


def test_a_field_of_the_wrong_shape_is_refused(coarse, fine):
    with pytest.raises(ValueError, match=r"one value per old cell and layer"):
        remap_field(_plane(fine)[None, :], coarse, fine)


def test_a_grid_in_another_crs_or_place_is_refused(coarse):
    with pytest.raises(ValueError, match="different coordinate systems"):
        remap_field(_plane(coarse)[None, :], coarse, _shifted(coarse, 0.0, crs="EPSG:2926"))
    with pytest.raises(ValueError, match="share NO area"):
        remap_field(_plane(coarse)[None, :], coarse, _shifted(coarse, 1.0e6))


class _FloPyLike:
    xcellcenters = None


class _ModelLike:
    vor = None


@pytest.mark.parametrize(
    "grid, names",
    [
        (lambda c: mf.GridSpec.voronoi(boundary=mf.GeoPackageSourceSpec(path="b.gpkg")),
         "spec.resolve"),
        (lambda c: mf.grid_ref("mesh"), "project.grid(key)"),
        (lambda c: _ModelLike(), "model.vor"),
        (lambda c: c.gdf_vorPolys, "no cell CENTRES"),
        (lambda c: _FloPyLike(), "load_mf6_run"),
    ],
)
def test_a_grid_that_is_not_built_is_refused_naming_the_spelling(coarse, grid, names):
    with pytest.raises(TypeError, match="BUILT grid") as error:
        remap_field(_plane(coarse)[None, :], coarse, grid(coarse))
    assert names in str(error.value)


def test_the_noun_defaults_are_the_engines():
    """``remap`` holds literal defaults (it imports the engine lazily); they must agree."""

    parameters = inspect.signature(DependentVariableFile.remap).parameters
    assert parameters["method"].default == DEFAULT_FIELD_METHOD
    assert parameters["uncovered"].default == DEFAULT_UNCOVERED
    assert parameters["per"].default == -1


# -- which saved time: the period's END, that period alone ----------------------------------


_ALL = [(k, p, (k + 1) * 10 / 3) for p in range(3) for k in range(3)]   # nstp=3, 10-day periods


@pytest.mark.parametrize("per, expected", [(0, (2, 0)), (2, (2, 2)), (-1, (2, 2)), (-3, (2, 0))])
def test_a_period_is_read_at_its_end(per, expected):
    assert period_end_record(_ALL, per, perlen=[10.0] * 3, file="the head file") == expected


def test_a_period_short_of_its_end_names_save_head_last():
    first_only = [item for item in _ALL if item[0] == 0]
    with pytest.raises(ValueError, match=r"END is missing.*SAVE HEAD LAST"):
        period_end_record(first_only, 0, perlen=[10.0] * 3, file="the head file")


def test_a_run_that_stopped_still_reads_the_periods_it_finished():
    """The file ends two steps into period 1 of 3 -- MF6 stopped there."""

    stopped = [item for item in _ALL if item[1] == 0] + [(0, 1, 10 / 3), (1, 1, 20 / 3)]
    assert period_end_record(stopped, 0, perlen=[10.0] * 3, file="the head file") == (2, 0)
    with pytest.raises(ValueError, match=r"run may have stopped inside it.*mfsim\.lst"):
        period_end_record(stopped, 1, perlen=[10.0] * 3, file="the head file")
    with pytest.raises(ValueError, match=r"saved no record in stress period 2.*mfsim\.lst"):
        period_end_record(stopped, -1, perlen=[10.0] * 3, file="the head file")


@pytest.mark.parametrize("per, error", [(3, ValueError), (-4, ValueError), (True, TypeError), (1.0, TypeError)])
def test_a_period_that_is_not_one_is_refused(per, error):
    with pytest.raises(error):
        period_end_record(_ALL, per, perlen=[10.0] * 3)


def test_without_a_tdis_the_last_saved_step_is_taken_unchecked():
    first_only = [item for item in _ALL if item[0] == 0]
    assert period_end_record(first_only, -1) == (0, 2)


# -- the reader, on what MF6 actually wrote -----------------------------------------------


def test_the_same_mesh_returns_the_run_s_own_heads(sentinel_run, triangulated_grid_finer):
    """``sentinel_run`` is built on ``triangulated_grid_finer``'s mesh: identity.

    Its dry corner has nothing below it (one layer) and its inactive block nothing at
    all, so both are filled from the nearest valid cell -- with a warning -- and the
    result holds no NaN.
    """

    expected = sentinel_run.hds.get(per=0)["head"].to_numpy()
    valid = np.isfinite(expected)
    assert (~valid).any(), "the fixture lost its dry and inactive cells"
    for method in ("linear", "area"):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            out = sentinel_run.hds.remap(triangulated_grid_finer, method=method)
        assert out.shape == (1, triangulated_grid_finer.ncpl)
        np.testing.assert_allclose(out[0][valid], expected[valid], atol=1e-6)
        assert np.isfinite(out).all()
        assert any("sat only over cells holding no head" in str(w.message) for w in caught)
        assert any("the end of stress period 0" in str(w.message) for w in caught)


def test_a_new_model_starts_from_the_remapped_heads_and_runs(
    sentinel_run, triangulated_grid, tmp_path
):
    """The user's line: ``mf.ic(strt=old.hds.remap(new_vor))`` builds and MF6 runs."""

    new = triangulated_grid
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        strt = sentinel_run.hds.remap(new)
    props = new.get_disv_gridprops()
    edges = new.gdf_vorPolys.geometry.bounds
    chd = [[(0, int(c)), 9_990.0] for c in np.flatnonzero(edges.minx <= 1e-6)]
    model = mf.gwf("m", grid=new, packages=[
        mf.disv(nlay=1, ncpl=props["ncpl"], nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=10_100.0,
                botm=9_900.0),
        mf.npf(k=10.0, icelltype=1), mf.ic(strt=strt),
        mf.chd(stress_period_data={0: chd}),
        mf.oc(head_filerecord="m.hds", saverecord=[("HEAD", "ALL")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]), mf.ims(models=["m"])])
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "Run 'r' built but its frozen snapshot", UserWarning)
        run = mf.Project(tmp_path / "p", name="p").run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])
    written = run.model().fl.ic.strt.array.ravel()
    np.testing.assert_allclose(written, strt.ravel(), rtol=1e-7)


def test_a_reader_with_no_model_remaps_on_the_file_s_own_clock(
    sentinel_run, triangulated_grid_finer
):
    """The recipe's precondition: a head file and a grid, no model object."""

    reader = HeadsPlus(hds_path=sentinel_run.hds.output_path, vor=triangulated_grid_finer)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        out = reader.remap(triangulated_grid_finer)
        expected = sentinel_run.hds.remap(triangulated_grid_finer)
    np.testing.assert_allclose(out, expected)


def test_per_kstpkper_and_the_refusals_reach_the_reader(sentinel_run, triangulated_grid):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        by_period = sentinel_run.hds.remap(triangulated_grid, per=0)
        by_step = sentinel_run.hds.remap(triangulated_grid, kstpkper=(0, 0))
    np.testing.assert_allclose(by_period, by_step)
    with pytest.raises(ValueError, match="is not a stress period"):
        sentinel_run.hds.remap(triangulated_grid, per=1)
    with pytest.raises(ValueError, match=r"kstpkper \(0, 3\) is unavailable"):
        sentinel_run.hds.remap(triangulated_grid, kstpkper=(0, 3))
    with pytest.raises(TypeError, match="BUILT grid"):
        sentinel_run.hds.remap(sentinel_run)


def test_a_disu_model_is_refused_by_name(sentinel_run, triangulated_grid, monkeypatch):
    monkeypatch.setattr(sentinel_run, "_grid_type_override", "disu", raising=False)
    with pytest.raises(ValueError, match="DISU"):
        sentinel_run.hds.remap(triangulated_grid)


@pytest.mark.slow
@pytest.mark.parametrize("kind, noun", [("gwt", "conc"), ("gwe", "temp")])
def test_a_concentration_and_a_temperature_remap_too(kind, noun, tmp_path, triangulated_grid):
    """A transport model's starting field, through the same verb.

    The coupled toy lives on an 800 ft square; ``triangulated_grid`` covers the 600 ft
    square inside it, so every new cell has old cells under it.
    """

    from tests.test_gwt_gwe_results import _coupled_run

    model = _coupled_run(tmp_path, kind).model("trans")
    reader = getattr(model, noun)
    out = reader.remap(triangulated_grid)
    assert out.shape == (1, triangulated_grid.ncpl)
    assert np.isfinite(out).all()
    expected = reader.get(per=0)[reader.value_name]
    assert expected.min() - 1e-9 <= out.min() and out.max() <= expected.max() + 1e-9


# -- the review's gaps (2026-10-15): each test names the mutant it kills --------------------


@pytest.mark.parametrize("missing", [(0, 1), (1, 2), (0, 2)])
def test_one_non_local_edge_is_enough_to_refuse_a_triangle(missing):
    """Any ONE edge joining cells that share no vertex makes the triangle unusable.

    At a notch's corner a bridging triangle has one such edge, not two as across a strip;
    checking only two of the three edges would let it through.
    """

    import scipy.sparse as sparse

    from myflopy.modflow.mf6.grid.remap import _Triangulation

    points = np.array([[0.0, 0.0], [10.0, 0.0], [0.0, 10.0]])
    touch = np.ones((3, 3), dtype=int)
    touch[missing[0], missing[1]] = touch[missing[1], missing[0]] = 0
    triangle = _Triangulation(
        points, np.arange(3), np.array([[2.0, 2.0]]), sparse.csr_matrix(touch)
    )
    assert not triangle.usable[0]
    whole = _Triangulation(
        points, np.arange(3), np.array([[2.0, 2.0]]), sparse.csr_matrix(np.ones((3, 3)))
    )
    assert whole.usable[0]


def test_rim_and_bridged_count_cell_layers_the_hull_and_the_local_test_name(coarse, fine):
    """``rim`` is every reached new centre outside the old hull -- counted independently."""

    from scipy.spatial import Delaunay

    one, report_one = remap_field(_plane(coarse)[None, :], coarse, fine)
    outside_hull = Delaunay(np.asarray(coarse.points)).find_simplex(np.asarray(fine.points)) < 0
    assert report_one["rim"] == int(outside_hull.sum()) == 38
    _, report_two = remap_field(np.vstack([_plane(coarse)] * 2), coarse, fine)
    assert report_two["rim"] == 2 * report_one["rim"]          # summed per old layer read
    assert report_two["bridged"] == 2 * report_one["bridged"]
    _, report_rows = remap_field(np.vstack([_plane(coarse)] * 2), coarse, fine, layers=[0, 0, 0])
    assert report_rows["rim"] == report_one["rim"]              # one old layer, read thrice


def test_outside_cells_take_the_nearest_valid_centre_and_the_warning_is_exact(coarse, grown):
    """The values, the area and the distance a fill reports -- computed independently."""

    from scipy.spatial import cKDTree

    data = _plane(coarse)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        out, report = remap_field(data[None, :], coarse, grown, layers=0)
    footprint = coarse.gdf_vorPolys.geometry.union_all()
    outside = grown.gdf_vorPolys.geometry.intersection(footprint).area.to_numpy() < 1e-6
    distance, nearest = cKDTree(np.asarray(coarse.points)).query(np.asarray(grown.points)[outside])
    np.testing.assert_allclose(out[outside], data[nearest])
    assert report["outside"] == int(outside.sum())
    np.testing.assert_allclose(
        report["outside_area"], grown.gdf_vorPolys.geometry.area.to_numpy()[outside].sum()
    )
    np.testing.assert_allclose(report["distance"], distance.max())
    assert report["no_value_under"] == {0: 0}
    (warning,) = [w for w in caught if "remap" in str(w.message)]
    assert f"{report['outside']:,} new cells ({report['outside_area']:,.0f}" in str(warning.message)
    assert f"up to {report['distance']:,.0f} away" in str(warning.message)


def test_uncovered_raise_names_each_layer_s_count(coarse, grown):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        _, report = remap_field(np.vstack([_plane(coarse)] * 2), coarse, grown)
    n = report["outside"]
    with pytest.raises(
        RemapRefused,
        match=rf"old layer 0: {n}; old layer 1: {n} -- of which {n} lie outside",
    ):
        remap_field(np.vstack([_plane(coarse)] * 2), coarse, grown, uncovered="raise")


def test_a_gap_is_filled_from_its_own_layer_not_another(coarse, fine):
    """The agreed rule: the nearest valid cell IN ITS LAYER (ledger 595's NOT-built list)."""

    from scipy.spatial import cKDTree

    x = np.asarray(coarse.points, dtype=float)[:, 0]
    top, bottom = _plane(coarse), _plane(coarse) + 100.0     # layer 1 far from layer 0
    top[x < 300.0] = 1.0e30
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        out, _ = remap_field(np.vstack([top, bottom]), coarse, fine)
    valid_area = coarse.gdf_vorPolys.geometry[x >= 300.0].union_all()
    gap = fine.gdf_vorPolys.geometry.intersection(valid_area).area.to_numpy() < 1e-6
    assert gap.any()
    ok = np.flatnonzero(x >= 300.0)
    _, nearest = cKDTree(np.asarray(coarse.points)[ok]).query(np.asarray(fine.points)[gap])
    np.testing.assert_allclose(out[0][gap], top[ok][nearest])


def test_a_row_does_not_depend_on_which_other_rows_are_read(coarse, fine):
    """Triangulations are cached by validity mask; a row reads only its own layer."""

    x = np.asarray(coarse.points, dtype=float)[:, 0]
    top, bottom = _plane(coarse), _plane(coarse, a=80.0, b=-0.03, c=0.01)
    top[x > 900.0] = -1.0e30                                  # dry: takes layer 1 below
    bottom[x < 200.0] = 1.0e30
    field = np.vstack([top, bottom])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        alone, _ = remap_field(field, coarse, fine, layers=[1])
        both, report = remap_field(field, coarse, fine, layers=[1, 0])
        reversed_, _ = remap_field(field, coarse, fine, layers=[0, 1])
    np.testing.assert_array_equal(alone[0], both[0])
    np.testing.assert_array_equal(both[0], reversed_[1])
    np.testing.assert_array_equal(both[1], reversed_[0])
    assert report["dry_filled"] == int((x > 900.0).sum() - ((x > 900.0) & (x < 200.0)).sum())


def test_a_remap_that_guessed_nothing_says_nothing(coarse):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        remap_field(np.vstack([_plane(coarse)] * 2), coarse, coarse)
    assert not [w for w in caught if "remap" in str(w.message)]


def test_active_keeps_a_cell_the_new_model_never_reads_out_of_the_counts(coarse, grown):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        _, report = remap_field(_plane(coarse)[None, :], coarse, grown)
    footprint = coarse.gdf_vorPolys.geometry.union_all()
    outside = grown.gdf_vorPolys.geometry.intersection(footprint).area.to_numpy() < 1e-6
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        out, quiet = remap_field(
            _plane(coarse)[None, :], coarse, grown, uncovered="raise", active=~outside[None, :]
        )
    assert quiet["outside"] == 0 and not [w for w in caught if "remap" in str(w.message)]
    assert np.isfinite(out).all()                    # filled all the same
    with pytest.raises(ValueError, match="active must say which new cells"):
        remap_field(_plane(coarse)[None, :], coarse, grown, active=np.ones((2, grown.ncpl)))


def test_a_disu_s_stacked_cells_are_refused_on_either_side(coarse, fine):
    """Two layers' cells on one polygon: centres repeat, and the cells overlap."""

    stacked = VoronoiGridPlus(
        verts=np.asarray(coarse.verts, dtype=float),
        iverts=list(coarse.iverts) + list(coarse.iverts),
        xcyc=np.vstack([coarse.points, coarse.points]),
        crs=coarse.crs, name="stacked",
    )
    with pytest.raises(RemapRefused, match="old grid share a centre"):
        remap_field(np.vstack([_plane(stacked)]), stacked, fine)
    with pytest.raises(RemapRefused, match="new grid share a centre"):
        remap_field(_plane(fine)[None, :], fine, stacked)


def test_with_grid_reads_a_starting_head_s_no_value_numbers_by_the_same_rule(coarse, fine):
    """``remap_arrays`` averaged a -1e30 into every new cell over it (the review measured
    magnitudes to 1e30); a dry starting head now takes the water table below it."""

    from myflopy.modflow.mf6.grid.remap import remap_arrays

    strt = np.vstack([_plane(fine), _plane(fine) - 5.0])
    x = np.asarray(fine.points, dtype=float)[:, 0]
    strt[0, x > 1000.0] = -1.0e30
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        moved, _ = remap_arrays(
            {"strt": strt}, old_polys=fine.gdf_vorPolys, new_polys=coarse.gdf_vorPolys,
            cellwise=("strt",), package="ic", old_grid=fine, new_grid=coarse,
        )
    assert np.all(np.abs(moved["strt"]) < 1.0e29)
    assert np.isfinite(moved["strt"]).all()
    lowest = (strt[1].min(), strt[0][strt[0] > -1e29].max())
    assert moved["strt"].min() >= lowest[0] - 1e-9 and moved["strt"].max() <= lowest[1] + 1e-9


def test_a_period_with_no_record_names_the_last_one_zero_based():
    stopped = [item for item in _ALL if item[1] == 0] + [(0, 1, 10 / 3)]
    with pytest.raises(ValueError, match=r"kstpkper=\(0, 1\)"):
        period_end_record(stopped, 2, perlen=[10.0] * 3, file="the head file")


def test_the_stopped_run_wording_appears_only_for_the_period_the_file_stops_in():
    first_only = [(0, p, 10 / 3) for p in range(3)]
    with pytest.raises(ValueError) as error:
        period_end_record(first_only, 0, perlen=[10.0] * 3, file="the head file")
    assert "stopped" not in str(error.value)
    with pytest.raises(ValueError, match="stopped inside it"):
        period_end_record(first_only, 2, perlen=[10.0] * 3, file="the head file")


def test_nonconverged_steps_reads_mf6_s_own_lines(tmp_path):
    from myflopy.run_record import nonconverged_steps

    (tmp_path / "mfsim.lst").write_text(
        " Solution 1 did not converge for stress period 2 and time step 3\n"
        " Solution Group 1 did not converge for stress period 2 and time step 3\n"
        " Solution Group 1 did not converge for stress period 5 and time step 1\n"
    )
    assert nonconverged_steps(tmp_path) == {(2, 1), (0, 4)}
    assert nonconverged_steps(tmp_path / "absent") == set()


# -- the reader, on a run with a clock ------------------------------------------------------


@pytest.fixture(scope="module")
def clock_run(tmp_path_factory):
    """Two layers, three steady periods of three steps, the east boundary rising per period.

    Output control saves every step of periods 0 and 2 and only the FIRST of period 1, so
    period 1's end is missing. Layer 0's north-east corner sits above every head: DRY.
    """

    import geopandas as gpd
    from shapely.geometry import box

    root = tmp_path_factory.mktemp("clock")
    boundary = root / "boundary.gpkg"
    gpd.GeoDataFrame({"n": ["d"]}, geometry=[box(0, 0, 600, 600)], crs=2927).to_file(boundary)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        vor = mf.GridSpec.voronoi(
            boundary=mf.GeoPackageSourceSpec(path=boundary), max_area=12_000.0
        ).resolve(workspace=root / "grid")
    x = np.asarray(vor.centroids_x, dtype=float)
    y = np.asarray(vor.centroids_y, dtype=float)
    raised = (x > 400) & (y > 400)
    edges = vor.gdf_vorPolys.geometry.bounds
    west = np.flatnonzero(edges.minx <= 1e-6)
    east = np.flatnonzero(edges.maxx >= 600 - 1e-6)
    chd = {
        per: [[(1, int(c)), 9_990.0] for c in west]
        + [[(1, int(c)), 10_010.0 + 10.0 * per] for c in east]
        for per in range(3)
    }
    props = vor.get_disv_gridprops()
    model = mf.gwf("m", grid=vor, packages=[
        mf.disv(nlay=2, ncpl=props["ncpl"], nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=10_100.0,
                botm=[np.where(raised, 10_050.0, 9_950.0), np.full(vor.ncpl, 9_800.0)]),
        mf.npf(k=10.0, icelltype=1), mf.ic(strt=10_000.0), mf.chd(stress_period_data=chd),
        mf.oc(head_filerecord="m.hds", saverecord={
            0: [("HEAD", "ALL")], 1: [("HEAD", "FIRST")], 2: [("HEAD", "ALL")]}),
    ])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=3, perioddata=[(1.0, 3, 1.0)] * 3), mf.ims(models=["m"])])
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "Run 'r' built but its frozen snapshot", UserWarning)
        run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])
    return run.model()


def _raw(model, key):
    return model.hds.get_data(kstpkper=key).reshape(model.hds.nlay, -1)


@pytest.mark.slow
def test_per_reads_the_record_that_ends_the_period(clock_run, triangulated_grid):
    hds = clock_run.hds
    assert is_no_value(_raw(clock_run, (2, 0))).any(), "the fixture lost its dry corner"
    for per, key in ((0, (2, 0)), (2, (2, 2)), (-1, (2, 2))):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            out = hds.remap(triangulated_grid, per=per, method="area")
            expected, _ = remap_field(_raw(clock_run, key), hds.vor, triangulated_grid,
                                      method="area")
        np.testing.assert_array_equal(out, expected)
    with pytest.raises(ValueError, match=r"END is missing.*SAVE HEAD LAST"):
        hds.remap(triangulated_grid, per=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        exact = hds.remap(triangulated_grid, kstpkper=(0, 1))
        expected, _ = remap_field(_raw(clock_run, (0, 1)), hds.vor, triangulated_grid)
    np.testing.assert_array_equal(exact, expected)


@pytest.mark.slow
@pytest.mark.parametrize(
    "kstpkper, error",
    [((0, 1.7), TypeError), ("01", TypeError), ((True, True), TypeError), (3, TypeError),
     ((0, 1, 2), TypeError), ((1, 1), ValueError)],
)
def test_kstpkper_must_be_two_saved_integers(clock_run, triangulated_grid, kstpkper, error):
    with pytest.raises(error, match="kstpkper"):
        clock_run.hds.remap(triangulated_grid, kstpkper=kstpkper)


@pytest.mark.slow
def test_the_noun_forwards_every_choice_to_the_engine(clock_run, triangulated_grid):
    """layers, method, uncovered and the dry rule reach the engine as given."""

    hds = clock_run.hds
    raw = _raw(clock_run, (2, 2))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        out = hds.remap(triangulated_grid, layers=[1, 0, 0], method="area")
        expected, report = remap_field(raw, hds.vor, triangulated_grid, layers=[1, 0, 0],
                                       method="area", water_table=True)
    np.testing.assert_array_equal(out, expected)
    assert report["dry_filled"] > 0                         # the dry corner took layer 1
    one = hds.remap(triangulated_grid, layers=1, uncovered="raise")
    assert one.shape == (triangulated_grid.ncpl,)


@pytest.mark.slow
def test_a_step_mf6_did_not_converge_on_is_refused(clock_run, triangulated_grid, tmp_path):
    """MF6 saves a failed step's heads; the listing names it, and the read refuses it."""

    import shutil

    from myflopy.modflow.mf6.headsplus import read_saved_field

    hds = clock_run.hds
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    shutil.copy(hds.output_path, run_dir / "m.hds")
    (run_dir / "mfsim.lst").write_text(
        "Solution Group 1 did not converge for stress period 3 and time step 3\n"
    )
    reader = HeadsPlus(hds_path=run_dir / "m.hds", vor=hds.vor)
    with pytest.raises(ValueError, match=r"did NOT converge on"):
        read_saved_field(reader, per=-1, perlen=[1.0] * 3, ncpl=hds.vor.ncpl, run_dir=run_dir)
    with pytest.raises(ValueError, match=r"did NOT converge on"):
        reader.remap(triangulated_grid)                     # finds mfsim.lst beside the file
    raw, key, _ = read_saved_field(reader, per=0, perlen=[1.0] * 3, ncpl=hds.vor.ncpl,
                                   run_dir=run_dir)
    assert key == (2, 0)


@pytest.mark.slow
def test_a_reader_with_no_model_says_its_record_is_unchecked(clock_run, triangulated_grid, grown):
    reader = HeadsPlus(hds_path=clock_run.hds.output_path, vor=clock_run.hds.vor)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        reader.remap(_grid(8, cell=100.0))                # 800 ft: cells outside the model
    (warning,) = [w for w in caught if "remap" in str(w.message)]
    assert "the last saved step of stress period 2 (no TDIS" in str(warning.message)


def _lattice(n, cell):
    """A true n x n lattice of square cells as a built grid: co-circular quads everywhere."""

    edges = np.arange(n + 1) * cell
    gx, gy = np.meshgrid(edges, edges)
    verts = np.column_stack((gx.ravel(), gy.ravel()))
    rings, centres = [], []
    for row in range(n):
        for col in range(n):
            a = row * (n + 1) + col
            rings.append([a, a + 1, a + n + 2, a + n + 1])
            centres.append(((col + 0.5) * cell, (row + 0.5) * cell))
    return VoronoiGridPlus(verts=verts, iverts=rings, xcyc=np.array(centres), crs="EPSG:2927",
                           name="lattice")


def test_on_a_lattice_no_triangle_is_refused_and_a_plane_is_exact():
    """Every Delaunay diagonal of a lattice joins cells sharing ONE vertex: local.

    A face-sharing test (two vertices) would refuse the co-circular quads' diagonals and
    hand every interior cell to the extension.
    """

    old, new = _lattice(8, 150.0), _grid(12)
    out, report = remap_field(_plane(old)[None, :], old, new)
    assert report["bridged"] == 0
    np.testing.assert_allclose(
        out[0], np.clip(_plane(new), _plane(old).min(), _plane(old).max()), atol=1e-9
    )


def test_an_old_grid_whose_cells_overlap_is_refused(coarse, fine):
    """Cells of two meshes in one grid: distinct centres, overlapping polygons."""

    offset = len(coarse.verts)
    both = VoronoiGridPlus(
        verts=np.vstack([coarse.verts, fine.verts]),
        iverts=list(coarse.iverts) + [[v + offset for v in ring] for ring in fine.iverts],
        xcyc=np.vstack([coarse.points, fine.points]),
        crs=coarse.crs, name="overlapping",
    )
    with pytest.raises(RemapRefused, match="cells OVERLAP one another"):
        remap_field(_plane(both)[None, :], both, fine)


def _columns(xedges, yedges, *, centre_x=None):
    """A built grid of rectangles on ``xedges`` x ``yedges``; ``centre_x`` moves column centres."""

    xedges, yedges = np.asarray(xedges, float), np.asarray(yedges, float)
    nx, ny = len(xedges) - 1, len(yedges) - 1
    gx, gy = np.meshgrid(xedges, yedges)
    verts = np.column_stack((gx.ravel(), gy.ravel()))
    middles = (xedges[:-1] + xedges[1:]) / 2 if centre_x is None else np.asarray(centre_x, float)
    rings, centres = [], []
    for row in range(ny):
        for col in range(nx):
            a = row * (nx + 1) + col
            rings.append([a, a + 1, a + nx + 2, a + nx + 1])
            centres.append((middles[col], (yedges[row] + yedges[row + 1]) / 2))
    return VoronoiGridPlus(verts=verts, iverts=rings, xcyc=np.array(centres), crs="EPSG:2927",
                           name="columns")


def test_the_extension_anchors_on_a_cell_under_the_new_one_not_across_the_gap():
    """The review's case: the nearest valid CENTRE is across the strip, the cells under it are not.

    Old columns of 100 ft with x 500-700 inactive; the column east of the strip has its MF6
    centre at x=710, inside its cell. A new column 480-690 sits over the west side's last
    20 ft and the strip: its centre (585) is nearer x=710 than the west centre at 450, so
    anchoring on the nearest centre read the EAST side's 20.0 where only west cells (10.0)
    lie under it.
    """

    xedges = np.arange(0.0, 1300.0, 100.0)
    centre_x = (xedges[:-1] + xedges[1:]) / 2
    centre_x[7] = 710.0                                        # [700, 800]'s centre
    old = _columns(xedges, np.arange(0.0, 1300.0, 100.0), centre_x=centre_x)
    x = np.asarray(old.points)[:, 0]
    field = np.where(x < 600.0, 10.0, 20.0)
    field[(x > 500.0) & (x < 700.0)] = 1.0e30
    new = _columns([0.0, 480.0, 690.0, 1200.0], np.arange(0.0, 1300.0, 100.0))
    out, report = remap_field(field[None, :], old, new)
    middle = np.isclose(np.asarray(new.points)[:, 0], 585.0)
    assert report["bridged"] + report["rim"] > 0, "the column did not take the extension"
    np.testing.assert_allclose(out[0][middle], 10.0)


@pytest.mark.parametrize("low, high, worst", [(300.0, 380.0, 4.0), (820.0, 900.0, 3.0)])
def test_a_gradient_across_a_chain_of_cells_is_not_fitted(fine, low, high, worst):
    """One band of valid cells: each cell's neighbours lie near ONE line.

    The first version fitted the slope across the band from their jitter -- 9.7 and 22.8
    worst against 2.9 and 2.0 now (``field_remap_rim.py``, part 4).
    """

    xy = np.asarray(fine.points)
    band = (xy[:, 1] > low) & (xy[:, 1] < high)
    field = np.full(fine.ncpl, 1.0e30)
    field[band] = 100.0 - 30.0 * (xy[band, 0] / 1200.0) ** 2
    new = _grid(24, cell=50.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        out, _ = remap_field(field[None, :], fine, new)
    under = fine.gdf_vorPolys.geometry[band].union_all()
    reached = new.gdf_vorPolys.geometry.intersection(under).area.to_numpy() > 1e-6
    truth = 100.0 - 30.0 * (np.asarray(new.points)[reached, 0] / 1200.0) ** 2
    assert np.abs(out[0][reached] - truth).max() < worst


def test_with_grid_carries_a_starting_head_linear_between_cell_centres(coarse, fine):
    """One rule for moving heads (ledger 600): ``with_grid(remap=True)``'s ``ic.strt`` is
    the field engine's linear carry -- a plane comes across exactly inside the old range,
    where the area mean it used would have stair-stepped it."""

    from myflopy.modflow.mf6.grid.remap import remap_arrays

    strt = np.vstack([_plane(coarse), _plane(coarse) - 2.0])
    moved, report = remap_arrays(
        {"strt": strt}, old_polys=coarse.gdf_vorPolys, new_polys=fine.gdf_vorPolys,
        cellwise=("strt",), package="ic", old_grid=coarse, new_grid=fine,
    )
    expected, _ = remap_field(strt, coarse, fine)
    np.testing.assert_array_equal(moved["strt"], expected)
    np.testing.assert_allclose(
        moved["strt"][0], np.clip(_plane(fine), strt[0].min(), strt[0].max()), atol=1e-9
    )
    area, _ = remap_field(strt, coarse, fine, method="area")
    assert np.abs(area[0] - _plane(fine)).max() > 1.0      # what the area mean gave
    assert report["filled"] == 0 and report.get("strt_filled", 0) == 0


def test_a_starting_head_on_a_removed_old_cell_is_a_gap(coarse, fine):
    """A cell the OLD domain removed holds an arbitrary starting head; it is never carried,
    and a dry cell above it does not take it for the water table."""

    from myflopy.modflow.mf6.grid.remap import remap_arrays

    x = np.asarray(coarse.points, dtype=float)[:, 0]
    strt = np.vstack([_plane(coarse), _plane(coarse) - 2.0])
    removed = x < 300.0
    strt[1, removed] = 1.0e6                               # arbitrary, as MF6 ignores it
    strt[0, removed] = -1.0e30                              # dry above the removed cell
    domain = np.ones((2, coarse.ncpl), dtype=int)
    domain[1, removed] = 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        moved, _ = remap_arrays(
            {"strt": strt}, old_polys=coarse.gdf_vorPolys, new_polys=fine.gdf_vorPolys,
            cellwise=("strt",), package="ic", domain=domain, old_grid=coarse, new_grid=fine,
        )
    assert moved["strt"].max() < 1.0e5

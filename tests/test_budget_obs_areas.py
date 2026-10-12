"""``Budget.plot_budget_obs``: the areas, the records and the clock (ledger 476).

Each test is one of the defects ledger 441 measured: an area file in another CRS read
0.0 for every area, a second record on a cell was dropped, the rows were the HEAD
file's saved times, a ``q`` the budget lacks returned an empty table, the y axis was
fixed at ``[0, 4]``, and a multi-part ``.gpkg`` area or a CRS-less file raised.
"""

from __future__ import annotations

import warnings

import geopandas as gpd
import pandas as pd
import pytest
import shapely

from myflopy.modflow.mf6.budget import Budget


@pytest.fixture(scope="module")
def drn_areas(canonical_run, tmp_path_factory):
    """Three drain cells' polygons as one area, written in several shapes."""

    model = canonical_run
    budget = Budget(model, "drn")
    cells = sorted(set(budget.df.index.get_level_values("node")))[:3]
    polys = model.vor.gdf_vorPolys
    root = tmp_path_factory.mktemp("budget_obs_areas")
    hull = polys.loc[cells].union_all().convex_hull.buffer(-1)
    multi = polys.loc[cells].union_all().buffer(-1)
    paths = {}
    frame = gpd.GeoDataFrame({"name": ["A"]}, geometry=[hull], crs=model.vor.crs)
    paths["grid"] = root / "areas.gpkg"
    frame.to_file(paths["grid"])
    paths["wgs"] = root / "areas_4326.gpkg"
    frame.to_crs(4326).to_file(paths["wgs"])
    paths["multi"] = root / "areas_multi.gpkg"
    gpd.GeoDataFrame({"name": ["M"]}, geometry=[multi], crs=model.vor.crs).to_file(paths["multi"])
    paths["naive"] = root / "areas_no_crs.gpkg"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gpd.GeoDataFrame({"name": ["A"]}, geometry=[hull]).to_file(paths["naive"])
    paths["off"] = root / "areas_off_grid.gpkg"
    xmin, ymin, xmax, ymax = polys.total_bounds
    far = shapely.box(xmax + 1e4, ymax + 1e4, xmax + 2e4, ymax + 2e4)
    gpd.GeoDataFrame({"name": ["far"]}, geometry=[far], crs=model.vor.crs).to_file(paths["off"])
    assert multi.geom_type == "MultiPolygon"
    return model, cells, paths


def _table(model, path, **options):
    return Budget(model, "drn").plot_budget_obs(shp_gpkg=path, plot_fig=False, **options)


def test_an_area_file_in_another_crs_is_reprojected(drn_areas):
    model, _cells, paths = drn_areas
    in_grid = _table(model, paths["grid"])
    in_wgs = _table(model, paths["wgs"])
    assert (in_grid["A"] > 0).all()
    pd.testing.assert_frame_equal(in_grid, in_wgs, rtol=1e-9)


def test_a_file_with_no_crs_is_taken_in_the_grids_and_a_str_path_reads(drn_areas):
    model, _cells, paths = drn_areas
    pd.testing.assert_frame_equal(_table(model, paths["grid"]), _table(model, str(paths["naive"])))


def test_a_multi_part_area_is_one_area(drn_areas):
    model, cells, paths = drn_areas
    table = _table(model, paths["multi"])
    flows = Budget(model, "drn").df["q"]
    expected = -flows[flows.index.get_level_values("node").isin(cells)].groupby(
        level="kstpkper", sort=False).sum() / 86400
    pd.testing.assert_series_equal(table["M"], expected, check_names=False, check_index=False)


def test_every_record_on_a_cell_is_added(drn_areas):
    model, _cells, paths = drn_areas
    single = _table(model, paths["grid"])
    budget = Budget(model, "drn")
    budget._df = pd.concat([budget.df, budget.df])  # a second, equal record on every cell
    doubled = budget.plot_budget_obs(shp_gpkg=paths["grid"], plot_fig=False)
    pd.testing.assert_frame_equal(doubled, single * 2)


def test_the_rows_are_the_budget_files_own_saved_times(drn_areas):
    model, _cells, paths = drn_areas
    full = _table(model, paths["grid"])
    budget = Budget(model, "drn")
    times = budget.df.index.get_level_values("kstpkper")
    assert list(full.index) == list(dict.fromkeys(times))
    dropped = full.index[2]
    budget._df = budget.df[[time != dropped for time in times]]
    thinned = budget.plot_budget_obs(shp_gpkg=paths["grid"], plot_fig=False)
    assert dropped not in list(thinned.index)
    pd.testing.assert_frame_equal(thinned, full.drop(index=[dropped]))


def test_a_q_the_budget_lacks_is_refused(drn_areas):
    model, _cells, paths = drn_areas
    with pytest.raises(KeyError, match="'nope' is not a column of this budget"):
        _table(model, paths["grid"], q="nope")


def test_an_area_on_no_cell_is_refused(drn_areas):
    model, _cells, paths = drn_areas
    with pytest.raises(ValueError, match="1 of 1 observation areas cover no grid cell"):
        _table(model, paths["off"])


def test_the_y_axis_fits_the_data_by_default(drn_areas):
    model, _cells, paths = drn_areas
    figure = _table(model, paths["grid"], return_fig=True)
    assert figure.layout.yaxis.range is None
    fixed = _table(model, paths["grid"], return_fig=True, y_range=[0, 1])
    assert tuple(fixed.layout.yaxis.range) == (0, 1)


def test_times_must_label_every_row(drn_areas):
    model, _cells, paths = drn_areas
    with pytest.raises(ValueError, match="times= gives 2 labels for 6 saved times"):
        _table(model, paths["grid"], return_fig=True, times=["2020-01-01", "2020-01-02"])


def test_sfr_plot_flows_draws_each_reachs_outflow_on_its_own_clock(canonical_run, monkeypatch):
    """It raised ``KeyError: 'rno'`` on every model (ledger 476)."""

    from myflopy.viz import Fig

    shown = []
    monkeypatch.setattr(Fig, "show", lambda self, *a, **k: shown.append(self))
    budget = canonical_run.outputs.sfr.bud
    flows = budget.get("flow")
    first = flows["kstpkper"].iloc[0]
    budget.plot_flows()
    (figure,) = shown
    (trace,) = figure.data
    leaving = flows[(flows["kstpkper"] == first) & (flows["q"] < 0)]
    assert list(trace.x) == (leaving["node"] - 1).tolist()
    assert list(trace.y) == pytest.approx((-leaving["q"]).tolist())
    last = flows["kstpkper"].iloc[-1]
    budget.plot_flows(kstpkper=last)
    assert len(shown[-1].data) == 1
    with pytest.raises(ValueError, match="is not a time this SFR budget file saved"):
        budget.plot_flows(kstpkper=[(9, 9)])


def test_an_area_sums_every_layer_of_its_columns(canonical_run, tmp_path):
    """A record names a MODEL-WIDE node; matching it as a cell number read only layer 0.

    The canonical GHB sits in all four layers: an area over every GHB column read a
    quarter of the package total until ledger 482.
    """

    model = canonical_run
    budget = Budget(model, "ghb")
    nodes = budget.df.index.get_level_values("node")
    layers = {model.node_to_lni[int(node)][0] for node in nodes}
    assert len(layers) > 1  # the premise: records below layer 0
    columns = sorted({model.node_to_lni[int(node)][1] for node in nodes})
    area = model.vor.gdf_vorPolys.loc[columns].union_all().buffer(1)
    path = tmp_path / "ghb_columns.gpkg"
    gpd.GeoDataFrame({"name": ["all"]}, geometry=[area], crs=model.vor.crs).to_file(path)
    table = budget.plot_budget_obs(shp_gpkg=path, plot_fig=False)
    pd.testing.assert_series_equal(table["all"], table["total"], check_names=False)

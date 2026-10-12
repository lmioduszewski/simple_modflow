"""Input-field views combine a cell's records by what the number IS, and pair them once.

Ledger 472. Ledger 297 made a single model's input MAP combine a cell's records by
the field's kind (``package_registry.agg_for_kind``); the series panels, the group maps
and the diff still summed every field, so ``drn.inputs.elev.plot()`` drew the SUM of the
drains' elevations and a group map of two ``bhead`` records 5 and 7 drew 12. And
``GroupPackageInputs.compare`` joined a cell's records on the cell, each with each, so
a change of 1.0 in each of two records drew 4.0 through the group and the diff.

The bench: two one-period members on one Voronoi mesh, each with two GHB records
(``bhead`` 5, 7 against 6, 8; ``cond`` 100 each) and two DRN records (``elev`` 3, 4
against 4, 5; ``cond`` 50 each) on one east-edge cell.
"""

from __future__ import annotations

import warnings

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import box

pytestmark = pytest.mark.slow


def _z(picture, ncpl):
    """The per-cell values a map picture draws."""

    for trace in picture.fig.data:
        values = getattr(trace, "z", None)
        if values is not None and len(values) == ncpl:
            return np.asarray(values, dtype=float)
    raise AssertionError("no per-cell trace")


@pytest.fixture(scope="module")
def bench(tmp_path_factory):
    import myflopy as mf

    root = tmp_path_factory.mktemp("by_kind")
    boundary = root / "boundary.gpkg"
    gpd.GeoDataFrame({"n": ["d"]}, geometry=[box(0, 0, 600, 600)], crs=2927).to_file(boundary)
    grid = mf.GridSpec.voronoi(
        boundary=mf.GeoPackageSourceSpec(path=boundary), max_area=40_000.0
    ).resolve(workspace=root / "grid")
    edges = grid.gdf_vorPolys.geometry.bounds
    west = [[(0, int(c)), 10.0] for c in np.flatnonzero(edges.minx <= 1e-6)]
    cell = int(np.flatnonzero(edges.maxx >= 600 - 1e-6)[0])
    props = grid.get_disv_gridprops()

    def member(name, shift):
        model = mf.gwf("m", grid=grid, packages=[
            mf.disv(nlay=1, ncpl=props["ncpl"], nvert=len(props["vertices"]),
                    vertices=props["vertices"], cell2d=props["cell2d"], top=20.0, botm=[0.0]),
            mf.npf(k=10.0), mf.ic(strt=8.0),
            mf.chd(stress_period_data={0: west}),
            mf.ghb(stress_period_data={0: [[(0, cell), 5.0 + shift, 100.0],
                                           [(0, cell), 7.0 + shift, 100.0]]}),
            mf.drn(stress_period_data={0: [[(0, cell), 3.0 + shift, 50.0],
                                           [(0, cell), 4.0 + shift, 50.0]]}),
            mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
                  saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
        ])
        simulation = mf.simulation(model, name="s", packages=[
            mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]), mf.ims(models=["m"])])
        with warnings.catch_warnings():
            # A hand-built grid cannot be frozen into the run; nothing here reopens it.
            warnings.filterwarnings("ignore", message=r"Run 'r' built but its frozen snapshot")
            run = mf.Project(root / name, name=name).run("r", simulation)
        assert run.success, "\n".join(run.report[-20:])
        return mf.load_mf6_run(run.model().workspace)

    a, b = member("pa", 0.0), member("pb", 1.0)
    group = mf.ModelGroup({"a": a, "b": b}, reference="a")
    return a, group, cell


def test_a_single_models_series_averages_an_elevation_and_adds_a_conductance(bench):
    a, _, _ = bench
    assert list(a.drn.inputs.elev.plot().data[0].y) == pytest.approx([3.5])
    assert list(a.drn.inputs.cond.plot().data[0].y) == pytest.approx([100.0])


def test_a_group_map_and_series_average_an_elevation(bench):
    a, group, cell = bench
    ncpl = a.vor.ncpl
    assert _z(group.packages.ghb.inputs.bhead.map(), ncpl)[cell] == pytest.approx(6.0)
    assert _z(group.packages.ghb.inputs.cond.map(), ncpl)[cell] == pytest.approx(200.0)
    series = [list(trace.y) for trace in group.packages.ghb.inputs.bhead.plot().data]
    assert series == [pytest.approx([6.0]), pytest.approx([7.0])]


def test_compare_pairs_records_not_their_cross_product(bench):
    _, group, cell = bench
    rows = group.packages.ghb.inputs.compare(cells=[cell])
    assert len(rows) == 2
    assert rows["bhead_diff"].tolist() == pytest.approx([1.0, 1.0])
    assert rows["cond_diff"].tolist() == pytest.approx([0.0, 0.0])


def test_compare_maps_and_the_diff_draw_one_records_change(bench):
    a, group, cell = bench
    ncpl = a.vor.ncpl
    drawn = group.packages.ghb.inputs.compare_map(model_name="b", value_column="bhead")
    assert _z(drawn, ncpl)[cell] == pytest.approx(1.0)
    diff = group.diff().packages
    assert _z(diff.ghb.inputs.bhead.map("b"), ncpl)[cell] == pytest.approx(1.0)
    assert [list(t.y) for t in diff.drn.inputs.elev.plot().data] == [pytest.approx([1.0])]
    assert [list(t.y) for t in diff.ghb.inputs.plot().data] == [pytest.approx([1.0])]


def test_the_diff_counts_changed_cells_not_rows(bench):
    _, group, _ = bench
    summary = group.diff().packages.ghb.inputs.summary()
    assert summary["value_cells_changed"].tolist() == [1]


def test_the_multiplier_scales_each_record_before_agg(bench):
    """Ledger 473: ``multiplier=-1, agg="max"`` draws the largest NEGATED record."""

    a, _, cell = bench
    drawn = _z(a.drn.inputs.elev.map(multiplier=-1.0, agg="max"), a.vor.ncpl)
    assert drawn[cell] == pytest.approx(-3.0)  # max(-3, -4); after agg it read -4


def test_an_undeclared_column_asks_for_agg(bench):
    from myflopy.modflow.mf6.package_registry import input_field_agg

    assert input_field_agg("drn", "elev") == "mean"
    assert input_field_agg("drn", "cond") == "sum"
    assert input_field_agg("drn", "iface", "max") == "max"
    with pytest.raises(ValueError, match="agg='sum'"):
        input_field_agg("drn", "iface")


def test_uzf_fields_combine_by_kind(canonical_run, monkeypatch):
    """The UZF series and map default to the field's kind, not ``"sum"``."""

    import myflopy.modflow.mf6.package_inputs as package_inputs

    extdp = canonical_run.uzf.inputs.extdp
    table = extdp.get()
    expected = table.groupby("per")["extdp"].mean().to_numpy(float)
    assert np.asarray(extdp.plot().data[0].y, float) == pytest.approx(expected)

    seen = {}
    real = package_inputs.build_cell_input_map_payload

    def spy(*args, **kwargs):
        seen["agg"] = kwargs["agg"]
        return real(*args, **kwargs)

    monkeypatch.setattr(package_inputs, "build_cell_input_map_payload", spy)
    extdp.map().fig
    assert seen["agg"] == "mean"
    canonical_run.uzf.inputs.finf.map().fig
    assert seen["agg"] == "sum"


def test_a_lakes_series_averages_its_connections_elevations(canonical_run):
    """``model.lak.inputs.<field>.plot()`` pools a lake's connections by the field's kind.

    One line per lake over every connection of it: an elevation, a leakance or the
    starting stage averages, the connection area adds. They all added, so the canonical
    lake's ``belev.plot()`` drew 6912.0 for 72 connections at 96.0, and ``strt.plot()``
    72 times its one starting stage (ledger 475).
    """

    for field, agg in (("belev", "mean"), ("strt", "mean"), ("bedleak", "mean"),
                       ("connection_area", "sum")):
        view = getattr(canonical_run.lak.inputs, field)
        table = view.get()
        assert table.groupby("per")["lake"].count().max() > 1, "a lake of several connections"
        expected = table.groupby("per")[field].agg(agg).to_numpy(float)
        assert np.asarray(view.plot().data[0].y, float) == pytest.approx(expected), field

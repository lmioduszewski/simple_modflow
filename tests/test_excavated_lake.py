"""``mf.lak`` on an excavated layer stack: the rule chosen from the RECORD (ledger 325).

The excavation design's C7 (sections 4.1-4.8). What these tests pin:

- the rule comes from the excavation record, found through the SAME resolver as the
  layer elevations (``attach=True``, ``surfaces=<result>``, the grid's own live table),
  never from geometry: a bare table carries no record;
- an excavated lake's VERTICAL connection goes to the top active cell, whose top IS
  the bed, and its walls go INTO the neighbouring columns with MF6's own DISV
  CONNLEN/CONNWIDTH; a rim cell connects at the ground;
- ``lake_bottom`` defaults from the record (F13); a declared one must BE the record's
  bed; a partial match, two excavations, ``rectangular_interior`` and a wrong
  ``only_layer`` are refused whatever ``unexcavated=`` says;
- the legacy rule writes EXACTLY what it wrote before, on five fixtures and the
  canonical lake (``tests/data/lak_legacy_golden.json``, written by the code BEFORE this
  change), except faces thinner than 1e-6 ft, which are dropped (F9);
- the ``unexcavated=`` guard: its numbers, its three words, the "looks excavated" case;
- the build-time check against the built DISV, from BOTH call sites (LAK built before
  and after its DISV), and the four build warnings;
- and MF6 (slow): a pit lake below the water table in uniform flow carries the
  through-flow of the same pit as a K = 1e4 inclusion, to 1%, with the GWF and LAK
  budgets closed -- the design measured 1.005 for this rule and 0.39 for the old one.
"""

from __future__ import annotations

import contextlib
import io
import json
import pathlib
import warnings

import flopy
import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import Point, box

import myflopy as mf
from myflopy import LAKBuilder, ModelContext
from myflopy.grid_identity import LAKE_CONNECTIONS, LAKE_RULES, lak_geometry_complaints
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig
from myflopy.modflow.mf6.grid.geometry import footprint_cells
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.mf6.lakes import FACE_TOL

GOLDEN = pathlib.Path(__file__).parent / "data" / "lak_legacy_golden.json"
GUARD = r"mf\.lak: lake\(s\) NOT excavated"


def _square_grid(n=4, top=None, tops=(12.0, 8.0, 0.0)) -> VoronoiGridPlus:
    """``n`` x ``n`` unit cells (cell ``r * n + c`` centred at ``(c + .5, r + .5)``)."""

    verts = np.array([[x, y] for y in range(n + 1) for x in range(n + 1)], float)
    iverts, centres = [], []
    for r in range(n):
        for c in range(n):
            v0 = r * (n + 1) + c
            iverts.append([v0, v0 + 1, v0 + n + 2, v0 + n + 1])
            centres.append([c + 0.5, r + 0.5])
    grid = VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=np.array(centres))
    columns = {0: np.full(n * n, tops[0]) if top is None else np.asarray(top, float)}
    for k, bottom in enumerate(tops[1:]):
        columns[k + 1] = np.full(n * n, bottom)
    grid.gdf_topbtm = gpd.GeoDataFrame(columns, geometry=grid.gdf_vorPolys.geometry, crs=grid.crs)
    grid.get_disu_connectivity(validate=True)
    return grid


POND = [5, 6, 9, 10]


def _pond(grid):
    return gpd.GeoDataFrame({"name": ["pond"]}, geometry=[box(1, 1, 3, 3)], crs=grid.crs)


def _stack():
    return mf.LayerStack(top=mf.Flat(12.0)).add("upper", bottom=mf.Flat(8.0)).add(
        "lower", bottom=mf.Flat(0.0)
    )


def _bowl_raster(path, grid, beds):
    """A 0.25-unit raster holding ``beds[cell]`` over each unit cell (else 12)."""

    n = 16
    values = np.full((n, n), 12.0)
    for cell, bed in beds.items():
        row, col = divmod(cell, 4)
        # raster row 0 is the TOP (y = 4); cell row r spans y in [r, r + 1]
        values[(3 - row) * 4:(4 - row) * 4, col * 4:(col + 1) * 4] = bed
    with rasterio.open(path, "w", driver="GTiff", height=n, width=n, count=1, dtype="float64",
                       crs=grid.crs, transform=from_origin(0.0, 4.0, 0.25, 0.25),
                       nodata=-9999.0) as out:
        out.write(values, 1)
    return mf.Raster(path)


def _lake(context, lakes, **kw):
    kw.setdefault("starting_stage", 11.0)
    kw.setdefault("bed_leakance", 0.1)
    return LAKBuilder(context=context, nper=1, lakes=lakes, lake_id_field="name", **kw)


@pytest.fixture
def flat(tmp_path):
    """The 2 x 2 pond cut to a flat 9.0 into a 12 / 8 / 0 stack, published with attach=True."""

    grid = _square_grid()
    result = _stack().excavate(bottom=9.0, within=_pond(grid), name="pond").build(grid, attach=True)
    return grid, result


# -- the rule ------------------------------------------------------------------------------


def test_an_excavated_lake_connects_at_its_bed_and_into_its_neighbours(flat):
    grid, result = flat
    lake = _lake(ModelContext(grid=grid, domain=result.idomain), _pond(grid))
    assert lake.lake_rules == {"pond": "excavated"}
    rows = lake.connectiondata
    vertical = sorted((tuple(r[2]), r[5], r[6]) for r in rows if r[3] == "VERTICAL")
    assert vertical == [((0, c), 9.0, 9.0) for c in POND], "top active cell, its top IS the bed"
    walls = [r for r in rows if r[3] == "HORIZONTAL"]
    # Each corner cell of the 2 x 2 pond has two LAND neighbours; every wall is written
    # INTO the neighbour, from the bed to the ground, one per active layer it spans.
    assert len(walls) == 8
    assert {r[2][1] for r in walls}.isdisjoint(POND), "walls live in the neighbours' cells"
    assert {(r[2][0], r[5], r[6]) for r in walls} == {(0, 9.0, 12.0)}
    # MF6's own DISV geometry: the neighbour's centre is 0.5 from the shared face, which is
    # 1.0 long -- NOT grid.cl12's centre-to-centre 1.0.
    assert {(r[7], r[8]) for r in walls} == {(0.5, 1.0)}
    roles = lake.connection_roles["pond"]
    assert [role for _, _, role in roles].count("wall") == 8
    owner = {tuple(r[2]): cell for r, (_, cell, role) in zip(rows, roles, strict=True) if role == "wall"}
    assert set(owner.values()) <= set(POND) and owner[(0, 1)] == 5 and owner[(0, 4)] == 5


def test_a_bowl_steps_into_its_own_higher_cells_and_a_rim_connects_at_the_ground(tmp_path):
    grid = _square_grid()
    beds = {5: 9.0, 6: 10.0, 9: 10.0, 10: 13.0}      # cell 10's bed is ABOVE the ground
    bed = _bowl_raster(tmp_path / "bowl.tif", grid, beds)
    result = _stack().excavate(bottom=bed, within=_pond(grid), name="pond").build(grid, attach=True)
    lake = _lake(ModelContext(grid=grid, domain=result.idomain), _pond(grid))
    rows, roles = lake.connectiondata, lake.connection_roles["pond"]
    vertical = {r[2][1]: r[5] for r in rows if r[3] == "VERTICAL"}
    assert vertical == {5: 9.0, 6: 10.0, 9: 10.0, 10: 12.0}, "the rim at the ground"
    assert [role for _, cell, role in roles if role != "wall"] == ["floor", "floor", "floor", "rim"]
    # Cell 5 (bed 9) steps into its lake neighbours 6 and 9 between 9 and their beds (10)
    into = {(r[2][1], r[5], r[6]) for r, (_, cell, role) in zip(rows, roles, strict=True)
            if role == "wall" and cell == 5}
    assert (6, 9.0, 10.0) in into and (9, 9.0, 10.0) in into
    assert (1, 9.0, 12.0) in into and (4, 9.0, 12.0) in into
    # The rim (bed == ground) writes no wall: nothing around it is higher.
    assert not [cell for _, cell, role in roles if role == "wall" and cell == 10]


def test_the_record_is_found_where_the_elevations_are(flat):
    """Through ``attach=True`` (the grid), ``surfaces=<result>``, and the grid's own LIVE
    table; a COPY of that table is a bare table and carries no record."""

    grid, result = flat
    lakes = _pond(grid)
    routes = {
        "attach": ModelContext(grid=grid, domain=result.idomain),
        "surfaces": ModelContext(grid=grid, domain=result.idomain, surfaces=result),
        "live table": ModelContext(grid=grid, domain=result.idomain, surfaces=grid.gdf_topbtm),
    }
    for route, context in routes.items():
        assert _lake(context, lakes).lake_rules == {"pond": "excavated"}, route
    bare = ModelContext(grid=grid, domain=result.idomain, surfaces=grid.gdf_topbtm.copy())
    assert _lake(bare, lakes, lake_bottom=9.0, unexcavated="ignore").lake_rules == {"pond": "legacy"}


def test_lake_bottom_defaults_from_the_record_and_a_declared_one_must_be_it(flat):
    grid, result = flat
    context = ModelContext(grid=grid, domain=result.idomain)
    spec = mf.lak(context=context, nper=1, lakes=_pond(grid), lake_id_field="name",
                  starting_stage=11.0, unexcavated="raise")
    assert spec.metadata[LAKE_RULES] == {"pond": "excavated"}
    assert _lake(context, _pond(grid), lake_bottom=9.0).connectiondata   # the same bed
    with pytest.raises(ValueError, match="not the bed excavation 'pond' was cut to") as error:
        _lake(context, _pond(grid), lake_bottom=9.5).connectiondata
    assert "4 cells differ" in str(error.value) and "lake_bottom 9.5, excavation 9" in str(error.value)


def test_a_lake_that_is_only_part_of_an_excavation_is_refused_whatever_unexcavated_says(flat):
    grid, result = flat
    context = ModelContext(grid=grid, domain=result.idomain)
    part = gpd.GeoDataFrame({"name": ["pond"]}, geometry=[box(1, 1, 3, 2)], crs=grid.crs)
    for word in ("warn", "raise", "ignore"):
        with pytest.raises(ValueError, match="cover different cells") as error:
            _lake(context, part, unexcavated=word).connectiondata
        assert "2 dug cell(s) no lake of this package covers (9, 10)" in str(error.value)
    wide = gpd.GeoDataFrame({"name": ["pond"]}, geometry=[box(1, 1, 4, 3)], crs=grid.crs)
    with pytest.raises(ValueError, match=r"2 lake cell\(s\) the excavation did not dig \(7, 11\)"):
        _lake(context, wide).connectiondata


def test_two_lakes_may_share_one_excavation_but_one_lake_may_not_span_two():
    grid = _square_grid()
    two = gpd.GeoDataFrame({"name": ["west", "east"]},
                           geometry=[box(0, 0, 1, 4), box(3, 0, 4, 4)], crs=grid.crs)
    result = _stack().excavate(bottom=9.0, within=two, name="both").build(grid, attach=True)
    context = ModelContext(grid=grid, domain=result.idomain)
    assert _lake(context, two).lake_rules == {"west": "excavated", "east": "excavated"}
    split = (_stack().excavate(bottom=9.0, within=two.iloc[[0]], name="w")
             .excavate(bottom=9.0, within=two.iloc[[1]], name="e").build(grid, attach=True))
    whole = gpd.GeoDataFrame({"name": ["all"]}, geometry=[box(0, 0, 4, 4)], crs=grid.crs)
    with pytest.raises(ValueError, match="covers cells of 2 excavations"):
        _lake(ModelContext(grid=grid, domain=split.idomain), whole).lake_rules


def test_rectangular_interior_and_a_wrong_only_layer_are_refused_on_an_excavated_lake(flat):
    grid, result = flat
    context = ModelContext(grid=grid, domain=result.idomain)
    rect = dict(connection_modes="rectangular", lake_top=12.0)
    rows = _lake(context, _pond(grid), **rect).connectiondata
    assert sum(r[3] == "HORIZONTAL" for r in rows) == 8, "edges only, by itself"
    with pytest.raises(ValueError, match="rectangular_interior=True on excavated lake"):
        _lake(context, _pond(grid), rectangular_interior=True, **rect).connectiondata
    with pytest.raises(ValueError, match="only_layer=1 for excavated lake 'pond'"):
        _lake(context, _pond(grid), only_layer=1).connectiondata
    assert _lake(context, _pond(grid), only_layer=0).connectiondata


def test_an_absorbed_sliver_puts_the_floor_on_the_cell_below(tmp_path):
    """A 0.05 ft remnant of the upper layer is absorbed (ledger 320): the connected cell is
    the LOWER layer, with its top at the bed, and the vertical goes there."""

    grid = _square_grid()
    result = _stack().excavate(bottom=8.05, within=_pond(grid), name="pond").build(grid, attach=True)
    assert (result.idomain[0, POND] == 0).all()
    for context in (ModelContext(grid=grid, domain=result.idomain), ModelContext(grid=grid)):
        # With NO domain every layer reads active, and the zero-thickness layer above the
        # bed "contains" it too (top == bottom == bed): the floor still goes to the cell
        # that has THICKNESS (the design's section 4.2).
        lake = _lake(context, _pond(grid), starting_stage=9.0)
        assert {tuple(r[2]) for r in lake.connectiondata if r[3] == "VERTICAL"} == {(1, c) for c in POND}
        assert {r[5] for r in lake.connectiondata if r[3] == "VERTICAL"} == {8.05}


# -- the legacy rule and its guard -------------------------------------------------------


def test_the_legacy_rule_writes_what_it_wrote_before():
    """Byte-identical to the rule BEFORE ledger 325 (``tests/data/lak_legacy_golden.json``,
    written by that code), on a bowl, a rectangular basin with and without interior faces,
    a bed equal to its model top (a USG export's lake) and a flat floor over sloping
    ground. The canonical lake and a 612-cell copy of a real pond were compared the same
    way when this landed (``docs/reproducers/excavated_lake_rule.py``)."""

    golden = json.loads(GOLDEN.read_text())
    grid = _square_grid()
    lakes = _pond(grid)
    context = ModelContext(grid=grid, domain=np.ones((2, 16), int))
    bed = pd.Series({c: 12.0 for c in range(16)})
    bed[5], bed[6], bed[9], bed[10] = 9.0, 10.0, 10.0, 11.0
    base = dict(starting_stage=11.5, bed_leakance=0.1, unexcavated="ignore")
    rect = dict(connection_modes="rectangular", lake_bottom=9.0, lake_top=11.0)
    grid3 = _square_grid(top=12.0 + 0.25 * (np.arange(16) % 4) + 0.1 * (np.arange(16) // 4))
    top3 = np.asarray(grid3.gdf_topbtm[0], float)
    grid4 = _square_grid(tops=(27.0, 20.0, 0.0), top=25.0 + 2.0 * (np.arange(16) % 4) / 3.0)
    built = {
        "toy_bathy": _lake(context, lakes, lake_bottom=bed, **base),
        "toy_rect": _lake(context, lakes, **rect, **base),
        "toy_rect_interior": _lake(context, lakes, rectangular_interior=True, **rect, **base),
        "bed_is_top": _lake(ModelContext(grid=grid3, domain=np.ones((2, 16), int)), lakes,
                            lake_bottom=pd.Series({c: top3[c] for c in range(16)}),
                            **{**base, "starting_stage": 13.0}),
        "flat_floor_sloping": _lake(ModelContext(grid=grid4, domain=np.ones((2, 16), int)), lakes,
                                    connection_modes="rectangular", lake_bottom=25.5, lake_top=28.0,
                                    **{**base, "starting_stage": 26.0}),
    }
    for name, builder in built.items():
        assert builder.lake_rules == {"pond": "legacy"}, name
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            rows = [[r[0], r[1], list(r[2]), r[3], r[4], r[5], r[6], r[7], r[8]]
                    for r in builder.connectiondata]
        assert rows == golden[name], name


@pytest.mark.canonical
def test_the_canonical_lake_writes_what_it_wrote_before(canonical_model, canonical_config):
    """The canonical rectangular lake (unexcavated, F8) -- 72 connections, identical."""

    if canonical_config != CanonicalModelConfig.testing():
        pytest.skip("the golden rows are the TESTING profile's")
    rows = [[int(r[0]), int(r[1]), [int(r[2][0]), int(r[2][1])], str(r[3]), float(r[4]),
             float(r[5]), float(r[6]), float(r[7]), float(r[8])]
            for r in canonical_model.gwf.get_package("lak").connectiondata.get_data().tolist()]
    assert rows == json.loads(GOLDEN.read_text())["canonical_testing"]


def test_a_face_thinner_than_the_tolerance_is_not_written():
    """F9: MF6 reads ``belev == telev`` as the WHOLE cell, and 8-decimal output makes a 1e-9
    face exactly that -- so a step thinner than 1e-6 writes no face (it did before)."""

    grid = _square_grid()
    context = ModelContext(grid=grid, domain=np.ones((2, 16), int))
    lakes = gpd.GeoDataFrame({"name": ["pond"]}, geometry=[box(1, 1, 2, 2)], crs=grid.crs)
    bed = pd.Series({c: 12.0 for c in range(16)})
    bed[5] = 11.0
    bed[6] = 11.0 + 0.4 * FACE_TOL      # a step of 4e-7 up to the east neighbour
    bed[4] = 11.0 + 2.0 * FACE_TOL      # 2e-6 to the west: kept
    rows = _lake(context, lakes, lake_bottom=bed, unexcavated="ignore").connectiondata
    neighbours = sorted(round(r[6] - r[5], 9) for r in rows if r[3] == "HORIZONTAL")
    assert neighbours == [2e-6, 1.0, 1.0], "the 4e-7 face to cell 6 is gone, the 2e-6 one kept"


def test_the_guard_reports_a_lake_mf6_counts_from_its_cell_tops(recwarn):
    grid = _square_grid()
    context = ModelContext(grid=grid, domain=np.ones((2, 16), int))
    kw = dict(connection_modes="rectangular", lake_bottom=9.0, lake_top=11.0, starting_stage=10.0)
    with pytest.warns(UserWarning, match=GUARD) as caught:
        mf.lak(context=context, nper=1, lakes=_pond(grid), lake_id_field="name", **kw)
    text = str(caught[0].message)
    assert "4 of its 4 VERTICAL connections sit in a cell whose TOP is above the bed, by up to 3" in text
    assert "(median 3); 4 of those tops are above its starting stage 10" in text
    assert "stack.excavate(bottom=<the same lake_bottom>" in text
    # ...and the route for a model with no stack (the design's F14).
    assert "dug = result.excavate(bottom=<the same lake_bottom>" in text
    assert "mf.disv(dug)" in text and "ModelContext(surfaces=dug, domain=dug.idomain)" in text
    with pytest.raises(ValueError, match=GUARD):
        mf.lak(context=context, nper=1, lakes=_pond(grid), lake_id_field="name",
               unexcavated="raise", **kw)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        mf.lak(context=context, nper=1, lakes=_pond(grid), lake_id_field="name",
               unexcavated="ignore", **kw)
    with pytest.raises(ValueError, match="unexcavated must be one of"):
        _lake(context, _pond(grid), unexcavated="nearest", **kw)


def test_a_lake_whose_bed_is_its_cell_tops_is_not_reported():
    """A USG export's lakes: the bed IS the model top, so MF6 counts from the bed."""

    grid = _square_grid(top=12.0 + 0.25 * (np.arange(16) % 4))
    top = np.asarray(grid.gdf_topbtm[0], float)
    context = ModelContext(grid=grid, domain=np.ones((2, 16), int))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        mf.lak(context=context, nper=1, lakes=_pond(grid), lake_id_field="name",
               starting_stage=13.0, lake_bottom=pd.Series({c: top[c] for c in range(16)}),
               unexcavated="raise")
    assert not [w for w in caught if "excavat" in str(w.message)]


def test_a_lake_that_looks_excavated_but_carries_no_record_is_reported(flat):
    grid, result = flat
    table = grid.gdf_topbtm.copy()               # the elevations WITHOUT their record
    absorbed = _stack().excavate(bottom=8.05, within=_pond(grid), name="pond").build(grid)
    table = gpd.GeoDataFrame({0: absorbed.top, 1: absorbed.botm[0], 2: absorbed.botm[1]},
                             geometry=grid.gdf_vorPolys.geometry, crs=grid.crs)
    context = ModelContext(grid=grid, domain=absorbed.idomain, surfaces=table)
    with pytest.raises(ValueError, match="LOOKS excavated -- 4 of its 4 cells") as caught:
        mf.lak(context=context, nper=1, lakes=_pond(grid), lake_id_field="name",
               starting_stage=9.0, lake_bottom=8.05, unexcavated="raise")
    # Hand-built elevations are exactly the case with no stack to follow the advice on:
    # it names the BUILT result's own excavate too (the design's F14).
    assert "result.excavate(bottom=<the same lake_bottom>" in str(caught.value)


# -- the metadata --------------------------------------------------------------------------


def test_the_rules_and_roles_ride_in_metadata_and_an_edit_of_the_connections_drops_them(flat):
    grid, result = flat
    spec = mf.lak(context=ModelContext(grid=grid, domain=result.idomain), nper=1,
                  lakes=_pond(grid), lake_id_field="name", starting_stage=11.0)
    assert spec.metadata[LAKE_RULES] == {"pond": "excavated"}
    roles = spec.metadata[LAKE_CONNECTIONS]["pond"]
    assert len(roles) == len(spec.options["connectiondata"]) == 12
    assert all(type(i) is int and type(c) is int and type(r) is str for i, c, r in roles)
    edited = spec.with_options(connectiondata=spec.options["connectiondata"][:4])
    assert LAKE_RULES not in edited.metadata and LAKE_CONNECTIONS not in edited.metadata
    assert LAKE_RULES in spec.with_options(print_stage=True).metadata
    assert mf.PackageSpec.from_dict(spec.to_dict()).metadata[LAKE_RULES] == {"pond": "excavated"}


# -- the build-time check and warnings -------------------------------------------------------


def _model(context, result, lake_spec, *, lak_first=False, newton=True, icelltype=1, rch=False):
    packages = [mf.disv(result), mf.npf(k=10.0, k33=1.0, icelltype=icelltype), mf.ic(strt=10.0),
                mf.chd(stress_period_data={0: [[(1, 0), 10.0]]})]
    if rch:
        packages.append(mf.rch.array(context=context, recharge=1e-3))
    packages.insert(0 if lak_first else len(packages), lake_spec)
    return mf.gwf("m", context=context, packages=packages,
                  **({"newtonoptions": "NEWTON"} if newton else {}))


def _build(model, workspace):
    simulation = mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    return simulation.build(workspace)


@pytest.mark.parametrize("lak_first", [False, True], ids=["after-disv", "before-disv"])
def test_a_lake_resolved_on_the_cut_is_refused_against_an_uncut_disv(flat, lak_first, tmp_path):
    """Both call sites: the LAK builder (DISV already built -- FloPy would otherwise stop
    naming only a cellid) and ``ModelSpec.build``'s own check (LAK declared first)."""

    grid, result = flat
    context = ModelContext(grid=grid, domain=result.idomain)
    lake = mf.lak(context=context, nper=1, lakes=_pond(grid), lake_id_field="name",
                  starting_stage=11.0)
    uncut = _stack().build(grid)
    with pytest.raises(ValueError, match="resolved against a different geometry") as error:
        _build(_model(ModelContext(grid=grid, domain=uncut.idomain), uncut, lake, lak_first=lak_first),
               tmp_path / "uncut")
    text = str(error.value)
    assert "excavated lake 'pond' was resolved with its VERTICAL connections at the excavated bed" in text
    site = "does not agree with the grid it was built" if lak_first else "does not fit the DISV it is being built with"
    assert site in text, "the check ran at the call site this order reaches"
    # ...and on the DISV it was resolved against, it builds.
    _build(_model(context, result, lake, lak_first=lak_first), tmp_path / "cut")


def test_a_lake_naming_a_dug_cell_is_refused_before_flopy_can_name_only_the_cellid(tmp_path):
    """The reverse: a lake resolved on the UNCUT stack, built with a DISV that absorbed the
    pond's upper cells. FloPy would stop at construction with an ``MFDataException`` naming
    a cellid; the builder's check stops first and names the cause."""

    grid = _square_grid()
    cut = _stack().excavate(bottom=8.05, within=_pond(grid), name="pond").build(grid)
    uncut = _stack().build(grid)
    lake = mf.lak(context=ModelContext(grid=grid, domain=uncut.idomain), nper=1, lakes=_pond(grid),
                  lake_id_field="name", starting_stage=11.0, lake_bottom=9.0, unexcavated="ignore")
    with pytest.raises(ValueError, match="does not fit the DISV it is being built with") as error:
        _build(_model(ModelContext(grid=grid, domain=cut.idomain), cut, lake), tmp_path / "m")
    assert "lake 'pond' connects to 4 cell(s) the built DISV has REMOVED" in str(error.value)


def test_the_geometry_check_names_a_removed_cell_and_a_face_outside_its_cell():
    idomain = np.array([[1, 0], [1, 1]])
    top, botm = np.array([12.0, 12.0]), np.array([[8.0, 8.0], [0.0, 0.0]])
    rows = [[0, 0, (0, 1), "VERTICAL", 1.0, 12.0, 12.0, 0.0, 0.0],
            [0, 1, (1, 0), "HORIZONTAL", 1.0, 7.0, 9.0, 0.5, 1.0],
            [0, 2, (1, 1), "HORIZONTAL", 1.0, 0.0, 0.0, 0.5, 1.0]]   # belev == telev: whole cell
    found = lak_geometry_complaints(idomain, top, botm, rows, {"pond": "legacy"})
    assert len(found) == 2
    assert "connects to 1 cell(s) the built DISV has REMOVED" in found[0] and "(0, 1)" in found[0]
    assert "1 HORIZONTAL connection(s) reaching outside their cell" in found[1]
    assert "written 7-9 in a cell 0-8" in found[1]


def test_the_build_warns_about_newton_icelltype_recharge_and_a_transparent_bed(flat, tmp_path):
    grid, result = flat
    context = ModelContext(grid=grid, domain=result.idomain)
    lake = mf.lak(context=context, nper=1, lakes=_pond(grid), lake_id_field="name",
                  starting_stage=11.0, bed_leakance=3000.0)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _build(_model(context, result, lake, newton=False, icelltype=0, rch=True), tmp_path / "a")
    text = "\n".join(str(w.message) for w in caught)
    assert "the model does not use NEWTON" in text
    assert "8 wall connection(s) of excavated lake(s) 'pond' enter cells with icelltype 0" in text
    assert "puts recharge on 4 of excavated lake 'pond''s 4 VERTICAL cell(s)" in text
    # none of them with NEWTON, icelltype 1, rainfall, and no recharge
    rained = mf.lak(context=context, nper=1, lakes=_pond(grid), lake_id_field="name",
                    starting_stage=11.0, bed_leakance=3000.0, rainfall=1e-3)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _build(_model(context, result, rained, rch=True), tmp_path / "b")
    assert not [w for w in caught if "excavated lake" in str(w.message)]


def test_a_transparent_bed_over_uneven_cells_warns(tmp_path):
    """The half-cell term 0.5*thk/K33 dominates 1/bedleak, and the connected cells run from
    0.2 to 3.9 ft thick: the mesh, not the bed, sets the rate."""

    grid = _square_grid()
    beds = {5: 8.2, 6: 11.9, 9: 11.9, 10: 11.9}
    bed = _bowl_raster(tmp_path / "bowl.tif", grid, beds)
    stack = mf.LayerStack(top=mf.Flat(12.0)).add("upper", bottom=mf.Flat(8.0)).add(
        "lower", bottom=mf.Flat(0.0))
    result = stack.excavate(bottom=bed, within=_pond(grid), name="pond").build(grid, attach=True)
    context = ModelContext(grid=grid, domain=result.idomain)
    lake = mf.lak(context=context, nper=1, lakes=_pond(grid), lake_id_field="name",
                  starting_stage=12.0, bed_leakance=3000.0)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _build(_model(context, result, lake), tmp_path / "t")
    text = "\n".join(str(w.message) for w in caught)
    assert "excavated lake 'pond''s bed is TRANSPARENT" in text
    assert "0.2-3.9 thick" in text


# -- MF6 ----------------------------------------------------------------------------------


def _pit_run(root, grid, result, k, lake=None):
    xy = np.array([(c[1], c[2]) for c in grid.get_disv_gridprops()["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    west, east = np.flatnonzero(xy[:, 0] < xmin + 20), np.flatnonzero(xy[:, 0] > xmax - 20)
    context = ModelContext(grid=grid, domain=result.idomain, surfaces=result)
    chd = [[(k_, int(c)), h] for cells, h in ((west, 97.0), (east, 93.0)) for c in cells
           for k_ in range(2) if result.idomain[k_, c] > 0]
    packages = [mf.disv(result), mf.npf(k=k, k33=k, icelltype=1), mf.ic(strt=95.0),
                mf.chd(stress_period_data={0: chd}),
                mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
                      saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])]
    if lake is not None:
        packages.append(lake(context))
    model = mf.gwf("m", context=context, newtonoptions="NEWTON", packages=packages)
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
        mf.ims(models=["m"], complexity="COMPLEX", outer_dvclose=1e-8, inner_dvclose=1e-10,
               outer_maximum=300, inner_maximum=300, linear_acceleration="BICGSTAB")])
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        run = mf.Project(root, name="p").run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])
    workspace = pathlib.Path(run.workspace)
    listing = (workspace / "m.lst").read_text()
    discrepancy = [abs(float(line.split()[-1])) for line in listing.splitlines()
                   if "PERCENT DISCREPANCY" in line]
    return workspace, discrepancy


@pytest.mark.slow
def test_mf6_a_pit_lake_carries_the_through_flow_of_a_conductive_inclusion(tmp_path):
    """A pit below the water table in uniform flow is, hydraulically, an (almost) infinitely
    conductive inclusion. Reference: no lake, K = 1e4 in the pit's layer-1 cells. The
    excavated rule carries 1.005 of the reference's through-flow (the design measured
    1.005-1.006 on two grids; the old rule on the un-excavated stack, 0.39-0.42), and the
    percent discrepancy is asserted for BOTH budgets -- Normal termination alone proves
    nothing (the design's section 5.5)."""

    grid = irregular_voronoi_grid(nrow=40, ncol=40, cell_size=20.0, refined=False)
    xmin, ymin, xmax, ymax = grid.gdf_vorPolys.total_bounds
    pond = gpd.GeoDataFrame({"name": ["pit"]}, crs=grid.crs,
                            geometry=[Point((xmin + xmax) / 2, (ymin + ymax) / 2).buffer(80.0, 64)])
    stack = (mf.LayerStack(top=mf.Flat(100.0)).add("upper", bottom=mf.Flat(82.0))
             .add("lower", bottom=mf.Flat(40.0)))
    plain = stack.build(grid)
    cells = footprint_cells(grid.gdf_vorPolys, pond.geometry)[0]
    k = np.full((2, grid.ncpl), 10.0)
    k_ref = k.copy()
    k_ref[0, cells] = 1.0e4

    workspace, discrepancy = _pit_run(tmp_path / "ref", grid, plain, k_ref)
    grb = flopy.mf6.utils.MfGrdFile(str(workspace / "m.disv.grb"))
    ia, ja = np.asarray(grb.ia), np.asarray(grb.ja)
    flows = np.asarray(flopy.utils.CellBudgetFile(str(workspace / "m.cbc"), precision="double")
                       .get_data(text="FLOW-JA-FACE")[-1]).ravel()
    inside = set(cells)
    reference = sum(flows[p] for n in cells for p in range(ia[n] + 1, ia[n + 1])
                    if ja[p] not in inside and flows[p] > 0)
    assert max(discrepancy) < 0.1

    pit = stack.excavate(bottom=82.0, within=pond, name="pit").build(grid)
    lake = lambda context: mf.lak(context=context, nper=1, lakes=pond, lake_id_field="name",  # noqa: E731
                                  starting_stage=95.0, bed_leakance=3000.0, unexcavated="raise")
    workspace, discrepancy = _pit_run(tmp_path / "lake", grid, pit, k, lake)
    budget = flopy.utils.CellBudgetFile(str(workspace / "m_lak_budget.lak"), precision="double")
    q = np.asarray(budget.get_data(text="GWF")[-1]["q"])
    assert len(discrepancy) == 2 and max(discrepancy) < 0.1, "GWF and LAK budgets both close"
    ratio = float(q[q > 0].sum()) / reference
    assert ratio == pytest.approx(1.0, abs=0.01)
    assert abs(float(q.sum())) < 1e-3 * float(q[q > 0].sum()), "a steady lake passes what it takes"

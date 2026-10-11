"""A region group is all or nothing, and a region's tuple cellid is the grid's
own -- ``(layer, row, col)`` on a structured grid (ledgers 375, 376).

The structured model is a real MF6 run (2 layers, 3 x 4 DIS) opened with
``load_mf6_run``, the route that brings a DIS model to the region registry.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

import myflopy as mf
from myflopy.modflow.mf6.simulation.regions import (
    ModelRegion,
    RegionRegistry,
    normalize_region_inputs,
)
from tests.test_group_same_name_runs import NCOL, _dis_run


@pytest.fixture(scope="module")
def dis_run(tmp_path_factory):
    return _dis_run(tmp_path_factory.mktemp("dis_region") / "base")


# --- 375: add_group is all or nothing ------------------------------------


@pytest.fixture
def registry():
    regions = RegionRegistry(model=None)
    for name in ("a", "b", "c"):
        regions.add(ModelRegion(name=name, cells=[0]))
    return regions


def test_a_refused_member_leaves_no_group_behind(registry):
    with pytest.raises(KeyError):
        registry.add_group("g", members=["a", "typo", "b"])
    assert "g" not in registry
    assert registry.list_groups() == []


def test_a_refused_overwrite_keeps_the_group_it_would_have_replaced(registry):
    kept = registry.add_group("h", members=["a", "b"])
    with pytest.raises(KeyError):
        registry.add_group("h", members=["c", "nope"], overwrite=True)
    assert registry.get_group("h") is kept
    assert kept.members == ["a", "b"]
    assert registry.resolve_cells("h") == [0]


def test_a_cycle_through_the_replaced_group_is_refused_and_rolled_back(registry):
    kept = registry.add_group("h", members=["a"])
    registry.add_group("outer", members=["h"])
    with pytest.raises(ValueError, match="cyclic"):
        registry.add_group("h", members=["b", "outer"], overwrite=True)
    assert registry.get_group("h") is kept and kept.members == ["a"]


def test_an_accepted_overwrite_replaces_the_group(registry):
    registry.add_group("h", members=["a"])
    replaced = registry.add_group("h", members=["b", "c"], overwrite=True)
    assert registry.get_group("h") is replaced and replaced.members == ["b", "c"]


# --- 376: a DIS cellid is (layer, row, col) -------------------------------


def test_a_dis_region_reads_layer_row_col_as_the_grids_own_cell(dis_run):
    model = mf.load_mf6_run(dis_run)
    region = model.add_region("south_east", cellids=[(0, 2, 3), (0, 1, 3)])

    # Row-major, as FloPy's grid, the cell polygons and the heads table number it.
    assert region.cells == [1 * NCOL + 3, 2 * NCOL + 3]
    assert region.layer == 0
    heads = model.region_heads("south_east")
    assert sorted(heads["cell"]) == [7, 11]
    # The east column is the CHD at 1.0, and cell 11's polygon is the grid's
    # south-east corner.
    assert heads["elev"].tolist() == pytest.approx([1.0, 1.0])
    corner = heads.loc[heads["cell"] == 11, "geometry"].iloc[0]
    assert corner.bounds == (30.0, 0.0, 40.0, 10.0)


def test_a_dis_region_refuses_a_two_item_cellid_naming_the_form(dis_run):
    model = mf.load_mf6_run(dis_run)
    with pytest.raises(ValueError, match=r"\(layer, row, col\)"):
        model.add_region("r", cellids=[(0, 3)])
    with pytest.raises(ValueError, match="off the grid"):
        model.add_region("r", cellids=[(0, 3, 0)])
    assert "r" not in model.regions
    # A bare in-layer cell number means the same on every grid.
    assert model.add_region("r", cellids=[11]).cells == [11]


def test_a_vertex_grid_region_refuses_a_three_item_cellid():
    with pytest.raises(ValueError, match=r"\(layer, cell\)"):
        normalize_region_inputs([(0, 3, 7)])
    model = SimpleNamespace(modelgrid=SimpleNamespace(grid_type="vertex"))
    regions = RegionRegistry(model)
    with pytest.raises(ValueError, match=r"\(layer, cell\)"):
        regions.add_from_cells("r", [(0, 3, 7)])
    assert regions.add_from_cells("r", [(1, 3), (1, 5)]).cells == [3, 5]


def test_a_region_of_no_cells_is_still_an_empty_region():
    # Reading the cellids for a tuple must not make ``None`` -- no cells, which
    # ``normalize_region_inputs`` has always accepted -- a TypeError.
    regions = RegionRegistry(model=None)
    assert regions.add_from_cells("empty", None).cells == []

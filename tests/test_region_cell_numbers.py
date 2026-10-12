"""A region's cells and layers are whole numbers, and the registry holds regions (ledger 478).

``add_from_cells`` passed a non-tuple cellid through ``int()`` -- ``5.7`` was cell 5,
``"5"`` was accepted -- ``RegionRegistry.add`` registered any object, a ``RegionGroup``
included, which then broke ``summary()``, and ``resolve_cells("unknown")`` said
"Group" whichever was meant.
"""

from __future__ import annotations

import numpy as np
import pytest

from myflopy.modflow.mf6.simulation.regions import ModelRegion, RegionGroup, RegionRegistry
from myflopy.modflow.utils.validators import whole_number


@pytest.mark.parametrize("cellid", [5.7, "5", True, float("nan"), (2.5,)])
def test_a_cell_that_is_not_a_whole_number_is_refused(cellid):
    regions = RegionRegistry(model=None)
    with pytest.raises(TypeError, match="A cell is a whole number"):
        regions.add_from_cells("r", [cellid])
    assert regions.list() == []


def test_an_integral_float_and_a_numpy_integer_are_cells():
    regions = RegionRegistry(model=None)
    region = regions.add_from_cells("r", [5.0, np.int64(7), (np.float64(9.0),)])
    assert region.cells == [5, 7, 9]
    assert all(type(cell) is int for cell in region.cells)


def test_a_layer_that_is_not_a_whole_number_is_refused():
    regions = RegionRegistry(model=None)
    with pytest.raises(TypeError, match="A layer is a whole number, not 1.7"):
        regions.add_from_cells("r", [5], layer=1.7)
    assert regions.add_from_cells("r", [5], layer=np.int64(1)).layer == 1


def test_add_takes_a_region_only():
    regions = RegionRegistry(model=None)
    with pytest.raises(TypeError, match="A group is add_group's"):
        regions.add(RegionGroup(name="g"))
    regions.add(ModelRegion(name="a", cells=[1]))
    assert list(regions.summary()["name"]) == ["a"]


def test_an_unknown_name_is_neither_a_region_nor_a_group():
    with pytest.raises(KeyError, match="neither a region nor a group"):
        RegionRegistry(model=None).resolve_cells("unknown")


def test_whole_number_is_the_one_rule():
    assert whole_number(np.int32(3), what="cell") == 3
    with pytest.raises(TypeError):
        whole_number(np.inf, what="layer")

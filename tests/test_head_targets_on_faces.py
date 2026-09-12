"""A head target on a cell face is ONE observation, in one cell.

``HeadTargets.match_to_model`` joined its targets to the grid with
``intersects`` and never reduced the result, so a well on a shared face became
two rows and one on a vertex four -- all under the same name. Every consumer
assumes one row per name, so ``compare`` returned five rows for two targets,
``to_flopy_obs`` wrote five MF6 records, and PEST disagreed with itself: the
build read the first cell, the forward run the last. Unmatched targets were
snapped to the nearest cell at any distance, silently, and a distance tie
crashed the whole call.
"""

from __future__ import annotations

import logging
from types import SimpleNamespace

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import Point, box

from myflopy.modflow.mf6.grid import geometry as grid_geometry
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.mf6.observations import heads as heads_module
from myflopy.modflow.mf6.observations.conc import ConcTargets
from myflopy.modflow.mf6.observations.heads import HeadTargets

LOGGER = "myflopy.modflow.mf6.observations.heads"


def _grid(n: int = 4, length: float = 100.0, order: list[int] | None = None) -> VoronoiGridPlus:
    """n x n square cells, row-major from the bottom left (cell 5 is 25..50, 25..50)."""

    d = length / n
    verts = np.array(
        [[i * d, j * d] for j in range(n + 1) for i in range(n + 1)], dtype=float
    )
    iverts = [
        [j * (n + 1) + i, (j + 1) * (n + 1) + i, (j + 1) * (n + 1) + i + 1, j * (n + 1) + i + 1]
        for j in range(n)
        for i in range(n)
    ]
    xcyc = [[(i + 0.5) * d, (j + 0.5) * d] for j in range(n) for i in range(n)]
    if order is not None:
        iverts = [iverts[i] for i in order]
        xcyc = [xcyc[i] for i in order]
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=np.array(xcyc))


def _model(grid, idomain=None):
    """A model-like object: the grid, and a flopy-shaped DISV idomain when given."""

    disv = None if idomain is None else SimpleNamespace(
        idomain=SimpleNamespace(array=np.asarray(idomain, dtype=int))
    )
    return SimpleNamespace(vor=grid, gwf=SimpleNamespace(disv=disv))


def _targets(named, *, crs="EPSG:2927", layers=None, cls=HeadTargets):
    """Point targets, one value per name so the tables line up."""

    names = [name for name, _ in named]
    frame = gpd.GeoDataFrame(
        {"name": names, "layer": layers if layers is not None else [0] * len(named)},
        geometry=[geometry for _, geometry in named],
        crs=crs,
    )
    values = pd.DataFrame({"per": [0], **{name: [1.0] for name in names}})
    return cls(locations=frame, values=values, time_column="per")


def _cell_targets(named, layers=None):
    """Targets given as explicit cell ids."""

    names = [name for name, _ in named]
    frame = pd.DataFrame(
        {
            "name": names,
            "cell": [cell for _, cell in named],
            "layer": layers if layers is not None else [0] * len(named),
        }
    )
    values = pd.DataFrame({"per": [0], **{name: [1.0] for name in names}})
    return HeadTargets(locations=frame, values=values, time_column="per")


@pytest.fixture
def grid():
    return _grid()


@pytest.mark.parametrize(
    "point, cell",
    [
        (Point(37.5, 37.5), 5),        # inside cell 5
        (Point(25.0, 37.5), 4),        # on the face between 4 and 5
        (Point(0.0, 37.5), 4),         # on the grid's outer edge of cell 4
        (Point(50.0, 50.0), 5),        # on the vertex shared by 5, 6, 9 and 10
        (Point(25.0, 25.0), 0),        # on the vertex shared by 0, 1, 4 and 5
    ],
)
def test_a_point_resolves_to_the_one_lowest_numbered_cell(grid, point, cell):
    matched = _targets([("W1", point)]).match_to_model(_model(grid))

    assert matched["cell"].tolist() == [cell]


def test_a_vertex_target_is_one_row_not_four(grid):
    """The defect: this returned five rows for two targets, and PEST then read
    cell 5 at build time and cell 10 in the forward run."""

    matched = _targets([("W1", Point(50.0, 50.0)), ("W2", Point(12.5, 12.5))]).match_to_model(
        _model(grid)
    )

    assert len(matched) == 2
    assert matched.set_index("name")["cell"].to_dict() == {"W1": 5, "W2": 0}


def test_the_lowest_cell_wins_whatever_order_the_grid_is_numbered_in():
    """The spatial index answers in tree order, so the tie must be decided by cell
    NUMBER. On a grid numbered back to front, the vertex cells are 10, 11, 14, 15."""

    scrambled = _grid(order=list(reversed(range(16))))

    matched = _targets([("W1", Point(25.0, 25.0))]).match_to_model(_model(scrambled))

    assert matched["cell"].tolist() == [10]


def test_activity_is_checked_in_the_targets_own_layer(grid):
    """The same vertex resolves differently per layer, so a domain cannot be read
    at layer 0 for every target."""

    idomain = np.ones((2, 16), dtype=int)
    idomain[0, 0] = 0                      # cell 0 inactive in layer 0 only
    model = _model(grid, idomain=idomain)
    vertex = Point(25.0, 25.0)             # touches cells 0, 1, 4 and 5

    assert _targets([("W1", vertex)], layers=[0]).match_to_model(model)["cell"].tolist() == [1]
    assert _targets([("W1", vertex)], layers=[1]).match_to_model(model)["cell"].tolist() == [0]


def test_a_passthrough_cell_is_not_an_active_cell(grid):
    """idomain -1 is a vertical passthrough: MF6 writes no head there."""

    idomain = np.ones((1, 16), dtype=int)
    idomain[0, 0] = -1
    model = _model(grid, idomain=idomain)

    assert _targets([("W1", Point(25.0, 25.0))]).match_to_model(model)["cell"].tolist() == [1]


def test_a_target_layer_beyond_the_domain_raises(grid):
    model = _model(grid, idomain=np.ones((1, 16), dtype=int))

    with pytest.raises(ValueError, match="outside the model domain"):
        _targets([("W1", Point(37.5, 37.5))], layers=[3]).match_to_model(model)


def test_a_model_with_no_domain_treats_every_cell_as_active(grid):
    """A model-like object with no idomain must not be asked for one: the property
    that builds a domain raises when there is nothing to build from."""

    assert _targets([("W1", Point(37.5, 37.5))]).match_to_model(_model(grid))["cell"].tolist() == [5]


def test_targets_are_reprojected_to_the_grid(grid):
    """Two wells whose cells are 5 and 10 both came back as cell 4, located by
    their raw degrees on a State Plane grid."""

    lonlat = gpd.GeoSeries([Point(37.5, 37.5), Point(62.5, 62.5)], crs=grid.crs).to_crs("EPSG:4326")
    targets = _targets([("W1", lonlat.iloc[0]), ("W2", lonlat.iloc[1])], crs="EPSG:4326")

    matched = targets.match_to_model(_model(grid))

    assert matched.set_index("name")["cell"].to_dict() == {"W1": 5, "W2": 10}
    assert matched["x"].iloc[0] == pytest.approx(37.5, abs=1e-6)


def test_a_polygon_target_is_refused(grid):
    """A head is measured at a point; a polygon brought the duplicates back."""

    with pytest.raises(ValueError, match="measured at a"):
        _targets([("W1", box(25.0, 25.0, 50.0, 50.0))]).match_to_model(_model(grid))


def test_every_target_missing_raises_naming_both_extents(grid):
    """They used to be snapped to the nearest cell at any distance: 1.45e6 ft for
    State Plane coordinates on a grid at the origin."""

    targets = _targets([("W1", Point(1_300_000.0, 650_000.0))])

    with pytest.raises(ValueError, match="No head target lands on an active cell"):
        targets.match_to_model(_model(grid))
    with pytest.raises(ValueError, match="The targets span"):
        targets.match_to_model(_model(grid))


def test_some_targets_missing_warns_by_name_and_keeps_the_rest(grid, caplog):
    """(150, 50) is equidistant from cells 7 and 11; the nearest-cell snap raised
    a pandas error there and lost the in-grid targets with it."""

    targets = _targets(
        [("W1", Point(37.5, 37.5)), ("W2", Point(62.5, 62.5)), ("OFF", Point(150.0, 50.0))]
    )

    with caplog.at_level(logging.WARNING, logger=LOGGER):
        matched = targets.match_to_model(_model(grid))

    assert matched["name"].tolist() == ["W1", "W2"]
    assert any("'OFF'" in record.getMessage() for record in caplog.records)


def test_the_same_missed_target_is_warned_about_once(grid, caplog):
    """Every consumer calls match_to_model again; five calls said it five times."""

    targets = _targets([("W1", Point(37.5, 37.5)), ("OFF", Point(150.0, 50.0))])
    model = _model(grid)

    with caplog.at_level(logging.WARNING, logger=LOGGER):
        targets.match_to_model(model)
        targets.match_to_model(model)

    assert len([r for r in caplog.records if "'OFF'" in r.getMessage()]) == 1


def test_an_inactive_target_alone_raises_naming_the_layer(grid):
    idomain = np.ones((1, 16), dtype=int)
    idomain[0, 5] = 0
    targets = _targets([("W1", Point(37.5, 37.5))])

    with pytest.raises(ValueError, match="on no active cell in layer 0"):
        targets.match_to_model(_model(grid, idomain=idomain))


def test_a_cell_id_the_grid_does_not_have_raises(grid):
    """It reached the value table as NaN and dropped out of stats silently."""

    with pytest.raises(ValueError, match="does not have"):
        _cell_targets([("W1", 99)]).match_to_model(_model(grid))


def test_a_cell_id_inactive_in_its_layer_raises(grid):
    idomain = np.ones((2, 16), dtype=int)
    idomain[1, 5] = 0

    with pytest.raises(ValueError, match="not active in their own layer"):
        _cell_targets([("W1", 5)], layers=[1]).match_to_model(_model(grid, idomain=idomain))


def test_a_valid_cell_id_still_gets_its_centroid(grid):
    matched = _cell_targets([("W1", 5)]).match_to_model(_model(grid))

    assert matched["cell"].tolist() == [5]
    assert matched["x"].iloc[0] == pytest.approx(37.5)


def test_concentration_targets_inherit_the_rule(grid):
    matched = _targets([("W1", Point(25.0, 37.5))], cls=ConcTargets).match_to_model(_model(grid))

    assert matched["cell"].tolist() == [4]


def test_the_rule_is_the_shared_one_not_a_copy():
    assert heads_module.cell_shares is grid_geometry.cell_shares
    assert heads_module.domain_active is grid_geometry.domain_active

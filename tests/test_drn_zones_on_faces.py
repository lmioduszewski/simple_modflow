"""A DRN seepage zone reaches the cells it COVERS, not every cell it touches.

``DrnFlowTargets.zone_definitions`` resolved a zone's cells with a spatial join
on ``intersects``, which is true of a shared edge and of a shared corner. On a
4x4 grid of 25 ft cells a zone drawn exactly on cell 5 claimed all nine cells
around it, a point on a face claimed both cells and one on a vertex four, a zone
off the grid came back empty and its flow as a silent zero -- and two zones
abutting on the column of faces at x=50 BOTH claimed that column, so every drain
in it was summed into two different observations and calibrated twice.
"""

from __future__ import annotations

import logging
import os
from types import SimpleNamespace

os.environ.setdefault("MPLBACKEND", "Agg")

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import LineString, Point, box

from myflopy.modflow.mf6.grid import geometry as grid_geometry
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.mf6.observations import _zones as zones_module
from myflopy.modflow.mf6.observations.drn import DrnFlowTargets

LOGGER = "myflopy.modflow.mf6.observations._zones"
OFF_GRID = box(1_000.0, 1_000.0, 1_050.0, 1_050.0)


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


def _model(grid, idomain=None, budget=None):
    """A model-like object: the grid, a flopy-shaped DISV idomain, an optional DRN budget.

    ``node_to_lni`` is MF6's USER numbering, ``layer * ncpl + cell`` -- the node a
    cell-budget record carries, which ``tests/test_drn_zone_layers.py`` pins on a
    real run with an inactive cell.
    """

    disv = None if idomain is None else SimpleNamespace(
        idomain=SimpleNamespace(array=np.asarray(idomain, dtype=int))
    )
    ncpl = len(grid.gdf_vorPolys)
    nlay = 1 if idomain is None or np.ndim(idomain) == 1 else np.shape(idomain)[0]
    model = SimpleNamespace(
        vor=grid,
        gwf=SimpleNamespace(disv=disv, drn="drn_owner"),
        node_to_lni={layer * ncpl + cell: (layer, cell)
                     for layer in range(nlay) for cell in range(ncpl)},
    )
    if budget is not None:
        attach_budget(model, budget)
    return model


def attach_budget(model, budget, *, package_type="drn", names=("drn",)):
    """Give a model double the cell-budget file a zone set reads, and its packages.

    ``budget`` holds ``kstpkper``, ``node`` (MF6's ZERO-based user node) and ``q``,
    as columns or index levels, and optionally ``pname``, the package that wrote
    each row (default: the first of ``names``). Every package of ``names`` writes
    one record per saved time, empty or not, as MF6 does for a package declared with
    ``save_flows``; ``names`` beyond those rows still write empty records. The
    zone reader matches records by ``paknam2`` (ledger 534), so the
    double answers ``model.package(<type>)`` with packages of those names.
    """

    frame = budget.reset_index() if "node" not in budget.columns else budget.copy()
    if "pname" not in frame.columns:
        frame["pname"] = names[0]
    times = list(dict.fromkeys(tuple(item) for item in frame["kstpkper"]))
    headers, records = [], []
    for kstpkper in times:
        for name in names:
            keep = [tuple(item) == kstpkper and str(owner).upper() == name.upper()
                    for item, owner in zip(frame["kstpkper"], frame["pname"], strict=True)]
            rows = frame.loc[keep]
            headers.append((kstpkper[0] + 1, kstpkper[1] + 1, package_type.upper().rjust(16),
                            name.upper().ljust(16)))
            records.append(np.array(
                [(int(node) + 1, 0, float(q)) for node, q in zip(rows["node"], rows["q"], strict=True)],
                dtype=[("node", "<i4"), ("node2", "<i4"), ("q", "<f8")]))
    reader = SimpleNamespace(
        recordarray=np.array(headers, dtype=[("kstp", "<i4"), ("kper", "<i4"),
                                             ("text", "S16"), ("paknam2", "S16")]),
        get_record=lambda index: records[index],
        get_kstpkper=lambda: times,
    )
    model._get_budget_reader = lambda: reader
    if not hasattr(model, "name"):
        model.name = "m"
    packages = [SimpleNamespace(package_name=name) for name in names]
    previous = getattr(model, "package", None)

    def package(kind):
        if str(kind).lower() == package_type:
            return packages[0] if len(packages) == 1 else packages
        if previous is None:
            raise AttributeError(f"Package {kind!r} is not attached to model 'fake'.")
        return previous(kind)

    model.package = package
    return model


def _zones(named, *, crs="EPSG:2927", layers=None, values=None, time_column="per"):
    """Zone definitions given as geometry, one row per (name, geometry) pair."""

    frame = gpd.GeoDataFrame(
        {
            "name": [name for name, _ in named],
            "layer": layers if layers is not None else [0] * len(named),
        },
        geometry=[geometry for _, geometry in named],
        crs=crs,
    )
    return DrnFlowTargets(
        locations=frame, values=values, time_column=time_column, value_column="flow"
    )


@pytest.fixture
def grid():
    return _grid()


@pytest.mark.parametrize(
    "label, geometry, cells",
    [
        # Was [0, 1, 2, 4, 5, 6, 8, 9, 10] -- the cell and all eight neighbours.
        ("polygon equal to cell 5", box(25.0, 25.0, 50.0, 50.0), [5]),
        ("polygon over cells 5 and 6", box(25.0, 25.0, 75.0, 50.0), [5, 6]),
        ("point on the 4|5 face", Point(25.0, 37.5), [4]),          # was [4, 5]
        ("point on the 5/6/9/10 vertex", Point(50.0, 50.0), [5]),   # was [5, 6, 9, 10]
        ("point inside cell 5", Point(37.5, 37.5), [5]),
        ("a drain trace through 4 and 5", LineString([(10.0, 37.5), (40.0, 37.5)]), [4, 5]),
    ],
)
def test_a_zone_reaches_the_cells_it_covers(grid, label, geometry, cells):
    zones = _zones([("seep", geometry)]).zone_definitions(_model(grid))

    assert zones["cells"].tolist() == [cells], label


def test_two_zones_abutting_on_a_face_no_longer_share_a_column(grid):
    """The defect that mattered: every drain in the column at x=50 was summed into
    BOTH observations, so the same water was calibrated against twice."""

    zones = _zones(
        [("west", box(0.0, 0.0, 50.0, 100.0)), ("east", box(50.0, 0.0, 100.0, 100.0))]
    ).zone_definitions(_model(grid))

    west, east = zones.set_index("name")["cells"].loc["west"], zones.set_index("name")["cells"].loc["east"]
    assert west == [0, 1, 4, 5, 8, 9, 12, 13]
    assert east == [2, 3, 6, 7, 10, 11, 14, 15]
    assert set(west) & set(east) == set()


def test_the_lowest_cell_wins_whatever_order_the_grid_is_numbered_in():
    """The spatial index answers in tree order, so a point on a vertex must be
    decided by cell NUMBER. Numbered back to front, those cells are 5, 6, 9, 10."""

    scrambled = _grid(order=list(reversed(range(16))))

    zones = _zones([("seep", Point(50.0, 50.0))]).zone_definitions(_model(scrambled))

    assert zones["cells"].tolist() == [[5]]


def test_only_cells_active_in_the_zones_own_layer_are_kept(grid):
    """MF6 writes no DRN budget entry and no DRN observation in a dead cell."""

    idomain = np.ones((2, 16), dtype=int)
    idomain[0, 5] = 0       # cell 5 inactive in layer 0 only
    idomain[1, 6] = -1      # a vertical passthrough is not active either
    model = _model(grid, idomain=idomain)
    patch = box(25.0, 25.0, 75.0, 50.0)     # covers cells 5 and 6

    assert _zones([("seep", patch)], layers=[0]).zone_definitions(model)["cells"].tolist() == [[6]]
    assert _zones([("seep", patch)], layers=[1]).zone_definitions(model)["cells"].tolist() == [[5]]


def test_a_zone_in_a_layer_the_model_does_not_have_raises(grid):
    """New: the layer was ignored outright before, so a zone declared in layer 3 of
    a two-layer model resolved against layer 0 and read as fine."""

    model = _model(grid, idomain=np.ones((2, 16), dtype=int))

    with pytest.raises(ValueError, match="Layer 3 is outside the model domain"):
        _zones([("seep", box(25.0, 25.0, 50.0, 50.0))], layers=[3]).zone_definitions(model)


def test_a_zone_landing_on_nothing_raises_naming_both_extents(grid):
    """It came back as ``[]`` and its simulated flow as a silent zero."""

    zones = _zones([("seep", OFF_GRID)])

    with pytest.raises(ValueError) as raised:
        zones.zone_definitions(_model(grid))

    message = str(raised.value)
    assert "No DRN zone lands on an active cell" in message
    assert "nothing to sum" in message                  # not "compare against"
    assert "The zones span (1,000.0, 1,000.0, 1,050.0, 1,050.0)" in message
    assert "the grid (0.0, 0.0, 100.0, 100.0)" in message


def test_the_raise_reports_the_zone_files_own_bounds_not_the_reprojected_ones(grid):
    """A location mismatch is TWO boxes: reprojecting the zones into the grid's
    feet first would print two boxes that agree and hide the mismatch."""

    zones = _zones([("seep", box(-122.3, 47.6, -122.2, 47.7))], crs="EPSG:4326")

    with pytest.raises(ValueError) as raised:
        zones.zone_definitions(_model(grid))

    message = str(raised.value)
    assert "The zones span (-122.3, 47.6, -122.2, 47.7)" in message
    assert "the grid (0.0, 0.0, 100.0, 100.0)" in message
    assert "reprojected from EPSG:4326 to the grid's EPSG:2927" in message


def test_the_miss_reason_names_the_zones_own_layer(grid):
    """A zone over live cells but in a dead LAYER is a different miss from one off
    the grid, and the message has to say which layer was looked in."""

    idomain = np.ones((2, 16), dtype=int)
    idomain[1, :] = 0
    model = _model(grid, idomain=idomain)

    with pytest.raises(ValueError, match="on no active cell in layer 1"):
        _zones([("seep", box(25.0, 25.0, 50.0, 50.0))], layers=[1]).zone_definitions(model)


def test_one_zone_of_several_missing_warns_by_name_and_keeps_the_rest(grid, caplog):
    zones = _zones([("seep", box(25.0, 25.0, 50.0, 50.0)), ("off", OFF_GRID)])

    with caplog.at_level(logging.WARNING, logger=LOGGER):
        resolved = zones.zone_definitions(_model(grid))

    assert resolved.set_index("name")["cells"].to_dict() == {"seep": [5], "off": []}
    assert any("'off'" in record.getMessage() for record in caplog.records)


def test_a_long_miss_listing_is_truncated_and_counted(grid, caplog):
    """Twelve of thirteen zones missing: ten named, the rest counted."""

    off = [(f"off{i}", box(1_000.0 + 10 * i, 1_000.0, 1_005.0 + 10 * i, 1_005.0)) for i in range(12)]
    zones = _zones([("seep", box(25.0, 25.0, 50.0, 50.0)), *off])

    with caplog.at_level(logging.WARNING, logger=LOGGER):
        zones.zone_definitions(_model(grid))

    message = next(r.getMessage() for r in caplog.records if "DRN zones" in r.getMessage())
    assert "12 of 13 DRN zones land on no active cell" in message
    assert "'off0'" in message and "'off9'" in message
    assert "'off10'" not in message and message.endswith("and 2 more")


def test_the_same_missed_zone_is_warned_about_once(grid, caplog):
    """Seven call sites in drn.py resolve the zones again, and PEST an eighth; the
    warning must not be said once per call."""

    zones = _zones([("seep", box(25.0, 25.0, 50.0, 50.0)), ("off", OFF_GRID)])
    model = _model(grid)

    with caplog.at_level(logging.WARNING, logger=LOGGER):
        zones.zone_definitions(model)
        zones.zone_definitions(model)

    assert len([r for r in caplog.records if "'off'" in r.getMessage()]) == 1


def test_a_zone_bordering_cells_without_covering_one_says_so(grid):
    """The area rule creates this miss (ledger 210), so it gets its own wording:
    a zone drawn just outside the grid shares the outer faces of cells 7 and 11
    and covers neither."""

    outside = box(100.0, 25.0, 150.0, 75.0)

    with pytest.raises(ValueError, match="only borders cells"):
        _zones([("seep", outside)]).zone_definitions(_model(grid))


def test_a_zone_given_as_several_rows_is_one_zone(grid):
    """A multipart zone arrives exploded, one row per part. The grouped join
    unioned them and so must this: both rows carry the whole zone."""

    zones = _zones(
        [("seep", box(25.0, 25.0, 50.0, 50.0)), ("seep", box(50.0, 50.0, 75.0, 75.0))]
    ).zone_definitions(_model(grid))

    assert zones["cells"].tolist() == [[5, 10], [5, 10]]
    assert zones["cells"].iloc[0] is not zones["cells"].iloc[1]     # no shared list


def test_a_multipart_zone_that_partly_lands_is_not_reported_missing(grid, caplog):
    """Counted by NAME. Accumulating the misses per ROW warned "1 of 2 DRN zones
    land on no active cell ... 'seep'" of a zone that had landed and was still in
    the frame -- and the off-grid part is first here, which is where a
    first-row-wins union loses the zone outright."""

    zones = _zones([("seep", OFF_GRID), ("seep", box(25.0, 25.0, 50.0, 50.0))])

    with caplog.at_level(logging.WARNING, logger=LOGGER):
        resolved = zones.zone_definitions(_model(grid))

    assert resolved["cells"].tolist() == [[5], [5]]
    assert [r.getMessage() for r in caplog.records if "DRN zone" in r.getMessage()] == []


def test_a_multipart_zone_missing_every_part_still_raises(grid):
    """Two rows, one name, neither landing: that is ONE zone of one missing, so it
    is a total miss and raises rather than warning about a partial one."""

    zones = _zones([("seep", OFF_GRID), ("seep", box(2_000.0, 2_000.0, 2_050.0, 2_050.0))])

    with pytest.raises(ValueError, match="No DRN zone lands on an active cell"):
        zones.zone_definitions(_model(grid))


def test_zones_are_reprojected_to_the_grid(grid):
    """A zone file in degrees over a State Plane grid was matched by raw degrees."""

    lonlat = gpd.GeoSeries([box(25.0, 25.0, 50.0, 50.0)], crs=grid.crs).to_crs("EPSG:4326")

    zones = _zones([("seep", lonlat.iloc[0])], crs="EPSG:4326").zone_definitions(_model(grid))

    assert zones["cells"].tolist() == [[5]]
    # The frame handed back is the located one -- PEST writes it as the zones
    # snapshot, so it must be in the CRS its cell ids were resolved in.
    assert zones.crs == grid.crs
    assert zones.geometry.iloc[0].bounds[0] == pytest.approx(25.0, abs=1e-3)


def test_explicit_cell_ids_are_still_taken_as_given(grid):
    """A cell list is an assertion about the model, not geometry to resolve."""

    targets = DrnFlowTargets(locations={"seep": [0, 15]})

    assert targets.zone_definitions(_model(grid))["cells"].tolist() == [[0, 15]]


def test_zone_membership_is_geometry_not_the_drn_package(grid):
    """Recorded, not changed (ledger 238): a covered cell with no drain stays in
    the zone and adds zero to the sum."""

    budget = pd.DataFrame(
        {"kstpkper": [(0, 0), (0, 0)], "node": [5, 9], "q": [-2.0, -3.0]}
    ).set_index(["kstpkper", "node"])
    model = _model(grid, budget=budget)
    zones = _zones([("seep", box(25.0, 25.0, 50.0, 75.0))])     # cells 5 and 9

    assert zones.zone_definitions(model)["cells"].tolist() == [[5, 9]]
    assert zones.simulated_series(model)["seep"].tolist() == [-5.0]


def test_a_drn_flow_comparison_still_runs_end_to_end(grid):
    """Two zones, a real DRN budget, residuals: the columns a calibration reads."""

    budget = pd.DataFrame(
        {
            "kstpkper": [(0, 0), (0, 0), (0, 1), (0, 1)],
            "node": [5, 6, 5, 6],
            "q": [-2.0, -1.0, -2.5, -1.5],
        }
    ).set_index(["kstpkper", "node"])
    zones = _zones(
        [("west", box(25.0, 25.0, 50.0, 50.0)), ("east", box(50.0, 25.0, 75.0, 50.0))],
        values={"per": [0, 1], "west": [-2.2, -2.4], "east": [-0.9, -1.4]},
    )

    compare = zones.compare(_model(grid, budget=budget))

    west_0 = compare.loc[(compare["name"] == "west") & (compare["per"] == 0)]
    assert west_0["sim_flow"].iloc[0] == pytest.approx(-2.0)
    assert west_0["residual"].iloc[0] == pytest.approx(0.2)
    assert compare.loc[compare["name"] == "east", "sim_flow"].tolist() == pytest.approx([-1.0, -1.5])
    assert zones.stats(_model(grid, budget=budget))["n"].iloc[0] == 4


def test_the_flopy_observation_records_follow_the_resolved_cells(grid):
    records = _zones([("seep", box(25.0, 25.0, 50.0, 50.0))]).to_flopy_obs(_model(grid))

    assert list(records.values())[0] == [("seep_c5", "DRN", (0, 5))]


def test_the_rule_is_the_shared_one_not_a_copy():
    assert zones_module.cell_shares is grid_geometry.cell_shares
    assert zones_module.domain_active is grid_geometry.domain_active
    assert zones_module.miss_reason is grid_geometry.miss_reason

"""Derived per-cell PRT maps: release groups, normalization, and the view verbs.

The normalization tests are pure-frame units (no MF6); the view tests run one
tiny two-cell GWF+PRT simulation once per module so the maps are exercised
against real track records rather than a hand-written CSV whose ``icell``
basing could drift from MF6's.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus  # noqa: E402
from myflopy.modflow.mf6.prt import (  # noqa: E402
    PRTProject,
    PRTReleasePoints,
    open_prt_run,
)
from myflopy.modflow.mf6.prt_maps import (  # noqa: E402
    PRT_COLORSCALE,
    PRTCaptureView,
    PRTEndpointsView,
    PRTPathlineView,
    PRTTravelTimeView,
    _join_groups,
    particle_endpoint_table,
    pathline_cell_table,
)
from myflopy.modflow.mf6.simulation.base import SimulationBase  # noqa: E402
from myflopy.modflow.mf6.simulation.discretization import (  # noqa: E402
    DisvGrid,
    TemporalDiscretization,
)
from myflopy.modflow.mf6.simulation.packages import (  # noqa: E402
    CHD,
    InitialConditions,
    KFlow,
    OutputControl,
    Storage,
)
from myflopy.modflow.utils.datatypes.hover import NBSP  # noqa: E402
from myflopy.viz import Fig  # noqa: E402


def _flow_model(name: str, workspace: Path) -> SimulationBase:
    """A two-cell Voronoi flow model with a head gradient left to right."""

    verts = np.array(
        [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [2.0, 0.0], [2.0, 1.0]],
        dtype=float,
    )
    vor = VoronoiGridPlus(
        verts=verts,
        iverts=[[0, 3, 2, 1], [1, 2, 5, 4]],
        xcyc=np.array([[0.5, 0.5], [1.5, 0.5]], dtype=float),
    )
    model = SimulationBase(name=name, mf_folder_path=workspace, vor=vor, nper=1)
    DisvGrid(vor=vor, model=model, top=[10.0, 10.0], bottom=[[0.0, 0.0]], nlay=1, idomain=[[1, 1]])
    TemporalDiscretization(model=model, per_len=1.0, num_steps=1, multiplier=1.0)
    InitialConditions(model=model, vor=vor, nlay=1, strt=[10.0, 9.0])
    KFlow(model=model, k=[10.0, 10.0], k33_vert=[1.0, 1.0], save_specific_discharge=True)
    Storage(model=model, sto_steady={0: True})
    CHD(model=model, stress_period_data={0: [((0, 0), 10.0), ((0, 1), 9.0)]})
    OutputControl(model=model)
    return model


@pytest.fixture(scope="module")
def prt_run(tmp_path_factory):
    """One completed GWF+PRT run whose two particles carry release-group boundnames."""

    workspace = tmp_path_factory.mktemp("prt_maps")
    model = _flow_model("prt_maps_flow", workspace)
    assert model.run_simulation()

    release = PRTReleasePoints.merge(
        PRTReleasePoints.from_cells(model, [0], group="west_wells"),
        PRTReleasePoints.from_cells(model, [1], group="east_wells"),
    )
    project = PRTProject(
        model, workspace=workspace / "prt", release_points=release, porosity=0.25
    )
    return project.run(silent=True)


def _track_frame(**overrides) -> pd.DataFrame:
    """A minimal synthetic track frame (two particles, one terminating record each)."""

    data = {
        "imdl": [1, 1, 1],
        "iprp": [1, 1, 1],
        "irpt": [1, 1, 2],
        "ilay": [1, 1, 2],
        "icell": [1, 2, 3],
        "ireason": [0, 3, 3],
        "trelease": [0.0, 0.0, 5.0],
        "t": [0.0, 12.0, 30.0],
        "name": ["WEST", "WEST", "EAST"],
    }
    data.update(overrides)
    return pd.DataFrame(data)


# -- normalization ---------------------------------------------------------


def test_pathline_cell_table_converts_one_based_icell_to_cell_and_layer():
    """``icell`` is a ONE-based whole-grid node number, not a zero-based cell.

    Everything downstream (maps, hover, capture zones) reads ``cell``/``layer``,
    so this conversion is the load-bearing step: on a 2-cell-per-layer grid
    ``icell=3`` is layer 1's cell 0, NOT cell 3.
    """

    frame = pathline_cell_table(_track_frame(), ncpl=2)

    assert frame["cell"].tolist() == [0, 1, 0]
    assert frame["layer"].tolist() == [0, 0, 1]
    # layer agrees with MF6's own one-based ilay
    assert frame["layer"].tolist() == (frame["ilay"] - 1).tolist()
    # travel time is measured from each particle's own release
    assert frame["travel_time"].tolist() == [0.0, 12.0, 25.0]
    assert frame["release_group"].tolist() == ["WEST", "WEST", "EAST"]
    # one stable label per particle, shared by that particle's records
    assert frame["particle"].tolist() == ["1-1-1-0", "1-1-1-0", "1-1-2-5"]


def test_pathline_cell_table_handles_missing_columns_and_empty_runs():
    with pytest.raises(KeyError, match="icell"):
        pathline_cell_table(pd.DataFrame({"t": [1.0], "trelease": [0.0]}), ncpl=2)

    empty = pathline_cell_table(_track_frame().iloc[0:0], ncpl=2)
    assert empty.empty
    assert {"cell", "layer", "travel_time", "release_group"} <= set(empty.columns)


def test_particle_endpoint_table_keeps_one_terminating_row_per_particle():
    endpoints = particle_endpoint_table(_track_frame(), ncpl=2)

    assert len(endpoints) == 2
    assert endpoints["particle"].tolist() == ["1-1-1-0", "1-1-2-5"]
    assert endpoints["travel_time"].tolist() == [12.0, 25.0]
    assert endpoints["cell"].tolist() == [1, 0]


def test_particle_endpoint_table_falls_back_to_the_last_record_in_transit():
    """Particles still moving when tracking stopped are placed where they got to."""

    in_transit = _track_frame(ireason=[0, 1, 1])
    endpoints = particle_endpoint_table(in_transit, ncpl=2)

    assert len(endpoints) == 2
    assert endpoints["travel_time"].tolist() == [12.0, 25.0]


# -- release groups --------------------------------------------------------


def test_release_points_carry_group_labels_and_merge_renumbers(tmp_path):
    model = _flow_model("release_groups", tmp_path)

    one_label = PRTReleasePoints.from_cells(model, [0, 1], group="west_wells")
    per_point = PRTReleasePoints.from_cells(model, [0, 1], group=["west", "east"])
    ungrouped = PRTReleasePoints.from_cells(model, [0, 1])

    assert one_label.packagedata[0] == (0, (0, 0), 0.5, 0.5, 0.5, "west_wells")
    assert one_label.groups == ("west_wells",)
    assert per_point.groups == ("west", "east")
    assert one_label.has_groups and per_point.has_groups
    assert not ungrouped.has_groups
    assert ungrouped.groups == ()
    assert len(ungrouped.packagedata[0]) == 5

    # from_points takes the same labels
    from_points = PRTReleasePoints.from_points(
        model, [(0.5, 0.5), (1.5, 0.5)], group="wells"
    )
    assert from_points.groups == ("wells",)

    with pytest.raises(ValueError, match="one label or one per release point"):
        PRTReleasePoints.from_cells(model, [0, 1], group=["only_one"])

    # merge renumbers irpt across sets so grouped sets become one PRP
    merged = PRTReleasePoints.merge(per_point, from_points)
    assert [row[0] for row in merged.packagedata] == [0, 1, 2, 3]
    assert merged.groups == ("west", "east", "wells")


def test_prt_project_enables_boundnames_only_for_grouped_release_points(tmp_path):
    """PRP rejects a sixth field without BOUNDNAMES, so the flag follows the rows."""

    model = _flow_model("boundnames_flow", tmp_path)
    grouped = PRTProject(
        model,
        workspace=tmp_path / "grouped",
        release_points=PRTReleasePoints.from_cells(model, [0], group="west_wells"),
    )
    plain = PRTProject(
        model,
        workspace=tmp_path / "plain",
        release_points=PRTReleasePoints.from_cells(model, [0]),
    )

    assert grouped.prp.boundnames.get_data() is True
    assert not plain.prp.boundnames.get_data()
    grouped.write()
    assert "BOUNDNAMES" in (tmp_path / "grouped" / f"{grouped.name}.prp").read_text()


# -- the view nouns --------------------------------------------------------


def test_travel_time_view_maps_time_of_travel_per_cell(prt_run):
    view = prt_run.travel_time
    assert isinstance(view, PRTTravelTimeView)

    frame = view.get()
    assert list(frame.columns) == [
        "layer",
        "cell",
        "travel_time",
        "particle_count",
        "min_time",
        "max_time",
        "release_groups",
    ]
    # both particles stop in the downgradient cell, so one row summarizes both
    assert frame["cell"].tolist() == [1]
    assert frame["particle_count"].tolist() == [2]
    assert frame["travel_time"].iloc[0] == pytest.approx(
        (frame["min_time"].iloc[0] + frame["max_time"].iloc[0]) / 2
    )
    assert view.get(stat="max")["travel_time"].iloc[0] == pytest.approx(
        frame["max_time"].iloc[0]
    )
    assert view.summary()["label"].iloc[0] == "prt.travel_time"

    choro = view.map()
    assert choro.colorscale == PRT_COLORSCALE == "earth"
    # cells no particle reached stay blank rather than reading as a zero travel time
    assert np.isnan(choro.zs[0])
    assert choro.zs[1] == pytest.approx(frame["travel_time"].iloc[0])

    template = choro.get_choropleth().hovertemplate
    assert "Travel time (median)" in template
    assert "particles" in template
    # the map integrates the whole run -- no stress period to footer
    assert "Period" not in template


def test_endpoints_view_counts_terminations_per_cell(prt_run):
    view = prt_run.endpoints
    assert isinstance(view, PRTEndpointsView)

    frame = view.get()
    assert frame["cell"].tolist() == [1]
    assert frame["particle_count"].tolist() == [2]
    assert sorted(frame["release_groups"].iloc[0].split(", ")) == [
        "EAST_WELLS",
        "WEST_WELLS",
    ]
    assert view.summary()["label"].iloc[0] == "prt.endpoints"

    choro = view.map()
    assert choro.colorscale == "earth"
    assert choro.zs[1] == pytest.approx(2.0)
    assert "Endpoints" in choro.get_choropleth().hovertemplate


def test_capture_view_splits_endpoints_by_release_group(prt_run):
    view = prt_run.capture
    assert isinstance(view, PRTCaptureView)

    # MF6 echoes PRP boundnames into the track CSV uppercased
    assert view.groups() == ["EAST_WELLS", "WEST_WELLS"]
    frame = view.get()
    assert frame["group"].tolist() == ["EAST_WELLS", "WEST_WELLS"]
    assert frame["particle_count"].tolist() == [1, 1]

    mosaic = view.map()
    assert isinstance(mosaic, Fig)
    assert [note.text for note in mosaic.layout.annotations] == [
        "EAST_WELLS",
        "WEST_WELLS",
    ]
    # one shared color scale across the panels is what makes them comparable
    assert len(mosaic.data) == 2
    assert all(trace.coloraxis == "coloraxis" for trace in mosaic.data)

    single = view.map(group="WEST_WELLS")
    assert single.zs[1] == pytest.approx(1.0)
    assert "Capture: WEST_WELLS" in single.get_choropleth().hovertemplate

    with pytest.raises(KeyError, match="NORTH_WELLS"):
        view.map(group="NORTH_WELLS")
    with pytest.raises(ValueError, match="by must be one of"):
        view.get(by="whatever")


def test_capture_view_without_boundnames_names_the_fix(prt_run, tmp_path):
    """An ungrouped run points at ``group=`` and the release-point fallback."""

    workspace = tmp_path / "ungrouped"
    workspace.mkdir()
    _track_frame(name=["", "", ""]).to_csv(workspace / "plain.trk.csv", index=False)
    plain = open_prt_run(prt_run.flow_model, workspace)

    with pytest.raises(ValueError, match=r"group='west_wells'.*by='release_point'"):
        plain.capture.get()

    by_point = plain.capture.get(by="release_point")
    assert by_point["group"].tolist() == ["1", "2"]


def test_join_groups_dedupes_sorts_and_drops_blanks():
    assert _join_groups(["WEST", "EAST", "WEST"]) == "EAST, WEST"
    assert _join_groups(["", "WEST", ""]) == "WEST"
    assert _join_groups(["", ""]) == ""


@pytest.fixture(scope="module")
def two_group_run(prt_run, tmp_path_factory):
    """A run whose two groups terminate in DIFFERENT cells, in two layers.

    The real fixture is degenerate on purpose (two cells, both particles stop in
    the same one), which would let a broken layer filter or an identical-panel
    capture mosaic pass. This synthetic track CSV separates every axis: WEST ends
    in cell 0, EAST in cell 1, one record sits in layer 1, and cell 0's travel
    times are chosen so a pooled statistic differs from a statistic-of-statistics.
    """

    workspace = tmp_path_factory.mktemp("prt_two_group")
    pd.DataFrame(
        {
            "imdl": [1, 1, 1, 1],
            "iprp": [1, 1, 1, 1],
            "irpt": [1, 2, 3, 4],
            "ilay": [1, 1, 2, 1],
            # layer 0 cell 0, layer 0 cell 1, layer 1 cell 0, layer 0 cell 0
            "icell": [1, 2, 3, 1],
            "ireason": [3, 3, 3, 3],
            "trelease": [0.0, 0.0, 0.0, 0.0],
            "t": [10.0, 100.0, 1000.0, 20.0],
            # inside the two-cell grid, so the pathline map can reproject them
            "x": [0.5, 1.5, 0.5, 0.5],
            "y": [0.5, 0.5, 0.5, 0.5],
            "z": [5.0, 5.0, 4.0, 5.0],
            "name": ["WEST", "EAST", "EAST", "WEST"],
        }
    ).to_csv(workspace / "two_group.trk.csv", index=False)
    return open_prt_run(prt_run.flow_model, workspace)


def test_layer_selector_filters_and_pooling_recomputes_the_statistic(two_group_run):
    view = two_group_run.travel_time

    every_layer = view.get()
    assert every_layer[["layer", "cell"]].values.tolist() == [[0, 0], [0, 1], [1, 0]]
    assert every_layer["particle_count"].tolist() == [2, 1, 1]
    assert view.get(layer=1)["cell"].tolist() == [0]
    assert view.get(layer=1)["travel_time"].tolist() == [1000.0]
    assert view.get(layer=0)["travel_time"].tolist() == [15.0, 100.0]

    # layer=None pools layers per plan-view cell and RECOMPUTES the statistic over
    # the pooled particles: cell 0 holds 10 and 20 in layer 0 plus 1000 in layer 1,
    # so the pooled median is 20. A median of the per-layer medians (15, 1000)
    # would be 507.5 -- this assertion is what tells the two apart.
    assert view.map().zs[0] == pytest.approx(20.0)
    assert view.map(layer=0).zs[0] == pytest.approx(15.0)


def test_logscale_transforms_the_mapped_values_only(two_group_run):
    """``logscale`` was silently dropped on custom-zs maps before this (ledger 60)."""

    view = two_group_run.travel_time
    linear = view.map(layer=0, stat="min")
    logged = view.map(layer=0, stat="min", logscale=True)

    assert linear.zs == pytest.approx([10.0, 100.0])
    assert logged.zs == pytest.approx([1.0, 2.0])
    # the hover keeps real travel times, so a log map is still readable in units
    assert any("10" in str(value) for value in logged.get_choropleth().customdata[0])


def test_logscale_blanks_non_positive_values(prt_run):
    """A particle that terminates at release has travel time 0 -- log10(0) is not -inf here."""

    zs = prt_run.travel_time.map(stat="min", logscale=True).zs
    assert np.isnan(zs[1])  # min travel time in that cell is exactly 0.0


def test_capture_panels_differ_per_group(two_group_run):
    view = two_group_run.capture

    assert view.groups() == ["EAST", "WEST"]
    west = view.map(group="WEST", layer=0)
    east = view.map(group="EAST", layer=0)
    assert west.zs[0] == pytest.approx(2.0) and np.isnan(west.zs[1])
    assert np.isnan(east.zs[0]) and east.zs[1] == pytest.approx(1.0)

    mosaic = view.map(layer=0)
    assert [note.text for note in mosaic.layout.annotations] == ["EAST", "WEST"]
    # the panels carry genuinely different data, not the same map twice
    assert list(mosaic.data[0].z) != list(mosaic.data[1].z)


def test_maps_render_through_the_matplotlib_backend(two_group_run):
    from matplotlib.figure import Figure

    assert isinstance(two_group_run.travel_time.map(layer=0, backend="mpl"), Figure)
    assert isinstance(two_group_run.endpoints.map(layer=0, backend="mpl"), Figure)
    assert isinstance(
        two_group_run.capture.map(group="WEST", layer=0, backend="mpl"), Figure
    )
    assert isinstance(two_group_run.travel_time.plot(backend="mpl"), Figure)
    assert isinstance(two_group_run.endpoints.plot(backend="mpl"), Figure)
    assert isinstance(two_group_run.capture.plot(backend="mpl"), Figure)
    # a mosaic is a Plotly composition -- say so rather than half-render it
    with pytest.raises(ValueError, match="Plotly composition"):
        two_group_run.capture.map(backend="mpl")


def test_time_integrated_maps_reject_period_and_agg_selectors(two_group_run):
    """A swallowed per=/agg= would silently answer a different question."""

    with pytest.raises(TypeError, match="time-integrated"):
        two_group_run.travel_time.map(per=1)
    with pytest.raises(TypeError, match="stat="):
        two_group_run.travel_time.map(agg="max")
    with pytest.raises(NotImplementedError, match="no period"):
        two_group_run.endpoints.animate()


def test_travel_time_hover_reads_the_models_time_unit(prt_run):
    """Travel times are in TDIS units; the hover reads them off the model, not 'd'."""

    assert prt_run.flow_model.sim.tdis.time_units.get_data().upper() == "DAYS"
    trace = prt_run.travel_time.map().get_choropleth()

    # units are rendered into the customdata values (format_number), not the template
    rendered = [str(value) for value in trace.customdata[1]]
    assert any(value.endswith(f"{NBSP}d") for value in rendered), rendered
    # min/max are elapsed travel times, so they are not labelled as arrival clock times
    assert "fastest" in trace.hovertemplate and "slowest" in trace.hovertemplate
    assert "arrival" not in trace.hovertemplate

    # an explicit unit still wins, and an unreadable TDIS degrades to no unit
    explicit = prt_run.travel_time.map(units="yr").get_choropleth()
    assert any(str(value).endswith(f"{NBSP}yr") for value in explicit.customdata[1])


def test_derived_views_answer_the_shared_verbs(prt_run):
    """Every noun answers get/summary/plot/map; spatial composers come from SpatialView."""

    for view in (prt_run.travel_time, prt_run.endpoints, prt_run.capture):
        assert isinstance(view.plot(), Fig)
        assert not view.get().empty
        assert len(view.summary()) == 1
        # PRT has no stress-period axis: the series verb says so instead of
        # failing deep inside a groupby on a missing 'per' column
        with pytest.raises(NotImplementedError, match="no stress-period axis"):
            view._series_table()

    assert isinstance(prt_run.endpoints.mosaic(), Fig)
    assert isinstance(prt_run.travel_time.mosaic(), Fig)
    # a capture facet needs one named group -- the mosaic verb is map()'s job
    assert isinstance(prt_run.capture.mosaic(group="WEST_WELLS"), Fig)
    with pytest.raises(ValueError, match="needs one release group"):
        prt_run.capture.mosaic()


# -- pathlines: the trajectories themselves --------------------------------


def test_pathlines_is_a_view_over_the_raw_track_records(prt_run):
    """``results.pathlines`` is a view; the untouched CSV stays on ``track_records``."""

    view = prt_run.pathlines
    assert isinstance(view, PRTPathlineView)

    raw = prt_run.track_records
    records = view.get()
    # the normalized table is a superset: every raw column survives
    assert set(raw.columns) <= set(records.columns)
    assert {"cell", "layer", "travel_time", "release_group", "particle"} <= set(records.columns)
    assert len(records) == len(raw)
    # and the raw frame is still cached/refreshable for FloPy and PyVista
    assert prt_run.track_records is raw

    # selectors narrow the same table
    assert set(view.get(group="WEST_WELLS")["release_group"]) == {"WEST_WELLS"}
    assert set(view.get(layer=0)["layer"]) == {0}


def test_pathline_summary_is_one_row_per_particle(prt_run):
    digest = prt_run.pathlines.summary()

    assert len(digest) == prt_run.pathlines.get()["particle"].nunique()
    assert list(digest.columns) == [
        "particle",
        "release_group",
        "records",
        "start_cell",
        "start_layer",
        "end_cell",
        "end_layer",
        "travel_time",
        "terminated",
    ]
    assert digest["terminated"].all()
    assert (digest["records"] >= 1).all()
    assert set(digest["release_group"]) == {"WEST_WELLS", "EAST_WELLS"}


def test_pathline_map_draws_one_polyline_per_particle_over_the_base(prt_run):
    """The map is a Choro carrying one Scattermap line per particle."""

    view = prt_run.pathlines
    choro = view.map()
    overlays = choro.overlay_traces()

    assert len(overlays) == view.get()["particle"].nunique()
    assert {trace.type for trace in overlays} == {"scattermap"}
    # vertices are reprojected onto the grid's own lat/lon frame, not raw model x/y
    west, south, east, north = prt_run.flow_model.vor.gdf_latlon.total_bounds
    for trace in overlays:
        assert all(west <= lon <= east for lon in trace.lon)
        assert all(south <= lat <= north for lat in trace.lat)

    # the assembled figure is the cells plus every path
    fig = choro.fig
    assert isinstance(fig, Fig)
    assert [trace.type for trace in fig.data] == ["choroplethmap", "scattermap", "scattermap"]
    assert fig.layout.title.text == "Particle pathlines"


def test_pathline_hover_names_the_particle_and_reads_the_time_unit(prt_run):
    """Hover is per *vertex* of one particle -- not per grid cell."""

    trace = prt_run.pathlines.map().overlay_traces()[0]

    assert "Particle" in trace.hovertemplate
    assert "cell" in trace.hovertemplate and "elapsed" in trace.hovertemplate
    # time-integrated: no stress period to report
    assert "Period" not in trace.hovertemplate
    # one customdata row per vertex, with the TDIS unit baked into the value
    records = prt_run.pathlines.get()
    first = records.loc[records["particle"] == records["particle"].iloc[0]]
    assert len(trace.customdata) == len(first) == len(trace.lon)
    assert any(str(value).endswith(f"{NBSP}d") for row in trace.customdata for value in row)


def test_pathline_colors_are_the_shared_category_colors(prt_run):
    """A release group keeps one color across every PRT figure."""

    from myflopy.viz import PALETTE

    map_colors = {
        trace.name: trace.line.color for trace in prt_run.pathlines.map().overlay_traces()
    }
    assert set(map_colors) == {"WEST_WELLS", "EAST_WELLS"}
    assert set(map_colors.values()) <= set(PALETTE.categorical)
    assert len(set(map_colors.values())) == 2  # distinct groups, distinct colors

    curve_colors = {trace.name: trace.line.color for trace in prt_run.travel_time.plot().data}
    assert curve_colors == map_colors


def test_pathline_map_base_can_be_heads_nothing_or_another_map(prt_run):
    view = prt_run.pathlines

    blank = view.map(base=None)
    assert all(np.isnan(value) for value in blank.zs)  # framed, but no cell values
    assert len(blank.overlay_traces()) == 2

    # paths ride on any map you already built -- here the capture zones
    borrowed = view.map(base=prt_run.travel_time.map())
    assert len(borrowed.overlay_traces()) == 2
    assert borrowed.get_choropleth().type == "choroplethmap"

    with pytest.raises(ValueError, match="not a pathline base map"):
        view.map(base="concentration")
    with pytest.raises(TypeError, match="Cannot draw pathlines over"):
        view.map(base=42)


def test_pathline_map_caps_particles_and_says_which_are_drawn(two_group_run):
    """A capped figure samples across the run and never hides the cap."""

    view = two_group_run.pathlines
    assert view.get()["particle"].nunique() == 4

    with pytest.warns(UserWarning, match="2 of 4 particles"):
        choro = view.map(max_particles=2)

    assert len(choro.overlay_traces()) == 2
    assert "2 of 4 particles" in choro.fig.layout.title.text
    # sampled across the sorted particles, so both groups survive the cap
    assert {trace.name for trace in choro.overlay_traces()} == {"WEST", "EAST"}

    uncapped = view.map(max_particles=None)
    assert len(uncapped.overlay_traces()) == 4


def test_pathline_mosaic_gives_one_panel_per_release_group(two_group_run):
    fig = two_group_run.pathlines.mosaic()

    # each panel keeps its cells AND both of its paths -- overlays are not dropped
    assert [trace.type for trace in fig.data] == [
        "choroplethmap",
        "scattermap",
        "scattermap",
        "choroplethmap",
        "scattermap",
        "scattermap",
    ]
    assert [annotation.text for annotation in fig.layout.annotations] == ["EAST", "WEST"]

    with pytest.raises(ValueError, match="only by='release_group'"):
        two_group_run.pathlines.mosaic(by="layer")


@pytest.fixture(scope="module")
def ungrouped_run(tmp_path_factory):
    """A REAL run whose release points carry no boundnames at all.

    Confirms against MF6 what the ungrouped path actually sees: the track CSV
    still has a ``name`` column, but every value is blank -- so there are no
    groups to facet by and coloring falls back to per-particle.
    """

    workspace = tmp_path_factory.mktemp("prt_ungrouped")
    model = _flow_model("ungrouped_flow", workspace)
    assert model.run_simulation()

    project = PRTProject(
        model,
        workspace=workspace / "prt",
        release_points=PRTReleasePoints.from_cells(model, [0, 1]),
        porosity=0.25,
    )
    return project.run(silent=True)


def test_pathline_mosaic_without_groups_names_the_fix(ungrouped_run):
    """An ungrouped run says how to get groups instead of drawing one blank panel."""

    view = ungrouped_run.pathlines

    assert view.get()["release_group"].eq("").all()
    assert view.groups() == []
    with pytest.raises(ValueError, match="no release groups"):
        view.mosaic()

    # ...but the single-panel map still draws, colored per particle
    traces = view.map().overlay_traces()
    assert len(traces) == view.get()["particle"].nunique() == 2
    assert {trace.name for trace in traces} == set(view.get()["particle"])


def test_pathline_plot_draws_elevation_against_travel_time(prt_run):
    figure = prt_run.pathlines.plot()

    assert isinstance(figure, Fig)
    assert len(figure.data) == prt_run.pathlines.get()["particle"].nunique()
    assert figure.layout.xaxis.title.text == "Travel time"
    assert figure.layout.yaxis.title.text == "Elevation"


def test_pathline_map_and_plot_render_through_matplotlib(prt_run):
    """``backend="mpl"`` reuses the existing FloPy plan view rather than a blank map.

    Returns a bare ``Figure``, not FloPy's ``(fig, ax)``. Every other
    ``backend="mpl"`` in the grammar returns one figure, and a lone verb handing
    back a tuple means a caller who writes one loop over several pictures gets a
    ``TypeError`` from the odd one out (changed 2026-09-02 with the noun-tier
    signature pass; the axes are still reachable as ``figure.axes[0]``).
    """

    from matplotlib.figure import Figure

    figure = prt_run.pathlines.map(backend="mpl")
    assert isinstance(figure, Figure)
    assert figure.axes[0].get_title() == "Particle pathlines"
    assert isinstance(prt_run.pathlines.plot(backend="mpl"), Figure)


def test_pathline_view_refuses_the_verbs_it_cannot_answer(prt_run):
    with pytest.raises(NotImplementedError, match="cross-section"):
        prt_run.pathlines.section()
    with pytest.raises(NotImplementedError, match="no period"):
        prt_run.pathlines.animate()
    with pytest.raises(ValueError, match="not a pathline coloring"):
        prt_run.pathlines.map(color="layer")


@pytest.fixture(scope="module")
def stratified_run(prt_run, tmp_path_factory):
    """Six WEST particles surrounding a single EAST one, in sorted-id order.

    Built so a cap that samples the *whole sorted list* -- either "the first N" or
    an even spread across it -- misses EAST entirely, while a per-group quota
    cannot. Without this arrangement a flat sample happens to pick both groups and
    the stratification is untested.
    """

    workspace = tmp_path_factory.mktemp("prt_stratified")
    releases = [1, 2, 3, 4, 5, 6, 7]
    pd.DataFrame(
        {
            "imdl": [1] * 7,
            "iprp": [1] * 7,
            "irpt": releases,
            "ilay": [1] * 7,
            "icell": [1] * 7,
            "ireason": [3] * 7,
            "trelease": [0.0] * 7,
            "t": [float(irpt) for irpt in releases],
            "x": [0.5] * 7,
            "y": [0.5] * 7,
            "z": [5.0] * 7,
            # EAST sits in the MIDDLE of the sorted particle ids
            "name": ["WEST", "WEST", "WEST", "EAST", "WEST", "WEST", "WEST"],
        }
    ).to_csv(workspace / "stratified.trk.csv", index=False)
    return open_prt_run(prt_run.flow_model, workspace)


@pytest.fixture(scope="module")
def summary_run(prt_run, tmp_path_factory):
    """Two particles with several records each; one never terminates."""

    workspace = tmp_path_factory.mktemp("prt_summary")
    pd.DataFrame(
        {
            "imdl": [1] * 5,
            "iprp": [1] * 5,
            "irpt": [1, 1, 1, 2, 2],
            "ilay": [1] * 5,
            # rows deliberately out of time order for particle 1, so every column
            # below reads chronologically as: cell 0 -> cell 1 -> cell 1 (terminates)
            "icell": [2, 2, 1, 1, 1],
            "ireason": [3, 1, 0, 0, 1],  # particle 2 is still in transit
            "trelease": [0.0] * 5,
            "t": [10.0, 5.0, 0.0, 0.0, 3.0],
            "x": [1.5, 0.8, 0.2, 0.5, 0.6],
            "y": [0.5] * 5,
            "z": [3.0, 4.0, 5.0, 5.0, 4.5],
            "name": ["WEST", "WEST", "WEST", "EAST", "EAST"],
        }
    ).to_csv(workspace / "summary.trk.csv", index=False)
    return open_prt_run(prt_run.flow_model, workspace)


def test_particle_cap_is_stratified_across_release_groups(stratified_run):
    """A cap must not delete a whole capture zone -- that changes what the map says."""

    view = stratified_run.pathlines
    frame = view.get()
    assert frame["particle"].nunique() == 7

    chosen, dropped = view._selected_particles(frame, 2)
    assert dropped == 5
    # one per group, evenly spaced inside each: EAST's only particle + WEST's first
    assert chosen == ["1-1-1-0", "1-1-4-0"]
    # the arrangement really would defeat a whole-list sample
    assert sorted(frame["particle"].unique())[:2] == ["1-1-1-0", "1-1-2-0"]

    with pytest.warns(UserWarning, match="2 of 7 particles"):
        choro = view.map(max_particles=2)
    assert {trace.name for trace in choro.overlay_traces()} == {"WEST", "EAST"}

    # a cap at or above the population changes nothing
    assert view._selected_particles(frame, 7) == (sorted(frame["particle"].unique()), 0)
    assert view._selected_particles(frame, None)[1] == 0


def test_pathline_summary_reports_each_particle_start_end_and_fate(summary_run):
    digest = summary_run.pathlines.summary().set_index("particle")

    first = digest.loc["1-1-1-0"]
    assert first["records"] == 3
    # rows arrive out of order in the CSV; the digest is by time, not by row
    assert (first["start_cell"], first["end_cell"]) == (0, 1)
    assert first["travel_time"] == 10.0
    assert bool(first["terminated"]) is True

    second = digest.loc["1-1-2-0"]
    assert second["records"] == 2
    assert (second["start_cell"], second["end_cell"]) == (0, 0)
    assert second["travel_time"] == 3.0
    # still moving when tracking stopped -- not a termination
    assert bool(second["terminated"]) is False


def test_pathline_vertices_follow_time_not_csv_row_order(summary_run):
    """A track drawn in file order is a visibly wrong figure."""

    traces = {trace.name: trace for trace in summary_run.pathlines.map().overlay_traces()}
    west = traces["WEST"]

    # x rises with t in the fixture, so the drawn line must run west to east
    assert list(west.lon) == sorted(west.lon)
    assert len(west.lon) == 3


def test_pathline_map_styling_and_legend_group_one_entry_per_group(prt_run):
    view = prt_run.pathlines

    thick = view.map(width=7.5).overlay_traces()
    assert {trace.line.width for trace in thick} == {7.5}

    # one legend entry per group, with the whole group toggling together
    traces = view.map().overlay_traces()
    assert [trace.legendgroup for trace in traces] == [trace.name for trace in traces]
    for name in {trace.name for trace in traces}:
        shown = [trace.showlegend for trace in traces if trace.name == name]
        assert sum(bool(flag) for flag in shown) == 1

    # color="particle" names traces by particle instead of by group
    per_particle = view.map(color="particle").overlay_traces()
    assert {trace.name for trace in per_particle} == set(view.get()["particle"])


def test_pathline_map_on_a_borrowed_base_adds_to_that_map(prt_run):
    """`base=<Choro>` draws onto the map you passed -- and says so, rather than copying."""

    base = prt_run.travel_time.map()
    base.fig.update_layout(title="time of travel")
    returned = prt_run.pathlines.map(base=base)

    assert returned is base
    assert len(base.overlay_traces()) == 2
    # a borrowed map keeps its own title unless the caller names one
    assert base.fig.layout.title.text == "time of travel"
    assert prt_run.pathlines.map(base=base, title="both").fig.layout.title.text == "both"
    # ...and drawing twice really does draw twice (documented, not silently deduped)
    assert len(base.overlay_traces()) == 4


def test_pathline_map_rejects_an_unknown_release_group(prt_run):
    with pytest.raises(KeyError, match="WEST_WELLS"):
        prt_run.pathlines.map(group="nowhere")


def test_pathline_mosaic_passes_layout_options_through(two_group_run):
    fig = two_group_run.pathlines.mosaic(ncols=1, title="by group", sync_views=False)

    assert fig.layout.title.text == "by group"
    # ncols=1 stacks the panels: two rows, so two distinct map subplots
    assert {trace.subplot for trace in fig.data} == {"map", "map2"}
    assert not any("plotly_relayout" in script for script in getattr(fig, "_post_scripts", []))


def test_pathline_views_survive_a_run_with_no_particles(prt_run, tmp_path):
    """An empty track CSV yields empty tables and empty figures, not exceptions."""

    columns = ["imdl", "iprp", "irpt", "ilay", "icell", "ireason", "trelease", "t", "x", "y", "z"]
    pd.DataFrame(columns=columns).to_csv(tmp_path / "empty.trk.csv", index=False)
    empty = open_prt_run(prt_run.flow_model, tmp_path)

    view = empty.pathlines
    assert view.get().empty
    assert view.groups() == []
    assert len(view.summary()) == 0
    assert view.map().overlay_traces() == []
    assert len(view.plot().data) == 0


@pytest.fixture(scope="module")
def mixed_group_run(tmp_path_factory):
    """A REAL run merging a named release set with an un-named one.

    Verified against MF6 rather than assumed: merging makes the PRP carry
    BOUNDNAMES, and MF6 then **synthesizes** a label for the points that had none
    (``PRP000000002``) rather than leaving the name blank. So a mixed release
    yields two groups, one of them auto-named -- which is what shows up in the
    legend and in ``capture``.
    """

    workspace = tmp_path_factory.mktemp("prt_mixed")
    model = _flow_model("mixed_flow", workspace)
    assert model.run_simulation()

    release = PRTReleasePoints.merge(
        PRTReleasePoints.from_cells(model, [0], group="west_wells"),
        PRTReleasePoints.from_cells(model, [1]),  # deliberately un-named
    )
    project = PRTProject(
        model, workspace=workspace / "prt", release_points=release, porosity=0.25
    )
    return project.run(silent=True)


def test_a_mixed_named_and_unnamed_release_gets_an_mf6_synthesized_group(mixed_group_run):
    """The un-named half is not blank -- MF6 names it, and it colors like any group."""

    groups = mixed_group_run.pathlines.groups()
    assert groups == ["PRP000000002", "WEST_WELLS"]

    traces = mixed_group_run.pathlines.map().overlay_traces()
    assert {trace.name for trace in traces} == {"WEST_WELLS", "PRP000000002"}
    assert len({trace.line.color for trace in traces}) == 2
    # the same synthesized label reaches the capture zones
    assert set(mixed_group_run.capture.groups()) == {"WEST_WELLS", "PRP000000002"}


def test_blank_group_labels_from_an_external_track_csv_are_drawn_not_dropped(prt_run, tmp_path):
    """Defensive: a track CSV myflopy did not write may mix labels with blanks.

    MF6 itself does not produce this (it synthesizes a name -- see the test
    above), but ``open_prt_run`` accepts any MODFLOW-6 track CSV, and a blank
    label must not crash the color lookup.
    """

    pd.DataFrame(
        {
            "imdl": [1, 1],
            "iprp": [1, 2],
            "irpt": [1, 1],
            "ilay": [1, 1],
            "icell": [1, 2],
            "ireason": [3, 3],
            "trelease": [0.0, 0.0],
            "t": [4.0, 6.0],
            "x": [0.5, 1.5],
            "y": [0.5, 0.5],
            "z": [5.0, 5.0],
            "name": ["WEST", None],
        }
    ).to_csv(tmp_path / "external.trk.csv", index=False)
    external = open_prt_run(prt_run.flow_model, tmp_path)

    traces = external.pathlines.map().overlay_traces()

    assert {trace.name for trace in traces} == {"WEST", "particles"}
    assert len({trace.line.color for trace in traces}) == 2


def test_per_particle_coloring_does_not_pollute_the_shared_memo(prt_run, isolated_category_colors):
    """Particle ids are not a recurring category, so they stay out of the memo.

    Otherwise one 250-particle map would register 250 names and shift the colors
    every later figure's release groups receive.
    """

    prt_run.pathlines.map(color="particle")
    assert isolated_category_colors == {}

    prt_run.pathlines.map()  # group coloring is remembered
    assert set(isolated_category_colors) == {"WEST_WELLS", "EAST_WELLS"}


def test_pathline_map_refuses_options_it_cannot_honor(prt_run):
    """Selectors that would be silently dropped raise instead."""

    view = prt_run.pathlines

    # the matplotlib plan view takes no base-map options
    with pytest.raises(TypeError, match="no base-map options"):
        view.map(backend="mpl", contours=True)

    # an already-built base cannot be reshaped by this call
    base = prt_run.travel_time.map()
    with pytest.raises(TypeError, match=r"\['layer'\]"):
        view.map(base=base, layer=0)
    with pytest.raises(TypeError, match="contours"):
        view.map(base=base, contours=True)

    # base=None owns the options that blank the cells
    with pytest.raises(TypeError, match="showscale"):
        view.map(base=None, showscale=True)

    # a bad color is caught even when the selection is empty
    empty = view.get().iloc[0:0]
    assert empty.empty
    with pytest.raises(ValueError, match="not a pathline coloring"):
        view.map(color="layer")


def test_pathline_mosaic_refuses_a_shared_base_and_the_mpl_backend(two_group_run, prt_run):
    view = two_group_run.pathlines

    with pytest.raises(ValueError, match="cannot back several panels"):
        view.mosaic(base=prt_run.travel_time.map())
    with pytest.raises(ValueError, match="Plotly composition"):
        view.mosaic(backend="mpl")


def test_naming_the_base_map_options_did_not_empty_the_checks_that_read_them(prt_run):
    """`_explicit_options` must reconstruct what the `**base_kwargs` tail carried.

    The 8.8 conversion named ~21 base-map parameters that used to arrive through
    `**base_kwargs`, and `map` READ that tail three times: to refuse base-map
    options on the Matplotlib branch, to refuse them on an already-built base,
    and to forward them. Naming them emptied it, so all three would have gone on
    passing silently -- the general hazard of this conversion (ledger 176), and
    the reason it gets a test rather than a comment.

    `zmin` is the forwarding probe, NOT `contours`: this fixture's grid is two
    cells, and two points produce no contour traces at all, so a trace count
    would read as "never forwarded" when the option had arrived correctly.
    (Measured against the canonical grid the same call goes 1 trace to 15, and
    instrumenting `_explicit_options` here shows
    `{'zmin': 1.0, 'contours': True, 'contour_levels': 5}` reaching the base map.)
    """

    import matplotlib

    matplotlib.use("Agg")
    from matplotlib.figure import Figure

    view = prt_run.pathlines

    # forwarded: the option must reach the base map's choropleth trace
    pinned = view.map(zmin=100.0, zmax=125.0)
    limits = [getattr(trace, "zmin", None) for trace in pinned.fig.data]
    assert 100.0 in limits, f"zmin= never reached the base map (saw {limits})"

    # refused: FloPy's plan view takes no base-map options, and an already-built
    # base cannot be reshaped after the fact
    with pytest.raises(TypeError, match="no base-map options"):
        view.map(backend="mpl", contours=True)
    with pytest.raises(TypeError, match="already-built map"):
        view.map(base=view.map(), zmin=1.0)

    # and the plain Matplotlib branch still returns ONE bare figure
    assert isinstance(view.map(backend="mpl"), Figure)

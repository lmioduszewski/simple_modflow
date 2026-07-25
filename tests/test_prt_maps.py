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
    assert any(value.endswith("\xa0d") for value in rendered), rendered
    # min/max are elapsed travel times, so they are not labelled as arrival clock times
    assert "fastest" in trace.hovertemplate and "slowest" in trace.hovertemplate
    assert "arrival" not in trace.hovertemplate

    # an explicit unit still wins, and an unreadable TDIS degrades to no unit
    explicit = prt_run.travel_time.map(units="yr").get_choropleth()
    assert any(str(value).endswith("\xa0yr") for value in explicit.customdata[1])


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

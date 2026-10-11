"""Horizontal flow barriers (ledger 162).

A barrier sits on the FACE between two cells, so MODFLOW 6 addresses it as a pair
of cells -- and FloPy will take any pair at all. It validates each cellid against
`idomain` and never checks that the two are neighbours, so a wrong pair is
discovered by MODFLOW 6 aborting mid-run. What a modeller actually has is a line.
Both halves of that gap are what `mf.hfb` exists to close.
"""

from __future__ import annotations

import pytest
import shapely as shp

import myflopy as mf
from myflopy.modflow.mf6.grid.barriers import (
    barrier_faces,
    enclosed_faces,
    is_watertight,
    shared_face_segment,
    validate_barrier_pairs,
)
from myflopy.specs import ModelContext


@pytest.fixture
def ctx(canonical_run):
    return ModelContext(grid=canonical_run.vor)


@pytest.fixture
def crossing_line(canonical_run):
    """A wall straight across the middle of the grid."""

    bounds = canonical_run.vor.gdf_vorPolys.total_bounds
    middle = (bounds[1] + bounds[3]) / 2
    return shp.LineString([(bounds[0], middle), (bounds[2], middle)])


def _records(spec):
    return spec.options["stress_period_data"][0]


# --- the geometry -------------------------------------------------------------


def test_a_line_resolves_to_the_connections_it_cuts(canonical_run, crossing_line):
    faces = barrier_faces(canonical_run.vor, crossing_line)
    assert faces
    assert all(a < b for a, b in faces), "pairs come back canonically ordered"
    assert len(faces) == len(set(faces)), "and unique"


def test_every_cut_is_a_real_connection(canonical_run, crossing_line):
    """Every face the resolver returns must be one MODFLOW 6 agrees exists.

    `vor.adjacent_cells_idx` was verified to be exactly MF6's DISV lateral
    connectivity, so this is the property that makes the resolver trustworthy.
    """

    adjacency = canonical_run.vor.adjacent_cells_idx
    for cell_a, cell_b in barrier_faces(canonical_run.vor, crossing_line):
        assert cell_b in {int(n) for n in adjacency[cell_a]}
        assert shared_face_segment(canonical_run.vor, cell_a, cell_b) is not None


def test_a_stub_inside_one_cell_cuts_nothing(canonical_run):
    centroid = canonical_run.vor.gdf_vorPolys.geometry.iloc[0].centroid
    stub = shp.LineString([(centroid.x, centroid.y), (centroid.x + 0.01, centroid.y)])
    assert barrier_faces(canonical_run.vor, stub) == []


def test_cells_that_only_touch_share_no_face(canonical_run):
    """A zero-length intersection is a corner touch, not a connection."""

    vor = canonical_run.vor
    far = shared_face_segment(vor, 0, len(vor.gdf_vorPolys) - 1)
    assert far is None


# --- the closed-wall trap -----------------------------------------------------


def test_a_ring_seals_its_interior_like_the_enclosing_cut(canonical_run):
    """A ring through cells used to leave every cell connected: it took the faces the
    ring CROSSES, which run along it. It cuts the connections across it now (ledger
    216) and seals like `.enclose` -- which stays, because it decides inside by
    centroid and refuses a cut that leaks.
    """

    vor = canonical_run.vor
    bounds = vor.gdf_vorPolys.total_bounds
    cx, cy = (bounds[0] + bounds[2]) / 2, (bounds[1] + bounds[3]) / 2
    polygon = shp.box(cx - 540, cy - 540, cx + 540, cy + 540)

    cut, interior = enclosed_faces(vor, polygon)
    assert interior
    assert is_watertight(vor, cut, interior[0])

    ring = barrier_faces(vor, polygon.exterior)
    assert is_watertight(vor, ring, interior[0])


# --- validation, which is the other half of the value -------------------------


def test_a_non_adjacent_pair_is_refused(canonical_run, ctx):
    """FloPy passes it; MODFLOW 6 aborts mid-run. Catch it at the call."""

    last = len(canonical_run.vor.gdf_vorPolys) - 1
    with pytest.raises(ValueError, match="unconnected"):
        mf.hfb(pairs=[(0, last)], hydchr=1e-6, context=ctx)


def test_a_self_pair_is_refused(ctx):
    with pytest.raises(ValueError, match="itself"):
        mf.hfb(pairs=[(5, 5)], hydchr=1e-6, context=ctx)


def test_a_duplicated_face_is_collapsed(canonical_run, ctx, crossing_line):
    """The most dangerous input there is, and MODFLOW 6 says nothing about it.

    A repeated face makes MF6 apply its series formula twice, and `condsat_reset`
    then restores the ALREADY-MODIFIED conductance -- so the face stays wrong for
    the rest of the run, even after an empty period removes the barrier.
    """

    a, b = barrier_faces(canonical_run.vor, crossing_line)[0]
    spec = mf.hfb(pairs=[(a, b), (b, a), (a, b)], hydchr=1e-6, context=ctx)
    assert len(_records(spec)) == 1


def test_validation_is_skipped_without_a_grid():
    """No context means no grid means nothing to validate against -- not a crash."""

    spec = mf.hfb(pairs=[(0, 999999)], hydchr=1e-6)
    assert len(_records(spec)) == 1


def test_a_vertical_pair_is_accepted(canonical_run):
    """Legal on DISV since MODFLOW 6 6.7.0; FloPy 3.10's dfn still says otherwise."""

    pairs = [((0, 5), (1, 5))]
    assert validate_barrier_pairs(canonical_run.vor, pairs) == pairs


def test_a_layer_skipping_vertical_pair_is_refused(canonical_run):
    with pytest.raises(ValueError, match="unconnected"):
        validate_barrier_pairs(canonical_run.vor, [((0, 5), (2, 5))])


# --- the helper ---------------------------------------------------------------


def test_pairs_expand_across_layers(canonical_run, ctx, crossing_line):
    faces = barrier_faces(canonical_run.vor, crossing_line)
    spec = mf.hfb(pairs=faces, hydchr=1e-6, layers=[0, 1], context=ctx)
    assert len(_records(spec)) == 2 * len(faces)
    assert {record[0][0] for record in _records(spec)} == {0, 1}


def test_qualified_pairs_reject_a_layers_argument(ctx):
    with pytest.raises(ValueError, match="ambiguous"):
        mf.hfb(pairs=[((0, 5), (0, 6))], hydchr=1e-6, layers=[0], context=ctx)


def test_one_hydchr_per_pair(canonical_run, ctx, crossing_line):
    faces = barrier_faces(canonical_run.vor, crossing_line)
    spec = mf.hfb(pairs=faces, hydchr=[1e-6] * len(faces), context=ctx)
    assert [record[2] for record in _records(spec)] == [1e-6] * len(faces)


def test_a_mismatched_hydchr_count_says_so(canonical_run, ctx, crossing_line):
    faces = barrier_faces(canonical_run.vor, crossing_line)
    with pytest.raises(ValueError, match="one value or one per pair"):
        mf.hfb(pairs=faces, hydchr=[1e-6, 2e-6], context=ctx)


def test_line_builds_a_spec_directly(canonical_run, ctx, crossing_line):
    spec = mf.hfb.line(crossing_line, context=ctx, hydchr=1e-6)
    assert _records(spec)
    assert spec.name == "hfb"


def test_enclose_seals_the_interior(canonical_run, ctx):
    bounds = canonical_run.vor.gdf_vorPolys.total_bounds
    cx, cy = (bounds[0] + bounds[2]) / 2, (bounds[1] + bounds[3]) / 2
    spec = mf.hfb.enclose(shp.box(cx - 540, cy - 540, cx + 540, cy + 540),
                          context=ctx, hydchr=1e-8)
    assert _records(spec)


def test_enclose_refuses_a_polygon_that_holds_no_cell(canonical_run, ctx):
    """A cell counts as inside when its CENTROID is; a polygon off the grid holds none."""

    bounds = canonical_run.vor.gdf_vorPolys.total_bounds
    away = shp.box(bounds[2] + 1000, bounds[3] + 1000, bounds[2] + 1100, bounds[3] + 1100)
    with pytest.raises(ValueError, match="encloses no cell"):
        mf.hfb.enclose(away, context=ctx, hydchr=1e-8)


def test_enclose_can_seal_a_single_cell(canonical_run, ctx):
    """Sealing one cell is legitimate, not an error -- it is the degenerate cut."""

    centroid = canonical_run.vor.gdf_vorPolys.geometry.iloc[0].centroid
    spec = mf.hfb.enclose(
        shp.Point(centroid.x, centroid.y).buffer(0.01), context=ctx, hydchr=1e-8
    )
    assert _records(spec)


def test_geometry_verbs_need_a_grid(crossing_line):
    with pytest.raises(ValueError, match="needs the grid"):
        mf.hfb.line(crossing_line, context=None, hydchr=1e-6)


def test_flopy_is_the_unvalidated_escape_hatch():
    """Deliberately no checks -- the pairs go to MODFLOW 6 as given."""

    spec = mf.hfb.flopy(stress_period_data={0: [[(0, 0), (0, 999999), 1e-6]]})
    assert _records(spec)


def test_maxhfb_is_not_exposed(canonical_run, ctx, crossing_line):
    """FloPy computes it from the data; a hand-set value that disagrees is a bug."""

    spec = mf.hfb.line(crossing_line, context=ctx, hydchr=1e-6)
    assert "maxhfb" not in spec.options


def test_hfb_is_not_a_registry_package():
    """It is face-indexed, has no cellid and no budget record, so it stays a plain factory."""

    from myflopy.modflow.mf6.package_registry import get_package_explorer_spec

    assert get_package_explorer_spec("hfb") is None


def test_hfb_is_still_discovered_on_reopen():
    """The one thing registry membership would have bought, for one line."""

    from myflopy.project.run_model import _NON_REGISTRY_SUFFIXES

    assert _NON_REGISTRY_SUFFIXES["hfb"] == "HFB"


@pytest.mark.slow
def test_a_resolved_barrier_runs_in_mf6(canonical_run_fresh, crossing_line):
    """The claim that matters: MODFLOW 6 accepts every pair the resolver produced."""

    import flopy

    model = canonical_run_fresh
    faces = barrier_faces(model.vor, crossing_line)
    records = [[[0, a], [0, b], 1e-6] for a, b in faces]
    flopy.mf6.ModflowGwfhfb(model.gwf, stress_period_data={0: records})
    model.sim.write_simulation(silent=True)
    success, report = model.run_simulation()
    assert success, "\n".join(report[-25:])


# --- the explorer -------------------------------------------------------------


@pytest.fixture
def barrier_run(canonical_run_fresh, crossing_line):
    """The canonical model with a wall of barriers straight across it."""

    import flopy

    model = canonical_run_fresh
    faces = barrier_faces(model.vor, crossing_line)
    flopy.mf6.ModflowGwfhfb(
        model.gwf,
        stress_period_data={0: [[[0, a], [0, b], 1e-8] for a, b in faces]},
    )
    model.sim.write_simulation(silent=True)
    success, report = model.run_simulation()
    assert success, "\n".join(report[-25:])
    return model


@pytest.mark.slow
def test_the_package_is_reachable_from_the_model(barrier_run):
    """`model.hfb` -- bespoke, because the generic explorer is cell-keyed."""

    assert barrier_run.hfb is not None


@pytest.mark.slow
def test_get_lists_every_barrier(barrier_run, crossing_line):
    table = barrier_run.hfb.get()
    assert len(table) == len(barrier_faces(barrier_run.vor, crossing_line))
    assert {"per", "layer", "cell1", "cell2", "hydchr", "length"} <= set(table.columns)
    # `length` is the shared face's own length -- what MF6 multiplies hydchr by
    assert (table["length"] > 0).all()
    assert not table["vertical"].any()


@pytest.mark.slow
def test_summary_is_one_row_per_period(barrier_run):
    summary = barrier_run.hfb.summary()
    assert len(summary) == 1
    assert summary.iloc[0]["barriers"] == len(barrier_run.hfb.get())


@pytest.mark.slow
def test_segments_are_the_faces_themselves(barrier_run):
    """A barrier's geometry is the shared EDGE, and an edge survives a re-grid."""

    frame = barrier_run.hfb.segments()
    assert len(frame) == len(barrier_run.hfb.get())
    assert frame.crs is not None
    assert (frame.geometry.geom_type == "LineString").all()


@pytest.mark.slow
def test_the_map_draws_barriers_as_lines(barrier_run):
    """A barrier is an edge, so it is drawn over the cell field, not instead of it."""

    types = [trace.type for trace in barrier_run.hfb.map().fig.data]
    assert types == ["choroplethmap", "scattermap"]


# --- the results tier, which exists after all ---------------------------------


@pytest.mark.slow
def test_flow_across_each_barrier_is_recoverable(barrier_run):
    """MODFLOW 6 writes no HFB budget record -- but the flow is not lost.

    A barrier sits ON a cell-to-cell connection, and FLOW-JA-FACE carries the
    flow on every connection, so a barrier's flow is simply that entry.
    """

    flows = barrier_run.hfb.results.q.get()
    assert not flows.empty
    assert "q_cell1_to_cell2" in flows.columns
    assert flows["q_cell1_to_cell2"].notna().all()
    assert (flows["q_cell1_to_cell2"] != 0).any()


@pytest.mark.slow
def test_there_is_one_flow_per_barrier_per_step(barrier_run):
    flows = barrier_run.hfb.results.q.get()
    barriers = len(barrier_run.hfb.get())
    assert len(flows) % barriers == 0
    for _, block in flows.groupby(["per", "kstp"]):
        assert len(block) == barriers


@pytest.mark.slow
def test_the_flow_column_names_its_frame(barrier_run):
    """MODFLOW 6's own sign is kept; the COLUMN says which way positive runs.

    Negating to make a number read "intuitively" is how sign bugs get built in.
    """

    flows = barrier_run.hfb.results.q.get()
    assert "q_cell1_to_cell2" in flows.columns
    assert "q" not in flows.columns


@pytest.mark.slow
def test_a_barrier_across_the_grid_stops_the_flow_across_it(canonical_run, crossing_line, tmp_path):
    """What a barrier is FOR. This used to sum |q| on the barrier's own faces -- near
    zero at a tiny hydchr whatever faces were picked -- and passed on a barrier that
    let 95% of the flow across (ledger 216). Here CHD holds the southern cells high
    and the northern ones low, so every unit of CHD inflow crosses the trace."""

    import flopy
    import numpy as np

    vor = canonical_run.vor
    y = np.asarray(vor.points)[:, 1]
    low, high = np.quantile(y, [0.08, 0.92])
    south = [int(c) for c in np.flatnonzero(y <= low)]
    north = [int(c) for c in np.flatnonzero(y >= high)]

    def inflow(pairs, label):
        sim = flopy.mf6.MFSimulation(sim_name=label, sim_ws=str(tmp_path / label), exe_name="mf6")
        flopy.mf6.ModflowTdis(sim)
        flopy.mf6.ModflowIms(sim, complexity="SIMPLE")
        gwf = flopy.mf6.ModflowGwf(sim, modelname=label, save_flows=True)
        flopy.mf6.ModflowGwfdisv(gwf, nlay=1, top=100.0, botm=[0.0], **vor.get_disv_gridprops())
        flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
        flopy.mf6.ModflowGwfic(gwf, strt=50.0)
        flopy.mf6.ModflowGwfchd(
            gwf,
            stress_period_data={
                0: [[(0, c), 90.0] for c in south] + [[(0, c), 10.0] for c in north]
            },
        )
        if pairs:
            flopy.mf6.ModflowGwfhfb(
                gwf, stress_period_data={0: [[(0, a), (0, b), 1e-10] for a, b in pairs]}
            )
        flopy.mf6.ModflowGwfoc(gwf, budget_filerecord=f"{label}.cbc", saverecord=[("BUDGET", "ALL")])
        sim.write_simulation(silent=True)
        success, report = sim.run_simulation(silent=True)
        assert success, "\n".join(report[-20:])
        chd = gwf.output.budget().get_data(text="CHD")[-1]
        return float(chd["q"][chd["q"] > 0].sum())

    open_flow = inflow([], "open")
    assert open_flow > 0
    assert inflow(barrier_faces(vor, crossing_line), "sealed") < 1e-4 * open_flow


@pytest.mark.slow
def test_results_summary_totals_the_barriers(barrier_run):
    summary = barrier_run.hfb.results.q.summary()
    assert {"per", "kstp", "barriers", "q_net", "q_abs_total", "q_max_abs"} <= set(summary.columns)
    assert (summary["q_abs_total"] >= summary["q_net"].abs()).all()


@pytest.mark.slow
def test_the_results_map_draws_one_line_per_barrier(barrier_run):
    picture = barrier_run.hfb.results.q.map()
    scatter = [t for t in picture.fig.data if t.type == "scattermap"]
    assert len(scatter) == len(barrier_run.hfb.get())


@pytest.mark.slow
def test_the_matplotlib_backend_draws_the_barriers_rather_than_dropping_them(barrier_run):
    """`backend="mpl"` on an HFB map must carry the barriers, not just the cells.

    `Choro.plot_mpl` reads the cell values and NOTHING else -- overlays are a
    Plotly-side concept it has never drawn. So the ordinary `_apply_backend`
    route, which is right for every other noun, would have returned a picture of
    the grid with the barriers silently missing: the one thing the picture is of.
    Both HFB verbs draw the segments onto the axes themselves instead.

    Counting Line2D artists is the falsifiable form. A `plot_mpl` passthrough
    gives zero of them while still returning a perfectly valid Figure, so
    asserting only on the type would pass over the defect.
    """

    import matplotlib

    matplotlib.use("Agg")
    from matplotlib.figure import Figure

    expected = len(barrier_run.hfb.segments())
    assert expected, "the fixture built no barriers, so this proves nothing"

    for picture in (
        barrier_run.hfb.inputs.map(backend="mpl"),
        barrier_run.hfb.results.q.map(backend="mpl"),
    ):
        assert isinstance(picture, Figure)
        assert len(picture.axes[0].lines) == expected


@pytest.mark.slow
def test_an_hfb_map_refuses_a_field_override_and_an_unknown_backend(barrier_run):
    """The noun rules reach HFB too, with one deliberate exception.

    `hfb.results.q` draws the flow ACROSS the barrier, so `values=` would repaint
    it while the hover went on reporting q -- refused, like every other noun.
    `hfb.inputs` is the exception: its subject is the barriers and the cells are
    a backdrop, so choosing what the backdrop shows contradicts nothing, exactly
    as on `vor.plot.map`.
    """

    with pytest.raises(TypeError, match="values"):
        barrier_run.hfb.results.q.map(values=[1.0] * barrier_run.vor.ncpl)

    # ... but the input map's `values` is its documented backdrop, and works.
    assert barrier_run.hfb.inputs.map(values=[1.0] * barrier_run.vor.ncpl)

    for verb in (barrier_run.hfb.inputs.map,
                 barrier_run.hfb.results.q.map):
        with pytest.raises(ValueError, match="backend must be"):
            verb(backend="nope")


# --- any discretization, ordered by TYPE (ledger 378) ---------------------


def _row_of_three(name="dis"):
    """Three 10 ft cells in a row: the smallest model a barrier can act on."""

    return mf.dis(
        nlay=1, nrow=1, ncol=3, delr=10.0, delc=10.0, top=10.0, botm=[0.0], name=name
    )


def _flow_packages():
    return [
        mf.npf(k=1.0),
        mf.ic(strt=5.0),
        mf.chd(stress_period_data={0: [[(0, 0, 0), 6.0], [(0, 0, 2), 4.0]]}),
        mf.oc(head_filerecord="m.hds", saverecord=[("HEAD", "ALL")]),
    ]


def _two_cell_disv(name):
    """Two square DISV cells sharing one face, under a name other than ``disv``."""

    vertices = [[0, 0.0, 10.0], [1, 10.0, 10.0], [2, 20.0, 10.0],
                [3, 0.0, 0.0], [4, 10.0, 0.0], [5, 20.0, 0.0]]
    cell2d = [[0, 5.0, 5.0, 4, 0, 1, 4, 3], [1, 15.0, 5.0, 4, 1, 2, 5, 4]]
    return mf.disv(
        nlay=1, ncpl=2, nvert=6, vertices=vertices, cell2d=cell2d,
        top=10.0, botm=[0.0], name=name,
    )


def test_a_barrier_names_no_spec_it_requires():
    """`requires=("disv",)` named a SPEC, so a DIS, a DISU or a renamed DISV was refused."""

    assert mf.hfb(pairs=[(0, 1)], hydchr=1.0).requires == ()
    assert mf.hfb.flopy(stress_period_data={0: []}).requires == ()


@pytest.mark.slow
def test_a_barrier_runs_on_a_structured_model(tmp_path):
    """Built AND run on DIS: the barrier moves the middle head from 5.0 to the far end's."""

    import flopy

    heads = {}
    for label, extra in (
        ("open", []),
        ("barrier", [mf.hfb(pairs=[((0, 0, 0), (0, 0, 1))], hydchr=1e-6)]),
    ):
        model = mf.gwf("m", packages=[_row_of_three("grid"), *extra, *_flow_packages()])
        run = mf.simulation(model, name="s").build(workspace=tmp_path / label)
        run.execute()
        heads[label] = flopy.utils.HeadFile(tmp_path / label / "m.hds").get_data().ravel()
    assert heads["open"][1] == pytest.approx(5.0)
    assert heads["barrier"][1] == pytest.approx(4.0, abs=1e-3)


def test_a_barrier_builds_beside_a_renamed_disv(tmp_path):
    model = mf.gwf(
        "m",
        packages=[
            _two_cell_disv("mesh"),
            mf.hfb(pairs=[(0, 1)], hydchr=1e-6),
            mf.npf(k=1.0),
            mf.ic(strt=5.0),
        ],
    )
    built = mf.simulation(model, name="s").build(workspace=tmp_path)
    assert built.built.models["m"].packages["hfb"].package_type == "hfb"


def test_a_barrier_builds_on_a_disu_model(tmp_path):
    disu = mf.disu(
        nodes=2, nja=4, top=10.0, bot=0.0, area=100.0, iac=[2, 2], ja=[0, 1, 1, 0],
        ihc=[1, 1, 1, 1], cl12=[0.0, 5.0, 0.0, 5.0], hwva=[0.0, 10.0, 0.0, 10.0],
    )
    model = mf.gwf(
        "m",
        packages=[disu, mf.hfb(pairs=[((0,), (1,))], hydchr=1e-6), mf.npf(k=1.0), mf.ic(strt=5.0)],
    )
    built = mf.simulation(model, name="s").build(workspace=tmp_path)
    assert built.built.models["m"].packages["hfb"].package_type == "hfb"


def test_a_barrier_ahead_of_its_discretization_is_refused_by_name(tmp_path):
    """FloPy builds it and fails only at WRITE, naming neither the package nor the fix."""

    model = mf.gwf(
        "m",
        packages=[mf.hfb(pairs=[((0, 0, 0), (0, 0, 1))], hydchr=1e-6), _row_of_three(), *_flow_packages()],
    )
    with pytest.raises(ValueError, match=r"'hfb'.*declared BEFORE it.*mf\.dis"):
        mf.simulation(model, name="s").build(workspace=tmp_path)


def test_a_plan_pair_on_a_structured_model_is_refused_naming_the_spelling(tmp_path):
    model = mf.gwf(
        "m", packages=[_row_of_three(), mf.hfb(pairs=[(0, 1)], hydchr=1e-6), *_flow_packages()]
    )
    with pytest.raises(ValueError, match=r"DIS, whose cellid is \(layer, row, column\)"):
        mf.simulation(model, name="s").build(workspace=tmp_path)


def test_plan_and_full_pairs_do_not_mix():
    with pytest.raises(ValueError, match="not a mixture"):
        mf.hfb(pairs=[(0, 1), ((0, 1), (0, 2))], hydchr=1e-6)
    with pytest.raises(ValueError, match="one shape"):
        mf.hfb(pairs=[((0, 1), (0, 2)), ((0, 0, 1), (0, 0, 2))], hydchr=1e-6)


def test_a_structured_pair_against_a_disv_context_is_refused(ctx):
    with pytest.raises(ValueError, match=r"DISV mesh, whose cellid is \(layer, cell\)"):
        mf.hfb(pairs=[((0, 0, 0), (0, 0, 1))], hydchr=1e-6, context=ctx)

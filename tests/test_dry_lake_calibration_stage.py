"""A DRY lake is calibrated against its BOTTOM, on both routes (ledger 351).

The user's rule. MF6's two outputs disagree at a dry lake: its STAGE observation
reports the stage the solver clamped at the lake's bottom (``gwf-lak.f90``,
``lak_bd_obs`` writes ``xnewpak``, clamped at ``lakebot``), while the stage file
writes -1e30 (``lak_ot_dv``). ``LakeStageTargets`` read the first on its observation
route and ``NaN`` from the second on its stage-file route, which PEST's files then
wrote as 1e30 -- one lake, observed one way, phi 2.9 on one route and 3e60 on the
other, and PESTPP-IES stalled on the 1e30 (``docs/reproducers/dry_lake_bottom.py``).
Now the stage-file route puts the lake's bottom in place of -1e30, computed by MF6's
own rule (``package_tables.lake_bottoms``), and an INACTIVE lake (+1e30, MF6's
``STATUS INACTIVE``) stays no-value: ``NaN`` in review, 1e30 in the PEST files. The
stage VIEWS keep a dry day as a gap (ledger 329; ``test_lak_stage_sentinels``).

Four real MF6 runs, each six 10-day periods of one step, each declaring a LAK
STAGE observation for every lake (upper-case labels ``ONE``, ``TWO`` ...), so a
target named for it reads MF6's own number and a target named anything else reads
the stage file -- the two routes on one run:

``vertical``    a 6 x 6 irregular Voronoi grid, one confined layer (top 100 ft, one
                cell at 99.123456789123 -- more digits than FloPy writes), CHD 50 on
                the west and east edges. ``one``: ONE vertical connection (bottom
                100); ``two``: TWO, under tops 100 and 99.123456789123 (MF6's bottom
                is the LOWER, as FloPy wrote it, 99.12345679); both drawn down by a
                withdrawal, dry in periods 1-3 (``two`` in 3), refilled by rain in 4.
                ``off``: STATUS INACTIVE in every period.
``horizontal``  the same grid in two confined layers (100 / 60 / 40), CHD 30 in both.
                Two lakes with ONE HORIZONTAL connection each into a layer-1 cell --
                ``hneq`` with BELEV 45 < TELEV 58 (bottom 45) and ``heq`` with BELEV =
                TELEV = 55, which MF6 reads as the whole face (bottom 40, the cell's
                bottom, not the 55 written). A horizontal connection has no area, so
                both hold no water: dry in EVERY record.
``excavated``   a 12 x 12 grid of 20 ft cells, a 100 / 80 / 40 stack with a round pit
                EXCAVATED to 82.3456789012345 ft (``LayerStack.excavate``) and
                ``mf.lak`` on it -- ledger 325's rule: VERTICAL connections under the
                bed, HORIZONTAL walls into the neighbours' cells from the bed up. NEWTON,
                CHD 60 in the lower layer; drawn down dry in periods 1-3.
``implicit``    the excavated run with ``mf.lak(implicit=True)`` -- MF6's LAK
                ``IMPLICIT`` (ledger 326), the formulation the user's pond runs, where
                the stage is an unknown of the flow matrix. The SAME clamp holds it at
                the bottom: ``lak_cq`` re-runs ``lak_solve(update=.false.)``, which
                raises ``xnewpak`` to ``lakebot`` in the solution vector the implicit
                stage lives in, before either output is written. (``lak_nur`` acts
                only under ``NEWTON UNDER_RELAXATION``; this name file reads
                ``NEWTON``.)

What each file holds is asserted before it is used, not assumed.
"""

from __future__ import annotations

import warnings
from types import SimpleNamespace

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import Point

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.package_tables import lake_bottoms
from myflopy.modflow.mf6.pest.forward_run import write_named_series_targets
from myflopy.modflow.mf6.pest.observations import prepare_lake_stage_observations
from myflopy.modflow.mf6.pest.specs import LakeStageObservationSpec

NPER = 6
#: The top of one cell under lake ``two``, and what FloPy WRITES for it (8 decimals).
ODD_TOP, ODD_TOP_WRITTEN = 99.123456789123, 99.12345679
PIT_BED, PIT_BED_WRITTEN = 82.3456789012345, 82.3456789
#: MF6's bottom of each lake, as MF6 read the files: ``{run: {lake: bottom}}``.
BOTTOMS = {
    "vertical": {0: 100.0, 1: ODD_TOP_WRITTEN, 2: 100.0},
    "horizontal": {0: 45.0, 1: 40.0},
    "excavated": {0: PIT_BED_WRITTEN},
    "implicit": {0: PIT_BED_WRITTEN},
}
#: The records MF6 wrote -1e30 in the stage file for: ``{run: {lake: [record, ...]}}``.
DRY = {
    "vertical": {0: [1, 2, 3], 1: [3]},
    "horizontal": {0: list(range(NPER)), 1: list(range(NPER))},
    "excavated": {0: [1, 2, 3]},
    "implicit": {0: [1, 2, 3]},
}
NAMES = {
    "vertical": {"one": 0, "two": 1, "off": 2},
    "horizontal": {"hneq": 0, "heq": 1},
    "excavated": {"pit": 0},
    "implicit": {"pit": 0},
}
LABELS = ["vertical", "horizontal", "excavated", "implicit"]


def _grid():
    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
    interior = [c for c in range(grid.ncpl) if c not in set(edge.tolist())]
    return grid, props, [int(c) for c in edge], interior


def _observations(label):
    return mf.LakeStageTargets(locations=NAMES[label]).to_flopy_obs(csv_name="lake_stage_targets.csv")


def _run(root, model, packages=()):
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=NPER, perioddata=[(10.0, 1, 1.0)] * NPER), *packages])
    with warnings.catch_warnings():
        # A built grid object cannot be frozen into the run's snapshot, and these runs
        # are reopened from their FILES (mf.load_mf6_run), never from it.
        warnings.simplefilter("ignore", UserWarning)
        run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success, "\n".join(run.report[-20:])
    return run


def _vertical(root):
    grid, props, edge, interior = _grid()
    c0, c1, c2, c3 = interior[0], interior[5], interior[6], interior[10]
    top = np.full(grid.ncpl, 100.0)
    top[c2] = ODD_TOP
    lake = mf.lak.flopy(
        nlakes=3,
        packagedata=[[0, 100.5, 1, "one"], [1, 100.5, 2, "two"], [2, 104.0, 1, "off"]],
        connectiondata=[
            [0, 0, (0, c0), "VERTICAL", 0.5, 0.0, 0.0, 0.0, 0.0],
            [1, 0, (0, c1), "VERTICAL", 0.5, 0.0, 0.0, 0.0, 0.0],
            [1, 1, (0, c2), "VERTICAL", 0.05, 0.0, 0.0, 0.0, 0.0],
            [2, 0, (0, c3), "VERTICAL", 0.01, 0.0, 0.0, 0.0, 0.0],
        ],
        perioddata={
            0: [[0, "WITHDRAWAL", 300.0], [1, "WITHDRAWAL", 300.0], [2, "STATUS", "INACTIVE"]],
            4: [[0, "RAINFALL", 0.5], [1, "RAINFALL", 0.5]],
            5: [[0, "RAINFALL", 0.0], [1, "RAINFALL", 0.0]],
        },
        boundnames=True, stage_filerecord="m.lak.stg", budget_filerecord="m.lak.bud",
        observations=_observations("vertical"),
    )
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=top, botm=[0.0]),
        mf.npf(k=10.0, icelltype=0), mf.sto(ss=1e-5, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=50.0), mf.chd(stress_period_data={0: [[(0, c), 50.0] for c in edge]}),
        lake,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    return _run(root, model)


def _horizontal(root):
    grid, props, edge, interior = _grid()
    lake = mf.lak.flopy(
        nlakes=2,
        packagedata=[[0, 57.0, 1, "hneq"], [1, 57.0, 1, "heq"]],
        connectiondata=[
            [0, 0, (1, interior[1]), "HORIZONTAL", 0.5, 45.0, 58.0, 50.0, 100.0],
            [1, 0, (1, interior[9]), "HORIZONTAL", 0.5, 55.0, 55.0, 50.0, 100.0],
        ],
        perioddata={0: [[0, "WITHDRAWAL", 3000.0], [1, "WITHDRAWAL", 3000.0]]},
        boundnames=True, stage_filerecord="m.lak.stg", budget_filerecord="m.lak.bud",
        observations=_observations("horizontal"),
    )
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=2, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=100.0, botm=[60.0, 40.0]),
        mf.npf(k=10.0, icelltype=0), mf.sto(ss=1e-5, sy=0.2, iconvert=0, transient={0: True}),
        mf.ic(strt=30.0),
        mf.chd(stress_period_data={0: [[(k, c), 30.0] for c in edge for k in (0, 1)]}),
        lake,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    return _run(root, model)


def _excavated(root, *, implicit=False):
    grid = irregular_voronoi_grid(nrow=12, ncol=12, cell_size=20.0, refined=False)
    xmin, ymin, xmax, ymax = grid.gdf_vorPolys.total_bounds
    pond = gpd.GeoDataFrame({"name": ["pit"]}, crs=grid.crs,
                            geometry=[Point((xmin + xmax) / 2, (ymin + ymax) / 2).buffer(30.0, 64)])
    stack = (mf.LayerStack(top=mf.Flat(100.0)).add("upper", bottom=mf.Flat(80.0))
             .add("lower", bottom=mf.Flat(40.0)))
    result = stack.excavate(bottom=mf.Flat(PIT_BED), within=pond, name="pit").build(grid)
    xy = np.array([(c[1], c[2]) for c in grid.get_disv_gridprops()["cell2d"]])
    edge = np.flatnonzero((xy[:, 0] < xmin + 20) | (xy[:, 0] > xmax - 20))
    context = mf.ModelContext(grid=grid, domain=result.idomain, surfaces=result)
    lake = mf.lak(context=context, nper=NPER, lakes=pond, lake_id_field="name",
                  starting_stage=95.0, bed_leakance=0.5, unexcavated="raise", implicit=implicit,
                  withdrawals=[5000.0, 5000.0, 5000.0, 0.0, 0.0, 0.0],
                  rainfall=[0.0, 0.0, 0.0, 0.0, 2.0, 0.0],
                  stage_filerecord="m.lak.stg", budget_filerecord="m.lak.bud",
                  observations=_observations("implicit" if implicit else "excavated"))
    model = mf.gwf("m", context=context, newtonoptions="NEWTON", packages=[
        mf.disv(result), mf.npf(k=10.0, k33=10.0, icelltype=1),
        mf.sto(ss=1e-5, sy=0.2, iconvert=1, transient={0: True}),
        mf.ic(strt=60.0),
        mf.chd(stress_period_data={0: [[(1, int(c)), 60.0] for c in edge if result.idomain[1, c] > 0]}),
        lake,
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    return _run(root, model, [mf.ims(models=["m"], complexity="COMPLEX",
                                     linear_acceleration="BICGSTAB")])


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    """``{label: run}`` for the four runs of the module docstring."""

    return {label: build(tmp_path_factory.mktemp(label)) for label, build in (
        ("vertical", _vertical), ("horizontal", _horizontal), ("excavated", _excavated),
        ("implicit", lambda root: _excavated(root, implicit=True)))}


def _loaded(run):
    """The model as MF6 read it -- from its files, as PEST's forward run reloads it."""

    return mf.load_mf6_run(run.workspace, model_name="m", verbosity_level=0)


def _files(label, run):
    """``(raw stage file, MF6's STAGE observations)``, asserting the dry records."""

    model = run.model()
    raw = model.outputs.lak.stage.get(masked=False)
    for lake, records in DRY[label].items():
        np.testing.assert_array_equal(np.flatnonzero(raw[:, lake] == -1e30), records)
    observed = pd.read_csv(run.workspace / "lake_stage_targets.csv", float_precision="round_trip")
    assert len(observed) == NPER
    return raw, observed


def _file_route(label):
    """The same lakes under names no OBS declares: the stage-file route."""

    return {f"f_{name}": lake for name, lake in NAMES[label].items()}


@pytest.mark.slow
@pytest.mark.parametrize("label", LABELS)
def test_lake_bottoms_is_mf6s_own_stage_at_every_dry_record(runs, label):
    """MF6 writes its lake bottom nowhere but here: the STAGE observation of a dry
    lake. The rule is replicated -- VERTICAL takes the cell's top (the BELEV written,
    0.0 here, is ignored), HORIZONTAL its BELEV or, where BELEV equals TELEV, the
    cell's bottom, the lowest over the connections -- so it is held to MF6's number
    on every connection type, bit for bit."""

    raw, observed = _files(label, runs[label])
    bottoms = lake_bottoms(_loaded(runs[label]))
    assert bottoms.to_dict() == BOTTOMS[label]
    for name, lake in NAMES[label].items():
        for record in DRY[label].get(lake, []):
            assert observed[name.upper()].iloc[record] == bottoms[lake], (name, record)
    if label == "horizontal":  # what was WRITTEN for heq, which MF6 does not use
        written = pd.DataFrame(runs[label].model().fl.lak.connectiondata.get_data())
        assert list(written["belev"]) == [45.0, 55.0] and list(written["telev"]) == [58.0, 55.0]
    options = [line.strip() for line in (runs[label].workspace / "m.lak").read_text().splitlines()]
    assert ("IMPLICIT" in options) == (label == "implicit")


@pytest.mark.slow
@pytest.mark.parametrize("label", LABELS)
def test_both_routes_read_one_number_at_every_record(runs, label):
    """The stage-file route read ``NaN`` where the observation route read the
    bottom. They are one number now, dry or wet, record for record; an INACTIVE lake
    is ``NaN`` on both."""

    raw, _ = _files(label, runs[label])
    model = _loaded(runs[label])
    observed = mf.LakeStageTargets(locations=NAMES[label]).simulated_series(model)
    from_file = mf.LakeStageTargets(locations=_file_route(label)).simulated_series(model)

    assert list(observed["time"]) == list(from_file["time"]) == list(range(NPER))
    for name, lake in NAMES[label].items():
        np.testing.assert_array_equal(observed[name], from_file[f"f_{name}"])
        dry = raw[:, lake] == -1e30
        np.testing.assert_array_equal(from_file[f"f_{name}"][dry], BOTTOMS[label][lake])
        np.testing.assert_array_equal(from_file[f"f_{name}"].isna(), raw[:, lake] == 1e30)
    if label == "vertical":
        assert from_file["f_off"].isna().all()  # inactive in every record: its column stays
    if label == "horizontal":  # dry in every record: its column holds the bottom
        np.testing.assert_array_equal(from_file["f_heq"], 40.0)


@pytest.mark.slow
def test_a_model_still_in_memory_is_off_mf6s_bottom_by_flopys_rounding_only(runs):
    """``run.model()`` holds the elevations the model was BUILT with; MF6 read what
    FloPy wrote, at 8 decimals. So there a dry lake is the unrounded bottom -- within
    5e-9 of MF6's, never further (half of FloPy's last written digit)."""

    raw, observed = _files("vertical", runs["vertical"])
    built = runs["vertical"].model()
    assert lake_bottoms(built)[1] == ODD_TOP
    series = mf.LakeStageTargets(locations={"f_two": 1}).simulated_series(built)
    assert series["f_two"].iloc[3] == ODD_TOP
    assert 0.0 < abs(series["f_two"].iloc[3] - observed["TWO"].iloc[3]) <= 5e-9


@pytest.mark.slow
def test_compare_and_stats_count_a_dry_period_on_both_routes(runs):
    """A gauged lake that never ran dry, against a model that does in periods 1-3:
    the stage-file route dropped those periods (``n`` 3 where the observation route
    counted 6) -- one run, two statistics. Both compare the bottom now."""

    _files("vertical", runs["vertical"])
    model = _loaded(runs["vertical"])
    gauged = [100.45, 100.30, 100.20, 100.10, 100.35, 100.40]
    results = {}
    for name in ("one", "f_one"):
        targets = mf.LakeStageTargets(
            locations={name: 0}, values=pd.DataFrame({"time": range(NPER), name: gauged}))
        results[name] = (targets.compare(model).sort_values("per"), targets.stats(model).iloc[0])
    (obs_compare, obs_stats), (file_compare, file_stats) = results["one"], results["f_one"]

    np.testing.assert_array_equal(obs_compare["sim_stage"], file_compare["sim_stage"])
    np.testing.assert_array_equal(file_compare["residual"].iloc[[1, 2, 3]],
                                  100.0 - np.array(gauged)[[1, 2, 3]])
    assert obs_stats["n"] == file_stats["n"] == NPER
    assert obs_stats["rmse"] == file_stats["rmse"]


@pytest.mark.slow
@pytest.mark.parametrize("route", ["observation", "stage-file"])
def test_the_pest_files_hold_the_bottom_and_1e30_only_for_an_inactive_lake(runs, tmp_path, route):
    """The base run's CSV and the forward run's, on each route: a dry record is the
    bottom (it was 1e30 on the stage-file route, which stalled PESTPP-IES), and 1e30
    is left for the lake that has no stage at all -- ``off``, inactive throughout."""

    raw, observed = _files("vertical", runs["vertical"])
    model = _loaded(runs["vertical"])
    locations = NAMES["vertical"] if route == "observation" else _file_route("vertical")
    targets = mf.LakeStageTargets(locations=locations)
    project = SimpleNamespace(model=model, template_workspace=tmp_path,
                              pf=SimpleNamespace(add_observations=lambda *a, **k: None))
    prepared = prepare_lake_stage_observations(
        project, LakeStageObservationSpec(targets=targets, prefix="lk"))
    config = prepared["named_series_forward_run_config"]
    base = tmp_path / config["output_csv"]
    forward = tmp_path / "forward.csv"
    write_named_series_targets(
        "lake_stage", tmp_path / config["locations_file"], forward, sim_ws=model.workspace)

    assert forward.read_text() == base.read_text()
    written = pd.read_csv(base, float_precision="round_trip")
    prefix = "" if route == "observation" else "f_"
    for name, lake in NAMES["vertical"].items():
        column = written[prefix + name]
        if name == "off":
            np.testing.assert_array_equal(column, 1e30)
            continue
        np.testing.assert_array_equal(column, observed[name.upper()])  # MF6's own numbers
        np.testing.assert_array_equal(column[raw[:, lake] == -1e30], BOTTOMS["vertical"][lake])
    assert (written.drop(columns="time").to_numpy() == 1e30).sum() == NPER


def _embedded_model(tmp_path):
    """A lake of one VERTICAL connection beside one connected EMBEDDEDH, built and not run."""

    grid, props, edge, interior = _grid()
    lake = mf.lak.flopy(
        nlakes=2,
        packagedata=[[0, 100.5, 1, "plain"], [1, 100.5, 1, "embedded"]],
        connectiondata=[[0, 0, (0, interior[0]), "VERTICAL", 0.5, 0.0, 0.0, 0.0, 0.0],
                        [1, 0, (0, interior[5]), "EMBEDDEDH", 0.5, 0.0, 0.0, 10.0, 10.0]],
        perioddata={0: [[0, "WITHDRAWAL", 0.0]]},
        boundnames=True, stage_filerecord="m.lak.stg",
    )
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=100.0, botm=[0.0]),
        mf.npf(k=10.0), mf.ic(strt=50.0), lake])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=2, perioddata=[(10.0, 1, 1.0)] * 2)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        return mf.Project(tmp_path / "p", name="p").prepare_run("r", simulation).model()


def test_an_embedded_lake_is_refused_by_name_and_only_when_it_is_asked_for(tmp_path, monkeypatch):
    """MF6 takes an embedded lake's bottom from its lake TABLE's leading zero rows,
    which is not read: refused, naming the lake, the connection and the way round it
    -- only for the lakes asked about, so a dry plain lake beside it is read."""

    from myflopy.modflow.mf6.observations import lake as lake_module

    model = _embedded_model(tmp_path)
    with pytest.raises(ValueError, match=r"Lake 1 of model 'm' connects EMBEDDEDH: .*TABLE.*"
                                         r"STAGE observation"):
        lake_bottoms(model)
    assert lake_bottoms(model, lakes=[0]).to_dict() == {0: 100.0}
    with pytest.raises(ValueError, match=r"Lake\(s\) \[2\] are not lakes of model 'm'"):
        lake_bottoms(model, lakes=[2])

    # The table a target reads: MF6's -1e30 for whichever lake is dry, a stage for
    # the other. The bottoms are asked for the DRY lakes only.
    def raw_table(dry_lake):
        stage = {0: [101.0, 101.0], 1: [101.0, 101.0]}
        stage[dry_lake][1] = -1e30
        return pd.DataFrame({"per": [0, 1, 0, 1], "lake": [0, 0, 1, 1], "cell": [0, 0, 5, 5],
                             "stage": stage[0] + stage[1]})

    for dry_lake, refused in ((0, False), (1, True)):
        monkeypatch.setattr(lake_module, "build_lak_stage_result_table",
                            lambda model, *, masked=True, t=raw_table(dry_lake): t.copy())
        if refused:
            with pytest.raises(ValueError, match=r"Lake 1 of model 'm' connects EMBEDDEDH"):
                lake_module._calibration_stage_table(model)
        else:
            table = lake_module._calibration_stage_table(model)
            np.testing.assert_array_equal(table["stage"], [101.0, 100.0, 101.0, 101.0])


def test_lake_bottoms_reads_a_structured_and_an_unstructured_grid(tmp_path):
    """The rule by cellid on DIS ``(layer, row, column)`` and DISU ``(node,)`` too:
    a layer's cell top is the bottom of the layer above, and BELEV = TELEV takes the
    cell's bottom. Built, not run: nothing here reads an output."""

    top = np.full((4, 5), 100.0)
    top[2, 3] = 97.5
    botm = np.stack([np.full((4, 5), 60.0), np.full((4, 5), 40.0)])
    botm[0, 1, 1] = 62.5
    cases = {
        "dis": (mf.dis(nlay=2, nrow=4, ncol=5, delr=10.0, delc=10.0, top=top, botm=botm), [
            [0, 0, (0, 2, 3), "VERTICAL", 0.5, 0.0, 0.0, 0.0, 0.0],
            [0, 1, (0, 0, 0), "VERTICAL", 0.5, 0.0, 0.0, 0.0, 0.0],
            [1, 0, (1, 1, 1), "HORIZONTAL", 0.5, 55.0, 58.0, 10.0, 10.0],
            [2, 0, (1, 1, 1), "HORIZONTAL", 0.5, 50.0, 50.0, 10.0, 10.0]],
            {0: 97.5, 1: 55.0, 2: 40.0}),
        "disu": (mf.disu(
            nodes=3, nja=7, top=[100.0, 90.0, 80.0], bot=[50.0, 45.0, 20.0], area=100.0,
            iac=[2, 3, 2], ja=[0, 1, 1, 0, 2, 2, 1], ihc=[1] * 7,
            cl12=[0.0, 5.0, 0.0, 5.0, 5.0, 0.0, 5.0], hwva=[0.0, 10.0, 0.0, 10.0, 10.0, 0.0, 10.0]), [
            [0, 0, (1,), "VERTICAL", 0.5, 0.0, 0.0, 0.0, 0.0],
            [1, 0, (2,), "HORIZONTAL", 0.5, 30.0, 35.0, 10.0, 10.0],
            [2, 0, (2,), "HORIZONTAL", 0.5, 30.0, 30.0, 10.0, 10.0]],
            {0: 90.0, 1: 30.0, 2: 20.0}),
    }
    for label, (dis, connections, expected) in cases.items():
        counts = pd.Series([row[0] for row in connections]).value_counts()
        lake = mf.lak.flopy(
            nlakes=3, packagedata=[[i, 101.0, int(counts[i]), f"l{i}"] for i in range(3)],
            connectiondata=connections, perioddata={0: [[0, "WITHDRAWAL", 0.0]]}, boundnames=True)
        model = mf.gwf("m", packages=[dis, mf.npf(k=1.0), mf.ic(strt=90.0), lake])
        simulation = mf.simulation(model, name="s", packages=[
            mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            built = mf.Project(tmp_path / label, name="p").prepare_run("r", simulation).model()
        assert lake_bottoms(built).to_dict() == expected, label

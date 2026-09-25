"""A package declared against a context from BEFORE the cut names cells that are gone.

The excavation design's C5 (ledgers 321 and 322). `LayerStack.excavate` removes the cells
above a pond's bed from the DISV (idomain 0). Every package resolved from a domain that
predates the cut still addresses them, and nothing looked:

- an array recharge/ET whose ``irch``/``ievt`` names a removed cell loses that column's
  flux -- MF6 maps the index to node 0 and skips it WITHOUT a word;
- a list boundary or an observation on a removed cell is accepted by FloPy at
  construction (measured, a CHD at a dug cell built clean) and refused by MF6 at read.

``ModelSpec.build`` now refuses all three, naming the excavation when a record reaches the
model. ``irch=None`` (one layer for every column) is the documented opt-out, so the same
loss is WARNED, with the rate it costs -- and the slow test holds that rate against MF6's
own budget.
"""

from __future__ import annotations

import contextlib
import io
import pathlib
import re
import warnings
from types import SimpleNamespace

import flopy
import geopandas as gpd
import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy import ModelContext
from myflopy.grid_identity import active_cell_warnings, excavation_cells
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.observations.heads import HeadTargets

RATE = 1e-3
POND = (5, 6)


def _stack():
    return (mf.LayerStack(top=mf.Flat(100.0))
            .add("fill", thickness=5.0)
            .add("sand", thickness=30.0))


def _setup(grid, cells=POND):
    """The same stack built before and after a pond is cut through the fill at ``cells``."""

    before = _stack().build(grid)
    pond = mf.Area(grid.gdf_vorPolys.geometry.loc[list(cells)].union_all())
    after = _stack().excavate(bottom=93.0, within=pond, name="pond").build(grid, attach=True)
    assert (after.idomain[0, list(cells)] == 0).all(), "the premise: layer 1 is dug out"
    return SimpleNamespace(
        before=before, after=after,
        stale=ModelContext(grid=grid, domain=before.idomain),
        fresh=ModelContext(grid=grid, domain=after.idomain),
    )


@pytest.fixture(params=["rect", "triangulated"])
def grid(request, triangulated_grid):
    """A hand-built grid, and a real ``GridSpec.voronoi`` one (int64 cell2d, ledger 265)."""

    if request.param == "rect":
        return irregular_voronoi_grid(nrow=4, ncol=4, cell_size=50.0)
    return triangulated_grid


def _simulation(context, result, *packages, nper=1):
    model = mf.gwf("m", context=context, packages=[
        mf.disv(result), mf.npf(k=10.0), mf.ic(strt=99.0), *packages,
    ])
    return mf.simulation(model, name="s", packages=[
        mf.tdis(nper=nper, perioddata=[(1.0, 1, 1.0)] * nper),
    ])


def _build(tmp_path, spec, tag="ws"):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        built = spec.build_flopy(tmp_path / tag)
    return built, [str(w.message) for w in caught if "REMOVED" in str(w.message)]


def _pond_area(grid, cells=POND):
    return float(grid.gdf_vorPolys.geometry.loc[list(cells)].area.sum())


def _loss(message):
    """``(lost, total)`` as the message states them."""

    lost, total = re.search(r": (?:up to )?([0-9.e+-]+) of ([0-9.e+-]+) ft3/d", message).groups()
    return float(lost), float(total)


# -- refused: a per-column index, a list boundary, an observation --------------------


def test_a_stale_array_recharge_is_refused_naming_the_excavation_and_the_loss(tmp_path, grid):
    """The shared-library case: ONE `rch.array` declared against the pre-cut context
    builds beside the baseline and is refused beside the excavated DISV."""

    s = _setup(grid)
    rch = mf.rch.array(context=s.stale, recharge=RATE)
    _build(tmp_path, _simulation(s.stale, s.before, rch), "baseline")    # fine before the cut

    with pytest.raises(ValueError) as raised:
        _build(tmp_path, _simulation(s.fresh, s.after, rch), "excavated")
    message = str(raised.value)
    assert "package 'rch' (recharge as arrays)" in message
    assert "2 inside excavation 'pond'" in message
    assert "cell 5 (layer 1), cell 6 (layer 1)" in message
    assert "declared against a context from BEFORE the cut" in message
    assert "declare the package without context=" in message, "the eager package's own context"
    assert "up to" not in message, "a recharge loss is exact"
    lost, total = _loss(message)
    assert lost == pytest.approx(RATE * _pond_area(grid), rel=1e-5)
    assert total == pytest.approx(RATE * grid.gdf_vorPolys.area.sum(), rel=1e-5)


def test_a_lazy_array_resolved_against_a_stale_context_domain_is_refused(tmp_path, grid):
    """`irch='top_active'` without `context=` resolves at BUILD against the MODEL'S
    context -- which, declared before the cut, still says layer 1 is there. The refusal
    names THAT context: telling the user to declare the package without context= -- what
    they already did -- got them the same refusal again."""

    s = _setup(grid)
    with pytest.raises(ValueError, match=r"package 'rch' \(recharge as arrays\)") as raised:
        _build(tmp_path, _simulation(s.stale, s.after, mf.rch.array(recharge=RATE)))
    message = str(raised.value)
    assert "the model's context is from BEFORE the cut" in message
    assert "irch was left to the build" in message and "leave domain= out" in message
    assert "declare the package without context=" not in message
    # ...and the advice works: without domain= the build reads the DISV's own idomain
    built, warned = _build(tmp_path, _simulation(ModelContext(grid=grid), s.after,
                                                 mf.rch.array(recharge=RATE)), "fixed")
    assert warned == []
    irch = built.built_model("m").model.get_package("rcha").irch.get_data()[0]
    assert (irch[list(POND)] == 1).all()


def test_a_lazy_list_boundary_on_a_stale_model_context_names_that_context(tmp_path, grid):
    """A `.gpkg` boundary declared without `context=` resolves its cells against the
    MODEL's context domain too; re-declaring it would resolve them the same way."""

    s = _setup(grid)
    source = tmp_path / "chd.gpkg"
    gpd.GeoDataFrame({"head": [98.0], "layer": [1], "name": ["c"]},
                     geometry=[grid.gdf_vorPolys.geometry.loc[5].representative_point()],
                     crs=grid.crs).to_file(source)
    with pytest.raises(ValueError) as raised:
        _build(tmp_path, _simulation(s.stale, s.after, mf.chd.gpkg(source, head="head", nper=1)))
    message = str(raised.value)
    assert "package 'chd' names 1 cell(s) the built DISV has REMOVED" in message
    assert "the model's context is from BEFORE the cut" in message
    assert "leave domain= out" in message and "Re-declare it against" not in message


def test_a_bare_package_spec_is_checked_like_a_named_one(tmp_path, grid):
    """Paired with its spec by the OBJECT it built, not by name: FloPy names a
    `PackageSpec` without `pname` `chd_0`, and a name lookup skipped it -- a CHD on a dug
    cell built clean, and a stale per-column irch was warned as 'irch left out'."""

    s = _setup(grid)
    chd = mf.PackageSpec("chd", flopy.mf6.ModflowGwfchd,
                         {"stress_period_data": {0: [[(0, 5), 98.0]]}})
    with pytest.raises(ValueError, match=r"package 'chd' names 1 cell\(s\) the built DISV has REMOVED"):
        _build(tmp_path, _simulation(s.fresh, s.after, chd), "chd")
    rcha = mf.PackageSpec("rch", flopy.mf6.ModflowGwfrcha,
                          {"recharge": RATE, "irch": np.zeros(grid.ncpl, dtype=int)})
    with pytest.raises(ValueError, match=r"'rcha_0' \(recharge as arrays\) addresses its recharge"):
        _build(tmp_path, _simulation(s.fresh, s.after, rcha), "rcha")


def test_a_stale_array_evt_is_refused(tmp_path, grid):
    s = _setup(grid)
    evt = mf.evt.array(context=s.stale, rate=RATE, depth=5.0)
    with pytest.raises(ValueError, match=r"package 'evt' \(evapotranspiration as arrays\)"):
        _build(tmp_path, _simulation(s.fresh, s.after, evt))


def test_a_stale_typed_ievt_beside_an_omitted_surface_is_refused(tmp_path, grid):
    """A typed ``ievt`` array with the surface left to the build was DROPPED, so a stale
    one was warned as "ievt left out" -- downgraded from the refusal the DECLARED value
    calls for (ledger 324). It is carried now, and refused like any per-column index."""

    from myflopy.modflow.mf6.grid.geometry import top_active_layer

    s = _setup(grid)
    evt = mf.evt.array(rate=RATE, depth=5.0, ievt=top_active_layer(s.before.idomain))
    assert sorted(evt.inputs) == ["surface"]
    with pytest.raises(ValueError, match=r"package 'evt' \(evapotranspiration as arrays\) addresses"):
        _build(tmp_path, _simulation(s.fresh, s.after, evt))


def test_an_array_resolved_against_the_excavated_domain_builds_silently(tmp_path, grid):
    s = _setup(grid)
    for tag, rch in (("eager", mf.rch.array(context=s.fresh, recharge=RATE)),
                     ("lazy", mf.rch.array(recharge=RATE))):
        built, warned = _build(tmp_path, _simulation(s.fresh, s.after, rch), tag)
        assert warned == []
        irch = built.built_model("m").model.get_package("rcha").irch.get_data()[0]
        assert (irch[list(POND)] == 1).all(), "the connected cell under the bed"


def test_a_list_boundary_on_a_dug_cell_is_refused(tmp_path, grid):
    """FloPy accepted it at construction; MF6 refuses it at read. The guard says why."""

    s = _setup(grid)
    chd = mf.chd(stress_period_data={0: [[(1, 0), 98.0], [(0, 5), 98.0]]})
    with pytest.raises(ValueError) as raised:
        _build(tmp_path, _simulation(s.fresh, s.after, chd))
    message = str(raised.value)
    assert "package 'chd' names 1 cell(s) the built DISV has REMOVED" in message
    assert "(0, 5) -- 1 inside excavation 'pond'" in message


@pytest.mark.parametrize("where", ["package", "model"])
def test_an_observation_in_a_dug_cell_is_refused_naming_the_excavation(tmp_path, grid, where):
    s = _setup(grid)
    if where == "package":
        obs = mf.chd(stress_period_data={0: [[(1, 0), 98.0]]},
                     observations={"c.csv": [("c6", "chd", (0, 6))]})
        owner = "the observations of package 'chd'"
    else:
        obs = mf.PackageSpec("obs", flopy.mf6.ModflowUtlobs, {"continuous": {
            "h.csv": [("h5", "HEAD", (0, 5)), ("deep", "HEAD", (1, 5))]}})
        owner = "observation file"
    with pytest.raises(ValueError) as raised:
        _build(tmp_path, _simulation(s.fresh, s.after, obs))
    message = str(raised.value)
    assert owner in message and "inside excavation 'pond'" in message
    assert "'deep'" not in message, "the cell under the bed is active"


def test_a_removed_cell_with_no_excavation_record_is_refused_without_one(tmp_path, grid):
    """The guard is about ACTIVITY, not about ponds: a hand-edited idomain is caught too,
    and no excavation is blamed."""

    s = _setup(grid)
    idomain = np.array(s.before.idomain)
    idomain[0, 3] = 0
    disv = s.before.to_disv(grid).with_options(idomain=idomain)
    model = mf.gwf("m", context=ModelContext(grid=grid, domain=idomain), packages=[
        disv, mf.npf(k=10.0), mf.ic(strt=99.0),
        mf.chd(stress_period_data={0: [[(0, 3), 98.0]]}),
    ])
    spec = mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    grid.excavations = {}      # the fixture grid is shared; no record reaches this model
    with pytest.raises(ValueError) as raised:
        _build(tmp_path, spec)
    assert "(0, 3)" in str(raised.value) and "excavation" not in str(raised.value)


# -- warned: the documented one-layer opt-out ----------------------------------------


@pytest.mark.parametrize("make", [
    lambda: mf.rch.array(recharge=RATE, irch=None),
    lambda: mf.rch.array(recharge={0: RATE, 1: 2 * RATE}, irch=0),
    lambda: mf.evt.array(rate=RATE, depth=5.0, surface=100.0, ievt=None),
], ids=["irch=None", "irch=0 per period", "ievt=None"])
def test_one_layer_for_every_column_warns_with_the_rate_it_loses(tmp_path, grid, make):
    s = _setup(grid)
    built, warned = _build(tmp_path, _simulation(s.fresh, s.after, make(), nper=2))
    assert built is not None, "legal: warned, not refused"
    assert len(warned) == 1
    message = warned[0]
    assert "puts" in message and "one layer in every column" in message
    assert "inside excavation 'pond'" in message
    assert "BEFORE the cut" not in message, "a typed layer is not a stale context"
    lost, total = _loss(message)
    assert lost == pytest.approx(RATE * _pond_area(grid), rel=1e-5)
    assert total == pytest.approx(RATE * grid.gdf_vorPolys.area.sum(), rel=1e-5)
    # ET's figure is the MAXIMUM rate over the lost area -- MF6 scales ET by the head's
    # depth below the surface -- so it is stated as a bound; recharge's is exact.
    evt = "evapotranspiration" in message
    assert ("up to" in message) is evt
    assert ("at the maximum ET rate" in message) is evt


def test_the_worst_period_is_named_when_it_is_not_the_first(tmp_path, grid):
    s = _setup(grid)
    rch = mf.rch.array(recharge={0: RATE, 1: 3 * RATE}, irch=None)
    _, warned = _build(tmp_path, _simulation(s.fresh, s.after, rch, nper=2))
    assert "in stress period 1" in warned[0]
    worst = re.search(r"and ([0-9.e+-]+) \(([0-9.]+)%\) in stress period 2", warned[0])
    assert worst and float(worst.group(1)) == pytest.approx(3 * RATE * _pond_area(grid), rel=1e-5)


def test_the_opt_out_is_silent_where_layer_one_is_everywhere(tmp_path, grid):
    s = _setup(grid)
    _, warned = _build(tmp_path, _simulation(s.stale, s.before, mf.rch.array(recharge=RATE, irch=None)))
    assert warned == []


def test_a_barrier_touching_a_removed_cell_is_not_this_guards_to_refuse():
    """MF6 WARNS about an HFB touching an inactive cell and drops it (gwf-hfb.f90), so
    refusing one would refuse a model MF6 runs. FloPy 3.10 refuses the cellid itself at
    construction ("Cellid (0, 5) is outside of the active model grid"), so this is pinned
    on the rule, with the package TYPE deciding: the same records as a CHD are refused."""

    from myflopy.grid_identity import _inactive_cellid_complaints

    idomain = np.ones((2, 16), dtype=int)
    idomain[0, 5] = 0
    spec = mf.hfb.flopy(stress_period_data={0: [[(0, 4), (0, 5), 1e-3]]})

    def model(kind):
        return SimpleNamespace(packagelist=[SimpleNamespace(name=[spec.name], package_type=kind)])

    assert _inactive_cellid_complaints(model("hfb"), [spec], idomain, {}) == []
    [complaint] = _inactive_cellid_complaints(model("chd"), [spec], idomain, {})
    assert "(0, 5)" in complaint


def test_the_helpers_answer_nothing_without_a_disv():
    assert active_cell_warnings(object(), None, []) == []
    assert excavation_cells(None) == {}


# -- head targets --------------------------------------------------------------------


def _target_model(grid, result):
    disv = SimpleNamespace(idomain=SimpleNamespace(array=np.asarray(result.idomain)))
    return SimpleNamespace(vor=grid, gwf=SimpleNamespace(disv=disv))


def test_a_head_target_in_a_dug_cell_names_the_excavation(grid, caplog):
    s = _setup(grid)
    model = _target_model(grid, s.after)
    by_cell = HeadTargets(
        locations=pd.DataFrame({"name": ["W5"], "layer": [0], "cell": [5]}),
        values=pd.DataFrame({"per": [0], "W5": [1.0]}), time_column="per",
    )
    with pytest.raises(ValueError, match=r"cell 5, layer 0, inside excavation 'pond'.*under its bed"):
        by_cell.match_to_model(model)

    inside = grid.gdf_vorPolys.geometry.loc[5].representative_point()
    frame = gpd.GeoDataFrame({"name": ["W5", "W0"], "layer": [0, 0]},
                             geometry=[inside, grid.gdf_vorPolys.geometry.loc[0].representative_point()],
                             crs=grid.crs)
    by_point = HeadTargets(locations=frame, values=pd.DataFrame({"per": [0], "W5": [1.0], "W0": [1.0]}),
                           time_column="per")
    # A partial miss is reported through the LOGGER only -- ledger 288's warnings.warn
    # rule has not reached the target families (recorded in ledger 321), so read it there.
    with caplog.at_level("WARNING", logger="myflopy.modflow.mf6.observations.heads"):
        matched = by_point.match_to_model(model)
    assert matched["name"].tolist() == ["W0"]
    assert any("'W5' (on no active cell in layer 0, inside excavation 'pond'" in r.getMessage()
               for r in caplog.records)


# -- MF6 -----------------------------------------------------------------------------


@pytest.mark.slow
def test_the_warned_loss_is_what_mf6_drops(tmp_path, triangulated_grid):
    """Hold the warning's number against MF6's own budget: with ``irch=None`` over the
    dug pond, the RCH inflow MF6 reports is the declared total LESS what the warning
    said is lost -- under Normal termination, which is why it is warned at all. The
    budget is asserted too (the percent discrepancy), not Normal termination alone."""

    grid = triangulated_grid
    s = _setup(grid)
    east = [int(c) for c in np.flatnonzero(grid.gdf_vorPolys.geometry.bounds.maxx >= 600 - 1e-6)]
    model = mf.gwf("m", context=s.fresh, newtonoptions="NEWTON", packages=[
        mf.disv(s.after), mf.npf(k=10.0, k33=1.0, icelltype=1), mf.ic(strt=99.0),
        mf.chd(stress_period_data={0: [[(1, c), 98.0] for c in east]}),
        mf.rch.array(recharge=RATE, irch=None),
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    spec = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
        mf.ims(models=["m"], complexity="COMPLEX", outer_dvclose=1e-6, inner_dvclose=1e-8,
               rcloserecord=[1e-4, "STRICT"]),
    ])
    with warnings.catch_warnings(record=True) as caught, contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("always")
        run = mf.Project(tmp_path / "p", name="p").run("r", spec)
    assert run.success, "\n".join(run.report[-20:])
    [message] = [str(w.message) for w in caught if "one layer in every column" in str(w.message)]
    lost, total = _loss(message)

    workspace = pathlib.Path(run.workspace)
    listing = (workspace / "m.lst").read_text()
    discrepancy = float([line for line in listing.splitlines()
                         if "PERCENT DISCREPANCY" in line][-1].split()[-1])
    assert abs(discrepancy) < 0.1
    budget = flopy.utils.CellBudgetFile(str(workspace / "m.cbc"))
    applied = float(np.sum(budget.get_data(text="RCH")[-1]["q"]))
    assert applied == pytest.approx(total - lost, rel=1e-6)
    assert lost == pytest.approx(RATE * _pond_area(grid), rel=1e-5) and lost > 0

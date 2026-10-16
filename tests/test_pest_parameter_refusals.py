"""Parameter declarations that could not act, or would break a realization, are refused.

The leftovers of the PEST-preparation parameter lanes:

* a ``riv.stage`` offset moves each record's river bottom with its stage (ledger 519),
  and a lower bound below a record's clearance carried the bottom through its CELL's
  bottom, which MF6 refuses at ``icelltype != 0`` and over-drains at ``0`` (ledger 290)
  -- nothing looked (ledger 568);
* ``zones=`` was resolved on the CALLER's model object, not the copy of the run the
  build stages (ledger 356), and explicit ``pp_points=`` were located in its cells;
* ``zones=`` with ``style="constant"`` on a list target was accepted and did nothing:
  pyEMU reads it only for a ``"zone"`` or ``"grid"`` list parameter (ledger
  569);
* a correlated list ``"grid"`` whose records pyEMU places at one point either made its
  prior draw fail late (``Matrix is not positive definite``) or drew it uncorrelated
  -- refused for ``sfr.rhk`` alone (ledger 558; ledger 570);
* an ``**extra`` key the facade sets itself silently overrode it (ledger
  571), and an ``initial_value`` there could start a parameter outside its bounds
  (ledger 518's rule).

The models are FloPy-written and opened with ``mf.load_mf6_run``: the two-RIV DISV/DIS/
DISU model of ``tests/test_pest_several_packages.py`` (cell bottoms 50 and 0 ft), and
here a two-layer 4 x 5 model whose CHD records stack in a column.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import flopy
import numpy as np
import pytest
from flopy.utils.gridutil import get_disu_kwargs, get_disv_kwargs

import myflopy as mf
from myflopy.modflow.mf6.pest import native_parameters as np_mod
from tests.test_pest_several_packages import NCPL, _heads, _run

pyemu = pytest.importorskip("pyemu")

NLAY, NROW, NCOL = 2, 4, 5


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    root = tmp_path_factory.mktemp("parameter_refusals")
    return {grid: _run(root / grid, grid) for grid in ("disv", "disu")}


def _cal(run: Path, tmp_path: Path, name: str = "pr"):
    model = mf.load_mf6_run(run)
    cal = model.pest(name, workspace=tmp_path / name, start_datetime="2024-01-01")
    cal.observe(_heads())
    return cal


def _build(cal, name: str = "pr"):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return cal.build(f"{name}.pst")


# -- riv.stage: a lower bound may not carry a river bottom through its cell -------------


def test_a_stage_lower_bound_past_a_records_clearance_is_refused_naming_the_safe_one(
        runs, tmp_path):
    """``crisp``'s record sits at rbot 96 over a 50 ft cell bottom, ``riv``'s at 98 and
    98.5: an offset of -47 carries crisp's bottom to 49 ft, inside the cell below. The
    refusal names the count, the worst record and the most negative safe lower bound --
    46 ft less the clearance that survives FloPy's writing (ledger 290), rounded toward
    zero -- and that bound builds."""

    cal = _cal(runs["disv"], tmp_path, "deep")
    cal.parameterize("riv.stage", package=("riv", "crisp"), bounds=(-47.0, 2.0))
    with pytest.raises(ValueError, match=(
            r"lower bound -47 carries the river bottom of 1 record\(s\) below their cell's "
            r"bottom -- worst, row 1 of m\.crisp_stress_period_data_1\.txt \(package "
            r"'crisp'\), cellid \(0, 12\) \(zero-based\): rbot 96 over a cell bottom of 50.*"
            r"bounds=\(-45\.999, 2\)")):
        _build(cal, "deep")

    cal = _cal(runs["disv"], tmp_path, "safe")
    cal.parameterize("riv.stage", package=("riv", "crisp"), bounds=(-45.999, 2.0))
    _build(cal, "safe")


def test_the_lower_bounds_of_stage_offsets_on_one_file_add(runs, tmp_path):
    """Two ``riv.stage`` parameters on one file are two offsets the forward run ADDS:
    -30 and -20 carry crisp's bottom 50 ft down, through its cell, though each alone
    would not."""

    cal = _cal(runs["disv"], tmp_path, "stacked")
    cal.parameterize("riv.stage", package="crisp", bounds=(-30.0, 2.0), name="broad")
    cal.parameterize("riv.stage", package="crisp", bounds=(-20.0, 2.0), name="fine",
                     style="grid")
    with pytest.raises(ValueError, match=(
            r"'fine': its lower bound -20 carries .* with the -30 the riv\.stage parameters "
            r"already on its file add.*bounds=\(-15\.999, 2\)")):
        _build(cal, "stacked")


def test_the_cell_bottom_is_the_records_own_on_disu(runs, tmp_path):
    """On DISU a record names its NODE: crisp's is node 12 of the upper layer, bottom 50."""

    cal = _cal(runs["disu"], tmp_path, "node")
    cal.parameterize("riv.stage", package="crisp", bounds=(-47.0, 2.0))
    with pytest.raises(ValueError, match=r"cellid \(12,\) \(zero-based\): rbot 96 over a cell "
                                         r"bottom of 50"):
        _build(cal, "node")


# -- zones= is shaped on the staged model ------------------------------------------------


def test_zones_are_resolved_on_the_model_the_build_staged(runs, tmp_path, monkeypatch):
    """Everything the build reads off a FloPy model comes from the copy of the run it
    staged (ledger 356); ``zones=`` was shaped on the caller's object."""

    seen = []
    real = np_mod.resolve_zone_array

    def spy(zones, *, family, model):
        seen.append(model)
        return real(zones, family=family, model=model)

    monkeypatch.setattr(np_mod, "resolve_zone_array", spy)
    cal = _cal(runs["disv"], tmp_path, "zoned")
    zones = np.ones(NCPL, dtype=int)
    zones[12] = 2
    cal.parameterize("riv.cond", style="zone", zones=zones, package="crisp", bounds=(0.1, 10.0))
    cal.parameterize("k", style="zone", zones=zones, bounds=(0.5, 2.0))
    staged = []
    real_staged = np_mod.staged_model
    monkeypatch.setattr(np_mod, "staged_model",
                        lambda project, name: staged.append(real_staged(project, name))
                        or staged[-1])
    _build(cal, "zoned")
    assert len(seen) == 2
    for model in seen:
        assert isinstance(model, flopy.mf6.MFModel), type(model)
        assert any(model is copy for copy in staged), "zones= shaped on the caller's model"


@pytest.mark.parametrize("grid", ["disv", "disu"])
def test_explicit_pilot_points_are_located_on_the_model_the_build_staged(
        runs, tmp_path, monkeypatch, grid):
    """``pp_points=`` were located in the cell polygons of the CALLER's model object --
    the last grid the pilot-point build still read there -- while the net, the
    interpolation and ``zones=`` read the copy of the run it staged (ledger 356)."""

    from myflopy.modflow.mf6.pest import pilot_points as pp_mod

    located = []
    real = pp_mod.model_cell_polygons
    monkeypatch.setattr(pp_mod, "model_cell_polygons",
                        lambda model: located.append(model) or real(model))
    staged = []
    real_staged = pp_mod.staged_model
    monkeypatch.setattr(pp_mod, "staged_model",
                        lambda project, name: staged.append(real_staged(project, name))
                        or staged[-1])
    cal = _cal(runs[grid], tmp_path, f"pp{grid}")
    modelgrid = cal.model.gwf.modelgrid
    cells = [0, 7]
    points = np.c_[np.asarray(modelgrid.xcellcenters).reshape(-1)[cells],
                   np.asarray(modelgrid.ycellcenters).reshape(-1)[cells]]
    spec = cal.parameterize("k", style="pilotpoints", pp_points=points, bounds=(0.5, 2.0))
    _build(cal, f"pp{grid}")
    assert located, "the explicit points were never located"
    for model in located:
        assert isinstance(model, flopy.mf6.MFModel), type(model)
        assert any(model is copy for copy in staged), "pp_points= located on the caller's model"
    assert [item["points"] for item in spec.pilot_point_sets] == [2] * len(spec.pilot_point_sets)


# -- zones= with style='constant' on a list target ----------------------------------------


def test_zones_with_a_list_constant_are_refused_naming_zone(runs, tmp_path):
    """pyEMU scaled every CHD record with one number whatever the labels (ledger 559)."""

    zones = np.repeat([1, 2], NCPL // 2)
    with pytest.raises(ValueError, match=(
            r"zones= with style='constant' on list target 'riv\.cond'.*one multiplier would "
            r"scale every record.*Use style='zone'")):
        np_mod.NativeParameterSpec(target="riv.cond", style="constant", zones=zones)
    with pytest.raises(ValueError, match=r"one offset would scale every record"):
        np_mod.NativeParameterSpec(target="chd", style="constant", zones=zones)
    # 'zone' and 'grid' read the labels; an ARRAY constant masks by them.
    np_mod.NativeParameterSpec(target="riv.cond", style="zone", zones=zones)
    np_mod.NativeParameterSpec(target="riv.cond", style="grid", zones=zones)
    np_mod.NativeParameterSpec(target="k", style="constant", zones=zones)


def test_recharge_is_asked_once_its_form_is_known(runs, tmp_path):
    """``recharge`` is a list target on a list RCH (the DISU run) and an array target on
    an array RCH (the DISV run), where a constant masks by its zones."""

    zones = np.repeat([1, 2], NCPL // 2)
    cal = _cal(runs["disu"], tmp_path, "listrch")
    with pytest.raises(ValueError, match=r"zones= with style='constant' on list target 'recharge'"):
        cal.parameterize("recharge", style="constant", zones=np.tile(zones, 2))
    cal = _cal(runs["disv"], tmp_path, "arrayrch")
    cal.parameterize("recharge", style="constant", zones=zones)


# -- a correlated list grid pyEMU cannot correlate ----------------------------------------


def _stacked_run(workspace: Path, grid: str, chd: list[tuple[int, int]],
                 drn: tuple[tuple[int, int], ...] = ()) -> Path:
    """A 2 x 4 x 5 model whose CHD records sit at ``chd`` (layer, cell-of-layer) pairs,
    and DRN records at ``drn`` (MF6 refuses two CHD in one cell, not two drains)."""

    sim = flopy.mf6.MFSimulation("s", sim_ws=str(workspace), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=2, perioddata=[(1.0, 1, 1.0)] * 2)
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m", save_flows=True)
    delr, delc = np.full(NCOL, 10.0), np.full(NROW, 10.0)
    if grid == "disv":
        flopy.mf6.ModflowGwfdisv(gwf, **get_disv_kwargs(NLAY, NROW, NCOL, delr, delc, 100.0,
                                                        [50.0, 0.0]))

        def cell(layer, number):
            return (layer, number)
    elif grid == "dis":
        flopy.mf6.ModflowGwfdis(gwf, nlay=NLAY, nrow=NROW, ncol=NCOL, delr=10.0, delc=10.0,
                                top=100.0, botm=[50.0, 0.0])

        def cell(layer, number):
            return (layer, number // NCOL, number % NCOL)
    else:
        flopy.mf6.ModflowGwfdisu(gwf, **get_disu_kwargs(NLAY, NROW, NCOL, delr, delc, 100.0,
                                                        [50.0, 0.0], return_vertices=True))

        def cell(layer, number):
            return (layer * NCPL + number,)
    flopy.mf6.ModflowGwfnpf(gwf, k=5.0, icelltype=0)
    flopy.mf6.ModflowGwfic(gwf, strt=100.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[
        [cell(layer, number), 100.0 - 0.1 * number] for layer, number in chd])
    if drn:
        flopy.mf6.ModflowGwfdrn(gwf, stress_period_data=[
            [cell(layer, number), 96.0, 1.0 + index] for index, (layer, number) in enumerate(drn)])
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord="m.hds", saverecord=[("HEAD", "ALL")])
    sim.write_simulation(silent=True)
    assert sim.run_simulation(silent=True)[0]
    return workspace


#: Three columns, each with a CHD record in both layers.
COLUMNS = [(layer, number) for number in (0, 1, 2) for layer in (0, 1)]


@pytest.fixture(scope="module")
def stacked(tmp_path_factory):
    root = tmp_path_factory.mktemp("stacked_chd")
    return {grid: _stacked_run(root / grid, grid, COLUMNS) for grid in ("disv", "dis", "disu")}


@pytest.mark.parametrize("grid, together", [("disv", r"\[\(0, 0\), \(1, 0\)\]"),
                                            ("disu", r"\[\(0,\), \(20,\)\]")])
def test_a_columns_records_on_disv_or_disu_are_refused_as_a_singular_prior(
        stacked, tmp_path, grid, together):
    """pyEMU keeps a DISV/DISU set in ONE covariance group whatever lies at one point, so
    a CHD in both layers of a column has correlation 1 with itself: the build passed and
    the prior draw died with ``LinAlgError: Matrix is not positive definite`` (measured,
    ``docs/reproducers/correlated_list_grid_shared_points.py``)."""

    cal = _cal(stacked[grid], tmp_path, "col")
    cal.parameterize("chd", style="grid", correlation=30.0)
    with pytest.raises(ValueError, match=(
            rf"records {together} \(zero-based\) lie at one point -- one model column -- in "
            r"one covariance group \(6 selected records.*not positive definite.*Drop "
            r"correlation=")):
        _build(cal, "col")


def test_a_nugget_does_not_let_a_columns_records_through_and_the_refusal_says_why(
        stacked, tmp_path):
    """A nugget makes the stacked set drawable -- measured from ``1e-10`` up -- but
    correlates a column's records at ``1 / (1 + nugget)``, a structure along the column
    the map-plane variogram never described: refused whatever the nugget, as ledger 452
    refuses a DISU array's stacked nodes. The refusal said the covariance was singular,
    which a nugget makes false (review of ledger 570)."""

    cal = _cal(stacked["disv"], tmp_path, "nug")
    cal.parameterize("chd", style="grid", correlation=30.0, nugget=0.25)
    with pytest.raises(ValueError, match=(
            r"lie at one point .* so nugget=0\.25 correlates them at 0\.8 -- drawable, but "
            r"a structure along the column .*ledger 452")) as caught:
        _build(cal, "nug")
    assert "singular" not in str(caught.value)


def test_a_columns_records_on_dis_are_correlated_within_each_layer(stacked, tmp_path):
    """On DIS pyEMU splits a set whose records share a point by their LAYER: each layer's
    records are correlated with each other and not across layers. Let through."""

    cal = _cal(stacked["dis"], tmp_path, "dis")
    cal.parameterize("chd", style="grid", correlation=30.0)
    _build(cal, "dis")
    frame = cal._native_parameter_frames["chd"]
    assert sorted(set(frame["covgp"])) == ["chd_cov0", "chd_cov1"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        np.random.seed(7)
        drawn = cal.draw_prior(400)._df
    names = list(frame["parnme"])
    corr = drawn[names].corr()
    layer = {name: int(name.split("idx0:")[1].split("_")[0]) for name in names}
    within = [corr.loc[a, b] for a in names for b in names if a < b and layer[a] == layer[b]]
    across = [corr.loc[a, b] for a in names for b in names if layer[a] != layer[b]]
    assert min(within) > 0.4 and max(abs(value) for value in across) < 0.3


def test_a_dis_split_into_one_record_per_group_is_refused(tmp_path):
    """On DIS, one record per layer of one column, and nothing else, is split one group
    per record: an uncorrelated prior."""

    one_each = _stacked_run(tmp_path / "run_one_each", "dis", [(0, 0), (1, 0)])
    cal = _cal(one_each, tmp_path, "one")
    cal.parameterize("chd", style="grid", correlation=30.0)
    with pytest.raises(ValueError, match=r"draws EVERY one's prior uncorrelated"):
        _build(cal, "one")

    # Two drains in one cell and layer are ONE parameter (pyEMU's rule: one per distinct
    # cellid), not two points; the column's other layer is a group with a neighbour.
    twice = _stacked_run(tmp_path / "run_twice", "dis", [(0, 0)],
                         drn=((0, 7), (0, 7), (1, 7), (0, 9), (1, 9)))
    cal = _cal(twice, tmp_path, "two")
    cal.parameterize("drn.cond", style="grid", correlation=30.0)
    _build(cal, "two")
    assert len(cal._native_parameter_frames["drncond"]) == 4


def test_records_apart_still_build_and_draw_correlated(tmp_path):
    apart = _stacked_run(tmp_path / "run_apart", "disv", [(0, 0), (0, 1), (0, 2)])
    cal = _cal(apart, tmp_path, "apart")
    cal.parameterize("chd", style="grid", correlation=30.0)
    _build(cal, "apart")


# -- **extra may not override what the facade sets ----------------------------------------


@pytest.mark.parametrize("key, through", [
    # Set by features=/style='feature' alone: package= picks FILES, every row of them.
    ("use_rows", "features= or style='feature' on lak.bedleak/sfr.rhk"),
    ("index_cols", "the target"),
    ("use_cols", "the target"),
    ("par_type", "style="),
    ("pargp", "name="),
    ("lower_bound", "bounds="),
    ("ult_lbound", "physical="),
    ("par_style", "additive="),
    ("zone_array", "zones="),
    ("geostruct", "correlation="),
])
def test_an_extra_keyword_the_facade_sets_is_refused_naming_its_argument(key, through):
    with pytest.raises(TypeError, match=rf"sets {key}= through {through}"):
        np_mod.NativeParameterSpec(target="riv.cond", extra={key: [0]})


def test_an_extra_keyword_the_facade_leaves_alone_still_passes(runs, tmp_path):
    cal = _cal(runs["disv"], tmp_path, "fmt")
    spec = cal.parameterize("riv.cond", package="crisp", bounds=(0.1, 10.0),
                            mfile_fmt="free")
    assert spec.extra == {"mfile_fmt": "free"}
    with pytest.raises(TypeError, match=r"sets pargp= through name="):
        cal.parameterize("riv.cond", package="riv", pargp="mine")


@pytest.mark.parametrize("target, bounds, value", [
    ("riv.cond", None, 3.0),          # a multiplier, default (0.5, 2.0)
    ("riv.cond", (0.5, 2.0), 0.25),
    ("riv.stage", None, -2.5),        # an offset, default (-2.0, 2.0)
    ("k", (0.5, 2.0), 2.000001),
])
def test_an_extra_initial_value_outside_the_bounds_is_refused(target, bounds, value):
    """pyEMU writes ``initial_value`` as every parameter's start, so one outside the
    bounds began PEST outside them -- the start ledger 518 refuses when the bounds
    exclude the facade's own -- and it passed as a keyword the facade leaves alone
    (ledger 571)."""

    with pytest.raises(ValueError, match=r"initial_value=.* outside its bounds=.*ledger 518"):
        np_mod.NativeParameterSpec(target=target, bounds=bounds,
                                   extra={"initial_value": value})


@pytest.mark.parametrize("value", [float("nan"), float("inf"), "high", True, [1.0]])
def test_an_extra_initial_value_that_is_not_one_number_is_refused(value):
    with pytest.raises(ValueError, match=r"initial_value=.*give one finite number"):
        np_mod.NativeParameterSpec(target="riv.cond", extra={"initial_value": value})


def test_an_extra_initial_value_inside_the_bounds_still_reaches_pyemu(runs, tmp_path):
    cal = _cal(runs["disv"], tmp_path, "start")
    for value in (None, 0.5, 1.5, 2.0):
        spec = np_mod.NativeParameterSpec(target="riv.cond", bounds=(0.5, 2.0),
                                          extra={"initial_value": value})
        assert spec.extra == {"initial_value": value}
    spec = cal.parameterize("riv.cond", package="crisp", bounds=(0.5, 2.0),
                            initial_value=1.5)
    pst = _build(cal, "start")
    names = [name for name in pst.parameter_data.index if name.startswith("pname:rivcond")]
    assert names and (pst.parameter_data.loc[names, "parval1"] == 1.5).all()
    assert spec.extra == {"initial_value": 1.5}

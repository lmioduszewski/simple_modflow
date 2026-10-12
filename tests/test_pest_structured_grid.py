"""PEST on a structured (DIS) model (ledger 420).

Every array parameter on a DIS model failed ``cal.build()`` -- ``style="constant"``
on ``k`` included -- in pyEMU's ``write_array_tpl`` (``IndexError: index 5 is out of
bounds for axis 0 with size 5``): the build wrote each externalized array one value
per line, and pyEMU, handed the model's ``StructuredGrid``, reads an array file as
``(nrow, ncol)``. A list parameter built, onto the wrong column: the recipes name
DISV's columns (``layer cell value``), and on DIS (``layer row col value``) a CHD
multiplier scaled each record's COLUMN NUMBER. And a target located on a model
built from ``mf.dis(...)`` -- whose ``model.vor`` is ``None`` -- died with
``AttributeError: 'NoneType' object has no attribute 'gdf_vorPolys'``.

The model here is FloPy-written and opened with ``mf.load_mf6_run``, as a DIS model
from outside myflopy would be; the ``mf.dis`` route is the second fixture.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
import warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest

import myflopy as mf
from myflopy.modflow.mf6.grid.geometry import model_cell_polygons
from myflopy.modflow.mf6.pest.ies import IesResults
from myflopy.modflow.mf6.pest.native_parameters import _RECIPES, list_columns

pyemu = pytest.importorskip("pyemu")
flopy = pytest.importorskip("flopy")

NLAY, NROW, NCOL = 2, 5, 6
K = np.linspace(1.0, 5.0, NLAY * NROW * NCOL).reshape(NLAY, NROW, NCOL)
ZONES = np.repeat([1, 2], NROW * NCOL // 2)


def _dis_run(workspace: Path) -> Path:
    sim = flopy.mf6.MFSimulation("s", sim_ws=str(workspace), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=3, perioddata=[(1.0, 1, 1.0)] * 3)
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m", save_flows=True)
    flopy.mf6.ModflowGwfdis(gwf, nlay=NLAY, nrow=NROW, ncol=NCOL, delr=10.0, delc=10.0,
                            top=100.0, botm=[50.0, 0.0])
    flopy.mf6.ModflowGwfnpf(gwf, k=K, k33=0.5, icelltype=0)
    flopy.mf6.ModflowGwfic(gwf, strt=100.0)
    flopy.mf6.ModflowGwfsto(gwf, steady_state={0: True}, transient={1: True}, ss=1e-5, sy=0.1)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=(
        [[(0, row, 0), 100.0] for row in range(NROW)]
        + [[(0, row, NCOL - 1), 95.0] for row in range(NROW)]))
    flopy.mf6.ModflowGwfrcha(gwf, recharge=1e-3)
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord="m.hds", budget_filerecord="m.cbc",
                           saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])
    sim.write_simulation(silent=True)
    success, _ = sim.run_simulation(silent=True)
    assert success
    return workspace


def _heads():
    points = gpd.GeoDataFrame({"name": ["W", "V"], "layer": [0, 1]},
                              geometry=gpd.points_from_xy([25.0, 45.0], [25.0, 15.0]))
    return mf.HeadTargets(locations=points, values=pd.DataFrame(
        {"time": [0, 1, 2], "W": [98.0, 98.0, 98.0], "V": [97.0, 97.0, 97.0]}))


def _build(cal, name: str):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return cal.build(f"{name}.pst", noptmax=0)


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    root = tmp_path_factory.mktemp("dis_pest")
    model = mf.load_mf6_run(_dis_run(root / "run"))
    cal = model.pest("dis", workspace=root / "pest", start_datetime="2024-01-01")
    cal.parameterize("k", style="grid", correlation=30.0, bounds=(0.5, 2.0))
    cal.parameterize("k33", style="zone", zones=ZONES, bounds=(0.5, 2.0))
    cal.parameterize("chd", style="constant", bounds=(0.9, 1.1))
    cal.observe(_heads())
    pst = _build(cal, "dis")
    return model, cal, pst


def test_list_columns_follow_the_cellid_width():
    """DISV is what the recipes are written for; DIS's cellid is one column wider,
    DISU's one narrower, and UZF's leading ``iuzno`` stays where it is."""

    chd, uzf = _RECIPES["chd"], _RECIPES["uzf.vks"]
    assert list_columns(chd, "vertex") == ([0, 1], 2)
    assert list_columns(chd, "structured") == ([0, 1, 2], 3)
    assert list_columns(chd, "unstructured") == ([0], 1)
    assert list_columns(uzf, "vertex") == ([1, 2], 6)
    assert list_columns(uzf, "structured") == ([1, 2, 3], 7)


def test_array_and_list_parameters_on_a_dis_model_scale_what_mf6_reads(built, tmp_path):
    """The build writes each array as ``(nrow, ncol)``; the forward run's own
    parameter application then scales every K cell and every CHD HEAD -- never a
    cellid -- and the file stays in the layout pyEMU and MF6 both read. Run in a copy
    of the template, as a PEST worker does."""

    _model, cal, pst = built
    template = tmp_path / "worker"
    shutil.copytree(cal.template_workspace, template)
    k_file = template / "m.npf_k_layer1.txt"
    assert np.loadtxt(k_file).shape == (NROW, NCOL)
    assert np.allclose(np.loadtxt(k_file), K[0])

    par = pst.parameter_data.copy()
    par.loc[par.index[par["pargp"].str.startswith("kl")], "parval1"] = 2.0
    par.loc[par.index[par["pargp"] == "chd"], "parval1"] = 1.1
    zone2 = [name for name in par.index if name.startswith("pname:k33") and name.endswith("zone:2")]
    par.loc[zone2, "parval1"] = 0.5
    for tpl, inp in zip(pst.template_files, pst.input_files, strict=True):
        pyemu.pst_utils.write_to_template(par["parval1"], str(template / tpl), str(template / inp))
    subprocess.run([sys.executable, "forward_run.py"], cwd=template, check=True,
                   capture_output=True, text=True)

    assert np.allclose(np.loadtxt(template / "m.npf_k_layer1.txt"), 2.0 * K[0])
    assert np.allclose(np.loadtxt(template / "m.npf_k_layer2.txt"), 2.0 * K[1])
    k33 = np.loadtxt(template / "m.npf_k33_layer1.txt").reshape(-1)
    assert np.allclose(k33, np.where(ZONES == 2, 0.25, 0.5))
    chd = pd.read_csv(template / "m.chd_stress_period_data_1.txt", sep=r"\s+", header=None)
    assert chd.iloc[:, :3].to_numpy().tolist() == (
        [[1, row + 1, 1] for row in range(NROW)] + [[1, row + 1, NCOL] for row in range(NROW)])
    assert np.allclose(chd.iloc[:, 3], [110.0] * NROW + [104.5] * NROW)


def test_head_targets_on_a_dis_model_are_its_cells(built):
    """A point is located on the DIS cell holding it, numbered ``row * ncol + col``."""

    model, _cal, _pst = built
    matched = _heads().match_to_model(model).set_index("name")
    assert matched.loc["W", "cell"] == 2 * NCOL + 2 and matched.loc["V", "cell"] == 3 * NCOL + 4


@pytest.mark.slow
def test_pestpp_ies_runs_a_dis_template(built):
    """PESTPP-IES's prior (NOPTMAX=-1) runs the DIS template; its BASE realization
    reads the heads the build read from the run."""

    _model, cal, _pst = built
    template = cal.template_workspace
    expected = {}
    for ins in sorted(template.glob("*.csv.ins")):
        read = pyemu.pst_utils.InstructionFile(str(ins)).read_output_file(
            str(template / ins.name[: -len(".ins")]))
        expected.update(read["obsval"].astype(float).to_dict())
    assert len(expected) == 6
    results = cal.prior(reals=3)
    obs = pd.read_csv(Path(results.workspace) / "dis.0.obs.csv", index_col=0)
    assert len(obs) == 3
    base = obs.loc[[label for label in obs.index if str(label) == "base"]].iloc[0]
    names = sorted(expected)
    assert np.allclose(base[names].to_numpy(dtype=float), [expected[name] for name in names],
                       rtol=1e-5, atol=0.0)


def test_a_model_built_from_mf_dis_locates_targets_and_pilot_points(tmp_path):
    """``mf.dis(...)`` leaves ``model.vor`` ``None``: its cells are FloPy's grid's,
    and head targets, list zones and pilot points are placed on them."""

    project = mf.Project(tmp_path / "p", name="p")
    flow = mf.gwf("m", packages=[
        mf.dis(nlay=NLAY, nrow=NROW, ncol=NCOL, delr=10.0, delc=10.0, top=100.0,
               botm=[50.0, 0.0]),
        mf.ic(strt=100.0), mf.npf(k=1.0, k33=0.5),
        mf.chd(stress_period_data=[[(0, row, 0), 100.0] for row in range(NROW)]
               + [[(0, row, NCOL - 1), 95.0] for row in range(NROW)]),
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    project.add_simulation(mf.simulation(flow, name="s", packages=[
        mf.tdis(nper=2, perioddata=[(1.0, 1, 1.0)] * 2)]))
    run = project.prepare_run("r", "s")
    run.execute()
    model = run.model("m")
    assert model.vor is None
    polygons = model_cell_polygons(model)
    assert len(polygons) == NROW * NCOL
    assert polygons.geometry.iloc[NCOL + 1].centroid.coords[0] == pytest.approx((15.0, 35.0))

    cal = model.pest("dis", workspace=tmp_path / "pest", start_datetime="2024-01-01")
    cal.parameterize("k", style="pilotpoints", pp_space=2, bounds=(0.5, 2.0))
    cal.parameterize("chd", style="zone", zones=ZONES, bounds=(0.9, 1.1))
    points = gpd.GeoDataFrame({"name": ["W"], "layer": [0]},
                              geometry=gpd.points_from_xy([25.0], [25.0]))
    cal.observe(mf.HeadTargets(locations=points,
                               values=pd.DataFrame({"time": [0, 1], "W": [98.0, 98.0]})))
    pst = _build(cal, "dis")
    names = pst.parameter_data.index
    assert {name for name in names if name.startswith("pname:chd")} == {
        "pname:chd_inst:0_ptype:zn_usecol:3_pstyle:m_zone:1",
        "pname:chd_inst:0_ptype:zn_usecol:3_pstyle:m_zone:2"}
    assert sum(name.startswith("kl0pp") for name in names) > 1
    assert list(pst.observation_data.index) == [
        "oname:hds_otype:lst_usecol:w_per:0", "oname:hds_otype:lst_usecol:w_per:1"]


def test_spectral_simulation_draws_on_a_dis_model(built):
    """Ledger 384's regular-grid branch, reachable now that a DIS model builds: pyEMU's
    spectral simulation draws every grid-style K parameter of the template."""

    _model, cal, pst = built
    ensemble = cal.draw_prior(4, use_specsim=True)
    grid = [name for name in pst.parameter_data.index if "ptype:gr" in name]
    assert len(grid) == NLAY * NROW * NCOL
    assert ensemble.shape[0] == 4 and set(grid) <= set(ensemble.columns)
    assert ensemble.loc[:, grid].to_numpy().std() > 0


@pytest.mark.slow
@pytest.mark.parametrize("style", ["constant", "pilotpoints"])
def test_a_captured_dis_field_reviews_on_its_own_cells(tmp_path, style):
    """A captured array on DIS is written a grid row per line, so its observations are
    named ``i:<row>_j:<col>``; the review read ``i`` alone as the cell, putting every
    value of a row on cell ``row``. A pilot-point array is written back one value per
    line, so there the line IS the cell, and reading it as a grid row of ``ncol`` put a
    layer's 30 values on cells 0..174 (review of ledger 420). The layout is read
    from the observation names, so each value is on its own cell either way."""

    model = mf.load_mf6_run(_dis_run(tmp_path / "run"))
    cal = model.pest("cap", workspace=tmp_path / "pest", start_datetime="2024-01-01")
    extra = {"pp_space": 2} if style == "pilotpoints" else {}
    cal.parameterize("k", style=style, bounds=(0.5, 2.0), capture=True, **extra)
    cal.observe(_heads())
    _build(cal, "cap")
    results = cal.prior(reals=3)
    for layer in range(NLAY):
        field = results.field("k", layer=layer)
        assert field["cell"].tolist() == list(range(NROW * NCOL))
        # PESTPP-IES writes its ensembles to six significant digits.
        assert np.allclose(field["base"], K[layer].reshape(-1), rtol=1e-5)


@pytest.mark.slow
def test_the_review_maps_of_a_model_built_from_mf_dis_draw_its_cells(tmp_path):
    """``IesResults``' maps read ``model.vor``, which a model declared with
    ``mf.dis(...)`` does not have, so ``plot_field`` and ``plot_obs_residuals`` raised
    ``AttributeError: 'NoneType' object has no attribute 'ncpl'`` on both backends
    (ledger 451). They draw on FloPy's cells now: each captured K on its own cell,
    and the target on the point the same run opened with ``load_mf6_run`` places it."""

    project = mf.Project(tmp_path / "p", name="p")
    flow = mf.gwf("m", packages=[
        mf.dis(nlay=NLAY, nrow=NROW, ncol=NCOL, delr=10.0, delc=10.0, top=100.0,
               botm=[50.0, 0.0]),
        mf.ic(strt=100.0), mf.npf(k=K, k33=0.5),
        mf.chd(stress_period_data=[[(0, row, 0), 100.0] for row in range(NROW)]
               + [[(0, row, NCOL - 1), 95.0] for row in range(NROW)]),
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    project.add_simulation(mf.simulation(flow, name="s", packages=[
        mf.tdis(nper=2, perioddata=[(1.0, 1, 1.0)] * 2)]))
    run = project.prepare_run("r", "s")
    run.execute()
    model = run.model("m")
    assert model.vor is None
    cal = model.pest("dis", workspace=tmp_path / "pest", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), capture=True)
    points = gpd.GeoDataFrame({"name": ["W"], "layer": [0]},
                              geometry=gpd.points_from_xy([25.0], [25.0]))
    cal.observe(mf.HeadTargets(locations=points,
                               values=pd.DataFrame({"time": [0, 1], "W": [98.0, 98.0]})))
    _build(cal, "dis")
    results = cal.prior(reals=3)
    assert results.model is model

    for layer in range(NLAY):
        cells = results.plot_field("k", stat="base", layer=layer).data[0]
        assert list(cells.locations) == list(range(NROW * NCOL))
        # The map colours log10(K); PESTPP-IES writes six significant digits.
        assert np.allclose(10.0 ** np.asarray(cells.z, dtype=float), K[layer].reshape(-1),
                           rtol=1e-5)
    assert results.plot_field("k", stat="base", backend="mpl").axes

    reopened = mf.load_mf6_run(run.workspace)
    marker = results.plot_obs_residuals().data[-1]
    expected = IesResults(results.workspace, model=reopened).plot_obs_residuals().data[-1]
    assert np.allclose(marker.lon, expected.lon) and np.allclose(marker.lat, expected.lat)
    assert results.plot_obs_residuals(backend="mpl").axes


@pytest.mark.parametrize("nlay", [1, 2])
def test_an_array_target_the_package_does_not_declare_is_refused_by_name(tmp_path, nlay):
    """``parameterize("k33")`` on a model whose NPF declares no K33 died in the
    relayering with ``cannot reshape array of size 1 into shape (2,newaxis)`` on two
    layers and ``No external input file found`` on one, naming neither the target nor
    the cause; it is refused naming both, and what MF6 uses instead.
    A target the package declares still builds."""

    sim = flopy.mf6.MFSimulation("s", sim_ws=str(tmp_path / "run"), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=1)
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
    flopy.mf6.ModflowGwfdis(gwf, nlay=nlay, nrow=3, ncol=4, delr=10.0, delc=10.0,
                            top=100.0, botm=[50.0, 0.0][:nlay])
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0, icelltype=0)
    flopy.mf6.ModflowGwfic(gwf, strt=100.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 0), 100.0], [(0, 2, 3), 95.0]])
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord="m.hds", saverecord=[("HEAD", "ALL")])
    sim.write_simulation(silent=True)
    assert sim.run_simulation(silent=True)[0]
    model = mf.load_mf6_run(tmp_path / "run")

    cal = model.pest("v", workspace=tmp_path / "pest", start_datetime="2024-01-01")
    cal.parameterize("k33", style="constant", bounds=(0.5, 2.0))
    with pytest.raises(ValueError, match=r"(?s)Target 'k33'.*NPF of model 'm'.*does not "
                                         r"declare it.*K33 equal to K.*Declare `k33=`"):
        _build(cal, "v")

    cal = model.pest("k", workspace=tmp_path / "pest_k", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0))
    assert len(_build(cal, "k").parameter_data) == 1

"""A parameter finds its package by TYPE, and several of one type are never picked in silence.

Ledger 517. Every list recipe globbed ``{model}.<type>_stress_period_data_*.txt``
for its files, but FloPy names an external file after the package's FILE
(``m.crisp_stress_period_data_1.txt``), so a GHB named ``edge`` raised
``FileNotFoundError`` and a second DRN beside ``drn`` was left unscaled with no word.
The files now come from FloPy's own records of the staged model, by type
(``package_registry.packages_of_type``), and ``package=`` names the one -- or the
several sharing one set -- when the model carries more than one (ledger 352).

The models here are FloPy-written and opened with ``mf.load_mf6_run``: DISV with
two RIV (``riv`` and ``crisp``, as Ten Trails carries), two DRN (``drn`` and
``seeps``) and a lone GHB named ``edge``; and DIS / DISU variants for the column
rule (ledger 420) of ``riv.stage``'s pair.
"""

from __future__ import annotations

import json
import re
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
from myflopy.modflow.mf6.pest import native_parameters as np_mod

pyemu = pytest.importorskip("pyemu")
flopy = pytest.importorskip("flopy")
from flopy.utils.gridutil import get_disu_kwargs, get_disv_kwargs  # noqa: E402

NLAY, NROW, NCOL = 2, 4, 5
NCPL = NROW * NCOL
#: (layer, cell) -> (stage, cond, rbot) of each river package.
RIV = {(0, 6): (99.0, 10.0, 98.0), (0, 7): (99.0, 20.0, 98.5)}
CRISP = {(0, 12): (98.0, 5.0, 96.0)}


def _flow_packages(gwf, cell, *, array_recharge=True):
    """Every package but the discretization, on a cellid of the grid's form. MF6 reads
    array recharge on DIS and DISV only."""

    flopy.mf6.ModflowGwfnpf(gwf, k=5.0, k33=0.5, icelltype=0)
    flopy.mf6.ModflowGwfic(gwf, strt=100.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[cell(0, 0), 100.0], [cell(0, 4), 99.0]])
    flopy.mf6.ModflowGwfriv(gwf, pname="riv", filename="m.riv", stress_period_data=[
        [cell(*key), *values] for key, values in RIV.items()])
    flopy.mf6.ModflowGwfriv(gwf, pname="crisp", filename="m.crisp", stress_period_data=[
        [cell(*key), *values] for key, values in CRISP.items()])
    flopy.mf6.ModflowGwfghb(gwf, pname="edge", filename="m.edge",
                            stress_period_data=[[cell(0, 19), 97.0, 3.0]])
    flopy.mf6.ModflowGwfdrn(gwf, pname="drn", filename="m.drn",
                            stress_period_data=[[cell(0, 15), 96.0, 1.0]])
    flopy.mf6.ModflowGwfdrn(gwf, pname="seeps", filename="m.seeps",
                            stress_period_data=[[cell(0, 16), 96.5, 2.0]])
    if array_recharge:
        flopy.mf6.ModflowGwfrcha(gwf, recharge=1e-3)
    else:
        flopy.mf6.ModflowGwfrch(gwf, stress_period_data=[[cell(0, n), 1e-3] for n in range(NCPL)])
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord="m.hds", budget_filerecord="m.cbc",
                           saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])


def _run(workspace: Path, grid: str) -> Path:
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
    _flow_packages(gwf, cell, array_recharge=grid != "disu")
    sim.write_simulation(silent=True)
    assert sim.run_simulation(silent=True)[0]
    return workspace


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    root = tmp_path_factory.mktemp("several_packages")
    return {grid: _run(root / grid, grid) for grid in ("disv", "dis", "disu")}


def _heads():
    points = gpd.GeoDataFrame({"name": ["W"], "layer": [0]},
                              geometry=gpd.points_from_xy([25.0], [25.0]))
    return mf.HeadTargets(locations=points, values=pd.DataFrame(
        {"time": [0, 1], "W": [99.0, 99.0]}))


def _cal(runs, tmp_path, grid="disv", name="sp"):
    model = mf.load_mf6_run(runs[grid])
    cal = model.pest(name, workspace=tmp_path / name, start_datetime="2024-01-01")
    cal.observe(_heads())
    return cal


def _build(cal, name="sp"):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return cal.build(f"{name}.pst", noptmax=0)


def _forward_run(cal, pst, tmp_path, values: dict[str, float]) -> Path:
    """Run the template's forward run in a copy, as a PEST worker does, with each
    parameter group in ``values`` set to its value."""

    worker = tmp_path / "worker"
    shutil.rmtree(worker, ignore_errors=True)
    shutil.copytree(cal.template_workspace, worker)
    par = pst.parameter_data.copy()
    for group, value in values.items():
        names = par.index[par["pargp"] == group]
        assert len(names), f"no parameter of group {group!r}: {sorted(set(par['pargp']))}"
        par.loc[names, "parval1"] = value
    for tpl, inp in zip(pst.template_files, pst.input_files, strict=True):
        pyemu.pst_utils.write_to_template(par["parval1"], str(worker / tpl), str(worker / inp))
    done = subprocess.run([sys.executable, "forward_run.py"], cwd=worker, capture_output=True,
                          text=True)
    assert done.returncode == 0, done.stdout[-2000:] + done.stderr[-2000:]
    assert "Normal termination" in (worker / "mfsim.lst").read_text()
    return worker


def _records(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=r"\s+", header=None)


# -- finding the package -----------------------------------------------------------


def test_several_packages_of_the_type_are_refused_naming_each_and_the_spelling(runs, tmp_path):
    """``riv`` and ``crisp`` are two creeks; scaling one, or both with one number, is a
    choice the caller makes. Refused at the call, naming both and every spelling."""

    cal = _cal(runs, tmp_path)
    with pytest.raises(mf.SeveralPackagesError, match=(
            r"(?s)2 RIV packages \('riv', 'crisp'\).*package='riv'.*package='crisp'.*"
            r"package=\('riv', 'crisp'\)")):
        cal.parameterize("riv.cond")
    with pytest.raises(mf.SeveralPackagesError, match=r"'drn', 'seeps'"):
        cal.parameterize("drn.cond")
    with pytest.raises(ValueError, match=r"package='creek' names no RIV package.*'riv', 'crisp'"):
        cal.parameterize("riv.cond", package="creek")
    assert cal._native_parameter_specs == []


def test_the_build_refuses_several_packages_on_the_staged_model(runs, tmp_path):
    """The call's check asks the model object; the build's asks the run's FILES, so a
    spec that reached it unnamed is refused there too."""

    cal = _cal(runs, tmp_path)
    cal._native_parameter_specs.append(np_mod.NativeParameterSpec(target="riv.cond"))
    with pytest.raises(mf.SeveralPackagesError, match=r"'riv', 'crisp'"):
        _build(cal)


def test_package_picks_one_river_and_the_forward_run_scales_only_it(runs, tmp_path):
    """``package='crisp'`` scales crisp's files and nothing of ``riv``'s; before, a RIV
    recipe written like the others would have globbed ``m.riv_*`` and scaled ``riv``."""

    cal = _cal(runs, tmp_path)
    spec = cal.parameterize("riv.cond", package="crisp", bounds=(0.1, 10.0))
    assert spec.name == "rivcond_crisp"
    pst = _build(cal)
    assert spec.resolved_files == ["m.crisp_stress_period_data_1.txt"]
    assert [target.package for target in spec.resolved_targets] == ["crisp"]
    info = pd.read_csv(cal.template_workspace / "mult2model_info.csv")
    assert set(info["model_file"]) == {"m.crisp_stress_period_data_1.txt"}

    worker = _forward_run(cal, pst, tmp_path, {"rivcond_crisp": 3.0})
    crisp = _records(worker / "m.crisp_stress_period_data_1.txt")
    assert np.allclose(crisp.iloc[:, 3], [3.0 * values[1] for values in CRISP.values()])
    riv = _records(worker / "m.riv_stress_period_data_1.txt")
    assert np.allclose(riv.iloc[:, 3], [values[1] for values in RIV.values()])


def test_a_sequence_of_packages_shares_one_parameter(runs, tmp_path):
    cal = _cal(runs, tmp_path)
    spec = cal.parameterize("riv.cond", package=["RIV", "crisp"], bounds=(0.1, 10.0))
    assert spec.package == ("riv", "crisp") and spec.name == "rivcond_riv_crisp"
    pst = _build(cal)
    assert int((pst.parameter_data["pargp"] == "rivcond_riv_crisp").sum()) == 1
    worker = _forward_run(cal, pst, tmp_path, {"rivcond_riv_crisp": 3.0})
    for filename, table in (("m.riv_stress_period_data_1.txt", RIV),
                            ("m.crisp_stress_period_data_1.txt", CRISP)):
        assert np.allclose(_records(worker / filename).iloc[:, 3],
                           [3.0 * values[1] for values in table.values()])
    with pytest.raises(ValueError, match="names 'riv' twice"):
        cal.parameterize("riv.stage", package=("riv", "riv"))


def test_a_lone_package_not_named_after_its_type_is_found(runs, tmp_path):
    """The model's only GHB is named ``edge``: the glob looked for ``m.ghb_*`` and raised
    ``FileNotFoundError``; FloPy's records name ``m.edge_stress_period_data_1.txt``."""

    cal = _cal(runs, tmp_path)
    cond = cal.parameterize("ghb.cond", bounds=(0.1, 10.0))
    head = cal.parameterize("ghb.bhead")
    pst = _build(cal)
    assert cond.resolved_files == head.resolved_files == ["m.edge_stress_period_data_1.txt"]
    worker = _forward_run(cal, pst, tmp_path, {"ghbcond": 2.0, "ghbbhead": -1.5})
    edge = _records(worker / "m.edge_stress_period_data_1.txt")
    assert np.allclose(edge.iloc[:, 2:4].to_numpy(), [[95.5, 6.0]])


def test_array_recharge_is_refused_naming_its_form(runs, tmp_path):
    """``recharge`` scales a list RCH's records; on an array RCH the build used to die
    blaming ``set_all_data_external``."""

    cal = _cal(runs, tmp_path)
    cal.parameterize("recharge")
    with pytest.raises(ValueError, match=r"(?s)package 'rcha_0'.*ARRAY form \(READASARRAYS\)"):
        _build(cal)


# -- riv.stage moves the pair ------------------------------------------------------


@pytest.mark.parametrize("grid", ["disv", "dis", "disu"])
def test_river_stage_moves_its_bottom_with_it_and_mf6_runs(runs, tmp_path, grid):
    """An offset of +1.5 on both rivers moves every stage AND every bottom by 1.5,
    leaves conductance alone, and MF6 terminates normally -- on each grid, whose cellid
    puts the columns in different places (ledger 420)."""

    cal = _cal(runs, tmp_path, grid)
    spec = cal.parameterize("riv.stage", package=("riv", "crisp"))
    assert spec.additive and spec.bounds == (-2.0, 2.0) and spec.resolved_transform == "none"
    pst = _build(cal)
    worker = _forward_run(cal, pst, tmp_path, {"rivstage_riv_crisp": 1.5})
    width = {"disv": 2, "dis": 3, "disu": 1}[grid]
    for filename, table in (("m.riv_stress_period_data_1.txt", RIV),
                            ("m.crisp_stress_period_data_1.txt", CRISP)):
        rows = _records(worker / filename)
        stage, cond, rbot = (rows.iloc[:, width + offset].to_numpy() for offset in range(3))
        assert np.allclose(stage, [values[0] + 1.5 for values in table.values()])
        assert np.allclose(cond, [values[1] for values in table.values()])
        assert np.allclose(rbot, [values[2] + 1.5 for values in table.values()])


def test_a_river_stage_below_its_bottom_is_refused_at_build(runs, tmp_path):
    """MF6 refuses such a record at any icelltype; with the pair moving together every
    realization would carry it."""

    run = tmp_path / "low"
    shutil.copytree(runs["disv"], run)
    path = run / "m.crisp"
    path.write_text(path.read_text().replace("9.80000000E+01", "9.50000000E+01"))
    model = mf.load_mf6_run(run)
    cal = model.pest("low", workspace=tmp_path / "pest_low", start_datetime="2024-01-01")
    cal.observe(_heads())
    cal.parameterize("riv.stage", package="crisp")
    with pytest.raises(ValueError, match=r"stage 95 < rbot 96"):
        _build(cal, "low")


# -- what an argument means, and the ones that cannot act --------------------------


@pytest.mark.parametrize("target", ["chd", "ghb.bhead", "drn.elev", "riv.stage"])
def test_an_elevation_is_an_additive_offset_and_never_a_multiplier(target):
    """Ledger 518: ``FieldSpec.kind == "elevation"`` makes the offset the
    default; ``chd`` and ``ghb.bhead`` were log multipliers with bounds (0.5, 2), so a
    400 ft head ranged 200-800 ft. A multiplier there is refused."""

    spec = np_mod.NativeParameterSpec(target=target)
    assert spec.additive is True and spec.bounds == (-2.0, 2.0)
    assert spec.resolved_transform == "none"
    with pytest.raises(ValueError, match="ELEVATION.*drop additive=False"):
        np_mod.NativeParameterSpec(target=target, additive=False)
    for target_kind in ("ghb.cond", "riv.cond", "drn.cond", "wel", "recharge", "k"):
        multiplier = np_mod.NativeParameterSpec(target=target_kind)
        assert multiplier.additive is False and multiplier.bounds == (0.5, 2.0)


@pytest.mark.parametrize("kwargs, start", [
    ({"target": "drn.elev", "bounds": (0.5, 2.0)}, "0"),
    ({"target": "k", "bounds": (2.0, 5.0)}, "1"),
    ({"target": "k", "style": "pilotpoints", "bounds": (2.0, 5.0)}, "1"),
    ({"target": "ghb.bhead", "bounds": (-3.0, -1.0)}, "0"),
])
def test_bounds_that_exclude_the_starting_value_are_refused(kwargs, start):
    """PEST would start every such parameter outside its own bounds (params-5); pilot
    points start each point at 1.0, so they are held to the same rule."""

    with pytest.raises(ValueError, match=rf"exclude {start}, the value every"):
        np_mod.NativeParameterSpec(**kwargs)
    with pytest.raises(ValueError, match="lower < upper"):
        np_mod.NativeParameterSpec(target="k", bounds=(2.0, 0.5))
    with pytest.raises(ValueError, match="must be above 0"):
        np_mod.NativeParameterSpec(target="k", bounds=(0.0, 2.0))
    assert np_mod.NativeParameterSpec(target="k", bounds=(1.0, 5.0)).bounds == (1.0, 5.0)


@pytest.mark.parametrize("kwargs, message", [
    ({"target": "ghb.cond", "layers": (2,)}, r"layers= on list target 'ghb.cond'.*zones="),
    ({"target": "k", "style": "pilotpoints", "additive": True}, r"pilot points MULTIPLY"),
    ({"target": "k", "style": "pilotpoints", "extra": {"initial_value": 2.0}},
     r"\['initial_value'\].*never call"),
    ({"target": "k", "temporal": 30.0}, r"temporal=30.0: nothing reads"),
    ({"target": "k", "style": "constant", "correlation": 500.0}, r"constant parameter has nothing"),
    ({"target": "k", "style": "zone", "correlation": 500.0}, r"zone parameter has nothing"),
    ({"target": "k", "style": "grid", "anisotropy": 3.0}, r"anisotropy=3.0 without correlation="),
    ({"target": "k", "style": "grid", "nugget": 0.5}, r"nugget=0.5 without correlation="),
    ({"target": "k", "style": "grid", "correlation": 500.0, "bearing": 30.0},
     r"bearing=30.0 with anisotropy=1.0"),
    ({"target": "k", "style": "pp"}, r"'pp' is not one of.*myflopy's is style='pilotpoints'"),
    ({"target": "k", "style": "pilot_points"}, r"style='pilotpoints'"),
    ({"target": "k", "style": "kl"}, r"style='kl' is not one of"),
    # Review of R1-params: pyEMU's own pilot points reject a list file on every grid.
    ({"target": "ghb.cond", "style": "pilotpoints"}, r"pilotpoints' on list target.*'grid'"),
    ({"target": "recharge", "style": "PilotPoints"}, r"pilotpoints' on list target"),
])
def test_an_argument_that_cannot_act_is_refused(kwargs, message):
    """params-10 and params-12: each of these was accepted and did nothing (or reached
    pyEMU's structured-grid pilot points and died on DISV)."""

    with pytest.raises(ValueError, match=message):
        np_mod.NativeParameterSpec(**kwargs)


def test_a_repeated_parameter_name_is_refused(runs, tmp_path):
    cal = _cal(runs, tmp_path)
    cal.parameterize("k", style="constant")
    with pytest.raises(ValueError, match=r"Parameter name 'K' is already declared.*name="):
        cal.parameterize("k33", name="K")
    with pytest.raises(ValueError, match=r"'k' is already declared"):
        cal.parameterize("k", style="grid")
    assert [spec.name for spec in cal._native_parameter_specs] == ["k"]


# -- layers, factors, metadata, size -----------------------------------------------


def test_grid_groups_are_named_by_their_zero_based_layer(runs, tmp_path):
    """params-15: groups were numbered by file position from 1, so ``kl1`` was layer 0
    beside a pilot-point ``kl1`` of layer 1."""

    cal = _cal(runs, tmp_path)
    cal.parameterize("k", style="grid", layers=[0, 1])
    cal.parameterize("k33", style="grid", name="kv")
    cal.parameterize("k", style="grid", layers=[1], name="one")
    pst = _build(cal)
    assert set(pst.parameter_data["pargp"]) == {"kl0", "kl1", "kvl0", "kvl1", "one"}


def test_a_whole_grid_array_keeps_its_factor(tmp_path):
    """params-7: a run whose NPF K reads ``INTERNAL FACTOR 2.0`` was re-stored one array
    per layer from the RAW values, so the template's MF6 read K = 3 where the run
    used 6. The template's per-layer files now hold the values MF6 used."""

    sim = flopy.mf6.MFSimulation("s", sim_ws=str(tmp_path / "run"), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=1)
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
    flopy.mf6.ModflowGwfdis(gwf, nlay=2, nrow=2, ncol=3, top=10.0, botm=[5.0, 0.0])
    flopy.mf6.ModflowGwfnpf(gwf, k={"data": np.full((2, 2, 3), 3.0), "factor": 2.0, "iprn": 1})
    flopy.mf6.ModflowGwfic(gwf, strt=8.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 0), 9.0], [(0, 1, 2), 8.0]])
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord="m.hds", saverecord=[("HEAD", "ALL")])
    sim.write_simulation(silent=True)
    assert re.search(r"FACTOR\s+2\.0", (tmp_path / "run" / "m.npf").read_text())
    assert sim.run_simulation(silent=True)[0]
    model = mf.load_mf6_run(tmp_path / "run")
    cal = model.pest("f", workspace=tmp_path / "pest", start_datetime="2024-01-01")
    cal.observe(mf.HeadTargets(
        locations=pd.DataFrame({"name": ["W"], "layer": [0], "cell": [1]}),
        values=pd.DataFrame({"time": [0], "W": [8.5]})))
    spec = cal.parameterize("k")
    _build(cal, "f")
    assert [target.layer for target in spec.resolved_targets] == [0, 1]
    for filename in spec.resolved_files:
        assert np.allclose(np.loadtxt(cal.template_workspace / filename), 6.0)


def test_the_build_metadata_records_every_parameter_call(runs, tmp_path):
    """build-6: ``parameter_sets`` was always ``[]`` -- its loop read a list nothing
    filled -- so nothing on disk said what a template's parameters scale."""

    cal = _cal(runs, tmp_path)
    cal.parameterize("riv.stage", package="crisp")
    cal.parameterize("k", style="grid", layers=[0], bounds=(0.2, 5.0), physical=(0.1, 50.0))
    cal.parameterize("ghb.cond", style="zone", zones=np.repeat([1, 2], NCPL // 2),
                     bounds=(0.1, 10.0))
    _build(cal)
    metadata = json.loads((cal.template_workspace / "myflopy_pest_metadata.json").read_text())
    assert "parameter_sets" not in metadata
    records = {record["name"]: record for record in metadata["parameters"]}
    assert set(records) == {"rivstage_crisp", "k", "ghbcond"}
    stage = records["rivstage_crisp"]
    assert stage["target"] == "riv.stage" and stage["package"] == "crisp"
    assert stage["additive"] is True and stage["transform"] == "none"
    assert stage["bounds"] == [-2.0, 2.0] and stage["pargps"] == ["rivstage_crisp"]
    assert stage["files"] == [{"filename": "m.crisp_stress_period_data_1.txt",
                               "package": "crisp", "layer": None, "period": 0}]
    k = records["k"]
    assert k["style"] == "grid" and k["layers"] == [0] and k["physical"] == [0.1, 50.0]
    assert k["files"] == [{"filename": "m.npf_k_layer1.txt", "package": "npf",
                           "layer": 0, "period": None}]
    assert k["pargps"] == ["k"]
    assert records["ghbcond"]["zones"] == [1, 2]


def test_a_correlated_grid_block_above_the_limit_is_refused(runs, tmp_path, monkeypatch):
    """build-11: one correlated block's prior is a dense n x n covariance; above the
    limit the build refuses, naming its size and pilot points. The limit is lowered
    here so a 20-cell layer crosses it; the real one is measured by
    ``docs/reproducers/correlated_grid_prior_size.py``."""

    monkeypatch.setattr(np_mod, "CORRELATED_GRID_LIMIT", 10)
    cal = _cal(runs, tmp_path)
    cal.parameterize("k", style="grid", layers=[0], correlation=30.0)
    with pytest.raises(ValueError, match=r"(?s)puts 20 parameters in one correlated block "
                                         r"\('k'\).*style='pilotpoints'"):
        _build(cal)
    uncorrelated = _cal(runs, tmp_path, name="u")
    uncorrelated.parameterize("k", style="grid", layers=[0])
    _build(uncorrelated, "u")


# -- review of R1-params ------------------------------------------------------------


def test_a_rebuild_keeps_the_river_bottom_moving_with_its_stage(runs, tmp_path):
    """The files whose forward run carries rbot were remembered on the project ACROSS
    builds, so a second ``build()`` -- a re-run notebook cell, an observation added --
    skipped every one: its template moved the stage alone and the declared depth
    changed in every realization. Each build makes its own forward run."""

    cal = _cal(runs, tmp_path)
    cal.parameterize("riv.stage", package="crisp")
    _build(cal)
    pst = _build(cal)
    assert "carry_river_bottom(" in (cal.template_workspace / "forward_run.py").read_text()
    worker = _forward_run(cal, pst, tmp_path, {"rivstage_crisp": -1.0})
    rows = _records(worker / "m.crisp_stress_period_data_1.txt")
    assert np.allclose(rows.iloc[:, 2], [values[0] - 1.0 for values in CRISP.values()])
    assert np.allclose(rows.iloc[:, 4], [values[2] - 1.0 for values in CRISP.values()])


def test_numpy_layers_and_correlation_reach_the_metadata_as_python_numbers(runs, tmp_path):
    """``layers=np.arange(2)`` and ``correlation=np.int64(30)`` built every parameter and
    then crashed writing the metadata: JSON refuses ``numpy.int64``."""

    cal = _cal(runs, tmp_path)
    spec = cal.parameterize("k", style="grid", layers=np.arange(2), correlation=np.int64(30))
    assert all(type(layer) is int for layer in spec.layers) and spec.layers == (0, 1)
    assert type(spec.correlation) is float
    assert np_mod.NativeParameterSpec(target="k", layers=np.int64(1)).layers == (1,)
    _build(cal)
    metadata = json.loads((cal.template_workspace / "myflopy_pest_metadata.json").read_text())
    (record,) = metadata["parameters"]
    assert record["layers"] == [0, 1] and record["correlation"] == 30.0
    for layers in ([1.5], [-1], [True], ["0"], []):
        with pytest.raises(ValueError, match="layers="):
            np_mod.NativeParameterSpec(target="k", layers=layers)

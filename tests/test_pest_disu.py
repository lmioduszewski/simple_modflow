"""PEST on an unstructured (DISU) model (ledger 452).

Ledger 420 gave DISU a one-column cellid by MF6's file format, but no PEST build of a
DISU model had been run. Run, every array parameter died in pyEMU's
``write_array_tpl`` (``IndexError: too many indices for array``: handed FloPy's
``UnstructuredGrid`` it indexed the one-dimensional centre array at ``[i, j]``), a
zone array was the wrong shape for pyEMU (``passed shape (24, 1) != zone_array.shape
(24,)``), a pilot-point cell table came out half the length of its centres (FloPy's
IHC-inferred ``nlay`` split the one node array), and a head target in layer 1 matched
a node and was written to a head-file record that does not exist (``IndexError``).

The model is FloPy-written with ``get_disu_kwargs`` -- two layers of 3 x 4 -- and
opened with ``mf.load_mf6_run``: MF6 numbers its 24 cells by node, one layer of
every node (ledger 419).
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

pyemu = pytest.importorskip("pyemu")
flopy = pytest.importorskip("flopy")
from flopy.utils.gridutil import get_disu_kwargs  # noqa: E402

NLAY, NROW, NCOL = 2, 3, 4
NODES = NLAY * NROW * NCOL
K = np.linspace(1.0, 5.0, NODES)
ZONES = np.repeat([1, 2], NODES // 2)  # zone 2 is the lower layer's nodes
CHD_NODES = [row * NCOL for row in range(NROW)] + [row * NCOL + NCOL - 1 for row in range(NROW)]


def _disu_run(workspace: Path) -> Path:
    sim = flopy.mf6.MFSimulation("s", sim_ws=str(workspace), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=2, perioddata=[(1.0, 1, 1.0)] * 2)
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
    flopy.mf6.ModflowGwfdisu(gwf, **get_disu_kwargs(
        NLAY, NROW, NCOL, np.full(NCOL, 10.0), np.full(NROW, 10.0), 100.0, [50.0, 0.0],
        return_vertices=True))
    flopy.mf6.ModflowGwfnpf(gwf, k=K, k33=0.5, icelltype=0)
    flopy.mf6.ModflowGwfic(gwf, strt=100.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=(
        [[(node,), 100.0] for node in CHD_NODES[:NROW]]
        + [[(node,), 95.0] for node in CHD_NODES[NROW:]]))
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord="m.hds", saverecord=[("HEAD", "ALL")])
    sim.write_simulation(silent=True)
    assert sim.run_simulation(silent=True)[0]
    return workspace


def _heads(layer: int = 0):
    """W by its point; V, a node of the lower layer, by its cell (node 12 + 6)."""

    points = gpd.GeoDataFrame({"name": ["W"], "layer": [layer]},
                              geometry=gpd.points_from_xy([15.0], [15.0]))
    by_point = mf.HeadTargets(locations=points, values=pd.DataFrame(
        {"time": [0, 1], "W": [98.0, 98.0]}))
    by_cell = mf.HeadTargets(
        locations=pd.DataFrame({"name": ["V"], "layer": [0], "cell": [18]}),
        values=pd.DataFrame({"time": [0, 1], "V": [97.0, 97.0]}))
    return by_point, by_cell


def _observe(cal):
    by_point, by_cell = _heads()
    cal.observe(by_point)
    cal.observe(by_cell, prefix="deep")


def _build(cal, name: str):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return cal.build(f"{name}.pst", noptmax=0)


@pytest.fixture(scope="module")
def run_dir(tmp_path_factory):
    return _disu_run(tmp_path_factory.mktemp("disu_pest") / "run")


@pytest.mark.parametrize("style", ["constant", "grid", "zone", "pilotpoints"])
def test_array_and_list_parameters_on_a_disu_model_scale_what_mf6_reads(run_dir, tmp_path, style):
    """Each array style builds and the forward run, in a copy of the template as a
    PEST worker runs it, scales every K node and every CHD head -- never a node
    number -- and MF6 runs on the result. The base heads are the run's."""

    model = mf.load_mf6_run(run_dir)
    cal = model.pest("u", workspace=tmp_path / "pest", start_datetime="2024-01-01")
    extra = {"zone": {"zones": ZONES},
             "pilotpoints": {"pp_space": 2}}.get(style, {})
    cal.parameterize("k", style=style, bounds=(0.5, 2.0), **extra)
    cal.parameterize("chd", style="constant", bounds=(-1.0, 1.0))
    _observe(cal)
    pst = _build(cal, f"u_{style}")
    heads = flopy.utils.HeadFile(str(run_dir / "m.hds")).get_data(kstpkper=(0, 1)).reshape(-1)
    observed = pst.observation_data["obsval"]
    assert {name.split("_usecol:")[1].split("_")[0] for name in observed.index} == {"w", "v"}

    template = tmp_path / "worker"
    shutil.copytree(cal.template_workspace, template)
    par = pst.parameter_data.copy()
    k_names = par.index[par["pargp"].str.startswith("k")]
    par.loc[k_names, "parval1"] = 2.0
    if style == "zone":
        par.loc[[name for name in k_names if name.endswith("zone:2")], "parval1"] = 0.5
    par.loc[par.index[par["pargp"] == "chd"], "parval1"] = 0.5
    for tpl, inp in zip(pst.template_files, pst.input_files, strict=True):
        pyemu.pst_utils.write_to_template(par["parval1"], str(template / tpl), str(template / inp))
    subprocess.run([sys.executable, "forward_run.py"], cwd=template, check=True,
                   capture_output=True, text=True)

    scale = np.where(ZONES == 2, 0.5, 2.0) if style == "zone" else 2.0
    assert np.allclose(np.loadtxt(template / "m.npf_k.txt").reshape(-1), scale * K)
    chd = pd.read_csv(template / "m.chd_stress_period_data_1.txt", sep=r"\s+", header=None)
    assert chd.iloc[:, 0].tolist() == [node + 1 for node in CHD_NODES]
    assert np.allclose(chd.iloc[:, 1], [100.5] * NROW + [95.5] * NROW)
    simulated = flopy.utils.HeadFile(str(template / "m.hds")).get_data(kstpkper=(0, 1))
    assert simulated.reshape(-1).size == NODES and np.isfinite(simulated).all()
    # The base heads PEST holds are the run's: W at node 5 (row 1, column 1 of the
    # top layer, the lowest-numbered node holding its point), V at node 18.
    base = pd.read_csv(cal.template_workspace / "hds_simulated_heads.csv")
    assert np.allclose(base["W"], heads[5]) and "V" not in base
    deep = pd.read_csv(cal.template_workspace / "deep_simulated_heads.csv")
    assert np.allclose(deep["V"], heads[18])


def test_a_head_target_below_layer_0_on_disu_is_refused_naming_the_node(run_dir, tmp_path):
    """A DISU model is one layer of nodes: a target in layer 1 matched the node holding
    its point and died reading head-file layer 1 (``IndexError: index 1 is out of
    bounds for axis 0 with size 1``). It is refused, naming the spelling that reaches
    a deeper node."""

    model = mf.load_mf6_run(run_dir)
    by_point, _ = _heads(layer=1)
    with pytest.raises(ValueError, match=r"(?s)Layer\(s\) \[1\] are off the grid.*"
                                         r"one layer.*by its `cell`"):
        by_point.match_to_model(model)


@pytest.mark.slow
def test_pestpp_ies_runs_a_disu_template(run_dir, tmp_path):
    """PESTPP-IES's prior runs the DISU template, and its base realization reads the
    heads the build read from the run."""

    model = mf.load_mf6_run(run_dir)
    cal = model.pest("u", workspace=tmp_path / "pest", start_datetime="2024-01-01")
    cal.parameterize("k", style="pilotpoints", pp_space=2, correlation=30.0, bounds=(0.5, 2.0))
    cal.parameterize("chd", style="zone", zones=ZONES, bounds=(-1.0, 1.0))
    _observe(cal)
    pst = _build(cal, "u")
    expected = pst.observation_data["obsval"]
    template = cal.template_workspace
    built = {}
    for ins in sorted(template.glob("*.csv.ins")):
        read = pyemu.pst_utils.InstructionFile(str(ins)).read_output_file(
            str(template / ins.name[: -len(".ins")]))
        built.update(read["obsval"].astype(float).to_dict())
    assert len(built) == len(expected) == 4
    results = cal.prior(reals=3)
    obs = pd.read_csv(Path(results.workspace) / "u.0.obs.csv", index_col=0)
    base = obs.loc[[label for label in obs.index if str(label) == "base"]].iloc[0]
    names = sorted(built)
    assert np.allclose(base[names].to_numpy(dtype=float), [built[name] for name in names],
                       rtol=1e-5, atol=0.0)


def test_a_correlated_grid_parameter_on_stacked_disu_nodes_is_refused(run_dir, tmp_path):
    """A layered DISU's nodes in one column share a centre, so one correlated set over
    every node has a singular covariance, and PESTPP-IES's prior died in pyEMU's draw
    with ``LinAlgError: Matrix is not positive definite``. A nugget makes it drawable
    since ledger 489, but would still correlate a column's nodes as one location,
    so it is refused at build whatever the nugget, naming a node, the node above or
    below it, and the spellings that draw."""

    model = mf.load_mf6_run(run_dir)
    for nugget in (0.0, 0.5):
        cal = model.pest("u", workspace=tmp_path / f"pest{nugget}", start_datetime="2024-01-01")
        cal.parameterize("k", style="grid", correlation=30.0, nugget=nugget, bounds=(0.5, 2.0))
        _observe(cal)
        with pytest.raises(ValueError, match=r"(?s)24 of its 24 nodes share a centre.*node 0 "
                                             r"at \(5, 25\) with node\(s\) \[12\].*style='zone'"):
            _build(cal, "u")

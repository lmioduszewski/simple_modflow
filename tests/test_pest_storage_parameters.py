"""``sto.ss`` and ``sto.sy`` are array calibration targets (ledger 548).

``resolve_target("sy")`` raised ``KeyError``: storage could not be calibrated at all,
though only transient data informs it and the Ten Trails recalibration fits a
transient cycle. Both are STO griddata, scaled as K is -- ``layers=``, ``zones=`` and
every style, pilot points included -- and refused where they cannot act: a run whose
every period is steady state (MF6 reads no storage there), and ``sto.sy`` on cells
that are all confined (``ICONVERT`` 0; MF6 reads SY only in a convertible cell).

The model is FloPy-written and opened with ``mf.load_mf6_run``: DISV, 2 layers, a
steady period then two transient ones, ``sy`` LAYERED (0.1, 0.2), ``ss`` one scalar
(so one whole-grid file the build relayers), ``iconvert`` 1 in layer 0 and 0 in layer
1; and a twin whose STO is steady state in every period.
"""

from __future__ import annotations

import json
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
from flopy.utils.gridutil import get_disv_kwargs  # noqa: E402

NLAY, NROW, NCOL = 2, 4, 5
NCPL = NROW * NCOL
SY = (0.1, 0.2)
SS = 1e-5


def _run(workspace: Path, *, transient: bool, iconvert=(1, 0)) -> Path:
    sim = flopy.mf6.MFSimulation("s", sim_ws=str(workspace), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=3, perioddata=[(1.0, 1, 1.0)] + [(10.0, 2, 1.0)] * 2,
                          time_units="days")
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m", save_flows=True)
    flopy.mf6.ModflowGwfdisv(gwf, **get_disv_kwargs(
        NLAY, NROW, NCOL, np.full(NCOL, 10.0), np.full(NROW, 10.0), 100.0, [50.0, 0.0]))
    flopy.mf6.ModflowGwfnpf(gwf, k=5.0, icelltype=[1, 0])
    flopy.mf6.ModflowGwfic(gwf, strt=99.0)
    flopy.mf6.ModflowGwfsto(gwf, sy=list(SY), ss=SS, iconvert=iconvert,
                            steady_state={0: True},
                            transient={1: True} if transient else None)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0), 100.0], [(0, NCPL - 1), 95.0]])
    flopy.mf6.ModflowGwfwel(gwf, stress_period_data={1: [[(1, 12), -5.0]]})
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord="m.hds", budget_filerecord="m.cbc",
                           saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])
    sim.write_simulation(silent=True)
    assert sim.run_simulation(silent=True)[0]
    return workspace


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    root = tmp_path_factory.mktemp("storage")
    # "half": layer 0 convertible in its first ten cells only, layer 1 confined.
    half = np.zeros((NLAY, NCPL), dtype=int)
    half[0, :NCPL // 2] = 1
    return {"transient": _run(root / "transient", transient=True),
            "steady": _run(root / "steady", transient=False),
            "half": _run(root / "half", transient=True, iconvert=half)}


def _cal(runs, tmp_path, kind="transient", name="st"):
    model = mf.load_mf6_run(runs[kind])
    cal = model.pest(name, workspace=tmp_path / name, start_datetime="2024-01-01")
    points = gpd.GeoDataFrame({"name": ["W"], "layer": [0]},
                              geometry=gpd.points_from_xy([25.0], [25.0]))
    cal.observe(mf.HeadTargets(locations=points, values=pd.DataFrame(
        {"time": [0, 1, 2], "W": [99.0, 99.0, 99.0]})))
    return cal


def _build(cal, name="st"):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return cal.build(f"{name}.pst", noptmax=0)


def _forward_run(cal, pst, tmp_path, values) -> Path:
    worker = tmp_path / "worker"
    shutil.rmtree(worker, ignore_errors=True)
    shutil.copytree(cal.template_workspace, worker)
    par = pst.parameter_data.copy()
    for key, value in values.items():
        names = par.index[par["pargp"] == key] if key in set(par["pargp"]) else [key]
        par.loc[names, "parval1"] = value
    for tpl, inp in zip(pst.template_files, pst.input_files, strict=True):
        pyemu.pst_utils.write_to_template(par["parval1"], str(worker / tpl), str(worker / inp))
    done = subprocess.run([sys.executable, "forward_run.py"], cwd=worker, capture_output=True,
                          text=True)
    assert done.returncode == 0, done.stdout[-2000:] + done.stderr[-2000:]
    assert "Normal termination" in (worker / "mfsim.lst").read_text()
    return worker


def _values(path: Path) -> np.ndarray:
    return np.loadtxt(path).reshape(-1)


@pytest.mark.parametrize("alias, canonical", [
    ("ss", "sto.ss"), ("specific_storage", "sto.ss"), ("sto.ss", "sto.ss"),
    ("sy", "sto.sy"), ("Specific_Yield", "sto.sy"), ("sto.sy", "sto.sy"),
])
def test_the_storage_targets_and_their_aliases_resolve(alias, canonical):
    recipe = np_mod.resolve_target(alias)
    assert (recipe.canonical, recipe.family, recipe.package, recipe.variable) == (
        canonical, "array", "sto", canonical.split(".")[1])
    spec = np_mod.NativeParameterSpec(target=alias)
    assert spec.additive is False and spec.bounds == np_mod.MULTIPLIER_BOUNDS


def test_sy_on_layer_0_scales_only_that_layer_s_file(runs, tmp_path):
    cal = _cal(runs, tmp_path)
    spec = cal.parameterize("sto.sy", layers=[0], bounds=(0.5, 4.0))
    pst = _build(cal)
    assert spec.resolved_files == ["m.sto_sy_layer1.txt"]
    worker = _forward_run(cal, pst, tmp_path, {"stosy": 3.0})
    assert np.allclose(_values(worker / "m.sto_sy_layer1.txt"), 3.0 * SY[0])
    assert np.allclose(_values(worker / "m.sto_sy_layer2.txt"), SY[1])


def test_a_scalar_ss_is_relayered_and_each_layer_is_its_own_grid_set(runs, tmp_path):
    """``ss=1e-5`` is one whole-grid file; the build splits it per layer (its FACTOR
    kept) so a grid set per layer can index it."""

    cal = _cal(runs, tmp_path)
    spec = cal.parameterize("ss", style="grid", bounds=(0.5, 4.0))
    pst = _build(cal)
    assert [(t.filename, t.layer) for t in spec.resolved_targets] == [
        ("m.sto_ss_layer1.txt", 0), ("m.sto_ss_layer2.txt", 1)]
    groups = pst.parameter_data.groupby("pargp").size()
    assert groups.to_dict() == {"stossl0": NCPL, "stossl1": NCPL}
    # The build metadata names the groups the parameters are in, per layer
    # (review of R2-arrays: it named one group, "stoss", which held none of them).
    metadata = json.loads((cal.template_workspace / "myflopy_pest_metadata.json").read_text())
    (record,) = metadata["parameters"]
    assert [group["pargps"] for group in record["period_groups"]] == [["stossl0", "stossl1"]]
    worker = _forward_run(cal, pst, tmp_path, {"stossl1": 2.0})
    assert np.allclose(_values(worker / "m.sto_ss_layer1.txt"), SS)
    assert np.allclose(_values(worker / "m.sto_ss_layer2.txt"), 2.0 * SS)


def test_pilot_points_on_sy_multiply_sy_not_k(runs, tmp_path):
    """The pilot-point base is the template's own SY file, found by type -- never another
    array (k33 once read K)."""

    cal = _cal(runs, tmp_path)
    cal.parameterize("sy", style="pilotpoints", layers=[0], pp_spacing=20.0, bounds=(0.5, 2.0))
    pst = _build(cal)
    base = cal.template_workspace / "m.sto_sy_layer1.txt.pp_base"
    assert np.allclose(_values(base), SY[0])
    worker = _forward_run(cal, pst, tmp_path, {"stosyl0": 1.5})
    assert np.allclose(_values(worker / "m.sto_sy_layer1.txt"), 1.5 * SY[0])
    assert np.allclose(_values(worker / "m.npf_k.txt" if (worker / "m.npf_k.txt").exists()
                               else worker / "m.npf_k_layer1.txt"), 5.0)


# -- refusals ---------------------------------------------------------------------------


@pytest.mark.parametrize("target", ["sto.ss", "sto.sy"])
def test_storage_on_a_run_that_is_steady_throughout_is_refused(runs, tmp_path, target):
    cal = _cal(runs, tmp_path, kind="steady")
    cal.parameterize(target)
    with pytest.raises(ValueError, match=r"every one of the run's 3 stress period\(s\) is steady"):
        _build(cal)


def test_sy_on_confined_cells_alone_is_refused(runs, tmp_path):
    cal = _cal(runs, tmp_path)
    cal.parameterize("sto.sy", layers=[1])
    with pytest.raises(ValueError, match=r"(?s)every cell of the layers and zones selected has "
                                         r"ICONVERT 0.*convertible"):
        _build(cal)


def test_a_grid_sy_set_on_a_confined_layer_is_refused_naming_the_layer(runs, tmp_path):
    """Each layer file of a grid style is a set of its own: a set on the confined layer
    is noise, even though layer 0's set acts."""

    cal = _cal(runs, tmp_path)
    cal.parameterize("sto.sy", style="grid")
    with pytest.raises(ValueError, match=r"every cell of layer 1 \(m.sto_sy_layer2.txt\) selected"):
        _build(cal)


def test_a_zone_label_only_on_confined_cells_is_refused_naming_it(runs, tmp_path):
    cal = _cal(runs, tmp_path)
    zones = np.ones(NCPL, dtype=int)
    cal.parameterize("sto.sy", style="zone", zones=zones, layers=[1])
    with pytest.raises(ValueError, match=r"every cell of zone 1 on the layers selected"):
        _build(cal)


def test_a_pilot_point_sy_set_on_confined_cells_of_one_zone_is_refused_naming_it(runs, tmp_path):
    """Pilot points under ``zones=`` are one set per layer file AND label (ledger 527): a
    label whose cells in the layer are all confined is a set that cannot act, though
    the layer's other label acts. The check read the whole layer, so it passed (ledger
    548, restated); a zone over convertible cells still builds."""

    zones = np.where(np.arange(NCPL) < NCPL // 2, 1, 2)
    cal = _cal(runs, tmp_path, kind="half")
    cal.parameterize("sto.sy", style="pilotpoints", layers=[0], zones=zones, pp_spacing=10.0)
    with pytest.raises(ValueError, match=r"every cell of zone 2 of layer 0 \(m.sto_sy_layer1.txt\) "
                                         r"has ICONVERT 0"):
        _build(cal)

    acting = _cal(runs, tmp_path, kind="half", name="ok")
    acting.parameterize("sto.sy", style="pilotpoints", layers=[0], pp_spacing=10.0,
                        zones=np.where(np.arange(NCPL) < NCPL // 2, 1, 0))
    _build(acting, name="ok")


def test_a_constant_sy_over_both_layers_acts_through_the_convertible_one(runs, tmp_path):
    cal = _cal(runs, tmp_path)
    spec = cal.parameterize("sto.sy")
    _build(cal)
    assert spec.resolved_files == ["m.sto_sy_layer1.txt", "m.sto_sy_layer2.txt"]

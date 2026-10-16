"""``recharge`` on an ARRAY RCH (READASARRAYS) scales its per-period arrays (ledger 547).

Ten Trails writes its recharge with ``mf.rch.array`` -- one ``recharge`` array per
period block, never a list -- and ``parameterize("recharge")`` had only a list recipe,
so the build refused it (findings params-9, build-7). The recipe now reads the staged
package's FORM: an array RCH's ``recharge`` dataset, one file per period BLOCK. FloPy
numbers those files by the block's ORDINAL, not its period (blocks at periods 0 and
2 write ``_1`` and ``_2``), so each file's period comes from FloPy's records; ``irch``
is never touched. The period files share ONE parameter set (pyEMU's own rule) and pilot
points interpolate ONE field per run, multiplied into every period file against that
file's own base.

The model here is FloPy-written and opened with ``mf.load_mf6_run``: DISV, 2 layers,
3 periods, an array RCH with blocks at periods 0 and 2 only (period 1 repeats period
0's), and a twin with a LIST RCH. Every forward run is a real MF6 run in a copy of the
template, as a PEST worker runs it.
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
#: The recharge of each block, per cell: period 0 (also served to period 1) and 2.
RECHARGE = {0: 1e-3 * (1.0 + np.arange(NCPL) / NCPL), 2: 2e-3 * (1.0 + np.arange(NCPL) / NCPL)}


def _run(workspace: Path, *, array: bool) -> Path:
    sim = flopy.mf6.MFSimulation("s", sim_ws=str(workspace), exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=3, perioddata=[(1.0, 1, 1.0)] * 3, time_units="days",
                          start_date_time="2024-01-01")
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m", save_flows=True)
    flopy.mf6.ModflowGwfdisv(gwf, **get_disv_kwargs(
        NLAY, NROW, NCOL, np.full(NCOL, 10.0), np.full(NROW, 10.0), 100.0, [50.0, 0.0]))
    flopy.mf6.ModflowGwfnpf(gwf, k=5.0, icelltype=0)
    flopy.mf6.ModflowGwfic(gwf, strt=100.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0), 100.0], [(0, NCPL - 1), 99.0]])
    if array:
        flopy.mf6.ModflowGwfrcha(gwf, recharge=dict(RECHARGE))
    else:
        flopy.mf6.ModflowGwfrch(gwf, stress_period_data={
            per: [[(0, cell), float(values[cell])] for cell in range(NCPL)]
            for per, values in RECHARGE.items()})
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord="m.hds", budget_filerecord="m.cbc",
                           saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])
    sim.write_simulation(silent=True)
    assert sim.run_simulation(silent=True)[0]
    return workspace


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    root = tmp_path_factory.mktemp("array_recharge")
    return {"array": _run(root / "array", array=True), "list": _run(root / "list", array=False)}


def _cal(runs, tmp_path, form="array", name="ar"):
    model = mf.load_mf6_run(runs[form])
    cal = model.pest(name, workspace=tmp_path / name, start_datetime="2024-01-01")
    points = gpd.GeoDataFrame({"name": ["W"], "layer": [0]},
                              geometry=gpd.points_from_xy([25.0], [25.0]))
    cal.observe(mf.HeadTargets(locations=points, values=pd.DataFrame(
        {"time": [0, 1, 2], "W": [99.0, 99.0, 99.0]})))
    return cal


def _build(cal, name="ar"):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return cal.build(f"{name}.pst")


def _forward_run(cal, pst, tmp_path, values) -> Path:
    """Run the template's forward run in a copy, as a PEST worker does. ``values`` maps
    a parameter group to one value, or a parameter name to its value."""

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


def _recharge_in(worker: Path) -> list[float]:
    """MF6's RCH budget IN of each period, from the run's own cell-budget file."""

    budget = flopy.utils.CellBudgetFile(str(worker / "m.cbc"), precision="double")
    totals = []
    for kstpkper in budget.get_kstpkper():
        records = budget.get_data(kstpkper=kstpkper, text="RCH")
        totals.append(float(sum(record["q"].clip(min=0).sum() for record in records)))
    return totals


def _values(path: Path) -> np.ndarray:
    return np.loadtxt(path).reshape(-1)


# -- the files and their periods ------------------------------------------------------


def test_recharge_on_an_array_rch_builds_on_its_period_files_and_never_irch(runs, tmp_path):
    """The build used to refuse it (no array recipe). Exactly the two block files are
    scaled, each with the period FloPy's records give it -- the period-2 block is
    ``_2.txt``, numbered by ordinal -- and no ``irch`` file."""

    cal = _cal(runs, tmp_path)
    spec = cal.parameterize("recharge", name="rch")
    _build(cal)
    assert spec.recipe.family == "array" and spec.recipe.variable == "recharge"
    assert [(t.filename, t.period) for t in spec.resolved_targets] == [
        ("m.rcha_recharge_1.txt", 0), ("m.rcha_recharge_2.txt", 2)]
    info = pd.read_csv(cal.template_workspace / "mult2model_info.csv")
    assert set(info["model_file"]) == {"m.rcha_recharge_1.txt", "m.rcha_recharge_2.txt"}
    assert not any("irch" in name for name in info["model_file"])
    # Without periods= one set serves every period, period 1's repeat of 0 included.
    assert [(group.label, group.periods) for group in spec.period_groups] == [(None, (0, 1, 2))]


def test_a_multiplier_scales_every_period_file_and_mf6_s_recharge(runs, tmp_path):
    """x2 on the constant doubles both block files and MF6's RCH budget IN in every
    period, the period that repeats a block included."""

    cal = _cal(runs, tmp_path)
    cal.parameterize("recharge", name="rch", bounds=(0.5, 4.0))
    pst = _build(cal)
    base = _forward_run(cal, pst, tmp_path / "one", {"rch": 1.0})
    doubled = _forward_run(cal, pst, tmp_path / "two", {"rch": 2.0})
    for name, period in (("m.rcha_recharge_1.txt", 0), ("m.rcha_recharge_2.txt", 2)):
        assert np.allclose(_values(base / name), RECHARGE[period])
        assert np.allclose(_values(doubled / name), 2.0 * RECHARGE[period])
    assert np.allclose(_recharge_in(doubled), 2.0 * np.asarray(_recharge_in(base)))


def test_a_grid_parameter_is_one_set_per_cell_shared_by_the_period_files(runs, tmp_path):
    """``style="grid"`` splits per LAYER file, never per period file: ncpl parameters,
    not nperiods x ncpl (pyEMU shares one per-cell set across the files of a call)."""

    cal = _cal(runs, tmp_path)
    cal.parameterize("recharge", style="grid", name="rchg")
    pst = _build(cal)
    names = pst.parameter_data.index[pst.parameter_data["pargp"] == "rchg"]
    assert len(names) == NCPL
    cell = 7
    by_cell = pst.parameter_data.loc[names]
    target = [name for name in by_cell.index if f"i:{cell}_" in name or name.endswith(f"i:{cell}")]
    assert len(target) == 1, sorted(by_cell.index)[:3]
    worker = _forward_run(cal, pst, tmp_path, {target[0]: 3.0})
    for name, period in (("m.rcha_recharge_1.txt", 0), ("m.rcha_recharge_2.txt", 2)):
        expected = RECHARGE[period].copy()
        expected[cell] *= 3.0
        assert np.allclose(_values(worker / name), expected)


def test_pilot_points_on_recharge_interpolate_one_field_into_every_period_file(runs, tmp_path):
    """One set, one forward-run command for both period files, one factor field; each
    file is its own base times that field -- and it stacks on a constant."""

    cal = _cal(runs, tmp_path)
    cal.parameterize("recharge", name="rch", bounds=(0.5, 4.0))
    spec = cal.parameterize("recharge", style="pilotpoints", name="rchpp", pp_spacing=20.0,
                            bounds=(0.2, 5.0))
    pst = _build(cal)
    assert [item["name"] for item in spec.pilot_point_sets] == ["rchpp"]
    manifest = json.loads((cal.template_workspace / "pilot_point_sets.json").read_text())
    (record,) = manifest["sets"]
    assert record["array_files"] == ["m.rcha_recharge_1.txt", "m.rcha_recharge_2.txt"]
    assert record["layer"] is None and record["periods"] == [0, 1, 2]
    # The points CSV and the set's frame carry the manifest's layer: blank / None, not a
    # -1 a reader would take for a layer number (integration of R2-arrays).
    (item,) = cal._pilot_point_sets
    assert item["frame"]["layer"].isna().all()
    points_text = (cal.template_workspace / record["points_csv"]).read_text().splitlines()
    assert points_text[0] == "parnme,x,y,layer,zone"
    assert all(line.split(",")[3] == "" for line in points_text[1:])
    forward =(cal.template_workspace / "forward_run.py").read_text()
    assert forward.count("apply_pilotpoints_to_array(") == 2  # its definition and ONE call

    points = pst.parameter_data.index[pst.parameter_data["pargp"] == "rchpp"]
    values = dict(zip(points, np.linspace(0.5, 3.0, len(points)), strict=True))
    worker = _forward_run(cal, pst, tmp_path, {"rch": 2.0, **values})
    first = _values(worker / "m.rcha_recharge_1.txt") / (2.0 * RECHARGE[0])
    second = _values(worker / "m.rcha_recharge_2.txt") / (2.0 * RECHARGE[2])
    assert np.allclose(first, second)                 # one field in both files
    assert first.min() < 0.9 and first.max() > 1.5   # and it is the points', not 1.0
    stored = np.load(cal.template_workspace / record["factors"])
    pv = np.array([values[name] for name in pd.read_csv(
        cal.template_workspace / record["points_csv"])["parnme"]])
    expected = 10.0 ** np.sum(stored["weight"] * np.log10(pv)[stored["index"]], axis=1)
    assert np.allclose(first[stored["cells"]], expected)


# -- refusals ---------------------------------------------------------------------------


def test_layers_on_recharge_is_refused_in_either_form():
    with pytest.raises(ValueError, match=r"layers= on target 'recharge'.*no layers to restrict"):
        np_mod.NativeParameterSpec(target="recharge", layers=[0])
    with pytest.raises(ValueError, match=r"layers= on target 'recharge'"):
        np_mod.NativeParameterSpec(target="recharge", style="pilotpoints", pp_spacing=20.0,
                                   layers=[0])


def test_pilot_points_on_a_list_rch_are_refused_at_the_call(runs, tmp_path):
    """Pilot points interpolate onto an array; a list RCH has records. Refused when it is
    written, naming the array form and the styles that act."""

    cal = _cal(runs, tmp_path, form="list")
    with pytest.raises(ValueError, match=r"(?s)hold it as LIST records.*mf.rch.array.*'grid'"):
        cal.parameterize("recharge", style="pilotpoints", pp_spacing=20.0)
    assert cal._native_parameter_specs == []
    spec = cal.parameterize("recharge", name="rch")
    assert spec.recipe.family == "list"


def test_capture_on_several_period_files_is_refused(runs, tmp_path):
    """A capture records one field under one prefix; two period files would collide."""

    cal = _cal(runs, tmp_path)
    cal.parameterize("recharge", name="rch", capture=True)
    with pytest.raises(NotImplementedError, match=r"capture=True on 'recharge'.*2 per-period files"):
        _build(cal)


def test_the_list_form_still_scales_its_records(runs, tmp_path):
    """A list RCH keeps the list recipe: its block files' records, its value column."""

    cal = _cal(runs, tmp_path, form="list")
    spec = cal.parameterize("recharge", name="rch", bounds=(0.5, 4.0))
    pst = _build(cal)
    assert spec.recipe.family == "list"
    assert [t.period for t in spec.resolved_targets] == [0, 2]
    worker = _forward_run(cal, pst, tmp_path, {"rch": 2.0})
    records = pd.read_csv(worker / spec.resolved_files[1], sep=r"\s+", header=None)
    assert np.allclose(records.iloc[:, 2], 2.0 * RECHARGE[2])

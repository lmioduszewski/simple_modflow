"""``cal.run_template``: the one run of a built template, scored as PEST++ scores it (ledger 605).

The canonical TESTING model (21x21, 6 periods) calibrated on ``kl0``/``k33l0`` against
its head targets ``hd``, two one-sided bound sets -- ``land`` (``at_most``, weight 2, three
of six cells violated at the start) and ``floor`` (``at_least``, weight 3, two violated)
-- and a zero-weight forecast. Pinned here:

- the ORACLE: phi by group from ``run_template`` equals what PESTPP-IES itself writes in
  ``<case>.phi.group.csv`` for the same parameter values (its NOPTMAX 0 run in a launch
  copy of the template), and every simulated value equals its ``.base.rei`` -- at the
  start and at chosen values;
- ``outputs_from=`` reads a completed run's outputs without running MODFLOW 6, and
  refuses a run that is not a normal, current run of the template's model, clock and
  mesh (another mesh was read in silence, 20.9 ft median);
- the template, the run directory and ``cal.pst`` are left as they were;
- the refusals of ``values=``, ``workspace=`` and ``overwrite=``;
- ``run_ies(iterations=0)`` is refused naming ``run_template``.
"""

from __future__ import annotations

import contextlib
import io
import os
import re
import shutil
import subprocess
import time
import warnings
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from myflopy.modflow.mf6.canonical_calibration import build_canonical_calibration_demo
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig
from myflopy.modflow.mf6.observations import HeadTargets
from myflopy.modflow.mf6.pest.project import TEMPLATE_RUN_MARKER

pyemu = pytest.importorskip("pyemu")

NEEDS_PESTPP = pytest.mark.skipif(shutil.which("pestpp-ies") is None,
                                  reason="needs pestpp-ies on PATH")
NEEDS_MF6 = pytest.mark.skipif(shutil.which("mf6") is None, reason="needs mf6 on PATH")
pytestmark = [pytest.mark.slow, NEEDS_MF6]


def _demo(root: Path, config: CanonicalModelConfig):
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return build_canonical_calibration_demo(root, config=config, n_head_wells=8)


@pytest.fixture(scope="module")
def preflight(tmp_path_factory):
    """The calibration the module docstring describes, built; its demo, cal and root."""

    root = tmp_path_factory.mktemp("run_template")
    demo = _demo(root / "model", CanonicalModelConfig.testing())
    model = demo.model
    cells = [int(c) for c in demo.head_targets.locations_gdf["cell"].iloc[:6]]
    probe = HeadTargets.from_cells(cells=cells, values=np.full((model.nper, 6), np.nan),
                                   layers=0, names=[f"b{i}" for i in range(6)])
    start = probe.simulated_heads(model)[[f"b{i}" for i in range(6)]].to_numpy(dtype=float)
    land = HeadTargets.from_cells(cells=cells, values=start + np.array([-1, -1, -1, 1, 1, 1.0]),
                                  layers=0, names=[f"b{i}" for i in range(6)], weights=2.0,
                                  inequality="at_most")
    floor = HeadTargets.from_cells(cells=cells, values=start + np.array([1, 1, -1, -1, -1, -1.0]),
                                   layers=0, names=[f"g{i}" for i in range(6)], weights=3.0,
                                   inequality="at_least")
    cal = model.pest("t", start_datetime="2020-01-01", workspace=root / "pest_t")
    cal.parameterize("k", style="constant", layers=[0], name="kl0")
    cal.parameterize("k33", style="constant", layers=[0], name="k33l0")
    cal.observe(demo.head_targets, prefix="hd")
    cal.observe(land, prefix="land")
    cal.observe(floor, prefix="floor")
    cal.forecast(demo.forecast_targets, prefix="fore")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build("t.pst")
    return SimpleNamespace(demo=demo, model=model, cal=cal, root=root)


@pytest.fixture(scope="module")
def start(preflight):
    """``run_template`` at the starting values, MODFLOW 6 run in the copy."""

    return preflight.cal.run_template(preflight.root / "start")


def _kl0(cal) -> str:
    return next(name for name in cal.pst.par_names if "kl0" in name and "k33" not in name)


def _pestpp_once(cal, directory: Path, *, values: dict | None = None) -> Path:
    """PESTPP-IES at NOPTMAX 0 in a launch copy of the template, ``values`` written as parval1."""

    cal._launch_copy(directory, method="test", announce=False)
    if values:
        pst = pyemu.Pst(str(directory / "t.pst"))
        for name, value in values.items():
            pst.parameter_data.loc[name, "parval1"] = value
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pst.write(str(directory / "t.pst"), version=2)
    done = subprocess.run(["pestpp-ies", "t.pst"], cwd=directory, capture_output=True,
                          text=True, check=False)
    assert done.returncode == 0, done.stdout[-2000:] + done.stderr[-2000:]
    return directory


def _assert_matches_pestpp(table: pd.DataFrame, directory: Path) -> None:
    groups = pd.read_csv(directory / "t.phi.group.csv")
    base = groups.loc[groups["par_realization"].astype(str) == "base"].iloc[0]
    columns = [c for c in groups.columns
               if c not in ("iteration", "total_runs", "obs_realization", "par_realization")]
    assert {"hd", "l_land", "g_floor"} <= set(columns)
    mine = table.groupby("group")["phi"].sum()
    for column in columns:
        assert float(mine.get(column, 0.0)) == pytest.approx(float(base[column]), rel=1e-5,
                                                             abs=1e-9), column
    rei = pyemu.pst_utils.read_resfile(str(directory / "t.base.rei"))
    modelled = rei["modelled"].astype(float)
    modelled.index = modelled.index.str.lower()
    assert float((table["simulated"] - modelled.reindex(table.index)).abs().max()) <= 1e-8


# -- the oracle: PESTPP-IES's own phi.group.csv ------------------------------------------


@NEEDS_PESTPP
def test_phi_by_group_is_what_pestpp_ies_writes_for_the_same_run(preflight, start, tmp_path):
    _assert_matches_pestpp(start, _pestpp_once(preflight.cal, tmp_path / "pestpp"))
    # Both bound sets bite at the start, so the oracle is not checking zeros.
    phi = start.groupby("group")["phi"].sum()
    assert phi["l_land"] > 0 and phi["g_floor"] > 0 and phi["hd"] > 0


@NEEDS_PESTPP
def test_chosen_values_are_what_pestpp_ies_runs_at_them(preflight, start, tmp_path):
    cal = preflight.cal
    name = _kl0(cal)
    table = cal.run_template(tmp_path / "half", values={name: 0.5, **{
        other: 1.0 for other in cal.pst.adj_par_names if other != name}})
    _assert_matches_pestpp(table, _pestpp_once(cal, tmp_path / "pestpp", values={name: 0.5}))
    assert table["phi"].sum() != pytest.approx(start["phi"].sum())


def test_a_bound_scores_zero_where_it_holds_and_its_misfit_where_violated(start):
    land = start.loc[start["group"] == "l_land"]
    floor = start.loc[start["group"] == "g_floor"]
    # Three of six land cells violated (simulated above), two of six floor cells (below),
    # in each of the six periods.
    assert (land["raw_residual"] > 0).sum() == 3 * 6 and len(land) == 6 * 6
    assert (floor["raw_residual"] < 0).sum() == 2 * 6 and len(floor) == 6 * 6
    np.testing.assert_allclose(land["residual"], land["raw_residual"].clip(lower=0.0))
    np.testing.assert_allclose(floor["residual"], floor["raw_residual"].clip(upper=0.0))
    assert (land["residual"] == 0).any() and (land["residual"] > 0).any()
    np.testing.assert_allclose(start["phi"], (start["weight"] * start["residual"]) ** 2)
    np.testing.assert_allclose(start["raw_residual"], start["simulated"] - start["obsval"])
    fore = start.loc[start["group"] == "fore"]
    assert len(fore) and (fore["weight"] == 0).all() and (fore["phi"] == 0).all()
    assert list(start.columns) == ["group", "obsval", "weight", "simulated", "residual",
                                   "raw_residual", "phi"]
    assert start["simulated"].notna().all()


# -- outputs_from: a completed run, no MODFLOW 6 -----------------------------------------


def _snapshot(directory: Path) -> dict:
    return {path.relative_to(directory).as_posix(): (path.stat().st_size, path.stat().st_mtime_ns)
            for path in directory.rglob("*") if path.is_file()}


def test_outputs_from_reads_the_run_and_runs_no_mf6(preflight, start, tmp_path):
    cal, model = preflight.cal, preflight.model
    run_dir = Path(model.workspace)
    before = (_snapshot(cal.template_workspace), _snapshot(run_dir),
              cal.pst.parameter_data.copy(), int(cal.pst.control_data.noptmax))
    table = cal.run_template(tmp_path / "outputs", outputs_from=model)
    np.testing.assert_allclose(table["simulated"], start["simulated"], rtol=0, atol=1e-8)
    # The listing was COPIED, not written by a run in the copy.
    assert (tmp_path / "outputs" / "mfsim.lst").read_bytes() == (run_dir / "mfsim.lst").read_bytes()
    # values= with outputs_from reach the result only through a reader of an INPUT:
    # this calibration has no screened well, so nothing moves.
    moved = cal.run_template(tmp_path / "moved", outputs_from=run_dir,
                             values={_kl0(cal): 0.5, **{
                                 n: 1.0 for n in cal.pst.adj_par_names if _kl0(cal) != n}})
    np.testing.assert_array_equal(moved["simulated"], table["simulated"])
    # The template, the run and cal.pst are left as they were.
    after = (_snapshot(cal.template_workspace), _snapshot(run_dir),
             cal.pst.parameter_data, int(cal.pst.control_data.noptmax))
    assert after[0] == before[0] and after[1] == before[1]
    pd.testing.assert_frame_equal(after[2], before[2])
    assert after[3] == before[3] == 0


def test_outputs_from_another_mesh_is_refused_naming_both_fingerprints(preflight, tmp_path):
    """Before, a 484-cell run against this 441-cell template read 48 heads 20.9 ft off (median)."""

    other = _demo(tmp_path / "big", CanonicalModelConfig(
        nrow=22, ncol=22, nlay=4, nper=6, cell_size=100.0, steps_per_period=1)).model
    assert other.name == preflight.model.name
    with pytest.raises(ValueError, match=r"another mesh") as raised:
        preflight.cal.run_template(tmp_path / "w", outputs_from=other)
    found = re.search(r"its mesh is DISV, 4 layer\(s\), fingerprint (\w{16}) against the "
                      r"template's DISV, 4 layer\(s\), fingerprint (\w{16})", str(raised.value))
    assert found is not None and found.group(1) != found.group(2), str(raised.value)
    assert not (tmp_path / "w").exists()


def _copied_run(preflight, directory: Path) -> Path:
    shutil.copytree(preflight.model.workspace, directory)
    return directory


def _age(path: Path, seconds: float = 3600.0) -> None:
    then = time.time() - seconds
    os.utime(path, (then, then))


def test_outputs_from_another_clock_or_model_is_refused(preflight, tmp_path):
    from tests.test_pest_handed_targets import run_toy

    run = _copied_run(preflight, tmp_path / "clock")
    tdis = next(run.glob("*.tdis"))
    lines = tdis.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.strip().upper().startswith("BEGIN PERIODDATA"))
    fields = lines[start + 2].split()
    lines[start + 2] = "  " + "  ".join([str(float(fields[0]) + 1.0), *fields[1:]])
    tdis.write_text("\n".join(lines) + "\n")
    _age(tdis)
    with pytest.raises(ValueError, match=r"stress period 1: "):
        preflight.cal.run_template(tmp_path / "w1", outputs_from=run)
    toy = run_toy(tmp_path / "toy")
    with pytest.raises(ValueError, match=rf"it has no model {preflight.model.name!r}"):
        preflight.cal.run_template(tmp_path / "w2", outputs_from=toy)


def test_outputs_from_a_run_that_did_not_finish_or_moved_on_is_refused(preflight, tmp_path):
    unfinished = _copied_run(preflight, tmp_path / "unfinished")
    listing = unfinished / "mfsim.lst"
    listing.write_text(listing.read_text().replace("Normal termination", "Still going"))
    with pytest.raises(ValueError, match=r"reads 'incomplete', not a normal termination"):
        preflight.cal.run_template(tmp_path / "w1", outputs_from=unfinished)
    moved_on = _copied_run(preflight, tmp_path / "moved_on")
    _age(moved_on / "mfsim.lst")
    with pytest.raises(ValueError, match=r"written after its mfsim.lst"):
        preflight.cal.run_template(tmp_path / "w2", outputs_from=moved_on)
    with pytest.raises(FileNotFoundError, match=r"holds no mfsim.nam"):
        preflight.cal.run_template(tmp_path / "w3", outputs_from=tmp_path)
    with pytest.raises(TypeError, match=r"outputs_from= takes a completed run"):
        preflight.cal.run_template(tmp_path / "w4", outputs_from=42)


def test_outputs_from_the_template_is_refused(preflight, tmp_path):
    with pytest.raises(ValueError, match=r"IS this calibration's template"):
        preflight.cal.run_template(tmp_path / "w", outputs_from=preflight.cal.template_workspace)



def test_outputs_from_a_pest_working_directory_is_refused(preflight, start, tmp_path):
    """A template COPY is not a model run (ledger 607).

    Its top-level forward run, target maps and instruction files replaced this
    template's in the copy: after the wells moved and the calibration was rebuilt, an
    earlier ``run_template`` folder returned the OLD wells' heads, 58 ft off, with no
    error (``docs/reproducers/run_template_is_pestpp_once.py``, part D).
    """

    cal = preflight.cal
    earlier = preflight.root / "start"  # the `start` fixture's folder: MF6 ran there
    assert (earlier / "mfsim.lst").is_file() and (earlier / TEMPLATE_RUN_MARKER).is_file()
    with pytest.raises(ValueError, match=r"is a PEST working directory, not a model run: it "
                                         r"holds myflopy_pest_build.json, "
                                         r"myflopy_template_run.json, forward_run.py"):
        cal.run_template(tmp_path / "w1", outputs_from=earlier)
    # A model run that a PEST step wrote into (a forward run beside it) is one too.
    run = _copied_run(preflight, tmp_path / "with_forward_run")
    shutil.copy2(cal.template_workspace / "forward_run.py", run / "forward_run.py")
    _age(run / "forward_run.py")
    with pytest.raises(ValueError, match=r"it holds forward_run.py"):
        cal.run_template(tmp_path / "w2", outputs_from=run)
    assert not (tmp_path / "w1").exists() and not (tmp_path / "w2").exists()


# -- workspace= and overwrite= -----------------------------------------------------------


def test_a_workspace_holding_anything_is_refused_and_kept(preflight, tmp_path):
    folder = tmp_path / "mine"
    folder.mkdir()
    (folder / "notes.txt").write_text("keep me")
    for overwrite in (False, True):
        with pytest.raises(FileExistsError, match=r"already holds 1 entry \(notes.txt\); pass a new"):
            preflight.cal.run_template(folder, overwrite=overwrite)
    assert (folder / "notes.txt").read_text() == "keep me"
    (tmp_path / "file").write_text("x")
    with pytest.raises(FileExistsError, match=r"is a FILE"):
        preflight.cal.run_template(tmp_path / "file")


def test_overwrite_deletes_only_an_earlier_run_template_of_this_calibration(preflight, tmp_path):
    cal = preflight.cal
    hand = cal.template_workspace / "by_hand.txt"
    hand.write_text("not the build's")
    try:
        with pytest.warns(UserWarning, match=r"run_template\(\): 1 entry .* is not in the copy: "
                                             r"by_hand.txt"):
            first = cal.run_template(tmp_path / "again")
    finally:
        hand.unlink()
    assert (tmp_path / "again" / TEMPLATE_RUN_MARKER).is_file()
    with pytest.raises(FileExistsError, match=r"earlier run_template\(\) of this calibration: "
                                              r"pass overwrite=True"):
        cal.run_template(tmp_path / "again")
    (tmp_path / "again" / "stale.txt").write_text("from the first run")
    second = cal.run_template(tmp_path / "again", overwrite=True)
    assert not (tmp_path / "again" / "stale.txt").exists()
    pd.testing.assert_frame_equal(first, second)
    # A launch copy carries the BUILD marker naming the template -- an IES master or a
    # worker does -- and is never deleted by overwrite=True.
    master = tmp_path / "master_like"
    cal._launch_copy(master, method="test", announce=False)
    held = sorted(path.name for path in master.iterdir())
    with pytest.raises(FileExistsError, match=r"never one a build or a PESTPP-IES launch made"):
        cal.run_template(master, overwrite=True)
    assert sorted(path.name for path in master.iterdir()) == held


@pytest.mark.parametrize("where", ["template", "template_child", "run", "run_child"])
def test_a_workspace_in_the_template_or_the_run_is_refused(preflight, where):
    template = preflight.cal.template_workspace.resolve()
    run = Path(preflight.model.workspace).resolve()
    target, relation = {
        "template": (template, "IS this calibration's template"),
        "template_child": (template / "x", "sits INSIDE this calibration's template"),
        "run": (run, "IS the model's run directory"),
        "run_child": (run / "x", "sits INSIDE the model's run directory"),
    }[where]
    with pytest.raises(ValueError, match=relation):
        preflight.cal.run_template(target)
    assert not (template / "x").exists() and not (run / "x").exists()


# -- values= -----------------------------------------------------------------------------


def test_values_must_name_every_adjustable_parameter(preflight, tmp_path):
    with pytest.raises(ValueError, match=r"Adjustable parameter\(s\) .* are not in values=: "
                                         r"name every adjustable parameter"):
        preflight.cal.run_template(tmp_path / "w", values={_kl0(preflight.cal): 0.5})
    assert not (tmp_path / "w").exists()


def test_a_name_the_control_file_lacks_is_refused(preflight, tmp_path):
    values = dict.fromkeys(preflight.cal.pst.adj_par_names, 1.0)
    values["kl0_typo"] = 2.0
    with pytest.raises(ValueError, match=r"names 1 parameter\(s\) the control file does not "
                                         r"hold: \['kl0_typo'\]"):
        preflight.cal.run_template(tmp_path / "w", values=values)



@pytest.mark.parametrize("bad", [float("nan"), float("inf"), np.float64("-inf")])
def test_a_value_that_is_not_finite_is_refused(preflight, tmp_path, bad):
    values = dict.fromkeys(preflight.cal.pst.adj_par_names, 1.0)
    values[_kl0(preflight.cal)] = bad
    with pytest.raises(ValueError, match=r"a parameter's value is a finite number"):
        preflight.cal.run_template(tmp_path / "w", values=values)
    assert not (tmp_path / "w").exists()


def test_an_observation_with_no_finite_simulated_value_is_refused(preflight, tmp_path):
    """PESTPP-IES fails a run whose output reads ``nan`` ("run failed"), where the phi
    by group summed over the table skipped the row, as if the template ran (ledger
    607; ``docs/reproducers/run_template_is_pestpp_once.py``, part E)."""

    cal = preflight.cal
    forward_run = cal.template_workspace / "forward_run.py"
    original = forward_run.read_text()
    output = str(next(name for name in cal.pst.output_files if "hd" in str(name)))
    hook = (f"\n    import pandas as _pd\n    _f = _pd.read_csv({output!r})\n"
            f"    _f.iloc[0, 1] = float('nan')\n    _f.to_csv({output!r}, index=False, "
            "na_rep='nan')\n")
    cut = original.index("\nif __name__")
    forward_run.write_text(original[:cut].rstrip("\n") + hook + "\n" + original[cut:])
    try:
        with pytest.raises(RuntimeError, match=r"1 of \d+ observation\(s\) read no finite "
                                               r"simulated value .* PESTPP-IES fails such a run"):
            cal.run_template(tmp_path / "w")
    finally:
        forward_run.write_text(original)
    assert (tmp_path / "w" / output).is_file()  # left for inspection


def test_one_parameter_named_twice_is_refused(preflight, tmp_path):
    name = _kl0(preflight.cal)
    values = {name: 1.0, name.upper(): 2.0}
    with pytest.raises(ValueError, match=r"names parameter .* twice"):
        preflight.cal.run_template(tmp_path / "w", values=values)


def test_values_of_another_type_are_refused(preflight, tmp_path):
    frame = pd.DataFrame({"parval1": [1.0]}, index=[_kl0(preflight.cal)])
    with pytest.raises(TypeError, match=r"values= takes None .* not DataFrame"):
        preflight.cal.run_template(tmp_path / "w", values=frame)


def test_before_build_is_refused(preflight, tmp_path):
    cal = preflight.model.pest("unbuilt", start_datetime="2020-01-01",
                               workspace=tmp_path / "unbuilt")
    with pytest.raises(ValueError, match=r"Call build\(\) before run_template\(\)"):
        cal.run_template(tmp_path / "w")


# -- one spelling: run_ies(iterations=0) is refused --------------------------------------


@pytest.mark.parametrize(("iterations", "error", "match"), [
    (0, ValueError, r"run_ies\(iterations=0\) .* cal\.run_template\("),
    (-1, ValueError, r"NOPTMAX -1 is a prior Monte Carlo, cal\.prior"),
    (2.5, TypeError, r"A number of iterations is a whole number"),
])
def test_run_ies_refuses_what_is_not_a_number_of_iterations(preflight, iterations, error, match):
    control = preflight.cal.template_workspace / "t.pst"
    before = control.stat().st_mtime_ns
    with pytest.raises(error, match=match):
        preflight.cal.run_ies(reals=4, iterations=iterations)
    assert control.stat().st_mtime_ns == before
    assert not (preflight.cal.template_workspace / "t.phi.actual.csv").exists()

"""What the ``ModelDiff`` results tier reads, measured while documenting its parameters.

Ledger 436.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/diff_results_tier.py

Needs ``mf6`` on PATH (measured with 6.8.0.dev0). About 20 s. Found while writing the
numpydoc entries of ``project/model_results_diff.py`` (2026-10-11); the docstrings state
these as current behaviour and none is fixed here.

THE SETUP. Three runs of the canonical model on its 21 x 21 testing grid (six 30-day
periods): ``a`` and ``b`` at one time step per period, ``c`` at two. And a 3 x 3 DIS
toy run twice: CHD on the west column, two drains (``drn_a`` at cell (0, 1, 2),
``drn_b`` at (0, 0, 2)) and a pumping well at (0, 2, 2), all three moving their water
by MVR into a one-reach SFR; the second run has ``drn_a``'s conductance 30 for 10.

WHAT IT MEASURES
  1. MVR. ``diff.packages.mvr.results`` reads records named exactly ``FROM-MVR`` /
     ``TO-MVR`` from the GROUNDWATER budget file. MF6 names a drain's and a well's
     mover term ``DRN-TO-MVR`` / ``WEL-TO-MVR`` there, and writes SFR's and LAK's
     ``FROM-MVR``/``TO-MVR`` into their own package budget files -- so every
     ``get(package=, direction=)`` and ``summary()`` is empty although the drains
     moved 69.5 and 122.8 volume units per day in the two runs (``drn_a`` 35.0
     against 91.3 of it) and the canonical model's SFR reach 55 feeds its lake.
  2. Heads. ``diff.hds`` reads each member's raw ``all_heads``: MF6's dry sentinel
     (-1e30, ledger 315) put in ONE cell of ``b`` -- a stand-in for a cell that dried,
     since the diff reads the table as it is -- makes ``summary()`` read
     ``max_abs_diff`` 1e30, ``mean_abs_diff`` ~9.4e25 and ``rmse`` ~9.7e27.
  3. LAK exchange. ``diff.packages.lak.results.q`` is joined on ``kstpkper``, not
     reduced to each member's period end: at ``per=1`` the ``nstp=2`` run is read at
     ``(0, 1)``, the FIRST step of period 1, against the reference's only step. GHB's
     ``q`` is reduced (the comparison carries no ``kstpkper``). The LAK join also
     keys on MF6's wetted ``FLOW-AREA``: a connection dry on one side only has no
     partner, and 2 of the 72 connections drop out at period 3 between ``a`` and
     ``c``.
  4. Budget. ``diff.bud`` differences every numeric column of FloPy's incremental
     listing budget, including ``time_step``, ``stress_period`` and ``tslen``: between
     ``a`` and ``c`` ``tslen`` reads -50% and ``time_step`` out of tolerance.
  5. Saturation. ``diff.packages.uzf.results.sat`` raises ``ValueError``: NPF's
     ``DATA-SAT`` record carries an always-zero ``q`` beside ``SAT``, and the table
     read renames ``q`` to ``sat`` next to the ``sat`` already there.
  6. Package names. A 3 x 3 toy run twice, with a DRN named ``seeps`` (two records
     on one cell, conductances 10/20 against 30/60) and a GHB named ``lake_edge``:
     ``diff.package_names`` lists only ``chd`` and ``connection_package_names``
     ``lak``, because both find a package by its NAME prefix, so ``report()`` reads
     "identical to reference"; ``packages.drn.inputs.summary()`` counts 4 changed
     value rows for the one cell.
"""

from __future__ import annotations

import contextlib
import io
import pathlib
import shutil
import tempfile
import warnings

import flopy

import myflopy as mf
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, build_canonical_model
from myflopy.modflow.mf6.package_budget import build_budget_result_table

warnings.simplefilter("ignore")
WORK = pathlib.Path(tempfile.mkdtemp(prefix="diff_results_tier_"))


def _records(path, precision="auto"):
    """``(text, paknam2)`` of every record in a budget file."""

    budget = flopy.utils.CellBudgetFile(str(path), precision=precision)
    return sorted({(r["text"].decode().strip(), r["paknam2"].decode().strip())
                   for r in budget.recordarray})


def _mover_toy(name, drn_a_cond):
    """The 3 x 3 toy with two drains and a well moving water into one SFR reach."""

    ws = WORK / name
    sim = flopy.mf6.MFSimulation(sim_name="mv", sim_ws=str(ws), exe_name="mf6",
                                 verbosity_level=0)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, complexity="SIMPLE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="mv", save_flows=True)
    flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=3, ncol=3, delr=10.0, delc=10.0,
                            top=10.0, botm=0.0)
    flopy.mf6.ModflowGwfnpf(gwf, k=10.0)
    flopy.mf6.ModflowGwfic(gwf, strt=9.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, r, 0), 9.0] for r in range(3)])
    flopy.mf6.ModflowGwfdrn(gwf, mover=True, pname="drn_a",
                            stress_period_data=[[(0, 1, 2), 5.0, drn_a_cond]])
    flopy.mf6.ModflowGwfdrn(gwf, mover=True, pname="drn_b", filename="mv_b.drn",
                            stress_period_data=[[(0, 0, 2), 5.0, 10.0]])
    flopy.mf6.ModflowGwfwel(gwf, mover=True, pname="wel",
                            stress_period_data=[[(0, 2, 2), -1.0]])
    flopy.mf6.ModflowGwfsfr(
        gwf, mover=True, pname="sfr", nreaches=1,
        packagedata=[[0, (0, 1, 1), 10.0, 1.0, 0.001, 8.0, 1.0, 0.0, 0.03, 0, 1.0, 0]],
        connectiondata=[[0]], perioddata=[[0, "INFLOW", 0.0]],
    )
    flopy.mf6.ModflowGwfmvr(
        gwf, maxmvr=3, maxpackages=4, packages=[["drn_a"], ["drn_b"], ["wel"], ["sfr"]],
        perioddata={0: [["drn_a", 0, "sfr", 0, "FACTOR", 1.0],
                        ["drn_b", 0, "sfr", 0, "FACTOR", 1.0],
                        ["wel", 0, "sfr", 0, "FACTOR", 1.0]]},
    )
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord="mv.cbc", head_filerecord="mv.hds",
                           saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])
    sim.write_simulation()
    ok, _ = sim.run_simulation(silent=True)
    assert ok, f"{name} did not run"
    return ws


def _canonical(name, nstp):
    """A run canonical model on the testing grid at ``nstp`` steps per period."""

    config = CanonicalModelConfig(nrow=21, ncol=21, nlay=4, nper=6, cell_size=100.0,
                                  steps_per_period=nstp)
    with contextlib.redirect_stdout(io.StringIO()):     # MF6's progress lines
        model = build_canonical_model(WORK / name, config=config)
        ok, _ = model.run_simulation()
    assert ok, f"{name} did not run"
    return mf.load_mf6_run(model.workspace)


def part_1_mover():
    print("PART 1 -- MVR")
    first, second = _mover_toy("mvr_a", 10.0), _mover_toy("mvr_b", 30.0)
    print("  groundwater budget records:", _records(second / "mv.cbc"))
    a, b = mf.load_mf6_run(first), mf.load_mf6_run(second)
    mvr = mf.ModelGroup({"a": a, "b": b}, reference="a").diff().packages.mvr.results
    for package in ("drn", "wel", "sfr"):
        for direction in ("from_mvr", "to_mvr"):
            rows = len(mvr.get(package=package, direction=direction))
            print(f"  mvr.get(package={package!r}, direction={direction!r}): {rows} rows")
    print("  mvr.summary():", len(mvr.summary()), "rows")
    for run, model in (("a", a), ("b", b)):
        moved = build_budget_result_table(model, budget_text="DRN-TO-MVR", package_name="drn")
        print(f"  run {run}: DRN-TO-MVR in the .cbc:",
              dict(zip(moved["pname"], moved["q"].round(3), strict=True)))


def part_2_to_4(a, b, c):
    print("PART 1b -- the canonical SFR-to-lake mover")
    ws = pathlib.Path(a.workspace)
    print("  groundwater budget:", [r for r in _records(next(ws.glob("*.cbc"))) if "MVR" in r[0]])
    for budget in sorted(ws.glob("*_budget.sfr")) + sorted(ws.glob("*_budget.lak")):
        print(f"  {budget.name}:", [r for r in _records(budget) if "MVR" in r[0]])

    print("PART 2 -- a dry cell's sentinel in diff.hds")
    table = b.all_heads.copy()
    table.iloc[5, 0] = -1e30

    class Dried:
        """``b`` with one saved head replaced by what MF6 writes for a dry cell."""

        def __init__(self, model):
            self._model = model

        def __getattr__(self, name):
            return getattr(self._model, name)

        @property
        def all_heads(self):
            return table

    summary = mf.ModelGroup({"a": a, "b": Dried(b)}, reference="a").diff().hds.summary()
    print(summary[["model", "max_abs_diff", "mean_abs_diff", "rmse", "within_tolerance"]]
          .to_string(index=False))

    diff = mf.ModelGroup({"a": a, "c": c}, reference="a").diff()
    print("PART 3 -- LAK exchange at nstp=2 against nstp=1, per=1")
    lak = diff.packages.lak.results.q.get(per=1)
    print("  kstpkper of the compared rows:", lak["kstpkper"].value_counts().to_dict())
    ghb = diff.packages.ghb.results.q.get(per=1)
    print("  GHB comparison carries kstpkper:", "kstpkper" in ghb.columns, f"({len(ghb)} rows)")
    paired = diff.packages.lak.results.q.get()["per"].value_counts().sort_index().to_dict()
    print("  LAK connections paired per period (72 connections):", paired)

    print("PART 4 -- FloPy's bookkeeping columns as budget terms")
    budget = diff.bud.summary()
    print(budget[budget["term"].isin(["stress_period", "time_step", "tslen"])]
          .to_string(index=False))

    print("PART 5 -- diff.packages.uzf.results.sat")
    try:
        diff.packages.uzf.results.sat.get()
        print("  get() returned")
    except ValueError as error:
        print(f"  get() raises ValueError: {error}")


def part_6_package_names():
    print("PART 6 -- packages found by NAME")

    def toy(name, cond):
        ws = WORK / name
        sim = flopy.mf6.MFSimulation(sim_name="t", sim_ws=str(ws), exe_name="mf6",
                                     verbosity_level=0)
        flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
        flopy.mf6.ModflowIms(sim, complexity="SIMPLE")
        gwf = flopy.mf6.ModflowGwf(sim, modelname="t", save_flows=True)
        flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=3, ncol=3, delr=10.0, delc=10.0,
                                top=10.0, botm=0.0)
        flopy.mf6.ModflowGwfnpf(gwf, k=10.0)
        flopy.mf6.ModflowGwfic(gwf, strt=9.0)
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, r, 0), 9.0] for r in range(3)])
        flopy.mf6.ModflowGwfdrn(gwf, pname="seeps", filename="t.drn", stress_period_data=[
            [(0, 1, 2), 5.0, cond], [(0, 1, 2), 5.0, 2.0 * cond]])
        flopy.mf6.ModflowGwfghb(gwf, pname="lake_edge", filename="t.ghb",
                                stress_period_data=[[(0, 0, 2), 8.0, 1.0]])
        flopy.mf6.ModflowGwfoc(gwf, budget_filerecord="t.cbc", head_filerecord="t.hds",
                               saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")])
        sim.write_simulation()
        ok, _ = sim.run_simulation(silent=True)
        assert ok, f"{name} did not run"
        return mf.load_mf6_run(ws)

    diff = mf.ModelGroup({"a": toy("names_a", 10.0), "b": toy("names_b", 30.0)},
                         reference="a").diff()
    print("  package_names:", diff.package_names,
          "connection_package_names:", diff.connection_package_names)
    print("  report():", diff.report().splitlines()[2])
    drn = diff.packages.drn.inputs.summary()
    print("  drn.inputs.summary() value_cells_changed:", int(drn["value_cells_changed"].iloc[0]),
          "for", int(drn["cells_shared"].iloc[0]), "shared cell")


def main():
    try:
        part_1_mover()
        a, b, c = _canonical("a", 1), _canonical("b", 1), _canonical("c", 2)
        part_2_to_4(a, b, c)
        part_6_package_names()
    finally:
        shutil.rmtree(WORK, ignore_errors=True)


if __name__ == "__main__":
    main()

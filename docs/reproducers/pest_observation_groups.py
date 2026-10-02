"""A target's group is its PEST observation group: measured through PEST++ itself.

Ledger 358. Every target class carried a ``group`` label that reached no PEST file:
pyEMU named each simulated column's observations a group of their own, so a group was
one well and PEST++'s phi by group could not see the user's grouping. This script
measures what the control file holds and what PEST++ reports, the PEST++ rules a group
name meets, and what the IES review draws.

    /home/lukem/python/envs/gw/bin/python docs/reproducers/pest_observation_groups.py

Sections A and C use only the public API, so they measure whichever tree is on
``PYTHONPATH`` -- the "before" column below is this script run against the parent
commit's ``src`` (``git archive 592d430 src | tar -x -C <dir>``, then
``PYTHONPATH=<dir>/src``). Section B edits a control file directly and measures PEST++
alone, so it reads the same on both. Needs ``mf6``, ``pestpp-glm`` and ``pestpp-ies``.

THE SETUP. ``tests/test_pest_handed_targets.py``'s model: a 6 x 6 irregular Voronoi grid
of ~100 ft cells, one confined layer, CHD at 102 ft west and east, heads filling from
90 ft, a four-reach stream, a three-cell DRN zone and a two-cell lake; three 10-day
periods of three steps, ``SAVE HEAD LAST`` + ``SAVE BUDGET LAST``. One parameter, ``k``.
Heads W2, W1 (labelled ``North``) and W10 (``South``) measured at 91-97 ft, W2 not in
period 2 and W1 not in period 1; the DRN zone labelled ``Seeps``; the lake unlabelled;
a forecast P labelled ``North``, with a placeholder value.

  A. The build, then ``pestpp-glm`` NOPTMAX=0: the groups the control file holds and
     the phi PEST++ reports for each (its ``.rec``).
  B. PEST++ 5.2.16 on the heads alone, every observation put in one group named as
     listed, ``pestpp-ies`` NOPTMAX=-1 (3 reals) unless named: the base realization's
     phi for that group. The simulated heads (~92-102 ft) sit ABOVE the targets, so a
     "simulated at least" group scores 0; with every target set to 120 ft (above them)
     a "simulated at most" group does. Then the characters and lengths a group name
     may hold, and how ``ies_phi_factor_file`` reads groups (NOPTMAX=1, 4 reals).
  C. ``cal.prior(reals=3)`` and the review's ``plot_vs_obs``: its panels, and the
     measured markers of each; then the same run with its groups set as this tree
     writes them (W2 and W1 ``north``, W10 ``south``), reviewed by the tree on
     ``PYTHONPATH``.
  D. (Added with ledger 359.) The heads, a lake and an SFR stage set, all UNLABELLED, so
     their groups are the default prefixes ``hds``, ``stage`` and ``sfr_stage``; a
     phi-factor file tagging each (``pestpp-ies`` NOPTMAX=1, 4 reals). Then the same
     with the lake set registered under ``prefix="lake"``.
  E. (Added with ledger 359.) W2/W1 labelled ``Regulatory wells`` -- PEST names a
     regularization group ``regul...`` -- through pestpp-glm NOPTMAX=0 and pestpp-ies
     NOPTMAX=-1; then ``cal.write()``, which lets pyEMU pick the file version, holding a
     219-character group and one with ``.`` and ``:``, through pestpp-glm.

MEASURED (mf6 6.8.0.dev0, PEST++ 5.2.16, pyEMU 1.4.0, 2026-10-02)
                                                 before (592d430)      after
  A  groups of the measured observations      6, one per column:     5: north (4 obs),
                                              w2, w1, w10, zone,     south (3), seeps (1),
                                              pond, a (3 obs each)   stage (1), sfr_flow (1)
     of those, also holding zero-weight rows  5 of 6                 0 of 5
     groups at weight 0 only                  fore1's p              hds_, drn_flow_,
                                                                     stage_, sfr_flow_,
                                                                     fore1_unmeasured; fore1
     phi by group, pestpp-glm .rec            w2 201.134, w1         north 324.021 (=
                                              122.887, w10 100.843,  201.134 + 122.887),
                                              zone 741302., pond     south 100.843, seeps
                                              121.000, a 5495.09     741302., stage 121.000,
                                                                     sfr_flow 5495.09
  B  (PEST++ alone, the same on both trees) pestpp-ies, heads ABOVE their targets:
     north 424.864; greater_n, greaterville, g_n 0 -- read as "simulated at least the
     target" -- while gt_n is 424.864. Every target at 120 ft: north 2377.88; less_n,
     lessard, l_n 0 -- "at most" -- while ln, lt_n, l-n, north_less are 2377.88.
     pestpp-glm NOPTMAX=0 reads greater_n and less_n as ordinary groups (424.864 each).
     A comma or a double quote in a group name: pestpp-glm exit 1, "Error condition
     prevents further execution"; a space, ``/``, ``#``, ``:``, ``'``, ``!`` and 5000
     characters: read, exit 0. ``ies_phi_factor_file``: tags "orth"/"outh" took north
     and south (a tag matches anywhere in a name); "north" alone: exit 1, "group 'SOUTH'
     not identified with any tags"; "north"/"south"/"nor": exit 1, "group 'NORTH' mapped
     to multiple tags: NOR NORTH"; north/south tagged beside an untagged zero-weight
     hds_unmeasured: exit 0.
  C  plot_vs_obs panels: their measured markers
     as built by the tree                     3, one per column:     3, one per well:
                                              w2 91, 92, 101.962;    W2 91, 92;
                                              w1 93, 101.847, 94;    W1 93, 94;
                                              w10 95, 96, 97         W10 95, 96, 97
     the same run grouped north/south         2: north 93, 91, 92,   the same 3 as above
                                              94 (one line through
                                              two wells); south 95,
                                              96, 97

  D  (this tree, ledger 359) tags hds/stage/sfr_stage on the default groups: exit 1,
     "group 'SFR_STAGE' mapped to multiple tags: SFR_STAGE STAGE" -- every tag that
     matches ``stage`` is inside ``sfr_stage`` too, so the lake's group cannot be given
     a fraction of its own; with the lake set under prefix "lake", tags hds/lake/sfr_stage:
     exit 0.
  E  (this tree) regulatory_wells 324.021 and south 100.843 in pestpp-glm's .rec and in
     pestpp-ies's base phi.group.csv -- an ordinary group in estimation mode, the same
     phi as `north` in A. cal.write() wrote VERSION 1; pestpp-glm exit 0, the 219-character
     group 100.843 and north.wells:a 324.021.

"before" in C draws every observation of a group as a red "measured" triangle: 101.962
and 101.847 are the base run's own heads in the period W2 and W1 were NOT measured,
weight 0 since ledger 340.
Why each is so, and what was decided: ledger 358.
"""

from __future__ import annotations

import pathlib
import re
import shutil
import subprocess
import tempfile
import warnings

import numpy as np
import pandas as pd
import pyemu

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.modflow.mf6.pest.ies import IesResults

NPER, NSTP, PERLEN = 3, 3, 10.0
REACHES = [13, 14, 15, 16]
DRN_CELLS = [20, 21, 27]
LAKE_CELLS = [8, 9]
HOME = pathlib.Path.home() / ".local/bin"
GLM = shutil.which("pestpp-glm") or str(HOME / "pestpp-glm")
IES = shutil.which("pestpp-ies") or str(HOME / "pestpp-ies")
#: pestpp-glm's ``.rec``: 'Contribution to phi from observation group "north" : 324.021'.
CONTRIBUTION = re.compile(r'observation group\s+"(.+?)"\s*:\s*(\S+)')
HEADS = pd.DataFrame(
    {"time": [0, 1, 2], "W2": [91.0, 92.0, None], "W1": [93.0, None, 94.0], "W10": [95.0, 96.0, 97.0]}
)


def _toy(root):
    grid = irregular_voronoi_grid(nrow=6, ncol=6, cell_size=100.0, refined=False)
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    xmin, _, xmax, _ = grid.gdf_vorPolys.total_bounds
    edge = np.flatnonzero((xy[:, 0] < xmin + 100) | (xy[:, 0] > xmax - 100))
    last = len(REACHES) - 1
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]),
                vertices=props["vertices"], cell2d=props["cell2d"], top=110.0, botm=[0.0]),
        mf.npf(k=10.0, icelltype=0),
        mf.sto(ss=1e-3, sy=0.2, iconvert=0, transient={0: True}), mf.ic(strt=90.0),
        mf.chd(stress_period_data={0: [[(0, int(c)), 102.0] for c in edge]}),
        mf.sfr.flopy(
            packagedata=[[i, (0, c), 100.0, 5.0, 0.001, 95.0 - 0.1 * i, 1.0, 1.0, 0.03,
                          1 if i in (0, last) else 2, 1.0, 0] for i, c in enumerate(REACHES)],
            connectiondata=[[i, *([i - 1] if i > 0 else []), *([-(i + 1)] if i < last else [])]
                            for i in range(len(REACHES))],
            perioddata={0: [[0, "INFLOW", 50.0]], 1: [[0, "INFLOW", 200.0]],
                        2: [[0, "INFLOW", 20.0]]},
            stage_filerecord="m.sfr.stg", budget_filerecord="m.sfr.bud"),
        mf.drn(stress_period_data={0: [[(0, c), 96.0, 50.0] for c in DRN_CELLS]}),
        mf.lak.flopy(
            nlakes=1, packagedata=[[0, 104.0, len(LAKE_CELLS), "pond"]],
            connectiondata=[[0, i, (0, c), "VERTICAL", 0.5, 0.0, 0.0, 0.0, 0.0]
                            for i, c in enumerate(LAKE_CELLS)],
            perioddata={0: [[0, "WITHDRAWAL", 10.0]]}, boundnames=True,
            stage_filerecord="m.lak.stg"),
        mf.oc(head_filerecord="m.hds", budget_filerecord="m.cbc",
              saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")])])
    simulation = mf.simulation(model, name="s", packages=[
        mf.tdis(nper=NPER, perioddata=[(PERLEN, NSTP, 1.0)] * NPER)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = mf.Project(root / "p", name="p").run("r", simulation)
    assert run.success
    return run.model()


def _heads(groups):
    return mf.HeadTargets(
        locations=pd.DataFrame({"name": ["W2", "W1", "W10"], "cell": [1, 2, 34], "layer": 0,
                                "group": groups}),
        values=HEADS)


def _calibration(model, workspace, name):
    cal = model.pest(name, workspace=workspace, start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
    return cal


def _built(cal, name):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return cal.build(f"{name}.pst", noptmax=0)


def section_a(model, root):
    cal = _calibration(model, root / "a", "a")
    cal.observe(_heads(["North", "North", "South"]))
    cal.observe(mf.DrnFlowTargets(
        locations=pd.DataFrame({"name": ["zone"], "cells": [DRN_CELLS], "group": ["Seeps"]}),
        values=pd.DataFrame({"time": [2.0], "zone": [-3.0]})))
    cal.observe(mf.LakeStageTargets(locations={"pond": 0},
                                    values=pd.DataFrame({"time": [1.0], "pond": [99.0]})))
    cal.observe(mf.SfrFlowTargets(locations={"a": 0},
                                  values=pd.DataFrame({"time": [1.0], "a": [-7.0]})))
    cal.forecast(mf.HeadTargets(
        locations=pd.DataFrame({"name": ["P"], "cell": [23], "layer": 0, "group": ["North"]}),
        values=pd.DataFrame({"time": [2], "P": [100.0]})))
    obs = _built(cal, "a").observation_data
    weighted = obs.loc[obs["weight"] > 0]
    measured = sorted(weighted["obgnme"].unique())
    zero_only = sorted(set(obs["obgnme"]) - set(measured))
    mixed = [g for g in measured if (obs.loc[obs["obgnme"] == g, "weight"] == 0).any()]
    print(f"A  groups of the measured observations ({len(measured)}):")
    for group in measured:
        print(f"     {group}: {len(obs.loc[obs['obgnme'] == group])} observation(s)")
    print(f"   groups at weight 0 only: {zero_only}")
    print(f"   measured groups also holding zero-weight rows: {len(mixed)} of {len(measured)}")
    done = subprocess.run([GLM, "a.pst"], cwd=cal.template_workspace, capture_output=True,
                          text=True, timeout=600)
    rec = (cal.template_workspace / "a.rec").read_text()
    contributions = CONTRIBUTION.findall(rec)
    print(f"   pestpp-glm exit {done.returncode}; phi by group:")
    for group, phi in contributions:
        print(f"     {group}: {phi}")


def _phi_of(template, label, groups, *, obsval=None, exe=IES, options=None, noptmax=None):
    """Copy ``template``, set every group (and optionally every target), run PEST++."""

    work = template.parent / f"b_{label}"
    if work.exists():
        shutil.rmtree(work)
    shutil.copytree(template, work)
    pst = pyemu.Pst(str(work / "b.pst"))
    pst.observation_data["obgnme"] = groups
    if obsval is not None:
        pst.observation_data["obsval"] = obsval
    pst.control_data.noptmax = noptmax if noptmax is not None else (-1 if exe == IES else 0)
    if exe == IES:
        pst.pestpp_options["ies_num_reals"] = 3
    pst.pestpp_options.update(options or {})
    pst.write(str(work / "b.pst"), version=2)
    done = subprocess.run([exe, "b.pst"], cwd=work, capture_output=True, text=True, timeout=600)
    errors = [line.strip() for line in (done.stdout + done.stderr).splitlines()
              if "rror" in line and "tmp file" not in line]
    return work, done.returncode, errors


def section_b(model, root):
    cal = _calibration(model, root / "b", "b")
    cal.observe(_heads([None, None, None]))
    _built(cal, "b")
    template = pathlib.Path(cal.template_workspace)

    def ies_base(work, group):
        table = pd.read_csv(work / "b.phi.group.csv")
        return table.loc[table["par_realization"].astype(str) == "base", group].iloc[0]

    print("B  pestpp-ies, base realization's phi, every observation in one group:")
    for group in ("north", "greater_n", "greaterville", "g_n", "gt_n"):
        work, code, _ = _phi_of(template, group, group)
        print(f"     {group}: {ies_base(work, group)}")
    print("   ... with every target at 120 ft (above every simulated head):")
    for group in ("north", "less_n", "lessard", "l_n", "ln", "lt_n", "l-n", "north_less"):
        work, code, _ = _phi_of(template, f"{group}_120", group, obsval=120.0)
        print(f"     {group}: {ies_base(work, group)}")
    print("   pestpp-glm NOPTMAX=0, phi of greater_n / less_n (targets as built):")
    for group in ("greater_n", "less_n"):
        work, code, _ = _phi_of(template, f"{group}_glm", group, exe=GLM)
        rec = (work / "b.rec").read_text()
        print(f"     {group}: {CONTRIBUTION.findall(rec)}")
    print("   characters and lengths (pestpp-glm NOPTMAX=0):")
    for label, group in (("comma", "north,wells"), ("dquote", 'north"wells'), ("space", "north wells"),
                         ("slash", "north/wells"), ("hash", "north#wells"), ("colon", "north:wells"),
                         ("quote", "north'wells"), ("bang", "north!wells"), ("long", "a" * 5000)):
        _work, code, errors = _phi_of(template, label, group, exe=GLM)
        print(f"     {label}: exit {code} {errors[:1]}")
    pst = pyemu.Pst(str(template / "b.pst"))
    split = ["north" if re.search(r"usecol:w(1|2)_", name) else "south"
             for name in pst.observation_data.index]
    zero = [g if w > 0 else "hds_unmeasured"
            for g, w in zip(split, pst.observation_data["weight"], strict=True)]
    print("   ies_phi_factor_file (NOPTMAX=1, 4 reals):")
    for label, text, groups in (("substring", "orth,0.8\nouth,0.2\n", split),
                                ("one_tag", "north,1.0\n", split),
                                ("two_tags", "north,0.5\nsouth,0.3\nnor,0.2\n", split),
                                ("zero_untagged", "north,0.8\nsouth,0.2\n", zero)):
        (template / "phifac.csv").write_text(text)
        work, code, _errors = _phi_of(
            template, f"phi_{label}", groups, noptmax=1,
            options={"ies_phi_factor_file": "phifac.csv", "ies_num_reals": 4})
        rec = (work / "b.rec").read_text()
        said = sorted({line.strip() for line in rec.splitlines()
                       if "tag" in line.lower() and ("rror" in line or "adjusted" in line)})
        print(f"     {label}: exit {code} {said[:2]}")


def _panels(results):
    figure = results.plot_vs_obs()
    titles = [annotation.text for annotation in figure.layout.annotations]
    measured = [trace for trace in figure.data if trace.name == "measured"]
    return [(title, [round(float(y), 3) for y in trace.y]) for title, trace in zip(titles, measured, strict=True)]


def section_c(model, root):
    cal = _calibration(model, root / "c", "c")
    cal.observe(_heads(["North", "North", "South"]))
    _built(cal, "c")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        results = cal.prior(reals=3)
    print("C  plot_vs_obs panels, and each one's measured markers:")
    for title, values in _panels(results):
        print(f"     {title}: {values}")
    grouped = pathlib.Path(results.workspace).parent / "c_grouped"
    if grouped.exists():
        shutil.rmtree(grouped)
    shutil.copytree(results.workspace, grouped)
    pst = pyemu.Pst(str(grouped / "c.pst"))
    obs = pst.observation_data
    obs["obgnme"] = [
        ("north" if re.search(r"usecol:w(1|2)_", name) else "south") if weight > 0 else "hds_unmeasured"
        for name, weight in zip(obs.index, obs["weight"], strict=True)
    ]
    pst.write(str(grouped / "c.pst"), version=2)
    print("   the same run grouped north/south, reviewed:")
    for title, values in _panels(IesResults(grouped, case_name="c")):
        print(f"     {title}: {values}")


def section_d(model, root):
    """The DEFAULT prefixes as phi-factor groups: a lake set's ``stage`` is inside an SFR
    set's ``sfr_stage``, so no tag can match the lake's group alone (ledger 359)."""

    print("D  ies_phi_factor_file on the default prefixes (pestpp-ies NOPTMAX=1, 4 reals):")
    for label, lake_prefix, text in (
        ("defaults", None, "hds,0.5\nstage,0.25\nsfr_stage,0.25\n"),
        ("lake_renamed", "lake", "hds,0.5\nlake,0.25\nsfr_stage,0.25\n"),
    ):
        cal = _calibration(model, root / f"d_{label}", "d")
        cal.observe(_heads([None, None, None]))
        cal.observe(mf.LakeStageTargets(locations={"pond": 0}, values=pd.DataFrame(
            {"time": [1.0, 2.0], "pond": [99.0, 98.0]})), prefix=lake_prefix)
        cal.observe(mf.SfrStageTargets(locations={"r1": 1}, values=pd.DataFrame(
            {"time": [0.0, 1.0], "r1": [95.0, 96.0]})))
        pst = _built(cal, "d")
        work = pathlib.Path(cal.template_workspace)
        (work / "phifac.csv").write_text(text)
        pst.control_data.noptmax = 1
        pst.pestpp_options.update({"ies_num_reals": 4, "ies_phi_factor_file": "phifac.csv"})
        pst.write(str(work / "d.pst"), version=2)
        done = subprocess.run([IES, "d.pst"], cwd=work, capture_output=True, text=True, timeout=900)
        rec = (work / "d.rec").read_text() if (work / "d.rec").exists() else ""
        said = sorted({line.strip() for line in (rec + done.stdout).splitlines()
                       if "tags" in line and "rror" in line})
        groups = sorted(set(pst.observation_data.loc[pst.observation_data["weight"] > 0, "obgnme"]))
        print(f"     {label}: measured groups {groups}; exit {done.returncode} {said[:1]}")


def section_e(model, root):
    """Two more names the build writes, through PEST++ (ledger 359): a label starting
    PEST's ``regul`` prefix, and ``cal.write()``'s version-1 file, which pyEMU picks for a
    small problem, holding a long group with ``.`` and ``:``."""

    cal = _calibration(model, root / "e_regul", "e")
    cal.observe(_heads(["Regulatory wells", "Regulatory wells", "South"]))
    _built(cal, "e")
    work = pathlib.Path(cal.template_workspace)
    print("E  a 'regul...' group (estimation mode), base phi by group:")
    for exe in (GLM, IES):
        pst = pyemu.Pst(str(work / "e.pst"))
        pst.control_data.noptmax = -1 if exe == IES else 0
        pst.pestpp_options["ies_num_reals"] = 3
        pst.write(str(work / "e.pst"), version=2)
        done = subprocess.run([exe, "e.pst"], cwd=work, capture_output=True, text=True, timeout=600)
        if exe == IES:
            table = pd.read_csv(work / "e.phi.group.csv")
            base = table.loc[table["par_realization"].astype(str) == "base"].iloc[0]
            found = [(group, round(float(base[group]), 3)) for group in ("regulatory_wells", "south")]
        else:
            found = CONTRIBUTION.findall((work / "e.rec").read_text())
        print(f"     {pathlib.Path(exe).name}: exit {done.returncode} {found}")
    cal = _calibration(model, root / "e_v1", "v")
    cal.observe(_heads(["North.Wells:A", "North.Wells:A", "Very Long Group Label " * 10]))
    _built(cal, "v")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        written = cal.write("v1.pst")
    version = "version 2" if "* observation data external" in written.read_text() else "version 1"
    done = subprocess.run([GLM, "v1.pst"], cwd=cal.template_workspace, capture_output=True,
                          text=True, timeout=600)
    found = CONTRIBUTION.findall((pathlib.Path(cal.template_workspace) / "v1.rec").read_text())
    print(f"   cal.write() -> {version}; pestpp-glm exit {done.returncode}; "
          f"{[(group if len(group) < 30 else f'{len(group)} chars', phi) for group, phi in found]}")


def main():
    root = pathlib.Path(tempfile.mkdtemp(prefix="pest_groups_"))
    try:
        model = _toy(root)
        section_a(model, root)
        section_b(model, root)
        section_c(model, root)
        section_d(model, root)
        section_e(model, root)
    finally:
        shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    main()

"""What MF6 does with a DISV ``LENGTH_UNITS`` word (ledger 404).

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/mf6_length_units.py

Needs `mf6` on PATH (measured with 6.8.0.dev0, FloPy 3.10). About 3 s.

PART 1 -- MF6 itself. A two-cell FloPy DISV model (one CHD, steady state) written with
each word and run; the listing's "Model length unit [0=UND, 1=FEET, 2=METERS,
3=CENTIMETERS] set as N" line is read back. Measured:

    FEET 1   feet 1   METERS 2   UNKNOWN 0
    FT 0     ft 0     M 0        FOO 0

Every run terminates normally and nothing in stdout or the listing names the word: MF6
reads any word it does not know as 0, undefined.

PART 2 -- what myflopy writes for a stack declared in an abbreviation. Before the fix
``LayerStack(length_units="ft").to_disv()`` wrote ``FT`` (and ``"m"`` wrote ``M``),
so the run reported unit 0; it writes ``FEET``/``METERS`` now and the run reports 1/2.

PART 3 -- does MF6 COMPUTE anything from the unit? (ledger 432) The canonical
testing model (every package, five OBS families) is written once, copied, its DISV's
``LENGTH_UNITS feet`` rewritten to ``FT`` in the copy, and both run. Measured: unit 1
against unit 0, and all 21 output files -- heads, the cell budget, the LAK/SFR/UZF/MVR
budgets and stages, every OBS CSV, the solver CSVs and the ``.grb`` -- BYTE-IDENTICAL.
MF6's source agrees: ``lenuni`` is read by the listing line and the NetCDF exporters
alone. So the unit is a label, and a wrong word moves no number.

PART 4 -- a spec that HOLDS an abbreviation (a saved project, a frozen run, or one
edited by ``with_options`` before the edit routes were checked) is now written as MF6's
word by the DIS builder: ``PackageSpec(..., {"length_units": "ft"})`` is handed to FloPy
as ``FEET`` and written ``LENGTH_UNITS feet`` (FloPy writes its enum in lower case), so a
rebuilt run reports unit 1 where it ran under 0, with the same numbers (PART 3).
About 15 s in all (PART 3 builds the canonical model).
"""

from __future__ import annotations

import pathlib
import re
import shutil
import subprocess
import tempfile

import flopy

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid

UNIT_LINE = re.compile(r"Model length unit \[0=UND, 1=FEET, 2=METERS, 3=CENTIMETERS\] set as\s+(\d)")


def _run(workspace: pathlib.Path) -> int:
    run = subprocess.run(["mf6"], cwd=workspace, capture_output=True, text=True, check=False)
    assert "Normal termination" in run.stdout, run.stdout[-1500:]
    return int(UNIT_LINE.search((workspace / "m.lst").read_text()).group(1))


def part_1(root: pathlib.Path) -> None:
    for word in ["FEET", "feet", "METERS", "UNKNOWN", "FT", "ft", "M", "FOO"]:
        workspace = root / f"flopy_{word}"
        sim = flopy.mf6.MFSimulation(sim_ws=str(workspace), exe_name="mf6")
        flopy.mf6.ModflowTdis(sim)
        flopy.mf6.ModflowIms(sim)
        gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
        flopy.mf6.ModflowGwfdisv(
            gwf, length_units=word, nlay=1, ncpl=2, nvert=6, top=10.0, botm=0.0,
            vertices=[[0, 0.0, 0.0], [1, 1.0, 0.0], [2, 2.0, 0.0],
                      [3, 0.0, 1.0], [4, 1.0, 1.0], [5, 2.0, 1.0]],
            cell2d=[[0, 0.5, 0.5, 4, 0, 3, 4, 1], [1, 1.5, 0.5, 4, 1, 4, 5, 2]],
        )
        flopy.mf6.ModflowGwfnpf(gwf)
        flopy.mf6.ModflowGwfic(gwf, strt=5.0)
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0), 5.0]])
        sim.write_simulation(silent=True)
        print(f"PART 1  LENGTH_UNITS {word:<8} -> unit {_run(workspace)}")


def part_2(root: pathlib.Path) -> None:
    vor = irregular_voronoi_grid(nrow=2, ncol=2, cell_size=100.0)
    for units in ["ft", "m", "cm"]:
        disv = (mf.LayerStack(top=100.0, length_units=units).add("a", thickness=10.0)
                .to_disv(vor, attach=False))
        model = mf.gwf("m", packages=[disv, mf.npf(k=1.0), mf.ic(strt=95.0),
                                      mf.chd(stress_period_data={0: [[(0, 0), 95.0]]})])
        sim = mf.simulation(model, name="s",
                            packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
        workspace = root / f"stack_{units}"
        sim.build_flopy(workspace).simulation.write_simulation(silent=True)
        print(f"PART 2  LayerStack(length_units={units!r}) writes "
              f"{disv.options['length_units']!r} -> unit {_run(workspace)}")


def part_3(root: pathlib.Path) -> None:
    import filecmp
    import warnings

    from myflopy.modflow.mf6.canonical_example import (
        CanonicalModelConfig,
        build_canonical_model,
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = build_canonical_model(root / "canonical", config=CanonicalModelConfig.testing())
        model.sim.write_simulation(silent=True)
    source = pathlib.Path(model.sim.sim_path)
    inputs = {path.name for path in source.iterdir()}
    runs = {}
    for word in ["feet", "FT"]:
        workspace = root / f"canonical_{word}"
        shutil.copytree(source, workspace)
        for disv in workspace.glob("*.disv"):
            disv.write_text(re.sub(r"LENGTH_UNITS\s+\S+", f"LENGTH_UNITS  {word}",
                                   disv.read_text()))
        run = subprocess.run(["mf6"], cwd=workspace, capture_output=True, text=True, check=False)
        assert "Normal termination" in run.stdout, run.stdout[-1500:]
        listing = next(path for path in workspace.glob("*.lst") if path.name != "mfsim.lst")
        unit = UNIT_LINE.search(listing.read_text()).group(1)
        print(f"PART 3  canonical model under LENGTH_UNITS {word:<4} -> unit {unit}")
        runs[word] = workspace
    outputs = sorted(path.name for path in runs["feet"].iterdir()
                     if path.is_file() and path.suffix != ".lst" and path.name not in inputs)
    differ = [name for name in outputs
              if not filecmp.cmp(runs["feet"] / name, runs["FT"] / name, shallow=False)]
    print(f"PART 3  {len(outputs)} output files compared, {len(differ)} differ: {differ}")


def part_4(root: pathlib.Path) -> None:
    vor = irregular_voronoi_grid(nrow=2, ncol=2, cell_size=100.0)
    declared = mf.LayerStack(top=100.0).add("a", thickness=10.0).to_disv(vor, attach=False)
    held = mf.PackageSpec(declared.name, declared.builder,
                          {**declared.options, "length_units": "ft"})
    model = mf.gwf("m", packages=[held, mf.npf(k=1.0), mf.ic(strt=95.0),
                                  mf.chd(stress_period_data={0: [[(0, 0), 95.0]]})])
    sim = mf.simulation(model, name="s",
                        packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    workspace = root / "held_ft"
    sim.build_flopy(workspace).simulation.write_simulation(silent=True)
    written = re.search(r"LENGTH_UNITS\s+(\S+)", (workspace / "m.disv").read_text()).group(1)
    print(f"PART 4  a spec holding 'ft' writes {written!r} -> unit {_run(workspace)}")


if __name__ == "__main__":
    root = pathlib.Path(tempfile.mkdtemp(prefix="mf6_length_units_"))
    try:
        part_1(root)
        part_2(root)
        part_3(root)
        part_4(root)
    finally:
        shutil.rmtree(root, ignore_errors=True)

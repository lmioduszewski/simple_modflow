"""What MF6 does with a TDIS ``TIME_UNITS`` word, and what myflopy reads from one (ledger 442).

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/tdis_time_units.py

Needs `mf6` on PATH (measured with 6.8.0.dev0, FloPy 3.10). About 3 s.

PART 1 -- MF6 itself. A two-cell FloPy DIS model (one CHD, steady state) written with
each TDIS word and run. Measured:

    unknown  UNKNOWN  seconds  DAYS  Hours  years   -> Normal termination
                                                       ("SIMULATION TIME UNIT IS ...")
    undefined  d  day  none  foo                    -> MF6 STOPS: "Unrecognized input
                                                       value for TIME_UNITS option"

FloPy writes every one of them as given (it lower-cases it). So unlike LENGTH_UNITS,
which MF6 reads as undefined whatever the word (ledger 404,
``mf6_length_units.py``), a TIME_UNITS word outside MF6's six stops the run.

PART 2 -- the one reader, ``hydraulics.declared_time_unit``, and the labels made from
it, for each case a FloPy TDIS can hold (FloPy itself refuses a blank one with an
``IndexError``): unset, ``unknown`` and a word MF6 refuses all read as no unit; each of
the five units reads as itself, in any case.

PART 3 -- the deprecated ``mp3du.legacy_prt.PRT`` copies FloPy's
``modeltime.time_units`` into its TDIS, and FloPy answers ``"unknown"`` for a flow TDIS
that names no unit, so it writes ``TIME_UNITS unknown`` -- which PART 1 shows MF6
reads. It writes nothing MF6 refuses for such a model, so that tier is left as it is.
"""

from __future__ import annotations

import pathlib
import shutil
import subprocess
import tempfile
import types
import warnings

import flopy

from myflopy.modflow.mf6.hydraulics import declared_time_unit, unit_labels
from myflopy.modflow.mf6.package_budget import budget_value_units


def _flow_sim(workspace: pathlib.Path, **tdis) -> flopy.mf6.MFSimulation:
    sim = flopy.mf6.MFSimulation(sim_ws=str(workspace), exe_name="mf6", verbosity_level=0)
    flopy.mf6.ModflowTdis(sim, **tdis)
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m")
    flopy.mf6.ModflowGwfdis(gwf, nrow=1, ncol=2, length_units="feet")
    flopy.mf6.ModflowGwfnpf(gwf)
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[((0, 0, 0), 1.0)])
    return sim


def _mf6_verdict(workspace: pathlib.Path) -> str:
    run = subprocess.run(["mf6"], cwd=workspace, capture_output=True, text=True, check=False)
    listing = (workspace / "mfsim.lst").read_text()
    if "Normal termination" in listing:
        unit = [line.strip() for line in listing.splitlines() if "SIMULATION TIME UNIT" in line]
        return f"normal ({unit[0] if unit else 'no unit line'})"
    refused = [line.strip() for line in run.stdout.splitlines() if "TIME_UNITS" in line]
    return f"STOPPED ({refused[0] if refused else run.stdout[-200:]})"


def _tdis_line(workspace: pathlib.Path) -> str:
    lines = next(workspace.glob("*.tdis")).read_text().splitlines()
    return next((line.strip() for line in lines if "TIME_UNITS" in line.upper()), "(none)")


def part_1(root: pathlib.Path) -> None:
    for word in ["unknown", "UNKNOWN", "seconds", "DAYS", "Hours", "years",
                 "undefined", "d", "day", "none", "foo"]:
        workspace = root / f"part1_{word}"
        _flow_sim(workspace, time_units=word).write_simulation(silent=True)
        print(f"PART 1  TIME_UNITS {word:<9} writes {_tdis_line(workspace)!r:<24} "
              f"-> {_mf6_verdict(workspace)}")


def part_2(root: pathlib.Path) -> None:
    for word in [None, "unknown", "d", "DAYS", " Hours ", "years"]:
        sim = _flow_sim(root / "part2", **({} if word is None else {"time_units": word}))
        model = types.SimpleNamespace(sim=sim, gwf=sim.get_model("m"), model_type="gwf6")
        print(f"PART 2  TDIS {str(word)!r:<10} declared_time_unit={declared_time_unit(model)!r:<10} "
              f"unit_labels={unit_labels(model)}  budget={budget_value_units(model)}")


def part_3(root: pathlib.Path) -> None:
    from myflopy.modflow.mp3du.legacy_prt import PRT

    for word in [None, "days"]:
        flow = _flow_sim(root / "part3_flow", **({} if word is None else {"time_units": word}))
        model = types.SimpleNamespace(
            name="m", gwf=flow.get_model("m"), model_output_folder_path=root / "part3_flow",
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            legacy = PRT(model=model, name="legacy", mf_folder_path=root)
        legacy.tdis.write()
        written = _tdis_line(legacy.model_output_folder_path)
        # The legacy PRT simulation itself is not runnable (it pairs PRT with IMS), so
        # its TDIS is checked under a flow model, where MF6 parses TDIS the same way.
        check = root / f"part3_check_{word}"
        sim = _flow_sim(check)
        sim.write_simulation(silent=True)
        shutil.copy(next(legacy.model_output_folder_path.glob("*.tdis")), next(check.glob("*.tdis")))
        print(f"PART 3  legacy PRT, flow TDIS {str(word)!r:<7} writes {written!r:<24} "
              f"-> {_mf6_verdict(check)}")


if __name__ == "__main__":
    root = pathlib.Path(tempfile.mkdtemp(prefix="tdis_time_units_"))
    try:
        part_1(root)
        part_2(root)
        part_3(root)
    finally:
        shutil.rmtree(root, ignore_errors=True)

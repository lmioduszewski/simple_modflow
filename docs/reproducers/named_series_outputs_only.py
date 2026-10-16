"""The PEST forward run's zone sets read the run's OUTPUTS alone: no model load per realization.

Ledger 582. ``write_named_series_targets`` read every DRN/RIV zone set
through ``mf.load_mf6_run`` -- one full FloPy load per realization (ledger 540:
8.3-9.0 s on Ten Trails ``ss`` beside a 17-20 s MF6 run). Given ``outputs=`` (the
build's ``named_series_outputs.json``) the zone sets read the run's budget file and
OBS output through ``pest.run_outputs.RunOutputs`` instead.

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/named_series_outputs_only.py <run copy> [repeats]

``<run copy>`` is a COPY of a completed Ten Trails run directory (the user's runs are
read only). The script builds a PEST template beside it (``<run copy>.pest/outputs``)
observing the sets of ledger 540's reproducer -- a RivFlowTargets of two zones (every
river cell of ``riv`` and of ``crisp``, ``package=None``), a ``crisp``-only
RivFlowTargets and a DRN zone of every drain cell -- copies the run's OUTPUT files
into the template, then times the forward run's named-series step ``repeats`` times
(default 3) each way, in this process: ``load`` is the call a template built before
writes (no ``outputs``: one ``load_mf6_run`` for the sets), ``outputs`` the one a
build writes now. Both must write identical tables, which it checks.

MEASURED 2026-10-12, a copy of Ten Trails ``ss`` (1 period, 5 x 20,323 DISV cells,
NEWTON; budget file ``tentrails.cbb`` from its OC), load average ~12.5-13.5 on 20
cores (other sessions' suites running), 3 repeats each:

  load     (the call a template built before writes):  min 10.55 s  (10.55 / 14.69 / 13.75 s)
  outputs  (the call a build writes now):              min  0.08 s  ( 0.08 /  0.08 /  0.08 s)
  every table identical: True

  The record itself: 0.05 s at build (the model already loaded for the base run).
  Ledger 540 measured the load at 8.3-9.0 s at load averages of 5-15.
"""

from __future__ import annotations

import ast
import json
import os
import re
import shutil
import sys
import time
import warnings
from pathlib import Path

import pandas as pd

import myflopy as mf
from myflopy.modflow.mf6.pest.forward_run import write_named_series_targets
from myflopy.modflow.mf6.pest.observations import snapshot_run_outputs
from myflopy.modflow.mf6.pest.run_outputs import OUTPUTS_SNAPSHOT

OUTPUT_SUFFIXES = {".cbb", ".cbc", ".hds", ".lst", ".csv", ".grb"}


def _cells_of(model, package_type: str, name: str | None = None) -> list[dict]:
    """Every (layer, cell) a package's period-0 records name, as zone rows by layer."""

    rows: dict[int, list[int]] = {}
    for package in model._packages_of_type(package_type):
        if name is not None and package.package_name.lower() != name:
            continue
        data = package.stress_period_data.get_data(0)
        for cellid in data["cellid"]:
            rows.setdefault(int(cellid[0]), []).append(int(cellid[-1]))
    return [{"layer": layer, "cells": sorted(set(cells))} for layer, cells in sorted(rows.items())]


def _zone(name: str, rows: list[dict]) -> list[dict]:
    return [{"name": name, **row} for row in rows]


def main():
    run = Path(sys.argv[1]).resolve()
    repeats = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    model = mf.load_mf6_run(run, verbosity_level=0)
    model.load_all()
    sets = [
        ("riv2", mf.RivFlowTargets(locations=_zone("riv", _cells_of(model, "riv", "riv"))
                                   + _zone("crisp", _cells_of(model, "riv", "crisp")))),
        ("crisp", mf.RivFlowTargets(locations=_zone("crisp", _cells_of(model, "riv", "crisp")),
                                    package="crisp")),
        ("drn", mf.DrnFlowTargets(locations=_zone("drains", _cells_of(model, "drn")))),
    ]
    cal = model.pest("outputs", start_datetime="2024-01-01")
    cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 1e4))
    for prefix, targets in sets:
        cal.observe(targets, prefix=prefix)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build("outputs.pst")
    template = Path(cal.template_workspace)
    for path in run.iterdir():
        if path.is_file() and path.suffix.lower() in OUTPUT_SUFFIXES:
            shutil.copy(path, template / path.name)
    call = re.search(r"write_named_series_targets\(sets=(\[.*\]), model_name='([^']+)'",
                     (template / "forward_run.py").read_text())
    entries, name = ast.literal_eval(call.group(1)), call.group(2)
    os.chdir(template)

    start = time.perf_counter()
    snapshot_run_outputs(model, ["riv", "drn"])
    print(f"snapshot (once, at build): {time.perf_counter() - start:.2f} s")
    tables = {}
    for mode in ("load", "outputs"):
        times = []
        for _ in range(repeats):
            written = [dict(entry, output_csv=f"{mode}_{entry['output_csv']}") for entry in entries]
            start = time.perf_counter()
            write_named_series_targets(written, sim_ws=template, model_name=name,
                                       outputs=OUTPUTS_SNAPSHOT if mode == "outputs" else None)
            times.append(time.perf_counter() - start)
        tables[mode] = [pd.read_csv(entry["output_csv"]) for entry in written]
        shown = " / ".join(f"{value:.2f}" for value in times)
        print(f"{mode:8s} {len(entries)} sets: min {min(times):.2f} s  ({shown} s)")
    same = all(a.equals(b) for a, b in zip(tables["load"], tables["outputs"], strict=True))
    print(f"every table identical: {same}")
    print(json.dumps({k: v for k, v in json.loads(Path(OUTPUTS_SNAPSHOT).read_text()).items()
                      if k != "packages"}))


if __name__ == "__main__":
    main()

"""The PEST forward run's named-series step: one model load per SET against one per RUN.

Ledger 540. ``write_named_series_targets`` took ONE lake/SFR/DRN set and ran
``mf.load_mf6_run`` + ``load_all()`` for it; the build wrote one call per set, so every
realization loaded the whole simulation once per set. It takes every set now and the
build writes one call.

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/named_series_one_load.py <template> before|after

``<template>`` is a built PEST template (``cal.build(...)``) holding the run's outputs --
run ``mf6`` in it, or copy the run's output files in. The sets are the ones the build
wrote into its ``forward_run.py``. ``before`` repeats what the old writer did per set
(a fresh ``load_mf6_run`` + ``load_all`` + the family's ``simulated_series``); ``after``
is one ``write_named_series_targets`` call. Wall clock, so read it at a known load.

TEMPLATES MEASURED (2026-10-12):

  * Ten Trails ``ss`` (5 x 20,323 DISV cells, NEWTON; a COPY of the read-only run),
    template externalized by the build: ``cal.observe`` of a two-zone RivFlowTargets
    (every river cell of ``riv`` and of ``crisp``, ``package=None``), a ``crisp``-only
    RivFlowTargets and a DRN polygon zone; outputs copied in.
  * The canonical model's testing profile (``build_canonical_model(...,
    CanonicalModelConfig.testing())``, run), observing its registered ``drn_flow``,
    ``lake_stage``, ``sfr_flow`` and ``sfr_stage`` sets; ``mf6`` run in the template.

MEASURED (load average 5-15 on 20 cores; MF6 alone takes 17-20 s on Ten Trails ``ss``):

  Ten Trails ``ss``, 3 sets:  before 30.2-30.4 s   after 8.3-9.0 s
  canonical testing, 4 sets:  before 1.04-1.10 s   after 0.37-0.44 s

  One set alone costs about what three do after (7.4 s on Ten Trails): the one full
  load. Loading only ``disv``, ``oc``, ``riv`` and ``drn`` (FloPy's ``load_only``) took
  4.5 s with identical output -- NOT done, see the ledger entry.
"""

from __future__ import annotations

import ast
import os
import re
import sys
import time
from pathlib import Path

import myflopy as mf
from myflopy.modflow.mf6.pest.forward_run import read_target_table, write_named_series_targets

BUILDERS = {
    "lake_stage": mf.LakeStageTargets,
    "sfr_stage": mf.SfrStageTargets,
    "sfr_flow": mf.SfrFlowTargets,
    "drn_flow": mf.DrnFlowTargets,
    "riv_flow": mf.RivFlowTargets,
}


def _old_way(template: Path, sets: list[dict], name: str) -> None:
    """What the writer did before, once PER SET."""

    for entry in sets:
        model = mf.load_mf6_run(template, model_name=name, verbosity_level=0)
        model.load_all()
        path = entry["locations_file"]
        if path.endswith(".gpkg"):
            import geopandas as gpd

            locations = gpd.read_file(path)
        else:
            locations = read_target_table(path)
        if "cells" in locations.columns:
            locations["cells"] = locations["cells"].fillna("").apply(
                lambda text: [int(v) for v in str(text).split(",") if v.strip()])
        BUILDERS[entry["kind"]](locations=locations, values=None,
                                package=entry.get("package")).simulated_series(model)


def main():
    template, mode = Path(sys.argv[1]).resolve(), sys.argv[2]
    call = re.search(r"write_named_series_targets\(sets=(\[.*\]), model_name='([^']+)'\)",
                     (template / "forward_run.py").read_text())
    sets, name = ast.literal_eval(call.group(1)), call.group(2)
    os.chdir(template)
    start = time.perf_counter()
    if mode == "before":
        _old_way(template, sets, name)
    else:
        write_named_series_targets(sets, sim_ws=template, model_name=name)
    print(f"{mode}: {len(sets)} sets, {time.perf_counter() - start:.2f} s")


if __name__ == "__main__":
    main()

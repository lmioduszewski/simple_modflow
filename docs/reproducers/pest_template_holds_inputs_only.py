"""What a PEST template holds of the run it was built from (ledger 383).

The build's staging copy (ledgers 356-357) copied the WHOLE run directory, and
``PstFrom`` copies the staging directory into the template -- and pyEMU copies the
template into every IES worker. It now copies what ``mfsim.nam`` leads MF6 to read
(``project._simulation_input_files``). This script builds one calibration on the
canonical model's run and sorts what the template holds, on whichever tree is on
``PYTHONPATH``:

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python \\
        docs/reproducers/pest_template_holds_inputs_only.py [--canonical full|validation|testing]

The "before" column is the parent commit's ``src`` (``git archive <parent> src | tar -x
-C <dir>``, then ``PYTHONPATH=<dir>/src``); the classes are read with THIS checkout's
``_simulation_input_files``, taken from its source file, so both columns sort by one
rule. It writes only under a temporary directory, which it removes.

THE SETUP. ``build_canonical_calibration_demo`` (4 head wells) writes and runs the
canonical model in its workspace, as every canonical notebook does; a 1 MB file in
``<run>/pest/`` stands in for a PEST workspace built inside the run before ledger 107.
The calibration: ``k`` on layer 0, ``style="grid"``, ``correlation=3000``; the head
targets. Each template file is put in ONE class: ``build`` (a file the build wrote --
pst, csv, ins, tpl, ``mult/``, ``org/``, ``forward_run.py``, logs, markers), ``input``
(a file ``mfsim.nam`` leads MF6 to read in the template), ``old pest`` (under
``pest/``), or ``run's other`` (anything else: outputs, run metadata).

MEASURED (FloPy 3.10.0, pyEMU 1.4.0, 2026-10-11; "before" is f71e442)   before          after
  testing (21x21x4), run directory 50 files / 3.89 MB
    template files / MB                                          180 / 4.18      154 / 0.73
      input                                                      134 / 0.52      134 / 0.52
      build                                                       20 / 0.22       20 / 0.22
      old pest                                                     2 / 1.00        0
      run's other (heads, budgets, listings, OBS CSVs, a pickle)  24 / 2.45        0
  full (100x100x4), run directory 50 files / 102.84 MB
    template files / MB                                          180 / 106.82    154 / 13.77
      input / build                                     134 / 11.33, 20 / 2.45   the same
      old pest / run's other                              2 / 1.00, 24 / 92.04   0 / 0
Every input and build file is byte-identical between the two columns but for FloPy's
date stamp and the paths in the metadata (compared on the testing profile). pyEMU
copies the template into every IES worker, so each worker carried the 92 MB too.
"""

from __future__ import annotations

import argparse
import ast
import os
import re
import shutil
import tempfile
import warnings
from pathlib import Path

#: This checkout's PEST project module, whose input-file rule sorts both columns.
_PROJECT = Path(__file__).resolve().parents[2] / "src/myflopy/modflow/mf6/pest/project.py"
_RULE = {"_NAME_FILE_BLOCKS", "_MF6_FILE_TYPE", "_MF6_FIELD", "_MF6_INPUT_REFERENCE",
         "_mf6_field", "_mf6_text", "_simulation_input_files"}


def _input_rule():
    """``_simulation_input_files`` as this checkout defines it, whatever is imported."""

    tree = ast.parse(_PROJECT.read_text(encoding="utf-8"))
    nodes = [node for node in tree.body
             if (isinstance(node, ast.FunctionDef) and node.name in _RULE)
             or (isinstance(node, ast.Assign) and {getattr(t, "id", None) for t in node.targets} & _RULE)]
    namespace = {"re": re, "os": os, "Path": Path}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(_PROJECT), "exec"), namespace)
    return namespace["_simulation_input_files"]


def _classify(template: Path, inputs: set[str]) -> dict[str, list[int]]:
    classes = {"input": [0, 0], "build": [0, 0], "old pest": [0, 0], "run's other": [0, 0]}
    build_suffixes = (".pst", ".csv", ".ins", ".tpl", ".py", ".sh", ".log", ".json", ".gpkg")
    for path in sorted(template.rglob("*")):
        if not path.is_file():
            continue
        name = str(path.relative_to(template))
        if name in inputs:
            kind = "input"
        elif name.startswith("pest/"):
            kind = "old pest"
        elif name.split("/")[0] in ("mult", "org") or (
            "/" not in name and name.endswith(build_suffixes)
            and not name.startswith(("ims_", "master_", "viz_prt_master"))
        ):
            kind = "build"
        else:
            kind = "run's other"
        classes[kind][0] += 1
        classes[kind][1] += path.stat().st_size
    return classes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--canonical", default="testing", choices=("full", "validation", "testing"))
    args = parser.parse_args()

    from myflopy.modflow.mf6.canonical_calibration import build_canonical_calibration_demo
    from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig

    _simulation_input_files = _input_rule()
    warnings.simplefilter("ignore")
    config = {"full": CanonicalModelConfig, "validation": CanonicalModelConfig.validation,
              "testing": CanonicalModelConfig.testing}[args.canonical]()
    root = Path(tempfile.mkdtemp(prefix="pest_inputs_only_"))
    try:
        demo = build_canonical_calibration_demo(root / "model", config=config, n_head_wells=4)
        run = Path(demo.model.workspace)
        (run / "pest").mkdir()
        (run / "pest" / "old.pst").write_text("pcf\n")
        (run / "pest" / "old.bin").write_bytes(b"\0" * 1_000_000)
        run_files = [path for path in run.rglob("*") if path.is_file()]
        cal = demo.model.pest("inputs", start_datetime="2024-01-01")
        cal.parameterize("k", style="grid", correlation=3000, layers=[0], bounds=(0.5, 2.0))
        cal.observe(demo.head_targets)
        cal.build("inputs.pst", noptmax=0)
        template = cal.template_workspace
        inputs = {str(path) for path in _simulation_input_files(template)}
        classes = _classify(template, inputs)
        total = [sum(c[0] for c in classes.values()), sum(c[1] for c in classes.values())]
        print(f"run directory: {len(run_files)} files, "
              f"{sum(p.stat().st_size for p in run_files) / 1e6:.2f} MB")
        print(f"template: {total[0]} files / {total[1] / 1e6:.2f} MB")
        for kind, (count, size) in classes.items():
            print(f"  {kind:12s} {count:5d} / {size / 1e6:.2f} MB")
    finally:
        shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    main()

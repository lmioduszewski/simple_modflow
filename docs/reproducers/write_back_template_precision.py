"""A write-back fills templates at full precision; pyEMU's filler does not (ledger 560).

``IesResults.calibrated()`` rebuilds what a realization's MODFLOW 6 read by filling the
calibration's templates with that realization's values and replaying its forward run.
pyEMU's own filler, ``pyemu.pst_utils.write_to_template`` (what ``Pst.write_input_files``
calls), fits each value to the WIDTH of its marker span. The two template shapes a
myflopy build writes are both narrow:

  * a pilot-point set's ``<set>_pp.csv.tpl`` -- ``kl0pp0000,~  kl0pp0000  ~`` (15 wide);
  * a pyEMU multiplier file's ``<name>_inst0_constant.csv.tpl`` -- one span per cell.

This fills both shapes with the same values through pyEMU's filler and through the
write-back's (``calibrated._fill_templates``, ``repr`` of the value), and reads them
back. No model and no PEST++ run:

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/write_back_template_precision.py

MEASURED (pyEMU 1.4.0, 2026-10-12): the worst relative error over 200 lognormal values
(sigma 1 decade) and the base realization's 4.681308128714071

                         pyEMU write_to_template      write-back (repr)
  pilot-point span       4.53e-04                     0.00e+00
  multiplier span        4.11e-07                     0.00e+00

The pilot-point span is the one that matters: 15 characters with two blanks each side
hold 4 significant digits. (Read back with ``float_precision="round_trip"``: pandas'
default C reader is itself off by up to ~5e-13 relative, which is the reader, not the
file.)

(The design measured the same thing end to end on a pilot-point K field: 3.4e-4 against
the realization's captured field with pyEMU's filler, 5e-6 -- the captured field's own 6
digits -- with the replay.)
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pyemu

from myflopy.modflow.mf6.pest.calibrated import _fill_templates

rng = np.random.default_rng(2026)
values = np.concatenate([[4.681308128714071], 10 ** rng.normal(0.0, 1.0, 199)])
names = [f"kl0pp{index:04d}" for index in range(values.size)]

SHAPES = {
    "pilot-point span": "ptf ~\nparnme,value\n" + "".join(
        f"{name},~  {name}  ~\n" for name in names),
    "multiplier span": "ptf ~\n" + "".join(f" ~   {name}    ~\n" for name in names),
}


def _read(path: Path, shape: str) -> np.ndarray:
    if shape == "pilot-point span":
        return pd.read_csv(path, float_precision="round_trip")["value"].to_numpy(dtype=float)
    return np.loadtxt(path, dtype=float)


with tempfile.TemporaryDirectory() as scratch:
    root = Path(scratch)
    print(f"{'':24s}{'pyEMU write_to_template':>26s}{'write-back (repr)':>22s}")
    for shape, text in SHAPES.items():
        (root / "t.tpl").write_text(text)
        series = pd.Series(values, index=names)
        pyemu.pst_utils.write_to_template(series, str(root / "t.tpl"), str(root / "pyemu.in"))
        pst = SimpleNamespace(template_files=["t.tpl"], input_files=["ours.in"])
        _fill_templates(pst, dict(zip(names, values, strict=True)), root, root)
        worst = [
            float(np.max(np.abs(_read(root / f, shape) / values - 1.0)))
            for f in ("pyemu.in", "ours.in")
        ]
        print(f"{shape:24s}{worst[0]:>26.2e}{worst[1]:>22.2e}")

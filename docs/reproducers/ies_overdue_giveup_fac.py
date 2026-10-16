"""What ``overdue_giveup_fac`` a parallel PESTPP-IES run actually applies (ledger 510).

PESTPP-IES kills a parallel model run that takes longer than ``overdue_giveup_fac``
times the average run and counts its realization failed. Left unset, it resets the
factor to 2.0 when the run starts. ``run_ies``/``prior`` now set 5.0 unless told
otherwise (``overdue_giveup_fac=None`` removes the option). This script launches a
parallel prior Monte Carlo three ways and reads, from each master's ``.rec`` record,
the factor PESTPP-IES states:

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/ies_overdue_giveup_fac.py

Needs ``mf6``, ``pestpp-ies`` and pyEMU; writes only under a temporary directory,
which it removes.

THE SETUP. The 6 x 6 irregular Voronoi toy of ``tests/test_pest_handed_targets.py``
(one layer, CHD east and west, a stream, a drain zone, a lake, three 10-day periods),
run once; the calibration is ``k`` constant and three head targets;
``prior(reals=3, workers=2)`` with ``overdue_giveup_fac`` left at its default, set to
12, and ``None``.

MEASURED (PEST++ 5.2.16, pyEMU 1.4.0, 2026-10-12)
  overdue_giveup_fac=     the .rec's lines naming it
  (default)               overdue_giveup_fac: 5
  12                      overdue_giveup_fac: 12
  None                    overdue_giveup_fac: 100
                          ...resetting overdue_giveup_fac to 2.0
The run is not slow enough to kill anything; what this pins is the factor in force. The
kill itself was seen once (scale-7): on a 3-realization Ten Trails steady run under a
load average of 7-16, realization 0 ran 1.83 min against a 0.91 min average and was
killed at the 2.0 PESTPP-IES had reset, and dropped.
"""

from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _absolute_pythonpath import make_pythonpath_absolute  # noqa: E402

make_pythonpath_absolute()

import tempfile  # noqa: E402
import warnings  # noqa: E402

import pandas as pd  # noqa: E402

import myflopy as mf  # noqa: E402

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from tests.test_pest_handed_targets import run_toy  # noqa: E402

DEFAULT = object()


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="overdue_") as directory:
        root = pathlib.Path(directory)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = run_toy(root)
            cal = model.pest("overdue", workspace=root / "overdue", start_datetime="2024-01-01")
            cal.parameterize("k", style="constant", bounds=(0.5, 2.0), physical=(1e-3, 100.0))
            cal.observe(mf.HeadTargets(
                locations=pd.DataFrame({"name": ["W2", "W1", "W10"], "cell": [1, 2, 34],
                                        "layer": 0}),
                values=pd.DataFrame({"time": [2], "W2": [92.0], "W1": [93.0], "W10": [96.0]}),
            ))
            cal.build("overdue.pst")
        for value in (DEFAULT, 12, None):
            options = {} if value is DEFAULT else {"overdue_giveup_fac": value}
            # Each run replaces the last one's master (ledger 532).
            results = cal.prior(reals=3, workers=2, overwrite=True, **options)
            record = (pathlib.Path(results.workspace) / "overdue.rec").read_text(errors="replace")
            lines = [line.strip() for line in record.splitlines() if "overdue_giveup_fac" in line]
            shown = "(default)" if value is DEFAULT else repr(value)
            print(f"overdue_giveup_fac={shown}:")
            for line in lines:
                print(f"    {line}")


if __name__ == "__main__":
    main()

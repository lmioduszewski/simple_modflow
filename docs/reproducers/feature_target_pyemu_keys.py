"""What pyEMU does with LAK connections and SFR reaches, by how a row is keyed.

Ledger 558 (``lak.bedleak`` and ``sfr.rhk``). Each section below calls
pyEMU's ``PstFrom.add_parameters`` directly on the external files of a small MF6 model
-- ``tests/pest_feature_models.py``'s: lake 0 ``Keevie`` with three connections, two of
them to one cell, lake 1 ``marjorie``; reaches 0-2 ``creek_a`` and 3-4 ``Creek_B``,
reaches 1 and 2 in one cell -- and applies the multipliers as a forward run does.

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/feature_target_pyemu_keys.py

About 30 s; it writes only under a temporary directory it removes.

MEASURED (pyEMU 1.4.0, FloPy 3.10.0, MF6 6.8.0.dev0, 2026-10-12)
  A. LAK keyed by ``lakeno`` alone, ``use_rows`` one lake, x10: APPLIED CORRECTLY --
     lake 1 -> 1.0 2.0 2.5 0.3 0.35, lake 2 -> 0.1 0.2 0.25 3.0 3.5. The design pass
     reported an ``AssertionError`` (``Probable miss-alignment in tpl indices``) here;
     it does not reproduce on pyEMU 1.4.0, so the library's ``(lakeno, iconn)`` key is
     kept because it is UNIQUE per row, not because the lake-only key was seen to fail.
  B. Keyed by ``(lakeno, iconn)``, ``use_rows`` lake 1's three connections, x10:
     bedleak 0.1 0.2 0.25 0.3 0.35 -> 1.0 2.0 2.5 0.3 0.35 (Keevie alone; claktype kept).
  C. SFR, reach 4's cellid ``NONE``, every reach x2 on column 8 (rhk on DISV): reaches
     0-3 rhk doubled, and reach 4 -- one column short -- had its Manning's n doubled,
     0.034 -> 0.068, its rhk untouched. pyEMU says nothing.
  D. ``use_rows`` naming no row of the file: pyEMU warns ``None of request rows found``
     and sets up a parameter on EVERY row (5 of 5 scaled).
  E. ``style="grid"`` with a variogram (range 50) over all five reaches, keyed
     ``(ifno, layer, cell)``: reaches 1 and 2 sit at one cell centre, and pyEMU gives
     EVERY reach a covariance group of its own (``sfrrhk_cov0`` .. ``_cov4``); 2,000
     prior draws correlate the reaches within 0.03 of 0. Without reach 2 (no shared
     cell) the same draw correlates neighbours 10 ft apart at 0.82.
"""

from __future__ import annotations

import contextlib
import io
import os
import shutil
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))


def _external(run: Path, workspace: Path):
    import flopy

    sim = flopy.mf6.MFSimulation.load(sim_ws=str(run), verbosity_level=0)
    sim.set_sim_path(str(workspace))
    sim.set_all_data_external()
    sim.write_simulation(silent=True)


def _apply(template: Path, pst, values: dict[str, float]) -> None:
    import pyemu

    par = pst.parameter_data
    for fragment, value in values.items():
        par.loc[[name for name in par.index if fragment in name], "parval1"] = value
    pst.write_input_files(pst_path=str(template))
    here = Path.cwd()
    os.chdir(template)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            pyemu.helpers.apply_list_and_array_pars(arr_par_file="mult2model_info.csv")
    finally:
        os.chdir(here)


def _pf(workspace: Path, template: Path):
    import pyemu

    return pyemu.utils.PstFrom(str(workspace), str(template), remove_existing=True,
                               longnames=True, zero_based=False)


def _column(path: Path, column: int) -> list[str]:
    return [line.split()[column] for line in path.read_text().splitlines() if line.strip()]


def main() -> None:
    from tests.pest_feature_models import build_model

    root = Path(tempfile.mkdtemp(prefix="feature_keys_"))
    try:
        lak = "m.lak_connectiondata.txt"
        sfr = "m.sfr_packagedata.txt"
        run = build_model(root / "run")
        _external(run, root / "ext")
        none = build_model(root / "none", none_reach=True)
        _external(none, root / "ext_none")
        for lake in (1, 2):
            with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
                warnings.simplefilter("ignore")
                pf = _pf(root / "ext", root / "a")
                pf.add_parameters(lak, par_type="constant", index_cols=[0], use_cols=[5],
                                  use_rows=[(lake,)], par_name_base="k", pargp="k")
                pst = pf.build_pst()
            try:
                _apply(root / "a", pst, {"k": 10.0})
                print(f"A. keyed by lakeno alone, lake {lake}: applied --",
                      _column(root / "a" / lak, 5))
            except AssertionError as error:
                print(f"A. keyed by lakeno alone, lake {lake}: AssertionError -- "
                      f"{str(error)[:60]}")

        with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
            warnings.simplefilter("ignore")
            pf = _pf(root / "ext", root / "b")
            pf.add_parameters(lak, par_type="constant", index_cols=[0, 1], use_cols=[5],
                              use_rows=[(1, 1), (1, 2), (1, 3)], par_name_base="k", pargp="k")
            pst = pf.build_pst()
        _apply(root / "b", pst, {"k": 10.0})
        print("B. keyed by (lakeno, iconn):", _column(root / "b" / lak, 5),
              _column(root / "b" / lak, 4))

        with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
            warnings.simplefilter("ignore")
            pf = _pf(root / "ext_none", root / "c")
            pf.add_parameters(sfr, par_type="constant", index_cols=[0], use_cols=[8],
                              par_name_base="r", pargp="r")
            pst = pf.build_pst()
        _apply(root / "c", pst, {"r": 2.0})
        rows = [line.split() for line in (root / "c" / sfr).read_text().splitlines()]
        print("C. NONE cellid, x2 on column 8: rows", [row[7:10] for row in rows])

        with warnings.catch_warnings(record=True) as caught, \
                contextlib.redirect_stdout(io.StringIO()):
            warnings.simplefilter("always")
            pf = _pf(root / "ext", root / "d")
            frame = pf.add_parameters(lak, par_type="constant", index_cols=[0, 1],
                                      use_cols=[5], use_rows=[(9, 9)],
                                      par_name_base="k", pargp="k")
        tpl = pd.read_csv(root / "d" / Path(frame["tpl_filename"].iloc[0]).name, skiprows=1)
        said = [str(item.message)[:50] for item in caught if "request rows" in str(item.message)]
        print(f"D. use_rows naming no row: {len(tpl)} of 5 rows parameterized; pyEMU: {said}")

        from myflopy.modflow.mf6.pest import native_parameters as np_mod
        from tests.pest_feature_models import build, calibration

        for label, features in (("all five reaches", None), ("reaches 0, 1, 3, 4", [0, 1, 3, 4])):
            cal = calibration(run, root / f"e_{len(features or [])}", name="e")
            spec = cal.parameterize("rhk", style="grid", correlation=50.0, bounds=(0.5, 2.0),
                                    features=features)
            # Past the build's refusal of the shared cell, to see what pyEMU does.
            original = np_mod._refuse_an_uncorrelatable_reach_grid
            np_mod._refuse_an_uncorrelatable_reach_grid = lambda *args: None
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    build(cal, "e")
            finally:
                np_mod._refuse_an_uncorrelatable_reach_grid = original
            frame = cal._native_parameter_frames[spec.name]
            with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
                warnings.simplefilter("ignore")
                drawn = cal.draw_prior(2000, use_specsim=False)
            values = drawn._df if hasattr(drawn, "_df") else drawn
            logs = np.log10(values[[name for name in values.columns if "sfrrhk" in name]])
            corr = logs.corr().to_numpy()
            off = corr[~np.eye(len(corr), dtype=bool)]
            print(f"E. {label}: covgp {sorted(set(frame['covgp']))}; off-diagonal "
                  f"correlation {off.min():.2f} .. {off.max():.2f}")
    finally:
        shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    main()

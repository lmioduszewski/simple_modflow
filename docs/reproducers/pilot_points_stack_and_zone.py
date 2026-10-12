"""Pilot points stack with other multipliers and honour zones: measured with real forward runs.

Ledger 525 and 527. A pilot-point set used to REBUILD its layer
file from the base captured at build, in a forward-run command that runs after pyEMU's
``apply_list_and_array_pars``: whatever pyEMU had just written to that file -- a
constant, zone or grid multiplier -- and whatever an earlier pilot-point set had
written, was thrown away, while the forward run exited 0. ``zones=`` was accepted and
ignored. Each case below builds a template on the canonical testing model, writes one
realization's parameter values into it, runs its ``forward_run.py`` (pyEMU's apply,
the pilot-point step, MF6) and reads the K file MF6 read.

    PYTHONPATH=src /home/lukem/python/envs/gw/bin/python docs/reproducers/pilot_points_stack_and_zone.py

Uses only the public API, so it measures whichever tree is on ``PYTHONPATH`` (made
absolute at start, ledger 425): the "before" column is the parent commit's ``src``,
where the net spacing was ``pp_space=`` in median cell widths (6 and 3 on these
~100 m cells) rather than ``pp_spacing=`` (600.0 and 300.0). Needs ``mf6``. About a
minute.

THE SETUP. ``CanonicalModelConfig.testing()`` (21 x 21 Voronoi cells of ~100 m, 4
layers), four head targets. Layer 0's K file, written / base:

  A. ``k`` constant on layer 0 at x2, plus pilot points on layer 0, every point 1.0.
  B. two pilot-point calls on layer 0 (``ka`` every point 2.0, ``kb`` every point 1.5).
  C. ``zones=`` west half 1, east half 0, every pilot point 2.0.
  D. two calls under the default name, ``layers=[0]`` and ``layers=[1]``.
  E. ``k`` constant on layer 0 at 1.0 with ``physical=`` capped at 1.2x the layer's
     largest K, plus pilot points on layer 0 with no ``physical=``, every point 5.0
     (ledger 528; "before" there is the R1-pilot commit 09f3a1b9, the parent
     having no stacking at all).

MEASURED (mf6 6.8.0.dev0, pyEMU 1.4.0, 2026-10-12; "before" is claude/pestprep 3a2d6f9)
                                                    before          after
  A  written / base on layer 0                      1.000000        2.000000
  B  written / base on layer 0                      1.500000        3.000000
     (before: the LAST set alone -- each set rebuilt the file from the base)
  C  east (zone 0) cells moved, of 221              221             0
     west (zone 1) written / base                   2.000000        2.000000
  D  pilot-point groups in the control file         kl1             kl0, kl1
  E  largest written / the constant's cap           4.166667        1.000000
"""

from __future__ import annotations

import inspect
import pathlib
import shutil
import subprocess
import sys
import tempfile
import warnings

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _absolute_pythonpath import make_pythonpath_absolute  # noqa: E402

# The forward run starts in the template, where a relative PYTHONPATH would import
# another myflopy: made absolute here, before anything changes directory (ledger 425).
make_pythonpath_absolute()

from myflopy.modflow.mf6.canonical_calibration import (  # noqa: E402
    build_canonical_calibration_demo,
)
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig  # noqa: E402
from myflopy.modflow.mf6.pest.project import PestProject  # noqa: E402

_NEW_API = "pp_spacing" in inspect.signature(PestProject.parameterize).parameters


def _spacing(cells: float) -> dict:
    """The net spacing, ``cells`` cell widths of ~100 m, in whichever spelling the tree has."""

    return {"pp_spacing": 100.0 * cells} if _NEW_API else {"pp_space": int(cells)}


def _write_inputs(cal):
    """``Pst.write_input_files`` serially (pyEMU's opens a process pool)."""

    from pyemu.pst.pst_utils import write_to_template

    par = cal.pst.parameter_data
    values = par["parval1"] * par["scale"] + par["offset"]
    values.index = values.index.str.lower()
    for tpl, model_file in zip(cal.pst.template_files, cal.pst.input_files, strict=True):
        write_to_template(values, str(cal.template_workspace / tpl),
                          str(cal.template_workspace / model_file))


def _case(demo, root, label, declarations, values):
    cal = demo.model.pest(label, workspace=root / label, start_datetime="2024-01-01")
    for kwargs in declarations:
        cal.parameterize("k", **kwargs)
    cal.observe(demo.head_targets)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cal.build(f"{label}.pst", noptmax=0)
    par = cal.pst.parameter_data
    for group, value in values.items():
        par.loc[par["pargp"] == group, "parval1"] = value
    _write_inputs(cal)
    result = subprocess.run([sys.executable, "forward_run.py"], cwd=cal.template_workspace,
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stdout[-1500:] + result.stderr[-1500:]
    written = np.loadtxt(cal.template_workspace / f"{demo.model.name}.npf_k_layer1.txt")
    return cal, written / np.asarray(demo.model.gwf.npf.k.get_data(), dtype=float)[0]


def main():
    root = pathlib.Path(tempfile.mkdtemp(prefix="pp_stack_"))
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            demo = build_canonical_calibration_demo(
                root / "model", config=CanonicalModelConfig.testing(), n_head_wells=4)
        bounds = {"bounds": (0.1, 10.0)}

        _, ratio = _case(demo, root, "a", [
            dict(style="constant", layers=[0], name="kc", **bounds),
            dict(style="pilotpoints", layers=[0], name="kp", **_spacing(6), **bounds),
        ], {"kc": 2.0})
        print(f"A  constant x2 under pilot points at 1.0: written/base {np.median(ratio):.6f}")

        _, ratio = _case(demo, root, "b", [
            dict(style="pilotpoints", layers=[0], name="ka", **_spacing(6), **bounds),
            dict(style="pilotpoints", layers=[0], name="kb", **_spacing(3), **bounds),
        ], {"kal0": 2.0, "kbl0": 1.5})
        print(f"B  two sets at 2.0 and 1.5: written/base {np.median(ratio):.6f}")

        xc = np.asarray(demo.model.gwf.modelgrid.xcellcenters, dtype=float).reshape(-1)
        zones = np.where(xc < np.median(xc), 1, 0)
        cal, ratio = _case(demo, root, "c", [
            dict(style="pilotpoints", layers=[0], name="kp", zones=zones, **_spacing(6),
                 **bounds),
        ], {"kpl0z1": 2.0, "kpl0": 2.0})
        moved = int((~np.isclose(ratio[zones == 0], 1.0, rtol=1e-9)).sum())
        print(f"C  zones= west 1 / east 0, points at 2.0: {moved} of {int((zones == 0).sum())} "
              f"east cells moved; west written/base {np.median(ratio[zones == 1]):.6f}")

        cap = 1.2 * float(np.asarray(demo.model.gwf.npf.k.get_data(), dtype=float)[0].max())
        _, ratio = _case(demo, root, "e", [
            dict(style="constant", layers=[0], name="kc", bounds=(0.1, 10.0),
                 physical=(1e-6, cap)),
            dict(style="pilotpoints", layers=[0], name="kp", **_spacing(6), **bounds),
        ], {"kpl0": 5.0})
        base0 = np.asarray(demo.model.gwf.npf.k.get_data(), dtype=float)[0]
        print(f"E  constant capped by physical=, pilot points at 5.0: largest written / cap "
              f"{(ratio * base0).max() / cap:.6f}")

        cal = demo.model.pest("d", workspace=root / "d", start_datetime="2024-01-01")
        cal.parameterize("k", style="pilotpoints", layers=[0], **_spacing(6), **bounds)
        cal.parameterize("k", style="pilotpoints", layers=[1], **_spacing(3), **bounds)
        cal.observe(demo.head_targets)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cal.build("d.pst", noptmax=0)
        groups = sorted(set(cal.pst.parameter_data["pargp"]))
        print(f"D  two calls under one name on layers 0 and 1: groups {groups}")
    finally:
        shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    main()

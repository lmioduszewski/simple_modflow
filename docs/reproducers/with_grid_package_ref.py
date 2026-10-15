"""A per-cell array behind `mf.ref`, or on a deferred grid, moved as `unchanged`.

`with_grid`/`move_report` classified a `PackageRef` rather than the library package it
names. A ref carries no options, inputs or metadata, so `package_move` answered
`unchanged` -- "names no cells" -- and `with_grid` never stamped it. Separately, a model
on a DEFERRED grid (`mf.grid_ref`) had no mesh in hand at the call, so even an INLINE
array had nothing to be measured against (the limit ledgers 271/295 recorded).

Setup: two DIFFERENT meshes with the SAME cell count (4x7 and 7x4
`rectangular_voronoi`, ncpl 28 each -- the shape no range check sees), a lazy
`mf.disv`, `strt` a `(2, 28)` array, `npf(k=1.0)` as the control that genuinely names no
cells. Each row: `move_report` before the move, after it, and what `prepare_run` does.

Measured on `cb8073e5` (before; no `library=` exists there)::

    context grid  inline  unverifiable  unverifiable  REFUSED
    context grid  ref     unchanged     unchanged     BUILT on the old strt
    grid_ref      inline  unchanged     unchanged     BUILT on the old strt
    grid_ref      ref     unchanged     unchanged     BUILT on the old strt
    remap=True    ref     unchanged  (inline twin: remapped)

After (ledger 597, `library=project`)::

    context grid  inline  unverifiable  unverifiable  REFUSED
    context grid  ref     unverifiable  unverifiable  REFUSED
    grid_ref      inline  unverifiable  unverifiable  REFUSED
    grid_ref      ref     unverifiable  unverifiable  REFUSED
    remap=True    ref     remapped   (an inline copy; the library package unchanged)

`npf` reads `unchanged` in every row, before and after.

Run: `PYTHONPATH=src python docs/reproducers/with_grid_package_ref.py`
"""

from __future__ import annotations

import inspect
import tempfile
import warnings
from pathlib import Path

import numpy as np

import myflopy as mf
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, rectangular_voronoi
from myflopy.specs import ModelContext
from myflopy.workspace import Project

warnings.simplefilter("ignore")
HAS_LIBRARY = "library" in inspect.signature(mf.SimulationSpec.with_grid).parameters
a = rectangular_voronoi(CanonicalModelConfig(nrow=4, ncol=7, nlay=2, nper=1))
b = rectangular_voronoi(CanonicalModelConfig(nrow=7, ncol=4, nlay=2, nper=1))
assert a.ncpl == b.ncpl == 28
strt = np.linspace(50.0, 90.0, 2 * a.ncpl).reshape(2, a.ncpl)


def _stack():
    stack = mf.LayerStack(top=mf.Surface.flat(100.0))
    stack.add("upper", thickness=20.0)
    stack.add("lower", thickness=30.0)
    return stack


def _sim(ic, npf, **gwf):
    model = mf.gwf(
        "m", packages=[mf.disv(_stack()), npf, ic, mf.oc(head_filerecord="m.hds")], **gwf
    )
    return mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])


def _lib(project):
    return {"library": project} if HAS_LIBRARY else {}


def _outcome(project, moved):
    try:
        run = project.prepare_run("r", moved, overwrite=True)
    except ValueError:
        return "REFUSED"
    built = run.built.built_model("m").model.get_package("ic").strt.array
    return "BUILT on the old strt" if np.array_equal(built, strt) else "BUILT"


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="with_grid_package_ref-") as scratch:
        _measure(Path(scratch))


def _measure(root: Path) -> None:
    on_a = {"context": ModelContext(grid=a, domain=np.ones((2, a.ncpl), dtype=int))}
    for label, gwf in (("context grid", on_a), ("grid_ref", {"grid": mf.grid_ref("mesh")})):
        for form in ("inline", "ref"):
            project = Project(root / f"{label[:4]}_{form}", name="p")
            project.add_grid("mesh", mf.GridSpec.from_object(a))
            project.add_package("ic_regrid", mf.ic(strt=strt))
            project.add_package("npf_scalar", mf.npf(k=1.0))
            if form == "ref":
                sim = _sim(mf.ref("ic_regrid"), mf.ref("npf_scalar"), **gwf)
                ic, npf = "ic_regrid", "npf_scalar"
            else:
                sim = _sim(mf.ic(strt=strt), mf.npf(k=1.0), **gwf)
                ic, npf = "ic", "npf"
            before = {r["package"]: r["how"] for r in sim.move_report(**_lib(project))}
            moved = sim.with_grid(b, **_lib(project))
            after = {r["package"]: r["how"] for r in moved.move_report(**_lib(project))}
            print(
                f"{label:13s} {form:6s}  {before[ic]:12s}  {after[ic]:12s}  "
                f"{_outcome(project, moved):22s}  (npf: {after[npf]})"
            )

    project = Project(root / "remap", name="p")
    project.add_package("ic_regrid", mf.ic(strt=strt))
    by_ref = _sim(mf.ref("ic_regrid"), mf.npf(k=1.0), **on_a).with_grid(
        b, remap=True, **_lib(project)
    )
    inline = _sim(mf.ic(strt=strt), mf.npf(k=1.0), **on_a).with_grid(b, remap=True)
    ref_row = next(r for r in by_ref.move_report(**_lib(project)) if r["package"].startswith("ic"))
    inline_row = next(r for r in inline.move_report() if r["package"] == "ic")
    print(
        f"remap=True    ref     {ref_row['how']}  (inline twin: {inline_row['how']}; "
        f"{ref_row['why'].split('; ')[-1]})"
    )


if __name__ == "__main__":
    main()

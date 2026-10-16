"""A calibrated K behind `mf.ref`, beside an eager DISV, could not be remapped onto a new mesh.

The route a calibration's write-back takes (`calibrated().to_library` registers an `npf`
built by `substitute_values`, and the model names it by `mf.ref`) beside an EAGER
`mf.disv(<LayerBuildResult>)` -- the shape of every Ten Trails run. Ledger 597 made a
move read the ref; five more rows went wrong on it (ledger 610).

Setup: `rectangular_voronoi` meshes, 4x7 and 7x4 (ncpl 28 each) and 9x9 (ncpl 81), two
layers, `k` a `(2, 28)` array; no MF6 run. Each row prints what the call did.

Measured on 55947911 (before)::

    eager DISV inline, remap=True                  RAISED RemapRefused
    eager DISV by ref + npf by ref, remap=True     RAISED RemapRefused
    grid_free=('npf',) under remap=True            k carried: True
    stamped npf moved 28 -> 81, built              RAISED MFDataException
    reloaded, simulation-tier with_grid            moved models: []
    reloaded, model-tier with_grid, built          BUILT, k == the OLD mesh's K: True
    remap onto project.grid('new') off the root    RAISED FileNotFoundError
    remap onto mf.grid_ref('new'), library=        RAISED TypeError

After (ledger 610)::

    eager DISV inline, remap=True                  disv stale, npf remapped
    eager DISV by ref + npf by ref, remap=True     dis/old stale, npf remapped
    grid_free=('npf',) under remap=True            k carried: False
    stamped npf moved 28 -> 81, built              RAISED ValueError (names the move and remap=True)
    reloaded, simulation-tier with_grid            RAISED ValueError (names the recorded mesh)
    reloaded, model-tier with_grid, built          RAISED ValueError (names the recorded mesh)
    remap onto project.grid('new') off the root    dis/old stale, npf remapped
    remap onto mf.grid_ref('new'), library=        dis/old stale, npf remapped

Run: `PYTHONPATH=src python docs/reproducers/remap_a_calibrated_ref.py`
"""

from __future__ import annotations

import os
import tempfile
import warnings
from pathlib import Path

import numpy as np

import myflopy as mf
from myflopy.modflow.mf6.canonical_example import CanonicalModelConfig, rectangular_voronoi
from myflopy.specs import ModelContext, substitute_values
from myflopy.workspace import Project

warnings.simplefilter("ignore")


def _mesh(nrow, ncol):
    return rectangular_voronoi(CanonicalModelConfig(nrow=nrow, ncol=ncol, nlay=2, nper=1))


a, b, c = _mesh(4, 7), _mesh(7, 4), _mesh(9, 9)
k_old = np.linspace(1.0, 50.0, 2 * a.ncpl).reshape(2, a.ncpl)


def _layering(grid):
    stack = mf.LayerStack(top=mf.Surface.flat(100.0))
    stack.add("upper", thickness=20.0)
    stack.add("lower", thickness=30.0)
    return stack.build(grid)


def _calibrated():
    return substitute_values(
        mf.npf(k=np.ones_like(k_old), icelltype=1), {"k": k_old},
        provenance={"pest": "cal", "realization": "base", "source_key": "npf/base"},
    )


def _model(grid, *packages):
    context = ModelContext(grid=grid, domain=np.ones((2, grid.ncpl), dtype=int))
    return mf.gwf("m", context=context, packages=[*packages, mf.ic(strt=90.0)])


def _sim(model):
    return mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])


def _how(model, library=None):
    rows = model.move_report(library=library)
    return ", ".join(f"{row['package']} {row['how']}" for row in rows if row["package"] != "ic")


def row(label, call):
    try:
        said = call()
    except Exception as error:  # noqa: BLE001 -- a reproducer reports whatever was raised
        said = f"RAISED {type(error).__name__}"
        text = str(error)
        if isinstance(error, ValueError) and "was declared on mesh" in text:
            said += " (names the recorded mesh)"
        elif isinstance(error, ValueError) and "remap=True" in text:
            said += " (names the move and remap=True)"
    print(f"{label:<47}{said}")


with tempfile.TemporaryDirectory() as scratch:
    scratch = Path(scratch)
    project = Project(scratch / "p", name="p")
    project.add_package("npf/base/cal", _calibrated())
    project.add_package("dis/old", mf.disv(_layering(a)))

    row("eager DISV inline, remap=True",
        lambda: _how(_model(a, mf.disv(_layering(a)), _calibrated()).with_grid(b, remap=True)))
    row("eager DISV by ref + npf by ref, remap=True",
        lambda: _how(_model(a, mf.ref("dis/old"), mf.ref("npf/base/cal"))
                     .with_grid(b, remap=True, library=project), project))

    def _grid_free():
        stack = mf.LayerStack(top=mf.Surface.flat(100.0))
        stack.add("upper", thickness=20.0)
        stack.add("lower", thickness=30.0)
        moved = _model(a, mf.disv(stack), _calibrated()).with_grid(
            b, remap=True, grid_free=("npf",))
        npf = next(p for p in moved.packages if getattr(p, "name", None) == "npf")
        return f"k carried: {not np.array_equal(npf.options['k'], k_old)}"

    row("grid_free=('npf',) under remap=True", _grid_free)

    def _stamped_to_81():
        stack = mf.LayerStack(top=mf.Surface.flat(100.0))
        stack.add("upper", thickness=20.0)
        stack.add("lower", thickness=30.0)
        moved = _model(a, mf.disv(stack), mf.npf(k=k_old, icelltype=1)).with_grid(c)
        mf.simulation(moved, name="s81", packages=[mf.tdis(nper=1)]).build_flopy(scratch / "b81")
        return "BUILT"

    row("stamped npf moved 28 -> 81, built", _stamped_to_81)

    root = scratch / "saved"
    saved = Project(root, name="saved")
    saved.add_package("npf/cal", mf.npf(k=k_old, icelltype=1))
    saved.add_package("dis/old", mf.disv(_layering(a)))
    saved.add_package("dis/same", mf.disv(_layering(b)))
    saved.add_grid("old", mf.GridSpec.from_object(a))
    saved.add_grid("new", mf.GridSpec.from_object(b))
    saved.add_simulation(_sim(_model(a, mf.ref("dis/old"), mf.ref("npf/cal"))))
    saved.add_simulation(_sim(mf.gwf("m", grid=mf.grid_ref("old"), packages=[
        mf.ref("dis/old"), mf.ref("npf/cal"), mf.ic(strt=90.0)])).derive("s_ref"))
    saved.save()
    elsewhere = scratch / "elsewhere"
    elsewhere.mkdir()
    here = os.getcwd()
    os.chdir(elsewhere)
    try:
        q = Project.load(root)
        row("reloaded, simulation-tier with_grid",
            lambda: f"moved models: {q.simulations['s'].with_grid(b, library=q).lineage[-1]['models']}")

        def _model_tier():
            moved = q.simulations["s"].model("m").with_grid(b, library=q)
            moved = moved.with_packages(mf.ref("dis/same"), slot="dis/old")
            run = q.prepare_run("r", q.simulations["s"].with_model(moved), overwrite=True)
            k = run.built.built_model("m").model.get_package("npf").k.array
            return f"BUILT, k == the OLD mesh's K: {np.array_equal(k, k_old)}"

        row("reloaded, model-tier with_grid, built", _model_tier)
        on_ref = q.simulations["s_ref"].model("m")
        row("remap onto project.grid('new') off the root",
            lambda: _how(on_ref.with_grid(q.grid("new"), remap=True, library=q), q))
        row("remap onto mf.grid_ref('new'), library=",
            lambda: _how(on_ref.with_grid(mf.grid_ref("new"), remap=True, library=q), q))
    finally:
        os.chdir(here)

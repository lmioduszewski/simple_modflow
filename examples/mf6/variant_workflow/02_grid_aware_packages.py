"""Pattern 2 -- a grid-aware package that recomputes per mesh.

When the *mesh* varies you cannot reuse a pickled array (it is the wrong size).
Instead the package is a *builder function* that computes from the grid at build
time, so it recomputes for whatever mesh the model is built on. This is the
right tool for mesh-refinement variants.

Run:
    python 02_grid_aware_packages.py
"""

from __future__ import annotations

from pathlib import Path

import flopy
import numpy as np

import myflopy as mf

HERE = Path(__file__).resolve().parent


def grid_aware_npf(model):
    """Compute K from the model's grid at build time (re-runs for each mesh).

    The package builder receives the FloPy model, whose `modelgrid` reflects the
    DIS that was built just before it. In a real model this would call your
    `get_ks(vor)` against the resolved grid.
    """

    grid = model.modelgrid
    col = np.tile(np.arange(grid.ncol), (grid.nrow, 1))
    k = (1.0 + col).reshape(grid.nlay, grid.nrow, grid.ncol)
    return flopy.mf6.ModflowGwfnpf(model, k=k, save_flows=True)


def model(name: str, ncol: int) -> mf.SimulationSpec:
    """A flow model at a chosen resolution; K is computed for that resolution."""

    dis = mf.PackageSpec(
        "dis",
        flopy.mf6.ModflowGwfdis,
        {
            "nlay": 1,
            "nrow": 1,
            "ncol": ncol,
            "delr": 100.0,
            "delc": 100.0,
            "top": 50.0,
            "botm": 0.0,
        },
    )
    chd = mf.chd(
        stress_period_data={0: [[(0, 0, 0), 50.0], [(0, 0, ncol - 1), 40.0]]}
    )
    flow = mf.gwf(
        "flow",
        packages=[
            dis,
            mf.ic(strt=45.0),
            mf.PackageSpec("npf", grid_aware_npf),   # builder runs against the resolved grid
            chd,
            mf.oc(head_filerecord="flow.hds", saverecord=[("HEAD", "ALL")]),
        ],
    )
    return mf.SimulationSpec(
        name,
        models=[flow],
        packages=[
            mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
            mf.ims(models=["flow"], complexity="SIMPLE"),
        ],
    )


project = mf.Project(HERE / "runs" / "grid_aware_demo", name="grid_aware_demo")

# Same npf "recipe", two meshes. K is recomputed for each -- nothing pickled.
coarse = project.run("coarse", model("coarse", ncol=10))
fine = project.run("fine", model("fine", ncol=40))


def ncells(run):
    return run.built.built_model("flow").package("npf").k.array.size


print(f"coarse success={coarse.success}  npf cells={ncells(coarse)}")
print(f"fine   success={fine.success}  npf cells={ncells(fine)}")
print("-> the same builder produced correctly-sized K for each mesh")

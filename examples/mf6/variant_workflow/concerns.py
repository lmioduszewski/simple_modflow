"""Reusable building blocks for the variant-workflow examples.

The "messy" array logic lives here, in small reusable functions, so the
model/variant files stay short and readable. In a real model this is where
your ``get_ks(vor)`` / raster-sampling / zonal-mask code would live.

A simple structured grid stands in for your Voronoi mesh; the patterns are
identical for unstructured grids.
"""

from __future__ import annotations

import numpy as np

import myflopy as mf

# One shared grid definition (stand-in for a built Voronoi mesh).
NLAY, NROW, NCOL = 1, 20, 20
DELR = DELC = 100.0


def dis(top: float = 50.0, botm: float = 0.0) -> mf.PackageSpec:
    """Discretization package.

    Layer top/bottom live in DIS, so *varying layer elevations is a DIS-package
    variant* on a fixed mesh -- not a new grid.
    """

    return mf.dis(
        nlay=NLAY,
        nrow=NROW,
        ncol=NCOL,
        delr=DELR,
        delc=DELC,
        top=top,
        botm=botm,
    )


def _k_array(scale: float) -> np.ndarray:
    """A *computed* K field (stand-in for raster sampling + zonal masks).

    K increases west-to-east. This is the kind of array that cannot be
    JSON-serialized, so the project pickles it when you call ``save()``.
    """

    col = np.tile(np.arange(NCOL), (NROW, 1))
    return ((1.0 + col) * scale).reshape(NLAY, NROW, NCOL)


def npf(scale: float = 1.0) -> mf.PackageSpec:
    """NPF carrying a computed K field -> array-bearing, pickled on save()."""

    return mf.npf(k=_k_array(scale))


def chd(west_head: float = 50.0, east_head: float = 40.0) -> mf.PackageSpec:
    """Constant heads on the west and east edges (drives flow across the grid)."""

    spd = []
    for row in range(NROW):
        spd.append([(0, row, 0), west_head])
        spd.append([(0, row, NCOL - 1), east_head])
    return mf.chd(stress_period_data={0: spd})


def oc(name: str = "flow") -> mf.PackageSpec:
    """Output control saving heads."""

    return mf.oc(head_filerecord=f"{name}.hds", saverecord=[("HEAD", "ALL")])


def k0(run, model: str = "flow") -> float:
    """Helper: first K value of a built model, to show which npf was resolved."""

    npf_pkg = run.built.built_model(model).package("npf")
    return float(np.asarray(npf_pkg.k.array).reshape(-1)[0])


def top0(run, model: str = "flow") -> float:
    """Helper: first DIS top value of a built model, to show a dis swap took."""

    dis_pkg = run.built.built_model(model).package("dis")
    return float(np.asarray(dis_pkg.top.array).reshape(-1)[0])

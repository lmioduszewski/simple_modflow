"""Pattern 3 -- a coupled GWF-GWT variant built from the package library.

The exchange is cheap assembly that lives in the variant function; the reusable
FLOW packages come from the library. Each model writes into its own
subdirectory automatically (runs/.../gwf/ and runs/.../gwt/).

Run:
    python 03_coupled_gwf_gwt.py
"""

from __future__ import annotations

from pathlib import Path

import concerns
import flopy

import myflopy as mf

HERE = Path(__file__).resolve().parent
GRID = {
    "nlay": concerns.NLAY,
    "nrow": concerns.NROW,
    "ncol": concerns.NCOL,
    "delr": concerns.DELR,
    "delc": concerns.DELC,
    "top": 50.0,
    "botm": 0.0,
}

project = mf.Project(HERE / "runs" / "coupled_demo", name="coupled_demo")

# Reusable FLOW packages (computed once, pickled on save()).
project.add_package("dis/base", concerns.dis())
project.add_package("npf/base", concerns.npf(scale=1.0))
project.save()


def chd_with_concentration(inflow_conc: float = 1.0) -> mf.PackageSpec:
    """CHD carrying an auxiliary concentration for transport (coupling needs it)."""

    spd = []
    for row in range(concerns.NROW):
        spd.append([(0, row, 0), 50.0, inflow_conc])              # west: inflow
        spd.append([(0, row, concerns.NCOL - 1), 40.0, 0.0])      # east: outflow
    return mf.chd(stress_period_data={0: spd}, auxiliary="concentration")


def transport_model() -> mf.ModelSpec:
    """A minimal GWT model on the same grid (assembled inline)."""

    return mf.gwt(
        "gwt",
        packages=[
            # ``dis`` (structured) stays a raw PackageSpec -- there is no mf.dis
            # helper (myflopy is Voronoi/DISV-first). ic/oc/adv/mst/ssm use the
            # package-first helpers: mf.ic/mf.oc dispatch on the GWT model kind.
            mf.PackageSpec("dis", flopy.mf6.ModflowGwtdis, dict(GRID)),
            mf.ic(strt=0.0),
            mf.adv(scheme="UPSTREAM"),
            mf.mst(porosity=0.25),
            mf.ssm(sources=[["chd", "AUX", "concentration"]]),
            mf.oc(
                concentration_filerecord="gwt.ucn",
                saverecord=[("CONCENTRATION", "ALL")],
            ),
        ],
    )


def coupled(name: str) -> mf.SimulationSpec:
    """Flow reuses the library; transport + exchange are assembled in code."""

    flow = mf.gwf(
        "gwf",
        packages=[
            mf.ref("dis/base"),
            mf.ic(strt=45.0),
            mf.ref("npf/base"),
            chd_with_concentration(),
            concerns.oc("gwf"),
        ],
    )
    return mf.SimulationSpec(
        name,
        models=[flow, transport_model()],
        packages=[
            mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
            mf.ims(name="ims_gwf", models=["gwf"], complexity="SIMPLE"),
            # transport matrix is asymmetric (advection) -> BICGSTAB, not CG
            mf.ims(
                name="ims_gwt",
                models=["gwt"],
                complexity="SIMPLE",
                linear_acceleration="BICGSTAB",
            ),
        ],
        exchanges=[
            mf.ExchangeSpec("gwfgwt", mf.build_gwf_gwt_exchange, models=("gwf", "gwt")),
        ],
    )


# Build + write, then run. Each model lands in its own subdir.
run = project.prepare_run("coupled_base", coupled("coupled_base"), overwrite=True)
run.write()
print("wrote coupled simulation to", run.workspace)
for sub in ("gwf", "gwt"):
    files = sorted(p.name for p in (run.workspace / sub).glob("*"))
    print(f"  {sub}/: {files}")

success, _ = run.execute()
print("run success:", success)

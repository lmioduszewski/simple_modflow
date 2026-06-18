"""Pattern 1 -- package library + stress/property variants (fixed mesh).

Compute the expensive packages once, save them to disk (array packages are
pickled), then define variants that reuse most packages and swap one. Each
variant is a single readable line: which package it references.

Run:
    python 01_package_library_and_variants.py
"""

from __future__ import annotations

from pathlib import Path

import myflopy as mf

import concerns

HERE = Path(__file__).resolve().parent
project = mf.Project(HERE / "runs" / "library_demo", name="library_demo")

# --- Build the library once. Each package is keyed by a string. ---
project.add_package("dis/base", concerns.dis(top=50.0))      # discretization
project.add_package("npf/base", concerns.npf(scale=1.0))     # computed K -> pickled
project.add_package("npf/high_k", concerns.npf(scale=10.0))  # a K variant
project.add_package("chd/base", concerns.chd())

# Persist the library. npf/* are pickled (arrays); chd/dis serialize to JSON.
project.save()
print("saved library to", project.layout.package_specs_dir)


def flow_model(name: str, npf_key: str) -> mf.SimulationSpec:
    """A flow model assembled from library references. Only `npf_key` varies."""

    flow = mf.gwf(
        "flow",
        packages=[
            mf.ref("dis/base"),
            mf.ic(strt=45.0),
            mf.ref(npf_key),          # <-- the only difference between variants
            mf.ref("chd/base"),
            concerns.oc("flow"),
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


# Each variant is one line. Runs land in runs/library_demo/runs/<name>/.
base = project.run("base", flow_model("base", "npf/base"))
high = project.run("high_k", flow_model("high_k", "npf/high_k"))

print(f"base   success={base.success}  k[0]={concerns.k0(base):.1f}")
print(f"high_k success={high.success}  k[0]={concerns.k0(high):.1f}")

# Reload in a fresh project: the pickled packages come back, no recompute.
reloaded = mf.Project.load(project.root)
print("reloaded packages:", sorted(reloaded.packages))
again = reloaded.run("base_again", flow_model("base_again", "npf/base"))
print(f"base_again success={again.success}  k[0]={concerns.k0(again):.1f}")

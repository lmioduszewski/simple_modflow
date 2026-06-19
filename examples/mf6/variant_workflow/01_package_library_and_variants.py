"""Pattern 1 -- build a model, then derive variants from it, serially.

This mirrors the normal, exploratory workflow: build one complete model and run
it, then -- when you decide something needs to change -- derive a variant FROM
the finished simulation and replace whatever package turns out to be the issue.
You are not locked to one axis (it is not "the npf variant"); you replace the
package you care about at that moment. Re-run a variant, overwriting, after a
fix. Everything is serial.

`derive`/`replace_package` are immutable: each returns a new SimulationSpec, so
the base stays intact and a variant can itself be tweaked and re-run.

Run:
    python 01_package_library_and_variants.py
"""

from __future__ import annotations

from pathlib import Path

import myflopy as mf

import concerns

HERE = Path(__file__).resolve().parent
project = mf.Project(HERE / "runs" / "library_demo", name="library_demo")


# --- 1. Build one complete model and run it. -------------------------------
flow = mf.gwf(
    "flow",
    packages=[
        concerns.dis(top=50.0),
        mf.ic(strt=45.0),
        concerns.npf(scale=1.0),     # computed K field
        concerns.chd(),
        concerns.oc("flow"),
    ],
)
sim = mf.simulation(flow)            # steady TDIS + one IMS for `flow`

# prepare_run builds the FloPy simulation in memory -- no MF6 files yet -- so you
# can inspect it before committing to a run. (project.run(...) is the one-shot
# shortcut that builds + writes + executes in a single call.)
base = project.prepare_run("base", sim, overwrite=True)
flow_view = base.model("flow")       # ModelView: raw .gwf/.sim + your plotting/explorer tools
print("inspect:", list(flow_view.gwf.package_dict), "k[0] =", concerns.k0(base))
#   ^ open your choropleths / package_explorer here, then run when satisfied
base.execute()
print(f"base    success={base.success}  k[0]={concerns.k0(base):.1f}")


# --- 2. Derive a variant from the finished model, inspect it, then run it.
#        Same loop as the base: prepare (build in memory) -> look -> execute.
high = sim.derive("high_k").replace_package("flow", concerns.npf(scale=10.0))
high_run = project.prepare_run("high_k", high, overwrite=True)
print("inspect high_k: k[0] =", concerns.k0(high_run))   # plot/check before running
high_run.execute()
print(f"high_k  success={high_run.success}  k[0]={concerns.k0(high_run):.1f}")


# --- 3. That variant needs a tweak: adjust the same package, re-run (overwrite).
high = high.replace_package("flow", concerns.npf(scale=25.0))
high_run = project.run("high_k", high, overwrite=True)
print(f"high_k  success={high_run.success}  k[0]={concerns.k0(high_run):.1f}  (re-run)")


# --- 4. Not locked to npf: derive another variant changing a DIFFERENT package.
shallow = sim.derive("shallow").replace_package("flow", concerns.dis(top=30.0))
shallow_run = project.run("shallow", shallow, overwrite=True)
print(f"shallow success={shallow_run.success}  top[0]={concerns.top0(shallow_run):.1f}")


# --- Reuse: keep expensive computed packages in the library, swap by reference.
# Compute once, save (array packages are pickled), reuse across sessions/variants
# via mf.ref -- no recompute. The authoring stays serial: derive, swap, run.
project.add_package("npf/high_k", concerns.npf(scale=10.0))
project.save()

reloaded = mf.Project.load(project.root)
from_lib = sim.derive("from_library").replace_package("flow", mf.ref("npf/high_k"))
lib_run = reloaded.run("from_library", from_lib, overwrite=True)
print(f"from_library success={lib_run.success}  k[0]={concerns.k0(lib_run):.1f}  (from pickle)")

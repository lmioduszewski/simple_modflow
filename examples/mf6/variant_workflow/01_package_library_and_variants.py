"""Pattern 1 -- build a model, then derive variants from it, serially.

This mirrors the normal, exploratory workflow: build one complete model and run
it, then -- when you decide something needs to change -- derive a variant FROM
the finished simulation and replace whatever package turns out to be the issue.
You are not locked to one axis (it is not "the npf variant"); you replace the
package you care about at that moment. Re-run a variant, overwriting, after a
fix. Everything is serial.

`derive`/`replace_package` are immutable: each returns a new SimulationSpec, so
the base stays intact and a variant can itself be tweaked and re-run.

Section 5 does the same work through the LIBRARY verbs (`packages.fork`/`edit`,
`simulations.fork(swap=)`), which register what they make and record where it
came from -- and shows `project.runs.drift()` naming the runs an edit
invalidated.

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


# --- 5. The same work through the library verbs (2026-09, ledger 186/187). ---
# `derive`/`replace_package` above are the SPEC-level, immutable spelling: they
# hand you a new object and leave registering it to you. The library verbs are
# the registering, mutating spelling -- same mechanism underneath, and they
# record where each variant came from.
#
# Fork a package to a new key, edit one in place, and fork a simulation by
# swapping a library key for another:
project.packages.fork("npf/high_k", "npf/very_high_k", note="second pass")

# Register the ref-based variant, then fork IT by swapping one library key for
# another. `swap` matches the LABEL an entry carries, which for `mf.ref` is its
# library key -- so this targets the npf reference and leaves everything else.
project.add_simulation(from_lib)
v2 = project.simulations.fork(
    "from_library", "from_library_v2", swap={"npf/high_k": "npf/very_high_k"}
)
print()
print("from_library_v2 packages:",
      [getattr(e, "key", None) or e.name for e in v2.model("flow").packages])
print("derived_from:", v2.derived_from, "| lineage:", [s["operation"] for s in v2.lineage])

# `edit` changes what an EXISTING key means -- the operation that can invalidate
# a run that already used it. Provenance alone does NOT count as a change: the
# fingerprint covers what a package BUILDS, so stamping a note leaves every run
# reading `unchanged`.
project.packages.edit("npf/high_k", note="reviewed 2026-09")
print()
print("after a note-only edit:",
      project.refresh_runs().drift("from_library").iloc[0]["state"])

# Changing a VALUE does count, and names the runs it reaches.
project.packages.edit("npf/high_k", k=[[42.0]], note="post-calibration K")

# What is defined, and what each one is:
print()
print(project)
print(project.packages.summary()[["key", "package", "builder"]].to_string(index=False))

# What has been run, what MF6 said about it, and what an edit invalidated.
# `status` is myflopy's verb; `outcome` is MF6's verdict from mfsim.lst.
print()
print(project.runs.summary()[["key", "status", "outcome", "drift"]].to_string(index=False))
print()
print(project.runs.drift().to_string(index=False))

# attic

Retired and experimental code kept for reference only. **Nothing here is part of
the importable `myflopy` package or the built wheel** — packaging discovers
modules under `src/` only (`[tool.setuptools.packages.find] where = ["src"]`),
and nothing in `src/` imports from this directory.

Moved out of `src/myflopy/modflow/mf6/` during the framework-hardening pass so
the importable source tree contains only shipped code.

- `archive/` — earlier package implementations superseded by the current code.
- `maybe_junk/` — one-off scripts, notebooks, and spreadsheets (pump run-time
  calcs, a Bokeh drawing experiment, an old heads-plotting module). Kept in case
  a snippet is still wanted; not maintained.

If something here turns out to be dead for good, delete it — git history
preserves it.

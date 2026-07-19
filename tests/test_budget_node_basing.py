"""Budget node ids are zero-based on EVERY path, for every cell-based BC.

Regression guard for the coupled off-by-one found by the 2026-07-18 review:

* ``budget_tables._zero_base_budget_frame`` normalized a hardcoded
  ``{"drn", "ghb", "rch"}`` only, so ``model.bud("chd"|"riv"|"wel"|"evt").df``
  handed back MF6's raw 1-based node ids.
* ``GroupBudget._normalize_nodes`` then re-normalized with an
  ``if node.min() >= 1: subtract 1`` heuristic. That heuristic is correct on RAW
  MF6 records (whose nodes are always >= 1) but was being applied to output that
  ``budget_df`` had ALREADY zero-based -- so ``group.bud("drn"|"ghb")`` came back
  one cell low, while ``rch`` escaped only because it covers node 0.

The two bugs cancelled for chd/riv/wel/evt, so fixing either alone traded one
off-by-one for another. These tests pin both paths against ground truth taken
from the packages' own stress-period data.
"""

from __future__ import annotations

import pytest

from myflopy.modflow.mf6.package_registry import _PACKAGE_EXPLORER_SPECS
from myflopy.project.model_group import ModelGroup

CELL_BC_PACKAGES = sorted(
    name for name, spec in _PACKAGE_EXPLORER_SPECS.items() if spec.kind == "cell_stress"
)


def _truth_nodes(model, package: str) -> set[int]:
    """Zero-based flat node ids a package actually occupies, from its own inputs."""

    ncpl = int(model.vor.ncpl)
    records = model.gwf.get_package(package).stress_period_data.get_data(key=0)
    return {int(row[0][0]) * ncpl + int(row[0][1]) for row in records}


def test_cell_bc_registry_covers_the_packages_under_test():
    # if a new cell BC is added, it is covered here automatically -- and if the
    # registry ever loses one, this fails rather than silently shrinking coverage
    assert set(CELL_BC_PACKAGES) == {"chd", "drn", "evt", "ghb", "rch", "riv", "wel"}


@pytest.mark.parametrize("package", CELL_BC_PACKAGES)
def test_model_budget_nodes_are_zero_based(canonical_run, package):
    """``model.bud(pkg).df`` node ids must be zero-based model cells."""

    frame = canonical_run.bud(package).df.reset_index()
    column = "node" if "node" in frame.columns else frame.columns[0]
    nodes = {int(value) for value in frame[column].dropna().unique()}

    assert nodes, f"{package} budget produced no nodes"
    assert nodes <= _truth_nodes(canonical_run, package), (
        f"{package}: model.bud() node ids are not the package's zero-based cells "
        f"(smallest returned {min(nodes)}, smallest true {min(_truth_nodes(canonical_run, package))})"
    )


@pytest.mark.parametrize("package", CELL_BC_PACKAGES)
def test_group_budget_nodes_match_the_single_model_path(canonical_run, package):
    """``group.bud(pkg)`` must agree with ``model.bud(pkg)`` -- no re-normalization."""

    group = ModelGroup({"a": canonical_run}, reference="a")
    grouped = {int(v) for v in group.bud(package).get()["node"].dropna().unique()}

    assert grouped, f"{package} group budget produced no nodes"
    assert grouped <= _truth_nodes(canonical_run, package), (
        f"{package}: group.bud() node ids drifted from the package's zero-based "
        f"cells (smallest returned {min(grouped)})"
    )

    frame = canonical_run.bud(package).df.reset_index()
    column = "node" if "node" in frame.columns else frame.columns[0]
    single = {int(value) for value in frame[column].dropna().unique()}
    assert grouped == single, (
        f"{package}: the group and single-model budget paths disagree "
        f"(group-only {sorted(grouped - single)[:5]}, model-only {sorted(single - grouped)[:5]})"
    )


@pytest.mark.parametrize("package", CELL_BC_PACKAGES)
def test_modern_explorer_path_agrees_with_the_legacy_budget_path(canonical_run, package):
    """The registry-backed explorer and the legacy ``bud()`` path must not disagree."""

    ncpl = int(canonical_run.vor.ncpl)
    results = getattr(canonical_run.packages, package).results.q.get()
    explorer = {
        int(row.layer) * ncpl + int(row.cell) for row in results.itertuples()
    }

    frame = canonical_run.bud(package).df.reset_index()
    column = "node" if "node" in frame.columns else frame.columns[0]
    legacy = {int(value) for value in frame[column].dropna().unique()}

    assert explorer and legacy
    assert explorer == legacy, (
        f"{package}: explorer and legacy budget paths disagree by "
        f"{sorted(explorer ^ legacy)[:5]}"
    )

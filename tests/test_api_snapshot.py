"""The package-API golden snapshot (plan 4.7.0) -- the 4.7 refactor's safety net.

``tests/api_snapshot.json`` records the public shape of the package API BEFORE
the 4.7 consolidation: helper signatures, ``*_spec`` outputs from fixed inputs,
and the package set visible at every surface that encodes package knowledge.
The consolidation collapses 11 spec factories, 11 helper classes and 7
GeoPackage resolvers into parameterized implementations; without this file,
"behaviour-preserving" would be a claim rather than something anyone can check.

If a test here fails:

* **unintended** -- the refactor changed the public surface; fix the code.
* **intended** -- rerun ``python scripts/derive_api_snapshot.py`` and REVIEW THE
  DIFF in the pull request. The diff is the change log for the public API.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SNAPSHOT_PATH = REPO_ROOT / "tests" / "api_snapshot.json"
DERIVE_SCRIPT = REPO_ROOT / "scripts" / "derive_api_snapshot.py"


@pytest.fixture(scope="module")
def snapshot() -> dict:
    return json.loads(SNAPSHOT_PATH.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def rebuilt() -> dict:
    sys.path.insert(0, str(DERIVE_SCRIPT.parent))
    from derive_api_snapshot import build_snapshot

    return build_snapshot()


def test_snapshot_matches_the_current_api(snapshot, rebuilt):
    """The checked-in snapshot must equal what the code produces right now."""

    for dimension in ("signatures", "specs", "surfaces"):
        current = rebuilt[dimension]
        stored = snapshot[dimension]

        added = sorted(set(current) - set(stored))
        removed = sorted(set(stored) - set(current))
        assert not added and not removed, (
            f"{dimension}: entries appeared/disappeared "
            f"(added={added}, removed={removed}). If intended, rerun "
            "scripts/derive_api_snapshot.py and review the diff."
        )

        changed = [key for key in stored if stored[key] != current[key]]
        assert not changed, (
            f"{dimension}: {changed} changed shape.\n"
            f"  stored : {json.dumps({k: stored[k] for k in changed[:2]}, indent=2)[:900]}\n"
            f"  current: {json.dumps({k: current[k] for k in changed[:2]}, indent=2)[:900]}\n"
            "If intended, rerun scripts/derive_api_snapshot.py and review the diff."
        )


def test_derive_script_is_idempotent():
    """``--check`` must pass, so CI catches a snapshot nobody regenerated."""

    result = subprocess.run(
        [sys.executable, str(DERIVE_SCRIPT), "--check"],
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_snapshot_covers_the_whole_public_package_surface(snapshot):
    """Guard against the net silently SHRINKING.

    A snapshot that quietly stops covering something still passes the equality
    test above, which would make it worthless exactly when it matters. Assert it
    spans every public export instead.
    """

    from myflopy import advanced, package_api

    missing_helpers = sorted(set(package_api.__all__) - set(snapshot["signatures"]))
    assert not missing_helpers, f"helpers absent from the snapshot: {missing_helpers}"

    missing_specs = sorted(set(advanced.__all__) - set(snapshot["specs"]))
    assert not missing_specs, f"spec factories absent from the snapshot: {missing_specs}"

    # every helper records at least its call signature
    for helper, methods in snapshot["signatures"].items():
        assert methods, f"{helper} has no recorded methods"

    # every list BC records all three documented entry points
    for helper in ("chd", "ghb", "drn", "riv", "wel", "rch", "evt"):
        assert {"__call__", "flopy", "gpkg"} <= set(snapshot["signatures"][helper]), helper


def test_surfaces_dimension_names_the_known_inconsistencies(snapshot):
    """The surfaces block is 4.7's site inventory -- it must stay comprehensive.

    These entries deliberately DISAGREE with one another today (that is the
    problem 4.7 removes). Pinning the disagreement means the refactor cannot
    quietly change a surface, in either direction.
    """

    surfaces = snapshot["surfaces"]
    required = {
        "registry.cell_stress",
        "package_api.__all__",
        "myflopy.__preferred__",
        "ModelPackages.properties",
        "GroupPackages.properties",
        "PackageDiffNamespace.properties",
        "model_diff._DIFF_PACKAGES",
        "results_diff._CELL_BUDGET_PACKAGES",
        "run_model._PACKAGE_SUFFIX_TO_TYPE",
        "budget_tables.cell_based",
        "components.SUPPORTED_ARTIFACT_TYPES",
    }
    assert required <= set(surfaces), f"missing surfaces: {sorted(required - set(surfaces))}"

    cell_bcs = set(surfaces["registry.cell_stress"])
    # surfaces that 4.7.1 already brought in line with the registry stay that way
    for name in (
        "model_diff._DIFF_PACKAGES",
        "budget_tables.cell_based",
        "components.SUPPORTED_ARTIFACT_TYPES",
    ):
        assert cell_bcs <= set(surfaces[name]), (
            f"{name} regressed: missing {sorted(cell_bcs - set(surfaces[name]))}"
        )

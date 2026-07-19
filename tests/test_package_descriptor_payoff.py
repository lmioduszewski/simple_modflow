"""The 4.7 payoff test: what does adding ONE registry entry actually get you?

Phase 4.7 exists because `mf.riv`/`mf.evt` followed a documented four-piece
checklist, under two adversarial reviews, and still missed six sites. The fix
was to make per-package knowledge derive from one descriptor. This test is how
that claim stops being a claim.

It registers a synthetic package descriptor and reports, surface by surface,
where it shows up on its own and where a human still has to type something. Two
assertions follow, and the second is the interesting one:

* every surface listed as AUTOMATIC must stay automatic -- a regression there
  re-opens the exact defect class 4.7 closed;
* the surfaces listed as MANUAL must be EXACTLY that set -- so a new
  hand-written site fails immediately, and closing one *also* fails, forcing
  the win to be recorded rather than absorbed silently.

**Why a subprocess.** The derived lists are module-level constants evaluated at
import. Injecting a descriptor into an already-imported process would not
propagate, and `importlib.reload` would leave other modules holding stale
references -- dangerous under `-n 10` with a session-scoped canonical model. A
fresh interpreter that injects BEFORE importing any consumer is the only honest
way to test "adding a descriptor entry makes it appear".
"""

from __future__ import annotations

import json
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]

#: Surfaces that a registry entry reaches with no other edit. Derived in 4.7.3.
#: This set may GROW; it must never shrink.
AUTOMATIC_SURFACES = frozenset(
    {
        "model_diff._DIFF_PACKAGES",
        "results_diff._CELL_BUDGET_PACKAGES",
        "results_diff._MOVER_PACKAGES",
        "budget_tables.cell_based",
        "run_model.suffix_map",
        "components.LIST_BC",
        "components.SUPPORTED",
        "components.APPLY_ORDER",
        "registry.budget_term",
    }
)

#: Surfaces that still require hand-written code. Two distinct reasons:
#:
#: * the first three are DELIBERATE (ledger 38) -- named factories, helper
#:   classes and resolvers carry per-package prose and explicit signatures, and
#:   generating them would cost documentation quality and IDE completion;
#: * the last four are namespace PROPERTIES, hand-written on purpose (4.7.7).
#:   Generating them would need a `.pyi` stub to stay visible to type checkers
#:   -- and a stub replaces the WHOLE module, so it would mean hand-maintaining
#:   stubs for every other public name in four large modules. Completeness is
#:   enforced by assertion instead: see
#:   `test_every_namespace_exposes_every_registry_package`.
MANUAL_SURFACES = frozenset(
    {
        "advanced.zzz_spec",
        "package_api.zzz",
        "geopackage.zzz_resolver",
        "ModelPackages.zzz",
        "SimulationBase.zzz",
        "GroupPackages.zzz",
        "PackageDiffNamespace.zzz",
    }
)

_PROBE = textwrap.dedent(
    '''
    import json

    # Inject a synthetic descriptor BEFORE importing any consumer, so the
    # module-level derivations see it exactly as they would a real addition.
    from myflopy.modflow.mf6.package_registry import (
        _PACKAGE_EXPLORER_SPECS, PackageExplorerSpec, PackageCapabilities,
        PackageTiers, FieldSpec, ResultSpec,
    )

    _PACKAGE_EXPLORER_SPECS["zzz"] = PackageExplorerSpec(
        name="zzz", kind="cell_stress", default_input="elev", colorscale="earth",
        inputs={"elev": FieldSpec("elev", label="Synthetic", colorscale="earth")},
        results={"q": ResultSpec("q", budget_text="ZZZ", value_name="q", colorscale="RdBu")},
        # a real FloPy class: the descriptor is exercised, not FloPy itself
        flopy_class="ModflowGwfdrn",
        record_fields=("elev", "cond"),
        gpkg_defaults={"elevation": "elevation", "conductance": "conductance"},
        capabilities=PackageCapabilities(mover=True, edges_only=True),
        tiers=PackageTiers(
            diffable=True, results_diffable=True, artifact_serializable=True,
            artifact_apply_order=55, model_accessor=False,
        ),
        file_suffix="zzz", zero_base_budget_nodes=True,
        blurb="A synthetic package that exists only inside this test.",
    )

    found = {}

    from myflopy.project.model_diff import _DIFF_PACKAGES, _PackageDiffNamespace
    found["model_diff._DIFF_PACKAGES"] = "zzz" in _DIFF_PACKAGES

    from myflopy.project.model_results_diff import _CELL_BUDGET_PACKAGES, MvrResultDiff
    found["results_diff._CELL_BUDGET_PACKAGES"] = "zzz" in _CELL_BUDGET_PACKAGES
    found["results_diff._MOVER_PACKAGES"] = "zzz" in MvrResultDiff._MOVER_PACKAGES

    from myflopy.modflow.mf6.budget_tables import _cell_based_budget_packages
    found["budget_tables.cell_based"] = "zzz" in _cell_based_budget_packages()

    from myflopy.project.run_model import _PACKAGE_SUFFIX_TO_TYPE
    found["run_model.suffix_map"] = "zzz" in _PACKAGE_SUFFIX_TO_TYPE

    from myflopy.project.components import (
        LIST_BC_ARTIFACT_TYPES, SUPPORTED_PACKAGE_ARTIFACT_TYPES,
        PACKAGE_ARTIFACT_APPLY_ORDER,
    )
    found["components.LIST_BC"] = "zzz" in LIST_BC_ARTIFACT_TYPES
    found["components.SUPPORTED"] = "zzz" in SUPPORTED_PACKAGE_ARTIFACT_TYPES
    found["components.APPLY_ORDER"] = "zzz" in PACKAGE_ARTIFACT_APPLY_ORDER

    from myflopy.modflow.mf6.package_registry import get_default_budget_term
    found["registry.budget_term"] = get_default_budget_term("zzz") is not None

    # ...and the surfaces a human still has to touch
    from myflopy import advanced, package_api
    found["advanced.zzz_spec"] = hasattr(advanced, "zzz_spec")
    found["package_api.zzz"] = hasattr(package_api, "zzz")

    from myflopy.geopackage import GeoPackageSource
    found["geopackage.zzz_resolver"] = hasattr(GeoPackageSource, "zzz")

    from myflopy.modflow.mf6.package_model import ModelPackages
    from myflopy.modflow.mf6.simulation.base import SimulationBase
    from myflopy.project.group.packages import GroupPackages

    def is_property(cls):
        return isinstance(getattr(cls, "zzz", None), property)

    found["ModelPackages.zzz"] = is_property(ModelPackages)
    found["SimulationBase.zzz"] = is_property(SimulationBase)
    found["GroupPackages.zzz"] = is_property(GroupPackages)
    found["PackageDiffNamespace.zzz"] = is_property(_PackageDiffNamespace)

    print("PAYOFF_JSON:" + json.dumps(found))
    '''
)


@pytest.fixture(scope="module")
def reach() -> dict[str, bool]:
    """Run the probe in a fresh interpreter and return ``{surface: reached}``."""

    result = subprocess.run(
        [sys.executable, "-c", _PROBE],
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
        timeout=300,
    )
    assert result.returncode == 0, (
        f"probe failed:\n{result.stdout[-2000:]}\n{result.stderr[-3000:]}"
    )
    line = next(
        line for line in result.stdout.splitlines() if line.startswith("PAYOFF_JSON:")
    )
    return json.loads(line[len("PAYOFF_JSON:") :])


def test_one_descriptor_entry_reaches_every_derived_surface(reach):
    """The 4.7 guarantee: these surfaces need no second edit.

    If this fails, a subsystem went back to carrying its own package list and
    the riv/evt six-site miss is possible again.
    """

    missing = sorted(name for name in AUTOMATIC_SURFACES if not reach.get(name))
    assert not missing, (
        "these surfaces no longer pick up a new registry entry on their own: "
        f"{missing} -- a derivation was replaced by a hardcoded list"
    )


def test_the_hand_written_surfaces_are_exactly_the_known_set(reach):
    """A ratchet in both directions.

    Fails if a NEW hand-written site appears (the defect class returning), and
    equally if one is closed without updating this list -- so progress has to be
    recorded rather than absorbed silently. Closing one is a good failure: fix
    it by moving the name out of ``MANUAL_SURFACES``.
    """

    still_manual = {name for name, reached in reach.items() if not reached}

    unexpected = sorted(still_manual - MANUAL_SURFACES)
    assert not unexpected, (
        f"new hand-written surfaces appeared: {unexpected} -- adding a package "
        "now needs edits that 4.7 was supposed to have removed"
    )

    now_automatic = sorted(MANUAL_SURFACES - still_manual)
    assert not now_automatic, (
        f"these are automatic now: {now_automatic}. Good -- remove them from "
        "MANUAL_SURFACES and note it in the plan/ledger."
    )


def test_the_probe_actually_covers_both_kinds_of_surface(reach):
    """Guard the guard: a probe that checked nothing would pass both tests."""

    assert set(reach) == AUTOMATIC_SURFACES | MANUAL_SURFACES
    assert len(AUTOMATIC_SURFACES) >= 9, "the automatic set shrank"
    assert any(reach.values()) and not all(reach.values()), (
        "the probe should find a mix; all-true or all-false means it is broken"
    )

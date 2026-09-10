"""Regenerate the package-API golden snapshot (plan 4.7.0).

The snapshot is the safety net for the 4.7 consolidation: it records the public
shape of the package API *before* the refactor so that "behaviour-preserving"
can be proven rather than asserted. Three dimensions:

1. **signatures** -- every ``mf.<pkg>`` helper's public methods and their full
   parameter lists (names, kinds, defaults). A collapsed helper class must keep
   these byte-identical.
2. **specs** -- every ``*_spec`` factory's output from FIXED inputs: package
   name, builder identity, and a compact option summary (reusing
   ``PackageSpec.option_summary``, so data payloads are summarized rather than
   dumped).
3. **surfaces** -- the package set visible at each place that encodes package
   knowledge. This doubles as the machine-readable version of 4.7's site
   inventory, and makes the CURRENT inconsistencies explicit.

Usage::

    python scripts/derive_api_snapshot.py          # rewrite tests/api_snapshot.json
    python scripts/derive_api_snapshot.py --check   # exit 1 if it would change

Intentional API changes are landed by rerunning this and REVIEWING THE DIFF.
"""

from __future__ import annotations

import argparse
import inspect
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SNAPSHOT_PATH = REPO_ROOT / "tests" / "api_snapshot.json"

sys.path.insert(0, str(REPO_ROOT / "src"))


def _default_repr(value) -> str:
    """Stable text for a parameter default (reprs of sentinels stay readable)."""

    if value is inspect.Parameter.empty:
        return "<required>"
    return repr(value)


def _signature(func) -> list[str]:
    """``name:kind=default`` for each parameter, in declaration order."""

    try:
        signature = inspect.signature(func)
    except (TypeError, ValueError):  # pragma: no cover - builtins
        return ["<unavailable>"]
    return [
        f"{name}:{parameter.kind.name}={_default_repr(parameter.default)}"
        for name, parameter in signature.parameters.items()
        if name != "self"
    ]


def collect_signatures() -> dict[str, dict[str, list[str]]]:
    """Public method signatures for every package-first helper singleton."""

    import myflopy as mf
    from myflopy import package_api

    snapshot: dict[str, dict[str, list[str]]] = {}
    for name in sorted(package_api.__all__):
        helper = getattr(mf, name)
        if inspect.isfunction(helper) or inspect.isbuiltin(helper):
            snapshot[name] = {"__call__": _signature(helper)}
            continue
        methods: dict[str, list[str]] = {}
        for attribute in sorted(dir(type(helper))):
            if attribute.startswith("_") and attribute != "__call__":
                continue
            member = getattr(type(helper), attribute, None)
            if callable(member):
                methods[attribute] = _signature(member)
        snapshot[name] = methods
    return snapshot


# Fixed inputs -- deliberately tiny and literal so the snapshot is stable.
_CELL = (0, 0)
_SPEC_INPUTS: dict[str, tuple[tuple, dict]] = {
    "chd_spec": (({0: [[_CELL, 1.0]]},), {}),
    "drn_spec": (({0: [[_CELL, 1.0, 2.0]]},), {}),
    "ghb_spec": (({0: [[_CELL, 1.0, 2.0]]},), {}),
    "riv_spec": (({0: [[_CELL, 1.0, 2.0, 0.5]]},), {}),
    "wel_spec": (({0: [[_CELL, -1.0]]},), {}),
    "rch_spec": (({0: [[_CELL, 1.0e-4]]},), {}),
    "evt_spec": (({0: [[_CELL, 10.0, 1.0e-4, 2.0]]},), {}),
    "uzf_spec": (([[0, _CELL, 1, -1, 0.0, 0.1, 0.05, 0.3, 0.1, 3.5]], {0: [[0, 0.001]]}), {}),
    "lak_spec": (
        ([[0, 9.0, 1]], [[0, 0, _CELL, "VERTICAL", 1.0e-5, 0.0, 0.0, 0.0, 0.0]], {0: [[0, "STATUS", "ACTIVE"]]}),
        {},
    ),
    "sfr_spec": (
        ([[0, _CELL, 1.0, 1.0, 0.001, 9.0, 1.0, 1.0, 0.03, 0, 1.0, 0]], [[0]], {0: [[0, "STATUS", "ACTIVE"]]}),
        {},
    ),
    "mvr_spec": ((([["sfr"], ["lak"]]), {0: [["sfr", 0, "lak", 0, "FACTOR", 1.0]]}), {}),
}


def collect_specs() -> dict[str, dict]:
    """Each ``*_spec`` factory's output for a fixed, literal input."""

    from myflopy import advanced
    from myflopy.specs import _builder_label

    snapshot: dict[str, dict] = {}
    for factory_name in sorted(advanced.__all__):
        factory = getattr(advanced, factory_name)
        args, kwargs = _SPEC_INPUTS[factory_name]
        spec = factory(*args, **kwargs)
        builder = spec.builder
        snapshot[factory_name] = {
            "name": spec.name,
            "builder": _builder_label(builder),
            "builder_target": getattr(
                getattr(builder, "args", [None])[0], "__name__", "<none>"
            ),
            "requires": list(spec.requires),
            "options": dict(sorted(spec.option_summary.items())),
        }
    return snapshot


def collect_surfaces() -> dict[str, list[str]]:
    """The package set visible at every place that encodes package knowledge.

    This is the machine-readable form of the 4.7 site inventory. Entries here
    are EXPECTED to disagree with one another today -- that disagreement is the
    thing 4.7 removes, and pinning it means the refactor cannot change a surface
    without the diff showing it.
    """

    import myflopy as mf
    from myflopy import package_api
    from myflopy.modflow.mf6.budget_tables import _cell_based_budget_packages
    from myflopy.modflow.mf6.package_model import ModelPackages
    from myflopy.modflow.mf6.package_registry import _PACKAGE_EXPLORER_SPECS
    from myflopy.modflow.mf6.simulation.base import SimulationBase
    from myflopy.project.components import (
        PACKAGE_ARTIFACT_APPLY_ORDER,
        SUPPORTED_PACKAGE_ARTIFACT_TYPES,
    )
    from myflopy.project.group.packages import GroupPackages
    from myflopy.project.model_diff import _DIFF_PACKAGES, _PackageDiffNamespace
    from myflopy.project.model_results_diff import _CELL_BUDGET_PACKAGES, MvrResultDiff
    from myflopy.project.run_model import _PACKAGE_SUFFIX_TO_TYPE

    def properties(cls) -> list[str]:
        return sorted(
            name
            for name, value in vars(cls).items()
            if isinstance(value, property) and not name.startswith("_")
        )

    def members(cls) -> list[str]:
        """Every public member a class DEFINES -- methods and properties alike.

        The project/run surface was covered by nothing until 2026-09-06, which is
        how ``ProjectSpec.package_spec(key)`` disappeared in the 2026-06-17 fold
        into ``Project`` without a single test failing (ledger 182). Signatures
        are deliberately not captured here: these are hand-written verbs whose
        parameters change often, and the value is in noticing a verb APPEARING or
        VANISHING.
        """

        return sorted(
            name
            for name, value in vars(cls).items()
            if not name.startswith("_")
            and (callable(value) or isinstance(value, (property, classmethod, staticmethod)))
        )

    def spec_fields(cls) -> list[str]:
        """The dataclass FIELDS a spec declares, which ``members`` cannot see.

        A field vanishing is at least as load-bearing as a verb vanishing, and
        nothing covered it. ``PackageSpec.enabled`` is the worked example:
        it participates in ``package_digest``, so removing it would make every
        historical run manifest report ``state='edited'`` -- and the vocabulary
        work has ``disabled()`` (the METHOD) on its cut list, one careless step
        away from the field.
        """

        import dataclasses

        return sorted(field.name for field in dataclasses.fields(cls))

    from myflopy.library import (
        GridLibrary,
        PackageLibrary,
        SimulationLibrary,
        SpecLibrary,
    )
    from myflopy.specs import (
        GridSpec,
        ModelContext,
        ModelSpec,
        PackageSpec,
        SimulationSpec,
    )
    from myflopy.workspace import Project, Run, RunLibrary, SaveReport

    return {
        "registry.cell_stress": sorted(
            name for name, spec in _PACKAGE_EXPLORER_SPECS.items()
            if spec.kind == "cell_stress"
        ),
        "registry.all": sorted(_PACKAGE_EXPLORER_SPECS),
        "package_api.__all__": sorted(package_api.__all__),
        "myflopy.__preferred__": sorted(mf.__preferred__),
        "myflopy.__engine__": sorted(mf.__engine__),
        "advanced.__all__": sorted(__import__("myflopy.advanced", fromlist=["x"]).__all__),
        "ModelPackages.properties": properties(ModelPackages),
        "GroupPackages.properties": properties(GroupPackages),
        "PackageDiffNamespace.properties": properties(_PackageDiffNamespace),
        "SimulationBase.package_properties": properties(SimulationBase),
        "model_diff._DIFF_PACKAGES": sorted(_DIFF_PACKAGES),
        "results_diff._CELL_BUDGET_PACKAGES": sorted(_CELL_BUDGET_PACKAGES),
        "results_diff._MOVER_PACKAGES": sorted(MvrResultDiff._MOVER_PACKAGES),
        "run_model._PACKAGE_SUFFIX_TO_TYPE": sorted(_PACKAGE_SUFFIX_TO_TYPE),
        "budget_tables.cell_based": sorted(_cell_based_budget_packages()),
        "components.SUPPORTED_ARTIFACT_TYPES": sorted(SUPPORTED_PACKAGE_ARTIFACT_TYPES),
        "components.APPLY_ORDER": sorted(PACKAGE_ARTIFACT_APPLY_ORDER),
        # The project/run surface. Covered by nothing before 2026-09-06.
        "Project.members": members(Project),
        "Run.members": members(Run),
        "SaveReport.members": members(SaveReport),
        "SpecLibrary.members": members(SpecLibrary),
        "PackageLibrary.members": members(PackageLibrary),
        "GridLibrary.members": members(GridLibrary),
        "SimulationLibrary.members": members(SimulationLibrary),
        "RunLibrary.members": members(RunLibrary),
        # The SPEC surface. Covered by nothing before 2026-09-10: `SimulationSpec`
        # appeared in this file exactly once, as a bare name in `mf.__all__`, so a
        # verb or field could appear or vanish on the declarative layer without a
        # single test noticing. Verbs and fields both, because a field going is at
        # least as consequential as a method going.
        "PackageSpec.members": members(PackageSpec),
        "PackageSpec.fields": spec_fields(PackageSpec),
        "ModelSpec.members": members(ModelSpec),
        "ModelSpec.fields": spec_fields(ModelSpec),
        "SimulationSpec.members": members(SimulationSpec),
        "SimulationSpec.fields": spec_fields(SimulationSpec),
        "ModelContext.members": members(ModelContext),
        "ModelContext.fields": spec_fields(ModelContext),
        "GridSpec.members": members(GridSpec),
        "GridSpec.fields": spec_fields(GridSpec),
    }


def build_snapshot() -> dict:
    """Assemble the full golden snapshot."""

    return {
        "_comment": (
            "Golden snapshot of the package API (plan 4.7.0). Regenerate with "
            "scripts/derive_api_snapshot.py and REVIEW THE DIFF -- an unexpected "
            "change here means a refactor altered the public surface."
        ),
        "signatures": collect_signatures(),
        "specs": collect_specs(),
        "surfaces": collect_surfaces(),
    }


def main() -> int:
    """Write (or check) the snapshot file."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="exit 1 if the snapshot would change")
    arguments = parser.parse_args()

    snapshot = build_snapshot()
    text = json.dumps(snapshot, indent=2, sort_keys=True) + "\n"

    if arguments.check:
        current = SNAPSHOT_PATH.read_text(encoding="utf-8") if SNAPSHOT_PATH.exists() else ""
        if current != text:
            print("api snapshot is STALE - rerun scripts/derive_api_snapshot.py")
            return 1
        print("api snapshot is up to date")
        return 0

    SNAPSHOT_PATH.write_text(text, encoding="utf-8")
    counts = {key: len(value) for key, value in snapshot.items() if isinstance(value, dict)}
    print(f"wrote {SNAPSHOT_PATH.relative_to(REPO_ROOT)} - {counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

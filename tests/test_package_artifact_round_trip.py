"""Every cell-based list BC survives an artifact capture -> apply round trip.

Regression guard for the 2026-07-18 sweep finding: the artifact subsystem
(`project/components.py`) carried its own hardcoded package list, so

* ``wel`` was unsupported for its entire life -- ``Project.add_package`` could
  never round-trip a well, and nobody noticed;
* ``riv``/``evt`` were unsupported from birth;
* ``chd``/``drn`` CAPTURED ``auxiliary`` but silently dropped it on restore,
  because each package had a hand-written restore block.

The set is now derived from the package registry, so these tests read it from
the registry too -- a new cell BC is covered automatically.
"""

from __future__ import annotations

import numpy as np

from myflopy.project.components import (
    LIST_BC_ARTIFACT_TYPES,
    SUPPORTED_PACKAGE_ARTIFACT_TYPES,
    apply_package_artifact,
    build_package_artifact,
    package_artifact_apply_order,
)

LIST_BCS = sorted(LIST_BC_ARTIFACT_TYPES)


def _flatten_text(value) -> set[str]:
    """Lowercased leaf strings from FloPy's nested list/tuple/recarray shapes."""

    if value is None:
        return set()
    if isinstance(value, (str, bytes)):
        return {str(value).lower()}
    if isinstance(value, np.void):  # a numpy record row (np.record subclasses this)
        return {leaf for item in tuple(value) for leaf in _flatten_text(item)}
    if isinstance(value, (list, tuple, np.ndarray)):
        return {leaf for item in value for leaf in _flatten_text(item)}
    return {str(value).lower()}


def test_artifact_support_covers_every_registry_cell_bc():
    from myflopy.modflow.mf6.package_registry import _PACKAGE_EXPLORER_SPECS

    registry_bcs = {
        name for name, spec in _PACKAGE_EXPLORER_SPECS.items() if spec.kind == "cell_stress"
    }
    assert registry_bcs == LIST_BC_ARTIFACT_TYPES
    assert registry_bcs <= SUPPORTED_PACKAGE_ARTIFACT_TYPES
    # every supported type needs a deterministic apply slot, or ordering silently
    # falls back to 999 and dependency-sensitive packages can be applied too early
    for package_type in SUPPORTED_PACKAGE_ARTIFACT_TYPES:
        assert package_artifact_apply_order(package_type) < 999, package_type


def test_every_list_bc_artifact_round_trips(canonical_model, canonical_run_fresh):
    """Capture each list BC from the shared model, re-apply to a private copy.

    One test rather than seven parametrized ones: ``canonical_run_fresh`` builds
    and runs its own model, so looping here costs one build instead of seven.
    """

    target = canonical_run_fresh
    for package in LIST_BCS:
        artifact = build_package_artifact(
            canonical_model,
            artifact_id=f"{package}_artifact",
            package_name=package,
            description=f"round-trip check for {package}",
        )
        assert artifact.package_type == package, package
        assert artifact.package_data["stress_period_data"], package

        target.gwf.remove_package(package)
        apply_package_artifact(target, artifact, validate=True)

        restored = target.gwf.get_package(package)
        assert restored is not None, f"{package} was not re-attached"

        original = canonical_model.gwf.get_package(package).stress_period_data.get_data(key=0)
        recovered = restored.stress_period_data.get_data(key=0)
        assert len(recovered) == len(original), package

        # cellids and the leading numeric field must survive intact
        for before, after in zip(original[:20], recovered[:20], strict=False):
            assert tuple(before[0]) == tuple(after[0]), package
            assert np.isclose(float(before[1]), float(after[1])), package


def test_captured_auxiliary_survives_the_round_trip(canonical_run_fresh):
    """chd/drn used to capture ``auxiliary`` and then drop it on restore."""

    import myflopy as mf

    target = canonical_run_fresh
    target.gwf.remove_package("chd")
    mf.chd.flopy(
        stress_period_data={0: [[(0, 0), 10.0, 1.5]]}, auxiliary=["aux1"]
    ).build(target.gwf)

    artifact = build_package_artifact(target, artifact_id="chd_aux", package_name="chd")
    # FloPy represents this as [["auxiliary", "aux1"]] / [("auxiliary", "aux1")]
    # (keyword + value), so assert on the flattened text rather than a shape
    assert "aux1" in _flatten_text(artifact.package_data.get("auxiliary"))

    target.gwf.remove_package("chd")
    apply_package_artifact(target, artifact, validate=True)

    aux = target.gwf.get_package("chd").auxiliary.get_data()
    assert "aux1" in _flatten_text(aux), f"auxiliary lost on restore: {aux}"


def test_evt_nseg_survives_the_round_trip(canonical_model):
    """``nseg`` is part of EVT's record shape, so it must round-trip."""

    artifact = build_package_artifact(
        canonical_model, artifact_id="evt_nseg", package_name="evt"
    )
    assert artifact.package_data.get("nseg") == 1


def test_mover_survives_the_round_trip_for_every_mover_capable_list_bc(
    canonical_run_fresh,
):
    """A restored list BC must keep MOVER, or any mvr record naming it breaks.

    Found by the 4.7.3 adversarial review (2026-07-18). Capture, restore AND
    dependency inference each hardcoded the mover set to ``{uzf, lak, sfr}``,
    so ``mf.drn(spd, mover=True)`` captured, restored, and came back with MOVER
    off — silently. Nothing caught it: this file asserted only
    stress_period_data/auxiliary/boundnames, and ``_artifact_dependencies``
    reported no ``mvr`` dependency, so ordering validation could not notice
    either.

    The set is now the registry's ``capabilities.mover``, which is itself
    pinned against the live FloPy constructor, so this covers whichever list
    BCs FloPy accepts ``mover`` on rather than a list someone remembered.
    """

    import myflopy as mf
    from myflopy.modflow.mf6.package_registry import _PACKAGE_EXPLORER_SPECS
    from myflopy.project.components import _artifact_dependencies

    mover_capable = sorted(
        name
        for name in LIST_BC_ARTIFACT_TYPES
        if _PACKAGE_EXPLORER_SPECS[name].capabilities.mover
    )
    assert mover_capable, "no mover-capable list BC -- the registry regressed"

    target = canonical_run_fresh
    for package_type in mover_capable:
        record = {
            "drn": [[(0, 0), 10.0, 1.5]],
            "ghb": [[(0, 0), 10.0, 1.5]],
            "riv": [[(0, 0), 10.0, 1.5, 8.0]],
            "wel": [[(0, 0), -1.0]],
        }[package_type]

        target.gwf.remove_package(package_type)
        getattr(mf, package_type).flopy(
            stress_period_data={0: record}, mover=True
        ).build(target.gwf)

        artifact = build_package_artifact(
            target, artifact_id=f"{package_type}_mover", package_name=package_type
        )
        assert artifact.package_data.get("mover") is True, (
            f"{package_type}: mover not captured"
        )
        assert "mvr" in _artifact_dependencies(artifact), (
            f"{package_type}: a mover-bearing artifact must depend on mvr"
        )

        target.gwf.remove_package(package_type)
        apply_package_artifact(target, artifact, validate=True)

        restored = target.gwf.get_package(package_type).mover.get_data()
        assert restored, f"{package_type}: mover lost on restore (got {restored!r})"


def test_artifacts_without_mover_stay_clean(canonical_model):
    """Packages that never had MOVER must not gain it, or acquire a false dep."""

    from myflopy.project.components import _artifact_dependencies

    artifact = build_package_artifact(
        canonical_model, artifact_id="drn_plain", package_name="drn"
    )
    assert artifact.package_data.get("mover") is False
    assert "mvr" not in _artifact_dependencies(artifact)


def test_non_mover_capable_packages_store_no_mover_key(canonical_model):
    """chd/rch/evt have no MF6 mover; storing the key would imply otherwise."""

    from myflopy.modflow.mf6.package_registry import _PACKAGE_EXPLORER_SPECS

    for package_type in sorted(LIST_BC_ARTIFACT_TYPES):
        if _PACKAGE_EXPLORER_SPECS[package_type].capabilities.mover:
            continue
        artifact = build_package_artifact(
            canonical_model,
            artifact_id=f"{package_type}_nomover",
            package_name=package_type,
        )
        assert "mover" not in artifact.package_data, package_type

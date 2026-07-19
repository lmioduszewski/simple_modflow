"""The package descriptor (plan 4.7.2) must equal the code it will replace.

4.7.2 moved per-package knowledge that was restated at ~92 sites into
``package_registry.py``. Nothing consumes it yet — 4.7.3 deletes the hardcoded
lists one at a time. The whole value of that sequencing depends on the
descriptor being *provably* equal to today's behaviour first, so:

**Every assertion here reads the ORIGINAL source, never a literal.** Capability
flags are checked against ``inspect.signature`` of the live FloPy class;
``gpkg_defaults`` against the live ``GeoPackageSource`` signatures; tier flags
against the actual tuples in ``model_diff`` / ``model_results_diff`` /
``components``; node basing against ``budget_tables``. A test that merely
restated the registry's own values would pass while the descriptor was wrong,
which is exactly the failure mode this file exists to prevent.

When 4.7.3 deletes a hardcoded list, its assertion here becomes circular and
must be deleted in the same commit — otherwise it silently stops testing
anything. Each such test is marked ``RETIRE WITH 4.7.3``.
"""

from __future__ import annotations

import inspect

import pytest

from myflopy.modflow.mf6.package_registry import (
    _PACKAGE_EXPLORER_SPECS,
    PackageExplorerSpec,
)

ALL_PACKAGES = sorted(_PACKAGE_EXPLORER_SPECS)
LIST_BCS = ("chd", "drn", "ghb", "riv", "wel", "rch", "evt")
ADVANCED = ("uzf", "sfr", "lak")


def spec(name: str) -> PackageExplorerSpec:
    return _PACKAGE_EXPLORER_SPECS[name]


# ---------------------------------------------------------------------------
# completeness — a descriptor with holes is worse than none
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("package", ALL_PACKAGES)
def test_every_package_carries_the_whole_descriptor(package):
    """No package may be half-described once 4.7.3 starts reading it."""

    entry = spec(package)
    assert entry.flopy_class, package
    assert entry.file_suffix, package
    assert entry.blurb.strip(), package
    assert entry.blurb.strip().endswith("."), f"{package}: blurb is not a sentence"
    assert entry.capabilities is not None
    assert entry.tiers is not None


@pytest.mark.parametrize("package", LIST_BCS)
def test_list_bcs_declare_their_record_and_gpkg_shape(package):
    """The list BCs are what 4.7.4 collapses; they need the full shape."""

    entry = spec(package)
    assert entry.record_fields, package
    assert entry.gpkg_defaults, package


@pytest.mark.parametrize("package", ADVANCED)
def test_advanced_packages_declare_no_flat_record(package):
    """UZF/SFR/LAK take packagedata, not a flat stress-period record.

    Recording this as empty (rather than omitting it) is what stops 4.7.4 from
    trying to generate a list-BC record path for them.
    """

    assert spec(package).record_fields == ()
    assert spec(package).gpkg_defaults == {}


# ---------------------------------------------------------------------------
# capabilities — checked against LIVE FloPy, not against literals
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("package", ALL_PACKAGES)
def test_flopy_class_name_resolves(package):
    import flopy

    entry = spec(package)
    assert hasattr(flopy.mf6, entry.flopy_class), (
        f"{package}: flopy.mf6 has no {entry.flopy_class}"
    )


@pytest.mark.parametrize("package", ALL_PACKAGES)
def test_capabilities_match_the_live_flopy_signature(package):
    """A FloPy upgrade that changes a package's options must fail here."""

    import flopy

    entry = spec(package)
    parameters = set(
        inspect.signature(getattr(flopy.mf6, entry.flopy_class).__init__).parameters
    )
    capabilities = entry.capabilities

    assert capabilities.mover == ("mover" in parameters), package
    assert capabilities.auxiliary == ("auxiliary" in parameters), package
    assert capabilities.boundnames == ("boundnames" in parameters), package
    assert capabilities.observations == ("observations" in parameters), package


# RETIRED IN 4.7.3: ``test_the_mover_list_equals_the_mover_capability``.
# ``MvrResultDiff._MOVER_PACKAGES`` is now DERIVED from
# ``capabilities.mover``, so asserting the two agree would compare the registry
# with itself and pass unconditionally. Coverage did not disappear, it moved
# up the chain: ``test_capabilities_match_the_live_flopy_signature`` checks
# ``capabilities.mover`` against the live FloPy constructor, and the derived
# tuple inherits that guarantee. Ledger entry 14 (``rch`` wrongly listed) is
# closed and the error class is now unrepresentable.


def _record_fields_from_dfn(package: str) -> tuple[str, ...]:
    """Read a package's stress-period record order out of FloPy's own dfn.

    The dfn lists ``cellid`` then the record fields in order, terminated by the
    generic ``aux`` / ``boundname`` entries. This is the authoritative order --
    deriving it here is what stops ``record_fields`` from being the one
    descriptor value that only agrees with itself.
    """

    import flopy

    names = [
        line.split()[-1]
        for block in getattr(flopy.mf6, f"ModflowGwf{package}").dfn
        if isinstance(block, list)
        for line in block
        if isinstance(line, str) and line.startswith("name ")
    ]
    after_cellid = names[names.index("cellid") + 1 :]
    fields: list[str] = []
    for name in after_cellid:
        if name in {"aux", "boundname"}:
            break
        fields.append(name)
    return tuple(fields)


@pytest.mark.parametrize("package", [p for p in LIST_BCS if p != "evt"])
def test_record_fields_match_the_flopy_dfn(package):
    """The record order must be MF6's, not one we remembered."""

    assert spec(package).record_fields == _record_fields_from_dfn(package)


def test_evt_records_the_single_segment_shape():
    """EVT is the one list BC whose record length depends on an option.

    The dfn carries the fully segmented form (``pxdp``/``petm``/``petm0``).
    myflopy's ``.gpkg`` path supports ``nseg=1`` only and rejects anything else
    (``geopackage.py``), so the descriptor records that prefix. Asserting it IS
    a prefix — rather than hardcoding three names — keeps this honest if MF6
    ever reorders the segmented tail.
    """

    full = _record_fields_from_dfn("evt")
    declared = spec("evt").record_fields

    assert declared == full[: len(declared)], f"{declared} is not a prefix of {full}"
    assert declared == ("surface", "rate", "depth")
    assert len(full) > len(declared), "expected a segmented tail beyond nseg=1"


@pytest.mark.parametrize("package", LIST_BCS)
def test_gpkg_defaults_match_the_live_resolver_signature(package):
    """Parameter names AND their default column names, read off the method."""

    from myflopy.geopackage import GeoPackageSource

    resolver = getattr(GeoPackageSource, package)
    parameters = inspect.signature(resolver).parameters
    declared = spec(package).gpkg_defaults

    for name, default_column in declared.items():
        assert name in parameters, f"{package}.gpkg has no parameter {name!r}"
        assert parameters[name].default == default_column, (
            f"{package}.{name} defaults to {parameters[name].default!r}, "
            f"descriptor says {default_column!r}"
        )


@pytest.mark.parametrize("package", LIST_BCS)
def test_edges_only_capability_matches_the_resolver(package):
    """``edges_only`` is a myflopy capability — only some resolvers expose it."""

    from myflopy.geopackage import GeoPackageSource

    parameters = inspect.signature(getattr(GeoPackageSource, package)).parameters
    assert spec(package).capabilities.edges_only == ("edges_only" in parameters), package


# ---------------------------------------------------------------------------
# the tier flags — checked against the real lists (RETIRE WITH 4.7.3)
# ---------------------------------------------------------------------------


# RETIRED IN 4.7.3: ``test_diffable_flags_match_model_diff``.
# ``_DIFF_PACKAGES`` and ``_CONNECTION_PACKAGES`` are now derived from
# ``tiers.diffable`` / ``tiers.connection_diffable``. The real invariant --
# that the diff tier matches the ModelGroup's still-hardcoded
# ``GroupPackageInputs`` accessors -- lives in
# ``tests/test_model_group_symmetry.py`` and is now a genuine cross-check
# between the registry and the group, rather than between two hand-typed lists.


def test_diffable_and_connection_tiers_are_disjoint():
    """A package is diffed row-by-row OR by connection geometry, never both.

    Not circular: it constrains the descriptor's own shape. The two diff paths
    build different frames, so a package flagged for both would silently take
    whichever branch ran first.
    """

    for package in ALL_PACKAGES:
        tiers = spec(package).tiers
        assert not (tiers.diffable and tiers.connection_diffable), package


# RETIRED IN 4.7.3: ``test_results_diffable_flags_match_model_results_diff``.
# ``_CELL_BUDGET_PACKAGES`` is now derived from ``tiers.results_diffable``.
# Unlike ``mover`` there is no independent source to check the flag against --
# which packages are worth results-diffing is a judgment, and the descriptor is
# now where that judgment lives. Behavioural coverage is in
# ``tests/test_model_results_diff.py``, which exercises the diff end to end.


def test_the_results_diff_namespace_still_covers_every_flagged_package():
    """Guards the derivation itself: every flagged package must be reachable.

    Not circular -- it checks the *consumer* can actually build a namespace for
    each package the registry claims, which is the thing that would break if a
    descriptor flag and the subsystem's real capability diverged.
    """

    from myflopy.project.model_results_diff import _CELL_BUDGET_PACKAGES

    flagged = {name for name in ALL_PACKAGES if spec(name).tiers.results_diffable}
    assert flagged, "no package is flagged results_diffable"
    assert set(_CELL_BUDGET_PACKAGES) == flagged
    # uzf is deliberately absent -- it has no per-cell q term (ledger 4.7 gaps)
    assert "uzf" not in flagged


def test_artifact_flags_and_order_match_components():
    """RETIRE WITH 4.7.3."""

    from myflopy.project.components import (
        PACKAGE_ARTIFACT_APPLY_ORDER,
        SUPPORTED_PACKAGE_ARTIFACT_TYPES,
    )

    for package in ALL_PACKAGES:
        tiers = spec(package).tiers
        assert tiers.artifact_serializable == (
            package in SUPPORTED_PACKAGE_ARTIFACT_TYPES
        ), package
        assert tiers.artifact_apply_order == PACKAGE_ARTIFACT_APPLY_ORDER.get(
            package
        ), package


def test_model_accessor_flag_matches_simulation_base():
    """RETIRE WITH 4.7.3.

    Several packages deliberately have no ``model.<pkg>`` accessor today
    (CLAUDE.md notes chd/drn/ghb/wel are missing while rch/uzf/sfr/lak exist).
    The descriptor records the gap rather than papering over it.
    """

    from myflopy.modflow.mf6.simulation.base import SimulationBase

    for package in ALL_PACKAGES:
        actual = isinstance(getattr(SimulationBase, package, None), property)
        assert spec(package).tiers.model_accessor == actual, package


# ---------------------------------------------------------------------------
# suffix + budget node basing (RETIRE WITH 4.7.3)
# ---------------------------------------------------------------------------


def test_file_suffix_matches_run_model():
    """RETIRE WITH 4.7.3."""

    from myflopy.project.run_model import _PACKAGE_SUFFIX_TO_TYPE

    for package in ALL_PACKAGES:
        suffix = spec(package).file_suffix
        assert suffix in _PACKAGE_SUFFIX_TO_TYPE, f"{package}: {suffix} unmapped"
        assert _PACKAGE_SUFFIX_TO_TYPE[suffix].lower() == package, package


def test_zero_base_budget_nodes_matches_budget_tables():
    """RETIRE WITH 4.7.3.

    This is the flag behind the budget off-by-one fixed in 4.7.1 — MF6 reports
    1-based node numbers for the cell-stress list BCs. Getting it wrong was a
    live bug for years, so it is pinned against the deriving function itself.
    """

    from myflopy.modflow.mf6.budget_tables import _cell_based_budget_packages

    cell_based = set(_cell_based_budget_packages())
    for package in ALL_PACKAGES:
        assert spec(package).zero_base_budget_nodes == (package in cell_based), package


def test_node_basing_splits_cell_stress_from_advanced():
    """The flag must track ``kind``, which is the physical reason behind it."""

    for package in ALL_PACKAGES:
        entry = spec(package)
        assert entry.zero_base_budget_nodes == (entry.kind == "cell_stress"), package


# ---------------------------------------------------------------------------
# blurbs are prose, not generated
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("package", ALL_PACKAGES)
def test_blurbs_describe_hydrology_not_implementation(package):
    """A blurb naming a field or a FloPy class has drifted into being a spec."""

    entry = spec(package)
    blurb = entry.blurb
    assert 40 < len(blurb) < 220, f"{package}: blurb length {len(blurb)}"
    assert entry.flopy_class not in blurb, package
    for field_name in entry.record_fields:
        assert f"``{field_name}``" not in blurb, package


def test_blurbs_are_distinct():
    """A copy-pasted blurb is worse than none — it misdescribes a package."""

    blurbs = [spec(name).blurb for name in ALL_PACKAGES]
    assert len(set(blurbs)) == len(blurbs)

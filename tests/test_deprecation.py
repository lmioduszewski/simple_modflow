"""The single deprecation mechanism (implementation plan 3.1/3.2, D12).

Contracts pinned here:

1. Every deprecated name works with a ``DeprecationWarning`` naming its
   replacement, and resolves to the same object the replacement names.
2. D12 hiding — deprecated names are absent from ``__all__`` and ``dir()``
   at every facade (module- and class-level).
3. ``myflopy.__compatibility__`` is the authoritative registry: it matches
   exactly what the helpers registered once all deprecating modules import.
4. The helper primitives behave (warn text, fallback chaining, unknown-name
   AttributeError, module-alias targets).
"""

from __future__ import annotations

import warnings

import pytest

import myflopy
from myflopy._deprecation import (
    deprecated_instance_getattr,
    deprecated_module_getattr,
    registered_deprecations,
    warn_deprecated,
)

# ---------------------------------------------------------------------------
# helper primitives
# ---------------------------------------------------------------------------


@pytest.fixture
def registry_snapshot():
    """Keep test-only registrations out of the process-global registry."""

    from myflopy import _deprecation

    before = dict(_deprecation._REGISTRY)
    yield
    _deprecation._REGISTRY.clear()
    _deprecation._REGISTRY.update(before)


def test_warn_deprecated_message_names_old_new_and_since():
    with pytest.warns(DeprecationWarning) as record:
        warn_deprecated("myflopy.test.old_name", "new_name", since="0.9.0", stacklevel=2)
    message = str(record[0].message)
    assert message.startswith("myflopy.test.old_name is deprecated since myflopy 0.9.0")
    assert "use new_name instead" in message
    assert "deprecation_policy" in message


def test_module_getattr_resolves_warns_and_hides(registry_snapshot):
    mapping = {"OldName": ("myflopy._deprecation:warn_deprecated", "warn_deprecated", "0.9.0")}
    getattr_fn, dir_fn = deprecated_module_getattr(mapping, "myflopy._deprecation")

    with pytest.warns(DeprecationWarning, match="OldName is deprecated"):
        assert getattr_fn("OldName") is warn_deprecated
    assert "OldName" not in dir_fn()
    assert "warn_deprecated" in dir_fn()  # real module contents still advertised

    with pytest.raises(AttributeError, match="myflopy._deprecation"):
        getattr_fn("nonsense")


def test_module_getattr_supports_module_targets_and_fallbacks(registry_snapshot):
    mapping = {"old_mod": ("myflopy._deprecation", "myflopy._deprecation", "0.9.0")}
    fallback_hits = []

    def fallback(name):
        fallback_hits.append(name)
        return "fallback-value"

    getattr_fn, dir_fn = deprecated_module_getattr(
        mapping,
        "myflopy._deprecation",
        fallback_getattr=fallback,
        fallback_dir=lambda: ["from_fallback", "old_mod"],
    )

    import myflopy._deprecation as dep_module

    with pytest.warns(DeprecationWarning):
        assert getattr_fn("old_mod") is dep_module
    assert getattr_fn("anything_else") == "fallback-value"
    assert fallback_hits == ["anything_else"]
    # fallback_dir feeds __dir__, deprecated names are still subtracted
    assert dir_fn() == ["from_fallback"]


def test_instance_getattr_warns_hides_and_raises_for_unknown(registry_snapshot):
    class Holder:
        def __init__(self):
            self._value = 42

        __getattr__ = deprecated_instance_getattr(
            {"value": ("_value", "holder.new_value", "0.9.0")},
            "myflopy.tests.Holder",
        )

    holder = Holder()
    with pytest.warns(DeprecationWarning, match="myflopy.tests.Holder.value is deprecated"):
        assert holder.value == 42
    assert "value" not in dir(holder)
    with pytest.raises(AttributeError, match="no attribute 'missing'"):
        holder.missing


# ---------------------------------------------------------------------------
# GHB/DRN vector-builder rename (3.2)
# ---------------------------------------------------------------------------


def test_ghb_module_alias_warns_resolves_and_hides():
    import myflopy.modflow.mf6.ghb as ghb_module

    with pytest.warns(DeprecationWarning, match="use GHBFromVector instead"):
        assert ghb_module.GHB is ghb_module.GHBFromVector
    assert ghb_module.GHBFromVector.__name__ == "GHBFromVector"
    assert "GHB" not in dir(ghb_module)


def test_drn_module_alias_warns_resolves_and_hides():
    import myflopy.modflow.mf6.drn as drn_module

    with pytest.warns(DeprecationWarning, match="use DRNFromVector instead"):
        assert drn_module.DRN is drn_module.DRNFromVector
    assert drn_module.DRNFromVector.__name__ == "DRNFromVector"
    assert "DRN" not in dir(drn_module)


def test_mf6_facade_aliases_warn_resolve_and_hide():
    import myflopy.modflow.mf6 as mf6
    from myflopy.modflow.mf6.drn import DRNFromVector
    from myflopy.modflow.mf6.ghb import GHBFromVector

    with pytest.warns(DeprecationWarning):
        assert mf6.GHB is GHBFromVector
    with pytest.warns(DeprecationWarning):
        assert mf6.DRN is DRNFromVector

    for old, new in (("GHB", "GHBFromVector"), ("DRN", "DRNFromVector")):
        assert old not in mf6.__all__
        assert old not in dir(mf6)
        assert new in mf6.__all__
        assert new in dir(mf6)


def test_simulation_layer_ghb_is_a_different_undeprecated_class():
    # The bare GHB at the simulation layer is the OO package wrapper, not the
    # vector builder — it keeps its name and must not warn.
    from myflopy.modflow.mf6.ghb import GHBFromVector
    from myflopy.modflow.mf6.simulation.packages import GHB as SimulationGHB

    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        from myflopy.modflow.mf6.simulation import GHB as ReExportedGHB
    assert ReExportedGHB is SimulationGHB
    assert SimulationGHB is not GHBFromVector


def test_legacy_utils_surfaces_path_warns_resolves_and_hides():
    # Plan 4.3: InterpolatedSurface moved to mf6/grid/interpolated_surface.
    import myflopy.modflow.utils.surfaces as legacy_surfaces
    from myflopy.modflow.mf6.grid.interpolated_surface import InterpolatedSurface

    with pytest.warns(DeprecationWarning, match="interpolated_surface"):
        assert legacy_surfaces.InterpolatedSurface is InterpolatedSurface
    assert "InterpolatedSurface" not in dir(legacy_surfaces)


# ---------------------------------------------------------------------------
# legacy PRT aliases (mp3du) — absorbed ad-hoc __getattr__
# ---------------------------------------------------------------------------


def test_legacy_prt_names_warn_resolve_and_hide():
    import myflopy.modflow.mp3du.legacy_prt as legacy_prt
    import myflopy.modflow.mp3du.particles as particles

    with pytest.warns(DeprecationWarning, match="legacy_prt.PRT instead"):
        assert particles.PRT is legacy_prt.PRT

    for name in ("PRT", "PrtMip", "PrtOc", "PrtPrp", "PrtDisv", "PrtFmi"):
        assert name not in dir(particles)
        assert name not in particles.__all__


# ---------------------------------------------------------------------------
# the authoritative registry
# ---------------------------------------------------------------------------


def test_compatibility_registry_matches_the_helpers_exactly():
    # Import every module that registers deprecated names, then compare.
    import myflopy.modflow.mf6  # noqa: F401
    import myflopy.modflow.mf6.drn  # noqa: F401
    import myflopy.modflow.mf6.ghb  # noqa: F401
    import myflopy.modflow.mp3du.particles  # noqa: F401
    import myflopy.modflow.utils.surfaces  # noqa: F401
    import myflopy.project.model_group  # noqa: F401

    registered = set(registered_deprecations())
    declared = set(myflopy.__compatibility__)
    assert registered == declared, (
        "myflopy.__compatibility__ is out of sync with the deprecation "
        f"helpers; only-registered={sorted(registered - declared)} "
        f"only-declared={sorted(declared - registered)}"
    )
    # Every entry documents a replacement and a since-version.
    for old, (replacement, since) in registered_deprecations().items():
        assert old.startswith("myflopy"), old
        assert replacement, old
        assert since.count(".") == 2, (old, since)


def test_engine_tier_is_unwarned_and_importable():
    assert myflopy.__engine__  # the second-tier marker survives the repurpose
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        for name in myflopy.__engine__:
            assert getattr(myflopy, name) is not None
            assert name not in myflopy.__all__

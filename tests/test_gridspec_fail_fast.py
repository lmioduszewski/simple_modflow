"""GridSpec.structured / from_geopackage must fail fast, not exist-then-explode.

myflopy resolves only generated Voronoi (``GridSpec.voronoi``) and Python
(``GridSpec.python``) grid specs. The structured and existing/GeoPackage
constructors are declared but unwired; rather than let a caller build a spec
that only explodes later at ``resolve()``, they raise immediately with a pointer
to the supported alternatives.
"""

from __future__ import annotations

import geopandas as gpd
import pytest
import shapely as shp

from myflopy.sources import ShapeSource
from myflopy.specs import _VORONOI_OPTION_KEYS, GridSpec


def test_structured_constructor_fails_fast_with_pointer():
    with pytest.raises(NotImplementedError) as excinfo:
        GridSpec.structured(
            nlay=1, nrow=2, ncol=2, delr=10.0, delc=10.0, top=1.0, botm=0.0
        )
    message = str(excinfo.value)
    assert "not implemented yet" in message
    assert "GridSpec.voronoi" in message
    assert "GridSpec.from_object" in message


def test_from_geopackage_constructor_fails_fast_with_pointer():
    with pytest.raises(NotImplementedError) as excinfo:
        GridSpec.from_geopackage("some_grid.gpkg", layer="cells")
    message = str(excinfo.value)
    assert "not implemented yet" in message
    assert "GridSpec.from_object" in message


def test_resolve_rejects_non_voronoi_method_reconstructed_via_from_dict():
    """A serialized structured spec bypasses the fail-fast constructor;
    ``resolve()`` must still refuse it with an actionable message."""

    payload = {
        "kind": "GridSpec",
        "name": "g",
        "grid_type": "dis",
        "method": "structured",
        "options": {},
    }
    spec = GridSpec.from_dict(payload)
    with pytest.raises(NotImplementedError) as excinfo:
        spec.resolve()
    assert "GridSpec.from_object" in str(excinfo.value)


# --- GridSpec.voronoi rejects what it used to swallow ----------------------- #
#
# `**engine_options` accepted ANY keyword and stored it unread, so a misspelled
# option silently meshed at the default -- measured, `breakline_bufer=40` gave
# corridors a quarter the intended width under a clean run. Naming the options
# (so an editor completes them) and rejecting the rest are the two halves.

def _boundary(tmp_path):
    path = tmp_path / "domain.gpkg"
    gpd.GeoDataFrame(
        {"n": [1]}, geometry=[shp.box(0, 0, 100, 100)], crs="EPSG:2927"
    ).to_file(path)
    return ShapeSource(path)


def _voronoi(tmp_path, **kwargs):
    return GridSpec.voronoi(
        boundary=_boundary(tmp_path), crs="EPSG:2927",
        boundary_max_area=1000.0, **kwargs
    )


def test_the_sizing_options_are_named_parameters_not_kwargs():
    """PyCharm and Pylance read the `def` line and never run the module, so a
    forwarded option only exists for a user if it is in the signature."""

    import inspect

    named = set(inspect.signature(GridSpec.voronoi).parameters)
    for option in (
        "boundary_max_area", "refinement_max_area", "refinement_priority",
        "refinement_buffer", "breakline_max_area", "breakline_priority",
        "breakline_buffer", "region_point_tolerance",
    ):
        assert option in named, f"{option} is reachable but not discoverable"


def test_a_misspelled_option_raises_and_suggests_the_real_one(tmp_path):
    with pytest.raises(TypeError) as excinfo:
        _voronoi(tmp_path, breakline_bufer=40)
    message = str(excinfo.value)
    assert "breakline_bufer" in message
    assert "Did you mean 'breakline_buffer'?" in message


def test_an_unrecognisable_option_still_raises(tmp_path):
    """No near match, so no suggestion -- but silence is not an option."""

    with pytest.raises(TypeError, match="unknown option 'compleletly_made_up'"):
        _voronoi(tmp_path, compleletly_made_up=123)


def test_legacy_option_aliases_still_resolve(tmp_path):
    """`line_buffer`/`max_area`/`default_*_area` are read by the resolver, so the
    allowlist has to admit them even though they are not named parameters."""

    spec = _voronoi(tmp_path, line_buffer=40, default_refinement_area=100.0)
    assert spec.options["line_buffer"] == 40
    assert spec.options["default_refinement_area"] == 100.0


def test_named_options_are_absent_when_not_passed(tmp_path):
    """`None` must mean "not given", or a present-but-None `boundary_max_area`
    would shadow its own `default_cell_area` alias in `_option`."""

    spec = _voronoi(tmp_path)
    for key in ("refinement_buffer", "breakline_priority", "boundary_label"):
        assert key not in spec.options


def test_the_option_allowlist_matches_what_the_resolver_reads():
    """Ratchets BOTH ways: a new `_option(options, ...)` read in the resolver
    that is not admitted here would be rejected at the call, and an allowlist
    entry nothing reads is dead weight that silently accepts a typo."""

    import re
    from pathlib import Path as _Path

    import myflopy.grid_spec_resolver as resolver

    source = _Path(resolver.__file__).read_text()
    read: set[str] = set()
    for match in re.finditer(
        r'_option\(\s*(?:build_)?options,\s*((?:"[a-z_]+"\s*,?\s*)+)', source
    ):
        read |= set(re.findall(r'"([a-z_]+)"', match.group(1)))
    read |= set(re.findall(r'options\.get\(\s*"([a-z_]+)"', source))

    missing = read - _VORONOI_OPTION_KEYS
    assert not missing, f"resolver reads options the allowlist rejects: {sorted(missing)}"
    extra = _VORONOI_OPTION_KEYS - read
    assert not extra, f"allowlist admits options nothing reads: {sorted(extra)}"


@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({"refinement": ["a", "b"]}, "takes ONE source"),
        ({"breaklines": ShapeSource("x.gpkg")}, "takes a SEQUENCE"),
        ({"points": ShapeSource("x.gpkg")}, "takes a SEQUENCE"),
    ],
)
def test_grid_sources_reject_the_wrong_shape(tmp_path, kwargs, match):
    """Each of these used to reach the resolver and die as `AttributeError:
    'X' object has no attribute 'path'`, naming neither argument nor fix."""

    with pytest.raises(TypeError, match=match):
        _voronoi(tmp_path, **kwargs)


def test_an_in_memory_geodataframe_is_refused_with_the_way_out(tmp_path):
    """A GridSpec is a serializable recipe, so a live frame cannot be a source --
    but the error has to say so and name the one-line fix."""

    frame = gpd.GeoDataFrame(
        {"n": [1]}, geometry=[shp.box(0, 0, 10, 10)], crs="EPSG:2927"
    )
    with pytest.raises(TypeError) as excinfo:
        _voronoi(tmp_path, refinement=frame)
    message = str(excinfo.value)
    assert "GeoDataFrame" in message
    assert "to_file" in message and "GeoPackageSourceSpec" in message

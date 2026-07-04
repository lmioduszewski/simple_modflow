"""GridSpec.structured / from_geopackage must fail fast, not exist-then-explode.

myflopy resolves only generated Voronoi (``GridSpec.voronoi``) and Python
(``GridSpec.python``) grid specs. The structured and existing/GeoPackage
constructors are declared but unwired; rather than let a caller build a spec
that only explodes later at ``resolve()``, they raise immediately with a pointer
to the supported alternatives.
"""

from __future__ import annotations

import pytest

from myflopy.specs import GridSpec


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

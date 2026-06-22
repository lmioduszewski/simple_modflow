from __future__ import annotations

import warnings

import pytest

from myflopy.modflow.utils.derived_raster import DerivedRaster
from myflopy.surfaces import Surface


def test_derived_raster_lifecycle(tmp_path):
    src = tmp_path / "src.txt"
    src.write_text("v1")
    out = tmp_path / "out.tif"
    calls = {"n": 0}

    def produce():
        calls["n"] += 1
        out.write_text("raster")

    dr = DerivedRaster(out, [src], {"resolution": 4}, produce)
    assert dr.status() == "missing"

    dr.ensure()  # produce on first use
    assert calls["n"] == 1
    assert out.exists() and dr.sidecar.exists()
    assert dr.status() == "fresh"

    dr.ensure()  # fresh -> no rebuild
    assert calls["n"] == 1

    src.write_text("v2")  # source changed -> stale
    assert dr.status() == "stale"
    with pytest.warns(UserWarning, match="stale"):
        dr.ensure()
    assert calls["n"] == 1  # stale is reused, NOT rebuilt

    dr.ensure(refresh=True)  # explicit rebuild
    assert calls["n"] == 2
    assert dr.status() == "fresh"


def test_derived_raster_params_change_is_stale(tmp_path):
    src = tmp_path / "s.txt"
    src.write_text("x")
    out = tmp_path / "o.tif"
    produce = lambda: out.write_text("r")  # noqa: E731

    DerivedRaster(out, [src], {"resolution": 4}, produce).ensure()
    # Same source, different params -> stale.
    assert DerivedRaster(out, [src], {"resolution": 8}, produce).status() == "stale"


def test_derived_raster_warn_false_is_silent(tmp_path):
    src = tmp_path / "s.txt"
    src.write_text("x")
    out = tmp_path / "o.tif"
    calls = {"n": 0}

    def produce():
        calls["n"] += 1
        out.write_text("r")

    dr = DerivedRaster(out, [src], {}, produce)
    dr.ensure()
    src.write_text("y")  # now stale
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # any warning would fail the test
        dr.ensure(warn=False)
    assert calls["n"] == 1  # reused silently


def test_derived_raster_missing_when_output_deleted(tmp_path):
    src = tmp_path / "s.txt"
    src.write_text("x")
    out = tmp_path / "o.tif"
    dr = DerivedRaster(out, [src], {}, lambda: out.write_text("r"))
    dr.ensure()
    assert dr.status() == "fresh"
    out.unlink()
    assert dr.status() == "missing"


def test_contour_surface_is_derived_and_cache_status(tmp_path):
    # No GRASS is invoked: we only inspect derived/cache metadata.
    surface = Surface.from_contours(
        tmp_path / "c.gpkg", z="elev", out=tmp_path / "c.interp.tif"
    )
    assert surface.is_derived is True
    assert surface.cache_status() == "missing"


def test_plain_raster_surface_is_not_derived(tmp_path):
    surface = Surface.raster(tmp_path / "r.tif")
    assert surface.is_derived is False
    assert surface.cache_status() is None

"""A bad ``gaps=``/``min_above_bottom=``/``snap_within=`` on the EAGER
``mf.<pkg>.gpkg(context=...)`` route names ``mf.<pkg>.gpkg``, as the lazy route does.

``GeoPackageSource.__post_init__`` checked them under "A GeoPackage source" before the
helper could tell the source which call built it.
"""

from __future__ import annotations

import geopandas as gpd
import pytest
from shapely.geometry import Point

import myflopy as mf
from tests.test_gpkg_on_cell_faces import _context

_VALUES = {
    "chd": {"head": 9.0},
    "ghb": {"head": 9.0, "cond": 1.0},
    "drn": {"elevation": 9.0, "cond": 1.0},
    "riv": {"stage": 9.0, "cond": 1.0, "rbot": 8.0},
    "wel": {"rate": -1.0},
    "rch": {"recharge": 1e-3},
    "evt": {"surface": 9.0, "rate": 1e-3, "depth": 1.0},
}


@pytest.fixture
def source(tmp_path):
    path = tmp_path / "points.gpkg"
    gpd.GeoDataFrame({"name": ["a"]}, geometry=[Point(12.5, 12.5)], crs="EPSG:2927").to_file(
        path, layer="bc", driver="GPKG"
    )
    return path


#: ``min_above_bottom=`` exists on the four helpers with a below-the-bottom elevation.
_CASES = [
    (package, option, value, message)
    for package in sorted(_VALUES)
    for option, value, message in (
        ("gaps", "nearst", "unknown gaps mode"),
        ("min_above_bottom", -1.0, "min_above_bottom"),
        ("snap_within", 0.0, "snap_within"),
    )
    if option != "min_above_bottom" or package in {"chd", "ghb", "drn", "riv"}
]


@pytest.mark.parametrize(("package", "option", "value", "message"), _CASES)
@pytest.mark.parametrize("eager", [True, False], ids=["eager", "lazy"])
def test_the_refusal_names_the_helper(source, package, option, value, message, eager):
    helper = getattr(mf, package).gpkg
    context = {"context": _context()} if eager else {}
    with pytest.raises(ValueError, match=rf"^mf\.{package}\.gpkg: .*{message}"):
        helper(source, nper=1, layer_field=None, **context, **_VALUES[package], **{option: value})

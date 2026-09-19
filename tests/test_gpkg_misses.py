"""A GeoPackage that lands nowhere on the grid is refused; a feature that misses is named.

A feature that resolved to no active cell was dropped without a word. When EVERY
feature missed -- a source in another location, a file with no CRS, a model placed
at the origin -- the package came back empty: MF6 rejects an empty list package at
read time with a message about the file rather than the cause, and with
``maxbound=`` runs without the boundary at all, under normal termination.
``k_array`` returned its defaults everywhere. Both now raise, naming where the
features are, where the grid is, and both CRSs. A partial miss can be deliberate,
so it warns -- by name.
"""

from __future__ import annotations

import logging
import warnings

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import box

import myflopy as mf
from myflopy.geopackage import GeoPackageSource, Spread
from myflopy.specs import ModelContext
from tests.test_gpkg_on_cell_faces import _VALUES, CELL_5, LIST_BCS, _cells, _context

# Real State Plane coordinates; every test grid here sits at the origin.
FAR = box(1_300_000.0, 650_000.0, 1_300_100.0, 650_100.0)


def _write(tmp_path, *geometries, crs="EPSG:2927", names=None):
    path = tmp_path / "bc.gpkg"
    n = len(geometries)
    gpd.GeoDataFrame(
        {
            "name": names or [f"f{i}" for i in range(n)],
            "layer": [1] * n,
            **{k: [v] * n for k, v in _VALUES.items()},
        },
        geometry=list(geometries),
        crs=crs,
    ).to_file(path, layer="bc", driver="GPKG")
    return path


# --- nothing lands: refused ------------------------------------------------- #


@pytest.mark.parametrize("package", LIST_BCS)
def test_a_source_that_lands_nowhere_is_refused_by_every_list_bc(tmp_path, package):
    with pytest.raises(ValueError, match="lands on an active cell"):
        getattr(mf, package).gpkg(_write(tmp_path, FAR), layer="bc", context=_context(), nper=1)


def test_the_refusal_says_where_the_features_and_the_grid_are(tmp_path):
    with pytest.raises(ValueError) as caught:
        mf.drn.gpkg(_write(tmp_path, FAR), layer="bc", context=_context(), nper=1)
    message = str(caught.value)
    assert "'bc.gpkg' (layer 'bc')" in message
    assert "the package would have no records" in message
    assert "'f0' (outside the grid)" in message
    assert "(1,300,000.0, 650,000.0, 1,300,100.0, 650,100.0)" in message
    assert "the grid (0.0, 0.0, 100.0, 100.0)" in message
    assert "both are in EPSG:2927" in message


def test_a_file_with_no_crs_is_named_as_the_likely_cause(tmp_path):
    lonlat = box(-122.5, 47.0, -122.4, 47.1)
    with pytest.raises(ValueError, match="the file has NO CRS"):
        mf.drn.gpkg(_write(tmp_path, lonlat, crs=None), layer="bc", context=_context(), nper=1)


def test_a_reprojected_source_names_both_crss(tmp_path):
    lonlat = box(-122.5, 47.0, -122.4, 47.1)
    with pytest.raises(ValueError, match="reprojected from EPSG:4326 to the grid's EPSG:2927"):
        mf.drn.gpkg(_write(tmp_path, lonlat, crs="EPSG:4326"), layer="bc", context=_context(), nper=1)


def test_features_on_inactive_cells_only_are_refused_as_such(tmp_path):
    grid = _context().grid
    dead = ModelContext(grid=grid, domain=np.zeros((1, 16), dtype=int))
    with pytest.raises(ValueError, match=r"'f0' \(on no active cell\)"):
        mf.drn.gpkg(_write(tmp_path, CELL_5), layer="bc", context=dead, nper=1)


def test_edges_only_says_it_was_looking_at_the_edge(tmp_path):
    """Cell 5 is interior: with edges_only it has nothing to land on."""

    with pytest.raises(ValueError, match="grid-edge cell"):
        mf.chd.gpkg(_write(tmp_path, CELL_5), layer="bc", context=_context(), nper=1, edges_only=True)


def test_a_full_miss_under_spread_is_refused_not_just_logged(tmp_path):
    with pytest.raises(ValueError, match="lands on an active cell"):
        mf.drn.gpkg(
            _write(tmp_path, FAR), layer="bc", context=_context(), nper=1,
            conductance=Spread("conductance"),
        )


def test_an_empty_layer_is_refused(tmp_path):
    path = _write(tmp_path, CELL_5)
    frame = gpd.read_file(path, layer="bc").iloc[0:0]
    frame.to_file(path, layer="bc", driver="GPKG")
    with pytest.raises(ValueError, match="has no features"):
        mf.drn.gpkg(path, layer="bc", context=_context(), nper=1)


def test_k_array_refuses_zones_that_land_nowhere(tmp_path):
    source = GeoPackageSource(_write(tmp_path, FAR), _context(), 1, layer="bc")
    with pytest.raises(ValueError, match="defaults everywhere"):
        source.k_array(value="k", nlay=1, defaults=1.0)


def test_a_polygon_that_only_borders_the_grid_is_told_to_overlap_it(tmp_path):
    """An ocean polygon on the coastline reached the edge cells only through the
    face-touch bug (ledger 210). It covers none of them, and now says so."""

    ocean = box(-50.0, 0.0, 0.0, 100.0)
    with pytest.raises(ValueError, match="only borders cells"):
        mf.chd.gpkg(
            _write(tmp_path, ocean), layer="bc", context=_context(), nper=1, edges_only=True
        )


# --- some land: warned, by name --------------------------------------------- #


def test_a_feature_that_misses_is_named_through_the_channel_a_caller_HEARS(tmp_path):
    """`caplog` is not evidence that anyone is told.

    `import myflopy` installs a `NullHandler` on the `myflopy` logger and root has
    none, so `logging.callHandlers` counts a handler as found and never falls back to
    `logging.lastResort` -- the partial-miss line produced **zero bytes on stderr** in a
    bare session and in a real ipykernel. These tests passed only because pytest
    attaches the handler the user does not have, which is the "a test can certify
    instead of check" pattern. Asserted on `pytest.warns` now, which pytest does not
    supply for free (ledger 288). The logger still emits, for a configured application.
    """

    path = _write(tmp_path, CELL_5, FAR, names=["inside", "away"])

    with pytest.warns(UserWarning, match="1 of 2 features") as caught:
        spec = mf.drn.gpkg(path, layer="bc", context=_context(), nper=1)

    assert _cells(spec) == [5]
    assert "'away' (outside the grid)" in str(caught[0].message)


def test_the_message_names_the_PACKAGE_so_four_of_them_are_four_warnings(tmp_path):
    """`warnings.warn` dedupes on (text, category, module, lineno).

    One GeoPackage driving four packages emitted four BYTE-IDENTICAL lines, so a
    default caller saw one of four -- in exactly the shape ledger 211 names as the
    legitimate reason to warn rather than refuse, one layer serving several packages.
    The package name is in the message now, which makes each event distinct AND is the
    accounting; one change buys both.
    """

    path = _write(tmp_path, CELL_5, FAR, names=["inside", "away"])
    context = _context()

    with pytest.warns(UserWarning) as caught:
        for helper in (mf.drn, mf.ghb, mf.riv, mf.chd):
            helper.gpkg(path, layer="bc", context=context, nper=1)

    subjects = {str(record.message).split(" -> ")[1].split(":")[0] for record in caught}
    assert subjects == {"drn", "ghb", "riv", "chd"}


def test_a_partial_miss_under_spread_is_reported_once(tmp_path):
    """Spread used to log its own 'covers no active cell' line on top."""

    path = _write(tmp_path, CELL_5, FAR, names=["inside", "away"])

    with pytest.warns(UserWarning) as caught:
        mf.drn.gpkg(
            path, layer="bc", context=_context(), nper=1, conductance=Spread("conductance")
        )

    assert len(caught) == 1
    assert "'away' (outside the grid)" in str(caught[0].message)


def test_k_array_names_a_zone_that_misses_and_keeps_the_rest(tmp_path, caplog):
    path = _write(tmp_path, CELL_5, FAR, names=["clay", "elsewhere"])
    with caplog.at_level(logging.WARNING, logger="myflopy.geopackage"):
        k = GeoPackageSource(path, _context(), 1, layer="bc").k_array(value="k", nlay=1, defaults=1.0)
    assert np.flatnonzero(k[0] != 1.0).tolist() == [5]
    assert "'elsewhere' (outside the grid)" in caplog.text


def test_a_source_that_lands_entirely_says_nothing(tmp_path, caplog):
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # any warning at all fails this
        with caplog.at_level(logging.WARNING, logger="myflopy.geopackage"):
            mf.drn.gpkg(_write(tmp_path, CELL_5), layer="bc", context=_context(), nper=1)
    assert caplog.text == ""

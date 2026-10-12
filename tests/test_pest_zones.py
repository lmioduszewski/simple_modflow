"""The zone front door for ``style="zone"`` (plan §5.8 item 2).

Before this, ``zones=`` was an unvalidated pass-through: whatever the caller
handed in went straight to pyEMU as ``zone_array``. That is a problem because
pyEMU's contract is ASYMMETRIC between the two target families and its
rejections name neither the target nor the shape — so the tests here are mostly
about the normalization, not about the raster reader the plan item advertised.
"""

from __future__ import annotations

import tracemalloc
import warnings

import geopandas as gpd
import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import box

import myflopy as mf
from myflopy.modflow.mf6.canonical_example import (
    CanonicalModelConfig,
    build_canonical_model,
)
from myflopy.modflow.mf6.pest.zones import ZoneSpec, resolve_zone_array


@pytest.fixture(scope="module")
def zone_model(tmp_path_factory):
    """A canonical model to resolve zones against (grid only; never run)."""

    return build_canonical_model(
        tmp_path_factory.mktemp("zones") / "m", config=CanonicalModelConfig.testing()
    )


def _halves(model):
    """One zone definition: west half zone 1, east half zone 2."""

    xc = np.asarray(model.gwf.modelgrid.xcellcenters).reshape(-1)
    return np.where(xc < np.median(xc), 1, 2)


def test_one_definition_resolves_to_each_family_s_own_shape(zone_model):
    """The whole point of the front door. Measured pyEMU contract (2026-07-30):
    array targets take a per-cell array and REJECT ``(nlay, ncpl)``; list targets
    take ``(nlay, ncpl)`` and reject per-cell with
    ``IndexError: index 1 is out of bounds for axis 1``. One definition has to
    serve both, because a zonation is a property of the map, not of the package.
    """

    ncpl = int(zone_model.vor.ncpl)
    nlay = int(zone_model.gwf.modelgrid.nlay)
    spec = ZoneSpec.from_array(_halves(zone_model))

    array_zones = resolve_zone_array(spec, family="array", model=zone_model)
    list_zones = resolve_zone_array(spec, family="list", model=zone_model)

    assert array_zones.shape == (ncpl,)
    assert list_zones.shape == (nlay, ncpl)
    # Same labels either way -- the shape changes, the zonation does not.
    assert np.array_equal(list_zones[0], array_zones)
    assert np.array_equal(list_zones[-1], array_zones)


def test_a_single_layer_list_target_is_padded_past_pyemus_reshape(zone_model):
    """pyEMU cannot zone a list target on a one-layer vertex grid at all: the
    CORRECT ``(1, ncpl)`` shape hits ``checker2`` (pst_from.py:2121-2134), which
    assumes a ``(1, n)`` array on a vertex grid is an idomain array and reshapes
    it to ``(n, 1)``. Every lookup then fails with the same cryptic IndexError as
    a wrong shape. Padding to two rows sidesteps it; the extra row addresses no
    real cell."""

    from types import SimpleNamespace

    ncpl = int(zone_model.vor.ncpl)
    one_layer = SimpleNamespace(
        vor=zone_model.vor,
        gwf=SimpleNamespace(modelgrid=SimpleNamespace(nlay=1)),
    )
    resolved = resolve_zone_array(
        ZoneSpec.from_array(_halves(zone_model)), family="list", model=one_layer
    )

    assert resolved.shape == (2, ncpl), "a (1, ncpl) array would be reshaped by pyEMU"
    assert np.array_equal(resolved[0], resolved[1])


def test_a_categorical_raster_is_sampled_by_majority_not_mean(zone_model, tmp_path):
    """Zone ids are labels, not quantities. The area-weighted MEAN this repo uses
    for surfaces would turn a cell straddling zones 1 and 3 into zone 2 -- a zone
    the cell does not touch and which may not exist at all."""

    cells = zone_model.vor.gdf_vorPolys
    minx, miny, maxx, maxy = cells.total_bounds
    width = height = 60
    resolution = (maxx - minx) / width
    # West half labelled 1, east half labelled 3. Their mean is 2.
    band = np.repeat(
        np.where(np.arange(width)[None, :] < width // 2, 1, 3).astype("int32"),
        height, axis=0,
    )[:height]
    path = tmp_path / "hsu.tif"
    with rasterio.open(
        path, "w", driver="GTiff", height=height, width=width, count=1,
        dtype="int32", crs=cells.crs,
        transform=from_origin(minx, maxy, resolution, resolution), nodata=-9999,
    ) as dst:
        dst.write(band, 1)

    spec = ZoneSpec.from_raster(path, zone_model)

    assert spec.ids == [1, 3], f"expected the raster's own labels, got {spec.ids}"
    assert 2 not in spec.ids, "a mean-of-categories invented a zone"
    xc = np.asarray(zone_model.gwf.modelgrid.xcellcenters).reshape(-1)
    middle = minx + (maxx - minx) / 2
    assert np.all(spec.values[xc < middle - resolution] == 1)
    assert np.all(spec.values[xc > middle + resolution] == 3)


def test_polygon_zones_take_their_label_from_a_column(zone_model):
    """`docs/myflopy_context.md` claimed polygon zones already existed for
    `parameterize`; they did not -- `zone_array` had exactly one occurrence in
    all of src/, the raw pass-through."""

    cells = zone_model.vor.gdf_vorPolys
    minx, miny, maxx, maxy = cells.total_bounds
    middle = minx + (maxx - minx) / 2
    polygons = gpd.GeoDataFrame(
        {"hsu": [7, 9]},
        geometry=[box(minx, miny, middle, maxy), box(middle, miny, maxx, maxy)],
        crs=cells.crs,
    )

    spec = ZoneSpec.from_polygons(polygons, zone_model, column="hsu")
    assert spec.ids == [7, 9]

    xc = np.asarray(zone_model.gwf.modelgrid.xcellcenters).reshape(-1)
    assert np.all(spec.values[xc < middle] == 7)

    with pytest.raises(KeyError, match="not in the polygon layer"):
        ZoneSpec.from_polygons(polygons, zone_model, column="nope")


def test_a_wrong_length_definition_says_so_instead_of_reaching_pyemu(zone_model):
    """pyEMU's own message for this is `write_array_tpl() error: passed shape
    (441, 1) != zone_array.shape (4, 441)`, which names no target and no
    parameter."""

    with pytest.raises(ValueError, match="one integer label per cell"):
        resolve_zone_array(np.ones(7, dtype=int), family="array", model=zone_model)

    with pytest.raises(ValueError, match="contains NaN"):
        ZoneSpec.from_array(np.array([1.0, np.nan, 2.0]))


def test_per_layer_zones_that_disagree_are_refused(zone_model):
    """Accepting them would mean silently taking layer 0 and discarding the
    rest, since the resolved array is one label per cell."""

    ncpl = int(zone_model.vor.ncpl)
    nlay = int(zone_model.gwf.modelgrid.nlay)
    stacked = np.ones((nlay, ncpl), dtype=int)
    stacked[1] = 2

    with pytest.raises(ValueError, match="differs between layers"):
        resolve_zone_array(stacked, family="array", model=zone_model)

    # Identical layers are fine -- that is just a per-cell definition spelled out.
    stacked[1] = 1
    assert resolve_zone_array(stacked, family="array", model=zone_model).shape == (ncpl,)


def test_zone_zero_warns_because_the_two_families_disagree_about_it(zone_model):
    """pyEMU skips zone ids below 1 for ARRAY targets but makes a real adjustable
    parameter for zone 0 on LIST targets. Normalizing the shape cannot normalize
    that, so it is called out rather than papered over."""

    values = _halves(zone_model).copy()
    values[:5] = 0

    with pytest.warns(UserWarning, match="Zone id 0"):
        resolve_zone_array(values, family="list", model=zone_model)

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        resolve_zone_array(values, family="array", model=zone_model)


@pytest.mark.slow
@pytest.mark.parametrize("target", ["k", "recharge"])
def test_parameterize_accepts_one_per_cell_definition_for_either_family(tmp_path, target):
    """The wiring, not the resolver. A per-cell definition on a LIST target is
    exactly what pyEMU rejects with `IndexError: index 1 is out of bounds for
    axis 1 with size 1`, so this fails if `parameterize` hands `zones=` straight
    through — which is what it did before this item."""

    pytest.importorskip("pyemu")
    from myflopy.modflow.mf6.canonical_calibration import (
        build_canonical_calibration_demo,
    )

    demo = build_canonical_calibration_demo(
        tmp_path / "model", config=CanonicalModelConfig.testing(), n_head_wells=4
    )
    cal = demo.model.pest("zoned", start_datetime="2024-01-01")
    cal.parameterize(target, style="zone", zones=_halves(demo.model),
                     bounds=(0.5, 2.0))
    cal.observe(demo.head_targets)
    pst = cal.build("zoned.pst", noptmax=0)

    assert pst.npar_adj >= 1
    zone_ids = {
        name.split("_zone:")[-1].split("_")[0]
        for name in pst.parameter_data.index if "_zone:" in name
    }
    assert zone_ids, "no zone-typed parameters were created"


def test_zone_spec_is_on_the_public_surface():
    assert mf.ZoneSpec is ZoneSpec
    assert "ZoneSpec" in set(mf.__preferred__)


# --- cells too small to hold a pixel centre (ledger 213) -------------------- #


def _raster(path, band, *, bounds, crs, nodata=-9999):
    minx, miny, maxx, maxy = bounds
    height, width = band.shape
    with rasterio.open(
        path, "w", driver="GTiff", height=height, width=width, count=1, dtype="int32",
        crs=crs, transform=from_origin(minx, maxy, (maxx - minx) / width, (maxy - miny) / height),
        nodata=nodata,
    ) as dst:
        dst.write(band.astype("int32"), 1)
    return path


def test_a_cell_too_small_for_a_pixel_centre_takes_the_label_under_its_centroid(zone_model, tmp_path):
    """rasterize burns a pixel only where its CENTRE lands, so a cell finer than the
    raster held none and took zone 0 -- 24% of a refined grid, frozen."""

    cells = zone_model.vor.gdf_vorPolys
    band = np.array([[1, 3], [5, 7]])  # four pixels, each far larger than a cell
    path = _raster(tmp_path / "coarse.tif", band, bounds=cells.total_bounds, crs=cells.crs)

    spec = ZoneSpec.from_raster(path, zone_model)

    centroids = cells.geometry.centroid
    with rasterio.open(path) as src:
        expected = np.array([band[src.index(x, y)] for x, y in zip(centroids.x, centroids.y, strict=True)])
    assert 0 not in spec.ids
    assert np.array_equal(spec.values, expected)


def test_a_centroid_off_the_raster_keeps_nodata_zone_without_a_declared_nodata(zone_model, tmp_path):
    """rasterio samples its fill -- 0 -- off a raster that declares no nodata."""

    cells = zone_model.vor.gdf_vorPolys
    minx, miny, maxx, maxy = cells.total_bounds
    middle = (minx + maxx) / 2
    path = _raster(tmp_path / "west.tif", np.full((1, 1), 4), bounds=(minx, miny, middle, maxy),
                   crs=cells.crs, nodata=None)

    spec = ZoneSpec.from_raster(path, zone_model, nodata_zone=-1)

    xc = cells.geometry.centroid.x.to_numpy()
    assert set(spec.values[xc > middle + 1.0].tolist()) == {-1}
    assert set(spec.values[xc < middle - 1.0].tolist()) == {4}


def test_a_centroid_on_nodata_keeps_nodata_zone(zone_model, tmp_path):
    cells = zone_model.vor.gdf_vorPolys
    minx, _, maxx, _ = cells.total_bounds
    band = np.array([[-9999, 3], [-9999, 3]])
    path = _raster(tmp_path / "half.tif", band, bounds=cells.total_bounds, crs=cells.crs)

    spec = ZoneSpec.from_raster(path, zone_model, nodata_zone=-1)

    xc = cells.geometry.centroid.x.to_numpy()
    middle = (minx + maxx) / 2
    assert set(spec.values[xc < middle - 1.0].tolist()) == {-1}
    assert set(spec.values[xc > middle + 1.0].tolist()) == {3}


def test_the_majority_still_decides_a_cell_that_holds_pixel_centres(zone_model, tmp_path):
    """One odd pixel under a centroid does not outvote the rest of the cell."""

    cells = zone_model.vor.gdf_vorPolys
    minx, miny, maxx, maxy = cells.total_bounds
    band = np.ones((60, 60), dtype=int)
    resolution = (maxx - minx) / 60, (maxy - miny) / 60
    target = len(cells) // 2
    centre = cells.geometry.iloc[target].centroid
    band[int((maxy - centre.y) / resolution[1]), int((centre.x - minx) / resolution[0])] = 9
    path = _raster(tmp_path / "fine.tif", band, bounds=(minx, miny, maxx, maxy), crs=cells.crs)

    spec = ZoneSpec.from_raster(path, zone_model)

    assert spec.values[target] == 1
    assert 9 not in spec.ids


# --- zone rasters: memory by labels, NaN pixels never vote ------------------------------- #


def _zone_raster(path, band, cells, *, nodata, dtype):
    import rasterio
    from rasterio.transform import from_origin

    minx, miny, maxx, maxy = cells.total_bounds
    height, width = band.shape
    with rasterio.open(
        path, "w", driver="GTiff", height=height, width=width, count=1, dtype=dtype,
        crs=cells.crs, nodata=nodata,
        transform=from_origin(minx, maxy, (maxx - minx) / width, (maxy - miny) / height),
    ) as dst:
        dst.write(band.astype(dtype), 1)
    return path


def test_a_sparse_label_raster_costs_memory_by_its_labels_not_their_span(zone_model, tmp_path):
    """Labels {1, 2, 65535} (an undeclared uint16 nodata): the per-cell count table was
    ncells x 65,535 int64 -- 463 MB peak at 441 cells."""

    cells = zone_model.vor.gdf_vorPolys
    band = np.ones((60, 60), dtype=int)
    band[:, 30:] = 2
    band[:5, :5] = 65535
    path = _zone_raster(tmp_path / "sparse.tif", band, cells, nodata=None, dtype="uint16")
    tracemalloc.start()
    spec = ZoneSpec.from_raster(path, zone_model)
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    assert peak < 50e6, f"{peak / 1e6:.0f} MB"
    assert set(spec.ids) <= {1, 2, 65535} and {1, 2} <= set(spec.ids)


def test_a_float_zone_raster_with_nan_nodata_is_read(zone_model, tmp_path):
    """NaN never equals itself, so NaN pixels voted and became int64-min: OverflowError."""

    cells = zone_model.vor.gdf_vorPolys
    band = np.full((60, 60), 4.0)
    band[:, 30:] = np.nan
    path = _zone_raster(tmp_path / "nan.tif", band, cells, nodata=float("nan"), dtype="float32")
    spec = ZoneSpec.from_raster(path, zone_model, nodata_zone=0)
    assert set(spec.ids) == {0, 4}
    xc = cells.geometry.centroid.x.to_numpy()
    minx, _, maxx, _ = cells.total_bounds
    assert set(spec.values[xc < (minx + maxx) / 2 - 200.0].tolist()) == {4}
    assert set(spec.values[xc > (minx + maxx) / 2 + 200.0].tolist()) == {0}

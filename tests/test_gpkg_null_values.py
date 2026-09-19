"""A GeoPackage NULL must not become a non-finite number in a written MF6 package.

Two paths, one hole. ``k_array`` did ``float(row[value])``, which turns a null
into NaN; the boundary record kept the ``None``, and numpy casts that to NaN when
FloPy builds the recarray. Either way MF6 was written a non-finite value and
answered a NaN budget under "Normal termination" -- measured: a NULL ``k`` gave
cell 5 = nan where the legacy KFromVector's ``fillna`` gave the layer default,
and a NULL head/conductance gave the record ``[(0, 0), None, None, 'a']``.

Substituting a default for a value the user never supplied is its own silent
guess, so the row is named instead.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import box

from myflopy.geopackage import CellSurfaceOffset, GeoPackageSource, Spread
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.specs import ModelContext

N = 3
D = 25.0
CRS = 2927


@pytest.fixture
def vor() -> VoronoiGridPlus:
    """A 3 x 3 grid of 25 ft square cells, row-major from the bottom left."""

    verts = np.array([[i * D, j * D] for j in range(N + 1) for i in range(N + 1)], dtype=float)
    iverts = [
        [j * (N + 1) + i, (j + 1) * (N + 1) + i, (j + 1) * (N + 1) + i + 1, j * (N + 1) + i + 1]
        for j in range(N)
        for i in range(N)
    ]
    xcyc = np.array([[(i + 0.5) * D, (j + 0.5) * D] for j in range(N) for i in range(N)])
    return VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc, crs=CRS)


def _source(tmp_path, vor, columns: dict, *, geometry=None, nper=1, **kwargs) -> GeoPackageSource:
    """A one-feature source over cell 0 unless ``geometry`` says otherwise."""

    gdf = gpd.GeoDataFrame(
        {"name": ["a"], "layer": [1], **columns},
        geometry=[box(0, 0, D, D) if geometry is None else geometry],
        crs=CRS,
    )
    path = tmp_path / "bc.gpkg"
    gdf.to_file(path, layer="bc", driver="GPKG")
    return GeoPackageSource(
        path=path, context=ModelContext(grid=vor, domain=None), nper=nper, layer="bc", **kwargs
    )


def test_a_null_head_names_the_feature_and_the_field(tmp_path, vor):
    """The error a user gets instead of a NaN CHD record."""

    source = _source(tmp_path, vor, {"head": [None]})
    with pytest.raises(ValueError, match=r"feature 'a' has no usable number for 'head'"):
        source.chd(head="head")


def test_a_null_says_the_attribute_was_empty(tmp_path, vor):
    """A NULL is described as empty, not as the word None."""

    source = _source(tmp_path, vor, {"head": [None]})
    with pytest.raises(ValueError, match=r"empty \(a GeoPackage NULL\)"):
        source.chd(head="head")


def test_a_null_conductance_names_conductance_not_the_first_field(tmp_path, vor):
    """Each value field is checked on its own, so the message points at the bad one."""

    source = _source(tmp_path, vor, {"head": [9.0], "conductance": [None]})
    with pytest.raises(ValueError, match=r"no usable number for 'conductance'"):
        source.ghb(head="head", conductance="conductance")


def test_a_nan_in_a_float_column_is_refused(tmp_path, vor):
    """A null read back into a float column arrives as NaN rather than None."""

    source = _source(tmp_path, vor, {"head": [float("nan")]})
    with pytest.raises(ValueError, match=r"feature 'a' has nan for 'head'"):
        source.chd(head="head")


def test_an_infinite_value_is_refused(tmp_path, vor):
    """MF6 cannot use an infinity either."""

    source = _source(tmp_path, vor, {"head": [float("inf")]})
    with pytest.raises(ValueError, match=r"has inf for 'head'"):
        source.chd(head="head")


def test_a_text_value_is_refused(tmp_path, vor):
    """A column of the wrong type fails here, naming the field, not later inside FloPy."""

    source = _source(tmp_path, vor, {"head": ["high"]})
    with pytest.raises(ValueError, match=r"no usable number for 'head' -- it is 'high'"):
        source.chd(head="head")


def test_a_null_in_one_periods_column_names_the_per_period_fields(tmp_path, vor):
    """A per-period list is only as good as its worst column."""

    source = _source(tmp_path, vor, {"h0": [9.0], "h1": [None]}, nper=2)
    with pytest.raises(ValueError, match=r"its per-period fields 'h0', 'h1'"):
        source.chd(head=["h0", "h1"])


def test_a_null_in_a_spread_field_names_it_as_a_spread(tmp_path, vor):
    """Spread multiplies by a share; a null would spread NaN over every cell."""

    source = _source(
        tmp_path, vor, {"elevation": [5.0], "conductance": [None]},
        geometry=box(0, 0, 2 * D, D),
    )
    with pytest.raises(ValueError, match=r"'conductance' \(spread\)"):
        source.drn(elevation="elevation", conductance=Spread("conductance"))


def test_a_null_offset_on_a_surface_reference_names_the_column(tmp_path, vor):
    """The OFFSET column is read through the same funnel now, so it is caught AT the read.

    It used to be caught downstream, once the SUM went non-finite -- which named the
    fields but not the feature and not which of `offset=`/`minimum=` was empty. The
    `minimum` column had no such luck in either direction: a NULL read back as `nan`
    SILENTLY DISABLED the floor (`max(value, nan)` returns `value`, and the result is
    finite so this guard passed), and a TEXT NULL raised a bare `TypeError` naming
    nothing. Both columns go through `_finite` now (ledger 289), which is what makes
    ledger 247's own claim -- that it names the offset and minimum columns -- true.
    """

    vor.gdf_topbtm = gpd.GeoDataFrame(
        {0: [10.0] * (N * N), 1: [0.0] * (N * N)},
        geometry=list(vor.gdf_vorPolys.geometry),
        crs=CRS,
    )
    source = _source(tmp_path, vor, {"drop": [float("nan")], "conductance": [1.0]})
    with pytest.raises(ValueError, match=r"feature 'a' \(offset=\) has nan for 'drop'"):
        source.drn(
            elevation=CellSurfaceOffset("cell_top", offset="drop"), conductance="conductance"
        )


def test_a_null_minimum_is_refused_rather_than_silently_disabling_the_floor(tmp_path, vor):
    """`max(value, nan)` returns `value`, and the result is FINITE -- so this guard passed.

    A caller who asked for a floor got none, and the written output was byte-identical to
    declaring no floor at all: measured, 11 of 11 records below their own cell bottom
    with zero warnings (ledger 289). That is the worst shape of a silent input defect --
    the safeguard the user reached for is the thing that vanished.
    """

    vor.gdf_topbtm = gpd.GeoDataFrame(
        {0: [10.0] * (N * N), 1: [0.0] * (N * N)},
        geometry=list(vor.gdf_vorPolys.geometry),
        crs=CRS,
    )
    source = _source(tmp_path, vor, {"floor": [float("nan")], "conductance": [1.0]})

    with pytest.raises(ValueError, match=r"\(minimum=\) has nan for 'floor'"):
        source.drn(
            elevation=CellSurfaceOffset("cell_top", offset=-20.0, minimum="floor"),
            conductance="conductance",
        )


def test_a_real_minimum_still_floors(tmp_path, vor):
    """The refusing direction: the guard must not break a floor that is actually given."""

    vor.gdf_topbtm = gpd.GeoDataFrame(
        {0: [10.0] * (N * N), 1: [0.0] * (N * N)},
        geometry=list(vor.gdf_vorPolys.geometry),
        crs=CRS,
    )
    source = _source(tmp_path, vor, {"floor": [4.0], "conductance": [1.0]})

    spec = source.drn(
        elevation=CellSurfaceOffset("cell_top", offset=-20.0, minimum="floor"),
        conductance="conductance",
    )

    assert {record[1] for record in spec.options["stress_period_data"][0]} == {4.0}


def test_a_feature_that_lands_nowhere_reports_the_MISS_not_the_null(tmp_path, vor):
    """Location is the more fundamental complaint, and keeps ledger 211's message."""

    source = _source(tmp_path, vor, {"head": [None]}, geometry=box(-999, -999, -998, -998))
    with pytest.raises(ValueError, match=r"lands on an active cell"):
        source.chd(head="head")


def test_a_k_zone_that_lands_nowhere_reports_the_MISS_not_the_null(tmp_path, vor):
    """k_array has its own loop, so it needs its own precedence test.

    A mutant that dropped the ``continue`` here survived the CHD version of this
    test: ``_boundary_data`` already skipped a missed feature, ``k_array`` did not.
    """

    source = _source(tmp_path, vor, {"k": [None]}, geometry=box(-999, -999, -998, -998))
    with pytest.raises(ValueError, match=r"lands on an active cell"):
        source.k_array(value="k", nlay=1, defaults=1.0)


def test_a_null_k_attribute_names_the_zone(tmp_path, vor):
    """k_array's own conversion had the same hole; legacy filled it with the default."""

    source = _source(tmp_path, vor, {"k": [None]})
    with pytest.raises(ValueError, match=r"feature 'a' has no usable number for 'k'"):
        source.k_array(value="k", nlay=1, defaults=1.0)


def test_a_usable_k_still_paints_and_the_defaults_still_fill(tmp_path, vor):
    """The control: the guard refuses bad numbers and passes good ones through."""

    source = _source(tmp_path, vor, {"k": [50.0]})
    array = source.k_array(value="k", nlay=1, defaults=1.0)
    assert array[0, 0] == 50.0
    assert np.isfinite(array).all()
    assert sorted(set(array[0])) == [1.0, 50.0]


def test_an_integer_value_is_accepted(tmp_path, vor):
    """An int column is a number; the guard must not insist on a float."""

    source = _source(tmp_path, vor, {"head": [9]})
    records = source.chd(head="head").options["stress_period_data"][0]
    assert [r[1] for r in records] == [9]


def test_a_scalar_passed_in_code_is_still_accepted(tmp_path, vor):
    """A constant given at the call site never came from the file, but is still checked."""

    source = _source(tmp_path, vor, {"head": [9.0]})
    records = source.chd(head=4.5).options["stress_period_data"][0]
    assert [r[1] for r in records] == [4.5]


def test_a_scalar_nan_passed_in_code_is_refused(tmp_path, vor):
    """...and a constant NaN is refused the same way a column one is."""

    source = _source(tmp_path, vor, {"head": [9.0]})
    with pytest.raises(ValueError, match=r"has nan for"):
        source.chd(head=float("nan"))


def test_a_boolean_is_not_a_number(tmp_path, vor):
    """True is not an elevation; bool is an int subclass, so it is excluded on purpose."""

    source = _source(tmp_path, vor, {"head": [True]})
    with pytest.raises(ValueError, match=r"no usable number for 'head'"):
        source.chd(head="head")

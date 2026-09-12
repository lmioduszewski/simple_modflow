from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely as shp
from shapely.geometry import LineString, Polygon

from myflopy._defaults import DEFAULT_CRS
from myflopy.modflow.utils.datatypes.readers import read_shp_gpkg


def get_voronoi_polygons(vor) -> list[Polygon]:
    """
    Build Shapely polygons for each Voronoi cell from vertex arrays.
    """
    polygons = []
    for region in vor.iverts:
        vertices = [vor.verts[j] for j in region]
        polygons.append(Polygon(vertices))
    return polygons


def generate_grid_coordinates(vor, grid_spacing: int) -> tuple[list[int], list[int]]:
    """
    Generate evenly spaced coordinates spanning the grid extent.
    """
    if hasattr(vor, 'x_vor') and hasattr(vor, 'y_vor'):
        xs = np.asarray(vor.x_vor)
        ys = np.asarray(vor.y_vor)
    else:
        verts = np.asarray(vor.verts)
        xs = verts[:, 0]
        ys = verts[:, 1]

    x_min = xs.min()
    x_max = xs.max()
    y_min = ys.min()
    y_max = ys.max()
    x_coords = []
    y_coords = []
    for x in range(int(x_min), int(x_max) + grid_spacing, grid_spacing):
        for y in range(int(y_min), int(y_max) + grid_spacing, grid_spacing):
            x_coords.append(x)
            y_coords.append(y)
    return x_coords, y_coords


def get_centroids(vor) -> tuple[list[float], list[float]]:
    """
    Return centroid coordinates for each Voronoi polygon.
    """
    if vor._gdf_vorPolys is None:
        vor.get_gdf_vorPolys()

    centroids = vor._gdf_vorPolys.geometry.centroid
    return centroids.x.tolist(), centroids.y.tolist()


def find_adjacent_polygons(gdf: gpd.GeoDataFrame) -> list[list[int]]:
    """
    Determine which polygons in the grid share an edge.
    """
    adjacent = [[] for _ in range(len(gdf))]
    gdf_sindex = gdf.sindex

    for i, poly in gdf.geometry.items():
        candidates = list(gdf_sindex.intersection(poly.bounds))

        for j in candidates:
            if i == j:
                continue

            shared_edge = poly.intersection(gdf.iloc[j].geometry)
            if isinstance(shared_edge, LineString):
                adjacent[i].append(j)

    return adjacent


def find_adjacent_cells(vor, cell_id: int) -> list[int]:
    """
    Return neighboring cell ids from the DISU connectivity vectors.
    """
    cell_id = int(cell_id)  # tolerate numpy ints / float-typed ids from spatial joins
    start_index = int(vor.ia[cell_id])  # O(1) via cached row pointer, not sum(iac[:cell])
    num_connections = int(vor.iac[cell_id])
    connections = vor.ja[start_index:start_index + num_connections]
    return [int(cell) for cell in connections if int(cell) != cell_id]


def calculate_distance(gdf: gpd.GeoDataFrame, poly_idx1: int, poly_idx2: int) -> float:
    """
    Calculate centroid-to-shared-face distance for two adjacent polygons.
    """
    poly1 = gdf.iloc[poly_idx1].geometry
    poly2 = gdf.iloc[poly_idx2].geometry
    shared_edge = poly1.intersection(poly2)
    centroid = poly1.centroid
    return centroid.distance(shared_edge)


def shared_face_length(poly1, poly2) -> float:
    """
    Calculate the length of the shared face between adjacent polygons.
    """
    intersection = poly1.intersection(poly2)
    if not isinstance(intersection, LineString):
        return 0
    return intersection.length


def get_gdf_vor_polys(vor, crs=None) -> gpd.GeoDataFrame:
    """
    Build and cache the GeoDataFrame of Voronoi polygons.
    """
    crs = vor.crs if crs is None else crs
    vertices_by_cells = []
    xvertices_by_cells = []
    yvertices_by_cells = []
    polygons = []

    for cell in vor.iverts:
        thiscell = []
        thiscellx = []
        thiscelly = []
        for vidx in cell:
            thisvert = vor.verts[vidx].tolist()
            thiscell.append(thisvert)
            thiscellx.append(thisvert[0])
            thiscelly.append(thisvert[1])
        xvertices_by_cells.append(thiscellx)
        yvertices_by_cells.append(thiscelly)
        vertices_by_cells.append(thiscell)

    for cell in vertices_by_cells:
        polygons.append(shp.Polygon(cell))

    gdf = gpd.GeoDataFrame(geometry=polygons, crs=crs)
    vor.gdf_vorPolys = gdf
    vor.x_coords_by_node = xvertices_by_cells
    vor.y_coords_by_node = yvertices_by_cells
    return gdf


def get_gdf_latlon(vor):
    """
    Return the Voronoi polygons in latitude/longitude coordinates.
    """
    if vor._gdf_latlon is None:
        vor._gdf_latlon = vor.gdf_vorPolys.to_crs(vor.crs_latlon)
    return vor._gdf_latlon


def get_latlon(vor):
    """
    Return GeoJSON-like lat/lon geometry for the Voronoi grid.
    """
    if vor._latlon is None:
        vor._latlon = json.loads(vor.gdf_latlon["geometry"].to_json())
    return vor._latlon


def points_to_latlon(vor, x, y):
    """Reproject model-coordinate points onto the map frame, as ``(lon, lat)`` arrays.

    The bridge between model results carrying raw ``x``/``y`` -- PRT track
    records, observation locations -- and the WGS84 space the choropleth lives
    in. The grid's own vertices are reprojected the same way
    (:func:`get_gdf_latlon`), so points land where the cells are.
    """

    crs = getattr(vor, "crs", None)
    if crs is None:
        raise ValueError(
            "Cannot place points on the map: the grid has no CRS. Build the grid "
            "with VoronoiGridPlus(..., crs='EPSG:XXXX') so model coordinates can "
            "be reprojected to lat/lon."
        )
    points = gpd.GeoSeries(gpd.points_from_xy(x, y), crs=crs).to_crs(vor.crs_latlon)
    return points.x.to_numpy(), points.y.to_numpy()


def get_grid_centroid(vor):
    """
    Return the centroid of the Voronoi domain in lat/lon coordinates.
    """
    return vor.gdf_latlon.union_all().centroid


def get_overlapping_area(vor, shp_gpkg=None, cell_list=None):
    """
    Return the area of cells overlapping a geometry or explicit cell ids.
    """
    if shp_gpkg is not None:
        cells = vor.get_vor_cells_as_series(shp_gpkg).to_list()
    elif cell_list is not None:
        cells = cell_list
    else:
        raise ValueError('You must provide a shp_gpkg or cell_list')
    return vor.gdf_vorPolys.loc[cells].union_all().area


def get_vor_idx_from_geometry(
    vor,
    shp_to_query: shp = None,
    gdf_to_query: gpd.GeoDataFrame = None,
    crs: str = DEFAULT_CRS,
    crs_latlon: str = "EPSG:4326",
    name_col: str = None,
    predicate: str = "intersects",
) -> dict:
    """
    Determine intersecting Voronoi cells for each geometry in a query layer.
    """
    if shp_to_query:
        gdf_query = read_shp_gpkg(shp_to_query).to_crs(crs)
    elif gdf_to_query is not None:
        gdf_query = gdf_to_query.to_crs(crs)
    else:
        raise ValueError('Provide shp_to_query or gdf_to_query')

    vor_idx_dict = {}
    latlon_vor = vor.gdf_vorPolys["geometry"].to_crs(crs_latlon)
    latlon_query = gdf_query.to_crs(crs_latlon)

    for idx in latlon_query.index:
        cell_intersections = (
            gpd.GeoSeries(latlon_query.iloc[idx]["geometry"])
            .sindex.query(latlon_vor, predicate=predicate)[0]
            .tolist()
        )
        key = gdf_query.iloc[idx][name_col] if name_col is not None else idx
        vor_idx_dict[key] = cell_intersections

    return vor_idx_dict


def get_vor_idx_from_geometry_idx(
    vor,
    gdf_to_query: gpd.GeoDataFrame = None,
    idx: int = 0,
    predicate: str = "intersects",
) -> list[int]:
    """
    Determine intersecting Voronoi cells for a single geometry in a GeoDataFrame.
    """
    return (
        vor.gdf_vorPolys["geometry"]
        .sindex.query(gdf_to_query['geometry'].iloc[idx], predicate=predicate, sort=True)
        .tolist()
    )


def set_k_vor(vor, k_dict: dict = None, k_default=100) -> list:
    """
    Set hydraulic conductivity values on the Voronoi polygon GeoDataFrame.
    """
    vor.gdf_vorPolys["Kh"] = k_default
    for key in k_dict.keys():
        vor.gdf_vorPolys.loc[vor.gdf_vorPolys.index.isin(k_dict[key]), ["Kh"]] = key
    return vor.gdf_vorPolys["Kh"].to_list()


def generate_grid_around_point(center_point: shp.Point, spacing: float, size: int, crs: str) -> gpd.GeoSeries:
    """
    Generate a GeoSeries of grid points around a center point.
    """
    minx = center_point.x - size / 2 * spacing
    miny = center_point.y - size / 2 * spacing
    maxx = center_point.x + size / 2 * spacing
    maxy = center_point.y + size / 2 * spacing

    x_coords = list(range(int(minx), int(maxx) + 1, spacing))
    y_coords = list(range(int(miny), int(maxy) + 1, spacing))
    points = [shp.Point(x, y) for x in x_coords for y in y_coords]
    return gpd.GeoSeries(points).set_crs(crs)


def generate_grid_polygons(center_point, spacing, size, gap):
    """
    Generate square polygons in a grid pattern around a center point.
    """
    minx = center_point.x - size / 2 * spacing
    miny = center_point.y - size / 2 * spacing
    maxx = center_point.x + size / 2 * spacing
    maxy = center_point.y + size / 2 * spacing

    x_coords = list(range(int(minx), int(maxx) + 1, spacing))
    y_coords = list(range(int(miny), int(maxy) + 1, spacing))

    polygons = []
    for x in x_coords[:-1]:
        for y in y_coords[:-1]:
            polygons.append(shp.box(x + gap / 2, y + gap / 2, x + spacing - gap / 2, y + spacing - gap / 2))

    return shp.MultiPolygon(polygons)


def voronoi_refine_by_point(point: shp.Point, spacing: int, tri) -> gpd.GeoDataFrame:
    """
    Add a small refinement stencil around a point to a Triangle grid.
    """
    polypoints = [
        (point.x, point.y),
        (point.x + spacing, point.y),
        (point.x + spacing, point.y - spacing),
        (point.x, point.y - spacing),
    ]
    poly_main = shp.Polygon(polypoints)
    transform_dist = spacing * 2

    poly_e = shp.transform(poly_main, lambda x: x + [transform_dist, 0])
    poly_w = shp.transform(poly_main, lambda x: x - [transform_dist, 0])
    poly_n = shp.transform(poly_main, lambda x: x + [0, transform_dist])
    poly_s = shp.transform(poly_main, lambda x: x - [0, transform_dist])
    poly_ne = shp.transform(poly_main, lambda x: x + [transform_dist, transform_dist])
    poly_nw = shp.transform(poly_main, lambda x: x + [-transform_dist, transform_dist])
    poly_se = shp.transform(poly_main, lambda x: x + [transform_dist, -transform_dist])
    poly_sw = shp.transform(poly_main, lambda x: x + [-transform_dist, -transform_dist])
    all_polys = [poly_main, poly_e, poly_w, poly_n, poly_s, poly_ne, poly_nw, poly_se, poly_sw]

    gdf_all_polys = gpd.GeoDataFrame(geometry=all_polys)
    for poly_idx in range(len(gdf_all_polys)):
        tri.add_polygon(gdf_all_polys.loc[poly_idx, 'geometry'])

    return gdf_all_polys


def geometry_first(gdf):
    """Return ``gdf`` with its geometry column(s) first, elevations in order.

    ``gdf_topbtm`` is read POSITIONALLY by several consumers -- ``iloc[:, 1]`` for
    the model top, ``iloc[:, 2:]`` for the bottoms -- and nothing enforced where
    the geometry column sat. :meth:`LayerSurfaces.sample` with ``reconcile=True``
    (the default) drops geometry, reconciles, then re-wraps with
    ``gpd.GeoDataFrame(df, geometry=...)``, which appends geometry at the **end**
    and undoes the explicit geometry-first ordering three lines above it. So a
    keyword about physics moved the column every positional reader depends on.

    Measured on a three-layer frame: the layer-elevation hover then read the first
    bottom as the model top (5.0 where the truth was 10.0) and handed a Polygon
    back as the last bottom, raising nothing. Integer labels do not protect
    against it, which is why unifying the column NAMES would not have helped.

    Geometry is found by DTYPE rather than by the name ``"geometry"``: a frame
    whose geometry column is called ``geom`` was otherwise swept in as a layer
    bottom. Idempotent, and returns the frame unchanged when it is already
    ordered, so callers can apply it on every read without copying.
    """

    if gdf is None or getattr(gdf, "columns", None) is None:
        return gdf
    geo = [c for c in gdf.columns if str(getattr(gdf[c], "dtype", "")) == "geometry"]
    if not geo or list(gdf.columns)[: len(geo)] == geo:
        return gdf
    return gdf[geo + [c for c in gdf.columns if c not in geo]]


# An areal piece narrower than this, in the grid's length unit, is dropped. Almost
# always it is a face-aligned edge whose coordinates were perturbed. A
# perturbation is a WIDTH, so the test is the piece's mean width, not its share of
# the cell: a share rule caught a NAD83 round trip on 500 ft cells and missed it on
# 300 ft ones. Measured: 3e-7 after a double-precision reprojection, 3.4e-4 after
# NAD83 <-> NAD83(HARN), up to 5e-3 from vertices FloPy wrote at 9 significant
# digits. It also drops REAL corner clips that thin, deliberately, as
# Spread.min_share does -- on real overlays at most 1.5e-3 ft2, 8.3e-6 of a cell.
# Coordinates ROUNDED on export (GeoJSON at 7 decimal degrees: 0.003 to 0.018)
# straddle it and are not reliably covered. Ledger 210 -- and SFR drops a reach
# SHORTER than it, whatever min_reach_length says (ledger 214).
_SLIVER_WIDTH = 0.01

# ...and a sliver is also under this share of its cell. Noise is width / cell width
# of a cell, under 1% on any cell wider than half a unit, so this changes nothing on
# a real grid; without it a lab-scale grid (0.01-unit cells) lost ALL its coverage,
# every piece being narrower than the tolerance.
_SLIVER_SHARE = 0.01


def _mean_width(geometry) -> float:
    """Twice the area over the perimeter of the AREAL parts: a strip's width, half a square's side.

    Areal parts only: a piece that also touches its cell along a face comes back as
    a GeometryCollection holding that face as a line, and counted in the perimeter
    it made a 1.2 ft square a 0.00998 "sliver".
    """

    if geometry.geom_type == "GeometryCollection":
        geometry = shp.union_all([part for part in shp.get_parts(geometry) if part.area > 0])
    perimeter = float(geometry.length)
    return 2.0 * float(geometry.area) / perimeter if perimeter > 0 else 0.0


def _extent(geometry) -> float:
    """A geometry's own measure: area for a polygon, length for a line, 0 for a point."""

    if geometry is None or geometry.is_empty:
        return 0.0
    area = float(getattr(geometry, "area", 0.0) or 0.0)
    return area if area > 0 else float(getattr(geometry, "length", 0.0) or 0.0)


def domain_active(domain, layer: int, cell: int) -> bool:
    """Whether ``(layer, cell)`` can carry a boundary, per an MF6 idomain.

    MF6 idomain is three-valued, not a flag: ``> 0`` active, ``0`` inactive, and
    ``< 0`` **vertical passthrough** -- the cell is removed and flow passes
    straight through it. A passthrough cell holds no boundary; MF6 refuses one
    with ``Cell is outside active grid domain``.

    So the test is ``> 0``, not truthiness. ``bool(-1)`` is ``True``, which put
    347 DRN records into passthrough cells on a real model and failed the run at
    read time -- a ``LayerStack`` using ``pinch="passthrough"`` produces these in
    quantity (13,797 on the model in question).

    A ``domain`` of ``None`` is a wholly active grid, and a one-dimensional one is
    the same activity in every layer. A layer beyond a ``(nlay, ncpl)`` domain
    raises rather than resolving against the last one.
    """

    if domain is None:
        return True
    values = np.asarray(domain)
    if values.ndim == 1:
        return bool(values[cell] > 0)
    if layer >= values.shape[0]:
        raise ValueError(
            f"Layer {layer} is outside the model domain, which has "
            f"{values.shape[0]} layers."
        )
    return bool(values[layer, cell] > 0)


def cell_shares(polys, geometry, *, usable=None) -> list[tuple[int, float]]:
    """The cells a geometry reaches, each with its fraction of the whole feature.

    ``polys`` is the grid's cell frame (``vor.gdf_vorPolys``), indexed by cell.
    ``usable(cells)`` is handed the cells a candidate touches and returns the
    eligible ones: it is how a caller keeps one layer's active cells, or only the
    grid edge, computing that set ONCE and closing over it rather than deciding
    per candidate. The default keeps every candidate. Either way the survivors are
    SORTED before they are measured and the pairs come back in cell order, which
    is load-bearing: the spatial index answers in tree order and a tie goes to the
    lowest-numbered cell, so a caller handing its cells back as a set must not be
    able to cost it that.

    The fraction is by length for a line, by area for a polygon, and 1.0 for a
    point or anything with no extent -- there is nothing to divide. A tie goes
    ONE way (decided 2026-09-11): a stretch of line lying on a face is in both
    cells' pieces, so it is split evenly between them and a line's shares sum to
    one rather than two; a point on a face or vertex goes whole to the first
    usable cell it touches, not whole to each. Fractions are of the WHOLE
    feature, so they sum to less than one where part of it lies outside the
    usable cells; ``myflopy.geopackage.Spread`` decides what that means.

    The FEATURE picks the measure, never the piece. A polygon meets each
    face-adjacent cell along the edge they share -- a line, with no area -- and
    measuring that piece by its length admitted every neighbour of a polygon
    drawn on cell faces, at the full value for any field not wrapped in
    ``Spread``: a K zone that was exactly one cell painted five. Pieces
    narrower than ``_SLIVER_WIDTH`` and under ``_SLIVER_SHARE`` of their cell are
    dropped for the same reason, unless every piece of the feature is one.

    Candidates come from the spatial index rather than an elementwise
    ``intersects`` over every cell, which was O(features x ncpl).
    """

    if geometry is None or geometry.is_empty:
        return []
    measure = _extent(geometry)
    if measure <= 0:
        chosen = set()
        for part in shp.get_parts(geometry):
            hits = _usable_candidates(polys, part, usable)
            if hits:
                chosen.add(hits[0])
        return [(cell, 1.0) for cell in sorted(chosen)]

    # A GeometryCollection is measured by its highest dimension: the lines in a
    # polygon collection add no area. In practice they are what gpd.clip and
    # gpd.overlay leave (keep_geom_type=False) where a zone ran along the clip
    # edge, so refusing the collection refused their default output.
    areal = float(getattr(geometry, "area", 0.0) or 0.0) > 0
    # A feature that is itself this thin is what was drawn, not a sliver of it.
    thin_feature = areal and _mean_width(geometry) <= _SLIVER_WIDTH
    sizes: dict[int, float] = {}
    slivers: set[int] = set()
    on_faces: dict[int, Any] = {}
    for cell in _usable_candidates(polys, geometry, usable):
        polygon = polys.geometry.loc[cell]
        piece = geometry.intersection(polygon)
        if areal:
            size = float(piece.area)
            if size <= 0:
                continue
            if (
                not thin_feature
                and _mean_width(piece) < _SLIVER_WIDTH
                and size < _SLIVER_SHARE * float(polygon.area)
            ):
                slivers.add(cell)
        else:
            size = float(piece.length)
            if size <= 0:
                continue
            along = piece.intersection(polygon.boundary)
            if along.length > 0:
                on_faces[cell] = along
        sizes[cell] = size
    # Slivers go -- unless they are all there is, when they ARE the feature: a
    # strip 0.015 wide centred on a face is two 0.0075 halves, and lost both.
    if len(slivers) < len(sizes):
        for cell in slivers:
            del sizes[cell]
    # Pairwise, but only over cells whose piece runs along their own boundary:
    # usually none, so a line off every face pays nothing for this.
    faced = sorted(on_faces)
    for i, first in enumerate(faced):
        for second in faced[i + 1:]:
            shared = float(on_faces[first].intersection(on_faces[second]).length)
            if shared > 0:
                sizes[first] -= shared / 2.0
                sizes[second] -= shared / 2.0
    return [(cell, size / measure) for cell, size in sizes.items() if size > 0]


def _usable_candidates(polys, geometry, usable) -> list[int]:
    """The cells the spatial index offers for a geometry, as ``usable`` allows them.

    Sorted whoever chose them, so that the tie rule rests on this function rather
    than on three callers each remembering it.
    """

    hits = polys.sindex.query(geometry, predicate="intersects")
    cells = [int(cell) for cell in polys.index[hits]]
    return sorted(cells if usable is None else usable(cells))


def miss_reason(polys, geometry, *, target: str = "active cell") -> str:
    """Why a feature reached no usable cell, worded once for every caller.

    ``target`` names what was being looked for, so a search restricted to the
    grid edge says ``"on no active grid-edge cell"``. ``"only borders cells"`` is
    its own reason because the area rule created it (ledger 210): a polygon drawn
    along the outer faces touches the edge cells and covers none of them.
    """

    if geometry is None or geometry.is_empty:
        return "no geometry"
    hits = polys.index[polys.sindex.query(geometry, predicate="intersects")]
    if len(hits) == 0:
        return "outside the grid"
    if geometry.area > 0 and all(
        geometry.intersection(polys.geometry.loc[cell]).area <= 0 for cell in hits
    ):
        return "only borders cells -- a polygon must overlap a cell to reach it"
    return f"on no {target}"


def reconcile_surfaces(vor, df: pd.DataFrame = None, min_sep=0.1, trigger_sep=1, which='bottom'):
    """Adjust stacked surface elevations to preserve layer ordering.

    ``min_sep`` is the separation opened where two surfaces conflict, and
    ``trigger_sep`` the separation below which they count as conflicting. Both accept
    a scalar applied to every contact, or **one value per layer** -- which is what
    lets a stack declare a per-unit minimum thickness rather than one number for the
    whole model. A per-layer sequence has one entry per BOTTOM surface, so it is one
    shorter than the frame's column count (the model top is never moved by
    ``which="bottom"``).

    Corrections cascade downward because each contact is compared against the
    already-adjusted surface above it.
    """
    df = vor.gdf_topbtm.copy() if df is None else df
    if isinstance(df, gpd.GeoDataFrame):
        df = df.drop(columns='geometry').map(lambda x: pd.to_numeric(x, errors='coerce'))
    else:
        df = df.map(lambda x: pd.to_numeric(x, errors='coerce'))

    labels = list(df.columns)
    df = df.loc[:, labels]
    n_bottoms = max(len(labels) - 1, 0)
    seps = _per_contact(min_sep, n_bottoms, 'min_sep')
    triggers = _per_contact(trigger_sep, n_bottoms, 'trigger_sep')
    for i, label in enumerate(labels):
        if i == 0:
            continue
        sep, trigger = seps[i - 1], triggers[i - 1]
        diffs = df.diff(axis=1)
        diff_list = list(diffs[diffs[label] >= -trigger].index)
        if which == 'bottom':
            df.iloc[diff_list, i] = df.iloc[diff_list, (i - 1)] - sep
        elif which == 'top':
            df.iloc[diff_list, (i - 1)] = df.iloc[diff_list, i] + sep
        else:
            raise ValueError(f'which arg {which} is not valid. Must be "top" or "bottom"')

    return df


def _per_contact(value, n: int, name: str) -> list[float]:
    """Normalize a scalar or per-layer separation to one value per bottom contact."""

    if np.isscalar(value) or value is None:
        return [0.1 if value is None else float(value)] * n
    values = [float(v) for v in value]
    if len(values) != n:
        raise ValueError(
            f"{name} has {len(values)} value(s) but the stack has {n} bottom "
            "surface(s); pass a scalar or one value per layer."
        )
    return values


def adjust_cells_by_id(
    vor,
    cell_ids: list,
    adjustment: int | float,
    df: pd.DataFrame = None,
    layer: int = 0,
    reconcile: bool = True,
):
    """
    Adjust elevations for selected cells and optionally reconcile surfaces.
    """
    df = vor.gdf_topbtm.copy() if df is None else df.copy()
    df.loc[cell_ids, layer] = df.loc[cell_ids, layer] + adjustment
    return vor.reconcile_surfaces(df) if reconcile else df


def adjust_top_btm_overlaps(
    vor,
    elev_df=None,
    shp_path: Path = None,
    buffer=1,
    layer_bottom_name=1,
    min_sep=None,
) -> pd.DataFrame:
    """
    Lower selected bottoms so they overlap adjacent cell tops consistently.
    """
    elev_df = vor.reconcile_surfaces(min_sep=min_sep) if elev_df is None else elev_df
    cells_to_adjust = vor.get_vor_cells_as_series(shp_path)
    # get_vor_cells_as_series returns one row per geometry, each holding a LIST of
    # intersected cell ids; flatten to unique individual cells.
    cell_ids: list[int] = []
    for value in cells_to_adjust:
        if isinstance(value, (list, tuple, set, np.ndarray, pd.Series)):
            cell_ids.extend(int(cell) for cell in value)
        else:
            cell_ids.append(int(value))
    cell_ids = sorted(dict.fromkeys(cell_ids))
    new_bottoms = {}

    for cell_id in cell_ids:
        adjacent_cells = vor.find_adjacent_cells(cell_id)
        ja_cell_tops = elev_df[0].loc[adjacent_cells]
        ja_min = ja_cell_tops.min()
        new_bottoms[cell_id] = ja_min - buffer

    new_bottoms = pd.Series(new_bottoms, name=layer_bottom_name)
    # Update on a copy then assign back: chained `df[col].update(...)` operates on
    # a copy under copy-on-write (pandas 3.0) and would silently drop the change.
    column = elev_df[layer_bottom_name].copy()
    column.update(new_bottoms)
    elev_df[layer_bottom_name] = column
    return vor.reconcile_surfaces(df=elev_df)

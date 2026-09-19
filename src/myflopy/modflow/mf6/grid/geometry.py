from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely as shp
from shapely.geometry import LineString, Polygon

from myflopy._logging import get_logger
from myflopy.modflow.utils.datatypes.readers import read_shp_gpkg

logger = get_logger(__name__)


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
    Return the area of the cells a geometry reaches, or of the cells named.

    ``shp_gpkg`` goes through the shared selector, whose answer is one entry per
    FEATURE, each a LIST of cells; :func:`selected_cells` flattens it. Handing
    that list of lists to ``.loc`` asked pandas for the label ``(5,)``, so this
    raised ``KeyError: None of [Index([(5,)], ...)] are in the [index]`` for
    every geometry. A geometry reaching no cell covers no area, so the answer is
    0.0, which is what an empty ``.loc`` already unions to -- and the selector
    has already named the feature in a warning.
    """
    if shp_gpkg is not None:
        cells = selected_cells(vor.get_vor_cells_as_series(shp_gpkg))
    elif cell_list is not None:
        cells = cell_list
    else:
        raise ValueError('You must provide a shp_gpkg or cell_list')
    return vor.gdf_vorPolys.loc[cells].union_all().area


#: The frame ``get_vor_idx_from_geometry`` used to compare in, kept only so its
#: ``crs_latlon`` default still means "unchanged". Cells are found in the grid's
#: own frame now.
_LATLON_FRAME = "EPSG:4326"


def get_vor_idx_from_geometry(
    vor,
    shp_to_query: shp = None,
    gdf_to_query: gpd.GeoDataFrame = None,
    crs: str = None,
    crs_latlon: str = _LATLON_FRAME,
    name_col: str = None,
    predicate: str | None = None,
) -> dict:
    """The cells each geometry in a query layer reaches, keyed by feature.

    ``predicate=None`` -- the default -- selects by MEASURE, through the shared
    :func:`cell_shares`: a polygon reaches the cells it COVERS by area, a line
    those it runs through, and a point the ONE cell holding it (the lowest-numbered
    where it sits on a face or a vertex). Passing a geopandas predicate
    (``'intersects'``, ``'contains'``, ...) keeps the raw spatial-index answer,
    which takes any touch; it is the same escape hatch
    :func:`~myflopy.modflow.mf6.grid.selection.get_vor_cells_as_series` offers, and
    asked the way that selector asks it -- of the CELL first, so ``'contains'``
    means "cells CONTAINING the feature". This asked the index the other way round
    from its own ``_idx`` twin, so one name covered two questions and only a
    symmetric predicate hid it.

    The bare ``'intersects'`` this defaulted to answered about cells a feature only
    borders. On a 4x4 grid of 25 ft cells: a polygon that was exactly cell 5 came
    back as 9 cells, and a point on the vertex shared by 5, 6, 9 and 10 as all
    four.

    **The grid decides where a cell is.** This forced ``crs=DEFAULT_CRS`` and then
    compared in EPSG:4326, so it answered about the wrong place for any grid in
    another frame -- on a UTM grid the cell-5 polygon came back as cells 4, 5, 8
    and 9 -- while a query already in the grid's frame made two reprojections it
    never needed. A query carrying its own CRS is reprojected to the grid's now, as
    the selector does (:func:`to_grid_crs`); ``crs`` names the frame a query
    carrying NONE is read in, defaulting to the grid's, where a naive geometry used
    to raise ``Cannot transform naive geometries``. Passed alongside a query that
    carries its own frame it cannot act, and says so rather than being dropped.

    ``crs_latlon`` was that comparison frame and has no job left -- measuring in
    degrees over a grid in feet is the defect, not a setting. It is still accepted,
    so an existing call keeps working, and REFUSES any other value rather than
    quietly answering about the grid's frame instead of the one asked for. It is
    read as a CRS, not as a string, so ``4326`` and ``"epsg:4326"`` are still the
    default value spelled differently.

    Keys are ``name_col``'s values, or each feature's own index label, taken by
    POSITION: a FILTERED frame works now, where ``gdf.iloc[[2]]`` carries the label
    2, was read as row 2 of a one-row frame, and raised ``IndexError: single
    positional indexer is out-of-bounds`` on every subset. A plain
    ``pandas.DataFrame`` with a ``geometry`` column is coerced, because that
    reached the ``_idx`` twin and worked.

    A feature reaching no cell keeps its key with an empty list -- the shape this
    always returned -- and is named in a warning with its reason; when EVERY
    feature misses, the warning carries the extents and both CRSs
    (:func:`where_they_are`), which is what the historical failure here -- a
    cell-sized polygon in degrees over a grid in feet -- is diagnosed from. The
    selector leaves such a feature out instead; that is the one difference left
    between them.
    """

    if crs_latlon is not None and not _is_latlon(crs_latlon):
        raise ValueError(
            f"crs_latlon={crs_latlon!r} cannot act: cells are found in the grid's own "
            f"frame ({getattr(vor, 'crs', None)}), never in a lat/lon one -- comparing "
            "in EPSG:4326 is what answered about the wrong place. Drop crs_latlon; pass "
            "crs= to name the frame of a query that carries none."
        )
    if shp_to_query is not None:
        gdf_query = read_shp_gpkg(shp_to_query)
    elif gdf_to_query is not None:
        gdf_query = gdf_to_query
    else:
        raise ValueError('Provide shp_to_query or gdf_to_query')

    polys = vor.gdf_vorPolys
    source = _query_geometries(gdf_query)
    source_crs = crs if source.crs is None else source.crs
    geometries = to_grid_crs(source, polys, assume_crs=crs)
    keys = list(gdf_query[name_col]) if name_col is not None else list(geometries.index)

    vor_idx_dict = {}
    missed = []
    for key, geometry in zip(keys, geometries, strict=True):
        cells = _reached_cells(polys, geometry, predicate)
        vor_idx_dict[key] = cells
        if not cells:
            missed.append((key, geometry))
    if missed:
        _report_misses(polys, source, keys, missed, source_crs)
    return vor_idx_dict


def get_vor_idx_from_geometry_idx(
    vor,
    gdf_to_query: gpd.GeoDataFrame = None,
    idx: int = 0,
    predicate: str | None = None,
) -> list[int]:
    """The cells ONE row of a query frame reaches: flat, ascending, each once.

    ``idx`` is positional. The rule is :func:`get_vor_idx_from_geometry`'s -- which
    this asks for that one row, so the pair cannot drift apart -- flattened by
    :func:`selected_cells`.

    This kept a bare ``'intersects'`` of its own: the cell-5 polygon came back as 9
    cells, the point on a vertex as 4, a point on the face between 4 and 5 as BOTH
    of them. It also read the query in the GRID's frame while its sibling
    reprojected to EPSG:2927 and compared in EPSG:4326, so the two disagreed --
    ``[4, 5]`` here against ``[5]`` there for that face point, and for a query in
    any other CRS, EMPTY here (degrees measured against feet) against six cells
    there. Under an explicit ``predicate`` they asked the spatial index opposite
    ways round as well, so ``'contains'`` meant "cells containing the feature"
    there and "cells inside it" here; both now ask it the way the shared selector
    does.
    """

    answer = get_vor_idx_from_geometry(
        vor, gdf_to_query=gdf_to_query.iloc[[idx]], predicate=predicate,
    )
    return selected_cells(list(answer.values()))


def _is_latlon(crs) -> bool:
    """Whether ``crs_latlon`` still names EPSG:4326, however it was spelled.

    Read as a CRS rather than compared as a string: ``4326`` is an int geopandas
    accepts and a ``pyproj.CRS`` for it does not ``str()`` to the literal, so a
    string test would refuse the default value under another spelling. Normalised
    through a geopandas object rather than by importing pyproj, which is
    geopandas' dependency and not one this project declares; an unreadable value
    raises pyproj's own ``CRSError``, which names it better than a guess here
    could.
    """

    return gpd.GeoSeries([], crs=crs).crs == gpd.GeoSeries([], crs=_LATLON_FRAME).crs


def _query_geometries(gdf_query) -> gpd.GeoSeries:
    """The query's geometries as a GeoSeries, whatever frame they arrived in.

    A plain ``pandas.DataFrame`` carrying a ``geometry`` column of shapely objects
    reached ``get_vor_idx_from_geometry_idx`` and worked, so it is coerced rather
    than broken: a bare ``Series`` has no CRS to reconcile and failed with
    ``AttributeError: 'Series' object has no attribute 'set_crs'`` several frames
    away from the call. A query with no geometry at all is named here, the way the
    shared selector names a bad ``overlapping_geometry``.
    """

    geometries = getattr(gdf_query, "geometry", None)
    if geometries is None:
        raise ValueError(
            "the query must carry geometry -- a GeoDataFrame, a GeoSeries, or a frame "
            "with a 'geometry' column of shapely objects -- got a "
            f"{type(gdf_query).__name__} with none."
        )
    if isinstance(geometries, gpd.GeoSeries):
        return geometries
    return gpd.GeoSeries(geometries)


def _reached_cells(polys, geometry, predicate: str | None) -> list[int]:
    """The cells one geometry reaches: by measure, or by a caller's own predicate.

    One spelling for both public selectors, and the SELECTOR's ORIENTATION for the
    predicate branch. ``get_vor_cells_as_series`` joins the cells to the query, so
    a predicate there reads of the cell first: ``'contains'`` means "cells
    CONTAINING the feature". The pair handed the spatial index its two sides
    opposite ways -- one of them asking ``contains(feature, cell)`` -- so one name
    asked two different questions, and only a symmetric predicate hid it.

    Cells come back sorted whatever order the index answered in, because a grid
    frame need not be held in cell order.
    """

    if predicate is None:
        return [cell for cell, _share in cell_shares(polys, geometry)]
    hits = gpd.GeoSeries([geometry]).sindex.query(polys.geometry, predicate=predicate)[0]
    return sorted(int(cell) for cell in polys.index[hits])


def _report_misses(polys, source, keys, missed, source_crs) -> None:
    """Name the features that reached nothing -- and where they are when NONE did.

    Ledger 211's policy for a warning rather than a raise, worded as
    :func:`miss_reason` and :func:`where_they_are` word it everywhere else. The
    extents-and-CRS sentence is added only when every feature missed, which is
    where it is the diagnosis: the historical failure at this site was a CRS
    mismatch, a cell-sized polygon in degrees over a grid in feet. Its bounds are
    the SOURCE frame's own, as every other caller passes them, because a mismatch
    is two boxes and reprojecting one of them first hides it.
    """

    reasons = ', '.join(
        f"{key!r} ({miss_reason(polys, geometry, target='cell')})"
        for key, geometry in missed[:10]
    )
    if len(missed) > 10:
        reasons += f", and {len(missed) - 10} more"
    message = '%d of %d features reach no grid cell and come back empty: %s'
    args = [len(missed), len(keys), reasons]
    if len(missed) == len(keys):
        message += ' %s'
        args.append(where_they_are(source, polys, source_crs=source_crs))
    logger.warning(message, *args)


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


# An areal piece narrower than this, in the grid's length unit, is dropped -- and a
# LINE piece SHORTER than it, the same comparison SFR makes and the one branch that
# had no floor until ledger 237. Almost
# always it is a face-aligned edge whose coordinates were perturbed. A
# perturbation is a DISTANCE, so the test is the piece's own measure -- an areal
# piece's mean width, a line piece's length -- and not its share of
# the cell: a share rule caught a NAD83 round trip on 500 ft cells and missed it on
# 300 ft ones. Measured: 3e-7 after a double-precision reprojection, 3.4e-4 after
# NAD83 <-> NAD83(HARN), up to 5e-3 from vertices FloPy wrote at 9 significant
# digits. It also drops REAL corner clips that thin, deliberately, as
# Spread.min_share does -- on real overlays at most 1.5e-3 ft2, 8.3e-6 of a cell.
# Coordinates ROUNDED on export (GeoJSON at 7 decimal degrees: 0.003 to 0.018)
# straddle it and are not reliably covered. Ledger 210 -- and SFR drops a reach
# SHORTER than it, whatever min_reach_length says (ledger 214).
_SLIVER_WIDTH = 0.01

# ...and a sliver is also under this share of its cell: an areal piece's width, or a
# line piece's own LENGTH, against the cell's mean width. Noise is width / cell width
# of a cell, under 1% on any cell wider than half a unit, so this changes nothing on
# a real grid; without it a lab-scale grid (0.01-unit cells) lost ALL its coverage,
# every piece being narrower -- or, for a line, shorter -- than the tolerance.
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


def to_grid_crs(geometries, polys, *, assume_crs=None):
    """Put a query in the GRID's frame, because the grid decides where a cell is.

    The one copy of the reconcile every selector needs. A query carrying its own
    CRS is reprojected to the grid's, saying so; a query carrying NONE is read as
    the grid's -- or as ``assume_crs`` where a caller offers one, which is the only
    surviving job of ``get_vor_idx_from_geometry``'s ``crs=``. Left in another
    frame a query answers about another place: a cell-sized polygon handed over in
    degrees against a grid in feet reached NO cell at all, and one measured in
    EPSG:2927 over a UTM grid reached cells 4, 5, 8 and 9 where 5 was the whole
    truth.

    ``assume_crs`` over a query that carries its OWN frame cannot act, so it is
    ignored with a warning rather than silently: a parameter that reads as the
    frame a query will be measured in, and is not, is the defect this function
    exists to close.

    Returns a NEW object rather than stamping a frame on the caller's:
    ``geometries.crs = ...`` wrote through the ``.geometry`` view it was handed.
    Neither side georeferenced is left alone -- raw coordinates against raw
    coordinates is what that caller has, and it is not this function's place to
    refuse it.
    """

    grid_crs = getattr(polys, "crs", None)
    source_crs = getattr(geometries, "crs", None)
    if source_crs is None:
        assumed = grid_crs if assume_crs is None else assume_crs
        if assumed is None:
            return geometries
        geometries = geometries.set_crs(assumed, allow_override=True)
        source_crs = geometries.crs
    elif assume_crs is not None and source_crs != assume_crs:
        logger.warning(
            'the geometries carry their own CRS (%s), so %s -- the frame named for '
            'a query that carries none -- is ignored', source_crs, assume_crs,
        )
    if grid_crs is None or source_crs == grid_crs:
        return geometries
    logger.warning(
        'geometries are in %s but the grid is in %s; reprojecting them to '
        'match the grid', source_crs, grid_crs,
    )
    return geometries.to_crs(grid_crs)


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


def centroid_zones(polys, zones) -> np.ndarray:
    """Which of ``zones`` holds each cell's CENTROID, as a row position; -1 for none.

    The centroid rule, in the one place both its callers can read it. A ZONE labels
    a WHOLE cell -- a K zone for PEST, a layer's extent for a surface -- so it is
    decided by the cell's centre, deliberately, and not by area (ledgers 118, 213):
    ``cell_shares`` is the area rule and is NOT what this does.

    A zone is CLOSED: a centroid lying exactly ON its edge is in it. ``within``
    excludes a boundary, and a boundary is drawn along a row of centroids more
    often than anywhere else, since a grid's centroids are where its data is. Two
    zones split at y=37.5 on a 4x4 grid of 25 ft cells -- the line through the
    row-1 centroids -- left that row in NO zone, while a split at y=37.4 gave it
    the upper one.

    Where two zones share that edge, or overlap anywhere, the FIRST ROW wins. That
    is what the callers' docstrings always claimed and what neither delivered: the
    winner came out of the spatial index's own tree order, so on a two-zone frame
    whose rows were reversed the same zone still won.

    ``polys`` is the grid's cell frame (``vor.gdf_vorPolys``) and the answer is
    positional over its rows, matching the per-cell arrays both callers build
    beside it. ``zones`` is a GeoDataFrame, a GeoSeries, or a plain sequence of
    geometries, in the grid's CRS -- reprojecting is the caller's, since only it
    knows what its row labels mean. A row whose geometry is MISSING or empty never
    matches, which a caller may lean on: blanking a row's geometry takes it out of
    the running while keeping every later row's position.
    """

    geometries = np.asarray(
        list(zones.geometry if hasattr(zones, "geometry") else zones), dtype=object
    )
    rows = np.full(len(polys), -1, dtype=int)
    if geometries.size == 0:
        return rows
    # `covers` is evaluated input-to-tree -- zone.covers(centroid) -- so the CELLS go
    # in the tree. The index keeps this off an O(ncpl x nzones) elementwise test, and
    # it never matches a row whose geometry is missing or empty, so neither needs
    # filtering out first.
    zone_hits, cell_hits = polys.geometry.centroid.sindex.query(
        geometries, predicate="covers"
    )
    # The lowest matching row per cell, in one pass: the index answers in tree order,
    # so a tie must not be left to whichever pair came back first.
    winner = np.full(len(polys), geometries.size, dtype=int)
    np.minimum.at(winner, cell_hits, zone_hits)
    found = winner < geometries.size
    rows[found] = winner[found]
    return rows


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
    ``Spread``: a K zone that was exactly one cell painted five. A piece that is
    noise is dropped for the same reason -- an areal piece NARROWER than
    ``_SLIVER_WIDTH``, a line piece SHORTER than it, each also under
    ``_SLIVER_SHARE`` of its cell's area or mean width -- unless every piece of the
    feature is one, when they are the feature. The line floor is the one SFR has
    applied to a reach since ledger 214: without it a line overshooting a face by
    1e-6 wrote a full second record for the cell beyond, at the same elevation and
    conductance, so a ``.gpkg`` DRN and an SFR built from one geometry disagreed.

    Candidates come from the spatial index rather than an elementwise
    ``intersects`` over every cell, which was O(features x ncpl).
    """

    if geometry is None or geometry.is_empty:
        return []
    measure = _extent(geometry)
    if measure <= 0:
        # A point has no extent to divide, so the FEATURE is divided instead: each of
        # its parts carries an equal share of it. Giving every part 1.0 -- what this did
        # -- broke the promise three paragraphs up, that fractions are of the WHOLE
        # feature. Measured: a MultiPoint of 3 wells declared at -1000 wrote -3000 under
        # `mf.Spread`, which exists to conserve a feature's total, and a part 50,000 ft
        # off the grid simply vanished: 2 cells at 1.0 each, summing to 2.0, so nothing
        # downstream could tell a complete feature from a truncated one (ledger 287).
        # The extent-based branch below was always right -- a MultiLineString with one
        # part off-grid already summed to 0.5 -- which is why only points were affected.
        parts = [part for part in shp.get_parts(geometry) if not part.is_empty]
        if not parts:
            return []
        each = 1.0 / len(parts)
        # Accumulated, not a set: two parts landing on ONE cell owe it both their
        # shares, and collapsing them would lose the second.
        sizes: dict[int, float] = {}
        for part in parts:
            hits = _usable_candidates(polys, part, usable)
            if hits:
                sizes[hits[0]] = sizes.get(hits[0], 0.0) + each
        return sorted(sizes.items())

    # A GeometryCollection is measured by its highest dimension: the lines in a
    # polygon collection add no area. In practice they are what gpd.clip and
    # gpd.overlay leave (keep_geom_type=False) where a zone ran along the clip
    # edge, so refusing the collection refused their default output.
    areal = float(getattr(geometry, "area", 0.0) or 0.0) > 0
    # A feature that is itself this thin -- this SHORT, for a line -- is what was
    # drawn, not a sliver of something larger.
    thin_feature = (_mean_width(geometry) if areal else measure) <= _SLIVER_WIDTH
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
            # The floor the areal branch has had since ledger 210, worded for a
            # length: a line ending ON a face that overshoots it by 1e-6 -- a CRS
            # round trip moves an end by ~3e-7 -- gave the cell beyond a piece of
            # 1e-6, and every field not wrapped in ``Spread`` is written whole, so
            # that was a second DRN record at the same elevation and conductance.
            # A length is measured against the CELL's own mean width -- half the
            # side of a square cell -- never against its area or its perimeter, so
            # a sub-unit grid, whose every piece is shorter than the tolerance,
            # keeps the pieces that are a real part of a cell.
            if (
                not thin_feature
                and size < _SLIVER_WIDTH
                and size < _SLIVER_SHARE * _mean_width(polygon)
            ):
                slivers.add(cell)
            along = piece.intersection(polygon.boundary)
            if along.length > 0:
                on_faces[cell] = along
        sizes[cell] = size
    # Slivers go -- unless they are all there is, when they ARE the feature: a
    # strip 0.015 wide centred on a face is two 0.0075 halves, and lost both.
    if len(slivers) < len(sizes):
        for cell in slivers:
            del sizes[cell]
            # ...and out of the face bookkeeping below, which would otherwise ask a
            # dropped cell to hand back half of a stretch it no longer holds --
            # KeyError on a line running a few thousandths along a face before
            # turning into one of the two cells. Clearing the bookkeeping instead
            # leaves an untouched pair sharing a face with the WHOLE stretch each,
            # which is the doubling ledger 212 removed.
            on_faces.pop(cell, None)
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


def selected_cells(selection) -> list[int]:
    """The flat, sorted, unique cells a selector answer names.

    :func:`myflopy.modflow.mf6.grid.selection.get_vor_cells_as_series` answers one
    entry PER FEATURE, each a LIST of cells -- the shape ledger 192 kept, so that
    a caller can still tell which feature reached what. A caller wanting only
    "which cells" has to flatten it, and three of them did not: they handed the
    list of lists straight to ``.loc``, which asked pandas for the label ``(5,)``
    and raised on every input. Written once here rather than at each of the six
    callers that want the flat answer.

    Accepts the Series, a list of lists, or a flat list of cell ids, so a caller
    that already holds ids can pass them through unchanged.
    """

    cells: list[int] = []
    values = selection.tolist() if hasattr(selection, "tolist") else selection
    for value in values:
        if isinstance(value, (list, tuple, set, np.ndarray, pd.Series)):
            cells.extend(int(cell) for cell in np.asarray(value).ravel())
        elif value is not None:
            cells.append(int(value))
    return sorted(dict.fromkeys(cells))


def _bounding_box(obj):
    """A frame's ``total_bounds`` or a single geometry's ``bounds``, as one tuple."""

    bounds = getattr(obj, "total_bounds", None)
    return obj.bounds if bounds is None else bounds


def where_they_are(
    features,
    polys,
    *,
    noun: str = "features",
    source_crs=None,
    grid_crs=None,
    no_crs_subject: str | None = None,
) -> str:
    """Where the features are and where the grid is -- what an empty result is acted on with.

    Ledger 211 refuses a source landing on NO cell, and the two bounding boxes
    plus the two CRSs are what that is diagnosed from: a model placed at the
    origin, a file read in degrees over a State Plane grid, a file carrying no CRS
    at all. The answer is the same wherever the miss was found, so it is worded
    once. ``features`` is anything with ``total_bounds`` -- a GeoDataFrame, a
    GeoSeries -- or a single geometry, which has ``bounds``.

    ``noun`` names the subject, plural, because each caller's differs
    ("features", "targets", "clip bounds"); ``no_crs_subject`` replaces only the
    opening of the no-CRS clause, for a caller whose subject there is the FILE
    rather than what was read out of it. ``grid_crs`` defaults to the polygons'
    own CRS, which is the grid's: a caller reading it off a grid OBJECT gets
    ``None`` whenever that object carries no ``.crs``, and the sentence then said
    "both are in ..." of a file that had in fact been reprojected.
    """

    if grid_crs is None:
        grid_crs = getattr(polys, "crs", None)
    if source_crs is None:
        subject = no_crs_subject or f"the {noun} carry NO CRS, so they"
        crs = f"{subject} were read as the grid's ({grid_crs})"
    elif grid_crs is not None and source_crs != grid_crs:
        crs = f"the {noun} were reprojected from {source_crs} to the grid's {grid_crs}"
    else:
        crs = f"both are in {source_crs}"
    span = ", ".join(f"{v:,.1f}" for v in _bounding_box(features))
    grid = ", ".join(f"{v:,.1f}" for v in _bounding_box(polys))
    return f"The {noun} span ({span}) and the grid ({grid}); {crs}."


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
    # One entry per geometry, each a LIST of cells: `selected_cells` is the one
    # copy of that flattening.
    cell_ids = selected_cells(vor.get_vor_cells_as_series(shp_path))
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

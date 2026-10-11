from __future__ import annotations

import json
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely as shp
from shapely.geometry import LineString, Polygon

from myflopy._logging import get_logger
from myflopy.modflow.utils.datatypes.readers import read_shp_gpkg, read_vector_layer

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
    """Distance from one polygon's CENTROID to the nearest point it shares with another.

    Plain planar geometry on the polygons -- NOT MODFLOW 6's connection length. MF6
    measures a connection (a DISU's ``CL12``, a LAK wall's ``CONNLEN``) from the cell's
    DISV centre -- ``cell2d``'s ``xc``/``yc``, which for a Voronoi cell is its generator
    point, not its centroid -- perpendicular to the LINE through the two vertices the
    cells share. This measures from the centroid to the nearest point of the shared
    boundary, so the two part wherever the centre is not the centroid, or the foot of the
    perpendicular falls off the shared segment. For MF6's number read ``grid.cl12``
    (after ``grid.get_disu_connectivity()``), or
    :meth:`~myflopy.modflow.mf6.grid.connectivity.DisvFaces.length_and_width` for one
    face (ledger 328). Nothing in myflopy calls this.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        The polygons, e.g. ``grid.gdf_vorPolys``. Rows are read by POSITION
        (``iloc``), never by index label.
    poly_idx1 : int
        Position of the polygon measured FROM: the distance starts at its centroid.
    poly_idx2 : int
        Position of the other polygon.

    Returns
    -------
    float
        In the polygons' length units. Swapping the two positions measures from the
        other centroid, so the answer is per side. ``nan`` when the polygons share no
        point; the distance to the corner when they touch only there.

    Raises
    ------
    IndexError
        A position outside ``gdf``.

    See Also
    --------
    myflopy.modflow.mf6.grid.connectivity.DisvFaces.length_and_width : MF6's length for one face.
    myflopy.modflow.mf6.grid.voronoi.VoronoiGridPlus.cl12 : MF6's length for every face.

    Examples
    --------
    >>> import geopandas as gpd
    >>> from shapely.geometry import box
    >>> from myflopy.modflow.mf6.grid.geometry import calculate_distance
    >>> cells = gpd.GeoDataFrame(geometry=[box(0, 0, 1, 1), box(1, 0, 3, 1)])
    >>> calculate_distance(cells, 0, 1), calculate_distance(cells, 1, 0)
    (0.5, 1.0)

    MF6 measures the same face from the DISV centre instead -- here ``(0.25, 0.5)``:

    >>> from myflopy.modflow.mf6.grid.connectivity import DisvFaces
    >>> props = {"vertices": [[0, 0.0, 0.0], [1, 1.0, 0.0], [2, 0.0, 1.0], [3, 1.0, 1.0],
    ...                       [4, 3.0, 0.0], [5, 3.0, 1.0]],
    ...          "cell2d": [[0, 0.25, 0.5, 4, 0, 2, 3, 1], [1, 2.0, 0.5, 4, 1, 3, 5, 4]]}
    >>> DisvFaces(props).length_and_width(0, 1, toward=0)
    (0.75, 1.0)
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

    The union is taken of the REPROJECTED cells. A cell with an edge a few ulps long can
    cross itself once reprojected -- 17 of 156,582 cells at (5e6, 5e6) ft in EPSG:2927 did,
    each valid where it was drawn, and GEOS refused the union, so the grid did not build
    (ledger 316; ``docs/reproducers/voronoi_centroid_reprojection.py``). The tessellation
    no longer writes such edges; should GEOS refuse anyway, the centroid is taken in the
    model's own coordinates and only that point is reprojected.
    """
    try:
        return vor.gdf_latlon.union_all().centroid
    except shp.errors.GEOSException as exc:
        logger.debug("grid centroid taken in model coordinates: the lat/lon cells do not union (%s)", exc)
        centroid = vor.gdf_vorPolys.union_all().centroid
        return gpd.GeoSeries([centroid], crs=vor.gdf_vorPolys.crs).to_crs(vor.crs_latlon).iloc[0]


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

    # Checked before ``.iloc``: a missing frame reached it as ``None`` and died with
    # ``AttributeError: 'NoneType' object has no attribute 'iloc'`` instead of the
    # refusal every other bad query gets.
    _query_geometries(gdf_to_query)
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

    if gdf_query is None:
        raise ValueError(
            "no query was given: pass the features as gdf_to_query -- a GeoDataFrame, "
            "a GeoSeries, or a frame with a 'geometry' column of shapely objects."
        )
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


def set_k_vor(vor, k_dict: dict = None, k_default: float = 100.0) -> list[float]:
    """
    Set hydraulic conductivity values on the Voronoi polygon GeoDataFrame.

    Documented on :meth:`VoronoiGridPlus.set_k_vor`. No mapping gives every cell
    ``k_default``. The column is FLOAT whatever was passed: seeded with the int
    ``100`` it took a float key only through pandas' incompatible-dtype upcast,
    which warns today and is deprecated to raise.
    """
    frame = vor.gdf_vorPolys
    frame["Kh"] = float(k_default)
    for key, cells in ({} if k_dict is None else k_dict).items():
        frame.loc[frame.index.isin(cells), "Kh"] = float(key)
    return frame["Kh"].to_list()


def _lattice_axis(centre: float, spacing: float, size: float) -> np.ndarray:
    """One axis of the square lattice: ``size * spacing`` wide, centred on ``centre``.

    ``centre - size / 2 * spacing`` plus whole multiples of ``spacing`` up to the
    other end (within 1e-9 of a step, so a span the step divides keeps its end).
    Computed as ``centre + spacing * (arange(n) - size / 2)`` rather than
    ``arange(start, stop, spacing)``, whose float stop decides the count by
    rounding, and so that an even ``size`` puts a point exactly on the centre. ``range`` refused a
    float step, and the old ``int()`` of each end truncated toward zero, which
    moved a non-integer end by up to one unit, toward the origin.
    """

    spacing = float(spacing)
    size = float(size)
    if not np.isfinite(spacing) or spacing <= 0:
        raise ValueError(f"spacing must be a finite number above zero, got {spacing!r}.")
    if not np.isfinite(size) or size < 0:
        raise ValueError(f"size must be a finite number of steps, zero or more, got {size!r}.")
    count = int(np.floor(size + 1e-9)) + 1
    return float(centre) + spacing * (np.arange(count) - size / 2)


def generate_grid_around_point(
    center_point: shp.Point, spacing: float, size: float, crs: str
) -> gpd.GeoSeries:
    """
    Generate a GeoSeries of grid points around a center point.

    Documented on :meth:`VoronoiGridPlus.generate_grid_around_point`.
    """
    x_coords = _lattice_axis(center_point.x, spacing, size)
    y_coords = _lattice_axis(center_point.y, spacing, size)
    points = [shp.Point(float(x), float(y)) for x in x_coords for y in y_coords]
    return gpd.GeoSeries(points).set_crs(crs)


def generate_grid_polygons(
    center_point: shp.Point, spacing: float, size: float, gap: float
) -> shp.MultiPolygon:
    """
    Generate square polygons in a grid pattern around a center point.

    Documented on :meth:`VoronoiGridPlus.generate_grid_polygons`.
    """
    x_coords = _lattice_axis(center_point.x, spacing, size)
    y_coords = _lattice_axis(center_point.y, spacing, size)
    spacing = float(spacing)

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


def top_active_layer(domain) -> np.ndarray:
    """The zero-based layer of the highest ACTIVE cell in each column, ``-1`` where none.

    The one statement of "the top active layer", which array-form recharge and ET
    (``irch``/``ievt="top_active"``), the areal RCH/EVT builders, UZF, SFR's default
    ``reach_layer`` and the MODFLOW-USG importer all ask. Active means ``idomain > 0``,
    the rule :func:`domain_active` states: ``0`` is inactive and ``-1`` is a
    **vertical pass-through** cell, removed from the model while flow passes through it,
    so neither can carry a flux.

    Parameters
    ----------
    domain : array_like of int, shape (nlay, ncpl) or (ncpl,)
        An MF6 IDOMAIN -- a numpy array, a nested list, or anything else
        :func:`numpy.asarray` turns into an integer array. A one-dimensional domain
        is one layer, the reading :func:`domain_active` gives it (the same activity
        in every layer), so its answer is ``0`` where a cell is active and ``-1``
        where it is not. Values are compared with ``> 0``; floats are accepted and
        read the same way.

    Returns
    -------
    numpy.ndarray of int, shape (ncpl,)
        Zero-based -- FloPy adds one when it writes MODFLOW's 1-based ``IRCH``/
        ``IEVT`` -- and ``-1`` for a column with no active cell in any layer. The
        ``-1`` is a flag, not a layer: every caller decides what an empty column
        means (array recharge writes layer 0, whose cell is inactive either way;
        the areal builders and UZF skip the column; SFR refuses the reach), so this
        function does not decide it for them.

    Raises
    ------
    ValueError
        If ``domain`` has neither one nor two dimensions -- a structured
        ``(nlay, nrow, ncol)`` IDOMAIN is not a DISV domain; flatten each layer to
        ``(nlay, nrow * ncol)`` first.

    Notes
    -----
    It replaces five copies of the rule with two predicates between them. The
    array-form helpers used ``!= 0``, so a pass-through cell counted as active: the
    column's recharge was written to that cell, which MF6 has removed, and MF6 skips
    such a record without a word. Measured on a 92-column, three-layer model with 40
    pass-through cells in layer 1, those columns received 0.0 of their 152.02 ft3/d
    -- 42.2% of the model's recharge -- under Normal termination (ledger 319,
    ``docs/reproducers/top_active_passthrough.py``).

    Examples
    --------
    >>> import numpy as np
    >>> from myflopy.modflow.mf6.grid.geometry import top_active_layer
    >>> idomain = np.array([[1, -1, 0, 0],
    ...                     [1,  1, 1, 0]])
    >>> top_active_layer(idomain).tolist()   # the -1 column's flux goes to layer 1
    [0, 1, 1, -1]
    >>> top_active_layer([1, 0, 1]).tolist()
    [0, -1, 0]
    """

    values = np.asarray(domain)
    if values.ndim == 1:
        values = values[np.newaxis, :]
    if values.ndim != 2:
        raise ValueError(
            "top_active_layer: expected a (nlay, ncpl) or (ncpl,) IDOMAIN, got shape "
            f"{values.shape}. A structured (nlay, nrow, ncol) domain must be reshaped "
            "to (nlay, nrow * ncol) first."
        )
    active = values > 0
    return np.where(active.any(axis=0), np.argmax(active, axis=0), -1).astype(int)


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


def polygon_cell_values(cell_polys, zones, values, *, rule: str = "majority"):
    """Per-cell values painted from a polygon map; ``(values, reached)``.

    The three polygon-to-cell rules in one place, because they DISAGREE and the
    disagreement is visible in what MF6 is written. The split is ``mean`` against the
    other two: on a 64-cell grid with one geology contact it wrote 12 distinct values
    where the map holds 2, and over 19 contact positions its area-integral ran 0.0% to
    60.5% from ``majority``'s. ``majority`` and ``centroid`` differed on at most ONE
    cell of 64 over the same sweep -- close, but not one rule wearing two names. A
    caller that inherits a rule instead of naming one is choosing by accident
    (ledger 299).

    ``majority`` gives a cell whichever polygon covers most of it. It is the default
    for :class:`~myflopy.sources.Polygons` because it invents nothing: an area-weighted
    ``mean`` across a contact writes a value no mapped unit has. ``mean`` is that
    area-weighted mean, right for a continuous field. ``centroid`` defers to
    :func:`centroid_zones` -- the same rule, not a second implementation of it.

    **The overlay, not** :func:`cell_shares`. That function normalizes by the whole
    FEATURE (share of the polygon falling in each cell, which is what ``mf.Spread``
    needs to divide a total) and applies a sliver rule doing it. Filling an array needs
    share of the CELL -- the same intersection, a different denominator -- and must
    keep the slivers, since a sliver is still the only thing covering its own sliver of
    cell. ``gpd.overlay`` conserves to 5.6e-16 where ``cell_shares`` loses pieces at
    scale (ledger 298).

    Ties go to the LOWEST row of ``zones``, matching :func:`centroid_zones`, so two
    polygons sharing an edge resolve the same way under every rule. Returns ``NaN`` and
    ``False`` for a cell no polygon reaches: what an unreached cell MEANS is the
    caller's, and filling it here would hide it.

    ``cell_polys`` is the grid's cell frame (``vor.gdf_vorPolys``); the answer is
    positional over its rows. ``zones`` is a GeoDataFrame, GeoSeries or sequence of
    geometries in the grid's CRS, and ``values`` the number each carries, positionally.
    """

    import geopandas as gpd

    numbers = np.asarray(values, dtype=float)
    painted = np.full(len(cell_polys), np.nan, dtype=float)
    reached = np.zeros(len(cell_polys), dtype=bool)
    if numbers.size == 0 or len(cell_polys) == 0:
        return painted, reached

    if rule == "centroid":
        rows = centroid_zones(cell_polys, zones)
        found = rows >= 0
        painted[found] = numbers[rows[found]]
        return painted, found

    if rule not in ("majority", "mean"):
        raise ValueError(
            f"polygon_cell_values(rule={rule!r}) is not a rule. Use 'majority', "
            "'mean' or 'centroid'."
        )

    geometries = list(zones.geometry if hasattr(zones, "geometry") else zones)
    # Positions assigned explicitly rather than taken from `reset_index`: the answer is
    # positional over BOTH frames, and a cell frame that has been filtered or reindexed
    # carries an index that is no longer 0..n-1.
    cells = cell_polys[["geometry"]].copy()
    cells["_cell"] = np.arange(len(cells))
    zone_frame = gpd.GeoDataFrame(
        {"_zone": np.arange(len(geometries))},
        geometry=geometries,
        crs=getattr(cell_polys, "crs", None),
    )
    pieces = gpd.overlay(cells, zone_frame, how="intersection", keep_geom_type=True)
    if pieces.empty:
        return painted, reached

    cell_index = pieces["_cell"].to_numpy()
    zone_index = pieces["_zone"].to_numpy()
    area = pieces.geometry.area.to_numpy()

    if rule == "mean":
        weight = np.bincount(cell_index, weights=area, minlength=len(cells))
        total = np.bincount(
            cell_index, weights=area * numbers[zone_index], minlength=len(cells)
        )
        reached = weight > 0
        painted[reached] = total[reached] / weight[reached]
        return painted, reached

    # Majority: per cell, the largest intersection, ties to the lowest zone row.
    #
    # No aggregation step, deliberately. `gpd.overlay` returns ONE row per (cell, zone)
    # and keeps a multi-part intersection WHOLE -- measured: a two-lobe unit covering
    # 60% of a cell comes back as a single MultiPolygon of area 60, not two pieces of
    # 30. So the rows ARE the per-unit totals, and a groupby-sum over them could not
    # change an outcome; a mutation swapping it for `max` passed every test. Ledger
    # 298's rule -- unfalsifiable code is the same defect as a certifying assertion.
    # What holds the behaviour is
    # `test_a_multi_part_polygon_is_judged_by_its_TOTAL_area_in_a_cell`, which fails
    # loudly if that ever stops being true.
    order = np.lexsort((zone_index, -area, cell_index))
    ordered = cell_index[order]
    leading = np.empty(len(order), dtype=bool)
    leading[0] = True
    leading[1:] = ordered[1:] != ordered[:-1]
    winners = order[leading]
    painted[cell_index[winners]] = numbers[zone_index[winners]]
    reached[cell_index[winners]] = True
    return painted, reached


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


def _capitalized(noun: str) -> str:
    """``noun`` with its first letter upper-cased and the rest left as written."""

    return noun[:1].upper() + noun[1:]


def read_footprints(
    source,
    *,
    layer: str | int | None = None,
    crs=None,
    id_field: str | None = None,
    noun: str = "footprint",
) -> gpd.GeoDataFrame:
    """Read polygon footprints -- lakes, excavations -- into one frame in the grid's CRS.

    The one reader behind ``mf.lak(lakes=)``, ``LayerStack.excavate(within=)`` and
    ``LakeTableBuilder(footprint=)``, so they can never read different features, a
    different layer or a different frame out of one file (the excavation design,
    section 3.4; ledger 318 for the table, which reprojects to its RASTER's CRS). It reads, reprojects, labels
    each feature, and refuses what a footprint cannot be: a feature that is not a
    polygon, and one with no area -- a footprint covers cells by AREA, and the shared
    cell rule falls back to LENGTH for a geometry with none, so a collapsed outline
    would come back as the cells along it rather than as the miss it is.

    Parameters
    ----------
    source : str, pathlib.Path, sequence of str or pathlib.Path, or geopandas.GeoDataFrame
        A vector file, several files (concatenated in order; files in different
        CRSs raise ``ValueError: Cannot determine common CRS``), or a frame already
        in memory, which is copied, never modified. Each file's features are labelled
        from its STEM when every file holds one feature (so ``["west.gpkg",
        "east.gpkg"]`` gives ``"west"`` and ``"east"``); a frame has no stem.
    layer : str or int, optional
        The layer to read from each FILE, by name or zero-based index (see
        :func:`~myflopy.modflow.utils.datatypes.readers.read_vector_layer`). ``None``
        reads the first. With a GeoDataFrame it cannot act -- a frame has no layers
        -- and RAISES rather than being ignored.
    crs : pyproj.CRS, str, int or None, optional
        The grid's CRS. Features in another CRS are reprojected to it; features
        carrying NO CRS are taken to be in it already (a model built at the origin
        on purpose). ``None`` leaves the features in their own frame. The frame the
        features arrived in is kept as ``frame.attrs["source_crs"]``, because a miss
        is diagnosed from the two frames and after reprojection only the grid's is
        on hand.
    id_field : str, optional
        A column whose values label the features (converted with ``str``). Without
        it: the file stems when they are unique, else ``"<noun>_<i>"`` by row.
    noun : str, default "footprint"
        What a feature is, singular -- ``"lake"`` for ``mf.lak`` -- used in the
        fallback labels and the refusals.

    Returns
    -------
    geopandas.GeoDataFrame
        One row per feature, indexed by its string label, carrying the source's own
        columns plus ``_source_name`` (the file stem) for features read from a file.
        Its CRS is ``crs`` when given and the features carry one.

    Raises
    ------
    ValueError
        ``layer`` given with a GeoDataFrame; an empty sequence of files; ``id_field``
        not a column; two features with one label; a feature that is not a Polygon
        or MultiPolygon; a feature with no area -- each naming the labels at fault.
        A ``layer`` naming no layer of a file raises through
        :func:`~myflopy.modflow.utils.datatypes.readers.read_vector_layer`, listing
        the layers present.
    pyogrio.errors.DataSourceError
        A file that does not exist or is not a vector file.

    See Also
    --------
    footprint_cells : The cells each footprint covers.
    refuse_footprint_misses : The refusal for footprints that land on no cell.

    Examples
    --------
    >>> import geopandas as gpd, shapely
    >>> from myflopy.modflow.mf6.grid.geometry import read_footprints
    >>> frame = gpd.GeoDataFrame({"name": ["pond"]},
    ...                          geometry=[shapely.Point(0, 0).buffer(10)], crs=2927)
    >>> read_footprints(frame, id_field="name", noun="lake").index.tolist()
    ['pond']
    >>> read_footprints(frame, crs=4326).crs.to_epsg()
    4326
    >>> read_footprints(gpd.GeoDataFrame(geometry=[shapely.LineString([(0, 0), (1, 1)])]),
    ...                 noun="lake")
    Traceback (most recent call last):
    ...
    ValueError: Lake geometries must be Polygon or MultiPolygon.
    """

    title = _capitalized(noun)
    if isinstance(source, gpd.GeoDataFrame):
        if layer is not None:
            raise ValueError(
                f"layer={layer!r} selects a layer of a FILE, and the {noun}s were "
                "given as a GeoDataFrame, which has none -- drop layer=, or pass the path."
            )
        records = source.copy()
    else:
        paths = [source] if isinstance(source, (str, Path)) else list(source)
        frames = []
        for path_value in paths:
            path = Path(path_value)
            frame = read_vector_layer(path, layer)
            frame["_source_name"] = path.stem
            frames.append(frame)
        if not frames:
            raise ValueError(f"At least one {noun} geometry is required.")
        records = gpd.GeoDataFrame(pd.concat(frames, ignore_index=True), crs=frames[0].crs)

    source_crs = records.crs
    if crs is not None and records.crs is not None and records.crs != crs:
        records = records.to_crs(crs)
    if id_field is not None:
        if id_field not in records:
            raise ValueError(f"{title} ID field is missing: {id_field}")
        ids = records[id_field].astype(str)
    elif "_source_name" in records and not records["_source_name"].duplicated().any():
        ids = records["_source_name"].astype(str)
    else:
        ids = pd.Series([f"{noun}_{index}" for index in range(len(records))], index=records.index)
    if ids.duplicated().any():
        duplicates = sorted(ids[ids.duplicated(keep=False)].unique())
        raise ValueError(f"{title} IDs must be unique: {', '.join(duplicates)}")
    if any(geometry.geom_type not in {"Polygon", "MultiPolygon"} for geometry in records.geometry):
        raise ValueError(f"{title} geometries must be Polygon or MultiPolygon.")
    flat = [
        str(label)
        for label, geometry in zip(ids, records.geometry, strict=True)
        if geometry.area <= 0.0
    ]
    if flat:
        raise ValueError(
            f"{title} geometry has no area: {', '.join(flat)}. A {noun} covers cells "
            "by area, so a collapsed or zero-width outline connects to nothing."
        )
    records = records.copy()
    records.index = pd.Index(ids.to_numpy(dtype=object))
    records.attrs["source_crs"] = source_crs
    return records


def footprint_cells(polys, geometries) -> list[list[int]]:
    """The cells each footprint COVERS, by area -- one sorted list per feature.

    The one cell rule for a footprint: a lake's cells (``LAKBuilder.lake_cells``) and
    an excavation's (``LayerStack.excavate``) are this, so the two agree cell for
    cell by construction. It is :func:`cell_shares` per feature with the shares
    dropped -- a footprint covers a cell or it does not -- which makes a polygon
    drawn along cell faces reach only the cells inside it, and drops a piece that is
    a sliver of its cell (a CRS round trip leaves these) unless every piece of the
    feature is one (ledgers 210, 241). The rule is PER FEATURE because of that last
    clause: taken over the union of several features, one feature made wholly of
    slivers would lose them to its neighbours.

    Parameters
    ----------
    polys : geopandas.GeoDataFrame
        The grid's cell polygons (``vor.gdf_vorPolys``), indexed by cell.
    geometries : iterable of shapely geometry
        The footprints -- a GeoSeries, a GeoDataFrame's ``.geometry``, or a list --
        already in the grid's frame (:func:`read_footprints` puts them there).

    Returns
    -------
    list of list of int
        One sorted list of python ``int`` cell numbers per footprint, in the order
        given; an empty list for a footprint covering no cell (this function does
        not refuse it -- :func:`refuse_footprint_misses` does).

    Raises
    ------
    AttributeError
        An item that is not a geometry.

    See Also
    --------
    cell_shares : The shared rule, with the shares.
    refuse_footprint_misses : Refuse or report footprints on no cell.

    Examples
    --------
    >>> import shapely
    >>> from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
    >>> from myflopy.modflow.mf6.grid.geometry import footprint_cells
    >>> vor = irregular_voronoi_grid(nrow=3, ncol=3, cell_size=100.0)
    >>> polys = vor.gdf_vorPolys
    >>> footprint_cells(polys, [polys.geometry.loc[4], shapely.box(-900, -900, -800, -800)])
    [[4], []]
    """

    return [[int(cell) for cell, _share in cell_shares(polys, geometry)] for geometry in geometries]


def refuse_footprint_misses(
    polys,
    footprints,
    cells,
    *,
    noun: str = "footprints",
    consequence: str = "",
    partial: str = "raise",
    source_crs=None,
    grid_crs=None,
) -> list[str]:
    """Refuse footprints covering no grid cell, saying why and where (ledger 211).

    The miss report :func:`footprint_cells` leaves to its callers, worded once. Each
    missing footprint is named with :func:`miss_reason`, and a refusal carries both
    extents and both CRSs (:func:`where_they_are`), because a location mismatch is two
    boxes. Only two of ``miss_reason``'s answers can come back for an areal footprint
    -- "outside the grid" and "only borders cells" -- because nothing filters the
    cells: any overlap by area is kept.

    Parameters
    ----------
    polys : geopandas.GeoDataFrame
        The grid's cell polygons, indexed by cell.
    footprints : geopandas.GeoDataFrame
        The footprints, indexed by label (what :func:`read_footprints` returns).
    cells : sequence of list of int
        :func:`footprint_cells` for ``footprints``, in the same order.
    noun : str, default "footprints"
        What the features are, PLURAL, as the message reads: ``"lakes"``.
    consequence : str, default ""
        Appended to "cover no grid cell", to say what a miss costs -- ``mf.lak``
        passes ``", so they would carry no LAK connection"``.
    partial : {"raise", "warn"}, default "raise"
        ``"raise"``: ANY footprint missing raises -- a lake with no cells carries no
        LAK connection, which MF6 refuses. ``"warn"``: only ALL missing raises, and
        a partial miss WARNS by name (``warnings.warn`` and the logger, ledger 288),
        which is an excavation's rule -- one file can draw ponds over several grids.
    source_crs : pyproj.CRS, str or None, optional
        The CRS the footprints arrived in, before any reprojection
        (``frame.attrs["source_crs"]``); ``None`` says they carried none.
    grid_crs : pyproj.CRS, str or None, optional
        The grid's CRS; ``None`` reads it off ``polys``.

    Returns
    -------
    list of str
        The labels of the footprints that were left out (empty when none were).

    Raises
    ------
    ValueError
        Under ``partial="raise"`` any footprint on no cell; under ``"warn"`` every
        footprint on no cell. Also an unknown ``partial``, or ``cells`` of another
        length than ``footprints``.

    Warns
    -----
    UserWarning
        Under ``partial="warn"``, some but not all footprints on no cell.

    Examples
    --------
    >>> import geopandas as gpd, shapely
    >>> from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
    >>> from myflopy.modflow.mf6.grid.geometry import (
    ...     footprint_cells, read_footprints, refuse_footprint_misses)
    >>> vor = irregular_voronoi_grid(nrow=3, ncol=3, cell_size=100.0)
    >>> ponds = read_footprints(gpd.GeoDataFrame(
    ...     {"name": ["in", "out"]},
    ...     geometry=[vor.gdf_vorPolys.geometry.loc[4], shapely.box(-900, -900, -800, -800)]),
    ...     id_field="name")
    >>> cells = footprint_cells(vor.gdf_vorPolys, ponds.geometry)
    >>> refuse_footprint_misses(vor.gdf_vorPolys, ponds, cells, noun="ponds")
    Traceback (most recent call last):
    ...
    ValueError: 1 of 2 ponds cover no grid cell: 'out' (outside the grid). ...
    """

    if partial not in ("raise", "warn"):
        raise ValueError(f"partial must be 'raise' or 'warn', not {partial!r}.")
    if len(cells) != len(footprints):
        raise ValueError(
            f"{len(cells)} cell lists for {len(footprints)} {noun}; pass footprint_cells "
            "of the same frame."
        )
    missing = [
        (str(label), geometry)
        for label, geometry, found in zip(footprints.index, footprints.geometry, cells, strict=True)
        if not found
    ]
    if not missing:
        return []
    reasons = ", ".join(f"'{label}' ({miss_reason(polys, geometry)})" for label, geometry in missing)
    if partial == "raise" or len(missing) == len(footprints):
        raise ValueError(
            f"{len(missing)} of {len(footprints)} {noun} cover no grid cell{consequence}: "
            f"{reasons}. "
            + where_they_are(
                footprints.geometry,
                polys,
                noun=noun,
                source_crs=source_crs,
                grid_crs=grid_crs,
            )
        )
    message = (
        f"{len(missing)} of {len(footprints)} {noun} cover no grid cell and are left "
        f"out{consequence}: {reasons}."
    )
    logger.warning(message)
    warnings.warn(message, UserWarning, stacklevel=3)
    return [label for label, _geometry in missing]


def reconcile_surfaces(vor, df: pd.DataFrame = None, min_sep=0.1, trigger_sep=1, which='bottom'):
    """Adjust stacked surface elevations to preserve layer ordering.

    ``min_sep`` is the separation SET where two surfaces conflict, and
    ``trigger_sep`` the separation at or below which they count as conflicting -- so
    a contact between ``min_sep`` and ``trigger_sep`` apart is moved CLOSER, which
    thins a layer as well as opening one. Both accept a scalar applied to every
    contact, or **one value per layer** -- which is what lets a stack declare a
    per-unit minimum thickness rather than one number for the whole model. A
    per-layer sequence has one entry per BOTTOM surface, so it is one shorter than
    the frame's column count.

    Rows are taken by POSITION, whatever the frame's index (ledger 377): the
    index is carried to the result untouched. ``which="bottom"`` walks the contacts
    top-down, moving each lower surface against the already-adjusted one above it;
    ``which="top"`` walks them bottom-up, moving each upper surface against the
    already-adjusted one below it. Either way a contact is final once handled, so
    every contact ends exactly ``min_sep`` apart or farther than ``trigger_sep``.
    """
    if which not in ('bottom', 'top'):
        raise ValueError(f'which arg {which} is not valid. Must be "top" or "bottom"')
    df = vor.gdf_topbtm.copy() if df is None else df
    if isinstance(df, gpd.GeoDataFrame):
        df = df.drop(columns='geometry')
    df = df.map(lambda x: pd.to_numeric(x, errors='coerce'))

    values = df.to_numpy(dtype=float, copy=True)
    n_bottoms = max(values.shape[1] - 1, 0)
    seps = _per_contact(min_sep, n_bottoms, 'min_sep')
    triggers = _per_contact(trigger_sep, n_bottoms, 'trigger_sep')
    contacts = range(1, n_bottoms + 1) if which == 'bottom' else range(n_bottoms, 0, -1)
    for i in contacts:
        upper, lower = values[:, i - 1], values[:, i]
        # A NaN on either side compares False, so it is never a conflict.
        conflict = lower - upper >= -triggers[i - 1]
        if which == 'bottom':
            values[conflict, i] = upper[conflict] - seps[i - 1]
        else:
            values[conflict, i - 1] = lower[conflict] + seps[i - 1]

    return pd.DataFrame(values, index=df.index, columns=df.columns)


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
    trigger_sep=1,
) -> pd.DataFrame:
    """
    Lower selected bottoms so they overlap adjacent cell tops consistently.

    Both reconciles -- the one building a default ``elev_df`` and the final one --
    use the caller's ``min_sep`` and ``trigger_sep``, and a frame passed in is
    copied, never written (ledger 377).
    """
    elev_df = (
        vor.reconcile_surfaces(min_sep=min_sep, trigger_sep=trigger_sep)
        if elev_df is None
        else elev_df.copy()
    )
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
    return vor.reconcile_surfaces(df=elev_df, min_sep=min_sep, trigger_sep=trigger_sep)


#: Below this FRACTION of a cell's own area, an uncovered sliver is arithmetic noise
#: rather than a gap. Polygon intersection is exact to ~1e-15 relative, so on a mesh
#: that nests perfectly the residual measured 4e-15 of a cell -- about 1e-11 ft2 on a
#: 3,800 ft2 cell. Relative rather than absolute so it means the same on a 4 ft2 cell
#: and a 40,000 ft2 one, and five orders of magnitude above the noise so it can only
#: ever hide something already invisible (ledger 298).
_COVERAGE_NOISE = 1e-9


@dataclass(frozen=True)
class CoverageReport:
    """How much of a NEW mesh has no data under it on an OLD one.

    Reported as AREA IN THE MODEL'S OWN UNITS (user decision, 2026-09-22), not as a
    coverage fraction. A fraction needs a tolerance before it means anything -- on a
    perfectly nesting 237 -> 2,352 refinement, comparing ``coverage < 1`` flags 765
    cells (32.5%) that are short by 4e-15 of themselves, and zero at 1e-12. Area
    dissolves that: the same residual is ~1e-11 ft2, which reads as the nothing it is,
    and a threshold a modeller might actually want is then expressible in ft2 rather
    than as an epsilon.
    """

    #: Total area of the new mesh, in the model's units squared.
    total_area: float
    #: Area of the new mesh with no old cell beneath it.
    uncovered_area: float
    #: Per new cell, the area not covered -- only cells above the noise floor.
    uncovered: dict[int, float]
    #: True when the two meshes share no area at all, which is a FRAME mismatch
    #: rather than a coverage gap: `grid_pin` is connectivity-only, so the same mesh
    #: in another frame pins IDENTICALLY while overlapping nothing.
    disjoint: bool
    #: Both bounding boxes and both CRSs, worded once by `where_they_are`.
    where: str

    @property
    def fraction_uncovered(self) -> float:
        """Uncovered area as a fraction of the new mesh, for a one-number summary."""

        return self.uncovered_area / self.total_area if self.total_area else 0.0

    def summary(self) -> str:
        """One sentence naming the area, the cells, and the worst of them."""

        if self.disjoint:
            return f"the two meshes share NO area -- {self.where}"
        if not self.uncovered:
            return f"every cell of the new mesh has data under it ({self.total_area:,.0f} covered)"
        worst = max(self.uncovered, key=self.uncovered.get)
        return (
            f"{self.uncovered_area:,.1f} of {self.total_area:,.1f} "
            f"({100 * self.fraction_uncovered:.2f}%) of the new mesh has no old cell "
            f"under it, across {len(self.uncovered)} cells; the worst is cell {worst} "
            f"at {self.uncovered[worst]:,.1f}"
        )


#: How several OLD cells combine into one NEW cell. See :meth:`CellOverlay.reduce`.
REMAP_RULES: tuple[str, ...] = ("mean", "geometric", "harmonic", "majority")


@dataclass(frozen=True)
class CellOverlay:
    """The ``(old cell, new cell, intersection area)`` triples between two meshes.

    Computed ONCE and reused. Both :func:`cell_coverage` and every remapped array of
    every layer ride on one of these: the geometry does not change between layers or
    between parameters, and re-overlaying per call measured a **24x** tax (0.661 s for
    six array parameters across four layers, against 0.027 s for one shared overlay).

    Positions, not index labels. ``old_index`` and ``new_index`` are positional over the
    two frames' ROWS. ``cell_coverage`` used to bincount by index LABEL and read the
    result by POSITION, which is correct only for a ``RangeIndex``; a row-permuted new
    frame reported 63 uncovered cells where the truth was 0, and a sliced one raised a
    broadcast error. Unreachable in practice because ``gdf_vorPolys`` is always a
    ``RangeIndex`` -- so it was a latent bug, and sharing this makes it reachable.
    """

    old_index: np.ndarray
    new_index: np.ndarray
    area: np.ndarray
    n_old: int
    n_new: int

    @property
    def disjoint(self) -> bool:
        """Whether the two meshes intersect in nothing at all."""

        return self.area.size == 0

    def reduce(self, values, *, rule: str = "mean", active=None):
        """Carry ``values`` from the old mesh onto the new one; ``(values, reached)``.

        INTENSIVE only -- a conductivity, a head, a storage coefficient: a quantity that
        does not scale with the size of the cell holding it. **The denominator is the
        new cell's COVERED area, not its own area**, which is what makes a partially
        covered cell come back with the value of what IS under it rather than that value
        diluted toward zero. Measured on a constant 7.0 field over a mesh grown past its
        source: the covered-area denominator returns 7.000000 for every partially
        covered cell, and the cell's own area returns 0.0066 to 6.9878.

        ``active`` is a boolean mask over the OLD cells, and it masks the WEIGHTS rather
        than the values. Passing an inactive old cell's value through as ``NaN`` poisons
        the new cells above it and ``reached`` does not say so -- measured on a pinching
        wedge, ``rule="mean"`` left ten reached new cells ``NaN``, THREE of them active
        on the new mesh. The active set genuinely differs between meshes, so this is a
        per-LAYER call, not a one-off.

        ``reached`` is ``False`` for a new cell with NO piece under it, and that is the
        FILL predicate. It is a different set from ``CoverageReport.uncovered``, which
        names cells that are merely SHORT -- one measured at 98.9% covered. Filling from
        that set threw away 15 correct area-weighted values out of 51 named.

        **An EXTENSIVE quantity is not this function.** A rate or a total needs the OLD
        cell's area as the denominator (the share of each old cell going out, rather than
        the share of each new cell coming in), and the two are not scalings of one
        another. C6 deliberately does not serve them -- see ledger 300.
        """

        if rule not in REMAP_RULES:
            raise ValueError(
                f"CellOverlay.reduce(rule={rule!r}) is not a rule. Use one of "
                + ", ".join(repr(name) for name in REMAP_RULES)
                + "."
            )
        source = np.asarray(values, dtype=float)
        if source.shape != (self.n_old,):
            raise ValueError(
                f"reduce() takes one value per OLD cell -- {self.n_old} of them -- and "
                f"was given an array of shape {source.shape}. Reduce one layer at a "
                "time; an (nlay, ncpl) array has to be split by the caller, which is "
                "also where the per-layer active mask belongs."
            )

        keep = np.ones(self.area.size, dtype=bool)
        if active is not None:
            mask = np.asarray(active, dtype=bool)
            if mask.shape != (self.n_old,):
                raise ValueError(
                    f"active is a mask over the {self.n_old} OLD cells, not "
                    f"{mask.shape}."
                )
            keep = mask[self.old_index]
        old_at = self.old_index[keep]
        new_at = self.new_index[keep]
        weight = self.area[keep]

        painted = np.full(self.n_new, np.nan, dtype=float)
        reached = np.zeros(self.n_new, dtype=bool)
        if weight.size == 0:
            return painted, reached
        total = np.bincount(new_at, weights=weight, minlength=self.n_new)
        reached = total > 0

        if rule == "majority":
            # Largest contributing area wins, ties to the lowest OLD cell -- the same
            # tie rule `polygon_cell_values` and `centroid_zones` use, so the three
            # agree about a shared edge.
            order = np.lexsort((old_at, -weight, new_at))
            leading = np.empty(order.size, dtype=bool)
            leading[0] = True
            leading[1:] = new_at[order][1:] != new_at[order][:-1]
            winners = order[leading]
            painted[new_at[winners]] = source[old_at[winners]]
            return painted, reached

        taken = source[old_at]
        if rule == "mean":
            weighted = np.bincount(new_at, weights=weight * taken, minlength=self.n_new)
            painted[reached] = weighted[reached] / total[reached]
            return painted, reached

        # Geometric and harmonic are undefined at or below zero -- a log and a
        # reciprocal. Refused by name rather than returning a NaN or an inf that MF6
        # would accept and report success on (ledger 247's failure mode). A HEAD can
        # legitimately be negative, which is why these two are offered for conductivity
        # and the elevation-like arrays use the mean.
        if not np.all(taken > 0.0):
            bad = int(np.count_nonzero(~(taken > 0.0)))
            raise ValueError(
                f"rule={rule!r} needs every contributing value to be positive -- it "
                f"takes a {'logarithm' if rule == 'geometric' else 'reciprocal'} of "
                f"each -- and {bad} of {taken.size} contributions are zero or negative "
                f"(smallest {float(taken.min()):g}). Use rule='mean' for a quantity "
                "that can be zero or negative, such as a head or an elevation."
            )
        if rule == "geometric":
            logs = np.bincount(new_at, weights=weight * np.log(taken), minlength=self.n_new)
            painted[reached] = np.exp(logs[reached] / total[reached])
            return painted, reached
        inverses = np.bincount(new_at, weights=weight / taken, minlength=self.n_new)
        painted[reached] = total[reached] / inverses[reached]
        return painted, reached

    def rule_spread(self, values, *, active=None):
        """How far apart the candidate rules are, per new cell: ``(ratio, worst_cell)``.

        The number that makes a conductivity default honest. For K the CORRECT average
        depends on the FLOW DIRECTION -- a property of the solution, not of the field --
        so no default and no per-field declaration can be right everywhere. Measured in
        MF6 on one 2:1 coarsening: with the contrast ACROSS the flow the arithmetic mean
        gave 264x the fine-mesh boundary flux and the harmonic mean was nearly exact;
        with the SAME contrast rotated 90 degrees the arithmetic mean was exact and the
        harmonic wrong by about 250x. The ORDERING REVERSES.

        So instead of choosing, this reports arithmetic/harmonic per cell, and the
        caller warns where they diverge. Where the ratio is ~1 no choice matters; where
        it is large the reader is told which cells to look at.
        """

        high, _ = self.reduce(values, rule="mean", active=active)
        low, _ = self.reduce(values, rule="harmonic", active=active)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(low > 0, high / low, np.nan)
        if not np.any(np.isfinite(ratio)):
            return np.array([]), None
        worst = int(np.nanargmax(ratio))
        return ratio, worst


def cell_overlay(old_polys, new_polys) -> CellOverlay:
    """Intersect two cell frames and return the pieces as positional triples.

    One ``gpd.overlay(how="intersection", keep_geom_type=True)``, which is what
    :func:`cell_coverage` already used and what a remap needs -- the same call, with the
    OLD cell kept rather than discarded. ``cell_coverage`` selected ``_old`` into the
    overlay and then never read it, so these triples existed for one statement and were
    thrown away with the local.

    NOT :func:`cell_shares`. That normalizes by the whole FEATURE, applies a sliver rule
    and loops in Python: used as a remap it measured 22 of 2,500 old cells losing more
    than 1e-9 of themselves and cost 2x. Note the denominator objection is only half the
    reason -- when the "feature" IS an old cell polygon, its per-feature share is
    exactly the EXTENSIVE weight; the sliver rule and the cost are the rest.
    """

    old_frame = old_polys[["geometry"]].copy()
    old_frame["_old"] = np.arange(len(old_frame))
    new_frame = new_polys[["geometry"]].copy()
    new_frame["_new"] = np.arange(len(new_frame))
    pieces = gpd.overlay(old_frame, new_frame, how="intersection", keep_geom_type=True)
    if pieces.empty:
        empty = np.array([], dtype=int)
        return CellOverlay(empty, empty, np.array([], dtype=float), len(old_frame), len(new_frame))
    return CellOverlay(
        old_index=pieces["_old"].to_numpy(),
        new_index=pieces["_new"].to_numpy(),
        area=pieces.geometry.area.to_numpy(),
        n_old=len(old_frame),
        n_new=len(new_frame),
    )


def refuse_a_different_frame(old_polys, new_polys, *, verb: str) -> None:
    """Refuse two meshes that are not in the same frame, BEFORE anything is carried.

    **The CRS test comes first and it is not optional.** ``CoverageReport.disjoint``
    catches only the zero-pieces case -- two meshes that overlap in nothing, which is
    what a translation produces. It does NOT catch a mesh in the wrong frame that still
    overlaps: measured, the same 1,200 ft footprint relabelled EPSG:2927 -> EPSG:2926
    comes back ``disjoint=False, uncovered_area=0.0`` and ``summary()`` says "every cell
    of the new mesh has data under it". ``gpd.overlay`` warns and proceeds. So neither
    test is sufficient alone and both are here.

    That matters most under a nearest-neighbour fill, where the mirror image of the same
    error fills most of the new mesh from one edge of the old one.

    **THE BLIND SPOT, stated out loud because no tolerance can see it.** A mesh rescaled
    into other UNITS while keeping the same declared CRS passes both tests: measured, the
    same mesh with its coordinates divided by 3.28084 and still labelled EPSG:2927 sits
    entirely inside the old footprint, every new cell is covered, and the remap is
    geometrically self-consistent -- it simply describes a tenth of the domain the
    modeller meant. Nothing geometric distinguishes that from a legitimate SUBDOMAIN
    refinement, which is a real and supported thing to do, so a test that fired on it
    would fire on correct work. This is the same blind spot ledger 211 records for a
    ``.gpkg`` source ("a unit/scale CRS error leaves near features ON the grid in the
    WRONG cells... no tolerance can see that"); the coverage report's own numbers are
    what a reader has to check.

    A grid is refused rather than reprojected, unlike every ``.gpkg`` SOURCE (which does
    ``to_crs`` silently). A source file legitimately arrives in its own frame; two model
    meshes are both the modeller's own, and one of them being in another frame is far
    more likely a declaration error than an intention. Ledger 298 took the same
    direction for a model with no grid: refuse rather than answer a question that cannot
    be asked.
    """

    old_crs = getattr(old_polys, "crs", None)
    new_crs = getattr(new_polys, "crs", None)
    if old_crs is not None and new_crs is not None and old_crs != new_crs:
        raise ValueError(
            f"{verb}: the two meshes are in different coordinate systems -- the old one "
            f"in {old_crs.to_string() if hasattr(old_crs, 'to_string') else old_crs} and "
            f"the new one in "
            f"{new_crs.to_string() if hasattr(new_crs, 'to_string') else new_crs}. "
            "Nothing can be carried between them until they agree: an overlay across two "
            "frames still finds pieces, so this would report full coverage while every "
            "value landed in the wrong cell. Reproject the new grid before moving to it."
        )
    overlay = cell_overlay(old_polys, new_polys)
    if overlay.disjoint:
        raise ValueError(
            f"{verb}: the two meshes share NO area, so there is nothing to carry. "
            + where_they_are(
                old_polys, new_polys, noun="old mesh cells",
                source_crs=old_crs, grid_crs=new_crs,
            )
            + " A mesh in the same frame but another PLACE pins identically under "
            "`grid_pin`, which hashes connectivity only -- so this is exactly the case a "
            "nearest-cell fill would paper over, filling every cell from a source a long "
            "way away."
        )


def cell_coverage(old_polys, new_polys) -> CoverageReport:
    """How much of ``new_polys`` has no ``old_polys`` beneath it.

    The measurement a grid move should be read BEFORE it happens: what a remap would
    have no data for. It answers about GEOMETRY only and knows nothing about packages,
    so it is the same number whatever is being moved.

    ``gpd.overlay(..., how="intersection")`` rather than :func:`cell_shares`: measured
    on a 237 -> 2,352 refinement the overlay conserves to 5.6e-16 while `cell_shares`,
    which exists to divide ONE feature among cells and applies a sliver rule doing it,
    loses pieces at scale. Both are right for their own job.

    **A disjoint pair is reported, not returned as "100% uncovered".** Two meshes in
    different frames overlap in nothing, and that is a different fact from a mesh grown
    past its data -- `grid_pin` cannot tell them apart because it hashes connectivity
    only, so the same mesh shifted to a local origin pins IDENTICALLY (ledger 298).
    """


    new_area = new_polys.geometry.area.to_numpy()
    total = float(new_area.sum())
    # `where_they_are` supplies "The {noun} span ..." itself, so the noun carries no
    # article -- and it words the two boxes and the two CRSs the same way every other
    # miss in this library does (ledger 211), which is the point of borrowing it.
    where = where_they_are(
        old_polys, new_polys, noun="old mesh cells",
        source_crs=getattr(old_polys, "crs", None), grid_crs=getattr(new_polys, "crs", None),
    )
    overlay = cell_overlay(old_polys, new_polys)
    if overlay.disjoint:
        return CoverageReport(total, total, {}, disjoint=True, where=where)

    covered = np.bincount(
        overlay.new_index, weights=overlay.area, minlength=len(new_polys)
    )
    # NOT clipped at zero. An intersection can overshoot its own cell by ~1e-15 of it,
    # so `short` can be negative -- and the threshold below excludes a negative for the
    # same reason it excludes a tiny positive. A `np.maximum(..., 0.0)` here read as
    # careful and was unfalsifiable: a mutation removing it changed no test, because it
    # could not change an outcome (ledger 298).
    short = new_area - covered
    uncovered = {
        int(cell): float(short[position])
        for position, cell in enumerate(new_polys.index)
        if short[position] > _COVERAGE_NOISE * new_area[position]
    }
    return CoverageReport(
        total_area=total,
        uncovered_area=float(sum(uncovered.values())),
        uncovered=uncovered,
        disjoint=False,
        where=where,
    )

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely as shp

from myflopy._deprecation import warn_deprecated
from myflopy._logging import get_logger
from myflopy.modflow.mf6.grid.geometry import (
    cell_shares,
    miss_reason,
    selected_cells,
    to_grid_crs,
    where_they_are,
)
from myflopy.modflow.utils.datatypes.readers import read_shp_gpkg, read_vector_layer

logger = get_logger(__name__)

#: A column name no file carries, so a caller that names no field keeps the
#: feature's own index as its key.
_NO_NAME_FIELD = "__myflopy_unnamed__"


def _measured_cells(gdf_vor_polys, geometries, names) -> pd.Series:
    """The cells each feature actually reaches, by the shared coverage rule.

    A polygon reaches the cells it COVERS by area, a line those it runs through,
    and a point the one cell holding it, with the tie and sliver rules of ledgers
    210-212 and 237. The spatial join this replaced took any touch, so a polygon drawn on
    cell faces selected its neighbours as well: 9 cells for a polygon that was
    exactly cell 5, and 4 for a point on a vertex.

    A feature that reaches nothing is left out, as it always was, but it is named
    in a warning rather than silently absent. This does not RAISE: an empty
    selection is documented as a warning where it is drawn (``map(select=...)``),
    and the callers that build model input raise for themselves.
    """

    cells = {}
    missed = []
    for index, geometry in geometries.items():
        hits = [cell for cell, _share in cell_shares(gdf_vor_polys, geometry)]
        if hits:
            cells[index] = hits
        else:
            missed.append((index, geometry))
    if missed:
        labelled = not isinstance(names, list)
        reasons = ", ".join(
            f"{(names.loc[index] if labelled else index)!r} "
            f"({miss_reason(gdf_vor_polys, geometry, target='cell')})"
            for index, geometry in missed[:10]
        )
        if len(missed) > 10:
            reasons += f", and {len(missed) - 10} more"
        logger.warning(
            "%d of %d features reach no grid cell and are left out: %s",
            len(missed), len(geometries), reasons,
        )
    return pd.Series(cells, name="cells", dtype=object)


def get_vor_cells_as_series(
    gdf_vor_polys,
    overlapping_geometry: shp.Polygon | shp.Point | gpd.GeoSeries | Path = None,
    predicate: str | None = None,
    return_dict: bool = False,
    name_field: str = 'ExploName',
    layer: str | int | None = None,
) -> pd.Series | dict:
    """
    Identify the Voronoi cells each geometry reaches.

    ``predicate=None`` -- the default -- selects by MEASURE: a polygon reaches the
    cells it covers by area, a line those it runs through, and a point the one
    cell holding it, under the tie and sliver rules of ledgers 210-212 and 237. Passing a
    geopandas predicate (``'intersects'``, ``'contains'``, ...) keeps the raw
    spatial join, which takes any touch, including a shared face or a single
    vertex; it is the escape hatch for a caller that wants topology rather than
    coverage.

    ``layer`` names which layer to read when ``overlapping_geometry`` is a path to
    a multi-layer file (a GeoPackage); ``None`` keeps GDAL's default of the first
    one. It is refused for an in-memory geometry, which has no layers to choose
    between -- a parameter that cannot act is as much a defect as a missing one.
    """
    names = None
    if layer is not None and not isinstance(overlapping_geometry, (str, Path)):
        raise ValueError(
            "layer= selects a layer inside a vector FILE, but overlapping_geometry is "
            f"a {type(overlapping_geometry).__name__}. Drop layer=, or pass the file path."
        )

    # Any shapely geometry, tested on the base class rather than a list of
    # concrete types -- the list left out MultiLineString, so a multi-part fault
    # trace fell through to the error branch below.
    if isinstance(overlapping_geometry, shp.geometry.base.BaseGeometry):
        geometries = gpd.GeoSeries([overlapping_geometry])
    elif isinstance(overlapping_geometry, gpd.GeoDataFrame):
        geometries = overlapping_geometry.geometry
        if name_field in overlapping_geometry.columns:
            names = overlapping_geometry[name_field]
    elif isinstance(overlapping_geometry, gpd.GeoSeries):
        # Separately from GeoDataFrame: a GeoSeries has no `.columns`, so the
        # shared branch raised AttributeError on every bare series.
        geometries = overlapping_geometry
    elif isinstance(overlapping_geometry, (str, Path)):
        geoms = read_vector_layer(overlapping_geometry, layer)
        geometries = geoms.geometry
        if name_field in geoms.columns:
            names = geoms[name_field]
    else:
        # RAISE, not return. This returned the exception as DATA, so every caller
        # got a ValueError object where it expected a Series and failed one line
        # later with something unrelated -- `TypeError: 'ValueError' object is not
        # iterable` being the usual disguise.
        raise ValueError(
            "overlapping_geometry must be a shapely geometry, a GeoDataFrame, a "
            f"GeoSeries, or a path to a vector file -- got {type(overlapping_geometry).__name__}."
        )

    if names is None:
        names = list(range(len(geometries)))
    assert len(names) == len(geometries), 'length of names must match length of geometries'

    # One copy of the reconcile, shared with the two get_vor_idx_* selectors, which
    # forced a CRS of their own and so answered about the wrong place.
    geometries = to_grid_crs(geometries, gdf_vor_polys)

    if predicate is None:
        group_join = _measured_cells(gdf_vor_polys, geometries, names)
    else:
        joined = gpd.sjoin(gdf_vor_polys, geometries.to_frame('geometry'), predicate=predicate, how='inner')
        group_join = joined.groupby('index_right').apply(lambda frame: frame.index.to_list(), include_groups=False)
        group_join.name = 'cells'
    if not isinstance(names, list):
        names_idx = names.loc[group_join.index]
        group_join.index = names_idx

    if return_dict:
        return group_join.to_dict()
    return group_join


def get_vor_cells_as_dict(
    vor,
    locs: Path,
    crs: str = None,
    predicate: str | None = None,
    loc_name_field: str = None,
    return_gdf: bool = False,
    layer: str | int | None = None,
) -> dict:
    """
    Provide a shapefile (locs) and get back the voronoi cells that contain them.

    ``layer`` names which layer to read out of a multi-layer file (a GeoPackage);
    ``None`` keeps GDAL's default of the first one, which is what this did before
    the parameter existed.
    """
    crs = vor.crs if crs is None else crs
    gdf_locs = read_vector_layer(locs, layer).to_crs(crs)

    # One call for the whole frame. Per feature, each miss was reported twice and
    # the join's overhead was paid once each; the selector now names the misses.
    cells = get_vor_cells_as_series(
        vor.gdf_vorPolys,
        gdf_locs,
        predicate,
        name_field=_NO_NAME_FIELD if loc_name_field is None else loc_name_field,
    )
    # The {name: [[cells]]} shape ledger 192 left alone.
    loc_vor_cell_dict = {key: [value] for key, value in cells.items()}

    if return_gdf:
        return loc_vor_cell_dict, gdf_locs
    return loc_vor_cell_dict


def get_model_boundary_polygons(gdf_vor_polys) -> dict:
    """
    Return a dict of polygons that form the model domain boundary.
    """
    warn_deprecated(
        "myflopy.modflow.mf6.grid.selection.get_model_boundary_polygons",
        "get_grid_edge",
        since="0.1",
    )
    grid = shp.MultiPolygon(gdf_vor_polys.geometry.to_list())
    polygons = gdf_vor_polys.geometry.to_list()
    convex_hull_boundary = grid.convex_hull.boundary

    boundary_polygons = []
    boundary_polygons_idx = []
    for idx, polygon in enumerate(polygons):
        if convex_hull_boundary.dwithin(polygon, 0.01):
            boundary_polygons.append(polygon)
            boundary_polygons_idx.append(idx)
    return dict(zip(boundary_polygons_idx, boundary_polygons, strict=False))


def _inactive_cells(vor, idomain, idomain_path) -> set[int]:
    """The cells to treat as inactive, from a list of IDS or from a vector file.

    An explicit argument wins, and ``idomain=[]`` means nothing is inactive: both
    used to fall through to the grid's own ``idomain_path`` and quietly answer
    about that instead.
    """

    if idomain is not None and idomain_path is not None:
        raise ValueError(
            "Give idomain= (the INACTIVE cell ids) or idomain_path= (a vector file "
            "covering them), not both."
        )
    ncpl = len(vor.adjacent_cells_idx)
    if idomain is None and idomain_path is None:
        idomain_path = getattr(vor, "idomain_path", None)
    if idomain_path is not None:
        return _inactive_from_file(vor, Path(idomain_path))
    if idomain is None:
        return set()

    values = np.asarray(list(idomain))
    if values.dtype == bool:
        raise ValueError(
            "idomain= is a list of INACTIVE cell ids, not a mask. Pass "
            "np.flatnonzero(domain <= 0)."
        )
    if values.size == ncpl and set(np.unique(values)).issubset({-1, 0, 1}):
        # The trap this parameter sets: MF6's idomain is an ARRAY of flags, and
        # this is neither that nor a list of ACTIVE ids.
        raise ValueError(
            f"idomain= is a list of INACTIVE cell ids, but this is {ncpl} values of "
            "-1/0/1, which reads as an MF6 idomain array. Pass "
            "np.flatnonzero(domain <= 0)."
        )
    cells = {int(value) for value in values.reshape(-1)}
    unknown = sorted(cell for cell in cells if not 0 <= cell < ncpl)
    if unknown:
        raise ValueError(
            f"idomain= names cells this grid does not have (it has {ncpl}, numbered "
            f"from 0): {unknown}."
        )
    return cells


def _inactive_from_file(vor, path: Path) -> set[int]:
    """The cells an idomain file's features COVER, by area (ledger 210).

    Read as a frame, so the file's own CRS reaches the selector and is reprojected
    to the grid's. The union this replaced dropped the CRS, and a file in another
    one then raised ``KeyError: 0`` rather than saying where its features were.
    A file with a layer column is read whole: the edge unions every layer.
    """

    features = read_shp_gpkg(path)
    series = get_vor_cells_as_series(vor.gdf_vorPolys, features)
    cells = set(selected_cells(series))
    if not cells:
        grid = vor.gdf_vorPolys
        raise ValueError(
            f"No feature in '{path.name}' covers a cell of this grid, so nothing "
            "would be inactive. "
            + where_they_are(
                features,
                grid,
                source_crs=features.crs,
                no_crs_subject="the file has NO CRS, so its coordinates",
            )
        )
    return cells


def _active_components(adjacency, active, inactive) -> list[set[int]]:
    """The face-connected groups of active cells."""

    seen: set[int] = set()
    components = []
    for seed in active:
        if seed in seen:
            continue
        members = {seed}
        seen.add(seed)
        stack = [seed]
        while stack:
            cell = stack.pop()
            for neighbour in adjacency[cell]:
                if neighbour not in inactive and neighbour not in members:
                    members.add(neighbour)
                    seen.add(neighbour)
                    stack.append(neighbour)
        components.append(members)
    return components


def _outer_edge_only(adjacency, outside, active, inactive) -> list[int]:
    """Each active component's OUTER edge, dropping the cells that ring a hole.

    A hole is what the outside cannot reach: per component, flood the rest of the
    mesh inward from the cells bordering the outside, and a member whose only
    inactive neighbours were never reached is ringing a hole of its own component.
    An island inside a ring keeps its edge, because that edge is the island's own
    outer boundary.
    """

    edge: list[int] = []
    for members in _active_components(adjacency, active, inactive):
        reached = {cell for cell in range(len(adjacency)) if cell not in members and outside[cell]}
        stack = list(reached)
        while stack:
            cell = stack.pop()
            for neighbour in adjacency[cell]:
                if neighbour not in members and neighbour not in reached:
                    reached.add(neighbour)
                    stack.append(neighbour)
        edge += [
            cell for cell in members
            if outside[cell] or any(neighbour in reached for neighbour in adjacency[cell])
        ]
    return sorted(edge)


def get_grid_edge_cells(
    vor,
    idomain: list = None,
    idomain_path: Path = None,
    include_interiors: bool = True,
) -> list:
    """The ACTIVE cells on the model's edge, ascending and each once.

    A cell is on the edge when it borders the outside of the mesh -- it has more
    faces than it has neighbours -- or when one of its face-neighbours is
    inactive. That is the same topological test this always used when given no
    idomain; the branch that took one instead outlined the active area and asked
    which cells its RING touched, which is where four defects came from:

    - the inactive set came from a touch join, so a cell merely BORDERING the
      inactive polygon counted as inactive. On a 6x6 grid a face-aligned strip
      returned 14 cells against a true 16, having lost the true edge column and
      gained four interior cells.
    - a ring runs through vertices, so cells touching it at a single corner were
      kept: an inactive corner cell added its diagonal neighbour.
    - outer-ring and hole-ring hits were concatenated without de-duplication, so a
      cell on both came back twice -- 36 entries for the 27 cells around a
      one-cell hole, which the simple model's drain mode wrote as 36 DRN records
      at 1.5x the intended conductance.
    - ``vor.idomain_path`` overrode an explicit ``idomain=``, and ``idomain=[]``
      fell through to it.

    ``idomain`` is the list of INACTIVE cell ids -- not a mask, and not MF6's
    array; both are refused, naming the conversion. ``idomain_path`` is a vector
    file whose features cover the inactive cells, by AREA like every other source
    (ledger 210), reprojected to the grid. Give one or the other, never both.
    ``include_interiors=False`` keeps only each active component's outer edge,
    dropping the cells that ring a hole inside it.
    """

    inactive = _inactive_cells(vor, idomain, idomain_path)
    adjacency = vor.adjacent_cells_idx
    counts = vor.gdf_vorPolys.geometry.exterior.apply(lambda geom: len(geom.coords) - 1).to_numpy()
    outside = counts > np.fromiter((len(neighbours) for neighbours in adjacency), dtype=int)
    active = [cell for cell in range(len(adjacency)) if cell not in inactive]

    if not include_interiors:
        logger.info('collecting outer grid edge cells')
        return _outer_edge_only(adjacency, outside, active, inactive)

    logger.info('collecting grid edge cells')
    return [
        cell for cell in active
        if outside[cell] or any(neighbour in inactive for neighbour in adjacency[cell])
    ]

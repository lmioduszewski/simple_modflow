from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely as shp
from shapely.geometry import Polygon

from myflopy._deprecation import warn_deprecated
from myflopy._logging import get_logger
from myflopy.modflow.mf6.grid.geometry import cell_shares, miss_reason
from myflopy.modflow.utils.datatypes.readers import read_shp_gpkg

logger = get_logger(__name__)

#: A column name no file carries, so a caller that names no field keeps the
#: feature's own index as its key.
_NO_NAME_FIELD = "__myflopy_unnamed__"


def _measured_cells(gdf_vor_polys, geometries, names) -> pd.Series:
    """The cells each feature actually reaches, by the shared coverage rule.

    A polygon reaches the cells it COVERS by area, a line those it runs through,
    and a point the one cell holding it, with the tie and sliver rules of ledgers
    210-212. The spatial join this replaced took any touch, so a polygon drawn on
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


def _layer_names(path: Path) -> list[str] | None:
    """The layer names in a vector file, or ``None`` when they cannot be listed.

    ``geopandas.list_layers`` is 1.0+, and the declared floor is 0.14, so on an
    older geopandas the caller skips validation and GDAL reports the failure
    instead -- less helpful, never wrong.
    """

    lister = getattr(gpd, "list_layers", None)
    if lister is None:
        return None
    return lister(path)["name"].tolist()


def _read_locs(path: str | Path, layer: str | int | None = None) -> gpd.GeoDataFrame:
    """Read one vector file, optionally naming which layer inside it to use.

    ``layer=None`` is GDAL's own default -- the FIRST layer in the file -- which
    is what this module did before the parameter existed, so existing calls are
    unaffected. That default is the reason the parameter is needed: a GeoPackage
    holding several layers reads its first and only emits a ``UserWarning``, so
    reading the wrong table is silent in any code that does not surface warnings.

    Raises ``ValueError`` naming the layers actually present, because GDAL's own
    message (``Layer 'x' could not be opened``) does not say what IS in the file.
    """

    path = Path(path)
    if layer is None:
        return gpd.read_file(path)

    available = _layer_names(path)
    if available is not None:
        if isinstance(layer, int):
            if not 0 <= layer < len(available):
                raise ValueError(
                    f"{path.name} has {len(available)} layer(s), so index {layer} is out "
                    f"of range. Layers: {', '.join(map(repr, available))}."
                )
        elif layer not in available:
            raise ValueError(
                f"{path.name} has no layer {layer!r}. "
                f"Layers: {', '.join(map(repr, available))}."
            )
    return gpd.read_file(path, layer=layer)


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
    cell holding it, under the tie and sliver rules of ledgers 210-212. Passing a
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
        geoms = _read_locs(overlapping_geometry, layer)
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

    if geometries.crs is None:
        geometries.crs = gdf_vor_polys.crs
    elif geometries.crs != gdf_vor_polys.crs:
        logger.warning(
            'geometries are in %s but the grid is in %s; reprojecting them to '
            'match the grid', geometries.crs, gdf_vor_polys.crs,
        )
        geometries = geometries.to_crs(gdf_vor_polys.crs)

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
    gdf_locs = _read_locs(locs, layer).to_crs(crs)

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


def get_grid_edge_cells(vor, idomain: list = None, idomain_path: Path = None, include_interiors: bool = True) -> list:
    """
    Get edge cells for the voronoi grid, optionally excluding inactive cells.
    """
    df = vor.gdf_vorPolys
    idomain_path = vor.idomain_path if idomain_path is None else idomain_path

    adj_lengths = np.fromiter((len(vor.adjacent_cells_idx[idx]) for idx in df.index), dtype=int)
    face_counts = df.geometry.exterior.apply(lambda geom: len(geom.coords) - 1).to_numpy()
    edge_mask = face_counts > adj_lengths
    edges = df.index[edge_mask].tolist()

    if idomain_path or idomain:
        if idomain_path:
            idomain_geom = read_shp_gpkg(idomain_path).union_all()
            icells = get_vor_cells_as_series(vor.gdf_vorPolys, idomain_geom)[0]
        else:
            icells = idomain

        icell_set = set(icells)
        remaining_cells = df.loc[~df.index.isin(icell_set)]
        union_geom = remaining_cells.geometry.union_all()

        if isinstance(union_geom, Polygon):
            edge_geoms = [union_geom.exterior]
            interiors = union_geom.interiors if include_interiors else []
        elif isinstance(union_geom, shp.MultiPolygon):
            edge_geoms = [poly.exterior for poly in union_geom.geoms]
            interiors = [ring for poly in union_geom.geoms for ring in poly.interiors] if include_interiors else []
        else:
            edge_geoms, interiors = [], []

        edge_cells = [pd.Series(get_vor_cells_as_series(vor.gdf_vorPolys, geom)[0]) for geom in edge_geoms]
        if include_interiors and interiors:
            edge_cells += [pd.Series(get_vor_cells_as_series(vor.gdf_vorPolys, ring)[0]) for ring in interiors]

        flat = pd.concat(edge_cells)
        flat.index = flat.values

        logger.info('collecting grid edge cells, excluding inactive cells')
        return flat.loc[~flat.index.isin(icell_set)].index.tolist()

    logger.info('collecting all grid edge cells')
    return edges

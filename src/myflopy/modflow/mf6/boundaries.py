"""Shared GIS/grid helpers used by boundary-specific MF6 setup classes."""

from __future__ import annotations

from typing import TYPE_CHECKING

from geopandas import GeoDataFrame

if TYPE_CHECKING:
    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus as Vor
    from myflopy.modflow.mf6.simulation.base import SimulationBase

from pathlib import Path

import geopandas as gpd
import pandas as pd
import shapely as shp

from myflopy._logging import get_logger
from myflopy.areas import _same_crs
from myflopy.modflow.mf6.boundary_support import filter_inactive_cells
from myflopy.modflow.utils.datatypes.readers import read_shp_gpkg

logger = get_logger(__name__)

# Conversion factors
inches_to_feet = 1 / 12


def remove_duplicates(lst: list, seen: set = None):
    """Removes duplicates from a list"""
    seen = set() if seen is None else seen
    new_lst = []
    for num in lst:
        if num not in seen:
            new_lst.append(num)
            seen.add(num)
    return new_lst


class Boundaries:
    """Base helper for GIS-driven boundary workflows.

    Subclasses such as ``DRN``, ``GHB``, and recharge helpers use this class to
    intersect geometry with model cells, filter inactive cells, and optionally
    register model regions.
    """

    def __init__(
            self,
            model: SimulationBase = None,
            vor: Vor = None,
            shp_gpkg: Path = None,
            uid: str = None,
            crs: int = None,
            bound_type: str = None,
            idomain: list[int] | pd.Series = None,
            idomain_path: Path = None,
    ):
        """Parameters
        ----------
        model
            Parent model to which the boundary applies.
        vor
            Grid helper used for cell selection and geometry operations.
        shp_gpkg
            Path to the geometry source describing the boundary.
        uid
            Unique-id field in the geometry attributes.
        crs
            EPSG code or CRS string the geometry source is reprojected to when there
            is no grid CRS to follow (no ``vor``; a ``VoronoiGridPlus`` always declares
            one). The source is intersected against the grid's cells, which are in the
            grid's CRS, so it is reprojected to THAT, and a ``crs`` naming another
            raises ``ValueError`` -- here, and on any later read once a relabel makes
            the two differ; it moved the source out of the cells' frame (ledger
            466). Another spelling of the grid's own CRS (``2927`` for
            ``"EPSG:2927"``) is accepted. Default None: the grid's ``crs``, read LIVE,
            so relabelling the grid after this object is built reaches it (ledger 431).
        bound_type
            Short identifier such as ``"drn"`` or ``"rch"``.
        idomain, idomain_path
            Optional active-domain definition used to filter inactive cells.
        """

        self.model = model
        # Resolve the grid independently of the model: an explicit ``vor`` must
        # survive even when there is no model (grid-only use), and a missing
        # ``nper`` must not wipe out the grid.
        self.vor = vor if vor is not None else getattr(model, "vor", None)
        # The grid's CRS is read live (the `crs` property): a copy taken here stayed
        # on the old label after `grid.crs = ...` (ledger 431).
        self._crs = crs
        self._gdf_crs = None
        self.nper = getattr(model, "nper", None)
        self.bound_type = bound_type
        self.uid = uid
        self._gdf = None
        self._shp_gpkg = shp_gpkg
        self._intersections = None
        self._edge_intersections = None
        self._intersections_no_duplicates = None
        self._vor_bound_polys = None
        self._rch_scale = None
        self._sorted_cells_along_line = None
        self._inactive_cells = None
        self.idomain = idomain
        self.idomain_path = idomain_path
        self._boundary_dict = None
        self.limit_to_k33 = True
        self.limit_to_k33_by = 0.1
        self.verbose = False
        self.crs  # a `crs=` naming another CRS than the grid's is refused here

    #: The lazy caches computed in :attr:`crs`: the source frame and every cell set,
    #: area or value derived from it. A relabel of the grid drops them all (ledger
    #: 431). A subclass that caches more from them names its own here, and its
    #: property calls :meth:`_follow_the_grid_crs` before reading its cache -- a
    #: subclass cache left out stayed on the OLD label beside a refreshed source
    #: (Review of U6, ledger 431).
    _CRS_DERIVED: tuple[str, ...] = (
        "_gdf",
        "_intersections",
        "_edge_intersections",
        "_intersections_no_duplicates",
        "_vor_bound_polys",
        "_rch_scale",
    )

    @property
    def crs(self):
        """The CRS the geometry source is reprojected to: the grid's, else the ``crs`` given.

        The grid's is read live, so a grid relabelled after this object was built
        is followed -- and every frame and cell set cached LAZILY from the source is
        dropped then, because they were reprojected to the OLD label (ledger
        431). A :attr:`boundary_dict` a builder already assembled is not: it is
        an explicit build's output, so rebuild it after relabelling the grid.

        Raises
        ------
        ValueError
            The ``crs`` given names another CRS than the grid declares: the source
            would be reprojected away from the frame of the cells it is intersected
            against (ledger 466).
        """

        grid_crs = getattr(self.vor, "crs", None)
        if self._crs is None:
            return grid_crs
        if grid_crs is not None and not _same_crs(self._crs, grid_crs):
            raise ValueError(
                f"{type(self).__name__}(crs={self._crs!r}) names another CRS than its "
                f"grid's, {str(grid_crs)!r}. The source is intersected against the "
                "grid's cells, which are in the grid's CRS, so it is reprojected to the "
                "grid's: leave crs= out. If the grid's label is the wrong one, relabel "
                f"it with grid.crs = {self._crs!r}."
            )
        return self._crs

    def _follow_the_grid_crs(self):
        """Drop every cache built in another CRS than :attr:`crs` now names."""

        if self._gdf is None or self._gdf_crs == self.crs:
            return
        if self._gdf_crs is None:
            # A frame set from outside (`obj._gdf = frame`) is taken to be in the CRS
            # the object names when it is first read.
            self._gdf_crs = self.crs
            return
        logger.debug("grid relabelled %s -> %s; rebuilding the boundary source",
                     self._gdf_crs, self.crs)
        for name in self._CRS_DERIVED:
            setattr(self, name, None)

    @property
    def gdf(self):
        """gets a GeoDataFrame of the shapefile polygons"""
        self._follow_the_grid_crs()
        if self._gdf is None:
            if self._shp_gpkg is not None:
                gdf: GeoDataFrame = gpd.read_file(self._shp_gpkg)
                if self.uid is not None:
                    gdf = gdf.set_index(self.uid)
                if self.crs is not None:
                    if isinstance(self.crs, str):
                        gdf.to_crs(inplace=True, crs=self.crs)
                    else:
                        gdf.to_crs(inplace=True, epsg=self.crs)
                self._gdf = gdf
                self._gdf_crs = self.crs
        return self._gdf

    @property
    def inactive_cells(self):
        """The grid cells excluded by idomain (from ``idomain``/``idomain_path``), or ``None`` (cached)."""

        if self._inactive_cells is None:
            if self.idomain is None and self.idomain_path is None:
                return None
            elif self.idomain is not None and self.vor is not None:
                assert len(self.idomain) == self.vor.ncpl, 'idomain length must be equal to num cells in vor grid'
                if isinstance(self.idomain,list):
                    self.idomain = pd.Series(self.idomain)
                assert isinstance(self.idomain,pd.Series), 'error: could not make idomain a pd.Series'
                inactive_cells = self.idomain[self.idomain == 0].index.tolist()
                self._inactive_cells = inactive_cells
            elif self.idomain_path is not None and self.vor is not None:
                idomain = read_shp_gpkg(self.idomain_path)
                icells = self.vor.get_vor_cells_as_series(idomain.union_all())[0]
                self._inactive_cells = icells
            else:
                logger.warning(
                    'no grid or idomain given, so inactive cells cannot be '
                    'determined; treating every cell as active'
                )
        return self._inactive_cells

    @property
    def intersections(self):
        """gets a DataFrame with unique ids (uid) for each shapefile polygon and the associated
        intersecting voronoi grid cells"""
        self._follow_the_grid_crs()
        if self._intersections is None:
            vor_polys = self.vor.gdf_vorPolys
            df_intersect = self.gdf.geometry.apply(
                lambda geom: vor_polys[vor_polys.intersects(geom)].index.tolist())
            df_intersect.name = 'intersect'
            self._intersections = df_intersect
        return self._intersections

    @property
    def edge_intersections(self):
        """gets a DataFrame with unique ids (uid) for each shapefile polygon and the associated
        intersecting voronoi grid cells, but then filters for only those on a grid edge"""
        self._follow_the_grid_crs()
        if self._edge_intersections is None:
            if self.idomain_path is not None:
                edge_cells = self.vor.get_grid_edge(idomain_path=self.idomain_path)
            elif self.idomain is not None:
                edge_cells = self.vor.get_grid_edge(idomain=self.inactive_cells)
            else:
                edge_cells = self.vor.get_grid_edge()
            intersections = self.intersections
            filtered = intersections.apply(lambda x: [cell for cell in x if cell in edge_cells])
            self._edge_intersections = filtered
        return self._edge_intersections

    @property
    def intersections_no_duplicates(self):
        """gets a DataFrame of intersecting cells with duplicate cells removed"""
        self._follow_the_grid_crs()
        if self._intersections_no_duplicates is None:
            seen = set()
            no_dups = self.intersections.copy()
            lens = no_dups.apply(lambda x: len(x)).sort_values()
            lens.name = 'len'
            no_dups = pd.concat([lens, no_dups], axis=1)
            no_dups['no_dup'] = no_dups.loc[:, 'intersect'].apply(
                lambda x: remove_duplicates(x, seen))
            no_dups.drop(['len', 'intersect'], inplace=True, axis='columns')
            self._intersections_no_duplicates = no_dups
        return self._intersections_no_duplicates

    @property
    def vor_bound_polys(self):
        """gets the intersecting voronoi polygons equivalent to the shapefile polygons"""
        self._follow_the_grid_crs()
        if self._vor_bound_polys is None:
            vor_polys = self.intersections_no_duplicates.copy()
            vor_polys['geometry'] = vor_polys['no_dup'].apply(lambda x: self.vor.gdf_vorPolys.loc[x].union_all())
            self._vor_bound_polys = gpd.GeoDataFrame(vor_polys, geometry='geometry').drop(columns='no_dup')
        return self._vor_bound_polys

    @property
    def shp_to_vor_poly_scale(self):
        """gets a DataFrame giving the scaling between the areas of the shapefile vs. voronoi polys"""
        self._follow_the_grid_crs()
        if self._rch_scale is None:
            rch_scale = self.gdf.area / self.vor_bound_polys.area
            self._rch_scale = rch_scale
        return self._rch_scale

    def sorted_cells_along_line(self, idx=0):
        """
        method to get a list of voronoi model cells along the length of a line in order from
        the beginning to end of the linestring.
        :param idx: index of the line, if there was more than one provided in the shapefile
        :return: dataframe of intersecting cells, sorted with distance along line and centroids
        """

        line_gdf = self.gdf
        line = line_gdf.geometry[idx]
        assert isinstance(line, shp.geometry.linestring.LineString), 'geometry not a LineString!'
        intersecting_cells = self.vor.gdf_vorPolys.loc[self.intersections.to_list()[idx], :]
        intersecting_cells['centroid'] = intersecting_cells.centroid

        # Calculate the distance along the line for each cell centroid
        intersecting_cells['distance_along_line'] = intersecting_cells.centroid.apply(
            lambda point: line.project(point)
        )
        # Sort the cells based on the distance along the line
        sorted_cells = intersecting_cells.sort_values('distance_along_line')
        sorted_cells.index.name = 'cell'

        return sorted_cells

    def _require_vor(self):
        """Return the active Voronoi helper or raise when geometry-to-cell mapping is unavailable."""

        if self.vor is None:
            raise ValueError("A Voronoi grid is required for this boundary workflow")
        return self.vor

    def _indexed_gdf(self, name_field: str) -> GeoDataFrame:
        """Return the boundary GeoDataFrame indexed by the requested name field."""

        if self.gdf is None:
            raise ValueError("No boundary geometry is loaded")
        return self.gdf if self.gdf.index.name == name_field else self.gdf.set_index(name_field)

    def _candidate_cells_by_name(self, *, edges_only: bool = False) -> dict:
        """Return intersecting model cells keyed by boundary-feature name."""

        intersections = self.edge_intersections if edges_only else self.intersections
        return intersections.to_dict()

    def iter_polygon_boundary_features(
        self,
        *,
        name_field: str,
        edges_only: bool = False,
    ):
        """Yield ``(name, row, active_cells)`` tuples for polygon-driven boundary builders."""

        self._require_vor()
        gdf = self._indexed_gdf(name_field)
        for name, cell_nums in self._candidate_cells_by_name(edges_only=edges_only).items():
            yield name, gdf.loc[name], filter_inactive_cells(cell_nums, self.inactive_cells)

    @property
    def boundary_dict(self):
        """The assembled MF6 stress-period boundary data dict for this boundary type."""

        return self._boundary_dict

    def _register_region(
        self,
        name: str,
        *,
        cellids,
        layer: int | list[int] | None = None,
        geometry=None,
        tags: list[str] | None = None,
        metadata: dict | None = None,
        overwrite: bool = False,
    ):
        """Register a named model region for these boundary cells (no-op without a bound model)."""

        if self.model is None:
            return None
        return self.model.add_region_from_cells(
            name,
            cellids=cellids,
            layer=layer,
            category="boundary",
            package=self.bound_type,
            tags=tags,
            geometry=geometry,
            metadata=metadata,
            overwrite=overwrite,
        )

    def _register_boundary_groups(
        self,
        *,
        cellids_by_name: dict[str, list],
        layers_by_name: dict[str, int | list[int]] | None = None,
        geometries_by_name: dict[str, object] | None = None,
        region_name_prefix: str | None = None,
        combined_region_name: str | None = None,
        tags: list[str] | None = None,
        metadata_by_name: dict[str, dict] | None = None,
        overwrite: bool = False,
    ):
        """Register one model region per named boundary group (plus an optional combined region).

        Returns a ``{group_name: region}`` map (with a ``"__combined__"`` entry
        when ``combined_region_name`` is given); a no-op without a bound model.
        """

        if self.model is None:
            return {}

        registered = {}
        layers_by_name = {} if layers_by_name is None else layers_by_name
        geometries_by_name = {} if geometries_by_name is None else geometries_by_name
        metadata_by_name = {} if metadata_by_name is None else metadata_by_name

        combined_cells = []
        for group_name, cellids in cellids_by_name.items():
            region_name = (
                f"{region_name_prefix}_{group_name}"
                if region_name_prefix is not None
                else str(group_name)
            )
            region = self._register_region(
                region_name,
                cellids=cellids,
                layer=layers_by_name.get(group_name),
                geometry=geometries_by_name.get(group_name),
                tags=tags,
                metadata=metadata_by_name.get(group_name),
                overwrite=overwrite,
            )
            registered[group_name] = region
            combined_cells.extend(cellids)

        if combined_region_name is not None and combined_cells:
            registered["__combined__"] = self._register_region(
                combined_region_name,
                cellids=combined_cells,
                geometry=None,
                tags=tags,
                metadata={"groups": list(cellids_by_name)},
                overwrite=overwrite,
            )

        return registered


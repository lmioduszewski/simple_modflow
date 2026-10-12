"""Voronoi grid wrapper and convenience helpers for MODFLOW workflows.

``VoronoiGridPlus`` extends FloPy's Voronoi grid representation with cached
GeoDataFrame views, DISU connectivity helpers, plotting shortcuts, raster-based
surface helpers, and selection utilities used throughout ``myflopy``.
"""

from __future__ import annotations

import dataclasses
import numbers
import types
import warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely as shp
from flopy.utils.cvfdutil import get_disv_gridprops
from flopy.utils.triangle import Triangle
from flopy.utils.voronoi import VoronoiGrid

from myflopy._defaults import DEFAULT_CRS
from myflopy._logging import get_logger
from myflopy.grid_identity import grid_pin
from myflopy.modflow.mf6.grid.connectivity import build_disu_connectivity
from myflopy.modflow.mf6.grid.geometry import (
    adjust_cells_by_id as geometry_adjust_cells_by_id,
)
from myflopy.modflow.mf6.grid.geometry import (
    adjust_top_btm_overlaps as geometry_adjust_top_btm_overlaps,
)
from myflopy.modflow.mf6.grid.geometry import (
    calculate_distance as geometry_calculate_distance,
)
from myflopy.modflow.mf6.grid.geometry import (
    find_adjacent_cells as geometry_find_adjacent_cells,
)
from myflopy.modflow.mf6.grid.geometry import (
    find_adjacent_polygons as geometry_find_adjacent_polygons,
)
from myflopy.modflow.mf6.grid.geometry import frame_in_grid_crs, geometry_first
from myflopy.modflow.mf6.grid.geometry import (
    generate_grid_around_point as geometry_generate_grid_around_point,
)
from myflopy.modflow.mf6.grid.geometry import (
    generate_grid_coordinates as geometry_generate_grid_coordinates,
)
from myflopy.modflow.mf6.grid.geometry import (
    generate_grid_polygons as geometry_generate_grid_polygons,
)
from myflopy.modflow.mf6.grid.geometry import (
    get_centroids as geometry_get_centroids,
)
from myflopy.modflow.mf6.grid.geometry import (
    get_gdf_latlon as geometry_get_gdf_latlon,
)
from myflopy.modflow.mf6.grid.geometry import (
    get_gdf_vor_polys as geometry_get_gdf_vor_polys,
)
from myflopy.modflow.mf6.grid.geometry import (
    get_grid_centroid as geometry_get_grid_centroid,
)
from myflopy.modflow.mf6.grid.geometry import (
    get_latlon as geometry_get_latlon,
)
from myflopy.modflow.mf6.grid.geometry import (
    get_overlapping_area as geometry_get_overlapping_area,
)
from myflopy.modflow.mf6.grid.geometry import (
    get_vor_idx_from_geometry as geometry_get_vor_idx_from_geometry,
)
from myflopy.modflow.mf6.grid.geometry import (
    get_vor_idx_from_geometry_idx as geometry_get_vor_idx_from_geometry_idx,
)
from myflopy.modflow.mf6.grid.geometry import (
    get_voronoi_polygons as geometry_get_voronoi_polygons,
)
from myflopy.modflow.mf6.grid.geometry import (
    points_to_latlon as geometry_points_to_latlon,
)
from myflopy.modflow.mf6.grid.geometry import (
    reconcile_surfaces as geometry_reconcile_surfaces,
)
from myflopy.modflow.mf6.grid.geometry import (
    set_k_vor as geometry_set_k_vor,
)
from myflopy.modflow.mf6.grid.geometry import (
    shared_face_length as geometry_shared_face_length,
)
from myflopy.modflow.mf6.grid.geometry import (
    voronoi_refine_by_point as geometry_voronoi_refine_by_point,
)
from myflopy.modflow.mf6.grid.helpers import (
    get_griddata_from_disu,
    get_griddata_from_gsf,
)
from myflopy.modflow.mf6.grid.plotting import GridPlots
from myflopy.modflow.mf6.grid.selection import (
    get_grid_edge_cells,
    get_model_boundary_polygons,
)
from myflopy.modflow.mf6.grid.selection import (
    get_vor_cells_as_dict as selection_get_vor_cells_as_dict,
)
from myflopy.modflow.mf6.grid.selection import (
    get_vor_cells_as_series as selection_get_vor_cells_as_series,
)
from myflopy.modflow.mf6.grid.surfaces import (
    get_cell_areas as surface_get_cell_areas,
)
from myflopy.modflow.mf6.grid.surfaces import (
    get_domain as surface_get_domain,
)
from myflopy.modflow.mf6.grid.surfaces import (
    get_gdf_topbtm_multilyr as surface_get_gdf_topbtm_multilyr,
)
from myflopy.modflow.mf6.grid.surfaces import (
    get_normal_from_strike_and_dip as surface_get_normal_from_strike_and_dip,
)
from myflopy.modflow.mf6.grid.surfaces import (
    get_origin_xy as surface_get_origin_xy,
)
from myflopy.modflow.mf6.grid.surfaces import (
    get_raster_from_strike_dip as surface_get_raster_from_strike_dip,
)
from myflopy.modflow.mf6.grid.surfaces import (
    get_raster_vals_at_centroids as surface_get_raster_vals_at_centroids,
)
from myflopy.modflow.mf6.grid.tessellation import tessellate
from myflopy.modflow.utils.datatypes.readers import as_path

logger = get_logger(__name__)


def normalize_crs(crs):
    """Accept a bare EPSG code as well as a full CRS identifier.

    ``2927``, any numpy integer of that value and ``"2927"`` become
    ``"EPSG:2927"``; anything else passes through untouched for pyproj/rasterio
    to interpret.

    Worth doing because the failure was silent until it was baffling. pyproj --
    and therefore geopandas ``to_crs`` -- accepts a bare ``"2927"``, so a grid
    built with it reprojects, clips and plots perfectly well. Only rasterio's
    ``warp_transform`` rejects it, so the first thing to break is raster
    sampling, several steps later, with ``CRSError: The WKT could not be parsed``
    -- which names WKT at someone who typed a number. Reported 2026-08-27.

    A bare integer string is never valid WKT or PROJ, so reading it as an EPSG
    code cannot be wrong.
    """

    # `numbers.Integral`, not `int`: `numpy.int64` subclasses neither `int` nor
    # anything rasterio reads, and slipped through as itself (ledger 370).
    if isinstance(crs, numbers.Integral) and not isinstance(crs, (bool, np.bool_)):
        return f"EPSG:{int(crs)}"
    if isinstance(crs, str) and crs.strip().isdigit():
        return f"EPSG:{crs.strip()}"
    return crs


def _readable_crs(crs, *, spelling: str):
    """``crs`` through :func:`normalize_crs`, or a ``ValueError`` in ``spelling``.

    The one check behind the constructor's ``crs=`` and the :attr:`VoronoiGridPlus.crs`
    setter (ledgers 370, 389): ``None`` and anything pyproj cannot read are
    refused, each message opening with the caller's own spelling -- ``spelling``
    holds one ``{}`` for the value, e.g. ``"grid.crs = {}"``.
    """

    crs = normalize_crs(crs)
    if crs is None:
        raise ValueError(
            f"{spelling.format(None)}: a grid needs the CRS its vertices are in -- its lat/lon "
            "frames and grid_centroid are reprojections of the cells. Give the real "
            "one, e.g. 'EPSG:2926'; a grid built with no crs= takes the default "
            "'EPSG:2927'."
        )
    import pyproj

    try:
        pyproj.CRS.from_user_input(crs)
    except pyproj.exceptions.CRSError as error:
        raise ValueError(f"{spelling.format(repr(crs))}: pyproj cannot read it as a CRS ({error}).") from error
    return crs


class VoronoiGridPlus(VoronoiGrid):
    """An unstructured Voronoi (DISV) grid with MODFLOW-modeling conveniences.

    Extends FloPy's ``VoronoiGrid`` with the geometry and helpers myflopy needs to
    build and post-process unstructured models: DISV grid properties
    (:meth:`get_disv_gridprops` -> ``ncpl``/``vertices``/``cell2d``), DISU
    connectivity (``iac``/``ja``/``cl12``/``hwva``), cell-centroid and adjacency
    lookups, CRS handling, and the per-cell ``gdf_vorPolys`` / ``gdf_topbtm``
    GeoDataFrames the surface-aware builders read.

    Build one from a :class:`~myflopy.modflow.mf6.grid.triangle.TriangleGrid`
    triangulation, or declaratively from GeoPackages via
    ``mf.GridSpec.voronoi(...).resolve(ws)``. Put the result in a
    :class:`~myflopy.specs.ModelContext` so the package-first GIS helpers can map
    features onto its cells, and pass ``get_disv_gridprops()`` to ``mf.disv``.

    The constructor's parameters are documented on :meth:`__init__`. After a build
    from a triangulation, ``tessellation_report`` says what the repair changed.
    """

    def __init__(
        self,
        tri: Triangle = None,
        crs: str = DEFAULT_CRS,
        rasters: list | Path = None,
        name: str = 'voronoi_grid',
        qhull_options: str = None,
        idomain_path: Path = None,
        verts: np.ndarray = None,
        iverts: list[list] = None,
        xcyc: np.ndarray = None,
        tessellation: str = "dual",
        **kwargs,
    ):
        """Create a Voronoi view from Triangle output or DISU-style arrays.

        Two forms: from a built triangulation (``tri``), or from explicit vertex-grid
        arrays (``verts`` + ``iverts`` + ``xcyc``, all three). ``tri`` wins when both
        are given; neither raises ``ValueError``.

        Parameters
        ----------
        tri : flopy.utils.triangle.Triangle, optional
            A BUILT triangulation -- FloPy's ``Triangle`` or myflopy's
            :class:`~myflopy.modflow.mf6.grid.triangle.TriangleGrid` (a subclass),
            after ``.build()``. The cells are its Voronoi dual, made by FloPy's
            ``tri2vor`` and, under ``tessellation="dual"``, repaired against the dual
            of Triangle's own mesh. FloPy raises ``TypeError`` for any other type.
            Default None: the array form.
        crs : str, int or pyproj.CRS, optional
            The CRS the grid's coordinates are in, stamped on ``gdf_vorPolys`` and
            every frame derived from it. Vector features in another CRS are
            reprojected to it; a raster in another CRS is sampled by reprojecting the
            cells to the raster's. A bare EPSG code (``2927``, a numpy integer, or
            ``"2927"``) is stored as ``"EPSG:2927"``; anything else is stored as
            given. Default ``"EPSG:2927"`` (the package-wide default, which places a
            location-agnostic model at that CRS's origin on purpose). ``None``, and
            anything pyproj cannot read, is REFUSED (``ValueError`` naming
            ``crs=``) before any cell is built: the constructor reprojects the cells
            to lat/lon for ``grid_centroid``, which an unreferenced frame cannot do
            (ledgers 370, 389). Assigning :attr:`crs` later RELABELS the grid
            and every frame it caches; it never reprojects.
        rasters : list of (str or os.PathLike or float), optional
            Surfaces to sample onto the cells the first time :attr:`gdf_topbtm` is
            read, in order: the model top, then one bottom per layer (column labels
            ``0``, ``1``, ...; ``nlay`` is their count less one). Each item is a
            raster path, ``str`` or ``os.PathLike`` (area-weighted mean of the
            pixels in each cell, ``gaps="warn"``), or a number (a constant column);
            a bare path is a one-item list, the model top alone. Any other item
            raises ``TypeError`` naming it at that first read (see
            :meth:`get_raster_vals_at_centroids`). Default None: no surfaces until
            one is assigned to :attr:`gdf_topbtm`.
        name : str, optional
            A label stored as ``grid.name``; nothing in the grid's own methods reads
            it. Default ``"voronoi_grid"``.
        qhull_options : str, optional
            Passed as ``qhull_options`` to ``scipy.spatial.Voronoi`` through FloPy's
            ``tri2vor``. The ``tri`` form only: given with the array form it raises
            ``TypeError``, since nothing there runs Qhull. Default None: SciPy's own.
        idomain_path : str or os.PathLike, optional
            A vector file (``.shp`` or ``.gpkg``) whose features cover the INACTIVE
            cells, by area. Stored as :attr:`idomain_path` and read by
            :meth:`get_grid_edge` when that call is given neither ``idomain`` nor
            ``idomain_path`` (so also by the ``edges_only=`` boundary helpers when
            their context carries no ``domain``). A ``str`` or ``os.PathLike``, stored
            as a ``pathlib.Path``; another type raises ``TypeError`` (ledger
            457). Default None.
        verts : numpy.ndarray, optional
            The array form: vertex coordinates, shape ``(nverts, 2)``. Needs
            ``iverts`` and ``xcyc`` too.
        iverts : list of list of int, optional
            The array form: each cell's ring as indices into ``verts``, one list per
            cell (``ncpl = len(iverts)``).
        xcyc : numpy.ndarray, optional
            The array form: each cell's centre, shape ``(ncpl, 2)`` -- the DISV
            ``cell2d`` ``xc``/``yc`` points MODFLOW 6 connects. Stored as
            ``grid.points``. The three arrays are taken exactly as given -- never
            repaired, so a grid rebuilt from a saved run or snapshot is the grid that
            ran.
        tessellation : {"dual", "flopy"}, optional
            How cells are made from ``tri``; ignored for ``verts``/``iverts``/``xcyc``
            (but validated either way: any other value raises ``ValueError``).
            ``"dual"`` (the default) runs FloPy's ``tri2vor`` and then repairs every cell
            that does not match the dual of Triangle's own mesh -- a reflex corner's
            cell drawn out of order, a vertex on the boundary dropped, a duplicate vertex
            id MF6 cannot match (ledger 316), a domain part lost (ledger 317). Cells
            FloPy got right keep its rings exactly, so a grid with nothing wrong is
            byte-identical to FloPy's. ``"flopy"`` keeps ``tri2vor``'s cells as they
            are, defects included -- the escape hatch for reproducing a grid made before
            the repair. ``tessellation_report`` records what was changed and the
            pre-repair ``grid_pin`` (``legacy_pin``).
        **kwargs
            The ``tri`` form only: forwarded to FloPy's ``VoronoiGrid``/``tri2vor``
            and on to ``scipy.spatial.Voronoi`` as keyword arguments. Given with the
            array form they raise ``TypeError`` naming them (they were dropped in
            silence). ``idomain=`` raises ``TypeError`` in either form: a grid takes
            no active domain (ledger 388) -- give the model one with
            ``ModelContext(domain=...)``, or exclude cells from the grid edge with
            :meth:`get_grid_edge`'s ``idomain=`` or the ``idomain_path`` above.

        Raises
        ------
        ValueError
            ``crs`` is ``None`` or unreadable; ``tessellation`` is not ``"dual"``
            or ``"flopy"``; neither ``tri`` nor all three of
            ``verts``/``iverts``/``xcyc`` was given; or the ``"dual"`` repair found
            a mesh whose dual is not a tiling
            (:class:`~myflopy.modflow.mf6.grid.tessellation.NotConformingDelaunay`).
        TypeError
            ``idomain=`` was passed, or the array form was given ``qhull_options``
            or other keyword arguments it cannot use.
        """
        if "idomain" in kwargs:
            raise TypeError(
                "VoronoiGridPlus takes no idomain=: it was accepted and dropped, and no "
                "grid ever read it (ledger 388). Give the model its active domain "
                "with mf.ModelContext(grid=..., domain=...) -- the boundary helpers' "
                "edges_only= reads it per layer -- or exclude cells from the grid edge "
                "with grid.get_grid_edge(idomain=[inactive ids]) or idomain_path=."
            )
        # Checked FIRST: the setter below refuses the same values (ledger 370), but
        # only after the triangulation's cells were built, and naming `grid.crs =`
        # at a caller who typed `crs=` (ledger 389).
        crs = _readable_crs(crs, spelling="VoronoiGridPlus(crs={})")
        if not tri and (kwargs or qhull_options is not None):
            unused = sorted(kwargs) + (["qhull_options"] if qhull_options is not None else [])
            raise TypeError(
                f"VoronoiGridPlus got {', '.join(unused)} with the array form "
                "(verts/iverts/xcyc): they configure the Voronoi diagram of a "
                "triangulation and have nothing to act on here, so they were dropped "
                "in silence. Pass tri= to build from a triangulation, or leave them out."
            )
        logger.info('building the Voronoi grid')
        if tessellation not in ("dual", "flopy"):
            raise ValueError(
                f"tessellation={tessellation!r} is not a tessellation. Use 'dual' (repair "
                "FloPy's cells against Triangle's mesh, the default) or 'flopy' (keep "
                "tri2vor's cells as they are)."
            )
        self.tessellation_report = None
        #: "dual" or "flopy" for a grid made from a triangulation, None for one rebuilt
        #: from arrays -- what a run records about the grid it was built on.
        self.tessellation = tessellation if tri else None
        if tri:
            # FloPy's tri2vor helper currently emits harmless RuntimeWarnings from
            # point_in_polygon when domain edges are horizontal. Suppress only that
            # known warning pattern so realistic Triangle-built grids stay quiet.
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message="invalid value encountered in divide",
                    category=RuntimeWarning,
                    module=r"flopy\.utils\.geometry",
                )
                warnings.filterwarnings(
                    "ignore",
                    message="divide by zero encountered in divide",
                    category=RuntimeWarning,
                    module=r"flopy\.utils\.geometry",
                )
                # tri2vor averages an EMPTY vertex list for each cell it is about to pop --
                # every cell of a boundary's second part, say -- before the repair below
                # puts those cells back.
                warnings.filterwarnings("ignore", message="Mean of empty slice", category=RuntimeWarning)
                warnings.filterwarnings(
                    "ignore",
                    message="invalid value encountered in divide",
                    category=RuntimeWarning,
                    module=r"numpy\._core\._methods",
                )
                super().__init__(tri, qhull_options=qhull_options, **kwargs)
            if tessellation == "dual":
                self._repair_tessellation(tri)
        elif all(v is not None for v in [verts, iverts, xcyc]):
            self.verts = verts
            self.iverts = iverts
            self.points = xcyc
            self.ncpl = len(iverts)
            self.nverts = verts.shape[0]
        else:
            raise ValueError('Either tri or (verts, iverts, xcyc) must be provided.')

        self.tri = tri
        self.crs_latlon = "EPSG:4326"
        self._centroids = None
        self._gdf_latlon = None
        self._nlay = None
        self._latlon = None
        self.rasters = rasters
        self.crs = normalize_crs(crs)
        self.name = name
        self.x_coords_by_node = []
        self.y_coords_by_node = []
        self._gdf_vorPolys = None
        self._iac = None
        self._ia = None
        self._nja = None
        self._idomain_path = None
        self._gdf_topbtm = None
        self._excavations = {}
        self._area_list = None
        self._adjacent_cells_idx = None
        self._ja, self._cl12, self._hwva = None, None, None
        self.centroids_x, self.centroids_y = self.get_centroids()
        self.idomain_path = idomain_path
        self.config = {'scrollZoom': True}
        self.scatt_layout = {'height': 1000, 'width': 1000, 'dragmode': 'pan'}
        self.grid_centroid = self.get_grid_centroid()
        logger.info('Voronoi grid ready: %d cells', getattr(self, 'ncpl', -1))

    # Attributes dropped from pickles: the Triangle builder plus caches that
    # recompute from the stored mesh geometry. Keeping them would bloat the
    # pickle and tie it more tightly to dependency versions.
    _PICKLE_VOLATILE = (
        "tri",
        "_centroids",
        "_gdf_latlon",
        "_latlon",
        "_gdf_vorPolys",
        "_iac",
        "_ia",
        "_nja",
        "_ja",
        "_cl12",
        "_hwva",
        "_area_list",
        "_adjacent_cells_idx",
    )

    def __getstate__(self):
        """Return a lean picklable state, dropping ``tri`` and recomputable caches."""

        return {
            key: value
            for key, value in self.__dict__.items()
            if key not in self._PICKLE_VOLATILE
        }

    def __setstate__(self, state):
        """Restore state, leaving dropped caches empty for lazy recomputation.

        Never repairs: a pickled grid is the grid something was built on. One saved
        before ledger 316 has no ``tessellation_report`` and reads ``None``. One saved
        while :attr:`crs` was a plain attribute carries it under ``"crs"``, which the
        property would shadow; it is moved to ``"_crs"`` (ledger 368).
        """

        state = dict(state)
        if "crs" in state:
            state.setdefault("_crs", state.pop("crs"))
        self.__dict__.update(state)
        for key in self._PICKLE_VOLATILE:
            self.__dict__.setdefault(key, None)
        self.__dict__.setdefault("tessellation_report", None)
        self.__dict__.setdefault("tessellation", None)

    def _repair_tessellation(self, tri) -> None:
        """Replace FloPy's cells by :func:`~myflopy.modflow.mf6.grid.tessellation.tessellate`'s.

        Runs inside ``__init__`` before any cache is built, so nothing needs
        invalidating. Records the pre-repair ``grid_pin`` when anything changed, which
        is what lets a refusal say a package was resolved on the unrepaired grid.
        """

        old_verts, old_iverts, old_points = self.verts, self.iverts, self.points
        verts, iverts, points, report = tessellate(tri, old_verts, old_iverts, old_points)
        if report.changed:
            legacy = types.SimpleNamespace(
                get_disv_gridprops=lambda: get_disv_gridprops(
                    np.asarray(old_verts), old_iverts, xcyc=np.asarray(old_points)
                )
            )
            report = dataclasses.replace(report, legacy_pin=grid_pin(legacy))
            self.verts = np.asarray(verts)
            self.iverts = iverts
            self.points = np.asarray(points)
            self.ncpl = len(iverts)
            self.nverts = self.verts.shape[0]
        self.tessellation_report = report

    @property
    def vor_list(self):
        """Return the Voronoi polygons as a plain geometry list."""
        return self.gdf_vorPolys.geometry.to_list()

    @property
    def cell_list(self):
        """Return the cell indices for the Voronoi grid."""
        return self.gdf_vorPolys.index.to_list()

    @property
    def area_list(self):
        """Return cached polygon areas for all Voronoi cells."""
        if self._area_list is None:
            self._area_list = [cell.area for cell in self.vor_list]
        return self._area_list

    @property
    def idomain_path(self):
        """Path used by some edge helpers to read external idomain data."""
        return self._idomain_path

    @idomain_path.setter
    def idomain_path(self, value):
        """Set the optional idomain path used by edge-detection helpers.

        A ``str`` or ``os.PathLike`` is stored as a ``pathlib.Path``; anything else
        but ``None`` raises ``TypeError`` naming it (ledger 457) -- a ``str``
        used to fail an ``assert``.
        """
        self._idomain_path = as_path(value, name="idomain_path")

    @property
    def nlay(self):
        """Infer the number of layers from ``gdf_topbtm``."""
        if self._nlay is None:
            self._nlay = len(self.gdf_topbtm.drop('geometry', axis=1).columns) - 1
        return self._nlay

    @property
    def gdf_topbtm(self):
        """Return cached top/bottom surfaces, building them from rasters if needed.

        The column order is normalized on READ, not on write, because two routes
        bypass the setter entirely. ``__setstate__`` writes ``__dict__`` directly,
        so a reloaded grid otherwise keeps whatever order it was saved in; and an
        in-place edit (``vor.gdf_topbtm[2] = ...``, which the canonical notebook
        does) calls this getter and mutates the frame without the setter ever
        running.

        The normalized frame is STORED and returned -- never a normalized copy.
        Returning a copy would hand that in-place edit a temporary, and the write
        would be silently lost: measured, a bottom set to 99.0 read back as 5.0.
        That would be a worse bug than the one this fixes.
        """
        if self._gdf_topbtm is None and self.rasters is not None:
            self._gdf_topbtm = self.get_gdf_topbtm_multilyr(rasters=self.rasters)
        self._gdf_topbtm = geometry_first(self._gdf_topbtm)
        return self._gdf_topbtm

    @gdf_topbtm.setter
    def gdf_topbtm(self, value):
        """Set the cached top/bottom surface GeoDataFrame.

        ``_nlay`` is dropped as well: it caches on first read and nothing
        invalidated it, so re-attaching a frame with a different layer count left
        the old count in place.

        A GeoDataFrame must be in the grid's :attr:`crs`, as for
        :attr:`gdf_vorPolys`: with no CRS it takes the grid's label, and in another
        CRS it is refused naming ``grid.crs = ...`` and ``frame.to_crs(grid.crs)``
        (ledger 431).
        """
        self._gdf_topbtm = geometry_first(
            frame_in_grid_crs(self, value, spelling="grid.gdf_topbtm = frame")
        )
        self._nlay = None
        # The excavation records ride WITH these elevations: whoever writes a new
        # table drops the old records, and `LayerBuildResult.attach_to_grid` sets its
        # own after the table (the excavation design, section 2.2).
        self._excavations = {}

    @property
    def excavations(self) -> dict:
        """The excavation records published with :attr:`gdf_topbtm`, by name.

        Set by ``LayerBuildResult.attach_to_grid`` (``LayerStack.build(attach=True)``)
        to the result's ``excavations`` -- :class:`~myflopy.surfaces.ExcavationRecord`
        objects saying which cells were dug and the bed at each -- and CLEARED by any
        other write of ``gdf_topbtm``, so it never describes elevations it did not
        come with. What ``mf.lak`` reads to tell an excavated lake from one nobody
        dug. ``{}`` when nothing was published (and on a grid pickled before it
        existed).
        """

        return self.__dict__.get("_excavations") or {}

    @excavations.setter
    def excavations(self, value) -> None:
        """Publish records (a mapping of name -> ExcavationRecord); ``None`` clears them."""

        self._excavations = dict(value or {})

    def get_disu_connectivity(self, *, tol: float = 0.0, validate: bool = True):
        """Compute and cache the DISU connection table MF6 would build from this grid's DISV.

        ``iac``/``ja``/``cl12``/``hwva``/``nja``, from :meth:`get_disv_gridprops`
        through :func:`~myflopy.modflow.mf6.grid.connectivity.build_disu_connectivity`
        -- MF6's own DISV geometry (``DisvGeom.f90``), so a DISU written from them
        reproduces the DISV to rounding. ``cl12`` is per SIDE: at ``(n, m)`` the
        distance from cell ``n``'s DISV centre to the shared face. It held the full
        centre-to-centre spacing from 2025-07-02 until ledger 328, which made a
        DISU's horizontal conductance about half the DISV's and a legacy lake
        wall's ``CONNLEN`` about twice MF6's.

        Parameters
        ----------
        tol : float, default 0.0
            A face whose width is ``<= tol`` (model length units) is not connected;
            ``0.0`` drops only zero-width faces.
        validate : bool, default True
            Check the table before caching it (rows headed by their own cell, every
            connection in both rows, ``hwva`` symmetric, ``cl12`` finite and
            non-negative); ``False`` skips the checks.

        Returns
        -------
        tuple
            ``iac`` (int64), ``ja`` (int64, zero-based), ``cl12`` (float64),
            ``hwva`` (float64), ``nja`` (int), also cached on the grid and read back
            by the :attr:`iac`/:attr:`ja`/:attr:`cl12`/:attr:`hwva`/:attr:`nja`
            properties; the :attr:`ia` row pointer is dropped and recomputed.

        Raises
        ------
        TypeError
            ``tol`` is not a number.
        ValueError
            ``validate`` found a defect in the table (the message names the cells).

        Examples
        --------
        >>> import numpy as np
        >>> from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
        >>> verts = np.array([[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]], float)
        >>> grid = VoronoiGridPlus(verts=verts, iverts=[[0, 3, 4, 1], [1, 4, 5, 2]],
        ...                        xcyc=np.array([[0.5, 0.5], [1.5, 0.5]]))
        >>> iac, ja, cl12, hwva, nja = grid.get_disu_connectivity()
        >>> ja.tolist(), cl12.tolist(), hwva.tolist()
        ([0, 1, 1, 0], [0.0, 0.5, 0.0, 0.5], [0.0, 1.0, 0.0, 1.0])
        """
        iac, ja, cl12, hwva, nja = build_disu_connectivity(
            self.get_disv_gridprops(), tol=tol, validate=validate
        )
        self._iac = iac
        self._ia = None  # invalidate the row-pointer cache; recomputed from iac on demand
        self._ja = ja
        self._cl12 = cl12
        self._hwva = hwva
        self._nja = nja
        return iac, ja, cl12, hwva, nja

    @property
    def iac(self):
        """Return the cached DISU ``iac`` array, computing it if needed."""
        if self._iac is None:
            self.get_disu_connectivity()
        return self._iac

    @property
    def ia(self):
        """Return the CRS row pointer ``ia = [0, cumsum(iac)]`` for the ``ja`` arrays.

        ``ia[cell]`` is the offset of ``cell``'s block in ``ja``/``cl12``/``hwva``
        (the diagonal/self entry is first). Cached so per-cell connection lookups are
        O(1) instead of recomputing ``sum(iac[:cell])`` -- an O(N) scan -- every call,
        which is O(N^2) when sweeping every cell of a large lake footprint.
        """
        if self._ia is None:
            iac = self.iac
            ia = np.zeros(len(iac) + 1, dtype=np.int64)
            np.cumsum(iac, out=ia[1:])
            self._ia = ia
        return self._ia

    @property
    def ja(self):
        """Return the cached DISU ``ja`` array, computing it if needed."""
        if self._ja is None:
            self.get_disu_connectivity()
        return self._ja

    @property
    def cl12(self):
        """Per-side connection lengths: at ``(n, m)``, cell ``n``'s DISV centre to the shared face.

        MF6's DISU ``CL12`` -- "the center of cell n [to] the shared face with each
        adjacent m cell" -- computed as MF6 computes a DISV face (the perpendicular
        from ``cell2d``'s ``xc``/``yc`` to the line through the two shared vertices),
        indexed like :attr:`ja`, zero on each row's own entry. NOT symmetric: the two
        rows of one connection hold the two sides. Until ledger 328 this held the full
        centre-to-centre spacing on both -- 1.74-2.35x (median 2.0) the lake cell's own
        side over the legacy lake walls of two real runs. Computed by
        :meth:`get_disu_connectivity` on first read.

        Examples
        --------
        >>> import numpy as np
        >>> from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
        >>> verts = np.array([[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]], float)
        >>> grid = VoronoiGridPlus(verts=verts, iverts=[[0, 3, 4, 1], [1, 4, 5, 2]],
        ...                        xcyc=np.array([[0.5, 0.5], [1.25, 0.5]]))
        >>> grid.cl12.tolist()   # row 0: cell 0's side (0.5); row 1: cell 1's (0.25)
        [0.0, 0.5, 0.0, 0.25]
        """
        if self._cl12 is None:
            self.get_disu_connectivity()
        return self._cl12

    @property
    def hwva(self):
        """Shared-face widths, indexed like :attr:`ja`: the distance between the two shared vertices.

        MF6's DISV face width (``DisvGeom.f90``'s ``distance``), the same in both rows
        of a connection and zero on each row's own entry. Computed by
        :meth:`get_disu_connectivity` on first read.

        Examples
        --------
        >>> import numpy as np
        >>> from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
        >>> verts = np.array([[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]], float)
        >>> grid = VoronoiGridPlus(verts=verts, iverts=[[0, 3, 4, 1], [1, 4, 5, 2]],
        ...                        xcyc=np.array([[0.5, 0.5], [1.5, 0.5]]))
        >>> grid.hwva.tolist()
        [0.0, 1.0, 0.0, 1.0]
        """
        if self._hwva is None:
            self.get_disu_connectivity()
        return self._hwva

    @property
    def nja(self):
        """Return the cached number of ``ja`` entries."""
        if self._nja is None:
            self.get_disu_connectivity()
        return self._nja

    @property
    def crs(self):
        """The CRS the grid's coordinates are in, and the label on every frame it caches.

        Setting it RELABELS, never reprojects: it declares that the vertices already
        are in that CRS. The cached :attr:`gdf_vorPolys` and :attr:`gdf_topbtm` take
        the new label in place (the same objects, so a reference held elsewhere agrees
        too), and the lat/lon frames (:attr:`gdf_latlon`, :attr:`latlon`) and
        ``grid_centroid``, which are REPROJECTIONS of the cells, are rebuilt from it.
        Assigning only the attribute left one grid reporting two CRSs -- the grid's
        own and its frames' -- and the readers split between them (ledger 368).
        To move the coordinates instead, reproject a frame:
        ``grid.gdf_vorPolys.to_crs(...)``.

        A bare EPSG code (``2927``, ``numpy.int64(2927)`` or ``"2927"``) is stored as
        ``"EPSG:2927"`` (:func:`normalize_crs`); anything else as given. ``None`` and
        a value pyproj cannot read are REFUSED (``ValueError``) before anything
        changes: a grid with no CRS cannot build its lat/lon frames, and the
        refusal used to come from geopandas halfway through the relabel, leaving
        the grid and its frames on ``None`` (ledger 370).
        """
        return self.__dict__.get("_crs")

    @crs.setter
    def crs(self, value):
        value = _readable_crs(value, spelling="grid.crs = {}")
        changed = "_crs" in self.__dict__ and self.__dict__["_crs"] != value
        self.__dict__["_crs"] = value
        if not changed:
            return
        for key in ("_gdf_vorPolys", "_gdf_topbtm"):
            frame = self.__dict__.get(key)
            if frame is not None and hasattr(frame, "set_crs"):
                frame.set_crs(value, allow_override=True, inplace=True)
        # Reprojections of the cells: computed from the OLD label, so they are wrong.
        self._gdf_latlon = None
        self._latlon = None
        if "grid_centroid" in self.__dict__:
            self.grid_centroid = self.get_grid_centroid()
        logger.debug("grid relabelled as %s; lat/lon frames rebuilt", value)

    @property
    def gdf_vorPolys(self):
        """Return the Voronoi polygons as a GeoDataFrame."""
        if self._gdf_vorPolys is None:
            self._gdf_vorPolys = self.get_gdf_vorPolys(crs=self.crs)
        return self._gdf_vorPolys

    @gdf_vorPolys.setter
    def gdf_vorPolys(self, value):
        """Set the cached Voronoi polygon GeoDataFrame.

        The frame must be in the grid's :attr:`crs`: one with no CRS takes the
        grid's label, and one in another CRS is REFUSED (``ValueError``) naming
        ``grid.crs = ...`` to relabel the grid and ``frame.to_crs(grid.crs)`` to
        reproject the frame -- it was stored as given, leaving one grid with two
        CRSs (ledger 431). ``None`` drops the cache, so the next read rebuilds
        it from the vertices.
        """
        self._gdf_vorPolys = frame_in_grid_crs(
            self, value, spelling="grid.gdf_vorPolys = frame"
        )

    @property
    def adjacent_cells_idx(self):
        """Return cached adjacency lists for each Voronoi cell."""
        if self._adjacent_cells_idx is None:
            self._adjacent_cells_idx = self.find_adjacent_polygons(self.gdf_vorPolys)
        return self._adjacent_cells_idx

    def get_disu_props(self):
        """Populate and cache the standard DISU connectivity properties."""
        self.iac
        self.ja
        self.nja
        self.cl12
        self.hwva

    def get_voronoi_polygons(self):
        """Return the Voronoi polygons as shapely geometry objects."""
        return geometry_get_voronoi_polygons(self)

    @property
    def plot(self) -> GridPlots:
        """The plotting verbs for this grid: ``map``, ``section``, ``grid``.

        ``vor.plot.map(values=...)`` for a basemap choropleth,
        ``vor.plot.section(line)`` for the grid's layers in cross-section,
        ``vor.plot.grid()`` -- or just ``vor.plot()`` -- for the bare mesh.
        Everything returned is a :class:`~myflopy.viz.Picture`, except
        ``section``'s default, a Matplotlib ``Figure`` (``backend="plotly"`` for
        the picture).

        **This deliberately shadows FloPy's inherited ``VoronoiGrid.plot()``.**
        That renderer is not lost: ``vor.plot.grid().plot_mpl()`` is the same
        call, and ``vor.plot()`` still draws the grid, so the shorthand people
        actually used keeps working. Replaces eleven plotting aliases.
        """

        return GridPlots(self)

    def get_vor_cells_as_series(
        self,
        overlapping_geometry: shp.Polygon | shp.Point | gpd.GeoSeries | Path = None,
        predicate: str | None = None,
        return_dict: bool = False,
        name_field: str = 'ExploName',
        layer: str | int | None = None,
    ) -> pd.Series | dict:
        """Map overlapping geometry to intersected Voronoi cell IDs.

        ``layer`` picks one layer out of a multi-layer vector file when
        ``overlapping_geometry`` is a path; ``None`` reads the first, as before.
        Passing it alongside an in-memory geometry raises, since there is then no
        file to pick a layer from.

        :func:`~myflopy.modflow.mf6.grid.selection.get_vor_cells_as_series` on this
        grid's :attr:`gdf_vorPolys`. A query carrying its own CRS is reprojected to
        the grid's; one carrying none is read in the grid's. A feature that reaches
        no cell is left out of the answer -- named in a logged warning under the
        default ``predicate=None``, silently under a geopandas predicate.

        Parameters
        ----------
        overlapping_geometry : shapely geometry, GeoDataFrame, GeoSeries, str or Path
            The features to look up: any shapely geometry (one feature), a
            GeoDataFrame or GeoSeries (one feature per row), or a path to a vector
            file. Required -- the ``None`` default, or any other type, raises
            ``ValueError``.
        predicate : str, optional
            Default None: select by MEASURE -- a polygon reaches the cells it
            covers by area, a line those it runs through, a point the one cell
            holding it (the lowest-numbered on a shared face or vertex). A geopandas
            predicate (``"intersects"``, ``"contains"``, ...) instead keeps the raw
            spatial join of cells to features, which takes any touch; it is asked of
            the CELL first, so ``"contains"`` means cells containing the feature.
        return_dict : bool, optional
            Default False: return a Series. True: the same mapping as a ``dict``.
        name_field : str, optional
            The column of a GeoDataFrame or file whose values key the answer.
            Default ``"ExploName"``. When the column is absent -- and always for a
            bare geometry or GeoSeries -- each feature is keyed by its own index
            label (``0`` for a single geometry). A misspelled name therefore falls
            back to index labels without complaint.
        layer : str or int, optional
            The layer to read (name, or zero-based position) when
            ``overlapping_geometry`` is a path to a multi-layer file. Default None:
            the file's first layer. Given with anything but a path, raises
            ``ValueError``.

        Returns
        -------
        pandas.Series or dict
            Feature key -> list of cell ids (a Series named ``"cells"`` unless
            ``return_dict``).
        """
        return selection_get_vor_cells_as_series(
            self.gdf_vorPolys,
            overlapping_geometry=overlapping_geometry,
            predicate=predicate,
            return_dict=return_dict,
            name_field=name_field,
            layer=layer,
        )

    def get_vor_cells_as_dict(
        self,
        locs: Path,
        crs: str = None,
        predicate: str | None = None,
        loc_name_field: str = None,
        return_gdf: bool = False,
        layer: str | int | None = None,
    ) -> dict:
        """Return intersected Voronoi cell IDs keyed by feature name.

        ``layer`` names which layer to read out of a multi-layer file such as a
        GeoPackage. ``None`` keeps GDAL's default -- the first layer -- so a call
        written before this parameter existed behaves identically. A wrong name
        raises listing the layers the file actually holds.

        The cells are chosen by :meth:`get_vor_cells_as_series`'s rule; a feature
        reaching no cell is left out (named in a logged warning only under the
        default ``predicate=None``).

        Parameters
        ----------
        locs : str or pathlib.Path
            A vector file of features to look up.
        crs : str, int or pyproj.CRS, optional
            The CRS the file is reprojected to before the lookup, and the CRS of the
            frame ``return_gdf=True`` returns. The lookup itself is made in the
            grid's CRS whatever this is. Default None: the grid's :attr:`crs`. The
            file must carry a CRS: geopandas refuses to reproject one that does not.
        predicate : str, optional
            Default None: select by measure; a geopandas predicate keeps the raw
            any-touch join. See :meth:`get_vor_cells_as_series`.
        loc_name_field : str, optional
            The column whose values key the answer. Default None: each feature's
            row index (``0``, ``1``, ...). A name the file does not carry falls back
            to the row index without complaint.
        return_gdf : bool, optional
            Default False: return the mapping alone. True: return
            ``(mapping, features)``, the features as read and reprojected to
            ``crs``.
        layer : str or int, optional
            The layer to read (name, or zero-based position) out of a multi-layer
            file. Default None: the first layer.

        Returns
        -------
        dict or tuple
            Feature key -> ``[[cell, ...]]`` -- a list holding ONE list of cell ids
            (the shape ledger 192 kept) -- or that and the features, with
            ``return_gdf``.
        """
        return selection_get_vor_cells_as_dict(
            self,
            locs=locs,
            crs=crs,
            predicate=predicate,
            loc_name_field=loc_name_field,
            return_gdf=return_gdf,
            layer=layer,
        )

    def get_model_boundary_polygons(self) -> dict:
        """Return polygon geometries for the model boundary components."""
        return get_model_boundary_polygons(self.gdf_vorPolys)

    def get_grid_edge(self, idomain: list = None, idomain_path: Path = None, include_interiors: bool = True) -> list:
        """The active cells on the model's edge, ascending and each once.

        ``idomain`` is the list of INACTIVE cell ids (a mask or an MF6 array is
        refused, naming the conversion); ``idomain_path`` a vector file whose
        features cover them, by area. Give one or the other; with neither, the
        grid's own ``idomain_path`` is used when it has one, and otherwise every
        cell is active. ``include_interiors=False`` drops the cells ringing a hole.

        A cell is on the edge when it borders the outside of the mesh or has an
        inactive face-neighbour. See
        :func:`~myflopy.modflow.mf6.grid.selection.get_grid_edge_cells`.

        Parameters
        ----------
        idomain : list of int, optional
            The INACTIVE cell ids (any iterable of integers). ``[]`` means none is
            inactive. Refused with ``ValueError``: a boolean mask, an ``ncpl``-long
            array of only -1/0/1 (an MF6 idomain -- pass
            ``np.flatnonzero(domain <= 0)``), and ids outside ``0..ncpl-1``. A grid
            holds no idomain of its own to fall back on (ledger 388).
        idomain_path : str or pathlib.Path, optional
            A ``.shp``/``.gpkg`` whose features cover the inactive cells, by area,
            reprojected to the grid's CRS; a file covering no cell raises
            ``ValueError``. Default None.
        include_interiors : bool, optional
            Default True: every active cell bordering the outside OR any inactive
            cell, holes included. False: only each active component's OUTER edge --
            the cells ringing a hole inside it are dropped, while an island inside a
            ring keeps its own edge.

        Returns
        -------
        list of int
            The edge cells, ascending, each once.

        Raises
        ------
        ValueError
            Both ``idomain`` and ``idomain_path`` given, or either is refused as
            above.
        """
        return get_grid_edge_cells(
            self,
            idomain=idomain,
            idomain_path=idomain_path,
            include_interiors=include_interiors,
        )

    def generate_grid_coordinates(self, grid_spacing: float) -> list:
        """Generate regularly spaced coordinates over the grid extent.

        Parameters
        ----------
        grid_spacing : float
            The step between points in x and in y, in model length units: any finite
            number above zero (``2.5`` works); zero, a negative or a non-finite value
            raises ``ValueError``. The points start at the vertex extent's minimum,
            exactly, and continue until they reach or pass its maximum (within 1e-9
            of a step), so the last row and column can lie up to one step outside
            the grid.

        Returns
        -------
        tuple of (list of float, list of float)
            ``(xs, ys)``, one entry per point, x-major.
        """
        return geometry_generate_grid_coordinates(self, grid_spacing)

    @property
    def centroids(self):
        """Return cached x/y centroid coordinate arrays."""
        if self._centroids is None:
            self._centroids = self.get_centroids()
        return self._centroids

    def get_centroids(self) -> tuple[list[float], list[float]]:
        """Compute x and y centroid coordinates for all Voronoi cells."""
        return geometry_get_centroids(self)

    def find_adjacent_polygons(self, gdf: gpd.GeoDataFrame) -> list:
        """Return adjacency lists for the polygons in ``gdf``.

        Two polygons are adjacent when their intersection is a single
        ``LineString`` -- a shared edge. A corner touch, and a pair sharing two
        separate segments (a ``MultiLineString``), are not. The grid itself is not
        read; :attr:`adjacent_cells_idx` is this over :attr:`gdf_vorPolys`, cached.

        Parameters
        ----------
        gdf : geopandas.GeoDataFrame
            The polygons, e.g. :attr:`gdf_vorPolys`. Must carry a default
            ``RangeIndex`` (``0..n-1``): index labels are used as list positions.

        Returns
        -------
        list of list of int
            Per polygon (by position), the positions of its neighbours, in
            spatial-index order rather than sorted.
        """
        return geometry_find_adjacent_polygons(gdf)

    def find_adjacent_cells(self, cell_id):
        """Return neighboring cell IDs for one Voronoi cell.

        Read from the DISU connection table (:attr:`ia`/:attr:`iac`/:attr:`ja`,
        computed by :meth:`get_disu_connectivity` on first use), so a neighbour is a
        cell sharing a face of non-zero width -- MF6's rule -- and a corner touch is
        not one.

        Parameters
        ----------
        cell_id : int
            The zero-based cell. Any integer-valued number is accepted (a numpy
            integer, or a float id from a spatial join), via ``int()``.

        Returns
        -------
        list of int
            The neighbouring cells, in the connection table's order (ascending).
        """
        return geometry_find_adjacent_cells(self, cell_id)

    def calculate_distance(self, gdf, poly_idx1, poly_idx2):
        """Distance from one polygon's CENTROID to the nearest point it shares with another.

        :func:`~myflopy.modflow.mf6.grid.geometry.calculate_distance` on ``gdf``; the
        grid itself is not read. NOT MODFLOW 6's connection length, which is measured
        from the DISV centre perpendicular to the line of the shared face: that is
        :attr:`cl12` (ledger 328), and the two part wherever the centre is not the
        centroid.

        Parameters
        ----------
        gdf : geopandas.GeoDataFrame
            The polygons, e.g. :attr:`gdf_vorPolys`. Rows are read by POSITION
            (``iloc``), never by index label.
        poly_idx1 : int
            Position of the polygon measured FROM: the distance starts at its centroid.
        poly_idx2 : int
            Position of the other polygon.

        Returns
        -------
        float
            In the polygons' length units, per side; ``nan`` when the polygons share
            no point.

        Raises
        ------
        IndexError
            A position outside ``gdf``.

        See Also
        --------
        cl12 : MF6's connection length, per side, for every face.

        Examples
        --------
        >>> import numpy as np
        >>> from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
        >>> verts = np.array([[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]], float)
        >>> grid = VoronoiGridPlus(verts=verts, iverts=[[0, 3, 4, 1], [1, 4, 5, 2]],
        ...                        xcyc=np.array([[0.25, 0.5], [1.5, 0.5]]))
        >>> grid.calculate_distance(grid.gdf_vorPolys, 0, 1)   # from the centroid
        0.5
        >>> float(grid.get_disu_connectivity()[2][1])          # MF6's: from (0.25, 0.5)
        0.75
        """
        return geometry_calculate_distance(gdf, poly_idx1, poly_idx2)

    def shared_face_length(self, poly1, poly2):
        """Return the length of the shared face between two polygons.

        Plain geometry on the two polygons given; the grid itself is not read.

        Parameters
        ----------
        poly1 : shapely.Polygon
            One polygon, e.g. ``grid.gdf_vorPolys.geometry[5]``.
        poly2 : shapely.Polygon
            The other.

        Returns
        -------
        float or int
            The length of their intersection, in the polygons' length units, when
            that intersection is a single ``LineString``; the integer ``0``
            otherwise -- no contact, a corner touch, or two separate shared segments
            (a ``MultiLineString``).
        """
        return geometry_shared_face_length(poly1, poly2)

    def shared_face(self, cell_a: int, cell_b: int):
        """The ``LineString`` two cells share, or ``None`` when they only touch.

        The segment itself, where :meth:`shared_face_length` gives only its
        length -- a horizontal flow barrier needs the geometry.

        Parameters
        ----------
        cell_a : int
            One zero-based cell, by position in :attr:`gdf_vorPolys`.
        cell_b : int
            The other.

        Returns
        -------
        shapely.LineString or None
            The shared segment (the longest, when the two boundaries meet in
            several); ``None`` when the cells share no edge longer than 1e-9 model
            length units -- a corner touch, no contact, or a degenerate ring GEOS
            refuses.
        """

        from myflopy.modflow.mf6.grid.barriers import shared_face_segment

        return shared_face_segment(self, cell_a, cell_b)

    def barrier_faces(self, barrier):
        """The cell pairs a barrier trace cuts -- the connections it blocks.

        See :func:`myflopy.modflow.mf6.grid.barriers.barrier_faces`. This is what
        turns a digitised fault trace into the pairs MODFLOW 6's HFB wants.

        A connection is cut when the segment joining its two cell centres (the DISV
        ``cell2d`` points MODFLOW 6 connects) crosses the trace; a centre lying ON
        the trace counts as on its left, looking along the trace as drawn.

        Parameters
        ----------
        barrier : shapely.LineString, shapely.MultiLineString, GeoSeries or GeoDataFrame
            The trace. A GeoSeries/GeoDataFrame is unioned into one geometry and
            reprojected to the grid's CRS when both declare one; a bare geometry is
            read in the grid's CRS. An empty geometry cuts nothing (``[]``, with a
            logged warning). Not a file path.

        Returns
        -------
        list of tuple of (int, int)
            ``(low, high)`` cell pairs, sorted and unique.
        """

        from myflopy.modflow.mf6.grid.barriers import barrier_faces

        return barrier_faces(self, barrier)

    def enclosed_faces(self, polygon):
        """The watertight cut sealing the cells inside ``polygon`` from the rest.

        See :func:`myflopy.modflow.mf6.grid.barriers.enclosed_faces` -- and note
        this is NOT the faces the polygon's ring crosses, which is not watertight.

        Parameters
        ----------
        polygon : shapely.Polygon, shapely.MultiPolygon, GeoSeries or GeoDataFrame
            The region to seal. A cell is inside when its polygon CENTROID lies
            strictly inside it. A GeoSeries/GeoDataFrame is unioned and reprojected
            to the grid's CRS when both declare one; a bare geometry is read in the
            grid's CRS. Not a file path.

        Returns
        -------
        tuple of (list of tuple of (int, int), list of int)
            ``(faces, interior)``: every ``(low, high)`` pair with one cell inside
            and one outside, sorted; and the enclosed cells, sorted. ``([], [])``,
            with a logged warning, when no centroid falls inside.
        """

        from myflopy.modflow.mf6.grid.barriers import enclosed_faces

        return enclosed_faces(self, polygon)

    def get_gdf_vorPolys(self, crs=None):
        """Build a GeoDataFrame of Voronoi polygons.

        One polygon per cell, from ``verts``/``iverts``, with a default
        ``RangeIndex``. **Side effect:** the frame REPLACES the cached
        :attr:`gdf_vorPolys` (and ``x_coords_by_node``/``y_coords_by_node``).

        Parameters
        ----------
        crs : str, int or pyproj.CRS, optional
            Default None: the grid's :attr:`crs`, the label the frame carries. Any
            other value must NAME the grid's CRS (``2927`` for ``"EPSG:2927"``
            does): since the frame becomes the grid's cache, another CRS raises
            ``ValueError`` naming ``grid.crs = ...`` to relabel the grid and
            ``grid.gdf_vorPolys.to_crs(...)`` to reproject -- it used to leave the
            cache labelled differently from the grid (ledger 431).

        Returns
        -------
        geopandas.GeoDataFrame
            A ``geometry`` column of cell polygons, rows in cell order.
        """
        return geometry_get_gdf_vor_polys(self, crs=crs)

    def get_gdf_topbtm(self, rasters: list, labels: list = None):
        """Backward-compatible alias for ``get_gdf_topbtm_multilyr``.

        Parameters
        ----------
        rasters : list of (str or os.PathLike or float), or str or os.PathLike
            As for :meth:`get_gdf_topbtm_multilyr`.
        labels : list of str or int, optional
            As for :meth:`get_gdf_topbtm_multilyr`.
        """
        return self.get_gdf_topbtm_multilyr(rasters=rasters, labels=labels)

    def get_raster_vals_at_centroids(
        self,
        raster_files: str | Path | list[str | Path | float] = None,
        labels: str | list[str] = None,
        *,
        method: str = "area",
        nodata: float | None = None,
        gaps: str = "warn",
    ) -> gpd.GeoDataFrame | None:
        """Sample raster values at Voronoi cell centroids.

        ``method``, ``nodata`` and ``gaps`` are
        :func:`~myflopy.modflow.mf6.grid.surfaces.get_raster_vals_at_centroids`'s;
        the first two were unreachable from here, which also put ``gaps`` --
        what to do about a cell the raster gave no value for -- out of reach
        of the scope a user samples from.

        Parameters
        ----------
        raster_files : str or os.PathLike, or list of (str or os.PathLike or float)
            What to sample: one raster path, or a list whose items are raster paths
            or numbers (a number, numpy's included, becomes a constant column). A
            path is a ``str`` or any ``os.PathLike`` such as ``pathlib.Path``. A
            raster in another CRS is sampled by reprojecting the cells to it. Any
            other item -- a ``bool``, ``None``, a nested list -- raises
            ``TypeError`` naming it (ledger 414). Default None: logs a warning
            and returns ``None``.
        labels : str, or list of str or int, optional
            Column names, one per item of ``raster_files`` (a bare ``str`` names the
            one column). Default None: ``0``, ``1``, ... A length mismatch raises
            ``ValueError``.
        method : {"area", "centroid"}, optional
            Default ``"area"``: the mean of the raster pixels inside each cell
            polygon, a cell holding no valid pixel falling back to the pixel at its
            centroid (warned as slow above 50,000 cells). ``"centroid"``: the one
            pixel at each cell centroid. Anything else raises ``ValueError``.
        nodata : float, optional
            A value to treat as no-data, overriding each raster's declared one, for
            files that use a sentinel without recording it; applies to every raster
            in the call. Default None: each file's own declared value.
        gaps : {"warn", "nearest", "raise", "ignore"}, optional
            What a cell the raster gave no value for becomes. ``"warn"`` (default):
            NaN, with a warning naming the cells. ``"nearest"``: the value of the
            nearest valued cell, by centroid distance. ``"raise"``: ``ValueError``.
            ``"ignore"``: NaN, silently. Unless ``"ignore"``, a raster that gives NO
            cell a value while not overlapping the grid at all raises
            ``ValueError`` printing both extents and both CRSs.

        Returns
        -------
        geopandas.GeoDataFrame or None
            ``geometry`` (the cell CENTROIDS, as points) then one column per label,
            one row per cell; ``None`` when ``raster_files`` is None.
        """
        return surface_get_raster_vals_at_centroids(
            self, raster_files=raster_files, labels=labels,
            method=method, nodata=nodata, gaps=gaps,
        )

    def get_gdf_topbtm_multilyr(self, rasters: list, labels: list = None):
        """Build top/bottom surfaces for one or more layers from rasters.

        :meth:`get_raster_vals_at_centroids` with its defaults (area-weighted,
        ``gaps="warn"``). Returns a new frame; it does not set :attr:`gdf_topbtm`.

        Parameters
        ----------
        rasters : list of (str or os.PathLike or float), or str or os.PathLike
            The model top, then one bottom per layer: raster paths (``str`` or
            ``os.PathLike``) or constant numbers. A bare path is a one-item list.
            Any other item raises ``TypeError`` naming it, as for
            :meth:`get_raster_vals_at_centroids`.
        labels : list of str or int, optional
            One column name per item. Default None: ``0``, ``1``, ... -- the same
            integer labels ``LayerBuildResult.attach_to_grid`` writes.

        Returns
        -------
        geopandas.GeoDataFrame
            ``geometry`` (cell centroids) then one column per surface.
        """
        return surface_get_gdf_topbtm_multilyr(self, rasters=rasters, labels=labels)

    def get_cell_areas(self):
        """Return per-cell polygon areas."""
        return surface_get_cell_areas(self)

    def get_origin_xy(self):
        """Return an origin point suitable for plotting or raster work."""
        return surface_get_origin_xy(self)

    @property
    def gdf_latlon(self):
        """Return the Voronoi polygons reprojected to latitude/longitude."""
        return geometry_get_gdf_latlon(self)

    @property
    def latlon(self):
        """Return centroid coordinates in latitude/longitude."""
        return geometry_get_latlon(self)

    def points_to_latlon(self, x, y):
        """Reproject model-coordinate points to ``(lon, lat)`` for map overlays.

        Parameters
        ----------
        x : array_like of float
            The points' x coordinates, in the grid's :attr:`crs`.
        y : array_like of float
            Their y coordinates, same length as ``x``.

        Returns
        -------
        tuple of (numpy.ndarray, numpy.ndarray)
            ``(lon, lat)`` in ``grid.crs_latlon`` (EPSG:4326).

        Raises
        ------
        ValueError
            The grid has no CRS, so there is nothing to reproject from.
        """
        return geometry_points_to_latlon(self, x, y)

    def get_grid_centroid(self):
        """Return the centroid of the overall grid extent."""
        return geometry_get_grid_centroid(self)

    def get_overlapping_area(self, shp_gpkg=None, cell_list=None):
        """Area of the WHOLE cells a geometry reaches, or of the cells named.

        Not the area of the geometry's intersection with the cells: every cell the
        geometry reaches counts in full (a small polygon inside one cell returns
        that cell's whole area). The cells are unioned, so the result is in model
        length units squared.

        Parameters
        ----------
        shp_gpkg : shapely geometry, GeoDataFrame, GeoSeries, str or Path, optional
            Features whose cells to measure, chosen by
            :meth:`get_vor_cells_as_series`'s default rule (by coverage, not any
            touch). Wins over ``cell_list`` when both are given. A geometry reaching
            no cell gives ``0.0``.
        cell_list : list of int, optional
            Cell ids (index labels of :attr:`gdf_vorPolys`) to measure instead.

        Returns
        -------
        float
            The area of the union of those cells.

        Raises
        ------
        ValueError
            Neither ``shp_gpkg`` nor ``cell_list`` was given.
        """
        return geometry_get_overlapping_area(self, shp_gpkg=shp_gpkg, cell_list=cell_list)

    def get_vor_idx_from_geometry(
        self,
        shp_to_query: shp = None,
        gdf_to_query: gpd.GeoDataFrame = None,
        crs: str = None,
        crs_latlon: str = "EPSG:4326",
        name_col: str = None,
        predicate: str | None = None,
    ) -> dict:
        """The cells each feature of a query layer COVERS, keyed by feature.

        By measure, through the shared rule: a polygon reaches the cells it covers
        by area, a line those it runs through, a point the one cell holding it. An
        explicit ``predicate=`` keeps the raw any-touch answer this used to
        default to, asked of the cell first as the shared selector asks it. Cells
        are found in THIS grid's frame -- ``crs`` names only the frame of a query
        that carries none -- and ``crs_latlon`` has no job left; see
        :func:`~myflopy.modflow.mf6.grid.geometry.get_vor_idx_from_geometry`.

        Parameters
        ----------
        shp_to_query : str or os.PathLike, or list of them, optional
            Vector files to read (``.shp``, or ``.gpkg`` -- first layer only),
            concatenated; they must share one CRS. A path is a ``str`` or any
            ``os.PathLike``. Wins over ``gdf_to_query``.
        gdf_to_query : geopandas.GeoDataFrame, optional
            The features, when not read from files. A GeoSeries, or a plain
            ``pandas.DataFrame`` with a ``geometry`` column, also works. One of the
            two is required, else ``ValueError``.
        crs : str, int or pyproj.CRS, optional
            The CRS to read a query that carries NONE in. Default None: the grid's
            :attr:`crs`. Ignored, with a logged warning, when the query carries its
            own CRS (it is then reprojected to the grid's).
        crs_latlon : str, int or pyproj.CRS, optional
            Retired; kept so existing calls still run. Any spelling of EPSG:4326 (the
            default) or ``None`` is accepted and has no effect; any other CRS raises
            ``ValueError``, and a value pyproj cannot read raises pyproj's
            ``CRSError`` (a ``RuntimeError``, measured).
        name_col : str, optional
            The query column whose values key the answer. Default None: each
            feature's own index label.
        predicate : str, optional
            Default None: select by measure (a polygon reaches the cells it covers by
            area, a line those it runs through, a point the one cell holding it). A
            geopandas predicate keeps the raw any-touch answer, asked of the cell
            first (``"contains"`` = cells containing the feature).

        Returns
        -------
        dict
            Feature key -> sorted list of cell ids. A feature reaching no cell keeps
            its key with an empty list and is named in a logged warning.
        """
        return geometry_get_vor_idx_from_geometry(
            self,
            shp_to_query=shp_to_query,
            gdf_to_query=gdf_to_query,
            crs=crs,
            crs_latlon=crs_latlon,
            name_col=name_col,
            predicate=predicate,
        )

    def get_vor_idx_from_geometry_idx(
        self,
        gdf_to_query: gpd.GeoDataFrame = None,
        idx: int = 0,
        predicate: str | None = None,
    ) -> list:
        """The cells ONE row of a query frame covers, flat and ascending.

        The same measured rule as :meth:`get_vor_idx_from_geometry`, which this
        asks for that single row; ``idx`` is positional.

        Parameters
        ----------
        gdf_to_query : geopandas.GeoDataFrame
            The features. A GeoSeries, or a plain ``pandas.DataFrame`` with a
            ``geometry`` column, also works. A query carrying its own CRS is
            reprojected to the grid's; one carrying none is read in the grid's.
            Required: the ``None`` default, and a frame with no geometry, raise
            ``ValueError``.
        idx : int, optional
            The POSITION of the row to look up (``gdf_to_query.iloc[[idx]]``), not
            its index label. Default ``0``, the first row. Out of range raises
            ``IndexError``.
        predicate : str, optional
            Default None: select by measure (a polygon reaches the cells it covers by
            area, a line those it runs through, a point the one cell holding it). A
            geopandas predicate keeps the raw any-touch answer, asked of the cell
            first (``"contains"`` = cells containing the feature).

        Returns
        -------
        list of int
            The cells that row reaches, ascending, each once; empty (with a logged
            warning) when it reaches none.
        """
        return geometry_get_vor_idx_from_geometry_idx(
            self,
            gdf_to_query=gdf_to_query,
            idx=idx,
            predicate=predicate,
        )

    def set_k_vor(self, k_dict: dict = None, k_default: float = 100.0) -> list[float]:
        """Build a per-cell hydraulic-conductivity list from region mappings.

        **Side effect:** writes the values into a ``"Kh"`` column of the cached
        :attr:`gdf_vorPolys` (replacing one already there), so they also reach
        anything that frame is written to afterwards, such as :meth:`to_shapefile`.

        Parameters
        ----------
        k_dict : dict of float to list of int, optional
            ``{K value: cell ids}``: every cell id in a value's list (row labels of
            :attr:`gdf_vorPolys`) gets that key as its K. Each value must be
            list-like. A cell named under several keys takes the LAST one in the
            dict's order. Each key is read as a ``float``. An id the grid does not
            have raises ``ValueError`` naming it, before anything is written.
            Default None: no cell is named, so every cell takes ``k_default``.
        k_default : float, optional
            The K of every cell no list names, read as a ``float``. Default
            ``100.0``. No unit is assumed: use the model's length/time units.

        Returns
        -------
        list of float
            One K per cell, in :attr:`gdf_vorPolys` row order. Always floats: the
            ``"Kh"`` column is float whatever was passed, so a float K never lands
            in an int column (pandas warns on that upcast today and is deprecated to
            raise).
        """
        return geometry_set_k_vor(self, k_dict=k_dict, k_default=k_default)

    def get_raster_from_strike_dip(
        self,
        strike: float,
        dip: float,
        known_point: tuple,
        pixel_size: float = 1,
        output_filename: str | Path | None = None,
    ):
        """Write a strike/dip plane as a raster over the grid, then sample it onto the cells.

        The GeoTIFF covers the grid's domain snapped outward to whole pixels (one
        float32 band, in the grid's :attr:`crs`) and is then read back by
        :meth:`get_raster_vals_at_centroids` with its defaults (area-weighted).

        **The plane is the geological strike/dip plane** (ledger 362).
        Azimuths are degrees clockwise from north (+y), and the right-hand rule
        holds: the plane dips toward ``strike + 90``. It passes through
        ``known_point`` at its elevation, is level along the strike line, and
        descends at ``tan(dip)`` per unit of horizontal distance in the dip
        direction -- ``strike=0, dip=45`` falls 1 ft per ft to the east. It is
        built on :meth:`get_normal_from_strike_and_dip`. Each pixel holds the plane
        at the pixel's CENTRE (ledger 314's registration), so the sampled cells
        carry no half-pixel offset.

        Parameters
        ----------
        strike : float
            Strike azimuth, in degrees clockwise from north (+y): the horizontal
            line along which the plane is level. Any finite value (taken modulo
            360); a non-number raises ``TypeError``.
        dip : float
            Dip, in degrees below the horizontal, toward ``strike + 90``. Must lie
            in ``[0, 90)``: a negative dip raises ``ValueError`` naming the strike
            turned 180 degrees that spells it, and ``90`` or more raises (a vertical
            plane has no elevation as a function of x and y). ``0`` is level at
            the known elevation.
        known_point : tuple of float
            ``(x, y, elevation)`` the plane passes through: ``x``/``y`` in the grid's
            coordinates, ``elevation`` in the units the raster should carry. Must
            unpack to exactly three values.
        pixel_size : float, optional
            The raster's pixel side, in model length units; the domain's bounds are
            snapped outward to multiples of it. Default ``1``. The plane is
            evaluated at each pixel's centre. Anything but a positive finite number
            raises ``ValueError`` naming it.
        output_filename : str or pathlib.Path, optional
            Where the GeoTIFF is written, overwriting any file there. Default
            ``None``: ``raster.tif`` in the working directory at the time of the
            CALL.

        Returns
        -------
        geopandas.GeoDataFrame
            ``geometry`` (the cell centroids, as points) and ``elev``, one row per
            cell.
        """
        return surface_get_raster_from_strike_dip(
            self,
            strike=strike,
            dip=dip,
            known_point=known_point,
            pixel_size=pixel_size,
            output_filename=output_filename,
        )

    def to_shapefile(self, filepath: str | Path = 'vor_shp.shp'):
        """Write the Voronoi polygons to a shapefile or other GeoPandas target.

        Writes :attr:`gdf_vorPolys` as it stands -- the cell polygons, in the grid's
        CRS, plus any column added to that cached frame (``"Kh"`` after
        :meth:`set_k_vor`). For per-cell values or layer surfaces, see
        :meth:`to_gpkg`.

        Parameters
        ----------
        filepath : str or pathlib.Path, optional
            The destination, passed to ``GeoDataFrame.to_file``, which picks the
            format from the extension (``.shp``, ``.gpkg``, ...). Default
            ``"vor_shp.shp"``, relative to the working directory.

        Returns
        -------
        None
        """
        return self.gdf_vorPolys.to_file(filepath)

    def to_gpkg(
        self,
        path: str | Path = "grid.gpkg",
        *,
        data=None,
        surfaces: bool = False,
        layer: str = "cells",
    ):
        """Write the grid's cells to a GeoPackage, ready for QGIS.

        The vector counterpart to :meth:`to_xugrid`. One ``cells`` layer of flat
        polygons carrying whatever per-cell values you pass, plus -- with
        ``surfaces=True`` -- one layer PER LAYER SURFACE whose geometry carries
        that surface's elevation as its **Z**.

        The Z is the point. A QGIS elevation profile reads a vector layer's own
        Z when its Elevation tab is set to ``Clamping: Absolute``, and that is
        the only setting needed: no expression, no data-defined override. The
        alternative, driving the elevation from an attribute, needs three
        settings and an expression per layer. Measured against QGIS 4.2.2: both
        produce the same profile, and a MESH layer produces a flat line at zero
        however it is configured, because mesh elevation comes from vertex Z and
        a UGRID file carries none (ledger 261).

        One GeoPackage layer per surface, rather than one layer with a surface
        column, because elevation in QGIS is a per-LAYER property: several
        surfaces on one profile means several map layers.

        Parameters
        ----------
        path : str or pathlib.Path, optional
            Destination GeoPackage; missing parent directories are created. Created
            if absent, added to if present: a layer of the same name is replaced
            and the file's other layers are kept (measured). Default
            ``"grid.gpkg"``, relative to the working directory.
        data : dict of str to array-like, optional
            ``{name: per-cell values}`` written as float attribute columns on the
            flat cell layer, beside a ``cell`` id column. Each array is flattened
            and must then hold exactly ``ncpl`` values (so a ``(nlay, ncpl)``
            array is refused), else ``ValueError``. With ``surfaces=True`` a name
            matching a surface (``top``, ``botm_1``, ...) is overwritten on that
            layer by the surface. Default None: the ``cell`` column alone.
        surfaces : bool, optional
            Also write ``top`` and ``botm_1`` ... ``botm_n`` as columns of the
            flat layer AND as one Z-valued layer each, from
            :meth:`layer_surface_fields`. If the grid carries no layer surfaces
            that raises ``ValueError`` before the GeoPackage is written. Default
            False.
        layer : str, optional
            Name for the flat cell layer (default ``"cells"``). Converted with
            ``str()``.

        Returns
        -------
        pathlib.Path
            The GeoPackage written.

        See Also
        --------
        to_xugrid : the same data as a UGRID mesh, for NetCDF interchange.
        layer_surface_fields : the elevations this writes, as plain arrays.

        Examples
        --------
        >>> vor.to_gpkg("grid.gpkg", surfaces=True)
        >>> # then in QGIS: add the `top` / `botm_1` layers, and set
        >>> # Properties > Elevation > Clamping: Absolute on each.

        >>> vor.to_gpkg("grid.gpkg", data={"k": npf_k_layer0,
        ...                                "head": model.hds.array(layer=0)},
        ...             surfaces=True)
        """

        from shapely import force_3d

        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        ncpl = int(self.ncpl)
        geometry = self.gdf_vorPolys.geometry
        columns = {"cell": np.arange(ncpl, dtype=int)}
        for name, values in (data or {}).items():
            array = np.asarray(values, dtype=float).ravel()
            if array.shape[0] != ncpl:
                raise ValueError(
                    f"to_gpkg(): '{name}' has {array.shape[0]} value(s), "
                    f"expected ncpl={ncpl}."
                )
            columns[str(name)] = array

        # Resolved BEFORE the first write, so a grid with no surfaces fails with
        # nothing written rather than leaving a half-built GeoPackage behind.
        elevations = self.layer_surface_fields() if surfaces else {}
        if surfaces:
            columns.update(elevations)

        gpd.GeoDataFrame(columns, geometry=geometry.to_numpy(), crs=self.crs).to_file(
            path, layer=str(layer), driver="GPKG"
        )
        for name, values in elevations.items():
            gpd.GeoDataFrame(
                {"cell": np.arange(ncpl, dtype=int), name: values},
                geometry=[force_3d(geom, float(z)) for geom, z in zip(geometry, values, strict=True)],
                crs=self.crs,
            ).to_file(path, layer=name, driver="GPKG")
        return path

    def ugrid2d(self):
        """Return this grid's topology as an :class:`xugrid.Ugrid2d` (no data).

        Builds the UGRID 2-D unstructured mesh -- nodes (cell vertices), faces
        (cells), and the face-node connectivity -- from :meth:`get_disv_gridprops`,
        carrying the grid :attr:`crs`. This is the shared topology builder behind
        :meth:`to_xugrid` and ``HeadsPlus.to_xugrid`` / ``model.to_xugrid``; call
        it directly when you want to attach your own ``xarray`` data on the face
        dimension (``grid.face_dimension``).

        Returns
        -------
        xugrid.Ugrid2d
            The mesh topology, with ``n_face == ncpl`` and ``n_node == nvert``.

        Raises
        ------
        ImportError
            If the optional ``xugrid`` package is not installed.
        """

        from myflopy._optional import require

        xu = require("xugrid", feature="ugrid2d() unstructured-grid export")

        gp = self.get_disv_gridprops()
        nvert = int(gp["nvert"])
        ncpl = int(gp["ncpl"])

        # Node coordinates, placed at their vertex id so connectivity lines up.
        node_x = np.full(nvert, np.nan, dtype=float)
        node_y = np.full(nvert, np.nan, dtype=float)
        for iv, x, y in gp["vertices"]:
            node_x[int(iv)] = float(x)
            node_y[int(iv)] = float(y)

        # Ragged face -> node connectivity, padded to a rectangular array.
        fill_value = -1
        rows = [
            [int(v) for v in cell[4:4 + int(cell[3])]]
            for cell in gp["cell2d"]
        ]
        max_nodes = max(len(row) for row in rows)
        face_node_connectivity = np.full((ncpl, max_nodes), fill_value, dtype=np.int64)
        for i, row in enumerate(rows):
            face_node_connectivity[i, : len(row)] = row

        return xu.Ugrid2d(
            node_x,
            node_y,
            fill_value,
            face_node_connectivity,
            name="mesh2d",
            is_projected=True,
            crs=self.crs,
        )

    def layer_surface_fields(self) -> dict[str, np.ndarray]:
        """The layer frame as FLAT per-cell arrays: ``top``, ``botm_1`` ... ``botm_n``.

        One variable PER LAYER, deliberately, rather than a single
        ``(nlay, ncpl)`` array. A mesh dataset in QGIS is one value per face
        (optionally per time); a variable carrying an extra ``layer`` dimension
        has no place in that model and QGIS drops it **silently** -- measured
        against QGIS 4.2.2, which listed no elevation datasets at all for a
        layered array and every one of them once flattened. Naming them this way
        is what makes :meth:`to_xugrid` ``surfaces=True`` usable there.

        Read positionally from :attr:`gdf_topbtm`, which is geometry-first by
        contract (ledger 196): column 1 is the model top and column ``2 + n`` the
        bottom of layer ``n``. The same convention the hover and the cross
        sections use, and it survives the two different column NAMINGS that
        ``attach_to_grid()`` and ``LayerStack.to_disv()`` publish.

        Returns
        -------
        dict of str to ndarray
            ``{"top": (ncpl,), "botm_1": (ncpl,), ...}``, one entry per layer.

        Raises
        ------
        ValueError
            If this grid carries no layer surfaces, rather than returning an
            empty mapping -- a caller asking for surfaces that are not there has
            a problem to fix, not a quietly emptier file.

        Examples
        --------
        >>> vor.layer_surface_fields()["botm_1"]        # layer 1's floor, per cell
        >>> vor.to_xugrid(surfaces=True)                # the same, on the mesh
        """

        frame = self.gdf_topbtm
        if frame is None:
            raise ValueError(
                "this grid carries no layer surfaces (gdf_topbtm is None). Attach "
                "them with LayerBuildResult.attach_to_grid(), or build the grid "
                "with rasters= so they are sampled."
            )
        nlay = int(frame.shape[1]) - 2
        if nlay < 1:
            raise ValueError(
                f"gdf_topbtm has {frame.shape[1]} column(s); a geometry, a top and "
                "at least one layer bottom are needed."
            )
        fields = {"top": np.asarray(frame.iloc[:, 1], dtype=float)}
        for index in range(nlay):
            fields[f"botm_{index + 1}"] = np.asarray(frame.iloc[:, 2 + index], dtype=float)
        return fields

    def to_xugrid(
        self,
        data=None,
        *,
        name: str = "data",
        layer_dim: str = "layer",
        surfaces: bool = False,
    ):
        """Export this DISV grid (and optional per-cell data) as an xugrid object.

        Builds a UGRID 2-D unstructured-mesh representation of the Voronoi grid --
        nodes (cell vertices), faces (cells), and the face-node connectivity --
        wrapped in an `xugrid <https://deltares.github.io/xugrid/>`_ object. xugrid
        is xarray extended for unstructured grids, so the result plugs directly
        into xarray-style analysis, unstructured plotting/cross-sections, and a
        UGRID-NetCDF export that QGIS, ParaView, and other tools can open -- which
        makes DISV results far easier to share than raw MODFLOW binaries.

        The grid topology is taken from :meth:`get_disv_gridprops` and the grid's
        :attr:`crs` is carried onto the mesh.

        Parameters
        ----------
        data : array-like or dict of str to array-like, optional
            Optional per-cell values to attach on the face dimension, read as
            floats. One of:

            - ``None`` (default) -- return just the mesh topology;
            - a 1-D array of length ``ncpl`` -- one field, e.g. a head layer
              (``model.hds.array(layer=0)``) or a K array;
            - a 2-D array shaped ``(nlay, ncpl)`` -- a layered field, given the
              extra ``layer_dim`` dimension;
            - a ``dict`` of ``{variable_name: array}`` -- several fields at once
              (each 1-D or 2-D as above).

            Any other shape, or a cell axis other than ``ncpl``, raises
            ``ValueError``.
        name : str, optional
            Variable name used when ``data`` is a single array (ignored for a
            ``dict`` and for ``data=None``). Defaults to ``"data"``.
        layer_dim : str, optional
            Dimension name given to the first axis of 2-D ``(nlay, ncpl)`` inputs;
            unused when every input is 1-D. Default ``"layer"``.
        surfaces : bool, optional
            Default False. True: also attach this grid's layer frame -- ``top`` and ``botm_1`` ...
            ``botm_n``, one FLAT variable per layer
            (:meth:`layer_surface_fields`). This is the form QGIS can read: a
            mesh dataset there is one value per face, so a layered array is
            dropped silently. A ``data`` entry of the same name wins, so you can
            override a surface with your own. Forces a ``UgridDataset`` return,
            since there is more than one variable. Raises if the grid has no
            surfaces attached.

        Returns
        -------
        xugrid.UgridDataArray | xugrid.UgridDataset
            A ``UgridDataArray`` when ``data`` is a single array and
            ``surfaces`` is False, otherwise a ``UgridDataset`` (also for
            ``data=None`` -- topology only).

        Raises
        ------
        ImportError
            If the optional ``xugrid`` / ``xarray`` packages are not installed.
        ValueError
            If an input array's cell axis does not match ``ncpl``.

        Notes
        -----
        ``.ugrid.plot()`` draws one value per cell, so plot only a 1-D
        ``(ncpl,)`` field. For a 2-D ``(nlay, ncpl)`` input, select a layer first
        (``uda.isel(layer=0).ugrid.plot()``); xugrid does not facet over the layer
        dimension. The 2-D array still writes to NetCDF whole.

        Examples
        --------
        >>> uda = vor.to_xugrid(model.hds.array(layer=0), name="head")
        >>> uda.ugrid.plot()                      # unstructured choropleth (1-D field)
        >>> uda.ugrid.to_netcdf("heads.nc")       # UGRID NetCDF for QGIS/ParaView

        >>> # several fields, including a layered (nlay, ncpl) one
        >>> k_2d = np.vstack([k_layer0, k_layer1])
        >>> uds = vor.to_xugrid({"head": heads_2d, "k": k_2d})
        >>> uds.ugrid.to_netcdf("model.nc")       # share; isel a layer to plot

        >>> # the layer frame, flat, so QGIS lists every elevation as a dataset
        >>> vor.to_xugrid(surfaces=True).ugrid.to_netcdf("grid.nc")
        >>> vor.to_xugrid(model.hds.array(layer=0), name="head", surfaces=True)
        """

        try:
            import xarray as xr
            import xugrid as xu
        except ImportError as err:  # pragma: no cover - optional dependency
            raise ImportError(
                "to_xugrid() requires the optional 'xugrid' and 'xarray' "
                "packages. Install them with `pip install xugrid xarray`."
            ) from err

        grid = self.ugrid2d()
        ncpl = int(self.ncpl)
        face_dim = grid.face_dimension

        def _face_dataarray(values, varname):
            """Wrap a per-cell (face) array as an xugrid ``DataArray`` (validating its length)."""

            array = np.asarray(values, dtype=float)
            if array.ndim == 1:
                if array.shape[0] != ncpl:
                    raise ValueError(
                        f"to_xugrid(): '{varname}' has length {array.shape[0]}, "
                        f"expected ncpl={ncpl}."
                    )
                dims = (face_dim,)
            elif array.ndim == 2:
                if array.shape[-1] != ncpl:
                    raise ValueError(
                        f"to_xugrid(): '{varname}' last axis is {array.shape[-1]}, "
                        f"expected ncpl={ncpl}."
                    )
                dims = (layer_dim, face_dim)
            else:
                raise ValueError(
                    f"to_xugrid(): '{varname}' must be 1-D (ncpl,) or 2-D "
                    f"(nlay, ncpl); got {array.ndim} dimensions."
                )
            return xr.DataArray(array, dims=dims, name=varname)

        # Surfaces FIRST, so an explicit `data` entry of the same name wins: what
        # the caller handed over beats what the grid derived.
        fields: dict = dict(self.layer_surface_fields()) if surfaces else {}
        if isinstance(data, dict):
            fields.update(data)
        elif data is not None:
            fields[name] = data

        if not fields:
            return xu.UgridDataset(grids=[grid])
        if surfaces or isinstance(data, dict):
            dataset = xr.Dataset(
                {key: _face_dataarray(values, key) for key, values in fields.items()}
            )
            return xu.UgridDataset(dataset, grids=[grid])
        return xu.UgridDataArray(_face_dataarray(data, name), grid)

    def get_domain(self):
        """Return the overall polygonal domain of the Voronoi grid."""
        return surface_get_domain(self)

    @staticmethod
    def get_normal_from_strike_and_dip(strike: float, dip: float) -> np.ndarray:
        """Return the upward unit normal of the plane with this strike and dip.

        Convention (ledger 362): azimuths in degrees clockwise from north
        (+y); the plane dips toward ``strike + 90`` (the right-hand rule) and is
        level along the strike. The normal is
        ``[sin(dip) * cos(strike), -sin(dip) * sin(strike), cos(dip)]`` (x, y, z):
        perpendicular to the strike line and to the dip line, unit length, with
        ``z >= 0`` -- a level plane (``dip=0``) gives ``[0, 0, 1]``.

        Parameters
        ----------
        strike : float
            Strike azimuth in degrees clockwise from north (+y). Any finite value
            (taken modulo 360); a non-number raises ``TypeError``.
        dip : float
            Dip in degrees below the horizontal, toward ``strike + 90``, in
            ``[0, 90)``; anything else raises ``ValueError`` naming it (90 is a
            vertical plane; a negative dip is the strike turned 180 degrees).

        Returns
        -------
        numpy.ndarray
            Shape ``(3,)``, unit length, upward (``z > 0``).
        """
        return surface_get_normal_from_strike_and_dip(strike, dip)

    @staticmethod
    def generate_grid_around_point(
        center_point: shp.Point, spacing: float, size: float, crs: str
    ) -> gpd.GeoSeries:
        """Generate a square lattice of points around a center point.

        The x values are ``cx + spacing * (i - size / 2)`` for ``i = 0, 1, ...`` up
        to ``cx + size / 2 * spacing`` -- reached when ``size`` is a whole number
        (to 1e-9 of a step), and otherwise the last whole step short of it; y
        likewise. One point per (x, y) pair, so a whole ``size`` gives
        ``(size + 1) ** 2`` points, and an even one puts a point on the centre.

        Parameters
        ----------
        center_point : shapely.geometry.Point
            The lattice's centre; only its ``x``/``y`` are read.
        spacing : float
            The step between points in x and in y, in the point's length units. Any
            finite number above zero (``2.5`` works). Zero, a negative or a
            non-finite value raises ``ValueError``.
        size : float
            The lattice's width in steps: it spans ``size * spacing`` along each
            axis, half on either side of the centre. Zero gives the centre alone; a
            negative or non-finite value raises ``ValueError``.
        crs : str, int or pyproj.CRS
            The CRS stamped on the result (``GeoSeries.set_crs``); the coordinates
            are not reprojected. Required.

        Returns
        -------
        geopandas.GeoSeries
            The points, x-major, with ``crs`` set.
        """
        return geometry_generate_grid_around_point(center_point, spacing, size, crs)

    @staticmethod
    def generate_grid_polygons(
        center_point: shp.Point, spacing: float, size: float, gap: float
    ) -> shp.MultiPolygon:
        """Generate a square lattice of boxes around a center point.

        The lattice is :meth:`generate_grid_around_point`'s; one box is made per
        lattice cell, each ``spacing - gap`` wide and inset by ``gap / 2`` from its
        cell's edges, so neighbouring boxes are ``gap`` apart.

        Parameters
        ----------
        center_point : shapely.geometry.Point
            The lattice's centre; only its ``x``/``y`` are read.
        spacing : float
            The lattice step in x and in y, in the point's length units: any finite
            number above zero, as for :meth:`generate_grid_around_point`.
        size : float
            The lattice's width in steps, as for :meth:`generate_grid_around_point`:
            a whole ``size`` gives ``size ** 2`` boxes. A size under ``1`` holds no
            box and raises ``ValueError``.
        gap : float
            The space left between neighbouring boxes, in the same units. ``0``
            makes them touch. A negative or non-finite gap, or one of ``spacing``
            or more (no box would be left), raises ``ValueError``.

        Returns
        -------
        shapely.geometry.MultiPolygon
            The boxes, with no CRS.
        """
        return geometry_generate_grid_polygons(center_point, spacing, size, gap)

    @staticmethod
    def voronoi_refine_by_point(point: shp.Point, spacing: int, tri: Triangle) -> gpd.GeoDataFrame:
        """Add a 3x3 pattern of squares around a point to a triangulation.

        Nine squares of side ``spacing``: one with ``point`` at its upper-left
        corner (it extends +x and -y from it), and eight copies of it shifted by
        ``2 * spacing`` east, west, north, south and diagonally, so neighbouring
        squares are ``spacing`` apart. **Side effect:** each is passed to
        ``tri.add_polygon``.

        Parameters
        ----------
        point : shapely.geometry.Point
            The upper-left corner of the central square; only ``x``/``y`` are read.
        spacing : float
            The squares' side, and the gap between them, in model length units.
        tri : flopy.utils.triangle.Triangle
            The triangulation to add the squares to. On a
            :class:`~myflopy.modflow.mf6.grid.triangle.TriangleGrid` with a domain
            already set, each becomes a refinement region with that method's
            defaults (maximum triangle area ``1`` square length unit, clipped to the
            domain); with no domain set, the FIRST square becomes the domain.

        Returns
        -------
        geopandas.GeoDataFrame
            The nine squares, with no CRS: centre, E, W, N, S, NE, NW, SE, SW.
        """
        return geometry_voronoi_refine_by_point(point, spacing, tri)

    def reconcile_surfaces(self, df: pd.DataFrame = None, min_sep=0.1, trigger_sep=1, which='bottom'):
        """Reconcile top/bottom surfaces so adjacent layers do not overlap.

        Surfaces are read POSITIONALLY, top first: column 0 is the model top and
        each later column the next bottom down; rows are read by position too, so
        any index works and is carried to the result. At a contact, every row whose
        lower surface is at most ``trigger_sep`` below the upper one -- or above
        it -- is SET to exactly ``min_sep`` apart. That can THIN a layer as well as
        open one: measured with the defaults, a bottom 0.5 below the top was moved
        UP to 0.1 below it (``LayerStack`` never lets that reach a unit at or above
        its own ``min_thickness``; this method does what it says to whatever it is
        given). Each contact is handled once, against a neighbour that is already
        final, so every contact ends exactly ``min_sep`` apart or farther than
        ``trigger_sep`` (ledger 377). Returns a new frame; neither ``df`` nor
        the grid's :attr:`gdf_topbtm` is changed.

        Parameters
        ----------
        df : pandas.DataFrame or geopandas.GeoDataFrame, optional
            The surfaces, one column each. A GeoDataFrame's ``geometry`` column is
            dropped; every other value is coerced to a number, anything that is not
            one becoming NaN (a NaN is never treated as conflicting). Any index:
            rows are taken by position and the index is carried to the result.
            Default None: a copy of :attr:`gdf_topbtm`.
        min_sep : float or sequence of float, optional
            The separation set at a conflicting contact, in elevation units. A
            scalar applies to every contact; a sequence gives one per bottom surface
            (one fewer than the frame's columns), else ``ValueError``. ``None``
            means ``0.1``. Default ``0.1``.
        trigger_sep : float or sequence of float, optional
            The separation at or below which a contact counts as conflicting
            (measured: a bottom exactly ``1`` below the top was moved under the
            defaults), shaped as ``min_sep`` (``None`` also means ``0.1``). Default
            ``1``.
        which : {"bottom", "top"}, optional
            Which surface of a conflicting pair moves -- to exactly ``min_sep``
            from the other, which can move it either way (measured with the
            defaults and a gap of 0.5: ``"bottom"`` moved the bottom UP to 0.1
            below the top, ``"top"`` moved the top DOWN to 0.1 above the bottom).
            ``"bottom"`` (default) moves the lower one, so the model top never
            moves; contacts are handled top-down. ``"top"`` moves the upper one,
            model top included, and handles the contacts bottom-up, so a raised
            surface is checked again against the one above it (measured:
            ``[10.0, 9.5, 9.45]`` gives ``[9.65, 9.55, 9.45]``). Anything else
            raises ``ValueError``.

        Returns
        -------
        pandas.DataFrame
            The reconciled surfaces, same columns and index, without geometry.
        """
        return geometry_reconcile_surfaces(
            self,
            df=df,
            min_sep=min_sep,
            trigger_sep=trigger_sep,
            which=which,
        )

    def adjust_cells_by_id(
        self,
        cell_ids: list,
        adjustment: int | float,
        df: pd.DataFrame = None,
        layer: int = 0,
        reconcile: bool = True,
    ):
        """Adjust elevations for selected cells by a fixed amount.

        Works on a copy: neither ``df`` nor the grid's :attr:`gdf_topbtm` is
        changed.

        Parameters
        ----------
        cell_ids : list of int
            The rows to adjust, by index LABEL (cell ids, for :attr:`gdf_topbtm`).
        adjustment : float
            Added to each selected value, in elevation units; negative lowers.
        df : pandas.DataFrame or geopandas.GeoDataFrame, optional
            The surfaces to adjust, copied first. Default None: a copy of
            :attr:`gdf_topbtm`.
        layer : int or str, optional
            The COLUMN LABEL of the surface to adjust, not a layer position. Under
            the integer labels ``rasters=`` and ``LayerBuildResult.attach_to_grid``
            write, ``0`` is the MODEL TOP and ``n`` the bottom of the n-th layer.
            Default ``0``. A label the frame lacks raises ``KeyError``.
        reconcile : bool, optional
            Default True: pass the result through :meth:`reconcile_surfaces` with
            its defaults (``min_sep=0.1``, ``trigger_sep=1``, ``which="bottom"``),
            which returns a plain DataFrame without geometry, keeping the index.
            False: return the adjusted copy as it is.

        Returns
        -------
        pandas.DataFrame or geopandas.GeoDataFrame
            The adjusted surfaces.
        """
        return geometry_adjust_cells_by_id(
            self,
            cell_ids=cell_ids,
            adjustment=adjustment,
            df=df,
            layer=layer,
            reconcile=reconcile,
        )

    def adjust_top_btm_overlaps(
        self,
        elev_df=None,
        shp: Path = None,
        buffer=1,
        layer_bottom_name=1,
        min_sep=None,
        trigger_sep=1,
    ) -> pd.DataFrame:
        """Drop one bottom under selected cells below their neighbours' tops.

        For each cell ``shp`` reaches, the ``layer_bottom_name`` surface is set to
        the LOWEST model top (column label ``0``) among its face neighbours, less
        ``buffer``; the frame is then passed through :meth:`reconcile_surfaces`
        with ``min_sep`` and ``trigger_sep``.

        Parameters
        ----------
        elev_df : pandas.DataFrame, optional
            The surfaces, with the model top under column label ``0`` (else
            ``KeyError``) and indexed by cell id. Copied: the frame you pass, and
            :attr:`gdf_topbtm`, are never written. Default None:
            :meth:`reconcile_surfaces` of :attr:`gdf_topbtm` with ``min_sep`` and
            ``trigger_sep``.
        shp : shapely geometry, GeoDataFrame, GeoSeries, str or pathlib.Path
            The features whose cells to lower, chosen by
            :meth:`get_vor_cells_as_series`'s default rule (by coverage). Required:
            the ``None`` default raises ``ValueError``.
        buffer : float, optional
            How far below the neighbours' lowest top the bottom is set, in elevation
            units. Default ``1``.
        layer_bottom_name : int or str, optional
            The column label of the surface to set. Under integer labels ``1`` is the
            first layer's bottom (the default).
        min_sep : float or sequence of float, optional
            The separation both reconciles set at a conflicting contact -- the one
            building the default ``elev_df`` and the final one -- as
            :meth:`reconcile_surfaces` takes it (``None``: ``0.1``). Default None.
        trigger_sep : float or sequence of float, optional
            The separation at or below which both reconciles count a contact as
            conflicting, as :meth:`reconcile_surfaces` takes it (``None``: ``0.1``).
            Default ``1``.

        Returns
        -------
        pandas.DataFrame
            The reconciled surfaces, without geometry.
        """
        return geometry_adjust_top_btm_overlaps(
            self,
            elev_df=elev_df,
            shp_path=shp,
            buffer=buffer,
            layer_bottom_name=layer_bottom_name,
            min_sep=min_sep,
            trigger_sep=trigger_sep,
        )

    @classmethod
    def vor_from_disu(
        cls,
        disu_path: Path,
        crs: str = DEFAULT_CRS,
        rasters: list | Path = None,
        name: str = 'voronoi_grid',
        idomain_path: Path = None,
    ):
        """Reconstruct a Voronoi-style grid view from a DISU input file.

        Reads an **MF6** ``DISU`` -- its ``VERTICES`` and ``CELL2D`` blocks, through
        FloPy -- and builds the grid by :class:`VoronoiGridPlus`'s array form, so the
        cells are taken exactly as written (never repaired) and each ring is put in
        clockwise order. A MODFLOW-USG ``DISU`` has no coordinates; use
        :meth:`from_gsf` for one.

        Parameters
        ----------
        disu_path : pathlib.Path
            The MF6 ``DISU`` file. Must be a ``pathlib.Path`` (its ``parent`` and
            ``name`` are read): a ``str`` raises ``AttributeError``.
        crs : str, int or pyproj.CRS, optional
            The CRS the file's coordinates are in, stamped on the grid. A bare EPSG
            code (``2927`` or ``"2927"``) is stored as ``"EPSG:2927"``. Default
            ``"EPSG:2927"``; ``None`` raises ``ValueError`` in the constructor. As
            for :class:`VoronoiGridPlus`.
        rasters : list of (str or os.PathLike or float), optional
            Surfaces sampled onto the cells the first time :attr:`gdf_topbtm` is read:
            the model top, then one bottom per layer. Default None. As for
            :class:`VoronoiGridPlus`.
        name : str, optional
            A label stored as ``grid.name``. Default ``"voronoi_grid"``.
        idomain_path : str or os.PathLike, optional
            A vector file whose features cover the INACTIVE cells, stored as
            :attr:`idomain_path` and read by :meth:`get_grid_edge`. A ``str`` or
            ``os.PathLike``, stored as a ``pathlib.Path``; another type raises
            ``TypeError`` (ledger 457). Default None.

        Returns
        -------
        VoronoiGridPlus
            The grid, one cell per ``CELL2D`` row.
        """
        verts, iverts, xcyc = get_griddata_from_disu(disu_path)
        return cls(
            verts=verts,
            iverts=iverts,
            xcyc=xcyc,
            crs=crs,
            name=name,
            rasters=rasters,
            idomain_path=idomain_path,
        )

    @classmethod
    def from_gsf(
        cls,
        gsf_path: Path,
        layer: int = 0,
        crs: str = DEFAULT_CRS,
        rasters: list | Path = None,
        name: str = 'voronoi_grid',
        idomain_path: Path = None,
    ):
        """Build a grid view from a MODFLOW-USG ``.gsf`` grid specification file.

        MODFLOW-USG's ``DISU`` carries no coordinates -- only connectivity and
        geometric measures -- so the ``.gsf`` is where a USG model's cell
        outlines actually live. Use this for a USG model;
        :meth:`vor_from_disu` reads an **MF6** ``DISU``, whose ``VERTICES`` and
        ``CELL2D`` blocks a USG ``DISU`` does not have.

        A ``.gsf`` describes 3-D cells and lists every layer separately, so this
        collapses each cell to its distinct plan corners and takes one layer.
        Passing the file's raw records to :class:`VoronoiGridPlus` instead builds
        polygons that trace their outline twice -- invalid geometry with double
        the true area, which still *renders* but silently loses map hover -- and
        one stacked copy of the grid per layer.

        Parameters
        ----------
        gsf_path : str or pathlib.Path
            Path to the ``.gsf`` file. A file whose first line is not
            ``UNSTRUCTURED``, or whose counts disagree with its records, raises
            ``ValueError``.
        layer : int, optional
            Zero-based layer to take the plan view from -- a position among the
            layer numbers the file uses, sorted. Every layer shares one plan
            geometry, so this selects the cell centers and the node numbers the
            cells correspond to. Default ``0``. Out of range raises ``ValueError``.
        crs : str, int or pyproj.CRS, optional
            The CRS the file's coordinates are in, stamped on the grid. A bare EPSG
            code (``2927`` or ``"2927"``) is stored as ``"EPSG:2927"``. Default
            ``"EPSG:2927"``; ``None`` raises ``ValueError`` in the constructor. As
            for :class:`VoronoiGridPlus`.
        rasters : list of (str or os.PathLike or float), optional
            Surfaces sampled onto the cells the first time :attr:`gdf_topbtm` is read:
            the model top, then one bottom per layer. Default None. As for
            :class:`VoronoiGridPlus`.
        name : str, optional
            A label stored as ``grid.name``. Default ``"voronoi_grid"``.
        idomain_path : str or os.PathLike, optional
            A vector file whose features cover the INACTIVE cells, stored as
            :attr:`idomain_path` and read by :meth:`get_grid_edge`. A ``str`` or
            ``os.PathLike``, stored as a ``pathlib.Path``; another type raises
            ``TypeError`` (ledger 457). Default None.

        Returns
        -------
        VoronoiGridPlus
            A grid whose cells are the plan polygons of the chosen layer.

        Examples
        --------
        >>> vor = mf.VoronoiGridPlus.from_gsf("flow-tt01.gsf", crs="EPSG:2927")
        >>> vor.plot.map().show()
        """
        verts, iverts, xcyc = get_griddata_from_gsf(gsf_path, layer=layer)
        return cls(
            verts=verts,
            iverts=iverts,
            xcyc=xcyc,
            crs=crs,
            name=name,
            rasters=rasters,
            idomain_path=idomain_path,
        )


def _top_per_cell(top, ncpl: int) -> np.ndarray:
    """Normalize top elevation data to a one-dimensional ``ncpl``-length array."""

    values = np.asarray(top, dtype=float).reshape(-1)
    if values.size == 1:
        values = np.repeat(values, ncpl)
    if values.size != ncpl:
        raise ValueError(f"Top array size {values.size} does not match ncpl {ncpl}.")
    return values


def _bottoms_per_layer(botm, ncpl: int) -> np.ndarray:
    """Normalize bottom elevation data to a ``(nlay, ncpl)`` array."""

    values = np.asarray(botm, dtype=float)
    if values.ndim == 0:
        return np.repeat(values.reshape(1), ncpl).reshape(1, ncpl)
    if values.ndim == 1:
        if values.size == ncpl:
            return values.reshape(1, ncpl)
        if values.size % ncpl == 0:
            return values.reshape(values.size // ncpl, ncpl)
    if values.ndim >= 2:
        trailing = int(np.prod(values.shape[1:]))
        if trailing == ncpl:
            return values.reshape(values.shape[0], ncpl)
    raise ValueError(f"Bottom array shape {values.shape} cannot be reshaped to ncpl {ncpl}.")


def grid_view_from_modelgrid(modelgrid, *, crs, name: str) -> VoronoiGridPlus:
    """A :class:`VoronoiGridPlus` view of a FloPy ``modelgrid``, with its layer elevations.

    The view :func:`myflopy.load_mf6_run` gives any MF6 run -- DIS, DISV or DISU --
    and :func:`model_grid_view` gives a live model that has none: one polygon per
    cell of a layer from ``modelgrid.verts``/``iverts`` in FloPy's order (on DIS,
    row ``r`` column ``c`` is cell ``r * ncol + c``), the cell centres from
    ``xcellcenters``/``ycellcenters``, and ``gdf_topbtm`` from ``top``/``botm``.

    Parameters
    ----------
    modelgrid : flopy.discretization.Grid
        A FloPy model grid: ``StructuredGrid``, ``VertexGrid`` or
        ``UnstructuredGrid``. Its ``verts``, ``iverts``, ``xcellcenters``,
        ``ycellcenters``, ``top`` and ``botm`` are read.
    crs : str, int, pyproj.CRS or None
        The view's CRS, as :class:`VoronoiGridPlus` takes it. Nothing is
        reprojected: the coordinates are the grid's own.
    name : str
        A label stored as the view's ``name``.

    Returns
    -------
    VoronoiGridPlus
        The view, its ``gdf_topbtm`` set (geometry first, then the top, then one
        bottom per layer).

    Raises
    ------
    ValueError
        When ``top`` or ``botm`` cannot be read as one value per cell of a layer.
    """

    verts = np.asarray(modelgrid.verts, dtype=float)
    iverts = [[int(vertex) for vertex in cell] for cell in modelgrid.iverts]
    xcyc = np.column_stack(
        [
            np.asarray(modelgrid.xcellcenters, dtype=float).reshape(-1),
            np.asarray(modelgrid.ycellcenters, dtype=float).reshape(-1),
        ]
    )
    view = VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc, crs=crs, name=name)
    ncpl = view.ncpl
    data: dict[object, object] = {
        "geometry": view.gdf_vorPolys.geometry,
        0: _top_per_cell(modelgrid.top, ncpl),
    }
    for idx, layer_botm in enumerate(_bottoms_per_layer(modelgrid.botm, ncpl), start=1):
        data[idx] = layer_botm
    view.gdf_topbtm = gpd.GeoDataFrame(data, geometry="geometry", crs=view.crs)
    return view


def model_grid_view(model) -> VoronoiGridPlus:
    """A model's grid view: ``model.vor``, else one built from FloPy's ``modelgrid``.

    A model declared with ``mf.dis(...)`` carries no grid view (``model.vor`` is
    ``None``), while the same run opened with :func:`myflopy.load_mf6_run` has one.
    A view that only draws a model's cells -- the PEST review's maps -- reads it
    here, so both draw: the live model gets the view :func:`grid_view_from_modelgrid`
    builds (the one ``load_mf6_run`` builds), in FloPy's ``modelgrid.crs``, else
    :data:`~myflopy._defaults.DEFAULT_CRS` (``load_mf6_run``'s default). Its
    cells are the cells :func:`~myflopy.modflow.mf6.grid.geometry.model_cell_polygons`
    returns, in the same order (ledger 451).

    Parameters
    ----------
    model : SimulationBase or LoadedMf6Run
        The model. Read: its ``vor`` (may be ``None``), else the ``modelgrid`` of
        its ``gwf``, else its own ``modelgrid``.

    Returns
    -------
    VoronoiGridPlus
        ``model.vor`` when set; otherwise the built view, cached on the model for
        as long as its FloPy grid object is the same one.

    Raises
    ------
    ValueError
        When the model has neither a grid view nor a FloPy ``modelgrid``.
    """

    vor = getattr(model, "vor", None)
    if vor is not None:
        return vor
    gwf = getattr(model, "gwf", None)
    grid = getattr(gwf, "modelgrid", None) if gwf is not None else None
    if grid is None:
        grid = getattr(model, "modelgrid", None)
    if grid is None:
        raise ValueError(
            f"Model {getattr(model, 'name', model)!r} has no grid to draw: neither a "
            "grid view (`model.vor`) nor a FloPy `modelgrid`."
        )
    cached = getattr(model, "_grid_view_from_modelgrid", None)
    if cached is not None and cached[0] is grid:
        return cached[1]
    view = grid_view_from_modelgrid(
        grid, crs=getattr(grid, "crs", None) or DEFAULT_CRS, name=str(getattr(model, "name", "model"))
    )
    try:
        model._grid_view_from_modelgrid = (grid, view)
    except AttributeError:
        logger.debug("Grid view of %r not cached: the model takes no attribute", model)
    return view

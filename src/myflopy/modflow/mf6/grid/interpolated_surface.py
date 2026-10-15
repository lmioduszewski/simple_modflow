"""`InterpolatedSurface` — model-aware interpolated surfaces + plotly traces.

Moved from `modflow/utils/surfaces.py` (implementation plan 4.3); that path
remains as a warned facade.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus as Vor
    from myflopy.modflow.mf6.simulation.base import SimulationBase

from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objs as go
import rasterio
import shapely as shp
from pandas import IndexSlice as idxx
from rasterio.features import geometry_mask
from rasterio.io import MemoryFile
from rasterio.mask import mask
from rasterio.transform import from_origin
from scipy.interpolate import RBFInterpolator, griddata
from scipy.spatial import QhullError
from shapely.geometry import MultiPolygon, Polygon, mapping

from myflopy import viz as f
from myflopy._logging import get_logger
from myflopy.modflow.mf6.grid.geometry import (
    miss_reason,
    selected_cells,
    where_they_are,
)
from myflopy.modflow.mf6.headsplus import DependentVariableFile, _is_sentinel
from myflopy.modflow.mf6.headsplus import HeadsPlus as Hp
from myflopy.modflow.mf6.package_explorer_utils import check_view_layers
from myflopy.modflow.mf6.period_end import last_saved_step
from myflopy.modflow.utils.datatypes.readers import read_shp_gpkg
from myflopy.viz import Picture

logger = get_logger(__name__)


#: What ``InterpolatedSurface(surf_type=...)`` interpolates: the head field, or a
#: layer's top elevation.
SURFACE_TYPES = ("hds", "lyr")


class InterpolatedSurface(Picture):

    #: Whether `.fig` has been assembled. A CLASS attribute, matching the other
    #: pictures, so instances built via `object.__new__` still answer.
    _assembled = False
    _fig = None


    def __init__(
            self,
            xs: np.array = None,
            ys: np.array = None,
            zs: np.array = None,
            vor: Vor = None,
            hds: Hp = None,
            model: SimulationBase = None,
            layer: int = 0,
            per: int = None,
            kstpkper: tuple = None,
            resolution: int = 1000,
            use_rbf: bool = False,
            surf_type: str = 'hds',
            clip: shp.Polygon | Path = None,
            interpolator: str = None,
            crs = None
    ):
        """
        Base class for interpolated surfaces.
        :param xs: optional, array of x coordinates. If not provided will get from vor
        :param ys: optional, array of y coordinates. If not provided will get from vor
        :param zs: optional, array of z coordinates. If not provided will get from hds obj
        :param vor: voronoi grid object corresponding to model, will get from model if not provided
        :param hds: optional dependent-variable reader of the model -- ``HeadsPlus``,
            ``ConcResults`` or ``TempResults`` (any ``DependentVariableFile``);
            without one, the model's own reader (``model.hds``/``.conc``/``.temp``
            by its kind) is used
        :param model: mf6 SimulationBase model
        :param layer: zero-based layer, defaults to 0. With ``surf_type='hds'``
            and no ``zs``, a layer the model's views do not number -- ONE on a
            DISU -- raises ``ValueError`` naming the count (ledger 496); it
            raised pandas' bare ``KeyError`` when the surface was drawn.
        :param per: stress period, zero-based: its LAST saved step is drawn
            (ledger 337, the one rule -- ``period_end.last_saved_step``); it
            indexed the saved records before, so ``per=1`` at ``nstp=3`` drew
            ``(1, 0)``. Give ``per`` or ``kstpkper``, not both.
        :param kstpkper: tuple of time step and period for surface
        :param resolution: defaults to 1000
        :param use_rbf: defaults to False
        :param surf_type: ``'hds'`` (default) interpolates the model's own field
            -- heads, or concentration/temperature on a GWT/GWE model -- at
            ``kstpkper``; ``'lyr'`` the TOP of layer ``layer`` -- the model top for
            ``layer=0``, the bottom of layer ``k - 1`` for ``layer=k`` -- read from
            the model's own discretization when a model is given, else from the
            grid's ``gdf_topbtm``. Anything else raises.
        :param clip: optional region -- a (Multi)Polygon or a vector file path -- to
            clip the surface to. Without one, a surface that knows its grid is
            masked to the grid's footprint (see ``region``).
        """
        self.model = model
        try:
            self.vor = self.model.vor if vor is None else vor
        except Exception:  # noqa: BLE001 - see below; the tuple cannot be closed
            # Deliberately broad. `model.vor` is a plain attribute on a live
            # SimulationBase, but a LAZY PROPERTY on the file-backed
            # LoadedMf6Run that loads the flopy simulation and rebuilds the
            # Voronoi grid -- so it reaches flopy (MFDataException),
            # geopandas/pyproj (CRSError, a RuntimeError), shapely (ValueError),
            # and the filesystem, and the caller-supplied `crs=` means a bad
            # projection string is a user input, not a bug. A surface plotted
            # without a CRS is the documented fallback; refusing to build the
            # object is not.
            logger.debug("could not resolve a grid from the model", exc_info=True)
            self.vor = None
        self._interpolators = ['griddata', 'rbf', 'linearND', 'cloughTocher2D']
        self._interpolator = None
        self.interpolator = interpolator
        self.crs = crs
        if surf_type not in SURFACE_TYPES:
            # Anything else left `zs` unset, and the refusal said every cell was
            # NaN rather than naming the word.
            raise ValueError(
                f"surf_type must be one of {SURFACE_TYPES}, not {surf_type!r}."
            )
        self.surf_type = surf_type
        self._xs = xs
        self._ys = ys
        self._zs = zs
        self._xys = None
        self.use_rbf = use_rbf
        self._griddata_interp = None
        self._rbf_interp = None
        self._linearND_interp = None
        self._cloughTocher2D_interp = None
        self._meshgrid = None
        self._kstpkper = None
        self._hds = hds
        if per is not None and kstpkper is not None:
            raise ValueError("per= and kstpkper= both choose the time; pass one.")
        if per is not None:
            self.kstpkper = last_saved_step(self.model.kstpkper, per)
        else:
            self.kstpkper = kstpkper
        if surf_type == 'hds' and zs is None:
            # A layer the model's views do not number is refused here, naming the
            # count -- ONE on a DISU -- rather than at `.fig` with pandas' bare
            # `KeyError` (`check_view_layers`, ledger 496).
            gwf = getattr(model, "gwf", None) if model is not None else None
            if gwf is not None:
                check_view_layers(gwf, layer)
            elif hds is not None:
                check_view_layers(int(hds.nlay), layer)
        self.layer = layer
        self.resolution = resolution
        self.neighbors = 10
        self.colorscale = 'Earth_r'  # reverse earth so blue is low nums
        self._clip = None
        self._clipped_cells = None
        self._domain = None

        self.clip = clip

    @property
    def clip(self):
        """clip is a shapely (Multi)Polygon to clip the interpolated surface to"""
        return self._clip

    @clip.setter
    def clip(self, val):
        """Set the clip region from ``None``, a vector file path, or a (Multi)Polygon.

        Anything else RAISES. A ``str`` path and a ``MultiPolygon`` used to fall
        through a setter with no branch for them, so the surface came back
        unclipped with nothing said (ledger 314).
        """

        if val is None:
            self._clip = None
        elif isinstance(val, (str, Path)):
            clip = read_shp_gpkg(Path(val))
            assert len(clip) == 1, 'clip must be a single polygon'
            clip = clip.union_all()
            assert isinstance(clip, (Polygon, MultiPolygon)), 'clip must be a shapely Polygon'
            self._clip = clip
        elif isinstance(val, (Polygon, MultiPolygon)):
            self._clip = val
        else:
            raise TypeError(
                "clip must be a Polygon, a MultiPolygon or the path of a vector "
                f"file, not {type(val).__name__}."
            )


    @property
    def clipped_cells(self) -> list[int] | None:
        """The cells the ``clip=`` region covers, or ``None`` when there is no clip.

        The selector answers one entry per FEATURE, each a LIST of cells, so
        ``.to_list()`` handed ``.loc`` the label ``(5,)`` and EVERY clipped
        surface raised ``KeyError``; :func:`selected_cells` flattens it.

        Resolved here, once, rather than inside ``xs`` and ``ys`` so that ``zs``
        reads the same cells whichever is asked for first -- reading ``zs`` on its
        own took every cell in the grid while ``xs`` took the clipped ones, and the
        two arrays then disagreed in length. Memoised, so the selector runs once
        per surface rather than once per ``xs``/``ys``/``zs``.

        A clip region that covers no cell RAISES: an empty point set reaches the
        interpolator, which fails for a reason nothing connects back to the clip.
        """

        if not self.clip:
            return None
        if self._clipped_cells is None:
            cells = selected_cells(self.vor.get_vor_cells_as_series(self.clip))
            if not cells:
                polys = self.vor.gdf_vorPolys
                raise ValueError(
                    "The clip region covers no grid cell "
                    f"({miss_reason(polys, self.clip, target='grid cell')}), so there "
                    "is nothing to interpolate. "
                    + where_they_are(self.clip, polys, noun="clip bounds")
                )
            self._clipped_cells = cells
        return self._clipped_cells

    @property
    def region(self):
        """What a drawn surface is masked to: ``clip=``, else the grid's footprint.

        ``None`` for a bare point cloud, which has no footprint to mask to.

        The interpolator is NOT a mask. ``griddata`` fills the convex hull of the
        points it is given and is NaN outside it, so a surface left to it LOOKS
        clipped while spanning every concavity of a domain that is not convex --
        and ``.show()`` drew exactly that, while only :meth:`plot_heatmap` masked
        (ledger 314). The footprint is the union of every cell, active or not,
        and is computed once.
        """

        if self.clip:
            return self.clip
        if self.vor is None:
            return None
        if self._domain is None:
            self._domain = self.vor.gdf_vorPolys.union_all()
        return self._domain

    def _required_region(self):
        """:attr:`region`, or a ``ValueError`` saying there is nothing to clip to."""

        region = self.region
        if region is None:
            raise ValueError(
                "There is nothing to clip to: no clip= was given, and this surface "
                "has no grid to take a footprint from. Pass clip=, or vor=."
            )
        return region

    @property
    def interpolator(self):
        """The interpolation method used to build the surface raster."""

        return self._interpolator

    @interpolator.setter
    def interpolator(self, val):
        """Set the interpolation method (ignored with a warning if unsupported)."""

        if val is not None:
            if val not in self._interpolators:
                logger.warning(
                    'interpolator %r is not one of %s; ignoring it',
                    val, self._interpolators,
                )
            else:
                self._interpolator = val

    @property
    def kstpkper(self):
        """The ``(timestep, period)`` sampled for a head surface (defaults to the model's first)."""

        if self._kstpkper is None:
            if self.model:
                self._kstpkper = self.model.kstpkper[0]
            else:
                raise ValueError(
                    "no model was given, so there is no default (timestep, period) "
                    "to sample; pass kstpkper=(kstp, kper) or model=."
                )
        return self._kstpkper

    @kstpkper.setter
    def kstpkper(self, val):
        """Set the sampled ``(timestep, period)`` tuple (rejects a non-length-2 value)."""

        if val is not None:
            if len(val) == 2:
                assert isinstance(val, tuple), (
                    'kstpkper must be a tuple of (timestep, period)'
                )
                self._kstpkper = val
            else:
                logger.warning(
                    'kstpkper must be a (timestep, period) pair, got %r; '
                    'keeping the previous value', val,
                )

    @property
    def xs(self):
        """Cell-centroid x-coordinates (clipped to the clip region when set; cached)."""

        if self._xs is None:
            if self.clipped_cells is not None:
                xs = pd.Series(self.vor.centroids_x).loc[self.clipped_cells].to_numpy()
            else:
                xs = np.array(self.vor.centroids_x)
            self._xs = xs
        return self._xs

    @property
    def ys(self):
        """Cell-centroid y-coordinates (clipped to the clip region when set; cached)."""

        if self._ys is None:
            if self.clipped_cells is not None:
                ys = pd.Series(self.vor.centroids_y).loc[self.clipped_cells].to_numpy()
            else:
                ys = np.array(self.vor.centroids_y)
            self._ys = ys
        return self._ys

    @property
    def zs(self):
        """Per-cell surface values to interpolate: heads (``hds``) or a layer elevation (``lyr``).

        Head surfaces read the selected ``kstpkper``/``layer``, with MODFLOW's
        no-value sentinels -- +1e30 inactive, -1e30 DRY -- as NaN and every real
        head kept, whatever its size; layer surfaces read ``gdf_topbtm``. Cached.

        The sentinels were cut at ``> 10_000``, which kept a dry cell's -1e30 as a
        sample and erased real heads above 10,000 ft. On a 207-cell run with 25
        dry cells, every finite ``griddata`` pixel fell outside the real 30-50 ft
        range; with heads 9,950-10,050 ft, 115 cells were erased and the surface
        stopped at 9,998 (ledger 315).

        Values the caller supplied are returned as given, without asking which
        cells ``clip=`` covers: that question needs a grid, and ran the selector
        -- and could refuse the clip -- for points that never came from one.
        """

        if self._zs is None:
            cells = slice(None) if self.clipped_cells is None else self.clipped_cells
            if self.surf_type == 'hds':
                """zs of the self.hds HeadPlus oject for a given layer at a certain stress-and-time period"""
                zs = np.asarray(
                    self.hds.all_values.loc[idxx[self.kstpkper, self.layer, cells], :].values,
                    dtype=float,
                )
                zs[_is_sentinel(zs)] = np.nan  # inactive (+1e30) and dry (-1e30) cells
                self._zs = zs
            if self.surf_type == 'lyr':
                self._zs = self._layer_top_table().loc[cells, self.layer].values
        return self._zs

    def _layer_top_table(self) -> pd.DataFrame:
        """Per-cell layer elevations, columns ``0`` (model top) to ``nlay`` (last bottom).

        A model's OWN discretization -- what MF6 ran, and what the section's model
        top and layer bottoms read -- when there is one. The grid's ``gdf_topbtm``
        was the only source, and a grid built for a hand-written DISV has none, so
        ``surf_type='lyr'`` died with ``'NoneType' object has no attribute 'loc'``
        on such a model while working on the same run reloaded (ledger 337). A
        bare grid still reads its own table, and without either this says so.
        """

        flopy_model = getattr(self.model, "gwf", None)
        modelgrid = getattr(flopy_model, "modelgrid", None)
        if modelgrid is not None and modelgrid.top is not None:
            top = np.asarray(modelgrid.top, dtype=float).ravel()
            botm = np.asarray(modelgrid.botm, dtype=float).reshape(-1, top.size)
            return pd.DataFrame({0: top, **{k + 1: row for k, row in enumerate(botm)}})
        table = getattr(self.vor, "gdf_topbtm", None) if self.vor is not None else None
        if table is None:
            raise ValueError(
                "surf_type='lyr' needs layer elevations: pass model= (its "
                "discretization is read), or a grid whose gdf_topbtm is set."
            )
        return table

    @zs.setter
    def zs(self, val):
        """Override the per-cell surface values with a caller-supplied array."""

        self._zs = val

    @property
    def xys(self):
        """The cell-centroid points as a list of ``(x, y)`` tuples (cached)."""

        if self._xys is None:
            xys = list(zip(self.xs, self.ys, strict=False))
            self._xys = xys
        return self._xys

    @property
    def xy_meshgrid(self):
        """A ``(grid_x, grid_y)`` regular meshgrid over the data extent at ``resolution`` (cached)."""

        if self._meshgrid is None:
            xs = self.xs
            ys = self.ys
            grid_x, grid_y = np.meshgrid(
                np.linspace(xs.min(), xs.max(), self.resolution),
                np.linspace(ys.max(), ys.min(), self.resolution)
            )
            self._meshgrid = (grid_x, grid_y)
        return self._meshgrid

    @property
    def hds(self):
        """The model's dependent-variable reader -- heads, concentration or temperature (cached).

        A model's OWN reader, ``model.hds``/``.conc``/``.temp`` by its kind (the
        file its OC names). It opened ``<model name>.hds`` whatever the kind, so
        the interpolated profile of a GWT model raised ``FileNotFoundError`` --
        the model writes ``<name>.ucn`` (ledger 436). A duck-typed model with no
        ``_field_reader`` reads the head file its OC declares, by the same rule
        (ledger 598).
        """

        if self._hds is None:
            if self.model:
                reader = getattr(self.model, "_field_reader", None)
                self._hds = reader if reader is not None else Hp(
                    hds_path=Hp._declared_path(self.model),
                    vor=self.vor
                )
            else:
                raise ValueError(
                    "no heads file was given and no model was given to find one on; "
                    "pass hds= or model=."
                )
        assert isinstance(self._hds, DependentVariableFile), (
            "the hds property must be a dependent-variable reader (HeadsPlus, "
            "ConcResults or TempResults)"
        )
        return self._hds

    @property
    def _finite_samples(self):
        """``(points, values)`` with the no-value cells dropped.

        You cannot interpolate FROM a cell that has no value, and scipy does not
        drop them for you: `griddata` triangulates every point it is given, so a
        single NaN vertex poisons each output cell of every triangle touching it.
        Scattered NaNs -- a DEM masked to a boundary, inactive cells in a head
        array -- therefore produce an ALL-NaN surface that renders as nothing at
        all, which reads as "plotting is broken" rather than "some cells have no
        data". Measured on a real stack: 300 valid cells of 582 gave 0 finite
        pixels of 14,400.

        Every interpolator must take its samples from here. ``linearND`` and
        ``cloughTocher2D`` were missed when this was added (ledger 152) and kept
        reading the raw arrays: 300 points with 40% NaN left them 292 and 0
        finite pixels of 1,600, where ``griddata`` kept 1,357.
        """

        xs = np.asarray(self.xs, dtype=float)
        ys = np.asarray(self.ys, dtype=float)
        zs = np.asarray(self.zs, dtype=float).ravel()
        keep = np.isfinite(xs) & np.isfinite(ys) & np.isfinite(zs)
        if not keep.any():
            raise ValueError(
                "Nothing to interpolate: every cell in this surface is NaN. A "
                "layer with no coverage, or a raster whose no-data value swallowed "
                "the whole grid -- check `mf.Raster(..., nodata=)` and the extent."
            )
        return np.column_stack((xs[keep], ys[keep])), zs[keep]

    @property
    def griddata_interp(self):
        """interpolated surface using scipy griddata"""
        if self._griddata_interp is None:
            points, values = self._finite_samples
            zis = griddata(
                points=points,
                values=values,
                xi=self.xy_meshgrid,
                method='cubic')
            self._griddata_interp = zis.squeeze()
        return self._griddata_interp

    @property
    def rbf_interp(self):
        """interpolated surface using scipy RBFInterpolator"""
        if self._rbf_interp is None:
            coords, values = self._finite_samples
            xis = self.xy_meshgrid[0].ravel()
            yis = self.xy_meshgrid[1].ravel()
            xyis = np.column_stack((xis, yis))
            interpolator = RBFInterpolator(
                coords,
                values,
                neighbors=self.neighbors
            )
            grid_z = interpolator(xyis).reshape(self.xy_meshgrid[0].shape)
            self._rbf_interp = grid_z
        return self._rbf_interp

    @property
    def linearND_interp(self):
        """interpolated surface using scipy LinearNDInterpolator"""
        if self._linearND_interp is None:
            from scipy.interpolate import LinearNDInterpolator
            coords, values = self._finite_samples
            interpolator = LinearNDInterpolator(coords, values)
            grid_z = interpolator(self.xy_meshgrid).reshape(self.xy_meshgrid[0].shape)
            self._linearND_interp = grid_z
        return self._linearND_interp

    @property
    def cloughTocher2D_interp(self):
        """interpolated surface using scipy CloughTocher2DInterpolator"""
        if self._cloughTocher2D_interp is None:
            from scipy.interpolate import CloughTocher2DInterpolator
            coords, values = self._finite_samples
            interpolator = CloughTocher2DInterpolator(coords, values)
            grid_z = interpolator(self.xy_meshgrid).reshape(self.xy_meshgrid[0].shape)
            self._cloughTocher2D_interp = grid_z
        return self._cloughTocher2D_interp

    @property
    def transform(self):
        """The raster's affine transform, registered to :attr:`xy_meshgrid`.

        ``xy_meshgrid`` puts its first and last pixel CENTRES on the edges of the
        points' extent, and ``from_origin`` takes a pixel's CORNER, so the origin
        is half a pixel up and to the left of the first point. It sat ON the first
        point, which put every pixel half a pixel east and south of where the
        surface draws it -- where a clip was tested, a saved raster was
        georeferenced and a cross-section sampled (ledger 314).
        """
        xs = self.xs
        ys = self.ys

        xmin, xmax = np.min(xs), np.max(xs)
        ymin, ymax = np.min(ys), np.max(ys)
        xres = (xmax - xmin) / (self.resolution - 1)
        yres = (ymax - ymin) / (self.resolution - 1)

        return from_origin(xmin - xres / 2, ymax + yres / 2, xres, yres)

    @property
    def memfile(self):
        """An in-memory single-band GeoTIFF (:class:`MemoryFile`) of the interpolated surface.

        Unmasked -- NaN only outside the convex hull of the points -- and it
        declares NaN as its no-data value, so a clip of it fills with NaN.
        """

        grid_z = np.asarray(self.surface, dtype=float)
        try:
            crs = self.vor.crs
        except AttributeError:
            # `vor` is None (see __init__) or is a grid object that carries no
            # crs; fall back to the one the caller passed.
            crs = self.crs

        memfile = MemoryFile()
        with memfile.open(
                driver='GTiff',
                height=self.resolution,
                width=self.resolution,
                count=1,
                dtype=grid_z.dtype,
                crs=crs,
                transform=self.transform,
                nodata=np.nan) as dst:
            dst.write(grid_z, 1)

        return memfile

    @property
    def projected_xys(self):
        """Gets projected x and y coordinates for the surface. Used in plotting."""
        with self.memfile.open() as dataset:

            data = dataset.read(1)
            transform = dataset.transform
            nrows, ncols = data.shape
            xs, ys = np.meshgrid(np.arange(ncols), np.arange(nrows))
            projected_x, projected_y = rasterio.transform.xy(transform, ys, xs, offset='center')
            projected_x = np.array(projected_x)
            projected_y = np.array(projected_y)

            return projected_x, projected_y

    @property
    def surface(self):
        """returns griddata interpolation first, if that fails then return the rbf interpolation"""
        if self.interpolator is not None:
            if self.interpolator == 'griddata':
                return self.griddata_interp
            elif self.interpolator == 'rbf':
                return self.rbf_interp
            elif self.interpolator == 'linearND':
                return self.linearND_interp
            elif self.interpolator == 'cloughTocher2D':
                return self.cloughTocher2D_interp
            else:
                # `self.interpolator`, not `val`: this is the `interp` property,
                # not the setter that has a `val` parameter. The name survived a
                # print -> logger conversion (017a69a) and would raise NameError
                # if a caller set `_interpolator` directly, past the setter's
                # own check.
                logger.warning(
                    'interpolator %r is not one of %s; ignoring it',
                    self.interpolator, self._interpolators,
                )
        if self.use_rbf is False:
            try:
                return self.griddata_interp
            except (QhullError, ValueError):
                # QhullError (a RuntimeError) is what scipy raises for a
                # degenerate point set -- collinear points, or too few of them
                # to triangulate. ValueError covers NaN/inf in the inputs. Both
                # mean "griddata cannot triangulate this", which is exactly when
                # the radial-basis interpolator is the right second choice.
                logger.debug(
                    "griddata interpolation failed; falling back to rbf",
                    exc_info=True,
                )
                return self.rbf_interp
        elif self.use_rbf is True:
            return self.rbf_interp

    @property
    def clipped_surface(self):
        """:attr:`surface` with every pixel outside :attr:`region` set to NaN -- what is drawn.

        The same values as :attr:`surface` when there is no region.
        """

        region = self.region
        if region is None:
            return np.array(self.surface, dtype=float)
        return self._masked(self.surface, region)

    def _inside(self, region) -> np.ndarray:
        """``(resolution, resolution)`` booleans: True where a pixel's CENTRE is in ``region``.

        The test ``rasterio.mask.mask`` applies to :attr:`memfile` --
        ``geometry_mask`` on :attr:`transform`, ``all_touched=False`` -- so the
        drawn surface, the heatmap and a saved raster keep the same pixels.
        """

        return geometry_mask(
            [mapping(region)],
            out_shape=(self.resolution, self.resolution),
            transform=self.transform,
            invert=True,
        )

    def _masked(self, surface, region):
        """A float COPY of ``surface`` with every pixel outside ``region`` set to NaN.

        A copy, because :attr:`surface` is the cached interpolation and masking it
        in place would mask every later read of it. Masking that leaves nothing to
        draw RAISES: the points and the region do not meet, which is a location or
        a CRS, and an empty picture says neither.
        """

        surface = np.array(surface, dtype=float)
        shape = (self.resolution, self.resolution)
        if surface.shape != shape:
            raise ValueError(
                f"surface= must be on this surface's own pixel grid, shape {shape}; "
                f"got {surface.shape}."
            )
        had_values = bool(np.isfinite(surface).any())
        surface[~self._inside(region)] = np.nan
        if had_values and not np.isfinite(surface).any():
            points = (np.nanmin(self.xs), np.nanmin(self.ys), np.nanmax(self.xs), np.nanmax(self.ys))
            span = ", ".join(f"{v:,.1f}" for v in points)
            where = ", ".join(f"{v:,.1f}" for v in region.bounds)
            raise ValueError(
                "Masking this surface to its region leaves nothing to draw: no "
                f"interpolated pixel falls inside it. The points span ({span}) and "
                f"the region ({where}); when those do not overlap it is almost "
                "always a location or a CRS."
            )
        return surface

    def clip_raster_with_polygon(self, polygon: Polygon = None, buffer=None):
        """Clip the in-memory raster with a polygon and return the clipped raster data.

        ``polygon`` defaults to :attr:`region` -- ``clip=``, else the grid's
        footprint. Pixels outside it come back NaN. They came back 0, a value a
        head or an elevation can really take, and callers then blanked every real
        0 along with them (ledger 314). The raster is not cropped: the result
        keeps the full ``resolution`` x ``resolution`` pixel grid. A pixel is
        kept when its CENTRE is inside the polygon.

        Parameters
        ----------
        polygon : shapely.geometry.Polygon or shapely.geometry.MultiPolygon, optional
            The polygon to clip to, in the surface's own x/y coordinates (the
            grid's CRS and length units). The default ``None`` uses
            :attr:`region`, and raises ``ValueError`` when there is none (a bare
            point cloud given no ``clip=``). A file path is NOT read here, unlike
            ``clip=``. It fails with ``AttributeError``.
        buffer : float, optional
            Distance, in the surface's length units, to grow the polygon by
            (negative shrinks it) with shapely's ``buffer`` before clipping. It
            applies to the default region as well. The default ``None``, like
            ``0``, leaves the polygon as it is.

        Returns
        -------
        tuple of (numpy.ndarray, affine.Affine, dict)
            The clipped image, shape ``(1, resolution, resolution)`` with NaN
            outside the polygon; its transform; and the in-memory raster's
            metadata (``nodata`` NaN, the grid's CRS or ``crs=``).
        """
        polygon = self._required_region() if polygon is None else polygon
        polygon = polygon.buffer(buffer) if buffer else polygon
        with self.memfile.open() as dataset:
            shapes = [mapping(polygon)]

            # Clip the raster with the polygon
            clipped_image, clipped_transform = mask(dataset, shapes, crop=False, nodata=np.nan)

            return clipped_image, clipped_transform, dataset.meta

    def save_raster(
            self,
            output_tif='clipped_raster.tif',
            clipped_image=None,
            clipped_transform=None,
            meta=None,
            polygon_clip: Polygon = None,
            buffer=0
    ):
        """Save clipped raster data to a GeoTIFF file.

        Pass all three of ``clipped_image``, ``clipped_transform`` and ``meta``
        (what :meth:`clip_raster_with_polygon` returns), or none of them to clip
        here. The written file declares NaN as no-data, and a real 0 is written
        as 0 -- it was written as NaN, with the clipped-away pixels.

        Parameters
        ----------
        output_tif : str or os.PathLike, default "clipped_raster.tif"
            The GeoTIFF to write. An existing file is overwritten, and its
            directory must exist. The default is a RELATIVE name, so it lands
            in the current working directory.
        clipped_image : numpy.ndarray, optional
            The pixels to write, shaped ``(bands, rows, cols)`` as
            :meth:`clip_raster_with_polygon` returns them. Written as given.
        clipped_transform : affine.Affine, optional
            The affine transform of ``clipped_image``.
        meta : dict, optional
            Raster metadata (rasterio's ``dataset.meta``). It is copied, never
            edited. ``driver``, ``height``, ``width`` and ``transform`` are set
            from the image and transform, and every other key (``dtype``,
            ``count``, ``crs``, ``nodata``) is written as given.
        polygon_clip : shapely.geometry.Polygon or shapely.geometry.MultiPolygon, optional
            What to clip to when the three pieces above are not given, in the
            surface's own coordinates. The default ``None`` uses ``clip=``,
            else the grid's footprint, and raises ``ValueError`` when the
            surface has neither. Given WITH the three pieces it raises
            ``ValueError``: nothing is clipped then (it was ignored in silence
            until ledger 481).
        buffer : float, default 0
            Distance, in the surface's length units, to grow (negative: shrink)
            the clip polygon by. ``0`` leaves it as it is. Anything else with the
            three pieces raises ``ValueError``, as ``polygon_clip`` does.

        Returns
        -------
        None

        Raises
        ------
        ValueError
            When some but not all of ``clipped_image``, ``clipped_transform``
            and ``meta`` are given, or ``polygon_clip``/``buffer`` with all three.
        """

        pieces = (clipped_image, clipped_transform, meta)
        if all(piece is None for piece in pieces):
            polygon_clip = self.clip if polygon_clip is None else polygon_clip
            clipped_image, clipped_transform, meta = self.clip_raster_with_polygon(polygon=polygon_clip, buffer=buffer)
        elif any(piece is None for piece in pieces):
            raise ValueError(
                "Pass all three of clipped_image, clipped_transform and meta, or "
                "none of them."
            )
        elif polygon_clip is not None or buffer != 0:
            # Nothing is clipped when the pieces are given, so these could not act
            # (ledger 481: they were ignored in silence).
            raise ValueError(
                "polygon_clip=/buffer= say how to CLIP, and the clipped image given is "
                "written as it is; clip it with them (clip_raster_with_polygon), or "
                "drop the image and let save_raster clip."
            )

        # Update the metadata with the new transform and dimensions -- a copy, so
        # the caller's own `meta` is left as it was.
        meta = {
            **meta,
            "driver": "GTiff",
            "height": clipped_image.shape[1],
            "width": clipped_image.shape[2],
            "transform": clipped_transform
        }

        # Write the clipped raster to a GeoTIFF file
        with rasterio.open(output_tif, "w", **meta) as dst:
            dst.write(clipped_image)

    def plot_heatmap(self):
        """Show the surface as a Plotly heatmap, north up, masked like :attr:`fig`.

        Drawn on the surface's own x/y at one unit to one unit. Row 0 of the
        raster is its NORTH edge, and a heatmap given no ``y`` draws row 0 at the
        BOTTOM, so this was mirrored north-south; the clipped-away pixels were
        drawn as the value 0, which also dragged the colour scale down to 0
        (ledger 314).
        """

        grid_x, grid_y = self.xy_meshgrid
        z = self.clipped_surface
        fig = f.Fig()
        # The raster's rows run north to south; hand plotly ascending y, and the
        # rows reversed to match, so north is up.
        fig.add_heatmap(x=grid_x[0], y=grid_y[::-1, 0], z=z[::-1])
        fig.update_yaxes(scaleanchor="x", scaleratio=1)
        fig.show()


    def surface_trace(self, *, surface=None, colorscale=None, **kwargs):
        """Return a Plotly ``go.Surface`` trace for this interpolated surface.

        The single builder shared by :attr:`fig` and higher-level multi-surface
        views (e.g. ``myflopy.layers.LayerBuildResult._surface_fig``). ``surface``
        overrides the z-array, which defaults to :attr:`clipped_surface` -- what
        :attr:`fig` draws; any extra keyword arguments pass straight through to
        ``go.Surface``.

        Parameters
        ----------
        surface : numpy.ndarray, optional
            The z values to draw, one per pixel of this surface's grid, shape
            ``(resolution, resolution)``, row 0 at the NORTH edge, in the
            surface's own units (a head or an elevation, model length units).
            NaN pixels are gaps. It is not checked against the grid here. The
            default ``None`` is :attr:`clipped_surface`, masked to
            :attr:`region` (unmasked when there is none, a bare point cloud).
        colorscale : str or list, optional
            A Plotly colorscale name or ``[[position, color], ...]`` stops. The
            default ``None`` uses this surface's ``colorscale`` attribute
            (``"Earth_r"`` as constructed).
        **kwargs
            Passed on to ``plotly.graph_objects.Surface`` (``opacity``,
            ``showscale``, ``name``, ``cmin``/``cmax``, ...). ``x``, ``y`` and
            ``z`` are set here, so passing one raises ``TypeError``.

        Returns
        -------
        plotly.graph_objects.Surface
            The trace, not yet added to any figure.
        """
        surface = self.clipped_surface if surface is None else surface
        xg, yg = self.xy_meshgrid
        return go.Surface(
            x=xg,
            y=yg,
            z=surface,
            colorscale=self.colorscale if colorscale is None else colorscale,
            **kwargs,
        )

    @property
    def fig(self):
        """The 3-D surface figure, masked to :attr:`region`.

        That is ``clip=`` when given, else the grid's footprint; a bare point cloud
        has neither and is drawn as interpolated. It was never masked at all: the
        edge it showed was the convex hull of the points, which spans every
        concavity of a domain that is not convex (ledger 314).

        Replaces a `plot()` that built a BARE `go.Figure` and immediately called
        `.show(renderer="browser")` -- so it returned None, could not be modified
        or embedded, forced a browser window, and silently dropped the house
        template, `scrollZoom` and `dragmode="pan"` that every other myflopy
        figure carries. Use `.clipped_fig()` for a variant -- that one is a METHOD
        precisely because it takes arguments, so it rebuilds.

        Assembled once and cached: `Picture` requires repeated access to return
        the SAME figure, so `surface.fig.update_layout(...)` then
        `surface.show()` acts on one figure.
        """

        if not self._assembled:
            self._fig = self.clipped_fig()
            self._assembled = True
        return self._fig

    def clipped_fig(self, *, surface=None, clip: bool | None = None):
        """The surface figure, masked to :attr:`region` -- ``clip=`` or the grid's footprint.

        ``clip=None`` masks whenever there is a region, which is what :attr:`fig`
        draws. ``clip=False`` draws the bare interpolation, whose only edge is the
        convex hull of the points. ``clip=True`` insists, and raises when there is
        nothing to mask to. ``surface=`` draws another array on this surface's
        pixel grid, masked the same way -- it used to be discarded whenever
        ``clip=True``.

        A new figure on every call. :attr:`fig` is the cached one.

        Parameters
        ----------
        surface : numpy.ndarray, optional
            The z values to draw on this surface's pixel grid, shape
            ``(resolution, resolution)``, row 0 at the north edge. When masking,
            any other shape raises ``ValueError``. With ``clip=False`` it is
            drawn unchecked. The default ``None`` is :attr:`surface`, the
            unmasked interpolation, which the mask then clips.
        clip : bool or None, default None
            ``None`` masks to :attr:`region` when there is one and draws
            unmasked when there is not (a bare point cloud). ``True`` masks and
            raises ``ValueError`` when there is no region. ``False`` never
            masks. Masking that removes EVERY finite pixel raises
            ``ValueError`` naming both extents. A surface with no finite pixel
            to begin with is drawn empty, without an error.

        Returns
        -------
        myflopy.viz.Fig
            A figure holding the one surface trace.
        """

        surface = self.surface if surface is None else surface
        if clip is None:
            clip = self.region is not None
        if clip:
            surface = self._masked(surface, self._required_region())

        figure = f.Fig()
        figure.add_trace(self.surface_trace(surface=surface))
        return figure

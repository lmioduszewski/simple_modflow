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
from rasterio.io import MemoryFile
from rasterio.mask import mask
from rasterio.transform import from_origin
from scipy.interpolate import RBFInterpolator, griddata
from scipy.spatial import QhullError
from shapely.geometry import Polygon, mapping

from myflopy import viz as f
from myflopy._logging import get_logger
from myflopy.viz import Picture
from myflopy.modflow.mf6.headsplus import HeadsPlus as Hp
from myflopy.modflow.utils.datatypes.readers import read_shp_gpkg

logger = get_logger(__name__)


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
        :param hds: Optional, HeadsPlus object corresponding to model, will get from model if not provided
        :param model: mf6 SimulationBase model
        :param layer: defaults to 0
        :param per: stress period number, O-based index; will take precedence over kstpkper if provided
        :param kstpkper: tuple of time step and period for surface
        :param resolution: defaults to 1000
        :param use_rbf: defaults to False
        :param surf_type: defaults to 'hds', can be 'lyr' or 'hds'. 'lyr' returns just model surfaces
        :param clip: optional shapely Polygon to clip interpolated surface to
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
        if per is not None:
            self.kstpkper = self.model.kstpkper[per]
        else:
            self.kstpkper = kstpkper
        self.layer = layer
        self.resolution = resolution
        self.neighbors = 10
        self.colorscale = 'Earth_r'  # reverse earth so blue is low nums
        self._clip = None
        self._clipped_cells = None

        self.clip = clip

    @property
    def clip(self):
        """clip is a shapely Polygon or Path to clip interpolated surface to"""
        return self._clip

    @clip.setter
    def clip(self, val):
        """Set the clip region from ``None``, a shapefile/GeoPackage ``Path``, or a Polygon."""

        if val is None:
            self._clip = None
        elif isinstance(val, Path):
            clip = read_shp_gpkg(val)
            assert len(clip) == 1, 'clip must be a single polygon'
            clip = clip.union_all()
            assert isinstance(clip, shp.geometry.Polygon), 'clip must be a shapely Polygon'
            self._clip = clip
        elif isinstance(val, Polygon):
            self._clip = val


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
            if self.clip:
                self._clipped_cells = self.vor.get_vor_cells_as_series(self.clip).to_list()
                xs = pd.Series(self.vor.centroids_x).loc[self._clipped_cells].to_numpy()
            else:
                xs = np.array(self.vor.centroids_x)
            self._xs = xs
        return self._xs

    @property
    def ys(self):
        """Cell-centroid y-coordinates (clipped to the clip region when set; cached)."""

        if self._ys is None:
            if self.clip:
                self._clipped_cells = self.vor.get_vor_cells_as_series(self.clip).to_list()
                ys = pd.Series(self.vor.centroids_y).loc[self._clipped_cells].to_numpy()
            else:
                ys = np.array(self.vor.centroids_y)
            self._ys = ys
        return self._ys

    @property
    def zs(self):
        """Per-cell surface values to interpolate: heads (``hds``) or a layer elevation (``lyr``).

        Head surfaces read the selected ``kstpkper``/``layer`` (inactive-cell
        sentinels become NaN); layer surfaces read ``gdf_topbtm``. Cached.
        """

        cells = slice(None) if self._clipped_cells is None else self._clipped_cells
        if self._zs is None:
            if self.surf_type == 'hds':
                """zs of the self.hds HeadPlus oject for a given layer at a certain stress-and-time period"""
                zs = self.hds.all_heads.loc[idxx[self.kstpkper, self.layer, cells], :].values
                zs[zs > 10_000] = np.nan # remove large zs, which would be inactive cells
                self._zs = zs
            if self.surf_type == 'lyr':
                zs = self.vor.gdf_topbtm.loc[cells, self.layer].values
                self._zs = zs
        return self._zs

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
        """The model's :class:`HeadsPlus` reader, opened from its ``.hds`` output (cached)."""

        if self._hds is None:
            if self.model:
                hds = Hp(
                    hds_path=self.model.model_output_folder_path / f'{self.model.name}.hds',
                    vor=self.vor
                )
                self._hds = hds
            else:
                raise ValueError(
                    "no heads file was given and no model was given to find one on; "
                    "pass hds= or model=."
                )
        assert isinstance(self._hds, Hp), 'the hds property must be a HeadsPlus instance'
        return self._hds

    @property
    def griddata_interp(self):
        """interpolated surface using scipy griddata"""
        if self._griddata_interp is None:
            zis = griddata(
                points=self.xys,
                values=self.zs,
                xi=self.xy_meshgrid,
                method='cubic')
            self._griddata_interp = zis.squeeze()
        return self._griddata_interp

    @property
    def rbf_interp(self):
        """interpolated surface using scipy RBFInterpolator"""
        if self._rbf_interp is None:
            coords = np.column_stack((self.xs, self.ys))
            xis = self.xy_meshgrid[0].ravel()
            yis = self.xy_meshgrid[1].ravel()
            xyis = np.column_stack((xis, yis))
            interpolator = RBFInterpolator(
                coords,
                self.zs,
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
            coords = np.column_stack((self.xs, self.ys))
            interpolator = LinearNDInterpolator(coords, self.zs)
            grid_z = interpolator(self.xy_meshgrid).reshape(self.xy_meshgrid[0].shape)
            self._linearND_interp = grid_z
        return self._linearND_interp

    @property
    def cloughTocher2D_interp(self):
        """interpolated surface using scipy CloughTocher2DInterpolator"""
        if self._cloughTocher2D_interp is None:
            from scipy.interpolate import CloughTocher2DInterpolator
            coords = np.column_stack((self.xs, self.ys))
            interpolator = CloughTocher2DInterpolator(coords, self.zs)
            grid_z = interpolator(self.xy_meshgrid).reshape(self.xy_meshgrid[0].shape)
            self._cloughTocher2D_interp = grid_z
        return self._cloughTocher2D_interp

    @property
    def transform(self):
        """
        Returns an affine transformation matrix for the given surface
        """
        xs = self.xs
        ys = self.ys

        xmin, xmax = np.min(xs), np.max(xs)
        ymin, ymax = np.min(ys), np.max(ys)
        xres = (xmax - xmin) / (self.resolution - 1)
        yres = (ymax - ymin) / (self.resolution - 1)

        return from_origin(xmin, ymax, xres, yres)

    @property
    def memfile(self):
        """An in-memory single-band GeoTIFF (:class:`MemoryFile`) of the interpolated surface."""

        grid_z = self.surface
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
                transform=self.transform) as dst:
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
                logger.warning(
                    'interpolator %r is not one of %s; ignoring it',
                    val, self._interpolators,
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

    def clip_raster_with_polygon(self, polygon: Polygon = None, buffer=None):
        """
        Clips the in-memory raster with a vector polygon and returns the clipped raster data.

        Parameters:
        memfile (MemoryFile): In-memory GeoTIFF raster.
        polygon (shapely.geometry.Polygon): Polygon to use for clipping.

        Returns:
        tuple: Clipped raster data array and the updated transform.
        """
        if self.clip and polygon is None:
            polygon = self.clip
        polygon = self.vor.gdf_vorPolys.union_all() if polygon is None else polygon
        polygon = polygon.buffer(buffer) if buffer else polygon
        with self.memfile.open() as dataset:
            shapes = [mapping(polygon)]

            # Clip the raster with the polygon
            clipped_image, clipped_transform = mask(dataset, shapes, crop=False)

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
        """
        Saves the clipped raster data to a GeoTIFF file.

        Parameters:
        clipped_image (numpy.ndarray): Clipped raster data array.
        clipped_transform (Affine): Transform for the clipped raster.
        meta (dict): Metadata of the original raster dataset.
        output_tif (str): Path to the output GeoTIFF file.

        Returns:
        None
        """

        polygon_clip = self.clip if polygon_clip is None else polygon_clip
        if not any([clipped_image, clipped_transform, meta]):
            clipped_image, clipped_transform, meta = self.clip_raster_with_polygon(polygon=polygon_clip, buffer=buffer)
            clip = pd.DataFrame(clipped_image[0]).replace(0.0, np.nan)
            clipped_image = clip.to_numpy().reshape((1, self.resolution, self.resolution))

        # Update the metadata with the new transform and dimensions
        meta.update({
            "driver": "GTiff",
            "height": clipped_image.shape[1],
            "width": clipped_image.shape[2],
            "transform": clipped_transform
        })

        # Write the clipped raster to a GeoTIFF file
        with rasterio.open(output_tif, "w", **meta) as dst:
            dst.write(clipped_image)

    def plot_heatmap(self):
        """Show the clipped interpolated surface as a Plotly heatmap."""

        clipped_image = self.clip_raster_with_polygon(polygon=self.clip)[0]
        fig = f.Fig()
        fig.add_heatmap(z=clipped_image[0])
        fig.show()


    def surface_trace(self, *, surface=None, colorscale=None, **kwargs):
        """Return a Plotly ``go.Surface`` trace for this interpolated surface.

        The single builder shared by :meth:`plot` and higher-level multi-surface
        views (e.g. ``myflopy.layers.LayerBuildResult.surface_3d``). ``surface``
        overrides the interpolated z-array (defaults to :attr:`surface`); any
        extra keyword arguments pass straight through to ``go.Surface``.
        """
        surface = self.surface if surface is None else surface
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
        """The 3-D surface figure.

        Replaces a `plot()` that built a BARE `go.Figure` and immediately called
        `.show(renderer="browser")` -- so it returned None, could not be modified
        or embedded, forced a browser window, and silently dropped the house
        template, `scrollZoom` and `dragmode="pan"` that every other myflopy
        figure carries. Use `.clipped_fig()` for the clipped variant -- that one
        is a METHOD precisely because it takes arguments, so it rebuilds.

        Assembled once and cached: `Picture` requires repeated access to return
        the SAME figure, so `surface.fig.update_layout(...)` then
        `surface.show()` acts on one figure.
        """

        if not self._assembled:
            self._fig = self.clipped_fig(clip=False)
            self._assembled = True
        return self._fig

    def clipped_fig(self, *, surface=None, clip: bool = False):
        """The surface figure, optionally clipped to the configured polygon."""

        if clip:
            surface = self.clip_raster_with_polygon()[0].squeeze()
            surface = surface.astype(float)
            surface[surface == 0] = np.nan

        figure = f.Fig()
        figure.add_trace(self.surface_trace(surface=surface))
        return figure

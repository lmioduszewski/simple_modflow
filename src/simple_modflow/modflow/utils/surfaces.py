from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.grid.voronoi import VoronoiGridPlus as Vor
    from simple_modflow.modflow.mf6.simulation.base import SimulationBase

import numpy as np
from scipy.interpolate import griddata, RBFInterpolator
import figs as f
import pandas as pd
import rasterio
from pathlib import Path
import plotly.graph_objs as go
from shapely.geometry import Polygon
from rasterio.io import MemoryFile
from rasterio.transform import from_origin
from rasterio.mask import mask
from pandas import IndexSlice as idxx
from shapely.geometry import mapping
from simple_modflow.modflow.mf6.headsplus import HeadsPlus as Hp
import shapely as shp
from simple_modflow.modflow.utils.datatypes.readers import read_shp_gpkg
from scipy.spatial import cKDTree


class InterpolatedSurface:

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
        except:
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
        return self._interpolator

    @interpolator.setter
    def interpolator(self, val):
        if val is not None:
            if val not in self._interpolators:
                print(f'interpolator must be one of {self._interpolators}')
            else:
                self._interpolator = val

    @property
    def kstpkper(self):
        if self._kstpkper is None:
            if self.model:
                self._kstpkper = self.model.kstpkper[0]
            else:
                print('no model defined')
                raise ValueError
        return self._kstpkper

    @kstpkper.setter
    def kstpkper(self, val):
        if val is not None:
            if len(val) == 2:
                assert isinstance(val, tuple), (
                    'kstpkper must be a tuple of (timestep, period)'
                )
                self._kstpkper = val
            else:
                print('kstpkper must be a tuple of (timestep, period)')

    @property
    def xs(self):
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
        self._zs = val

    @property
    def xys(self):
        if self._xys is None:
            xys = list(zip(self.xs, self.ys))
            self._xys = xys
        return self._xys

    @property
    def xy_meshgrid(self):
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
        if self._hds is None:
            if self.model:
                hds = Hp(
                    hds_path=self.model.model_output_folder_path / f'{self.model.name}.hds',
                    vor=self.vor
                )
                self._hds = hds
            else:
                print('no heads file defined')
                raise ValueError
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

        grid_z = self.surface
        try:
            crs = self.vor.crs
        except:
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
                print(f'interpolator must be one of {self._interpolators}')
        if self.use_rbf is False:
            try:
                return self.griddata_interp
            except:
                print('error with griddata interpolation, using rbf')
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

        clipped_image = self.clip_raster_with_polygon(polygon=self.clip)[0]
        fig = f.Fig()
        fig.add_heatmap(z=clipped_image[0])
        fig.show()


    def plot(self, surface=None, clip=False, renderer='browser'):
        surface = self.surface if surface is None else surface
        """
        plot surface using plotly, defaults to griddata_interp
        :param surface: surface to plot, ex. self.griddata_interp or self.rbf_interp
        :return: plots surface to browser
        """
        if clip:
            surface = self.clip_raster_with_polygon()[0].squeeze()
            surface = surface.astype(float)
            surface[surface == 0] = np.nan

        xg, yg = self.xy_meshgrid

        fig = go.Figure()
        fig.add_surface(
            x=xg,
            y=yg,
            z=surface,
            colorscale=self.colorscale
        )
        fig.show(renderer=renderer)

def rasterize_points(
    df: pd.DataFrame,
    xcol="x",
    ycol="y",
    zcol="z",
    cellsize=20.0,
    crs="EPSG:4326",            # set to your projected CRS
    bounds=None,                 # (xmin, ymin, xmax, ymax); if None uses data extent
    method="griddata",                # "idw" | "griddata" | "rbf"
    idw_power=2.0,
    idw_k=16,                    # neighbors for IDW
    idw_radius=None,             # optional search radius (same units as coords)
    rbf_kernel="linear",         # e.g., "linear", "thin_plate_spline"
    rbf_smoothing=0.0,           # >0 smooths
    nodata=np.nan,
    out_path=None
):
    """
    Interpolates a raster from a DataFrame of x, y, z - points.
    Returns (Z, transform) if out_path is None, else saves GeoTIFF and returns out_path.
    df must have columns x, y, z. Coordinates should be in a projected CRS (meters/feet).

    Example usage:

        rasterize_points(
            df,
            cellsize=20.0,
            crs="EPSG:2927",
            method="griddata",   # or "rbf"
            out_path="surface_z.tif"
        )
    """
    # --- prep points ---
    x = df[xcol].to_numpy(dtype=float)
    y = df[ycol].to_numpy(dtype=float)
    z = df[zcol].to_numpy(dtype=float)

    if bounds is None:
        xmin, xmax = np.nanmin(x), np.nanmax(x)
        ymin, ymax = np.nanmin(y), np.nanmax(y)
        # tiny padding to fully include edges
        pad = 1e-9
        xmin -= pad; ymin -= pad; xmax += pad; ymax += pad
    else:
        xmin, ymin, xmax, ymax = bounds

    # build grid (top-left origin)
    xi = np.arange(xmin, xmax + cellsize, cellsize)
    yi = np.arange(ymax, ymin - cellsize, -cellsize)  # descending so row 0 is top
    cols = xi.size
    rows = yi.size
    XI, YI = np.meshgrid(xi, yi)
    grid_pts = np.column_stack([XI.ravel(), YI.ravel()])
    transform = from_origin(xmin, ymax, cellsize, cellsize)

    # --- interpolation methods ---
    if method.lower() == "idw":
        tree = cKDTree(np.column_stack([x, y]))
        # query k neighbors (fast, vectorized). SciPy >=1.6 supports workers=-1
        dist, idx = tree.query(grid_pts, k=min(idw_k, len(x)), workers=-1)

        # handle 1 neighbor vs many neighbors consistently
        if dist.ndim == 1:
            dist = dist[:, None]
            idx = idx[:, None]

        if idw_radius is not None:
            mask = dist > idw_radius
        else:
            mask = np.zeros_like(dist, dtype=bool)

        with np.errstate(divide="ignore"):
            w = 1.0 / np.power(dist, idw_power)
        w[mask] = 0.0
        # exact hits (dist=0) -> set weight=inf so they dominate
        w[np.isinf(w)] = 1e12

        vals = z[idx]
        wsum = w.sum(axis=1)
        out = np.full(grid_pts.shape[0], np.nan, dtype=float)
        good = wsum > 0
        out[good] = (w[good] * vals[good]).sum(axis=1) / wsum[good]

        Z = out.reshape(rows, cols)

    elif method.lower() == "griddata":
        pts = np.column_stack([x, y])
        # linear inside hull
        Z = griddata(pts, z, (XI, YI), method="linear")
        # nearest fill for holes / outside hull
        Zn = griddata(pts, z, (XI, YI), method="nearest")
        if Z is None:
            Z = Zn
        else:
            fill = np.isnan(Z)
            Z[fill] = Zn[fill]

    elif method.lower() == "rbf":
        pts = np.column_stack([x, y])
        rbf = RBFInterpolator(pts, z, kernel=rbf_kernel, smoothing=rbf_smoothing)
        Z = rbf(grid_pts).reshape(rows, cols)
    else:
        raise ValueError("method must be 'idw', 'griddata', or 'rbf'")

    # --- write GeoTIFF or return ---
    if out_path:
        profile = {
            "driver": "GTiff",
            "height": rows,
            "width": cols,
            "count": 1,
            "dtype": "float32",
            "crs": crs,
            "transform": transform,
            "nodata": nodata if np.isnan(nodata) else float(nodata),
            "compress": "deflate",
            "tiled": True,
            "predictor": 3,
        }
        with rasterio.open(out_path, "w", **profile) as dst:
            dst.write(Z.astype("float32"), 1)
        return out_path
    else:
        return Z, transform


if __name__ == '__main__':
    import pickle
    model_path_v7b_et = Path(r"C:\Users\lukem\mf6\cum7bET\cum7bET.model")
    with open(model_path_v7b_et, 'rb') as file:
        model7b: SimulationBase = pickle.load(file)

    surf = InterpolatedSurface(
        model=model7b,
        per=75,
        interpolator='linearND',
        resolution=200
    )
    surf.plot(clip=True)

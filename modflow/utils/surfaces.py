from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
    from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor

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
            kstpkper: tuple = None,
            resolution: int = 1000,
            use_rbf: bool = False,
            surf_type: str = 'hds'
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
        :param kstpkper: tuple of time step and period for surface
        :param resolution: defaults to 1000
        :param use_rbf: defaults to False
        :param surf_type: defaults to 'hds', can be 'lyr' or 'hds'. 'lyr' returns just model surfaces
        """
        self.model = model
        self.vor = self.model.vor if vor is None else vor
        self.surf_type = surf_type
        self._xs = xs
        self._ys = ys
        self._zs = zs
        self._xys = None
        self.use_rbf = use_rbf
        self._griddata_interp = None
        self._rbf_interp = None
        self._meshgrid = None
        self._hds = hds
        self.kstpkper = kstpkper
        self.layer = layer
        self.resolution = resolution
        self.neighbors = 10
        self.colorscale = 'Earth_r' # reverse earth so blue is low nums

    @property
    def xs(self):
        if self._xs is None:
            xs = np.array(self.vor.centroids_x)
            self._xs = xs
        return self._xs

    @property
    def ys(self):
        if self._ys is None:
            ys = np.array(self.vor.centroids_y)
            self._ys = ys
        return self._ys

    @property
    def zs(self):
        if self.surf_type == 'hds':
            """zs of the self.hds HeadPlus oject for a given layer at a certain stress-and-time period"""
            zs = self.hds.all_heads.loc[idxx[self.kstpkper, self.layer, :], :].values
            self._zs = zs
        if self.surf_type == 'lyr':
            zs = self.vor.gdf_topbtm.loc[:, self.layer].values
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
                np.linspace(ys.min(), ys.max(), self.resolution))
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
    def transform(self):
        """
        Returns an affine transformation matrix for the given surface
        """
        xs = self.xs
        ys = self.ys

        xmin, ymin = np.min(xs), np.min(ys)
        xmax, ymax = np.max(xs), np.max(ys)
        xres = (xmax - xmin) / self.resolution
        yres = (ymax - ymin) / self.resolution
        transform = from_origin(xmin, ymin, xres, -yres)

        return transform

    @property
    def memfile(self):

        if self.use_rbf:
            grid_z = self.rbf_interp
        else:
            grid_z = self.griddata_interp

        memfile = MemoryFile()
        with memfile.open(
                driver='GTiff',
                height=self.resolution,
                width=self.resolution,
                count=1,
                dtype=grid_z.dtype,
                crs=self.vor.crs,
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
        try:
            return self.griddata_interp
        except:
            return self.rbf_interp

    def clip_raster_with_polygon(self, polygon: Polygon = None):
        """
        Clips the in-memory raster with a vector polygon and returns the clipped raster data.

        Parameters:
        memfile (MemoryFile): In-memory GeoTIFF raster.
        polygon (shapely.geometry.Polygon): Polygon to use for clipping.

        Returns:
        tuple: Clipped raster data array and the updated transform.
        """
        polygon = self.vor.gdf_vorPolys.union_all() if polygon is None else polygon
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
            polygon_clip: Polygon = None
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

        if not any([clipped_image, clipped_transform, meta]):
            clipped_image, clipped_transform, meta = self.clip_raster_with_polygon(polygon=polygon_clip)
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

        clipped_image = self.clip_raster_with_polygon()[0]
        fig = f.Fig()
        fig.add_heatmap(z=clipped_image[0])
        fig.show()

    def plot(self, surface=None):
        """
        plot surface using plotly, defaults to griddata_interp
        :param surface: surface to plot, ex. self.griddata_interp or self.rbf_interp
        :return: plots surface to browser
        """
        surface = self.surface if surface is None else surface
        fig = go.Figure()
        fig.add_surface(
            z=surface,
            x=self.projected_xys[0],
            y=self.projected_xys[1],
            colorscale=self.colorscale
        )
        fig.show(renderer='browser')

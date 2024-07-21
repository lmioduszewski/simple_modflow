import numpy as np
from scipy.interpolate import griddata, RBFInterpolator
import figs as f
import pandas as pd
import rasterio
from pathlib import Path
from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor
from simple_modflow.modflow.mf6.headsplus import HeadsPlus as Hp
from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
import plotly.graph_objs as go
from pandas import IndexSlice as idxx

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
            resolution: int = 1000
    ):
        """
        Base class for interpolated surfaces.
        :param xs: optional, array of x coordinates. If not provided will get from vor
        :param ys: optional, array of y coordinates. If not provided will get from vor
        :param zs:
        :param vor: voronoi grid object corresponding to model
        :param hds: Optional, HeadsPlus object corresponding to model, will get from model if not provided
        :param model: mf6 SimulationBase model
        :param layer: defaults to 0
        :param kstpkper: tuple of time step and period for surface
        :param resolution: defaults to 1000
        """
        self.model = model
        self.vor = vor
        self._xs = xs
        self._ys = ys
        self._zs = zs
        self._xys = None
        self._griddata_interp = None
        self._rbf_interp = None
        self._meshgrid = None
        self._hds = hds
        self.kstpkper = kstpkper
        self.layer = layer
        self.resolution = resolution
        self.neighbors = 10

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
        zs = self.hds.all_heads.loc[idxx[self.kstpkper, self.layer, :], :].values
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

    def plot(self, surface=None):
        if surface is None:
            surface = self.griddata_interp
        fig = go.Figure()
        fig.add_surface(z=surface)
        fig.show(renderer='browser')


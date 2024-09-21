from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.mfsimbase import SimulationBase

from simple_modflow.modflow.utils.surfaces import InterpolatedSurface
from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor
import figs as f
import pandas as pd
import geopandas as gpd


class ModelSurface:

    def __init__(
            self,
            model: SimulationBase = None,
    ):
        self._model = model
        self._surfaces = {}
        self.nper = model.gwf.modeltime.nper
        self.nstp = model.gwf.modeltime.nstp

    def hds(self, layer=0, kstpkper: tuple = None, plot: bool = False, **kwargs):
        kstpkper = (self.nstp[0] - 1, 0) if kstpkper is None else kstpkper
        surf = InterpolatedSurface(model=self.model, layer=layer, kstpkper=kstpkper, **kwargs)
        if plot:
            surf.plot()
        return surf

    def lyr(self, layer=0, plot: bool = False, **kwargs):
        surf = InterpolatedSurface(model=self.model, layer=layer, surf_type='lyr', **kwargs)
        if plot:
            surf.plot()
        return surf

    @property
    def model(self):
        return self._model


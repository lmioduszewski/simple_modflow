from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase


from myflopy.modflow.mf6.grid.interpolated_surface import InterpolatedSurface


class ModelSurface:

    def __init__(
            self,
            model: SimulationBase = None,
    ):
        """Bind an interpolated-surface factory to ``model`` (caching its timing)."""

        self._model = model
        self._surfaces = {}
        self.nper = model.gwf.modeltime.nper
        self.nstp = model.gwf.modeltime.nstp
        self.kstpkper = model.kstpkper

    def hds(self, per=None, layer=0, kstpkper: tuple = None, plot: bool = False, **kwargs) -> InterpolatedSurface:
        """
        define a head-based surface for the given model for a specific stress period and layer
        :param per: integer ranging from 0 (first stress period) to the last stress period
        :param layer: defaults to 0, which is Layer 1. It's a zero-based indexing system...0, 1, 2, 3, and so on...
        :param kstpkper: tuple of (timestep, stress period). Or can just define per. kstpkper will then
        befined based on per as an index to model.kstpkper. ex. mondel.kstpkper[0] for the first stress period
        :param plot: defaults to False, if true hds surface will automatically plot
        :param kwargs: other args to pass to InterpolatedSurface class
        :return: an InterpolatedSurface object
        """
        kstpkper = (self.nstp[0] - 1, 0) if kstpkper is None else kstpkper
        if per is not None:
            assert isinstance(per, int) and per >= 0, 'per arg must be an integer greater than or equal to 0'
            kstpkper = self.kstpkper[per]
        surf = InterpolatedSurface(model=self.model, layer=layer, kstpkper=kstpkper, **kwargs)
        if plot:
            surf.plot()
        return surf

    def lyr(self, layer=0, plot: bool = False, **kwargs):
        """An interpolated layer-elevation surface for ``layer`` (optionally plotted)."""

        surf = InterpolatedSurface(model=self.model, layer=layer, surf_type='lyr', **kwargs)
        if plot:
            surf.plot()
        return surf

    @property
    def model(self):
        """The bound MODFLOW simulation."""

        return self._model


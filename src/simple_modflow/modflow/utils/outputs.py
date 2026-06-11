from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.simulation.base import SimulationBase

from simple_modflow.modflow.mf6.budget import LakBudget, SFRBudget
import pandas as pd


class LakStageOutput:
    """Lake stage output accessor hosted under ``model.outputs.lak``."""

    def __init__(self, model: "SimulationBase"):
        self.model = model

    def get(self):
        """Return lake stage output reshaped to ``(nper, nlakes)``."""

        nlakes = self.model.lak.nlakes.data
        stages = self.model.lak.output.stage()
        return stages.get_alldata().reshape(-1, nlakes)


class SFRStageOutput:
    """Stream stage output accessor hosted under ``model.outputs.sfr``."""

    def __init__(self, model: "SimulationBase"):
        self.model = model

    def get(self):
        """Return a reach-by-saved-timestep dataframe of stream stage output."""

        nreaches = self.model.sfr.nreaches.data
        reader = self.model.sfr.output.stage()
        kstpkpers = [tuple(value) for value in reader.get_kstpkper()]
        stages = reader.get_alldata()
        stage_data = stages.reshape(len(kstpkpers), nreaches).transpose()
        frame = pd.DataFrame(stage_data, columns=kstpkpers)
        frame.index.name = "reaches"
        frame.columns.name = "kstpkper"
        return frame


class LakOutputData:
    """Namespace for LAK-specific output helpers."""

    def __init__(self, model: SimulationBase):
        self.model = model

    @property
    def bud(self):
        """Return the LAK package-output budget accessor."""

        return LakBudget(self.model)

    @property
    def stage(self):
        """Return the LAK stage output accessor.

        The returned helper exposes ``get()`` so existing code can follow the
        same pattern as other output helpers while keeping stage logic under the
        ``outputs`` namespace.
        """

        return LakStageOutput(self.model)


class SFROutputData:
    """Namespace for SFR-specific output helpers."""

    def __init__(self, model: SimulationBase):
        self.model = model

    @property
    def bud(self):
        """Return the SFR package-output budget accessor."""

        return SFRBudget(self.model)

    @property
    def stage(self):
        """Return the SFR stage output accessor.

        The returned helper exposes ``get()`` so existing code can follow the
        same pattern as other output helpers while keeping stage logic under the
        ``outputs`` namespace.
        """

        return SFRStageOutput(self.model)

class UzfOutputData:
    """Namespace for UZF-specific output helpers."""

    def __init__(self, model: SimulationBase):
        self.model = model
        self._ifno_to_cellid = None

    @property
    def ifno_to_cellid(self):
        if self._ifno_to_cellid is None:
            if self.model.modelgrid.grid_type == 'vertex':
                ifno_to_cellid = pd.DataFrame(self.model.uzf.packagedata.get_data().cellid.tolist())
                ifno_to_cellid.index.name = 'ifno'
                ifno_to_cellid.columns = ['layer', 'cellid']
                self._ifno_to_cellid = ifno_to_cellid
            else:
                raise ValueError("Model grid type is not vertex. You'll have to figure this out manually...")
        return self._ifno_to_cellid

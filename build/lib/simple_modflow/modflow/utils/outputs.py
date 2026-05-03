from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
    from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor

from simple_modflow.modflow.mf6.budget import LakBudget, LakStage, SFRBudget, SFRStage
import pandas as pd


class LakOutputData:
    """base class for all kinds of lak output data"""

    def __init__(self, model: SimulationBase):
        self.model = model

    @property
    def bud(self):
        return LakBudget(self.model)

    @property
    def stage(self):
        return LakStage(self.model)


class SFROutputData:
    """base class for all kinds of sfr output data"""

    def __init__(self, model: SimulationBase):
        self.model = model

    @property
    def bud(self):
        return SFRBudget(self.model)

    @property
    def stage(self):
        return SFRStage(self.model)

class UzfOutputData:
    """base class for all kinds of uzf output data"""

    def __init__(self, model: SimulationBase):
        self.model = model
        self._ifno_to_cellid = None

    @property
    def ifno_to_cellid(self):
        if self._ifno_to_cellid is None:
            if self.model.modelgrid.grid_type == 'vertex':
                ifno_to_cellid = pd.DataFrame(self.model.gwf.uzf.packagedata.get_data().cellid.tolist())
                ifno_to_cellid.index.name = 'ifno'
                ifno_to_cellid.columns = ['layer', 'cellid']
                self._ifno_to_cellid = ifno_to_cellid
            else:
                raise ValueError("Model grid type is not vertex. You'll have to figure this out manually...")
        return self._ifno_to_cellid
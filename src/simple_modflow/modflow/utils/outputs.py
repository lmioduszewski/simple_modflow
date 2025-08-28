from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
    from simple_modflow.modflow.mf6.voronoiplus import VoronoiGridPlus as Vor

from simple_modflow.modflow.mf6.budget import LakBudget, LakStage, SFRBudget, SFRStage


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


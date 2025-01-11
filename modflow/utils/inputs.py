from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.mfsimbase import SimulationBase

import pandas as pd
from simple_modflow.modflow.utils.datatypes.choros import Choro


class RchInput:

    def __init__(
            self,
            model: SimulationBase,
    ):
        self.model = model

    def plot_input(
            self,
            per: tuple = None,
            to_in_per_yr: bool = False
    ):
        model = self.model
        vor = model.vor

        df = pd.DataFrame(model.gwf.rch.stress_period_data.data[per])
        df['cellid'] = df['cellid'].apply(lambda x: x[1])  # extract cell numbers
        c_hov = df.set_index('cellid').reindex(list(range(vor.ncpl)), fill_value=0)  # insert missing cell nums
        if to_in_per_yr:
            c_hov = c_hov * 365 * 12
        c_hov = c_hov.to_dict(orient='list')
        Choro(vor=vor, custom_hover=c_hov, custom_zs=c_hov[list(c_hov.keys())[0]]).plot()


class DrnInput:

    def __init__(
            self,
            model: SimulationBase,
    ):

        self.model = model

    def plot_input(
            self,
            per: tuple = None,
    ):
        model = self.model
        vor = model.vor

        drn_cells = pd.DataFrame(model.gwf.drn.stress_period_data.data[per])['cellid'].apply(lambda x: x[1]).to_list()
        vor.show_selected_cells(drn_cells)

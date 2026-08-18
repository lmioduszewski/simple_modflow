from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase

import pandas as pd

from myflopy.modflow.utils.datatypes.choros import Choro

from myflopy._logging import get_logger

logger = get_logger(__name__)


class Inputs:
    def __init__(self, model: SimulationBase):
        """Bind the legacy input-explorer namespace to ``model``."""

        self.model = model

    @property
    def rch(self):
        """Recharge input helper for the bound model."""

        return RchInput(self.model)

    @property
    def uzf(self):
        """UZF input helper for the bound model."""

        return UzfInput(self.model)


class UzfInput:
    def __init__(self, model: SimulationBase):
        """Bind the UZF input helper to ``model``."""

        self.model = model

    def finf(self, per: int | str = 0, multiplier: float = 12 * 30):
        """
        Calculate the modified infiltration rate for UZF (Unsaturated Zone Flow) package
        based on the given period and multiplier. The output is adjusted to include all
        model cells even if they have zero infiltration.

        :param per: The stress period index for which the infiltration rates need
            to be calculated. Defaults to 0.
        :type per: int
        :param multiplier: A scaling factor applied to the infiltration rate.
            Defaults to 12 * 30 - converts ft/day to in/month.
        :type multiplier: float
        :return: A pandas Series containing the adjusted infiltration rates for all
            model cells, indexed by cell ID.
        :rtype: pd.Series
        """
        uzf = self.model.gwf.uzf
        ncpl = self.model.modelgrid.ncpl
        packagedata = pd.DataFrame(uzf.packagedata.get_data()).set_index('ifno')

        if isinstance(per, str):
            if per.lower() == 'all':
                dfs = []
                for k, v in uzf.perioddata.data.items():
                    df = pd.DataFrame(v).loc[:, ['ifno', 'finf']]
                    df.columns = ['ifno', k]
                    df = df.set_index('ifno')
                    dfs.append(df)

                perioddata = pd.concat(dfs, axis=1).reset_index()
            else:
                raise ValueError(f'Invalid period value: {per}. "all" is only accepted for this argument.')

        elif isinstance(per, int):
            perioddata = pd.DataFrame(uzf.perioddata.data[per]).loc[:, ['ifno', 'finf']]

        # create 'cell' column - converts ifno to cellid
        perioddata['cell'] = perioddata.ifno.apply(lambda x: packagedata.loc[x].cellid[1])
        perioddata = perioddata.drop(columns='ifno')

        # create dataframe with just cellid and uzf infiltration data
        # then reindex so all model cells are included even if zero finf
        finf = perioddata.set_index('cell').reindex(list(range(ncpl)), fill_value=0)
        logger.info('scaling UZF infiltration rates by %s', multiplier)
        finf = finf * multiplier

        return finf

    def plot(self, per: int = 0, multiplier: float = 12 * 30, **kwargs):
        """Plot the UZF infiltration (``finf``) choropleth for one period."""

        finf = self.finf(per, multiplier).finf
        return self.model.cor(per=per, custom_zs=finf.to_list(), **kwargs)

class RchInput:

    def __init__(
            self,
            model: SimulationBase,
    ):
        """Bind the RCH input helper to ``model`` (default period 0)."""

        self.model = model
        self._df = None
        self.per = 0  # default stress period

    def df(self, per: int|str = None):
        """The RCH stress-period-data table for ``per`` indexed by cell number."""

        per = self.per if per is None else per
        df = pd.DataFrame(self.model.gwf.rch.stress_period_data.data[per])
        df['cellid'] = df['cellid'].apply(lambda x: x[1])  # extract cell numbers
        df = df.set_index('cellid')
        return df

    def plot(
            self,
            per: tuple = None,
            multiplier: float | int = 1
    ):
        """Plot the recharge choropleth for one period (filling zero-recharge cells)."""

        model = self.model
        vor = model.vor

        df = pd.DataFrame(model.gwf.rch.stress_period_data.data[per])
        df['cellid'] = df['cellid'].apply(lambda x: x[1])  # extract cell numbers
        df = df.set_index('cellid')
        c_hov = df.reindex(list(range(vor.ncpl)), fill_value=0)  # insert missing cell nums
        c_hov = c_hov * multiplier
        c_hov = c_hov.to_dict(orient='list')
        return Choro(vor=vor, custom_hover=c_hov, custom_zs=c_hov[list(c_hov.keys())[0]])


class DrnInput:

    def __init__(
            self,
            model: SimulationBase,
    ):
        """Bind the DRN input helper to ``model``."""

        self.model = model

    def plot(
            self,
            per: tuple = None,
    ):
        """Highlight the DRN cells on the grid for one stress period."""

        model = self.model
        vor = model.vor

        drn_cells = pd.DataFrame(model.gwf.drn.stress_period_data.data[per])['cellid'].apply(lambda x: x[1]).to_list()
        vor.show_selected_cells(drn_cells)

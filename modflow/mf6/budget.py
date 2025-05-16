from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from simple_modflow.modflow.mf6.mfsimbase import SimulationBase
    from .voronoiplus import VoronoiGridPlus as Vor
    from flopy.utils.binaryfile import CellBudgetFile, HeadFile

import flopy
import pandas as pd
import numpy as np
from figs import Fig
from pandas import IndexSlice as idxx


class Budget:

    def __init__(self, model: SimulationBase, gwf_package: str = None):
        """
        Class for basic gwf budget operations. Advanced packages, like SFR and LAK, have their own classes.
        Though you can still access their GWF budget data from this class.
        :param model: model file
        :param gwf_package: package name for budget to retrieve, ex. DRN, RCH, GHB
        """
        self.model = model
        self._gwf_package = None
        self._df = None

        self.gwf_package = gwf_package

    @property
    def gwf_package(self):
        """package name for budget to retrieve, ex. DRN, RCH, GHB"""
        return self._gwf_package

    @gwf_package.setter
    def gwf_package(self, pck: str):
        if pck is not None:
            assert any([pck.upper() in pack for pack in self.types]), f'package: [{pck}] not included in budget file'
        self._gwf_package = pck

    @property
    def types(self):
        """returns types of records in the model budget file"""
        types = self.model.gwf.output.budget().get_unique_record_names()
        types = [pack.astype(str) for pack in types]
        return types

    @property
    def budget(self):
        """gets the entire model budget if package is not provided. If package is set, returns
        raw budget data for that package. Call the df property to get a concatenated Dataframe
        of all data for a package"""
        bud: CellBudgetFile = self.model.gwf.output.budget()
        if self.gwf_package is None:
            return bud
        else:
            return bud.get_data(text=self.gwf_package)

    @property
    def df(self):
        """Gets a concatenated Dataframe of all budget data for a package for all stress periods"""
        if self.gwf_package is None:
            return 'Must define a package first.'
        if self._df is None:

            buds = self.budget
            kstpkper_all = self.model.kstpkper
            dfs = []
            # assert len(buds) == self.model.nper, 'number of budgets does not match number of periods'
            for per, bud in enumerate(buds):
                kstpkper = kstpkper_all[per]
                df = pd.DataFrame(bud)
                df['kstpkper'] = [kstpkper for _ in range(len(df))]
                dfs.append(df)

            df = pd.concat(dfs)
            df = df.set_index([df.columns[0], 'kstpkper'])
            self._df = df

        return self._df


class LakBudget:

    def __init__(self, model: SimulationBase):
        self.model = model

    @property
    def types(self):
        types = self.model.gwf.lak.output.budget().get_unique_record_names()
        types = [bud.astype(str) for bud in types]
        return types

    @property
    def budget(self):
        bud: CellBudgetFile = self.model.gwf.lak.output.budget()
        return bud

    def get(self, bud_type: str = None, return_df: bool = True):
        """returns the budget data for a given budget type. Default behavior is to return
        a concatenated Dataframe of all budget data for a package for all stress periods.
        Can return just the list of numpy record arrays by setting return_df=False."""
        if bud_type is None:
            print('must define type')
            return self.types
        bud = self.budget.get_data(text=bud_type)
        if return_df is False:
            return bud
        elif return_df is True:
            if len(bud) != len(self.model.kstpkper):
                print('Warning: length of periods and budget array list do not match')
                print('cannot automatically assign periods to each array')
                return bud
            else:
                for idx, per in enumerate(self.model.kstpkper):
                    df = pd.DataFrame(bud[idx])
                    df['kstpkper'] = [per for _ in range(len(df))]
                    bud[idx] = df
                return pd.concat(bud)


class LakStage:
    """class to work with lak stage output data"""

    def __init__(self, model: SimulationBase):
        self.model = model
        self.nlakes = self.model.gwf.lak.nlakes.data

    def get(self):
        """returns lake stage output data"""
        stages: HeadFile = self.model.gwf.lak.output.stage()
        stg_data = stages.get_alldata().reshape(-1, self.nlakes)
        return stg_data


class SFRBudget:
    """class to work with sfr budget data"""

    def __init__(self, model: SimulationBase):
        self.model = model
        self.nreaches = self.model.gwf.sfr.nreaches.data

    @property
    def types(self):
        types = self.model.gwf.sfr.output.budget().get_unique_record_names()
        types = [bud.astype(str) for bud in types]
        return types

    @property
    def budget(self):
        bud: CellBudgetFile = self.model.gwf.sfr.output.budget()
        return bud

    def get(self, bud_type: str = None, return_df: bool = True):
        """returns the budget data for a given budget type. Default behavior is to return
        a concatenated Dataframe of all budget data for a package for all stress periods.
        Can return just the list of numpy record arrays by setting return_df=False."""
        if bud_type is None:
            print('must define type')
            return self.types
        bud = self.budget.get_data(text=bud_type)
        if return_df is False:
            return bud
        elif return_df is True:
            if len(bud) != len(self.model.kstpkper):
                print('Warning: length of periods and budget array list do not match')
                print('cannot automatically assign periods to each array')
                return bud
            else:
                for idx, per in enumerate(self.model.kstpkper):
                    df = pd.DataFrame(bud[idx])
                    df['kstpkper'] = [per for _ in range(len(df))]
                    bud[idx] = df
                return pd.concat(bud)

    def plot_flows(
            self,
            kstpkper: tuple = None,
            cfd_to_cfs: bool = False,
            cfd_to_gpm: bool = False,
            html: str = None
    ):
        """
        plot flows for all streams in given stress period. Can convert to cfs or gpm using boolean arguments.
        Defaults to first stress period if none is specified.
        :param kstpkper: stress period and timestep to plot
        :param cfd_to_cfs: convert to cfs
        :param cfd_to_gpm: convert to gpm
        :param html: if you want to save to html, provide html name
        :return:
        """
        if kstpkper is None:
            kstpkper = [self.model.kstpkper[0]]
        gpm = 448.8 if cfd_to_gpm else 1
        cfs = (24 * 60 * 60) if cfd_to_cfs or cfd_to_gpm else 1
        flows: pd.DataFrame = self.get('flow')
        flows = flows.set_index(['kstpkper', 'node'], drop=True)
        stream_flows = []
        for riv in self.model.sfr_input.stream_reaches:
            for stpper in kstpkper:
                riv_flows = flows.loc[idxx[stpper, riv], :]
                riv_flows = riv_flows[riv_flows.q < 0].q / cfs
                riv_flows: pd.Series = riv_flows * -1 * gpm
                stream_flows.append(riv_flows)
        fig = Fig()
        for i, flows in enumerate(stream_flows):
            flows = flows.reset_index().drop('kstpkper', axis=1)
            fig.add_scattergl(x=flows.node, y=flows.q, name=f'stream {i}')
        if html is not None:
            fig.write_html(file=html)
        fig.show()


class SFRStage:
    """class to work with sfr stage output data"""

    def __init__(self, model: SimulationBase):
        self.model = model
        self.nreaches = self.model.gwf.sfr.nreaches.data

    def get(self):
        """gets a Dataframe of all stages for each reach for each stress period"""
        stages: np.ndarray = self.model.gwf.sfr.output.stage().get_alldata()
        stg_data = stages.reshape(self.model.nper, self.nreaches).transpose()
        stg_data = pd.DataFrame(stg_data)
        stg_data.index.name = 'reaches'
        return stg_data


class DRNBudget:

    def __init__(self, model: SimulationBase = None):
        self.model = model

    def plot_choro(self, per: int = 0, zmax=None):

        model = self.model
        kstpkper = model.kstpkper[per]
        drn_df = model.bud('drn').df
        drn_flows = drn_df.loc[idxx[:, kstpkper], :].q.droplevel(1) * -1
        if zmax is None:
            zmax = drn_flows.max()
        node = drn_flows.reset_index().drop_duplicates('node')
        node['node'] = node['node'] - 1
        node.set_index('node', inplace=True)
        full_idx = range(model.vor.ncpl)
        # Set index to include all cells, and fill empty values with zero
        drn_flows = node.reindex(full_idx, fill_value=0)
        model.choro(per=per, custom_zs=drn_flows.q.to_list(), zmin=0, zmax=zmax).plot()


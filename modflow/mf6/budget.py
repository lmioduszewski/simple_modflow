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
from pathlib import Path
from simple_modflow.modflow.utils.datatypes.readers import read_shp_gpkg
from simple_modflow.modflow.mf6.boundaries import Boundaries
from plotly.colors import DEFAULT_PLOTLY_COLORS as colors
from itertools import cycle


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
    def gwf_package(self, pkg: str):
        if pkg is not None:
            assert any([pkg.upper() in pack for pack in self.types]), f'package: [{pkg}] not included in budget file'
        self._gwf_package = pkg

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

    def plot_budget_obs(
            self,
            model: SimulationBase = None,
            shp_gpkg: Path = None,
            q: str = 'q',
            plot_fig=True,
            name_field: str = 'name'
    ):
        """
        Plots the water budget observation values for specified areas (polygons) defined in a
        geospatial file, computes the total observed inflow or outflow, and calculates
        percent contributions of all obs areas to the total flow.

        :param model: The simulation model. Defaults to None. If None, uses self.model.
        :type model: SimulationBase
        :param shp_gpkg: The path to the shapefile or GeoPackage file that defines
            the observation areas.
        :type shp_gpkg: Path
        :param q: The column name representing the flow parameter of interest.
            Defaults to 'q'.
        :type q: str
        :param plot_fig: Flag to indicate whether to plot the results as a figure.
            Defaults to True.
        :type plot_fig: bool
        :param name_field: text string corresponding to name field in shp_gpkg. Defaults to 'name'.
        :return: A DataFrame containing flow values and percentage contributions for
            observation areas and stress periods.
        :rtype: pd.DataFrame
        """
        model = self.model if model is None else model

        full_idx = range(model.vor.ncpl)  # valid cell ids
        # remove duplicates
        df = self.df[~self.df.index.duplicated(keep='first')]
        per0 = model.kstpkper[0]

        # get valid cell ids for the package budget
        pkg_cells = self.df.loc[idxx[:, per0], :].index.get_level_values(0).to_list()
        # subtract one to transform to 0-based since node values are 1-based in the budget file
        pkg_cells = [i - 1 for i in pkg_cells]
        b_obs = read_shp_gpkg(shp_gpkg)
        crs = b_obs.crs.to_epsg()
        # find cells intersecting each observation area in the shp_gpkg
        b_obs['cells'] = Boundaries(model, model.vor, shp_gpkg=shp_gpkg, crs=crs).intersections

        def drop_non_pkg_cells(x):
            new_cells = []
            for c in x:
                if c in pkg_cells:
                    new_cells.append(c)
            return new_cells

        # only keep valid package cell ids for each observation area
        b_obs['cells'] = b_obs['cells'].apply(lambda x: drop_non_pkg_cells(x))

        bdict = {}
        bdict['total'] = {}

        # create dictionary to hold info for each observation area defined by shp_gpkg
        for obs in b_obs.iterrows():
            bdict[obs[1].loc[name_field]] = {}
            cs = obs[1].loc['cells']
            # get flows for each observation area for each stress period
            for kstpkper in model.kstpkper:
                try:
                    # get all flows for the given self.gwf_package in this stress period
                    flows = df.loc[idxx[:, kstpkper], :][q].droplevel(1) * -1
                    tot_flow = flows.sum() / (24 * 60 * 60)
                    bdict['total'][kstpkper] = tot_flow
                    flows = flows.reindex(full_idx, fill_value=0)
                    # add up flows for this observation area for this stress period
                    bdict[obs[1].loc[name_field]][kstpkper] = flows.loc[cs].sum() / (24 * 60 * 60)
                except:
                    # skip stress period if it doesn't exist
                    continue

        df = pd.DataFrame.from_dict(bdict)
        # calculate percent of total package flows for all given obs areas for each stress period
        df['percent_of_tot'] = df.iloc[:, 1:].reset_index(drop=True).sum(axis=1).values / df['total']
        df = df.reset_index(drop=True)
        if plot_fig:
            fig = Fig()
            for col in df.columns:
                fig.add_scattergl(
                    x=df.index,
                    y=df[col],
                    name=col,
                )
            fig.update_yaxes(range=[0, 4])
            fig.show()
        return df

    @classmethod
    def multimodel_plot_budget_obs(
            cls,
            models: list[SimulationBase],
            model_package: str = None,
            shp_gpkg: Path = None,
            q: str = 'q',
            name_field: str = 'name',
            plot_fig: bool = True
    ):
        """
        Generate and plot the budget observations from multiple simulation models.

        This method aggregates budget observation data generated by different
        simulation models' specified packages. These observations are plotted
        with different colors for each model. The plotted figure can be displayed
        or returned for further manipulation.

        :param models: List of simulation models that are instances of `SimulationBase`
        :type models: list[SimulationBase]
        :param model_package: Specific model package identifier from which budget
            observations are to be extracted for each model. Defaults to `'drn'`
            if not explicitly provided.
        :type model_package: str, optional
        :param shp_gpkg: File path to a shapefile or geopackage containing
            spatial reference information for the budget observations. Used to
            link geographic data with the model's budget data.
        :type shp_gpkg: Path, optional
        :param q: Column name in the budget dataframe for identifying the flow
            rate or budget quantity to analyze. Defaults to `'q'`.
        :type q: str, optional
        :param name_field: Column name in the shapefile or geopackage used to
            uniquely identify geographic entities related to the budget data.
            Defaults to 'name'.
        :type name_field: str, optional
        :param plot_fig: Flag indicating whether or not the method should display
            the generated plot figure. Defaults to 'True'.
        :type plot_fig: bool, optional
        :return: A plotly Figure object containing the budget observation plot.
        :rtype: Figure
        """

        model_package = 'drn' if model_package is None else model_package
        obs_dict = {}
        color_cycle = cycle(colors)
        fig = Fig()

        # get DataFrames of budget obs for each model and add to dict
        for i, model in enumerate(models):
            obs_dict[model.name] = model.bud(model_package).plot_budget_obs(
                shp_gpkg=shp_gpkg, q=q, name_field=name_field, plot_fig=False)

        # add budget obs traces for each model to figure
        for modelname, df in obs_dict.items():
            color = next(color_cycle)
            for col in df.columns:
                fig.add_scattergl(
                    x=df.index,
                    y=df[col],
                    name=f'{col}-{modelname}',
                    line=dict(color=color)
                )

        if plot_fig:
            fig.show()

        return fig


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
    """
    Represents a module for DRN Budget management and plot visualization.

    This class is designed to handle budget data related to Drain (DRN) packages,
    and provide a visualization of the data using choropleth plots.

    :ivar model: The simulation model object that contains budget data and plotting methods.
    :type model: SimulationBase
    """

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

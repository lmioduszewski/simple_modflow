"""Budget access wrappers for groundwater and package-level MF6 outputs.

This module keeps the long-lived public budget API in one place while
delegating newer table-building and plotting responsibilities to helper modules:

- :mod:`myflopy.modflow.mf6.budget_tables`
- :mod:`myflopy.modflow.mf6.budget_plotting`

The public classes here are intentionally compatibility-oriented wrappers. They
retain the established `myflopy` calling style while sharing common
logic with the newer lazy-loading and grouped-model layers.
"""

from __future__ import annotations
from typing import TYPE_CHECKING
import warnings

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase
    from flopy.utils.binaryfile import CellBudgetFile, HeadFile

import pandas as pd
import numpy as np
from pandas import IndexSlice as idxx
from pathlib import Path

from myflopy.modflow.mf6.budget_plotting import (
    multimodel_plot_budget_obs as _multimodel_plot_budget_obs,
    plot_budget_obs as _plot_budget_obs,
    plot_drn_choropleth as _plot_drn_choropleth,
)
from myflopy.modflow.mf6.budget_tables import (
    budget_df as _budget_df,
    budget_types as _budget_types,
    package_output_budget as _package_output_budget,
    package_output_df as _package_output_df,
    package_output_types as _package_output_types,
    raw_budget as _raw_budget,
)


class Budget:
    """Access groundwater-model budget records through a model-aware wrapper.

    Parameters
    ----------
    model
        Parent simulation model providing access to the MF6 budget file and the
        grid/period metadata needed to shape results.
    gwf_package
        Optional groundwater-package filter such as ``"drn"``, ``"ghb"``, or
        ``"rch"``. When provided, :attr:`df` and :attr:`budget` return package-
        specific data; when omitted, :attr:`budget` returns the raw budget file
        reader.

    Notes
    -----
    The dataframe view exposed by :attr:`df` uses zero-based cell indexing for
    supported groundwater packages so the returned node ids line up with FloPy's
    grid conventions and the rest of ``myflopy``.
    """
    def __init__(self, model: SimulationBase, gwf_package: str = None):
        """Create one groundwater-budget accessor for ``model``."""
        self.model = model
        self._gwf_package = None
        self._df = None

        self.gwf_package = gwf_package

    @property
    def gwf_package(self):
        """Package name filter used when retrieving one groundwater budget."""
        return self._gwf_package

    @gwf_package.setter
    def gwf_package(self, pkg: str):
        """Validate and store the optional groundwater package filter."""

        if pkg is not None:
            assert any([pkg.upper() in pack for pack in self.types]), f'package: [{pkg}] not included in budget file'
        self._gwf_package = pkg

    @property
    def types(self):
        """Return record types present in the model's groundwater budget file."""
        return _budget_types(self.model)

    @property
    def budget(self):
        """Return the raw budget reader or one package's raw record list.

        Returns
        -------
        CellBudgetFile or list
            The full MF6 budget reader when :attr:`gwf_package` is unset, or
            the raw list of record arrays for the selected package when
            :attr:`gwf_package` is set.
        """
        return _raw_budget(self.model, self.gwf_package)

    @property
    def df(self):
        """Return one concatenated dataframe for the selected groundwater package.

        The dataframe is built lazily and cached on first access. For supported
        packages such as ``drn``, ``ghb``, and ``rch``, the node identifiers are
        normalized to zero-based indexing exactly once.
        """
        if self.gwf_package is None:
            raise ValueError(
                "gwf_package is required before requesting Budget.df. "
                "Use model.bud('drn').df or a similar package-specific accessor."
            )
        if self._df is None:
            self._df = _budget_df(self.model, self.gwf_package)

        return self._df

    def plot_budget_obs(
            self,
            model: SimulationBase = None,
            shp_gpkg: Path = None,
            q: str = 'q',
            plot_fig=True,
            return_fig=False,
            name_field: str = 'name',
            times = None,
            multiplier = 24 * 60 * 60,
            y_range = [0, 4]
    ):
        """Aggregate package budget flows by observation polygon.

        Parameters
        ----------
        model
            Optional model override. Defaults to the accessor's parent model.
        shp_gpkg
            Polygon layer defining the observation areas to summarize.
        q
            Budget value column to aggregate, typically ``"q"``.
        plot_fig
            If ``True``, display a line plot of the resulting observation-area
            totals.
        return_fig
            If ``True``, return the constructed figure object instead of only
            returning the dataframe.
        name_field
            Attribute field in ``shp_gpkg`` used to label each observation area.
        times
            Optional explicit x-axis datetime labels.
        multiplier
            Unit-conversion divisor applied after summing flows. The historical
            default converts cubic feet per day to cubic feet per second.
        y_range
            Optional plot y-axis range used when a figure is built.

        Returns
        -------
        pandas.DataFrame or figs.Fig
            Observation-area totals by stress period, or a figure when
            ``return_fig=True``.
        """
        model = self.model if model is None else model
        return _plot_budget_obs(
            self,
            model=model,
            shp_gpkg=shp_gpkg,
            q=q,
            plot_fig=plot_fig,
            return_fig=return_fig,
            name_field=name_field,
            times=times,
            multiplier=multiplier,
            y_range=y_range,
        )

    @classmethod
    def multimodel_plot_budget_obs(
            cls,
            models: list[SimulationBase],
            model_package: str = 'drn',
            shp_gpkg: Path = None,
            q: str = 'q',
            name_field: str = 'name',
            times: pd.DatetimeIndex = None,
            plot_fig: bool = True
    ):
        """Plot observation-area groundwater budgets for several models.

        Parameters
        ----------
        models
            Models whose budget observations should be plotted together.
        model_package
            Groundwater budget package to summarize, such as ``"drn"``.
        shp_gpkg
            Polygon layer defining the observation areas.
        q
            Budget value column to aggregate.
        name_field
            Attribute field used to label observation areas.
        times
            Optional datetime labels for the x-axis.
        plot_fig
            If ``True``, display the combined figure before returning it.

        Returns
        -------
        figs.Fig
            Combined line plot containing traces for each model/observation-area
            combination.
        """

        return _multimodel_plot_budget_obs(
            models,
            model_package='drn' if model_package is None else model_package,
            shp_gpkg=shp_gpkg,
            q=q,
            name_field=name_field,
            times=times,
            plot_fig=plot_fig,
        )


class LakBudget:
    """Access the LAK package's package-output budget file for one model."""

    def __init__(self, model: SimulationBase):
        """Create one LAK budget accessor for ``model``."""
        self.model = model

    @property
    def types(self):
        """Return record types available in the LAK budget output file."""
        return _package_output_types(self.model.lak)

    @property
    def budget(self):
        """Return the raw FloPy LAK package-output budget reader."""
        return _package_output_budget(self.model.lak)

    def get(self, bud_type: str = None, return_df: bool = True):
        """Return LAK budget data for one output record type.

        Parameters
        ----------
        bud_type
            Output record type to retrieve.
        return_df
            If ``True``, return one concatenated dataframe across all periods.
            If ``False``, return the raw list of record arrays.
        """
        if bud_type is None:
            raise ValueError(
                "bud_type is required when requesting LAK budget data. "
                f"Available types: {self.types}"
            )
        bud = self.budget.get_data(text=bud_type)
        if return_df is False:
            return bud
        elif return_df is True:
            try:
                return _package_output_df(self.model, self.model.lak, bud_type=bud_type)
            except ValueError:
                warnings.warn(
                    "Length of LAK budget periods and record arrays do not match; "
                    "returning the raw record list instead of a concatenated dataframe.",
                    stacklevel=2,
                )
                return bud


class LakStage:
    """Compatibility wrapper for lake stage output access.

    The canonical home for this behavior is now ``model.outputs.lak.stage``.
    This class remains for backward compatibility with older call sites.
    """

    def __init__(self, model: SimulationBase):
        """Create one LAK stage accessor for ``model``."""
        self.model = model
        self.nlakes = self.model.lak.nlakes.data

    def get(self):
        """Return lake stage output reshaped to ``(nper, nlakes)``."""
        return self.model.outputs.lak.stage.get()


class SFRBudget:
    """Access the SFR package's package-output budget file for one model."""

    def __init__(self, model: SimulationBase):
        """Create one SFR budget accessor for ``model``."""
        self.model = model
        self.nreaches = self.model.sfr.nreaches.data

    @property
    def types(self):
        """Return record types available in the SFR budget output file."""
        return _package_output_types(self.model.sfr)

    @property
    def budget(self):
        """Return the raw FloPy SFR package-output budget reader."""
        return _package_output_budget(self.model.sfr)

    def get(self, bud_type: str = None, return_df: bool = True):
        """Return SFR budget data for one output record type.

        Parameters
        ----------
        bud_type
            Output record type to retrieve.
        return_df
            If ``True``, return one concatenated dataframe across all periods.
            If ``False``, return the raw list of record arrays.
        """
        if bud_type is None:
            raise ValueError(
                "bud_type is required when requesting SFR budget data. "
                f"Available types: {self.types}"
            )
        bud = self.budget.get_data(text=bud_type)
        if return_df is False:
            return bud
        elif return_df is True:
            try:
                return _package_output_df(self.model, self.model.sfr, bud_type=bud_type)
            except ValueError:
                warnings.warn(
                    "Length of SFR budget periods and record arrays do not match; "
                    "returning the raw record list instead of a concatenated dataframe.",
                    stacklevel=2,
                )
                return bud

    def plot_flows(
            self,
            kstpkper: tuple = None,
            cfd_to_cfs: bool = False,
            cfd_to_gpm: bool = False,
            html: str = None
    ):
        """Plot flow-by-reach values for one or more stress periods.

        Parameters
        ----------
        kstpkper
            One or more explicit ``(kstp, kper)`` selectors. Defaults to the
            first stress period.
        cfd_to_cfs
            If ``True``, convert cubic feet per day to cubic feet per second.
        cfd_to_gpm
            If ``True``, convert cubic feet per day to gallons per minute.
            This implies the cubic-feet-per-second conversion as part of the
            existing historical workflow.
        html
            Optional output HTML path for the figure.
        """
        if kstpkper is None:
            kstpkper = [self.model.kstpkper[0]]
        gpm = 448.8 if cfd_to_gpm else 1
        cfs = (24 * 60 * 60) if cfd_to_cfs or cfd_to_gpm else 1
        flows: pd.DataFrame = self.get('flow')
        flows = flows.set_index(['kstpkper', 'node'], drop=True)
        stream_flows = []
        reaches = pd.DataFrame(self.model.gwf.sfr.packagedata.get_data())["rno"].tolist()
        for stpper in kstpkper:
            riv_flows = flows.loc[idxx[stpper, reaches], :]
            riv_flows = riv_flows[riv_flows.q < 0].q / cfs
            riv_flows: pd.Series = riv_flows * -1 * gpm
            stream_flows.append(riv_flows)
        from myflopy.viz import Fig

        fig = Fig()
        for i, flows in enumerate(stream_flows):
            flows = flows.reset_index().drop('kstpkper', axis=1)
            fig.add_scattergl(x=flows.node, y=flows.q, name=f'stream {i}')
        if html is not None:
            fig.write_html(file=html)
        fig.show()


class SFRStage:
    """Compatibility wrapper for stream stage output access.

    The canonical home for this behavior is now ``model.outputs.sfr.stage``.
    This class remains for backward compatibility with older call sites.
    """

    def __init__(self, model: SimulationBase):
        """Create one SFR stage accessor for ``model``."""
        self.model = model
        self.nreaches = self.model.sfr.nreaches.data

    def get(self):
        """Return a reach-by-period dataframe of stream stage output."""
        return self.model.outputs.sfr.stage.get()


class DRNBudget:
    """Convenience wrapper for DRN-specific budget visualization.

    This is a small compatibility wrapper around the older DRN choropleth
    plotting entry point.
    """

    def __init__(self, model: SimulationBase = None):
        """Create one DRN budget plotting helper for ``model``."""
        self.model = model

    def plot_choro(self, per: int = 0, zmax=None):
        """Plot DRN flows for one zero-based stress-period index."""
        _plot_drn_choropleth(self.model, per=per, zmax=zmax)

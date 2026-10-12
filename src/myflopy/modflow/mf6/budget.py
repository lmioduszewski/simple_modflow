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

if TYPE_CHECKING:

    from myflopy.modflow.mf6.simulation.base import SimulationBase

from pathlib import Path

import numpy as np
import pandas as pd

from myflopy.modflow.mf6.budget_plotting import (
    multimodel_plot_budget_obs as _multimodel_plot_budget_obs,
)
from myflopy.modflow.mf6.budget_plotting import (
    plot_budget_obs as _plot_budget_obs,
)
from myflopy.modflow.mf6.budget_tables import (
    budget_df as _budget_df,
)
from myflopy.modflow.mf6.budget_tables import (
    budget_types as _budget_types,
)
from myflopy.modflow.mf6.budget_tables import (
    package_output_budget as _package_output_budget,
)
from myflopy.modflow.mf6.budget_tables import (
    package_output_df as _package_output_df,
)
from myflopy.modflow.mf6.budget_tables import (
    package_output_types as _package_output_types,
)
from myflopy.modflow.mf6.budget_tables import (
    raw_budget as _raw_budget,
)
from myflopy.modflow.mf6.package_registry import one_package

#: Model-budget records MF6 names for the process rather than for the package
#: that wrote them, so a package-name lookup can never find them.
_BUDGET_RECORD_ALIASES = {
    "ssm": "SOURCE-SINK MIX",
}


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
        """Resolve and store the optional groundwater/transport package filter."""

        if pkg is None:
            self._gwf_package = None
            return

        available = [str(record).strip() for record in self.types]
        wanted = str(pkg).strip()
        if not any(wanted.upper() in record.upper() for record in available):
            # MF6 names a few model-budget records for the PROCESS rather than
            # the package that wrote them, so the package's own name never
            # appears in the file. SSM is the one users reach for: its record is
            # "SOURCE-SINK MIX", which left model.bud("ssm") unreachable on
            # EVERY transport model (measured 2026-07-27).
            alias = _BUDGET_RECORD_ALIASES.get(wanted.lower())
            if alias is None or not any(
                alias.upper() in record.upper() for record in available
            ):
                raise ValueError(
                    f"package: [{pkg}] not included in budget file. "
                    f"Available records: {available}"
                )
            wanted = alias
        self._gwf_package = wanted

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
            y_range = None
    ):
        """Aggregate package budget flows by observation polygon.

        Reads :attr:`df`, so the accessor needs a package (``model.bud("drn")``;
        a bare ``model.bud()`` raises ``ValueError``). Each row is one saved time
        of THIS budget file, labelled by its own zero-based ``kstpkper`` (ledger
        301; it read the HEAD file's clock until ledger 476). Flows are
        NEGATED (water leaving the aquifer, a drain's outflow, reads positive),
        every record on a cell is ADDED (a flow is extensive, ledger 297; it kept
        the first until ledger 476), and the sum is divided by
        ``multiplier``.

        Parameters
        ----------
        model : SimulationBase, optional
            The model whose grid the areas are matched against. Default ``None``:
            the accessor's own model. The budget values and their clock always
            come from the ACCESSOR's model (:attr:`df`).
        shp_gpkg : str or pathlib.Path
            A polygon file (its FIRST layer) whose rows are the observation areas,
            read through the one footprint reader (``mf.lak(lakes=)`` and
            ``LayerStack.excavate`` read through it too). Required despite the
            ``None`` default. The areas are REPROJECTED to the grid's CRS; a file
            carrying no CRS is taken to be in it. A multi-part polygon is one area.
            A cell belongs to an area that COVERS it by area (the footprint rule:
            an area drawn along cell faces reaches only the cells inside it), in
            EVERY layer: a record below layer 0 counts under the area over its
            column (it was left out, read as another cell number, until ledger
            482). An area covering no cell RAISES ``ValueError`` naming both extents and
            both CRSs -- before ledger 476 a file in another CRS read 0.0 for
            every area in silence. A feature that is not a polygon, or one with no
            area, raises ``ValueError``.
        q : str, default 'q'
            The column of :attr:`df` to sum: ``"q"``, the flow (``df`` also holds
            ``node2`` and whatever auxiliary columns the record carries). A name
            ``df`` lacks raises ``KeyError`` naming the columns it has.
        plot_fig : bool, default True
            Draw a line per area (and ``total``) and ``show()`` it.
        return_fig : bool, default False
            Return the figure instead of the table; it is built even when
            ``plot_fig`` is False. With both False no figure is built.
        name_field : str, default 'name'
            The attribute column of ``shp_gpkg`` naming each area; its values
            (as text) become the table's columns, after ``total``. A missing
            column, or two areas of one name, raises ``ValueError``.
        times : sequence of datetime-like or pandas.DatetimeIndex, optional
            The figure's x values, one per table row, in place of each row's
            ``"kstp, kper"``; converted with ``pandas.DatetimeIndex``. A length
            other than the table's raises ``ValueError``. The table itself keeps
            its ``kstpkper`` index. Default ``None``.
        multiplier : float, default 86400
            The DIVISOR applied to every sum (``24 * 60 * 60``: ft3/d to ft3/s in a
            feet-and-days model). Pass ``1`` to keep the model's own flow units.
        y_range : list of float, optional
            The figure's y-axis range. Default ``None`` lets plotly fit the data
            (it was a fixed ``[0, 4]`` until ledger 476).

        Returns
        -------
        pandas.DataFrame or Fig
            The table -- indexed by the budget file's ``kstpkper``, a ``total``
            column over every record of the package, then one column per area --
            or, with ``return_fig=True``, the figure (a :class:`myflopy.viz.Fig`).
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

        Each model is summarized by :meth:`plot_budget_obs` at its defaults
        (``multiplier=86400``, flows negated so outflow reads positive), and
        every model's ``total`` and area lines are drawn in one colour per model.

        Parameters
        ----------
        models : list of SimulationBase
            The models to compare. Each must have run, carry a budget record
            matching ``model_package`` and a ``name``, which suffixes its traces
            (``"<area>-<model name>"``).
        model_package : str, default 'drn'
            The budget record to summarize, read through ``model.bud(...)`` -- a
            case-insensitive substring of a record name. ``None`` also means
            ``"drn"``.
        shp_gpkg : str or pathlib.Path
            The observation-area polygons, read as :meth:`plot_budget_obs` reads
            them (reprojected to each model's grid CRS; an area covering no cell
            raises). Required despite the ``None`` default.
        q : str, default 'q'
            The budget column to sum, as in :meth:`plot_budget_obs`.
        name_field : str, default 'name'
            The attribute column naming each area.
        times : sequence of datetime-like or pandas.DatetimeIndex, optional
            The x values for every model's lines, one per row of each model's
            table (another length raises ``ValueError``); default ``None`` draws
            each row at its budget file's ``"kstp, kper"``.
        plot_fig : bool, default True
            ``show()`` the figure before returning it.

        Returns
        -------
        Fig
            One :class:`myflopy.viz.Fig` holding every model's traces; returned
            whatever ``plot_fig`` says.
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
    """Access one LAK package's package-output budget file for one model."""

    def __init__(self, model: SimulationBase, package: str | None = None):
        """Create one LAK budget accessor for ``model``'s LAK package ``package``.

        ``package`` is the package's NAME, or ``None`` for the model's only LAK
        package; with several, reading raises ``SeveralPackagesError`` naming
        ``model.outputs.lak(package=...)`` (ledger 352).
        """
        self.model = model
        self.package = package

    def _package(self):
        """The LAK package this budget reads (``package_registry.one_package``)."""
        return one_package(
            self.model, "lak", self.package, spelling="model.outputs.lak(package={name!r})"
        )

    @property
    def types(self):
        """Return record types available in the LAK budget output file."""
        return _package_output_types(self._package())

    @property
    def budget(self):
        """Return the raw FloPy LAK package-output budget reader."""
        return _package_output_budget(self._package())

    def get(self, bud_type: str = None, return_df: bool = True):
        """Return LAK budget data for one output record type.

        Parameters
        ----------
        bud_type : str
            The record to read from the package's own budget file -- one of
            :attr:`types` (``"GWF"``, ``"STORAGE"``, ``"FROM-MVR"``, ...). Matched
            by FloPy case-insensitively as a SUBSTRING of the record names, so
            ``"gwf"`` and ``"stor"`` work and the first record containing it is
            read. Required despite the ``None`` default: ``None`` raises
            ``ValueError`` listing :attr:`types`, and a string no record contains
            raises FloPy's ``ValueError``.
        return_df : bool, default True
            ``True``: one concatenated dataframe -- always a dataframe, one row
            per record entry, columns ``node``/``node2``/``q`` plus the record's
            auxiliary (``FLOW-AREA``, ``VOLUME``, ...) exactly as MF6 wrote them
            (1-based ids, not shifted), and a ``kstpkper`` column holding the
            zero-based ``(kstp, kper)`` the LAK budget file itself recorded.
            Only the literal ``False`` returns the raw list of FloPy record
            arrays instead; any other value (``0``, ``None``) is read as ``True``.

        Notes
        -----
        Until ledger 301 this warned and returned the RAW LIST when the record
        count disagreed with ``model.kstpkper``, silently breaking the dataframe
        contract above; the comparison was against the HEAD file's save clock,
        which an advanced package's budget file has never been written on. The
        clock is now read from that file, so there is nothing left to disagree.
        Ask for the raw list with ``return_df=False``.
        """
        if bud_type is None:
            raise ValueError(
                "bud_type is required when requesting LAK budget data. "
                f"Available types: {self.types}"
            )
        if return_df is False:
            return self.budget.get_data(text=bud_type)
        return _package_output_df(self.model, self._package(), bud_type=bud_type)


class LakStage:
    """Compatibility wrapper for lake stage output access.

    The canonical home for this behavior is now ``model.outputs.lak.stage``.
    This class remains for backward compatibility with older call sites.
    """

    def __init__(self, model: SimulationBase):
        """Create one LAK stage accessor for ``model`` (its only LAK package)."""
        self.model = model
        self.nlakes = one_package(self.model, "lak").nlakes.data

    def get(self):
        """Return ``model.outputs.lak.stage.get()``: one row per saved record, one
        column per lake, ``NaN`` where a lake is dry or inactive (ledger 329)."""
        return self.model.outputs.lak.stage.get()


class SFRBudget:
    """Access one SFR package's package-output budget file for one model."""

    def __init__(self, model: SimulationBase, package: str | None = None):
        """Create one SFR budget accessor for ``model``'s SFR package ``package``.

        ``package`` is the package's NAME, or ``None`` for the model's only SFR
        package; with several, reading raises ``SeveralPackagesError`` naming
        ``model.outputs.sfr(package=...)`` (ledger 352).
        """
        self.model = model
        self.package = package

    def _package(self):
        """The SFR package this budget reads (``package_registry.one_package``)."""
        return one_package(
            self.model, "sfr", self.package, spelling="model.outputs.sfr(package={name!r})"
        )

    @property
    def nreaches(self):
        """The package's reach count (``NREACHES``)."""
        return self._package().nreaches.data

    @property
    def types(self):
        """Return record types available in the SFR budget output file."""
        return _package_output_types(self._package())

    @property
    def budget(self):
        """Return the raw FloPy SFR package-output budget reader."""
        return _package_output_budget(self._package())

    def get(self, bud_type: str = None, return_df: bool = True):
        """Return SFR budget data for one output record type.

        Parameters
        ----------
        bud_type : str
            The record to read from the package's own budget file -- one of
            :attr:`types` (``"GWF"``, ``"STORAGE"``, ``"FROM-MVR"``, ...). Matched
            by FloPy case-insensitively as a SUBSTRING of the record names, so
            ``"gwf"`` and ``"stor"`` work and the first record containing it is
            read. Required despite the ``None`` default: ``None`` raises
            ``ValueError`` listing :attr:`types`, and a string no record contains
            raises FloPy's ``ValueError``.
        return_df : bool, default True
            ``True``: one concatenated dataframe -- always a dataframe, one row
            per record entry, columns ``node``/``node2``/``q`` plus the record's
            auxiliary (``FLOW-AREA``, ``VOLUME``, ...) exactly as MF6 wrote them
            (1-based ids, not shifted), and a ``kstpkper`` column holding the
            zero-based ``(kstp, kper)`` the SFR budget file itself recorded.
            Only the literal ``False`` returns the raw list of FloPy record
            arrays instead; any other value (``0``, ``None``) is read as ``True``.

        Notes
        -----
        Until ledger 301 this warned and returned the RAW LIST when the record
        count disagreed with ``model.kstpkper``, silently breaking the dataframe
        contract above; the comparison was against the HEAD file's save clock,
        which an advanced package's budget file has never been written on. The
        clock is now read from that file, so there is nothing left to disagree.
        Ask for the raw list with ``return_df=False``.
        """
        if bud_type is None:
            raise ValueError(
                "bud_type is required when requesting SFR budget data. "
                f"Available types: {self.types}"
            )
        if return_df is False:
            return self.budget.get_data(text=bud_type)
        return _package_output_df(self.model, self._package(), bud_type=bud_type)

    def plot_flows(
            self,
            kstpkper: tuple = None,
            cfd_to_cfs: bool = False,
            cfd_to_gpm: bool = False,
            html: str = None
    ):
        """Plot flow-by-reach values for one or more stress periods.

        Draws, per selected time, each reach's NEGATIVE ``FLOW-JA-FACE`` entries
        (its flow to the next reach) sign-flipped to positive against the
        reach's zero-based ``ifno``, and ``show()``s the figure. Returns
        ``None``. It raised ``KeyError: 'rno'`` on every model until ledger
        476 (FloPy 3.10 names the column ``ifno``); it no longer reads
        ``packagedata`` at all, since every node of the SFR budget is a reach.

        Parameters
        ----------
        kstpkper : tuple of (int, int) or list of tuple of (int, int), optional
            The zero-based ``(kstp, kper)`` times to draw, one trace each, on the
            SFR budget file's OWN clock (ledger 301): a time that file did not save
            raises ``ValueError`` naming the times it did. A bare tuple is one time.
            Default ``None``: the first time the SFR budget file saved (it was the
            HEAD file's first time until ledger 476).
        cfd_to_cfs : bool, default False
            Divide every flow by 86,400 (ft3/d to ft3/s in a feet-and-days model).
        cfd_to_gpm : bool, default False
            Draw gallons per minute: divide by 86,400 and multiply by 448.8, so it
            implies ``cfd_to_cfs`` and wins over it.
        html : str or pathlib.Path, optional
            Also write the figure to this HTML file (``Fig.write_html``). Default
            ``None``: nothing is written.
        """
        flows: pd.DataFrame = self.get('flow')
        saved = list(dict.fromkeys(flows['kstpkper']))
        if kstpkper is None:
            kstpkper = [saved[0]]
        elif isinstance(kstpkper, tuple) and len(kstpkper) == 2 and all(
            isinstance(part, (int, np.integer)) for part in kstpkper
        ):
            kstpkper = [kstpkper]
        times = [tuple(int(part) for part in stpper) for stpper in kstpkper]
        unsaved = [time for time in times if time not in saved]
        if unsaved:
            raise ValueError(
                f"kstpkper {unsaved} is not a time this SFR budget file saved; it saved "
                f"{saved[:6]}{' ...' if len(saved) > 6 else ''}."
            )
        gpm = 448.8 if cfd_to_gpm else 1
        cfs = (24 * 60 * 60) if cfd_to_cfs or cfd_to_gpm else 1
        from myflopy.viz import Fig

        fig = Fig()
        for time in times:
            # A reach's NEGATIVE FLOW-JA-FACE entries are its flow to the reaches
            # downstream. MF6 writes the reach number 1-based; FloPy's packagedata
            # ``ifno`` (the rno this read before FloPy 3.10 renamed it) is
            # zero-based, so the x axis is that ``ifno`` (ledger 476).
            at_time = flows[(flows['kstpkper'] == time) & (flows['q'] < 0)]
            fig.add_scattergl(
                x=at_time['node'] - 1,
                y=at_time['q'] * -1 / cfs * gpm,
                name=f'{time}',
            )
        fig.update_xaxes(title_text='reach (ifno, zero-based)')
        if html is not None:
            fig.write_html(file=html)
        fig.show()


class SFRStage:
    """Compatibility wrapper for stream stage output access.

    The canonical home for this behavior is now ``model.outputs.sfr.stage``.
    This class remains for backward compatibility with older call sites.
    """

    def __init__(self, model: SimulationBase):
        """Create one SFR stage accessor for ``model`` (its only SFR package)."""
        self.model = model
        self.nreaches = one_package(self.model, "sfr").nreaches.data

    def get(self):
        """Return ``model.outputs.sfr.stage.get()``: one row per reach, one column per
        saved record, ``NaN`` where a reach is dry or inactive (ledger 353)."""
        return self.model.outputs.sfr.stage.get()


class DRNBudget:
    """Convenience wrapper for DRN-specific budget visualization.

    This is a small compatibility wrapper around the older DRN choropleth
    plotting entry point.
    """

    def __init__(self, model: SimulationBase = None):
        """Create one DRN budget plotting helper for ``model``."""
        self.model = model


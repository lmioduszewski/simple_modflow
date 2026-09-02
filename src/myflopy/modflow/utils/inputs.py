from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase

import pandas as pd

from myflopy._logging import get_logger
from myflopy.modflow.mf6.package_plotting import (
    _apply_backend,
    refuse_noun_parameters,
    resolve_noun_hover,
)
from myflopy.modflow.utils.datatypes.choros import Choro

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

    def map(
            self,
            per: int = 0,
            multiplier: float = 12 * 30,
            *,
            # -- colour ----------------------------------------------------------
            zmin: float | None = None,
            zmax: float | None = None,
            colorscale: str | list | tuple | None = None,
            logscale: bool = False,
            # -- contours --------------------------------------------------------
            contours: bool | str = False,
            contour_values=None,
            contour_levels: int | float | list = 10,
            contour_color: str = "black",
            contour_width: float = 1.5,
            contour_name: str | None = None,
            contour_clip: bool = True,
            contour_resolution: int = 150,
            contour_method: str = "linear",
            # -- highlighting ----------------------------------------------------
            select=None,
            select_style: str = "outline",
            select_color: str | None = None,
            # -- overlays and framing --------------------------------------------
            locs=None,
            hillshade_path=None,
            fit_bounds: bool = True,
            bounds_padding: float = 0.05,
            # -- hover -----------------------------------------------------------
            hover=None,
            # -- renderer --------------------------------------------------------
            backend: str = "plotly",
            **trace_kwargs,
    ):
        """UZF infiltration (``finf``) for one period, as a choropleth.

        A RECORD noun, and a legacy one: ``model.inputs.uzf`` predates the
        ``model.packages.<pkg>.inputs.<field>`` grammar and survives because UZF
        is indexed by ``ifno`` rather than by cellid, so its period table has to
        be joined back to the package data before it can be coloured. That join
        is :meth:`finf`, and it is why this signature carries ``multiplier`` on
        top of the shared drawing parameters -- the free verb never converts
        units, because it is handed the values directly. Cells with no UZF object
        are filled with zero rather than left blank, so the picture is of the
        whole grid.

        Every parameter is named rather than swept into ``**kwargs``: PyCharm and
        Pylance read the ``def`` line and never run the module, so a parameter
        that arrives through a tail is one no editor can ever offer (plan 8.8).

        ``per`` and ``multiplier`` stay positional-or-keyword, because that is how
        this method has always been called (``uzf.map(3, 1.0)``); everything added
        below is keyword-only.

        The per-layer arguments the free verb has -- ``kstpkper``,
        ``per_timestep``, ``bgs``, ``hover_layers``, ``hover_surfaces``,
        ``show_layer_elevs``, ``show_mounding`` -- are deliberately absent. Each
        was measured to leave a record noun's figure byte-identical, and
        ``show_mounding`` is worse than inert: it injects a head-derived row into
        an elevation tooltip on every cell. ``values``/``type`` and the legacy
        hover trio raise instead of being forwarded: infiltration is what this
        noun draws, and overriding the numbers under a tooltip that still reports
        ``finf`` produces a figure whose cells and whose labels disagree.

        Parameters
        ----------
        per : int, default 0
            Stress period to read, zero-based -- both the UZF period whose
            ``finf`` is drawn and the period the underlying map reads for its
            tooltip. Unlike :meth:`finf`, ``"all"`` is not accepted here: a
            choropleth colours one number per cell, so pick a period (or build
            several maps and hand them to ``model.plot.animate``).
        multiplier : float, default 360
            Scale every infiltration rate before drawing. The default ``12 * 30``
            converts MODFLOW's ft/day into inches/month, which is the unit
            recharge is usually reviewed in; pass ``1.0`` to see the raw model
            units, or ``365.25`` for ft/year. Applied to every cell, including
            the zero-filled ones, so it never changes which cells are blank.

        Returns
        -------
        Choro or matplotlib.figure.Figure
            A :class:`~myflopy.viz.Picture` under the default backend --
            ``.show()``, ``.save(path)``, ``.html(path)`` -- or a bare
            Matplotlib figure with ``backend="mpl"``.

        See Also
        --------
        finf : the per-cell infiltration series behind the picture.

        Examples
        --------
        >>> model.inputs.uzf.map()                       # in/month, period 0
        >>> model.inputs.uzf.map(per=2)
        >>> model.inputs.uzf.map(3, 1.0)                 # period 3, raw ft/day
        >>> model.inputs.uzf.map(logscale=True, colorscale="Blues")
        >>> model.inputs.uzf.map(contours=True, contour_levels=6, zmax=4)
        >>> model.inputs.uzf.map(backend="mpl").savefig("finf.png")
        """

        refuse_noun_parameters("inputs.uzf", "finf", trace_kwargs)
        hover = resolve_noun_hover(hover, trace_kwargs)
        finf = self.finf(per, multiplier).finf
        choro = self.model.plot.map(
            per=per,
            custom_zs=finf.to_list(),
            zmin=zmin,
            zmax=zmax,
            colorscale=colorscale,
            logscale=logscale,
            contours=contours,
            contour_values=contour_values,
            contour_levels=contour_levels,
            contour_color=contour_color,
            contour_width=contour_width,
            contour_name=contour_name,
            contour_clip=contour_clip,
            contour_resolution=contour_resolution,
            contour_method=contour_method,
            select=select,
            select_style=select_style,
            select_color=select_color,
            locs=locs,
            hillshade_path=hillshade_path,
            fit_bounds=fit_bounds,
            bounds_padding=bounds_padding,
            hover=hover,
            **trace_kwargs,
        )
        return _apply_backend(choro, backend)


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

    def map(
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

    def map(
            self,
            per: tuple = None,
    ):
        """Highlight the DRN cells on the grid for one stress period."""

        model = self.model
        vor = model.vor

        drn_cells = pd.DataFrame(model.gwf.drn.stress_period_data.data[per])['cellid'].apply(lambda x: x[1]).to_list()
        return vor.plot.map(select=drn_cells)

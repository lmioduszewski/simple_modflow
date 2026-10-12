"""Shared plotting helpers for bound named-series targets."""

from __future__ import annotations


class BoundNamedSeriesTargetPlots:
    """Shared plotting facade for model-bound stage/flow target sets.

    ``model.targets.<family>.plot`` for ``lake_stage``, ``sfr_stage``,
    ``sfr_flow`` and ``drn_flow``. Every picture is drawn from the bound set's
    ``compare()`` on the bound model, so it reads what that reads (a period's
    value at its END, a flow in MF6's sign -- zero or negative, water leaving --
    a dry lake at its bottom); a ``baseline`` model adds the same comparison on
    a second run. Values are in the family's units: model length units for a
    stage, length cubed per time for a flow.
    """

    def __init__(
        self,
        bound_targets,
        *,
        target_column: str,
        simulated_column: str,
        title: str,
        yaxis_title: str,
    ):
        """Configure the shared plotting facade with the target/simulated column names and labels."""

        self.bound_targets = bound_targets
        self.target_column = target_column
        self.simulated_column = simulated_column
        self.title = title
        self.yaxis_title = yaxis_title

    def obs_vs_sim(self, *, baseline=None, backend: str = "plotly"):
        """Return an observed-vs-simulated cross plot (optionally with a ``baseline`` model).

        One point per comparison row holding both a target and a simulated
        value, with a 1:1 line spanning every value drawn, titled for the family
        ("Observed vs simulated lake stage", "... SFR stage", "... SFR flow",
        "... DRN seepage"). A simulated value that is MODFLOW's no-value marker
        (``|value| >= 1e29``: a dry or inactive SFR reach on the stage-file
        route) is not drawn (ledger 355); ``compare()`` keeps it. When no row of
        the bound model's comparison holds both values, its series is empty, the
        figure carries a note naming the missing side, and a ``UserWarning`` says
        so -- for a set given no values (every named-series set on the canonical
        model) that is all it draws: no point and no 1:1 line.

        Parameters
        ----------
        baseline : SimulationBase or None, optional
            A second run model: the SAME targets compared against it (the
            wrapped set's ``compare(baseline)``; DRN zones drawn as geometry are
            resolved against ITS grid) are drawn as a "Baseline" series beside the
            bound model's "Model" series. Default ``None`` draws the bound model
            alone. Its comparison raises as the bound model's does.
        backend : str, optional
            ``"plotly"`` (default; also ``"interactive"``) returns a
            :class:`~myflopy.modflow.calcs.calibration.CalibrationPlot`;
            ``"matplotlib"`` (also ``"mpl"``, ``"static"``) a static
            ``matplotlib.figure.Figure``. Case and surrounding space are
            ignored; anything else raises ``ValueError``, after the comparisons
            have run.

        Returns
        -------
        CalibrationPlot or matplotlib.figure.Figure
            The cross plot, per ``backend``.
        """

        from myflopy.modflow.calcs.calibration import CalibrationPlot

        current = self.bound_targets.compare()
        baseline_frame = None if baseline is None else self.bound_targets.targets.compare(baseline)
        return CalibrationPlot.from_obs_vs_sim(
            current,
            baseline_compare=baseline_frame,
            target_column=self.target_column,
            simulated_column=self.simulated_column,
            title=self.title,
            yaxis_title=self.yaxis_title,
            backend=backend,
        )

    def timeseries(self, name: str | None = None, *, baseline=None, backend: str = "plotly"):
        """Return an observed/simulated time-series plot for one target ``name``.

        The target's values ("Target") and the bound model's ("Model") against
        the target's ``time`` -- the zero-based stress period -- titled with the
        target's name. A simulated no-value marker (``|value| >= 1e29``) is a
        gap in the line (ledger 355).

        Parameters
        ----------
        name : str or None, optional
            The target (lake, reach or zone name) to draw, matched against the
            comparison's ``name`` column ignoring case (and the argument's
            surrounding space). Default ``None`` draws the set's only target,
            and raises ``ValueError`` when the set holds more than one -- it
            never draws them all; a name matching no row raises ``ValueError``
            too.
        baseline : SimulationBase or None, optional
            A second run model whose simulated values for the same target are
            drawn as a "Baseline" line (the wrapped set's ``compare(baseline)``).
            Default ``None`` draws the bound model alone.
        backend : str, optional
            ``"plotly"`` (default; also ``"interactive"``) or ``"matplotlib"`` (also
            ``"mpl"``, ``"static"``), as in :meth:`obs_vs_sim`;
            anything else raises ``ValueError``.

        Returns
        -------
        CalibrationPlot or matplotlib.figure.Figure
            The time series, per ``backend``.
        """

        from myflopy.modflow.calcs.calibration import CalibrationPlot

        current = self.bound_targets.compare()
        baseline_frame = None if baseline is None else self.bound_targets.targets.compare(baseline)
        return CalibrationPlot.from_timeseries(
            current,
            name=name,
            baseline_compare=baseline_frame,
            target_column=self.target_column,
            simulated_column=self.simulated_column,
            yaxis_title=self.yaxis_title,
            title=None,
            backend=backend,
        )

    def residuals_by_period(self, *, baseline=None, backend: str = "plotly"):
        """Return a by-period residual summary plot (optionally with a ``baseline`` model).

        The mean absolute residual of the measured rows in each ``time`` (the
        zero-based stress period), "Model MAE", titled "Residual MAE by period"
        for every family; a period with no measured row is a gap, and a residual
        against MODFLOW's no-value marker is left out of its mean (ledger 355).

        Parameters
        ----------
        baseline : SimulationBase or None, optional
            A second run model whose per-period MAE for the same targets is drawn
            as "Baseline MAE" (the wrapped set's ``compare(baseline)``). Default
            ``None`` draws the bound model alone.
        backend : str, optional
            ``"plotly"`` (default; also ``"interactive"``) or ``"matplotlib"`` (also
            ``"mpl"``, ``"static"``), as in :meth:`obs_vs_sim`;
            anything else raises ``ValueError``.

        Returns
        -------
        CalibrationPlot or matplotlib.figure.Figure
            The period summary, per ``backend``.
        """

        from myflopy.modflow.calcs.calibration import CalibrationPlot

        current = self.bound_targets.compare()
        baseline_frame = None if baseline is None else self.bound_targets.targets.compare(baseline)
        return CalibrationPlot.from_residuals_by_period(current, baseline_compare=baseline_frame, backend=backend)

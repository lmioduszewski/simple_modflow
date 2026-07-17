"""Shared plotting helpers for bound named-series targets."""

from __future__ import annotations


class BoundNamedSeriesTargetPlots:
    """Shared plotting facade for model-bound stage/flow target sets."""

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
        """Return an observed-vs-simulated cross plot (optionally with a ``baseline`` model)."""

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
        """Return an observed/simulated time-series plot for one target ``name`` (or all)."""

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
        """Return a by-period residual summary plot (optionally with a ``baseline`` model)."""

        from myflopy.modflow.calcs.calibration import CalibrationPlot

        current = self.bound_targets.compare()
        baseline_frame = None if baseline is None else self.bound_targets.targets.compare(baseline)
        return CalibrationPlot.from_residuals_by_period(current, baseline_compare=baseline_frame, backend=backend)



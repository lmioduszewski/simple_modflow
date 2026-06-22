"""Read and visualize a completed PESTPP-IES run.

:class:`IesResults` opens the output of an iterative ensemble smoother run
(launched via :meth:`PestProject.run_ies`, or any pestpp-ies run) and turns the
raw output files into the assessment plots and ensembles you actually want:
phi convergence, prior-vs-posterior ensembles against the observations, and
posterior forecast distributions.

The headline plots are interactive Plotly figures (myflopy house style); the
underlying ``pyemu.ObservationEnsemble`` / ``pyemu.ParameterEnsemble`` objects
remain available for anything bespoke.

PESTPP-IES output files this reads (``<case>`` is the control-file stem):

- ``<case>.phi.actual.csv`` / ``<case>.phi.meas.csv`` -- objective function by
  iteration (one column per realization, plus a ``base`` column);
- ``<case>.<iter>.obs.csv`` -- the observation ensemble at each iteration
  (iteration ``0`` is the prior; the highest iteration is the posterior);
- ``<case>.<iter>.par.csv`` -- the parameter ensemble at each iteration;
- ``<case>.obs+noise.csv`` -- the measurement-noise ensemble.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots


# Plotly colors used consistently across the IES plots.
_PRIOR_COLOR = "rgba(150,150,150,0.45)"
_POST_COLOR = "rgba(31,119,180,0.55)"
_NOISE_COLOR = "rgba(214,39,40,0.45)"
_MEAS_COLOR = "rgb(214,39,40)"
_TRUTH_COLOR = "rgb(214,39,40)"

# Matplotlib/seaborn equivalents (used when backend="matplotlib").
_MPL_PRIOR = "0.6"
_MPL_POST = "#1f77b4"
_MPL_MEAS = "crimson"
_MPL_TRUTH = "crimson"


def _normalize_backend(backend: str) -> str:
    """Map friendly backend names to ``"plotly"`` or ``"matplotlib"``."""

    value = str(backend).strip().lower()
    if value in ("plotly", "go"):
        return "plotly"
    if value in ("matplotlib", "mpl", "seaborn", "sns"):
        return "matplotlib"
    raise ValueError(f"backend must be 'plotly' or 'matplotlib', got {backend!r}.")


def _new_mpl_axes(**kwargs):
    """Create a seaborn-styled matplotlib figure/axes without global side effects."""

    import matplotlib.pyplot as plt
    import seaborn as sns

    with sns.axes_style("whitegrid"):
        fig, ax = plt.subplots(**kwargs)
    return fig, ax

# Trailing ``:<value>`` time token in a pyEMU long observation name.
_TIME_RE = re.compile(r":([0-9eE.+\-]+)$")


def _import_pyemu():
    try:
        import pyemu
    except ModuleNotFoundError as exc:  # pragma: no cover - environment guard
        raise ModuleNotFoundError(
            "pyemu is required to read PESTPP-IES results."
        ) from exc
    return pyemu


def _parse_time(obsnme: str) -> float | None:
    """Parse the trailing time index out of a pyEMU long observation name."""

    match = _TIME_RE.search(str(obsnme))
    if not match:
        return None
    try:
        return float(match.group(1))
    except ValueError:
        return None


@dataclass
class IesSettings:
    """Printable summary of how an IES run was configured and what it produced."""

    case: str
    workspace: str
    num_reals: int | None
    noptmax: int | None
    iterations_available: list[int]
    nnz_obs: int
    n_forecasts: int
    has_noise: bool
    pestpp_options: dict

    def __str__(self) -> str:
        ies_opts = {k: v for k, v in self.pestpp_options.items() if str(k).startswith("ies_")}
        lines = [
            f"PESTPP-IES run: {self.case}",
            f"  workspace   : {self.workspace}",
            f"  realizations: {self.num_reals}",
            f"  noptmax     : {self.noptmax}   iterations on disk: {self.iterations_available}",
            f"  observations: {self.nnz_obs} nonzero-weight   forecasts: {self.n_forecasts}",
            f"  noise ensemble: {'yes' if self.has_noise else 'no'}",
            f"  ies options : {ies_opts if ies_opts else '(defaults)'}",
        ]
        return "\n".join(lines)

    def __repr__(self) -> str:
        return self.__str__()


@dataclass
class IesForecast:
    """One forecast's prior and posterior ensemble, with summary and plot.

    Returned by :meth:`IesResults.forecast`. Wraps the prior and posterior
    samples of a single prediction of interest so its uncertainty can be
    summarized or plotted in one line.
    """

    name: str
    prior: np.ndarray
    posterior: np.ndarray
    truth: float | None = None

    def summary(self) -> pd.Series:
        """Return prior/posterior mean, std, and 5/50/95 percentiles."""

        def _stats(values, tag):
            values = np.asarray(values, dtype=float)
            return {
                f"{tag}_mean": float(np.mean(values)),
                f"{tag}_std": float(np.std(values)),
                f"{tag}_p05": float(np.percentile(values, 5)),
                f"{tag}_p50": float(np.percentile(values, 50)),
                f"{tag}_p95": float(np.percentile(values, 95)),
            }

        data = {"forecast": self.name}
        data.update(_stats(self.prior, "prior"))
        data.update(_stats(self.posterior, "posterior"))
        if self.truth is not None:
            data["truth"] = float(self.truth)
        # Fraction by which posterior standard deviation shrank vs the prior.
        prior_std = np.std(self.prior)
        if prior_std > 0:
            data["uncertainty_reduction"] = float(1.0 - np.std(self.posterior) / prior_std)
        return pd.Series(data)

    def plot(self, *, bins: int = 20, title: str | None = None, backend: str = "plotly"):
        """Overlay prior and posterior histograms with the truth/target line.

        Grey is the prior forecast distribution, blue the posterior. A vertical
        red line marks the target/known value when one is available.

        ``backend`` selects ``"plotly"`` (interactive, default) or
        ``"matplotlib"`` (static matplotlib/seaborn).
        """

        prior = np.asarray(self.prior, dtype=float)
        posterior = np.asarray(self.posterior, dtype=float)
        edges = np.histogram_bin_edges(np.concatenate([prior, posterior]), bins=bins)

        if _normalize_backend(backend) == "matplotlib":
            import seaborn as sns

            fig, ax = _new_mpl_axes(figsize=(6, 4))
            ax.hist(prior, bins=edges, density=True, color=_MPL_PRIOR, alpha=0.55, label="prior")
            ax.hist(posterior, bins=edges, density=True, color=_MPL_POST, alpha=0.6, label="posterior")
            if self.truth is not None:
                ax.axvline(float(self.truth), color=_MPL_TRUTH, linestyle="--")
            ax.set_xlabel("forecast value")
            ax.set_ylabel("probability density")
            ax.set_title(title or f"Forecast: {self.name}")
            ax.legend()
            sns.despine(fig)
            return fig

        fig = go.Figure()
        fig.add_histogram(x=prior, xbins=dict(start=edges[0], end=edges[-1],
                          size=(edges[-1] - edges[0]) / bins), name="prior",
                          marker_color=_PRIOR_COLOR, histnorm="probability density")
        fig.add_histogram(x=posterior, xbins=dict(start=edges[0], end=edges[-1],
                          size=(edges[-1] - edges[0]) / bins), name="posterior",
                          marker_color=_POST_COLOR, histnorm="probability density")
        fig.update_layout(barmode="overlay", title=title or f"Forecast: {self.name}",
                          xaxis_title="forecast value", yaxis_title="probability density",
                          template="plotly_white")
        if self.truth is not None:
            fig.add_vline(x=float(self.truth), line_color=_TRUTH_COLOR, line_dash="dash")
        return fig


class IesResults:
    """Open and assess a completed PESTPP-IES run.

    Parameters
    ----------
    workspace
        Directory holding the pestpp-ies output (the master directory for a
        parallel run, or the template workspace for a serial run).
    case_name
        Control-file stem (e.g. ``"calib"`` for ``calib.pst``). Auto-discovered
        from the single ``.pst`` in ``workspace`` when omitted.
    model
        Optional myflopy model (carrying the grid) used for spatial parameter
        maps (:meth:`plot_field`). Populated automatically when the run is
        launched via :meth:`PestProject.run_ies`.

    Notes
    -----
    Iteration ``0`` is the prior ensemble; the highest iteration written to disk
    is treated as the posterior (pestpp-ies may stop before ``NOPTMAX`` if it
    converges). Use :meth:`obs_ensemble` / :meth:`par_ensemble` to reach any
    specific iteration.
    """

    def __init__(self, workspace, *, case_name: str | None = None, model=None):
        self.workspace = Path(workspace)
        self.pyemu = _import_pyemu()
        self.case = case_name or self._discover_case()
        self.pst = self.pyemu.Pst(str(self.workspace / f"{self.case}.pst"))
        self.model = model
        try:
            self.pst.try_parse_name_metadata()
        except Exception:  # pragma: no cover - metadata parsing is best-effort
            pass

    @cached_property
    def _metadata(self) -> dict:
        """The myflopy PEST metadata written at build time, if present."""

        path = self.workspace / "myflopy_pest_metadata.json"
        if not path.exists():
            return {}
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except (ValueError, OSError):  # pragma: no cover - best effort
            return {}

    @property
    def capture_fields(self) -> list[dict]:
        """Captured parameter-field definitions available for spatial maps."""

        return list(self._metadata.get("capture_fields", []))

    # -- discovery --------------------------------------------------------

    def _discover_case(self) -> str:
        pst_files = sorted(self.workspace.glob("*.pst"))
        if not pst_files:
            raise FileNotFoundError(f"No .pst control file found in {self.workspace}.")
        # Prefer one with companion phi output (a finished ies run).
        for candidate in pst_files:
            if (self.workspace / f"{candidate.stem}.phi.actual.csv").exists():
                return candidate.stem
        return pst_files[0].stem

    @cached_property
    def iterations(self) -> list[int]:
        """Sorted iteration numbers with an observation ensemble on disk."""

        nums = []
        for path in self.workspace.glob(f"{self.case}.*.obs.csv"):
            middle = path.name[len(self.case) + 1 : -len(".obs.csv")]
            if middle.isdigit():
                nums.append(int(middle))
        return sorted(nums)

    @property
    def prior_iteration(self) -> int:
        """Iteration treated as the prior (the lowest available, normally 0)."""

        if not self.iterations:
            raise FileNotFoundError("No observation-ensemble files were found for this run.")
        return self.iterations[0]

    @property
    def posterior_iteration(self) -> int:
        """Iteration treated as the posterior (the highest available)."""

        if not self.iterations:
            raise FileNotFoundError("No observation-ensemble files were found for this run.")
        return self.iterations[-1]

    # -- ensembles --------------------------------------------------------

    def obs_ensemble(self, iteration: int):
        """Return the observation ensemble at ``iteration`` as a pyEMU object."""

        path = self.workspace / f"{self.case}.{int(iteration)}.obs.csv"
        if not path.exists():
            raise FileNotFoundError(f"No observation ensemble at iteration {iteration}: {path}")
        return self.pyemu.ObservationEnsemble.from_csv(pst=self.pst, filename=str(path))

    def par_ensemble(self, iteration: int):
        """Return the parameter ensemble at ``iteration`` as a pyEMU object."""

        path = self.workspace / f"{self.case}.{int(iteration)}.par.csv"
        if not path.exists():
            raise FileNotFoundError(f"No parameter ensemble at iteration {iteration}: {path}")
        return self.pyemu.ParameterEnsemble.from_csv(pst=self.pst, filename=str(path))

    @property
    def prior(self):
        """The prior observation ensemble (iteration 0)."""

        return self.obs_ensemble(self.prior_iteration)

    @property
    def posterior(self):
        """The posterior observation ensemble (highest available iteration)."""

        return self.obs_ensemble(self.posterior_iteration)

    @property
    def prior_parameters(self):
        """The prior parameter ensemble."""

        return self.par_ensemble(self.prior_iteration)

    @property
    def posterior_parameters(self):
        """The posterior parameter ensemble."""

        return self.par_ensemble(self.posterior_iteration)

    @cached_property
    def noise(self):
        """The observation+noise ensemble, or ``None`` if the run used no noise."""

        path = self.workspace / f"{self.case}.obs+noise.csv"
        if not path.exists():
            return None
        return self.pyemu.ObservationEnsemble.from_csv(pst=self.pst, filename=str(path))

    @cached_property
    def phi(self) -> pd.DataFrame:
        """The 'actual' objective-function progression (one column per realization)."""

        return pd.read_csv(self.workspace / f"{self.case}.phi.actual.csv")

    @cached_property
    def phi_measured(self) -> pd.DataFrame:
        """The 'measured+noise' objective-function progression."""

        return pd.read_csv(self.workspace / f"{self.case}.phi.meas.csv")

    @property
    def forecast_names(self) -> list[str]:
        """Observation names registered as forecasts (predictions of interest)."""

        return list(self.pst.forecast_names)

    # -- plots ------------------------------------------------------------

    def plot_phi(self, *, log: bool = True, measured: bool = False, backend: str = "plotly"):
        """Plot objective-function (phi) convergence across iterations.

        One faint line per realization shows how its misfit dropped each
        iteration; the bold line is the ensemble mean. This is the first plot to
        look at -- it tells you whether history matching actually improved the
        fit and whether the ensemble collapsed.

        Parameters
        ----------
        log
            Use a log-scaled phi axis (default ``True``; recommended).
        measured
            Plot the 'measured+noise' phi instead of the 'actual' phi.
        backend
            ``"plotly"`` (interactive, default) or ``"matplotlib"`` (static
            matplotlib/seaborn).
        """

        frame = self.phi_measured if measured else self.phi
        realization_cols = list(frame.columns[6:])
        x = frame["iteration"].to_numpy()
        title = "Phi convergence" + (" (measured+noise)" if measured else "")

        if _normalize_backend(backend) == "matplotlib":
            import seaborn as sns

            fig, ax = _new_mpl_axes(figsize=(7, 4))
            for col in realization_cols:
                ax.plot(x, frame[col], color="0.5", lw=0.8, alpha=0.4)
            ax.plot(x, frame["mean"], color=_MPL_POST, lw=2.5, label="mean phi")
            if log:
                ax.set_yscale("log")
            ax.set_xlabel("iteration")
            ax.set_ylabel("phi")
            ax.set_title(title)
            ax.legend()
            sns.despine(fig)
            return fig

        fig = go.Figure()
        for col in realization_cols:
            fig.add_scatter(x=x, y=frame[col], mode="lines", line=dict(color="rgba(80,80,80,0.35)", width=1),
                            name=str(col), showlegend=False, hoverinfo="skip")
        fig.add_scatter(x=x, y=frame["mean"], mode="lines+markers",
                        line=dict(color=_POST_COLOR.replace("0.55", "1.0"), width=3), name="mean phi")
        fig.update_layout(title=title, xaxis_title="iteration", yaxis_title="phi",
                          template="plotly_white")
        if log:
            fig.update_yaxes(type="log")
        return fig

    def _observation_metadata(self) -> pd.DataFrame:
        """Return observation rows annotated with parsed time and group."""

        obs = self.pst.observation_data.copy()
        obs["_time"] = [_parse_time(name) for name in obs.index]
        return obs

    def plot_vs_obs(self, *, groups: list[str] | None = None, max_groups: int = 12,
                    iteration: int | None = None, backend: str = "plotly"):
        """Plot prior and posterior ensembles against the measured observations.

        For every nonzero-weight observation group (one timeseries each), grey
        lines are the prior ensemble, blue lines the posterior, and red markers
        the measured values (with the noise ensemble shaded when available).
        A good fit means the blue lines bracket the red markers.

        Parameters
        ----------
        groups
            Observation group names to plot. Defaults to all nonzero-weight
            groups (capped at ``max_groups``).
        max_groups
            Maximum number of groups to draw when ``groups`` is not given.
        iteration
            Posterior iteration to plot (default: the highest available).
        backend
            ``"plotly"`` (interactive, default) or ``"matplotlib"`` (static
            matplotlib/seaborn).
        """

        obs = self._observation_metadata()
        nnz_groups = list(self.pst.nnz_obs_groups)
        chosen = groups if groups is not None else nnz_groups[:max_groups]
        if not chosen:
            raise ValueError("No nonzero-weight observation groups to plot.")

        prior = self.prior._df
        posterior = self.obs_ensemble(iteration if iteration is not None else self.posterior_iteration)._df

        if _normalize_backend(backend) == "matplotlib":
            import matplotlib.pyplot as plt
            import seaborn as sns

            with sns.axes_style("whitegrid"):
                fig, axes = plt.subplots(len(chosen), 1, figsize=(8, 2.6 * len(chosen)), squeeze=False)
            for ax, group in zip(axes[:, 0], chosen):
                group_obs = obs.loc[obs["obgnme"] == group].sort_values("_time")
                names = group_obs.index.tolist()
                times = group_obs["_time"].to_numpy()
                for real in prior.index:
                    ax.plot(times, prior.loc[real, names].to_numpy(dtype=float), color="0.6", lw=0.8, alpha=0.4)
                for real in posterior.index:
                    ax.plot(times, posterior.loc[real, names].to_numpy(dtype=float), color=_MPL_POST, lw=0.8, alpha=0.5)
                ax.plot(times, group_obs["obsval"].to_numpy(dtype=float), "^", color=_MPL_MEAS, ms=7)
                ax.set_title(group, loc="left", fontsize=9)
                ax.set_ylabel("value")
            axes[-1, 0].set_xlabel("time")
            fig.suptitle("Simulated ensemble vs measured observations")
            fig.tight_layout()
            return fig

        fig = make_subplots(rows=len(chosen), cols=1, subplot_titles=chosen, shared_xaxes=False)
        for row, group in enumerate(chosen, start=1):
            group_obs = obs.loc[obs["obgnme"] == group].copy()
            group_obs = group_obs.sort_values("_time")
            names = group_obs.index.tolist()
            times = group_obs["_time"].to_numpy()
            first = True
            for real in prior.index:
                fig.add_scatter(x=times, y=prior.loc[real, names].to_numpy(dtype=float), mode="lines",
                                line=dict(color=_PRIOR_COLOR, width=1), row=row, col=1,
                                name="prior", legendgroup="prior", showlegend=first and row == 1, hoverinfo="skip")
                first = False
            first = True
            for real in posterior.index:
                fig.add_scatter(x=times, y=posterior.loc[real, names].to_numpy(dtype=float), mode="lines",
                                line=dict(color=_POST_COLOR, width=1), row=row, col=1,
                                name="posterior", legendgroup="posterior", showlegend=first and row == 1, hoverinfo="skip")
                first = False
            fig.add_scatter(x=times, y=group_obs["obsval"].to_numpy(dtype=float), mode="markers",
                            marker=dict(color=_MEAS_COLOR, size=8, symbol="triangle-up"), row=row, col=1,
                            name="measured", legendgroup="measured", showlegend=row == 1)
        fig.update_layout(title="Simulated ensemble vs measured observations",
                          template="plotly_white", height=260 * len(chosen))
        return fig

    def forecast(self, name) -> IesForecast:
        """Return one forecast's prior+posterior ensemble for summary/plotting.

        Parameters
        ----------
        name
            A forecast observation name (see :attr:`forecast_names`), or a
            :class:`HeadTargets`-style object / spec that was registered with
            :meth:`PestProject.forecast` (its single observation is resolved).

        Returns
        -------
        IesForecast
            With ``.summary()`` and ``.plot()``.
        """

        obsnme = self._resolve_forecast_name(name)
        prior = self.prior._df.loc[:, obsnme].to_numpy(dtype=float)
        posterior = self.posterior._df.loc[:, obsnme].to_numpy(dtype=float)
        truth = self.pst.observation_data.loc[obsnme, "obsval"]
        truth = float(truth) if pd.notna(truth) else None
        return IesForecast(name=obsnme, prior=prior, posterior=posterior, truth=truth)

    def _resolve_forecast_name(self, name) -> str:
        names = self.forecast_names
        if isinstance(name, str):
            if name in names:
                return name
            matches = [n for n in names if name.lower() in n.lower()]
            if len(matches) == 1:
                return matches[0]
            if not matches:
                raise KeyError(f"No forecast matching {name!r}. Available: {names}")
            raise KeyError(f"Ambiguous forecast {name!r}; matches {matches}.")
        raise TypeError(
            "forecast() expects a forecast observation name (a string from "
            "`forecast_names`)."
        )

    def forecasts(self) -> pd.DataFrame:
        """Return a prior/posterior uncertainty summary for every forecast."""

        rows = [self.forecast(name).summary() for name in self.forecast_names]
        if not rows:
            return pd.DataFrame()
        return pd.DataFrame(rows).set_index("forecast")

    # -- phi & weight diagnostics ----------------------------------------

    def plot_phi_distribution(self, *, measured: bool = False, bins: int = 25,
                              backend: str = "plotly"):
        """Histogram of objective-function (phi) values: prior vs posterior.

        Purpose
        -------
        Shows the *spread* of misfit across the ensemble before and after history
        matching (grey = prior / iteration 0, blue = posterior / final
        iteration), on a log scale. This is the "did we learn, and did we learn
        too much?" plot.

        What to look for
        ----------------
        - A clear shift to **lower** phi from prior to posterior means the
          ensemble assimilated information from the data.
        - A posterior that **collapses to a narrow spike at very low phi** is a
          warning sign of over-fitting (especially if forecast distributions also
          collapse) -- with an imperfect model you should not expect to drive phi
          to zero. As a rule of thumb the achievable phi is around the number of
          nonzero-weight observations (``self.pst.nnz_obs``).

        Parameters
        ----------
        measured
            Use the 'measured+noise' phi (each realization vs its noisy data copy)
            instead of the 'actual' phi (vs the raw observed values).
        backend
            ``"plotly"`` (interactive, default) or ``"matplotlib"``.
        """

        frame = self.phi_measured if measured else self.phi
        cols = list(frame.columns[6:])
        prior = np.log10(np.asarray(frame.iloc[0][cols], dtype=float))
        posterior = np.log10(np.asarray(frame.iloc[-1][cols], dtype=float))
        edges = np.histogram_bin_edges(np.concatenate([prior, posterior]), bins=bins)
        title = "Phi distribution: prior vs posterior" + (" (measured+noise)" if measured else "")

        if _normalize_backend(backend) == "matplotlib":
            import seaborn as sns

            fig, ax = _new_mpl_axes(figsize=(6.5, 4))
            ax.hist(prior, bins=edges, color=_MPL_PRIOR, alpha=0.55, label="prior")
            ax.hist(posterior, bins=edges, color=_MPL_POST, alpha=0.6, label="posterior")
            ax.set_xlabel(r"$\log_{10}\phi$")
            ax.set_ylabel("realizations")
            ax.set_title(title)
            ax.legend()
            sns.despine(fig)
            return fig

        size = (edges[-1] - edges[0]) / bins
        fig = go.Figure()
        fig.add_histogram(x=prior, xbins=dict(start=edges[0], end=edges[-1], size=size),
                          name="prior", marker_color=_PRIOR_COLOR)
        fig.add_histogram(x=posterior, xbins=dict(start=edges[0], end=edges[-1], size=size),
                          name="posterior", marker_color=_POST_COLOR)
        fig.update_layout(barmode="overlay", title=title, template="plotly_white",
                          xaxis_title="log10(phi)", yaxis_title="realizations")
        return fig

    def parameters_at_bounds(self, *, which: str = "posterior", tol: float = 0.01) -> pd.DataFrame:
        """Return, by parameter group, how many parameters are pinned at a bound.

        A parameter counts as "at a bound" when its ensemble-mean value sits
        within ``tol`` (a fraction of the bound range, measured in the
        parameter's transform space) of its lower or upper bound.

        Returns
        -------
        pandas.DataFrame
            Indexed by parameter group with columns ``n_parameters``,
            ``n_at_lower``, ``n_at_upper`` and ``pct_at_bound``.
        """

        ensemble = (self.posterior_parameters if which == "posterior" else self.prior_parameters)._df
        par = self.pst.parameter_data.loc[[c for c in ensemble.columns if c in self.pst.parameter_data.index]]
        means = ensemble[par.index].mean()
        lower = par["parlbnd"].astype(float)
        upper = par["parubnd"].astype(float)
        is_log = par["partrans"].astype(str).str.lower().eq("log")

        value, low, high = means.copy(), lower.copy(), upper.copy()
        log_mask = is_log & (lower > 0) & (upper > 0) & (means > 0)
        value[log_mask] = np.log10(means[log_mask])
        low[log_mask] = np.log10(lower[log_mask])
        high[log_mask] = np.log10(upper[log_mask])
        span = (high - low).replace(0, np.nan)
        position = (value - low) / span
        at_lower = position <= tol
        at_upper = position >= (1.0 - tol)

        summary = pd.DataFrame({
            "pargp": par["pargp"].to_numpy(),
            "at_lower": at_lower.to_numpy(),
            "at_upper": at_upper.to_numpy(),
        })
        grouped = summary.groupby("pargp").agg(
            n_parameters=("pargp", "size"),
            n_at_lower=("at_lower", "sum"),
            n_at_upper=("at_upper", "sum"),
        )
        grouped["pct_at_bound"] = 100.0 * (grouped["n_at_lower"] + grouped["n_at_upper"]) / grouped["n_parameters"]
        return grouped.sort_values("pct_at_bound", ascending=False)

    def plot_parameters_at_bounds(self, *, which: str = "posterior", tol: float = 0.01,
                                  backend: str = "plotly"):
        """Bar chart of the percentage of parameters at their bounds, by group.

        Purpose
        -------
        Parameters pinned at their bounds are trying to move further than you
        allowed them to. A quick health check on the prior.

        What to look for
        ----------------
        A large fraction of a group at its bounds usually means the **prior is
        too tight** (widen the bounds) or those parameters are **compensating for
        structural error** elsewhere in the model. A few at bounds is normal; a
        whole group at bounds deserves investigation.

        Parameters
        ----------
        which
            ``"posterior"`` (default) or ``"prior"``.
        tol
            Closeness to a bound (fraction of the bound range) that counts as
            "at" the bound.
        backend
            ``"plotly"`` (interactive, default) or ``"matplotlib"``.
        """

        table = self.parameters_at_bounds(which=which, tol=tol)
        groups = table.index.tolist()
        pct = table["pct_at_bound"].to_numpy(dtype=float)
        title = f"Parameters at bounds ({which})"

        if _normalize_backend(backend) == "matplotlib":
            import seaborn as sns

            fig, ax = _new_mpl_axes(figsize=(6.5, max(2.5, 0.5 * len(groups) + 1)))
            ax.barh(groups, pct, color=_MPL_POST)
            ax.set_xlabel("% of parameters at a bound")
            ax.set_title(title)
            ax.invert_yaxis()
            sns.despine(fig)
            return fig

        fig = go.Figure(go.Bar(x=pct, y=groups, orientation="h", marker_color=_POST_COLOR))
        fig.update_layout(title=title, template="plotly_white",
                          xaxis_title="% of parameters at a bound", yaxis_title="parameter group")
        return fig

    def phi_contributions(self, *, realization: str = "base") -> pd.Series:
        """Return the phi (misfit) contributed by each observation group.

        Computed for one realization (default the ``base`` / minimum-error-
        variance one) as the sum of squared weighted residuals within each
        nonzero-weight observation group.

        Returns
        -------
        pandas.Series
            Indexed by observation group, sorted from largest contribution down.
        """

        obs = self.pst.observation_data
        obs = obs.loc[obs["weight"].astype(float) > 0]
        ensemble = self.posterior._df
        if realization in ensemble.index:
            simulated = ensemble.loc[realization]
        else:
            simulated = ensemble.mean()
        names = [name for name in obs.index if name in simulated.index]
        obs = obs.loc[names]
        residual = obs["weight"].astype(float) * (simulated[names].astype(float) - obs["obsval"].astype(float))
        contribution = pd.Series((residual.to_numpy() ** 2), index=names)
        by_group = contribution.groupby(obs["obgnme"].to_numpy()).sum()
        return by_group.sort_values(ascending=False)

    def plot_phi_contributions(self, *, realization: str = "base", kind: str = "bar",
                               max_groups: int = 20, backend: str = "plotly"):
        """Plot the misfit (phi) contributed by each observation group.

        Purpose
        -------
        Shows *which observations the objective function is actually made of* --
        the basis for "visibility weighting". Groups that dominate phi dominate
        the calibration.

        What to look for
        ----------------
        - One or two groups **dominating** total phi means those observations
          (often high-magnitude or densely sampled) are steering the fit; consider
          re-weighting so groups contribute more evenly toward your model purpose.
        - A group with stubbornly high contribution is something the model
          **cannot fit** -- worth understanding before you trust forecasts that
          depend on it.

        Parameters
        ----------
        realization
            Realization to evaluate (default the ``base`` realization).
        kind
            ``"bar"`` (default) or ``"pie"``.
        max_groups
            Cap the number of groups shown (largest contributors first).
        backend
            ``"plotly"`` (interactive, default) or ``"matplotlib"``.
        """

        contributions = self.phi_contributions(realization=realization).head(max_groups)
        labels = [str(name) for name in contributions.index]
        values = contributions.to_numpy(dtype=float)
        title = "Phi contribution by observation group"

        if _normalize_backend(backend) == "matplotlib":
            import seaborn as sns

            if str(kind).lower() == "pie":
                fig, ax = _new_mpl_axes(figsize=(5.5, 5.5))
                ax.pie(values, labels=labels, autopct="%1.0f%%", textprops={"fontsize": 8})
                ax.set_title(title)
                return fig
            fig, ax = _new_mpl_axes(figsize=(6.5, max(2.5, 0.4 * len(labels) + 1)))
            ax.barh(labels, values, color=_MPL_POST)
            ax.set_xlabel("phi contribution")
            ax.set_title(title)
            ax.invert_yaxis()
            sns.despine(fig)
            return fig

        if str(kind).lower() == "pie":
            fig = go.Figure(go.Pie(labels=labels, values=values))
            fig.update_layout(title=title, template="plotly_white")
            return fig
        fig = go.Figure(go.Bar(x=values, y=labels, orientation="h", marker_color=_POST_COLOR))
        fig.update_layout(title=title, template="plotly_white",
                          xaxis_title="phi contribution", yaxis_title="observation group")
        return fig

    # -- prior Monte Carlo & prior-data conflict -------------------------

    def conflict(self, *, coverage: float = 1.0) -> pd.DataFrame:
        """Detect prior-data conflict: observations the prior ensemble can't reach.

        For every nonzero-weight observation, compares the measured value against
        the range the *prior* simulated ensemble can produce. An observation is
        "in conflict" when the measured value falls outside that range -- meaning
        no plausible parameter set reproduces it, so history matching cannot
        either (the problem is in the model, bounds, or weights, not the
        parameters).

        Parameters
        ----------
        coverage
            Central fraction of the prior ensemble treated as "reachable".
            ``1.0`` (default) uses the full min-max range; ``0.95`` uses the
            2.5-97.5 percentile band (more robust to a few extreme realizations).

        Returns
        -------
        pandas.DataFrame
            Indexed by observation name with columns ``obgnme``, ``time``,
            ``measured``, ``prior_lo``, ``prior_hi``, ``prior_mean`` and
            ``in_conflict``.
        """

        obs = self.pst.observation_data
        obs = obs.loc[obs["weight"].astype(float) > 0]
        prior = self.prior._df
        names = [name for name in obs.index if name in prior.columns]
        obs = obs.loc[names]
        sims = prior[names].astype(float)
        alpha = (1.0 - float(coverage)) / 2.0
        lo = sims.quantile(alpha)
        hi = sims.quantile(1.0 - alpha)
        measured = obs["obsval"].astype(float)
        meta = self._observation_metadata().loc[names]
        return pd.DataFrame(
            {
                "obgnme": obs["obgnme"].to_numpy(),
                "time": meta["_time"].to_numpy(),
                "measured": measured.to_numpy(),
                "prior_lo": lo.to_numpy(),
                "prior_hi": hi.to_numpy(),
                "prior_mean": sims.mean().to_numpy(),
                "in_conflict": ((measured < lo) | (measured > hi)).to_numpy(),
            },
            index=names,
        )

    def plot_prior_vs_obs(self, *, groups: list[str] | None = None, max_groups: int = 12,
                          show_conflict: bool = True, coverage: float = 1.0, backend: str = "plotly"):
        """Plot the prior ensemble (grey spaghetti) against the measured data.

        Purpose
        -------
        The prior Monte Carlo check, done *before* history matching: each grey
        line is one prior realization's simulated equivalent through time; red
        markers are the measured observations. It answers "does the prior even
        contain the data?" -- the cheap go/no-go gate before calibrating.

        What to look for
        ----------------
        - The grey envelope should **bracket** the red markers. Then there is a
          plausible parameter set that reproduces the data, and history matching
          has something to find.
        - Any red marker **outside** the grey envelope is **prior-data conflict**
          (highlighted in orange when ``show_conflict``): no plausible parameter
          set reaches it. Fix the model / bounds / weights before history
          matching -- don't try to calibrate your way out of it.

        Parameters
        ----------
        groups
            Observation groups to plot (default: all nonzero-weight groups,
            capped at ``max_groups``).
        max_groups
            Maximum number of groups when ``groups`` is not given.
        show_conflict
            Highlight observations in prior-data conflict (default ``True``).
        coverage
            Passed to :meth:`conflict` for the conflict test.
        backend
            ``"plotly"`` (interactive, default) or ``"matplotlib"``.
        """

        obs = self._observation_metadata()
        nnz_groups = list(self.pst.nnz_obs_groups)
        chosen = groups if groups is not None else nnz_groups[:max_groups]
        if not chosen:
            raise ValueError("No nonzero-weight observation groups to plot.")
        prior = self.prior._df
        conflict = self.conflict(coverage=coverage) if show_conflict else None

        def _group_data(group):
            group_obs = obs.loc[obs["obgnme"] == group].sort_values("_time")
            names = group_obs.index.tolist()
            times = group_obs["_time"].to_numpy()
            measured = group_obs["obsval"].to_numpy(dtype=float)
            if conflict is not None:
                flags = conflict["in_conflict"].reindex(names).fillna(False).to_numpy(dtype=bool)
            else:
                flags = np.zeros(len(names), dtype=bool)
            return names, times, measured, flags

        if _normalize_backend(backend) == "matplotlib":
            import matplotlib.pyplot as plt
            import seaborn as sns

            with sns.axes_style("whitegrid"):
                fig, axes = plt.subplots(len(chosen), 1, figsize=(8, 2.6 * len(chosen)), squeeze=False)
            for ax, group in zip(axes[:, 0], chosen):
                names, times, measured, flags = _group_data(group)
                for real in prior.index:
                    ax.plot(times, prior.loc[real, names].to_numpy(dtype=float), color="0.6", lw=0.8, alpha=0.4)
                ax.plot(times[~flags], measured[~flags], "^", color=_MPL_MEAS, ms=7, label="measured")
                if flags.any():
                    ax.plot(times[flags], measured[flags], "X", color="darkorange", ms=10, label="prior-data conflict")
                ax.set_title(group, loc="left", fontsize=9)
                ax.set_ylabel("value")
            axes[-1, 0].set_xlabel("time")
            fig.suptitle("Prior ensemble vs measured observations")
            fig.tight_layout()
            return fig

        fig = make_subplots(rows=len(chosen), cols=1, subplot_titles=chosen, shared_xaxes=False)
        for row, group in enumerate(chosen, start=1):
            names, times, measured, flags = _group_data(group)
            first = True
            for real in prior.index:
                fig.add_scatter(x=times, y=prior.loc[real, names].to_numpy(dtype=float), mode="lines",
                                line=dict(color=_PRIOR_COLOR, width=1), row=row, col=1,
                                name="prior", legendgroup="prior", showlegend=first and row == 1, hoverinfo="skip")
                first = False
            fig.add_scatter(x=times[~flags], y=measured[~flags], mode="markers",
                            marker=dict(color=_MEAS_COLOR, size=8, symbol="triangle-up"),
                            row=row, col=1, name="measured", legendgroup="measured", showlegend=row == 1)
            if flags.any():
                fig.add_scatter(x=times[flags], y=measured[flags], mode="markers",
                                marker=dict(color="darkorange", size=11, symbol="x"),
                                row=row, col=1, name="prior-data conflict", legendgroup="conflict", showlegend=row == 1)
        fig.update_layout(title="Prior ensemble vs measured observations",
                          template="plotly_white", height=260 * len(chosen))
        return fig

    def plot_conflict(self, *, coverage: float = 1.0, backend: str = "plotly"):
        """Bar chart of the percentage of observations in prior-data conflict, by group.

        Purpose
        -------
        Summarizes :meth:`conflict` so you can see *where* the prior fails to
        bracket the data.

        What to look for
        ----------------
        Groups with a high conflict percentage are where the model (or its prior /
        weights) most needs revisiting before history matching. Zero across the
        board is the green light to proceed.

        Parameters
        ----------
        coverage
            Passed to :meth:`conflict`.
        backend
            ``"plotly"`` (interactive, default) or ``"matplotlib"``.
        """

        table = self.conflict(coverage=coverage)
        grouped = table.groupby("obgnme")["in_conflict"].agg(["sum", "size"])
        grouped["pct"] = 100.0 * grouped["sum"] / grouped["size"]
        grouped = grouped.sort_values("pct", ascending=False)
        labels = [str(name) for name in grouped.index]
        pct = grouped["pct"].to_numpy(dtype=float)
        title = "Prior-data conflict by observation group"

        if _normalize_backend(backend) == "matplotlib":
            import seaborn as sns

            fig, ax = _new_mpl_axes(figsize=(6.5, max(2.5, 0.4 * len(labels) + 1)))
            ax.barh(labels, pct, color="darkorange")
            ax.set_xlabel("% of observations in prior-data conflict")
            ax.set_title(title)
            ax.invert_yaxis()
            sns.despine(fig)
            return fig

        fig = go.Figure(go.Bar(x=pct, y=labels, orientation="h", marker_color="darkorange"))
        fig.update_layout(title=title, template="plotly_white",
                          xaxis_title="% in prior-data conflict", yaxis_title="observation group")
        return fig

    # -- spatial parameter maps ------------------------------------------

    _ARR_RE = re.compile(r"arr_i:(\d+)_j:(\d+)")

    def _capture_info(self, target: str) -> dict:
        """Resolve a captured-field definition by target or parameter name."""

        fields = self.capture_fields
        if not fields:
            raise ValueError(
                "No captured parameter fields are available. Re-run the setup with "
                "cal.parameterize(..., capture=True) to record per-cell fields."
            )
        key = str(target).strip().lower()
        for info in fields:
            if key in (str(info.get("target", "")).lower(), str(info.get("name", "")).lower()):
                return info
        available = [f.get("target") for f in fields]
        raise KeyError(f"No captured field for {target!r}. Available: {available}")

    def _field_cell_index(self, prefix: str, layer: int) -> dict[str, int]:
        """Map captured-field observation names to cell ids for one layer."""

        if self.model is None:
            raise ValueError(
                "Spatial maps need the model grid. Open the run via "
                "cal.run_ies(...) (which passes the model), or pass model=... to IesResults."
            )
        ncpl = int(self.model.vor.ncpl)
        token = f"oname:{prefix.lower()}_"
        mapping: dict[str, int] = {}
        for name in self.pst.observation_data.index:
            if token not in name:
                continue
            match = self._ARR_RE.search(name)
            if not match:
                continue
            flat = int(match.group(1))
            if flat // ncpl == int(layer):
                mapping[name] = flat % ncpl
        if not mapping:
            raise ValueError(f"No captured field cells found for layer {layer}.")
        return mapping

    def field(self, target: str, *, layer: int = 0) -> pd.DataFrame:
        """Return per-cell prior/posterior statistics for a captured field.

        Requires the parameter to have been declared with
        ``cal.parameterize(..., capture=True)``.

        Returns
        -------
        pandas.DataFrame
            Indexed by ``cell`` with columns ``prior_mean``, ``prior_std``,
            ``posterior_mean``, ``posterior_std``, ``change``
            (``posterior_mean / prior_mean``), and ``base`` (the
            minimum-error-variance realization) when available.
        """

        info = self._capture_info(target)
        if info.get("family") != "array":
            raise NotImplementedError(
                f"Spatial maps currently support array fields (K, K33); "
                f"{info.get('target')!r} is a list field."
            )
        cells = self._field_cell_index(info["prefix"], layer)
        names = list(cells.keys())
        prior = self.prior._df.loc[:, names]
        posterior = self.posterior._df.loc[:, names]
        frame = pd.DataFrame({
            "cell": [cells[name] for name in names],
            "prior_mean": prior.mean().to_numpy(),
            "prior_std": prior.std().to_numpy(),
            "posterior_mean": posterior.mean().to_numpy(),
            "posterior_std": posterior.std().to_numpy(),
        })
        if "base" in posterior.index:
            frame["base"] = posterior.loc["base"].to_numpy()
        frame["change"] = frame["posterior_mean"] / frame["prior_mean"]
        return frame.sort_values("cell").reset_index(drop=True)

    def plot_field(self, target: str, *, stat: str = "mean", which: str = "posterior",
                   layer: int = 0, backend: str = "plotly", **choropleth_kwargs):
        """Map a captured parameter field on the model grid (Voronoi choropleth).

        This answers "property patterns -- plausible or laughable?": it shows the
        spatial pattern of a calibrated property and how history matching changed
        it. Requires ``cal.parameterize(..., capture=True)``.

        Parameters
        ----------
        target
            A captured target name, e.g. ``"k"``.
        stat
            ``"mean"`` (ensemble mean field), ``"std"`` (posterior spread --
            where the property is still uncertain), ``"base"`` (the
            minimum-error-variance realization), or ``"change"``
            (posterior_mean / prior_mean -- where calibration moved the property).
        which
            ``"prior"`` or ``"posterior"`` for ``stat`` in {``"mean"``, ``"std"``}.
        layer
            Model layer to map (default 0).
        backend
            ``"plotly"`` (interactive map, default) or ``"matplotlib"`` (static
            matplotlib choropleth of the Voronoi cells).
        **choropleth_kwargs
            Forwarded to the Plotly Voronoi choropleth builder (``backend="plotly"``).
        """

        frame = self.field(target, layer=layer)
        column = {
            "mean": f"{which}_mean",
            "std": f"{which}_std",
            "base": "base",
            "change": "change",
        }.get(str(stat).lower())
        if column is None or column not in frame.columns:
            raise ValueError(
                f"Unknown stat {stat!r} (or unavailable). Use 'mean', 'std', 'base', or 'change'."
            )

        ncpl = int(self.model.vor.ncpl)
        values = np.full(ncpl, np.nan)
        values[frame["cell"].to_numpy(dtype=int)] = frame[column].to_numpy(dtype=float)
        label = f"{target} {stat}" + (f" ({which})" if stat in ("mean", "std") else "")

        if _normalize_backend(backend) == "matplotlib":
            import matplotlib.pyplot as plt

            gdf = self.model.vor.gdf_vorPolys.copy()
            gdf["value"] = values
            fig, ax = plt.subplots(figsize=(6, 6))
            cmap = "RdBu_r" if str(stat).lower() == "change" else "viridis"
            gdf.plot(column="value", ax=ax, legend=True, cmap=cmap)
            ax.set_title(label)
            ax.set_axis_off()
            return fig

        from myflopy.modflow.mf6.grid.plotting import build_choropleth

        choro = build_choropleth(
            self.model.vor, custom_zs=list(values), layer=layer, **choropleth_kwargs
        )
        return choro.plot()

    def best(self, *, criterion: str = "base") -> str:
        """Return the label of the single 'best' realization to carry forward.

        Parameters
        ----------
        criterion
            ``"base"`` (default, recommended) returns the ``base`` realization --
            the minimum-error-variance parameter set, the closest analogue to a
            deterministic GLM calibration. ``"min_phi"`` returns the lowest-phi
            realization; this is tempting but discouraged (it tends to be
            over-fit and noisy), so prefer ``"base"`` unless you have a reason.
        """

        criterion = str(criterion).lower()
        phi = self.phi
        realization_cols = list(phi.columns[6:])
        if criterion == "base":
            if "base" in realization_cols:
                return "base"
            criterion = "min_phi"  # fall through when no base realization exists
        if criterion == "min_phi":
            final = phi.iloc[-1]
            candidates = {col: float(final[col]) for col in realization_cols}
            return min(candidates, key=candidates.get)
        raise ValueError(f"Unknown criterion {criterion!r}; use 'base' or 'min_phi'.")

    @property
    def settings(self) -> IesSettings:
        """A printable summary of how this IES run was configured."""

        options = dict(self.pst.pestpp_options)
        num_reals = options.get("ies_num_reals")
        return IesSettings(
            case=self.case,
            workspace=str(self.workspace),
            num_reals=int(num_reals) if num_reals is not None else None,
            noptmax=int(self.pst.control_data.noptmax),
            iterations_available=self.iterations,
            nnz_obs=int(self.pst.nnz_obs),
            n_forecasts=len(self.forecast_names),
            has_noise=self.noise is not None,
            pestpp_options=options,
        )

    def report(self, html_path, *, max_groups: int = 6) -> Path:
        """Write a single self-contained HTML report of the headline IES plots.

        Bundles phi convergence, the ensemble-vs-observation comparison, and one
        histogram per forecast into one file -- the ensemble analogue of the
        deterministic ``PestRunResults.export_review``.

        Returns
        -------
        pathlib.Path
            The written HTML path.
        """

        html_path = Path(html_path)
        figures = [self.plot_phi(), self.plot_phi_distribution()]
        try:
            figures.append(self.plot_vs_obs(max_groups=max_groups))
        except ValueError:
            pass
        for diagnostic in (self.plot_phi_contributions, self.plot_parameters_at_bounds):
            try:
                figures.append(diagnostic())
            except Exception:  # pragma: no cover - diagnostics are best-effort in the bundle
                pass
        for name in self.forecast_names:
            figures.append(self.forecast(name).plot())
        if self.model is not None:
            for info in self.capture_fields:
                if info.get("family") != "array":
                    continue
                try:
                    figures.append(self.plot_field(info["target"], stat="mean", which="posterior"))
                    figures.append(self.plot_field(info["target"], stat="change"))
                except Exception:  # pragma: no cover - field maps are best-effort in the bundle
                    pass

        parts = ["<html><head><meta charset='utf-8'><title>IES review: "
                 f"{self.case}</title></head><body>", f"<h1>PESTPP-IES review: {self.case}</h1>",
                 f"<pre>{self.settings}</pre>"]
        for index, fig in enumerate(figures):
            parts.append(fig.to_html(full_html=False, include_plotlyjs="cdn" if index == 0 else False))
        parts.append("</body></html>")
        html_path.write_text("\n".join(parts), encoding="utf-8")
        return html_path


def open_ies_run(workspace, *, case_name: str | None = None, model=None) -> IesResults:
    """Open a completed PESTPP-IES run directory for assessment.

    Convenience wrapper around :class:`IesResults` mirroring ``open_pest_run``.
    Pass ``model=`` (a myflopy model carrying the grid) to enable spatial
    parameter maps via :meth:`IesResults.plot_field`.
    """

    return IesResults(workspace, case_name=case_name, model=model)

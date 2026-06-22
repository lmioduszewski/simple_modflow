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

    def plot(self, *, bins: int = 20, title: str | None = None) -> go.Figure:
        """Overlay prior and posterior histograms with the truth/target line.

        Grey is the prior forecast distribution, blue the posterior. A vertical
        red line marks the target/known value when one is available.
        """

        prior = np.asarray(self.prior, dtype=float)
        posterior = np.asarray(self.posterior, dtype=float)
        edges = np.histogram_bin_edges(np.concatenate([prior, posterior]), bins=bins)
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

    def plot_phi(self, *, log: bool = True, measured: bool = False) -> go.Figure:
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
        """

        frame = self.phi_measured if measured else self.phi
        realization_cols = list(frame.columns[6:])
        x = frame["iteration"].to_numpy()
        fig = go.Figure()
        for col in realization_cols:
            fig.add_scatter(x=x, y=frame[col], mode="lines", line=dict(color="rgba(80,80,80,0.35)", width=1),
                            name=str(col), showlegend=False, hoverinfo="skip")
        fig.add_scatter(x=x, y=frame["mean"], mode="lines+markers",
                        line=dict(color=_POST_COLOR.replace("0.55", "1.0"), width=3), name="mean phi")
        fig.update_layout(title="Phi convergence" + (" (measured+noise)" if measured else ""),
                          xaxis_title="iteration", yaxis_title="phi", template="plotly_white")
        if log:
            fig.update_yaxes(type="log")
        return fig

    def _observation_metadata(self) -> pd.DataFrame:
        """Return observation rows annotated with parsed time and group."""

        obs = self.pst.observation_data.copy()
        obs["_time"] = [_parse_time(name) for name in obs.index]
        return obs

    def plot_vs_obs(self, *, groups: list[str] | None = None, max_groups: int = 12,
                    iteration: int | None = None) -> go.Figure:
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
        """

        obs = self._observation_metadata()
        nnz_groups = list(self.pst.nnz_obs_groups)
        chosen = groups if groups is not None else nnz_groups[:max_groups]
        if not chosen:
            raise ValueError("No nonzero-weight observation groups to plot.")

        prior = self.prior._df
        posterior = self.obs_ensemble(iteration if iteration is not None else self.posterior_iteration)._df
        noise = self.noise._df if self.noise is not None else None

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
                   layer: int = 0, **choropleth_kwargs) -> go.Figure:
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
        **choropleth_kwargs
            Forwarded to the Voronoi choropleth builder.
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
        figures = [self.plot_phi()]
        try:
            figures.append(self.plot_vs_obs(max_groups=max_groups))
        except ValueError:
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

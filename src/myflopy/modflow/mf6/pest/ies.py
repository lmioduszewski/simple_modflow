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
import warnings
from collections.abc import Sequence
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from myflopy import viz
from myflopy._logging import get_logger
from myflopy._optional import require
from myflopy.modflow.mf6.grid.plotting import build_choropleth
from myflopy.modflow.mf6.package_plotting import (
    _red_white_blue_diverging_colorscale,
    _symmetric_color_limit,
)
from myflopy.modflow.utils.datatypes.choros import mpl_colormap_for
from myflopy.modflow.utils.datatypes.hover import (
    parameter_field_hover,
    residual_hover,
)

# Plot colors come from the shared palette (one findable place: myflopy.viz).
_PRIOR_COLOR = viz.PALETTE.prior
_POST_COLOR = viz.PALETTE.posterior
_POST_SOLID = viz.PALETTE.posterior_solid
_ENSEMBLE_COLOR = viz.PALETTE.ensemble
_MEAS_COLOR = viz.PALETTE.measured
_TRUTH_COLOR = viz.PALETTE.truth
_CONFLICT_COLOR = viz.PALETTE.conflict
_MPL_PRIOR = viz.PALETTE.mpl_prior
_MPL_POST = viz.PALETTE.mpl_posterior
_MPL_MEAS = viz.PALETTE.mpl_measured
_MPL_TRUTH = viz.PALETTE.mpl_truth
_MPL_ENSEMBLE = viz.PALETTE.mpl_ensemble
_MPL_CONFLICT = viz.PALETTE.mpl_conflict


logger = get_logger(__name__)


def _normalize_backend(backend: str) -> str:
    """Map friendly backend names to ``"plotly"`` or ``"matplotlib"``."""

    value = str(backend).strip().lower()
    if value in ("plotly", "go"):
        return "plotly"
    if value in ("matplotlib", "mpl", "seaborn", "sns"):
        return "matplotlib"
    raise ValueError(f"backend must be 'plotly' or 'matplotlib', got {backend!r}.")


def _log_decade_colorbar(low: float, high: float, *, unit: str = "") -> dict:
    """Relabel a log10 colorbar in real units, with ticks at whole decades.

    Nothing in ``Choro`` sets ``tickvals``/``ticktext``, so a ``logscale`` map is
    otherwise read in log10 units -- a colorbar running -1 to 1 for a field that
    actually runs 0.1 to 10. This rides ``Choro``'s trace kwargs. It does NOT
    survive :func:`myflopy.viz.mosaic`, whose shared coloraxis carries no
    colorbar key (compromise ledger 71).
    """

    decades = [
        decade for decade in range(int(np.floor(low)), int(np.ceil(high)) + 1)
        if low <= decade <= high
    ]
    # Fewer than three whole decades in range gives a colorbar with one or two
    # ticks, which reads as broken; fall back to five evenly spaced labels.
    ticks = [float(d) for d in decades] if len(decades) >= 3 else [
        float(t) for t in np.linspace(low, high, 5)
    ]
    return {
        "tickvals": ticks,
        "ticktext": [f"{10.0 ** tick:.3g}{unit}" for tick in ticks],
    }


def _log_decade_colorbar_for_mosaic(low, high) -> dict:
    """``viz.mosaic`` colorbar callback: decades over the POOLED panel limits.

    A mosaic pools its panels onto one shared color axis, so the real-unit ticks
    have to be computed from the limits *it* resolved, not from any one panel's
    (compromise ledger 71). ``viz.mosaic`` leaves the limits unset when no panel
    has finite data, hence the guard.
    """

    if low is None or high is None:
        return {}
    return _log_decade_colorbar(float(low), float(high))


# Display names for stats whose column name reads badly as a figure label.
_STAT_LABELS = {"reduction": "uncertainty reduction"}

# Which divisor was zero when a derived stat comes back non-finite.
_UNDEFINED_REASONS = {
    "change": "undefined (prior mean 0)",
    "reduction": "undefined (prior sd 0)",
}


def _field_map_policy(stat: str, values) -> dict[str, object]:
    """House color policy for one captured-parameter-field map, as ``Choro`` kwargs.

    Registry-free by design: ``package_registry`` is keyed by package/field
    NAME, and a *statistic* of a captured array is not a package field. This
    mirrors ``prt_maps.PRT_COLORSCALE`` -- a derived map whose scale is pinned at
    its single source. Pinned by ``tests/test_colorscale_policy.py``.

    Returned as ``Choro`` kwargs rather than a bare colorscale name because for a
    ratio the scale and its LIMITS are one indivisible decision: a diverging
    scale whose ``zmin``/``zmax`` are not symmetric silently puts white somewhere
    other than 1.
    """

    key = str(stat).lower()
    if key == "change":
        # posterior_mean / prior_mean is a RATIO whose neutral value is 1, not 0,
        # so it is mapped as log10(ratio): a halving and a doubling then sit the
        # same distance either side of the middle. errstate because a zero prior
        # mean makes this warn -- Choro's own _logscaled is already guarded; this
        # second log10 is ours.
        with np.errstate(divide="ignore", invalid="ignore"):
            logged = np.log10(np.asarray(values, dtype=float))
        # _symmetric_color_limit drops non-finite values and returns 0.0 for an
        # empty or all-NaN input, so an inf cannot widen the range. `or 1.0` gives
        # a defined +/- one decade for the degenerate cases: a prior Monte-Carlo
        # run (prior == posterior, so change is identically 1.0) and an all-NaN
        # layer.
        half = _symmetric_color_limit(logged) or 1.0
        return {
            "colorscale": _red_white_blue_diverging_colorscale(),
            "logscale": True,
            # Symmetric zmin/zmax in LOG space is what actually centers this map:
            # plot_mpl reads self._zmin/_zmax and never consults the trace kwargs,
            # so zmid alone would center the plotly map and leave matplotlib
            # autoscaled. zmid is belt, not braces.
            "zmin": -half,
            "zmax": half,
            "zmid": 0.0,
            "colorbar": _log_decade_colorbar(-half, half, unit="x"),
        }
    if key in ("mean", "base"):
        # The only array-family targets field() can reach are k/k33 --
        # conductivities, which span decades.
        with np.errstate(divide="ignore", invalid="ignore"):
            logged = np.log10(np.asarray(values, dtype=float))
        finite = logged[np.isfinite(logged)]
        colorbar = (
            _log_decade_colorbar(float(finite.min()), float(finite.max()))
            if finite.size else None
        )
        return {"colorscale": "earth", "logscale": True, "colorbar": colorbar}
    if key == "std":
        # A spread in model units is legitimately ZERO where the ensemble
        # collapsed, and log10(0) blanks the cell -- which would erase exactly
        # the cells this map exists to show. Stay linear.
        return {"colorscale": "earth", "logscale": False}
    if key == "reduction":
        policy = {"colorscale": "earth", "logscale": False}
        finite = np.asarray(values, dtype=float)
        finite = finite[np.isfinite(finite)]
        # Unlike a spread in model units, variance reduction has an ABSOLUTE
        # frame: 0 = the data told you nothing, 1 = the ensemble collapsed.
        # Anchor to it so two layers (or two runs) are comparable, and so a
        # field that only ever reduces 0.95-1.0 does not autoscale into a
        # dramatic-looking map of a trivial range. A posterior spread CAN grow
        # (negative reduction); when it does, fall back to the data's own range
        # rather than clipping those cells to the bottom color, where they
        # would read as "no reduction" instead of "worse than the prior".
        if finite.size and float(finite.min()) >= 0.0 and float(finite.max()) <= 1.0:
            policy.update({"zmin": 0.0, "zmax": 1.0})
        return policy
    raise ValueError(
        f"No field-map color policy for stat {stat!r}; "
        "expected 'mean', 'std', 'base', 'change', or 'reduction'."
    )


def _relabel_log_colorbar(figure, colorbar: dict | None) -> None:
    """Put a log map's real-unit ticks on the matplotlib colorbar.

    ``Choro.plot_mpl`` builds its colorbar through geopandas and exposes no tick
    hook, so a log-scaled static map would label its colorbar ``-3 … 2`` while
    the interactive one reads ``0.001 … 100`` — the same field, two readings, and
    the static one silently off by orders of magnitude. Applied only when the
    effective ``logscale`` is on, so an explicit ``logscale=False`` override does
    not get decade labels over linear values.
    """

    tickvals = (colorbar or {}).get("tickvals")
    ticktext = (colorbar or {}).get("ticktext")
    # axes[-1] is the colorbar geopandas appends; with colorbar=False there is
    # only the map axes and nothing to relabel.
    if not tickvals or not ticktext or len(figure.axes) < 2:
        return
    axis = figure.axes[-1]
    axis.set_yticks(list(tickvals))
    axis.set_yticklabels([str(text) for text in ticktext])


def _parse_cell_list(value) -> list[int]:
    """Parse a snapshotted ``"0,1,2"`` cell list, tolerating an empty zone.

    `_resolve_drn_zone_cells` assigns ``[]`` to a zone that intersects nothing,
    which the CSV snapshot round-trips as NaN (``str(nan) == "nan"``) while the
    GeoPackage snapshot round-trips as ``""``. Both mean the same thing.
    """

    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    return [
        int(cell) for cell in str(value).split(",")
        if cell.strip() and cell.strip().lower() != "nan"
    ]


def _hover_column(values, *, reason: str = "undefined (prior mean 0)") -> list:
    """Raw per-cell hover values: NaN -> blank, +/-inf -> a named reason.

    ``change = posterior_mean / prior_mean`` is an unguarded divide, so a zero
    prior mean gives +/-inf and 0/0 gives NaN. Both draw as an identical gap on
    the map, and ``format_number`` renders both as ``""``, which would make an
    undefined ratio indistinguishable from a cell that was never captured. A
    string payload entry passes through ``format_number`` verbatim. ``reason``
    names the divisor that was zero -- ``reduction`` divides by the prior sd,
    not the prior mean.
    """

    return [
        None if np.isnan(value) else (
            reason if np.isinf(value) else float(value)
        )
        for value in np.asarray(values, dtype=float)
    ]


# Trailing ``:<value>`` time token in a pyEMU long observation name.
_TIME_RE = re.compile(r":([0-9eE.+\-]+)$")


def _import_pyemu():
    """Import pyEMU (a required optional dependency for reading PESTPP-IES results)."""

    return require("pyemu", feature="reading PESTPP-IES results")


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
        """A multi-line human-readable summary of the run's configuration and outputs."""

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
        """Same as :meth:`__str__` (the settings summary)."""

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
            """Mean/std/5-50-95 percentiles of an ensemble, keyed by ``<tag>_<stat>``."""

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

            fig, ax = viz.mpl_axes(figsize=(6, 4))
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

        fig = viz.Fig()
        fig.add_histogram(x=prior, xbins=dict(start=edges[0], end=edges[-1],
                          size=(edges[-1] - edges[0]) / bins), name="prior",
                          marker_color=_PRIOR_COLOR, histnorm="probability density")
        fig.add_histogram(x=posterior, xbins=dict(start=edges[0], end=edges[-1],
                          size=(edges[-1] - edges[0]) / bins), name="posterior",
                          marker_color=_POST_COLOR, histnorm="probability density")
        fig.update_layout(barmode="overlay", title=title or f"Forecast: {self.name}",
                          xaxis_title="forecast value", yaxis_title="probability density",
                          dragmode="pan")
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
        """Open a completed PESTPP-IES run in ``workspace``, loading its ``.pst`` control file.

        Discovers the case name if not given, imports pyEMU, and best-effort
        parses observation-name metadata.
        """

        self.workspace = Path(workspace)
        self.pyemu = _import_pyemu()
        self.case = case_name or self._discover_case()
        self.pst = self.pyemu.Pst(str(self.workspace / f"{self.case}.pst"))
        self.model = model
        # No `try_parse_name_metadata()` call here: `pyemu.Pst(filename)` loads
        # via `Pst.load`, which calls it already (verified against the installed
        # pyemu). The explicit re-call was redundant, and its `except Exception`
        # was guarding work that had already happened.

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

    @property
    def observation_sets(self) -> list[dict]:
        """Observation-set definitions recorded at build time.

        One entry per ``cal.observe(...)`` **and** ``cal.forecast(...)`` call —
        a forecast is registered as an ordinary observation set carrying weight
        0, and nothing in the persisted metadata distinguishes the two. Each
        entry has its ``kind`` (``head_targets``, ``drn_flow``, ``lake_stage``,
        ``sfr_stage``, ``sfr_flow``), observation-name ``prefix``, and the
        snapshot files written beside the control file. This is what lets a
        completed run be re-joined to where its observations are on the grid.
        """

        return list(self._metadata.get("observation_sets", []))

    # -- discovery --------------------------------------------------------

    def _discover_case(self) -> str:
        """The case name from the workspace's ``.pst`` files, preferring one with phi output."""

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
        """Observation names registered as forecasts (predictions of interest).

        Empty when the run registered none -- which is the ordinary case, not an
        exotic one, since ``cal.forecast(...)`` is optional and
        ``PestProject._apply_forecasts`` writes no ``forecasts`` option without
        it. pyEMU returns ``None`` there (and ``[""]`` for an empty option
        string), so both are normalized away here rather than at the four call
        sites downstream.
        """

        return [name for name in (self.pst.forecast_names or []) if name]

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

            fig, ax = viz.mpl_axes(figsize=(7, 4))
            for col in realization_cols:
                ax.plot(x, frame[col], color=_MPL_ENSEMBLE, lw=0.8, alpha=0.4)
            ax.plot(x, frame["mean"], color=_MPL_POST, lw=2.5, label="mean phi")
            if log:
                ax.set_yscale("log")
            ax.set_xlabel("iteration")
            ax.set_ylabel("phi")
            ax.set_title(title)
            ax.legend()
            sns.despine(fig)
            return fig

        fig = viz.Fig()
        for col in realization_cols:
            fig.add_scatter(x=x, y=frame[col], mode="lines", line=dict(color=_ENSEMBLE_COLOR, width=1),
                            name=str(col), showlegend=False, hoverinfo="skip")
        fig.add_scatter(x=x, y=frame["mean"], mode="lines+markers",
                        line=dict(color=_POST_SOLID, width=3), name="mean phi")
        fig.update_layout(title=title, xaxis_title="iteration", yaxis_title="phi",
                          dragmode="pan")
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
            fig, axes = viz.mpl_axes(len(chosen), 1, figsize=(8, 2.6 * len(chosen)), squeeze=False)
            for ax, group in zip(axes[:, 0], chosen, strict=False):
                group_obs = obs.loc[obs["obgnme"] == group].sort_values("_time")
                names = group_obs.index.tolist()
                times = group_obs["_time"].to_numpy()
                for real in prior.index:
                    ax.plot(times, prior.loc[real, names].to_numpy(dtype=float), color=_MPL_PRIOR, lw=0.8, alpha=0.4)
                for real in posterior.index:
                    ax.plot(times, posterior.loc[real, names].to_numpy(dtype=float), color=_MPL_POST, lw=0.8, alpha=0.5)
                ax.plot(times, group_obs["obsval"].to_numpy(dtype=float), "^", color=_MPL_MEAS, ms=7)
                ax.set_title(group, loc="left", fontsize=9)
                ax.set_ylabel("value")
            axes[-1, 0].set_xlabel("time")
            fig.suptitle("Simulated ensemble vs measured observations")
            fig.tight_layout()
            return fig

        fig = viz.subplots(rows=len(chosen), cols=1, subplot_titles=chosen, shared_xaxes=False)
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
                          dragmode="pan", height=260 * len(chosen))
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
        """Resolve a forecast selector (exact name or index) to a full observation name."""

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
        # Realizations that failed during an iteration carry a non-finite phi
        # (a normal IES occurrence); drop them so the histogram range is finite.
        prior = prior[np.isfinite(prior)]
        posterior = posterior[np.isfinite(posterior)]
        combined = np.concatenate([prior, posterior])
        if combined.size == 0:
            raise ValueError("No finite phi values available to plot the distribution.")
        edges = np.histogram_bin_edges(combined, bins=bins)
        title = "Phi distribution: prior vs posterior" + (" (measured+noise)" if measured else "")

        if _normalize_backend(backend) == "matplotlib":
            import seaborn as sns

            fig, ax = viz.mpl_axes(figsize=(6.5, 4))
            ax.hist(prior, bins=edges, color=_MPL_PRIOR, alpha=0.55, label="prior")
            ax.hist(posterior, bins=edges, color=_MPL_POST, alpha=0.6, label="posterior")
            ax.set_xlabel(r"$\log_{10}\phi$")
            ax.set_ylabel("realizations")
            ax.set_title(title)
            ax.legend()
            sns.despine(fig)
            return fig

        size = (edges[-1] - edges[0]) / bins
        fig = viz.Fig()
        fig.add_histogram(x=prior, xbins=dict(start=edges[0], end=edges[-1], size=size),
                          name="prior", marker_color=_PRIOR_COLOR)
        fig.add_histogram(x=posterior, xbins=dict(start=edges[0], end=edges[-1], size=size),
                          name="posterior", marker_color=_POST_COLOR)
        fig.update_layout(barmode="overlay", title=title, dragmode="pan",
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

    def sensitivity(self, *, forecast: str | None = None,
                    by_group: bool = True) -> pd.DataFrame:
        """Which parameters the data informed, and which drive a forecast.

        **Ensemble-based, so it needs no jacobian** — everything here is computed
        from the prior and posterior ensembles this run already wrote. That is
        also its limitation, and the limitation is the point:

        * ``learned`` is ``1 - posterior_sd / prior_sd`` per parameter, in the
          parameter's transform space. It answers *"did the data constrain
          this?"* — the per-parameter analogue of
          ``plot_field(stat="reduction")``.
        * ``forecast_corr`` (when ``forecast=`` is given) is the empirical
          correlation across realizations between the parameter and that
          forecast. It answers *"does this parameter drive that prediction?"*

        **This is not CSS, and the two can disagree.** Composite scaled
        sensitivity is a local derivative at one parameter set; these are global
        measures conditioned on the prior. A parameter can be highly sensitive
        locally and score near zero here because the prior never moved it far,
        or vice versa. Where a jacobian-based answer is wanted, that needs a
        PESTPP-GLM run, which myflopy does not currently launch (§5.8).

        **Sampling noise is real.** With ``n`` realizations a correlation of
        roughly ``1/sqrt(n)`` is indistinguishable from zero — at the common
        ``reals=50`` that is about 0.14, so small values should not be ranked
        against each other. ``n_reals`` is returned so the reader can judge.

        Parameters
        ----------
        forecast
            Forecast name to correlate parameters against. ``None`` (default)
            returns the ``learned`` column only.
        by_group
            Aggregate to parameter groups (default). ``False`` returns one row
            per parameter, which for a grid or pilot-point run is large.

        Returns
        -------
        pandas.DataFrame
            Sorted by the strongest signal available, most informative first.
        """

        prior = self.prior_parameters._df
        posterior = self.posterior_parameters._df
        shared = [name for name in posterior.columns if name in prior.columns]
        par = self.pst.parameter_data
        adjustable = [
            name for name in shared
            if name in par.index
            and str(par.loc[name, "partrans"]).lower() not in ("fixed", "tied")
        ]
        if not adjustable:
            return pd.DataFrame(
                columns=["n_parameters", "learned", "forecast_corr", "n_reals"]
            )

        prior_values = prior[adjustable].astype(float)
        posterior_values = posterior[adjustable].astype(float)
        # Compare spreads in the space the parameter is estimated in, or a log
        # multiplier's ratio spread would be read as a linear one.
        is_log = par.loc[adjustable, "partrans"].astype(str).str.lower().eq("log")
        for name in [n for n in adjustable if is_log[n]]:
            if (prior_values[name] > 0).all() and (posterior_values[name] > 0).all():
                prior_values[name] = np.log10(prior_values[name])
                posterior_values[name] = np.log10(posterior_values[name])

        prior_sd = prior_values.std()
        learned = 1.0 - (posterior_values.std() / prior_sd.replace(0.0, np.nan))

        table = pd.DataFrame({
            "pargp": par.loc[adjustable, "pargp"].astype(str).to_numpy(),
            "learned": learned.reindex(adjustable).to_numpy(dtype=float),
        }, index=adjustable)
        table["n_reals"] = int(len(posterior_values))

        if forecast is not None:
            values = self._forecast_realizations(forecast, posterior_values.index)
            # Correlate in the SAME space the spreads were measured in, so a log
            # parameter's relationship with the forecast is not distorted.
            table["forecast_corr"] = [
                float(posterior_values[name].corr(values)) for name in adjustable
            ]

        if not by_group:
            return table.sort_values(
                "forecast_corr" if forecast is not None else "learned",
                key=abs if forecast is not None else None,
                ascending=False,
            )

        aggregations = {"n_parameters": ("learned", "size"),
                        "learned": ("learned", "mean"),
                        "n_reals": ("n_reals", "max")}
        if forecast is not None:
            # Mean of the ABSOLUTE correlation: a group whose members push a
            # forecast in opposite directions is still influential, and averaging
            # signed values would cancel it to zero.
            table["abs_corr"] = table["forecast_corr"].abs()
            aggregations["forecast_corr"] = ("abs_corr", "mean")
        grouped = table.groupby("pargp").agg(**aggregations)
        return grouped.sort_values(
            "forecast_corr" if forecast is not None else "learned", ascending=False
        )

    def _forecast_realizations(self, forecast: str, index) -> pd.Series:
        """One forecast's posterior realizations, aligned to a parameter index."""

        ensemble = self.posterior._df
        if forecast not in ensemble.columns:
            available = sorted(self.forecast_names)
            raise ValueError(
                f"No forecast named {forecast!r} in this run. Available: {available}."
            )
        values = ensemble[forecast].astype(float)
        # Parameter and observation ensembles share realization labels, but a
        # realization that failed appears in one and not the other; aligning on
        # the index rather than on position keeps them paired.
        return values.reindex(index)

    def plot_sensitivity(self, *, forecast: str | None = None, top: int = 15,
                         backend: str = "plotly"):
        """Bar chart of what the data informed, or of what drives a forecast.

        Purpose
        -------
        Two questions with one figure. Without ``forecast``, it shows how much
        each parameter group's uncertainty the data removed — groups near zero
        were not constrained by the observations you have. With ``forecast``, it
        shows how strongly each group co-varies with that prediction across the
        ensemble, which is the "what would I need to measure better?" view.

        What to look for
        ----------------
        A group with **high forecast correlation and low learning** is the
        uncomfortable one: it matters for the prediction and the data did not
        pin it down. That combination is where forecast uncertainty comes from.

        Read small bars with care — see :meth:`sensitivity` on sampling noise
        (roughly ``1/sqrt(n_reals)``), and on why this is not CSS.

        Parameters
        ----------
        forecast
            Forecast name; ``None`` (default) plots uncertainty reduction.
        top
            Keep the strongest ``top`` groups.
        backend
            ``"plotly"`` (interactive, default) or ``"matplotlib"``.
        """

        table = self.sensitivity(forecast=forecast)
        if table.empty:
            raise ValueError(
                "This run has no adjustable parameters in both ensembles, so "
                "there is nothing to score."
            )
        column = "forecast_corr" if forecast is not None else "learned"
        table = table.head(int(top))
        groups = table.index.tolist()
        values = table[column].to_numpy(dtype=float)
        noise = 1.0 / np.sqrt(max(int(table["n_reals"].max()), 1))
        label = (
            f"|correlation| with {forecast}" if forecast is not None
            else "uncertainty reduction (1 - post sd / prior sd)"
        )
        title = (
            f"Parameter influence on {forecast}" if forecast is not None
            else "What the data informed"
        )

        if _normalize_backend(backend) == "matplotlib":
            import seaborn as sns

            fig, ax = viz.mpl_axes(figsize=(6.5, max(2.5, 0.5 * len(groups) + 1)))
            ax.barh(groups, values, color=_MPL_POST)
            ax.axvline(noise, color="grey", linestyle="--", linewidth=1)
            ax.set_xlabel(label)
            ax.set_title(title)
            ax.invert_yaxis()
            sns.despine(fig)
            return fig

        fig = viz.Fig().add_trace(
            go.Bar(x=values, y=groups, orientation="h", marker_color=_POST_COLOR)
        )
        # The noise floor drawn on the figure, so a bar below it is visibly
        # indistinguishable from zero rather than merely documented as such.
        fig.add_vline(
            x=noise, line_dash="dash", line_color="grey",
            annotation_text=f"noise floor ~1/sqrt({int(table['n_reals'].max())})",
        )
        fig.update_layout(title=title, dragmode="pan",
                          xaxis_title=label, yaxis_title="parameter group")
        return fig

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

            fig, ax = viz.mpl_axes(figsize=(6.5, max(2.5, 0.5 * len(groups) + 1)))
            ax.barh(groups, pct, color=_MPL_POST)
            ax.set_xlabel("% of parameters at a bound")
            ax.set_title(title)
            ax.invert_yaxis()
            sns.despine(fig)
            return fig

        fig = viz.Fig().add_trace(go.Bar(x=pct, y=groups, orientation="h", marker_color=_POST_COLOR))
        fig.update_layout(title=title, dragmode="pan",
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
                fig, ax = viz.mpl_axes(figsize=(5.5, 5.5))
                ax.pie(values, labels=labels, autopct="%1.0f%%", textprops={"fontsize": 8})
                ax.set_title(title)
                return fig
            fig, ax = viz.mpl_axes(figsize=(6.5, max(2.5, 0.4 * len(labels) + 1)))
            ax.barh(labels, values, color=_MPL_POST)
            ax.set_xlabel("phi contribution")
            ax.set_title(title)
            ax.invert_yaxis()
            sns.despine(fig)
            return fig

        if str(kind).lower() == "pie":
            fig = viz.Fig().add_trace(go.Pie(labels=labels, values=values))
            fig.update_layout(title=title, dragmode="pan")
            return fig
        fig = viz.Fig().add_trace(go.Bar(x=values, y=labels, orientation="h", marker_color=_POST_COLOR))
        fig.update_layout(title=title, dragmode="pan",
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
            """The time-ordered observation names/times/ensembles for one observation group."""

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
            fig, axes = viz.mpl_axes(len(chosen), 1, figsize=(8, 2.6 * len(chosen)), squeeze=False)
            for ax, group in zip(axes[:, 0], chosen, strict=False):
                names, times, measured, flags = _group_data(group)
                for real in prior.index:
                    ax.plot(times, prior.loc[real, names].to_numpy(dtype=float), color=_MPL_PRIOR, lw=0.8, alpha=0.4)
                ax.plot(times[~flags], measured[~flags], "^", color=_MPL_MEAS, ms=7, label="measured")
                if flags.any():
                    ax.plot(times[flags], measured[flags], "X", color=_MPL_CONFLICT, ms=10, label="prior-data conflict")
                ax.set_title(group, loc="left", fontsize=9)
                ax.set_ylabel("value")
            axes[-1, 0].set_xlabel("time")
            fig.suptitle("Prior ensemble vs measured observations")
            fig.tight_layout()
            return fig

        fig = viz.subplots(rows=len(chosen), cols=1, subplot_titles=chosen, shared_xaxes=False)
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
                                marker=dict(color=_CONFLICT_COLOR, size=11, symbol="x"),
                                row=row, col=1, name="prior-data conflict", legendgroup="conflict", showlegend=row == 1)
        fig.update_layout(title="Prior ensemble vs measured observations",
                          dragmode="pan", height=260 * len(chosen))
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

            fig, ax = viz.mpl_axes(figsize=(6.5, max(2.5, 0.4 * len(labels) + 1)))
            ax.barh(labels, pct, color=_MPL_CONFLICT)
            ax.set_xlabel("% of observations in prior-data conflict")
            ax.set_title(title)
            ax.invert_yaxis()
            sns.despine(fig)
            return fig

        fig = viz.Fig().add_trace(go.Bar(x=pct, y=labels, orientation="h", marker_color=_CONFLICT_COLOR))
        fig.update_layout(title=title, dragmode="pan",
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

    def _field_cell_index(self, info: dict, layer: int) -> dict[str, int]:
        """Map captured-field observation names to cell ids for one layer."""

        if self.model is None:
            raise ValueError(
                "Spatial maps need the model grid. Open the run via "
                "cal.run_ies(...) or model.pest_runs[i].review() (both pass the "
                "model), or name one yourself with review(model=...) / "
                "IesResults(..., model=...)."
            )
        ncpl = int(self.model.vor.ncpl)
        # JSON round-trips the int layer keys to strings.
        layer_prefixes = {int(k): v for k, v in (info.get("layer_prefixes") or {}).items()}
        mapping: dict[str, int] = {}
        if layer_prefixes:
            # Multi-layer DISV: one capture file per layer, each indexed
            # 0..ncpl-1, so the array index *is* the cell id.
            prefix = layer_prefixes.get(int(layer))
            if prefix is None:
                raise ValueError(
                    f"No captured field for layer {layer}. "
                    f"Captured layers: {sorted(layer_prefixes)}."
                )
            token = f"oname:{prefix.lower()}_"
            for name in self.pst.observation_data.index:
                if token not in name:
                    continue
                match = self._ARR_RE.search(name)
                if match:
                    mapping[name] = int(match.group(1))
        else:
            # A single whole-grid array covering every layer at once, which is
            # how MF6 writes griddata with no LAYERED keyword (MST porosity):
            # flat = layer * ncpl + cell.
            token = f"oname:{info['prefix'].lower()}_"
            for name in self.pst.observation_data.index:
                if token not in name:
                    continue
                match = self._ARR_RE.search(name)
                if match and int(match.group(1)) // ncpl == int(layer):
                    mapping[name] = int(match.group(1)) % ncpl
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
            One row per captured cell -- ``cell`` is a COLUMN, not the index
            (the frame carries a fresh RangeIndex) -- with columns
            ``prior_mean``, ``prior_std``, ``posterior_mean``,
            ``posterior_std``, ``change`` (``posterior_mean / prior_mean``),
            ``reduction`` (``1 - posterior_std / prior_std`` -- the share of the
            prior spread the data removed), and ``base`` (the
            minimum-error-variance realization) when available. Cells that were
            not captured are absent, so the frame is generally shorter than
            ``ncpl``.
        """

        info = self._capture_info(target)
        if info.get("family") != "array":
            raise NotImplementedError(
                f"Spatial maps currently support array fields (K, K33); "
                f"{info.get('target')!r} is a list field."
            )
        cells = self._field_cell_index(info, layer)
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
        # How much of the prior spread the data removed: 1 = the ensemble
        # collapsed (fully informed), 0 = the data said nothing here. NEGATIVE
        # where the posterior spread GREW, which is real and is not clipped away.
        # An unguarded divide, like `change`: a zero prior sd leaves NaN/+-inf,
        # which the map renders as a named reason rather than a silent gap.
        with np.errstate(divide="ignore", invalid="ignore"):
            frame["reduction"] = 1.0 - frame["posterior_std"] / frame["prior_std"]
        return frame.sort_values("cell").reset_index(drop=True)

    def _field_choro(self, target: str, *, stat: str, which: str, layer: int,
                     **choropleth_kwargs):
        """Build one field map as a ``Choro``, with its heading and resolved kwargs.

        The seam :meth:`plot_field` and :meth:`plot_field_mosaic` share: a mosaic
        needs the ``Choro`` itself (``viz.mosaic`` reads ``get_choropleth()``),
        while ``plot_field`` renders it to a figure. Returns
        ``(choro, heading, kwargs)`` -- the kwargs so the caller can see which
        policy keys actually took effect after ``setdefault``.
        """

        frame = self.field(target, layer=layer)
        column = {
            "mean": f"{which}_mean",
            "std": f"{which}_std",
            "base": "base",
            "change": "change",
            "reduction": "reduction",
        }.get(str(stat).lower())
        if column is None or column not in frame.columns:
            raise ValueError(
                f"Unknown stat {stat!r} (or unavailable). Use 'mean', 'std', "
                "'base', 'change', or 'reduction'."
            )

        stat_key = str(stat).lower()
        ncpl = int(self.model.vor.ncpl)
        cells = frame["cell"].to_numpy(dtype=int)

        def _scatter(name: str) -> np.ndarray:
            """One captured column spread onto the full grid, uncaptured cells NaN."""

            # field() has one row per CAPTURED cell; every hover column must be
            # exactly ncpl long or the hover assembler raises.
            spread = np.full(ncpl, np.nan)
            spread[cells] = frame[name].to_numpy(dtype=float)
            return spread

        values = _scatter(column)
        label = f"{target} {_STAT_LABELS.get(stat_key, stat_key)}" + (
            f" ({which})" if stat_key in ("mean", "std") else ""
        )

        # RAW values for the hover -- only custom_zs is log-transformed, so the
        # reader sees "2.5x", not "0.4".
        payload = {
            column: _hover_column(
                values, reason=_UNDEFINED_REASONS.get(stat_key, "undefined")
            )
        }
        for name in ("prior_mean", "posterior_mean", "posterior_std"):
            if name in frame.columns and name != column:
                payload[name] = _hover_column(_scatter(name))

        # Ensemble row counts actually on disk -- NOT settings.num_reals, which
        # is the REQUESTED count and is None for a run myflopy did not launch.
        n_prior = int(self.prior._df.shape[0])
        n_post = int(self.posterior._df.shape[0])
        if stat_key in ("change", "reduction"):
            # Both compare the two ensembles, so both span the two iterations.
            reals = f"{n_prior}" if n_prior == n_post else f"{n_prior}→{n_post}"
            context = (
                f"iterations {self.prior_iteration}→{self.posterior_iteration}, "
                f"{reals} realizations"
            )
        else:
            is_prior = stat_key in ("mean", "std") and which == "prior"
            iteration = self.prior_iteration if is_prior else self.posterior_iteration
            context = f"iteration {iteration}, {n_prior if is_prior else n_post} realizations"
        heading = f"{label} — {context}"

        # setdefault, never direct kwargs: these names now REACH Choro, so an
        # explicit plot_field(..., colorscale=..., logscale=False) must override
        # the policy rather than raise "multiple values for keyword".
        choropleth_kwargs.setdefault("custom_hover", payload)
        choropleth_kwargs.setdefault(
            "hover_spec", parameter_field_hover(column, title=heading)
        )
        choropleth_kwargs.setdefault("hover_heads", False)
        choropleth_kwargs.setdefault("hover_ks", False)
        policy = _field_map_policy(stat_key, values)
        # An explicit `logscale=False` must not keep the policy's decade
        # colorbar: its tickvals are positions in LOG space, so over linear data
        # every label crushes into the bottom of the bar (and a large K can
        # overflow `10 ** tick` outright).
        if choropleth_kwargs.get("logscale", policy.get("logscale")) is False:
            policy.pop("colorbar", None)
        for key, value in policy.items():
            if value is not None:
                choropleth_kwargs.setdefault(key, value)

        # `model` is deliberately not passed: Choro reads the model's head output
        # whenever model is not None, so a PEST model with missing or stale heads
        # would raise at construction -- and a parameter field is not a head map.
        choro = build_choropleth(
            self.model.vor, custom_zs=list(values), layer=layer, **choropleth_kwargs
        )
        return choro, heading, choropleth_kwargs

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
            ``"mean"`` (ensemble mean field), ``"std"`` (spread -- where the
            property is still uncertain), ``"base"`` (the minimum-error-variance
            realization), ``"change"`` (posterior_mean / prior_mean -- where
            calibration moved the property), or ``"reduction"``
            (``1 - posterior_std / prior_std`` -- how much of the prior spread
            the data removed, i.e. *did the data inform this region*).
        which
            ``"prior"`` or ``"posterior"`` for ``stat`` in {``"mean"``, ``"std"``}.
            Ignored by the others: ``change`` and ``reduction`` are already
            prior-to-posterior comparisons, and ``base`` is a single posterior
            realization.
        layer
            Model layer to map (default 0).
        backend
            ``"plotly"`` (interactive map, default) or ``"matplotlib"`` (static
            matplotlib choropleth of the Voronoi cells).
        **choropleth_kwargs
            Forwarded to the Voronoi choropleth builder; any of them overrides
            the house color policy for this map.

        See Also
        --------
        plot_field_mosaic : the same map for prior and posterior, side by side.
        """

        choro, heading, resolved = self._field_choro(
            target, stat=stat, which=which, layer=layer, **choropleth_kwargs
        )
        if _normalize_backend(backend) == "matplotlib":
            # No cmap=: plot_mpl derives it from the policy colorscale (stops
            # become a LinearSegmentedColormap, "earth" becomes "gist_earth") and
            # vmin/vmax default to the policy's zmin/zmax, so the colors and the
            # limits agree between backends. The tick LABELS do not come for
            # free -- plot_mpl has no colorbar hook, hence the relabel.
            figure = choro.plot_mpl(title=heading)
            if resolved.get("logscale"):
                _relabel_log_colorbar(figure, resolved.get("colorbar"))
            return figure
        # `title` is matplotlib-only and is not a Choroplethmap property, so it
        # must never ride the trace kwargs. plot() sets a 20px top margin; bump it
        # or the title clips.
        return choro.plot().update_layout(title=heading, margin={"t": 40})

    def plot_field_mosaic(self, target: str, *,
                          which: Sequence[str] = ("prior", "posterior"),
                          stat: str = "mean", layer: int = 0, ncols: int = 2,
                          title: str | None = None, sync_views: bool = True,
                          backend: str = "plotly", **choropleth_kwargs):
        """Compose one field map per ensemble, side by side on a shared color scale.

        The prior-vs-posterior figure: the same property, the same colors, the
        same limits, so the panels are actually comparable and the eye can only
        be reading a real difference.

        Parameters
        ----------
        target
            A captured target name, e.g. ``"k"``.
        which
            The ensembles to compose, one panel each (default prior then
            posterior).
        stat
            ``"mean"`` or ``"std"`` -- the only stats that HAVE a prior and a
            posterior form. ``"change"`` and ``"reduction"`` are already
            prior-to-posterior comparisons, and ``"base"`` is a single posterior
            realization with no prior counterpart; all three are refused here.
        layer, ncols, title, sync_views
            Grid layer, mosaic width, overall title, and whether the panels pan
            and zoom together (see :func:`myflopy.viz.mosaic`).
        backend
            ``"plotly"`` only -- ``viz.mosaic`` composes Plotly subplots. Use
            :meth:`plot_field` per panel for static figures.

        Notes
        -----
        A mosaic pools every panel onto ONE color axis, which discards each
        panel's own limits and colorbar (compromise ledger 71). That pooling is
        the point -- it is what makes the panels comparable -- but it means the
        real-unit ticks for a log-scaled ``mean`` mosaic have to be rebuilt from
        the pooled limits, which is what the ``colorbar`` callback does.
        """

        if _normalize_backend(backend) != "plotly":
            raise ValueError(
                "plot_field_mosaic composes Plotly subplots; there is no "
                "matplotlib mosaic. Call plot_field(..., backend='matplotlib') "
                "per panel instead."
            )
        stat_key = str(stat).lower()
        if stat_key not in ("mean", "std"):
            raise ValueError(
                f"stat={stat!r} has no separate prior and posterior form to "
                "compose: 'change' and 'reduction' are already prior-to-"
                "posterior comparisons, and 'base' is one posterior "
                f"realization. Use plot_field(..., stat={stat!r})."
            )
        sides = [str(side).lower() for side in which]
        if len(sides) < 2:
            raise ValueError(
                "A mosaic needs at least two panels, e.g. "
                "which=('prior', 'posterior')."
            )

        panels, resolved = [], {}
        for side in sides:
            choro, heading, resolved = self._field_choro(
                target, stat=stat_key, which=side, layer=layer, **choropleth_kwargs
            )
            panels.append((heading, choro))

        return viz.mosaic(
            panels,
            ncols=ncols,
            title=title or f"{target} {stat_key} — {' vs '.join(sides)}",
            sync_views=sync_views,
            # Keyed on the EFFECTIVE logscale, not on the stat: `mean` is mapped
            # in log space by default (without this the shared bar reads
            # "-3 ... 2" for a field running 0.001 to 100), but a caller passing
            # `logscale=False` gets linear data, where decade labels would be
            # wrong and `10 ** cmax` can overflow.
            colorbar=(
                _log_decade_colorbar_for_mosaic if resolved.get("logscale") else None
            ),
        )

    # Prefix, location name and time in a pyEMU list-style observation name:
    # "oname:hds_otype:lst_usecol:obs_00_per:0" -> ("hds", "obs_00", "0").
    # Parsed from the NAME rather than read off pst.try_parse_name_metadata()'s
    # `usecol` column, which truncates at the first underscore -- "obs_00"
    # arrives there as "obs", useless as a join key for any real location name.
    _OBS_NAME_RE = re.compile(
        r"oname:(?P<prefix>.+?)_otype:.*?usecol:(?P<location>.+)_(?:per|time):(?P<time>[^_]+)$"
    )

    _LOCATION_COLUMNS = (
        "prefix", "kind", "value_kind", "location", "layer", "cell", "x", "y",
        "cells",
    )

    #: Geometry for observation families built before the metadata carried a
    #: ``geometry`` key, so runs already on disk keep reviewing. New families
    #: declare it in their metadata instead of being added here.
    _LEGACY_GEOMETRY_BY_KIND = {
        "head_targets": "points",
        "conc_targets": "points",
        "drn_flow": "zones",
    }

    #: Same, for the value each family measures -- it names the residual's UNIT,
    #: which is what makes mixing families on one color scale detectable.
    _LEGACY_VALUE_BY_KIND = {
        "head_targets": "head",
        "conc_targets": "conc",
        "drn_flow": "flow_target",
    }

    def _simulated_realization(self, realization: str) -> pd.Series:
        """One realization's simulated values, falling back to the ensemble mean."""

        ensemble = self.posterior._df
        if realization in ensemble.index:
            return ensemble.loc[realization]
        return ensemble.mean()

    def _observation_locations(self) -> pd.DataFrame:
        """Grid locations for every observation set that recorded one.

        Dispatch is on the set's ``geometry``, not on its kind: **point**
        families (head, concentration) snapshot a point GeoPackage and a
        name -> cell map; **zone** families (DRN) snapshot the explicit cell
        list per zone. Anything that declares neither is skipped -- lake stage
        and SFR stage/flow record only a lake or reach NUMBER, which is not a
        place on the grid (compromise ledger 76).

        Runs written before the metadata carried ``geometry`` fall back to
        ``_LEGACY_GEOMETRY_BY_KIND`` so they keep reviewing.
        """

        import geopandas as gpd

        frames = []
        for entry in self.observation_sets:
            kind = str(entry.get("kind", ""))
            prefix = str(entry.get("prefix", ""))
            geometry = str(
                entry.get("geometry") or self._LEGACY_GEOMETRY_BY_KIND.get(kind, "")
            )
            if not prefix or geometry not in ("points", "zones"):
                continue
            locations_file = entry.get("locations_file")
            path = self.workspace / str(locations_file) if locations_file else None
            if path is None or not path.exists():
                continue
            table = (
                gpd.read_file(path) if path.suffix.lower() == ".gpkg"
                else pd.read_csv(path)
            )
            if table.empty:
                continue
            frame = pd.DataFrame({
                "prefix": prefix.lower(),
                "kind": kind,
                "value_kind": str(
                    entry.get("value_column")
                    or self._LEGACY_VALUE_BY_KIND.get(kind, kind)
                ),
                # The control file lowercases location names when it builds
                # observation names; the snapshot keeps the original case.
                "location": table["name"].astype(str).str.strip().str.lower(),
            })
            if geometry == "points":
                if hasattr(table, "geometry"):
                    frame["x"] = table.geometry.x.to_numpy()
                    frame["y"] = table.geometry.y.to_numpy()
                if "layer" in table.columns:
                    frame["layer"] = table["layer"].to_numpy()
                # Each family names its own map file (heads and concentration
                # do not share one); the default is the pre-`mapping_file`
                # layout, which only heads ever wrote.
                mapping_path = self.workspace / str(
                    entry.get("mapping_file") or f"{prefix}_head_target_map.csv"
                )
                if mapping_path.exists():
                    mapping = pd.read_csv(mapping_path)
                    mapping["location"] = (
                        mapping["name"].astype(str).str.strip().str.lower()
                    )
                    # One row per location before merging: a snapshot listing a
                    # name twice (duplicated source row, two screens) would
                    # otherwise cross-join into n**2 rows against a single pyEMU
                    # observation, and the map would draw n identical markers.
                    frame = frame.drop_duplicates(subset="location").merge(
                        mapping.drop_duplicates(subset="location")
                               .loc[:, ["location", "cell"]],
                        on="location", how="left",
                    )
            else:
                # Snapshotted as a comma-joined string, one zone per row -- but a
                # zone that resolved to no cells round-trips through CSV as NaN,
                # whose str() is "nan". Parsing that unguarded took out the whole
                # frame, every head target included, over one degenerate zone.
                frame["cells"] = [
                    _parse_cell_list(value) for value in table["cells"]
                ]
            frames.append(frame)

        if not frames:
            return pd.DataFrame(columns=list(self._LOCATION_COLUMNS))
        merged = pd.concat(frames, ignore_index=True)
        return merged.reindex(columns=list(self._LOCATION_COLUMNS))

    def _measured_observation_keys(self) -> set[tuple[str, str, str]] | None:
        """``(prefix, location, time)`` for observations that carry a MEASUREMENT.

        pyEMU creates one observation per row of the simulated output, not per
        row of the target table, and ``_assign_target_values`` only overwrites
        the ones it has a target for (`observations.py`, ``if obsnme not in
        obs.index: continue``). Every other row keeps pyEMU's defaults —
        **weight 1.0 and an ``obsval`` equal to the base model's own simulated
        head**. Averaging those into a location's residual pulls it toward zero
        and can flip its sign, so they are excluded here by joining back to the
        `` <prefix>_target_values.csv`` snapshot of what was actually measured.

        Returns ``None`` when no observation set recorded a values file, meaning
        "cannot tell" — the caller then does not filter, because dropping
        everything would be worse than the bias.
        """

        keys: set[tuple[str, str, str]] = set()
        found = False
        for entry in self.observation_sets:
            values_file = entry.get("values_file")
            prefix = str(entry.get("prefix", ""))
            if not values_file or not prefix:
                continue
            path = self.workspace / str(values_file)
            if not path.exists():
                continue
            table = pd.read_csv(path)
            if not {"time", "name"}.issubset(table.columns):
                continue
            found = True
            # The observation name's trailing token is str(time) verbatim --
            # `_build_index_row_labels` formats it as f"{index}:{value}" -- so
            # comparing the stringified time avoids having to know whether this
            # family indexed on `per` (heads) or `time` (named series).
            for name, time in zip(table["name"], table["time"], strict=True):
                keys.add((
                    prefix.lower(),
                    str(name).strip().lower(),
                    str(time).strip().lower(),
                ))
        return keys if found else None

    def _select_prefixes(self, locations: pd.DataFrame, prefix) -> pd.DataFrame:
        """Filter located observations to the requested set prefix(es)."""

        if prefix is None:
            return locations
        wanted = (
            {str(prefix).lower()} if isinstance(prefix, str)
            else {str(value).lower() for value in prefix}
        )
        selected = locations[locations["prefix"].isin(wanted)]
        if selected.empty:
            available = sorted(locations["prefix"].unique())
            raise ValueError(
                f"No locatable observation set matches prefix={prefix!r}. "
                f"This run recorded: {available}."
            )
        return selected

    def obs_residuals(self, *, realization: str = "base", prefix=None) -> pd.DataFrame:
        """Return per-location residuals for observations that can be placed on the grid.

        The residual is ``simulated - measured`` -- the sign
        :meth:`phi_contributions` already uses. (PEST's own ``.res`` file reports
        ``measured - modelled``; rather than sign one into the other, every frame
        and label here names which convention it is in.)

        Observations are summarized over TIME: a head target measured in six
        stress periods becomes one row whose ``residual`` is the mean of the six,
        with ``n`` recording how many went into it. A map needs one value per
        place, and a mean residual is the bias there.

        Only **measured, history-matched** observations count. pyEMU creates an
        observation for every row of the simulated output, so a target measured
        at one time out of six leaves five rows holding the base model's own
        output at weight 1.0; those are excluded by joining back to the
        ``*_target_values.csv`` snapshot. Zero-weight observations are excluded
        too, which is what keeps ``cal.forecast(...)`` sets — registered as
        ordinary observation sets — off a misfit figure.

        Parameters
        ----------
        realization
            Which posterior realization to score (default ``"base"``, the
            minimum-error-variance one); falls back to the ensemble mean when
            that label is absent, exactly as :meth:`phi_contributions` does.
        prefix
            One observation-set prefix, or several, to restrict the frame to.
            Default ``None`` keeps every locatable set — which mixes UNITS when
            a run history-matched more than one kind of measurement (heads in
            length, concentration in mass/volume, DRN seepage in volume/time).

        Returns
        -------
        pandas.DataFrame
            One row per (``prefix``, ``location``) with ``kind``,
            ``value_kind``, ``measured``, ``simulated``, ``residual``,
            ``weight``, ``n``, and whichever location the run recorded:
            ``cell``/``x``/``y`` for point targets (heads, concentration),
            ``cells`` (the zone's cell list) for DRN zones. Empty when the run
            recorded no locatable observation sets.
        """

        locations = self._observation_locations()
        if locations.empty:
            return pd.DataFrame(
                columns=["prefix", "location", "kind", "value_kind", "measured",
                         "simulated", "residual", "weight", "n", "cell", "x",
                         "y", "cells"]
            )
        locations = self._select_prefixes(locations, prefix)

        simulated = self._simulated_realization(realization)
        obs = self._observation_metadata()
        obs = obs.loc[[name for name in obs.index if name in simulated.index]].copy()
        obs["simulated"] = simulated.reindex(obs.index).astype(float).to_numpy()
        parsed = obs.index.to_series().str.extract(self._OBS_NAME_RE)
        obs["prefix"] = parsed["prefix"].str.lower().to_numpy()
        obs["location"] = parsed["location"].str.lower().to_numpy()
        obs["time"] = parsed["time"].str.lower().to_numpy()
        obs = obs.dropna(subset=["prefix", "location"])
        obs["obsval"] = obs["obsval"].astype(float)
        obs["weight"] = obs["weight"].astype(float)

        measured = self._measured_observation_keys()
        if measured is not None:
            obs = obs.loc[[
                key in measured for key in
                zip(obs["prefix"], obs["location"], obs["time"], strict=True)
            ]]
        # Zero weight means PEST did not history-match it: forecasts (which
        # `cal.forecast` registers as ordinary observation sets, so they arrive
        # here looking like head targets) and deliberately silenced targets. A
        # misfit map must not draw either as if the model had been fitted to it.
        obs = obs.loc[obs["weight"] > 0]

        summary = obs.groupby(["prefix", "location"], as_index=False).agg(
            measured=("obsval", "mean"),
            simulated=("simulated", "mean"),
            weight=("weight", "mean"),
            n=("obsval", "size"),
        )
        summary["residual"] = summary["simulated"] - summary["measured"]
        merged = summary.merge(locations, on=["prefix", "location"], how="inner")
        return merged.sort_values(["prefix", "location"]).reset_index(drop=True)

    def plot_obs_residuals(self, *, realization: str = "base", layer: int = 0,
                           backend: str = "plotly", title: str | None = None,
                           prefix=None, **choropleth_kwargs):
        """Map where the calibrated model is biased, and by how much.

        Point targets (heads, concentration) draw as markers at their
        coordinates; DRN zones color the cells they cover. Both read on ONE
        diverging scale centered on zero, so a point and the zone beneath it
        mean the same thing at the same color: **red is under-simulated**
        (simulated below measured), blue is over-simulated, white is on the
        money.

        That one shared scale is why ``prefix`` matters. A run that
        history-matched heads AND concentration produces residuals in different
        UNITS, and drawing them together puts the larger-magnitude family in
        charge of the color limit — the other one renders uniformly white and
        reads as a perfect fit. Passing ``prefix`` restricts the figure to one
        family; leaving it out draws them all and warns.

        Parameters
        ----------
        realization
            Posterior realization to score (default ``"base"``).
        prefix
            One observation-set prefix, or several, to draw. Default ``None``
            draws every locatable set.
        layer
            Grid layer to draw the cells of (default 0). Head targets are drawn
            at their map position whatever layer they were screened in --
            filtering them would hide data on a plan-view figure -- so read the
            layer off the hover rather than off the map.
        backend
            ``"plotly"`` (interactive, default) or ``"matplotlib"``.
        title
            Overrides the default heading.
        **choropleth_kwargs
            Forwarded to the Voronoi choropleth builder.

        See Also
        --------
        obs_residuals : the same numbers as a DataFrame.
        """

        frame = self.obs_residuals(realization=realization, prefix=prefix)
        if frame.empty:
            # Distinguish "nothing locatable was recorded" from "locations were
            # found but nothing joined to them" -- the second is a bug or a
            # weights-all-zero run, and blaming lake/SFR targets for it would
            # send the reader to the wrong place entirely.
            if self._observation_locations().empty:
                raise ValueError(
                    "No observation locations were recorded for this run. A "
                    "residual map needs the *_target_locations.gpkg / "
                    "*_head_target_map.csv snapshots written beside the control "
                    "file, and only point targets (heads, concentration) and "
                    "DRN zones record geometry -- lake and SFR targets record "
                    "only a lake or reach number."
                )
            raise ValueError(
                "Observation locations were found, but no measured, nonzero-"
                "weight observation joined to them. Check that this run has "
                "history-matched observations (forecasts carry weight 0) and "
                "that the *_target_values.csv snapshots are present."
            )

        # One color limit over families measured in different units hands the
        # scale to whichever has the larger numbers, and silently renders the
        # other as "no misfit anywhere". Warn rather than refuse: runs built
        # before `prefix` existed already draw heads and DRN seepage together,
        # and refusing would break reviewing them.
        drawn_kinds = sorted(frame["value_kind"].dropna().unique()) if "value_kind" in frame else []
        if len(drawn_kinds) > 1:
            warnings.warn(
                f"Drawing residuals of different kinds ({drawn_kinds}) on one "
                "color scale; they are in different units, so the "
                "larger-magnitude family sets the limit and the others read as "
                "near-zero misfit. Pass prefix=... to map one family at a time.",
                UserWarning,
                stacklevel=2,
            )

        points = frame[frame["x"].notna()] if "x" in frame else frame.iloc[:0]
        zones = frame[frame["cells"].notna()] if "cells" in frame else frame.iloc[:0]

        # ONE limit over every residual, points and zones together: two scales
        # on one figure would make a point and the cell under it different
        # colors for the same number.
        residuals = frame["residual"].to_numpy(dtype=float)
        half = _symmetric_color_limit(residuals) or 1.0
        colorscale = _red_white_blue_diverging_colorscale()

        ncpl = int(self.model.vor.ncpl)
        cell_row = np.full(ncpl, -1, dtype=int)
        for position, row in enumerate(zones.itertuples()):
            for cell in row.cells:
                if 0 <= int(cell) < ncpl:
                    cell_row[int(cell)] = position

        def _zone_column(name: str) -> np.ndarray:
            """One zone column spread over the cells it covers, elsewhere NaN."""

            spread = np.full(ncpl, np.nan)
            if not zones.empty:
                mask = cell_row >= 0
                spread[mask] = zones[name].to_numpy(dtype=float)[cell_row[mask]]
            return spread

        zone_values = _zone_column("residual")
        payload = {
            name: _hover_column(_zone_column(name))
            for name in ("residual", "measured", "simulated", "weight")
        }

        used = (
            realization if realization in self.posterior._df.index else "ensemble mean"
        )
        heading = title or (
            f"residuals (simulated − measured) — iteration "
            f"{self.posterior_iteration}, {used}, {len(frame)} locations"
        )

        choropleth_kwargs.setdefault("custom_hover", payload)
        choropleth_kwargs.setdefault("hover_spec", residual_hover(title=heading))
        choropleth_kwargs.setdefault("hover_heads", False)
        choropleth_kwargs.setdefault("hover_ks", False)
        choropleth_kwargs.setdefault("colorscale", colorscale)
        choropleth_kwargs.setdefault("zmin", -half)
        choropleth_kwargs.setdefault("zmax", half)
        choropleth_kwargs.setdefault("zmid", 0.0)
        # With no zones every cell is NaN, so the cells' own scale bar would be
        # an empty legend that reads as a broken figure. Hide it and let the
        # point markers carry the one scale instead -- `showscale`, not
        # `colorbar`, which on the trace is a dict of tick properties.
        if zones.empty:
            choropleth_kwargs.setdefault("showscale", False)

        choro = build_choropleth(
            self.model.vor, custom_zs=list(zone_values), layer=layer,
            **choropleth_kwargs
        )

        if _normalize_backend(backend) == "matplotlib":
            # Always a colorbar here, even with no zones and an all-NaN cell
            # column: geopandas draws a correct bar from the explicit
            # vmin/vmax, and it is the ONLY color key this backend has --
            # plot_mpl ignores overlays, so the scattered points cannot carry
            # one the way the Plotly markers do.
            figure = choro.plot_mpl(title=heading)
            if not points.empty:
                # axes[0] is the map (plot_mpl appends the colorbar after it),
                # drawn in MODEL coordinates -- no reprojection, unlike the
                # Plotly path. Same colormap as the cells, from the same stops.
                figure.axes[0].scatter(
                    points["x"].to_numpy(dtype=float),
                    points["y"].to_numpy(dtype=float),
                    c=points["residual"].to_numpy(dtype=float),
                    cmap=mpl_colormap_for(colorscale), vmin=-half, vmax=half,
                    s=45, edgecolor="black", linewidth=0.5, zorder=5,
                )
            return figure

        if not points.empty:
            lon, lat = self.model.vor.points_to_latlon(
                points["x"].to_numpy(dtype=float), points["y"].to_numpy(dtype=float)
            )
            choro.add_overlay(go.Scattermap(
                lon=lon, lat=lat, mode="markers", name="observations",
                marker={
                    "color": points["residual"].to_numpy(dtype=float),
                    "colorscale": colorscale, "cmin": -half, "cmax": half,
                    "size": 11, "opacity": 0.95,
                    # Exactly one scale bar on the figure: the cells' when there
                    # are zones, the markers' when there are not.
                    "showscale": bool(zones.empty),
                },
                customdata=[
                    [row.location, float(row.measured), float(row.simulated),
                     float(row.residual), int(row.n),
                     "" if pd.isna(row.layer) else int(row.layer)]
                    for row in points.itertuples()
                ],
                hovertemplate=(
                    "<b>%{customdata[0]}</b><br>"
                    "residual (sim − meas): %{customdata[3]:.4g}<br>"
                    "measured: %{customdata[1]:.4g}<br>"
                    "simulated: %{customdata[2]:.4g}<br>"
                    "observations: %{customdata[4]}<br>"
                    "layer: %{customdata[5]}<extra></extra>"
                ),
            ))
        return choro.plot().update_layout(title=heading, margin={"t": 40})

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

        Bundles, in order: phi convergence, the phi distribution, the
        ensemble-vs-observation comparison, phi contributions, parameters at
        bounds, one histogram per forecast, and -- when the run carries a model
        and captured fields -- the posterior mean and change maps of each.

        The bundle is the headline set, not every figure on the class -- among
        others it omits :meth:`plot_field` ``stat="reduction"``,
        :meth:`plot_field_mosaic`, :meth:`plot_obs_residuals`,
        :meth:`plot_prior_vs_obs` and :meth:`plot_conflict`. Call those directly.

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
            except (OSError, ValueError, TypeError, KeyError) as error:
                # pragma: no cover - diagnostics are best-effort in the bundle
                # A run killed between writing its obs and par ensembles leaves
                # the .par.csv missing (OSError); an object-dtype column makes
                # the mean a TypeError, not a ValueError.
                logger.debug(
                    "leaving %s out of the report bundle: %s",
                    getattr(diagnostic, "__name__", diagnostic), error,
                )
        for name in self.forecast_names:
            figures.append(self.forecast(name).plot())
        if self.model is not None:
            for info in self.capture_fields:
                if info.get("family") != "array":
                    continue
                try:
                    figures.append(self.plot_field(info["target"], stat="mean", which="posterior"))
                    figures.append(self.plot_field(info["target"], stat="change"))
                except (OSError, ValueError, TypeError, KeyError, AttributeError,
                        IndexError, NotImplementedError) as error:
                    # pragma: no cover - field maps are best-effort in the bundle
                    # NotImplementedError is reachable even though the loop
                    # skips non-array families: `plot_field` is handed
                    # `info["target"]` as a STRING and re-resolves it, so the
                    # family it lands on need not be the one filtered above.
                    logger.debug(
                        "no field map for %s: %s", info.get("target"), error,
                    )

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

    Convenience wrapper around :class:`IesResults`. Usually reached through
    ``model.pest_runs[i].review()`` rather than called directly.
    Pass ``model=`` (a myflopy model carrying the grid) to enable spatial
    parameter maps via :meth:`IesResults.plot_field`.
    """

    return IesResults(workspace, case_name=case_name, model=model)

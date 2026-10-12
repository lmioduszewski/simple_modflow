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
  (iteration ``0`` is the prior; the highest iteration the phi file records is the
  posterior);
- ``<case>.<iter>.par.csv`` -- the parameter ensemble at each iteration;
- ``<case>.obs+noise.csv`` -- the measurement-noise ensemble.

Each ensemble is read as ``.jcb`` too, the form ``ies_save_binary`` writes.
"""

from __future__ import annotations

import json
import re
import warnings
from collections.abc import Sequence
from dataclasses import dataclass
from functools import cached_property
from numbers import Integral
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from myflopy import viz
from myflopy._logging import get_logger
from myflopy._optional import require
from myflopy.modflow.mf6.grid.geometry import model_cell_polygons
from myflopy.modflow.mf6.grid.plotting import _choropleth_factory
from myflopy.modflow.mf6.grid.voronoi import model_grid_view
from myflopy.modflow.mf6.package_plotting import (
    _red_white_blue_diverging_colorscale,
    _symmetric_color_limit,
    normalize_backend,
)
from myflopy.modflow.mf6.pest.forward_run import read_target_table
from myflopy.modflow.utils.datatypes.choros import mpl_colormap_for
from myflopy.modflow.utils.datatypes.hover import (
    parameter_field_hover,
    residual_hover,
)
from myflopy.modflow.utils.validators import whole_number
from myflopy.viz import log_decade_colorbar

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

    # The one resolver (``package_plotting.normalize_backend``): this kept its own
    # alias list -- "go"/"sns" in, "interactive"/"static" out -- until ledger
    # 480. Surrounding spaces are stripped first, as they always were here.
    return "matplotlib" if normalize_backend(str(backend).strip()) == "mpl" else "plotly"


def _which_ensemble(which) -> str:
    """``"prior"`` or ``"posterior"`` (any case); anything else raises (ledger 480).

    ``parameters_at_bounds`` read the prior for ANY value but ``"posterior"`` -- a typo
    or ``"Posterior"`` included -- and ``plot_field`` raised "Unknown stat" for one.
    """

    side = str(which).strip().lower()
    if side not in ("prior", "posterior"):
        raise ValueError(f"which={which!r} is not an ensemble; use 'prior' or 'posterior'.")
    return side


def _histogram_bars(fig, values, edges, *, density: bool, name: str, color: str) -> None:
    """Draw ``values`` binned on ``edges`` as bars on a Plotly ``fig``.

    The bins are numpy's, so the Plotly figure takes whatever ``bins=`` numpy takes
    -- a count, a rule name, explicit (even unequal) edges -- exactly as the
    Matplotlib branch does. ``xbins`` divided the range by ``bins``, so a rule name
    raised ``TypeError`` and explicit edges drew a meaningless width (ledger
    480).
    """

    counts, _ = np.histogram(values, bins=edges, density=density)
    fig.add_bar(
        x=(edges[:-1] + edges[1:]) / 2.0, y=counts, width=np.diff(edges), name=name,
        marker_color=color, opacity=0.6,
    )


def _log_decade_colorbar(low: float, high: float, *, unit: str = "") -> dict:
    """Decade ticks in real units for a log10 colorbar; ``low``/``high`` in log10.

    Since 2026-09-13 ``Choro`` labels every log map this way itself
    (:func:`myflopy.viz.log_decade_colorbar`), so this is only still passed for
    the ``change`` map, whose ``unit="x"`` suffix the default cannot know. It is
    the same helper, so the two cannot label one field two ways.
    """

    return log_decade_colorbar(low, high, unit=unit)


def _log_decade_colorbar_for_mosaic(low, high) -> dict:
    """``viz.mosaic`` colorbar callback: decades over the POOLED panel limits.

    A mosaic pools its panels onto one shared color axis, so the real-unit ticks
    have to be computed from the limits *it* resolved, not from any one panel's
    (compromise ledger 71). ``viz.mosaic`` leaves the limits unset when no panel
    has finite data, hence the guard. ``viz.mosaic`` now does this by itself for
    a grid of log-scaled panels; kept as the explicit spelling.
    """

    if low is None or high is None:
        return {}
    return _log_decade_colorbar(float(low), float(high))


# Display names for stats whose column name reads badly as a figure label.
_STAT_LABELS = {"reduction": "uncertainty reduction"}

# The same, for a field whose statistics :meth:`IesResults.field` took in log10.
_LOG_STAT_LABELS = {
    "mean": "geometric mean",
    "std": "std of log10",
    "change": "change (geometric)",
    "reduction": "uncertainty reduction (log10)",
}

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
            # Limits in DATA units (a log map's zmin/zmax/zmid are, since
            # 2026-09-13): symmetric about 1 in ratio space, which Choro maps
            # to symmetric log10 limits about 0. Both backends read them, so
            # zmid alone would not do -- it centers the Plotly map only.
            "zmin": 10.0 ** -half,
            "zmax": 10.0 ** half,
            "zmid": 1.0,
            # Ticks at the positions Choro will actually send: log10 of the
            # data-unit limits, which round-trips ±half by an ulp for most
            # values -- a tick an ulp outside the range is a tick Plotly drops.
            "colorbar": _log_decade_colorbar(
                float(np.log10(10.0 ** -half)), float(np.log10(10.0 ** half)), unit="x"
            ),
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
    """Put a log map's ``tickvals``/``ticktext`` on the matplotlib colorbar.

    ``Choro.plot_mpl`` draws a log map through a ``LogNorm`` since 2026-09-13,
    so its colorbar is a real log axis in DATA units and already reads in
    decades; ``plot_mpl`` also applies caller ticks itself. This stays as the
    explicit spelling for the ``change`` map's ``x`` suffix: the ``tickvals``
    are the trace's log10 positions, so they are raised back to data units
    before being set. Applied only when the effective ``logscale`` is on, so an
    explicit ``logscale=False`` override does not get decade labels over linear
    values.
    """

    tickvals = (colorbar or {}).get("tickvals")
    ticktext = (colorbar or {}).get("ticktext")
    # axes[-1] is the colorbar geopandas appends; with colorbar=False there is
    # only the map axes and nothing to relabel.
    if not tickvals or not ticktext or len(figure.axes) < 2:
        return
    axis = figure.axes[-1]
    axis.set_yticks([10.0 ** float(tick) for tick in tickvals])
    axis.set_yticklabels([str(text) for text in ticktext])
    axis.minorticks_off()


def _parse_cell_list(value) -> list[int]:
    """Parse a snapshotted ``"0,1,2"`` cell list, tolerating an empty zone.

    `_resolve_drn_zone_cells` assigns ``[]`` to a zone that covers no active cell
    while other zones do (all of them missing raises, ledger 238), which the CSV
    snapshot round-trips as NaN (``str(nan) == "nan"``) while the GeoPackage
    snapshot round-trips as ``""``. Both mean the same thing.
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


#: The two forms PESTPP-IES writes an ensemble in: CSV, and PEST's binary ``.jcb``
#: under ``ies_save_binary`` (ledger 512).
_ENSEMBLE_SUFFIXES = (".csv", ".jcb")


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

        Parameters
        ----------
        bins : int, default 20
            The number of equal-width bins, shared by both histograms and spanning
            the pooled range of the prior and posterior samples, so the two are
            binned identically; each is normalized to a probability density. A
            positive integer: ``0`` and a float (``20.0`` included) raise from
            :func:`numpy.histogram_bin_edges`. Both branches also take anything
            else that function takes -- a rule name such as ``"auto"``, or
            explicit bin edges (the Plotly branch refused a rule name and drew
            edges with a meaningless width until ledger 480).
        title : str, optional
            The figure title. Default ``None``: ``"Forecast: <name>"``.
        backend : str, default 'plotly'
            ``"plotly"`` (or ``"interactive"``) returns an interactive :class:`viz.Fig`;
            ``"matplotlib"`` (or ``"mpl"``, ``"static"``) a static
            :class:`matplotlib.figure.Figure`. Case and surrounding spaces are
            ignored; anything else raises ``ValueError``.

        Returns
        -------
        viz.Fig or matplotlib.figure.Figure
            The overlaid histograms, with the truth line when :attr:`truth` is set.
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
        _histogram_bars(fig, prior, edges, density=True, name="prior", color=_PRIOR_COLOR)
        _histogram_bars(fig, posterior, edges, density=True, name="posterior", color=_POST_COLOR)
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
        Optional myflopy model used for spatial parameter maps (:meth:`plot_field`,
        :meth:`plot_obs_residuals`): its grid view (``model.vor``), else one built
        from its FloPy ``modelgrid`` -- a model declared with ``mf.dis(...)`` has no
        view (``grid.voronoi.model_grid_view``, ledger 451). Populated
        automatically when the run is launched via :meth:`PestProject.run_ies`.

    Notes
    -----
    Iteration ``0`` is the prior ensemble; the highest of THIS run's iterations --
    those its phi file records, :attr:`iterations` -- is treated as the posterior
    (pestpp-ies may stop before ``NOPTMAX`` if it converges). Use :meth:`obs_ensemble` / :meth:`par_ensemble` to reach any
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
        # Ensembles read so far, by (kind, path): (the file's (size, mtime_ns), ensemble).
        self._ensembles: dict[tuple[str, str], tuple[tuple[int, int], object]] = {}
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
        """Sorted iterations of THIS run with an observation ensemble on disk.

        An iteration counts when its ``<case>.<n>.obs.csv`` -- or ``.obs.jcb``, written
        under ``ies_save_binary`` -- exists AND the run's own phi file
        (``<case>.phi.actual.csv``, which every run rewrites from its first row) lists
        it. A file of an earlier run in the same directory, an iteration past this
        run's last, is not this run's (ledger 507): a rerun with fewer
        iterations reviewed the earlier run's last iteration as its posterior, and a
        prior Monte Carlo an IES posterior. Without a phi file every ensemble on
        disk counts.
        """

        nums = set()
        for suffix in _ENSEMBLE_SUFFIXES:
            for path in self.workspace.glob(f"{self.case}.*.obs{suffix}"):
                middle = path.name[len(self.case) + 1 : -len(f".obs{suffix}")]
                if middle.isdigit():
                    nums.add(int(middle))
        recorded = self._phi_iterations()
        if recorded is not None:
            nums &= recorded
        return sorted(nums)

    def _phi_iterations(self) -> set[int] | None:
        """The iterations the run's phi file records, or ``None`` when it has none."""

        path = self.workspace / f"{self.case}.phi.actual.csv"
        if not path.exists():
            return None
        try:
            column = pd.read_csv(path, usecols=["iteration"])["iteration"]
        except (OSError, ValueError) as error:
            # A phi file PESTPP-IES is still writing, or one with no `iteration`
            # column (pandas' EmptyDataError and ParserError are ValueErrors): no
            # record to read, so every ensemble on disk counts.
            logger.debug("phi file %s gives no iteration record: %s", path, error)
            return None
        return {int(value) for value in pd.to_numeric(column, errors="coerce").dropna()}

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

    def _option_is_set(self, name: str) -> bool:
        """Whether the run's control file sets the boolean PEST++ option ``name``."""

        value = self.pst.pestpp_options.get(name)
        return str(value).strip().lower() in ("true", "1", "yes", "t")

    def _ensemble_file(self, stem: str) -> Path | None:
        """The run's ensemble file ``<stem>.csv`` or ``<stem>.jcb``, or ``None``.

        Both forms are read (ledger 512); when both exist, the one the run's
        ``ies_save_binary`` option says it wrote.
        """

        suffixes = _ENSEMBLE_SUFFIXES
        if self._option_is_set("ies_save_binary"):
            suffixes = tuple(reversed(suffixes))
        for suffix in suffixes:
            path = self.workspace / f"{stem}{suffix}"
            if path.exists():
                return path
        return None

    def _read_ensemble(self, kind: str, path: Path):
        """Read one ensemble file, once while it is unchanged on disk.

        Cached by path and the file's size and modification time, so the review reads
        each file once however many views use it (ledger 513) -- a captured
        full-grid field is hundreds of MB, ~2 s a read -- and a file rewritten since
        (a run still going) is read again.
        """

        stat = path.stat()
        key = (kind, str(path))
        stamp = (stat.st_size, stat.st_mtime_ns)
        cached = self._ensembles.get(key)
        if cached is not None and cached[0] == stamp:
            return cached[1]
        cls = (
            self.pyemu.ObservationEnsemble if kind == "obs" else self.pyemu.ParameterEnsemble
        )
        if path.suffix == ".jcb":
            ensemble = cls.from_binary(pst=self.pst, filename=str(path))
        else:
            ensemble = cls.from_csv(pst=self.pst, filename=str(path))
        self._ensembles[key] = (stamp, ensemble)
        return ensemble

    def obs_ensemble(self, iteration: int):
        """Return the observation ensemble at ``iteration`` as a pyEMU object.

        Parameters
        ----------
        iteration : int
            The PESTPP-IES iteration whose ``<case>.<iteration>.obs.csv`` -- or
            ``.obs.jcb``, written under ``ies_save_binary`` (ledger 512) --
            to read: ``0`` is the prior, and :attr:`iterations` lists this run's. Passed
            through ``int()``, so ``"2"`` reads iteration 2 and a float is truncated
            toward zero (``2.9`` reads iteration 2). An iteration with no file raises
            ``FileNotFoundError`` naming the path. Not checked against
            :attr:`iterations`, so a file an earlier run left is read when asked for
            by number.

        Returns
        -------
        pyemu.ObservationEnsemble
            One row per realization (labelled as the file labels them: text when it
            holds a ``base`` row, integers when a CSV does not), one column per
            observation. The SAME object on every call while the file is unchanged
            (ledger 513): modify a ``.copy()``, never it.
        """

        stem = f"{self.case}.{int(iteration)}.obs"
        path = self._ensemble_file(stem)
        if path is None:
            raise FileNotFoundError(
                f"No observation ensemble at iteration {iteration}: "
                f"{self.workspace / stem}.csv (or .jcb)"
            )
        return self._read_ensemble("obs", path)

    def par_ensemble(self, iteration: int):
        """Return the parameter ensemble at ``iteration`` as a pyEMU object.

        Parameters
        ----------
        iteration : int
            The PESTPP-IES iteration whose ``<case>.<iteration>.par.csv`` -- or
            ``.par.jcb`` (ledger 512) -- to read: ``0`` is the prior. Passed
            through ``int()``, so ``"2"`` reads iteration 2 and a float is truncated
            toward zero. An iteration with no file raises ``FileNotFoundError`` naming
            the path -- including one :attr:`iterations` lists, when the run stopped
            between writing that iteration's observation and parameter ensembles.

        Returns
        -------
        pyemu.ParameterEnsemble
            One row per realization, one column per parameter, in the values
            PESTPP-IES wrote (not log-transformed). The SAME object on every call
            while the file is unchanged (ledger 513): modify a ``.copy()``.
        """

        stem = f"{self.case}.{int(iteration)}.par"
        path = self._ensemble_file(stem)
        if path is None:
            raise FileNotFoundError(
                f"No parameter ensemble at iteration {iteration}: "
                f"{self.workspace / stem}.csv (or .jcb)"
            )
        return self._read_ensemble("par", path)

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
        """The observation+noise ensemble, or ``None`` if the run used no noise.

        ``None`` whenever the run's control file sets ``ies_no_noise``, whatever an
        earlier run left in the directory (ledger 507); else
        ``<case>.obs+noise.csv`` or ``.jcb``.
        """

        if self._option_is_set("ies_no_noise"):
            return None
        path = self._ensemble_file(f"{self.case}.obs+noise")
        if path is None:
            return None
        return self._read_ensemble("obs", path)

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
        log : bool, default True
            Draw the phi axis on a log scale (recommended: phi typically falls by
            orders of magnitude). Tested for truth; ``False`` draws it linear.
        measured : bool, default False
            Plot the 'measured+noise' phi (:attr:`phi_measured`, each realization
            against its own noisy copy of the data, ``<case>.phi.meas.csv``)
            instead of the 'actual' phi (:attr:`phi`, against the observed values,
            ``<case>.phi.actual.csv``). The title says which.
        backend : str, default 'plotly'
            ``"plotly"`` (or ``"interactive"``) returns an interactive :class:`viz.Fig`;
            ``"matplotlib"`` (or ``"mpl"``, ``"static"``) a static
            :class:`matplotlib.figure.Figure`. Case and surrounding spaces are
            ignored; anything else raises ``ValueError``.

        Returns
        -------
        viz.Fig or matplotlib.figure.Figure
            One line per realization plus the bold ensemble mean, phi against
            iteration.
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

    # The simulated COLUMN an observation belongs to: pyEMU names a list-style
    # observation "<column base>_<index>:<value>" ("oname:hds_otype:lst_usecol:w1"
    # + "_per:0"), and that base is the group pyEMU gave every observation of the
    # column -- the only group a run built before ledger 358 has.
    _COLUMN_BASE_RE = re.compile(r"^(?P<base>.+)_[^_:]+:[^_:]+$")

    def _target_snapshot(self) -> pd.DataFrame:
        """Every target the run recorded by its observation: ``obsnme`` -> set, target, column.

        Read from each set's ``<prefix>_target_values.csv``, which names each
        target's observation since ledger 340 and its simulated ``column`` since
        ledger 358. Empty for a run older than both.
        """

        frames = []
        for entry in self.observation_sets:
            values_file = entry.get("values_file")
            path = self.workspace / str(values_file) if values_file else None
            if path is None or not path.exists():
                continue
            table = read_target_table(path)
            if "obsnme" not in table.columns:
                continue
            frames.append(pd.DataFrame({
                "obsnme": table["obsnme"].astype(str).str.strip().str.lower(),
                "prefix": str(entry.get("prefix", "")),
                "target": (
                    table["column"] if "column" in table.columns else table["name"]
                ).astype(str),
            }))
        if not frames:
            return pd.DataFrame(columns=["obsnme", "prefix", "target"])
        return pd.concat(frames, ignore_index=True).drop_duplicates("obsnme").set_index("obsnme")

    def _is_measured(self, obs: pd.DataFrame) -> pd.Series:
        """Which observations are measured, history-matched targets.

        Weight above zero, and -- when the run recorded what it targeted
        (:meth:`_measured_observation_keys`) -- one of those: a run built before
        ledger 340 left every row no target named at weight 1 and the base run's
        own value, which is no measurement.
        """

        weighted = obs["weight"].astype(float) > 0
        measured = self._measured_observation_keys()
        if measured is None:
            return weighted
        keys, names = measured
        parsed = obs.index.to_series().str.extract(self._OBS_NAME_RE)
        joined = [
            str(name).lower() in names
            or (str(prefix).lower(), str(location).lower(), str(time).lower()) in keys
            for name, prefix, location, time in zip(
                obs.index, parsed["prefix"], parsed["location"], parsed["time"], strict=True
            )
        ]
        return weighted & pd.Series(joined, index=obs.index)

    def _location_panels(self, groups, max_locations: int) -> list[dict]:
        """One time series per target LOCATION with a measured observation, in control-file order.

        A location is one simulated COLUMN (``_COLUMN_BASE_RE``): its every
        observation, measured or not, in time order, so the ensemble draws the
        whole series. Since ledger 358 a group holds every location measured
        under its label, so a group is no longer one location: drawn per group,
        a group of three wells was one line zig-zagging between them.

        Parameters
        ----------
        groups : sequence of str or None
            Observation groups to keep (any case): a location is kept when its
            measured observations are in one. ``None``: every location.
        max_locations : int
            The most panels returned, the first in control-file order.

        Returns
        -------
        list of dict
            ``title``, ``names`` (time-ordered), ``times`` and ``measured`` (a
            boolean array over ``names``) per location.
        """

        obs = self._observation_metadata()
        measured = self._is_measured(obs)
        base = obs.index.to_series().str.extract(self._COLUMN_BASE_RE)["base"]
        base = base.fillna(obs["obgnme"].astype(str))
        group = obs["obgnme"].astype(str).str.lower()
        wanted = None if groups is None else {str(name).strip().lower() for name in groups}
        snapshot = self._target_snapshot()
        panels = []
        for column in pd.unique(base[measured]):
            members = base.index[base == column]
            hits = members[measured.loc[members].to_numpy()]
            their_groups = list(dict.fromkeys(group.loc[hits]))
            if wanted is not None and not wanted.intersection(their_groups):
                continue
            known = snapshot.index.intersection(hits)
            if len(known):
                row = snapshot.loc[known[0]]
                title = f"{row['prefix']}: {row['target']}"
                if their_groups != [str(column).lower()]:
                    title += f" [{', '.join(their_groups)}]"
            else:
                title = str(column)
            frame = obs.loc[members].sort_values("_time", kind="stable")
            panels.append({
                "title": title,
                "names": frame.index.tolist(),
                "times": frame["_time"].to_numpy(),
                "measured": measured.loc[frame.index].to_numpy(dtype=bool),
            })
            if len(panels) >= int(max_locations):
                break
        return panels

    def plot_vs_obs(self, *, groups: list[str] | None = None, max_locations: int = 12,
                    iteration: int | None = None, backend: str = "plotly"):
        """Plot prior and posterior ensembles against the measured observations.

        One panel per target LOCATION -- a well, a reach, a lake, a DRN zone --
        that has a measured observation: grey lines are the prior ensemble and
        blue lines the posterior through every period the run saved, and red
        markers the measured values. A good fit means the blue lines bracket the
        red markers. Each panel is titled with its observation set's prefix, the
        target and, in brackets, its observation group.

        A panel was one observation GROUP until ledger 358, which was one
        location only because pyEMU gave every simulated column a group of its
        own; a target's group is now its PEST group, so a group of several wells
        drew one line through all of them. A run built before that -- one group
        per column -- draws the same panels, titled by target where its values
        snapshot names its observations (ledger 340) and by the group otherwise.
        Its markers are its MEASURED observations only: the panel drew every
        observation of its group as one, so each period nobody measured showed
        the base run's own value as data.

        Parameters
        ----------
        groups : list of str or None, optional
            Observation groups to draw (any case), named as the control file
            holds them -- ``north_wells`` for a label ``"North Wells"``, as
            :meth:`phi_contributions` lists them: the locations whose measured
            observations are in one of them. Default ``None``: every location
            with a measured observation, up to ``max_locations``.
        max_locations : int, optional
            The most panels drawn, the first locations in control-file order.
            Default 12.
        iteration : int or None, optional
            Posterior iteration to plot. Default ``None``: the highest available.
        backend : str, optional
            ``"plotly"`` (interactive, default) or ``"matplotlib"`` (static).

        Returns
        -------
        plotly.graph_objects.Figure or matplotlib.figure.Figure
            One row per location.

        Raises
        ------
        ValueError
            When no location has a measured observation (in ``groups``).
        """

        panels = self._location_panels(groups, max_locations)
        if not panels:
            raise ValueError(
                "No measured, nonzero-weight observation to plot"
                + (f" in groups {list(groups)}." if groups is not None else ".")
            )

        prior = self.prior._df
        posterior = self.obs_ensemble(iteration if iteration is not None else self.posterior_iteration)._df
        obsval = self.pst.observation_data["obsval"]
        titles = [panel["title"] for panel in panels]

        if _normalize_backend(backend) == "matplotlib":
            fig, axes = viz.mpl_axes(len(panels), 1, figsize=(8, 2.6 * len(panels)), squeeze=False)
            for ax, panel in zip(axes[:, 0], panels, strict=False):
                names, times, hit = panel["names"], panel["times"], panel["measured"]
                for real in prior.index:
                    ax.plot(times, prior.loc[real, names].to_numpy(dtype=float), color=_MPL_PRIOR, lw=0.8, alpha=0.4)
                for real in posterior.index:
                    ax.plot(times, posterior.loc[real, names].to_numpy(dtype=float), color=_MPL_POST, lw=0.8, alpha=0.5)
                ax.plot(times[hit], obsval.loc[names].to_numpy(dtype=float)[hit], "^", color=_MPL_MEAS, ms=7)
                ax.set_title(panel["title"], loc="left", fontsize=9)
                ax.set_ylabel("value")
            axes[-1, 0].set_xlabel("time")
            fig.suptitle("Simulated ensemble vs measured observations")
            fig.tight_layout()
            return fig

        fig = viz.subplots(rows=len(panels), cols=1, subplot_titles=titles, shared_xaxes=False)
        for row, panel in enumerate(panels, start=1):
            names, times, hit = panel["names"], panel["times"], panel["measured"]
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
            fig.add_scatter(x=times[hit], y=obsval.loc[names].to_numpy(dtype=float)[hit], mode="markers",
                            marker=dict(color=_MEAS_COLOR, size=8, symbol="triangle-up"), row=row, col=1,
                            name="measured", legendgroup="measured", showlegend=row == 1)
        fig.update_layout(title="Simulated ensemble vs measured observations",
                          dragmode="pan", height=260 * len(panels))
        return fig

    def forecast(self, name) -> IesForecast:
        """Return one forecast's prior+posterior ensemble for summary/plotting.

        Parameters
        ----------
        name : str or int
            A forecast observation name from :attr:`forecast_names`: matched
            exactly first, then as a case-insensitive SUBSTRING of those names,
            which must pick out exactly one (``"fore_1"`` finds
            ``"oname:fore1_otype:lst_usecol:fore_1_per:0"``). No match raises
            ``KeyError`` listing the forecasts, two or more raise ``KeyError``
            naming them. An integer is a POSITION in :attr:`forecast_names`
            (ledger 480; out of range raises ``IndexError``). Anything else
            -- the target object or spec that was passed to ``cal.forecast(...)``
            included -- raises ``TypeError``. Only names registered as forecasts
            are searched; a history-matched observation is not one.

        Returns
        -------
        IesForecast
            With ``.summary()`` and ``.plot()``; its ``name`` is the full
            observation name and its ``truth`` that observation's control-file
            value.
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
        if isinstance(name, Integral) and not isinstance(name, bool):
            if not -len(names) <= int(name) < len(names):
                raise IndexError(
                    f"forecast {name} is out of range: this run has {len(names)} forecasts."
                )
            return names[int(name)]
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
            "`forecast_names`) or its position in that list."
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
        measured : bool, default False
            Use the 'measured+noise' phi (:attr:`phi_measured`, each realization
            vs its noisy data copy) instead of the 'actual' phi (:attr:`phi`, vs
            the raw observed values).
        bins : int, default 25
            The number of equal-width bins over the pooled range of
            ``log10(phi)`` of both iterations, shared by the two histograms.
            Realizations whose phi is not finite (failed runs) are dropped first,
            and none left raises ``ValueError``. A positive integer: ``0`` and a
            float raise from :func:`numpy.histogram_bin_edges`. Both branches
            also take a rule name (``"auto"``) or explicit edges in
            ``log10(phi)`` (the Plotly branch refused a rule name until ledger
            480).
        backend : str, default 'plotly'
            ``"plotly"`` (or ``"interactive"``) returns an interactive :class:`viz.Fig`;
            ``"matplotlib"`` (or ``"mpl"``, ``"static"``) a static
            :class:`matplotlib.figure.Figure`. Case and surrounding spaces are
            ignored; anything else raises ``ValueError``.

        Returns
        -------
        viz.Fig or matplotlib.figure.Figure
            The prior (first row of the phi file) and posterior (last row)
            histograms of ``log10(phi)``, counted in realizations.
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

        fig = viz.Fig()
        _histogram_bars(fig, prior, edges, density=False, name="prior", color=_PRIOR_COLOR)
        _histogram_bars(fig, posterior, edges, density=False, name="posterior", color=_POST_COLOR)
        fig.update_layout(barmode="overlay", title=title, dragmode="pan",
                          xaxis_title="log10(phi)", yaxis_title="realizations")
        return fig

    def parameters_at_bounds(self, *, which: str = "posterior", tol: float = 0.01) -> pd.DataFrame:
        """Return, by parameter group, how many parameters are pinned at a bound.

        A parameter counts as "at a bound" when its ensemble mean, taken in the
        parameter's transform space, sits within ``tol`` (a fraction of the bound
        range, in that space) of its lower or upper bound. For a log-transformed
        parameter that is the mean of ``log10`` value -- a geometric mean -- as in
        :meth:`sensitivity`; until ledger 514 it was ``log10`` of the
        arithmetic mean, which a wide log-normal spread pulls toward the upper
        bound.

        Parameters
        ----------
        which : str, default 'posterior'
            Which parameter ensemble to average: ``"posterior"`` reads
            :attr:`posterior_parameters` (the highest of :attr:`iterations`),
            ``"prior"`` :attr:`prior_parameters`; any case. Anything else raises
            ``ValueError`` -- until ledger 480 it read the prior without a
            word.
        tol : float, default 0.01
            How close to a bound counts as "at" it, as a fraction of the
            parameter's bound range: a parameter whose ensemble mean sits at or
            below ``tol`` of the way from its lower bound counts at the lower, at
            or above ``1 - tol`` at the upper. Measured in ``log10`` for a
            ``partrans == "log"`` parameter whose bounds and every realization's
            value are positive, linearly otherwise. Not range-checked: ``tol >= 0.5`` counts
            every parameter at a bound, and a negative one counts none. A
            parameter whose bounds are equal is counted at neither.

        Returns
        -------
        pandas.DataFrame
            Indexed by parameter group with columns ``n_parameters``,
            ``n_at_lower``, ``n_at_upper`` and ``pct_at_bound``, most-pinned group
            first. Only parameters in both the ensemble and the control file count.
        """

        side = _which_ensemble(which)
        ensemble = (self.posterior_parameters if side == "posterior" else self.prior_parameters)._df
        par = self.pst.parameter_data.loc[[c for c in ensemble.columns if c in self.pst.parameter_data.index]]
        values = ensemble[par.index].astype(float)
        lower = par["parlbnd"].astype(float)
        upper = par["parubnd"].astype(float)
        is_log = par["partrans"].astype(str).str.lower().eq("log")

        value, low, high = values.mean(), lower.copy(), upper.copy()
        log_mask = is_log & (lower > 0) & (upper > 0) & (values > 0).all()
        # The mean in the space the parameter is estimated in: of log10 value for a
        # log parameter (ledger 514), as `sensitivity` measures it.
        value[log_mask] = np.log10(values.loc[:, log_mask]).mean()
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
        forecast : str or int, optional
            The observation to correlate every parameter against, across the
            POSTERIOR realizations (each paired by realization label, so a
            realization present in only one ensemble drops out). Any observation
            of the run by its full name (a column of the posterior observation
            ensemble, a history-matched one included), or a FORECAST as
            :meth:`forecast` reads it -- a case-insensitive substring picking out
            one, or a position (ledger 480; it took the exact full name
            only). A name matching neither raises ``KeyError`` listing
            :attr:`forecast_names`. Adds the ``forecast_corr`` column, and the
            table is then sorted by it. Default ``None``: the ``learned`` column
            only.
        by_group : bool, default True
            Aggregate to parameter groups (``pargp``): ``learned`` is then the
            group's mean, ``forecast_corr`` the mean of its members' ABSOLUTE
            correlations (so members pulling a forecast opposite ways do not
            cancel), and ``n_parameters`` its size. ``False`` returns one row per
            adjustable parameter, with its ``pargp`` and the SIGNED correlation,
            sorted by absolute correlation when ``forecast`` is given and by
            ``learned`` otherwise -- large for a grid or pilot-point run.

        Returns
        -------
        pandas.DataFrame
            Sorted by the strongest signal available, most informative first.
            Only adjustable parameters (not ``fixed`` or ``tied``) present in
            both the prior and the posterior ensemble count; with none the frame
            is empty.
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

    def _correlation_target(self, forecast) -> str:
        """An observation's full name: any column given exactly, else a forecast as
        :meth:`forecast` reads it -- the one resolver (it took the full name only
        until ledger 480)."""

        if isinstance(forecast, str) and forecast in self.posterior._df.columns:
            return forecast
        return self._resolve_forecast_name(forecast)

    def _forecast_realizations(self, forecast: str, index) -> pd.Series:
        """One forecast's posterior realizations, aligned to a parameter index."""

        values = self.posterior._df[self._correlation_target(forecast)].astype(float)
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
        forecast : str or int, optional
            An observation, read as :meth:`sensitivity` reads it (any full name,
            or a forecast's substring or position): the bars are each parameter
            group's mean absolute correlation with it. Default
            ``None``: the bars are each group's uncertainty reduction
            (``learned``).
        top : int, default 15
            Keep the strongest ``top`` groups, in the table's order: a whole
            number of at least 1. Anything else raises -- ``0`` raised pandas'
            ``cannot convert float NaN to integer`` and ``2.5`` was truncated
            until ledger 480.
        backend : str, default 'plotly'
            ``"plotly"`` (or ``"interactive"``) returns an interactive :class:`viz.Fig`;
            ``"matplotlib"`` (or ``"mpl"``, ``"static"``) a static
            :class:`matplotlib.figure.Figure`. Case and surrounding spaces are
            ignored; anything else raises ``ValueError``. Both draw the
            ``1/sqrt(n_reals)`` noise floor as a dashed line.

        Returns
        -------
        viz.Fig or matplotlib.figure.Figure
            One horizontal bar per parameter group. A run with no adjustable
            parameter in both ensembles raises ``ValueError``.
        """

        table = self.sensitivity(forecast=forecast)
        if table.empty:
            raise ValueError(
                "This run has no adjustable parameters in both ensembles, so "
                "there is nothing to score."
            )
        column = "forecast_corr" if forecast is not None else "learned"
        top = whole_number(top, what="group count")
        if top < 1:
            raise ValueError(f"top={top} keeps no parameter group; pass 1 or more.")
        if forecast is not None:
            forecast = self._correlation_target(forecast)
        table = table.head(top)
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
        which : str, default 'posterior'
            Passed to :meth:`parameters_at_bounds`: ``"posterior"`` reads
            :attr:`posterior_parameters`, ``"prior"`` :attr:`prior_parameters`;
            any case. Anything else raises ``ValueError`` (ledger 480).
        tol : float, default 0.01
            Passed to :meth:`parameters_at_bounds`: closeness to a bound, as a
            fraction of the bound range (``log10`` for log-transformed
            parameters), that counts as "at" the bound.
        backend : str, default 'plotly'
            ``"plotly"`` (or ``"interactive"``) returns an interactive :class:`viz.Fig`;
            ``"matplotlib"`` (or ``"mpl"``, ``"static"``) a static
            :class:`matplotlib.figure.Figure`. Case and surrounding spaces are
            ignored; anything else raises ``ValueError``.

        Returns
        -------
        viz.Fig or matplotlib.figure.Figure
            One horizontal bar per parameter group: the percentage of its
            parameters at either bound.
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

        Parameters
        ----------
        realization : str or int, default 'base'
            The label of a POSTERIOR realization (a row of :attr:`posterior`),
            compared as TEXT, so ``"3"`` and ``3`` are one label whichever the
            ensemble file uses. A label not there -- a typo, or ``"base"`` on a
            run with no ``base`` realization -- raises ``KeyError`` listing some
            labels; until ledger 480 it scored the ensemble's column MEAN
            in silence.

        Returns
        -------
        pandas.Series
            Indexed by observation group (``obgnme``), sorted from largest
            contribution down.
        """

        obs = self.pst.observation_data
        obs = obs.loc[obs["weight"].astype(float) > 0]
        simulated = self._posterior_row(realization)
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
        realization : str or int, default 'base'
            Passed to :meth:`phi_contributions`: the posterior realization to
            score; one not in the posterior ensemble raises ``KeyError``.
        kind : str, default 'bar'
            ``"bar"`` draws horizontal bars, ``"pie"`` a pie of the shares; any
            case. Anything else raises ``ValueError`` (a typo drew bars until
            ledger 480).
        max_groups : int, default 20
            The most groups shown, largest contributors first
            (:meth:`pandas.Series.head`): a float raises ``TypeError``, and a
            negative count keeps all but the last that many. On a pie the rest
            are left out, not pooled, so the shares are of the groups shown.
        backend : str, default 'plotly'
            ``"plotly"`` (or ``"interactive"``) returns an interactive :class:`viz.Fig`;
            ``"matplotlib"`` (or ``"mpl"``, ``"static"``) a static
            :class:`matplotlib.figure.Figure`. Case and surrounding spaces are
            ignored; anything else raises ``ValueError``.

        Returns
        -------
        viz.Fig or matplotlib.figure.Figure
            The phi contribution of each observation group shown.
        """

        kind = str(kind).strip().lower()
        if kind not in ("bar", "pie"):
            raise ValueError(f"kind={kind!r} is not a phi-contribution plot; use 'bar' or 'pie'.")
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
        coverage : float, default 1.0
            Central fraction of the prior ensemble treated as "reachable", per
            observation: the band runs from the ``(1 - coverage) / 2`` to the
            ``(1 + coverage) / 2`` quantile of the prior simulated values.
            ``1.0`` uses the full min-max range; ``0.95`` uses the 2.5-97.5
            percentile band (more robust to a few extreme realizations); ``0``
            collapses it to the median, so every observation not exactly on it
            is in conflict. Outside ``[0, 1]`` raises ``ValueError`` (a value from
            ``-1`` up to ``0`` swapped the band's ends, so every observation read
            in conflict, until ledger 480).

        Returns
        -------
        pandas.DataFrame
            Indexed by observation name with columns ``obgnme``, ``time``,
            ``measured``, ``prior_lo``, ``prior_hi``, ``prior_mean`` and
            ``in_conflict``.
        """

        if not 0.0 <= float(coverage) <= 1.0:
            raise ValueError(
                f"coverage={coverage!r} is the central FRACTION of the prior "
                "ensemble treated as reachable; it lies in [0, 1]."
            )
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

    def plot_prior_vs_obs(self, *, groups: list[str] | None = None, max_locations: int = 12,
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

        One panel per target LOCATION with a measured observation, as
        :meth:`plot_vs_obs` draws them (ledger 358); :meth:`plot_conflict`
        summarizes the same test by observation group.

        Parameters
        ----------
        groups : list of str or None, optional
            Observation groups to draw (any case), named as the control file
            holds them -- ``north_wells`` for a label ``"North Wells"``, as
            :meth:`phi_contributions` lists them: the locations whose measured
            observations are in one of them. Default ``None``: every location
            with a measured observation, up to ``max_locations``.
        max_locations : int, optional
            The most panels drawn, the first locations in control-file order.
            Default 12.
        show_conflict : bool, optional
            Mark the measured observations in prior-data conflict (orange ``x``).
            Default ``True``.
        coverage : float, optional
            Passed to :meth:`conflict`: the central fraction of the prior
            ensemble counted as reachable. Default ``1.0``, its full range.
        backend : str, optional
            ``"plotly"`` (interactive, default) or ``"matplotlib"`` (static).

        Returns
        -------
        plotly.graph_objects.Figure or matplotlib.figure.Figure
            One row per location.

        Raises
        ------
        ValueError
            When no location has a measured observation (in ``groups``).
        """

        panels = self._location_panels(groups, max_locations)
        if not panels:
            raise ValueError(
                "No measured, nonzero-weight observation to plot"
                + (f" in groups {list(groups)}." if groups is not None else ".")
            )
        prior = self.prior._df
        obsval = self.pst.observation_data["obsval"]
        conflict = self.conflict(coverage=coverage) if show_conflict else None

        def _panel_data(panel):
            """The time-ordered names/times, measured markers and conflict flags of one location."""

            names, hit = panel["names"], panel["measured"]
            measured = obsval.loc[names].to_numpy(dtype=float)[hit]
            times = panel["times"][hit]
            if conflict is not None:
                kept = [name for name, flag in zip(names, hit, strict=True) if flag]
                flags = conflict["in_conflict"].reindex(kept).fillna(False).to_numpy(dtype=bool)
            else:
                flags = np.zeros(len(times), dtype=bool)
            return names, panel["times"], times, measured, flags

        if _normalize_backend(backend) == "matplotlib":
            fig, axes = viz.mpl_axes(len(panels), 1, figsize=(8, 2.6 * len(panels)), squeeze=False)
            for ax, panel in zip(axes[:, 0], panels, strict=False):
                names, line_times, times, measured, flags = _panel_data(panel)
                for real in prior.index:
                    ax.plot(line_times, prior.loc[real, names].to_numpy(dtype=float), color=_MPL_PRIOR, lw=0.8, alpha=0.4)
                ax.plot(times[~flags], measured[~flags], "^", color=_MPL_MEAS, ms=7, label="measured")
                if flags.any():
                    ax.plot(times[flags], measured[flags], "X", color=_MPL_CONFLICT, ms=10, label="prior-data conflict")
                ax.set_title(panel["title"], loc="left", fontsize=9)
                ax.set_ylabel("value")
            axes[-1, 0].set_xlabel("time")
            fig.suptitle("Prior ensemble vs measured observations")
            fig.tight_layout()
            return fig

        fig = viz.subplots(rows=len(panels), cols=1, subplot_titles=[panel["title"] for panel in panels],
                           shared_xaxes=False)
        for row, panel in enumerate(panels, start=1):
            names, line_times, times, measured, flags = _panel_data(panel)
            first = True
            for real in prior.index:
                fig.add_scatter(x=line_times, y=prior.loc[real, names].to_numpy(dtype=float), mode="lines",
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
                          dragmode="pan", height=260 * len(panels))
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
        coverage : float, default 1.0
            Passed to :meth:`conflict`: the central fraction of the prior
            ensemble counted as reachable, ``1.0`` its full min-max range and
            ``0.95`` the 2.5-97.5 percentile band. Above ``1`` raises.
        backend : str, default 'plotly'
            ``"plotly"`` (or ``"interactive"``) returns an interactive :class:`viz.Fig`;
            ``"matplotlib"`` (or ``"mpl"``, ``"static"``) a static
            :class:`matplotlib.figure.Figure`. Case and surrounding spaces are
            ignored; anything else raises ``ValueError``.

        Returns
        -------
        viz.Fig or matplotlib.figure.Figure
            One horizontal bar per observation group: the percentage of its
            nonzero-weight observations in conflict.
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
        ncpl = len(model_cell_polygons(self.model))

        def captured(token: str) -> dict[str, tuple[int, int]]:
            """Each captured observation of ``token``'s file, by its ``(i, j)``."""

            found = {}
            for name in self.pst.observation_data.index:
                if token not in name:
                    continue
                match = self._ARR_RE.search(name)
                if match:
                    found[name] = (int(match.group(1)), int(match.group(2)))
            return found

        def flat_index(found: dict[str, tuple[int, int]]) -> dict[str, int]:
            """Each observation's position in its file, read row by row.

            pyEMU names a captured array value ``arr_i:<line>_j:<value on the
            line>``, so the position is ``i * <values per line> + j``, and the
            values per line are read from the names themselves: one per line on
            DISV, and on DIS either a grid ROW (``ncol``, as the build writes an
            array any parameter scales -- a pilot-point array too since ledger
            525) or one (a pilot-point array of a template built
            before it, whose interpolation wrote one value per line). Reading ``i`` alone put a DIS row's values on cell ``row``; the
            grid's ``ncol`` put a pilot-point layer's 30 values on cells 0..174
            (ledger 420).
            """

            per_line = 1 + max((j for _, j in found.values()), default=0)
            return {name: i * per_line + j for name, (i, j) in found.items()}

        # JSON round-trips the int layer keys to strings.
        layer_prefixes = {int(k): v for k, v in (info.get("layer_prefixes") or {}).items()}
        mapping: dict[str, int] = {}
        if layer_prefixes:
            # One capture file per layer, each holding that layer's cells in
            # order, so the position in the file *is* the cell id.
            prefix = layer_prefixes.get(int(layer))
            if prefix is None:
                raise ValueError(
                    f"No captured field for layer {layer}. "
                    f"Captured layers: {sorted(layer_prefixes)}."
                )
            mapping = flat_index(captured(f"oname:{prefix.lower()}_"))
        else:
            # A single whole-grid array covering every layer at once, which is
            # how MF6 writes griddata with no LAYERED keyword (MST porosity):
            # flat = layer * ncpl + cell.
            for name, flat in flat_index(captured(f"oname:{info['prefix'].lower()}_")).items():
                if flat // ncpl == int(layer):
                    mapping[name] = flat % ncpl
        if not mapping:
            raise ValueError(f"No captured field cells found for layer {layer}.")
        return mapping

    def field(self, target: str, *, layer: int = 0) -> pd.DataFrame:
        """Return per-cell prior/posterior statistics for a captured field.

        Requires the parameter to have been declared with
        ``cal.parameterize(..., capture=True)``.

        Parameters
        ----------
        target : str
            The captured field, matched case-insensitively (surrounding spaces
            ignored) against each entry of :attr:`capture_fields` by its
            ``target`` (``"k"``) or its ``name``. A run that captured nothing
            raises ``ValueError`` naming ``capture=True``, an unknown name
            ``KeyError`` listing the captured targets, and a captured LIST field
            ``NotImplementedError`` (only array fields -- K, K33, porosity -- map).
        layer : int, default 0
            The zero-based model layer whose cells to return. A field captured per
            layer must have captured this one, or ``ValueError`` names the layers
            it did; a field captured as one whole-grid array (MST porosity) is
            split by ``layer * ncpl + cell``. Needs :attr:`model` for the grid's
            ``ncpl`` -- ``ValueError`` without one -- and a layer with no captured
            cell raises ``ValueError``.

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
            ``ncpl``. ``frame.attrs["space"]`` names the space the statistics
            were taken in (ledger 514):

            * ``"log10"`` -- when every parameter of the field is log-transformed
              (``transform="log"``, :meth:`PestProject.parameterize`'s default) and
              every captured value of both ensembles is positive. The space
              PESTPP-IES estimates those parameters in, and :meth:`sensitivity`'s:
              each mean is GEOMETRIC (``10 ** mean(log10 value)``), each std is the
              std of ``log10`` value (in decades), so ``change`` is the ratio of
              geometric means and ``reduction`` agrees with :meth:`sensitivity`'s
              ``learned``. Arithmetic statistics of a log-normal K read an unchanged
              median as roughly halved wherever the data narrowed the spread.
            * ``"linear"`` -- otherwise (a ``transform="none"`` parameter of the
              field, no parameter of the field found in the control file, or a
              value at or below zero): arithmetic means and stds, in model units.

            The parameters of a field are those of its parameter set's group:
            ``pargp`` ``<name>`` or ``<name>l<layer>...``, ``<name>`` being the
            captured field's parameter-set name.

        Notes
        -----
        Capturing a field records EVERY cell of every layer it covers as a
        zero-weight observation: 101,615 observations for ``k33`` over Ten Trails'
        5 x 20,323 cells, which grew its template from 18 to 33 MB and puts that
        many columns in every observation ensemble each iteration writes (a 352 MB
        CSV for one 5-layer field at 200 realizations, ~2 s a read). A capture covers
        the layers its parameter does (``layers=``); save the ensembles in binary
        (``run_ies(ies_save_binary=True)``), which the review reads too.
        """

        info = self._capture_info(target)
        if info.get("family") != "array":
            raise NotImplementedError(
                f"Spatial maps currently support array fields (K, K33); "
                f"{info.get('target')!r} is a list field."
            )
        cells = self._field_cell_index(info, layer)
        names = list(cells.keys())
        prior = self.prior._df.loc[:, names].astype(float)
        posterior = self.posterior._df.loc[:, names].astype(float)
        space = self._field_space(info, prior, posterior)
        if space == "log10":
            # Statistics in the space the parameters are estimated in: a geometric
            # mean and a spread in decades (ledger 514).
            logged_prior, logged_posterior = np.log10(prior), np.log10(posterior)
            means = (10.0 ** logged_prior.mean(), 10.0 ** logged_posterior.mean())
            spreads = (logged_prior.std(), logged_posterior.std())
        else:
            means = (prior.mean(), posterior.mean())
            spreads = (prior.std(), posterior.std())
        frame = pd.DataFrame({
            "cell": [cells[name] for name in names],
            "prior_mean": means[0].to_numpy(),
            "prior_std": spreads[0].to_numpy(),
            "posterior_mean": means[1].to_numpy(),
            "posterior_std": spreads[1].to_numpy(),
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
        frame = frame.sort_values("cell").reset_index(drop=True)
        frame.attrs["space"] = space
        return frame

    def _field_space(self, info: dict, prior: pd.DataFrame, posterior: pd.DataFrame) -> str:
        """``"log10"`` or ``"linear"``: the space :meth:`field` takes a field's statistics in.

        ``"log10"`` when the control file holds parameters of the field -- group
        ``<name>`` or ``<name>l<layer>...`` -- all ``partrans == "log"``, and every
        value of both ensembles is positive; the rule :meth:`sensitivity` applies per
        parameter (ledger 514).
        """

        name = str(info.get("name") or info.get("target") or "")
        data = getattr(self.pst, "parameter_data", None)
        if not name or data is None or "pargp" not in data:
            return "linear"
        pattern = re.compile(rf"{re.escape(name.lower())}(l\d+.*)?")
        groups = data["pargp"].astype(str).str.lower()
        transforms = data.loc[groups.map(lambda group: bool(pattern.fullmatch(group))).to_numpy(),
                              "partrans"].astype(str).str.lower()
        if transforms.empty or not transforms.eq("log").all():
            return "linear"
        values = np.concatenate([prior.to_numpy(dtype=float).ravel(),
                                 posterior.to_numpy(dtype=float).ravel()])
        return "log10" if values.size and bool((values > 0).all()) else "linear"

    def _field_choro(self, target: str, *, stat: str, which: str, layer: int,
                     title: str | None = None, **choropleth_kwargs):
        """Build one field map as a ``Choro``, with its heading and resolved kwargs.

        The seam :meth:`plot_field` and :meth:`plot_field_mosaic` share: a mosaic
        needs the ``Choro`` itself (``viz.mosaic`` reads ``get_choropleth()``),
        while ``plot_field`` renders it to a figure. Returns
        ``(choro, heading, kwargs)`` -- the kwargs so the caller can see which
        policy keys actually took effect after ``setdefault``.
        """

        which = _which_ensemble(which)
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
        ncpl = int(model_grid_view(self.model).ncpl)
        cells = frame["cell"].to_numpy(dtype=int)

        def _scatter(name: str) -> np.ndarray:
            """One captured column spread onto the full grid, uncaptured cells NaN."""

            # field() has one row per CAPTURED cell; every hover column must be
            # exactly ncpl long or the hover assembler raises.
            spread = np.full(ncpl, np.nan)
            spread[cells] = frame[name].to_numpy(dtype=float)
            return spread

        values = _scatter(column)
        # Name the space the statistics were taken in (ledger 514).
        labels = _LOG_STAT_LABELS if frame.attrs.get("space") == "log10" else _STAT_LABELS
        label = f"{target} {labels.get(stat_key, stat_key)}" + (
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
        heading = f"{label} — {context}" if title is None else str(title)

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
        choro = _choropleth_factory(
            model_grid_view(self.model), custom_zs=list(values), layer=layer,
            **choropleth_kwargs
        )
        return choro, heading, choropleth_kwargs

    def plot_field(self, target: str, *, stat: str = "mean", which: str = "posterior",
                   layer: int = 0, backend: str = "plotly", title: str | None = None,
                   **choropleth_kwargs):
        """Map a captured parameter field on the model grid (a choropleth of its cells).

        This answers "property patterns -- plausible or laughable?": it shows the
        spatial pattern of a calibrated property and how history matching changed
        it. Requires ``cal.parameterize(..., capture=True)``.

        Parameters
        ----------
        target : str
            A captured target, e.g. ``"k"``, resolved as :meth:`field` resolves
            it (case-insensitive, by ``target`` or ``name``); it also heads the
            map's title.
        stat : str, default 'mean'
            ``"mean"`` (ensemble mean field), ``"std"`` (spread -- where the
            property is still uncertain), ``"base"`` (the minimum-error-variance
            realization), ``"change"`` (posterior_mean / prior_mean -- where
            calibration moved the property), or ``"reduction"``
            (``1 - posterior_std / prior_std`` -- how much of the prior spread
            the data removed, i.e. *did the data inform this region*). Any case.
            Each has its own colour policy: ``mean`` and ``base`` on a log scale,
            ``std`` linear (a collapsed spread is a real zero), ``change`` on a
            log, diverging scale symmetric about a ratio of 1, and ``reduction``
            linear, pinned to 0..1 whenever every value lies in it. Anything else,
            and ``"base"`` on a run with no ``base`` realization, raises
            ``ValueError``.
        which : str, default 'posterior'
            ``"prior"`` or ``"posterior"``: whose ensemble ``stat`` in
            {``"mean"``, ``"std"``} summarizes. Ignored by the others: ``change``
            and ``reduction`` are already prior-to-posterior comparisons, and
            ``base`` is a single posterior realization, so it changes nothing
            with them. Any case; anything else raises ``ValueError`` naming the
            two (``"Posterior"`` raised "Unknown stat" until ledger 480).
        layer : int, default 0
            The zero-based model layer to map: which cells :meth:`field` reads,
            and the layer the choropleth is built on. See :meth:`field` for the
            layers a capture offers.
        backend : str, default 'plotly'
            ``"plotly"`` (or ``"interactive"``) returns the interactive figure
            (``Choro.fig``, a :class:`viz.Fig`);
            ``"matplotlib"`` (or ``"mpl"``, ``"static"``) a static
            :class:`matplotlib.figure.Figure` of the model cells
            (``Choro.plot_mpl``), with the same colours, limits and log-decade
            labels. Case and surrounding spaces are ignored; anything else raises
            ``ValueError``.
        title : str, optional
            The figure heading, on both backends and in the hover. Default
            ``None``: ``"<target> <stat> (<which>) — iteration <n>, <count>
            realizations"``. Until ledger 480 it was not a parameter, so it
            reached the trace (Plotly refused it) or was ignored (Matplotlib).
        **choropleth_kwargs
            Forwarded to the choropleth builder
            (``grid.plotting._choropleth_factory``, then ``Choro``): its named
            options -- ``colorscale``, ``logscale``, ``zmin``/``zmax``,
            ``features``, ``select``, ``extent``, ``contours``, ``hover_spec``,
            ``custom_hover``, ... -- and, past them, properties of the
            ``go.Choroplethmap`` trace, which Plotly validates when the figure is
            built (an unknown name raises ``ValueError``). Each overrides the house
            colour policy and hover this method would set; ``logscale=False`` also
            drops the policy's decade colorbar. ``custom_zs`` is set here and
            raises ``TypeError`` when passed (``layer=`` is this method's own
            parameter, above, so it never reaches the tail).

        See Also
        --------
        plot_field_mosaic : the same map for prior and posterior, side by side.
        """

        choro, heading, resolved = self._field_choro(
            target, stat=stat, which=which, layer=layer, title=title, **choropleth_kwargs
        )
        if _normalize_backend(backend) == "matplotlib":
            # No cmap=: plot_mpl builds it from the policy colorscale (stops and
            # names alike become a LinearSegmentedColormap from Plotly's own
            # stops) and vmin/vmax default to the policy's zmin/zmax, so the
            # colors and the limits agree between backends. plot_mpl also puts
            # the policy's own ticks on its log colorbar; the relabel below is
            # the explicit spelling of that same step, and idempotent over it.
            figure = choro.plot_mpl(title=heading)
            if resolved.get("logscale"):
                _relabel_log_colorbar(figure, resolved.get("colorbar"))
            return figure
        # `title` is matplotlib-only and is not a Choroplethmap property, so it
        # must never ride the trace kwargs. plot() sets a 20px top margin; bump it
        # or the title clips.
        return choro.fig.update_layout(title=heading, margin={"t": 40})

    def plot_field_mosaic(self, target: str, *,
                          which: Sequence[str] = ("prior", "posterior"),
                          stat: str = "mean", layer: int = 0, ncols: int = 2,
                          title: str | None = None, sync_views: bool | None = None,
                          backend: str = "plotly", **choropleth_kwargs):
        """Compose one field map per ensemble, side by side on a shared color scale.

        The prior-vs-posterior figure: the same property, the same colors, the
        same limits, so the panels are actually comparable and the eye can only
        be reading a real difference.

        Parameters
        ----------
        target : str
            A captured target, e.g. ``"k"``, resolved as :meth:`field` resolves
            it (case-insensitive, by ``target`` or ``name``).
        which : sequence of str, default ('prior', 'posterior')
            The ensembles to compose, one panel each and in this order:
            ``"prior"`` and/or ``"posterior"``, any case, a repeat allowed. Fewer
            than two raises ``ValueError``, and an unknown entry ``ValueError``
            naming the two. A bare string raises ``TypeError`` (it was iterated
            character by character until ledger 480).
        stat : str, default 'mean'
            ``"mean"`` or ``"std"`` (any case) -- the only stats that HAVE a prior
            and a posterior form. ``"change"`` and ``"reduction"`` are already
            prior-to-posterior comparisons, and ``"base"`` is a single posterior
            realization with no prior counterpart; all three, and anything else,
            are refused here with ``ValueError`` naming :meth:`plot_field`.
        layer : int, default 0
            The zero-based model layer every panel maps, as in :meth:`plot_field`.
        ncols : int, default 2
            Panels per row of the mosaic (:func:`myflopy.viz.mosaic`), capped at
            the number of panels; ``0`` gives one column.
        title : str, optional
            The overall figure title. Default ``None``:
            ``"<target> <stat> — prior vs posterior"`` (the ``which`` sides joined
            by ``vs``). Each panel keeps its own heading naming its iteration and
            realization count.
        sync_views : bool, optional
            Whether the panels pan and zoom together, as
            :func:`myflopy.viz.mosaic` takes it: ``None`` (default) links panels
            that open on one window, ``False`` leaves each independent, and
            ``True`` with a LIST ``extent`` of different windows raises
            ``ValueError``.
        backend : str, default 'plotly'
            ``"plotly"`` (or ``"interactive"``; case and surrounding spaces ignored) only
            -- ``viz.mosaic`` composes Plotly subplots. A static spelling
            (``"matplotlib"``, ``"mpl"``, ``"static"``) raises
            ``ValueError`` naming :meth:`plot_field`, which draws one static
            panel; anything else raises ``ValueError`` naming the two backends.
        **choropleth_kwargs
            ``extent`` goes to the mosaic: one window for every panel, or a list
            with one per panel (``None`` for a panel's own), VIEW only. Everything
            else is passed to EACH panel's choropleth exactly as
            :meth:`plot_field`'s ``**choropleth_kwargs`` are -- the ``Choro``
            options, then ``go.Choroplethmap`` trace properties Plotly validates
            -- overriding the house colour policy per panel, with ``custom_zs``
            refused (``TypeError``; ``layer=`` is this method's own parameter).
            The panels' own colorbars are discarded: the mosaic pools them onto
            one colour axis over every panel's data, a ``zmin``/``zmax`` passed
            here holding its side.

        Notes
        -----
        A mosaic pools every panel onto ONE color axis, which discards each
        panel's own colorbar (compromise ledger 71) and keeps the WIDEST of any
        limits the panels set (ledger 306). That pooling is the point -- it is
        what makes the panels comparable -- but it means the real-unit ticks for
        a log-scaled ``mean`` mosaic have to be rebuilt from the pooled limits,
        which is what the ``colorbar`` callback does.
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
        if isinstance(which, str):
            raise TypeError(
                f"which={which!r} is one ensemble; a mosaic composes several -- "
                "pass a sequence, e.g. which=('prior', 'posterior')."
            )
        sides = [_which_ensemble(side) for side in which]
        if len(sides) < 2:
            raise ValueError(
                "A mosaic needs at least two panels, e.g. "
                "which=('prior', 'posterior')."
            )

        extent = choropleth_kwargs.pop("extent", None)  # the composer's (ledger 311)
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
            extent=extent,
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
        """One realization's simulated values (an absent label raises, ledger 480)."""

        return self._posterior_row(realization)

    def _posterior_row(self, realization) -> pd.Series:
        """The posterior ensemble's row ``realization``, its label compared as text.

        The one realization lookup: ``phi_contributions``, ``obs_residuals`` and the
        residual map each fell back to the ensemble MEAN for a label not there,
        saying nothing (ledger 480).
        """

        ensemble = self.posterior._df
        by_text = {str(label): label for label in ensemble.index}
        key = str(realization)
        if key not in by_text:
            shown = list(by_text)[:8]
            more = " ..." if len(by_text) > 8 else ""
            hint = (
                " This run has no 'base' realization: pass one of its labels, or "
                "best(criterion='min_phi')."
                if key == "base" else ""
            )
            raise KeyError(
                f"No posterior realization {realization!r}; it has {shown}{more}.{hint}"
            )
        return ensemble.loc[by_text[key]]

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
                else read_target_table(path)
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
                    mapping = read_target_table(mapping_path)
                    # The label when the mapping carries one (ledger 234): pyEMU's
                    # observation names come from that, not from the bare name.
                    label = "obsname" if "obsname" in mapping.columns else "name"
                    mapping["location"] = (
                        mapping[label].astype(str).str.strip().str.lower()
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

    def _measured_observation_keys(
        self,
    ) -> tuple[set[tuple[str, str, str]], set[str]] | None:
        """The observations that carry a TARGET, by name or by ``(prefix, location, time)``.

        pyEMU creates one observation per row of the simulated output, not per
        row of the target table. A run built before ledger 340 left every row no
        target named at pyEMU's defaults -- **weight 1.0 and an ``obsval`` equal
        to the base model's own simulated value** -- and averaging those into a
        location's residual pulls it toward zero and can flip its sign, so they
        are excluded here by joining back to the ``<prefix>_target_values.csv``
        snapshot of what was targeted.

        A run built since carries those rows at weight 0, and its snapshot names
        each target's OBSERVATION (``obsnme``), which is joined on directly. The
        older snapshots are joined on ``str(time)`` against the name's trailing
        token, which misses a key pyEMU spelled differently (``1.0`` against
        ``per:1``) -- the spelling the build itself used to miss.

        Returns
        -------
        tuple of (set, set) or None
            ``(keys, names)``: ``(prefix, location, time)`` triples from older
            snapshots and lower-cased observation names from newer ones.
            ``None`` when no observation set recorded a values file, meaning
            "cannot tell" -- the caller then does not filter, because dropping
            everything would be worse than the bias.
        """

        keys: set[tuple[str, str, str]] = set()
        names: set[str] = set()
        found = False
        for entry in self.observation_sets:
            values_file = entry.get("values_file")
            prefix = str(entry.get("prefix", ""))
            if not values_file or not prefix:
                continue
            path = self.workspace / str(values_file)
            if not path.exists():
                continue
            table = read_target_table(path)
            if "obsnme" in table.columns:
                found = True
                names.update(table["obsnme"].dropna().astype(str).str.strip().str.lower())
                continue
            if not {"time", "name"}.issubset(table.columns):
                continue
            found = True
            # An older snapshot, from a build that rebuilt each target's name as
            # f"{index}:{str(time)}" (before ledger 340): comparing the stringified
            # time against the name's trailing token avoids having to know whether
            # this family indexed on `per` (heads) or `time` (named series).
            for name, time in zip(table["name"], table["time"], strict=True):
                keys.add((
                    prefix.lower(),
                    str(name).strip().lower(),
                    str(time).strip().lower(),
                ))
        return (keys, names) if found else None

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
        output -- at weight 0 in a run built since ledger 340, at weight 1.0
        before it; those are excluded by joining back to the
        ``*_target_values.csv`` snapshot. Zero-weight observations are excluded
        too, which is what keeps ``cal.forecast(...)`` sets — registered as
        ordinary observation sets — off a misfit figure.

        Parameters
        ----------
        realization : str, default 'base'
            Which posterior realization to score (default ``"base"``, the
            minimum-error-variance one), its label compared as text, as
            :meth:`phi_contributions` reads it; an absent label raises
            ``KeyError`` (it scored the ensemble MEAN until ledger 480).
        prefix : str or iterable of str, optional
            One observation-set prefix (as given to ``cal.observe``), or several,
            to restrict the frame to; matched case-insensitively. A prefix that
            names no LOCATABLE set (lake and SFR sets record no location) raises
            ``ValueError`` listing the ones recorded -- unless the run recorded
            none, when the frame is simply empty. Default ``None`` keeps every
            locatable set — which mixes UNITS when a run history-matched more
            than one kind of measurement (heads in length, concentration in
            mass/volume, DRN seepage in volume/time).

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

        # Only the measured targets (`_is_measured`): a row no target named, and a
        # zero weight -- forecasts, which `cal.forecast` registers as ordinary
        # observation sets so they arrive here looking like head targets, and
        # deliberately silenced targets. A misfit map must not draw either as if
        # the model had been fitted to it.
        obs = obs.loc[self._is_measured(obs)]

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
        realization : str, default 'base'
            Passed to :meth:`obs_residuals`: the posterior realization to score;
            an absent label raises ``KeyError``.
        layer : int, default 0
            The zero-based layer the cell choropleth is built for (``Choro``'s
            ``layer``). Nothing is FILTERED by it: head targets
            are drawn at their map position whatever layer they were screened in
            -- filtering them would hide data on a plan-view figure -- and a DRN
            zone colours the cells its in-layer numbers name, whatever layer the
            zone lies in. A point target's hover names its layer; read it there
            rather than off the map.
        backend : str, default 'plotly'
            ``"plotly"`` (or ``"interactive"``) returns the interactive figure
            (``Choro.fig``, a :class:`viz.Fig`) with the point targets as map
            markers; ``"matplotlib"`` (or ``"mpl"``, ``"static"``) a
            static :class:`matplotlib.figure.Figure` in MODEL coordinates with
            the points scattered over the cells. Case and surrounding spaces are
            ignored; anything else raises ``ValueError``.
        title : str, optional
            The figure heading, on both backends and in the hover. Default
            ``None``: ``"residuals (simulated − measured) — iteration <n>,
            <realization>, <count> locations"``.
        prefix : str or iterable of str, optional
            Passed to :meth:`obs_residuals`: one observation-set prefix, or
            several, to draw (case-insensitive; one naming no locatable set
            raises ``ValueError``). Default ``None`` draws every locatable set.
            Whatever is drawn, the method WARNS when the sets measure different
            quantities, since they then share one colour scale.
        **choropleth_kwargs
            Forwarded to the choropleth builder
            (``grid.plotting._choropleth_factory``, then ``Choro``): its named
            options and, past them, ``go.Choroplethmap`` trace properties that
            Plotly validates when the figure is built. ``colorscale``, ``zmin``,
            ``zmax``, ``zmid``, ``custom_hover``, ``hover_spec``,
            ``hover_heads``, ``hover_ks`` -- and ``showscale`` when no zone is
            drawn -- override what this method sets for the CELLS, and the point
            markers take the same ``colorscale``/``zmin``/``zmax`` (they kept the
            default scale until ledger 480); ``custom_zs`` raises
            ``TypeError`` (``layer=`` is this method's own parameter, so it never
            reaches the tail).

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

        ncpl = int(model_grid_view(self.model).ncpl)
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

        heading = title or (
            f"residuals (simulated − measured) — iteration "
            f"{self.posterior_iteration}, {realization}, {len(frame)} locations"
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

        choro = _choropleth_factory(
            model_grid_view(self.model), custom_zs=list(zone_values), layer=layer,
            **choropleth_kwargs
        )
        # The markers take the scale the CELLS resolved -- a caller's colorscale/
        # zmin/zmax included -- so a point and the zone under it read alike (they
        # kept the default scale until ledger 480).
        colorscale = choropleth_kwargs["colorscale"]
        low, high = choropleth_kwargs["zmin"], choropleth_kwargs["zmax"]

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
                    cmap=mpl_colormap_for(colorscale), vmin=low, vmax=high,
                    s=45, edgecolor="black", linewidth=0.5, zorder=5,
                )
            return figure

        if not points.empty:
            lon, lat = model_grid_view(self.model).points_to_latlon(
                points["x"].to_numpy(dtype=float), points["y"].to_numpy(dtype=float)
            )
            choro.add_overlay(go.Scattermap(
                lon=lon, lat=lat, mode="markers", name="observations",
                marker={
                    "color": points["residual"].to_numpy(dtype=float),
                    "colorscale": colorscale, "cmin": low, "cmax": high,
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
        return choro.fig.update_layout(title=heading, margin={"t": 40})

    def best(self, *, criterion: str = "base") -> str:
        """Return the label of the single 'best' realization to carry forward.

        Parameters
        ----------
        criterion : str, default 'base'
            ``"base"`` (default, recommended) returns the ``base`` realization --
            the minimum-error-variance parameter set, the closest analogue to a
            deterministic GLM calibration. ``"min_phi"`` returns the lowest-phi
            realization in the LAST row of the 'actual' phi file among those
            PESTPP-IES kept: one with an empty phi cell (dropped or failed) or
            absent from the posterior observation or parameter ensemble is never
            returned (ledger 515), and ``ValueError`` is raised when no
            realization is left; this is tempting but discouraged (it tends to be over-fit and noisy), so
            prefer ``"base"`` unless you have a reason. Any case. A run with no
            ``base`` realization raises ``ValueError`` naming ``"min_phi"`` (it
            answered with the ``"min_phi"`` choice in silence until ledger
            480). Anything else raises ``ValueError``.

        Returns
        -------
        str
            The realization's label, as the phi file's column names it.
        """

        criterion = str(criterion).lower()
        phi = self.phi
        realization_cols = list(phi.columns[6:])
        if criterion == "base":
            if "base" in realization_cols:
                return "base"
            raise ValueError(
                "This run has no 'base' realization (PESTPP-IES writes one unless "
                "ies_include_base is off); pass criterion='min_phi' for the lowest-"
                "phi realization instead."
            )
        if criterion == "min_phi":
            # PESTPP-IES writes an EMPTY phi cell for a realization it dropped or that
            # failed, and no number compares below NaN, so `min` kept a leading one
            # and returned a realization with no posterior (ledger 515).
            final = pd.to_numeric(phi.iloc[-1][realization_cols], errors="coerce").dropna()
            final.index = final.index.astype(str)
            held = set(map(str, self.posterior._df.index)) & set(
                map(str, self.posterior_parameters._df.index)
            )
            final = final[final.index.isin(held)]
            if final.empty:
                raise ValueError(
                    "No realization of the last iteration has both a phi and a "
                    "posterior parameter set; PESTPP-IES dropped or failed every one."
                )
            return str(final.idxmin())
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

    def report(self, html_path, *, max_locations: int = 6) -> Path:
        """Write a single self-contained HTML report of the headline IES plots.

        Bundles, in order: phi convergence, the phi distribution, the
        ensemble-vs-observation comparison, phi contributions, parameters at
        bounds, one histogram per forecast, and -- when the run carries a model
        and captured fields -- the posterior mean and change maps of each.

        The bundle is the headline set, not every figure on the class -- among
        others it omits :meth:`plot_field` ``stat="reduction"``,
        :meth:`plot_field_mosaic`, :meth:`plot_obs_residuals`,
        :meth:`plot_prior_vs_obs` and :meth:`plot_conflict`. Call those directly.

        Parameters
        ----------
        html_path : str or pathlib.Path
            The HTML file to write, overwritten when it exists; its directory
            must exist. The plotly library is loaded from a CDN, so the page
            needs a network connection to draw.
        max_locations : int, optional
            The most target locations the ensemble-vs-observation figure draws
            (:meth:`plot_vs_obs`'s ``max_locations``). Default 6.

        Returns
        -------
        pathlib.Path
            The written HTML path.
        """

        html_path = Path(html_path)
        figures = [self.plot_phi(), self.plot_phi_distribution()]
        try:
            figures.append(self.plot_vs_obs(max_locations=max_locations))
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

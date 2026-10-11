from __future__ import annotations

from typing import TYPE_CHECKING

from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.viz import mpl_axes

if TYPE_CHECKING:
    import pandas as pd

import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from pandas import IndexSlice as idxx

from myflopy import viz as f
from myflopy._logging import get_logger
from myflopy.modflow.mf6.headsplus import _is_sentinel

logger = get_logger(__name__)


def calculate_calibration_statistics(observed, simulated):
    """
    Calculate common calibration statistics.

    Parameters:
    observed (np.array): Array of observed values.
    simulated (np.array): Array of simulated values.

    Returns:
    dict: Dictionary of calibration statistics.
    """
    if len(observed) != len(simulated):
        raise ValueError(f"Observed and simulated arrays must have the same length."
                         f"simulated length: {len(simulated)}, observed length: {len(observed)}")

    observed = np.array(observed, dtype=float).flatten()
    simulated = np.array(simulated, dtype=float).flatten()

    # Mean Error (ME)
    me = np.mean(simulated - observed)

    # Mean Absolute Error (MAE)
    mae = np.mean(np.abs(simulated - observed))

    # Root Mean Square Error (RMSE)
    rmse = np.sqrt(np.mean((simulated - observed) ** 2))

    # Normalized Root Mean Square Error (NRMSE)
    observed_range = np.max(observed) - np.min(observed)
    nrmse = rmse / observed_range

    # Nash-Sutcliffe Efficiency (NSE)
    nse = 1 - (np.sum((simulated - observed) ** 2) / np.sum((observed - np.mean(observed)) ** 2))

    # Coefficient of Determination (R^2)
    r_squared = np.corrcoef(observed, simulated)[0, 1] ** 2

    stats = {
        'ME': me,
        'MAE': mae,
        'RMSE': rmse,
        'NRMSE': nrmse,
        'NSE': nse,
        'R^2': r_squared
    }

    return stats


def _require_compare_columns(frame: pd.DataFrame, required: set[str], *, caller: str) -> pd.DataFrame:
    """Validate that a compare/residual frame includes the required columns."""

    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{caller} requires columns: {', '.join(sorted(missing))}")
    return frame.copy()


def _blank_no_value(frame: pd.DataFrame, columns) -> tuple[pd.DataFrame, int]:
    """``frame`` with MODFLOW's no-value markers in ``columns`` made ``NaN``, and how many.

    A compare table carries what a target READ, and on an SFR target's stage-file
    route that is the stage file raw: -1e30 for a DRY reach, +1e30 for an INACTIVE
    one (ledger 353 -- the user kept SFR targets as they were). A picture of it is
    a stage view, and a stage view leaves a reach with no stage blank (the user:
    "the plots for SFR should just not plot a stage if the reach is dry"), so the
    calibration pictures blank it by the head file's rule, ``headsplus._is_sentinel``
    (``|value| >= 1e29``) -- the cross plot's 1:1 line ran from -1e30 to +1e30, the
    time series drew -1e30 and the period MAE read 1e30 (ledger 355). Only the
    PICTURE changes: ``compare()`` and the files PEST reads keep the number.

    Parameters
    ----------
    frame : pandas.DataFrame
        A compare-style table; it is copied, never edited.
    columns : iterable of str
        The value columns to blank (a simulated value, an absolute residual);
        a name the frame does not carry is skipped.

    Returns
    -------
    tuple of (pandas.DataFrame, int)
        The copy, and the number of values blanked.
    """

    blanked = frame.copy()
    count = 0
    for column in columns:
        if column not in blanked.columns:
            continue
        values = pd.to_numeric(blanked[column], errors="coerce")
        no_value = _is_sentinel(values.to_numpy(dtype=float))
        count += int(no_value.sum())
        blanked[column] = values.mask(no_value)
    return blanked, count


def _describe_empty_pairing(
    frame: pd.DataFrame, target_column: str, simulated_column: str
) -> str:
    """Explain why a cross plot has no points, naming the side that is missing.

    A calibration plot with nothing on it is almost always a data problem --
    all-NaN observations, or a join that matched no rows -- and styled-but-empty
    axes hide that completely. Diagnosing it here means the caller is told which
    column failed instead of being left to guess from a blank figure.
    """

    observed = frame[target_column].notna().sum()
    simulated = frame[simulated_column].notna().sum()
    if not len(frame):
        return "the compare table is empty (no target rows at all)."
    if not observed and not simulated:
        return (
            f"both {target_column!r} and {simulated_column!r} are entirely NaN "
            f"across {len(frame)} rows."
        )
    if not observed:
        return (
            f"{target_column!r} is entirely NaN across {len(frame)} rows -- the "
            f"targets carry no measured values, so there is nothing to compare "
            f"the {simulated} simulated values against."
        )
    if not simulated:
        return (
            f"{simulated_column!r} is entirely NaN across {len(frame)} rows -- "
            f"the model produced no simulated values at the {observed} target "
            f"locations (was the model run?)."
        )
    return (
        f"no row has BOTH {target_column!r} and {simulated_column!r} present "
        f"({observed} observed and {simulated} simulated values, but they never "
        f"coincide) -- usually a name/time join that did not line up."
    )


def _warn_empty_cross_plot(reason: str, *, caller: str) -> str:
    """Warn about an empty cross plot and return the on-figure annotation text."""

    import warnings

    warnings.warn(
        f"{caller} has no points to plot: {reason}",
        UserWarning,
        stacklevel=3,
    )
    return f"No paired observed/simulated values<br><br>{reason}"


def _sort_compare_by_time(frame: pd.DataFrame) -> pd.DataFrame:
    """Sort a compare-style DataFrame by time/per when possible."""

    if "time" not in frame.columns:
        return frame.reset_index(drop=True)
    sort_values = pd.to_numeric(frame["time"], errors="coerce")
    if sort_values.notna().all():
        return frame.assign(_sort=sort_values).sort_values("_sort").drop(columns="_sort").reset_index(drop=True)
    return frame.sort_values("time").reset_index(drop=True)


def _multiindex_head_frames(compare: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build observed/simulated MultiIndex head frames from compare output."""

    frame = _require_compare_columns(
        compare,
        {"name", "layer", "head_target", "sim_head"},
        caller="CalibrationPlot.from_compare(type='heads')",
    )
    time_values = pd.to_numeric(frame["per"] if "per" in frame.columns else frame["time"], errors="coerce")
    if time_values.isna().any():
        raise ValueError(
            "CalibrationPlot.from_compare(type='heads') requires numeric per/time values."
        )
    tuples = list(
        zip(
            frame["name"].astype(str),
            pd.to_numeric(frame["layer"], errors="coerce").fillna(0).astype(int),
            [(0, int(value)) for value in time_values], strict=False,
        )
    )
    index = pd.MultiIndex.from_tuples(tuples, names=["locs", "layer", "kstpkper"])
    observed = pd.DataFrame({"elev": frame["head_target"].to_numpy()}, index=index)
    simulated = pd.DataFrame({"elev": frame["sim_head"].to_numpy()}, index=index)
    return observed, simulated


def _infer_compare_value_columns(
    frame: pd.DataFrame,
    *,
    target_column: str | None = None,
    simulated_column: str | None = None,
    caller: str,
) -> tuple[str, str]:
    """Resolve target/simulated value columns for compare-style frames."""

    if target_column is not None and simulated_column is not None:
        _require_compare_columns(
            frame,
            {target_column, simulated_column},
            caller=caller,
        )
        return target_column, simulated_column

    known_pairs = [
        ("head_target", "sim_head"),
        ("stage_target", "sim_stage"),
        ("flow_target", "sim_flow"),
        ("target_value", "sim_value"),
        ("value_target", "sim_value"),
    ]
    for candidate_target, candidate_simulated in known_pairs:
        if {candidate_target, candidate_simulated}.issubset(frame.columns):
            return candidate_target, candidate_simulated

    raise ValueError(
        f"{caller} could not infer target/simulated value columns. "
        "Pass target_column= and simulated_column= explicitly."
    )


# Matplotlib/seaborn colors mirroring the Plotly defaults, for backend="matplotlib".
_MPL_MODEL = "#1f77b4"
_MPL_BASELINE = "0.6"
_MPL_TARGET = "crimson"


def _normalize_backend(backend: str) -> str:
    """Map friendly backend names to ``"plotly"`` or ``"matplotlib"``."""

    value = str(backend).strip().lower()
    if value in ("plotly", "go"):
        return "plotly"
    if value in ("matplotlib", "mpl", "seaborn", "sns"):
        return "matplotlib"
    raise ValueError(f"backend must be 'plotly' or 'matplotlib', got {backend!r}.")


def _mpl_axes(**kwargs):
    """Create a seaborn-styled matplotlib figure/axes without global side effects."""

    import seaborn as sns

    with sns.axes_style("whitegrid"):
        fig, ax = mpl_axes(**kwargs)
    return fig, ax


class CalibrationPlot(f.Fig):

    @classmethod
    def from_compare(
        cls,
        compare: pd.DataFrame,
        *,
        type: str = "calibration",
        target_column: str | None = None,
        simulated_column: str | None = None,
        title: str | None = None,
        yaxis_title: str | None = None,
    ):
        """Build a calibration plot directly from a compare/residual DataFrame.

        Always a Plotly figure: this entry point has no ``backend=``.

        Parameters
        ----------
        compare : pandas.DataFrame
            A compare or residual table, one row per target row, as a target
            set's ``compare(model)`` returns it. It is copied, never edited. It
            must hold a target/simulated value pair: the two columns named
            below, else the first of ``("head_target", "sim_head")``,
            ``("stage_target", "sim_stage")``, ``("flow_target", "sim_flow")``,
            ``("target_value", "sim_value")`` and ``("value_target",
            "sim_value")`` it carries. With neither, ``ValueError``. That
            inference runs for every ``type``, including those that do not
            plot the pair.
        type : str, default "calibration"
            Which plot:

            - ``"calibration"``: target (x) against simulated (y), one point per
              row with fit statistics, titled "Calibration" by default;
            - ``"obs_vs_sim"``: the same cross plot without the statistics,
              titled "Observed vs simulated" by default;
            - ``"heads"``: target and simulated values against the period, one
              pair of lines per location and layer. A table without ``name``
              or ``layer`` raises ``ValueError``. The period is ``per``, else
              ``time``: a table with neither raises ``KeyError``, and a
              non-numeric or missing value raises ``ValueError``. Every
              location must have rows in EVERY layer the table names, or
              ``KeyError`` (a well in layer 0 beside one in layer 1 raises).
              The pair is read as heads, and ``title`` and ``yaxis_title`` are
              ignored;
            - ``"residuals_by_period"``: mean ``abs_residual`` per ``time``
              (:meth:`from_residuals_by_period`). ``title``, ``yaxis_title`` and
              the value pair are ignored.

            Any other string, a misspelling included, draws the
            ``"calibration"`` plot without a word.
        target_column : str, optional
            The column of target values (the x axis of a cross plot). It is used
            only when ``simulated_column`` is ALSO given. Given alone, it is
            ignored in silence and the pair is inferred as above. A named column
            the table lacks raises ``ValueError``.
        simulated_column : str, optional
            The column of simulated values (the y axis), under the same rule:
            it counts only together with ``target_column``.
        title : str, optional
            Title of a cross plot. The default ``None`` (or ``""``) gives that
            plot's own title, listed under ``type``.
        yaxis_title : str, optional
            Y-axis title of a cross plot. The default ``None`` (or ``""``) gives
            ``"Simulated"``. The x-axis is always titled ``"Observed"``.

        Returns
        -------
        CalibrationPlot
            A Plotly figure (a ``viz.Fig``).
        """

        resolved_target, resolved_simulated = _infer_compare_value_columns(
            compare,
            target_column=target_column,
            simulated_column=simulated_column,
            caller="CalibrationPlot.from_compare(...)",
        )
        frame = compare.copy()

        if type == "heads":
            if {resolved_target, resolved_simulated} != {"head_target", "sim_head"}:
                frame = frame.rename(columns={resolved_target: "head_target", resolved_simulated: "sim_head"})
            observed, simulated = _multiindex_head_frames(frame)
            return cls(observed=observed, simulated=simulated, type=type)

        if type == "obs_vs_sim":
            return cls.from_obs_vs_sim(
                frame,
                target_column=resolved_target,
                simulated_column=resolved_simulated,
                title=title or "Observed vs simulated",
                yaxis_title=yaxis_title or "Simulated",
            )

        if type == "residuals_by_period":
            return cls.from_residuals_by_period(frame)

        return cls.from_obs_vs_sim(
            frame,
            target_column=resolved_target,
            simulated_column=resolved_simulated,
            title=title or "Calibration",
            yaxis_title=yaxis_title or "Simulated",
            add_stats=True,
        )

    @classmethod
    def from_targets(cls, targets, *, model=None, type: str = "calibration"):
        """Build a calibration plot from ``HeadTargets`` or a bound target helper.

        The value columns are inferred from the comparison table
        (:meth:`from_compare`), and titles take their defaults. Use a target
        set's own ``plot`` for the family's titles.

        Parameters
        ----------
        targets : object
            Anything with a ``compare(model)`` method returning a compare table:
            a raw target set (``HeadTargets``, ``ConcTargets``,
            ``LakeStageTargets``, ``SfrStageTargets``, ``SfrFlowTargets``,
            ``DrnFlowTargets``) or a model-bound helper (``model.targets.heads``
            and its siblings). An object with a ``targets`` attribute is read as
            BOUND.
        model : SimulationBase, optional
            The run model to compare against. A raw target set needs it, and
            ``ValueError`` follows without it. For a bound helper, the default
            ``None`` uses the model it is bound to, and a model given here is
            compared INSTEAD.
        type : str, default "calibration"
            Which plot, passed to :meth:`from_compare`, whose entry lists the
            values. A misspelling draws the ``"calibration"`` plot. ``"heads"``
            needs a ``layer`` column, which lake, SFR and drain comparisons do
            not have (``ValueError``), and raises ``KeyError`` when a head
            target location has no row in a layer another location has.

        Returns
        -------
        CalibrationPlot
            A Plotly figure (a ``viz.Fig``).
        """

        if hasattr(targets, "targets") and model is None:
            compare = targets.compare()
        else:
            if model is None:
                raise ValueError(
                    "CalibrationPlot.from_targets(...) requires model= when given raw HeadTargets."
                )
            compare = targets.compare(model)
        return cls.from_compare(compare, type=type)

    @classmethod
    def from_obs_vs_sim(
        cls,
        compare: pd.DataFrame,
        *,
        baseline_compare: pd.DataFrame | None = None,
        target_column: str = "head_target",
        simulated_column: str = "sim_head",
        title: str = "Observed vs simulated heads",
        xaxis_title: str = "Observed",
        yaxis_title: str = "Simulated",
        add_stats: bool = False,
        backend: str = "plotly",
    ):
        """Build a target-vs-simulated cross plot from compare-style DataFrames.

        ``backend`` selects ``"plotly"`` (interactive, default) or
        ``"matplotlib"`` (static matplotlib/seaborn). A simulated value that is
        MODFLOW's no-value marker (``|value| >= 1e29``: a dry or inactive SFR
        reach on a stage target's stage-file route) is not drawn and is left out
        of the 1:1 line and the statistics (ledger 355); the compare table is not
        changed.

        A row missing either value is not drawn. When no row has both, the
        figure carries a note saying which side is missing and a
        ``UserWarning`` is issued.

        Parameters
        ----------
        compare : pandas.DataFrame
            A compare table holding ``target_column`` and ``simulated_column``
            (``ValueError`` naming a missing one). It is copied, never edited.
        baseline_compare : pandas.DataFrame, optional
            A second compare table with the same two columns, such as the run
            before a calibration, drawn as a "Baseline" series. Its values
            widen the 1:1 line and stay out of the statistics. The default
            ``None`` draws none.
        target_column : str, default "head_target"
            The column of target values, on the x axis.
        simulated_column : str, default "sim_head"
            The column of simulated values, on the y axis. Its no-value markers
            are blanked as described above.
        title : str, default "Observed vs simulated heads"
            Figure title.
        xaxis_title : str, default "Observed"
            X-axis title.
        yaxis_title : str, default "Simulated"
            Y-axis title.
        add_stats : bool, default False
            Annotate ME, MAE, RMSE, NRMSE, NSE and R^2 of the model's points
            (simulated minus target; ME, MAE and RMSE in the values' own units,
            the other three unitless), upper left: to 3
            decimals on Plotly (:meth:`add_stats`), to 3 significant figures on
            Matplotlib. Nothing is added when no row has both values.
        backend : str, default "plotly"
            ``"plotly"`` (or ``"go"``) returns this interactive Plotly figure.
            ``"matplotlib"`` (or ``"mpl"``, ``"seaborn"``, ``"sns"``) returns a
            static ``matplotlib.figure.Figure`` in seaborn's whitegrid style.
            Case and surrounding spaces are ignored, and any other value raises
            ``ValueError``.

        Returns
        -------
        CalibrationPlot or matplotlib.figure.Figure
            By ``backend``.
        """

        checked, blanked = _blank_no_value(
            _require_compare_columns(
                compare,
                {target_column, simulated_column},
                caller="CalibrationPlot.from_obs_vs_sim(...)",
            ),
            [simulated_column],
        )
        frame = checked.dropna(subset=[target_column, simulated_column])
        empty_note = None
        if frame.empty:
            reason = _describe_empty_pairing(checked, target_column, simulated_column)
            if blanked:
                reason += (
                    f" {blanked} simulated values were MODFLOW's no-value marker "
                    "(|value| >= 1e29: a dry or inactive feature) and are not drawn."
                )
            empty_note = _warn_empty_cross_plot(reason, caller="CalibrationPlot.from_obs_vs_sim(...)")

        baseline = None
        values = [frame[target_column], frame[simulated_column]]
        if baseline_compare is not None:
            baseline = _blank_no_value(
                _require_compare_columns(
                    baseline_compare,
                    {target_column, simulated_column},
                    caller="CalibrationPlot.from_obs_vs_sim(..., baseline_compare=...)",
                ),
                [simulated_column],
            )[0].dropna(subset=[target_column, simulated_column])
            values.extend([baseline[target_column], baseline[simulated_column]])
        merged_values = pd.concat(values, axis=0).dropna()
        one_to_one = (float(merged_values.min()), float(merged_values.max())) if not merged_values.empty else None

        if _normalize_backend(backend) == "matplotlib":
            import seaborn as sns

            fig, ax = _mpl_axes(figsize=(5.5, 5.5))
            ax.scatter(frame[target_column], frame[simulated_column], s=22, color=_MPL_MODEL, alpha=0.75, label="Model")
            if baseline is not None:
                ax.scatter(baseline[target_column], baseline[simulated_column], s=22, color=_MPL_BASELINE, alpha=0.75, label="Baseline")
            if one_to_one is not None:
                ax.plot(one_to_one, one_to_one, color=_MPL_TARGET, ls="--", lw=1, label="1:1 line")
                ax.set_aspect("equal", adjustable="box")
            ax.set_xlabel(xaxis_title)
            ax.set_ylabel(yaxis_title)
            ax.set_title(title)
            ax.legend()
            if empty_note is not None:
                ax.text(
                    0.5,
                    0.5,
                    empty_note.replace("<br>", "\n"),
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    fontsize=9,
                    wrap=True,
                    color="#8a2f2f",
                )
            if add_stats and not frame.empty:
                stats = calculate_calibration_statistics(frame[target_column], frame[simulated_column])
                text = "\n".join(f"{key} = {value:.3g}" for key, value in stats.items())
                ax.text(0.03, 0.97, text, transform=ax.transAxes, va="top", ha="left", fontsize=8,
                        bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))
            sns.despine(fig)
            return fig

        fig = cls()
        fig.add_scattergl(
            x=frame[target_column],
            y=frame[simulated_column],
            mode="markers",
            name="Model",
        )
        if baseline is not None:
            fig.add_scattergl(
                x=baseline[target_column],
                y=baseline[simulated_column],
                mode="markers",
                name="Baseline",
            )
        if one_to_one is not None:
            fig.add_scattergl(
                x=list(one_to_one),
                y=list(one_to_one),
                mode="lines",
                name="1:1 line",
            )
        fig.update_layout(
            title=title,
            xaxis_title=xaxis_title,
            yaxis_title=yaxis_title,
        )
        if empty_note is not None:
            fig.add_annotation(
                text=empty_note,
                showarrow=False,
                xref="paper",
                yref="paper",
                x=0.5,
                y=0.5,
                align="center",
                font={"size": 13, "color": "#8a2f2f"},
                bgcolor="rgba(255,255,255,0.85)",
                bordercolor="#8a2f2f",
                borderwidth=1,
                borderpad=10,
            )
        if add_stats and not frame.empty:
            fig.add_stats(observed=frame[target_column], simulated=frame[simulated_column])
        return fig

    @classmethod
    def from_timeseries(
        cls,
        compare: pd.DataFrame,
        *,
        name: str | None = None,
        baseline_compare: pd.DataFrame | None = None,
        target_column: str = "head_target",
        simulated_column: str = "sim_head",
        target_label: str = "Target",
        simulated_label: str = "Model",
        baseline_label: str = "Baseline",
        yaxis_title: str = "Head",
        title: str | None = None,
        backend: str = "plotly",
    ):
        """Build a time-series plot of targets and simulated heads for one target.

        ``backend`` selects ``"plotly"`` (interactive, default) or
        ``"matplotlib"`` (static matplotlib/seaborn). A simulated value that is
        MODFLOW's no-value marker (``|value| >= 1e29``: a dry or inactive SFR
        reach on a stage target's stage-file route) is a gap in the line, the
        model's and the baseline's (ledger 355); the compare table is not changed.

        Rows are drawn in ``time`` order, numerically when every ``time`` is a
        number. The x-axis is titled "Time / period".

        Parameters
        ----------
        compare : pandas.DataFrame
            A compare table holding ``name``, ``time``, ``target_column`` and
            ``simulated_column`` (``ValueError`` naming a missing one). It is
            copied, never edited.
        name : str, optional
            The target location to draw. It is matched against ``name``
            case-insensitively, after stripping spaces from the value given
            here. The default ``None`` is the table's one location, and it
            raises ``ValueError`` when the table holds more than one. A name
            no row carries raises ``ValueError``.
        baseline_compare : pandas.DataFrame, optional
            A second compare table, drawn as a third line from its rows for the
            same ``name``. It needs ``name``, ``time`` and ``simulated_column``
            only; its targets are not drawn. A baseline with no row for that
            name draws an empty line, with no error. The default ``None`` draws
            none.
        target_column : str, default "head_target"
            The column of target values.
        simulated_column : str, default "sim_head"
            The column of simulated values, read from ``compare`` and from
            ``baseline_compare``.
        target_label : str, default "Target"
            Legend name of the target line.
        simulated_label : str, default "Model"
            Legend name of the simulated line.
        baseline_label : str, default "Baseline"
            Legend name of the baseline line. Unused without
            ``baseline_compare``.
        yaxis_title : str, default "Head"
            Y-axis title. Give the quantity for a stage or flow table.
        title : str, optional
            Figure title. The default ``None`` (or ``""``) is the location's
            name as the earliest matching row (in ``time`` order) spells it.
        backend : str, default "plotly"
            ``"plotly"`` (or ``"go"``) returns this interactive Plotly figure.
            ``"matplotlib"`` (or ``"mpl"``, ``"seaborn"``, ``"sns"``) returns a
            static ``matplotlib.figure.Figure`` in seaborn's whitegrid style.
            Case and surrounding spaces are ignored, and any other value raises
            ``ValueError``.

        Returns
        -------
        CalibrationPlot or matplotlib.figure.Figure
            By ``backend``.
        """

        frame, _ = _blank_no_value(
            _require_compare_columns(
                compare,
                {"name", "time", target_column, simulated_column},
                caller="CalibrationPlot.from_timeseries(...)",
            ),
            [simulated_column],
        )
        names = frame["name"].astype(str)
        if name is None:
            unique_names = sorted(names.unique().tolist())
            if len(unique_names) != 1:
                raise ValueError(
                    "CalibrationPlot.from_timeseries(...) requires name= when more than one target location is present."
                )
            name = unique_names[0]
        key = str(name).strip().lower()
        frame = frame.loc[names.str.lower() == key].copy()
        if frame.empty:
            raise ValueError(f"No compare rows found for observation name {name!r}.")
        frame = _sort_compare_by_time(frame)

        baseline = None
        if baseline_compare is not None:
            baseline, _ = _blank_no_value(
                _require_compare_columns(
                    baseline_compare,
                    {"name", "time", simulated_column},
                    caller="CalibrationPlot.from_timeseries(..., baseline_compare=...)",
                ),
                [simulated_column],
            )
            baseline_names = baseline["name"].astype(str)
            baseline = _sort_compare_by_time(baseline.loc[baseline_names.str.lower() == key].copy())
        resolved_title = title or str(frame["name"].iloc[0])

        if _normalize_backend(backend) == "matplotlib":
            import seaborn as sns

            fig, ax = _mpl_axes(figsize=(8, 4))
            ax.plot(frame["time"], frame[target_column], "o-", color=_MPL_TARGET, label=target_label)
            ax.plot(frame["time"], frame[simulated_column], "o-", color=_MPL_MODEL, label=simulated_label)
            if baseline is not None:
                ax.plot(baseline["time"], baseline[simulated_column], "o-", color=_MPL_BASELINE, label=baseline_label)
            ax.set_xlabel("Time / period")
            ax.set_ylabel(yaxis_title)
            ax.set_title(resolved_title)
            ax.legend()
            sns.despine(fig)
            return fig

        fig = cls()
        fig.add_scattergl(
            x=frame["time"],
            y=frame[target_column],
            mode="lines+markers",
            name=target_label,
        )
        fig.add_scattergl(
            x=frame["time"],
            y=frame[simulated_column],
            mode="lines+markers",
            name=simulated_label,
        )
        if baseline is not None:
            fig.add_scattergl(
                x=baseline["time"],
                y=baseline[simulated_column],
                mode="lines+markers",
                name=baseline_label,
            )
        fig.update_layout(
            title=resolved_title,
            xaxis_title="Time / period",
            yaxis_title=yaxis_title,
        )
        return fig

    @classmethod
    def from_residuals_by_period(
        cls,
        compare: pd.DataFrame,
        *,
        baseline_compare: pd.DataFrame | None = None,
        title: str = "Residual MAE by period",
        yaxis_title: str = "Mean absolute error",
        backend: str = "plotly",
    ):
        """Build a by-period residual summary plot from compare-style data.

        ``backend`` selects ``"plotly"`` (interactive, default) or
        ``"matplotlib"`` (static matplotlib/seaborn). A residual against
        MODFLOW's no-value marker (``|abs_residual| >= 1e29``: a dry or inactive
        SFR reach on a stage target's stage-file route) is left out of its
        period's mean, and a period with nothing else is a gap (ledger 355); the
        compare table is not changed.

        Rows are grouped by their ``time`` value, not by ``per``. A row whose
        ``abs_residual`` is NaN (an unmeasured target row) is skipped by the
        mean. The x-axis is titled "Time / period".

        Parameters
        ----------
        compare : pandas.DataFrame
            A compare table holding ``time`` and ``abs_residual`` (in the target
            values' own units). ``ValueError`` names a missing one. It is never
            edited.
        baseline_compare : pandas.DataFrame, optional
            A second compare table with the same two columns, drawn as
            "Baseline MAE". The default ``None`` draws none.
        title : str, default "Residual MAE by period"
            Figure title.
        yaxis_title : str, default "Mean absolute error"
            Y-axis title.
        backend : str, default "plotly"
            ``"plotly"`` (or ``"go"``) returns this interactive Plotly figure.
            ``"matplotlib"`` (or ``"mpl"``, ``"seaborn"``, ``"sns"``) returns a
            static ``matplotlib.figure.Figure`` in seaborn's whitegrid style.
            Case and surrounding spaces are ignored, and any other value raises
            ``ValueError``.

        Returns
        -------
        CalibrationPlot or matplotlib.figure.Figure
            By ``backend``.
        """

        def _period_mae(data, caller):
            """Per-period mean absolute residual from a compare frame, sorted by time."""

            checked, _ = _blank_no_value(
                _require_compare_columns(data, {"time", "abs_residual"}, caller=caller),
                ["abs_residual"],
            )
            stats = (
                checked
                .groupby("time")
                .agg(mae=("abs_residual", "mean"))
                .reset_index()
            )
            return _sort_compare_by_time(stats.rename(columns={"mae": "value"}))

        current_stats = _period_mae(compare, "CalibrationPlot.from_residuals_by_period(...)")
        baseline_stats = None
        if baseline_compare is not None:
            baseline_stats = _period_mae(
                baseline_compare, "CalibrationPlot.from_residuals_by_period(..., baseline_compare=...)"
            )

        if _normalize_backend(backend) == "matplotlib":
            import seaborn as sns

            fig, ax = _mpl_axes(figsize=(7, 4))
            ax.plot(current_stats["time"], current_stats["value"], "o-", color=_MPL_MODEL, label="Model MAE")
            if baseline_stats is not None:
                ax.plot(baseline_stats["time"], baseline_stats["value"], "o-", color=_MPL_BASELINE, label="Baseline MAE")
            ax.set_xlabel("Time / period")
            ax.set_ylabel(yaxis_title)
            ax.set_title(title)
            ax.legend()
            sns.despine(fig)
            return fig

        fig = cls()
        fig.add_scattergl(
            x=current_stats["time"],
            y=current_stats["value"],
            mode="lines+markers",
            name="Model MAE",
        )
        if baseline_stats is not None:
            fig.add_scattergl(
                x=baseline_stats["time"],
                y=baseline_stats["value"],
                mode="lines+markers",
                name="Baseline MAE",
            )
        fig.update_layout(
            title=title,
            xaxis_title="Time / period",
            yaxis_title=yaxis_title,
        )
        return fig

    def __init__(
            self,
            observed: pd.Series | list | Path = None,
            simulated: pd.Series | list = None,
            model: SimulationBase = None,
            loc_cell_dict: dict = None,
            loc_shp_gpkg: Path = None,
            lak_obs_dict: dict = None,
            cal_range: list | tuple = None,
            exploration_name_field: str = 'ExploName',
            crs=None,
            obs_layers: int | list[int] = None,
            verbose: bool = False,
            type: str = 'calibration'
    ):
        """
        Initializes an object to plot observed and simulated calibration data.

        :param observed: Observed data, provided as a pandas Series, a list, or a Path
            object pointing to the data.
        :param simulated: Simulated model data, provided as a pandas Series or a list.
            if not provided, simulated data is generated using the model object and
            and spatial data provided in loc_shp_gpkg.
        :param model: Simulation model object adhering to the SimulationBase class.
        :param loc_cell_dict: Dictionary linking geospatial locations to cell mappings
            or identifiers.
        :param loc_shp_gpkg: Path to shapefile or GeoPackage used for geospatial
            location of point observations of heads (i.e. wells)
        :param lak_obs_dict: Dictionary linking lake names to lake number ids.
        :param cal_range: List or tuple specifying the calibration range, typically as
            time steps or indices for calibration analysis.
        :param exploration_name_field: Field name in the geospatial file (e.g.,
            shapefile or GeoPackage) representing exploration or location names.
        :param crs: Coordinate reference system used for geospatial data processing
            and mapping.
        :param obs_layers: Integer or list of integers representing observation layers
        :param type: can be 'calibration' for cross plot or 'heads' to show actual heads as a plot.
        """
        super().__init__()
        self._observed = None
        self._simulated = None
        self._model = None
        self._crs = crs
        self._cal_range = cal_range

        self._loc_cell_dict = loc_cell_dict
        self._loc_shp_gpkg = None
        self._lak_obs_dict = None
        self._exploration_name_field = exploration_name_field
        self._obs_layers = obs_layers
        self._verbose = verbose

        self.model = model
        self.lak_obs_dict = lak_obs_dict
        self.loc_shp_gpkg = loc_shp_gpkg
        self.simulated = simulated
        self.observed = observed

        if self.observed is not None and self.simulated is not None:
            if type == 'calibration':
                self.add_calib_scatter()
                self.add_stats()
                self.set_layout()
            if type == 'heads':
                self.add_heads()


    @property
    def model(self):
        """The bound MODFLOW simulation used to generate simulated values."""

        return self._model

    @model.setter
    def model(self, value):
        """Set the source model, asserting it is a :class:`SimulationBase` (or ``None``)."""

        if value is not None:
            assert isinstance(value, SimulationBase), 'model must be a SimulationBase object'
        self._model = value

    @property
    def loc_shp_gpkg(self):
        """Path to the observation-location shapefile/GeoPackage used to sample simulated heads."""

        return self._loc_shp_gpkg

    @loc_shp_gpkg.setter
    def loc_shp_gpkg(self, value):
        """Set the observation-location file, asserting it is a :class:`~pathlib.Path` (or ``None``)."""

        if value is not None:
            assert isinstance(value, Path), 'loc_shp_gpkg must be a Path object'
        self._loc_shp_gpkg = value

    @property
    def simulated(self):
        """gets simulated data from model object if no simulated data is provided"""
        if self._simulated is None and self.loc_shp_gpkg and self.model:

            if self._verbose:
                print('generating simulated data using model object and spatial data provided')
            simulated = self.model.hds.get_obs_heads(locs=self.loc_shp_gpkg, long_format=True)

            # add lake stage data to simulated data
            if self.lak_obs_dict is not None:
                lak_obs = {}
                # get lake stage data for each lake
                for k, v in self.lak_obs_dict.items():
                    lak_obs[k] = self.model.outputs.lak.stage.get()[:, v].tolist()
                # create a multi-index for lake stage data
                lake_names = list(lak_obs.keys())
                layers = [0]  # set layer to 0 for all lake observations
                kstpkper = self.model.gwf.output.head().get_kstpkper()
                assert len(lak_obs[lake_names[0]]) == len(kstpkper), \
                    'length of lake stage data does not match number of stress periods in model'
                lak_idx = pd.MultiIndex.from_product(
                    [lake_names, layers, kstpkper], names=['locs', 'layer', 'kstpkper'])
                lake_data = list(itertools.chain.from_iterable(lak_obs.values()))
                # create a DataFrame with lake stage data and append to simulated data
                lak_df = pd.DataFrame(lake_data, index=lak_idx, columns=['elev'])
                simulated = pd.concat([simulated, lak_df], axis=0)

            simulated = simulated[sorted(simulated.columns)]
            self._simulated = simulated

        return self._simulated

    @simulated.setter
    def simulated(self, value):
        """Set simulated data directly; ignored when ``loc_shp_gpkg`` is set (generated instead)."""

        if self.loc_shp_gpkg is not None:
            if self._verbose:
                print('path to shapefile or geopackage provided, '
                      'so provided simulated data is ignored, and'
                      'will generate observed data using model'
                      'object and spatial data provided')
        else:
            self._simulated = value

    @property
    def lak_obs_dict(self):
        """Mapping of lake observation name to lake id, folded into simulated data."""

        return self._lak_obs_dict

    @lak_obs_dict.setter
    def lak_obs_dict(self, value):
        """Set the lake-observation dict, asserting each lake id exists in the LAK package."""

        if value is not None:
            assert isinstance(value, dict), 'lak_obs_dict must be a dictionary'
            for lak_id in value.values():
                assert lak_id in range(self.model.outputs.lak.stage.nlakes), \
                    'lak_id not found in lak package'
        self._lak_obs_dict = value

    @property
    def observed(self):
        """The observed target data (long-format, indexed by kstpkper/layer)."""

        return self._observed

    @observed.setter
    def observed(self, value):
        """Set observed data from a DataFrame or an Excel ``Path`` (aligned to model kstpkper)."""

        if isinstance(value, Path):
            obs_data = pd.read_excel(value)
            assert len(obs_data) == self.model.nper, \
                'length of observed data does not match number of stress periods in model'
            if self._obs_layers is None:
                self._obs_layers = 0  # default to layer 1
            col1 = obs_data.columns[0]  # should be stress periods
            if len(obs_data) != len(self.model.kstpkper):
                logger.warning(
                    'the model has %d time steps but the observed data has %d rows; '
                    'trimming the observed data to match',
                    len(self.model.kstpkper), len(obs_data),
                )
                obs_data = obs_data.iloc[:len(self.model.kstpkper), :]
            obs_data['kstpkper'] = self.model.kstpkper
            # drop stress period column in favor of model kstpkper to match simulated data
            obs_data.drop(columns=col1, inplace=True)
            obs_data = obs_data[sorted(obs_data.columns)]
            obs_data['layer'] = self._obs_layers
            obs_data = obs_data.set_index(['kstpkper', 'layer'])
            obs_data = obs_data.melt(ignore_index=False, value_name='elev', var_name='locs')
            obs_data = obs_data.reset_index().set_index(['locs', 'layer', 'kstpkper'])
            value = obs_data
            """except ValueError as e:
                print('observed data path not readable, must be Excel file')
                traceback.print_exc()"""
        self._observed = value

    def add_stats(
            self,
            observed: list = None,
            simulated: list = None,
    ):
        """Annotate this figure with calibration statistics and return them as text.

        The statistics come from :func:`calculate_calibration_statistics`, on
        ``simulated - observed``: ME, MAE, RMSE, NRMSE (RMSE over the observed
        range), NSE and R^2, in the values' own units where they have one. They
        are written to 3 decimals in a box at the upper left.

        Parameters
        ----------
        observed : array-like or pandas.DataFrame, optional
            The target values. The default ``None`` uses ``self.observed``, set
            by the constructor or the ``observed`` property. On a figure built
            by a ``from_*`` method it is None, and leaving either argument out
            then raises ``TypeError``.
        simulated : array-like or pandas.DataFrame, optional
            The simulated values, defaulting to ``self.simulated`` the same way.
            When BOTH are DataFrames, each must have exactly one column
            (``ValueError`` otherwise). They are paired by INDEX, and a row
            with a NaN on either side is dropped. Otherwise both are flattened
            and paired by POSITION. Their lengths must match (``ValueError``),
            and a NaN is kept, which makes every statistic NaN.

        Returns
        -------
        str
            The statistics, one ``"name: value"`` per line, joined by ``<br>``.
        """

        observed = self.observed if observed is None else observed
        simulated = self.simulated if simulated is None else simulated

        if all(isinstance(x, pd.DataFrame) for x in [observed, simulated]):
            logger.debug('observed and simulated data were given as DataFrames')
            sim_obs = pd.concat([simulated, observed], axis=1).dropna()
            sim_obs.columns = ['simulated', 'observed']
            # drop nan rows so the calib stats are calculated correctly
            sim_obs.dropna()
            simulated = sim_obs.loc[:, 'simulated'].to_numpy()
            observed = sim_obs.loc[:, 'observed'].to_numpy()

        stats = calculate_calibration_statistics(observed, simulated)
        calib_stats = [f'{k}: {v:.3f}' for k, v in stats.items()]

        # Create a single string with line breaks
        stats_text = "<br>".join(calib_stats)

        # Add annotation in the upper-left area
        self.add_annotation(
            x=0.02,  # small offset from the left border
            y=0.98,  # near top
            xref="paper",  # 'paper' means figure reference (0=left, 1=right)
            yref="paper",  # 'paper' means figure reference (0=bottom, 1=top)
            text=stats_text,
            showarrow=False,
            align="left",  # text alignment
            bordercolor="black",
            borderwidth=1,
            borderpad=4,
            bgcolor="white",
            opacity=0.8
        )
        return stats_text

    def add_calib_scatter(self):
        """
        Adds scatter plot elements for calibration visualization in the plot. This method
        handles two scenarios: calibration using simple datasets (e.g., lists or
        pandas Series) and calibration using observation data linked to shapefiles
        or geopackages. If calibration data is provided with a shapefile or geopackage,
        the observed and calibrated data provided to class must be in DataFrame format.

        :param self: The instance of the class where this method is defined.
        :return: None
        """
        cal_range = self._cal_range if self._cal_range is not None else None

        # if observed and simulated data are simple lists or Series, add all data as one Scattergl trace
        if all([isinstance(x, pd.Series | list) for x in [self._observed, self.simulated]]):
            self.add_scattergl(x=self._observed, y=self._simulated, mode='markers', name='data')
            data_min = self._observed.min()
            data_max = self._observed.max()

        # if a shapefile or geopackage was provided, add each observation
        # location as a separate Scattergl trace'
        if self.loc_shp_gpkg is not None:
            data_min = None
            data_max = None
            if self._verbose:
                print('using observation data based on provided'
                      'shapefile or geopackage for calibration plot')
            assert isinstance(self.simulated, pd.DataFrame)
            loc_names = self.simulated.index.get_level_values(level=0).unique().to_list()
            layers = self.simulated.index.get_level_values(level=1).unique().to_list()

            for loc_name in loc_names:
                for layer in layers:

                    sim_data = pd.Series(self.simulated.loc[idxx[loc_name, layer, :], 'elev'].to_list())
                    obs_data = pd.Series(self.observed.loc[idxx[loc_name, layer, :], 'elev'].to_list())
                    sim_obs = pd.concat([sim_data, obs_data], axis=1)
                    sim_obs.columns = ['simulated', 'observed']
                    sim_obs = sim_obs.dropna()

                    # get max and min of obs data to draw cal line
                    if data_min is None:
                        data_min = obs_data.min()
                    if data_max is None:
                        data_max = obs_data.max()
                    if obs_data.min() < data_min:
                        data_min = obs_data.min()
                    if obs_data.max() > data_max:
                        data_max = obs_data.max()

                    self.add_scattergl(
                        x=sim_obs.loc[:, 'simulated'],
                        y=sim_obs.loc[:, 'observed'],
                        mode='markers',
                        name=f'{loc_name} {layer}')

        if cal_range is None:
            cal_range = [data_min, data_max]

        # add 1:1 calibration line
        self.add_scattergl(
            x=cal_range,
            y=cal_range,
            mode='lines',
            name='1:1 line'
        )

    def set_layout(self):
        """
        Updates the layout configuration for a plot with specific titles for axes, legend,
        and the plot itself.
        :return: None
        """
        self.update_layout(
            title='Calibration',
            xaxis_title='Simulated',
            yaxis_title='Observed'
        )

    def add_heads(self):
        """Add observed vs simulated head time-series traces, one per location and layer."""

        locs = list(self.observed.index.get_level_values(0).unique())
        layers = list(self.observed.index.get_level_values(1).unique())

        for loc in locs:
            for layer in layers:

                # get observed data for location and layer
                loc_obs = self.observed.loc[idxx[loc, layer, :], :].dropna().reset_index().drop(
                    ['locs', 'layer'], axis=1)
                loc_obs.kstpkper = loc_obs.kstpkper.apply(lambda x: x[1])
                loc_obs.set_index('kstpkper', inplace=True)
                if len(loc_obs) == 0:
                    continue

                self.add_scattergl(
                    name=f'{loc} observed',
                    x=loc_obs.index,
                    y=loc_obs.elev,
                    line_color='blue'
                )

                # get simulated data for location and layer
                sim_obs = self.simulated.loc[idxx[loc, layer, :], :].dropna().reset_index().drop(
                    ['locs', 'layer'], axis=1)
                sim_obs.kstpkper = sim_obs.kstpkper.apply(lambda x: x[1])
                sim_obs.set_index('kstpkper', inplace=True)

                self.add_scattergl(
                    name=f'{loc} simulated',
                    x=sim_obs.index,
                    y=sim_obs.elev,
                    line_color='red'
                )



if __name__ == '__main__':
    fig = CalibrationPlot()

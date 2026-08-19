"""Derived per-cell maps over a finished PRT (particle-tracking) run.

PRT results are trajectories, not a per-cell field, so they do not arrive in the
view grammar for free. This module turns a run's track records into the three
per-cell summaries that *are* grammar-compatible -- travel time, endpoint counts,
and capture zones -- and exposes each as a noun on
:class:`~myflopy.modflow.mf6.prt.PRTRunResults`::

    results.travel_time.get() / .summary() / .plot() / .map() / .mosaic()
    results.endpoints.get()   / ...
    results.capture.get(by="release_group") / .map()

The maps are **time-integrated**: they summarize a whole run rather than one
stress period, so they carry no period footer and reject ``per=`` outright
rather than accept a selector they would ignore.

``results.pathlines`` is the fourth noun and the odd one out: it keeps the
trajectories instead of collapsing them, drawing one map polyline per particle
over an optional base map (:class:`PRTPathlineView`). Its ``get()`` is also the
run's normalized record table.
"""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from flopy.mf6.mfbase import MFDataException

from myflopy import viz
from myflopy._logging import get_logger
from myflopy.modflow.mf6.package_plotting import (
    SpatialView,
    _apply_backend,
    build_cell_input_map_payload,
)
from myflopy.modflow.mf6.package_tables import summarize_input_table
from myflopy.modflow.utils.datatypes.hover import (
    HoverContext,
    pathline_hover,
    result_hover,
)

if TYPE_CHECKING:
    from myflopy.modflow.mf6.prt import PRTRunResults

logger = get_logger(__name__)

#: Track-CSV columns that together identify one released particle.
_PARTICLE_KEYS = ("imdl", "iprp", "irpt", "trelease")

#: ``ireason`` code MF6 writes on a particle's terminating record.
_TERMINATION_REASON = 3

#: Sequential scale for the derived PRT maps (counts and elapsed times are both
#: one-sided magnitudes -- never a signed difference, so never diverging).
PRT_COLORSCALE = "earth"

#: Choropleth options ``_blank_base`` sets itself (so ``base=None`` cannot take them).
_BLANK_BASE_KWARGS = frozenset(
    {"type", "custom_zs", "custom_hover", "hover_heads", "hover_ks", "showscale"}
)

#: Short unit labels for MF6's TDIS time units, for the travel-time hover.
_TIME_UNIT_LABELS = {
    "SECONDS": "s",
    "MINUTES": "min",
    "HOURS": "h",
    "DAYS": "d",
    "YEARS": "yr",
}


def _time_unit(model) -> str:
    """The flow model's TDIS time unit as a short label (``""`` when unknown).

    Travel times are in whatever unit the flow model's TDIS declares, so the
    hover reads it off the model rather than assuming days.
    """

    from myflopy.modflow.mf6.prt import _tdis_time_units

    try:
        declared = str(_tdis_time_units(model)).upper()
    except (AttributeError, MFDataException, OSError, ValueError):
        # A bare/partial model has no TDIS to read (AttributeError). The rest
        # come from `model.sim` being a LAZY property on a file-backed model:
        # the first access performs the whole MFSimulation.load, which raises
        # flopy's MFDataException, OSError on an unreadable workspace, or a
        # UnicodeDecodeError (a ValueError) on a binary mfsim.nam. An unlabelled
        # travel time is better than a map that will not draw.
        logger.debug("no TDIS time unit for the PRT hover", exc_info=True)
        return ""
    return _TIME_UNIT_LABELS.get(declared, "")


def pathline_cell_table(pathlines: pd.DataFrame, *, ncpl: int) -> pd.DataFrame:
    """Normalize raw PRT track records onto the flow model's grid.

    The track CSV numbers cells with ``icell``, a **one-based whole-grid node
    number**; every myflopy view works in zero-based per-layer cells. This is the
    single place that conversion happens::

        cell  = (icell - 1) % ncpl
        layer = (icell - 1) // ncpl      # == ilay - 1

    Parameters
    ----------
    pathlines
        Raw track records (``PRTRunResults.track_records``).
    ncpl
        Cells per layer in the flow model's grid.

    Returns
    -------
    pandas.DataFrame
        Every input column plus ``cell``, ``layer``, ``travel_time``
        (``t - trelease``), ``release_group`` (the PRP boundname MF6 echoes into
        ``name``, uppercased; empty when the run had no groups), and ``particle``
        -- a stable label joining the ``(imdl, iprp, irpt, trelease)`` key.
    """

    missing = [column for column in ("icell", "t", "trelease") if column not in pathlines.columns]
    if missing:
        raise KeyError(
            f"PRT track records are missing required column(s) {missing}; expected a "
            "MODFLOW 6 track CSV (kper, kstp, imdl, iprp, irpt, ilay, icell, ...)."
        )

    frame = pathlines.copy()
    if frame.empty:
        for column in ("cell", "layer"):
            frame[column] = pd.Series(dtype="int64")
        frame["travel_time"] = pd.Series(dtype=float)
        for column in ("release_group", "particle"):
            frame[column] = pd.Series(dtype=object)
        return frame

    node = frame["icell"].astype("int64") - 1
    frame["cell"] = node % int(ncpl)
    frame["layer"] = node // int(ncpl)
    frame["travel_time"] = frame["t"].astype(float) - frame["trelease"].astype(float)
    frame["release_group"] = (
        frame["name"].fillna("").astype(str) if "name" in frame.columns else ""
    )
    frame["particle"] = _particle_labels(frame)
    return frame


def _particle_labels(frame: pd.DataFrame) -> pd.Series:
    """A stable per-particle label from whichever id columns the CSV carries."""

    keys = [column for column in _PARTICLE_KEYS if column in frame.columns]
    if not keys:
        return frame.index.astype(str)
    label = None
    for column in keys:
        piece = frame[column].map(lambda value: format(value, "g"))
        label = piece if label is None else label + "-" + piece
    return label


def particle_endpoint_table(pathlines: pd.DataFrame, *, ncpl: int) -> pd.DataFrame:
    """One row per particle: where it stopped, how long it took, who released it.

    Endpoints are the ``ireason == 3`` (termination) records. When a run has no
    terminating records -- every particle still in transit when tracking stopped --
    each particle's **last** record stands in for its endpoint, so the maps show
    where particles had reached rather than nothing at all.
    """

    frame = pathline_cell_table(pathlines, ncpl=ncpl)
    if frame.empty:
        return frame
    if "ireason" in frame.columns:
        terminal = frame.loc[frame["ireason"] == _TERMINATION_REASON]
        if terminal.empty:
            terminal = frame
    else:
        terminal = frame
    return (
        terminal.sort_values("t", kind="stable")
        .drop_duplicates(subset="particle", keep="last")
        .reset_index(drop=True)
    )


def _aggregate_endpoints(frame: pd.DataFrame, keys: Sequence[str], stat: str) -> pd.DataFrame:
    """Collapse endpoint rows onto ``keys`` with travel-time stats and a particle count."""

    columns = ["travel_time", "particle_count", "min_time", "max_time", "release_groups"]
    if frame.empty:
        return pd.DataFrame(columns=[*keys, *columns])
    grouped = frame.groupby(list(keys), dropna=False)
    times = grouped["travel_time"]
    summary = pd.DataFrame(
        {
            "travel_time": times.agg(stat),
            "particle_count": grouped.size().astype(int),
            "min_time": times.min(),
            "max_time": times.max(),
            "release_groups": grouped["release_group"].agg(_join_groups),
        }
    ).reset_index()
    return summary


def _join_groups(values) -> str:
    """Join the distinct non-empty release-group labels behind one cell."""

    labels = sorted({str(value) for value in values if str(value)})
    return ", ".join(labels)


class _PRTDerivedView(SpatialView):
    """Shared plumbing for the views derived from a PRT run.

    Hosts a :class:`SpatialView` over a *derived* table (the
    ``SurfaceWaterExchangeResultsExplorer`` pattern): the view supplies the
    tabular ``get()`` and an atomic ``map()``, and inherits ``mosaic``/``animate``.
    The three per-cell subclasses use that inheritance as-is;
    :class:`PRTPathlineView` is not per-cell and overrides both.
    PRT runs have no stress-period axis, so ``per`` is fixed and the series verb
    is replaced by a purpose-built distribution figure on each subclass.
    """

    #: Extra per-cell columns offered to the map hover.
    _HOVER_FIELDS: tuple[str, ...] = ()

    def __init__(self, results: PRTRunResults):
        """Bind the view to a finished PRT run's results handle."""

        self.results = results
        self.model = results.flow_model

    @property
    def _ncpl(self) -> int:
        """Cells per layer of the flow model the particles were tracked through."""

        return int(self.model.vor.ncpl)

    def _endpoints(self, layer: int | None = None) -> pd.DataFrame:
        """The run's particle endpoints, optionally restricted to one layer."""

        frame = particle_endpoint_table(self.results.track_records, ncpl=self._ncpl)
        if layer is not None and not frame.empty:
            frame = frame.loc[frame["layer"] == int(layer)]
        return frame

    # -- SpatialView hooks -------------------------------------------------
    def _spatial_map(self, *, per, layer, model=None, **kwargs):
        """Facet panels by layer only -- these maps integrate over the whole run."""

        return self.map(layer=int(layer), **kwargs)

    def _spatial_periods(self) -> list[int]:
        """PRT maps have no period axis (one time-integrated panel)."""

        return [0]

    def _series_table(self) -> pd.DataFrame:
        """PRT results carry no per-period series; :meth:`plot` draws a distribution."""

        raise NotImplementedError(
            "PRT results have no stress-period axis, so there is no series to facet. "
            "Use plot() for the distribution across particles."
        )

    def animate(self, *args, **kwargs):
        """PRT maps integrate the whole run, so there is no axis to animate over."""

        raise NotImplementedError(
            "PRT maps summarize a whole tracking run, so there is no period (or "
            "model) axis to animate over -- a one-frame 'Period 0' animation would "
            "assert a period the values do not belong to. Use mosaic() to compare "
            "layers, capture.map() for per-group panels, or the pathline map for "
            "particle motion."
        )

    def _cell_map(
        self,
        frame: pd.DataFrame,
        *,
        value_column: str,
        layer: int | None,
        hover_spec,
        colorscale: str | None = None,
        logscale: bool = False,
        fill_value: float = float("nan"),
        backend: str = "plotly",
        **kwargs,
    ):
        """Render one pre-aggregated per-cell frame through the house choropleth path.

        ``logscale`` transforms the mapped values only -- the hover keeps the real
        travel times/counts, so a log-scaled map is still readable in real units.
        """

        # ``per``/``agg`` are set below; without this guard they would either
        # collide inside cor() with an opaque "multiple values" TypeError or be
        # swallowed by its **kwargs, silently ignoring what the caller asked for.
        if "per" in kwargs:
            raise TypeError(
                "PRT maps are time-integrated over the whole run, so they take no "
                "per=. Select with layer= (and stat= on travel_time) instead."
            )
        if "agg" in kwargs:
            raise TypeError(
                "PRT maps aggregate particles, not stress-period records: pass "
                "stat= (e.g. stat='max') rather than agg=."
            )
        values, hover = build_cell_input_map_payload(
            frame,
            ncpl=self._ncpl,
            value_column=value_column,
            per=None,
            layer=layer,
            fill_value=fill_value,
            agg="first",  # the frame is already one row per cell
        )
        kwargs.setdefault("hover_spec", hover_spec)
        choro = self.model.plot.map(
            per=0,
            layer=0 if layer is None else int(layer),
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or PRT_COLORSCALE,
            logscale=logscale,
            **kwargs,
        )
        return _apply_backend(choro, backend)

    def _distribution_figure(
        self, labels_values, *, title: str, x_title: str, y_title: str, colors=None
    ):
        """A bar figure of ``[(label, value)]`` in the house template.

        ``colors`` maps label -> color for bars that name a *category* (a release
        group); cells and other ad-hoc labels stay on the default bar color.
        """

        fig = viz.Fig()
        labels = [str(label) for label, _ in labels_values]
        values = [value for _, value in labels_values]
        marker = {"color": [colors[label] for label in labels]} if colors else None
        fig.add_bar(x=labels, y=values, name=y_title, marker=marker)
        fig.update_layout(title=title, xaxis_title=x_title, yaxis_title=y_title)
        return fig


class PRTTravelTimeView(_PRTDerivedView):
    """Particle travel time summarized onto the cells particles ended in.

    The classic time-of-travel figure: for every cell where particles stopped,
    how long they took to get there (``t - trelease``). ``stat`` picks the
    statistic across the particles sharing a cell.
    """

    value_name = "travel_time"
    _HOVER_FIELDS = ("particle_count", "min_time", "max_time", "release_groups")

    def get(self, *, stat: str = "median", layer: int | None = None) -> pd.DataFrame:
        """Travel-time rows, one per ``(layer, cell)`` particles terminated in.

        Columns: ``layer``, ``cell``, ``travel_time`` (the ``stat`` across that
        cell's particles), ``particle_count``, ``min_time``, ``max_time``,
        ``release_groups``. ``layer`` filters; ``None`` keeps every layer.
        """

        return _aggregate_endpoints(self._endpoints(layer), ["layer", "cell"], str(stat))

    def summary(self) -> pd.DataFrame:
        """Compact digest of the travel-time table."""

        return summarize_input_table(
            self.get(), label="prt.travel_time", value_columns=["travel_time"]
        )

    def map(
        self,
        *,
        stat: str = "median",
        layer: int | None = None,
        logscale: bool = False,
        colorscale: str | None = None,
        units: str | None = None,
        backend: str = "plotly",
        **kwargs,
    ):
        """Travel time per cell in plan view.

        ``layer=None`` (default) pools every layer into one plan-view panel --
        the statistic is recomputed over the pooled particles, never averaged
        from per-layer values. Cells no particle reached are left blank rather
        than drawn as zero. ``logscale`` helps when travel times span decades
        (it scales the colors only; the hover keeps real times). ``units`` labels
        the hover values and defaults to the flow model's TDIS time unit.
        """

        frame = _aggregate_endpoints(self._endpoints(layer), ["cell"], str(stat))
        unit = _time_unit(self.model) if units is None else str(units)
        hover = result_hover(
            "travel_time",
            title=f"Travel time ({stat})",
            extra_fields=self._HOVER_FIELDS,
            units={"travel_time": unit, "min_time": unit, "max_time": unit},
            labels={
                "travel_time": f"travel time ({stat})",
                "particle_count": "particles",
                "min_time": "fastest",
                "max_time": "slowest",
                "release_groups": "released from",
            },
            # time-integrated over the run: a "Period 0" footer would be a lie
            footer=(),
        )
        return self._cell_map(
            frame,
            value_column="travel_time",
            layer=layer,
            hover_spec=hover,
            colorscale=colorscale,
            logscale=logscale,
            backend=backend,
            **kwargs,
        )

    def plot(self, *, backend: str = "plotly", title: str | None = None):
        """Cumulative arrival curve: particles that have stopped by a given travel time.

        One step line per release group when the run has them, else one for the
        whole population. ``backend="mpl"`` returns a matplotlib figure.
        """

        frame = self._endpoints()
        heading = title or "Particle arrivals by travel time"
        series = []
        if not frame.empty:
            groups = (
                frame.groupby("release_group", dropna=False)
                if frame["release_group"].astype(bool).any()
                else [("particles", frame)]
            )
            for label, sub in groups:
                times = np.sort(sub["travel_time"].to_numpy(dtype=float))
                counts = np.arange(1, times.size + 1)
                series.append((str(label) or "particles", times, counts))

        # Shared policy colors, so a release group reads the same here as on the
        # pathline map and the capture panels.
        colors = viz.category_colors([label for label, _, _ in series])
        if self._normalize_backend(backend) == "plotly":
            fig = viz.Fig()
            for label, times, counts in series:
                fig.add_scatter(
                    x=times,
                    y=counts,
                    mode="lines+markers",
                    line_shape="hv",
                    name=label,
                    line_color=colors[label],
                    marker_color=colors[label],
                )
            fig.update_layout(
                title=heading, xaxis_title="Travel time", yaxis_title="Particles arrived"
            )
            return fig
        fig, axis = viz.mpl_axes(figsize=(8, 4))
        for label, times, counts in series:
            axis.step(
                times, counts, where="post", linewidth=2.0, label=label, color=colors[label]
            )
        axis.set_title(heading)
        axis.set_xlabel("Travel time")
        axis.set_ylabel("Particles arrived")
        if series:
            axis.legend()
        fig.tight_layout()
        return fig


class PRTPathlineView(_PRTDerivedView):
    """The trajectories themselves: particle tracks drawn over the model map.

    The other PRT nouns collapse a run onto cells; this one keeps the paths. Each
    particle becomes one map polyline, hovering a vertex reports where that
    particle was and how long it had been travelling, and release groups share
    one legend entry and one color with every other PRT figure.

    ``get()`` returns the normalized track records -- every raw track-CSV column
    plus ``cell``/``layer``/``travel_time``/``release_group``/``particle`` -- so
    it is also the way to reach the run's raw table.
    """

    value_name = "travel_time"

    #: Particles drawn before :meth:`map` starts sampling (see ``max_particles``).
    DEFAULT_MAX_PARTICLES = 250

    def get(self, *, group: str | None = None, layer: int | None = None) -> pd.DataFrame:
        """Normalized track records, optionally limited to one release group or layer.

        ``layer`` filters *records*, not particles: a particle that crosses layers
        contributes only the vertices inside the requested layer.
        """

        frame = pathline_cell_table(self.results.track_records, ncpl=self._ncpl)
        if group is not None and not frame.empty:
            frame = frame.loc[frame["release_group"] == str(group)]
        if layer is not None and not frame.empty:
            frame = frame.loc[frame["layer"] == int(layer)]
        return frame.reset_index(drop=True)

    def summary(self) -> pd.DataFrame:
        """One row per particle: where it started, where it stopped, how long it took."""

        frame = self.get()
        columns = [
            "particle",
            "release_group",
            "records",
            "start_cell",
            "start_layer",
            "end_cell",
            "end_layer",
            "travel_time",
            "terminated",
        ]
        if frame.empty:
            return pd.DataFrame(columns=columns)
        ordered = frame.sort_values("t", kind="stable")
        grouped = ordered.groupby("particle", sort=True)
        first, last = grouped.first(), grouped.last()
        digest = pd.DataFrame(
            {
                "release_group": last["release_group"],
                "records": grouped.size().astype(int),
                "start_cell": first["cell"],
                "start_layer": first["layer"],
                "end_cell": last["cell"],
                "end_layer": last["layer"],
                "travel_time": last["travel_time"],
                "terminated": (
                    grouped["ireason"].agg(lambda codes: bool((codes == _TERMINATION_REASON).any()))
                    if "ireason" in frame.columns
                    else False
                ),
            }
        ).reset_index()
        return digest[columns]

    def groups(self) -> list[str]:
        """The distinct release-group labels in the run (empty when it had none)."""

        frame = self.get()
        if frame.empty:
            return []
        return sorted({label for label in frame["release_group"].astype(str) if label})

    # -- the map ------------------------------------------------------------
    def _selected_particles(self, frame: pd.DataFrame, max_particles: int | None):
        """Particle labels to draw, and how many were left out.

        Sampling is **stratified by release group** and then evenly spaced inside
        each group, rather than "the first N": a cap that silently dropped a whole
        capture zone would change what the figure appears to say. A cap smaller
        than the number of groups cannot keep them all -- the quota is round-robin,
        so it keeps the first ``max_particles`` groups in sorted order.
        """

        labels = sorted(set(frame["particle"].tolist()))
        if max_particles is None or len(labels) <= int(max_particles):
            return labels, 0

        by_group: dict[str, set] = {}
        for particle, group in zip(frame["particle"], frame["release_group"], strict=False):
            by_group.setdefault(str(group), set()).add(particle)
        ordered = {group: sorted(members) for group, members in sorted(by_group.items())}

        keep = int(max_particles)
        quotas = dict.fromkeys(ordered, 0)
        while sum(quotas.values()) < keep and any(
            quotas[group] < len(members) for group, members in ordered.items()
        ):
            for group, members in ordered.items():
                if sum(quotas.values()) >= keep:
                    break
                if quotas[group] < len(members):
                    quotas[group] += 1

        chosen: list[str] = []
        for group, members in ordered.items():
            quota = quotas[group]
            if not quota:
                continue
            picks = np.linspace(0, len(members) - 1, quota).round().astype(int)
            chosen.extend(members[index] for index in dict.fromkeys(picks.tolist()))
        chosen.sort()
        return chosen, len(labels) - len(chosen)

    def _trace_color_key(self, color: str) -> str:
        """Resolve ``color=`` to the frame column that decides trace colors."""

        if color in ("release_group", "particle"):
            return color
        raise ValueError(
            f"color={color!r} is not a pathline coloring; use 'release_group' "
            "(default; falls back to per-particle when the run has no groups) or "
            "'particle'."
        )

    def _pathline_traces(self, frame: pd.DataFrame, *, color: str, width: float):
        """One :class:`plotly.graph_objects.Scattermap` polyline per particle."""

        key = self._trace_color_key(color)
        if key == "release_group" and not frame["release_group"].astype(bool).any():
            key = "particle"  # no boundnames in this run: every particle its own color
        unit = _time_unit(self.model)
        hover_fields = [
            name for name in ("layer", "z", "release_group") if name in frame.columns
        ]
        # Label first, THEN color: a run can mix named and un-named release points
        # (MF6 writes an empty name), and those particles are labelled "particles".
        labels = frame[key].astype(str).map(lambda value: value or "particles")
        # Particle ids are not a category that recurs across figures, so they are
        # colored per figure instead of filling the process-wide memo.
        colors = viz.category_colors(labels.unique().tolist(), memoize=key != "particle")

        traces, legended = [], set()
        for label, records in frame.groupby("particle", sort=True):
            records = records.sort_values("t", kind="stable")
            lon, lat = self.model.vor.points_to_latlon(records["x"], records["y"])
            legend_label = labels.loc[records.index[0]]
            spec = pathline_hover(particle=str(label), unit=unit or None, fields=hover_fields)
            context = HoverContext(
                ncpl=len(records),
                payload={name: records[name].tolist() for name in [*hover_fields, "travel_time"]},
                cells=records["cell"].tolist(),
            )
            customdata, template, hoverlabel = spec.render(context)
            traces.append(
                go.Scattermap(
                    mode="lines+markers",
                    lon=lon.tolist(),
                    lat=lat.tolist(),
                    line={"color": colors[legend_label], "width": float(width)},
                    marker={"size": 5, "color": colors[legend_label]},
                    name=legend_label,
                    legendgroup=legend_label,
                    showlegend=legend_label not in legended,
                    customdata=customdata,
                    hovertemplate=template,
                    hoverlabel=hoverlabel,
                )
            )
            legended.add(legend_label)
        return traces

    def _blank_base(self, **kwargs):
        """A map with the grid framed but no cell values -- context for ``base=None``."""

        return self.model.plot.map(
            type="custom",
            custom_zs=[float("nan")] * self._ncpl,
            custom_hover={},
            hover_heads=False,
            hover_ks=False,
            showscale=False,
            **kwargs,
        )

    def map(
        self,
        *,
        base: str | None = "heads",
        per: int | None = None,
        layer: int | None = None,
        group: str | None = None,
        color: str = "release_group",
        width: float = 2.0,
        max_particles: int | None = DEFAULT_MAX_PARTICLES,
        backend: str = "plotly",
        title: str | None = None,
        **base_kwargs,
    ):
        """Draw the particle tracks over a plan-view map of the flow model.

        Parameters
        ----------
        base
            What to draw beneath the paths: ``"heads"`` (the flow model's head
            choropleth at ``per``/``layer``), ``None`` (the grid framed on the
            basemap, no cell values), or an existing ``Choro`` -- so paths can
            ride on any map you have already built, including
            ``capture.map(group=...)`` or ``travel_time.map()``.
        per, layer
            Stress period and layer **of the base map**. The paths themselves are
            never filtered by layer: a particle's plan-view track is the whole
            three-dimensional trajectory, and hiding the parts that left one
            layer would draw a broken line.
        group
            Restrict to one release group (raises naming the available groups if
            it is not one of them).
        color
            ``"release_group"`` (default) or ``"particle"``. Colors come from
            :func:`myflopy.viz.category_colors`, so a group matches its arrival
            curve and capture bars.
        width
            Line width of each track, in pixels.
        max_particles
            Draw at most this many particles, sampled **stratified by release
            group** (a round-robin quota, then evenly spaced inside each group) so
            a cap at or above the group count cannot drop a whole capture zone;
            ``None`` draws every one. The
            cap is announced in the title and by a warning. One trace per particle
            is what makes per-particle hover and legend toggling work, so a run
            with thousands of them is capped rather than silently slow.
        backend
            ``"plotly"`` returns the ``Choro`` carrying the path overlays;
            ``"mpl"`` returns the matplotlib ``(fig, ax)`` plan view. The mpl path
            is the older FloPy plan view and honors only ``group``/``title``:
            ``base``, ``per``, ``layer``, ``color``, ``width`` and
            ``max_particles`` are plotly-side concerns and are ignored.
        title
            Figure title. On a caller-supplied ``base`` the existing title is left
            alone unless this is given.
        **base_kwargs
            Forwarded to the base map (``model.plot.map(...)``) -- e.g. ``contours=``,
            ``zmin=``/``zmax=``, ``locs=``.

        Returns
        -------
        Choro or tuple
            The ``Choro`` carrying one polyline overlay per particle (plotly), or
            matplotlib's ``(fig, ax)``. Passing a ``Choro`` as ``base`` **adds the
            overlays to that map** and returns it, so calling twice with the same
            base draws the paths twice.
        """

        self._trace_color_key(color)  # validate up front, not only if rows survive
        frame = self.get(group=group)
        if group is not None and frame.empty:
            available = self.groups()
            raise KeyError(
                f"No pathlines for release group {group!r}; this run has "
                f"{available or 'no groups (release with group=... to name them)'}."
            )
        if self._normalize_backend(backend) == "mpl":
            from myflopy.modflow.mf6.interactive_plotting import plot_particle_pathlines

            if base_kwargs:
                raise TypeError(
                    f"backend='mpl' draws FloPy's plan view, which takes no base-map "
                    f"options: {sorted(base_kwargs)}. Drop them, or use the default "
                    "plotly backend."
                )
            return plot_particle_pathlines(
                self.model,
                frame,
                title=title or "Particle pathlines",
            )

        kept, dropped = self._selected_particles(frame, max_particles)
        if dropped:
            drawn = frame.loc[frame["particle"].isin(kept)]
            warnings.warn(
                f"Drawing {len(kept)} of {len(kept) + dropped} particles; raise "
                "max_particles= (or None) to draw them all.",
                stacklevel=2,
            )
        else:
            drawn = frame

        borrowed = False
        if base is not None and not isinstance(base, str) and hasattr(base, "add_overlay"):
            # The map is already built; anything that would have shaped it is a
            # selector this call cannot honor, so say so rather than drop it.
            ignored = sorted(base_kwargs) + [
                name for name, value in (("per", per), ("layer", layer)) if value is not None
            ]
            if ignored:
                raise TypeError(
                    f"base= is an already-built map, so {ignored} cannot apply to it. "
                    "Set those when you build the base, or pass base='heads'/None."
                )
            choro, borrowed = base, True
        elif base is None or base == "heads":
            layer_index = 0 if layer is None else int(layer)
            if base is None:
                clashes = sorted(base_kwargs.keys() & _BLANK_BASE_KWARGS)
                if clashes:
                    raise TypeError(
                        f"base=None draws the grid with no cell values, so it sets "
                        f"{clashes} itself. Pass base='heads' (or a Choro) to control them."
                    )
                choro = self._blank_base(per=per, layer=layer_index, **base_kwargs)
            else:
                choro = self.model.plot.map(per=per, layer=layer_index, **base_kwargs)
        elif isinstance(base, str):
            raise ValueError(
                f"base={base!r} is not a pathline base map; use 'heads', None, "
                "or a Choro you built yourself (e.g. capture.map(group=...))."
            )
        else:
            raise TypeError(
                f"Cannot draw pathlines over a base of type {type(base).__name__}; "
                "pass 'heads', None, or a Choro."
            )

        if not drawn.empty:
            choro.add_overlay(*self._pathline_traces(drawn, color=color, width=width))
        # A caller-supplied base already says what it is; only name the figure when
        # asked to, or when this call built the map itself.
        heading = title if title is not None else (None if borrowed else "Particle pathlines")
        if dropped and heading is not None:
            heading = f"{heading} ({len(kept)} of {len(kept) + dropped} particles)"
        if heading is not None:
            choro.fig.update_layout(title=heading)
        return _apply_backend(choro, backend)

    def mosaic(
        self,
        *,
        by: str = "release_group",
        ncols: int = 2,
        sync_views: bool = True,
        backend: str = "plotly",
        title: str | None = None,
        **kwargs,
    ):
        """One map panel per release group, on a shared view.

        Only ``by="release_group"`` is offered: layer panels would repeat the same
        trajectories over different head layers, since a plan-view track is never
        layer-specific.
        """

        if by != "release_group":
            raise ValueError(
                f"by={by!r} is not a pathline facet; particles are not per-layer, "
                "so only by='release_group' faceting is meaningful."
            )
        if self._normalize_backend(backend) != "plotly":
            raise ValueError(
                "A pathline mosaic is a Plotly composition; pass group= to map() for "
                "a single matplotlib panel, or backend='plotly'."
            )
        borrowed_base = kwargs.get("base")
        if borrowed_base is not None and hasattr(borrowed_base, "add_overlay"):
            raise ValueError(
                "One already-built base cannot back several panels -- every group "
                "would be drawn onto the same map. Pass base='heads' or None, or "
                "build the panels yourself with map(group=..., base=...)."
            )
        labels = self.groups()
        if not labels:
            raise ValueError(
                "This run has no release groups to facet by. Release with "
                "PRTReleasePoints.from_cells(..., group=...) to name them, or call "
                "map() for the single-panel figure."
            )
        panels = [
            (label, self.map(group=label, backend=backend, title=label, **kwargs))
            for label in labels
        ]
        return viz.mosaic(panels, ncols=ncols, title=title, sync_views=sync_views)

    def plot(self, *, backend: str = "plotly", title: str | None = None):
        """Elevation against travel time -- the vertical half of each trajectory.

        The plan-view map says where particles went; this says how deep they got
        and when, one line per particle colored by release group.
        """

        frame = self.get()
        heading = title or "Particle elevation over travel time"
        if not frame.empty and "z" not in frame.columns:
            raise KeyError(
                "Track records carry no 'z' column, so there is no elevation to "
                "plot; use map() for the plan view."
            )
        series = []
        if not frame.empty:
            key = "release_group" if frame["release_group"].astype(bool).any() else "particle"
            for _label, records in frame.groupby("particle", sort=True):
                records = records.sort_values("t", kind="stable")
                series.append(
                    (
                        str(records[key].iloc[0]) or "particles",
                        records["travel_time"].to_numpy(dtype=float),
                        records["z"].to_numpy(dtype=float),
                    )
                )
        colors = viz.category_colors([label for label, _, _ in series])

        if self._normalize_backend(backend) == "plotly":
            fig = viz.Fig()
            legended = set()
            for label, times, elevations in series:
                fig.add_scatter(
                    x=times,
                    y=elevations,
                    mode="lines",
                    name=label,
                    legendgroup=label,
                    showlegend=label not in legended,
                    line_color=colors[label],
                )
                legended.add(label)
            fig.update_layout(
                title=heading, xaxis_title="Travel time", yaxis_title="Elevation"
            )
            return fig
        fig, axis = viz.mpl_axes(figsize=(8, 4))
        legended = set()
        for label, times, elevations in series:
            axis.plot(
                times,
                elevations,
                linewidth=1.5,
                color=colors[label],
                label=label if label not in legended else None,
            )
            legended.add(label)
        axis.set_title(heading)
        axis.set_xlabel("Travel time")
        axis.set_ylabel("Elevation")
        if legended:
            axis.legend()
        fig.tight_layout()
        return fig

    def section(self, *args, **kwargs):
        """Trajectories are not a per-cell field, so there is no section to slice."""

        raise NotImplementedError(
            "Pathlines are trajectories, not a per-cell field, so they cannot be "
            "sliced into a cross-section. plot() draws elevation against travel "
            "time, which is the vertical view of the same tracks."
        )

    def animate(self, *args, **kwargs):
        """One figure already holds the whole run, so there is no axis to animate."""

        raise NotImplementedError(
            "The pathline map already draws every particle's whole trajectory, so "
            "there is no period (or model) axis to animate over. Use mosaic() to "
            "compare release groups, or plot() for elevation against travel time. "
            "(Animating particle *position* over time would need a time axis this "
            "view does not build -- see the compromise ledger.)"
        )


class PRTEndpointsView(_PRTDerivedView):
    """Where particles ended: termination counts per cell."""

    value_name = "particle_count"
    _HOVER_FIELDS = ("release_groups", "min_time", "max_time")

    def get(self, *, layer: int | None = None) -> pd.DataFrame:
        """Endpoint counts, one row per ``(layer, cell)`` particles terminated in."""

        return _aggregate_endpoints(self._endpoints(layer), ["layer", "cell"], "median")

    def summary(self) -> pd.DataFrame:
        """Compact digest of the endpoint-count table."""

        return summarize_input_table(
            self.get(), label="prt.endpoints", value_columns=["particle_count"]
        )

    def map(
        self,
        *,
        layer: int | None = None,
        colorscale: str | None = None,
        logscale: bool = False,
        backend: str = "plotly",
        **kwargs,
    ):
        """Particle-termination counts per cell (blank where no particle stopped)."""

        frame = _aggregate_endpoints(self._endpoints(layer), ["cell"], "median")
        unit = _time_unit(self.model)
        hover = result_hover(
            "particle_count",
            title="Endpoints",
            extra_fields=self._HOVER_FIELDS,
            # min/max are elapsed travel times, not clock times -- "first/last
            # arrival" would read as absolute dates.
            units={"min_time": unit, "max_time": unit},
            labels={
                "particle_count": "particles",
                "release_groups": "released from",
                "min_time": "fastest",
                "max_time": "slowest",
            },
            footer=(),
        )
        return self._cell_map(
            frame,
            value_column="particle_count",
            layer=layer,
            hover_spec=hover,
            colorscale=colorscale,
            logscale=logscale,
            backend=backend,
            **kwargs,
        )

    def plot(self, *, top: int = 20, backend: str = "plotly", title: str | None = None):
        """Bar figure of the cells that caught the most particles (``top`` of them)."""

        frame = self.get().sort_values("particle_count", ascending=False).head(int(top))
        pairs = [
            (f"L{int(row.layer)} C{int(row.cell)}", int(row.particle_count))
            for row in frame.itertuples()
        ]
        heading = title or "Particle endpoints by cell"
        if self._normalize_backend(backend) == "plotly":
            return self._distribution_figure(
                pairs, title=heading, x_title="Cell", y_title="Particles"
            )
        fig, axis = viz.mpl_axes(figsize=(8, 4))
        axis.bar([label for label, _ in pairs], [value for _, value in pairs])
        axis.set_title(heading)
        axis.set_xlabel("Cell")
        axis.set_ylabel("Particles")
        fig.tight_layout()
        return fig


class PRTCaptureView(_PRTDerivedView):
    """Capture zones: where each release group's particles ended up.

    ``by="release_group"`` (default) groups particles by the PRP boundname their
    release point carried -- build them with
    ``PRTReleasePoints.from_cells(..., group="west_wells")``. ``by="release_point"``
    falls back to one group per release point id (``irpt``) for runs that were
    built without boundnames.
    """

    value_name = "particle_count"
    _BY_COLUMNS = {"release_group": "release_group", "release_point": "irpt"}

    def _group_column(self, by: str) -> str:
        """Resolve the ``by`` selector to the endpoint column that carries the group."""

        try:
            return self._BY_COLUMNS[str(by)]
        except KeyError:
            raise ValueError(
                f"by must be one of {sorted(self._BY_COLUMNS)}, got {by!r}."
            ) from None

    def _grouped_endpoints(self, by: str, layer: int | None = None) -> pd.DataFrame:
        """Endpoint rows with a ``group`` column resolved from the ``by`` selector."""

        column = self._group_column(by)
        frame = self._endpoints(layer)
        if frame.empty:
            return frame.assign(group=pd.Series(dtype=object))
        if column not in frame.columns:
            raise KeyError(f"PRT track records carry no {column!r} column.")
        labels = frame[column].astype(str)
        if column == "release_group" and not labels.str.strip().astype(bool).any():
            raise ValueError(
                "This PRT run has no release-group boundnames, so particles cannot be "
                "grouped by release group. Rebuild the release points with a group "
                "label -- PRTReleasePoints.from_cells(..., group='west_wells') -- or "
                "pass by='release_point' to group by release point id instead."
            )
        return frame.assign(group=labels)

    def get(self, *, by: str = "release_group", layer: int | None = None) -> pd.DataFrame:
        """Capture rows: ``group``, ``layer``, ``cell``, ``particle_count`` (+ travel times)."""

        frame = self._grouped_endpoints(by, layer)
        return _aggregate_endpoints(frame, ["group", "layer", "cell"], "median")

    def summary(self) -> pd.DataFrame:
        """Compact digest of the capture table."""

        return summarize_input_table(
            self.get(), label="prt.capture", value_columns=["particle_count"]
        )

    def groups(self, *, by: str = "release_group") -> list[str]:
        """The release groups present in this run, in sorted order."""

        frame = self._grouped_endpoints(by)
        return sorted({str(value) for value in frame["group"]}) if not frame.empty else []

    def map(
        self,
        *,
        by: str = "release_group",
        group: str | None = None,
        layer: int | None = None,
        colorscale: str | None = None,
        ncols: int = 2,
        sync_views: bool = True,
        backend: str = "plotly",
        title: str | None = None,
        **kwargs,
    ):
        """Capture-zone map(s).

        With ``group=None`` (default) this returns a **mosaic** -- one endpoint
        panel per release group on one shared color scale and (``sync_views``)
        linked pan/zoom, which is how capture zones are read side by side. Name a
        ``group`` to get that single group's ``Choro`` instead (the form the
        facet composers use); ``ncols``/``sync_views``/``title`` describe the
        mosaic and do not apply to that single panel.
        """

        frame = self._grouped_endpoints(by, layer)
        labels = sorted({str(value) for value in frame["group"]}) if not frame.empty else []
        if group is not None:
            return self._group_map(
                frame,
                group=str(group),
                layer=layer,
                colorscale=colorscale,
                backend=backend,
                **kwargs,
            )
        if self._normalize_backend(backend) != "plotly":
            raise ValueError(
                "A capture mosaic is a Plotly composition; pass group= for a single "
                "matplotlib panel, or backend='plotly'."
            )
        panels = [
            (
                label,
                self._group_map(
                    frame, group=label, layer=layer, colorscale=colorscale, **kwargs
                ),
            )
            for label in labels
        ]
        if not panels:
            raise ValueError("This PRT run has no particle endpoints to map.")
        return viz.mosaic(
            panels,
            ncols=int(ncols),
            title=title or f"Capture zones by {by.replace('_', ' ')}",
            sync_views=sync_views,
        )

    def _group_map(self, frame, *, group: str, layer, colorscale, backend="plotly", **kwargs):
        """One release group's endpoint choropleth."""

        selected = frame.loc[frame["group"].astype(str) == str(group)]
        if selected.empty:
            raise KeyError(f"No particles were released from group {group!r}.")
        unit = _time_unit(self.model)
        hover = result_hover(
            "particle_count",
            title=f"Capture: {group}",
            extra_fields=("release_groups", "min_time", "max_time"),
            # min/max are elapsed travel times, not clock times -- "first/last
            # arrival" would read as absolute dates.
            units={"min_time": unit, "max_time": unit},
            labels={
                "particle_count": "particles",
                "release_groups": "released from",
                "min_time": "fastest",
                "max_time": "slowest",
            },
            footer=(),
        )
        return self._cell_map(
            _aggregate_endpoints(selected, ["cell"], "median"),
            value_column="particle_count",
            layer=layer,
            hover_spec=hover,
            colorscale=colorscale,
            backend=backend,
            **kwargs,
        )

    def _spatial_map(self, *, per, layer, model=None, **kwargs):
        """Facet by layer for one named group (a facet needs a single scale)."""

        if kwargs.get("group") is None:
            available = self.groups(by=kwargs.get("by", "release_group"))
            if len(available) != 1:
                raise ValueError(
                    "Faceting capture panels needs one release group; pass "
                    f"group= (available: {available}) or call map() for the "
                    "per-group mosaic."
                )
            kwargs["group"] = available[0]
        return self.map(layer=int(layer), **kwargs)

    def plot(self, *, by: str = "release_group", backend: str = "plotly", title: str | None = None):
        """Bar figure of how many particles each release group contributed."""

        frame = self._grouped_endpoints(by)
        counts = (
            frame.groupby("group", dropna=False).size().sort_values(ascending=False)
            if not frame.empty
            else pd.Series(dtype=int)
        )
        pairs = [(str(label), int(value)) for label, value in counts.items()]
        heading = title or f"Particles by {by.replace('_', ' ')}"
        colors = viz.category_colors([label for label, _ in pairs])
        if self._normalize_backend(backend) == "plotly":
            return self._distribution_figure(
                pairs,
                title=heading,
                x_title="Release group",
                y_title="Particles",
                colors=colors,
            )
        fig, axis = viz.mpl_axes(figsize=(8, 4))
        axis.bar(
            [label for label, _ in pairs],
            [value for _, value in pairs],
            color=[colors[label] for label, _ in pairs],
        )
        axis.set_title(heading)
        axis.set_xlabel("Release group")
        axis.set_ylabel("Particles")
        fig.tight_layout()
        return fig

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
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from myflopy import viz
from myflopy.modflow.mf6.package_explorer_utils import _default_show_layer_elevs
from myflopy.modflow.mf6.package_plotting import (
    SpatialView,
    _apply_backend,
    build_cell_input_map_payload,
)
from myflopy.modflow.mf6.package_tables import summarize_input_table
from myflopy.modflow.utils.datatypes.hover import result_hover

if TYPE_CHECKING:
    from myflopy.modflow.mf6.prt import PRTRunResults

#: Track-CSV columns that together identify one released particle.
_PARTICLE_KEYS = ("imdl", "iprp", "irpt", "trelease")

#: ``ireason`` code MF6 writes on a particle's terminating record.
_TERMINATION_REASON = 3

#: Sequential scale for the derived PRT maps (counts and elapsed times are both
#: one-sided magnitudes -- never a signed difference, so never diverging).
PRT_COLORSCALE = "earth"

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
    except Exception:  # a bare/partial model has no TDIS to read
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
        Raw track records (``PRTRunResults.pathlines``).
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
    """Shared plumbing for the per-cell views derived from a PRT run.

    Hosts a :class:`SpatialView` over a *derived* table (the
    ``SurfaceWaterExchangeResultsExplorer`` pattern): the view supplies the
    tabular ``get()`` and an atomic ``map()``, and inherits ``mosaic``/``animate``.
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

        frame = particle_endpoint_table(self.results.pathlines, ncpl=self._ncpl)
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
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(self.model))
        kwargs.setdefault("hover_spec", hover_spec)
        choro = self.model.cor(
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

    def _distribution_figure(self, labels_values, *, title: str, x_title: str, y_title: str):
        """A bar figure of ``[(label, value)]`` in the house template."""

        fig = viz.Fig()
        labels = [str(label) for label, _ in labels_values]
        values = [value for _, value in labels_values]
        fig.add_bar(x=labels, y=values, name=y_title)
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

        if self._normalize_backend(backend) == "plotly":
            fig = viz.Fig()
            for label, times, counts in series:
                fig.add_scatter(
                    x=times, y=counts, mode="lines+markers", line_shape="hv", name=label
                )
            fig.update_layout(
                title=heading, xaxis_title="Travel time", yaxis_title="Particles arrived"
            )
            return fig
        fig, axis = viz.mpl_axes(figsize=(8, 4))
        for label, times, counts in series:
            axis.step(times, counts, where="post", linewidth=2.0, label=label)
        axis.set_title(heading)
        axis.set_xlabel("Travel time")
        axis.set_ylabel("Particles arrived")
        if series:
            axis.legend()
        fig.tight_layout()
        return fig


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
        if self._normalize_backend(backend) == "plotly":
            return self._distribution_figure(
                pairs, title=heading, x_title="Release group", y_title="Particles"
            )
        fig, axis = viz.mpl_axes(figsize=(8, 4))
        axis.bar([label for label, _ in pairs], [value for _, value in pairs])
        axis.set_title(heading)
        axis.set_xlabel("Release group")
        axis.set_ylabel("Particles")
        fig.tight_layout()
        return fig

"""Plot payload builders and plotting mixins for package explorers."""

from __future__ import annotations

import inspect
from collections.abc import Iterable
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from myflopy.modflow.utils.datatypes.hover import thousands_formatter
from myflopy.viz import (
    Fig,
    _build_frame_figure,
    is_extent_list,
    mpl_axes,
    mpl_colormap_for,
    subplots,
)
from myflopy.viz import mosaic as viz_mosaic

if TYPE_CHECKING:
    pass
from myflopy.modflow.mf6.package_explorer_utils import (
    _aggregate_hover_strings,
    views_nlay,
)
from myflopy.modflow.mf6.period_end import (
    at_last_saved_step,
    at_last_saved_steps,
    last_saved_steps,
)

#: Spellings each renderer answers to. One tuple pair, because the alias set was
#: written out three times (here, `SpatialView._normalize_backend`, and
#: `render_xsections`) and a fourth spelling added to one of them would have
#: reached only the surfaces that happened to share that copy.
PLOTLY_BACKENDS = ("plotly", "interactive")
MPL_BACKENDS = ("mpl", "matplotlib", "static")


def normalize_backend(backend: str) -> str:
    """Canonicalize a backend alias to ``"plotly"`` or ``"mpl"``.

    The single resolver for the renderer switch, at every scope: the free verbs
    in :mod:`myflopy.plot`, the bound namespaces, and the noun grammar all route
    through it, so ``backend="matplotlib"`` means the same thing everywhere and
    an unknown spelling raises the same message.

    Raises
    ------
    ValueError
        Naming both accepted values. A backend that is silently ignored is the
        defect this replaces -- ``plot.section(..., backend="mpl")`` used to
        vanish into a ``**kwargs`` tail and return a Plotly picture.
    """

    value = str(backend).lower()
    if value in PLOTLY_BACKENDS:
        return "plotly"
    if value in MPL_BACKENDS:
        return "mpl"
    raise ValueError(f"backend must be 'plotly' or 'mpl', got {backend!r}.")


def as_mpl_figure(drawn):
    """Normalize whatever a Matplotlib renderer returned to a bare ``Figure``.

    The renderers behind ``backend="mpl"`` disagree about what they hand back:
    ``Choro.plot_mpl`` returns a ``Figure``, the ``figs`` cross-section helper a
    ``(fig, ax)`` tuple, and FloPy's patch renderer (``GridMesh.plot_mpl``) an
    ``Axes``. The grammar contract is one figure at every scope, so the
    normalization lives here once instead of at each of the five call sites that
    were each about to grow their own copy.

    ``.plot_mpl()`` on the picture is left alone: it keeps returning what it
    always returned, so nothing that already calls it moves.
    """

    if isinstance(drawn, tuple):
        drawn = drawn[0]
    figure = getattr(drawn, "figure", None)
    # An `Axes` carries `.figure`; a `Figure` is its own, and matplotlib >= 3.10
    # answers `Figure.figure` with self, so this covers both without isinstance.
    return figure if figure is not None else drawn


def resolve_noun_hover(hover, trace_kwargs: dict):
    """The caller's hover spec, honouring the legacy ``hover_spec=`` spelling.

    ``hover_spec`` is what every noun's ``map()`` forwarded before plan 8.8 named
    ``hover``: the body ended in ``kwargs.setdefault("hover_spec", <default>)``,
    so a caller who passed ``hover_spec=`` kept it and everyone else got the
    default. Undocumented, but reachable, used, and asserted by
    ``test_package_map_hover_spec_is_overridable``.

    Naming ``hover=`` broke it SILENTLY rather than loudly. :class:`Choro` accepts
    both, stores ``hover=`` as the override, and prefers the override -- so the
    noun's own default, now arriving as an explicit ``hover=``, outranked the
    caller's ``hover_spec=``, which still rode the tail and still reached the
    constructor. No error, and no visible difference except the tooltip.

    ``hover_spec`` cannot simply join :data:`NOUN_REFUSED_PARAMS`:
    ``test_the_tiers_partition_the_free_verb`` asserts every declared name is a
    real ``plot.map`` parameter, and ``plot.map`` has only ``hover``. So it stays
    an unnamed legacy alias, resolved here and POPPED off the tail so it cannot
    reach the constructor and lose a second time.

    ``hover=`` wins when both are given -- it is the named, current spelling.
    """

    legacy = trace_kwargs.pop("hover_spec", None)
    return legacy if hover is None else hover


def _explicit_options(local_vars: dict, method, names) -> dict:
    """The named parameters in ``names`` a caller actually SET, by value.

    Lives here rather than with its first caller because the composers need it
    too: ``mosaic``/``animate`` name the same drawing parameters and must forward
    only what the caller actually set, or a default would override the noun's own
    colour policy on every panel.

    ``PRTPathlineView.map`` used to read its ``**base_kwargs`` tail three times
    -- to refuse base-map options on the Matplotlib branch, to refuse them on an
    already-built base, and to forward them -- so naming those parameters (plan
    8.8) would have emptied the tail and silently disabled all three checks. This
    reconstructs the same information from the named parameters instead, by
    comparing each against the default in the method's OWN signature, so the two
    cannot drift.

    A numpy array compared against ``None`` returns an elementwise array whose
    truth value raises; an array is never a default, so that is "supplied".
    """

    defaults = {
        name: parameter.default
        for name, parameter in inspect.signature(method).parameters.items()
    }
    supplied = {}
    for name in names:
        if name not in local_vars:
            continue
        value, default = local_vars[name], defaults.get(name, inspect.Parameter.empty)
        try:
            differs = value is not default and bool(value != default)
        except ValueError:
            differs = True
        if differs:
            supplied[name] = value
    return supplied


def _mpl_shared_scale(choros, *, diff: bool):
    """``(arrays, cmap, norm)`` for a static mosaic/animation of ``Choro`` panels.

    The per-panel values are the DATA-unit ones (``linear_zs``); the shared
    norm is symmetric about zero for a diff surface, a ``LogNorm`` over the
    pooled positive range when every panel is log-scaled, else linear over the
    pooled finite range. The colormap is the first panel's own colorscale
    (``mpl_colormap_for``, from Plotly's stops) so the static composer draws
    what the interactive one draws -- it used to hard-code ``viridis`` over the
    raw log10 exponents, so a static log mosaic read ``-3 ... 2`` under a label
    naming real units (ledger 252). Mixing log and linear panels raises, as
    ``viz.mosaic`` does: one bar cannot mean both.

    An explicit ``zmin``/``zmax`` on the panels WINS over the pooled range, on
    either side independently -- the composers name both and forwarded them to
    every panel, and this pooled over the data regardless, so
    ``hds.mosaic(zmin=100, zmax=120)`` drew 79.04..153.76 (ledger 306).
    """

    from matplotlib import colors as mcolors

    log_flags = {bool(getattr(choro, "logscale", False)) for choro in choros}
    if len(log_flags) > 1:
        raise ValueError(
            "cannot share one color scale between log-scaled and linear panels; "
            "give every panel the same logscale=."
        )
    logscale = log_flags == {True}
    arrays = [np.asarray(choro.linear_zs, dtype=float) for choro in choros]
    if logscale:
        arrays = [np.where(values > 0, values, np.nan) for values in arrays]
    finite_parts = [values[np.isfinite(values)] for values in arrays if np.isfinite(values).any()]
    finite = np.concatenate(finite_parts) if finite_parts else np.asarray([])
    lo, hi = _explicit_color_limits(choros)
    if diff:
        absmax = (float(np.nanmax(np.abs(finite))) if finite.size else 1.0) or 1.0
        return arrays, "RdBu", mcolors.Normalize(
            vmin=-absmax if lo is None else lo, vmax=absmax if hi is None else hi
        )
    cmap = mpl_colormap_for(getattr(choros[0], "colorscale", "earth"))
    if logscale:
        vmin = lo if lo is not None else (float(np.nanmin(finite)) if finite.size else 0.1)
        vmax = hi if hi is not None else (float(np.nanmax(finite)) if finite.size else 10.0)
        return arrays, cmap, mcolors.LogNorm(vmin=vmin, vmax=vmax)
    vmin = lo if lo is not None else (float(np.nanmin(finite)) if finite.size else 0.0)
    vmax = hi if hi is not None else (float(np.nanmax(finite)) if finite.size else 1.0)
    return arrays, cmap, mcolors.Normalize(vmin=vmin, vmax=vmax)


def _explicit_color_limits(choros) -> tuple[float | None, float | None]:
    """The widest caller-set ``(zmin, zmax)`` across ``choros``, in DATA units.

    Each side independently: a ``zmin`` alone leaves the top pooled. Widest,
    because panels a caller built separately may disagree, and a shared bar that
    clips one of them would misreport it.
    """

    lows = [float(c._zmin) for c in choros if getattr(c, "_zmin", None) is not None]
    highs = [float(c._zmax) for c in choros if getattr(c, "_zmax", None) is not None]
    return (min(lows) if lows else None), (max(highs) if highs else None)


def _apply_backend(choro, backend: str = "plotly"):
    """Return a ``Choro`` as an interactive Plotly figure or a static mpl one.

    Lets a single-model atomic ``map(..., backend=)`` honor the same
    ``backend="plotly"|"mpl"`` switch the :class:`SpatialView` grammar uses,
    without the map having to inherit the engine.
    """

    return choro if normalize_backend(backend) == "plotly" else choro.plot_mpl()


def _symmetric_color_limit(values: Iterable[float]) -> float:
    """Return the symmetric absolute color bound for a numeric sequence."""

    array = np.asarray(list(values), dtype=float)
    if array.size == 0:
        return 0.0
    finite = array[np.isfinite(array)]
    if finite.size == 0:
        return 0.0
    return float(np.max(np.abs(finite)))


def _blue_white_red_diverging_colorscale() -> list[list[object]]:
    """Return a blue-white-red diverging colorscale for signed maps.

    Plotly maps position 0.0 to ``zmin`` and 1.0 to ``zmax``, so this is
    blue at the NEGATIVE end and red at the POSITIVE end. Use
    :func:`_exchange_colorscale` for surface-water exchange rather than
    reaching for this directly -- which sign means "gaining" is a per-package
    fact, not a constant.
    """

    return [
        [0.0, "#1f77b4"],
        [0.5, "#ffffff"],
        [1.0, "#d62728"],
    ]


def _red_white_blue_diverging_colorscale() -> list[list[object]]:
    """Return the diff-map orientation as STOPS: red at ``zmin``, blue at ``zmax``.

    The same endpoint colors as :func:`_blue_white_red_diverging_colorscale`,
    reversed -- the "negative red, positive blue" rule
    ``get_default_group_compare_colorscale()`` declares for difference maps. It
    exists as stops rather than as a NAME because the helpers own the
    ORIENTATION -- which end is red is a policy decision a bare name cannot
    carry. (Until ledger 253 there was a second reason: a hand table sent the
    name ``'rdbu'`` to matplotlib's REVERSED colormap, so a named diverging
    scale rendered mirrored between backends. Names now render from Plotly's
    own stops on both, so only the orientation reason remains.)

    Derived from its sibling rather than restating the hexes, so the two
    orientations cannot drift apart.
    """

    colors = [color for _, color in reversed(_blue_white_red_diverging_colorscale())]
    return [[position, color] for position, color in zip((0.0, 0.5, 1.0), colors, strict=True)]


def _exchange_colorscale(frame: str) -> list[list[object]]:
    """Return the exchange scale oriented for ``frame``: gaining blue, losing red.

    ``frame`` is a ``ResultSpec.reference_frame`` value and MUST be supplied --
    which sign means "gaining" is a property of the source file, not a constant,
    because myflopy keeps MF6's raw sign and never normalizes it:

    * ``"gwf"`` -- an aquifer-referenced cell (.cbc) record: MF6 writes flow FROM
      the feature TO the cell, so a GAINING feature is NEGATIVE. Blue sits at the
      negative end -> ``[blue, white, red]``. (SFR and the list BCs.)
    * ``"feature"`` -- a feature-referenced package budget file (LAK's ``GWF``
      record) OR myflopy's normalized ``exchange_intensity``: a gaining feature
      is POSITIVE. Blue sits at the positive end -> ``[red, white, blue]``.

    Plotly maps position 0.0 to ``zmin`` and 1.0 to ``zmax``. Passing a stale
    frame is what inverted the SFR map on 2026-07-19; callers now derive it from
    the registry, which describes what MF6 wrote and so cannot drift from the
    data.
    """

    blue, white, red = "#1f77b4", "#ffffff", "#d62728"
    if frame == "gwf":
        return [[0.0, blue], [0.5, white], [1.0, red]]
    return [[0.0, red], [0.5, white], [1.0, blue]]

def _as_layer_cell_property(values, *, nlay: int, ncpl: int, label: str) -> np.ndarray:
    """Return static package arrays as a consistent ``(nlay, ncpl)`` array."""

    arr = np.asarray(values)
    arr = np.squeeze(arr)
    if arr.ndim == 0:
        return np.full((int(nlay), int(ncpl)), arr.item())
    if arr.ndim == 1:
        if arr.size == int(ncpl):
            return np.tile(arr.reshape(1, int(ncpl)), (int(nlay), 1))
        if arr.size == int(nlay) * int(ncpl):
            return arr.reshape(int(nlay), int(ncpl))
    if arr.ndim == 2:
        if arr.shape == (int(nlay), int(ncpl)):
            return arr
        if arr.shape == (int(ncpl), int(nlay)):
            return arr.T
        if arr.size == int(nlay) * int(ncpl):
            return arr.reshape(int(nlay), int(ncpl))
    if arr.ndim >= 3:
        if arr.shape[0] == int(nlay) and int(np.prod(arr.shape[1:])) == int(ncpl):
            return arr.reshape(int(nlay), int(ncpl))
        if arr.size == int(nlay) * int(ncpl):
            return arr.reshape(int(nlay), int(ncpl))
    raise ValueError(
        f"{label} must resolve to one value per layer/cell; got shape={arr.shape}."
    )


#: per-feature result columns passed through to LAK/SFR q-map hovers when a
#: caller has joined them onto the exchange frame (see ``join_*_stage``).
_FEATURE_HOVER_COLUMNS = ("stage", "depth")


def _rows_drawn_for_period(frame: pd.DataFrame, per: int | None) -> pd.DataFrame:
    """The rows a picture of stress period ``per`` draws: its LAST saved step (ledger 337).

    A results table labelled by its file's own clock (a ``kstpkper`` column) holds
    every step the period saved, and a map reduces the rows of one CELL to one value
    -- which ADDED a flow's steps (3.0x at ``nstp=3`` under ``SAVE BUDGET ALL``) and
    painted a stage's FIRST step (the default ``agg="first"``) until ledger 348. The
    step is chosen here, before any rows combine, through
    :func:`~myflopy.modflow.mf6.period_end.at_last_saved_step`. A table with no
    clock (an INPUT table: one record per period) or no period (``per=None``) is
    returned as given.

    Parameters
    ----------
    frame : pandas.DataFrame
        The selected rows, usually a noun's ``get(per=per, ...)``.
    per : int or None
        The stress period drawn, zero-based; ``None`` draws the rows as given.

    Returns
    -------
    pandas.DataFrame
        The rows of ``per``'s last saved step, or ``frame`` itself.

    Raises
    ------
    ValueError
        When ``frame`` has rows and a clock but none saved in ``per``.
    """

    if per is None or frame.empty or "kstpkper" not in frame.columns:
        return frame
    return at_last_saved_step(frame, per)


def _append_feature_hover_columns(hover, selected, full_index) -> None:
    """Add joined per-feature columns (stage, depth) to a q-map hover payload.

    Stage repeats on every cell a feature touches, so a per-cell mean is the
    faithful aggregate (matches ``StageResultsExplorer``'s series default).
    """

    for column in _FEATURE_HOVER_COLUMNS:
        if column not in selected.columns or column in hover:
            continue
        numeric = pd.to_numeric(selected[column], errors="coerce")
        hover[column] = (
            numeric.groupby(selected["cell"], dropna=False)
            .mean()
            .reindex(full_index, fill_value=np.nan)
            .tolist()
        )


def build_sfr_q_map_payload(
    frame: pd.DataFrame,
    *,
    ncpl: int,
    per: int | None,
    layer: int | None,
    multiplier: float = 1.0,
    fill_value: float = 0.0,
    value_column: str = "q_gwf",
) -> tuple[list[float], dict[str, list]]:
    """Convert SFR exchange rows into a cell choropleth normalized by reach length.

    The mapped value is ``sum(q) / sum(rlen)`` within each cell, which is more
    meaningful than raw exchange magnitude alone when reach lengths vary across
    cells. ``value_column`` is the frame's exchange column (``q_gwf``); it keeps
    MF6's raw aquifer-referenced sign, so gaining reaches are negative. The sums run
    over ONE time step, ``per``'s last saved (``_rows_drawn_for_period``, ledger
    348): they ran over every step the period saved, three flows added at
    ``nstp=3`` under ``SAVE BUDGET ALL``.
    """

    full_index = pd.Index(range(int(ncpl)), name="cell")
    if frame.empty:
        values = pd.Series(float(fill_value), index=full_index, name="q_per_length")
        hover = {
            "Package": ["" for _ in full_index],
            "Period": [per for _ in full_index],
            "Layer": [layer for _ in full_index],
            "Cell": full_index.to_list(),
            "Record Count": [0 for _ in full_index],
            "q_per_length": values.tolist(),
            "q": [0.0 for _ in full_index],
            "rlen": [0.0 for _ in full_index],
        }
        return values.tolist(), hover

    selected = _rows_drawn_for_period(frame, per).copy()
    selected[value_column] = pd.to_numeric(selected[value_column], errors="coerce")
    selected["rlen"] = pd.to_numeric(selected["rlen"], errors="coerce")
    grouped = selected.groupby("cell", dropna=False)
    q_sum = grouped[value_column].sum(min_count=1).reindex(full_index, fill_value=np.nan)
    rlen_sum = grouped["rlen"].sum(min_count=1).reindex(full_index, fill_value=np.nan)
    values = pd.Series(float(fill_value), index=full_index, dtype=float)
    valid = rlen_sum > 0.0
    # Scale the values, never the fill (ledger 473).
    values.loc[valid] = (q_sum.loc[valid] / rlen_sum.loc[valid]).astype(float) * float(multiplier)

    hover: dict[str, list] = {
        "Package": grouped["package"]
        .agg(_aggregate_hover_strings)
        .reindex(full_index, fill_value="")
        .tolist()
        if "package" in selected.columns
        else ["" for _ in full_index],
        "Period": [per for _ in full_index],
        "Layer": [layer for _ in full_index],
        "Cell": full_index.to_list(),
        "Record Count": grouped.size()
        .reindex(full_index, fill_value=0)
        .astype(int)
        .tolist(),
        "q_per_length": values.tolist(),
        "q": q_sum.fillna(0.0).astype(float).tolist(),
        "rlen": rlen_sum.fillna(0.0).astype(float).tolist(),
    }
    if "reach" in selected.columns:
        reach_strings = selected.assign(
            reach=selected["reach"].astype("string")
        ).groupby("cell", dropna=False)["reach"]
        hover["reach"] = (
            reach_strings.agg(_aggregate_hover_strings)
            .reindex(full_index, fill_value="")
            .tolist()
        )
    _append_feature_hover_columns(hover, selected, full_index)
    return values.tolist(), hover


def build_lak_q_map_payload(
    frame: pd.DataFrame,
    *,
    ncpl: int,
    per: int | None,
    layer: int | None,
    multiplier: float = 1.0,
    fill_value: float = 0.0,
    value_column: str = "q_lake",
) -> tuple[list[float], dict[str, list]]:
    """Convert LAK exchange rows into a cell choropleth normalized by area.

    The mapped value is ``sum(q) / sum(flow_area)`` within each cell, giving a
    lake-groundwater exchange intensity with units of length per time.
    ``value_column`` is the frame's exchange column (``q_lake``); it keeps MF6's
    raw feature-referenced sign, so a gaining lake is positive. The sums run over
    ONE time step, ``per``'s last saved (``_rows_drawn_for_period``, ledger 348):
    they ran over every step the period saved, so at ``nstp=3`` under ``SAVE BUDGET
    ALL`` a connection's three flows AND its area were each added, and the hover's
    ``q`` read -2793 ft3/d where the period ends at -587.
    """

    full_index = pd.Index(range(int(ncpl)), name="cell")
    if frame.empty:
        values = pd.Series(float(fill_value), index=full_index, name="q_per_area")
        hover = {
            "Package": ["" for _ in full_index],
            "Period": [per for _ in full_index],
            "Layer": [layer for _ in full_index],
            "Cell": full_index.to_list(),
            "Record Count": [0 for _ in full_index],
            "q_per_area": values.tolist(),
            "q": [0.0 for _ in full_index],
            "flow_area": [0.0 for _ in full_index],
        }
        return values.tolist(), hover

    selected = _rows_drawn_for_period(frame, per).copy()
    selected[value_column] = pd.to_numeric(selected[value_column], errors="coerce")
    selected["flow_area"] = pd.to_numeric(selected["flow_area"], errors="coerce")
    grouped = selected.groupby("cell", dropna=False)
    q_sum = grouped[value_column].sum(min_count=1).reindex(full_index, fill_value=np.nan)
    area_sum = (
        grouped["flow_area"].sum(min_count=1).reindex(full_index, fill_value=np.nan)
    )
    values = pd.Series(float(fill_value), index=full_index, dtype=float)
    valid = area_sum > 0.0
    # Scale the values, never the fill (ledger 473).
    values.loc[valid] = (q_sum.loc[valid] / area_sum.loc[valid]).astype(float) * float(multiplier)

    hover: dict[str, list] = {
        "Package": grouped["package"]
        .agg(_aggregate_hover_strings)
        .reindex(full_index, fill_value="")
        .tolist()
        if "package" in selected.columns
        else ["" for _ in full_index],
        "Period": [per for _ in full_index],
        "Layer": [layer for _ in full_index],
        "Cell": full_index.to_list(),
        "Record Count": grouped.size()
        .reindex(full_index, fill_value=0)
        .astype(int)
        .tolist(),
        "q_per_area": values.tolist(),
        "q": q_sum.fillna(0.0).astype(float).tolist(),
        "flow_area": area_sum.fillna(0.0).astype(float).tolist(),
    }
    if "lake" in selected.columns:
        lake_strings = selected.assign(lake=selected["lake"].astype("string")).groupby(
            "cell", dropna=False
        )["lake"]
        hover["lake"] = (
            lake_strings.agg(_aggregate_hover_strings)
            .reindex(full_index, fill_value="")
            .tolist()
        )
    if "claktype" in selected.columns:
        hover["claktype"] = (
            grouped["claktype"]
            .agg(_aggregate_hover_strings)
            .reindex(full_index, fill_value="")
            .tolist()
        )
    _append_feature_hover_columns(hover, selected, full_index)
    return values.tolist(), hover


def build_surface_water_q_map_payload(
    frame: pd.DataFrame,
    *,
    ncpl: int,
    per: int | None,
    layer: int | None,
    multiplier: float = 1.0,
    fill_value: float = 0.0,
) -> tuple[list[float], dict[str, list]]:
    """Convert combined SFR/LAK exchange rows into one shared cell map payload."""

    full_index = pd.Index(range(int(ncpl)), name="cell")
    if frame.empty:
        values = pd.Series(
            float(fill_value), index=full_index, name="exchange_intensity"
        )
        hover = {
            "Package": ["" for _ in full_index],
            "Period": [per for _ in full_index],
            "Layer": [layer for _ in full_index],
            "Cell": full_index.to_list(),
            "surface_water_exchange": values.tolist(),
            "sfr_exchange": [0.0 for _ in full_index],
            "lak_exchange": [0.0 for _ in full_index],
        }
        return values.tolist(), hover

    selected = frame.copy()
    selected["exchange_intensity"] = pd.to_numeric(
        selected["exchange_intensity"], errors="coerce"
    )
    grouped = selected.groupby("cell", dropna=False)
    values = (
        grouped["exchange_intensity"]
        .sum(min_count=1)
        .reindex(full_index, fill_value=np.nan)
    )
    # Scale the values, never the fill (ledger 473).
    values = (values.astype(float) * float(multiplier)).fillna(float(fill_value))

    sfr_component = selected.loc[selected["source"] == "sfr"].groupby(
        "cell", dropna=False
    )["exchange_intensity"].sum(min_count=1).reindex(full_index, fill_value=0.0).astype(
        float
    ) * float(multiplier)
    lak_component = selected.loc[selected["source"] == "lak"].groupby(
        "cell", dropna=False
    )["exchange_intensity"].sum(min_count=1).reindex(full_index, fill_value=0.0).astype(
        float
    ) * float(multiplier)

    hover: dict[str, list] = {
        "Package": grouped["package"]
        .agg(_aggregate_hover_strings)
        .reindex(full_index, fill_value="")
        .tolist()
        if "package" in selected.columns
        else ["" for _ in full_index],
        "Period": [per for _ in full_index],
        "Layer": [layer for _ in full_index],
        "Cell": full_index.to_list(),
        "Record Count": grouped.size()
        .reindex(full_index, fill_value=0)
        .astype(int)
        .tolist(),
        "surface_water_exchange": values.tolist(),
        "sfr_exchange": sfr_component.tolist(),
        "lak_exchange": lak_component.tolist(),
    }
    if "source" in selected.columns:
        hover["source"] = (
            grouped["source"]
            .agg(_aggregate_hover_strings)
            .reindex(full_index, fill_value="")
            .tolist()
        )
    return values.tolist(), hover


# --- what a NOUN's picture verb offers (plan 8.8, extended to the noun tier) ----
#
# `model.hds.map()` and `model.drn.inputs.elev.map()` are NOUNS: narrower
# than the verb, because the noun already fixes WHAT is drawn. The 8.8 rule --
# name every parameter, because editors read the `def` line and never run the
# module -- applies to them too, and until 2026-08-30 reached none of them.
#
# The split below is MEASURED, not stylistic. Every Tier 1 name changes the
# figure on both a results noun (`model.hds`) and a record noun
# (`packages.drn.inputs.elev`). Every Tier 2 name changes it only on a per-layer
# field: on a record noun the rows are period-indexed with one value per cell, so
# `kstpkper`/`bgs`/`hover_layers`/`hover_surfaces`/`show_layer_elevs` produce a
# byte-identical figure, and `show_mounding` is worse than inert -- it injects a
# head-derived row into an elevation tooltip on every cell. Naming a parameter
# that provably cannot act is the same defect as hiding one that can -- and so is
# FORWARDING one: a one-value noun's tail refuses Tier 2 by name
# (`refuse_noun_parameters`, ledger 465).

#: Drawing parameters every spatial noun offers, on top of its own selectors.
#:
#: `backend` is in this tier deliberately, and it is the one name here that is
#: not a drawing option: it chooses the RENDERER. It has always worked on the
#: nouns -- through their `**kwargs` tail -- so a conversion that named the other
#: twenty and forgot this one would DELETE the matplotlib backend from that noun,
#: silently, and no other test would notice. Listing it here makes
#: `test_noun_signatures` refuse the conversion instead.
NOUN_MAP_PARAMS = (
    "zmin", "zmax", "colorscale", "logscale",
    "contours", "contour_values", "contour_levels", "contour_color",
    "contour_width", "contour_name", "contour_clip", "contour_resolution",
    "contour_method",
    "select", "select_style", "select_color",
    "features", "hillshade_path",
    "fit_bounds", "bounds_padding", "extent",
    "hover", "show_layer_elevs", "hover_surfaces",
    "backend",
)

#: Extra parameters only a PER-LAYER field (heads, concentration, temperature)
#: can honour. Measured inert on record nouns -- see the note above.
#:
#: `show_layer_elevs`/`hover_surfaces` LEFT this tier on 2026-09-14 (ledger 258).
#: They read the GRID's layer frame, not the mapped field, so there was never a
#: reason they could not act on a K or conductance map -- and they did not, only
#: because `Choro._build_hover_context` fetched the elevations inside its heads
#: branch. With that fixed they are measured live on every spatial noun, so they
#: moved to Tier 1. `hover_layers` stayed: it needs a per-layer FIELD, which a
#: record noun does not carry, and it is still measured inert there.
LAYER_FIELD_MAP_PARAMS = (
    "kstpkper", "per_timestep", "bgs",
    "hover_layers", "show_mounding",
)

#: Parameters a noun must REFUSE rather than forward. Each contradicts what a
#: noun IS, and each was reachable and silently wrong:
#:   `values`/`type` -- the noun's whole premise is that it fixes the field.
#:     Measured: `hds.map(values=[1.0]*ncpl)` repainted every cell while the
#:     hover, title and colorscale went on reporting real heads.
#:   `custom_hover`/`hover_heads`/`hover_ks` -- the legacy hover path. A noun
#:     always passes `hover_spec=`, which beats them, so they were dead.
NOUN_REFUSED_PARAMS = ("values", "type", "custom_hover", "hover_heads", "hover_ks")

#: The per-PANEL drawing options the composers name and forward to each panel's
#: `map`. `backend` is excluded deliberately: on a composer it chooses the
#: renderer for the WHOLE figure, and a panel is always a `Choro`.
PANEL_DRAWING_PARAMS = tuple(name for name in NOUN_MAP_PARAMS if name != "backend")


def refuse_panel_drawing_options(kind: str, verb: str, options: dict) -> None:
    """Raise if map-drawing options were passed to a non-map composer kind.

    `mosaic(kind="plot", logscale=True)` used to drop `logscale` into a
    `**kwargs` tail that the series branch never reads -- accepted, ignored, and
    the figure came back linear with nothing said. Naming the parameters (8.8)
    makes the mistake reachable from an editor, so it has to be answered.

    Only the NAMED drawing options are refused; the `**kwargs` tail still rides
    to the section branch, which has drawing options of its own.
    """

    if not options:
        return
    raise TypeError(
        f"{verb}(kind={kind!r}) draws no choropleth, so it cannot honour "
        f"{sorted(options)}. Those are map options -- use kind='map', or drop them."
    )


def refuse_a_missing_line_profile(view, where: str) -> None:
    """Raise unless ``view`` has a line profile (a :class:`LineProfileView`).

    ``mosaic``/``animate`` are defined once, on :class:`SpatialView`, so every
    noun has ``kind="section"``; only a dependent-variable noun can draw it
    (ledger 462).
    """

    if isinstance(view, LineProfileView):
        return
    raise ValueError(
        f"{where} draws a LINE PROFILE, which only a dependent-variable noun has -- "
        "model.hds / .conc / .temp, group.hds, diff.hds. This noun has none; its "
        "pictures are map, plot, mosaic(kind='map'|'plot') and animate. The model's "
        "cross-section, cells coloured by layer or by any (nlay, ncpl) array, is "
        "model.plot.section(line=..., fill=...)."
    )


#: What a ``kind="plot"`` composer reads: the series panels take nothing else,
#: so any other keyword is refused rather than dropped (ledger 461).
SERIES_COMPOSER_PARAMS = {
    "mosaic": ("kind", "by", "per", "layer", "model", "cells", "agg", "ncols", "title",
               "backend"),
    "animate": ("kind", "over", "per", "layer", "model", "cells", "agg", "title",
                "backend"),
}


def refuse_series_composer_options(where: str, verb: str, unexpected: dict) -> None:
    """Raise a ``TypeError`` naming every keyword a ``kind="plot"`` composer did not take.

    The series panels read only :data:`SERIES_COMPOSER_PARAMS`; the composer's
    ``**kwargs`` tail used to swallow anything else and draw (ledger 461).
    """

    if not unexpected:
        return
    import difflib

    accepted = SERIES_COMPOSER_PARAMS[verb]
    described = []
    for name in sorted(unexpected):
        near = difflib.get_close_matches(name, list(accepted), n=1)
        described.append(f"{name!r}" + (f" (did you mean {near[0]!r}?)" if near else ""))
    raise TypeError(
        f"{where} got unexpected keyword argument(s) {', '.join(described)}: its series "
        f"panels take only {', '.join(accepted)}."
    )


#: Reachable on a noun but deliberately not offered:
#:   `hover_fields` -- NOT inert, measured (ledger 436): through a noun's tail it
#:     appends a context field the noun's hover carries (`drn.results.q.map(
#:     hover_fields=["Record Count"])` adds that line). Whether to name it on every
#:     noun or refuse it for `hover=` is an open decision (ledger 465), so it
#:     still rides the tail and acts.
#:   `zoom` -- inert under the default `fit_bounds=True`, which replaces it with
#:     a fitted zoom. (This line used to say it "needs a basemap the noun's
#:     choropleth does not draw"; measured, nouns draw `carto-voyager` and
#:     `fit_bounds=False, zoom=10` through the tail DOES act -- ledger 307.)
#:   `rch_scale` / `animation_kstpkpers` -- recharge- and animate-specific.
NOUN_INERT_PARAMS = ("hover_fields", "zoom", "rch_scale", "animation_kstpkpers")


def refuse_noun_parameters(noun: str, field: str, kwargs: dict, *, honours=()) -> None:
    """Raise if a noun was handed a parameter that contradicts what it draws.

    Called from a noun's ``map()`` before forwarding. Without it these arrive in
    the ``**trace_kwargs`` tail: `type=` collided with a `TypeError` naming
    `ModelPlots.map`, a class the caller never typed, and `values=` won outright,
    producing a figure whose cells and whose tooltips disagreed.

    The per-layer FIELD options (:data:`LAYER_FIELD_MAP_PARAMS`) are refused too,
    unless ``honours`` names them: a noun that draws one value per cell is not
    told apart by them. Measured on eight nouns (ledger 465): ``bgs``,
    ``kstpkper``, ``per_timestep`` and ``hover_layers`` left the figure
    byte-identical, and ``show_mounding`` kept the colours and added a head row
    to every tooltip -- except ``hover_layers`` on the two pictures that keep a
    per-layer hover (the transmissivity nouns and ``model.inputs.uzf``).

    Parameters
    ----------
    noun : str
        How the caller reached the noun, for the message (``"hfb.inputs"``).
    field : str
        What the noun draws, for the message (``"barrier flow"``).
    kwargs : dict
        The noun's ``**trace_kwargs`` tail, read and never modified.
    honours : iterable of str, default ()
        The per-layer field options this noun's picture DOES act on -- every one
        of :data:`LAYER_FIELD_MAP_PARAMS` on a dependent-variable noun or a map
        drawn over one, ``("hover_layers",)`` on a per-layer hover. The rest of
        that tier raises ``TypeError`` when it rides the tail.

    Raises
    ------
    TypeError
        Naming the first refused parameter found and what to call instead.
    """

    for name in LAYER_FIELD_MAP_PARAMS:
        if name in kwargs and name not in honours:
            raise TypeError(
                f"{noun}.map() draws {field}, one value per cell, so it takes no "
                f"{name}=: that option reads a per-layer field -- heads, "
                "concentration, temperature -- and changes nothing here"
                + (" but the tooltip, which it fills with a head-derived row"
                   if name == "show_mounding" else "")
                + f". Map the field it reads instead: model.hds.map({name}=...)."
            )
    for name in NOUN_REFUSED_PARAMS:
        if name not in kwargs:
            continue
        if name in ("values", "type"):
            raise TypeError(
                f"{noun}.map() draws {field}; it has no {name}= to override it. "
                f"To draw an arbitrary array on this grid use "
                f"`myflopy.plot.map(model, values=...)`."
            )
        raise TypeError(
            f"{noun}.map() does not accept {name}=: it always supplies its own "
            f"hover, so the legacy hover arguments have no effect here. Pass "
            f"`hover=` a HoverSpec instead."
        )


#: The line-profile options a field noun's ``section`` NAMES and forwards to
#: ``XSection``, each measured to change the picture
#: (``docs/reproducers/section_noun_kwargs_tail.py``, ledger 334) -- ``use_rbf``,
#: ``num_points`` and ``clip`` only when interpolating, and ``XSection`` refuses
#: them otherwise (ledger 337). Everything else the noun used to take through
#: ``**kwargs`` was either the cross-section's, inert on a noun
#: (``animation_kstpkpers`` -- a noun returns a figure, never ``.ani``; ``spacing``
#: and ``extrapolate_beyond_section_ends``, removed everywhere in ledger 337), or
#: ``section_name``, which crashed a group and a diff by colliding with the name
#: each member's line carries.
NOUN_SECTION_PARAMS = (
    "line", "cells", "per", "kstpkper", "layer", "interpolate", "interpolator",
    "use_rbf", "num_points", "clip", "x_or_y", "show_model_top", "show_model_btm",
)

#: The model's CROSS-SECTION's arguments. A field noun draws its field's line
#: profile, so these are refused BY NAME, pointing at the verb that draws the
#: cross-section -- the user's call (ledger 334), taken over routing them: a noun
#: fixes the field it draws (why map nouns refuse ``values=``), ``fill="layer"``
#: draws the geology and an array any field, and a group or a diff has no one
#: model's cells to fill.
NOUN_SECTION_REFUSED = (
    "fill", "fill_cmap", "fill_label", "layers", "head_layers", "layer_labels",
    "show_grid", "legend",
)

#: What a ``kind="section"`` composer forwards to each panel's sections: the
#: noun's options less the ones the composer owns (``per``, ``layer``, ``model``)
#: and ``kstpkper``, which each panel overwrites with its own time.
COMPOSER_SECTION_OPTIONS = tuple(
    name for name in NOUN_SECTION_PARAMS if name not in ("per", "kstpkper", "layer")
)


def refuse_section_arguments(where: str, unexpected: dict, accepted) -> None:
    """Raise a ``TypeError`` naming every keyword a noun's section did not take.

    ``where`` is how the caller spelled the call (``"model.hds.section()"``).
    The cross-section's arguments get the verb that draws them; ``surf_type``
    gets why a noun refuses it; anything else is named, with the nearest
    accepted spelling when there is one. Nothing is forwarded from here: a
    keyword that reaches this function was not taken.
    """

    if not unexpected:
        return
    names = sorted(unexpected)
    cross = [name for name in names if name in NOUN_SECTION_REFUSED]
    if cross:
        raise TypeError(
            f"{where} draws this field's LINE PROFILE, so it takes no "
            f"{', '.join(name + '=' for name in cross)} -- those draw the model's "
            f"CROSS-SECTION, the cells coloured by layer or by a field. That picture "
            f"is model.plot.section(line=..., fill=...); from a group or a diff, one "
            f"member's: group.models[name].plot.section(line=..., fill=...)."
        )
    if "surf_type" in unexpected:
        raise TypeError(
            f"{where} draws this field; surf_type= would section the layer "
            f"elevations instead. Section them from the model: "
            f"model.plot.section(line=..., fill=None, surf_type='lyr', interpolate=True)."
        )
    if "kstpkper" in unexpected and "per" not in accepted:
        raise TypeError(
            f"{where} takes no kstpkper=: each panel of a section composer reads "
            f"its own time. Choose the times with per=."
        )
    import difflib

    described = []
    for name in names:
        near = difflib.get_close_matches(name, list(accepted), n=1)
        described.append(f"{name!r}" + (f" (did you mean {near[0]!r}?)" if near else ""))
    raise TypeError(
        f"{where} got unexpected keyword argument(s) {', '.join(described)}. It "
        f"takes {', '.join(accepted)}."
    )


#: Aggregations whose 0 is a real answer for a group holding no value.
_COUNTING_AGGREGATIONS = ("count", "size", "nunique")


def _aggregate_leaving_gaps(by, aggregation):
    """``by.agg(aggregation)``, with a group that holds no value at all left ``NaN``.

    pandas sums an all-``NaN`` group to 0.0, so a lake with no stage (``NaN`` in
    the stage table, ledger 329; a reach's too since ledger 353) was drawn at 0 ft
    by any ``sum`` -- in a series and on a map, where it also stretched the colour
    range down to 0.
    ``mean``/``first``/``min``/``max`` already give ``NaN``; a counting
    aggregation's 0 is a real answer and is left alone. A tuple, not a set, so an
    unhashable aggregation (a list) still compares.

    ``by`` is a ``SeriesGroupBy``, grouped with or without ``as_index``.
    """

    result = by.agg(aggregation)
    if aggregation in _COUNTING_AGGREGATIONS:
        return result
    counts = by.count()
    if isinstance(counts, pd.DataFrame):
        counts = counts[by.obj.name]
    has_value = counts.to_numpy() > 0
    if isinstance(result, pd.DataFrame):
        result[by.obj.name] = result[by.obj.name].where(has_value)
        return result
    return result.where(has_value)


def build_cell_input_map_payload(
    frame: pd.DataFrame,
    *,
    ncpl: int,
    value_column: str,
    per: int | None,
    layer: int | None,
    multiplier: float = 1.0,
    fill_value: float = 0.0,
    agg: str = "sum",
) -> tuple[list[float], dict[str, list]]:
    """Convert a normalized cell-input table into choropleth values and hover.

    Parameters
    ----------
    frame
        Normalized table from :func:`build_cell_package_input_table` or a
        compatible filtered table.
    ncpl
        Number of cells per layer in the target grid.
    value_column
        Numeric column to render as choropleth values.
    per, layer
        Metadata added to the hover table. On a RESULTS table -- one carrying
        ``kstpkper``, its file's clock -- ``per`` also chooses the rows drawn: the
        period's LAST saved step (ledger 337's rule, ``_rows_drawn_for_period``),
        before ``agg`` combines anything. It combined every step the period saved
        until ledger 348: a flow's steps were added (``agg="sum"``) and a stage's
        FIRST step painted (``agg="first"``). An input table (no ``kstpkper``) or
        ``per=None`` is drawn as given.
    multiplier
        Optional scalar multiplier applied to the mapped values.
    fill_value
        Value used for cells not present in the input data for the selected
        period/layer. A cell that IS present but whose rows hold no value is
        ``NaN`` -- blank -- under every aggregation but a counting one, never
        ``fill_value`` and never 0.0 (ledger 329).
    agg
        Aggregation name passed to ``groupby().agg(...)`` for duplicate cells.
    """

    # Emptiness FIRST. A selection that matched no records comes back WITHOUT
    # its value column -- `wel.inputs.get(per=0, layer=0)` on a model whose
    # wells are in layers 1 and 3 returns
    # `['model', 'package', 'per', 'layer', 'cell']` -- so guarding the column
    # first turned "no records for this period/layer" into
    # `KeyError: Value column 'q' was not found`, naming a column the package
    # certainly has and blaming the wrong thing. It also made the empty branch
    # below, which exists precisely to draw an all-fill map, unreachable.
    full_index = pd.Index(range(int(ncpl)), name="cell")
    if frame.empty:
        values = pd.Series(float(fill_value), index=full_index, name=value_column)
        hover = {
            "Package": ["" for _ in full_index],
            "Period": [per for _ in full_index],
            "Layer": [layer for _ in full_index],
            "Cell": full_index.to_list(),
            "Record Count": [0 for _ in full_index],
            value_column: values.astype(float).tolist(),
        }
        return values.astype(float).tolist(), hover

    if value_column not in frame.columns:
        raise KeyError(
            f"Value column {value_column!r} was not found; the selection has "
            f"{sorted(frame.columns)}."
        )

    selected = _rows_drawn_for_period(frame, per).copy()
    # The multiplier scales each RECORD before ``agg`` combines them, and never the
    # fill: it acted after ``agg`` and on the ``fill_value`` cells, so a fill of 1.0
    # drew 2.0 at ``multiplier=2`` -- and 1.0 when the selection held no record at
    # all (ledger 473).
    scaled = pd.to_numeric(selected[value_column], errors="coerce") * float(multiplier)
    values = _aggregate_leaving_gaps(
        scaled.groupby(selected["cell"], dropna=False), agg
    ).reindex(full_index, fill_value=fill_value).astype(float)
    grouped = selected.groupby("cell", dropna=False)

    hover: dict[str, list] = {
        "Package": grouped["package"]
        .agg(_aggregate_hover_strings)
        .reindex(full_index, fill_value="")
        .tolist()
        if "package" in selected.columns
        else ["" for _ in full_index],
        "Period": [per for _ in full_index],
        "Layer": [layer for _ in full_index],
        "Cell": full_index.to_list(),
        "Record Count": grouped.size()
        .reindex(full_index, fill_value=0)
        .astype(int)
        .tolist(),
        value_column: values.tolist(),
    }

    excluded = {"model", "package", "per", "layer", "cell", value_column}
    for column in selected.columns:
        if column in excluded:
            continue
        column_group = grouped[column]
        if pd.api.types.is_numeric_dtype(selected[column]):
            hover[column] = (
                column_group.agg("first")
                .reindex(full_index, fill_value=np.nan)
                .tolist()
            )
        else:
            hover[column] = (
                column_group.agg(_aggregate_hover_strings)
                .reindex(full_index, fill_value="")
                .tolist()
            )
    return values.tolist(), hover


def build_group_input_compare_map_payload(
    frame: pd.DataFrame,
    *,
    ncpl: int,
    value_column: str,
    diff_column: str,
    model_name: str,
    reference_model: str,
    per: int | None,
    layer: int | None,
    multiplier: float = 1.0,
    fill_value: float = 0.0,
    agg: str = "sum",
) -> tuple[list[float], dict[str, list], float]:
    """Convert grouped comparison rows into diff-map values and hover data.

    Parameters
    ----------
    frame
        Filtered comparison table produced by a grouped ``compare()`` helper.
    ncpl
        Number of cells per layer in the reference grid.
    value_column
        Raw model value column, such as ``"recharge"`` or ``"finf"``.
    diff_column
        Difference column to map, such as ``"recharge_diff"``.
    model_name, reference_model
        Added to hover metadata.
    per, layer
        Added to hover metadata.
    multiplier
        Optional multiplier applied to the mapped difference values.
    fill_value
        Value used where a cell has no comparison row. A cell with rows but no
        difference in any of them is ``NaN``, not 0.0: an undefined difference
        is not "no difference" (ledger 329).
    agg
        Aggregation used when more than one row maps to the same cell.

    Returns
    -------
    tuple[list[float], dict[str, list], float]
        Choropleth z-values, hover payload, and the symmetric absolute color
        range suggested for diverging diff maps.
    """

    # Emptiness FIRST. A selection that matched no records comes back WITHOUT
    # its value column -- `wel.inputs.get(per=0, layer=0)` on a model whose
    # wells are in layers 1 and 3 returns
    # `['model', 'package', 'per', 'layer', 'cell']` -- so guarding the column
    # first turned "no records for this period/layer" into
    # `KeyError: Value column 'q' was not found`, naming a column the package
    # certainly has and blaming the wrong thing. It also made the empty branch
    # below, which exists precisely to draw an all-fill map, unreachable.
    full_index = pd.Index(range(int(ncpl)), name="cell")
    if frame.empty:
        values = pd.Series(float(fill_value), index=full_index, name=diff_column)
        hover = {
            "Model": [model_name for _ in full_index],
            "Reference Model": [reference_model for _ in full_index],
            "Period": [per for _ in full_index],
            "Layer": [layer for _ in full_index],
            "Cell": full_index.to_list(),
            "Record Count": [0 for _ in full_index],
            diff_column: values.astype(float).tolist(),
        }
        return values.astype(float).tolist(), hover, float(abs(fill_value))

    for column, label in ((diff_column, "Difference"), (value_column, "Value")):
        if column not in frame.columns:
            raise KeyError(
                f"{label} column {column!r} was not found; the selection has "
                f"{sorted(frame.columns)}."
            )

    grouped = frame.groupby("cell", dropna=False)
    # Each difference scaled BEFORE ``agg``, and the fill never (ledger 473).
    scaled = pd.to_numeric(frame[diff_column], errors="coerce") * float(multiplier)
    diff_values = _aggregate_leaving_gaps(
        scaled.groupby(frame["cell"], dropna=False), agg
    ).reindex(full_index, fill_value=fill_value).astype(float)
    current_values = (
        grouped[value_column].agg("first").reindex(full_index, fill_value=np.nan)
    )
    reference_column = f"reference_{value_column}"
    reference_values = (
        grouped[reference_column].agg("first").reindex(full_index, fill_value=np.nan)
        if reference_column in frame.columns
        else pd.Series(np.nan, index=full_index)
    )

    hover: dict[str, list] = {
        "Model": [model_name for _ in full_index],
        "Reference Model": [reference_model for _ in full_index],
        "Period": [per for _ in full_index],
        "Layer": [layer for _ in full_index],
        "Cell": full_index.to_list(),
        "Record Count": grouped.size()
        .reindex(full_index, fill_value=0)
        .astype(int)
        .tolist(),
        value_column: current_values.tolist(),
        reference_column: reference_values.tolist(),
        diff_column: diff_values.tolist(),
    }

    excluded = {
        "model",
        "reference_model",
        "package",
        "per",
        "layer",
        "cell",
        value_column,
        reference_column,
        diff_column,
    }
    for column in frame.columns:
        if column in excluded:
            continue
        column_group = grouped[column]
        if pd.api.types.is_numeric_dtype(frame[column]):
            hover[column] = (
                column_group.agg("first")
                .reindex(full_index, fill_value=np.nan)
                .tolist()
            )
        else:
            hover[column] = (
                column_group.agg(_aggregate_hover_strings)
                .reindex(full_index, fill_value="")
                .tolist()
            )

    absmax = (
        float(np.nanmax(np.abs(diff_values.to_numpy(dtype=float))))
        if len(diff_values)
        else 0.0
    )
    return diff_values.tolist(), hover, absmax


def _xy_panels_bounds(panels):
    """Global finite x/y bounds across ``[(label, [(name, x, y)])]`` panels."""

    xs_parts, ys_parts = [], []
    for _, lines in panels:
        for _, x_values, y_values in lines:
            x_array = np.asarray(x_values, dtype=float)
            y_array = np.asarray(y_values, dtype=float)
            xs_parts.append(x_array[np.isfinite(x_array)])
            ys_parts.append(y_array[np.isfinite(y_array)])
    xs_all = np.concatenate(xs_parts) if xs_parts else np.asarray([])
    ys_all = np.concatenate(ys_parts) if ys_parts else np.asarray([])
    x_bounds = (float(xs_all.min()), float(xs_all.max())) if xs_all.size else None
    y_bounds = (float(ys_all.min()), float(ys_all.max())) if ys_all.size else None
    return x_bounds, y_bounds


def _xy_mosaic_plotly(panels, *, ncols, title, xaxis_title, yaxis_title, markers):
    """Shared-scale grid of xy panels (series or section profiles)."""

    if not panels:
        raise ValueError("mosaic requires at least one panel.")
    labels = [label for label, _ in panels]
    ncols = min(int(ncols), len(panels)) or 1
    nrows = int(np.ceil(len(panels) / ncols))
    fig = subplots(
        nrows,
        ncols,
        shared_xaxes="all",
        shared_yaxes="all",
        subplot_titles=labels,
    )
    mode = "lines+markers" if markers else "lines"
    seen = set()
    for index, (_, lines) in enumerate(panels):
        row, col = index // ncols + 1, index % ncols + 1
        for name, x_values, y_values in lines:
            fig.add_scatter(
                x=x_values,
                y=y_values,
                mode=mode,
                name=str(name),
                legendgroup=str(name),
                showlegend=str(name) not in seen,
                row=row,
                col=col,
            )
            seen.add(str(name))
    fig.update_xaxes(title_text=xaxis_title, row=nrows)
    fig.update_yaxes(title_text=yaxis_title, col=1)
    fig.update_layout(title=title, uirevision="lock")
    return fig


def _xy_mosaic_mpl(panels, *, ncols, title, xlabel, ylabel, markers):
    """Matplotlib grid of xy panels with shared axes and one legend."""

    if not panels:
        raise ValueError("mosaic requires at least one panel.")
    labels = [label for label, _ in panels]
    ncols = min(int(ncols), len(panels)) or 1
    nrows = int(np.ceil(len(panels) / ncols))
    fig, axes = mpl_axes(
        nrows,
        ncols,
        figsize=(5.0 * ncols, 3.8 * nrows),
        squeeze=False,
        sharex=True,
        sharey=True,
    )
    marker = "o" if markers else None
    handles: dict[str, object] = {}
    for axis, label, (_, lines) in zip(axes.flat, labels, panels, strict=False):
        for name, x_values, y_values in lines:
            (handle,) = axis.plot(
                x_values, y_values, marker=marker, linewidth=1.8, label=str(name)
            )
            handles.setdefault(str(name), handle)
        axis.set_title(str(label))
    for axis in axes.flat[len(panels):]:
        axis.set_visible(False)
    for axis in axes[-1, :]:
        axis.set_xlabel(xlabel)
    for axis in axes[:, 0]:
        axis.set_ylabel(ylabel)
    if handles:
        fig.legend(handles.values(), handles.keys(), loc="upper right")
    fig.suptitle(title)
    return fig


def _xy_animation_plotly(frames, *, title, xaxis_title, yaxis_title, markers):
    """Plotly play/slider animation over ``[(label, [(name, x, y)])]`` frames."""

    if not frames:
        raise ValueError("animate requires at least one frame.")
    x_bounds, y_bounds = _xy_panels_bounds(frames)
    mode = "lines+markers" if markers else "lines"

    def _traces(lines):
        """Build one ``go.Scatter`` per ``(name, x, y)`` line for a single frame."""

        return [
            go.Scatter(x=x_values, y=y_values, mode=mode, name=str(name))
            for name, x_values, y_values in lines
        ]

    names = [str(label) for label, _ in frames]
    fig = Fig(data=_traces(frames[0][1]))
    fig.frames = [
        go.Frame(data=_traces(lines), name=str(label)) for label, lines in frames
    ]
    play = {"frame": {"duration": 600, "redraw": True}, "fromcurrent": True}
    pause = {"frame": {"duration": 0, "redraw": False}, "mode": "immediate"}
    fig.update_layout(
        title=title,
        xaxis_title=xaxis_title,
        yaxis_title=yaxis_title,
        uirevision="lock",
        updatemenus=[
            {
                "type": "buttons",
                "showactive": False,
                "buttons": [
                    {"label": "Play", "method": "animate", "args": [None, play]},
                    {"label": "Pause", "method": "animate", "args": [[None], pause]},
                ],
            }
        ],
        sliders=[
            {
                "active": 0,
                "steps": [
                    {
                        "method": "animate",
                        "args": [[name], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate"}],
                        "label": name,
                    }
                    for name in names
                ],
            }
        ],
    )
    # freeze axes so the animation does not rescale frame to frame
    if x_bounds is not None:
        span = (x_bounds[1] - x_bounds[0]) or 1.0
        fig.update_xaxes(range=[x_bounds[0] - 0.05 * span, x_bounds[1] + 0.05 * span])
    if y_bounds is not None:
        span = (y_bounds[1] - y_bounds[0]) or 1.0
        fig.update_yaxes(range=[y_bounds[0] - 0.05 * span, y_bounds[1] + 0.05 * span])
    return fig


def _xy_animation_mpl(frames, *, title, xlabel, ylabel, markers):
    """Matplotlib ``FuncAnimation`` over ``[(label, [(name, x, y)])]`` frames."""

    from matplotlib.animation import FuncAnimation

    if not frames:
        raise ValueError("animate requires at least one frame.")
    labels = [str(label) for label, _ in frames]
    x_bounds, y_bounds = _xy_panels_bounds(frames)
    marker = "o" if markers else None
    fig, axis = mpl_axes(1, 1, figsize=(8.0, 5.0))

    def _draw(index):
        """Redraw the axes for animation frame ``index`` with fixed shared limits."""

        axis.clear()
        for name, x_values, y_values in frames[index][1]:
            axis.plot(x_values, y_values, marker=marker, linewidth=1.8, label=str(name))
        axis.set_title(f"{title} - {labels[index]}")
        axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)
        if frames[index][1]:
            axis.legend(loc="best")
        if x_bounds is not None:
            span = (x_bounds[1] - x_bounds[0]) or 1.0
            axis.set_xlim(x_bounds[0] - 0.05 * span, x_bounds[1] + 0.05 * span)
        if y_bounds is not None:
            span = (y_bounds[1] - y_bounds[0]) or 1.0
            axis.set_ylim(y_bounds[0] - 0.05 * span, y_bounds[1] + 0.05 * span)

    _draw(0)
    return FuncAnimation(fig, _draw, frames=len(frames), interval=600, blit=False)


class SpatialView:
    """Unified view grammar for cell-mapped package fields.

    Panel verbs :meth:`map` (plan view) and :meth:`plot` (series view) plus
    composers :meth:`mosaic` and :meth:`animate`, exposed identically on
    single-model, model-group, and diff surfaces (inputs and results). A host
    supplies the atomic styled panel via :meth:`_spatial_map` (and optionally
    the series hooks) plus the dimensions it has; this mixin composes single
    panels, faceted mosaics (by ``layer`` or ``model``), and animations (over
    ``period`` or ``model``) for both backends.

    ``backend="plotly"`` (default) returns an interactive figure/``Choro``;
    ``backend="mpl"`` returns a matplotlib ``Figure`` (or, from :meth:`animate`,
    a ``FuncAnimation``). The mapped *content* -- raw value vs difference -- is
    decided by the node: a ``diff.*`` node maps deltas, an inputs/results node
    maps raw values.
    """

    _PLOTLY_BACKENDS = PLOTLY_BACKENDS
    _MPL_BACKENDS = MPL_BACKENDS

    # -- hooks the host provides ------------------------------------------
    def _spatial_map(self, *, per, layer, model=None, **kwargs):
        """Return one styled Plotly ``Choro`` for a (per, layer, model).

        Default: delegate to the host's atomic ``map(per=, layer=, ...)`` builder
        -- the single-model surfaces, whose ``map`` already returns a ``Choro``,
        so they gain :meth:`mosaic`/:meth:`animate` for free. Group/diff surfaces
        -- which need a ``model`` selector plus their own delta/colorscale policy
        -- override this hook and let :meth:`map` come from this mixin.
        """

        return self.map(per=int(per), layer=int(layer), **kwargs)

    def _spatial_reference_model(self):
        """Model whose grid/layers frame the view (single-model default)."""

        return self.model

    def _spatial_layers(self) -> list[int]:
        """Layers to offer for faceting: those present in the data, else every grid layer."""

        # Prefer the layers actually present in this field's data (many BC
        # packages live in a single layer); fall back to every grid layer.
        frame = self.get()
        columns = getattr(frame, "columns", [])
        if "layer" in columns and not frame.empty:
            present = sorted({int(value) for value in frame["layer"].dropna().tolist()})
            if present:
                return present
        return list(range(views_nlay(self._spatial_reference_model().gwf)))

    def _spatial_periods(self) -> list[int]:
        """Stress periods present in the data (``[0]`` when there is no period axis)."""

        frame = self.get()
        columns = getattr(frame, "columns", [])
        if "per" not in columns or frame.empty:
            return [0]
        return sorted({int(value) for value in frame["per"].dropna().tolist()})

    def _spatial_models(self) -> list[str] | None:
        """Model names for the model axis; ``None`` for a single-model surface."""

        return None  # single-model surfaces have no model axis

    def _spatial_default_facet(self) -> str:
        """Default mosaic facet: ``"model"`` for a group, else ``"layer"``."""

        return "model" if self._spatial_models() else "layer"

    def _spatial_value_label(self) -> str:
        """Human label for the mapped quantity (field/value/package name, else "value")."""

        for attribute in ("field_name", "value_name", "package_name"):
            value = getattr(self, attribute, None)
            if value is not None:
                return str(value)
        return "value"

    def _spatial_is_diff(self) -> bool:
        """Whether panels show a signed difference (drives diverging, zero-centered scales)."""

        return False

    # -- series hooks (the plot verb) ---------------------------------------
    #: entity id columns that get one line each in ``plot`` when present
    _SERIES_ENTITY_COLUMNS = ("lake", "reach", "well")

    def _series_table(self) -> pd.DataFrame:
        """Long table backing :meth:`plot`; defaults to the node's ``get()``."""

        get = getattr(self, "get", None)
        if get is None:
            raise NotImplementedError(
                "plot() requires a tabular get(); this node does not provide one."
            )
        return get()

    def _series_value_column(self, frame: pd.DataFrame) -> str:
        """Column of :meth:`_series_table` that :meth:`plot` draws."""

        value = getattr(self, "value_name", None)
        if value is not None and value in frame.columns:
            return str(value)
        label = self._spatial_value_label()
        if label in frame.columns:
            return label
        raise KeyError(
            "Could not infer the value column for plot(); available columns: "
            f"{list(frame.columns)}."
        )

    def _series_default_agg(self) -> str:
        """How :meth:`plot` collapses cells within a line (sum for fluxes)."""

        return "sum"

    # -- dimension resolution ---------------------------------------------
    def _resolve_layers(self, layer) -> list[int]:
        """Normalize the ``layer`` selector to a list: ``None`` -> all present, int -> one, iterable -> many."""

        if layer is None:
            return self._spatial_layers()
        if isinstance(layer, (int, np.integer)):
            return [int(layer)]
        resolved = [int(value) for value in layer]
        if not resolved:
            raise ValueError("layer must select at least one layer.")
        return resolved

    def _single_layer(self, layer) -> int:
        """Reduce the ``layer`` selector to one layer (0 by default; first of a set)."""

        if layer is None:
            return 0
        if isinstance(layer, (int, np.integer)):
            return int(layer)
        return self._resolve_layers(layer)[0]

    def _resolve_periods(self, per) -> list[int]:
        """Normalize the ``per`` selector to a list: ``None`` -> all present, int -> one, iterable -> many."""

        if per is None:
            return self._spatial_periods()
        if isinstance(per, (int, np.integer)):
            return [int(per)]
        resolved = [int(value) for value in per]
        if not resolved:
            raise ValueError("per must select at least one period.")
        return resolved

    def _resolve_models(self, model) -> list:
        """Normalize the ``model`` selector against this surface's model axis.

        Single-model surfaces accept only ``None``/``"*"`` (returning ``[None]``);
        group surfaces map ``None``/``"*"`` to all models and validate any explicit
        names against the group.
        """

        models = self._spatial_models()
        if not models:
            if model not in (None, "*"):
                raise ValueError(
                    "This surface has a single model; 'model' is not applicable."
                )
            return [None]
        if model in (None, "*"):
            return list(models)
        names = [str(model)] if isinstance(model, str) else [str(m) for m in model]
        missing = [name for name in names if name not in models]
        if missing:
            raise KeyError(f"Models not in this group: {missing!r}")
        return names

    @classmethod
    def _normalize_backend(cls, backend: str) -> str:
        """Canonicalize a backend alias to ``"plotly"`` or ``"mpl"`` (raises if unknown)."""

        return normalize_backend(backend)

    def map(
        self,
        model=None,
        *,
        per: int = 0,
        layer: int = 0,
        backend: str = "plotly",
        **map_kwargs,
    ):
        """One choropleth panel.

        The ``map`` of a GROUP or DIFF surface -- ``group.hds``,
        ``group.packages.drn.results.q``, ``diff.packages.drn.results.q``, ...
        A single model's nouns define their own ``map``, with its own parameters.
        ``backend="plotly"`` (default) returns an interactive ``Choro``;
        ``backend="mpl"`` returns a static matplotlib ``Figure``.

        Parameters
        ----------
        model : str, optional
            Which model to draw, by name; the one positional. On a group,
            ``None`` (the default) draws the group's reference model, and a name
            not in the group raises ``KeyError``. On a diff, the COMPARED model
            whose difference from the reference is drawn: ``None`` picks it when
            the group has exactly one non-reference model and raises
            ``ValueError`` otherwise, as does the reference's own name.
        per : int, default 0
            Stress period, zero-based. Passed through ``int()`` (so ``None``
            raises ``TypeError``) and handed to the panel builder with ``layer``.
        layer : int, default 0
            Layer, zero-based. Passed through ``int()`` the same way.
        backend : {'plotly', 'mpl'}, default 'plotly'
            ``'plotly'`` (or ``'interactive'``) returns the ``Choro`` picture;
            ``'mpl'`` (or ``'matplotlib'``, ``'static'``) draws that same picture
            with ``Choro.plot_mpl()`` and returns the bare ``Figure``. Anything
            else raises ``ValueError``.
        **map_kwargs
            Forwarded to this surface's panel builder. A group FIELD
            (``group.hds``/``.conc``/``.temp``) hands them to the member's own
            ``map(...)``, and a group STAGE view (``lak``/``sfr``
            ``results.stage``) to the member's stage noun's ``map`` (its ``agg``
            defaults to ``'first'``). Every other group view first takes its own
            -- ``multiplier`` (default 1.0), ``fill_value`` (0.0) and
            ``colorscale``; ``agg`` (``'sum'``; on an inputs view the field's
            kind decides, ledger 472), which the LAK and SFR exchange
            accept and do not use and the surface-water exchange does not take;
            ``value_column`` on an inputs view; ``connection_type`` on the LAK
            exchange; ``include``/``lak_connection_type`` on the surface-water
            exchange -- and hands the rest to the model's bound map verb -- the
            drawing options of :func:`myflopy.plot.map` (``zmin``/``zmax``,
            ``contours``, ``features``, ``extent``, ...). A diff hands them to the
            group's ``compare_map``, which ends at the reference model's bound map
            verb.

        Returns
        -------
        Choro or matplotlib.figure.Figure
        """

        choro = self._spatial_map(per=int(per), layer=int(layer), model=model, **map_kwargs)
        if self._normalize_backend(backend) == "plotly":
            return choro
        return choro.plot_mpl()

    @staticmethod
    def _series_line_label(keys, key) -> str:
        """Human label for one plotted line (model / entity / layer / cell)."""

        parts = []
        for column, value in zip(keys, key, strict=False):
            if pd.isna(value):
                continue
            if column == "model":
                parts.append(str(value))
            elif column == "layer":
                parts.append(f"L{int(value)}")
            elif column == "cell":
                parts.append(f"C{int(value)}")
            else:  # entity columns: lake, reach, well
                parts.append(f"{column.capitalize()} {int(value)}")
        return " / ".join(parts) if parts else "All"

    def _no_series_message(self) -> str:
        """Why ``plot()`` has nothing to draw here, and what to call instead.

        ``plot()`` is the node's NON-SPATIAL CHART, and for a per-cell field
        that chart is a series by stress period (see
        ``docs/view_layer_conventions.md``; the verb is not "time series" in
        general -- on other nodes it is a profile, a bar chart or a histogram).
        A static property (NPF ``k``, DISV ``top`` ...) has no stress-period
        dimension at all, so there is no series to draw -- and saying only
        "this node has no 'per' data" leaves the reader to guess that the
        spatial verb is ``map()``.
        """

        package = getattr(self, "package_name", None)
        field = getattr(self, "field_name", None)
        noun = f"{package}.{field}" if package and field else "this field"

        spatial = []
        for verb, call in (("map", "map(layer=0)"), ("section", "section(...)")):
            if callable(getattr(self, verb, None)):
                spatial.append(f".{call}")
        alternatives = " or ".join(spatial) if spatial else ".get()"

        return (
            f"plot() charts {noun} by stress period, but it is a static property "
            f"with no stress-period dimension -- there is nothing to put on the "
            f"time axis. Use {alternatives} for its spatial distribution, or "
            ".get() / .summary() for the values themselves."
        )

    def _series_selection(self, *, model=None, per=None, layer=None, cells=None):
        """Filter the series table; return ``(frame, value_column, line_keys)``."""

        selected_models = self._resolve_models(model)
        frame = self._series_table()
        if frame is None or "per" not in getattr(frame, "columns", []):
            raise ValueError(self._no_series_message())
        frame = frame.copy()

        # column-based filters (a filter is skipped when its column is absent)
        if model is not None and "model" in frame.columns:
            frame = frame[frame["model"].isin([str(name) for name in selected_models])]
        if per is not None:
            per_values = (
                [int(per)] if isinstance(per, (int, np.integer)) else [int(v) for v in per]
            )
            frame = frame[frame["per"].isin(per_values)]
        if layer is not None and "layer" in frame.columns:
            layer_values = (
                [int(layer)]
                if isinstance(layer, (int, np.integer))
                else [int(v) for v in layer]
            )
            frame = frame[frame["layer"].isin(layer_values)]
        if cells is not None and "cell" in frame.columns:
            cell_values = (
                [int(cells)]
                if isinstance(cells, (int, np.integer))
                else [int(v) for v in cells]
            )
            frame = frame[frame["cell"].isin(cell_values)]

        # One point per period: its LAST saved step, on the table's own clock and,
        # in a group table, on each MODEL's (ledger 337's rule). A results table
        # carries every step a period saved, and the point was the SUM of a flow's
        # steps (3.0x at nstp=3 under SAVE BUDGET ALL) or the MEAN of a stage's
        # until ledger 348. A table with no clock is one record per period already.
        if "kstpkper" in frame.columns:
            frame = at_last_saved_steps(
                frame, by="model" if "model" in frame.columns else None
            )

        value_column = self._series_value_column(frame)
        frame[value_column] = pd.to_numeric(frame[value_column], errors="coerce")

        # line keys: model axis, entity ids, multi-layer split, explicit cells.
        # A single-model surface may still carry a ``model`` column in its
        # table; only treat it as a line axis when the surface has one.
        keys = ["model"] if self._spatial_models() and "model" in frame.columns else []
        keys += [c for c in self._SERIES_ENTITY_COLUMNS if c in frame.columns]
        if "layer" in frame.columns and frame["layer"].nunique() > 1:
            keys.append("layer")
        if cells is not None and "cell" in frame.columns:
            keys.append("cell")
        return frame, value_column, keys

    def _series_lines(self, frame, keys, value_column, aggregation):
        """Aggregate a filtered series frame into ``[(label, x, y)]`` lines."""

        # A point with no value at all stays a gap: pandas sums an all-NaN group
        # to 0.0, so a dry lake's masked stage (ledger 329) was drawn at 0 ft.
        by = frame.groupby([*keys, "per"], dropna=False, as_index=False)[value_column]
        grouped = _aggregate_leaving_gaps(by, aggregation)
        lines = []
        if keys:
            for key, sub in grouped.groupby(keys, dropna=False):
                key = key if isinstance(key, tuple) else (key,)
                sub = sub.sort_values("per")
                lines.append(
                    (
                        self._series_line_label(keys, key),
                        sub["per"].astype(int).to_numpy(),
                        sub[value_column].astype(float).to_numpy(),
                    )
                )
        elif not grouped.empty:
            sub = grouped.sort_values("per")
            lines.append(
                (
                    self._spatial_value_label(),
                    sub["per"].astype(int).to_numpy(),
                    sub[value_column].astype(float).to_numpy(),
                )
            )
        return lines

    def plot(
        self,
        model=None,
        *,
        per=None,
        layer=None,
        cells=None,
        agg: str | None = None,
        backend: str = "plotly",
        title: str | None = None,
    ):
        """Series panel: this field by stress period.

        One line per model (groups/diffs), per entity (lake/reach) when the
        package has one, and per cell when ``cells`` selects specific cells;
        otherwise cells are aggregated (``agg``, default sum for fluxes).
        ``backend="plotly"`` (default) returns an interactive ``viz.Fig``;
        ``backend="mpl"`` a matplotlib ``Figure``. Stress periods are zero-based.
        A period's point is its LAST saved step (each model's own, in a group):
        ``agg`` combines the cells of that one step, never the steps of a period
        (ledgers 337, 348).

        A static property (NPF ``k``, DISV ``top``) has no stress periods, and
        raises ``ValueError`` naming ``map()`` instead.

        Parameters
        ----------
        model : str or sequence of str, optional
            On a group or a diff, the model(s) to draw, one line each -- on a diff,
            the COMPARED models: a name, a list of names, or ``'*'``/``None`` (the
            default) for every one. A name the surface does not draw raises
            ``KeyError`` (on a diff, the reference's own name included). On a
            single model's noun only ``None``/``'*'`` is accepted; anything else
            raises ``ValueError``.
        per : int or sequence of int, optional
            Keep only these stress periods (zero-based); ``None`` (the default)
            keeps every period in the table.
        layer : int or sequence of int, optional
            Keep only these layers (zero-based); ``None`` (the default) keeps all,
            and when more than one remains each layer is its own line. IGNORED,
            with nothing said, on a table that has no ``layer`` column (a lake's
            stage, for one): every row is kept.
        cells : int or sequence of int, optional
            Keep only these cells (zero-based, within a layer) and draw one line
            per cell rather than combining them. IGNORED, with nothing said, on a
            table that has no ``cell`` column.
        agg : str, optional
            How the cells in one line combine at a period: any pandas groupby
            reduction NAME -- ``'sum'``, ``'mean'``, ``'min'``, ``'max'``,
            ``'median'``, ``'first'``, ``'count'``, ... A callable is not accepted,
            nor is an unknown name: both raise ``AttributeError``. ``None`` (the
            default) is the surface's own: ``'mean'`` for heads, concentration and
            temperature (one model or a group), their difference, and a single
            model's LAK/SFR stage; a package INPUT field's own kind on one model,
            a group and a diff, a LAK/SFR input field's included
            (``package_registry.input_field_agg``: ``'sum'`` for a conductance, a
            rate or a connection area, ``'mean'`` for an elevation or a property
            -- ``model.drn.inputs.elev.plot()`` drew the SUM of the drains'
            elevations until ledger 472, and ``model.lak.inputs.belev.plot()``
            the sum over a lake's connections until ledger 475); ``'sum'`` everywhere else
            -- a flux and a difference of fluxes, and a group's stage view (one
            row per lake or reach and period, so the sum is the stage). A point
            whose rows hold no value at all stays a gap rather than summing to 0.
        backend : {'plotly', 'mpl'}, default 'plotly'
            ``'plotly'`` (or ``'interactive'``) returns a :class:`~myflopy.viz.Fig`;
            ``'mpl'`` (or ``'matplotlib'``, ``'static'``) an 8 x 4 inch Matplotlib
            ``Figure``. Anything else raises ``ValueError``.
        title : str, optional
            The figure title. Default ``'<field> by stress period'``, with
            ``' (model - reference)'`` appended on a diff.

        Returns
        -------
        myflopy.viz.Fig or matplotlib.figure.Figure
        """

        backend_kind = self._normalize_backend(backend)
        frame, value_column, keys = self._series_selection(
            model=model, per=per, layer=layer, cells=cells
        )
        aggregation = self._series_default_agg() if agg is None else str(agg)
        lines = self._series_lines(frame, keys, value_column, aggregation)
        heading = title or (
            f"{self._spatial_value_label()} by stress period"
            + (" (model - reference)" if self._spatial_is_diff() else "")
        )
        if backend_kind == "plotly":
            fig = Fig()
            for label, xs_values, ys_values in lines:
                fig.add_scatter(
                    x=xs_values, y=ys_values, mode="lines+markers", name=label
                )
            fig.update_layout(
                title=heading, xaxis_title="Stress Period", yaxis_title=value_column
            )
            return fig
        fig, ax = mpl_axes(figsize=(8, 4))
        for label, xs_values, ys_values in lines:
            ax.plot(xs_values, ys_values, marker="o", linewidth=2.0, label=label)
        ax.set_title(heading)
        ax.set_xlabel("Stress Period")
        ax.set_ylabel(value_column)
        if lines:
            ax.legend()
        fig.tight_layout()
        return fig

    def _noun_path(self) -> str:
        """How a caller reached this noun, for messages (``model.hds``, ``group.conc``)."""

        return type(self).__name__

    # -- the composers' own selectors, by kind (ledger 461) ------------
    def _panel_map_names_agg(self) -> bool | None:
        """Whether this surface's panel builder takes ``agg``: ``None`` when it cannot tell.

        A single model's noun builds each panel with its own ``map``, which either
        NAMES ``agg`` (a record noun) or hands its tail to Plotly as trace
        properties, which have no ``agg``. A group or a diff overrides
        :meth:`_spatial_map`; when that names ``agg`` it takes it, and when it only
        forwards a tail the answer lies further down, so it is forwarded.
        """

        if type(self)._spatial_map is SpatialView._spatial_map:
            return "agg" in inspect.signature(self.map).parameters
        parameters = inspect.signature(self._spatial_map).parameters
        if "agg" in parameters:
            return True
        if any(p.kind is inspect.Parameter.VAR_KEYWORD for p in parameters.values()):
            return None
        return False

    def _map_composer_options(self, where: str, cells, agg) -> dict:
        """``cells``/``agg`` as a ``kind="map"`` composer forwards them to each panel.

        ``agg`` reaches every panel whose ``map`` takes one -- a record noun's
        combines a cell's records by it -- and is refused where the panel has
        none. ``cells`` is refused: no map takes it (the highlight is ``select=``).
        Both were dropped in silence until ledger 461, so
        ``q.mosaic(agg="count")`` drew the sum.
        """

        if cells is not None:
            raise TypeError(
                f"{where} draws every cell of each panel, so it takes no cells=: "
                "cells= picks a series panel's lines (kind='plot') or a section's "
                "path (kind='section'). To highlight cells on a map, pass select=."
            )
        if agg is None:
            return {}
        if self._panel_map_names_agg() is False:
            raise TypeError(
                f"{where} got agg={agg!r}, but this noun's map has no agg: each "
                "panel draws one value per cell. agg= reduces a series panel "
                "(kind='plot') here."
            )
        return {"agg": agg}

    def _section_composer_options(self, where: str, cells, agg, kwargs: dict) -> dict:
        """The line-profile options a ``kind="section"`` composer hands each panel.

        ``cells`` is the section's path, as it is on :meth:`section` -- the
        composer's own named ``cells`` used to capture it, so the panel got no
        line (ledger 461). ``agg`` has nothing to reduce on a profile and is
        refused; any other keyword is checked by :func:`refuse_section_arguments`.
        """

        if agg is not None:
            raise TypeError(
                f"{where} draws a line profile, which has no agg= to reduce by; "
                "agg= reduces a series panel (kind='plot') or a record map "
                "(kind='map')."
            )
        refuse_section_arguments(
            where,
            {k: v for k, v in kwargs.items() if k not in COMPOSER_SECTION_OPTIONS},
            COMPOSER_SECTION_OPTIONS,
        )
        return kwargs if cells is None else {**kwargs, "cells": cells}

    def mosaic(
        self,
        *,
        # -- what is faceted, and how ------------------------------------------
        kind: str = "map",
        by: str | None = None,
        per=None,
        layer=None,
        model=None,
        cells=None,
        agg: str | None = None,
        ncols: int = 3,
        title: str | None = None,
        sync_views: bool | None = None,
        # -- colour of each panel ---------------------------------------------
        zmin: float | None = None,
        zmax: float | None = None,
        colorscale: str | list | tuple | None = None,
        logscale: bool = False,
        # -- contours ----------------------------------------------------------
        contours: bool | str = False,
        contour_values=None,
        contour_levels: int | float | list = 10,
        contour_color: str = "black",
        contour_width: float = 1.5,
        contour_name: str | None = None,
        contour_clip: bool = True,
        contour_resolution: int = 150,
        contour_method: str = "linear",
        # -- highlighting -------------------------------------------------------
        select=None,
        select_style: str = "outline",
        select_color: str | None = None,
        # -- overlays and framing -----------------------------------------------
        features=None,
        hillshade_path=None,
        fit_bounds: bool = True,
        bounds_padding: float = 0.05,
        extent=None,
        # -- hover ---------------------------------------------------------------
        hover=None,
        show_layer_elevs: bool | None = None,
        hover_surfaces: bool | None = None,
        # -- renderer -------------------------------------------------------------
        backend: str = "plotly",
        **kwargs,
    ):
        """Faceted grid of panels sharing one scale.

        ``kind`` picks the panel type: ``"map"`` (choropleths, default),
        ``"plot"`` (series panels), or ``"section"`` (cross sections, heads leaves).
        ``by`` picks the facet axis -- maps: ``"layer"``/``"model"``; series:
        ``"model"``/``"layer"``/an entity (``"lake"``, ``"reach"``); sections:
        ``"model"``/``"period"``. Defaults to the surface's natural axis.

        **Every drawing option of ``map`` is named here and reaches each panel**
        (8.8, 2026-09-14): a mosaic of K over layers is
        ``mosaic(logscale=True, features=wells, show_layer_elevs=True)``. Only the
        options a caller actually SET are forwarded, so a panel keeps its noun's
        own colour policy otherwise. They are refused, not ignored, on
        ``kind="plot"``/``"section"``, which draw no choropleth.

        Parameters
        ----------
        kind : {'map', 'plot', 'section'}, default 'map'
            What each panel is: a choropleth, a series plot, or a line profile.
            ``'section'`` exists only on a dependent-variable noun (``model.hds``,
            ``.conc``, ``.temp``, ``group.hds``, ``diff.hds``) and raises
            ``ValueError`` on any other.
        by : str, optional
            Facet axis. Maps: ``'layer'`` or ``'model'``; series: ``'model'``,
            ``'layer'`` or an entity (``'lake'``, ``'reach'``); sections:
            ``'model'`` or ``'period'``. Defaults to this surface's natural axis.
        per : int, optional
            Stress period held constant across the panels (0-based).
        layer : int or sequence of int, optional
            The layer faceted over (``by='layer'``), or held constant otherwise.
            Defaults to every layer this field has data for.
        model : str or sequence of str, optional
            The model faceted over (``by='model'``), or held constant otherwise.
            Groups and diffs only.
        cells : int or sequence of int, optional
            ``kind='plot'``: restrict the series panels to these cells (zero-based,
            within a layer), one line per cell. ``kind='section'``: each panel's
            path as cell numbers, as :meth:`section` takes it (give it or
            ``line``). ``kind='map'`` takes none -- a map draws every cell, and
            the highlight is ``select=`` -- and raises ``TypeError``.
        agg : str, optional
            ``kind='plot'``: how several rows per period are reduced for a series
            panel (``'sum'``/``'mean'``/...). ``kind='map'``: handed to every
            panel's map, which combines a cell's records by it (a record noun's
            ``agg``) -- on a noun whose map has no ``agg`` (heads, a static
            array) it raises ``TypeError``. ``kind='section'`` has nothing to
            reduce and raises ``TypeError``.
        ncols : int, default 3
            Grid width; rows grow as needed.
        title : str, optional
            Overall figure title. Defaults to a label naming the field and axis.
        sync_views : bool, optional
            Pan and zoom every map panel together. By default (``None``) panels on
            ONE window -- the grid, or one shared ``extent`` -- are linked and
            panels given a LIST of extents are not (the first pan would snap them
            all to one window); ``True`` then raises. Ignored for non-map kinds.
        backend : {'plotly', 'mpl'}, default 'plotly'
            Which renderer draws the GRID of panels -- not each panel, which is
            always built as a ``Choro`` first. ``'plotly'`` returns one house
            ``viz.Fig`` of linked subplots; ``'mpl'`` a single static
            :class:`matplotlib.figure.Figure` whose panels share one colorbar
            (a log-scaled mosaic gets a log axis in real units). It earns its
            keep here: a composer has no intermediate picture to call
            ``.plot_mpl()`` on.
        **kwargs
            Options for each panel's builder, by ``kind``. ``'map'``: forwarded,
            with the drawing options named here, to every panel's map -- a single
            model's noun ``map`` (``multiplier``/``fill_value`` on a record noun,
            ``per_timestep``/``bgs``/``hover_layers`` on heads, ..., and through
            its own tail Plotly trace properties), or on a group or a diff what
            :meth:`map` forwards as ``**map_kwargs``. ``per``, ``layer`` and
            ``model`` are this composer's and are set per panel; ``agg`` and
            ``cells`` are named above and never get here.
            ``'section'``: the line profile's own options -- ``line``,
            ``interpolate``, ``interpolator``, ``use_rbf``, ``num_points``,
            ``clip``, ``x_or_y``, ``show_model_top``, ``show_model_btm`` (see
            :meth:`section`) -- and any other keyword raises ``TypeError``: a
            cross-section argument (``fill``, ``layers``, ...) naming
            ``model.plot.section(fill=...)``, ``kstpkper`` naming ``per=``, and
            anything else listing the names accepted. ``'plot'``: none -- the
            series panels read only the parameters named above, so any keyword
            raises ``TypeError`` listing them.

            .. generated from myflopy.plot.map -- edit there, then run scripts/derive_docstrings.py
        zmin, zmax : float, optional
            Fixed color-scale limits. Set both to hold the scale steady across
            frames or panels; ``zmin >= zmax`` raises rather than rendering one flat
            color. With neither, the range comes from the data.
        colorscale : str or list of (float, str), optional
            Any Plotly colorscale name -- ``'earth'`` (the default), ``'Viridis'``,
            ``'Turbo'``, ``'RdBu'``, ...; append ``_r`` to reverse one -- or explicit
            ``[[position, color], ...]`` stops. The static backend builds its
            colormap from the same stops, so a name renders the same way on both.
            An unknown name warns and falls back to ``'earth'``.
        logscale : bool, default False
            Color on a log10 scale, for a field that spans orders of magnitude
            (hydraulic conductivity, recharge). The colorbar is labelled in real
            units at whole decades on both backends and the hover keeps the real
            values. ``zmin``/``zmax`` (and a ``zmid``) are then in DATA units --
            ``zmin=1e-4, zmax=10`` -- and a non-positive limit raises. Non-positive
            values draw as gaps.
        contours : bool or str, default False
            Overlay contour lines. ``True`` contours the mapped values; a string
            names a different field to contour instead.
        contour_values : sequence of float, optional
            Contour a supplied array rather than the mapped one.
        contour_levels : int or float or list of float, default 10
            A count of levels, a fixed interval, or explicit level values.
        contour_color : str, default 'black'
            Line colour for the contour trace.
        contour_width : float, default 1.5
            Line width for the contour trace.
        contour_name : str, optional
            Legend name for the contour trace.
        contour_clip : bool, default True
            Clip contours to the active domain instead of the full grid extent.
        contour_resolution : int, default 150
            Grid size used to build the contours, per axis. Higher is smoother and
            slower. **Only acts under ``contour_method="cubic"``**, which is the one
            that interpolates onto a square grid before contouring; the default
            linear method triangulates the cell centres directly, so there is no
            grid for this to size and passing it changes nothing (measured: 338
            contour points at 40, 150 and 400 alike).
        contour_method : {'linear', 'cubic'}, default 'linear'
            How the scattered cell values become contours. ``'linear'``
            triangulates the cell centres and contours the triangulation --
            fast, exact at the centres, and faceted. ``'cubic'`` interpolates onto
            a ``contour_resolution``-square grid first (Clough-Tocher) and contours
            that -- smoother, slower, and able to overshoot between cells.
            ``'tri'``/``'tricontour'`` and ``'clough'``/``'clough_tocher'``/
            ``'cloughtocher'`` are accepted as aliases. Anything else raises naming
            both; ``'nearest'`` in particular was documented here for a while and
            has never been implemented.
        select : sequence of int, ndarray, str, Path or geometry, optional
            Cells to highlight. Cell indices, a boolean mask of length ``ncpl``, the
            name of a registered model region (``"all_streams"``), a path to a
            vector file, or a shapely/GeoPandas geometry, by coverage. Highlights
            nothing (with a warning) when the selection is empty.
        select_style : {'outline', 'dim', 'both'}, default 'outline'
            How the highlight is drawn. ``'outline'`` traces the dissolved boundary
            of the selection and leaves every cell at full opacity. ``'dim'`` fades
            the *unselected* cells to 20% instead, which suits a bare grid where the
            selection is the subject, but on a field it costs a measured 6.9x of
            readable contrast everywhere you did not select -- and one box/lasso
            gesture in the browser overwrites it. ``'both'`` draws each.
        select_color : str, optional
            Highlight colour, defaulting to ``viz.PALETTE.highlight``. Worth setting
            when the default red collides with a red-blue diverging colorscale.
        features : path, GeoDataFrame, geometry, Feature, or list/dict of those, optional
            Custom geometry drawn over the cells -- wells, pond footprints, a
            boundary. A path to any vector file GeoPandas reads (name a GeoPackage
            layer with ``Feature(path, layer=...)``), a ``GeoDataFrame`` or
            ``GeoSeries``, a shapely geometry, or a :class:`myflopy.plot.Feature`
            carrying its own ``color``/``width``/``size``/``fill``/``label``/
            ``hover``. A list draws several layers; a dict names them. Points
            become markers, lines polylines, polygons outlines (or fills). Each
            layer gets a legend entry, a stable house color and a hover listing its
            attribute columns. A source without a CRS is taken to be in the grid's.
            Drawn by both backends and kept by ``mosaic`` and ``animate``.
        hillshade_path : Path, optional
            A hillshade GeoTIFF to draw beneath the cells for topographic context.
        fit_bounds : bool, default True
            Fit the initial view to the grid extent rather than using ``zoom``. The
            fit is for the figure's own size, or plotly's default 700x450 when it
            names none -- a larger display only adds margin around the grid -- and
            ``.save()`` to an image refits for the export's size.
        bounds_padding : float, default 0.05
            Fractional padding around the fitted bounds.
        extent : tuple, points, geometry, GeoDataFrame, area or str, optional
            The window the map OPENS on, on both backends -- the one way to frame a
            static (matplotlib) map, which cannot be zoomed. Four numbers
            ``(xmin, ymin, xmax, ymax)`` in the grid's CRS (shapely's order, NOT
            FloPy's ``(xmin, xmax, ymin, ymax)``); a run of ``(x, y)`` points (their
            bounding box); any shapely geometry -- a polygon, a transect line -- by its
            bounds; a GeoDataFrame, reprojected from its own CRS; or the name of an
            area of interest (``model.aoi``). It sets the VIEW only: no cell is
            cropped and the colour scale still spans the whole field, so pass
            ``zmin``/``zmax`` to tighten it for a close-up. ``bounds_padding`` pads it
            as it pads the grid. A box off the grid raises, printing both boxes and
            both CRSs. In a ``mosaic``, one extent frames every panel and a LIST of
            them frames each panel on its own.
        hover : HoverSpec, optional
            Replace the sectioned hover outright. See
            :mod:`myflopy.modflow.utils.datatypes.hover`.
        hover_surfaces : bool, optional
            Add the model-top / layer-bottom table to the sectioned hover.
        show_layer_elevs : bool, optional
            Add model-top and per-layer-bottom rows to the hover. Defaults to
            whether the grid actually carries layer elevations (``vor.gdf_topbtm``),
            because forcing it on a grid without them raises.

        Returns
        -------
        Choro or matplotlib.figure.Figure
            With ``backend='plotly'`` (the default), a
            :class:`~myflopy.viz.Picture`: it renders itself in Jupyter, and answers
            ``.fig``, ``.show()``, ``.save(path)`` and ``.html(path)``. It also
            carries ``.plot_mpl()`` for a static rendering and ``.ani`` for the
            animation over periods.

            With ``backend='mpl'``, a bare Matplotlib ``Figure`` -- not a Picture, so
            use ``.savefig(path)`` and ``.axes[0]`` rather than the picture verbs.

            .. end generated

        See Also
        --------
        map : one panel of this, drawn on its own.
        animate : the same panels in sequence rather than side by side.

        Examples
        --------
        >>> model.npf.k.mosaic()                     # one panel per layer
        >>> model.npf.k.mosaic(logscale=True, show_layer_elevs=True)
        >>> model.hds.mosaic(by="layer", contours=True, contour_levels=8)
        >>> model.drn.results.q.mosaic(features="wells.gpkg", ncols=2)
        >>> model.npf.k.mosaic(logscale=True, backend="mpl")
        """

        panel_options = _explicit_options(locals(), type(self).mosaic, PANEL_DRAWING_PARAMS)
        # A LIST of extents is one per PANEL (ledger 309), so it is dealt out
        # below instead of handed whole to every panel.
        extents = (
            list(panel_options.pop("extent"))
            if is_extent_list(panel_options.get("extent"))
            else None
        )
        kind_key = str(kind).lower()
        where = f"{self._noun_path()}.mosaic(kind={kind_key!r})"
        if kind_key == "map":
            axis = (by or self._spatial_default_facet()).lower()
            fixed_per = self._single_period(per)
            panels = self._facet_panels(
                axis, per=fixed_per, layer=layer, model=model, extents=extents,
                **panel_options, **self._map_composer_options(where, cells, agg),
                **kwargs,
            )
            heading = title or f"{self._spatial_value_label()} by {axis}"
            if self._normalize_backend(backend) == "plotly":
                return self._plotly_mosaic(
                    panels, ncols=int(ncols), title=heading, sync_views=sync_views
                )
            return self._mpl_mosaic(panels, ncols=int(ncols), title=heading)
        if extents is not None:
            panel_options["extent"] = extents
        refuse_panel_drawing_options(kind_key, "mosaic", panel_options)
        if kind_key == "plot":
            refuse_series_composer_options(where, "mosaic", kwargs)
            return self._plot_mosaic(
                by=by,
                per=per,
                layer=layer,
                model=model,
                cells=cells,
                agg=agg,
                backend=backend,
                ncols=ncols,
                title=title,
            )
        if kind_key == "section":
            refuse_a_missing_line_profile(self, where)
            return self._xs_mosaic(
                by=by,
                per=per,
                layer=layer,
                model=model,
                backend=backend,
                ncols=ncols,
                title=title,
                **self._section_composer_options(where, cells, agg, kwargs),
            )
        raise ValueError(f"kind must be 'map', 'plot', or 'section', got {kind!r}.")

    def _plot_mosaic(
        self, *, by, per, layer, model, cells, agg, backend, ncols, title
    ):
        """Grid of series panels faceted by model, layer, or an entity."""

        frame, value_column, keys = self._series_selection(
            model=model, per=per, layer=layer, cells=cells
        )
        default_axis = (
            "model"
            if (self._spatial_models() and "model" in frame.columns)
            else next(
                (c for c in self._SERIES_ENTITY_COLUMNS if c in frame.columns),
                "layer",
            )
        )
        axis = str(by or default_axis).lower()
        if axis == "period":
            raise ValueError(
                "Series panels already have period on the x-axis; facet by "
                "'model', 'layer', or an entity (e.g. 'lake') instead."
            )
        if axis not in frame.columns:
            raise ValueError(
                f"Cannot facet series by {axis!r}; available facet columns: "
                f"{[c for c in ('model', 'layer', *self._SERIES_ENTITY_COLUMNS) if c in frame.columns]}."
            )
        if axis == "model":
            in_frame = set(frame["model"])
            facet_values = [n for n in self._resolve_models(model) if n in in_frame]
        else:
            facet_values = sorted(frame[axis].dropna().unique().tolist())
        aggregation = self._series_default_agg() if agg is None else str(agg)
        panel_keys = [key for key in keys if key != axis]
        panels = [
            (
                self._series_line_label([axis], (value,)),
                self._series_lines(
                    frame[frame[axis] == value], panel_keys, value_column, aggregation
                ),
            )
            for value in facet_values
        ]
        heading = title or (
            f"{self._spatial_value_label()} by stress period per {axis}"
            + (" (model - reference)" if self._spatial_is_diff() else "")
        )
        if self._normalize_backend(backend) == "plotly":
            return _xy_mosaic_plotly(
                panels,
                ncols=int(ncols),
                title=heading,
                xaxis_title="Stress Period",
                yaxis_title=value_column,
                markers=True,
            )
        return _xy_mosaic_mpl(
            panels,
            ncols=int(ncols),
            title=heading,
            xlabel="Stress Period",
            ylabel=value_column,
            markers=True,
        )

    def _facet_panels(self, axis, *, per, layer, model, extents=None, **map_kwargs):
        """Build ``[(label, Choro), ...]`` mosaic panels faceted by ``"layer"`` or ``"model"``.

        A layer facet fixes the model (reference) and varies the layer; a model
        facet fixes the layer and varies the model. Both hold ``per`` constant.
        ``extents`` -- one window per panel, in facet order -- frames each on its
        own; its length must match the facet (ledger 309).
        """

        if axis == "layer":
            reference = None if self._spatial_models() is None else self._resolve_models(model)[0]
            facets = [
                (f"Layer {value + 1}", {"layer": value, "model": reference})
                for value in self._resolve_layers(layer)
            ]
        elif axis == "model":
            base_layer = self._single_layer(layer)
            facets = [
                (str(name), {"layer": base_layer, "model": name})
                for name in self._resolve_models(model)
            ]
        else:
            raise ValueError(f"by must be 'layer' or 'model', got {axis!r}.")
        if extents is None:
            windows = [{} for _ in facets]
        elif len(extents) != len(facets):
            raise ValueError(
                f"extent= lists {len(extents)} windows for {len(facets)} panels "
                f"({', '.join(label for label, _ in facets)}); give one per panel, "
                "or a single extent to frame them all."
            )
        else:
            windows = [{} if value is None else {"extent": value} for value in extents]
        return [
            (label, self._spatial_map(per=int(per), **dims, **window, **map_kwargs))
            for (label, dims), window in zip(facets, windows, strict=True)
        ]

    def _plotly_mosaic(self, panels, *, ncols, title, sync_views=None):
        """Compose ``Choro`` panels into one synchronized Plotly small-multiples figure."""

        return viz_mosaic(
            list(panels),
            ncols=int(ncols),
            title=title,
            diff=self._spatial_is_diff(),
            sync_views=sync_views,
        )

    def _mpl_mosaic(self, panels, *, ncols, title):
        """Render ``Choro`` panels as a matplotlib grid sharing one color scale + colorbar.

        Diff surfaces get a symmetric ``RdBu`` scale centered at zero; other
        surfaces take the first panel's own ``colorscale`` over the panels'
        pooled finite range -- through a ``LogNorm`` in DATA units when the
        panels are log-scaled, so the shared bar reads real units (ledger 252).
        """

        from matplotlib import cm

        labels = [label for label, _ in panels]
        choros = [choro for _, choro in panels]
        if not choros:
            raise ValueError("mosaic requires at least one panel.")
        ncols = min(int(ncols), len(choros)) or 1
        nrows = int(np.ceil(len(choros) / ncols))
        fig, axes = mpl_axes(nrows, ncols, figsize=(5.0 * ncols, 4.5 * nrows), squeeze=False)
        _arrays, cmap, norm = _mpl_shared_scale(choros, diff=self._spatial_is_diff())
        for axis, label, choro in zip(axes.flat, labels, choros, strict=False):
            # Through the panel's OWN static renderer, so what one `map(backend=
            # "mpl")` draws -- `features=`, the `select=` outline -- every panel
            # draws. This loop used to re-plot the cells itself and dropped both
            # (ledger 306).
            choro.plot_mpl(
                ax=axis,
                cmap=cmap,
                norm=norm,
                colorbar=False,
                title=str(label),
                linewidth=0.2,
                edgecolor="#666666",
            )
        for axis in axes.flat[len(choros):]:
            axis.set_visible(False)
        mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
        mappable.set_array([])
        colorbar = fig.colorbar(
            mappable,
            ax=list(axes.flat[: len(choros)]),
            shrink=0.85,
            label=self._spatial_value_label(),
        )
        if getattr(choros[0], "number_style", "sig") == "thousands":
            # The panels' own number, 112,000, on the one shared bar (ledger 313).
            from matplotlib.ticker import NullFormatter

            colorbar.ax.yaxis.set_major_formatter(thousands_formatter())
            colorbar.ax.yaxis.set_minor_formatter(NullFormatter())
        fig.suptitle(title)
        return fig

    def _single_period(self, per) -> int:
        """Reduce the ``per`` selector to one period (0 by default; first of a set)."""

        if per is None:
            return 0
        if isinstance(per, (int, np.integer)):
            return int(per)
        return self._resolve_periods(per)[0]

    def _frame_panels(self, axis, *, per, layer, model, **map_kwargs):
        """Build ``[(label, Choro), ...]`` animation frames over ``"period"`` or ``"model"``.

        A period animation fixes the model and layer and steps through periods; a
        model animation fixes the period and layer and steps through models.
        """

        base_layer = self._single_layer(layer)
        if axis == "period":
            reference = self._resolve_models(model)[0]
            return [
                (
                    f"Period {value}",
                    self._spatial_map(per=value, layer=base_layer, model=reference, **map_kwargs),
                )
                for value in self._resolve_periods(per)
            ]
        if axis == "model":
            fixed_per = self._single_period(per)
            return [
                (
                    str(name),
                    self._spatial_map(per=fixed_per, layer=base_layer, model=name, **map_kwargs),
                )
                for name in self._resolve_models(model)
            ]
        raise ValueError(f"over must be 'period' or 'model', got {axis!r}.")

    def animate(
        self,
        *,
        # -- what is swept, and how --------------------------------------------
        kind: str = "map",
        over: str | None = None,
        per=None,
        layer: int = 0,
        model=None,
        cells=None,
        agg: str | None = None,
        title: str | None = None,
        # -- colour of each panel ---------------------------------------------
        zmin: float | None = None,
        zmax: float | None = None,
        colorscale: str | list | tuple | None = None,
        logscale: bool = False,
        # -- contours ----------------------------------------------------------
        contours: bool | str = False,
        contour_values=None,
        contour_levels: int | float | list = 10,
        contour_color: str = "black",
        contour_width: float = 1.5,
        contour_name: str | None = None,
        contour_clip: bool = True,
        contour_resolution: int = 150,
        contour_method: str = "linear",
        # -- highlighting -------------------------------------------------------
        select=None,
        select_style: str = "outline",
        select_color: str | None = None,
        # -- overlays and framing -----------------------------------------------
        features=None,
        hillshade_path=None,
        fit_bounds: bool = True,
        bounds_padding: float = 0.05,
        extent=None,
        # -- hover ---------------------------------------------------------------
        hover=None,
        show_layer_elevs: bool | None = None,
        hover_surfaces: bool | None = None,
        # -- renderer -------------------------------------------------------------
        backend: str = "plotly",
        **kwargs,
    ):
        """Animate a panel across a dimension.

        ``kind`` picks the panel type: ``"map"`` (default), ``"section"`` (cross
        sections, heads leaves), or ``"plot"`` (series). ``over="period"``
        sweeps stress periods (maps and sections); ``over="model"`` sweeps the
        group's models (any kind). Both backends return an object for the caller
        to display or save -- nothing is written to disk.

        **Every drawing option of ``map`` is named here and reaches every frame**
        (8.8, 2026-09-14), and only the ones a caller actually SET are forwarded.
        They are refused, not ignored, on ``kind="plot"``/``"section"``.

        Parameters
        ----------
        kind : {'map', 'plot', 'section'}, default 'map'
            What each frame is: a choropleth, a series plot, or a line profile.
            ``'section'`` exists only on a dependent-variable noun (``model.hds``,
            ``.conc``, ``.temp``, ``group.hds``, ``diff.hds``) and raises
            ``ValueError`` on any other.
        over : {'period', 'model'}, optional
            The dimension swept. ``'period'`` (maps and sections) steps stress
            periods; ``'model'`` steps a group's models. Defaults to
            ``'period'`` for maps and sections, ``'model'`` for series.
        per : int or sequence of int, optional
            Periods to sweep (``over='period'``), or the one period held
            constant otherwise. Defaults to every period.
        layer : int, default 0
            Layer held constant across the frames (0-based).
        model : str or sequence of str, optional
            The models swept (``over='model'``), or the one held constant
            otherwise. Groups and diffs only.
        cells : int or sequence of int, optional
            ``kind='plot'``: restrict the series frames to these cells (zero-based,
            within a layer), one line per cell. ``kind='section'``: each frame's
            path as cell numbers, as :meth:`section` takes it (give it or
            ``line``). ``kind='map'`` takes none -- a map draws every cell, and
            the highlight is ``select=`` -- and raises ``TypeError``.
        agg : str, optional
            ``kind='plot'``: how several rows per period are reduced for a series
            frame. ``kind='map'``: handed to every frame's map, which combines a
            cell's records by it (a record noun's ``agg``) -- on a noun whose map
            has no ``agg`` (heads, a static array) it raises ``TypeError``.
            ``kind='section'`` has nothing to reduce and raises ``TypeError``.
        title : str, optional
            Figure title. Defaults to a label naming the field and the axis.
        backend : {'plotly', 'mpl'}, default 'plotly'
            Which renderer draws the animation. ``'plotly'`` returns an
            interactive house ``viz.Fig`` with play/pause and a slider;
            ``'mpl'`` a :class:`matplotlib.animation.FuncAnimation`. Both are
            objects to display or save -- nothing is written to disk.
        **kwargs
            Options for each frame's builder, by ``kind``. ``'map'``: forwarded,
            with the drawing options named here, to every frame's map -- a single
            model's noun ``map`` (``multiplier``/``fill_value`` on a record noun,
            ``per_timestep``/``bgs``/``hover_layers`` on heads, ..., and through
            its own tail Plotly trace properties), or on a group or a diff what
            :meth:`map` forwards as ``**map_kwargs``. ``per``, ``layer`` and
            ``model`` are this composer's and are set per frame; ``agg`` and
            ``cells`` are named above and never get here.
            ``'section'``: the line profile's own options -- ``line``,
            ``interpolate``, ``interpolator``, ``use_rbf``, ``num_points``,
            ``clip``, ``x_or_y``, ``show_model_top``, ``show_model_btm`` (see
            :meth:`section`) -- and any other keyword raises ``TypeError``: a
            cross-section argument (``fill``, ``layers``, ...) naming
            ``model.plot.section(fill=...)``, ``kstpkper`` naming ``per=``, and
            anything else listing the names accepted. ``'plot'``: none -- the
            series frames read only the parameters named above, so any keyword
            raises ``TypeError`` listing them.

            .. generated from myflopy.plot.map -- edit there, then run scripts/derive_docstrings.py
        zmin, zmax : float, optional
            Fixed color-scale limits. Set both to hold the scale steady across
            frames or panels; ``zmin >= zmax`` raises rather than rendering one flat
            color. With neither, the range comes from the data.
        colorscale : str or list of (float, str), optional
            Any Plotly colorscale name -- ``'earth'`` (the default), ``'Viridis'``,
            ``'Turbo'``, ``'RdBu'``, ...; append ``_r`` to reverse one -- or explicit
            ``[[position, color], ...]`` stops. The static backend builds its
            colormap from the same stops, so a name renders the same way on both.
            An unknown name warns and falls back to ``'earth'``.
        logscale : bool, default False
            Color on a log10 scale, for a field that spans orders of magnitude
            (hydraulic conductivity, recharge). The colorbar is labelled in real
            units at whole decades on both backends and the hover keeps the real
            values. ``zmin``/``zmax`` (and a ``zmid``) are then in DATA units --
            ``zmin=1e-4, zmax=10`` -- and a non-positive limit raises. Non-positive
            values draw as gaps.
        contours : bool or str, default False
            Overlay contour lines. ``True`` contours the mapped values; a string
            names a different field to contour instead.
        contour_values : sequence of float, optional
            Contour a supplied array rather than the mapped one.
        contour_levels : int or float or list of float, default 10
            A count of levels, a fixed interval, or explicit level values.
        contour_color : str, default 'black'
            Line colour for the contour trace.
        contour_width : float, default 1.5
            Line width for the contour trace.
        contour_name : str, optional
            Legend name for the contour trace.
        contour_clip : bool, default True
            Clip contours to the active domain instead of the full grid extent.
        contour_resolution : int, default 150
            Grid size used to build the contours, per axis. Higher is smoother and
            slower. **Only acts under ``contour_method="cubic"``**, which is the one
            that interpolates onto a square grid before contouring; the default
            linear method triangulates the cell centres directly, so there is no
            grid for this to size and passing it changes nothing (measured: 338
            contour points at 40, 150 and 400 alike).
        contour_method : {'linear', 'cubic'}, default 'linear'
            How the scattered cell values become contours. ``'linear'``
            triangulates the cell centres and contours the triangulation --
            fast, exact at the centres, and faceted. ``'cubic'`` interpolates onto
            a ``contour_resolution``-square grid first (Clough-Tocher) and contours
            that -- smoother, slower, and able to overshoot between cells.
            ``'tri'``/``'tricontour'`` and ``'clough'``/``'clough_tocher'``/
            ``'cloughtocher'`` are accepted as aliases. Anything else raises naming
            both; ``'nearest'`` in particular was documented here for a while and
            has never been implemented.
        select : sequence of int, ndarray, str, Path or geometry, optional
            Cells to highlight. Cell indices, a boolean mask of length ``ncpl``, the
            name of a registered model region (``"all_streams"``), a path to a
            vector file, or a shapely/GeoPandas geometry, by coverage. Highlights
            nothing (with a warning) when the selection is empty.
        select_style : {'outline', 'dim', 'both'}, default 'outline'
            How the highlight is drawn. ``'outline'`` traces the dissolved boundary
            of the selection and leaves every cell at full opacity. ``'dim'`` fades
            the *unselected* cells to 20% instead, which suits a bare grid where the
            selection is the subject, but on a field it costs a measured 6.9x of
            readable contrast everywhere you did not select -- and one box/lasso
            gesture in the browser overwrites it. ``'both'`` draws each.
        select_color : str, optional
            Highlight colour, defaulting to ``viz.PALETTE.highlight``. Worth setting
            when the default red collides with a red-blue diverging colorscale.
        features : path, GeoDataFrame, geometry, Feature, or list/dict of those, optional
            Custom geometry drawn over the cells -- wells, pond footprints, a
            boundary. A path to any vector file GeoPandas reads (name a GeoPackage
            layer with ``Feature(path, layer=...)``), a ``GeoDataFrame`` or
            ``GeoSeries``, a shapely geometry, or a :class:`myflopy.plot.Feature`
            carrying its own ``color``/``width``/``size``/``fill``/``label``/
            ``hover``. A list draws several layers; a dict names them. Points
            become markers, lines polylines, polygons outlines (or fills). Each
            layer gets a legend entry, a stable house color and a hover listing its
            attribute columns. A source without a CRS is taken to be in the grid's.
            Drawn by both backends and kept by ``mosaic`` and ``animate``.
        hillshade_path : Path, optional
            A hillshade GeoTIFF to draw beneath the cells for topographic context.
        fit_bounds : bool, default True
            Fit the initial view to the grid extent rather than using ``zoom``. The
            fit is for the figure's own size, or plotly's default 700x450 when it
            names none -- a larger display only adds margin around the grid -- and
            ``.save()`` to an image refits for the export's size.
        bounds_padding : float, default 0.05
            Fractional padding around the fitted bounds.
        extent : tuple, points, geometry, GeoDataFrame, area or str, optional
            The window the map OPENS on, on both backends -- the one way to frame a
            static (matplotlib) map, which cannot be zoomed. Four numbers
            ``(xmin, ymin, xmax, ymax)`` in the grid's CRS (shapely's order, NOT
            FloPy's ``(xmin, xmax, ymin, ymax)``); a run of ``(x, y)`` points (their
            bounding box); any shapely geometry -- a polygon, a transect line -- by its
            bounds; a GeoDataFrame, reprojected from its own CRS; or the name of an
            area of interest (``model.aoi``). It sets the VIEW only: no cell is
            cropped and the colour scale still spans the whole field, so pass
            ``zmin``/``zmax`` to tighten it for a close-up. ``bounds_padding`` pads it
            as it pads the grid. A box off the grid raises, printing both boxes and
            both CRSs. In a ``mosaic``, one extent frames every panel and a LIST of
            them frames each panel on its own.
        hover : HoverSpec, optional
            Replace the sectioned hover outright. See
            :mod:`myflopy.modflow.utils.datatypes.hover`.
        hover_surfaces : bool, optional
            Add the model-top / layer-bottom table to the sectioned hover.
        show_layer_elevs : bool, optional
            Add model-top and per-layer-bottom rows to the hover. Defaults to
            whether the grid actually carries layer elevations (``vor.gdf_topbtm``),
            because forcing it on a grid without them raises.

        Returns
        -------
        Choro or matplotlib.figure.Figure
            With ``backend='plotly'`` (the default), a
            :class:`~myflopy.viz.Picture`: it renders itself in Jupyter, and answers
            ``.fig``, ``.show()``, ``.save(path)`` and ``.html(path)``. It also
            carries ``.plot_mpl()`` for a static rendering and ``.ani`` for the
            animation over periods.

            With ``backend='mpl'``, a bare Matplotlib ``Figure`` -- not a Picture, so
            use ``.savefig(path)`` and ``.axes[0]`` rather than the picture verbs.

            .. end generated

        See Also
        --------
        map : one frame of this, drawn on its own.
        mosaic : the same panels side by side rather than in sequence.

        Examples
        --------
        >>> model.hds.animate()                               # heads over periods
        >>> model.hds.animate(logscale=True, features="wells.gpkg")
        >>> model.drn.results.q.animate(over="period", zmin=-50, zmax=50)
        >>> model.hds.animate(backend="mpl")
        """

        panel_options = _explicit_options(locals(), type(self).animate, PANEL_DRAWING_PARAMS)
        if is_extent_list(panel_options.get("extent")):
            raise TypeError(
                "animate() has ONE view, shared by every frame, so it takes one "
                "extent; a list of them is for a mosaic's panels."
            )
        kind_key = str(kind).lower()
        where = f"{self._noun_path()}.animate(kind={kind_key!r})"
        if kind_key == "map":
            axis = str(over or "period").lower()
            frames = self._frame_panels(
                axis, per=per, layer=layer, model=model, **panel_options,
                **self._map_composer_options(where, cells, agg), **kwargs,
            )
            heading = title or f"{self._spatial_value_label()} over {axis}"
            if self._normalize_backend(backend) == "plotly":
                return self._plotly_animation(frames, title=heading)
            return self._mpl_animation(frames, title=heading)
        refuse_panel_drawing_options(kind_key, "animate", panel_options)
        if kind_key == "plot":
            refuse_series_composer_options(where, "animate", kwargs)
            return self._plot_animation(
                over=over,
                per=per,
                layer=layer,
                model=model,
                cells=cells,
                agg=agg,
                backend=backend,
                title=title,
            )
        if kind_key == "section":
            refuse_a_missing_line_profile(self, where)
            return self._xs_animation(
                over=over,
                per=per,
                layer=layer,
                model=model,
                backend=backend,
                title=title,
                **self._section_composer_options(where, cells, agg, kwargs),
            )
        raise ValueError(f"kind must be 'map', 'plot', or 'section', got {kind!r}.")

    def _plot_animation(
        self, *, over, per, layer, model, cells, agg, backend, title
    ):
        """Flip series panels over the group's models."""

        axis = str(over or "model").lower()
        if axis != "model":
            raise ValueError(
                "Series panels already sweep periods on the x-axis; use "
                "over='model', or animate kind='map'/'section' over periods."
            )
        frame, value_column, keys = self._series_selection(
            model=model, per=per, layer=layer, cells=cells
        )
        if not (self._spatial_models() and "model" in frame.columns):
            raise ValueError("over='model' requires a model axis (groups/diffs).")
        aggregation = self._series_default_agg() if agg is None else str(agg)
        panel_keys = [key for key in keys if key != "model"]
        frames = [
            (
                str(name),
                self._series_lines(
                    frame[frame["model"] == name], panel_keys, value_column, aggregation
                ),
            )
            for name in self._resolve_models(model)
        ]
        heading = title or (
            f"{self._spatial_value_label()} by stress period over model"
            + (" (model - reference)" if self._spatial_is_diff() else "")
        )
        if self._normalize_backend(backend) == "plotly":
            return _xy_animation_plotly(
                frames,
                title=heading,
                xaxis_title="Stress Period",
                yaxis_title=value_column,
                markers=True,
            )
        return _xy_animation_mpl(
            frames,
            title=heading,
            xlabel="Stress Period",
            ylabel=value_column,
            markers=True,
        )

    def _plotly_animation(self, frames, *, title):
        """Delegate to the shared combinator -- one implementation, two spellings.

        This method WAS that implementation, privately. `plot.animate` promoted
        it (8.6a), and delegating rather than keeping a copy fixed a real bug on
        the way: this version passed only `choro.get_choropleth()`, so contours,
        location markers and pathlines drawn over a map silently vanished
        between the static picture and its animation.
        """

        return _build_frame_figure(frames, title=title)

    def _mpl_animation(self, frames, *, title):
        """Build a matplotlib ``FuncAnimation`` over ``[(label, Choro), ...]`` map frames.

        Uses one fixed color scale across all frames (symmetric ``RdBu`` for diff
        surfaces; otherwise the first frame's own ``colorscale`` over the pooled
        range, through a ``LogNorm`` in data units when the frames are
        log-scaled).
        """

        from matplotlib.animation import FuncAnimation

        if not frames:
            raise ValueError("animate requires at least one frame.")
        labels = [str(label) for label, _ in frames]
        choros = [choro for _, choro in frames]
        _arrays, cmap, norm = _mpl_shared_scale(choros, diff=self._spatial_is_diff())
        fig, axis = mpl_axes(1, 1, figsize=(7.0, 6.0))

        def _draw(index):
            """Redraw the map for animation frame ``index`` under the shared color scale."""

            axis.clear()
            # Through the frame's own static renderer, as `_mpl_mosaic` does.
            choros[index].plot_mpl(
                ax=axis,
                cmap=cmap,
                norm=norm,
                colorbar=False,
                title=f"{title} - {labels[index]}",
                linewidth=0.2,
                edgecolor="#666666",
            )

        _draw(0)
        return FuncAnimation(fig, _draw, frames=len(frames), interval=600, blit=False)


class LineProfileView(SpatialView):
    """A :class:`SpatialView` that has a LINE PROFILE: the ``section`` verb.

    Only a dependent-variable noun -- heads, concentration, temperature, on one
    model (``model.hds``), a group (``group.hds``) or a diff (``diff.hds``) -- has
    a value at every cell of every layer to profile along a line, so only these
    answer ``section`` and ``mosaic``/``animate(kind="section")``. Every other
    noun inherited a ``section`` that raised ``NotImplementedError`` on every
    field (50 of 50 on the canonical model, 14 of 14 on a group of it): a verb
    that can never act, now not offered (ledger 462). A host supplies
    :meth:`_sections`.
    """

    # -- the hook the host provides (the section verb + kind="section" composers) --
    def _sections(self, model=None, **kwargs):
        """Return ``{name: XSection}`` for this node; every host overrides it."""

        raise NotImplementedError(f"{type(self).__name__} must define _sections().")

    def _panel_map_names_agg(self) -> bool:
        """``False``: a dependent-variable panel draws one value per cell, on every tier.

        A group's or a diff's ``_spatial_map`` only passes its tail on, which
        :meth:`SpatialView._panel_map_names_agg` cannot see through, so
        ``group.hds.mosaic(agg=...)`` reached Plotly as a trace property and raised
        its ``ValueError`` where ``model.hds`` raised the ``TypeError`` the
        composers document (review of ledger 461).
        """

        return False

    @staticmethod
    def _profile_axis_labels(sections) -> tuple[str, str]:
        """``(xlabel, ylabel)`` of a profile drawn from ``sections``: the first one's.

        :meth:`XSection._axis_labels <myflopy.modflow.utils.datatypes.xsections.XSection._axis_labels>`
        -- the model's own length unit, and ``Concentration``/``Temperature`` on a
        transport field. The composers wrote ``Distance``/``Elevation`` whatever
        the field (review of ledger 464).
        """

        first = next(iter(sections.values()), None)
        return first._axis_labels() if first is not None else ("Distance", "Elevation")

    def section(
        self,
        model=None,
        *,
        line=None,
        cells=None,
        per: int | None = None,
        kstpkper: tuple[int, int] | None = None,
        layer=0,
        interpolate: bool = False,
        interpolator: str | None = None,
        use_rbf: bool = True,
        num_points: int = 100,
        clip=None,
        x_or_y: str = "x",
        show_model_top: bool = True,
        show_model_btm: bool = False,
        backend: str = "plotly",
        title: str | None = None,
        **refused,
    ):
        """This field's LINE PROFILE: its value against distance along a line.

        One line per ``layer``, plus the model top. On a group ``model=None``
        overlays every member; on a diff, the reference and each compared model.
        This is the picture ``model.plot.section(fill=None)`` draws. The model's
        CROSS-SECTION -- the cells the line cuts, coloured by layer or by a field --
        is ``model.plot.section(line=..., fill=...)``; its arguments (``fill``,
        ``layers``, ``head_layers``, ...) RAISE here and say so, rather than being
        taken and dropped (ledger 334).

        Parameters
        ----------
        model : str, optional
            On a group or a diff, one member's name, to draw that member alone.
            A single model's noun takes none.
        line : LineString, sequence of (x, y), or dict, optional
            The section line: a shapely ``LineString``, ``[(x, y), ...]``, or a
            flopy ``{"line": [...]}`` dict. Give ``line`` or ``cells``.
        cells : int or list of int, optional
            The section path as cell numbers, in order; the line runs through
            their centroids. One cell draws a straight line through it, east-west
            or north-south as ``x_or_y`` says.
        per : int, optional
            Stress period (zero-based) to draw, at its LAST saved step -- on each
            member's own clock on a group or a diff. With neither ``per`` nor
            ``kstpkper``, the first saved time.
        kstpkper : tuple of (int, int), optional
            The exact ``(time step, period)`` to draw. Giving it with ``per``
            raises, as does a time the output did not save.
        layer : int or list of int, default 0
            Whose values to draw, zero-based; a list overlays one line per layer.
            An interpolated profile draws one, so a list of more raises there.
        interpolate : bool, default False
            Sample an interpolated surface along the line, rather than stepping
            from cell centre to cell centre.
        interpolator : str, optional
            The interpolation method by name -- ``'griddata'``, ``'rbf'``,
            ``'linearND'`` or ``'cloughTocher2D'``. Giving one turns
            ``interpolate`` on.
        use_rbf : bool, default True
            With ``interpolate`` and no ``interpolator``: radial basis functions
            (``True``) or ``griddata`` (``False``). ``False`` without
            interpolation raises.
        num_points : int, default 100
            How many points to sample along the line when interpolating;
            changed without interpolation, it raises.
        clip : Path or geometry, optional
            Clip the interpolated surface to this region; raises without
            interpolation.
        x_or_y : {'x', 'y'}, default 'x'
            With a single cell in ``cells``, run the line east-west (``'x'``) or
            north-south (``'y'``) through it.
        show_model_top : bool, default True
            Draw the model top.
        show_model_btm : bool, default False
            Draw the bottoms of the layers in ``layer``.
        backend : {'plotly', 'mpl'}, default 'plotly'
            ``'plotly'`` returns an interactive :class:`~myflopy.viz.Fig`;
            ``'mpl'`` a Matplotlib ``Figure``.
        title : str, optional
            The figure title.
        **refused
            Nothing is accepted here. It exists so a keyword this noun does not
            take raises a ``TypeError`` that says why: the cross-section's
            arguments name ``model.plot.section(fill=...)``, and anything else is
            named with its nearest accepted spelling.

        Returns
        -------
        myflopy.viz.Fig or matplotlib.figure.Figure

        Raises
        ------
        TypeError
            For any keyword not named above -- the cross-section's pointing at
            ``model.plot.section(fill=...)``, a typo with its nearest spelling.
        ValueError
            ``use_rbf``/``num_points``/``clip`` without interpolation, naming
            ``interpolate=True``; a list of layers when interpolating; ``per``
            with ``kstpkper``; a period or time the output did not save.

        Examples
        --------
        >>> model.hds.section(line=line)
        >>> model.hds.section(line=line, layer=[0, 1], per=3)
        >>> model.hds.section(cells=[1653, 651, 1241], interpolator="griddata")
        >>> model.hds.section(line=line, backend="mpl").savefig("heads.png")
        """

        refuse_section_arguments(
            f"{self._noun_path()}.section()", refused, NOUN_SECTION_PARAMS
        )
        from myflopy.modflow.utils.datatypes.xsections import render_xsections

        sections = self._sections(
            model=model, line=line, cells=cells, per=per, kstpkper=kstpkper,
            layer=layer, interpolate=interpolate, interpolator=interpolator,
            use_rbf=use_rbf, num_points=num_points, clip=clip, x_or_y=x_or_y,
            show_model_top=show_model_top, show_model_btm=show_model_btm,
        )
        return render_xsections(sections, backend=backend, title=title)

    @staticmethod
    def _section_lines(sections) -> list:
        """Flatten sections to ``[(series_label, x, y)]`` at their current period."""

        from myflopy.modflow.utils.datatypes.xsections import combined_section_frame

        data = combined_section_frame(sections)
        return [
            (
                str(name),
                sub["distance"].to_numpy(dtype=float),
                sub["elevation"].to_numpy(dtype=float),
            )
            for name, sub in data.groupby("series", sort=False)
        ]

    @staticmethod
    def _period_end_kstpkper(model, period: int):
        """Return the last ``(kstp, kper)`` saved in ``period`` (or None if it saved none).

        The one rule, :func:`~myflopy.modflow.mf6.period_end.last_saved_steps`.
        """

        return last_saved_steps(model.kstpkper).get(int(period))

    def _xs_mosaic(
        self, *, by, per, layer, model, backend, ncols, title, **kwargs
    ):
        """Grid of cross-section panels faceted by model or period."""

        axis = str(by or ("model" if self._spatial_models() else "period")).lower()
        base_layer = 0 if layer is None else layer
        panels = []
        labels = ("Distance", "Elevation")
        if axis == "model":
            for name in self._resolve_models(model):
                sections = self._sections(
                    model=name, per=self._single_period(per), layer=base_layer, **kwargs
                )
                panels.append((str(name), self._section_lines(sections)))
                labels = self._profile_axis_labels(sections)
        elif axis == "period":
            sections = self._sections(model=model, layer=base_layer, **kwargs)
            for period in self._resolve_periods(per):
                for section in sections.values():
                    key = self._period_end_kstpkper(section.model, period)
                    if key is not None:
                        section.kstpkper = key
                panels.append((f"Period {period}", self._section_lines(sections)))
            labels = self._profile_axis_labels(sections)
        else:
            raise ValueError(f"by must be 'model' or 'period' for kind='section', got {axis!r}.")
        heading = title or f"{self._spatial_value_label()} cross sections by {axis}"
        if self._normalize_backend(backend) == "plotly":
            return _xy_mosaic_plotly(
                panels,
                ncols=int(ncols),
                title=heading,
                xaxis_title=labels[0],
                yaxis_title=labels[1],
                markers=False,
            )
        return _xy_mosaic_mpl(
            panels,
            ncols=int(ncols),
            title=heading,
            xlabel=labels[0],
            ylabel=labels[1],
            markers=False,
        )

    def _xs_animation(self, *, over, per, layer, model, backend, title, **kwargs):
        """Animate cross sections over stress periods or the group's models."""

        axis = str(over or "period").lower()
        base_layer = 0 if layer is None else layer
        frames = []
        labels = ("Distance", "Elevation")
        if axis == "period":
            sections = self._sections(model=model, layer=base_layer, **kwargs)
            for period in self._resolve_periods(per):
                for section in sections.values():
                    key = self._period_end_kstpkper(section.model, period)
                    if key is not None:
                        section.kstpkper = key
                frames.append((f"Period {period}", self._section_lines(sections)))
            labels = self._profile_axis_labels(sections)
        elif axis == "model":
            for name in self._resolve_models(model):
                sections = self._sections(
                    model=name, per=self._single_period(per), layer=base_layer, **kwargs
                )
                frames.append((str(name), self._section_lines(sections)))
                labels = self._profile_axis_labels(sections)
        else:
            raise ValueError(
                f"over must be 'period' or 'model' for kind='section', got {axis!r}."
            )
        heading = title or f"{self._spatial_value_label()} cross section over {axis}"
        if self._normalize_backend(backend) == "plotly":
            return _xy_animation_plotly(
                frames,
                title=heading,
                xaxis_title=labels[0],
                yaxis_title=labels[1],
                markers=False,
            )
        return _xy_animation_mpl(
            frames,
            title=heading,
            xlabel=labels[0],
            ylabel=labels[1],
            markers=False,
        )


class LeafFieldSugar:
    """``field=`` keyword and ``.field(name)`` nodes on an aggregator *leaf*.

    For hosts that are themselves a :class:`SpatialView` leaf (their verbs draw
    a default field) but carry several numeric input fields -- e.g. a grouped
    GHB input accessor with ``bhead`` and ``cond``. ``field=None`` keeps the
    leaf's own verb; ``field="cond"`` (or the ``inputs.cond`` attribute)
    dispatches to a field-pinned node from :meth:`_field_node`. Namespace-style
    hosts without their own verbs use :class:`FieldMappable` instead.

    Place this mixin *before* the SpatialView base in the MRO.
    """

    def _field_names(self) -> list[str]:
        """Return the selectable field names (override per host)."""

        raise NotImplementedError

    def field_names(self) -> list[str]:
        """List the mappable field names available on this node."""

        return list(self._field_names())

    def _field_node(self, name: str):
        """Return this leaf pinned to one field (override per host)."""

        raise NotImplementedError

    def field(self, name: str):
        """Return the field-pinned node for ``name`` (validated).

        The same node ``node.<name>`` reaches; ``map(field=name)`` and the other
        verbs dispatch through it.

        Parameters
        ----------
        name : str
            The field, by name -- one of :meth:`field_names`. Passed through
            ``str()`` and lower-cased first, so ``'COND'`` is ``'cond'``. A name
            not listed raises ``ValueError`` naming the choices; when the package
            lists no fields at all, the name is handed on unchecked.

        Returns
        -------
        SpatialView
            This leaf pinned to the one field.
        """

        key = str(name).lower()
        available = self._field_names()
        if available and key not in available:
            raise ValueError(
                f"Unknown field {name!r} for this package; choose from {available}."
            )
        return self._field_node(key)

    #: attributes _field_names itself may read -- never treat these as fields
    #: (prevents __getattr__ recursion when they are genuinely missing)
    _FIELD_SUGAR_GUARD = ("package_name", "group", "model", "field_name")

    def __getattr__(self, name: str):
        """Resolve ``node.<field>`` to that field's pinned node (else raise ``AttributeError``).

        Only names that are actual field names (and not private or guarded
        attributes) are dispatched, so a genuinely missing attribute still errors
        cleanly instead of recursing.
        """

        if (
            not name.startswith("_")
            and name not in self._FIELD_SUGAR_GUARD
            and name in self._field_names()
        ):
            return self._field_node(name)
        raise AttributeError(
            f"{type(self).__name__!s} has no attribute or input field {name!r}"
        )

    def map(self, *args, field=None, **kwargs):
        """Map this leaf's default field, or the one named by ``field=``.

        Parameters
        ----------
        *args
            Passed positionally to the map -- the leaf's own ``map`` or, with
            ``field=``, the pinned field's: on a group or a diff its one positional
            is ``model``.
        field : str, optional
            The input field to draw instead of the leaf's default, by name -- one
            of :meth:`field_names` (``'cond'``, ``'bhead'``, ...), matched
            case-insensitively; ``field="cond"`` is ``node.cond.map(...)``. A
            name not listed raises ``ValueError`` naming the choices. ``None`` (the
            default) draws the leaf's own default field.
        **kwargs
            Passed to the same ``map``, whose parameters they must be (see
            :meth:`SpatialView.map`).
        """

        if field is not None:
            return self.field(field).map(*args, **kwargs)
        return super().map(*args, **kwargs)

    def plot(self, *args, field=None, **kwargs):
        """Series-plot this leaf's default field, or the one named by ``field=``.

        Parameters
        ----------
        *args
            Passed positionally to the series plot -- the leaf's own ``plot`` or, with
            ``field=``, the pinned field's: on a group or a diff its one positional
            is ``model``.
        field : str, optional
            The input field to draw instead of the leaf's default, by name -- one
            of :meth:`field_names` (``'cond'``, ``'bhead'``, ...), matched
            case-insensitively; ``field="cond"`` is ``node.cond.plot(...)``. A
            name not listed raises ``ValueError`` naming the choices. ``None`` (the
            default) draws the leaf's own default field.
        **kwargs
            Passed to the same ``plot``, whose parameters they must be (see
            :meth:`SpatialView.plot`).
        """

        if field is not None:
            return self.field(field).plot(*args, **kwargs)
        return super().plot(*args, **kwargs)

    def mosaic(self, *args, field=None, **kwargs):
        """Mosaic this leaf's default field, or the one named by ``field=``.

        Parameters
        ----------
        *args
            Passed positionally to the leaf's own ``mosaic`` or, with ``field=``,
            the pinned field's. Neither takes a positional, so any raises
            ``TypeError``.
        field : str, optional
            The input field to draw instead of the leaf's default, by name -- one
            of :meth:`field_names` (``'cond'``, ``'bhead'``, ...), matched
            case-insensitively; ``field="cond"`` is ``node.cond.mosaic(...)``. A
            name not listed raises ``ValueError`` naming the choices. ``None`` (the
            default) draws the leaf's own default field.
        **kwargs
            Passed to the same ``mosaic``, whose parameters they must be (see
            :meth:`SpatialView.mosaic`).
        """

        if field is not None:
            return self.field(field).mosaic(*args, **kwargs)
        return super().mosaic(*args, **kwargs)

    def animate(self, *args, field=None, **kwargs):
        """Animate this leaf's default field, or the one named by ``field=``.

        Parameters
        ----------
        *args
            Passed positionally to the leaf's own ``animate`` or, with ``field=``,
            the pinned field's. Neither takes a positional, so any raises
            ``TypeError``.
        field : str, optional
            The input field to draw instead of the leaf's default, by name -- one
            of :meth:`field_names` (``'cond'``, ``'bhead'``, ...), matched
            case-insensitively; ``field="cond"`` is ``node.cond.animate(...)``. A
            name not listed raises ``ValueError`` naming the choices. ``None`` (the
            default) draws the leaf's own default field.
        **kwargs
            Passed to the same ``animate``, whose parameters they must be (see
            :meth:`SpatialView.animate`).
        """

        if field is not None:
            return self.field(field).animate(*args, **kwargs)
        return super().animate(*args, **kwargs)


class DiffSpatialView(SpatialView):
    """:class:`SpatialView` for a diff node: panels are per compared-model deltas.

    The mapped content is a difference (model - reference), so the shared color
    scale is diverging and centered at zero, and faceting defaults to one panel
    per compared model. The host supplies ``_spatial_map(per, layer, model)``
    that returns the delta ``Choro`` for one compared model plus the compared
    model names via ``_spatial_models``.
    """

    def _spatial_is_diff(self) -> bool:
        """Diff nodes always map signed deltas (diverging, zero-centered scale)."""

        return True

    def _spatial_default_facet(self) -> str:
        """Diff mosaics default to one panel per compared model."""

        return "model"


class FieldMappable:
    """Namespace-level ``map`` / ``mosaic`` / ``animate`` with a ``field=`` selector.

    A results/inputs namespace groups several mappable *fields* -- e.g. a lake's
    ``q`` (exchange), ``stage``, ``stage_change``; UZF's ``gwrch`` / ``sat``.
    Each field is a first-class accessor (``ns.stage``) with a :class:`SpatialView`
    surface. This mixin adds the same verbs at the namespace level with a
    ``field=`` selector that simply dispatches to the named accessor -- so
    ``ns.map(field="stage")`` is exactly ``ns.stage.map()``. ``field=None`` uses
    the package's default field; :meth:`field_names` lists the choices.
    """

    _default_field: str = "q"

    def _field_names(self) -> list[str]:
        """Return the selectable field names (override per namespace)."""

        raise NotImplementedError

    def field_names(self) -> list[str]:
        """List the mappable field names available on this namespace."""

        return list(self._field_names())

    def _field_accessor(self, field=None):
        """Return the accessor for ``field`` (or the default field), validating the name."""

        name = field if field is not None else self._default_field
        available = self._field_names()
        if available and name not in available:
            raise ValueError(
                f"Unknown field {field!r} for this package; choose from {available}."
            )
        return getattr(self, name)

    def map(self, *args, field=None, **kwargs):
        """Map one result/input field (``field=`` selects it; default otherwise).

        **The accepted arguments depend on ``field=``**, which is why this takes
        ``*args, **kwargs`` while every leaf verb is explicit. ``ns.map()`` is
        sugar for ``ns.<field>.map()``, and the accessors it dispatches to are
        genuinely different pictures with genuinely different knobs -- a lake's
        ``connections`` map takes ``lake=``/``value_column=`` and has no ``per``;
        a budget map takes ``per=``/``agg=``. Collapsing them into one signature
        would have to either drop arguments or invent ones that raise.

        So: **for completion and type checking, call the leaf** --
        ``ns.stage.map(`` shows its real signature, this does not. Use
        :meth:`field_names` to see the choices.

        Extra positional/keyword args pass through -- e.g. a group namespace
        forwards a positional model name: ``results.map("F9b")``.

        Parameters
        ----------
        *args
            Passed positionally to ``ns.<field>.map(...)`` -- on a group or a
            diff namespace, the model name (``results.map("F9b")``). A single
            model's field map takes no positional, so any raises ``TypeError``.
        field : str, optional
            The field, by name -- one of :meth:`field_names`, matched exactly
            (``'Stage'`` raises where ``'stage'`` works). A name not listed raises
            ``ValueError`` naming the choices. ``None`` (the default) is the
            namespace's default field: ``'q'`` on a package's results,
            ``'gwrch'`` on UZF's, ``'finf'`` on UZF's inputs, and on another
            package's inputs its default value column (``'connection_area'`` on
            LAK, ``'rhk'`` on SFR). On a diff's LAK and SFR results ``'stage'``
            is listed but its node has no verbs, so it raises ``AttributeError``.
        **kwargs
            Passed to ``ns.<field>.map(...)``, whose parameters they must be --
            which differ by field.

        Returns
        -------
        Choro or matplotlib.figure.Figure
            What ``ns.<field>.map(...)`` returns.

        Examples
        --------
        >>> ns.field_names()                 # ['q', 'stage', 'stage_change']
        >>> ns.map()                         # the default field
        >>> ns.map(field="stage")            # same as ns.stage.map()
        >>> ns.stage.map(per=3)              # <- the spelling an editor can help with
        """

        return self._field_accessor(field).map(*args, **kwargs)

    def plot(self, *args, field=None, **kwargs):
        """Series plot of one field by stress period (``field=`` selects it).

        Signature varies with ``field=`` -- see :meth:`map` for why, and call
        ``ns.<field>.plot(`` when you want an editor to complete the arguments.

        Parameters
        ----------
        *args
            Passed positionally to ``ns.<field>.plot(...)``: its one positional is
            ``model`` -- a name (or names) on a group or a diff namespace, and only
            ``None``/``'*'`` on a single model's, where anything else raises
            ``ValueError``.
        field : str, optional
            The field, by name -- one of :meth:`field_names`, matched exactly
            (``'Stage'`` raises where ``'stage'`` works). A name not listed raises
            ``ValueError`` naming the choices. ``None`` (the default) is the
            namespace's default field: ``'q'`` on a package's results,
            ``'gwrch'`` on UZF's, ``'finf'`` on UZF's inputs, and on another
            package's inputs its default value column (``'connection_area'`` on
            LAK, ``'rhk'`` on SFR). On a diff's LAK and SFR results ``'stage'``
            is listed but its node has no verbs, so it raises ``AttributeError``.
        **kwargs
            Passed to ``ns.<field>.plot(...)``, whose parameters they must be --
            which differ by field.
        """

        return self._field_accessor(field).plot(*args, **kwargs)

    def mosaic(self, *args, field=None, **kwargs):
        """Shared-scale mosaic of one field over layers/models.

        Signature varies with ``field=`` -- see :meth:`map` for why, and call
        ``ns.<field>.mosaic(`` when you want an editor to complete the arguments.

        Parameters
        ----------
        *args
            Passed positionally to ``ns.<field>.mosaic(...)``; the field's ``mosaic``
            takes none, so any raises ``TypeError``.
        field : str, optional
            The field, by name -- one of :meth:`field_names`, matched exactly
            (``'Stage'`` raises where ``'stage'`` works). A name not listed raises
            ``ValueError`` naming the choices. ``None`` (the default) is the
            namespace's default field: ``'q'`` on a package's results,
            ``'gwrch'`` on UZF's, ``'finf'`` on UZF's inputs, and on another
            package's inputs its default value column (``'connection_area'`` on
            LAK, ``'rhk'`` on SFR). On a diff's LAK and SFR results ``'stage'``
            is listed but its node has no verbs, so it raises ``AttributeError``.
        **kwargs
            Passed to ``ns.<field>.mosaic(...)``, whose parameters they must be --
            which differ by field.
        """

        return self._field_accessor(field).mosaic(*args, **kwargs)

    def animate(self, *args, field=None, **kwargs):
        """Animate one field over periods/models.

        Signature varies with ``field=`` -- see :meth:`map` for why, and call
        ``ns.<field>.animate(`` when you want an editor to complete the arguments.

        Parameters
        ----------
        *args
            Passed positionally to ``ns.<field>.animate(...)``; the field's ``animate``
            takes none, so any raises ``TypeError``.
        field : str, optional
            The field, by name -- one of :meth:`field_names`, matched exactly
            (``'Stage'`` raises where ``'stage'`` works). A name not listed raises
            ``ValueError`` naming the choices. ``None`` (the default) is the
            namespace's default field: ``'q'`` on a package's results,
            ``'gwrch'`` on UZF's, ``'finf'`` on UZF's inputs, and on another
            package's inputs its default value column (``'connection_area'`` on
            LAK, ``'rhk'`` on SFR). On a diff's LAK and SFR results ``'stage'``
            is listed but its node has no verbs, so it raises ``AttributeError``.
        **kwargs
            Passed to ``ns.<field>.animate(...)``, whose parameters they must be --
            which differ by field.
        """

        return self._field_accessor(field).animate(*args, **kwargs)


__all__ = [
    "_symmetric_color_limit",
    "_blue_white_red_diverging_colorscale",
    "_red_white_blue_diverging_colorscale",
    "_as_layer_cell_property",
    "build_sfr_q_map_payload",
    "build_lak_q_map_payload",
    "build_surface_water_q_map_payload",
    "build_cell_input_map_payload",
    "build_group_input_compare_map_payload",
    "SpatialView",
    "LineProfileView",
    "DiffSpatialView",
    "FieldMappable",
    "LeafFieldSugar",
]

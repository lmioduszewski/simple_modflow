"""Sectioned, styled hover specifications for choropleth maps.

Turns a flat per-cell field payload into a compact, sectioned Plotly hover in the
figs house style: a bold, colored primary value, labeled field blocks, an
optional per-layer / surface table, and a muted footer. The three axes are
independent:

- **blocks** -- which payload fields to show and how to group them
  (:class:`Fields` inline or stacked, :class:`LayerTable` for per-layer values).
- **layers** -- how the main investigated field displays vertically
  (``"active"`` | ``"active+strip"`` | ``"all"``).
- **surfaces** -- grid-derived elevations (model top + each layer bottom at the
  cell), independent of what is plotted; merges with the layer table on head maps.

Everything is precomputed into ``customdata`` strings in Python so per-cell logic
(precision, units, dry-cell marking) is exact; the returned ``hovertemplate`` is a
fixed skeleton of styling + ``%{customdata[i]}`` placeholders. See
:meth:`HoverSpec.render`.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from typing import Any

import numpy as np

# Non-breaking space: Plotly preserves it (regular spaces collapse), so it is the
# reliable way to right-align numbers inside a monospace hover table.
NBSP = " "
# Real minus sign reads cleaner than a hyphen for negative fluxes/heads.
MINUS = "−"


def _is_finite_number(value: Any) -> bool:
    """Return True if ``value`` coerces to a finite float (not NaN/inf/non-numeric)."""

    try:
        return np.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def format_number(value: Any, *, precision: int = 4, unit: str | None = None) -> str:
    """Format one value for a hover: significant digits, real minus, optional unit."""

    if value is None:
        return ""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not np.isfinite(number):
        return ""
    text = f"{number:.{precision}g}"
    if text.startswith("-"):
        text = MINUS + text[1:]
    if unit:
        text = f"{text}{NBSP}{unit}"
    return text


def _pad_right(strings: Sequence[str]) -> list[str]:
    """Left-pad each string with NBSP to the column's max width (right-align)."""

    width = max((len(s) for s in strings), default=0)
    return [NBSP * (width - len(s)) + s for s in strings]


@dataclass(frozen=True)
class HoverStyle:
    """The look of a hover box (figs house style by default)."""

    font_family: str = "Calibri, sans-serif"
    font_size: float = 11.5
    font_color: str = "#222222"
    bgcolor: str = "rgba(255,255,255,0.96)"
    bordercolor: str = "#b5b5b5"
    accent: str = "#0f6e56"
    muted: str = "#888780"
    align: str = "left"

    def to_hoverlabel(self) -> dict:
        """Return the Plotly ``hoverlabel`` mapping for the trace/layout."""

        return {
            "bgcolor": self.bgcolor,
            "bordercolor": self.bordercolor,
            "align": self.align,
            "font": {
                "family": self.font_family,
                "size": self.font_size,
                "color": self.font_color,
            },
        }

    def input(self) -> HoverStyle:
        """The same style with the blue 'input map' accent."""

        return replace(self, accent="#185FA5")


@dataclass
class HoverContext:
    """Per-cell data a :class:`HoverSpec` renders, as plain arrays (test-friendly).

    ``payload`` holds flat per-cell fields (the choropleth's own hover dict).
    ``layer_fields`` maps a field name to a list of per-layer arrays (e.g.
    ``{"head": [layer0, layer1, ...]}``) for :class:`LayerTable`. ``top`` and
    ``botm`` (per-layer arrays) drive the surfaces column. Scalars like
    ``period``/``date`` render as static footer text; list-valued footer sources
    (``area``, ``cells``) become per-cell customdata.
    """

    ncpl: int
    active_layer: int = 0
    payload: Mapping[str, Sequence[Any]] = field(default_factory=dict)
    layer_fields: Mapping[str, Sequence[Sequence[Any]]] = field(default_factory=dict)
    top: Sequence[Any] | None = None
    botm: Sequence[Sequence[Any]] | None = None
    area: Sequence[Any] | None = None
    cells: Sequence[Any] | None = None
    period: int | None = None
    step: int | None = None
    date: str | None = None
    model_name: str | None = None

    @property
    def nlay(self) -> int:
        """Number of layers implied by the layer fields (or the surfaces), else 0."""

        for arrays in self.layer_fields.values():
            return len(arrays)
        if self.botm is not None:
            return len(self.botm)
        return 0

    def has(self, name: str) -> bool:
        """True if ``name`` is available as a flat field or a per-layer field."""

        return name in self.payload or name in self.layer_fields

    def resolve(self, name: str) -> list:
        """Per-cell values for ``name``: a flat field, or a layer field's active layer."""

        if name in self.payload:
            return list(self.payload[name])
        if name in self.layer_fields:
            arrays = self.layer_fields[name]
            index = min(self.active_layer, len(arrays) - 1)
            return list(arrays[index])
        return ["" for _ in range(self.ncpl)]

    def formatted(self, name: str, *, precision: int, unit: str | None = None) -> list[str]:
        """Per-cell display strings for ``name`` (resolve + :func:`format_number`)."""

        return [format_number(v, precision=precision, unit=unit) for v in self.resolve(name)]


class HoverBlock:
    """A section of a hover. Subclasses append to an assembler in :meth:`build`."""

    def build(self, asm: _HoverAssembler, ctx: HoverContext, spec: HoverSpec) -> None:
        """Append this block's template fragments + customdata columns to ``asm``."""

        raise NotImplementedError


@dataclass(frozen=True)
class Fields(HoverBlock):
    """A labeled group of flat fields: inline ``a · b`` or a right-aligned grid."""

    title: str | None = None
    fields: tuple[str, ...] = ()
    inline: bool = False

    def build(self, asm, ctx, spec) -> None:
        """Emit the present fields as one inline row or a labeled right-aligned grid."""

        present = [f for f in self.fields if ctx.has(f)]
        if not present:
            return
        if self.title:
            asm.section_label(self.title)
        if self.inline:
            pieces = []
            for name in present:
                pieces.append(asm.inline_field(name, ctx, spec))
            asm.text(f"{NBSP}{MINUS}{NBSP}".join(pieces) + "<br>")
            return
        for name in present:
            label = spec.label(name)
            asm.text(f'<span style="color:{spec.style.muted}">{label}</span>{NBSP}{NBSP}')
            asm.value(ctx.formatted(name, precision=spec.precision, unit=spec.unit(name)))
            asm.text("<br>")


@dataclass(frozen=True)
class LayerTable(HoverBlock):
    """Per-layer values as a compact table; optionally merges the surfaces column.

    Rows are layers (plus a ``Top`` row when ``surfaces``). The active layer row
    is accent-bold. With ``mark_dry`` a head below its cell bottom gets a dagger.
    """

    fields: tuple[str, ...] = ("head",)
    layers: str = "all"
    surfaces: bool = False
    mark_dry: bool = True

    def build(self, asm, ctx, spec) -> None:
        """Emit the per-layer table (head column, optional bottom column + Top row).

        Builds NBSP-padded columns so numbers right-align, bolds the active layer,
        and appends a dagger to dry cells when ``mark_dry`` and surfaces are shown.
        Does nothing if the primary field has no per-layer arrays in ``ctx``.
        """

        primary = self.fields[0]
        arrays = ctx.layer_fields.get(primary)
        if not arrays:
            return
        nlay = len(arrays)
        unit = spec.unit(primary)
        show_bot = self.surfaces and ctx.botm is not None
        botm = ctx.botm if show_bot else None

        head_cols = [
            _pad_right([format_number(v, precision=spec.precision) for v in arrays[layer]])
            for layer in range(nlay)
        ]
        bot_cols = None
        top_col = None
        if show_bot:
            bot_cols = [
                _pad_right([format_number(v, precision=spec.precision) for v in botm[layer]])
                for layer in range(nlay)
            ]
            if ctx.top is not None:
                top_col = _pad_right([format_number(v, precision=spec.precision) for v in ctx.top])

        if show_bot:
            asm.text(
                f'<span style="color:{spec.style.muted};font-size:10px">'
                f'{NBSP * 4}{primary}{NBSP * 3}bot</span><br>'
            )
            if top_col is not None:
                asm.text(f'<span style="color:{spec.style.muted}">Top</span>{NBSP}{NBSP}')
                asm.text(NBSP * (max(len(s) for s in head_cols[0])) + f'{NBSP}{NBSP}')
                asm.value(top_col, muted=spec)
                asm.text("<br>")

        for layer in range(nlay):
            is_active = layer == ctx.active_layer
            label = f"L{layer + 1}"
            head_vals = list(head_cols[layer])
            if show_bot and self.mark_dry:
                head_vals = self._mark_dry(head_vals, arrays[layer], botm[layer])
            open_span = f'<span style="color:{spec.style.accent}"><b>' if is_active else ""
            close_span = "</b></span>" if is_active else ""
            label_color = spec.style.accent if is_active else spec.style.muted
            asm.text(f'<span style="color:{label_color}">{open_span}{label}{close_span}</span>{NBSP}{NBSP}')
            asm.value(head_vals, wrap=f"{open_span}{{}}{close_span}")
            if show_bot:
                asm.text(f"{NBSP}{NBSP}")
                asm.value(bot_cols[layer], muted=spec)
            asm.text("<br>")

        if unit:
            asm.text(f'<span style="color:{spec.style.muted};font-size:10px">{unit}</span><br>')

    @staticmethod
    def _mark_dry(formatted: list[str], heads, botm) -> list[str]:
        """Append a dagger to each formatted head that sits below its cell bottom."""

        out = []
        for text, head, bottom in zip(formatted, heads, botm, strict=False):
            dry = (
                _is_finite_number(head)
                and _is_finite_number(bottom)
                and float(head) < float(bottom)
            )
            out.append(text + ("†" if dry else ""))
        return out


@dataclass(frozen=True)
class HoverSpec:
    """A complete hover: header primary value, body blocks, and a footer."""

    primary: str | None = None
    title: str | None = None
    id_field: str | None = "cell"
    layers: str = "active"
    surfaces: bool | str = False
    blocks: tuple[HoverBlock, ...] = ()
    footer: tuple[str, ...] = ("period", "date")
    precision: int = 4
    units: Mapping[str, str] = field(default_factory=dict)
    labels: Mapping[str, str] = field(default_factory=dict)
    style: HoverStyle = field(default_factory=HoverStyle)

    def unit(self, name: str) -> str | None:
        """The display unit configured for field ``name`` (or None)."""

        return self.units.get(name)

    def label(self, name: str) -> str:
        """The display label for field ``name`` (falls back to the name itself)."""

        return self.labels.get(name, name)

    def with_style(self, **overrides: Any) -> HoverSpec:
        """A copy of this spec with the given :class:`HoverStyle` fields overridden."""

        return replace(self, style=replace(self.style, **overrides))

    def with_fields(self, *names: str, title: str | None = None) -> HoverSpec:
        """Append an extra inline :class:`Fields` block (call-site ``hover_fields=``)."""

        return replace(self, blocks=self.blocks + (Fields(title=title, fields=tuple(names), inline=True),))

    def _effective_blocks(self, ctx: HoverContext) -> list[HoverBlock]:
        """The blocks to render, prepending an auto layer/surface table when needed.

        When ``layers="all"`` (or surfaces are requested) and the primary field has
        per-layer data, a :class:`LayerTable` is synthesized ahead of the user's
        explicit ``blocks``; surface-only requests still get a Top/bottom table.
        """

        blocks: list[HoverBlock] = []
        primary = self.primary
        want_layers = primary in ctx.layer_fields and self.layers in ("active+strip", "all")
        want_surfaces = bool(self.surfaces)
        if primary in ctx.layer_fields and (self.layers == "all" or want_surfaces):
            blocks.append(
                # the head-below-bottom dagger is head physics: on only for heads,
                # meaningless for concentration/temperature.
                LayerTable(
                    fields=(primary,),
                    layers="all",
                    surfaces=want_surfaces,
                    mark_dry=(primary == "head"),
                )
            )
        elif want_surfaces and ctx.botm is not None:
            blocks.append(LayerTable(fields=("__surface__",), surfaces=True))
        blocks.extend(self.blocks)
        return blocks

    def render(self, ctx: HoverContext) -> tuple[list, str, dict]:
        """Render to ``(customdata, hovertemplate, hoverlabel)`` for a map trace."""

        asm = _HoverAssembler(ctx.ncpl, self.style)

        title = self.title or (self.label(self.primary) if self.primary else "")
        if title:
            asm.text(f'<span style="font-size:12.5px"><b>{title}</b></span>')
        if self.id_field == "cell" and ctx.cells is not None:
            if title:
                asm.text(f"{NBSP}{NBSP}")
            asm.text(f'<span style="color:{self.style.muted}">cell </span>')
            asm.value([str(c) for c in ctx.cells])
        asm.text("<br>")

        if self.primary and ctx.has(self.primary):
            unit = self.unit(self.primary)
            primary_vals = ctx.formatted(self.primary, precision=self.precision, unit=unit)
            if self.primary in ctx.layer_fields:
                left = f"Layer {ctx.active_layer + 1}"
            else:
                left = self.label(self.primary)
            asm.text(f'<span style="color:{self.style.muted}">{left}</span>{NBSP}{NBSP}')
            asm.value(
                primary_vals,
                wrap=f'<span style="font-size:15px;color:{self.style.accent}"><b>{{}}</b></span>',
            )
            asm.text("<br>")
            if self.layers == "active+strip" and self.primary in ctx.layer_fields:
                asm.layer_strip(self.primary, ctx, self)

        for block in self._effective_blocks(ctx):
            block.build(asm, ctx, self)

        self._render_footer(asm, ctx)
        return asm.finish()

    def _render_footer(self, asm: _HoverAssembler, ctx: HoverContext) -> None:
        """Append the muted footer line joining the configured period/date/area/model bits.

        Scalar sources (period, step, date, model) render as static text; list
        sources (area) become per-cell customdata via :meth:`_HoverAssembler.footer_ref`.
        Emits nothing when none of the footer keys resolve to data.
        """

        pieces: list[str] = []
        deferred_static: list[str] = []
        for key in self.footer:
            if key == "period" and ctx.period is not None:
                deferred_static.append(f"Period {ctx.period}")
            elif key == "step" and ctx.step is not None:
                deferred_static.append(f"step {ctx.step}")
            elif key == "date" and ctx.date is not None:
                deferred_static.append(str(ctx.date))
            elif key == "area" and ctx.area is not None:
                ref = asm.footer_ref([format_number(a, precision=max(self.precision, 5)) for a in ctx.area])
                pieces.append(("area", ref))
            elif key == "model" and ctx.model_name is not None:
                deferred_static.append(str(ctx.model_name))
        parts = list(deferred_static)
        for _, ref in pieces:
            parts.append(ref)
        if not parts:
            return
        joined = f"{NBSP}{MINUS}{NBSP}".join(parts)
        asm.footer_line(joined)


class _HoverAssembler:
    """Collects template fragments + customdata columns for one hover."""

    def __init__(self, ncpl: int, style: HoverStyle):
        """Start an empty assembler for a grid of ``ncpl`` cells with ``style``."""

        self.ncpl = int(ncpl)
        self.style = style
        self._cols: list[list] = []
        self._parts: list[str] = []
        self._footer_pending: list[str] = []

    def text(self, fragment: str) -> None:
        """Append a static (non-per-cell) template fragment; empty fragments are skipped."""

        if fragment:
            self._parts.append(fragment)

    def _register(self, values: Sequence[Any]) -> str:
        """Store a per-cell column and return its ``%{customdata[i]}`` placeholder.

        Raises ``ValueError`` if the column length does not match the grid's cell
        count, since Plotly indexes customdata positionally per point.
        """

        column = [("" if v is None else str(v)) for v in values]
        if len(column) != self.ncpl:
            raise ValueError(
                f"hover column has {len(column)} values but grid has {self.ncpl} cells"
            )
        index = len(self._cols)
        self._cols.append(column)
        return "%{customdata[" + str(index) + "]}"

    def value(self, values: Sequence[Any], *, wrap: str = "{}", muted: HoverSpec | None = None) -> None:
        """Append a per-cell value column, optionally wrapped in markup / muted color.

        ``wrap`` is a template whose ``{}`` is replaced by the customdata placeholder
        (e.g. a bold accent span); ``muted`` additionally wraps it in the spec's
        muted color.
        """

        ref = self._register(values)
        text = wrap.replace("{}", ref)
        if muted is not None:
            text = f'<span style="color:{muted.style.muted}">{text}</span>'
        self._parts.append(text)

    def inline_field(self, name: str, ctx: HoverContext, spec: HoverSpec) -> str:
        """Return a ``label value`` fragment for one inline field (no trailing break)."""

        label = spec.label(name)
        ref = self._register(ctx.formatted(name, precision=spec.precision, unit=spec.unit(name)))
        return f'<span style="color:{spec.style.muted}">{label}</span>{NBSP}{ref}'

    def section_label(self, title: str) -> None:
        """Append a small, muted, uppercase section heading followed by a line break."""

        self._parts.append(
            f'<span style="color:{self.style.muted};font-size:10px">{title.upper()}</span><br>'
        )

    def layer_strip(self, primary: str, ctx: HoverContext, spec: HoverSpec) -> None:
        """Append a one-line ``L2 v · L3 v`` strip of the non-active layers' values.

        Used by the ``layers="active+strip"`` mode to show neighbor layers compactly
        beneath the active-layer primary value.
        """

        arrays = ctx.layer_fields[primary]
        pieces = []
        for layer in range(len(arrays)):
            if layer == ctx.active_layer:
                continue
            ref = self._register(
                [format_number(v, precision=spec.precision) for v in arrays[layer]]
            )
            pieces.append(f"L{layer + 1}{NBSP}{ref}")
        if pieces:
            joined = f"{NBSP}{MINUS}{NBSP}".join(pieces)
            self._parts.append(f'<span style="color:{spec.style.muted};font-size:10.5px">{joined}</span><br>')

    def footer_ref(self, values: Sequence[Any]) -> str:
        """Register a per-cell footer column (e.g. cell area) as ``cell <ref> ft²``."""

        return "cell " + self._register(values) + " ft²"

    def footer_line(self, joined: str) -> None:
        """Append the assembled footer text as a small muted span (no line break)."""

        self._parts.append(
            f'<span style="color:{self.style.muted};font-size:10.5px">{joined}</span>'
        )

    def finish(self) -> tuple[list, str, dict]:
        """Close the template and return ``(customdata, hovertemplate, hoverlabel)``.

        Transposes the registered columns into one customdata row per cell, appends
        Plotly's ``<extra></extra>`` to suppress the default trace box, and returns
        the style's hoverlabel mapping.
        """

        template = "".join(self._parts) + "<extra></extra>"
        if self._cols:
            customdata = [list(row) for row in zip(*self._cols, strict=False)]
        else:
            customdata = [[] for _ in range(self.ncpl)]
        return customdata, template, self.style.to_hoverlabel()


# --- default specs per plot kind -------------------------------------------------

def head_hover(*, layers: str = "active+strip", surfaces: bool = False) -> HoverSpec:
    """Default heads hover: active layer primary, optional strip/table/surfaces."""

    return HoverSpec(
        primary="head",
        title="Head",
        layers=layers,
        surfaces=surfaces,
        footer=("period", "date"),
        units={"head": "ft"},
        labels={"head": "Head"},
    )


def conc_hover(*, unit: str = "mg/L", layers: str = "active+strip", surfaces: bool = False) -> HoverSpec:
    """Default GWT concentration hover (mirrors :func:`head_hover`).

    ``unit`` is model-dependent (mass/volume); the default ``"mg/L"`` is a common
    groundwater convention, not a physical law -- thread the model's real unit
    through when known.
    """

    return HoverSpec(
        primary="conc",
        title="Concentration",
        layers=layers,
        surfaces=surfaces,
        footer=("period", "date"),
        units={"conc": unit},
        labels={"conc": "Concentration"},
    )


def temp_hover(*, unit: str = "°C", layers: str = "active+strip", surfaces: bool = False) -> HoverSpec:
    """Default GWE temperature hover (mirrors :func:`head_hover`).

    ``unit`` defaults to ``"°C"`` by convention; thread the model's real unit
    through when known.
    """

    return HoverSpec(
        primary="temp",
        title="Temperature",
        layers=layers,
        surfaces=surfaces,
        footer=("period", "date"),
        units={"temp": unit},
        labels={"temp": "Temperature"},
    )


_RESIDUAL_LABELS = {
    # Named, not signed into agreement: PEST's own .res file reports
    # measured - modelled, myflopy reports simulated - measured (the sign
    # `phi_contributions` already uses). Spelling it out in the label is what
    # keeps a reader from guessing which frame a red cell is in.
    "residual": "residual (sim − meas)",
    "measured": "measured",
    "simulated": "simulated",
    "weight": "weight",
    "zone": "zone",
    "location": "location",
    "n": "observations",
}


def residual_hover(
    *,
    title: str | None = None,
    fields: Sequence[str] = ("measured", "simulated", "weight"),
    units: Mapping[str, str] | None = None,
    labels: Mapping[str, str] | None = None,
) -> HoverSpec:
    """Hover for an observation-residual map (PEST/IES ``plot_obs_residuals``).

    ``footer=()`` for the same reason :func:`parameter_field_hover` documents:
    a residual summarized over every time an observation was made does not
    belong to one stress period, so a "Period 0" line would assert one.
    """

    body = tuple(name for name in fields if name != "residual")
    return HoverSpec(
        primary="residual",
        title=title or "residual",
        blocks=(Fields(fields=body),) if body else (),
        footer=(),
        units=dict(units or {}),
        labels={**_RESIDUAL_LABELS, **(labels or {})},
    )


_PARAMETER_FIELD_LABELS = {
    "change": "posterior / prior",
    "prior_mean": "prior mean",
    "prior_std": "prior sd",
    "posterior_mean": "posterior mean",
    "posterior_std": "posterior sd",
    "base": "base realization",
}


def parameter_field_hover(
    stat_column: str,
    *,
    title: str | None = None,
    fields: Sequence[str] = ("prior_mean", "posterior_mean", "posterior_std"),
    units: Mapping[str, str] | None = None,
    labels: Mapping[str, str] | None = None,
) -> HoverSpec:
    """Hover for a calibrated parameter-field map (PEST/IES ``plot_field``).

    ``footer=()`` for the same reason :func:`result_hover` documents for PRT
    maps: a parameter field summarizes a whole calibration rather than being
    drawn at one stress period, so a "Period 0" line would assert a period the
    value does not belong to. Run identity (which iterations, how many
    realizations) does NOT go in the footer either -- :meth:`HoverSpec._render_footer`
    recognizes only period/step/date/area/model and silently drops anything
    else -- so ``plot_field`` puts it in the figure title, where the matplotlib
    backend can show it too.

    ``stat_column`` is dropped from ``fields`` so the plotted statistic is not
    printed twice, once as the primary line and again in the block.
    """

    body = tuple(name for name in fields if name != stat_column)
    return HoverSpec(
        primary=stat_column,
        title=title or stat_column,
        blocks=(Fields(fields=body),) if body else (),
        footer=(),
        units=dict(units or {}),
        labels={**_PARAMETER_FIELD_LABELS, **(labels or {})},
    )


def lak_hover() -> HoverSpec:
    """Default lake-exchange hover: q/area primary with a lake feature block."""

    return HoverSpec(
        primary="q_per_area",
        title="Lake exchange",
        blocks=(Fields(title=None, fields=("stage", "q", "flow_area", "claktype")),),
        footer=("period", "date", "area"),
        units={"q_per_area": "ft/d", "q": "ft³/d", "flow_area": "ft²", "stage": "ft"},
        labels={"q_per_area": "q / area", "flow_area": "area", "claktype": "type"},
    )


def sfr_hover() -> HoverSpec:
    """Default stream-exchange hover: q/length primary with a reach block."""

    return HoverSpec(
        primary="q_per_length",
        title="Stream exchange",
        blocks=(Fields(title=None, fields=("stage", "depth", "q", "rlen")),),
        footer=("period", "date"),
        units={"q_per_length": "ft²/d", "q": "ft³/d", "rlen": "ft", "stage": "ft", "depth": "ft"},
        labels={"q_per_length": "q / length"},
    )


def cell_input_hover(
    value_column: str,
    *,
    extra_fields: Sequence[str] = (),
    context_fields: Sequence[str] = ("Layer", "Record Count"),
) -> HoverSpec:
    """Default input-map hover (DRN/CHD/GHB/RCH): the input value + its fields.

    ``extra_fields`` are the package's OTHER record fields -- reading a GHB map
    means asking "what head, against what conductance", and a hover carrying only
    the coloured one answers half the question. Callers pass the siblings of
    ``value_column``; the untitled block matches :func:`lak_hover` / :func:`sfr_hover`,
    which list a feature's fields the same way.

    ``context_fields`` are cell facts rather than package fields, so they render
    as one muted inline row: which LAYER you are looking at (the commonest reason
    a boundary map comes back empty is that the package is not in the default
    layer 0), and how many records were aggregated into the cell, which is the
    only thing that explains a summed value. Both are skipped when the payload
    does not carry them.
    """

    blocks: list[HoverBlock] = []
    if extra_fields:
        blocks.append(Fields(title=None, fields=tuple(extra_fields)))
    if context_fields:
        blocks.append(Fields(title=None, fields=tuple(context_fields), inline=True))
    return HoverSpec(
        primary=value_column,
        title=value_column,
        blocks=tuple(blocks),
        footer=("period", "date"),
        style=HoverStyle().input(),
        labels={"Record Count": "records", "Layer": "layer"},
    )


def result_hover(
    value_column: str,
    *,
    title: str | None = None,
    extra_fields: Sequence[str] = (),
    units: Mapping[str, str] | None = None,
    labels: Mapping[str, str] | None = None,
    footer: Sequence[str] = ("period", "date"),
) -> HoverSpec:
    """Generic result-map hover (per-package cell budget ``q``, stage, UZF fields).

    ``footer`` keeps the period/date line by default. Pass ``footer=()`` for maps
    that are **time-integrated** rather than drawn at one stress period -- PRT
    travel-time / endpoint maps summarize a whole run, so a "Period 0" footer
    would assert a period the value does not belong to.
    """

    return HoverSpec(
        primary=value_column,
        title=title or value_column,
        blocks=(Fields(fields=tuple(extra_fields)),) if extra_fields else (),
        footer=tuple(footer),
        units=dict(units or {}),
        labels=dict(labels or {}),
    )


def pathline_hover(
    *,
    particle: str | None = None,
    unit: str | None = None,
    fields: Sequence[str] = ("layer", "z", "release_group"),
) -> HoverSpec:
    """Hover for one particle's track on the pathline map.

    Unlike every other spec here this renders **per vertex of a trajectory**, not
    per grid cell -- the hover context is built over the points of one particle's
    line, so ``cell`` is the cell that particle was in at that point and the
    elapsed time is its own. ``particle`` names the trace in the header (the id
    is fixed for the whole line); ``unit`` is the flow model's TDIS time unit.
    """

    return HoverSpec(
        primary="travel_time",
        title=f"Particle {particle}" if particle else "Particle",
        blocks=(Fields(fields=tuple(fields)),) if fields else (),
        footer=(),
        units={"travel_time": unit} if unit else {},
        labels={
            "travel_time": "elapsed",
            "release_group": "group",
            "z": "elevation",
            "t": "time",
        },
    )


def compare_hover(
    value_column: str,
    diff_column: str,
    *,
    title: str | None = None,
    units: Mapping[str, str] | None = None,
    labels: Mapping[str, str] | None = None,
) -> HoverSpec:
    """Default diff-map hover for group ``compare_map`` payloads.

    Compare payloads (:func:`build_group_input_compare_map_payload`) carry
    ``value_column``, ``reference_{value_column}``, ``diff_column``, plus
    per-cell ``Model`` / ``Reference Model`` strings. The primary is the mapped
    difference; the block shows the two values it came from. The purple accent
    signals "this map is a delta" next to teal result and blue input maps.
    """

    reference_column = f"reference_{value_column}"
    merged_labels = {
        diff_column: f"Δ {value_column}",
        value_column: "model",
        reference_column: "reference",
        "Model": "model",
        "Reference Model": "vs",
    }
    merged_labels.update(labels or {})
    return HoverSpec(
        primary=diff_column,
        title=title or f"Δ {value_column} vs reference",
        blocks=(
            Fields(fields=(value_column, reference_column)),
            Fields(fields=("Model", "Reference Model")),
        ),
        footer=("period", "date"),
        units=dict(units or {}),
        labels=merged_labels,
        style=HoverStyle(accent="#534AB7"),
    )


def surface_water_hover() -> HoverSpec:
    """Default combined SFR+LAK exchange hover with the per-source breakdown."""

    return HoverSpec(
        primary="surface_water_exchange",
        title="Surface-water exchange",
        blocks=(Fields(fields=("sfr_exchange", "lak_exchange", "source")),),
        footer=("period", "date"),
        labels={"sfr_exchange": "sfr", "lak_exchange": "lak"},
    )


def default_hover_from_payload(payload: Mapping[str, Sequence[Any]], *, primary: str | None = None) -> HoverSpec:
    """Fallback spec for an arbitrary payload: primary + the remaining fields inline."""

    keys = [k for k in payload if k not in {"Period", "Time Step", "Cell", "Layer"}]
    main = primary or (keys[0] if keys else None)
    rest = tuple(k for k in keys if k != main)
    return HoverSpec(
        primary=main,
        title=main or "Value",
        blocks=(Fields(fields=rest, inline=True),) if rest else (),
    )

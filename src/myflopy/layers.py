"""Layer-centric authoring facade for model discretization.

`LayerStack` is a thin, friendly front door over :class:`~myflopy.surfaces.LayerSurfaces`.
You declare a model top, then ``add`` named layers -- each carrying its own bottom
definition and (optionally) its own ``min_thickness`` / ``pinch`` policy::

    from myflopy.layers import LayerStack, Raster, Flat, Contours

    stack = LayerStack(vor, top=Raster("ground.tif"), length_units="feet")
    stack.add("alluvium", bottom=Raster("allu_bot.tif"))
    stack.add("clay",     thickness=30, min_thickness=1, pinch="passthrough")
    stack.add("bedrock",  bottom=Contours("bedrock.gpkg", z="elev"), fill="propagate")

    result = stack.build()          # result.top / .botm / .idomain / .report()
    disv   = stack.to_disv(vor)     # ready-to-use mf.disv spec

This module adds **no new geometry logic**: it translates "top + N named layers"
into the "N+1 surfaces" list and delegates sampling, reconcile, and pinch-out to
the existing, tested `LayerSurfaces` engine.
"""

from __future__ import annotations

import tempfile
from dataclasses import dataclass, field
from dataclasses import replace as _dc_replace
from pathlib import Path
from typing import Any

import numpy as np
import shapely.geometry as shp

from myflopy._logging import get_logger
from myflopy._optional import require
from myflopy.modflow.mf6.grid.interpolated_surface import InterpolatedSurface
from myflopy.modflow.mf6.grid.plotting import GridPlots, _as_linestring
from myflopy.modflow.mf6.grid.triangle import TriangleGrid
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus
from myflopy.modflow.utils.datatypes.readers import read_shp_gpkg
from myflopy.package_api import disv as disv_spec
from myflopy.surfaces import (
    _MF6_LENGTH_UNITS,
    LayerSurfaces,
    Surface,
    _depends_on_previous,
)
from myflopy.viz import Fig, MplPicture, Picture, VtkScene, mpl_axes

logger = get_logger(__name__)

# Convenience source aliases -- the existing Surface constructors under friendlier
# names for layer authoring (NOT new implementations).
Raster = Surface.raster
Flat = Surface.flat
Contours = Surface.from_contours
Points = Surface.from_points
Array = Surface.from_array
Isopach = Surface.isopach
# Fractional subdivision: the contact a share of the way down to a target.
Toward = Surface.toward
# Surface algebra (lower/upper envelopes, clamping, zone selection).
Min = Surface.minimum
Max = Surface.maximum
Clamp = Surface.clamp
Where = Surface.where

_UNSET = object()

#: What ``section(color_by=...)`` will actually key the fill to. Unknown values
#: used to fall through to the layer-coloured branch and draw a picture that
#: silently was not the one asked for.
_SECTION_COLOR_BY = ("layer", "thickness")

#: Per-cell scalars ``grid(backend="vtk")`` can colour by, each with the
#: colormap that suits it. ``layer`` is discrete -- one flat colour per layer,
#: named in the bar -- and the rest are continuous fields. ``elevation`` is the
#: odd one out and stays for compatibility: it is a POINT array of vertex z, so
#: it ramps smoothly within a cell, where ``top`` and ``botm`` are the flat
#: per-cell contact elevations.
_VTK_GRID_CMAPS = {
    "layer": "tab10",
    "thickness": "viridis",
    "top": "terrain",
    "botm": "terrain",
    "cellid": "tab20",
    "elevation": "terrain",
}


def _section_line_points(line) -> list[tuple[float, float]]:
    """Coerce a section line to a list of ``(x, y)`` points.

    Accepts everything the ``section`` docstrings advertise: a shapely
    ``LineString`` / ``MultiLineString``, a path to a vector file (``.shp`` /
    ``.gpkg``, read with the same reader ``vor.plot.section`` uses), or a plain
    sequence of points.

    Raises here, at the call, rather than deep inside FloPy: a section is a lazy
    Picture, so an unusable line used to survive ``stack.plot.section(...)``
    untouched and fail much later on ``.axes`` -- ``TypeError: 'LineString'
    object is not iterable`` -- pointing nowhere near the argument that caused
    it.
    """

    if isinstance(line, (str, Path)):
        # Not wrapped in a try: a broad handler here would eat the argument
        # validation below it, which is the exact bug already found twice in
        # `read_gpkg` and `contour_line_segments`.
        line = _as_linestring(read_shp_gpkg(Path(line)).union_all())
    if isinstance(line, (shp.LineString, shp.MultiLineString)):
        line = list(_as_linestring(line).coords)

    try:
        points = [tuple(float(v) for v in pt)[:2] for pt in line]
    except TypeError:
        raise TypeError(
            f"a section line must be a LineString, a path to a vector file, or "
            f"a sequence of (x, y) points; got {type(line).__name__}."
        ) from None
    if len(points) < 2 or any(len(pt) != 2 for pt in points):
        raise ValueError(
            f"a section line needs at least two (x, y) points; got {points!r}."
        )
    return points


def _check_layer_name(name, layers) -> None:
    """Reject a layer name that would collide with another surface.

    Names have to be unique because they IDENTIFY things: a contact in
    ``surface_names``, and -- since the 3-D volume draws one actor per layer --
    an actor in the PyVista scene, where a repeat silently REPLACES its
    predecessor rather than adding one. ``"top"`` is reserved for the model top
    for the same reason; without this it failed much later, inside pandas, as
    ``cannot reindex on an axis with duplicate labels``.
    """

    if name == "top":
        raise ValueError(
            "'top' is reserved for the model top surface; name the layer for "
            "the material it is, e.g. 'topsoil' or 'upper_sand'."
        )
    if any(layer.name == name for layer in layers):
        raise ValueError(f"layer {name!r} already exists.")


def _coerce_surface(x) -> Surface:
    """Accept a Surface, or coerce a path -> Raster and a number -> Flat."""

    if isinstance(x, Surface):
        return x
    if isinstance(x, (str, Path)):
        return Surface.raster(x)
    if isinstance(x, (int, float)):
        return Surface.flat(float(x))
    raise TypeError(f"Expected a Surface, path, or number; got {type(x).__name__}.")


def _make_surface(bottom, thickness, fill) -> Surface:
    """Build a layer's bottom surface from ``bottom=`` or ``thickness=``."""

    if (bottom is None) == (thickness is None):
        raise ValueError("Provide exactly one of bottom= or thickness=.")
    if thickness is not None:
        if isinstance(thickness, Surface):
            if thickness.kind not in _THICKNESS_KINDS:
                # Any OTHER Surface was passed through verbatim, and a Surface
                # is an ELEVATION -- so `thickness=Flat(20)` under a top of 100
                # set the bottom to 20 (an 80-thick layer), silently, while the
                # docstring promised a thickness. A raw ndarray at least failed
                # loudly. The three kinds above really do measure from the
                # surface above, so they mean here exactly what they say.
                raise TypeError(
                    f"thickness= got a {thickness.kind!r} surface, which is an "
                    f"ELEVATION -- it would set the layer BOTTOM to those values "
                    f"rather than cut that thickness below the surface above. "
                    f"For a per-cell thickness map write "
                    f"`bottom=Isopach(<map>)`; for a constant, `thickness=20.0`; "
                    f"to place the contact itself, `bottom=<surface>`."
                )
            surface = thickness
        else:
            surface = Surface.constant_thickness(float(thickness))
    else:
        surface = _coerce_surface(bottom)
    if fill is not None:
        surface = _dc_replace(surface, fill=fill)
    return surface


#: Surface kinds that genuinely express "this far below the surface above", and
#: so mean what they say when handed to ``thickness=``. Every other kind is an
#: elevation, and was silently used as the layer BOTTOM.
_THICKNESS_KINDS = ("isopach", "constant_thickness", "offset_below")


def _split_shares(split, name: str) -> list[float]:
    """Normalize a ``split=`` argument to shares of the unit, summing to 1.

    ``None`` -> one layer; an int ``N`` -> N equal shares; a sequence -> those
    shares verbatim.
    """

    if split is None:
        return [1.0]
    if isinstance(split, bool):  # bool is an int; `split=True` means nothing
        raise TypeError(
            f"layer {name!r}: split= takes a count or a list of shares, not a bool."
        )
    if isinstance(split, (int, np.integer)):
        n = int(split)
        if n < 1:
            raise ValueError(f"layer {name!r}: split= must be at least 1, got {n}.")
        return [1.0 / n] * n
    # `str` before the loop: it is iterable, so it would otherwise get as far as
    # float('t') and report a confusing conversion error rather than the type.
    if isinstance(split, str):
        raise TypeError(
            f"layer {name!r}: split= takes a count (split=3) or a sequence of "
            f"shares (split=[0.3, 0.7]); got the string {split!r}."
        )
    try:
        shares = [float(f) for f in split]
    except TypeError:
        raise TypeError(
            f"layer {name!r}: split= takes a count (split=3) or a sequence of "
            f"shares (split=[0.3, 0.7]); got {type(split).__name__}."
        ) from None
    if not shares:
        raise ValueError(f"layer {name!r}: split= got an empty sequence.")
    if any(s <= 0 for s in shares):
        raise ValueError(
            f"layer {name!r}: every split share must be positive, got {shares}."
        )
    total = sum(shares)
    if abs(total - 100.0) < 1e-6:
        raise ValueError(
            f"layer {name!r}: split shares must sum to 1, not 100 -- these look "
            f"like percentages. Write {[s / 100 for s in shares]} instead of "
            f"{shares}."
        )
    if abs(total - 1.0) > 1e-6:
        raise ValueError(
            f"layer {name!r}: split shares must sum to 1, got {total:g} from "
            f"{shares}."
        )
    return shares


def _split_names(names, shares: list[float], unit: str) -> list[str] | None:
    """Validate an explicit ``names=`` list for a split unit, else ``None``.

    ``None`` keeps the generated ``<unit>_1 .. <unit>_N``. An explicit list has
    to name EVERY slice: a partial list would leave the rest generated, so which
    MODFLOW layer a given name refers to would depend on where the caller
    stopped counting. Validated at declaration rather than at build, so a bad
    list is reported at the ``add()`` that wrote it.
    """

    if names is None:
        return None
    # `str` before the coercion: it is iterable, so a bare name would otherwise
    # be silently accepted as a list of one-character names.
    if isinstance(names, str):
        raise TypeError(
            f"layer {unit!r}: names= takes a sequence of names, one per split "
            f"layer; got the string {names!r}."
        )
    try:
        given = list(names)
    except TypeError:
        raise TypeError(
            f"layer {unit!r}: names= takes a sequence of names, one per split "
            f"layer; got {type(names).__name__}."
        ) from None
    if len(shares) == 1:
        raise ValueError(
            f"layer {unit!r}: names= needs a split of 2 or more to name -- an "
            f"unsplit unit is already called {unit!r}. Pass split= as well, or "
            f"rename the unit itself."
        )
    if len(given) != len(shares):
        raise ValueError(
            f"layer {unit!r}: names= must have one name per split layer, got "
            f"{len(given)} name(s) for {len(shares)} layers."
        )
    for candidate in given:
        if not isinstance(candidate, str) or not candidate.strip():
            raise ValueError(
                f"layer {unit!r}: every split name must be a non-empty string, "
                f"got {given!r}."
            )
        if candidate == "top":
            raise ValueError(
                f"layer {unit!r}: 'top' is reserved for the model top surface; "
                f"name the layer for the material it is."
            )
    repeats = sorted({n for n in given if given.count(n) > 1})
    if repeats:
        raise ValueError(
            f"layer {unit!r}: split names must be distinct -- a name identifies "
            f"a contact, a DataFrame column and a 3-D scene actor -- got "
            f"repeats {repeats}."
        )
    return given


def _sub_surface(surface: Surface, share: float, step: float, unit: str) -> Surface:
    """One sub-layer's bottom, for a unit being split.

    Two cases, because a unit is declared in one of two ways. A unit given a
    ``thickness=`` is RELATIVE, and splitting it just divides that thickness --
    no interpolation, nothing to interpolate toward. A unit given an absolute
    ``bottom=`` needs :meth:`~myflopy.surfaces.Surface.toward`, whose ``step``
    is the share of the REMAINING interval (see :func:`_split_steps`).
    """

    if surface.is_relative:
        return _dc_replace(surface, value=float(surface.value) * share)
    if _depends_on_previous(surface):
        raise TypeError(
            f"layer {unit!r} cannot be split: its bottom is a {surface.kind!r} "
            f"surface, which is measured from the layer above, so each cut "
            f"would measure from the previous cut and the contacts would drift "
            f"downward rather than divide the unit. Give the unit an absolute "
            f"bottom (a raster, contours, or points), or declare it with "
            f"thickness= and split that."
        )
    return Surface.toward(surface, step)


def _unit_thickness(thickness: np.ndarray, units: dict[str, list[int]]) -> np.ndarray:
    """Each layer's thickness replaced by its UNIT's total.

    Pinch policy is a statement about a geologic unit -- "where this clay is
    thinner than a foot, take it out" -- and a split is a numerical choice that
    should not change the answer. Evaluated per slice it does: a 2.5 ft unit
    split three ways comes back ``idomain [0, 1, 0]``, holes inside a unit that
    is fully present, and the zeros BLOCK vertical flow. Feeding the unit total
    in makes the verdict all-or-nothing.

    Note dividing ``min_thickness`` by N instead does NOT work -- reconcile has
    already crushed the outer slices by then, so the verdict is the same ragged
    one, and it breaks the ``min_sep < min_thickness`` invariant past N ~ 9.
    """

    if all(len(idx) == 1 for idx in units.values()):
        return thickness
    totals = np.array(thickness, dtype=float)
    for indices in units.values():
        if len(indices) > 1:
            totals[indices] = thickness[indices].sum(axis=0)
    return totals


def _split_steps(shares: list[float]) -> list[float]:
    """Turn shares of the WHOLE unit into shares of what is left at each cut.

    ``Surface.toward`` resolves against the contact immediately above it, not
    the unit top -- so feeding it cumulative fractions 1/3, 2/3, 1 yields
    thicknesses [20, 26.67, 13.33], not equal thirds. The conversion is
    ``g_i = share_i / (1 - sum(shares before i))``. The last step is forced to
    exactly 1.0: it must land ON the unit's base, and floating-point residue in
    the running remainder would otherwise leave it fractionally short.
    """

    steps, remaining = [], 1.0
    for share in shares:
        steps.append(share / remaining)
        remaining -= share
    steps[-1] = 1.0
    return steps


def _reconcile_args(reconcile) -> tuple[bool, str]:
    """Map the facade ``reconcile`` argument to ``LayerSurfaces.sample`` kwargs."""

    if reconcile in (False, None):
        return False, "bottom"
    if reconcile is True:
        return True, "bottom"
    if reconcile in ("bottom", "top"):
        return True, reconcile
    raise ValueError("reconcile must be 'bottom', 'top', or False.")


def modflow_surfaces(source, *, resample: bool = True):
    """Read an existing MODFLOW model's surfaces as :class:`Surface` objects.

    ``source`` is a flopy model (anything exposing ``.modelgrid``) or a flopy
    modelgrid directly. Returns ``(top_surface, [bottom_surface, ...])``. With
    ``resample=True`` (default) each surface is interpolated from the source cell
    centres, so it transfers onto a *different* grid; with ``resample=False`` the
    arrays are used verbatim (the target grid must match cell-for-cell).
    """
    mg = getattr(source, "modelgrid", source)
    top = np.asarray(mg.top, dtype=float).ravel()
    botm = np.asarray(mg.botm, dtype=float)
    botm = botm.reshape(botm.shape[0], -1)
    if resample:
        xc = np.asarray(mg.xcellcenters, dtype=float).ravel()
        yc = np.asarray(mg.ycellcenters, dtype=float).ravel()
        top_s = Surface.from_points(xc, yc, top)
        botm_s = [Surface.from_points(xc, yc, botm[k]) for k in range(botm.shape[0])]
    else:
        top_s = Surface.from_array(top)
        botm_s = [Surface.from_array(botm[k]) for k in range(botm.shape[0])]
    return top_s, botm_s


@dataclass
class _Layer:
    name: str
    surface: Surface
    min_thickness: float | None = None
    pinch: str | None = None
    #: How many model layers this geologic UNIT becomes. `None`/1 = itself; an
    #: int = that many equal shares; a sequence = those shares. The unit stays
    #: ONE record until `_expanded_layers` compiles it, so `min_thickness` and
    #: `pinch` keep meaning "this unit", not "each slice of it".
    split: Any = None
    #: Explicit names for a split unit's layers, one per share. `None` keeps the
    #: generated `<name>_1 .. <name>_N`. Validated where it is declared, so a
    #: stack never reaches `_expanded_layers` with a names/split length mismatch.
    names: Any = None


@dataclass
class LayerQCReport:
    """Diagnostics for a built layer stack -- the problems to fix before MF6 runs.

    The structured result of QC-ing a :class:`LayerBuildResult` (via
    ``result.qc()`` or :meth:`LayerStack.qc`). It counts the geometry pathologies
    that make a DISV grid fail or behave oddly -- cells with no source coverage
    (NaN top/bottom), active cells bounded by a NaN surface or with non-positive
    thickness (fatal), overly thin cells, pinched-out cells, isolated active cells
    with no connection, and how many disconnected active components exist. When
    produced by ``LayerStack.qc`` it also reports how much top-down reconciliation
    moved each surface. Check :attr:`ok` for a fatal/clean verdict, or ``str(report)``
    for a human-readable per-layer breakdown.

    The integer/list fields hold counts (per layer where noted); see the inline
    field comments. ``isolated_active`` lists ``(layer, cell)`` pairs.
    """

    nlay: int
    ncpl: int
    names: list[str]
    nan_top: int                       # cells with no top elevation (no coverage)
    nan_botm: list[int]                # per layer: cells with no bottom elevation
    nan_active_cells: int              # active cells bounded by a NaN surface (fatal)
    nonpositive_active: list[int]      # per layer: active cells with thickness <= 0
    thin: list[int]                    # per layer: cells thinner than min_thickness
    pinched: list[int]                 # per layer: cells removed (idomain != 1)
    isolated_active: list[tuple]       # (layer, cell) active cells with no connection
    n_active_components: int           # connected components of active cells
    component_sizes: list[int]         # sizes, largest first
    # reconcile diagnostics -- only filled by LayerStack.qc (needs the surfaces)
    reconcile_adjusted: list[int] | None = None    # per layer: cells reconcile moved
    reconcile_max_shift: list[float] | None = None  # per layer: largest move

    @property
    def ok(self) -> bool:
        """True when nothing fatal was found (no NaN-bounded or isolated active cells,
        no non-positive active thickness). Thin/pinched cells are expected, not fatal."""
        return (
            self.nan_active_cells == 0
            and sum(self.nonpositive_active) == 0
            and len(self.isolated_active) == 0
        )

    def __str__(self) -> str:
        """Render a human-readable multi-line QC summary (status, warnings, per-layer stats)."""

        head = "OK" if self.ok else "PROBLEMS FOUND"
        lines = [
            f"LayerStack QC [{head}]: {self.nlay} layers, {self.ncpl} cells, "
            f"{self.n_active_components} active component(s)"
        ]
        if self.nan_active_cells:
            lines.append(f"  ** {self.nan_active_cells} active cell(s) bounded by a NaN "
                         "surface (no source coverage) -- fix the source or fill='propagate'")
        if self.isolated_active:
            lines.append(f"  ** {len(self.isolated_active)} isolated active cell(s) with no "
                         "connection -- prune_isolated() deactivates them")
        if sum(self.nonpositive_active):
            lines.append(f"  ** {sum(self.nonpositive_active)} active cell(s) with thickness <= 0")
        if self.n_active_components > 1:
            shown = ", ".join(str(s) for s in self.component_sizes[:5])
            lines.append(f"  note: active cells split into {self.n_active_components} components "
                         f"(sizes: {shown}{'...' if self.n_active_components > 5 else ''})")
        for i, name in enumerate(self.names):
            extra = ""
            if self.reconcile_adjusted is not None:
                extra = (f"  reconcile_moved={self.reconcile_adjusted[i]} "
                         f"(max {self.reconcile_max_shift[i]:.2f})")
            lines.append(
                f"  [{i}] {name:<14} nan_bottom={self.nan_botm[i]} "
                f"thin={self.thin[i]} pinched={self.pinched[i]} "
                f"thickness<=0(active)={self.nonpositive_active[i]}{extra}"
            )
        return "\n".join(lines)


@dataclass
class LayerBuildResult:
    """The disv-ready arrays produced by :meth:`LayerStack.build`.

    Bundles the discretization arrays a layer stack resolves to -- ``top``
    ``(ncpl,)``, ``botm`` ``(nlay, ncpl)``, ``idomain`` ``(nlay, ncpl)`` (1 active,
    -1 pass-through, 0 inactive), and the derived ``thickness`` -- alongside the
    per-layer metadata (``names``, ``min_thickness``, ``pinch`` policy, units) and a
    back-reference to the ``vor`` grid. Feed ``.top``/``.botm``/``.idomain`` straight
    into ``mf.disv(...)``. Inspect quality with :meth:`report` (a per-layer
    thickness/pinch summary) or ``.qc()`` (a :class:`LayerQCReport`), and publish the
    surfaces onto the grid for the GIS-aware builders with :meth:`attach_to_grid`
    (or ``LayerStack.build(attach=True)``).

    The array shapes and idomain encoding are noted in the inline field comments.
    """

    top: np.ndarray          # (ncpl,)
    botm: np.ndarray         # (nlay, ncpl)
    idomain: np.ndarray      # (nlay, ncpl): 1 active, -1 pass-through, 0 inactive
    thickness: np.ndarray    # (nlay, ncpl)
    names: list[str]
    min_thickness: list[float]
    pinch: list[str]
    length_units: str
    time_units: str
    vor: Any = None          # grid the result was built on (for views)
    #: Declared unit name -> the model layers it became, e.g.
    #: ``{"fill": [0], "sand": [1, 2, 3], "clay": [4]}``. Every unit appears,
    #: split or not, so it is also the map from what you DECLARED to what MODFLOW
    #: got. See :meth:`per_layer`.
    units: dict[str, list[int]] = field(default_factory=dict)

    def per_layer(self, mapping, *, default=_UNSET) -> list:
        """Expand a ``{unit: value}`` mapping to one value per MODEL layer.

        The bridge between geology and discretization. Splitting a unit changes
        ``nlay``, and every per-layer argument downstream -- ``mf.npf(k=...)``,
        ``sto``, ``ic`` -- is a POSITIONAL list, so a hand-written one silently
        means something new the moment a split is added or removed::

            mf.npf(k=layers.per_layer({"fill": 30.0, "sand": 25.0, "clay": 0.05}))

        A list of the wrong length is at least rejected by FloPy; a stale list of
        the RIGHT length is accepted with new meaning, which is the failure this
        exists to prevent. Layer names work as keys too, so one slice of a unit
        can differ from its siblings -- a later key wins.

        Parameters
        ----------
        mapping : dict
            Unit names (or individual layer names) to values.
        default : optional
            Value for units the mapping does not name. Without it, an unnamed
            unit raises.
        """

        by_layer: list = [_UNSET] * self.nlay
        for unit, indices in self.units.items():
            if unit in mapping:
                for i in indices:
                    by_layer[i] = mapping[unit]
        for i, name in enumerate(self.names):  # a layer name overrides its unit
            if name in mapping:
                by_layer[i] = mapping[name]

        unknown = set(mapping) - set(self.units) - set(self.names)
        if unknown:
            raise KeyError(
                f"{sorted(unknown)} name neither a unit nor a layer; have units "
                f"{sorted(self.units)} and layers {self.names}."
            )
        missing = [self.names[i] for i, v in enumerate(by_layer) if v is _UNSET]
        if missing:
            if default is _UNSET:
                raise KeyError(
                    f"no value for {missing}; give one per unit, or pass "
                    f"default= to fill the rest."
                )
            by_layer = [default if v is _UNSET else v for v in by_layer]
        return by_layer

    @property
    def nlay(self) -> int:
        """Number of layers (rows of ``botm``)."""

        return self.botm.shape[0]

    @property
    def n_pinched(self) -> int:
        """Cells removed from the solution (idomain != 1)."""
        return int((self.idomain != 1).sum())

    def report(self) -> str:
        """Per-layer thickness + thin/pinched-cell summary."""
        n_units = len(self.units) or self.nlay
        split = "" if n_units == self.nlay else f" from {n_units} units"
        lines = [
            f"LayerStack: {self.nlay} layers{split}, {self.top.size} cells "
            f"[{self.length_units}/{self.time_units}], {self.n_pinched} pinched cells"
        ]
        # `thin` counts against the UNIT total, matching the pinch verdict: a
        # slice of a split unit is thinner than its unit by construction, and
        # reporting each slice as thin would flag every split stack.
        unit_thk = _unit_thickness(self.thickness, self.units)
        for i, name in enumerate(self.names):
            t = self.thickness[i]
            thin = int((unit_thk[i] < float(self.min_thickness[i])).sum())
            pinched = int((self.idomain[i] != 1).sum())
            lines.append(
                f"  [{i}] {name:<14} thk min={t.min():.2f} mean={t.mean():.2f} "
                f"max={t.max():.2f}  thin(<{self.min_thickness[i]})={thin} "
                f"pinched={pinched} ({self.pinch[i]})"
            )
        return "\n".join(lines)

    def to_disv(self, vor=None, *, name: str = "disv"):
        """Return the ``mf.disv`` spec for **this** build -- no second resolution.

        ``LayerStack.to_disv`` samples the stack again, so it is a different
        resolution from the one you are holding, and any argument the two calls do
        not share writes a model that is not the one you inspected. That drifted
        8.5 ft on a real stack, and because :meth:`attach_to_grid` publishes *this*
        build's surfaces, ``mf.CellSurfaceOffset("cell_bottom", ...)`` then placed
        drains against bottoms MODFLOW never received.

        Building the spec from the arrays already on this result removes the
        possibility by construction: there is one geometry, and it is the one you
        checked with :meth:`qc`.

        Parameters
        ----------
        vor : VoronoiGridPlus, optional
            The grid to take cell geometry from; defaults to the one built against.
        name : str, default "disv"
            Package name.

        Returns
        -------
        PackageSpec

        Examples
        --------
        >>> res = layer_stack.build(vor)
        >>> res.attach_to_grid()
        >>> disv = res.to_disv(vor)
        """


        grid = self.vor if vor is None else vor
        if grid is None:
            raise ValueError("No grid; pass vor= or build the stack with a grid.")
        props = grid.get_disv_gridprops()
        options = {}
        if self.length_units is not None:
            options["length_units"] = _MF6_LENGTH_UNITS.get(
                self.length_units.lower(), self.length_units.upper()
            )
        return disv_spec(
            nlay=self.nlay,
            ncpl=props["ncpl"],
            nvert=len(props["vertices"]),
            vertices=props["vertices"],
            cell2d=props["cell2d"],
            top=np.asarray(self.top, dtype=float),
            botm=np.asarray(self.botm, dtype=float),
            idomain=np.asarray(self.idomain),
            name=name,
            **options,
        )

    def attach_to_grid(self, vor=None):
        """Publish ``top``/``botm`` onto ``vor.gdf_topbtm`` for grid-aware builders.

        Writes a centroid GeoDataFrame with integer columns (``0`` = model top,
        ``1..nlay`` = layer bottoms) -- the format the surface-aware builders read
        (SFR reach tops, LAK lake-cell layering). Returns the grid.

        Usually you do not call this directly: pass ``attach=True`` to
        :meth:`LayerStack.build`.
        """
        import geopandas as gpd

        grid = self.vor if vor is None else vor
        if grid is None:
            raise ValueError("No grid to attach to; pass vor= or build the stack with a grid.")
        columns = {0: np.asarray(self.top, dtype=float)}
        for i in range(self.nlay):
            columns[i + 1] = np.asarray(self.botm[i], dtype=float)
        grid.gdf_topbtm = gpd.GeoDataFrame(
            {"geometry": grid.gdf_vorPolys.geometry, **columns},
            geometry="geometry", crs=grid.crs,
        )
        return grid

    # -- QC / validation -------------------------------------------------- #
    def _active_components(self):
        """Union-find over active cells (idomain == 1). Returns ``(labels, sizes)``:
        ``labels`` is an ``(nlay, ncpl)`` array of component roots (-1 where not
        active) and ``sizes`` maps a root to its cell count.

        Connections follow MF6: horizontal between active plan-neighbors, and
        vertical down a column where ``idomain == -1`` (pass-through) bridges
        active cells while ``idomain == 0`` (inactive) blocks them."""
        from collections import Counter

        nlay, ncpl = self.idomain.shape
        parent = list(range(nlay * ncpl))

        def find(a):
            """Union-find root of flat cell index ``a``, with path compression."""

            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        def union(a, b):
            """Merge the union-find sets containing flat cell indices ``a`` and ``b``."""

            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb

        active = self.idomain == 1
        adj = self.vor.adjacent_cells_idx
        for k in range(nlay):
            base = k * ncpl
            for i in range(ncpl):
                if not active[k, i]:
                    continue
                for j in adj[i]:
                    if 0 <= j < ncpl and active[k, j]:
                        union(base + i, k * ncpl + j)
        for i in range(ncpl):
            last = None
            for k in range(nlay):
                d = self.idomain[k, i]
                if d == 1:
                    if last is not None:
                        union(last * ncpl + i, k * ncpl + i)
                    last = k
                elif d == 0:
                    last = None      # inactive blocks vertical flow
                # d == -1: pass-through keeps ``last`` reachable
        labels = np.full((nlay, ncpl), -1, dtype=int)
        sizes = Counter()
        for k in range(nlay):
            for i in range(ncpl):
                if active[k, i]:
                    r = find(k * ncpl + i)
                    labels[k, i] = r
                    sizes[r] += 1
        return labels, sizes

    def qc(self) -> LayerQCReport:
        """Check the built stack for problems MODFLOW 6 would choke on.

        Catches **NaN-bounded active cells** (a surface had no source coverage),
        **non-positive active thickness**, and **isolated active cells** (no
        connection to any neighbour) -- plus thin/pinched counts and the number
        of connected active components. Returns a :class:`LayerQCReport` whose
        ``.ok`` is False if anything fatal was found. ``LayerStack.qc`` adds
        reconcile diagnostics on top of this."""
        nlay, ncpl = self.idomain.shape
        active = self.idomain == 1
        surf_top = np.vstack([self.top[None, :], self.botm[:-1]])  # top of each layer
        nan_bounded = np.isnan(surf_top) | np.isnan(self.botm)
        labels, sizes = self._active_components()
        isolated = [
            (k, i)
            for k in range(nlay)
            for i in range(ncpl)
            if active[k, i] and sizes[labels[k, i]] == 1
        ]
        return LayerQCReport(
            nlay=nlay,
            ncpl=ncpl,
            names=list(self.names),
            nan_top=int(np.isnan(self.top).sum()),
            nan_botm=[int(np.isnan(self.botm[k]).sum()) for k in range(nlay)],
            nan_active_cells=int((nan_bounded & active).sum()),
            nonpositive_active=[int(((self.thickness[k] <= 0) & active[k]).sum())
                                for k in range(nlay)],
            # Against the UNIT total, matching the pinch verdict: a slice of a
            # split unit is thinner than its unit by construction, so counting
            # per slice would flag every split stack as thin.
            thin=[int((_unit_thickness(self.thickness, self.units)[k]
                       < float(self.min_thickness[k])).sum())
                  for k in range(nlay)],
            pinched=[int((self.idomain[k] != 1).sum()) for k in range(nlay)],
            isolated_active=isolated,
            n_active_components=len(sizes),
            component_sizes=sorted(sizes.values(), reverse=True),
        )

    def validate(self) -> LayerBuildResult:
        """Raise ``ValueError`` if :meth:`qc` finds anything fatal; else return self."""
        report = self.qc()
        if not report.ok:
            raise ValueError("Layer stack failed QC:\n" + str(report))
        return self

    def prune_isolated(self) -> LayerBuildResult:
        """Deactivate (``idomain -> 0``) active cells that have no connection.

        Returns a new result. Isolated cells have no edges, so removing them
        never disconnects anything else -- a single pass is enough."""
        isolated = self.qc().isolated_active
        if not isolated:
            return self
        idomain = self.idomain.copy()
        for k, i in isolated:
            idomain[k, i] = 0
        return _dc_replace(self, idomain=idomain)

    # -- views (thin wrappers over flopy / the surface API) --------------- #
    def _resolve_line(self, line, x, y) -> dict:
        """Build a flopy cross-section line from an explicit line, x=, y=, or a
        default West-East line through the grid centre.

        ``line`` is anything :func:`_section_line_points` accepts -- geometry,
        vector file, or points."""
        if line is not None:
            return {"line": _section_line_points(line)}
        minx, miny, maxx, maxy = (float(v) for v in self.vor.gdf_vorPolys.total_bounds)
        if x is not None:
            return {"line": [(float(x), miny), (float(x), maxy)]}
        if y is not None:
            return {"line": [(minx, float(y)), (maxx, float(y))]}
        ymid = (miny + maxy) / 2.0  # default: West-East section through the centre
        return {"line": [(minx, ymid), (maxx, ymid)]}

    def vertex_grid(self):
        """A flopy ``VertexGrid`` carrying this result's geometry."""
        import flopy

        p = self.vor.get_disv_gridprops()
        return flopy.discretization.VertexGrid(
            vertices=p["vertices"], cell2d=p["cell2d"],
            top=self.top, botm=self.botm, idomain=self.idomain,
            nlay=self.nlay, ncpl=p["ncpl"], crs=str(getattr(self.vor, "crs", None)),
        )

    def _require_finite_geometry(self, verb: str) -> None:
        """Refuse to draw a stack whose elevations contain ``NaN``.

        A NaN top or bottom propagates into the renderer's axis limits, where
        Matplotlib rejects it as ``Axis limits cannot be NaN or Inf`` -- eight
        frames deep, naming neither the layer nor the cause. The cause is always
        the same: a source surface had no data over some cells. :meth:`qc`
        already counts them per layer, so point at it.
        """

        nan_top = int(np.isnan(self.top).sum())
        nan_botm = [int(np.isnan(self.botm[k]).sum()) for k in range(self.nlay)]
        if not nan_top and not any(nan_botm):
            return
        worst = [
            f"{self.names[k]!r} ({n} cell{'s' if n != 1 else ''})"
            for k, n in enumerate(nan_botm) if n
        ]
        detail = ""
        if nan_top:
            detail += f"the model top ({nan_top} cells)"
        if worst:
            detail += (", " if detail else "") + ", ".join(worst)
        raise ValueError(
            f"cannot draw {verb}: the stack has NaN elevations in {detail}. A "
            f"surface had no source data over those cells -- a raster that does "
            f"not cover the whole grid, contours interpolated inside a smaller "
            f"hull, or a nodata value read as elevation. Run .qc() for the "
            f"per-layer counts, then either extend the source, pass "
            f"fill='propagate' on that layer to inherit the surface above, or "
            f"deactivate the cells with idomain."
        )

    def _draw_cross_section(
        self, line=None, *, x=None, y=None, color_by="layer", ax=None,
        cmap="tab10", show_grid=True, legend=True, title=None,
    ):
        """Render the layer cross-section into ``ax`` and return it.

        Private: this is the RENDERER. ``stack.plot.section(y=300)`` is the verb,
        and it returns a Picture that answers ``.show()``/``.save()``/``.html()``
        like every other. Matplotlib-native by necessity -- see
        :class:`LayerSection`.

        Layer coloring delegates to the shared
        :func:`~myflopy.modflow.mf6.cross_section_plotting.plot_layered_cross_section`
        renderer (the same core behind the model-aware ``plot_model_cross_section``)."""
        import matplotlib.pyplot as plt

        from myflopy.modflow.mf6.cross_section_plotting import (
            ModelCrossSectionStyle,
            plot_layered_cross_section,
        )

        self._require_finite_geometry("a cross-section")
        vg = self.vertex_grid()
        line_spec = self._resolve_line(line, x, y)

        if color_by == "thickness":
            import flopy

            if ax is None:
                _, ax = mpl_axes(figsize=(9, 4))
            xsec = flopy.plot.PlotCrossSection(modelgrid=vg, line=line_spec, ax=ax)
            xsec.plot_array(self.thickness, cmap="viridis")
            if show_grid:
                xsec.plot_grid(lw=0.25, color="0.3")
            ax.set_title(title or "Layer cross-section")
            ax.set_xlabel(f"distance along section [{self.length_units}]")
            ax.set_ylabel(f"elevation [{self.length_units}]")
            return ax

        cm = plt.get_cmap(cmap)
        colors = [cm(i % 10) for i in range(self.nlay)]  # matplotlib RGBA tuples
        style = ModelCrossSectionStyle(
            figsize=(9, 4), grid_linewidth=0.25, grid_color="0.3",
            layer_alpha=0.75, title_fontsize=12, label_fontsize=10,
            legend_loc="upper right", legend_frameon=True, legend_fontsize=8,
        )
        _, ax = plot_layered_cross_section(
            vg, line_spec, ax=ax, style=style,
            layer_colors=colors, layer_labels=self.names,
            show_grid=show_grid, show_layers=True, show_head=False,
            show_legend=legend, title=title or "Layer cross-section",
            xlabel=f"distance along section [{self.length_units}]",
            ylabel=f"elevation [{self.length_units}]",
        )
        return ax

    @property
    def surface_names(self) -> list[str]:
        """Plottable surface names: the model top plus each layer's bottom contact."""
        return ["top"] + list(self.names)

    def _surface_z(self, name) -> np.ndarray:
        """Elevation array for a surface name (``"top"`` or a layer's bottom)."""
        if name == "top":
            return self.top
        return self.botm[self.names.index(name)]

    def _resolve_surface_names(self, layer) -> list[str]:
        """Normalize the ``layer`` selector to a list of valid surface names."""
        valid = self.surface_names
        if isinstance(layer, str) and layer == "all":
            return valid
        names = list(layer) if isinstance(layer, (list, tuple)) else [layer]
        bad = [n for n in names if n not in valid]
        if bad:
            raise KeyError(f"no surface named {bad!r}; choose from {valid} or 'all'.")
        # Deduplicate, order preserved: drawing one contact twice means nothing,
        # and in the VTK scene each name is an ACTOR name, where a repeat
        # silently replaces the first rather than adding a second.
        return list(dict.fromkeys(names))

    @staticmethod
    def _rgb(rgba) -> str:
        """Format a 0-1 RGBA tuple as a Plotly ``rgb(r,g,b)`` string (0-255)."""

        r, g, b = (int(round(255 * c)) for c in rgba[:3])
        return f"rgb({r},{g},{b})"

    def _surface_fig(
        self, layer="top", *, resolution=120, colorscale="Earth_r",
        color_by=None, opacity=None, height=None,
    ):
        """The 3-D surface figure. Reached as ``stack.plot.surface(...)``.

        Private because it is Layer 1 only: it BUILDS the figure and nothing
        else. Writing it out and opening it are `.html(path)` / `.show()` on the
        Picture that wraps it -- this used to take `html_path=`/`browser=` and do
        both itself.

        Interactive 3D surface(s) of one or more layers (plotly).

        ``layer`` selects which surface(s) to draw and may be:

        * ``"top"`` or any layer name -- that layer's *bottom* contact (default ``"top"``),
        * a list of names, e.g. ``["top", "clay", "bedrock"]``, overlaid in one scene,
        * ``"all"`` -- the model top plus every layer bottom.

        Discover the choices with ``result.surface_names``. A single surface with
        relief is shaded by elevation (``colorscale``); a flat surface or several
        surfaces get distinct solid colors. Override with
        ``color_by="elevation"`` / ``"surface"`` and tune ``opacity``.

        ``height`` sets the figure height in pixels (default ``None`` = fill the
        container / browser window)."""

        names = self._resolve_surface_names(layer)
        crs = str(getattr(self.vor, "crs", None))
        xs = np.asarray(self.vor.centroids[0])
        ys = np.asarray(self.vor.centroids[1])

        # Each surface is built by the shared InterpolatedSurface.surface_trace
        # so there is one go.Surface builder across the toolkit; we keep the
        # interpolated z-array to set a common elevation scale below.
        grids = []
        for name in names:
            isurf = InterpolatedSurface(
                xs=xs, ys=ys, zs=np.asarray(self._surface_z(name)),
                surf_type="lyr", resolution=resolution, crs=crs,
            )
            grids.append((name, isurf, isurf.surface))

        # Overall elevation span -> colour mode, shared colour range and z-axis.
        zmin = min(float(np.nanmin(zz)) for *_, zz in grids)
        zmax = max(float(np.nanmax(zz)) for *_, zz in grids)
        flat = zmax - zmin < 1e-6
        if color_by is None:
            color_by = "elevation" if (len(names) == 1 and not flat) else "surface"
        if opacity is None:
            opacity = 1.0 if len(names) == 1 else 0.85
        zaxis = dict(title=f"elev [{self.length_units}]")
        aspectmode = "manual"
        if flat:
            # A perfectly flat surface has zero vertical extent; a lone flat
            # go.Surface then fails to render in hardware-WebGL viewers (its faces
            # get a 0/0 colour and the trace has no z-depth), while multi-surface
            # scenes dodge this because their combined z-range is non-zero. Give
            # the sheet a faint, non-planar relief and a real z-axis so it always
            # draws; ``cube`` aspect keeps the near-flat sheet prominent.
            mid = 0.5 * (zmin + zmax)
            pad = max(1.0, abs(mid) * 0.01)
            gx0, gy0 = grids[0][1].xy_meshgrid
            xspan = float(np.nanmax(gx0) - np.nanmin(gx0)) or 1.0
            yspan = float(np.nanmax(gy0) - np.nanmin(gy0)) or 1.0
            ripple = (
                ((gx0 - np.nanmin(gx0)) / xspan - 0.5)
                + ((gy0 - np.nanmin(gy0)) / yspan - 0.5)
            ) * (0.5 * pad)                      # ~+/-pad/2 of imperceptible relief
            grids = [(n, s, zz + ripple) for (n, s, zz) in grids]
            zmin, zmax = mid - pad, mid + pad
            zaxis["range"] = [zmin, zmax]
            aspectmode = "cube"

        fig = Fig()
        if color_by == "elevation":
            for i, (name, isurf, zz) in enumerate(grids):
                fig.add_trace(isurf.surface_trace(
                    surface=zz, colorscale=colorscale, name=name,
                    cmin=zmin, cmax=zmax, opacity=opacity, showscale=(i == 0),
                    colorbar=dict(title=f"elev [{self.length_units}]"),
                ))
        else:  # one solid color per surface, distinguished by a legend
            import matplotlib.pyplot as plt

            cmap = plt.get_cmap("tab10")
            for i, (name, isurf, zz) in enumerate(grids):
                c = self._rgb(cmap(i % 10))
                fig.add_trace(isurf.surface_trace(
                    surface=zz, colorscale=[[0, c], [1, c]], name=name,
                    cmin=zmin, cmax=zmax, opacity=opacity,
                    showscale=False, showlegend=True,
                ))

        multi = len(names) > 1
        scene = dict(aspectmode=aspectmode, zaxis=zaxis)
        if aspectmode == "manual":
            scene["aspectratio"] = dict(x=1, y=0.6, z=0.45)
        fig.update_layout(
            title=("layer surfaces" if multi else f"{names[0]} surface"),
            autosize=True, height=height, margin=dict(l=0, r=0, t=40, b=0),
            showlegend=(color_by == "surface" and multi), scene=scene,
        )
        return fig

    def _resolve_layer_indices(self, layers) -> list[int]:
        """Normalize a layer selector to sorted, unique layer indices.

        Accepts ``None``/``"all"`` (every layer), a single layer name or integer
        index, a UNIT name (which expands to every layer it was split into), or a
        list mixing them."""
        if layers is None or (isinstance(layers, str) and layers == "all"):
            return list(range(self.nlay))
        items = list(layers) if isinstance(layers, (list, tuple)) else [layers]
        idx = []
        for it in items:
            if isinstance(it, (int, np.integer)) and not isinstance(it, bool):
                i = int(it)
                if not 0 <= i < self.nlay:
                    raise IndexError(f"layer index {i} out of range [0, {self.nlay}).")
                idx.append(i)
            elif it in self.names:
                idx.append(self.names.index(it))
            elif it in self.units:
                # A unit name selects every layer it became, so `grid("sand")`
                # still works once `sand` is split into `sand_1..sand_3`.
                idx.extend(self.units[it])
            else:
                choices = sorted(set(self.names) | set(self.units))
                raise KeyError(
                    f"no layer or unit named {it!r}; choose from {choices} or an index."
                )
        return sorted(set(idx))

    def _vtk_surface_plotter(
        self, layer="all", *, resolution=120, scale=8, cmap="tab10",
        opacity=1.0, width=900, height=580, show_edges=False,
    ):
        """Build the PyVista scene for one or more CONTACT SURFACES.

        Private, Layer 1 only: it builds the scene, and the
        :class:`~myflopy.viz.VtkScene` that wraps it does the showing and
        writing. Reached as ``stack.plot.surface(..., backend="vtk")``.

        The sibling of :meth:`_vtk_plotter`, and a different subject: that one
        renders the layered cell VOLUME, this renders the contacts as separate
        sheets. Sheets are what you want when the question is "where does this
        contact go" -- each is its own actor, so a viewer can hide the ones above
        to look underneath, which a single fused volume will not let you do.

        Each surface is interpolated by the shared
        :class:`~myflopy.modflow.mf6.grid.interpolated_surface.InterpolatedSurface`,
        so a VTK sheet and its Plotly counterpart are the same numbers.

        Returns the plotter and the sheet meshes, which the wrapping
        :class:`~myflopy.viz.VtkScene` carries as ``.meshes`` -- the handle for
        writing one out (``.save("contact.vtu")``) and opening it elsewhere.
        """

        pv = require("pyvista", feature="interactive 3-D scenes")
        import matplotlib as mpl

        names = self._resolve_surface_names(layer)
        xs = np.asarray(self.vor.centroids[0])
        ys = np.asarray(self.vor.centroids[1])
        colors = mpl.colormaps[cmap]

        # `off_screen`, matching `_vtk_plotter`: without it there is no render
        # window outside Jupyter, and `.save("sheets.png")` fails with "Nothing
        # to screenshot" in a plain script. A notebook forces it on anyway.
        plotter = pv.Plotter(off_screen=True, window_size=(width, height))
        sheets = []
        for i, name in enumerate(names):
            interp = InterpolatedSurface(
                xs=xs, ys=ys, zs=np.asarray(self._surface_z(name)),
                surf_type="lyr", resolution=resolution,
                crs=str(getattr(self.vor, "crs", None)),
            )
            gx, gy = interp.xy_meshgrid
            gz = np.asarray(interp.surface, dtype=float)
            grid = pv.StructuredGrid(
                np.asarray(gx, float), np.asarray(gy, float), gz * scale
            )
            # Cells with no data would otherwise render as a sheet pinned to z=0,
            # which reads as a real contact at sea level.
            grid["elevation"] = gz.ravel(order="F")
            sheet = grid.threshold(scalars="elevation")
            plotter.add_mesh(
                sheet,
                # PyVista wants floats or hex, not plotly's `rgb(r,g,b)` string,
                # so `_rgb` (which the Plotly path uses) is deliberately skipped.
                color=tuple(colors(i % colors.N)[:3]),
                opacity=opacity,
                show_edges=show_edges,
                label=str(name),
                # `name=` as well as `label=`: the label is what the LEGEND
                # says, the name is what `plotter.actors` is KEYED by. Without
                # it the keys are addresses -- `UnstructuredGrid(Addr=0x...)` --
                # so "hide the sheets above", which is the whole reason these
                # are separate actors, has no way to say which sheet it means.
                name=str(name),
                smooth_shading=True,
            )
            sheets.append(sheet)
        if len(names) > 1:
            plotter.add_legend(bcolor=None)
        plotter.add_axes()
        return plotter, sheets

    def _vtk_plotter(
        self, layers=None, *, color_by="layer", scale=8, cmap=None,
        width=900, height=580,
    ):
        """Build the PyVista plotter for the layered grid volume.

        Private: this is Layer 1 only -- it BUILDS the scene. Displaying it and
        writing it out belong to the :class:`~myflopy.viz.VtkScene` that wraps it.
        Reached as ``stack.plot.grid(backend="vtk")``.

        ``layers`` selects which layers to show: ``None``/``"all"`` (default) for
        every layer, a single layer name or index, or a list mixing them
        (e.g. ``["sand", "clay"]`` or ``[0, 2]``). Colors stay keyed to each
        layer's position, so a subset keeps the same colors it has in the full
        stack.

        Every cell carries ``layer``, ``thickness``, ``top``, ``botm`` and
        ``cellid`` whichever one ``color_by`` draws, so a picked cell can report
        all of them and a written ``.vtu`` keeps them all. ``cmap=None`` takes
        the default that suits the chosen scalar (:data:`_VTK_GRID_CMAPS`).

        Returns the plotter and, as the one mesh, the ASSEMBLED volume -- not
        the per-layer pieces the actors draw. Those are views of it, and it is
        the whole selection with every array on it, which is what makes
        ``scene.meshes[0].save("stack.vtu")`` the useful thing to hand to
        ParaView."""
        import tempfile
        from pathlib import Path

        import flopy
        from flopy.export.vtk import Vtk

        pv = require("pyvista", feature="interactive 3-D scenes")

        if color_by not in _VTK_GRID_CMAPS:
            raise ValueError(
                f"color_by must be one of {list(_VTK_GRID_CMAPS)}, not "
                f"{color_by!r}."
            )
        sel = self._resolve_layer_indices(layers)
        p = self.vor.get_disv_gridprops()
        ws = Path(tempfile.mkdtemp(prefix="layer_view_"))
        sim = flopy.mf6.MFSimulation(sim_name="layerview", sim_ws=str(ws))
        flopy.mf6.ModflowTdis(sim)
        flopy.mf6.ModflowIms(sim)
        gwf = flopy.mf6.ModflowGwf(sim, modelname="layers")
        flopy.mf6.ModflowGwfdisv(
            gwf, nlay=self.nlay, ncpl=p["ncpl"], nvert=len(p["vertices"]),
            vertices=p["vertices"], cell2d=p["cell2d"],
            top=self.top, botm=self.botm, idomain=self.idomain,
        )
        vtk = Vtk(model=gwf, vertical_exageration=1, binary=True, smooth=False)
        vtk.add_model(gwf)
        # Always attach a per-cell layer index so a subset can be selected.
        # Integer, so `add_array` leaves it alone -- see the float note below.
        vtk.add_array(np.repeat(np.arange(self.nlay)[:, None], p["ncpl"], axis=1), "layer")
        mesh = vtk.to_pyvista()
        if isinstance(mesh, pv.MultiBlock):
            mesh = mesh.combine()

        # Per-cell fields, attached HERE rather than through `Vtk.add_array`,
        # which NaN-masks FLOAT arrays wherever idomain == 0 -- i.e. on exactly
        # the cells a pinched-out layer creates, which are the ones worth
        # inspecting. (Measured: 138 of 843 cells on a stack with one pinching
        # layer; integer arrays come through intact.) `top` and `botm` are
        # already here from `add_model` and are OVERWRITTEN, not added: flopy's
        # own `top` is NaN for every layer below 0, so a per-layer top has to be
        # built from the contacts.
        mesh.cell_data["thickness"] = self.thickness.ravel()
        mesh.cell_data["top"] = np.vstack([self.top[None, :], self.botm[:-1]]).ravel()
        mesh.cell_data["botm"] = self.botm.ravel()
        mesh.cell_data["cellid"] = np.tile(np.arange(p["ncpl"]), self.nlay)

        if len(sel) != self.nlay:  # keep only the requested layers' cells
            layer_cell = np.asarray(mesh.cell_data["layer"]).astype(int)
            mesh = mesh.extract_cells(np.isin(layer_cell, sel))
        if color_by == "elevation":
            # The one POINT array, and the odd one out: vertex z, so it ramps
            # within a cell where `top`/`botm` are flat per-cell contacts.
            mesh["elevation"] = mesh.points[:, 2]

        mesh_kwargs = dict(
            show_edges=True,
            cmap=_VTK_GRID_CMAPS[color_by] if cmap is None else cmap,
        )
        if color_by == "layer":
            # Discrete colors keyed to each layer's global index; label only the
            # layers actually shown.
            mesh_kwargs.update(
                n_colors=self.nlay,
                clim=[-0.5, self.nlay - 0.5],
                annotations={float(i): self.names[i] for i in sel},
                scalar_bar_args=dict(title="layer", n_labels=0),
            )
        else:
            # An explicit clim over the whole selection, because each layer is
            # its own actor below and PyVista would otherwise scale each one to
            # its OWN range -- making a thin layer and a thick one look alike.
            finite = np.asarray(mesh[color_by], dtype=float)
            finite = finite[np.isfinite(finite)]
            mesh_kwargs.update(
                clim=[float(finite.min()), float(finite.max())] if finite.size else None,
                scalar_bar_args=dict(title=color_by),
            )

        plotter = pv.Plotter(off_screen=True, window_size=(width - 40, height - 40))
        # One actor per layer, named for it, rather than one fused volume. The
        # fused mesh drew the same picture, but a viewer could not take it
        # apart: `scene.scene.actors["clay"].visibility = False` needs a `clay`
        # actor to exist. Costs nlay draw calls instead of one.
        #
        # Each actor asks for the colour bar and the scene still gets exactly
        # one: PyVista keys bars by TITLE and reuses the existing one. That is
        # only correct because every actor shares the clim and cmap set above --
        # give them separate ranges and the single bar would describe one layer
        # while colouring all of them.
        layer_cell = np.asarray(mesh.cell_data["layer"]).astype(int)
        for k in sel:
            sheet = mesh.extract_cells(layer_cell == k)
            if sheet.n_cells == 0:  # selected but wholly absent from the mesh
                continue
            plotter.add_mesh(
                sheet, scalars=color_by, name=self.names[k], **mesh_kwargs,
            )
        plotter.set_scale(zscale=scale)
        plotter.add_axes()
        plotter.camera_position = "yz"
        return plotter, [mesh]

    def _thickness_values(self, layer=None):
        """Per-cell thickness (total, or a single named layer) and its label."""

        if layer is None:
            return self.thickness.sum(axis=0), "total thickness"
        return self.thickness[self.names.index(layer)], f"{layer!r} thickness"

    def _draw_thickness_map(self, *, layer=None, ax=None):
        """Render the thickness choropleth into ``ax`` (geopandas, no basemap)."""

        values, title = self._thickness_values(layer)
        gdf = self.vor.gdf_vorPolys.copy().assign(_thickness=values)
        if ax is None:
            _, ax = mpl_axes()
        gdf.plot(column="_thickness", ax=ax, legend=True)
        ax.set_title(f"{title} [{self.length_units}]")
        ax.set_aspect("equal")
        return ax

    @property
    def plot(self) -> StackPlots:
        """The plotting verbs for this layer geometry: ``map``, ``section``,
        ``surface``, ``grid``.

        Same verbs as ``model.plot`` and ``vor.plot``, and everything returned is
        a :class:`~myflopy.viz.Picture`. A stack has no results and no time, so
        there is no ``animate``; ``qc()`` stays a report rather than becoming a
        picture, because it is text you read.

        ``surface`` is a height field -- one contact as ``z(x, y)``. The layered
        cell VOLUME is ``grid(backend="vtk")``, a different shape entirely.

        Replaces ``thickness_map()``/``preview()`` (now ``plot.map()``),
        ``cross_section()`` (``plot.section()``), ``surface_3d()``
        (``plot.surface()``), ``vtk_3d()`` (``plot.grid()``) and ``views()``
        (compose what you want with ``myflopy.plot.mosaic``).

        Examples
        --------
        >>> stack.plot.map()                          # total thickness
        >>> stack.plot.map("sand", basemap=True)      # one layer, on a basemap
        >>> stack.plot.section(y=300)
        >>> stack.plot.surface("all")
        >>> stack.plot.grid(["sand", "clay"], scale=12)
        """

        return StackPlots(self)


class LayerSection(MplPicture):
    """A filled, layer-coloured cross-section through a built stack.

    Matplotlib by necessity, not by preference: the renderer is FloPy's
    ``PlotCrossSection``, and the Plotly ``GridSection`` draws cell outlines
    without layer fills. Rather than exempt it from the picture grammar, it is an
    :class:`~myflopy.viz.MplPicture` -- same verbs, Axes underneath.
    """

    def __init__(self, result, line=None, *, x=None, y=None, **kwargs):
        """Bind a section of ``result`` along ``line`` (or ``x=``/``y=``).

        The line and ``color_by`` are checked HERE, not at draw time. A Picture
        is lazy by design, but an argument that can never work is a caller
        mistake, and reporting it from ``.axes`` puts the traceback in the wrong
        place entirely."""

        self._result = result
        self._line = None if line is None else _section_line_points(line)
        self._x, self._y = x, y
        color_by = kwargs.get("color_by", "layer")
        if color_by not in _SECTION_COLOR_BY:
            raise ValueError(
                f"color_by must be one of {list(_SECTION_COLOR_BY)}, not "
                f"{color_by!r}."
            )
        self._kwargs = kwargs
        self.title = kwargs.get("title") or "Layer cross-section"

    def draw(self, ax=None, **kwargs):
        """Render the section into ``ax`` (or a new one) and return the Axes."""

        return self._result._draw_cross_section(
            self._line, x=self._x, y=self._y, ax=ax, **{**self._kwargs, **kwargs}
        )


class LayerThicknessMap(MplPicture):
    """Per-cell layer thickness, drawn on the grid without a basemap.

    Deliberately NOT a :class:`Choro`: a stack under construction is often on
    synthetic or local coordinates, and a web basemap would put it in the ocean.
    ``stack.plot.map(basemap=True)`` gives the georeferenced choropleth when the
    grid really is where it says it is.
    """

    def __init__(self, result, layer=None):
        """Bind a thickness map of ``result`` (total, or one named layer)."""

        self._result = result
        self._layer = layer
        _, self.title = result._thickness_values(layer)

    def draw(self, ax=None, **kwargs):
        """Render the thickness map into ``ax`` (or a new one)."""

        return self._result._draw_thickness_map(layer=self._layer, ax=ax, **kwargs)


class LayerSurface(Picture):
    """One or more layer contacts as an interactive 3-D Plotly surface."""

    def __init__(self, result, layer="top", **kwargs):
        """Bind a 3-D surface of ``result`` for the chosen contact(s)."""

        self._result = result
        self._layer = layer
        self._kwargs = kwargs
        self._built = None

    @property
    def fig(self) -> Fig:
        """The assembled 3-D figure (built once, then cached)."""

        if self._built is None:
            self._built = self._result._surface_fig(self._layer, **self._kwargs)
        return self._built


#: `_vtk_plotter`'s own defaults, so `grid(backend="plotly")` can tell a value
#: the caller chose from one it merely inherited. Mirrored, not imported, because
#: the signature is the public contract and a test pins the two equal.
_VTK_GRID_DEFAULTS = {
    "color_by": "layer", "scale": 8, "cmap": None, "width": 900, "height": 580,
}


class StackPlots:
    """The plotting verbs for layer geometry -- ``stack.plot.map()`` (plan 8.5a).

    Three verbs, because a layer stack can answer three questions: how thick is
    it (:meth:`map`), what does it look like in section (:meth:`section`), and
    what shape is a given contact (:meth:`surface`). No ``animate`` -- a stack has
    no stress periods. ``qc()`` stays a method on the stack, because a QC report
    is text you read, not a picture you look at.

    Everything returned is a :class:`~myflopy.viz.Picture`, so it renders inline
    and answers ``.show()`` / ``.save(path)`` / ``.html(path)``. Two of the three
    are Matplotlib underneath and answer those over ``.axes``; only
    :meth:`surface` has a Plotly ``.fig``.
    """

    def __init__(self, result):
        """Bind the plotting verbs to a built :class:`LayerBuildResult`."""

        self.result = result

    def __repr__(self):
        """Name the verbs, since tab-completion is how this gets found."""

        return f"StackPlots({self.result.nlay} layers: map, section, surface, grid)"

    def map(self, layer=None, *, basemap: bool = False, **kwargs):
        """Per-cell thickness -- total, or one named layer.

        Draws on the grid with no basemap by default, because a stack is often
        still on synthetic coordinates. ``basemap=True`` routes through the
        shared choropleth instead, for a grid that really is georeferenced.

        Parameters
        ----------
        layer : str or int, optional
            One layer by name or index. With none, total thickness.
        basemap : bool, default False
            Route through the shared choropleth (``vor.plot.map``) so the cells
            sit on a web basemap. Requires a real CRS.
        **kwargs
            Choropleth styling, forwarded to :meth:`~myflopy.modflow.mf6.grid
            .plotting.GridPlots.map`. **Only meaningful with ``basemap=True``**
            -- the default renderer is Matplotlib and takes none of them.

        Raises
        ------
        TypeError
            If styling arguments are given without ``basemap=True``. They used
            to be accepted and silently discarded, which quietly produced an
            unstyled picture.
        """

        if basemap:
            values, _ = self.result._thickness_values(layer)
            # `GridPlots`, not `myflopy.plot`: that module is a layer ABOVE this
            # one in the import graph, so reaching it would need a deferred
            # import and the exact-match ratchet only moves down. Same function
            # either way -- `vor.plot.map` is what the front door calls too.
            return GridPlots(self.result.vor).map(values=list(values), **kwargs)
        if kwargs:
            raise TypeError(
                f"{', '.join(sorted(kwargs))} style the choropleth, which is "
                f"only drawn with basemap=True; the default thickness map is "
                f"Matplotlib and ignores them."
            )
        return LayerThicknessMap(self.result, layer=layer)

    def section(
        self,
        line=None,
        *,
        x=None,
        y=None,
        color_by: str = "layer",
        cmap: str = "tab10",
        show_grid: bool = True,
        legend: bool = True,
        title: str | None = None,
        **kwargs,
    ) -> LayerSection:
        """A filled, layer-coloured cross-section: ``stack.plot.section(y=300)``.

        Parameters
        ----------
        line : LineString, MultiLineString, path, or sequence of points, optional
            The section line: a shapely geometry, a path to a ``.shp``/``.gpkg``
            to read it from, or the points themselves as ``[(x, y), ...]``.
            Alternatively give ``x=`` or ``y=`` for an axis-aligned slice. With
            none of them, a West-East line through the grid centre.
        x, y : float, optional
            Draw the section along a constant x or constant y.
        color_by : {'layer', 'thickness'}, default 'layer'
            Cell scalar the fill is keyed to: discrete layer colours, or a
            continuous viridis thickness field.
        cmap : str, default 'tab10'
            Colormap for that scalar.
        show_grid : bool, default True
            Draw cell edges over the fill.
        legend : bool, default True
            Include the layer legend.
        title : str, optional
            Plot title. Defaults to "Layer cross-section".
        **kwargs
            Forwarded to :class:`LayerSection`.

        Returns
        -------
        LayerSection
            A Matplotlib :class:`~myflopy.viz.Picture`; ``.axes`` rather than
            ``.fig``.

        Raises
        ------
        TypeError
            If ``line`` is not a geometry, a path, or a sequence of points.
        ValueError
            If ``line`` resolves to fewer than two points, or ``color_by`` is
            not one of the values above. Both are raised HERE, not later from
            ``.axes``.
        """

        return LayerSection(
            self.result, line, x=x, y=y, color_by=color_by, cmap=cmap,
            show_grid=show_grid, legend=legend, title=title, **kwargs,
        )

    def surface(
        self,
        layer="top",
        *,
        backend: str = "plotly",
        resolution: int = 120,
        colorscale: str = "Earth_r",
        color_by: str | None = None,
        opacity: float | None = None,
        height: int | None = None,
        scale: float = 8,
        cmap: str = "tab10",
        width: int = 900,
        show_edges: bool = False,
        **kwargs,
    ):
        """One or more contacts as an interactive 3-D surface.

        ``layer`` is a name, a list of names, or ``"all"``; discover them with
        ``result.surface_names``. That is a height field ``z(x, y)``; the layered
        grid VOLUME is :meth:`grid`, a different shape entirely.

        Parameters
        ----------
        layer : str or list of str, default 'top'
            Which contact(s) to draw. ``"all"`` draws every one.
        resolution : int, default 120
            Interpolation grid size per axis.
        colorscale : str, default 'Earth_r'
            Plotly colorscale for the height field.
        color_by : str, optional
            Colour by a scalar other than elevation.
        opacity : float, optional
            Surface opacity, useful when stacking several contacts.
        height : int, optional
            Figure height in pixels.
        backend : {'plotly', 'vtk'}, default 'plotly'
            ``'plotly'`` overlays the surfaces in one figure. ``'vtk'`` renders
            each contact as its own interactive sheet in a PyVista scene, which
            is the one to reach for when you need to look at surfaces
            INDIVIDUALLY -- every sheet is a separate actor, so you can hide the
            ones above and see underneath. Needs the ``viz3d`` extra.
        scale : float, default 8
            *(vtk only)* Vertical exaggeration.
        cmap : str, default 'tab10'
            *(vtk only)* Colormap the sheets are coloured from, one per surface.
        width : int, default 900
            *(vtk only)* Scene width in pixels.
        show_edges : bool, default False
            *(vtk only)* Draw the interpolation mesh on each sheet.
        **kwargs
            Forwarded to :class:`LayerSurface` (plotly) or the scene builder.

        Returns
        -------
        LayerSurface or VtkScene
            A :class:`~myflopy.viz.Picture` either way. The VTK scene exposes
            ``.scene`` (the PyVista ``Plotter``) rather than ``.fig``.

        Raises
        ------
        ValueError
            If ``backend`` is neither ``'plotly'`` nor ``'vtk'``.
        """

        if backend == "vtk":
            plotter, sheets = self.result._vtk_surface_plotter(
                layer, resolution=resolution, scale=scale, cmap=cmap,
                opacity=1.0 if opacity is None else opacity,
                width=width, height=height or 580, show_edges=show_edges,
                **kwargs,
            )
            return VtkScene(plotter, meshes=sheets, title="layer surfaces")
        if backend != "plotly":
            raise ValueError(f"backend must be 'plotly' or 'vtk', not {backend!r}.")
        return LayerSurface(
            self.result, layer, resolution=resolution, colorscale=colorscale,
            color_by=color_by, opacity=opacity, height=height, **kwargs,
        )

    def grid(
        self,
        layers=None,
        *,
        backend: str = "vtk",
        color_by: str = "layer",
        scale: float = 8,
        cmap: str | None = None,
        width: int = 900,
        height: int = 580,
    ):
        """The layered grid mesh itself.

        ``backend="vtk"`` (the default here) renders the cell VOLUME in 3-D,
        coloured by layer -- the picture the old ``vtk_3d()`` drew, minus its
        habit of writing an HTML file into the working directory on every call.
        ``backend="plotly"`` gives the flat 2-D mesh instead, the same picture as
        ``vor.plot.grid()``.

        A ``backend`` switch is honest here because both branches draw the SAME
        subject -- this grid -- and differ only in renderer. That is why the 3-D
        volume is `grid`, not `surface`: `surface` means a height field.

        Parameters
        ----------
        layers : str or int or list, optional
            *(vtk only)* Which layers to show: a name, an index, or a list mixing
            them. Colours stay keyed to each layer's position, so a subset looks
            the same as it does in the full stack.
        backend : {'vtk', 'plotly'}, default 'vtk'
            ``'vtk'`` renders the 3-D volume and needs the ``viz3d`` extra;
            ``'plotly'`` draws the flat 2-D mesh.
        color_by : {'layer', 'thickness', 'top', 'botm', 'cellid', 'elevation'}, default 'layer'
            *(vtk only)* Cell scalar to colour by. ``'layer'`` is discrete, with
            each layer's name in the colour bar; the rest are continuous fields.
            ``'top'`` and ``'botm'`` are the flat per-cell contacts, where
            ``'elevation'`` is vertex z and so ramps within a cell.

            Every one of these is attached to the mesh regardless of which is
            drawn, so ``scene.meshes[0]`` carries them all into a ``.vtu``.
        scale : float, default 8
            *(vtk only)* Vertical exaggeration.
        cmap : str, optional
            *(vtk only)* Colormap for ``color_by``. Defaults to the one that
            suits the chosen scalar -- ``tab10`` for layers, ``viridis`` for
            thickness, ``terrain`` for elevations.
        width, height : int, default 900, 580
            *(vtk only)* Scene size in pixels.

        Returns
        -------
        VtkScene or GridMesh
            A :class:`~myflopy.viz.Picture` either way. ``VtkScene`` exposes
            ``.scene`` instead of ``.fig``.

        Raises
        ------
        ValueError
            If ``backend`` is neither value, if ``color_by`` is not one of the
            values above, or if a vtk-only argument is given with
            ``backend="plotly"``.
        """

        scene_args = {"color_by": color_by, "scale": scale, "cmap": cmap,
                      "width": width, "height": height}
        if backend == "plotly":
            stray = [n for n, v in scene_args.items() if v != _VTK_GRID_DEFAULTS[n]]
            if layers is not None:
                stray.insert(0, "layers")
            if stray:
                raise ValueError(
                    f"{', '.join(stray)} configure the 3-D scene; the flat mesh "
                    f"has no use for them. Drop backend='plotly', or drop these."
                )
            return GridPlots(self.result.vor).grid()
        if backend != "vtk":
            raise ValueError(
                f"backend must be 'vtk' or 'plotly', not {backend!r}."
            )
        plotter, meshes = self.result._vtk_plotter(layers, **scene_args)
        return VtkScene(plotter, meshes=meshes, title="layered grid")


class LayerStack:
    """Author a model's vertical layering from a top surface + named layers.

    ``LayerStack`` is the **friendly facade** for building MODFLOW layer geometry.
    You declare the model top, then ``.add(...)`` one named layer at a time (by its
    *bottom* surface or its *thickness*), and ``.build()`` resolves the stack into
    DISV-ready ``top`` / ``botm`` / ``idomain`` arrays. It is a thin layer over the
    :class:`~myflopy.surfaces.LayerSurfaces` engine -- it compiles to it (see
    :meth:`_layer_surfaces`) and never re-implements the sampling / reconcile /
    pinch-out logic. Each surface is an atomic :class:`~myflopy.surfaces.Surface`
    (use the aliases :class:`Raster`, :class:`Contours`, :class:`Points`,
    :class:`Flat`, :class:`Array`).

    What :meth:`build` does for you: samples every surface onto the grid cells
    (area-weighted by default), resolves them **top-down** so a flat/relative
    surface stays flat where it fits and is lowered only where it would intrude on
    the surface above ("flat where possible, fit between"), and turns sub-minimum
    or inverted thicknesses into **pinch-outs** (idomain) per your per-layer policy.

    Parameters
    ----------
    vor
        The grid (a ``VoronoiGridPlus``/vertex grid) whose cells the surfaces are
        sampled onto.
    top
        The model-top surface -- a :class:`Surface` or any value the aliases accept
        (a raster path via :class:`Raster`, an array via :class:`Array`, a constant
        via :class:`Flat`, ...).
    length_units, time_units
        Model units (default feet / days); ``length_units`` flows to DISV and is
        used to convert any surface declaring different ``units=``.

    Examples
    --------
    >>> from myflopy import LayerStack
    >>> from myflopy.layers import Raster, Contours
    >>> stack = (
    ...     LayerStack(vor, top=Raster("ground.tif"))   # land surface from a DEM
    ...     .add("alluvium",   thickness=25)            # 25-ft upper aquifer
    ...     .add("aquitard",   thickness=8, pinch="inactive")   # pinches out where thin
    ...     .add("bedrock",    bottom=Contours("bedrock_top.shp"))
    ... )
    >>> layers = stack.build()          # -> layers.top, layers.botm, layers.idomain, layers.nlay
    >>> print(stack.qc())               # NaN/thickness/connectivity report
    >>> layers.cross_section(x=500)     # quick W-E section to eyeball it

    Feed the result straight into ``mf.disv`` and the model context::

        gp = vor.get_disv_gridprops()
        ctx = mf.ModelContext(grid=vor, domain=layers.idomain)
        flow = mf.gwf("flow", context=ctx, packages=[
            mf.disv(nlay=layers.nlay, ncpl=gp["ncpl"], nvert=len(gp["vertices"]),
                    vertices=gp["vertices"], cell2d=gp["cell2d"],
                    top=layers.top, botm=layers.botm, idomain=layers.idomain),
            ...,
        ])

    See Also
    --------
    from_modflow : Seed an editable stack from an existing model's top/botm.
    build : Resolve the stack into DISV arrays (the per-layer pinch options).
    qc : Geometry quality-control report (NaN, thin/pinched, connectivity).
    """

    def __init__(
        self,
        vor=None,
        top=None,
        *,
        length_units: str = "feet",
        time_units: str = "days",
    ):
        """Start an editable layer stack with the model-top surface ``top``.

        Layers are added below the top with :meth:`add`; ``length_units`` /
        ``time_units`` are carried onto the built :class:`LayerBuildResult`.

        **``vor`` is optional.** Surfaces are lazy and ``add`` resolves nothing, so
        the whole layering -- names, sources, per-layer ``min_thickness`` and
        ``pinch`` -- can be declared before a grid exists::

            stack = (mf.LayerStack(top=ground)          # no grid yet
                     .add("sand", thickness=40.0, pinch="inactive")
                     .add("clay", bottom=clay_base))
            layers = stack.build(vor)                   # grid arrives here

        Supply it later to :meth:`build`, :meth:`qc` or :meth:`to_disv`, or bind
        it once with :meth:`for_grid`. Passing it here still works and is right
        whenever the grid already exists -- the deferred form exists so that
        declaring the layering does not force the grid to be built first.

        Raises
        ------
        TypeError
            If ``top`` is missing, or if a :class:`~myflopy.surfaces.Surface` is
            passed as ``vor`` -- ``LayerStack(ground)`` reads as "the deferred
            form" but binds ``ground`` to the grid, which would fail much later
            and somewhere else.
        """

        if top is None:
            hint = (
                " It looks like you passed the top surface positionally: the "
                "first positional argument is the GRID, so write "
                "`LayerStack(top=...)` for the deferred form, or "
                "`LayerStack(vor, top)`."
            ) if isinstance(vor, Surface) else ""
            raise TypeError(f"LayerStack needs a `top` surface.{hint}")
        if isinstance(vor, Surface):
            raise TypeError(
                "LayerStack's first argument is the grid, not a surface. Write "
                "`LayerStack(top=...)` to defer the grid, or `LayerStack(vor, top)`."
            )

        self.vor = vor
        self._top = _coerce_surface(top)
        self._layers: list[_Layer] = []
        #: Cached draft grid for `.plot` on a gridless stack; see `draft_grid`.
        self._draft = None
        self.length_units = length_units
        self.time_units = time_units

    def _require_grid(self, vor, caller: str):
        """Resolve the grid for ``caller``, preferring an explicit one.

        One message for every consumer, because "NoneType has no attribute
        ncpl" three frames down is the failure a deferred stack would otherwise
        produce.
        """

        resolved = self.vor if vor is None else vor
        if resolved is None:
            raise ValueError(
                f"LayerStack.{caller}() needs a grid. This stack was declared "
                f"without one, so pass it now -- `stack.{caller}(vor)` -- or bind "
                f"it once with `stack.for_grid(vor)`."
            )
        return resolved

    def _georeferenced_surfaces(self) -> list[Surface]:
        """Every surface in the stack that carries its own extent, composites included.

        A raster or contour surface knows where it is; ``flat``/``array`` and the
        algebra kinds do not, but they are *built from* ones that might, so the
        walk descends through ``operands``.
        """

        found, seen = [], set()

        def walk(surface):
            if id(surface) in seen:
                return
            seen.add(id(surface))
            if surface.kind in ("raster", "contours"):
                found.append(surface)
            for operand in surface.operands or ():
                if isinstance(operand, Surface):
                    walk(operand)

        walk(self._top)
        for layer in self._layers:
            if isinstance(layer.surface, Surface):
                walk(layer.surface)
        return found

    def draft_grid(self, *, cells: int = 400, extent=None, crs: str | None = None):
        """Build a throwaway Voronoi grid covering this stack's own extent.

        For LOOKING at a stack before its real grid exists. The extent comes from
        the first raster- or contour-backed surface in the stack -- those are the
        only kinds that know where they are -- and the mesh is deliberately coarse.

        **Not a modelling grid.** It has no boundary polygon, no refinement, and
        an extent that is whatever your DEM happens to cover rather than your
        model domain. :meth:`build`, :meth:`qc` and :meth:`to_disv` will not use
        it and still demand a real grid; only :attr:`plot` falls back to it.

        Parameters
        ----------
        cells : int, default 400
            Roughly how many cells to aim for. Coarse on purpose -- a few hundred
            renders in well under a second and is plenty to see shape, ordering
            and where a layer pinches.
        extent : tuple of float, optional
            ``(xmin, ymin, xmax, ymax)`` to use instead of the surfaces' own --
            for previewing over your real domain, or a corner of it.
        crs : str, optional
            Overrides the CRS read from the source raster.

        Returns
        -------
        VoronoiGridPlus

        Raises
        ------
        ValueError
            If no surface carries an extent and none was given. A stack of
            ``Flat``/``Array`` surfaces describes thicknesses with no notion of
            where they are, so there is nothing to infer.
        """

        if extent is None:
            rasterio = require("rasterio", feature="reading a surface extent for draft_grid()")
            for surface in self._georeferenced_surfaces():
                source = surface.resolve_source()
                if isinstance(source, Path):
                    with rasterio.open(source) as handle:
                        extent = tuple(handle.bounds)
                        crs = crs or str(handle.crs)
                    break
            else:
                raise ValueError(
                    "draft_grid() needs an extent: no surface in this stack carries "
                    "one. `Raster` and `Contours` know where they are; `Flat`, "
                    "`Array` and the algebra built on them do not. Pass "
                    "extent=(xmin, ymin, xmax, ymax), or use a real grid."
                )

        xmin, ymin, xmax, ymax = (float(v) for v in extent)
        width, height = xmax - xmin, ymax - ymin
        if width <= 0 or height <= 0:
            raise ValueError(
                f"draft_grid() got an empty extent: {(xmin, ymin, xmax, ymax)}."
            )

        tri = TriangleGrid(model_ws=tempfile.mkdtemp(prefix="draft_grid_"), angle=30)
        tri.set_domain_rectangle(
            x_dist=width, y_dist=height, origin=(xmin, ymin),
            max_area=(width * height) / max(int(cells), 1),
        )
        tri.build()
        return VoronoiGridPlus(tri, crs=crs)

    def for_grid(self, vor) -> LayerStack:
        """Return a copy of this stack bound to ``vor``.

        The layering is the expensive thing to write and the grid is the thing
        you change, so binding returns a NEW stack rather than mutating this one:
        one declaration can serve a coarse test grid and a fine production grid
        without being restated.

            coarse = stack.for_grid(vor_coarse).build()
            fine   = stack.for_grid(vor_fine).build()

        Also what makes ``.plot`` reachable on a deferred stack, since a property
        cannot take a grid argument.
        """

        bound = LayerStack(
            vor,
            self._top,
            length_units=self.length_units,
            time_units=self.time_units,
        )
        bound._layers = list(self._layers)
        return bound

    @classmethod
    def from_modflow(
        cls,
        vor,
        source,
        *,
        names: list[str] | None = None,
        resample: bool = True,
        length_units: str = "feet",
        time_units: str = "days",
    ) -> LayerStack:
        """Build a stack on ``vor`` from an existing MODFLOW model's top/botm.

        ``source`` is a flopy model or modelgrid; its top and per-layer bottoms
        become this stack's surfaces (interpolated onto ``vor`` when
        ``resample=True``, used verbatim when ``False``). Name the layers with
        ``names`` (defaults to ``layer1..N``). Edit the returned stack like any
        other -- e.g. ``.replace(...)`` a bottom, then ``.build()`` / ``.to_disv()``.
        """
        top_s, botm_s = modflow_surfaces(source, resample=resample)
        if names is None:
            names = [f"layer{i + 1}" for i in range(len(botm_s))]
        if len(names) != len(botm_s):
            raise ValueError(f"{len(names)} names given for {len(botm_s)} layers.")
        stack = cls(vor, top=top_s, length_units=length_units, time_units=time_units)
        for name, surface in zip(names, botm_s, strict=False):
            stack.add(name, bottom=surface)
        return stack

    # -- authoring -------------------------------------------------------- #
    def add(
        self,
        name: str,
        *,
        bottom=None,
        thickness=None,
        min_thickness: float | None = None,
        pinch: str | None = None,
        fill: str | None = None,
        split=None,
        names=None,
    ) -> LayerStack:
        """Append a named geologic unit beneath the current bottom. Returns ``self``.

        Define the unit either by its **bottom** surface or its **thickness**
        (exactly one). Thickness is measured down from the surface above, so layers
        compose naturally as you stack them.

        ``split=`` discretizes the unit into several MODEL layers without
        changing its geometry -- see below.

        Parameters
        ----------
        name
            Layer name (used in QC, plots, and as the surface label).
        bottom
            The layer's bottom as a :class:`Surface` / alias (e.g.
            ``Contours("base.shp")``, ``Array(values)``, ``Flat(90.0)``). Mutually
            exclusive with ``thickness``.
        thickness
            Constant thickness below the surface above. Mutually exclusive with
            ``bottom``. For a per-cell thickness MAP write ``bottom=Isopach(map)``
            -- a bare Surface here is an elevation, not a thickness, and is
            refused for that reason.
        min_thickness
            Minimum thickness for the UNIT; thinner cells are handled per
            ``pinch``. Defaults to the stack-wide value passed to :meth:`build`.
            Not divided among a split unit's layers -- see ``split``.
        pinch
            What to do where the unit is thinner than ``min_thickness``:
            ``"passthrough"`` (keep the cell active, default), ``"inactive"``
            (idomain 0 -- a true pinch-out), or ``"floor"`` (clamp to the minimum).
        fill
            For raster/derived bottoms, how to fill cells with no source data
            (e.g. ``"propagate"`` to inherit the surface above -> pinch).
        split
            Discretize this unit into several model layers. ``split=3`` gives
            three equal shares; ``split=[0.3, 0.7]`` gives those shares, which
            must sum to 1. The layers are named ``<name>_1 .. <name>_N``, and
            ``<name>`` remains addressable through
            :attr:`LayerBuildResult.units`.

            The unit's GEOMETRY is unchanged: the last cut lands exactly on the
            declared bottom. ``min_thickness`` and ``pinch`` stay unit-scoped, so
            a unit either pinches out whole or not at all -- evaluating them per
            slice produces holes inside a unit that physically exists.

            A unit whose bottom is measured from the layer above (an
            ``Isopach``) cannot be split; give it an absolute bottom, or declare
            it with ``thickness=`` and split that.
        names
            Names for the split layers, top to bottom, replacing the generated
            ``<name>_1 .. <name>_N``. Requires ``split`` of 2 or more, and must
            give one name per layer -- naming only some of them would leave
            which layer a name refers to depending on where you stopped
            counting. The UNIT keeps ``name``, so
            :attr:`LayerBuildResult.units` and
            :meth:`LayerBuildResult.per_layer` are unaffected by the choice.

        Examples
        --------
        >>> stack.add("sand", thickness=20)                       # 20-ft layer
        >>> stack.add("clay", thickness=5, pinch="inactive")      # pinches out where thin
        >>> stack.add("bedrock", bottom=Contours("bedrock.shp"))  # bottom from contours
        >>> stack.add("sand", bottom=Raster("base.tif"), split=3) # -> sand_1..sand_3
        >>> stack.add("till", thickness=60, split=[0.25, 0.75])   # -> 15 ft, then 45 ft
        >>> stack.add("sand", bottom=Raster("base.tif"), split=3,  # -> your names,
        ...           names=["upper sand", "mid sand", "lower sand"])  # unit stays "sand"
        """
        _check_layer_name(name, self._layers)
        surface = _make_surface(bottom, thickness, fill)
        shares = _split_shares(split, name)  # reject a bad split here, not at build
        _split_names(names, shares, name)    # ...and a bad names= list with it
        self._layers.append(_Layer(name, surface, min_thickness, pinch, split, names))
        return self

    def insert_below(
        self, name: str, new_name: str, *, bottom=None, thickness=None,
        min_thickness: float | None = None, pinch: str | None = None,
        fill: str | None = None, split=None, names=None,
    ) -> LayerStack:
        """Insert a new unit directly below the existing unit ``name``."""
        _check_layer_name(new_name, self._layers)
        surface = _make_surface(bottom, thickness, fill)
        _split_names(names, _split_shares(split, new_name), new_name)
        self._layers.insert(
            self._index(name) + 1,
            _Layer(new_name, surface, min_thickness, pinch, split, names),
        )
        return self

    def replace(
        self, name: str, *, bottom=None, thickness=None,
        min_thickness=_UNSET, pinch=_UNSET, fill=None, split=_UNSET, names=_UNSET,
    ) -> LayerStack:
        """Update an existing unit in place; unspecified fields are kept.

        This is also how a unit is split after the fact --
        ``stack.replace("sand", split=3)``. There is deliberately no ``.split()``
        verb: this method already edits a declared unit by name, already returns
        ``self``, and already keeps what you do not mention.

        ``split`` and ``names`` are kept independently but validated TOGETHER,
        so changing one to a length the other cannot match raises here rather
        than at build. Dropping a split from a named unit therefore reads
        ``replace(name, split=None, names=None)`` -- the alternative, silently
        discarding names the caller wrote, is the kind of quiet renaming this
        module exists to prevent.
        """
        idx = self._index(name)
        layer = self._layers[idx]
        if bottom is not None or thickness is not None:
            surface = _make_surface(bottom, thickness, fill)
        else:
            surface = layer.surface if fill is None else _dc_replace(layer.surface, fill=fill)
        resolved_split = layer.split if split is _UNSET else split
        resolved_names = layer.names if names is _UNSET else names
        if split is not _UNSET or names is not _UNSET:
            _split_names(
                resolved_names, _split_shares(resolved_split, name), name
            )
        self._layers[idx] = _Layer(
            name,
            surface,
            layer.min_thickness if min_thickness is _UNSET else min_thickness,
            layer.pinch if pinch is _UNSET else pinch,
            resolved_split,
            resolved_names,
        )
        return self

    def remove(self, name: str) -> LayerStack:
        """Remove a layer by name."""
        del self._layers[self._index(name)]
        return self

    @property
    def names(self) -> list[str]:
        """The layer names, top-to-bottom (excluding the model top)."""

        return [layer.name for layer in self._layers]

    def _index(self, name: str) -> int:
        """The position of the layer named ``name`` (raises ``KeyError`` if absent)."""

        for i, layer in enumerate(self._layers):
            if layer.name == name:
                return i
        raise KeyError(f"no layer named {name!r} (have {self.names}).")

    # -- compilation ------------------------------------------------------ #
    def _expanded_layers(self) -> tuple[list[_Layer], dict[str, list[int]]]:
        """Expand each declared UNIT into the model layers it becomes.

        Returns the flat layer list and the ``unit name -> layer indices`` map
        that :attr:`LayerBuildResult.units` carries forward. A unit with no
        ``split`` passes through unchanged and keeps its own name, so declaring
        ``split=1`` -- or adding a split later -- never renames anything. A unit
        carrying ``names`` uses those instead of the generated ones; either way
        the UNIT key in the returned map is the declared ``name``.
        """

        out: list[_Layer] = []
        units: dict[str, list[int]] = {}
        for layer in self._layers:
            shares = _split_shares(layer.split, layer.name)
            sub_names = _split_names(layer.names, shares, layer.name)
            first = len(out)
            if len(shares) == 1:
                out.append(_dc_replace(layer, split=None, names=None))
            else:
                steps = _split_steps(shares)
                for i, (share, step) in enumerate(zip(shares, steps, strict=True), 1):
                    out.append(_dc_replace(
                        layer,
                        name=sub_names[i - 1] if sub_names else f"{layer.name}_{i}",
                        surface=_sub_surface(layer.surface, share, step, layer.name),
                        split=None,
                        names=None,
                    ))
            units[layer.name] = list(range(first, len(out)))

        seen: dict[str, str] = {}
        for unit, indices in units.items():
            for i in indices:
                clash = seen.get(out[i].name)
                if clash is not None:
                    raise ValueError(
                        f"splitting {unit!r} produces the layer name "
                        f"{out[i].name!r}, which unit {clash!r} already uses. "
                        f"Layer names identify a contact, a DataFrame column and "
                        f"a 3-D scene actor, so they have to be distinct -- "
                        f"rename one of the two units."
                    )
                seen[out[i].name] = unit
        return out, units

    def _layer_surfaces(self) -> LayerSurfaces:
        """Compile to the :class:`LayerSurfaces` engine: the top plus each layer bottom, labeled."""

        layers, _ = self._expanded_layers()
        surfaces = [self._top] + [layer.surface for layer in layers]
        labels = ["top"] + [layer.name for layer in layers]
        return LayerSurfaces(surfaces, labels=labels)

    def _per_layer_config(self, default_min_thickness, default_pinch):
        """Resolve each layer's ``(min_thickness, pinch)``, filling unset ones with the defaults.

        Per MODEL layer, but the values come from the UNIT: every slice of a
        split unit inherits the unit's policy undivided, which is what makes the
        pinch verdict a property of the unit rather than of an arbitrary
        discretization choice.
        """

        layers, _ = self._expanded_layers()
        min_thk = [
            default_min_thickness if layer.min_thickness is None else layer.min_thickness
            for layer in layers
        ]
        pinch = [
            default_pinch if layer.pinch is None else layer.pinch
            for layer in layers
        ]
        return min_thk, pinch

    @staticmethod
    def _reconcile_separations(min_thk, pinch, fallback: float, units=None):
        """Per-layer ``(min_sep, trigger_sep)`` for reconcile, from each layer's policy.

        A ``"floor"`` layer declares a real minimum thickness, so reconcile enforces
        **its own** ``min_thickness`` -- that is what makes ``min_thickness`` set
        geometry rather than only decide an idomain value.

        A ``"passthrough"`` or ``"inactive"`` layer must be allowed to come out too
        thin, because being too thin is the signal that it pinches out. Reconcile
        gives those only the ``fallback`` separation, just enough to keep the stack
        ordered, and ``min_thickness`` stays a pure threshold for them.

        ``min_thickness`` is declared for the UNIT, so a split unit's separation is
        divided among its slices: three slices of a 3 ft minimum get 1 ft each and the
        unit still comes out 3 ft. Handing every slice the unit's own figure would
        inflate the unit by its split factor -- the same failure ledger 158 records
        for the pinch threshold, which bites the geometry just as hard.
        """

        shares = [1] * len(min_thk)
        for indices in (units or {}).values():
            for index in indices:
                if 0 <= index < len(shares):
                    shares[index] = len(indices)

        seps, triggers = [], []
        for thickness, policy, n in zip(min_thk, pinch, shares, strict=False):
            if policy == "floor":
                seps.append(float(thickness) / n)
                triggers.append(float(thickness) / n)
            else:
                seps.append(float(fallback))
                triggers.append(None)      # caller fills the stack-wide trigger
        return seps, triggers

    def _max_split(self) -> int:
        """The largest number of model layers any one unit becomes."""

        return max(
            (len(_split_shares(layer.split, layer.name)) for layer in self._layers),
            default=1,
        )

    def _trigger_sep(self, trigger_sep):
        """Resolve ``trigger_sep=None`` to a value that suits the discretization.

        Reconcile rewrites any layer thinner than ``trigger_sep`` to exactly
        ``min_sep``. The engine's default of 1.0 length unit is sized for
        geologic units, and applying it unchanged to their SLICES is destructive:
        a 2.4 ft unit split three ways comes back ``[0.1, 1.5, 0.1]`` -- 1.7 ft
        of 2.4, with the base moved and the layer below silently absorbing the
        difference. Dividing by the largest split restores it to an exact
        ``[0.8, 0.8, 0.8]``, because the trigger then means the same fraction of
        a slice that 1.0 meant of a unit.

        Pass a number to override. An unsplit stack is unaffected: the maximum
        split is 1, so this returns the historical 1.0.
        """

        if trigger_sep is not None:
            return float(trigger_sep)
        return 1.0 / self._max_split()

    def refresh(self) -> LayerStack:
        """Rebuild the cached raster of every derived (contour) surface now."""
        for surface in [self._top] + [layer.surface for layer in self._layers]:
            if getattr(surface, "is_derived", False):
                surface.resolve_source(refresh=True, warn=False)
        return self

    def cache_status(self) -> dict[str, str]:
        """Map ``name -> "missing"/"fresh"/"stale"`` for each derived surface."""
        named = [("top", self._top)] + [(layer.name, layer.surface) for layer in self._layers]
        return {
            name: surface.cache_status()
            for name, surface in named
            if getattr(surface, "is_derived", False)
        }

    def build(
        self,
        vor=None,
        *,
        default_min_thickness: float = 1.0,
        default_pinch: str = "floor",
        reconcile="bottom",
        min_sep: float = 0.1,
        trigger_sep: float | None = None,
        method: str = "area",
        refresh: bool = False,
        attach: bool = False,
    ) -> LayerBuildResult:
        """Resolve, reconcile, and pinch out the stack into DISV-ready arrays.

        Samples every surface onto the grid, reconciles crossing surfaces, applies
        the per-layer pinch policy, and returns the ``top``/``botm``/``idomain``
        arrays ready for ``mf.disv(...)``. A stale derived (contour) surface is
        reused with a warning; pass ``refresh=True`` to rebuild its cache first.

        Parameters
        ----------
        vor : VoronoiGridPlus, optional
            The grid to sample onto. Required only if this stack was declared
            without one (``LayerStack(top=...)``); an explicit grid here wins
            over the stack's own. See :meth:`for_grid` to bind one instead.
        default_min_thickness : float, default 1.0
            Minimum layer thickness for layers that do not declare their own with
            ``.add(..., min_thickness=)``.
        default_pinch : str, default "floor"
            Thin-layer policy for layers that do not declare their own: ``"floor"``
            (hold the cell open at its minimum thickness), ``"inactive"`` (idomain
            0, a true pinch-out), or ``"passthrough"`` (idomain -1, flow passes
            vertically through and the cell can carry NO boundary condition).
        reconcile : {'bottom', 'top', True, False, None}, default 'bottom'
            What to do where two surfaces cross or come closer than 1 length unit
            of each other -- interpolated contacts routinely do, and a crossing
            means a negative layer thickness.

            ``'bottom'``
                Trust the surface ABOVE: push the lower one down to
                ``upper - min_sep``. **The usual choice**, and the only one that
                leaves your model top where you put it -- which matters when the
                top is a measured DEM. Corrections cascade downward, so one
                intruding contact pushes every contact below it clear.
            ``'top'``
                Trust the surface BELOW: raise the upper one to
                ``lower + min_sep``. Corrections propagate upward and CAN MOVE
                THE MODEL TOP -- with ``top=100, a=105`` it returns ``top=105.1``.
                Reach for it only when the deeper contact is the reliable one
                (say a well-picked bedrock surface under a coarse interpolated
                top), and check ``qc()`` afterwards.
            ``True``
                Same as ``'bottom'``.
            ``False`` or ``None``
                Do not reconcile. Crossing surfaces stay crossed, so layers get
                zero or negative thickness and the pinch policy is what saves
                you -- useful for seeing the raw sampled surfaces, not for
                building a model.

            Either way ``qc()`` reports, per layer, how many cells reconcile had
            to move and the largest move, which is how you tell "tidied two
            cells" from "rebuilt the geometry".
        min_sep : float, default 0.1
            **Fallback** separation, used only for layers that pinch
            (``"passthrough"`` / ``"inactive"``) -- enough to keep the stack ordered
            while still letting them come out thin enough to pinch. A ``"floor"``
            layer ignores it and uses its own ``min_thickness`` instead, which is
            how ``min_thickness`` comes to set geometry.
        trigger_sep : float, optional
            The separation below which reconcile ACTS -- surfaces closer than
            this are treated as conflicting even if they never actually cross,
            and the thinner layer is rewritten to ``min_sep``. Defaults to
            ``1.0 / <largest split>``: the historical 1 length unit for an
            unsplit stack, and proportionally smaller once a unit is cut into
            slices, so that splitting a unit does not shrink it. Pass a number to
            override.
        method : str, default "area"
            Raster sampling method: ``"area"`` (area-weighted) or ``"centroid"``.
        refresh : bool, default False
            Rebuild any stale derived (contour) surface caches before sampling.
        attach : bool, default False
            Also publish the result onto ``vor.gdf_topbtm`` (see
            :meth:`LayerBuildResult.attach_to_grid`) so surface-aware builders --
            ``mf.sfr`` reach tops, ``mf.lak`` lake-cell layering -- can read the
            elevations. Use it when your model has SFR/LAK on this grid.

        Returns
        -------
        LayerBuildResult
            Bundles ``top`` ``(ncpl,)``, ``botm`` / ``idomain`` ``(nlay, ncpl)``,
            derived ``thickness``, and per-layer metadata.

        Raises
        ------
        ValueError
            If no layers have been added with :meth:`add`.

        Examples
        --------
        >>> layers = (mf.LayerStack(vor, top=mf.Raster("ground.tif"))
        ...           .add("sand", thickness=20)
        ...           .add("clay", bottom=mf.Contours("base.shp"), pinch="inactive")
        ...           .build(attach=True))
        >>> mf.disv(nlay=layers.nlay, ..., top=layers.top, botm=layers.botm,
        ...         idomain=layers.idomain)
        """
        if not self._layers:
            raise ValueError("Add at least one layer with .add(...) before build().")
        vor = self._require_grid(vor, "build")
        layers, units = self._expanded_layers()
        ls = self._layer_surfaces()
        rec_on, which = _reconcile_args(reconcile)
        min_thk, pinch = self._per_layer_config(default_min_thickness, default_pinch)
        seps, triggers = self._reconcile_separations(min_thk, pinch, min_sep, units)
        stack_trigger = self._trigger_sep(trigger_sep)
        if trigger_sep is not None:
            triggers = [stack_trigger] * len(triggers)
        else:
            triggers = [stack_trigger if v is None else v for v in triggers]
        gdf = ls.sample(
            vor, reconcile=rec_on, which=which, min_sep=seps,
            trigger_sep=triggers,
            method=method, length_units=self.length_units, refresh=refresh,
        )
        top, botm = ls._split_top_botm(gdf)
        thickness = ls._thickness(gdf)
        ls._validate_pinch_invariant(
            min_thk, pinch, thickness.shape[0],
            {"reconcile": rec_on, "min_sep": min_sep},
        )
        idomain = ls._idomain_from_thickness(
            _unit_thickness(thickness, units), min_thk, pinch
        )
        result = LayerBuildResult(
            top=top, botm=botm, idomain=idomain, thickness=thickness,
            names=[layer.name for layer in layers],
            min_thickness=min_thk, pinch=pinch,
            length_units=self.length_units, time_units=self.time_units,
            vor=vor, units=units,
        )
        if attach:
            result.attach_to_grid(vor)
        return result

    def qc(
        self,
        vor=None,
        *,
        default_min_thickness: float = 1.0,
        default_pinch: str = "floor",
        reconcile="bottom",
        min_sep: float = 0.1,
        trigger_sep: float | None = None,
        method: str = "area",
        refresh: bool = False,
    ) -> LayerQCReport:
        """Build the stack and run QC, including reconcile diagnostics.

        On top of :meth:`LayerBuildResult.qc` (NaN coverage, isolated active
        cells, thickness), this samples the surfaces *without* reconciling and
        reports, per layer, how many cells reconcile had to move and the largest
        move -- showing where surfaces were crossing before reconcile fixed them.

        Read this before trusting a build. A large ``reconcile_moved`` count
        means the geometry you got is substantially not the geometry you
        described, which is worth knowing before it becomes an idomain.

        Parameters
        ----------
        vor : VoronoiGridPlus, optional
            Required only if the stack was declared without a grid.
        default_min_thickness, default_pinch, reconcile, min_sep, trigger_sep, method, refresh
            Passed straight through to :meth:`build`, so QC reports on the
            geometry you are actually going to build. Documented there -- in
            particular the ``reconcile`` options, which are what this report is
            diagnosing.

        Returns
        -------
        LayerQCReport
        """
        vor = self._require_grid(vor, "qc")
        result = self.build(
            vor,
            default_min_thickness=default_min_thickness, default_pinch=default_pinch,
            reconcile=reconcile, min_sep=min_sep, trigger_sep=trigger_sep,
            method=method, refresh=refresh,
        )
        report = result.qc()
        ls = self._layer_surfaces()
        raw = ls.sample(
            vor, reconcile=False, method=method,
            length_units=self.length_units, refresh=refresh,
        )
        _, raw_botm = ls._split_top_botm(raw)
        adjusted, max_shift = [], []
        for k in range(result.nlay):
            shift = np.abs(np.asarray(raw_botm[k], float) - np.asarray(result.botm[k], float))
            finite = shift[np.isfinite(shift)]
            adjusted.append(int((finite > 1e-6).sum()))
            max_shift.append(float(finite.max()) if finite.size else 0.0)
        report.reconcile_adjusted = adjusted
        report.reconcile_max_shift = max_shift
        return report

    # -- views: build with defaults, then view the result ----------------- #
    @property
    def plot(self) -> StackPlots:
        """The plotting verbs for this stack: ``map``, ``section``, ``surface``.

        Builds the stack with default options and returns the result's namespace,
        so ``stack.plot.map()`` is ``stack.build().plot.map()``. Build with
        non-default options first if you need them.

        **Works without a grid.** A stack declared with none draws on a coarse
        :meth:`draft_grid` covering the surfaces' own extent, so you can look at
        the layering before committing to a mesh. One INFO line records the
        extent and cell count, because a blocky preview should not be mistaken
        for your model. The draft is cached, so repeated ``.plot`` calls reuse it.

        Only pictures do this. :meth:`build`, :meth:`qc` and :meth:`to_disv`
        still require a real grid -- an approximate picture is useful, invented
        model geometry is not. Control the draft with
        ``stack.for_grid(stack.draft_grid(cells=2000)).plot``.
        """

        if self.vor is not None:
            return self.build().plot

        if self._draft is None:
            self._draft = self.draft_grid()
            xmin, ymin, xmax, ymax = self._draft.gdf_vorPolys.total_bounds
            logger.info(
                "LayerStack.plot: no grid on this stack, drawing on a draft grid "
                "of %d cells over (%.0f, %.0f)-(%.0f, %.0f). Pictures only -- "
                "build()/qc()/to_disv() still need a real grid.",
                int(self._draft.ncpl), xmin, ymin, xmax, ymax,
            )
        return self.build(self._draft).plot


    def to_disv(
        self,
        vor=None,
        *,
        name: str = "disv",
        attach: bool = True,
        default_min_thickness: float = 1.0,
        default_pinch: str = "floor",
        reconcile="bottom",
        min_sep: float = 0.1,
        trigger_sep: float | None = None,
        method: str = "area",
        refresh: bool = False,
    ):
        """Build a ready-to-use ``mf.disv`` spec (with pinch-out idomain).

        The one-call alternative to ``build()`` + hand-written ``mf.disv(...)``,
        and the form worth preferring: it returns a
        :class:`~myflopy.specs.PackageSpec`, so the stack can be registered on a
        project (``project.add_package("disv/base", spec)``) and referenced with
        ``mf.ref``, which is what keeps an array-heavy model persistable.

        Parameters
        ----------
        vor : VoronoiGridPlus, optional
            Required only if the stack was declared without a grid.
        name : str, default "disv"
            Package name on the built model.
        attach : bool, default True
            Publish the sampled elevations onto ``vor.gdf_topbtm`` -- note this
            defaults to True here and to False on :meth:`build`. Surface-aware
            builders (``mf.sfr`` reach tops, ``mf.lak`` lake-cell layering) and
            the layer-elevation hover rows read it.
        default_min_thickness, default_pinch, reconcile, min_sep, method, refresh
            Passed straight through to the sampling pass; documented on
            :meth:`build`, including the ``reconcile`` options.

        Returns
        -------
        PackageSpec
            A DISV spec with pinch-out idomain applied.
        """
        if not self._layers:
            raise ValueError("Add at least one layer before to_disv().")
        vor = self._require_grid(vor, "to_disv")
        ls = self._layer_surfaces()
        rec_on, which = _reconcile_args(reconcile)
        min_thk, pinch = self._per_layer_config(default_min_thickness, default_pinch)
        # `idomain=` explicitly, rather than `pinch_out=True`: the engine would
        # judge each layer on its own thickness, and a split unit has to be
        # judged whole (see `_unit_thickness`). Letting it decide would give a
        # DIFFERENT idomain here than `build()` produces from the same stack.
        result = self.build(
            vor,
            default_min_thickness=default_min_thickness, default_pinch=default_pinch,
            reconcile=reconcile, min_sep=min_sep, trigger_sep=trigger_sep,
            method=method, refresh=refresh,
        )
        # ...and the GEOMETRY has to be resolved the same way for the same reason.
        # `ls.to_disv` samples the stack a second time, so handing it the scalar
        # `min_sep` while `build()` used per-layer separations writes a DISV whose
        # bottoms disagree with `result.botm` -- measured at 8.5 ft on a real stack.
        # That is silent and dangerous: `attach_to_grid()` publishes `result`'s
        # surfaces, so `mf.CellSurfaceOffset("cell_bottom", ...)` places boundaries
        # against bottoms MF6 never sees, and drains land below their cell.
        _, units = self._expanded_layers()
        seps, triggers = self._reconcile_separations(min_thk, pinch, min_sep, units)
        stack_trigger = self._trigger_sep(trigger_sep)
        if trigger_sep is not None:
            triggers = [stack_trigger] * len(triggers)
        else:
            triggers = [stack_trigger if v is None else v for v in triggers]
        return ls.to_disv(
            vor,
            idomain=result.idomain,
            minimum_thickness=min_thk,
            pinch=pinch,
            length_units=self.length_units,
            name=name,
            attach=attach,
            reconcile=rec_on,
            which=which,
            min_sep=seps,
            trigger_sep=triggers,
            method=method,
            refresh=refresh,
        )


__all__ = [
    "LayerStack",
    "LayerBuildResult",
    "LayerQCReport",
    "modflow_surfaces",
    "Raster",
    "Flat",
    "Contours",
    "Points",
    "Array",
    "Isopach",
    "Toward",
    "Min",
    "Max",
    "Clamp",
    "Where",
]

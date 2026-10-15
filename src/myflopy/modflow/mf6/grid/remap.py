"""Carry a package's per-cell arrays from one mesh onto another.

B3's C6. The overlay engine: it intersects the OLD grid's own cell polygons with the
NEW grid's and uses the resulting areas to combine old values into new cells. That is
what lets a model whose arrays came from nowhere reproducible -- a PEST-calibrated K
field, an array read out of an imported model -- move to another mesh at all, where
``with_grid`` could previously only stamp it ``unverifiable`` and let the build refuse it.

**A remap is the WORSE route whenever a source exists, and that is measured.**
Re-declaring the same polygon map on the new mesh (``mf.Polygons``, ledger 299) held to a
maximum 7.1% and a median 2.1% boundary-flux error against the fine mesh, where the best
reducer's own worst case was 12.5% and the worst was four orders of magnitude. So this is
the route for an array with NO recipe, not a general-purpose grid changer, and the
warning says so.

**What it deliberately does not serve** (ledger 300, each with its measurement):
an EXTENSIVE quantity -- a rate or a total -- because the right denominator is the OLD
cell's area rather than the new cell's covered area, and because ``FieldSpec.kind``
cannot tell a rate from a total; and an INTEGER FLAG, because averaging a cell type or a
zero-based layer index produces a number MF6 has no meaning for.
"""

from __future__ import annotations

import warnings

import numpy as np

from myflopy._logging import get_logger
from myflopy._numbers import is_dry, is_no_value
from myflopy._warnings import stacklevel_outside_package
from myflopy.modflow.mf6.grid.geometry import (
    REMAP_RULES,
    refuse_a_different_frame,
)

logger = get_logger(__name__)

#: Per-cell arrays remapped with the plain area-weighted mean. Each is a quantity whose
#: value does not scale with the cell holding it and which may legitimately be zero or
#: negative -- a head especially -- so the geometric and harmonic rules, which take a
#: logarithm and a reciprocal, are not offered for them.
INTENSIVE_ARRAYS: frozenset[str] = frozenset({"strt", "ss", "sy"})

#: Conductivities. These are the reason ``rule=`` exists: the correct average depends on
#: the FLOW DIRECTION, which is a property of the solution and not of the field, so no
#: default can be right everywhere. Measured in MF6 on one 2:1 coarsening -- with the
#: contrast ACROSS the flow the arithmetic mean gave 264x the fine-mesh boundary flux and
#: the harmonic mean was nearly exact; with the SAME contrast rotated 90 degrees the
#: arithmetic mean was exact and the harmonic wrong by about 250x. The ordering REVERSES.
CONDUCTIVITY_ARRAYS: frozenset[str] = frozenset({"k", "k22", "k33"})

#: Integer flags, refused by name. Every rule here returns a float and two of the four
#: can return one that is in no source cell. Measured on an 8x8 -> 5x5 coarsening of an
#: ``icelltype`` contact: the mean gave 5 of 25 cells values like 0.3851 and 0.5459, and
#: even ``majority``, which invents nothing, moved the convertible AREA fraction from
#: 0.4949 to 0.4569. ``irch``/``ievt`` are worse again -- they are zero-based LAYER
#: INDICES, so an average of them is meaningless in a way nothing downstream catches.
FLAG_ARRAYS: frozenset[str] = frozenset({"icelltype", "iconvert", "irch", "ievt"})

#: What a new cell with NO old cell under it means.
UNCOVERED_MODES: tuple[str, ...] = ("nearest", "raise")

#: The conductivity rule and the uncovered mode a remap uses when the caller names none:
#: geometric is never best and never catastrophic (see ``CONDUCTIVITY_ARRAYS``).
DEFAULT_RULE = "geometric"
DEFAULT_UNCOVERED = "nearest"

#: Ratio between the arithmetic and harmonic means above which the choice of rule is
#: reported as material. 2.0 is a factor of two in a conductivity, which is the point at
#: which a modeller would want to know rather than a threshold with physics behind it.
_SPREAD_REPORTED_ABOVE = 2.0


class RemapRefused(ValueError):
    """A package's arrays cannot be carried onto the new mesh, and why."""


def _nearest_values(old_xy, new_xy, values, *, usable=None):
    """Each point of ``new_xy`` takes the value of the nearest USABLE old point.

    Returns ``(values, worst_distance)``, or ``(None, None)`` when no old point is usable.
    Point to point through a KD-tree, on whichever points the caller names: a grid move's
    package arrays pass POLYGON centroids (what ``with_grid(remap=True)`` has always
    used), a solved field passes MF6's own cell centres -- the point a head IS the value
    at, which on a Voronoi mesh sits up to 0.6 of a cell width from the centroid
    (measured 44.9 ft on Ten Trails). Polygon distance is the more defensible measure
    and is what ledger 288 chose for snapping a FEATURE to a cell -- the two measured
    3.07x apart on one real miss -- but a feature is one geometry against a grid, where
    this is thousands against thousands, and the answer only ever fills a cell that has
    NO data at all. The distance is RETURNED rather than logged, because it is the number
    the caller's warning has to name.
    """

    from scipy.spatial import cKDTree

    values = np.asarray(values, dtype=float)
    keep = np.isfinite(values)
    if usable is not None:
        keep &= np.asarray(usable, dtype=bool)
    if not keep.any():
        return None, None
    old_xy = np.asarray(old_xy, dtype=float)
    distance, index = cKDTree(old_xy[keep]).query(np.asarray(new_xy, dtype=float))
    return values[keep][index], float(np.max(distance))


def _centroids(polys, rows=None) -> np.ndarray:
    """``(n, 2)`` polygon centroids of ``polys``, optionally at positional ``rows``."""

    centres = polys.geometry.centroid
    xy = np.column_stack((centres.x.to_numpy(), centres.y.to_numpy()))
    return xy if rows is None else xy[rows]


def check_remap_choices(rule: str, uncovered: str) -> None:
    """Raise ``ValueError`` naming the choices when ``rule``/``uncovered`` is not one.

    The one check, run by :func:`remap_arrays` and at the ``with_grid`` CALL -- which
    must not wait for a package with an array to carry, or a typo is accepted
    whenever there is nothing to remap.
    """

    if rule not in REMAP_RULES:
        raise ValueError(
            f"rule={rule!r} is not a rule. Use one of "
            + ", ".join(repr(name) for name in REMAP_RULES) + "."
        )
    if uncovered not in UNCOVERED_MODES:
        raise ValueError(
            f"uncovered={uncovered!r} is not a mode. Use one of "
            + ", ".join(repr(name) for name in UNCOVERED_MODES) + "."
        )


def remap_arrays(
    options: dict,
    *,
    old_polys,
    new_polys,
    cellwise: tuple[str, ...],
    package: str,
    rule: str = DEFAULT_RULE,
    uncovered: str = DEFAULT_UNCOVERED,
    domain=None,
    verb: str = "with_grid(remap=True)",
) -> tuple[dict, dict]:
    """Return ``(new options, report)`` with every per-cell array on the new mesh.

    ``cellwise`` names the options that hold one value per cell, as
    :func:`myflopy.grid_identity.cellwise_options` decides it -- this does not re-derive
    that, because a per-LAYER list such as ``k=[10.0, 1.0]`` means the same thing on any
    mesh and must not be touched.

    ``domain`` is the OLD model's idomain, and it masks the WEIGHTS per layer rather than
    the values. Measured on a pinching wedge, handing an inactive old cell's value
    through as ``NaN`` left ten reached new cells ``NaN``, three of them ACTIVE on the new
    mesh -- and ``reached`` says nothing about it.
    """

    check_remap_choices(rule, uncovered)

    refused = sorted(set(cellwise) & FLAG_ARRAYS)
    if refused:
        raise RemapRefused(
            f"{package}: "
            + ", ".join(repr(name) for name in refused)
            + (" is an integer flag" if len(refused) == 1 else " are integer flags")
            + " held one value per cell, and combining several old cells into one new "
            "cell averages them -- there is no cell type 0.4, and irch/ievt are "
            "zero-based LAYER INDICES. Re-declare "
            + ("it" if len(refused) == 1 else "them")
            + " on the new mesh, or name this package in `grid_free=` if the value is "
            "uniform."
        )
    unknown = sorted(set(cellwise) - INTENSIVE_ARRAYS - CONDUCTIVITY_ARRAYS)
    if unknown:
        raise RemapRefused(
            f"{package}: "
            + ", ".join(repr(name) for name in unknown)
            + " holds one value per cell and this engine has no measured rule for "
            + ("it" if len(unknown) == 1 else "them")
            + ". It carries intensive arrays ("
            + ", ".join(sorted(INTENSIVE_ARRAYS))
            + ") and conductivities ("
            + ", ".join(sorted(CONDUCTIVITY_ARRAYS))
            + "). An EXTENSIVE quantity -- a rate or a total -- is deliberately not "
            "carried: its denominator is the OLD cell's area, and a rate and a total "
            "need opposite treatment while declaring the same kind (ledger 300)."
        )

    # The frame guard, before any value is touched and before `uncovered=` is consulted.
    # It computes the overlay for its disjoint test and hands it back: computing it again
    # here doubled the cost of every remap.
    overlay = refuse_a_different_frame(old_polys, new_polys, verb=verb)

    idomain = None if domain is None else np.asarray(domain)
    moved: dict = dict(options)
    report: dict = {"filled": 0, "filled_area": 0.0, "distance": None, "spread": {}}
    new_area = new_polys.geometry.area.to_numpy()

    for name in cellwise:
        source = np.asarray(options[name], dtype=float)
        layers = source if source.ndim == 2 else source.reshape(1, -1)
        if name == "strt":
            # A starting head may hold MF6's no-value numbers -- one copied from a run's
            # output does -- and an area mean averaged a -1e30 into every new cell over it.
            # The rule `remap_field` uses: a dry cell takes the water table below it, an
            # inactive one is a gap, and a gap is masked from the weights like an
            # inactive cell.
            layers, _ = no_value_to_nan(layers, water_table=True)
        chosen = rule if name in CONDUCTIVITY_ARRAYS else "mean"
        carried = []
        for index, layer_values in enumerate(layers):
            active = _active_mask(idomain, index, overlay.n_old)
            if name == "strt" and not np.isfinite(layer_values).all():
                finite = np.isfinite(layer_values)
                active = finite if active is None else (active & finite)
            painted, reached = overlay.reduce(layer_values, rule=chosen, active=active)
            if name in CONDUCTIVITY_ARRAYS:
                _note_spread(overlay, layer_values, active, report, name, index)
            if not reached.all():
                painted = _fill(
                    painted, reached, layer_values, overlay, old_polys, new_polys,
                    active=active, uncovered=uncovered, package=package, name=name,
                    new_area=new_area, report=report, verb=verb,
                )
            carried.append(painted)
        stacked = np.array(carried)
        moved[name] = stacked if source.ndim == 2 else stacked[0]

    _warn(report, package=package, cellwise=cellwise, rule=rule)
    return moved, report


def _active_mask(idomain, layer: int, n_old: int):
    """The old mesh's active cells in ``layer``, or ``None`` when it declares none."""

    if idomain is None:
        return None
    if idomain.ndim == 1:
        return idomain > 0
    if layer < idomain.shape[0]:
        return idomain[layer] > 0
    # Fewer idomain layers than array layers: the array is authoritative about its own
    # depth, and guessing a mask for a layer the domain does not describe would mask
    # cells that are in fact active. Left unmasked rather than assumed inactive.
    return None


def _note_spread(overlay, values, active, report, name, layer) -> None:
    """Record how far the candidate rules are apart, for the warning to name."""

    try:
        ratio, worst = overlay.rule_spread(values, active=active)
    except ValueError:
        # A non-positive contribution: the spread is undefined and the chosen rule has
        # already raised if it needed positivity. Nothing to report rather than a crash
        # inside a diagnostic.
        return
    if worst is None or ratio.size == 0:
        return
    material = int(np.count_nonzero(ratio > _SPREAD_REPORTED_ABOVE))
    if material:
        report["spread"][f"{name}[{layer}]"] = (material, float(np.nanmax(ratio)))


def _fill(painted, reached, layer_values, overlay, old_polys, new_polys, *,
          active, uncovered, package, name, new_area, report, verb):
    """Fill the new cells with NO old cell under them, per ``uncovered``."""

    targets = np.flatnonzero(~reached)
    short_area = float(new_area[targets].sum())
    if uncovered == "raise":
        raise RemapRefused(
            f"{verb}: {package}.{name} has {targets.size} new cells with no old cell "
            f"under them at all ({short_area:,.0f} in model units of area). Pass "
            "`uncovered='nearest'` to fill them from the nearest old cell, or move to a "
            "mesh the old one covers."
        )
    filled, distance = _nearest_values(
        _centroids(old_polys), _centroids(new_polys, targets), layer_values, usable=active
    )
    if filled is None:
        raise RemapRefused(
            f"{verb}: {package}.{name} has {targets.size} uncovered new cells and NO "
            "usable old cell to fill them from -- every old cell is inactive in this "
            "layer or holds a non-finite value."
        )
    painted = painted.copy()
    painted[targets] = filled
    report["filled"] += int(targets.size)
    report["filled_area"] += short_area
    report["distance"] = max(report["distance"] or 0.0, distance or 0.0)
    return painted


def _warn(report, *, package: str, cellwise, rule: str) -> None:
    """Say what was interpolated and what was guessed -- in a channel that reaches a human.

    ``warnings.warn`` as well as the logger, because ``import myflopy`` installs a
    ``NullHandler`` on the package root and nothing configures logging, so
    ``logging.lastResort`` never fires and a log line reaches NOBODY -- measured at zero
    bytes on stderr (ledger 288). The existing nearest-neighbour fill for rasters reports
    through ``logger.info`` and is the anti-precedent.

    The PACKAGE and the parameters are named in the text because ``warnings.warn``
    de-duplicates on ``(text, category, module, lineno)``, and one move remaps several
    arrays from one line of code -- without them the second package's warning would be
    swallowed as a repeat of the first.
    """

    pieces = []
    if report["filled"]:
        pieces.append(
            f"{report['filled']} new cells had no old cell under them and were filled "
            f"from the nearest one, up to {report['distance']:,.0f} away, covering "
            f"{report['filled_area']:,.0f} in model units of area"
        )
    for key, (count, worst) in sorted(report["spread"].items()):
        pieces.append(
            f"{key}: the arithmetic and harmonic means differ by more than "
            f"{_SPREAD_REPORTED_ABOVE:g}x on {count} cells (worst {worst:,.1f}x), so "
            f"rule={rule!r} is a real choice there and not a formality"
        )
    if not pieces:
        return
    message = (
        f"{package}: remapped "
        + ", ".join(sorted(cellwise))
        + " onto the new mesh by area overlay. "
        + "; ".join(pieces)
        + ". Values were INTERPOLATED, not re-derived: where this package has a source "
        "(a polygon map, a GeoPackage, a raster) re-declaring it on the new mesh is "
        "measurably more faithful -- 7.1% worst-case flux error against this engine's "
        "12.5% at best (ledger 300)."
    )
    logger.warning("%s", message)
    warnings.warn(message, UserWarning, stacklevel=stacklevel_outside_package())


# -- a SOLVED field (heads, concentration, temperature) onto another mesh -----------------

#: How a solved field crosses between meshes. ``"linear"`` interpolates between the old
#: cells' MF6 centres; ``"area"`` is :meth:`CellOverlay.reduce`'s area-weighted mean --
#: the rule ``with_grid(remap=True)`` carries ``ic.strt`` by.
FIELD_METHODS: tuple[str, ...] = ("linear", "area")

#: ``"linear"``, measured (``docs/reproducers/heads_remap_linear_vs_area.py``): an MF6
#: head is a value AT its cell's centre, and the area mean paints every new cell inside
#: an old one with that old cell's value -- flat inside it, then a jump of a whole old
#: cell's head drop at its face. Against the other mesh's own solved heads, away from a
#: pumping well, linear was 6.4x closer refining and 2.6x closer coarsening; on a field
#: with a known exact answer, 57x and 150x.
DEFAULT_FIELD_METHOD = "linear"


def check_field_choices(method: str, uncovered: str) -> None:
    """Raise ``ValueError`` naming the choices when ``method``/``uncovered`` is not one.

    Run at the call, before anything is read, so a typo is refused even on a call that
    would have had nothing to fill.

    Parameters
    ----------
    method : str
        One of :data:`FIELD_METHODS`.
    uncovered : str
        One of :data:`UNCOVERED_MODES`.

    Raises
    ------
    ValueError
        Either is not one of its choices.
    """

    if method not in FIELD_METHODS:
        raise ValueError(
            f"method={method!r} is not a method. Use one of "
            + ", ".join(repr(name) for name in FIELD_METHODS) + "."
        )
    if uncovered not in UNCOVERED_MODES:
        raise ValueError(
            f"uncovered={uncovered!r} is not a mode. Use one of "
            + ", ".join(repr(name) for name in UNCOVERED_MODES) + "."
        )


def _field_grid(grid, *, verb: str, side: str):
    """``(polygons, centres, iverts)`` of a BUILT grid, or a ``TypeError`` naming one.

    A built grid is what ``GridSpec.resolve()`` returns, what ``ModelContext(grid=)``
    takes and what ``model.vor`` / ``mf.load_mf6_run(ws).vor`` are: it carries its cell
    polygons, MF6's cell centres (``points`` -- the ``xc``/``yc`` its DISV writes) and
    each cell's vertex ring. Everything else is refused naming the spelling that gives
    one, rather than guessed at.
    """

    if all(hasattr(grid, name) for name in ("gdf_vorPolys", "points", "iverts")):
        polys = grid.gdf_vorPolys
        centres = np.asarray(grid.points, dtype=float)
        if centres.ndim != 2 or centres.shape != (len(polys), 2):
            raise TypeError(
                f"{verb}: the {side} grid carries {len(polys)} cell polygons but its "
                f"centres (`points`) have shape {centres.shape}; a built grid holds one "
                "(x, y) centre per cell."
            )
        return polys, centres, grid.iverts
    kind = type(grid).__name__
    if kind == "GridRef":
        hint = (
            "a project grid REFERENCE; build it first, "
            "project.grid(key).resolve(workspace=...)"
        )
    elif kind == "GridSpec" or hasattr(grid, "resolve"):
        hint = (
            "a grid RECIPE; build it first, spec.resolve(workspace=...), and give the "
            "new model that same grid -- the values are carried onto one exact mesh"
        )
    elif hasattr(grid, "vor"):
        hint = "a model; pass its grid, model.vor"
    elif hasattr(grid, "geometry"):
        hint = (
            "a table of polygons, which holds no cell CENTRES -- the points MF6 solved "
            "the field at; pass the built grid the polygons came from"
        )
    elif hasattr(grid, "xcellcenters"):
        hint = (
            "a FloPy grid; pass myflopy's view of it, mf.load_mf6_run(ws).vor, or build "
            "mf.VoronoiGridPlus(verts=..., iverts=..., xcyc=..., crs=...)"
        )
    else:
        hint = "not a grid"
    raise TypeError(
        f"{verb}: the {side} grid must be a BUILT grid (a VoronoiGridPlus -- its cell "
        f"polygons, centres and vertex rings); got {kind}, {hint}."
    )


def _field_rows(layers, nlay: int, *, verb: str) -> tuple[list[int], bool]:
    """``layers`` as a list of OLD layer indices, one per output row, and whether it was one int."""

    def _one(value):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
            raise TypeError(
                f"{verb}: layers must be an old layer index or a sequence of them "
                f"(zero-based integers), got {value!r}."
            )
        index = int(value)
        if not -nlay <= index < nlay:
            raise ValueError(
                f"{verb}: layer {index} is outside the old model's {nlay} layers "
                f"(zero-based 0..{nlay - 1}, or a negative counting up from the bottom)."
            )
        return index % nlay

    if layers is None:
        return list(range(nlay)), False
    if isinstance(layers, (int, np.integer, bool, np.bool_)):
        return [_one(layers)], True
    try:
        rows = [_one(value) for value in layers]
    except TypeError as error:
        if isinstance(layers, (str, bytes)) or not hasattr(layers, "__iter__"):
            raise TypeError(
                f"{verb}: layers must be an old layer index or a sequence of them, got "
                f"{layers!r}."
            ) from error
        raise
    if not rows:
        raise ValueError(f"{verb}: layers is empty; name at least one old layer.")
    return rows, False


def no_value_to_nan(values, *, water_table: bool) -> tuple[np.ndarray, int]:
    """A raw ``(nlay, ncpl)`` field with MODFLOW's no-value markers made ``NaN``.

    An INACTIVE cell (+1e30) has no value. A DRY one (-1e30) has none of its own -- but
    with ``water_table`` it takes the value of the first layer BELOW it in its column
    that holds one: the water table. That normally sits below the dry cell's bottom, so
    a model started from it starts the cell dry again, as the solved run had it (MF6's IC
    documentation: a head below the cell bottom starts the cell dry). Not always: without
    rewetting MF6 keeps a cell dry after the water below has risen past its bottom, and
    then the number is ABOVE that bottom and the new model starts the cell wet -- where the
    water now is. Either way it keeps the field continuous across the wet/dry edge, where
    a gap would make the interpolation step. A dry cell with nothing below it keeps no
    value. Under NEWTON MF6 writes a dry cell's computed head as an ordinary number, so
    nothing is filled.

    Parameters
    ----------
    values : array_like
        ``(nlay, ncpl)`` numbers as the binary file holds them.
    water_table : bool
        Fill a dry cell from the layer below -- right for a HEAD; a concentration or a
        temperature at a dry cell is simply absent.

    Returns
    -------
    tuple of (numpy.ndarray, int)
        The field as floats with ``NaN`` where no value is left, and how many dry cells
        took the water table.
    """

    raw = np.asarray(values, dtype=float)
    field = np.where(is_no_value(raw), np.nan, raw)
    if not water_table:
        return field, 0
    dry = is_dry(raw)
    below = np.full(raw.shape[1], np.nan)
    filled = 0
    for layer in range(raw.shape[0] - 1, -1, -1):
        take = dry[layer] & np.isfinite(below)
        field[layer, take] = below[take]
        filled += int(np.count_nonzero(take))
        below = np.where(np.isfinite(field[layer]), field[layer], below)
    return field, filled


def _touching(iverts, ncpl: int):
    """A sparse ``(ncpl, ncpl)`` relation: True where two cells share at least one vertex."""

    import scipy.sparse as sparse

    rings = [np.asarray(ring, dtype=int) for ring in iverts]
    rows = np.repeat(np.arange(ncpl), [ring.size for ring in rings])
    columns = np.concatenate(rings) if rings else np.array([], dtype=int)
    incidence = sparse.csr_matrix(
        (np.ones(rows.size, dtype=np.int32), (rows, columns)),
        shape=(ncpl, int(columns.max()) + 1 if columns.size else 0),
    )
    return (incidence @ incidence.T).tocsr()


class _Triangulation:
    """The linear reconstruction of one layer's field from its VALID old centres.

    Built once per distinct validity mask and reused by every layer that shares it --
    re-triangulating per layer measured 7x the cost (0.053 s against 0.0075 s at 20,537
    cells). Two parts:

    * **Inside**: a Delaunay triangulation of the valid old centres, interpolated
      barycentrically at each new centre lying in a LOCAL simplex -- one whose every edge
      joins two old cells sharing a vertex. A triangle failing that reaches past the
      cells' own neighbourhood: across a notch, a hole or an inactive strip, or along a
      sliver of the hull whose corners are not neighbours. Measured across an inactive
      strip between fields of 10 and 20, 6 of 36 new cells read 17.2-18.0 through such
      triangles (``docs/reproducers/field_remap_rim.py``). On a lattice the test flags
      nothing (0 of 12,482 simplices); a face-sharing test would have flagged its
      co-circular diagonals.
    * **Everywhere else** (see :meth:`extend`): a valid old cell's value moved along that
      cell's own least-squares gradient.
    """

    def __init__(self, old_xy, ids, new_xy, touch):
        from scipy.spatial import Delaunay, QhullError

        n_new = len(new_xy)
        self.ids = ids
        self.usable = np.zeros(n_new, dtype=bool)
        self.outside_hull = np.ones(n_new, dtype=bool)
        self.vertices = None
        self.weights = None
        self.pairs = None
        if ids.size == 0:
            return
        # Centred: Qhull and the plane fits work in the frame they are given, and State
        # Plane coordinates (~1.3e6) leave fewer digits than a local frame does.
        self.origin = old_xy[ids].mean(axis=0)
        self.local_xy = old_xy[ids] - self.origin
        self._plane_sums(ids, touch)
        if ids.size < 3:
            return
        try:
            tri = Delaunay(self.local_xy)
        except QhullError:
            # Three or more points on one line: no triangle to interpolate in, so every
            # reached cell of this layer takes the extension.
            logger.debug("valid old centres are collinear; no linear interpolation here")
            return
        local = np.ones(len(tri.simplices), dtype=bool)
        if touch is not None:
            corners = ids[tri.simplices]
            for a, b in ((0, 1), (1, 2), (0, 2)):
                local &= np.asarray(touch[corners[:, a], corners[:, b]]).ravel() > 0
        points = np.asarray(new_xy, dtype=float) - self.origin
        simplex = tri.find_simplex(points)
        found = simplex >= 0
        at = np.where(found, simplex, 0)
        self.outside_hull = ~found
        self.usable = found & local[at]
        transform = tri.transform[at]
        partial = np.einsum("nij,nj->ni", transform[:, :2], points - transform[:, 2])
        self.weights = np.column_stack((partial, 1.0 - partial.sum(axis=1)))
        self.vertices = tri.simplices[at]

    def _plane_sums(self, ids, touch) -> None:
        """The geometric half of each valid cell's least-squares gradient, shared by layers.

        The gradient solves the 2x2 normal equations of a plane through the cell and its
        valid vertex-touching neighbours. Where those neighbours lie close to ONE line
        through the cell, the across-line component is set by their jitter, not by the
        field; so the matrix is split into its two directions and a direction carrying
        less than :data:`_PLANE_SPREAD` of the other's spread is dropped -- the gradient is
        kept along the line the neighbours DO span, and is zero across it.
        """

        if touch is None:
            return
        position = np.full(touch.shape[0], -1)
        position[ids] = np.arange(ids.size)
        coo = touch[ids].tocoo()
        here = coo.row
        there = position[coo.col]
        keep = (there >= 0) & (there != here)
        here, there = here[keep], there[keep]
        dx = self.local_xy[there, 0] - self.local_xy[here, 0]
        dy = self.local_xy[there, 1] - self.local_xy[here, 1]
        n = ids.size
        sxx = np.bincount(here, weights=dx * dx, minlength=n)
        sxy = np.bincount(here, weights=dx * dy, minlength=n)
        syy = np.bincount(here, weights=dy * dy, minlength=n)
        half_gap = np.sqrt(((sxx - syy) / 2.0) ** 2 + sxy**2)
        big = (sxx + syy) / 2.0 + half_gap
        small = (sxx + syy) / 2.0 - half_gap
        # The unit direction of the larger spread. With no cross term the axes are the
        # directions; otherwise (big - syy, sxy) is its eigenvector.
        ex = np.where(sxy != 0.0, big - syy, np.where(sxx >= syy, 1.0, 0.0))
        ey = np.where(sxy != 0.0, sxy, np.where(sxx >= syy, 0.0, 1.0))
        norm = np.hypot(ex, ey)
        norm = np.where(norm > 0.0, norm, 1.0)
        full = small > _PLANE_SPREAD * big
        self.pairs = (here, there, dx, dy, sxx, sxy, syy, big, ex / norm, ey / norm, full)

    def interpolate(self, values_at_ids) -> np.ndarray:
        """The linear field at every new centre (``NaN`` where it is not usable)."""

        out = np.full(self.usable.size, np.nan)
        if self.weights is None:
            return out
        corners = np.asarray(values_at_ids, dtype=float)[self.vertices]
        out[self.usable] = np.einsum("ni,ni->n", corners, self.weights)[self.usable]
        return out

    def extend(self, values_at_ids, anchors, points) -> np.ndarray:
        """Each ``anchors`` cell's value moved along its own gradient to ``points``.

        ``anchors`` are positions in :attr:`ids` -- the valid old cell holding the largest
        share of each new cell, so the value comes from the cell UNDER it, never from one
        across the gap that made it unusable -- and ``points`` the new cells' centres. That
        serves a new centre outside the hull of the old ones (a band half a cell wide
        around a structured or jittered mesh, whose boundary cells are not centred ON the
        boundary: 38 of 144 cells on one 12x12 test mesh) and one in a non-local triangle.
        It is exact for a planar field wherever the plane stays inside the range the old
        layer holds -- the area mean it replaced was off by 2.79 on that mesh -- and is
        CLAMPED to that range, so it never invents a value the old field does not hold: a
        plane fitted across a sharp front extrapolates past it (measured, -0.055 from a
        field whose least value was 4.7e-9). The layer's range and not the fitted cells':
        at the rim a plane lies past its neighbours by definition, and clamping to them
        cost as much as the area mean (2.782 against 2.786) where the layer clamp clips
        only the extreme corners -- never the worse of the two on four fields.
        """

        values = np.asarray(values_at_ids, dtype=float)
        return np.clip(self.moved(values, anchors, points), values.min(), values.max())

    def moved(self, values_at_ids, anchors, points) -> np.ndarray:
        """:meth:`extend` before its clamp: the anchors' values along their gradients."""

        values = np.asarray(values_at_ids, dtype=float)
        out = values[anchors]
        if self.pairs is None:
            return out
        here, there, dx, dy, sxx, sxy, syy, big, ex, ey, full = self.pairs
        dh = values[there] - values[here]
        n = values.size
        sxh = np.bincount(here, weights=dx * dh, minlength=n)
        syh = np.bincount(here, weights=dy * dh, minlength=n)
        determinant = sxx * syy - sxy * sxy
        # `_PLANE_SPREAD` keeps a counted direction's determinant clear of zero; the
        # second test makes the division unable to see one whatever that constant is.
        full = full & (determinant > 0.0)
        safe = np.where(full, determinant, 1.0)
        gx = np.where(full, (syy * sxh - sxy * syh) / safe, 0.0)
        gy = np.where(full, (sxx * syh - sxy * sxh) / safe, 0.0)
        along = np.where(big > 0.0, (ex * sxh + ey * syh) / np.where(big > 0.0, big, 1.0), 0.0)
        gx = np.where(full, gx, along * ex)
        gy = np.where(full, gy, along * ey)
        offset = np.asarray(points, dtype=float) - self.origin - self.local_xy[anchors]
        return out + gx[anchors] * offset[:, 0] + gy[anchors] * offset[:, 1]


#: The least spread, as a fraction of the larger, that a cell's neighbours must have
#: ACROSS their main line for the gradient across it to be fitted: below it they lie too
#: close to one line, and the across-line slope would be their jitter, not the field.
#: 0.01 is a cross-line spread of a tenth of the along-line one.
_PLANE_SPREAD = 0.01


def _largest_valid_piece(overlay, ok) -> np.ndarray:
    """For each new cell, the VALID old cell covering most of it (``-1`` where none).

    Ties go to the lowest old cell, the tie rule :meth:`CellOverlay.reduce`'s majority
    uses.
    """

    keep = ok[overlay.old_index]
    old_at = overlay.old_index[keep]
    new_at = overlay.new_index[keep]
    area = overlay.area[keep]
    largest = np.full(overlay.n_new, -1)
    if area.size == 0:
        return largest
    order = np.lexsort((old_at, -area, new_at))
    leading = np.empty(order.size, dtype=bool)
    leading[0] = True
    leading[1:] = new_at[order][1:] != new_at[order][:-1]
    largest[new_at[order][leading]] = old_at[order][leading]
    return largest


def _refuse_overlapping_cells(overlay, new_area, *, verb: str) -> None:
    """Refuse an OLD grid whose cells overlap one another -- a DISU's stacked nodes.

    A plan-view mesh covers each point once, so the old pieces under a new cell can add
    up to its area at most. A DISU head file holds every node as one layer, and its view
    stacks a column's nodes on one polygon; carried as a layer, a new cell would average
    and triangulate nodes from DIFFERENT layers (measured: a 2-layer DISU read through a
    reader opened from its path returned a 1-row field without complaint). Checked on the
    overlay already computed, so it costs nothing.
    """

    covered = np.bincount(overlay.new_index, weights=overlay.area, minlength=overlay.n_new)
    over = covered > np.asarray(new_area, dtype=float) * (1.0 + 1e-6) + 1e-9
    if over.any():
        raise RemapRefused(
            f"{verb}: the old grid's cells OVERLAP one another -- {int(over.sum())} new "
            "cells lie under more than their own area of old cells -- so it is not a "
            "plan-view mesh: a DISU, whose head file holds every node as one layer, or a "
            "grid stacking several layers' cells. Carry a DIS or DISV model's field."
        )


def remap_field(
    values,
    old_grid,
    new_grid,
    *,
    layers=None,
    method: str = DEFAULT_FIELD_METHOD,
    uncovered: str = DEFAULT_UNCOVERED,
    water_table: bool = True,
    active=None,
    label: str = "head",
    verb: str = "model.hds.remap",
    when: str | None = None,
) -> tuple[np.ndarray, dict]:
    """Carry one saved time of a SOLVED field from ``old_grid`` onto ``new_grid``.

    The engine behind ``model.hds.remap`` (and ``model.conc``/``model.temp``) and the
    lazy recipe: it takes the field's own numbers at one time and the two built grids,
    and needs no model object. The caller picks the time.

    Each output row reads one OLD layer (``layers``). Per row: the no-value markers are
    made gaps (:func:`no_value_to_nan`, a DRY head taking the water table below it); the
    old cells still holding a value are overlaid on the new mesh ONCE, through the frame
    guard (:func:`~myflopy.modflow.mf6.grid.geometry.refuse_a_different_frame`), and a new
    cell with such a cell under it is REACHED. ``"linear"`` interpolates a reached cell's
    value at its MF6 centre from the old centres around it; where the centre falls outside
    those centres' hull, or in a triangle reaching past their neighbourhood
    (:class:`_Triangulation`), it takes the value of the valid old cell covering most of
    it, moved along that cell's own least-squares gradient and held inside the range the
    old layer holds: local, exact for a plane inside that range, and never a value the old
    field does not hold. A cell NOT reached -- outside the old model, or over cells holding
    no value -- is ``uncovered``: the nearest valid old centre's value with ONE warning, or
    a refusal. The result holds no ``NaN``: MF6 reads a ``NaN`` starting head, reports
    normal termination, and returns a ``NaN`` budget (ledger 247).

    Parameters
    ----------
    values : array_like
        ``(nlay, ncpl)`` numbers as the binary file holds them, ``ncpl`` the old grid's
        cell count, MODFLOW's no-value markers (+1e30 inactive, -1e30 dry) included; a
        ``NaN`` is a gap too.
    old_grid : VoronoiGridPlus
        The BUILT grid ``values`` were solved on: cell polygons (``gdf_vorPolys``), MF6's
        cell centres (``points``) and vertex rings (``iverts``). One whose cells overlap
        one another -- a DISU's stacked nodes -- is refused.
    new_grid : VoronoiGridPlus
        The built grid to carry the field onto, in the SAME CRS. A grid recipe, a grid
        reference, a model, a FloPy grid or a table of polygons is refused with a
        ``TypeError`` naming the spelling that gives a built grid.
    layers : int or sequence of int, optional
        The OLD layer each output row reads, zero-based, a negative counting up from the
        bottom; repeats are allowed, which is how a new model with MORE layers than the old
        one is filled (``[0, 0, 1, 2]``). An int returns one ``(new_ncpl,)`` row. Default
        ``None``: every old layer, in order.
    method : {"linear", "area"}, default "linear"
        ``"linear"`` interpolates between old cell centres; ``"area"`` takes the
        area-weighted mean of the valid old cells under each new cell -- the rule
        ``with_grid(remap=True)`` carries ``ic.strt`` by. See :data:`DEFAULT_FIELD_METHOD`
        for the measurement behind the default.
    uncovered : {"nearest", "raise"}, default "nearest"
        A new cell no valid old cell reaches: ``"nearest"`` takes the value of the nearest
        old centre holding one in that layer and warns once, naming the counts and the
        farthest distance; ``"raise"`` refuses with ``RemapRefused``.
    water_table : bool, default True
        A DRY old cell takes the value of the first layer below it holding one (heads);
        ``False`` leaves it a gap (concentration, temperature).
    active : array_like of bool, optional
        ``(len(layers), new_ncpl)``: which new cells each output row's model layer keeps
        (``idomain > 0``). An uncovered cell inactive in EVERY row reading its old layer is
        still filled -- the array stays finite -- but is neither counted, warned about nor
        refused, because the new model never reads it. Default ``None``: every new cell
        counts.
    label : str, default "head"
        What the field is, for the warning and refusals (``"concentration"``).
    verb : str, default "model.hds.remap"
        The spelling the caller was reached through, for the warning and refusals.
    when : str, optional
        Which saved time ``values`` is, for the warning (``"the end of stress period
        3"``). Default ``None``: not named.

    Returns
    -------
    tuple of (numpy.ndarray, dict)
        The field on ``new_grid`` -- ``(len(layers), new_ncpl)``, or ``(new_ncpl,)`` for an
        int ``layers`` -- finite everywhere; and a report: ``rows`` (the old layer of each
        row), ``dry_filled`` (dry cells of the layers read that took the water table),
        ``outside`` (counted new cells with no old cell under them at
        all) and ``outside_area``, ``no_value_under`` (per OLD layer, counted cells over
        only cells holding no value), ``distance`` (the farthest counted nearest-cell
        fill), ``rim`` and ``bridged`` (CELL-LAYERS -- summed over the distinct old layers
        read -- that ``"linear"`` could not interpolate and extended instead: ``rim`` outside
        the hull of the valid centres, which is every one when they cannot be triangulated,
        ``bridged`` in a non-local triangle).

    Raises
    ------
    TypeError
        ``old_grid``/``new_grid`` is not a built grid, or ``layers`` holds a non-integer.
    ValueError
        ``method``/``uncovered`` is not one of its choices; ``values`` is not one value
        per old cell and layer; ``active`` is not one row per output row; a layer is out
        of range; the two grids are in different CRSs or share no area.
        ``RemapRefused`` (a ``ValueError``) when the old grid's cells overlap, under
        ``uncovered="raise"`` when a counted cell is not reached, or when a requested old
        layer holds no value at all.

    Examples
    --------
    >>> strt, report = remap_field(raw, old.vor, new_vor)          # every layer
    >>> top, _ = remap_field(raw, old.vor, new_vor, layers=0, method="area")
    """

    check_field_choices(method, uncovered)
    old_polys, old_xy, old_iverts = _field_grid(old_grid, verb=verb, side="old")
    new_polys, new_xy, _ = _field_grid(new_grid, verb=verb, side="new")
    for side, centres in (("old", old_xy), ("new", new_xy)):
        if len(np.unique(centres, axis=0)) < len(centres):
            raise RemapRefused(
                f"{verb}: cells of the {side} grid share a centre, so it is not a plan-view "
                "mesh: a DISU, whose nodes stack a column's layers on one polygon, or a "
                "grid repeating its cells per layer. Carry between DIS or DISV grids."
            )
    raw = np.asarray(values, dtype=float)
    if raw.ndim != 2 or raw.shape[1] != len(old_polys):
        raise ValueError(
            f"{verb}: the field must hold one value per old cell and layer, "
            f"(nlay, {len(old_polys)}); got shape {raw.shape}."
        )
    rows, single = _field_rows(layers, raw.shape[0], verb=verb)
    counts = None
    if active is not None:
        counts = np.asarray(active, dtype=bool)
        if counts.ndim == 1:
            counts = counts[None, :]
        if counts.shape != (len(rows), len(new_polys)):
            raise ValueError(
                f"{verb}: active must say which new cells each output row keeps, "
                f"({len(rows)}, {len(new_polys)}); got shape {counts.shape}."
            )
    field, _ = no_value_to_nan(raw, water_table=water_table)
    valid = np.isfinite(field)
    read = sorted(set(rows))
    # Dry cells that took the water table, in the layers READ: the rest were filled only
    # so a dry cell above one of these could find its water.
    dry_filled = int(np.count_nonzero((is_dry(raw) & valid)[read])) if water_table else 0

    overlay = refuse_a_different_frame(old_polys, new_polys, verb=verb)
    new_area = new_polys.geometry.area.to_numpy()
    _refuse_overlapping_cells(overlay, new_area, verb=verb)
    anything_under = np.bincount(overlay.new_index, minlength=overlay.n_new) > 0
    touch = _touching(old_iverts, len(old_xy)) if method == "linear" else None

    triangulations: dict[bytes, _Triangulation] = {}
    painted_rows: dict[int, np.ndarray] = {}
    gaps_by_layer: dict[int, np.ndarray] = {}
    report = {
        "rows": list(rows), "dry_filled": dry_filled, "outside": 0, "outside_area": 0.0,
        "no_value_under": {}, "distance": None, "rim": 0, "bridged": 0,
    }
    for layer in sorted(set(rows)):
        ok = valid[layer]
        painted, reached = overlay.reduce(
            np.where(ok, field[layer], 0.0), rule="mean", active=ok
        )
        if method == "linear" and ok.any():
            key = ok.tobytes()
            if key not in triangulations:
                triangulations[key] = _Triangulation(old_xy, np.flatnonzero(ok), new_xy, touch)
            linear = triangulations[key]
            at_ids = field[layer, linear.ids]
            interpolated = linear.interpolate(at_ids)
            take = reached & np.isfinite(interpolated)
            painted = np.where(take, interpolated, painted)
            handed_back = np.flatnonzero(reached & ~take)
            if handed_back.size:
                position = np.full(len(old_xy), -1)
                position[linear.ids] = np.arange(linear.ids.size)
                anchors = position[_largest_valid_piece(overlay, ok)[handed_back]]
                painted[handed_back] = linear.extend(at_ids, anchors, new_xy[handed_back])
                rim = linear.outside_hull[handed_back]
                report["rim"] += int(np.count_nonzero(rim))
                report["bridged"] += int(np.count_nonzero(~rim))
        painted_rows[layer] = painted
        gaps_by_layer[layer] = ~reached

    # Which unreached cells COUNT: those some row reading that old layer keeps active.
    counted: dict[int, np.ndarray] = {}
    for layer, gaps in gaps_by_layer.items():
        if counts is None:
            counted[layer] = gaps
        else:
            reading = [row for row, old in enumerate(rows) if old == layer]
            counted[layer] = gaps & counts[reading].any(axis=0)
    asked = [layer for layer in sorted(gaps_by_layer) if counted[layer].any()]
    if asked and uncovered == "raise":
        outside = int(np.count_nonzero(np.logical_or.reduce(
            [counted[layer] & ~anything_under for layer in asked]
        )))
        detail = "; ".join(
            f"old layer {layer}: {int(np.count_nonzero(counted[layer]))}" for layer in asked
        )
        raise RemapRefused(
            f"{verb}: new cells with no {label} under them -- {detail} -- of which "
            f"{outside} lie outside the old model's cells altogether. Pass "
            "uncovered='nearest' to fill them from the nearest old cell holding one."
        )
    outside_cells = None
    for layer in sorted(gaps_by_layer):
        gaps = gaps_by_layer[layer]
        if not gaps.any():
            continue
        filled, distance = _nearest_values(
            old_xy, new_xy[gaps], field[layer], usable=valid[layer]
        )
        if filled is None:
            raise RemapRefused(
                f"{verb}: old layer {layer} holds no {label} at all -- every cell is "
                + ("inactive, or dry with no water below it" if water_table
                   else "inactive or dry")
                + " -- so there is nothing to carry into the rows reading it. Read "
                "another layer there: layers=[...] names the old layer of each row."
            )
        painted_rows[layer] = painted_rows[layer].copy()
        painted_rows[layer][gaps] = filled
        mine = counted[layer]
        if not mine.any():
            continue
        # The distance of the COUNTED fills only: the warning names it beside their count.
        _, counted_distance = _nearest_values(
            old_xy, new_xy[mine], field[layer], usable=valid[layer]
        )
        report["distance"] = max(report["distance"] or 0.0, counted_distance)
        report["no_value_under"][layer] = int(np.count_nonzero(mine & anything_under))
        here = mine & ~anything_under
        outside_cells = here if outside_cells is None else (outside_cells | here)
    if outside_cells is not None:
        report["outside"] = int(np.count_nonzero(outside_cells))
        report["outside_area"] = float(new_area[outside_cells].sum())

    out = np.vstack([painted_rows[layer] for layer in rows])
    if not np.isfinite(out).all():
        # A gap the passes above did not close is a defect here, not a property of the
        # data -- refuse rather than hand MF6 a NaN it would run on.
        raise RemapRefused(
            f"{verb}: {int(np.count_nonzero(~np.isfinite(out)))} values came back "
            "without a number; this is a defect in the remap, please report it."
        )
    if dry_filled:
        logger.info(
            "%s: %d dry cells took the %s of the first wet layer below them", verb,
            dry_filled, label,
        )
    if report["rim"] or report["bridged"]:
        logger.debug(
            "%s: %d cell-layers outside the hull of the old centres and %d in triangles "
            "reaching past their neighbourhood took the covering old cell's value along "
            "its own gradient", verb, report["rim"], report["bridged"],
        )
    if asked:
        _warn_field(report, verb=verb, label=label, when=when, method=method,
                    n_new=overlay.n_new, water_table=water_table)
    return (out[0] if single else out), report


def _warn_field(report, *, verb, label, when, method, n_new, water_table) -> None:
    """ONE warning naming what was GUESSED, through both channels (ledger 288).

    The time, the field and the counts are in the text: ``warnings.warn`` de-duplicates on
    the text and the calling line, and a re-run notebook cell is that same line.
    """

    pieces = []
    if report["outside"]:
        pieces.append(
            f"{report['outside']:,} new cells ({report['outside_area']:,.0f} in model "
            f"units of area) lie outside the old model's cells"
        )
    under = {layer: count for layer, count in report["no_value_under"].items() if count}
    if under:
        listing = ", ".join(f"{count:,} in old layer {layer}" for layer, count in under.items())
        why = "inactive, or dry with no water below" if water_table else "inactive or dry"
        pieces.append(f"{listing} sat only over cells holding no {label} ({why})")
    if not pieces:
        return
    at = f" at {when}" if when else ""
    how = "linear interpolation" if method == "linear" else "area-weighted mean"
    message = (
        f"{verb}: carried the {label}{at} (old layers {report['rows']}) onto {n_new:,} new "
        f"cells by {how}. "
        + "; ".join(pieces)
        + f" -- those took the {label} of the nearest old cell holding one in their "
        f"layer, up to {report['distance']:,.0f} away. Pass uncovered='raise' to refuse "
        "instead."
    )
    logger.warning("%s", message)
    warnings.warn(message, UserWarning, stacklevel=stacklevel_outside_package())

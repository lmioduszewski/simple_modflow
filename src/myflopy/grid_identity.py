"""Is this the same mesh, and which mesh was it?

Extracted from ``specs.py`` unchanged. It landed there because ``ModelSpec.build`` is
what asks the question, but it answers a question about GRIDS, not about specs: it
reads a FloPy package and a grid object and hashes connectivity. It referenced three
names from the rest of that module, two of which it owned -- so it was a cut, not an
untangling, and ``specs.py`` had grown 18.8% in four days against its own docstring
promising "intentionally small" spec classes.

Two fingerprints, on purpose, answering the same question from opposite sides:
:func:`grid_pin` reads a grid, :func:`disv_pin` reads a BUILT DISV package. They agree
on the same mesh, which is what lets a model be checked once its grid object is gone.
Both hash CONNECTIVITY only -- the cell count and each cell's canonical vertex ring --
so they are exact and need no tolerance, and are blind to a change that moves
coordinates without renumbering. The geometry comparison inside
:func:`_discretization_complaints` is what covers that, with a tolerance scaled to the
cell size, and the two are complementary rather than redundant.

Everything here is re-exported from ``myflopy.specs`` for the callers that already
import it there.
"""

from __future__ import annotations

import hashlib
from typing import Any

from myflopy._logging import get_logger

logger = get_logger(__name__)


#: Where a resolved package records the mesh its cell numbers refer to. It rides in
#: ``metadata`` rather than in a field of its own because it is provenance -- it says
#: where the records came from, not what they build -- which also keeps it out of
#: ``package_digest`` (the records already encode the grid, so digesting the pin as well
#: would report every stored run as drifted for no change in what it builds).
GRID_PIN = "grid_pin"

#: The option names a pin SPEAKS FOR, declared by the resolver that took it. An edit
#: that replaces one of these voids the pin; an edit that renames the file does not.
#: Declared at the pin site rather than guessed from a global list, because the same
#: option name is grid-bound on one package and not on another -- `packagedata` carries
#: a cellid on `sfr`/`uzf` and none on `lak` (ledger 275).
GRID_PIN_OPTIONS = "grid_pin_options"

#: What replaced a pin when an edit voided it: ``{"was": <pin>, "options": [names]}``.
#: The voided pin's VALUE is kept deliberately -- the build-time refusal reads it, so
#: voiding records what the records used to certify rather than deleting the evidence.
#: Measured: popping `GRID_PIN` instead removed the refusal outright, and a CHD pinned
#: to one mesh then built clean on another.
PIN_VOIDED = "grid_pin_voided"


#: Context fields that ``ModelSpec.to_dict`` cannot carry, and that something reads.
#: ``dates`` is deliberately absent: nothing in the library reads it, so losing it
#: costs nothing and warning about it would be noise.
DROPPED_CONTEXT_FIELDS = ("grid", "domain", "surfaces")


def _as_int(value):
    """``value`` as an int, or ``None`` when it is not one number.

    Every dimension this guard reads can legitimately be absent or non-scalar:
    FloPy computes ``NVERT`` at write time and ``nvert.get_data()`` is ``None``
    until then, and a FloPy ``UnstructuredGrid`` reports ``ncpl`` PER LAYER. An
    ``int()`` on either raises, which would make a guard the error it exists to
    prevent.
    """

    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None

def _canonical_ring(vertices) -> tuple[int, ...]:
    """One form for a cell's vertex ring, whichever way it was written down.

    ``[0, 1, 2, 3]``, ``[1, 2, 3, 0]`` and ``[0, 3, 2, 1]`` are the same polygon --
    rotated, and wound the other way. Hashing them raw made the fingerprint mean
    "this serialization of this mesh" rather than "this mesh", so two identical
    meshes could disagree and a correct model was refused after a reload. Rotating to
    the lowest index and taking the lesser of the two windings gives one answer.
    """

    ring = [int(part) for part in vertices]
    if not ring:
        return ()
    start = ring.index(min(ring))
    forward = tuple(ring[start:] + ring[:start])
    backward = (forward[0], *reversed(forward[1:]))
    return min(forward, backward)

def grid_pin(grid) -> str | None:
    """A fingerprint of a mesh's cell NUMBERING, or ``None`` if it cannot be taken.

    Hashes the connectivity -- the cell count and each cell's own canonical vertex
    ring, in cell order -- and nothing continuous. It does NOT hash the vertex count:
    FloPy computes ``NVERT`` at write time when a caller omits it, so folding it in
    made a correct package fingerprint differently from its own grid. That is what makes it
    exact where comparing coordinates needs a tolerance: the same mesh written to
    fewer decimal places, rendered float32, or shifted to a local origin all hash the
    same, while a renumbering does not. A renumbering is the thing that matters,
    because cell numbering is the whole meaning of a cellid.

    It deliberately says nothing about vertex POSITIONS, so a mesh stretched without
    being renumbered pins the same -- its cellids still name the same cells, and
    `_discretization_complaints`' geometry test is what notices the stretch.
    """

    if grid is None:
        return None
    try:
        props = grid.get_disv_gridprops()
        cell2d = props["cell2d"]
        # `nvert` is deliberately NOT hashed: FloPy computes it at write time when a
        # caller leaves it out, so `disv.nvert.get_data()` is None on a perfectly
        # correct package -- and folding it in made that model fingerprint differently
        # from its own grid. The vertex indices already carry the connectivity.
        payload = repr(
            (
                _as_int(props.get("ncpl")) or len(cell2d),
                tuple(_canonical_ring(row[4:]) for row in cell2d),
            )
        )
    except Exception:  # noqa: BLE001 - see below; the tuple cannot be closed
        # Deliberately broad, and it only ever weakens a guard. `grid` is whatever
        # object a caller put on a context -- a user's builder result, a grid
        # rebuilt from a run folder with no triangulation -- so reading its
        # connectivity reaches flopy, numpy and any property it computes lazily. A
        # pin that cannot be taken means "unpinned", never an error at the caller.
        logger.debug("could not fingerprint a grid for pinning", exc_info=True)
        return None
    return hashlib.sha256(payload.encode()).hexdigest()[:16]

def disv_pin(disv) -> str | None:
    """:func:`grid_pin`'s fingerprint, taken from a BUILT DISV package.

    The same connectivity hash, read off the package instead of the grid, so a model
    whose context no longer holds a grid object can still be checked against the mesh
    it was declared on. FloPy pads every ``cell2d`` row out to the widest cell, so the
    trailing ``None`` vertices have to go: ``ncvert`` says how many of the ``icvert``
    fields are real, and including the padding would make the answer depend on the
    widest cell in the grid rather than on this cell.
    """

    if disv is None:
        return None
    try:
        rows = []
        for row in disv.cell2d.get_data():
            count = int(row[3])
            rows.append(_canonical_ring(row[4 + i] for i in range(count)))
        payload = repr((_as_int(disv.ncpl.get_data()), tuple(rows)))
    except Exception:  # noqa: BLE001 - see below; the tuple cannot be closed
        # Deliberately broad, and it only ever weakens a guard. This reads a FloPy
        # package's data, which reaches flopy's own exception hierarchy and numpy; a
        # fingerprint that cannot be taken means "unpinned", never an error at a
        # caller who asked for a build.
        logger.debug("could not fingerprint a built DISV", exc_info=True)
        return None
    return hashlib.sha256(payload.encode()).hexdigest()[:16]

#: Metadata key recording a 16-character hash of the file a package was built FROM,
#: taken when it was declared. Excluded from `package_digest` like the rest of
#: `metadata`, so adding it moved no stored digest (ledger 276).
SOURCE_SHA = "source_sha"


def source_fingerprint(path) -> str | None:
    """A 16-character content hash of ``path``, or ``None`` if it cannot be read.

    ``None`` for a missing or unreadable file rather than raising: declaring a package
    must not fail because a bookkeeping column could not be filled, which is the same
    rule `package_digest` follows.

    Chunked, because a source is an arbitrary file on the user's disk. Measured on the
    real Ten Trails inputs: 3.2 ms for the largest GeoPackage (2.4 MB), against a
    declaration that already reads and spatially resolves the same file.

    **It hashes BYTES, not logical content.** Rewriting a GeoPackage without changing a
    feature still moves the hash -- SQLite rewrites pages, and the file carries its own
    timestamps. That is the safe direction: a false "edited" is a prompt to look, where
    a false "unchanged" is the lie this exists to remove.
    """

    import hashlib

    try:
        with open(path, "rb") as handle:
            digest = hashlib.sha256()
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                digest.update(chunk)
    except OSError:
        return None
    return digest.hexdigest()[:16]


def pin_metadata(context, *, certifies) -> dict[str, str]:
    """``{GRID_PIN: ..., GRID_PIN_OPTIONS: ...}`` for a pinnable grid, else ``{}``.

    ``certifies`` names the options this pin speaks for, and is REQUIRED -- every pin
    site must say what it resolved. Measured why it is not defaulted: an empty tuple
    writes an empty list, and ``void_pin`` intersects against it, so a site that
    forgot the keyword would be silently UN-VOIDABLE forever. An empty ``certifies``
    therefore omits the key entirely and lands in the conservative branch (void on any
    override), which is also what every pin serialized before ledger 275 does.

    The one place a builder stamps the mesh its records refer to, so the seven
    ``.gpkg`` forms and the advanced builders cannot drift apart about what a pin is.
    Spread ``**pin_metadata(self.context)`` into an existing ``with_metadata`` call.

    ``mf.hfb`` DOES carry one, since ledger 267. An earlier version of this docstring
    said it could not -- "face-indexed, no cellid column" -- which contradicted
    :func:`_cellids_in` eight lines below, whose own docstring says it reads "two for
    an HFB, which sits on the face between a pair". Measured: the range check finds 10
    cellids in a 5-segment barrier. What the range check cannot do is catch a barrier
    resolved on a SMALLER mesh, where every one of those numbers is still in range --
    which is the pin's job, not a reason to skip it.

    ``mf.mvr`` still carries none: it names packages and record indices rather than
    cells, and holds no context at all. It needs the symbolic form the plan describes,
    which is a different piece of work.
    """

    pin = grid_pin(getattr(context, "grid", None)) if context is not None else None
    if pin is None:
        return {}
    stamped: dict[str, str] = {GRID_PIN: pin}
    if certifies:
        stamped[GRID_PIN_OPTIONS] = sorted(certifies)
    return stamped


def void_pin(metadata, touched):
    """Return ``metadata`` with the pin voided if ``touched`` names an option it certifies.

    An edit that REPLACES a resolved option value makes the pin a lie: the records are
    now the caller's, and the package still certifies that the library resolved them on
    a named mesh. Measured before this -- an edited CHD kept pin
    ``639cf118e33c3df2`` and reported ``stale``, so a user who typed CORRECT cell
    numbers for the new grid was refused forever by a sentence about records that no
    longer existed.

    Keyed on the option NAME, not on a change of value: a value oracle means a deep
    walk of every option at every check, and ``fork_package`` already measures 808 ms
    copying a 654,480-record EVT. The cost is that re-supplying an identical dict voids
    a pin that was still true -- a spurious ``unverifiable`` rather than a silent stale
    pin, which is the direction that fails safe.

    Returns the SAME object when nothing voids: ``with_options`` shares its metadata
    dict with the parent, so this must never mutate in place.
    """

    pin = (metadata or {}).get(GRID_PIN)
    if pin is None:
        # `.get(...) is None` rather than falsy: the `.gpkg` site writes the key with a
        # None value for a grid that cannot be pinned.
        return metadata
    certified = metadata.get(GRID_PIN_OPTIONS)
    hit = (
        sorted(touched)
        if certified is None
        else sorted(set(touched) & set(certified))
    )
    if not hit:
        return metadata
    voided = {
        key: value
        for key, value in metadata.items()
        if key not in (GRID_PIN, GRID_PIN_OPTIONS)
    }
    # Plain JSON types: a run manifest is `json.dumps`'d with no `default=`.
    voided[PIN_VOIDED] = {"was": pin, "options": hit}
    return voided


# MF6 requires these three exchanges to join models on the SAME discretization --
# `exg-gwfgwt.f90:224` terminates on a node-count mismatch and `:239-255` on an
# idomain mismatch. It compares COUNTS and IDOMAIN only, never geometry, so two
# different meshes with the same `ncpl` and the same idomain pass both of its checks
# and it then maps node-for-node by index onto the wrong geometry.
#
# `gwfgwf` is deliberately absent: coupling a coarse model to a fine one is what that
# exchange is FOR.
SAME_MESH_EXCHANGES = frozenset(
    {"build_gwf_gwt_exchange", "build_gwf_gwe_exchange", "build_gwf_prt_exchange"}
)


#: Metadata key stamped by ``with_grid`` onto a package that is grid-bound but
#: carries no pin -- hand-written cellids, a ``.flopy`` list form, a per-cell array.
#: It records the mesh the package was moved AWAY from. The build refuses a package
#: carrying it, because nothing else can: its cell numbers are all still in range on
#: the new mesh and every one of them names a different cell.
MOVED_FROM = "myflopy_moved_from"

#: Metadata key listing the packages the caller asserted are mesh-independent, so a
#: stamp on them is cleared rather than refused.
GRID_FREE = "myflopy_grid_free"


def package_move(package, *, grid_free: tuple[str, ...] = ()) -> tuple[str, str]:
    """Classify one package for a grid move: ``(how, why)``.

    The four verdicts, and what each means for the build that follows:

    ``"rebuilt"``
        It carries a recipe rather than resolved records -- a lazy DISV today -- so it
        resolves itself against the new mesh. Nothing to check.
    ``"stale"``
        It carries a PIN naming another mesh. Refused at build by the pin check
        (ledger 250/267); re-declare it against the new grid.
    ``"unverifiable"``
        It names cells and carries no pin -- hand-written cellids, a ``.flopy`` form --
        OR its pin was VOIDED because an edit replaced the records the library
        resolved. **You cannot pin cell numbers the caller typed**, so nothing can tell
        a correct one from a stale one: every number is still in range on a larger mesh.
        Refused unless named in ``grid_free=``.
    ``"unchanged"``
        It names no cells at all -- a scalar ``npf``, a ``tdis`` -- or the caller
        asserted it is mesh-independent. It moves untouched.

    **The recipe test comes LAST, and that is the fix in ledger 275.** It used to come
    first, so a truthy ``inputs`` reported "nothing to check" before the hand-written
    DISV, the pin and the cellid test had been consulted -- measured, a pinned CHD with
    any ``inputs`` entry at all read ``rebuilt`` while its pin sat unread. Harmless
    while only a lazy DISV carried ``inputs``; B2 puts a source recipe in every lazy
    ``.gpkg``, which is exactly a package that carries BOTH a recipe and grid-bound
    content.
    """

    name = getattr(package, "name", "?")
    if name in grid_free:
        return "unchanged", "asserted mesh-independent by grid_free="
    options = getattr(package, "options", None) or {}
    inputs = getattr(package, "inputs", None) or {}
    metadata = getattr(package, "metadata", None) or {}
    if any(key in options for key in ("cell2d", "vertices", "ncpl")):
        # A hand-written DISV IS the old mesh -- it holds that grid's `cell2d` and
        # `vertices` verbatim. It carries no pin and names no cellid, so both of the
        # tests below miss it, and the first report called it "unchanged".
        return "stale", "it describes the old mesh itself (hand-written discretization)"
    pin = metadata.get(GRID_PIN)
    if pin is not None:
        certified = metadata.get(GRID_PIN_OPTIONS)
        if certified and set(certified) <= set(inputs):
            # Every option the pin speaks for is now supplied by the recipe, so the pin
            # describes nothing that survives to the new mesh. This is the contract
            # B2's lazy `.gpkg` has to meet: the source spec lands under the certified
            # name, and the stale pin stops masking it and being masked by it.
            return "rebuilt", "its pinned records are supplied by a recipe now"
        return "stale", f"its records are pinned to mesh {pin}"
    voided = metadata.get(PIN_VOIDED)
    if voided:
        return (
            "unverifiable",
            f"its records were edited after the library resolved them on mesh "
            f"{voided.get('was')}",
        )
    if any(True for _ in _cellids_in(options)):
        return "unverifiable", "names cells the library did not resolve, so it carries no pin"
    if inputs:
        return "rebuilt", "carries a recipe, so it resolves against the new mesh"
    return "unchanged", "names no cells"


def move_report(model, *, grid_free: tuple[str, ...] = ()) -> list[dict[str, str]]:
    """One row per package: what a grid move does to it, before anything is built."""

    return [
        {"package": getattr(package, "name", "?"), "how": how, "why": why}
        for package in getattr(model, "packages", ())
        for how, why in [package_move(package, grid_free=grid_free)]
    ]


def coupled_mesh_complaints(exchanges, built_models) -> list[str]:
    """Complaints for exchanges whose two models resolved to different meshes.

    Measured before this existed: a GWF on a 3x2 Voronoi and a GWT on a 2x3, both
    `ncpl=6`, joined by `mf.exchange` -- built clean, ran, zero warnings, and the
    concentration field was wrong. Nothing compared siblings; `_discretization_complaints`
    validates each model against its OWN context and never looks across.

    DISV only. A structured pair is not compared, because the pin is a DISV
    connectivity hash; MF6's node-count check is all that stands behind a DIS pair,
    and it has the same blind spot.
    """

    complaints: list[str] = []
    for exchange in exchanges:
        if getattr(exchange.builder, "__name__", "") not in SAME_MESH_EXCHANGES:
            continue
        pins: dict[str, str | None] = {}
        for name in exchange.models:
            built = built_models.get(name)
            model = getattr(built, "model", None)
            disv = model.get_package("disv") if model is not None else None
            pins[name] = disv_pin(disv) if disv is not None else None
        known = {name: pin for name, pin in pins.items() if pin is not None}
        if len(known) < 2 or len(set(known.values())) == 1:
            continue
        pairs = ", ".join(f"{name} is on mesh {pin}" for name, pin in known.items())
        complaints.append(
            f"exchange '{exchange.name}' joins models on different meshes: {pairs}. "
            "MODFLOW 6 requires this exchange's models to share a discretization, and it "
            "compares only node counts and idomain -- so a mesh of the same size passes "
            "its check and it maps node-for-node onto the wrong geometry."
        )
    return complaints


def moved_package_complaints(packages) -> list[str]:
    """Complaints for packages `with_grid` stamped and nobody re-declared.

    The stamp is the only thing standing behind an UNPINNED grid-bound package across
    a move: you cannot pin cell numbers the caller typed, and on a larger mesh every
    one of them is still in range and every one names a different cell. Measured
    before the stamp: a hand-written `wel` cellid and a `chd.flopy` list both built
    clean across a 28 -> 92 move (ledger 271).
    """

    complaints = []
    for package in packages:
        pin = (getattr(package, "metadata", None) or {}).get(MOVED_FROM)
        if pin is None:
            continue
        complaints.append(
            f"package {getattr(package, 'name', '?')!r} was moved to this grid by "
            f"`with_grid`, and its cell numbers were typed rather than resolved, so "
            f"they still mean mesh {pin}. Re-declare it against this grid, or name it "
            "in `grid_free=` if it really is mesh-independent."
        )
    return complaints


def _cellids_in(options: dict[str, Any]):
    """Every ``(layer, cell)``-shaped cellid in a package's stress-period data.

    A record's LEADING tuples of ints are its cellids: one for a list BC, two for
    an HFB, which sits on the face between a pair. Reading them all is still
    unambiguous -- a tuple of ints in front of a record cannot be anything else --
    and it is what lets a stale barrier be caught on either side.

    It stops at the first non-cellid element, so it never reaches the value fields,
    and it deliberately does NOT walk ``packagedata``/``connectiondata``, where the
    cellid's position differs per package (SFR's ``packagedata[1]`` is a cellid,
    LAK's is a stage, and SFR's ``connectiondata`` holds reach numbers that would
    read as cells). Guessing there would refuse correct models; advanced-package
    topology is what the grid PIN is for, not a range check.
    """

    import numpy as np

    data = options.get("stress_period_data")
    # A dict of periods is the usual shape, but a bare list of records is accepted
    # by the helpers and written by FloPy, and used to bypass this entirely.
    if isinstance(data, dict):
        blocks = data.values()
    elif isinstance(data, (list, tuple)):
        blocks = (data,)
    else:
        return

    def _cellid(value):
        """``value`` as a cellid tuple, or ``None`` if it is not one."""

        if not isinstance(value, (list, tuple)) or len(value) not in (2, 3):
            return None
        if any(
            isinstance(part, bool) or not isinstance(part, (int, np.integer))
            for part in value
        ):
            return None
        return tuple(int(part) for part in value)

    for records in blocks:
        if not isinstance(records, (list, tuple)):
            continue
        for record in records:
            if not isinstance(record, (list, tuple)):
                continue
            for element in record:
                cellid = _cellid(element)
                if cellid is None:
                    break
                yield cellid

def _mesh_scale(points) -> float:
    """A typical cell size for ``points``, for a tolerance that travels.

    An absolute tolerance cannot work here. The same mesh written to a hundredth of
    a foot, or rendered float32 at State Plane coordinates where the float32 step is
    already ~0.06 ft, differs from itself by more than any fixed epsilon -- while on
    a grid in DEGREES that same epsilon is ~110 m, and a genuinely displaced mesh
    passes. So the tolerance is a thousandth of a cell, and the cell comes from the
    data.
    """

    import numpy as np

    span = points.max(axis=0) - points.min(axis=0)
    area = float(span[0]) * float(span[1])
    if area <= 0.0 or len(points) == 0:
        return float(np.max(np.abs(points))) or 1.0
    return float(np.sqrt(area / len(points)))

def _mesh_residual(built, wanted) -> float:
    """How far the two cell-centre clouds disagree once their frames are aligned.

    Compared INDEX-ALIGNED, after removing each cloud's own centroid. That is exactly
    translation-invariant, which is what a converted model needs: ``to_mf6`` shifts to
    a local origin because MF6 loses DISV conductance precision on State Plane
    coordinates, so the DISV is written in the model's frame while the grid stays in
    the world's (ledgers 159, 218) -- 1.34e+06 length units apart on the synthetic USG
    fixture, and the same mesh.

    It is NOT invariant to cell order, and that is the point. An earlier version
    compared each point's distance from its own centroid, sorted, to survive rotation
    too -- and sorting threw away the one thing a cellid means. A renumbered mesh
    passed it while every package's cellid pointed somewhere else: measured on a
    10x10 grid, cell 5 moved 865.7 length units, the mean cell moved 517.6, and the
    build was clean. Mirroring and an x/y axis swap passed for the same reason.
    Rotation between ``cell2d`` and the grid is not a shape the toolkit produces --
    the local-origin shift only translates, and ``angrot`` lives in the DISV options
    rather than in ``cell2d`` -- so there is nothing to buy with that blind spot.
    """

    import numpy as np

    return float(
        np.max(np.abs((built - built.mean(axis=0)) - (wanted - wanted.mean(axis=0))))
    )

def _discretization_complaints(model, context, packages) -> list[str]:
    """Every way a built model's discretization disagrees with what it was declared against.

    The failure this exists to stop is silent. ``replace_grid`` swaps a grid and
    leaves every package's baked cellids, connection data and dimension counts on
    the old mesh; nothing looked, so MF6 ran the old cell numbers against the new
    geometry and reported Normal termination. Measured on that path: a lazy CHD
    raised ``IndexError`` on a finer grid and resolved to ZERO cells on a coarser
    one, and a lazy disv beside one eager ``drn.gpkg`` ran the DRN on the previous
    grid's cell numbers.

    A range check is not a correctness proof -- old cellids that all happen to fall
    inside a new, larger grid pass it. That gap is what the grid pin closes; this
    catches the part that is decidable from the built model alone, which is the part
    that is silent today.
    """

    import numpy as np

    complaints: list[str] = []
    grid = getattr(context, "grid", None)
    disv = model.get_package("disv") if hasattr(model, "get_package") else None
    ncpl = nlay = None
    if disv is not None:
        ncpl = _as_int(disv.ncpl.get_data())
        nlay = _as_int(disv.nlay.get_data())
        grid_ncpl = _as_int(getattr(grid, "ncpl", None)) if grid is not None else None
        if ncpl is not None and grid_ncpl is not None and grid_ncpl != ncpl:
            complaints.append(
                f"the DISV package has ncpl={ncpl} but the model's context grid has "
                f"{grid_ncpl} cells -- the discretization and the grid are not the "
                "same mesh"
            )
        elif grid is not None and ncpl is not None:
            try:
                props = grid.get_disv_gridprops()
            except Exception:  # noqa: BLE001 - see below; the tuple cannot be closed
                # Deliberately broad, and it only ever WEAKENS this guard. The grid is
                # whatever object a caller put on the context -- a user's own builder
                # result, a reconstructed grid with no triangulation -- so the failures
                # reach flopy, numpy, and any property it computes lazily. A guard that
                # cannot compare must not become the error the user sees.
                logger.debug("could not read gridprops for the consistency guard", exc_info=True)
                props = None
            built_nvert = _as_int(disv.nvert.get_data())
            grid_nvert = _as_int((props or {}).get("nvert"))
            if built_nvert is not None and grid_nvert is not None and built_nvert != grid_nvert:
                complaints.append(
                    f"the DISV package has nvert={built_nvert} but the context grid has "
                    f"{grid_nvert} vertices -- two meshes can share every cell centre and "
                    "still have different cell boundaries"
                )
            elif props is not None:
                # Same cell COUNT is not the same mesh: a 3x2 and a 2x3 grid both have 6.
                # Compared against the grid's own gridprops rather than its polygon
                # centroids, because a Voronoi cell's DISV centre is not its centroid --
                # measured 38.9 length units apart on one 4x4 grid, which would condemn
                # every model.
                built = np.array(
                    [[r[1], r[2]] for r in disv.cell2d.get_data()], dtype=float
                )
                wanted = np.array([[r[1], r[2]] for r in props["cell2d"]], dtype=float)
                if built.shape == wanted.shape and built.size:
                    residual = _mesh_residual(built, wanted)
                    tolerance = 1e-3 * _mesh_scale(wanted)
                    if residual > tolerance:
                        complaints.append(
                            "the DISV package and the context grid have the same cell "
                            f"count but different geometry -- once both frames are "
                            f"aligned their cell centres still differ by up to "
                            f"{residual:g} length units, against a tolerance of "
                            f"{tolerance:g}. A renumbering counts: cellid N means a "
                            "different cell in each mesh"
                        )

    domain = getattr(context, "domain", None)
    if domain is not None and ncpl is not None:
        shape = np.asarray(domain).shape
        if shape and shape[-1] != ncpl:
            complaints.append(
                f"the model's context domain is shaped {shape}, which is {shape[-1]} "
                f"cells, against the built DISV's ncpl={ncpl}"
            )
        elif len(shape) > 1 and nlay is not None and shape[0] != nlay:
            complaints.append(
                f"the model's context domain is shaped {shape}, which is {shape[0]} "
                f"layers, against the built DISV's nlay={nlay}. A two-dimensional "
                "domain answers per layer, so this one cannot answer for the layers it "
                "has no row for -- for one activity in every layer, give it as a flat "
                "(ncpl,) array"
            )

    # Three sources for "which mesh is this model on", strongest first. The order is
    # load-bearing: the recorded pin must be consulted BEFORE the DISV's own, because
    # comparing the DISV against the record is the whole check when a grid is gone.
    # Taking the DISV's fingerprint first makes it agree with itself and says nothing.
    model_pin = grid_pin(grid)
    if model_pin is None:
        # A reloaded model has no grid object: `ModelSpec.to_dict` cannot carry one,
        # so every check above is skipped and this guard used to evaporate on a save
        # and load -- a mismatched model that raised live built clean afterwards, and
        # a package pin survived the round trip and was then never compared. The
        # recorded pin is 16 characters and restores both (ledger 257).
        model_pin = (getattr(context, "metadata", None) or {}).get(GRID_PIN)
        built_pin = disv_pin(disv)
        if model_pin is not None and built_pin is not None and built_pin != model_pin:
            complaints.append(
                "the DISV package is not the mesh this model was declared against: it "
                f"fingerprints as {built_pin} and the model recorded {model_pin}. The "
                "grid object itself is not carried across a save, so the cell NUMBERING "
                "is all that could be compared here -- a mesh with the same "
                "connectivity and different coordinates would have passed. Declare the "
                "model with mf.gwf(..., grid=mf.grid_ref(key)) to get the full geometry "
                "check back"
            )
        if model_pin is None and built_pin is not None:
            # No grid and nothing recorded -- `mf.gwf(name, packages=[...])`, which
            # this repo's own examples use. The built DISV can still answer for the
            # mesh, which is what gives the package pins below something to compare
            # against; without it a package resolved on another mesh built CLEAN.
            model_pin = built_pin
    if model_pin is not None:
        for package in packages:
            metadata = package.metadata or {}
            pin = metadata.get(GRID_PIN)
            if pin is not None and pin != model_pin:
                complaints.append(
                    f"package {package.name!r} was resolved against a different mesh: "
                    f"its records are pinned to grid {pin}, and this model's grid is "
                    f"{model_pin}. Cell numbers do not carry between meshes, and a "
                    "number that happens to be in range is still the wrong cell"
                )
                continue
            voided = metadata.get(PIN_VOIDED)
            if voided and voided.get("was") not in (None, model_pin):
                # The sentence that answers a refusal the old one could not. Before
                # ledger 275 an edited package kept its pin and was told "re-declare it
                # against the new grid" -- so a user who had ALREADY typed correct
                # numbers for this mesh was refused forever, by a message about records
                # that no longer existed. Naming the edit is what makes the refusal
                # actionable, and the edit route can never re-pin (ledger 250: you
                # cannot pin cell numbers the caller typed).
                complaints.append(
                    f"package {package.name!r} holds records that were EDITED after the "
                    f"library resolved them on mesh {voided['was']}, replacing "
                    f"{', '.join(voided['options'])} -- and this model's grid is "
                    f"{model_pin}. Nothing verified the cell numbers that are there "
                    "now. If you typed them for THIS grid, declare the package fresh "
                    "rather than editing a resolved one; otherwise re-declare it "
                    "against this grid"
                )

    if ncpl is not None:
        for package in packages:
            stale = sorted(
                {
                    cellid
                    for cellid in _cellids_in(package.options)
                    if cellid[-1] >= ncpl or (nlay is not None and cellid[0] >= nlay)
                }
            )
            if stale:
                listing = ", ".join(str(cellid) for cellid in stale[:5])
                if len(stale) > 5:
                    listing += f", and {len(stale) - 5} more"
                complaints.append(
                    f"package {package.name!r} holds {len(stale)} cellid(s) this model's "
                    f"grid does not have (nlay={nlay}, ncpl={ncpl}): {listing}"
                )
    return complaints

def dropped_context_fields(context) -> tuple[str, ...]:
    """Which of a context's unserializable fields are set, and so will be lost."""

    return tuple(
        name
        for name in DROPPED_CONTEXT_FIELDS
        if getattr(context, name, None) is not None
    )

def _context_metadata(context) -> dict[str, Any]:
    """A context's metadata, plus the grid fingerprint that outlives its grid.

    The grid itself is not serialized -- inlining one was considered and refused
    (ledger 193) because it floods ``Project.validate`` and breaks ``to_yaml``, and it
    measures ~3.3 MB per model at 9,405 cells. Sixteen characters of it survive
    instead, which is enough for the build-time guard to tell one mesh from another.
    """

    metadata = dict(getattr(context, "metadata", None) or {})
    grid = getattr(context, "grid", None)
    if grid is not None:
        # A present grid IS the declaration, so it overwrites any pin carried in from
        # an earlier load. Refusing to refresh froze the first pin ever written: a
        # project loaded, re-gridded and saved again recorded the fingerprint of a
        # mesh it no longer contained, and was then refused forever.
        pin = grid_pin(grid)
        if pin is not None:
            metadata[GRID_PIN] = pin
    return metadata


__all__ = [
    "DROPPED_CONTEXT_FIELDS",
    "GRID_PIN",
    "_as_int",
    "_canonical_ring",
    "_cellids_in",
    "_context_metadata",
    "_discretization_complaints",
    "_mesh_residual",
    "_mesh_scale",
    "disv_pin",
    "dropped_context_fields",
    "grid_pin",
    "pin_metadata",
]

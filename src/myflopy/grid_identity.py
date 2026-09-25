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
from numbers import Real
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

#: The option names a package supplies ONE VALUE PER CELL under, when the library could
#: not pin it. The array-form boundary conditions address cells POSITIONALLY -- under
#: ``READASARRAYS`` "one value per column" is the format's own definition -- so they carry
#: no cellid for :func:`_cellids_in` to find and no pin when no grid was resolved against.
#: Without this marker such a package answered ``"unchanged"``, the verdict that means it
#: moves untouched (ledger 282).
CELLWISE_OPTIONS = "grid_cellwise_options"

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

#: Metadata key holding the lazy recipe an EAGER `.gpkg` declaration could be rebuilt
#: from. It is what lets `with_grid` re-lazify the largest `stale` class EXACTLY rather
#: than approximating it through an overlay (ledger 296). Excluded from
#: `package_digest` like the rest of `metadata`, so recording it moves no stored digest.
MOVE_RECIPE = "move_recipe"


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

#: Stamped on a package whose per-cell arrays were CARRIED onto the new mesh by the
#: overlay engine (ledger 300). It records the old mesh's pin and what the remap had to
#: guess, and it is what makes the fifth verdict answerable after the fact -- the moved
#: spec no longer holds the old grid, so nothing else could say it.
#: It is stamped INSTEAD of :data:`MOVED_FROM`, never beside it: that one is what
#: ``moved_package_complaints`` refuses the build on, so leaving it would refuse the very
#: package the remap just fixed.
REMAPPED_FROM = "myflopy_remapped_from"

#: Names the SIBLING packages whose features this one addresses by INDEX rather than by
#: cellid -- an MVR's providers and receivers. Stamped by ``mvr_spec``, so both
#: ``mf.mvr(moves=...)`` and ``mf.mvr.flopy(...)`` carry it.
#:
#: It exists because ``package_move`` sees ONE package and a mover's staleness lives
#: entirely in its siblings: nothing about the mover's own options changes when the mesh
#: does, so it answered "unchanged -- names no cells" beside a stale SFR (ledger 303).
#: Deliberately NOT ``requires=``, which already means "build this first" and is carried by
#: ``mf.hfb`` as ``("disv",)`` -- reading that as an index reference would refuse every
#: barrier whenever the discretization rebuilt, which is a different and correct thing.
INDEXED_PACKAGES = "grid_indexed_packages"

#: Metadata key listing the packages the caller asserted are mesh-independent, so a
#: stamp on them is cleared rather than refused.
GRID_FREE = "myflopy_grid_free"


def column_count(value) -> int | None:
    """How many columns a payload supplies a value for, or ``None`` if it is not per-column.

    A numpy array answers with ``shape[-1]``; a flat list or tuple of numbers answers
    with its length; a NESTED sequence (a per-layer list of per-cell lists, which FloPy
    accepts for ``npf(k=...)``) answers with its inner length when the rows agree.

    Reading ``.shape`` alone was the first gap: ``recharge=[1e-4] * ncpl`` is an
    ordinary way to write a per-column payload and was invisible here (ledger 282).
    Missing the nested case is the second, and it is what a per-cell ``k`` on a
    multi-layer model looks like written by hand. Strings are excluded because a
    time-series name is a scalar.

    Lives here, in import layer 0, because three callers need the SAME rule and they sit
    in layers 0, 2 and 5: ``package_move`` below, ``advanced._array_bc_spec``, and
    ``package_api``'s lazy areal pin.
    """

    shape = getattr(value, "shape", None)
    if shape is not None:
        return int(shape[-1]) if shape else None
    if not isinstance(value, (list, tuple)) or not value or isinstance(value, str):
        return None
    if all(isinstance(item, Real) for item in value):
        return len(value)
    inner = {column_count(item) for item in value}
    # One answer, and a real one: rows of different widths are not a per-cell payload,
    # and neither is a list holding a mix of arrays and scalars.
    return inner.pop() if len(inner) == 1 and None not in inner else None


def cellwise_options(options, ncpl: int | None = None) -> tuple[str, ...]:
    """WHICH options address cells: their NAMES, never a bool.

    With ``ncpl`` the test is EXACT -- a payload supplying one value per cell of that
    mesh. Without it, the test is "is this per-column at all", which is the
    conservative substitute used where the format itself defines a payload as
    per-column (``READASARRAYS``). The two are not equivalent and the caller must know
    which it is asking for: ``npf(k=[10.0, 1.0])`` is a per-LAYER list and would be a
    false positive under the ncpl-free form, which is why ``package_move`` only ever
    uses the exact one.

    The names rather than a bool because a pin declares exactly which options it speaks
    for, so an edit to a genuinely grid-free option beside a resolved one does not void
    it (ledger 275).
    """

    def per_cell(value) -> bool:
        if isinstance(value, dict):
            return any(per_cell(item) for item in value.values())
        count = column_count(value)
        return count is not None if ncpl is None else count == ncpl

    return tuple(key for key, value in options.items() if per_cell(value))


def package_move(
    package, *, grid_free: tuple[str, ...] = (), ncpl: int | None = None
) -> tuple[str, str]:
    """Classify one package for a grid move: ``(how, why)``.

    The five verdicts, and what each means for the build that follows:

    ``"remapped"``
        Its per-cell arrays were CARRIED onto the new mesh by the overlay engine
        (ledger 300), so they describe this mesh and the build accepts them. Distinct
        from ``"rebuilt"`` on purpose: a rebuilt package re-derives itself from a recipe
        and is exact, where this one was INTERPOLATED and may have had cells filled from
        a neighbour. Where a package has a source, re-declaring it is measurably more
        faithful than remapping it.
    ``"rebuilt"``
        It carries a recipe rather than resolved records -- a lazy DISV today -- so it
        resolves itself against the new mesh. Nothing to check.
    ``"stale"``
        It carries a PIN naming another mesh. Refused at build by the pin check
        (ledger 250/267); re-declare it against the new grid.
    ``"unverifiable"``
        It names cells and carries no pin -- hand-written cellids, a ``.flopy`` form --
        OR it supplies one value per cell POSITIONALLY and carries no pin (the array
        form, ledger 282) -- OR its pin was VOIDED because an edit replaced the records
        the library resolved. **You cannot pin cell numbers the caller typed**, so
        nothing can tell a correct one from a stale one: every number is still in range
        on a larger mesh, and every array is still the right length on a mesh of the
        same size. Refused unless named in ``grid_free=``.
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
    remapped = metadata.get(REMAPPED_FROM)
    if remapped is not None:
        # Before the pin and the cellwise branches, because a remapped package's values
        # ARE on this mesh -- it is the one case where a per-cell array can be trusted.
        # Reported rather than silently `unchanged`: the values were INTERPOLATED, and a
        # reader deciding whether to trust a run needs to know that and how much of it
        # was filled from a neighbour rather than derived.
        filled = int(remapped.get("filled") or 0)
        detail = f", {filled} cells filled from the nearest old cell" if filled else ""
        return (
            "remapped",
            f"its per-cell arrays were carried from mesh {remapped.get('was')} by area "
            f"overlay{detail}",
        )
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
    moved = metadata.get(MOVED_FROM)
    if moved is not None:
        # ALREADY DECIDED, by the `with_grid` that stamped it. Reporting it afresh would
        # need the OLD mesh, which the moved spec no longer carries -- so without this
        # the report read AFTER a move disagreed with the build that follows, saying
        # `unchanged` about the very package `moved_package_complaints` refuses. The
        # report exists to be read before `prepare_run`; it has to say the same thing.
        return (
            "unverifiable",
            f"`with_grid` could not verify it on this mesh; its values still mean "
            f"mesh {moved}",
        )
    voided = metadata.get(PIN_VOIDED)
    if voided:
        return (
            "unverifiable",
            f"its records were edited after the library resolved them on mesh "
            f"{voided.get('was')}",
        )
    if any(True for _ in _cellids_in(options)):
        return "unverifiable", "names cells the library did not resolve, so it carries no pin"
    cellwise = metadata.get(CELLWISE_OPTIONS)
    if cellwise and not set(cellwise) <= set(inputs):
        # An ARRAY-form package names its cells by POSITION, so it reaches here with no
        # cellid and no pin and used to fall through to "unchanged" -- measured, one
        # `mf.rch.array` declaration wrote mesh A's recharge pattern onto a different
        # mesh of the same size, clean and unwarned (ledger 282). The `<= inputs` test
        # mirrors the pin branch above: once a recipe supplies these options, the array
        # is rebuilt rather than carried.
        return (
            "unverifiable",
            "supplies one value per cell and carries no pin, so nothing can tell a "
            "payload shaped for this mesh from one shaped for another",
        )
    # NOTHING STAMPS npf/ic/sto, because they take no `context=` and so have no mesh
    # to measure a payload against at declaration. `ncpl` is the caller's -- `with_grid`
    # still holds the OLD grid when it classifies, since it does not mutate. Without it
    # this cannot fire, which is the honest limit: a model whose context carries no grid
    # has nothing to compare a length to.
    #
    # EXACT, never the conservative form: `npf(k=[10.0, 1.0])` is a per-LAYER list and
    # guessing from "it is a sequence" would refuse a correct move (ledger 295).
    if ncpl is not None:
        named = cellwise_options(options, ncpl)
        if named and not set(named) <= set(inputs):
            return (
                "unverifiable",
                f"supplies one value per cell ({', '.join(sorted(named))}) and carries "
                "no pin, so nothing can tell a payload shaped for this mesh from one "
                "shaped for another",
            )
    if inputs:
        return "rebuilt", "carries a recipe, so it resolves against the new mesh"
    return "unchanged", "names no cells"


def package_moves(
    packages, *, grid_free: tuple[str, ...] = (), ncpl: int | None = None
) -> dict[str, tuple[str, str]]:
    """Every package's verdict, with a sibling-INDEXING package corrected downward.

    Two passes, because :func:`package_move` is handed ONE package and one package is not
    enough to classify a MOVER. An MVR's records name a provider and a receiver by
    PACKAGE NAME plus an INTEGER -- a reach number, a lake number, a UZF object number --
    and none of its own options change when the mesh does. So it answered ``unchanged``,
    "names no cells", beside a stale SFR, while holding an index into the old reach
    network: measured, a mover pointing at reach 9 of 10 on one mesh, where the same
    declaration on a finer mesh resolves reach 20 of 21 (ledger 303).

    **A provider that is anything but ``unchanged`` invalidates the index, ``rebuilt``
    included** -- that is the whole point. A rebuilt SFR re-intersects its centreline and
    renumbers every reach, so a correct move of the provider is exactly what makes a
    hand-typed index wrong. Today the provider's own pin refuses the build first, which
    is what has been hiding this; the moment a provider becomes movable the mover is the
    only thing standing between a re-grid and water routed to the wrong reach.

    All three callers of :func:`package_move` go through this -- :func:`move_report`,
    ``ModelSpec.with_grid``'s stamping loop, and ``Run.resolution``. Updating two of the
    three would leave the run table and the pre-build report saying different things
    about the same mover, which is the disagreement ``move_report``'s own docstring
    exists to prevent.

    **DELIBERATELY CONSERVATIVE about WHICH provider invalidates the index.** A reach
    number is mesh-derived and a LAKE number is feature-derived, so a mover into a lake
    would usually survive a move that a mover into a stream would not -- and this names
    both. Two reasons. A lake number is only usually stable: a lake that lands on no
    active cell of the new mesh (ledger 211) is dropped and renumbers the rest. And the
    two errors are not symmetric -- a false refusal costs a re-declaration, where a false
    pass routes water to the wrong reach under Normal termination.

    The message is the SPECIFIC one only on the PRE-move report, which is the one the
    docs tell you to read. After ``with_grid`` the mover carries ``MOVED_FROM`` and
    reports through that branch instead, exactly as every other refused package does --
    the verdict was "ALREADY DECIDED, by the ``with_grid`` that stamped it" (ledger 271).
    """

    verdicts = {
        name: package_move(package, grid_free=grid_free, ncpl=ncpl)
        for package in packages
        if (name := getattr(package, "name", None)) is not None
    }
    for package in packages:
        name = getattr(package, "name", None)
        if name is None or name in grid_free:
            continue
        indexed = (getattr(package, "metadata", None) or {}).get(INDEXED_PACKAGES)
        if not indexed:
            continue
        moved = [
            other
            for other in indexed
            if verdicts.get(other, ("unchanged", ""))[0] != "unchanged"
        ]
        if not moved or verdicts.get(name, ("", ""))[0] != "unchanged":
            continue
        verdicts[name] = (
            "unverifiable",
            "it names features of "
            + ", ".join(sorted(moved))
            + " by INDEX, and "
            + ("that package" if len(moved) == 1 else "those packages")
            + " did not survive this move unchanged -- a reach or lake number means "
            "something else on the new mesh",
        )
    return verdicts


def move_report(
    model, *, grid_free: tuple[str, ...] = (), ncpl: int | None = None
) -> list[dict[str, str]]:
    """One row per package: what a grid move does to it, before anything is built.

    ``ncpl`` is the mesh the packages were DECLARED on, so the report can see a
    per-cell array. The report and the stamp must be handed the same one or they
    disagree about the same package -- which is the whole point of reading the report
    before ``prepare_run``.
    """

    packages = tuple(getattr(model, "packages", ()))
    verdicts = package_moves(packages, grid_free=grid_free, ncpl=ncpl)
    return [
        {"package": name, "how": verdicts[name][0], "why": verdicts[name][1]}
        for package in packages
        if (name := getattr(package, "name", None)) is not None
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
        # A numpy RECARRAY is neither a list nor a tuple, and it is the shape FloPy
        # itself produces -- so a package built from `.flopy(...)`, or rebuilt from a
        # written MF6 run, was INVISIBLE here and classified "unchanged", the verdict
        # that means MOVES UNTOUCHED. Measured: 118 cellids found through the `.gpkg`
        # path, 0 through the identical data as a recarray (ledger 281).
        if isinstance(records, np.ndarray):
            records = records.tolist()
        elif not isinstance(records, (list, tuple)):
            continue
        for record in records:
            if isinstance(record, np.void):
                record = record.tolist()
            if not isinstance(record, (list, tuple)):
                continue
            for element in record:
                cellid = _cellid(element)
                if cellid is None:
                    break
                yield cellid


#: The FloPy package types whose MF6 reader REFUSES a cell the discretization removed:
#: every list boundary built on ``BndExtType`` stores ``Cell is outside active grid
#: domain`` for a record whose reduced node is <= 0 -- idomain 0 AND -1 alike
#: (``BoundaryPackageExt.f90``). HFB is deliberately absent: MF6 WARNS about a barrier
#: touching an inactive cell and drops it (``gwf-hfb.f90``), so refusing one here would
#: refuse a model MF6 runs.
_LIST_BC_TYPES = frozenset(
    {"chd", "drn", "ghb", "riv", "wel", "rch", "evt", "cnc", "src", "ctp", "esl"}
)

#: The array-form packages: FloPy type -> (layer-index option, value option, what it
#: carries). MF6 maps an index naming a removed cell to node 0 and skips it, WITHOUT a
#: word (``Disv.f90`` ``nlarray_to_nodelist``, ``gwf-rch.f90`` ``rch_cf``).
_ARRAY_INDEX = {
    "rcha": ("irch", "recharge", "recharge"),
    "evta": ("ievt", "rate", "evapotranspiration"),
}

#: Short spellings of the length and time units ``ModelContext`` carries, for a rate.
_UNIT_ABBREVIATIONS = {
    "feet": "ft", "foot": "ft", "ft": "ft", "meters": "m", "meter": "m", "m": "m",
    "centimeters": "cm", "cm": "cm", "days": "d", "day": "d", "d": "d",
    "seconds": "s", "second": "s", "s": "s", "minutes": "min", "hours": "h",
    "years": "yr", "year": "yr",
}


def _built_idomain(disv, nlay: int, ncpl: int):
    """The built DISV's idomain as ``(nlay, ncpl)`` ints -- ``None`` in FloPy is all active.

    Read off the PACKAGE, never ``model.modelgrid``: flopy synthesizes an all-ones float
    array there where the package's own field says None (ledger 279).
    """

    import numpy as np

    field = getattr(disv, "idomain", None)
    data = field.get_data() if field is not None else None
    if data is None:
        return np.ones((nlay, ncpl), dtype=int)
    return np.asarray(data, dtype=int).reshape(nlay, ncpl)


def _disv_areas(disv, ncpl: int):
    """Each cell's plan area from the DISV's own vertices (the shoelace), or ``None``.

    The DISV rather than the context grid, because it is what MF6 multiplies a rate by.
    """

    import numpy as np

    try:
        vertices = disv.vertices.get_data()
        cells = disv.cell2d.get_data()
        xy = {int(row[0]): (float(row[1]), float(row[2])) for row in vertices}
        areas = np.zeros(ncpl, dtype=float)
        for row in cells:
            ring = [xy[int(index)] for index in list(row)[4 : 4 + int(row[3])]]
            x = np.array([p[0] for p in ring])
            y = np.array([p[1] for p in ring])
            areas[int(row[0])] = 0.5 * abs(float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))
    except (AttributeError, KeyError, IndexError, TypeError, ValueError) as error:
        logger.debug("could not read cell areas off the DISV; rates go unreported: %s", error)
        return None
    return areas


def excavation_cells(context) -> dict[int, str]:
    """``{cell: excavation name}`` for every cell an excavation record ``context`` reaches names.

    Read from what a model can carry a record on, in ``specs._layer_elevations``'
    order: ``context.surfaces`` -- a ``LayerBuildResult`` (its ``excavations``) or a
    table whose ``attrs["excavations"]`` a lazy DISV published -- then the grid's
    ``excavations`` attribute, which ``LayerBuildResult.attach_to_grid`` sets. Only used
    to NAME the cause of a refusal, so every source is read and the first to name a
    cell wins; a record is duck-typed (``.name``, ``.cells``) because this module is
    import layer 0 and ``ExcavationRecord`` lives in layer 2.

    Parameters
    ----------
    context : ModelContext or None
        The model's context. ``None``, or one carrying no record, gives ``{}``.

    Returns
    -------
    dict of int to str
        Cell number (zero-based, a python ``int``) to the name of the excavation that
        dug it.

    Examples
    --------
    >>> import myflopy as mf
    >>> from myflopy.grid_identity import excavation_cells
    >>> from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
    >>> vor = irregular_voronoi_grid(nrow=3, ncol=3, cell_size=100.0)
    >>> pond = mf.Area(vor.gdf_vorPolys.geometry.loc[4])
    >>> result = (mf.LayerStack(top=mf.Flat(100.0)).add("sand", thickness=20.0)
    ...           .excavate(bottom=90.0, within=pond, name="pond").build(vor, attach=True))
    >>> excavation_cells(mf.ModelContext(grid=vor))
    {4: 'pond'}
    """

    found: dict[int, str] = {}
    surfaces = getattr(context, "surfaces", None)
    sources = [getattr(surfaces, "excavations", None)]
    attrs = getattr(surfaces, "attrs", None)
    if isinstance(attrs, dict):
        sources.append(attrs.get("excavations"))
    sources.append(getattr(getattr(context, "grid", None), "excavations", None))
    for records in sources:
        if not isinstance(records, dict):
            continue
        for record in records.values():
            for cell in getattr(record, "cells", ()) or ():
                found.setdefault(int(cell), str(getattr(record, "name", "?")))
    return found


def _dug_note(cells, idomain, dug: dict[int, str], *, stale: bool = True) -> str:
    """Which of ``cells`` -- ``(layer, cell)`` pairs -- an excavation dug out, in words.

    A cell counts when its column is in a record AND it lies above the column's first
    active layer: that is the space the cut removed (or absorbed). A cell further down
    a dug column was removed by something else (a pinch), and is not blamed on the pond.
    ``stale`` adds the likely cause -- a context from before the cut -- which is true of
    a resolved index or cellid and not of a single layer typed for every column.
    """

    import numpy as np

    by_name: dict[str, int] = {}
    for layer, cell in cells:
        name = dug.get(int(cell))
        if name is None:
            continue
        active = np.flatnonzero(idomain[:, cell] > 0)
        if active.size and layer < active[0]:
            by_name[name] = by_name.get(name, 0) + 1
    if not by_name:
        return ""
    parts = [f"{count} inside excavation {name!r}" for name, count in sorted(by_name.items())]
    note = f" -- {', '.join(parts)}, dug out of the stack by `LayerStack.excavate`"
    if stale:
        note += ", so this was declared against a context from BEFORE the cut"
    return note


def _listing(items, limit: int = 5) -> str:
    """The first ``limit`` items, then how many more."""

    items = list(items)
    text = ", ".join(str(item) for item in items[:limit])
    if len(items) > limit:
        text += f", and {len(items) - limit} more"
    return text


def _flopy_types(model) -> dict[str, str]:
    """``{package name: FloPy package type}`` for every package built on ``model``."""

    return {
        _flopy_name(package): str(getattr(package, "package_type", ""))
        for package in getattr(model, "packagelist", ()) or ()
    }


def _inactive_cellid_complaints(model, packages, idomain, dug) -> list[str]:
    """A list boundary or an observation naming a cell the built DISV removed.

    MF6 REFUSES both at read (``Cell is outside active grid domain`` for a boundary,
    ``Error reading data from ID string`` for a head observation -- both measured in
    ``docs/reproducers/active_cell_guard.py``), and FloPy accepts both at construction
    -- measured, a CHD at a dug cell built clean. So this is the
    earlier, named form of an error MF6 would stop on anyway, and can refuse no model
    MF6 runs. What it adds is the CAUSE: a package declared against a context from
    before the stack was excavated (or re-pinched) names cells that no longer exist.
    """

    import numpy as np

    nlay, ncpl = idomain.shape
    complaints: list[str] = []
    types = _flopy_types(model)
    resolved = getattr(model, "_myflopy_resolved", None) or {}
    for package in packages:
        effective = resolved.get(package.name, package)
        if types.get(package.name) not in _LIST_BC_TYPES:
            continue
        removed = sorted(
            {
                cellid
                for cellid in _cellids_in(effective.options)
                if len(cellid) == 2
                and 0 <= cellid[0] < nlay
                and 0 <= cellid[1] < ncpl
                and idomain[cellid[0], cellid[1]] <= 0
            }
        )
        if removed:
            complaints.append(
                f"package {package.name!r} names {len(removed)} cell(s) the built DISV "
                f"has REMOVED (idomain <= 0), which MF6 refuses at read ('Cell is "
                f"outside active grid domain'): {_listing(removed)}"
                f"{_dug_note(removed, idomain, dug)}. Re-declare it against the context "
                "of the layer stack this DISV was built from"
            )

    for package in getattr(model, "packagelist", ()) or ():
        if getattr(package, "package_type", None) != "obs":
            continue
        try:
            blocks = package.continuous.get_data() or {}
        except AttributeError as error:
            logger.debug("an OBS package carried no continuous block; not checked: %s", error)
            continue
        removed = []
        for records in blocks.values():
            for record in records if records is not None else ():
                for field in ("id", "id2"):
                    value = record[field] if field in (record.dtype.names or ()) else None
                    if (
                        isinstance(value, tuple)
                        and len(value) == 2
                        and all(
                            isinstance(part, (int, np.integer)) and not isinstance(part, bool)
                            for part in value
                        )
                        and 0 <= int(value[0]) < nlay
                        and 0 <= int(value[1]) < ncpl
                        and idomain[int(value[0]), int(value[1])] <= 0
                    ):
                        removed.append((str(record["obsname"]), (int(value[0]), int(value[1]))))
        if removed:
            parent = getattr(package, "parent_file", None)
            owner = (
                f"the observations of package {_flopy_name(parent)!r}"
                if parent is not None and getattr(parent, "package_type", None) != "obs"
                else f"observation file {getattr(package, 'filename', None) or _flopy_name(package)!r}"
            )
            complaints.append(
                f"{owner} name {len(removed)} cell(s) the built DISV has REMOVED "
                "(idomain <= 0), which MF6 cannot read (a head observation there stops it "
                "with 'Error reading data from ID string'): "
                + _listing(f"{name!r} at {cellid}" for name, cellid in removed)
                + _dug_note([cellid for _name, cellid in removed], idomain, dug)
                + ". Observe an active cell -- in a dug column, the one under the bed"
            )
    return complaints


def _flopy_name(package) -> str:
    """A FloPy package's own name (FloPy keeps it in a one-item list)."""

    names = getattr(package, "name", None) or ["?"]
    return str(names[0] if isinstance(names, (list, tuple)) else names)


def _index_is_positional(value) -> bool:
    """Whether a declared ``irch``/``ievt`` names a layer PER COLUMN (an array, or a dict of them).

    ``None`` and a scalar -- also per period -- name one layer everywhere, which is the
    documented opt-out that :func:`active_cell_warnings` reports rather than refuses.
    """

    import numpy as np

    if value is None:
        return False
    if isinstance(value, dict):
        return any(_index_is_positional(item) for item in value.values())
    return np.ndim(value) >= 1


def _period_arrays(datum, ncpl: int) -> dict:
    """``{period: (ncpl,) array or None}`` from FloPy's per-period ``get_data()``.

    ``None`` for a period is FloPy's "nothing declared here" -- MF6 then keeps the
    value in force, which :func:`_latest` finds. A scalar is broadcast, as MF6's
    ``CONSTANT`` is. FloPy returning ``None`` outright (an index never given) is ``{}``.
    """

    import numpy as np

    if datum is None:
        return {}
    if not isinstance(datum, dict):
        datum = {0: datum}
    out: dict = {}
    for period, value in datum.items():
        values = None if value is None else np.asarray(value).reshape(-1)
        if values is not None and values.size in (1, ncpl):
            out[int(period)] = np.broadcast_to(values, (ncpl,)).copy()
        else:
            out[int(period)] = None
    return out


def _latest(by_period: dict, period: int, fill):
    """The value in force at ``period``: the last one declared at or before it, else ``fill``."""

    earlier = [p for p in sorted(by_period) if p <= period and by_period[p] is not None]
    return by_period[earlier[-1]] if earlier else fill


def _read_period_data(package, option: str, ncpl: int):
    """``_period_arrays`` of one FloPy option, or ``None`` when FloPy cannot read it.

    An array given as an ``OPEN/CLOSE`` file name is read from the run's workspace,
    which need not exist yet at build time.
    """

    from flopy.mf6.mfbase import MFDataException

    try:
        return _period_arrays(getattr(package, option).get_data(), ncpl)
    except (MFDataException, OSError, ValueError, TypeError) as error:
        logger.debug("could not read %s off %s; its loss goes unmeasured: %s",
                     option, _flopy_name(package), error)
        return None


def _flux_units(context) -> str:
    """``ft3/d`` for a feet/days model, else the unit names spelled out."""

    length = str(getattr(context, "length_units", None) or "length")
    time = str(getattr(context, "time_units", None) or "time")
    short_l = _UNIT_ABBREVIATIONS.get(length.lower())
    short_t = _UNIT_ABBREVIATIONS.get(time.lower())
    if short_l and short_t:
        return f"{short_l}3/{short_t}"
    return f"{length}^3/{time}"


def _array_index_findings(model, context, packages, idomain, dug):
    """``(complaints, warnings)`` for every array-form RCH/EVT built on the model.

    A per-column index (``irch``/``ievt`` as an array -- ``"top_active"`` resolved,
    eagerly or at build) naming a removed cell in a column that HAS an active cell is
    a COMPLAINT: it was resolved against another domain, and MF6 drops that column's
    flux without a word. ``None``, or one layer for every column, is the documented
    opt-out (ledger 163), so the same loss is a WARNING, with the rate it costs
    (ledger 322). A column with no active cell loses nothing that was not already
    outside the model, and is not counted.
    """

    import numpy as np

    nlay, ncpl = idomain.shape
    complaints: list[str] = []
    warned: list[str] = []
    resolved = getattr(model, "_myflopy_resolved", None) or {}
    declared = {package.name: resolved.get(package.name, package) for package in packages}
    has_active = (idomain > 0).any(axis=0)
    columns = np.arange(ncpl)
    areas = None
    units = _flux_units(context)
    for flopy_package in getattr(model, "packagelist", ()) or ():
        kind = getattr(flopy_package, "package_type", None)
        if kind not in _ARRAY_INDEX:
            continue
        index_name, value_name, noun = _ARRAY_INDEX[kind]
        name = _flopy_name(flopy_package)
        spec = declared.get(name)
        positional = spec is not None and (
            index_name in (spec.inputs or {})
            or _index_is_positional((spec.options or {}).get(index_name))
        )
        indices = _read_period_data(flopy_package, index_name, ncpl)
        if indices is None:
            continue
        rates = _read_period_data(flopy_package, value_name, ncpl) or {}
        top = np.zeros(ncpl, dtype=int)      # no index: MF6 uses layer 1 everywhere
        findings = []
        for period in sorted(set(indices) | set(rates) | {0}):
            index = np.asarray(_latest(indices, period, top), dtype=int)
            in_range = (index >= 0) & (index < nlay)
            removed = has_active & in_range & (idomain[np.where(in_range, index, 0), columns] <= 0)
            beyond = has_active & ~in_range
            if not (removed.any() or beyond.any()):
                continue
            rate = _latest(rates, period, None)
            lost = total = None
            if rate is not None:
                if areas is None:
                    areas = _disv_areas(model.get_package("disv"), ncpl)
                if areas is not None:
                    flux = np.abs(np.asarray(rate, dtype=float)) * areas
                    lost = float(flux[removed].sum())
                    total = float(flux[has_active].sum())
            findings.append((period, index, removed, beyond, lost, total))
        if not findings:
            continue

        period, index, removed, beyond, lost, total = findings[0]
        pairs = [(int(index[cell]), int(cell)) for cell in np.flatnonzero(removed)]
        where = _listing(f"cell {cell} (layer {layer + 1})" for layer, cell in pairs)
        if beyond.any():
            where = (where + "; " if where else "") + (
                f"{int(beyond.sum())} column(s) whose {index_name} names a layer this DISV "
                f"does not have (nlay={nlay}), which MF6 stops on"
            )
        cost = ""
        if lost is not None and total:
            cost = (
                f": {lost:.6g} of {total:.6g} {units} ({100.0 * lost / total:.1f}%) in "
                f"stress period {period + 1}"
            )
            worst = max(findings, key=lambda item: item[4] if item[4] is not None else -1.0)
            if worst[4] is not None and worst[4] > lost and worst[5]:
                cost += (
                    f", and {worst[4]:.6g} ({100.0 * worst[4] / worst[5]:.1f}%) in "
                    f"stress period {worst[0] + 1}"
                )
        count = int(removed.sum() + beyond.sum())
        note = _dug_note(pairs, idomain, dug)
        if positional:
            complaints.append(
                f"package {name!r} ({noun} as arrays) addresses its {noun} in {count} "
                f"column(s) to a cell the built DISV has REMOVED although the column has "
                f"an active cell below: {where}{note}. MF6 drops {noun} addressed to a "
                f"removed cell without a word{cost}. Its {index_name} was resolved "
                "against a domain that is not this DISV's: declare the package without "
                f"context= (so {index_name}='top_active' resolves at build against this "
                "DISV) or against the excavated result's idomain"
            )
        else:
            given = (spec.options or {}).get(index_name) if spec is not None else None
            said = f"{index_name} left out" if given is None else f"{index_name}={given!r}"
            note = _dug_note(pairs, idomain, dug, stale=False)
            warned.append(
                f"package {name!r} ({noun} as arrays): {said} puts {noun} in "
                f"one layer in every column, and the built DISV has REMOVED that cell in "
                f"{count} column(s) that have an active cell below: {where}{note}. MF6 "
                f"drops {noun} addressed to a removed cell without a word{cost}. Leave "
                f"{index_name} at its default, 'top_active', to put it on each column's "
                "highest active cell"
            )
    return complaints, warned


def _activity_inputs(model):
    """The built DISV's ``(nlay, ncpl)`` idomain, or ``None`` when there is nothing to check.

    ``None`` for a model with no DISV (the activity checks are DISV-only, as the mesh
    checks are), and for one whose every cell is active -- nothing can name a removed
    cell there, so the per-package walks are skipped entirely.
    """

    disv = model.get_package("disv") if hasattr(model, "get_package") else None
    if disv is None or isinstance(disv, list):
        return None
    ncpl = _as_int(disv.ncpl.get_data())
    nlay = _as_int(disv.nlay.get_data())
    if not ncpl or not nlay:
        return None
    try:
        idomain = _built_idomain(disv, nlay, ncpl)
    except ValueError as error:
        # An idomain that is not (nlay, ncpl) is not this check's to report.
        logger.debug("the built idomain could not be shaped; activity not checked: %s", error)
        return None
    return None if (idomain > 0).all() else idomain


def _active_cell_complaints(model, context, packages) -> list[str]:
    """What names a cell the BUILT DISV has removed, as refusals (ledger 321).

    Separate from the mesh checks because the question is different: not "is this the
    same mesh" but "does this package still address cells that exist on it". The route
    to a wrong answer is a context declared BEFORE the stack was excavated (or
    re-pinched): the DISV holds the cut, while every package resolved from the old
    domain still names the dug cells.
    """

    idomain = _activity_inputs(model)
    if idomain is None:
        return []
    dug = excavation_cells(context)
    complaints = _inactive_cellid_complaints(model, packages, idomain, dug)
    return complaints + _array_index_findings(model, context, packages, idomain, dug)[0]


def active_cell_warnings(model, context, packages) -> list[str]:
    """The WARNINGS of the build-time activity check: an ``irch``/``ievt`` opt-out that loses flux.

    ``irch=None`` (and ``ievt=None``, or one layer given for every column) is the
    documented way to put an array-form flux in one layer everywhere (ledger 163). It
    stays legal -- this does not refuse it -- but where that layer's cell is REMOVED in
    a column that has an active cell below, MF6 drops the column's flux without a word,
    and that is now said, with the rate it costs (ledger 322).

    Parameters
    ----------
    model : flopy.mf6.MFModel
        The model, after every package is built on it.
    context : ModelContext or None
        The model's context: its units label the rate, and any excavation record it
        reaches names the cells an excavation dug.
    packages : sequence of PackageSpec
        The packages built from specs, as ``ModelSpec.build`` passes them; each
        array-form package's declared ``irch``/``ievt`` decides warning against
        refusal.

    Returns
    -------
    list of str
        One message per package that loses flux; empty when none does, when the model
        has no DISV, or when every cell of it is active.

    Examples
    --------
    >>> from myflopy.grid_identity import active_cell_warnings
    >>> active_cell_warnings(object(), None, [])
    []
    """

    idomain = _activity_inputs(model)
    if idomain is None:
        return []
    return _array_index_findings(model, context, packages, idomain, excavation_cells(context))[1]


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
    # In range is not enough: a cell can exist on this mesh and still have been REMOVED
    # by the DISV's idomain -- what a context declared before an excavation leaves.
    complaints += _active_cell_complaints(model, context, packages)
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
    "active_cell_warnings",
    "_canonical_ring",
    "_cellids_in",
    "_context_metadata",
    "_discretization_complaints",
    "_mesh_residual",
    "_mesh_scale",
    "disv_pin",
    "dropped_context_fields",
    "excavation_cells",
    "grid_pin",
    "pin_metadata",
]

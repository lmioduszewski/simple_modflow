"""Copy the files a project's grids are built from into the project itself.

A saved project is a set of instructions, and some of those instructions name files
on disk. `Project.save(bundle_sources=True)` copies those files into
``<root>/inputs/`` and rewrites the stored paths to point there, so the project is an
archive rather than a set of references into whatever directory it was authored in.

**Grids AND lazy packages** (ledger 293). A ``GridSpec`` was always a live source --
``ModelSpec.build`` re-resolves it on every build -- and since ledgers 270/279/284/285
so is any package declared without a ``context=``: it holds a recipe rather than
records, and reads its file at build time.

An EAGER package still names nothing this copies, and the measurement behind that is
worth keeping: every `.gpkg` helper given a ``context=`` resolves at declaration, an
instrumented recorder saw zero file reads inside ``PackageSpec.build`` across 306
passing tests, and a project rebuilds with its GeoPackage deleted. So the package walk
finds nothing on an eager spec and returns it unchanged, rather than needing to ask
whether a package is lazy before looking.

**A lazy model has TWO carriers, and the easily-overlooked one fails FIRST.** The
boundary recipes are :class:`~myflopy.sources.DataSourceSpec` objects; the layer stack
behind ``mf.disv(<gridless LayerStack>)`` is not a dataclass, names no source spec, and
per ledger 270 is the ONLY shape a layered model saves in. Its raster is opened inside
the discretization, before any boundary condition is built, so bundling the recipes
alone produces an archive that reports success and cannot open.

Two of the grid's own paths are read LATER rather than at build: ``rasters`` (on the
first ``gdf_topbtm``) and ``idomain_path`` (on the first edge query). They are bundled
anyway, because a project that builds and then dies on its first cross-section is not
archived.
"""

from __future__ import annotations

import filecmp
import hashlib
import shutil
from dataclasses import replace
from pathlib import Path
from typing import Any

from myflopy._logging import get_logger
from myflopy.specs import DataSourceSpec, GridSpec, PackageSpec, repoint_source_paths

logger = get_logger(__name__)


#: Extensions whose format is MULTI-FILE: the named file is unusable without its
#: companions, which share its stem. Copying only the file the recipe names produced a
#: project that reported `bundled` and could no longer be OPENED -- measured, a shapefile
#: written with five sidecars (.cpg .dbf .prj .shp .shx) was bundled as the .shp alone and
#: the grid that resolved before the save then failed `DataSourceError` on a tree where
#: nothing had been deleted (ledger 292). A `.gdb` is a DIRECTORY and is already carried
#: whole by `copytree`; this is the other shape of the same mistake.
#: Raster sidecars are here because they are cheap and their absence is silent: a `.tif`
#: opens without its `.tfw` and lands in the wrong place.
_COMPANION_FORMATS = frozenset(
    {".shp", ".tab", ".mif", ".tif", ".tiff", ".asc", ".img", ".bil", ".bip", ".bsq"}
)


def _companions(source: Path) -> list[Path]:
    """Every sibling belonging to ``source``'s multi-file format, itself excluded.

    Matched by NAME PREFIX rather than by an enumerated suffix list: a shapefile's
    companion set is open (.shx .dbf .prj .cpg .qix .sbn .sbx .qpj ...) and a raster's
    sidecar is often a double extension (`ground.tif.aux.xml`), whose `Path.stem` is
    `ground.tif` rather than `ground`. A prefix match catches both without a table that
    silently omits the one format a user has.
    """

    if source.suffix.lower() not in _COMPANION_FORMATS:
        return []
    try:
        siblings = list(source.parent.iterdir())
    except OSError:
        return []
    prefix = f"{source.stem}."
    return [
        path
        for path in siblings
        if path != source and path.is_file() and path.name.startswith(prefix)
    ]


def _paths_of(spec: GridSpec) -> list[Path]:
    """Every file path a grid recipe names, in no particular order."""

    found: list[Any] = []
    for source in (spec.boundary, spec.refinement):
        if source is not None:
            found.append(source.path)
    for group in (spec.breaklines, spec.points, spec.inputs):
        for item in group:
            if isinstance(item, DataSourceSpec):
                found.append(item.path)
            elif isinstance(item, (str, Path)):
                found.append(item)
    if spec.script is not None:
        found.append(spec.script)
    if spec.pickle_path is not None:
        found.append(spec.pickle_path)
    options = spec.voronoi_options or {}
    rasters = options.get("rasters")
    if isinstance(rasters, (str, Path)):
        found.append(rasters)
    elif isinstance(rasters, (list, tuple)):
        found.extend(item for item in rasters if isinstance(item, (str, Path)))
    if options.get("idomain_path") is not None:
        found.append(options["idomain_path"])
    return [Path(item) for item in found]


def grid_sources(spec: GridSpec, *, project_root: Path | str | None = None) -> list[Path]:
    """The files a grid recipe reads, resolved the way the resolver resolves them.

    A relative path is anchored to the project root exactly as
    ``grid_spec_resolver._source_path`` does, so this reports the file that will
    actually be opened rather than the string that was stored.
    """

    resolved = []
    for path in _paths_of(spec):
        if path.is_absolute() or project_root is None:
            resolved.append(path)
        else:
            resolved.append(Path(project_root) / path)
    return resolved


def _inside(source: Path, directory: Path) -> bool:
    """Whether ``source`` already lives under ``directory``, resolving both."""

    try:
        return source.resolve().is_relative_to(directory.resolve())
    except (OSError, ValueError):
        return False


def _same_file(left: Path, right: Path) -> bool:
    """Whether two paths hold the same bytes -- or the same directory name.

    Used only to decide whether a name in the inputs directory is already claimed by
    a DIFFERENT source. A directory (a ``.gdb``) is compared by name: walking one to
    decide a filename would cost more than the collision it avoids, and two
    differently-named directories never collide in the first place.
    """

    try:
        if left.is_dir() or right.is_dir():
            return left.is_dir() and right.is_dir() and left.name == right.name
        return filecmp.cmp(left, right, shallow=False)
    except OSError:
        return False


class _Bundler:
    """Copies sources into one directory and remembers where each one went."""

    def __init__(self, inputs_dir: Path, root: Path, project_root: Path | None):
        self.inputs_dir = inputs_dir
        self.root = root
        self.project_root = project_root
        self.taken: dict[str, Path] = {}
        self.copied: list[dict[str, str]] = []
        self.missing: list[str] = []

    def _resolve(self, value: Any) -> Path:
        """Where a stored path actually points."""

        path = Path(value)
        if path.is_absolute() or self.project_root is None:
            return path
        return self.project_root / path

    def _destination(self, source: Path) -> Path:
        """A unique name for ``source`` inside the inputs directory.

        Two sources with the same basename from different directories would
        otherwise overwrite one another -- silently, and with the survivor answering
        for both. The second one keeps its name plus a short digest of its origin.

        The claim is checked against the DIRECTORY as well as against this call's own
        bookkeeping, because each grid is bundled separately and a project may be
        saved more than once: an in-memory map alone let two grids collide, and
        collided again on a re-save.
        """

        name = source.name
        claimed = self.taken.get(name)
        if claimed is None:
            landed = self.inputs_dir / name
            if landed.exists() and not _same_file(landed, source):
                claimed = landed        # taken by something else, on a previous pass
        if claimed is not None and claimed != source and not _same_file(claimed, source):
            digest = hashlib.sha256(str(source.parent).encode()).hexdigest()[:8]
            name = f"{source.stem}.{digest}{source.suffix}"
        self.taken.setdefault(name, source)
        return self.inputs_dir / name

    def take(self, value: Any, *, write_target: bool = False) -> Any:
        """Copy one source in and return the path to store, or the original.

        ``write_target=True`` marks a path the project WRITES rather than reads -- the
        GRASS interpolation cache behind a contour surface. An absent one is copied
        over in silence instead of being reported missing: it has simply not been
        produced yet, and listing it would put a fault in the report for a project in
        a perfectly ordinary state.
        """

        if value is None:
            return value
        source = self._resolve(value)
        if not source.exists():
            # Not this verb's complaint to make: the resolver already names a
            # missing source, and refusing to save a project because one file moved
            # would lose the rest of it.
            if not write_target:
                self.missing.append(str(source))
            return value
        # Against the PROJECT ROOT, not just `inputs/`. A source the project already
        # holds does not need copying anywhere -- and testing only `inputs/` made the
        # bundler copy the project's OWN artifact into itself: a `method="object"` grid
        # becomes a pickle on save, and that pickle was then duplicated into `inputs/`
        # and the recipe repointed at the copy, orphaning the original. On a real grid
        # the pickle is the largest thing a saved project holds (ledger 292). It is also
        # what `usg.export_gis(dir)` encourages -- writing sources INTO the project.
        already_inside = _inside(source, self.inputs_dir) or _inside(source, self.root)
        destination = source if already_inside else self._destination(source)
        if not already_inside:
            destination.parent.mkdir(parents=True, exist_ok=True)
            if source.is_dir():
                # A .gdb is a directory, and copying only its "file" would produce a
                # project that looks bundled and is not.
                if destination.exists():
                    shutil.rmtree(destination)
                shutil.copytree(source, destination)
            else:
                shutil.copy2(source, destination)
                for companion in _companions(source):
                    beside = destination.parent / (
                        destination.stem + companion.name[len(source.stem):]
                    )
                    shutil.copy2(companion, beside)
                    self.copied.append(
                        {"from": str(companion), "to": str(beside.relative_to(self.root))}
                    )
            self.copied.append(
                {"from": str(source), "to": str(destination.relative_to(self.root))}
            )
        return destination.relative_to(self.root).as_posix()

    def source(self, spec: DataSourceSpec | None) -> DataSourceSpec | None:
        """Copy a source spec's file in and repoint it, project-relative."""

        if spec is None:
            return None
        stored = self.take(spec.path)
        if stored == spec.path:
            return spec
        # `external` means "never anchor to the project root" -- true while the file
        # lived outside it, and false now that a copy is inside. Bundling an external
        # source is the point of asking for a bundle; it is reported, not silent.
        return replace(spec, path=stored, external=False)


def bundle_grid_sources(
    spec: GridSpec,
    *,
    root: Path | str,
    inputs_dir: Path | str,
    project_root: Path | str | None = None,
) -> tuple[GridSpec, dict[str, Any]]:
    """Return ``spec`` with its sources copied into ``inputs_dir`` and repointed.

    Paths are stored RELATIVE to the project root, which is what makes the result
    portable -- and is also why ``Run`` must know its project root when it rebuilds,
    or a relative path resolves against the current working directory instead.
    """

    root, inputs_dir = Path(root), Path(inputs_dir)
    bundler = _Bundler(inputs_dir, root, Path(project_root) if project_root else None)

    def _group(items):
        return tuple(
            bundler.source(item)
            if isinstance(item, DataSourceSpec)
            else (bundler.take(item) if isinstance(item, (str, Path)) else item)
            for item in items
        )

    options = dict(spec.voronoi_options or {})
    if "rasters" in options:
        rasters = options["rasters"]
        if isinstance(rasters, (str, Path)):
            options["rasters"] = bundler.take(rasters)
        elif isinstance(rasters, (list, tuple)):
            options["rasters"] = type(rasters)(
                bundler.take(item) if isinstance(item, (str, Path)) else item
                for item in rasters
            )
    if options.get("idomain_path") is not None:
        options["idomain_path"] = bundler.take(options["idomain_path"])

    bundled = replace(
        spec,
        boundary=bundler.source(spec.boundary),
        refinement=bundler.source(spec.refinement),
        breaklines=_group(spec.breaklines),
        points=_group(spec.points),
        inputs=_group(spec.inputs),
        script=bundler.take(spec.script) if spec.script is not None else None,
        pickle_path=(
            bundler.take(spec.pickle_path) if spec.pickle_path is not None else None
        ),
        voronoi_options=options,
    )
    report = {"copied": bundler.copied, "missing": bundler.missing}
    if bundler.missing:
        logger.warning(
            "bundling grid %r: %d source(s) could not be found and were left "
            "pointing where they were: %s",
            spec.name, len(bundler.missing), ", ".join(bundler.missing[:5]),
        )
    return bundled, report


def bundle_package_sources(
    spec: PackageSpec,
    *,
    root: Path | str,
    inputs_dir: Path | str,
    project_root: Path | str | None = None,
) -> tuple[PackageSpec, dict[str, Any]]:
    """Return ``spec`` with the files its RECIPE reads copied in and repointed.

    A package only has such files when it was declared LAZILY -- ``mf.chd.gpkg(...)``
    without a ``context=``, or ``mf.disv(<gridless LayerStack>)``. An eager package
    resolved at declaration and holds finished records, so it names no file at build
    time and nothing here finds one; that is the measurement this module's own header
    was written around, and it is still true of the eager form.

    **Both carriers, because the one that is easy to overlook fails first.** A model
    whose layer stack reads ``gis/ground.tif`` dies ``RasterioIOError`` inside the
    discretization, before any boundary condition is built -- so a bundler that carried
    only the ``.gpkg`` recipes would report a complete archive and produce one that
    cannot open. :func:`myflopy.specs.repoint_source_paths` walks both.

    ``report["notes"]`` names what was NOT carried, which is the half a report is
    usually missing: a stack holding per-cell arrays cannot be serialized (ledger 193)
    and is left exactly as it was.
    """

    root, inputs_dir = Path(root), Path(inputs_dir)
    bundler = _Bundler(inputs_dir, root, Path(project_root) if project_root else None)
    notes: list[tuple[str, str]] = []

    def _take(path, *, role, notes):
        """Copy one path in, whatever kind of thing named it."""

        return bundler.take(path, write_target=role == "write_target")

    inputs = repoint_source_paths(spec.inputs, _take, notes=notes) if spec.inputs else spec.inputs
    bundled = spec if inputs == spec.inputs else replace(spec, inputs=inputs)
    report = {"copied": bundler.copied, "missing": bundler.missing, "notes": notes}
    if bundler.missing:
        logger.warning(
            "bundling package %r: %d source(s) could not be found and were left "
            "pointing where they were: %s",
            spec.name, len(bundler.missing), ", ".join(bundler.missing[:5]),
        )
    return bundled, report


__all__ = ["bundle_grid_sources", "bundle_package_sources", "grid_sources"]

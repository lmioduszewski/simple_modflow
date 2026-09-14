"""Copy the files a project's grids are built from into the project itself.

A saved project is a set of instructions, and some of those instructions name files
on disk. `Project.save(bundle_sources=True)` copies those files into
``<root>/inputs/`` and rewrites the stored paths to point there, so the project is an
archive rather than a set of references into whatever directory it was authored in.

**Grids only, and that is a measurement rather than a scope cut.** Every `.gpkg`
package helper resolves EAGERLY at declaration, so a saved package holds finished
records: an instrumented recorder saw zero file reads inside ``PackageSpec.build``
across 306 passing tests, a project rebuilds with its GeoPackage deleted, and the
``source_path`` a package records is read by nothing. Bundling those today would copy
files nothing opens. A ``GridSpec`` is the opposite -- ``ModelSpec.build`` re-resolves
it on EVERY build -- so its sources are the live ones. The package half belongs with
the commit that makes package resolution lazy, when a recipe is what resolution
consults.

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
from myflopy.specs import DataSourceSpec, GridSpec

logger = get_logger(__name__)


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

    def take(self, value: Any) -> Any:
        """Copy one source in and return the path to store, or the original."""

        if value is None:
            return value
        source = self._resolve(value)
        if not source.exists():
            # Not this verb's complaint to make: the resolver already names a
            # missing source, and refusing to save a project because one file moved
            # would lose the rest of it.
            self.missing.append(str(source))
            return value
        try:
            already_inside = source.resolve().is_relative_to(self.inputs_dir.resolve())
        except (OSError, ValueError):
            already_inside = False
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


__all__ = ["bundle_grid_sources", "grid_sources"]

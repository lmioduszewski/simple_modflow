"""Serializable source specifications for project-backed model recipes."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, ClassVar


def _path_text(path: Path | str) -> str:
    """Serialize a path as a forward-slash string for portable JSON/YAML output."""

    return str(Path(path).as_posix())


@dataclass(frozen=True, slots=True)
class DataSourceSpec:
    """Durable, serializable reference to a file that feeds a grid or package.

    Source specs are how a :class:`~myflopy.specs.SimulationSpec` records *where*
    its input data lives without inlining the data itself, so a recipe can be
    written to JSON/YAML, version-controlled, and resolved again later. Each
    subclass adds the fields a given file type needs (a table name, a CRS, a
    raster band, ...); this base captures the parts they all share -- the file
    ``path``, whether it is ``external`` to the project, and free-form
    ``metadata``. ``to_dict``/``from_dict`` round-trip a spec through its
    ``kind`` tag, dispatching to the right subclass on the way back.

    You rarely build these by hand: the project library helpers
    (``project.add_package``/``add_grid``) and the GIS-aware package forms
    (``mf.ghb.gpkg`` etc.) construct the appropriate subclass for you. Reach for
    them directly only when authoring a serialized recipe.

    Parameters
    ----------
    path : pathlib.Path or str
        Filesystem path to the source file. Coerced to a :class:`~pathlib.Path`
        and written with forward slashes by ``to_dict``. A RELATIVE path is
        anchored to the project root by the code that reads it -- grid resolution
        (:meth:`~myflopy.specs.GridSpec.resolve`) and, for a lazy package's recipe
        in ``PackageSpec.inputs``, ``ModelSpec.build`` -- unless ``external`` is
        set. An absolute path is used as it is.
    external : bool, default False
        ``True`` means "never anchor": a relative ``path`` is not joined to the
        project root (so it resolves against the current working directory),
        and grid resolution skips its existence check for this source.
        ``project.save(bundle_sources=True)`` still copies the file into the
        project and clears the flag, because once the copy is inside the project
        the path is anchored like any other. Written by ``to_dict`` only when
        ``True``.
    metadata : dict, optional
        Free-form JSON-serializable annotations carried with the reference, and
        written by ``to_dict`` only when non-empty. Grid resolution reads two
        kinds of key from it: ``metadata["crs"]`` as a last-resort CRS for a file
        that declares none, and, for a :class:`TableSource`,
        ``metadata["x_field"]``/``metadata["y_field"]`` naming the coordinate
        columns. Anything else is carried, not read.

    See Also
    --------
    TableSource, ShapeSource, GeoPackageSourceSpec, RasterSource, LiteralSource
    """

    path: Path | str
    external: bool = False
    metadata: dict[str, Any] = field(default_factory=dict)

    kind: ClassVar[str] = "DataSourceSpec"

    def __post_init__(self) -> None:
        """Coerce ``path`` to a :class:`~pathlib.Path` (on this frozen dataclass)."""

        object.__setattr__(self, "path", Path(self.path))

    def to_dict(self) -> dict[str, Any]:
        """Serialize to a ``kind``-tagged dict, omitting default ``external``/``metadata``."""

        payload: dict[str, Any] = {
            "kind": self.kind,
            "path": _path_text(self.path),
        }
        if self.external:
            payload["external"] = True
        if self.metadata:
            payload["metadata"] = self.metadata
        return payload

    def digest_payload(self) -> dict[str, Any]:
        """What :func:`myflopy.specs.package_digest` fingerprints instead of the object.

        The serialized declaration, which for every source but one IS ``to_dict()``.
        The hook exists so a source that knows its own file's CONTENT can be
        fingerprinted by that rather than by where the file sits -- see
        :meth:`BoundarySource.digest_payload`. Separate from ``to_dict`` on purpose:
        ``to_dict`` is what gets WRITTEN and must still name a path a reader can open.
        """

        return self.to_dict()

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> DataSourceSpec:
        """Rebuild a source spec, dispatching on its ``kind`` tag to the right subclass.

        Parameters
        ----------
        data : dict
            A payload as ``to_dict`` writes it. Its ``"kind"`` value picks the
            class to rebuild -- ``data``'s kind, not the class this is called on, so
            ``TableSource.from_dict`` of a ``"ShapeSource"`` payload returns a
            :class:`ShapeSource` -- and may name any registered source, a
            :class:`LiteralSource` included. The rest of the keys are that class's
            fields; ``"path"`` is required by every file-backed kind (``KeyError``
            without it), and keys a class does not read are ignored.

        Returns
        -------
        DataSourceSpec or LiteralSource
            A new instance of the class the ``kind`` tag names.

        Raises
        ------
        ValueError
            ``data["kind"]`` is missing or names no registered source class.
        """

        kind = str(data.get("kind"))
        try:
            source_type = _SOURCE_TYPES[kind]
        except KeyError as error:
            raise ValueError(f"Unknown data source kind: {kind!r}") from error
        return source_type._from_dict(data)

    @classmethod
    def _from_dict(cls, data: dict[str, Any]) -> DataSourceSpec:
        """Construct this concrete class from an already-dispatched payload dict."""

        return cls(
            path=data["path"],
            external=bool(data.get("external", False)),
            metadata=dict(data.get("metadata", {})),
        )


@dataclass(frozen=True, slots=True)
class TableSource(DataSourceSpec):
    """Reference to a tabular file, as a serializable recipe entry.

    Meant for row-oriented inputs -- stress-period data, well lists, observation
    tables. **What reads one today is grid resolution only**: a ``TableSource``
    passed as a :meth:`~myflopy.specs.GridSpec.voronoi` source is read with
    :func:`pandas.read_csv` (so it must be a CSV there, whatever the suffix),
    taking geometry from a WKT ``geometry`` column or from x/y columns.
    Placed in ``PackageSpec.inputs`` it reaches the package builder unresolved --
    no myflopy builder reads one.

    Parameters
    ----------
    path : pathlib.Path or str
        The table file. Coerced to a :class:`~pathlib.Path` and anchored like any
        :class:`DataSourceSpec` path.
    external : bool, default False
        As :class:`DataSourceSpec`: ``True`` leaves a relative path unanchored.
    metadata : dict, optional
        As :class:`DataSourceSpec`. For grid resolution of a table with no
        ``geometry`` column, ``metadata["x_field"]``/``metadata["y_field"]`` name
        the coordinate columns (default ``"x"``/``"y"``; either missing from the
        table raises ``ValueError``). With or without that column,
        ``metadata["crs"]`` is the CRS the rows are taken to be in.
    table : str, optional
        A sheet or table name inside a multi-table workbook. Recorded and
        serialized only: no reader consults it today (grid resolution reads the
        whole file with ``pandas.read_csv``).

    Examples
    --------
    >>> TableSource("wells.csv")
    >>> TableSource("inputs.xlsx", table="ghb_stress")
    """

    table: str | None = None
    kind: ClassVar[str] = "TableSource"

    def to_dict(self) -> dict[str, Any]:
        """Serialize like the base spec, adding ``table`` when set."""

        payload = DataSourceSpec.to_dict(self)
        if self.table is not None:
            payload["table"] = self.table
        return payload

    @classmethod
    def _from_dict(cls, data: dict[str, Any]) -> TableSource:
        """Rebuild a :class:`TableSource` from a payload dict."""

        return cls(
            path=data["path"],
            table=data.get("table"),
            external=bool(data.get("external", False)),
            metadata=dict(data.get("metadata", {})),
        )


@dataclass(frozen=True, slots=True)
class ShapeSource(DataSourceSpec):
    """Reference to a vector file (shapefile/GeoJSON) used for grid or BC geometry.

    Use this for single-layer vector files -- a boundary polygon for grid
    generation, refinement lines, or features that drive a list boundary
    condition. For a multi-layer GeoPackage, use :class:`GeoPackageSourceSpec`
    instead (it can name a layer and carry a field mapping).

    Pass everything after ``path`` **by keyword**: the inherited fields come
    first in positional order -- ``(path, external, metadata, crs)`` -- so a
    second positional argument sets ``external``, not ``crs``.

    Parameters
    ----------
    path : Path or str
        Path to the ``.shp`` / ``.geojson`` file. Relative paths are anchored at
        resolve time to the caller's ``project_root`` unless ``external`` is set.
    external : bool, default False
        ``True`` marks the file as outside the project tree: referenced in
        place, and its relative path is not anchored to ``project_root``.
    metadata : dict, optional
        Free-form JSON-serializable annotations. ``metadata["crs"]`` doubles as
        a last-resort CRS fallback.
    crs : str or int, optional
        A CRS to **assume when the file declares none** -- a fallback, not an
        override, and not a reprojection target. A file carrying its own CRS
        keeps it and this is ignored; grid sources are reprojected to the CRS on
        the :class:`~myflopy.specs.GridSpec`. Either spelling works:
        ``"EPSG:2927"`` or ``2927``.

    Examples
    --------
    >>> ShapeSource("domain_boundary.shp")
    >>> ShapeSource("streams.geojson", crs="EPSG:2927")   # .geojson has no CRS block

    See Also
    --------
    GeoPackageSourceSpec : Multi-layer GeoPackage, with a layer name and field map.
    """

    crs: str | int | None = None
    kind: ClassVar[str] = "ShapeSource"

    def to_dict(self) -> dict[str, Any]:
        """Serialize like the base spec, adding ``crs`` when set."""

        payload = DataSourceSpec.to_dict(self)
        if self.crs is not None:
            payload["crs"] = self.crs
        return payload

    @classmethod
    def _from_dict(cls, data: dict[str, Any]) -> ShapeSource:
        """Rebuild a :class:`ShapeSource` from a payload dict."""

        return cls(
            path=data["path"],
            crs=data.get("crs"),
            external=bool(data.get("external", False)),
            metadata=dict(data.get("metadata", {})),
        )


@dataclass(frozen=True, slots=True)
class GeoPackageSourceSpec(DataSourceSpec):
    """Serializable reference to one layer inside a GeoPackage (``.gpkg``).

    A ``.gpkg`` can hold many layers, so a bare path is ambiguous in a way a
    shapefile path is not. This spec pins **which** layer, **which** of its
    columns carry the values a reader needs, and optionally a subset of its
    rows -- all as plain data, so a recipe can be written to JSON/YAML and
    resolved again later.

    Use it wherever a source may be a multi-layer GeoPackage or needs a column
    mapping; use :class:`ShapeSource` for a single-layer vector file that needs
    neither. Today the only reader is **grid resolution** -- the ``boundary`` /
    ``refinement`` / ``breaklines`` / ``points`` sources of
    :meth:`~myflopy.specs.GridSpec.voronoi`. The runtime
    :class:`~myflopy.geopackage.GeoPackageSource` behind ``mf.ghb.gpkg(...)``
    is a separate class and does not read any of these fields; it names its
    columns with ``name_field`` / ``layer_field`` / ``period_field``.

    .. warning::
       Pass everything after ``path`` **by keyword**. This is a dataclass
       extending :class:`DataSourceSpec`, so the inherited fields come first in
       positional order -- ``(path, external, metadata, layer, query, crs,
       fields)``. ``GeoPackageSourceSpec("x.gpkg", "areas")`` therefore sets
       ``external="areas"``, which is truthy, rather than naming a layer.

    Parameters
    ----------
    path : Path or str
        Path to the ``.gpkg`` file. Relative paths are anchored at resolve time
        to the ``project_root`` given to
        :meth:`~myflopy.specs.GridSpec.resolve`, unless ``external`` is set.
        Coerced to a :class:`~pathlib.Path`.
    external : bool, default False
        ``True`` marks the file as living outside the project tree: it is
        referenced in place and its relative path is **not** anchored to
        ``project_root``.
    metadata : dict, optional
        Free-form JSON-serializable annotations carried with the reference.
        ``metadata["crs"]`` doubles as a last-resort CRS fallback.
    layer : str, optional
        Which layer to read. ``None`` reads the file's first/default layer,
        which is fine for a single-layer GeoPackage and a silent source of
        surprise in a multi-layer one -- name it explicitly when the file has
        more than one. The same file may be referenced by several specs with
        different ``layer`` values, which is how one GeoPackage supplies both
        refinement polygons and breakline centerlines.
    query : str, optional
        Row filter, applied **after** the layer is read. This is a
        :meth:`pandas.DataFrame.query` expression, **not** SQL: write
        ``"active == 1"``, not ``"active = 1"`` (the latter raises
        ``ValueError: cannot assign without a target object``). Column names
        are bare identifiers; quote strings, e.g. ``"kind == 'stream'"``.
        Because it filters after reading, it subsets features but does not
        reduce IO.
    crs : str or int, optional
        A CRS to **assume when the file declares none** -- a fallback, not an
        override. A GeoPackage that carries its own CRS keeps it, and passing
        something different here does nothing rather than reinterpreting the
        coordinates. Accepts either spelling: ``"EPSG:2927"`` or ``2927``.

        This is not the reprojection target either. Grid sources are reprojected
        to the CRS declared on the :class:`~myflopy.specs.GridSpec`, so a layer
        in EPSG:4326 inside a spec built with ``crs=2927`` is converted for you
        and ``crs=`` on this spec is not what does it. ``metadata["crs"]`` is
        consulted as a further fallback.
    fields : dict of str to str, optional
        Mapping of **logical key -> the column in this layer that holds it**.
        The KEY is the name myflopy looks up; the VALUE is your column. So
        ``{"area": "max_area"}`` reads each feature's target cell area from a
        column called ``max_area``.

        Grid sources understand exactly three keys, all optional:

        ``"area"``
            Target cell area for that feature, in CRS units squared. Falls back
            to ``refinement_max_area`` / ``breakline_max_area``; one of the two
            is required or resolution raises naming both.
        ``"label"``
            Region name used in diagnostics and in
            ``tri._prepared_regions``. Defaults to ``refinement_<i>`` /
            ``breakline_<i>``.
        ``"priority"``
            Who claims the overlap where regions cross. Defaults to
            ``refinement_priority`` (0) or ``breakline_priority`` (1).

        An unrecognised key is **ignored silently**, so a typo degrades to the
        global default rather than raising -- check the spelling here first when
        a per-feature value seems not to apply. A row whose mapped column is
        null or missing also falls back.

    Notes
    -----
    Serialization writes ``kind`` as ``"GeoPackageSource"`` -- the tag, not the
    class name -- and omits ``layer`` / ``query`` / ``crs`` / ``fields`` when
    unset, so ``to_dict``/``from_dict`` round-trip a minimal payload.

    Frozen and slotted: build a modified copy with
    :func:`dataclasses.replace`, not by assignment.

    Examples
    --------
    Refinement polygons carrying their own cell sizes:

    >>> GeoPackageSourceSpec("grid_refinement.gpkg", layer="areas",
    ...                      fields={"area": "max_area", "label": "name"})

    Two layers of one file feeding two different grid channels:

    >>> refine = GeoPackageSourceSpec("refine.gpkg", layer="areas",
    ...                               fields={"area": "max_area"})
    >>> creeks = GeoPackageSourceSpec("refine.gpkg", layer="centerlines")

    A subset of one layer, filtered with pandas syntax:

    >>> GeoPackageSourceSpec("bcs.gpkg", layer="boundaries",
    ...                      query="kind == 'ghb' and active == 1")

    A layer written without a CRS, told what it is:

    >>> GeoPackageSourceSpec("no_crs.gpkg", crs=2927)

    See Also
    --------
    ShapeSource : Single-layer vector file; no layer name, no field map.
    DataSourceSpec : The shared ``path``/``external``/``metadata`` base.
    myflopy.specs.GridSpec.voronoi : Where these sources are consumed.
    """

    layer: str | None = None
    query: str | None = None
    crs: str | int | None = None
    fields: dict[str, str] = field(default_factory=dict)
    kind: ClassVar[str] = "GeoPackageSource"

    def to_dict(self) -> dict[str, Any]:
        """Serialize like the base spec, adding ``layer``/``query``/``crs``/``fields`` when set."""

        payload = DataSourceSpec.to_dict(self)
        if self.layer is not None:
            payload["layer"] = self.layer
        if self.query is not None:
            payload["query"] = self.query
        if self.crs is not None:
            payload["crs"] = self.crs
        if self.fields:
            payload["fields"] = self.fields
        return payload

    @classmethod
    def _from_dict(cls, data: dict[str, Any]) -> GeoPackageSourceSpec:
        """Rebuild a :class:`GeoPackageSourceSpec` from a payload dict."""

        return cls(
            path=data["path"],
            layer=data.get("layer"),
            query=data.get("query"),
            crs=data.get("crs"),
            fields=dict(data.get("fields", {})),
            external=bool(data.get("external", False)),
            metadata=dict(data.get("metadata", {})),
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class BoundarySource(DataSourceSpec):
    """How to resolve one list boundary condition from a GeoPackage, as plain data.

    This is the RECIPE behind a lazy ``mf.<pkg>.gpkg(...)``: everything needed to
    reproduce the records given a grid, and nothing that needs one. It rides in
    ``PackageSpec.inputs`` under the option it supplies (``stress_period_data``), which
    is the shape ledger 272 calls case (b) -- a source shadowing the option of the same
    name, and the only shape ``_check_input_collisions`` permits there.

    It is NOT :class:`~myflopy.geopackage.GeoPackageSource`, the runtime resolver. That
    class holds a live :class:`~myflopy.specs.ModelContext` (hence a grid) and a cached
    GeoDataFrame, has no ``to_dict``, and raises on a gridless context -- so it can
    never be serialized, and its digest moves the moment anything reads its cache. This
    carries the same fields minus the context, plus what the helper adds on top.

    It is also not :class:`GeoPackageSourceSpec`, whose ``fields`` maps logical keys to
    columns for GRID resolution. The value mapping here is ``values`` rather than
    ``fields`` deliberately: one name meaning two things is how that class's own
    docstring says a recipe goes silently inert.

    Nothing is validated when one is constructed. The ``mf.<pkg>.gpkg`` helpers check
    the file, ``nper``, ``gaps`` and ``snap_within`` at declaration; everything is
    checked again when the recipe is resolved at BUILD, where
    :class:`~myflopy.geopackage.GeoPackageSource` is constructed from these fields and
    its ``<package>`` method called with ``values``.

    Parameters
    ----------
    path : pathlib.Path or str
        The GeoPackage (or any vector file ``geopandas.read_file`` opens) holding the
        features. Coerced to a :class:`~pathlib.Path`. A lazy ``mf.<pkg>.gpkg``
        declaration (no ``context=``) stores it RESOLVED to an absolute path after
        checking it exists; the recipe an eager declaration records in its metadata
        keeps the path as given. A relative one (what
        ``project.save(bundle_sources=True)`` writes) is anchored to the project root
        at build unless ``external`` is set. A missing file raises
        ``FileNotFoundError`` at build.
    external : bool, default False
        As :class:`DataSourceSpec`: ``True`` leaves a relative ``path`` unanchored at
        build. The helpers never set it.
    metadata : dict, optional
        As :class:`DataSourceSpec`: free-form annotations, carried not read. (The
        package's provenance -- ``source_path``, the content hash -- goes in the
        ``PackageSpec``'s own ``metadata``, not here.)
    package : str
        Which boundary condition -- ``"chd"``, ``"ghb"``, ``"drn"``, ``"riv"``,
        ``"wel"``, ``"rch"``, ``"evt"``. It decides the record SCHEMA (the field names
        and their MF6 order), which the package registry owns, and which
        ``GeoPackageSource`` method resolves the features.
    nper : int
        Number of stress periods to resolve; every period ``0..nper-1`` gets records.
        Below 1 raises ``ValueError`` (at declaration through a helper, else at build).
    values : dict, optional
        Per-field value specifications, keyed by the helper's own argument names
        (``head=``, ``conductance=``, ``elevation=``...) and passed as keyword
        arguments to that ``GeoPackageSource`` method. A value may be a column name,
        a constant, a list of column names for per-period data, or a live
        :class:`~myflopy.geopackage.Spread` / :class:`~myflopy.geopackage.CellSurfaceOffset`.
        Empty (the default) leaves every field at the method's own default column.
    layer : str, optional
        The layer/table inside ``path`` to read. ``None`` (the default) reads the
        file's first layer.
    name_field : str or None, default "name"
        Attribute column naming each feature: its value identifies the feature in
        warnings and errors (the row index does when the column is absent or this is
        ``None``), and becomes each record's boundname when ``boundnames`` is
        ``True`` -- the column is then required, and a file without it raises
        ``ValueError`` naming it. ``None`` with ``boundnames=True`` writes the row
        index as the boundname.
    layer_field : str or None, default "layer"
        Attribute column giving each feature's model layer, numbered from
        ``layer_base``; the column is required when set. ``None`` puts every feature
        in the first model layer. A value that is not a whole number at or above
        ``layer_base`` (a NULL, a fraction, one below the base) raises ``ValueError``
        naming the feature.
    period_field : str or None, optional
        Attribute column giving the ONE stress period a feature applies to, numbered
        from ``period_base``. A value that is not a whole number naming a period
        ``0..nper-1`` (a NULL, a fraction, one below the base or past the last period)
        raises ``ValueError`` naming the feature, as ``layer_field`` does.
        ``None`` (the default) applies every feature to all ``nper`` periods.
    layer_base : int, default 1
        What the first model layer is called in ``layer_field``: ``1`` for the usual
        GIS numbering, ``0`` for a zero-based column.
    period_base : int, default 0
        What the first stress period is called in ``period_field``.
    gaps : {"warn", "nearest", "raise", "ignore"}, default "warn"
        What a feature that reaches no active cell means (ledger 288): ``"warn"``
        leaves it out and names it in a ``UserWarning`` -- the same message also
        reports features clipped at the grid's edge, elevations below their own cell
        bottom and records lifted by ``min_above_bottom``; ``"nearest"`` also moves a
        POINT that missed onto the nearest usable cell when it lies within
        ``snap_within``; ``"raise"`` raises ``ValueError`` wherever ``"warn"`` would
        report; ``"ignore"`` says nothing. Whatever the mode, a source with NO
        feature on an active cell raises ``ValueError``. Another word raises
        ``ValueError``.
    snap_within : float, optional
        How far, in the grid's own coordinate units, ``gaps="nearest"`` may move a
        point -- measured to the candidate cell's POLYGON. Must be positive. ``None``
        (the default) means one cell width, ``sqrt(area)`` of the cell the point would
        snap to, resolved at BUILD: a lazy declaration has no mesh, and a distance
        resolved at declaration would carry the OLD grid's cell width across
        ``with_grid``. Only read under ``gaps="nearest"``.
    min_above_bottom : float, optional
        Keep each resolved ELEVATION field at least this far above its own cell's
        bottom, in model length units (ledger 290) -- the head of ``chd``/``ghb``,
        the elevation of ``drn``, the surface of ``evt``, and on ``riv`` the ``rbot``,
        with ``stage`` raised by the same amount when it was below the floor too.
        ``0.0`` still lifts by a tiny nudge, because a value written exactly at the
        bottom can round below it. ``None`` (the default) changes nothing but still
        reports records below their cell bottom. Both need layer elevations: without
        any, the report is skipped and a value here raises ``ValueError`` at build
        (ledger 395). Negative raises ``ValueError`` at build. The
        ``wel``/``rch``/``evt`` helpers never set it; ``wel`` and ``rch`` have no
        elevation field.
    edges_only : bool, default False
        Keep only cells on the edge of the ACTIVE domain, per layer (ledger 246). Read
        only for ``chd``/``ghb``/``drn``/``riv``, the packages whose registry entry
        allows it; for ``wel``/``rch``/``evt`` it is recorded and ignored.
    boundnames : bool, default True
        Whether each record carries a boundname (from ``name_field``) and the package
        is written with ``BOUNDNAMES``.
    source_sha : str, optional
        16-character content hash of ``path`` taken at declaration -- the helpers fill
        it in, ``None`` when the file could not be read. The DIGEST reads it in place
        of the path (:meth:`digest_payload`); it is not compared with the file when the
        recipe is built.

    Notes
    -----
    **Keyword-only, and not merely as a courtesy.** ``nper`` has no default and the
    inherited fields do, so without ``kw_only`` the class fails to be created at all
    ("non-default argument follows default argument"). It also closes the positional
    trap :class:`GeoPackageSourceSpec`'s docstring warns about.

    **``to_dict()`` is not pure JSON.** A live ``Spread`` in ``values`` comes back as
    itself, because that is what ``_canonical`` must fingerprint and what
    ``from_dict`` must rebuild. :func:`myflopy.specs._json_value` is the converter --
    it recurses into this payload and writes the ``$spread`` / ``$cell_surface_offset``
    tokens. ``json.dumps(source.to_dict())`` raises; that is not a defect.

    **It revives its own ``values``, so the trip back works from ANYWHERE.**
    ``_spec_value(..., int_keys=True)`` recurses into a source payload on the way to
    ``options``/``inputs``, and for a long time that was the only path that rebuilt the
    nested objects -- a recipe reached through any other field came back with plain
    dicts where the ``Spread`` was, silently. That stopped being theoretical when an
    EAGER declaration began recording its recipe in ``metadata`` (ledger 296), which is
    a different path; ``_from_dict`` now revives ``values`` itself and is idempotent, so
    arriving pre-revived through ``inputs`` costs nothing.

    **``kind`` ends in "Source" on purpose.** ``_spec_value`` keeps an
    ``endswith("Source")`` arm alongside the registry so an unregistered kind RAISES
    rather than degrading to a dict -- the silent failure of ledger 269. A new source
    class that forgets :data:`_SOURCE_TYPES` fails loudly.

    **It carries no package NAME.** Identity belongs to the ``PackageSpec``, not to its
    recipe: ``fork_package`` renames through ``_renamed_options``, which rewrites
    ``options`` and never touches ``inputs``, so a name frozen in here would be
    silently discarded by a rename -- and two renamed forks would collapse into one MF6
    package under a FloPy "Replacing existing package" warning, which is the exact
    defect ``_renamed_options`` already exists to prevent.
    """

    package: str
    nper: int
    values: dict[str, Any] = field(default_factory=dict)
    layer: str | None = None
    name_field: str | None = "name"
    layer_field: str | None = "layer"
    period_field: str | None = None
    layer_base: int = 1
    period_base: int = 0
    gaps: str = "warn"
    snap_within: float | None = None
    min_above_bottom: float | None = None
    edges_only: bool = False
    boundnames: bool = True
    #: 16-character content hash of ``path`` when the recipe was declared, or ``None``
    #: if the file could not be read. It is what the DIGEST sees in place of the path
    #: (:meth:`digest_payload`), which is why it lives here and not in ``metadata``:
    #: ``package_digest`` excludes metadata by contract.
    source_sha: str | None = None

    kind: ClassVar[str] = "BoundarySource"

    def to_dict(self) -> dict[str, Any]:
        """Serialize like the base spec, adding every field that is not a default.

        Defaults are omitted so that adding a defaulted field later cannot move an
        existing package's digest -- ``_canonical`` fingerprints this payload rather
        than the object (ledger 273).
        """

        payload = DataSourceSpec.to_dict(self)
        payload["package"] = self.package
        payload["nper"] = self.nper
        if self.values:
            payload["values"] = dict(self.values)
        from myflopy.specs import _json_value

        for name, default in (
            ("layer", None),
            ("name_field", "name"),
            ("layer_field", "layer"),
            ("period_field", None),
            ("layer_base", 1),
            ("period_base", 0),
            ("gaps", "warn"),
            ("snap_within", None),
            ("min_above_bottom", None),
            ("edges_only", False),
            ("boundnames", True),
            ("source_sha", None),
        ):
            # NORMALIZED first, then compared. The `type(...) is not type(...)` clause
            # is what keeps a numpy-typed value from being mistaken for its plain
            # default -- and it also wrote that numpy value straight into the payload,
            # which then failed BOTH `json.dumps` and `yaml.safe_dump`. `snap_within`
            # makes that live: `snap_within=np.sqrt(vor.gdf_vorPolys.area.median())` is
            # the obvious spelling, and it produced a recipe that could not be saved,
            # frozen or round-tripped (ledger 288). Normalizing above the comparison
            # also means a `np.float64` equal to its default is correctly OMITTED, so
            # no stored digest moves.
            value = _json_value(getattr(self, name))
            if value != default or type(value) is not type(default):
                payload[name] = value
        return payload

    def digest_payload(self) -> dict[str, Any]:
        """The declaration with the file's CONTENT standing in for its LOCATION.

        A lazy recipe is the only source a package's digest reads that is also the
        thing the package is BUILT FROM, and fingerprinting it by path was wrong in
        both directions at once (ledger 293):

        * **it churned.** ``project.save(bundle_sources=True)`` copies the GeoPackage
          into ``inputs/`` and repoints the recipe, so a byte-identical archive moved
          every lazy package's digest and ``runs.drift()`` reported every stored run as
          edited on the day it was archived -- measured 4c3988e4 -> 113e32d2 for a file
          that had not changed;
        * **and it was blind.** A DIFFERENT file at the same path digested IDENTICALLY,
          which is the ``unchanged``-is-the-dangerous-verdict failure: the recipe says
          "read this file", and nothing about the digest noticed when the file stopped
          saying what it said.

        The content hash fixes both, and it is also what makes the lazy form agree with
        the EAGER one. An eager ``.gpkg`` resolves at declaration and stores RECORDS, so
        its digest is location-blind and content-sensitive by construction. Before this
        the two drifted in OPPOSITE directions from the same declaration.

        Falls back to the path when ``source_sha`` is ``None`` -- a file that could not
        be read at declaration. Dropping the path as well would make two unreadable
        sources digest alike, which is the blindness above with no hash to replace it.
        """

        payload = BoundarySource.to_dict(self)
        if payload.get("source_sha") is not None:
            payload.pop("path", None)
        return payload

    @classmethod
    def _from_dict(cls, data: dict[str, Any]) -> BoundarySource:
        """Rebuild a :class:`BoundarySource` from a payload dict.

        ``values`` is revived HERE rather than relied on having been revived on the way
        in. Through ``inputs`` it already has been and this is a no-op; through any
        other field -- ``metadata``, where an eager declaration records its recipe --
        nothing had recursed into it, and a ``Spread`` came back a plain dict with
        nothing raised (ledger 296). A source that can only be rebuilt correctly via one
        of its callers is a trap, and this class's own docstring had recorded it as one.
        """

        from myflopy.specs import _spec_options

        values = data.get("values", {})
        if values:
            values = _spec_options(dict(values))
        return cls(
            path=data["path"],
            package=data["package"],
            nper=data["nper"],
            values=dict(values),
            layer=data.get("layer"),
            name_field=data.get("name_field", "name"),
            layer_field=data.get("layer_field", "layer"),
            period_field=data.get("period_field"),
            layer_base=data.get("layer_base", 1),
            period_base=data.get("period_base", 0),
            gaps=data.get("gaps", "warn"),
            snap_within=data.get("snap_within"),
            min_above_bottom=data.get("min_above_bottom"),
            edges_only=bool(data.get("edges_only", False)),
            boundnames=bool(data.get("boundnames", True)),
            source_sha=data.get("source_sha"),
            external=bool(data.get("external", False)),
            metadata=dict(data.get("metadata", {})),
        )


@dataclass(frozen=True, slots=True)
class RasterSource(DataSourceSpec):
    """Reference to a raster (GeoTIFF/IMG), as a serializable recipe entry.

    It records a raster and how it was meant to be sampled -- band, target and
    method -- and round-trips through ``to_dict``/``from_dict``. **Nothing in
    myflopy reads its band, target, method or CRS today**: no builder samples it; as a
    :meth:`~myflopy.specs.GridSpec.voronoi` boundary/refinement source the grid
    resolver raises ``TypeError`` (it reads GeoPackage, shape and table sources only);
    in the grid's path-valued ``rasters=``/``idomain_path=`` options only its anchored
    ``path`` is taken; and in ``PackageSpec.inputs`` it reaches the builder
    unresolved. To fill a per-cell array from a raster, use :class:`RasterValues`
    (``mf.RasterValues``); for an elevation, :meth:`myflopy.surfaces.Surface.raster`
    (``mf.Raster``) in a ``LayerStack``; to sample one by hand,
    ``vor.get_raster_vals_at_centroids``.

    Parameters
    ----------
    path : pathlib.Path or str
        The raster file (``.tif`` / ``.img`` ...). Coerced to a
        :class:`~pathlib.Path`.
    external : bool, default False
        As :class:`DataSourceSpec`: ``True`` leaves a relative path unanchored.
    metadata : dict, optional
        As :class:`DataSourceSpec`: free-form annotations, carried not read.
    band : int, optional
        The band meant to be read, 1-based. Recorded and serialized only.
    map_to : str, optional
        The package input or surface this raster is meant to populate. Recorded
        and serialized only.
    method : str, optional
        The sampling method meant for it (``"nearest"``, ``"linear"``, ...).
        Recorded and serialized only; nothing validates the name.
    crs : str or int, optional
        A CRS meant to override the raster's own. Recorded and serialized only.

    Examples
    --------
    >>> RasterSource("ground_surface.tif", map_to="top")
    >>> RasterSource("recharge_mm_yr.tif", band=1, method="linear")
    """

    band: int | None = None
    map_to: str | None = None
    method: str | None = None
    crs: str | int | None = None
    kind: ClassVar[str] = "RasterSource"

    def to_dict(self) -> dict[str, Any]:
        """Serialize like the base spec, adding any set ``band``/``map_to``/``method``/``crs``."""

        payload = DataSourceSpec.to_dict(self)
        for key in ("band", "map_to", "method", "crs"):
            value = getattr(self, key)
            if value is not None:
                payload[key] = value
        return payload

    @classmethod
    def _from_dict(cls, data: dict[str, Any]) -> RasterSource:
        """Rebuild a :class:`RasterSource` from a payload dict."""

        return cls(
            path=data["path"],
            band=data.get("band"),
            map_to=data.get("map_to"),
            method=data.get("method"),
            crs=data.get("crs"),
            external=bool(data.get("external", False)),
            metadata=dict(data.get("metadata", {})),
        )


@dataclass(frozen=True, slots=True)
class LiteralSource:
    """Inline literal value for a package input that has no external file.

    The escape hatch among the source specs: when an input is a plain scalar,
    list, or small mapping that you simply want to embed in a serialized recipe
    (rather than point at a file), wrap it in ``LiteralSource``. The ``value`` is
    stored as-is, so it must be JSON-serializable to round-trip. Unlike the
    file-backed specs it has no ``path`` and is not a :class:`DataSourceSpec`
    subclass.

    Parameters
    ----------
    value : object
        The literal payload, stored as it is. It must be JSON-serializable (after
        myflopy's numpy/tuple normalization) for a spec holding it to be written.
        Nothing unwraps it on the way to a build: placed in ``PackageSpec.inputs`` it
        reaches the builder as the ``LiteralSource`` object, not as ``value`` --
        ``inputs={"k": LiteralSource(1e-4)}`` on ``ModflowGwfnpf`` fails inside
        FloPy with an ``MFDataException`` on ``k``.
    metadata : dict, optional
        Free-form JSON-serializable annotations, written by ``to_dict`` only when
        non-empty.

    Examples
    --------
    >>> LiteralSource(1.0e-4)                 # a constant K value
    >>> LiteralSource({0: True})              # steady-state flags per period
    """

    value: Any
    metadata: dict[str, Any] = field(default_factory=dict)

    kind: ClassVar[str] = "LiteralSource"

    def to_dict(self) -> dict[str, Any]:
        """Serialize the inline value to a ``kind``-tagged dict (omitting empty metadata)."""

        payload = {
            "kind": self.kind,
            "value": self.value,
        }
        if self.metadata:
            payload["metadata"] = self.metadata
        return payload

    def digest_payload(self) -> dict[str, Any]:
        """The serialized declaration -- a literal has no file, so this IS ``to_dict``."""

        return self.to_dict()

    @classmethod
    def _from_dict(cls, data: dict[str, Any]) -> LiteralSource:
        """Rebuild a :class:`LiteralSource` from a payload dict."""

        return cls(
            value=data.get("value"),
            metadata=dict(data.get("metadata", {})),
        )


#: The polygon-to-cell rules :class:`Polygons` accepts, in the order they are taught.
#: Declared HERE, the lowest layer that needs them, so a declaration can be refused at
#: the call rather than five frames into a build. The resolver in
#: ``modflow/mf6/grid/geometry.py`` re-states the three names in its own dispatch rather
#: than importing this -- that module sits a layer ABOVE this one, so the import would
#: have to run the other way, and grid geometry has no business knowing about source
#: specs. ``tests/test_polygon_arrays.py`` asserts the two agree, which is what keeps
#: the restatement from drifting.
POLYGON_RULES: tuple[str, ...] = ("majority", "mean", "centroid")


@dataclass(frozen=True, slots=True, kw_only=True)
class Polygons(DataSourceSpec):
    """A polygon map that fills a per-cell ARRAY parameter, as plain data.

    The recipe behind ``mf.npf(k=mf.Polygons("geology.gpkg", column="k"))``: a vector
    file whose features carry the value in an attribute column, resolved against
    whatever grid the build supplies. It rides in ``PackageSpec.inputs`` under the
    parameter it fills, which is the shape ``_check_input_collisions`` permits there.

    **What this adds that the library did not already have.**
    :meth:`myflopy.geopackage.GeoPackageSource.k_array` has painted K from polygons for
    a long time, and it still does; it is not deprecated by this. What it cannot be is
    LAZY or SAVED. It needs the grid at DECLARATION, so a model built with it cannot be
    moved to another mesh -- the whole point of Phase 2 -- and it hands back an
    ``ndarray``, which ``to_dict`` refuses (ledger 193), so ``Project.save()`` persisted
    **zero simulations** for such a model and ``load()`` returned an empty project.
    Measured, both: see ledger 299. This carries the RECIPE instead, so the declaration
    survives a save and resolves against the grid it is built on.

    For a BOUNDARY CONDITION -- a conductance, a stage, an invert -- this is the wrong
    tool and is not needed: the seven ``mf.<pkg>.gpkg(...)`` forms already read a
    polygon's attribute columns, and have been lazy since ledger 279.

    A RASTER of values is :class:`RasterValues`, accepted in the same slots.

    Parameters
    ----------
    path : pathlib.Path or str
        Vector file holding the polygons (GeoPackage, shapefile, GeoJSON...), read with
        ``geopandas.read_file`` at BUILD and reprojected to the grid's CRS when both
        declare one. Coerced to a :class:`~pathlib.Path`. The array helpers
        (``mf.npf``/``mf.ic``/``mf.sto``) check it exists and store it RESOLVED to an
        absolute path at declaration -- a missing file raises ``FileNotFoundError``
        there.
    external : bool, default False
        As :class:`DataSourceSpec`: ``True`` leaves a relative ``path`` unanchored at
        build. The helpers never set it.
    metadata : dict, optional
        As :class:`DataSourceSpec`: free-form annotations, carried not read.
    column : str, optional
        Attribute column holding the value to paint. A NULL in it is REFUSED naming
        the feature, never filled with a default -- substituting a number the user did
        not give is its own silent guess (ledger 247). An empty string raises
        ``ValueError`` here; a column the file lacks raises ``ValueError`` at build,
        naming the file's columns. Give exactly one of ``column`` and ``value``:
        both, or neither, raises ``ValueError`` here.
    value : float, optional
        One number painted into every cell the polygons take, in place of a
        column -- for a constant factor that is not worth adding to the file, most
        often a MULTIPLIER in a product: ``mf.RasterValues("k33_3.tif") *
        mf.Polygons("pit.gpkg", value=0.1, default=1.0)`` is the raster times 0.1
        in the cells the polygons cover MORE than half of and times 1.0 in the rest
        (ledger 541; ledger 551) -- the ``rule`` and ``default`` decide
        which cells the polygons take. Must be a finite real number (``bool``
        refused), else ``ValueError`` here. Exactly one of ``column`` and
        ``value``.
    layer : str, optional
        Layer/table name inside a multi-layer file. ``None`` (the default) reads the
        file's first layer, which is right for a single-layer one.
    rule : {"majority", "mean", "centroid"}, default "majority"
        How a cell overlapping several polygons picks its value. Another word raises
        ``ValueError`` here, at the declaration.

        ``majority`` gives the cell the VALUE whose polygons cover most of it -- the
        features carrying one value are united first, so a unit drawn as several
        touching features (30% + 30% of a cell) beats another unit's single 40%, as
        its dissolved self does, and a tie between values goes to the value of the
        earliest feature among them in that cell (ledger 553). A map whose units are
        each one feature is judged exactly as by feature.
        It is the default because it INVENTS NOTHING: on a geology map an
        area-weighted mean across a contact writes a conductivity that no mapped unit
        has, which is the arithmetic-versus-harmonic problem in miniature.
        ``mean`` is the area-weighted mean, right for a genuinely continuous field.
        ``centroid`` matches :func:`~myflopy.modflow.mf6.grid.geometry.centroid_zones`
        and :meth:`myflopy.pest.ZoneSpec.from_polygons`. With ``default`` given, the
        part of a cell no polygon covers competes under ``majority`` and is averaged
        in under ``mean`` -- see ``default``.

        The fork is ``mean`` against the other two, and it is not subtle. On a 64-cell
        grid with ONE contact, ``mean`` wrote 12 distinct conductivities where the map
        holds 2 -- ten values no mapped unit has -- and across 19 contact positions its
        area-integral ran between 0.0% and 60.5% away from ``majority``'s.
        ``majority`` and ``centroid`` are close but not the same rule: over those same
        19 positions they differed on at most ONE cell of 64, the cell whose centre and
        whose bulk fall on opposite sides of the contact. Which is why this is NAMED
        rather than inherited.
    default : float, optional
        Value for the ground no polygon covers. ``None`` (the default): a cell no
        polygon reaches is REFUSED at build naming the count -- an unfilled cell is a
        NaN, and MF6 answers a NaN budget while reporting normal termination (ledger
        247) -- and a PARTLY covered cell takes the polygon(s) covering it, the
        uncovered rest of it no competitor.

        Given, it fills every cell no polygon reaches AND the uncovered part of each
        cell competes as the default (ledger 551, the user's call):
        under ``majority`` the default is one more value, whose vote is the part of
        the cell no polygon covers PLUS the cover of any polygon whose value EQUALS
        ``default`` (ledger 553); the value covering most of the cell takes it only
        when it covers MORE than that vote, so a single-valued map -- a multiplier --
        takes exactly the cells it covers more than half of, however its features
        split it (ledger 552); a tie with the default's vote (to 1e-9 of the cell's
        area) leaves the default. Under
        ``mean`` the uncovered part is averaged in AT the default, by area.
        ``centroid`` is unchanged: a cell whose centre is in no polygon reads the
        default. A map that tiles the grid, none of whose values equals ``default``,
        paints exactly what it painted without ``default``. Under every rule, with
        or without ``default``, a feature that reaches the grid but sets no cell --
        outvoted in every cell it touches, holding no centre under ``centroid``, or
        meeting the cells only along an edge -- is WARNED by name, with why, at build
        (ledger 585).
    layer_field : str, optional
        Column naming the MODEL layer a feature paints, ``layer_base``-based. ``None``
        (the default) makes the map areal, broadcast to EVERY layer -- which is what a
        geology map digitized in plan view means. When set, a file without the column
        raises ``ValueError`` at build, and so does a feature whose value is not a
        whole number from ``layer_base`` to ``layer_base + nlay - 1`` -- below the
        base, past the model's last layer, a fraction or a NULL -- naming every such
        feature (ledger 393).
    layer_base : int, default 1
        What the first model layer is called in ``layer_field``.
    source_sha : str, optional
        16-character content hash of ``path`` when the recipe was declared. It is what
        the DIGEST sees in place of the path (:meth:`digest_payload`), which is why it
        lives here and not in ``metadata``: ``package_digest`` excludes metadata by
        contract (ledger 293). The array helpers stamp it at declaration, replacing
        any value passed here.

    Notes
    -----
    **``kind`` is ``"PolygonsSource"`` while the class is ``Polygons``.** The class name
    is the one taught -- ``mf.Polygons`` -- and the kind must end in ``"Source"``
    because ``_spec_value`` keeps an ``endswith("Source")`` arm beside the registry so
    an UNREGISTERED kind raises rather than degrading to a plain dict (the silent
    failure of ledger 269). The registry maps the tag to the class; nothing requires
    them to be spelled alike, and :class:`DataSourceSpec`'s own ``"DataSourceSpec"``
    already is not.

    **It carries no package name**, for the reason :class:`BoundarySource` records: a
    name frozen in a recipe survives ``fork_package``'s rename and two renamed forks
    then collapse into one MF6 package.

    **``*`` builds a product.** ``polygons * factor`` (or ``factor * polygons``) is an
    :class:`ArrayProduct`, whose factors may be :class:`RasterValues`, ``Polygons`` and
    numbers; see that class for what a factor may not be.

    Examples
    --------
    >>> mf.npf(k=mf.Polygons("geology.gpkg", column="k", layer="units", default=1.0))
    >>> mf.ic(strt=mf.Polygons("zones.shp", column="head", rule="mean", default=9.0))
    >>> mf.Polygons("pit.gpkg", value=0.1, default=1.0)     # a constant multiplier
    """

    #: Numpy must not claim ``ndarray * Polygons`` element by element before this
    #: class's ``__rmul__`` refuses the array (an ARRAY factor is position-bound).
    __array_ufunc__: ClassVar[None] = None

    column: str | None = None
    value: float | None = None
    layer: str | None = None
    rule: str = "majority"
    default: float | None = None
    layer_field: str | None = None
    layer_base: int = 1
    source_sha: str | None = None

    kind: ClassVar[str] = "PolygonsSource"

    def __post_init__(self) -> None:
        """Coerce ``path`` and refuse an unknown ``rule`` at the DECLARATION."""

        DataSourceSpec.__post_init__(self)
        if self.rule not in POLYGON_RULES:
            raise ValueError(
                f"Polygons(rule={self.rule!r}) is not a rule. Use one of "
                + ", ".join(repr(name) for name in POLYGON_RULES)
                + ". 'majority' gives a cell the value whose polygons cover most of it, "
                "'mean' the area-weighted mean of them, 'centroid' the one holding the "
                "cell's centre."
            )
        if self.column is not None and self.value is not None:
            raise ValueError(
                f"Polygons(column={self.column!r}, value={self.value!r}) gives two "
                "values to paint. Give `column=` to read each polygon's own value from "
                "the file, or `value=` to paint one number into every polygon -- not "
                "both."
            )
        if self.value is not None:
            object.__setattr__(
                self, "value", _finite_number(self.value, "Polygons(value=...)")
            )
            return
        if not self.column:
            raise ValueError(
                "Polygons(column=...) names the attribute column holding the value to "
                "paint, and cannot be empty -- or give `value=` to paint one number "
                "into every polygon."
            )

    def to_dict(self) -> dict[str, Any]:
        """Serialize like the base spec, adding every field that is not a default.

        Defaults are omitted so that adding a defaulted field later cannot move an
        existing package's digest -- ``_canonical`` fingerprints this payload rather
        than the object (ledger 273).
        """

        payload = DataSourceSpec.to_dict(self)
        # `column` stays written FIRST and unconditionally when set, exactly as before
        # `value=` existed, so no stored Polygons digest moves (ledger 273's rule).
        if self.column is not None:
            payload["column"] = self.column
        from myflopy.specs import _json_value

        for name, default in (
            ("value", None),
            ("layer", None),
            ("rule", "majority"),
            ("default", None),
            ("layer_field", None),
            ("layer_base", 1),
            ("source_sha", None),
        ):
            # NORMALIZED first, then compared -- the `type(...) is not type(...)` clause
            # is what keeps a numpy-typed value from being mistaken for its plain
            # default and written into a payload that then fails both `json.dumps` and
            # `yaml.safe_dump`. `default=np.float64(...)` is an obvious spelling here,
            # exactly as `snap_within=` was for `BoundarySource` (ledger 288).
            value = _json_value(getattr(self, name))
            if value != default or type(value) is not type(default):
                payload[name] = value
        return payload

    def digest_payload(self) -> dict[str, Any]:
        """The declaration with the file's CONTENT standing in for its LOCATION.

        Ledger 293's rule, for the same two reasons in the same two directions: keying
        on the path CHURNS (``project.save(bundle_sources=True)`` repoints the recipe,
        so a byte-identical archive moves the digest and ``runs.drift()`` calls every
        stored run edited) and it is BLIND (a different file at the same path digests
        identically). Falls back to the path when the file could not be read -- dropping
        it as well would make two unreadable sources digest alike.
        """

        payload = Polygons.to_dict(self)
        if payload.get("source_sha") is not None:
            payload.pop("path", None)
        return payload

    def __mul__(self, other: Any) -> ArrayProduct:
        """``polygons * factor``: an :class:`ArrayProduct` of the two."""

        return ArrayProduct.of(self, other)

    def __rmul__(self, other: Any) -> ArrayProduct:
        """``factor * polygons``: an :class:`ArrayProduct` of the two."""

        return ArrayProduct.of(other, self)

    @classmethod
    def _from_dict(cls, data: dict[str, Any]) -> Polygons:
        """Rebuild a :class:`Polygons` from a payload dict."""

        return cls(
            path=data["path"],
            column=data.get("column"),
            value=data.get("value"),
            layer=data.get("layer"),
            rule=data.get("rule", "majority"),
            default=data.get("default"),
            layer_field=data.get("layer_field"),
            layer_base=data.get("layer_base", 1),
            source_sha=data.get("source_sha"),
            external=bool(data.get("external", False)),
            metadata=dict(data.get("metadata", {})),
        )


#: How :class:`RasterValues` reads a cell's value off a raster, and how it averages the
#: pixels inside one. Stated here, the lowest layer that needs them, so a declaration is
#: refused at the call; the sampler in ``modflow/mf6/grid/surfaces.py`` re-states them
#: (it sits a layer above) and ``tests/test_raster_arrays.py`` holds the two together.
RASTER_METHODS: tuple[str, ...] = ("area", "centroid")
RASTER_MEANS: tuple[str, ...] = ("arithmetic", "geometric", "harmonic")
#: What a cell the raster gives no value for becomes. Narrower than the sampler's four
#: words on purpose: ``"warn"`` and ``"ignore"`` would hand MF6 a NaN, which it answers
#: with a NaN budget under Normal termination (ledger 247).
RASTER_GAPS: tuple[str, ...] = ("raise", "nearest")


def _finite_number(value: Any, what: str) -> float:
    """``value`` as a finite ``float``, or ``ValueError`` naming ``what``.

    ``bool`` is refused although it subclasses ``int``: ``value=True`` is a typo for a
    flag, not the number 1.
    """

    import math
    import numbers

    if isinstance(value, bool) or type(value).__name__ == "bool_":
        raise ValueError(f"{what} must be a number, not {value!r}.")
    if not isinstance(value, numbers.Real):
        raise ValueError(f"{what} must be a real number, not {value!r} ({type(value).__name__}).")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{what} must be finite, not {value!r}.")
    return number


@dataclass(frozen=True, slots=True, kw_only=True)
class RasterValues(DataSourceSpec):
    """A raster that fills a per-cell ARRAY parameter, as plain data.

    The recipe behind ``mf.npf(k=mf.RasterValues("k_1.tif"))``: one raster band,
    sampled onto whatever grid the build supplies -- the same sampler, the same
    numbers, as ``mf.Surface.raster("k_1.tif").values(vor)``, but held as a
    DECLARATION. So a model built from it moves to another mesh (``with_grid`` reads
    ``rebuilt``), saves with its project (an ``ndarray`` is refused by ``to_dict``,
    ledger 193, so a model holding sampled arrays persisted no simulation at all), and
    its run freezes the arrays the build sampled (ledger 251). Its digest reads the
    raster's CONTENT, not its location (ledger 293), and
    ``project.save(bundle_sources=True)`` copies the raster into the project (ledger
    256/293). See ledger 541.

    It is NOT :class:`myflopy.surfaces.Surface` (``mf.Raster``), which is an
    ELEVATION: it carries length units and the relative kinds that only mean something
    for a layer surface, and the array slots refuse it on purpose. ``RasterValues`` is
    a property field -- a conductivity, a storage coefficient, a starting head, a
    recharge rate -- in whatever units the raster holds, written as it is.

    One object is ONE file. A layered parameter takes a per-layer LIST --
    ``k=[mf.RasterValues("k_1.tif"), ..., mf.RasterValues("k_5.tif")]`` -- and one
    recipe in a layered slot is broadcast to every layer. ``*`` builds an
    :class:`ArrayProduct` with another recipe or a number.

    Where it is accepted: ``mf.npf(k=, k22=, k33=)``, ``mf.ic(strt=)``,
    ``mf.sto(ss=, sy=)`` -- a recipe, a per-layer list of recipes and numbers, or a
    product -- and, AREAL, ``mf.rch.array(recharge=)`` and ``mf.evt.array(rate=,
    surface=, depth=)``, as a recipe or a ``{period: recipe}`` dict. Never in an
    integer flag (``icelltype``, ``iconvert``, ``irch``, ``ievt``).

    Parameters
    ----------
    path : pathlib.Path or str
        The raster file (GeoTIFF, IMG, anything ``rasterio`` opens). Coerced to a
        :class:`~pathlib.Path`. The array helpers check it exists and store it
        RESOLVED to an absolute path at declaration -- a missing file raises
        ``FileNotFoundError`` there. A raster in another CRS than the grid is sampled
        by reprojecting the cells to it; with either CRS missing it is sampled
        without reprojection and warned about.
    external : bool, default False
        As :class:`DataSourceSpec`: ``True`` leaves a relative ``path`` unanchored at
        build. The helpers never set it.
    metadata : dict, optional
        As :class:`DataSourceSpec`: free-form annotations, carried not read.
    band : int, default 1
        The 1-based band to read. Below 1, or not a whole number, raises
        ``ValueError`` here; a band the file lacks raises ``ValueError`` at build,
        naming how many it has.
    method : {"area", "centroid"}, default "area"
        How a cell's value is read off the raster. ``"area"`` averages every pixel
        whose CENTRE falls inside the cell polygon (``mean`` says how), and a cell
        holding no such pixel -- smaller than a pixel, or on the raster's edge --
        falls back to the pixel at its centroid. ``"centroid"`` reads that one pixel
        for every cell. Another word raises ``ValueError`` here.
    mean : {"arithmetic", "geometric", "harmonic"}, default "arithmetic"
        How ``method="area"`` averages the pixels inside one cell. ``"arithmetic"``
        is the default because it reproduces, bit for bit, the arrays
        ``Surface.raster(path).values(vor)`` samples -- so moving a model onto
        recipes changes no number. ``"geometric"`` (the exponential of the mean
        log) and ``"harmonic"`` (the reciprocal of the mean reciprocal) are for a
        conductivity whose pixels are much finer than the cells. **Which one is
        right for K depends on the direction of FLOW, which is a property of the
        solution, not of the field** (ledger 300): with the contrast across the
        flow the harmonic is nearly exact and the arithmetic far too high (264x on
        one measured coarsening), and turned 90 degrees the ordering reverses. The
        geometric is never best and never catastrophic. Both refuse a pixel that is
        zero or negative, at build, naming the count. Not ``"arithmetic"`` with
        ``method="centroid"`` raises ``ValueError`` here: one pixel has no mean.
    nodata : float, optional
        A value to treat as no-data, overriding the file's own declared one, for a
        raster that uses a sentinel (often 0 or -9999) without recording it.
        ``None`` (the default) uses the file's.
    gaps : {"raise", "nearest"}, default "raise"
        What a cell the raster gives no value for means -- off the raster, or on its
        no-data. ``"raise"`` refuses at build naming how many cells and which, unless
        ``default`` is given; ``"nearest"`` takes each from its nearest valued cell,
        by centroid distance, and logs which. The sampler's other two words,
        ``"warn"`` and ``"ignore"``, are refused here (``ValueError``): both hand MF6 a
        NaN, which it answers with a NaN budget under Normal termination (ledger
        247). Whatever the mode, a raster giving NO cell a value raises at build,
        printing both extents and both CRSs.
    default : float, optional
        The value for a cell the raster gives no value for, in place of refusing
        it. Must be a finite number. Giving it with ``gaps="nearest"`` raises
        ``ValueError`` here: they are two answers to one question.
    source_sha : str, optional
        16-character content hash of ``path`` when the recipe was declared. It is
        what the DIGEST sees in place of the path (:meth:`digest_payload`), so
        bundling a project moves no digest and editing the raster does (ledger 293).
        The array helpers stamp it at declaration, replacing any value passed here.

    Notes
    -----
    **``kind`` is ``"RasterValuesSource"``**, for the reason :class:`Polygons`
    records: the kind must end in ``"Source"`` so that an unregistered one raises.

    It sits beside :class:`RasterSource`, which RECORDS a raster for a grid recipe and
    is sampled by nothing; this one is resolved.

    Examples
    --------
    >>> mf.npf(k=[mf.RasterValues(f"k_{n}.tif") for n in range(1, 6)])
    >>> mf.sto(ss=mf.RasterValues("ss.tif", gaps="nearest"))
    >>> mf.npf(k=mf.RasterValues("k.tif", mean="geometric", default=1.0))
    >>> mf.rch.array(recharge={0: mf.RasterValues("recharge_01.tif"),
    ...                        1: mf.RasterValues("recharge_02.tif")})
    """

    #: See :attr:`Polygons.__array_ufunc__`.
    __array_ufunc__: ClassVar[None] = None

    band: int = 1
    method: str = "area"
    mean: str = "arithmetic"
    nodata: float | None = None
    gaps: str = "raise"
    default: float | None = None
    source_sha: str | None = None

    kind: ClassVar[str] = "RasterValuesSource"

    def __post_init__(self) -> None:
        """Coerce ``path`` and refuse what cannot be resolved, at the DECLARATION."""

        DataSourceSpec.__post_init__(self)
        import numbers

        band = self.band
        if isinstance(band, bool) or not isinstance(band, numbers.Integral):
            raise ValueError(f"RasterValues(band={band!r}) must be a whole number from 1.")
        object.__setattr__(self, "band", int(band))
        if self.band < 1:
            raise ValueError(
                f"RasterValues(band={self.band}) is not a band: bands are numbered from 1."
            )
        for name, allowed in (
            ("method", RASTER_METHODS),
            ("mean", RASTER_MEANS),
            ("gaps", RASTER_GAPS),
        ):
            if getattr(self, name) not in allowed:
                extra = (
                    " 'warn' and 'ignore' would hand MODFLOW a NaN for those cells, "
                    "which it answers with a NaN budget under Normal termination; give "
                    "`default=` to say what an uncovered cell means."
                    if name == "gaps"
                    else ""
                )
                raise ValueError(
                    f"RasterValues({name}={getattr(self, name)!r}) is not one of "
                    + ", ".join(repr(word) for word in allowed)
                    + "."
                    + extra
                )
        if self.method == "centroid" and self.mean != "arithmetic":
            raise ValueError(
                f"RasterValues(method='centroid', mean={self.mean!r}): a centroid read is "
                "ONE pixel, which has no mean. Use method='area' to average the pixels "
                "inside each cell."
            )
        if self.default is not None:
            object.__setattr__(
                self, "default", _finite_number(self.default, "RasterValues(default=...)")
            )
            if self.gaps == "nearest":
                raise ValueError(
                    "RasterValues(gaps='nearest', default=...) gives two answers to one "
                    "question -- what a cell the raster gives no value for means. Give "
                    "one."
                )
        if self.nodata is not None:
            object.__setattr__(
                self, "nodata", _finite_number(self.nodata, "RasterValues(nodata=...)")
            )

    def to_dict(self) -> dict[str, Any]:
        """Serialize like the base spec, adding every field that is not a default.

        Defaults are omitted so that adding a defaulted field later cannot move an
        existing package's digest (ledger 273).
        """

        payload = DataSourceSpec.to_dict(self)
        # No `_json_value` pass, unlike `Polygons.to_dict`: `__post_init__` already made
        # every field a plain Python value (`band` an `int`, `nodata`/`default` a
        # `float`), so a numpy scalar cannot reach the payload (ledger 265).
        for name, default in (
            ("band", 1),
            ("method", "area"),
            ("mean", "arithmetic"),
            ("nodata", None),
            ("gaps", "raise"),
            ("default", None),
            ("source_sha", None),
        ):
            value = getattr(self, name)
            if value != default or type(value) is not type(default):
                payload[name] = value
        return payload

    def digest_payload(self) -> dict[str, Any]:
        """The declaration with the file's CONTENT standing in for its LOCATION.

        :meth:`Polygons.digest_payload`'s rule (ledger 293): falls back to the path
        when the file could not be read.
        """

        payload = RasterValues.to_dict(self)
        if payload.get("source_sha") is not None:
            payload.pop("path", None)
        return payload

    def __mul__(self, other: Any) -> ArrayProduct:
        """``raster * factor``: an :class:`ArrayProduct` of the two."""

        return ArrayProduct.of(self, other)

    def __rmul__(self, other: Any) -> ArrayProduct:
        """``factor * raster``: an :class:`ArrayProduct` of the two."""

        return ArrayProduct.of(other, self)

    @classmethod
    def _from_dict(cls, data: dict[str, Any]) -> RasterValues:
        """Rebuild a :class:`RasterValues` from a payload dict."""

        return cls(
            path=data["path"],
            band=data.get("band", 1),
            method=data.get("method", "area"),
            mean=data.get("mean", "arithmetic"),
            nodata=data.get("nodata"),
            gaps=data.get("gaps", "raise"),
            default=data.get("default"),
            source_sha=data.get("source_sha"),
            external=bool(data.get("external", False)),
            metadata=dict(data.get("metadata", {})),
        )


@dataclass(frozen=True, slots=True)
class ArrayProduct:
    """The cell-by-cell PRODUCT of array recipes and numbers, as plain data.

    What ``recipe * factor`` builds: ``mf.RasterValues("k33_3.tif") *
    mf.Polygons("pit.gpkg", value=0.1, default=1.0)`` is the raster's K33, times 0.1
    in the cells the pit polygons cover more than half of and unchanged elsewhere -- a
    polygon MULTIPLIER without a hand-edited copy of the raster (ledger 541; the
    reach is ledger 551's, see :class:`Polygons` ``default``). Every factor is resolved
    on the grid the build supplies and the results multiplied, so the product moves,
    saves, digests and bundles exactly as its factors do.

    Accepted wherever :class:`RasterValues` is, and inside a per-layer list.

    Parameters
    ----------
    factors : tuple of (RasterValues or Polygons or float)
        What to multiply, in order. A factor that is itself a product is flattened
        into this one. At least one must be a recipe (``ValueError`` otherwise: a
        product of numbers is a number, so write it). A number must be finite. An
        ``ndarray`` is REFUSED (``TypeError``): it names cells by POSITION, which is
        exactly what makes a package unsaveable and unmovable -- the reason the
        recipe exists; ``with_grid(remap=True)`` carries a per-cell array to a new
        mesh. A ``Surface``/``LayerStack`` (an ELEVATION) and anything else are
        refused too, naming the type. Build it with ``*`` or :meth:`of`; the
        constructor checks the same.

    Notes
    -----
    Serialized under its own ``"$array_product"`` token, each factor written as the
    factor writes itself, so the walkers that repoint, bundle and anchor a source
    reach a factor's path unchanged. Its digest is its factors' digests -- content
    keyed -- in order.

    Examples
    --------
    >>> pit = mf.Polygons("pit.gpkg", value=0.1, default=1.0)
    >>> mf.npf(k33=[mf.RasterValues("k33_1.tif"), mf.RasterValues("k33_2.tif") * pit])
    >>> mf.ic(strt=0.5 * mf.RasterValues("head.tif"))
    """

    #: See :attr:`Polygons.__array_ufunc__`.
    __array_ufunc__: ClassVar[None] = None

    factors: tuple[Any, ...]

    #: The key it is written under; :mod:`myflopy.specs` dispatches on it.
    TOKEN: ClassVar[str] = "$array_product"

    def __post_init__(self) -> None:
        """Flatten nested products and refuse a factor that is not a recipe or a number."""

        flat: list[Any] = []
        for factor in tuple(self.factors):
            if isinstance(factor, ArrayProduct):
                flat.extend(factor.factors)
            elif isinstance(factor, (RasterValues, Polygons)):
                flat.append(factor)
            else:
                flat.append(_product_number(factor))
        if not any(isinstance(factor, (RasterValues, Polygons)) for factor in flat):
            raise ValueError(
                "An ArrayProduct needs at least one RasterValues or Polygons factor; a "
                "product of numbers is a number, so pass the number."
            )
        object.__setattr__(self, "factors", tuple(flat))

    @classmethod
    def of(cls, *factors: Any) -> ArrayProduct:
        """The product of ``factors``, as ``*`` builds it.

        Parameters
        ----------
        *factors : RasterValues, Polygons, ArrayProduct or float
            What to multiply, in order; a product among them is flattened in. See
            the class's ``factors`` for what is refused.

        Returns
        -------
        ArrayProduct
        """

        return cls(factors=tuple(factors))

    def __mul__(self, other: Any) -> ArrayProduct:
        """``product * factor``: one longer product."""

        return ArrayProduct.of(self, other)

    def __rmul__(self, other: Any) -> ArrayProduct:
        """``factor * product``: one longer product."""

        return ArrayProduct.of(other, self)

    def to_dict(self) -> dict[str, Any]:
        """``{"factors": [...]}``, each recipe as its own ``to_dict()`` writes it."""

        return {
            "factors": [
                factor.to_dict() if isinstance(factor, DataSourceSpec) else factor
                for factor in self.factors
            ]
        }

    def digest_payload(self) -> dict[str, Any]:
        """The factors' content-keyed payloads, in order (ledger 293)."""

        return {
            "factors": [
                factor.digest_payload() if isinstance(factor, DataSourceSpec) else factor
                for factor in self.factors
            ]
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ArrayProduct:
        """Rebuild a product from what :meth:`to_dict` wrote.

        Parameters
        ----------
        data : dict
            ``{"factors": [...]}``; a dict factor is revived through
            :func:`source_from_dict` (or is already a live recipe), a number is kept.

        Returns
        -------
        ArrayProduct
        """

        return cls(
            factors=tuple(
                source_from_dict(factor) if isinstance(factor, dict) else factor
                for factor in data["factors"]
            )
        )


def _product_number(factor: Any) -> float:
    """A number factor of a product, or the refusal naming what the factor is."""

    if getattr(factor, "ndim", 0) and hasattr(factor, "shape"):
        raise TypeError(
            "An array cannot be a factor of an array recipe: it names cells by "
            "POSITION, so the package would again be unsaveable and could not move to "
            "another mesh -- which is what the recipe exists to avoid. Put the values "
            "in a raster (mf.RasterValues) or a polygon map (mf.Polygons(value=...)); "
            "to carry a calibrated per-cell array to a new mesh use "
            "sim.with_grid(grid, remap=True)."
        )
    token = getattr(type(factor), "_SPEC_TOKEN", None)
    if token in ("$surface", "$layer_stack"):
        raise TypeError(
            "A Surface or LayerStack is an ELEVATION and cannot be a factor of an array "
            "recipe. A raster of values is mf.RasterValues(path)."
        )
    if isinstance(factor, (int, float)) or type(factor).__module__ == "numpy":
        return _finite_number(factor, "A number factor of an ArrayProduct")
    raise TypeError(
        f"{factor!r} ({type(factor).__name__}) cannot be a factor of an array recipe: a "
        "factor is mf.RasterValues, mf.Polygons, or a number."
    )


_SOURCE_TYPES: dict[str, type[Any]] = {
    DataSourceSpec.kind: DataSourceSpec,
    TableSource.kind: TableSource,
    ShapeSource.kind: ShapeSource,
    GeoPackageSourceSpec.kind: GeoPackageSourceSpec,
    BoundarySource.kind: BoundarySource,
    Polygons.kind: Polygons,
    RasterValues.kind: RasterValues,
    RasterSource.kind: RasterSource,
    LiteralSource.kind: LiteralSource,
}


def is_source_payload(data: Any) -> bool:
    """Whether ``data`` is a serialized source dict :func:`source_from_dict` rebuilds.

    The registry is the authority. A ``kind.endswith("Source")`` test looks
    equivalent and is not: :class:`DataSourceSpec`'s own kind is
    ``"DataSourceSpec"``, which ends in "Spec", so such a test silently passed
    the public base class through as a plain ``dict``.

    Parameters
    ----------
    data : object
        Anything. Only a ``dict`` whose ``"kind"`` value names a registered source
        class answers ``True``; the rest of its keys are not checked.

    Returns
    -------
    bool
    """

    return isinstance(data, dict) and str(data.get("kind")) in _SOURCE_TYPES


def source_from_dict(data: dict[str, Any]) -> DataSourceSpec | LiteralSource:
    """Recreate a source specification from its serialized representation.

    The same dispatch as :meth:`DataSourceSpec.from_dict`.

    Parameters
    ----------
    data : dict
        A payload as a source's ``to_dict`` writes it; its ``"kind"`` value picks the
        class, and ``"path"`` is required by every file-backed kind (``KeyError``
        without it).

    Returns
    -------
    DataSourceSpec or LiteralSource

    Raises
    ------
    ValueError
        ``data["kind"]`` is missing or names no registered source class.
    """

    kind = str(data.get("kind"))
    try:
        source_type = _SOURCE_TYPES[kind]
    except KeyError as error:
        raise ValueError(f"Unknown source kind: {kind!r}") from error
    return source_type._from_dict(data)


__all__ = [
    "ArrayProduct",
    "BoundarySource",
    "POLYGON_RULES",
    "Polygons",
    "RASTER_GAPS",
    "RASTER_MEANS",
    "RASTER_METHODS",
    "RasterValues",
    "DataSourceSpec",
    "is_source_payload",
    "GeoPackageSourceSpec",
    "LiteralSource",
    "RasterSource",
    "ShapeSource",
    "TableSource",
    "source_from_dict",
]

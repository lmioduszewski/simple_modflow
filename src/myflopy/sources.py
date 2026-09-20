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
    path
        Filesystem path to the source file. Stored as a :class:`~pathlib.Path`.
    external
        ``True`` marks the file as living outside the project tree (it is
        referenced in place rather than copied into the run workspace).
    metadata
        Arbitrary JSON-serializable annotations carried alongside the reference.

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
        """Rebuild a source spec, dispatching on its ``kind`` tag to the right subclass."""

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
    """Reference to a tabular file (CSV/Excel/Parquet) that supplies package records.

    Use this for stress-period data, well lists, observation tables, and similar
    row-oriented inputs. ``table`` names the sheet/table within a multi-table
    workbook (an Excel sheet name, for example) and may be omitted for a single
    flat file such as a CSV.

    Parameters
    ----------
    path
        Path to the ``.csv`` / ``.xlsx`` / ``.parquet`` file.
    table
        Optional sheet or table name inside a multi-table workbook.

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
    crs : str or int, optional
        A CRS to **assume when the file declares none** -- a fallback, not an
        override, and not a reprojection target. A file carrying its own CRS
        keeps it and this is ignored; grid sources are reprojected to the CRS on
        the :class:`~myflopy.specs.GridSpec`. Either spelling works:
        ``"EPSG:2927"`` or ``2927``.
    external : bool, default False
        ``True`` marks the file as outside the project tree: referenced in
        place, and its relative path is not anchored to ``project_root``.
    metadata : dict, optional
        Free-form JSON-serializable annotations. ``metadata["crs"]`` doubles as
        a last-resort CRS fallback.

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
    external : bool, default False
        ``True`` marks the file as living outside the project tree: it is
        referenced in place and its relative path is **not** anchored to
        ``project_root``.
    metadata : dict, optional
        Free-form JSON-serializable annotations carried with the reference.
        ``metadata["crs"]`` doubles as a last-resort CRS fallback.

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

    Parameters
    ----------
    package : str
        Which boundary condition -- ``"chd"``, ``"ghb"``, ``"drn"``, ``"riv"``,
        ``"wel"``, ``"rch"``, ``"evt"``. It decides the record SCHEMA (the field names
        and their MF6 order), which the package registry owns.
    nper : int
        Number of stress periods to resolve.
    values : dict
        Per-field value specifications, keyed by the helper's own argument names
        (``head=``, ``conductance=``, ``elevation=``...). A value may be a column name,
        a constant, a list of column names for per-period data, or a live
        :class:`~myflopy.geopackage.Spread` / :class:`~myflopy.geopackage.CellSurfaceOffset`.
    layer, name_field, layer_field, period_field, layer_base, period_base
        As :class:`~myflopy.geopackage.GeoPackageSource` reads them.
    gaps, snap_within
        What a feature that reaches no cell means, and how far a POINT may be moved to
        reach one (ledger 288). ``snap_within`` stays ``None`` in the recipe when the
        caller gave no number: a lazy declaration has no mesh, and a distance resolved
        at declaration would carry the OLD grid's cell width across ``with_grid``.
    edges_only : bool
        Keep only cells on the edge of the ACTIVE domain, per layer (ledger 246).
    boundnames : bool
        Whether to write a boundname column.

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

    **The reverse trip works only through ``options``/``inputs``.**
    ``_spec_value(..., int_keys=True)`` is what recurses INTO a source payload and
    revives those tokens before this class is rebuilt. A recipe reached through any
    other field comes back with plain dicts where the objects were.

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

        ``values`` is taken as it arrives: by the time this runs, ``_spec_value`` has
        already revived any ``$spread`` / ``$cell_surface_offset`` inside it.
        """

        return cls(
            path=data["path"],
            package=data["package"],
            nper=data["nper"],
            values=dict(data.get("values", {})),
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
    """Reference to a raster (GeoTIFF/IMG) sampled onto the grid for a surface or array.

    Use this for elevation surfaces (top/botm), recharge grids, or any gridded
    package input that lives in a raster. The raster is sampled at cell
    locations using ``method`` and, optionally, assigned to a named target via
    ``map_to``.

    Parameters
    ----------
    path
        Path to the raster file (``.tif`` / ``.img`` ...).
    band
        1-based band index to read (default: the first band).
    map_to
        Optional name of the package input/surface this raster populates.
    method
        Sampling/resampling method (e.g. ``"nearest"``, ``"linear"``). ``None``
        uses the caller's default.
    crs
        Optional CRS override; defaults to the raster's stored CRS.

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
    value
        The literal payload (must be JSON-serializable to survive ``to_dict``).
    metadata
        Arbitrary JSON-serializable annotations.

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


_SOURCE_TYPES: dict[str, type[Any]] = {
    DataSourceSpec.kind: DataSourceSpec,
    TableSource.kind: TableSource,
    ShapeSource.kind: ShapeSource,
    GeoPackageSourceSpec.kind: GeoPackageSourceSpec,
    BoundarySource.kind: BoundarySource,
    RasterSource.kind: RasterSource,
    LiteralSource.kind: LiteralSource,
}


def is_source_payload(data: Any) -> bool:
    """Whether ``data`` is a serialized source dict :func:`source_from_dict` rebuilds.

    The registry is the authority. A ``kind.endswith("Source")`` test looks
    equivalent and is not: :class:`DataSourceSpec`'s own kind is
    ``"DataSourceSpec"``, which ends in "Spec", so such a test silently passed
    the public base class through as a plain ``dict``.
    """

    return isinstance(data, dict) and str(data.get("kind")) in _SOURCE_TYPES


def source_from_dict(data: dict[str, Any]) -> DataSourceSpec | LiteralSource:
    """Recreate a source specification from its serialized representation."""

    kind = str(data.get("kind"))
    try:
        source_type = _SOURCE_TYPES[kind]
    except KeyError as error:
        raise ValueError(f"Unknown source kind: {kind!r}") from error
    return source_type._from_dict(data)


__all__ = [
    "BoundarySource",
    "DataSourceSpec",
    "is_source_payload",
    "GeoPackageSourceSpec",
    "LiteralSource",
    "RasterSource",
    "ShapeSource",
    "TableSource",
    "source_from_dict",
]

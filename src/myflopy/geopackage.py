"""GeoPackage-first construction of reusable MODFLOW model inputs."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Protocol

import geopandas as gpd
import numpy as np

from myflopy._logging import get_logger
from myflopy.advanced import (
    chd_spec,
    drn_spec,
    evt_spec,
    ghb_spec,
    rch_spec,
    riv_spec,
    wel_spec,
)
from myflopy.modflow.mf6.package_registry import _PACKAGE_EXPLORER_SPECS
from myflopy.specs import ModelContext, PackageSpec

logger = get_logger(__name__)

SurfaceReference = str


def _extent(geometry) -> float:
    """A geometry's own measure: area for a polygon, length for a line, 0 for a point."""

    if geometry is None or geometry.is_empty:
        return 0.0
    area = float(getattr(geometry, "area", 0.0) or 0.0)
    return area if area > 0 else float(getattr(geometry, "length", 0.0) or 0.0)


def _period_value(row, field: str | list[str] | tuple[str, ...], period: int):
    """Read a constant or per-period value from one GeoPackage feature."""

    if isinstance(field, str):
        return row[field]
    if not field:
        raise ValueError("At least one value field is required.")
    return row[field[min(period, len(field) - 1)]]


@dataclass(frozen=True, slots=True)
class CellSurfaceOffset:
    """An elevation or head placed relative to a model surface, per cell.

    Use it wherever a ``.gpkg`` helper takes an elevation-like value -- a drain
    invert, a river stage, a GHB head -- when what you mean is "so far below the
    ground" or "just above the cell floor" rather than a fixed number. Each cell the
    feature covers resolves its own value::

        value = surface(reference) + offset          # then floored at `minimum`

    Parameters
    ----------
    reference : {"model_top", "cell_top", "cell_bottom"}
        Which surface to measure from. **There is no default -- name it.**

        ``"model_top"``
            The ground surface. The same elevation for every layer, so a record's
            ``layer_field`` does not change it.
        ``"cell_top"``
            Top of the cell in *that record's* layer. Layer 1's cell top is the
            model top; layer 2's is layer 1's bottom, and so on.
        ``"cell_bottom"``
            Bottom of the cell in that record's layer.
    offset : float or str, default 0.0
        **Added** to the surface, so below it is NEGATIVE and above it is positive.
        A drain 2 ft below ground is ``offset=-2.0``; a drain 2 ft above the cell
        floor is ``CellSurfaceOffset("cell_bottom", offset=2.0)``. A ``str`` names a
        GeoPackage column, read per feature, so each row carries its own depth.
    minimum : float or str, optional
        A floor applied afterwards (``max(value, minimum)``), as a constant or a
        column name. It is a plain number, not a surface -- it cannot express "keep
        this inside its cell", because it never sees the cell bottom.

    Examples
    --------
    Every cell a drain line crosses, 2 ft above that cell's own floor::

        mf.drn.gpkg(path, layer="drains", context=ctx, nper=nper,
                    elevation=mf.CellSurfaceOffset("cell_bottom", offset=2.0),
                    conductance=mf.Spread("conductance"))

    A per-feature depth below ground, carried in a ``depth`` column::

        elevation=mf.CellSurfaceOffset("model_top", offset="depth")

    Notes
    -----
    ``"cell_top"`` and ``"cell_bottom"`` below layer 1 read ``grid.gdf_topbtm`` by
    integer column, which is the layout ``LayerBuildResult.attach_to_grid()``
    publishes. ``LayerStack.to_disv()`` publishes *named* columns instead, and with
    only those in place these two raise a ``KeyError``; call ``attach_to_grid()``
    last. ``"model_top"`` works with either, since it also matches ``"top"``.
    """

    reference: SurfaceReference
    offset: int | float | str = 0.0
    minimum: int | float | str | None = None

    def __post_init__(self) -> None:
        """Lowercase and validate ``reference`` (one of model_top/cell_top/cell_bottom)."""

        reference = self.reference.lower()
        valid = {"model_top", "cell_top", "cell_bottom"}
        if reference not in valid:
            raise ValueError(f"reference must be one of: {', '.join(sorted(valid))}.")
        object.__setattr__(self, "reference", reference)

    @staticmethod
    def _row_or_value(row, value: int | float | str | None) -> float | None:
        """Resolve a numeric constant or a GeoPackage field name to a float (``None`` -> ``None``)."""

        if value is None:
            return None
        if isinstance(value, str):
            return float(row[value])
        return float(value)

    def required_fields(self) -> set[str]:
        """Return GeoPackage fields required to resolve this value."""

        fields = set()
        if isinstance(self.offset, str):
            fields.add(self.offset)
        if isinstance(self.minimum, str):
            fields.add(self.minimum)
        return fields

    def resolve(self, source: SupportsSurfaceValue, row, *, layer: int, cell: int) -> float:
        """Return the resolved elevation/head for one mapped cell.

        ``source`` is anything exposing ``surface_value`` -- a :class:`GeoPackageSource`
        (the file path) or a :class:`SurfaceResolver` (the file-less areal-builder path).
        """

        surface = source.surface_value(self.reference, layer=layer, cell=cell)
        value = surface + self._row_or_value(row, self.offset)
        minimum = self._row_or_value(row, self.minimum)
        if minimum is not None:
            value = max(value, minimum)
        return float(value)


class SupportsSurfaceValue(Protocol):
    """Anything that can resolve a model/cell surface elevation for a cell."""

    def surface_value(self, reference: SurfaceReference, *, layer: int, cell: int) -> float:
        ...


@dataclass(frozen=True, slots=True)
class SurfaceResolver:
    """Resolve model/cell surface elevations for :class:`CellSurfaceOffset`.

    Wraps a :class:`~myflopy.specs.ModelContext` and owns the single copy of the
    ``top``/``bottom`` column lookup. Both :class:`GeoPackageSource` (the file
    path) and the file-less areal builders (``mf.evt(context=...)``) resolve a
    ``cell_top`` / ``model_top`` / ``cell_bottom`` elevation through this one
    place, so the ``gdf_topbtm`` column layout is never re-encoded at a call site.
    """

    context: ModelContext

    @property
    def grid(self):
        """Grid helper carried by the model context."""

        return self.context.grid

    def _surfaces(self):
        """The surface source for :class:`CellSurfaceOffset`: the context's, else the grid's top/botm."""

        surfaces = self.context.surfaces
        if surfaces is None:
            surfaces = getattr(self.grid, "gdf_topbtm", None)
        if surfaces is None:
            raise ValueError(
                "ModelContext.surfaces or grid.gdf_topbtm is required for CellSurfaceOffset."
            )
        return surfaces

    @staticmethod
    def _cell_sequence(value: Any, cell: int) -> float:
        """One cell's value from a scalar (broadcast) or a per-cell array."""

        array = np.asarray(value, dtype=float)
        if array.ndim == 0:
            return float(array)
        return float(array.reshape(-1)[cell])

    def surface_value(self, reference: SurfaceReference, *, layer: int, cell: int) -> float:
        """Return a model/cell surface value for a mapped boundary cell."""

        surfaces = self._surfaces()
        if isinstance(surfaces, dict):
            if reference in {"model_top", "cell_top"} and layer == 0:
                return self._cell_sequence(surfaces["top"], cell)
            if reference == "model_top":
                return self._cell_sequence(surfaces["top"], cell)
            bottom = surfaces.get("bottom", surfaces.get("botm"))
            if bottom is None:
                raise ValueError("Surface dictionary requires 'bottom' or 'botm'.")
            bottom_array = np.asarray(bottom, dtype=float)
            if reference == "cell_top":
                if bottom_array.ndim == 1:
                    raise ValueError("cell_top for layers below 0 requires multilayer bottom surfaces.")
                return self._cell_sequence(bottom_array[layer - 1], cell)
            if bottom_array.ndim == 1:
                if layer != 0:
                    raise ValueError("cell_bottom for layers below 0 requires multilayer bottom surfaces.")
                return self._cell_sequence(bottom_array, cell)
            return self._cell_sequence(bottom_array[layer], cell)

        columns = getattr(surfaces, "columns", ())
        if reference == "model_top":
            candidates = (0, "top", "model_top")
        elif reference == "cell_top":
            candidates = ((0, "top", "model_top") if layer == 0 else (layer, f"layer_{layer}_top"))
        else:
            candidates = (layer + 1, "bottom" if layer == 0 else f"layer_{layer}_bottom")
        for column in candidates:
            if column in columns:
                return float(surfaces.loc[cell, column])
        raise ValueError(
            f"Could not resolve {reference!r} for layer {layer}; available surface columns are {list(columns)!r}."
        )


@dataclass(frozen=True, slots=True)
class Spread:
    """Split an EXTENSIVE field across the cells a feature covers.

    A conductance or a flux belongs to the *feature*, not to each cell the feature
    happens to touch. The default resolution writes a row's value to every
    intersected cell, which is correct for an elevation or a head and multiplies an
    extensive value by the cell count: measured on a real DRN drawn as one line per
    source cell, ``mf.drn.gpkg`` on the grid the values came from turned 205 records
    into 424 and 66,007.43 ft2/d into 133,631.25 -- a factor of 2.02, on an
    unchanged mesh. Wrap the field to split it instead::

        mf.drn.gpkg(path, layer="drn_lines", context=ctx, nper=1,
                    conductance=mf.Spread("conductance"))

    Each cell receives ``value * share``, where ``share`` is the cell's fraction of
    the feature -- length for a line, area for a polygon. A point has no extent to
    split, so it keeps the whole value and one cell, which is already correct.

    ``mode`` decides what happens to the part of a feature that falls outside the
    active grid:

    ``"clip"`` (default)
        Shares are fractions of the WHOLE feature, so the part outside is simply not
        applied and the total drops. Honest: that length is not in the model.
    ``"retained"``
        Shares are renormalized over the cells actually used, preserving the total by
        concentrating it on the part that remains. A different model -- choose it
        deliberately.

    Either way the resolved fraction is logged at DEBUG, and a feature losing more
    than ``warn_below`` of itself is logged at WARNING.

    ``min_share`` drops cells a feature barely grazes. A line drawn between cell
    centres clips the corners of its neighbours, which adds cells the feature does
    not meaningfully occupy: measured on a real DRN, 11 of 79 cells were corner
    clips holding 0.90% of the conductance between them, none above 0.37%
    individually. Dropped shares are removed before ``mode`` is applied, so
    ``"retained"`` still preserves the total exactly.
    """

    field: str
    mode: str = "clip"
    min_share: float = 0.01
    warn_below: float = 0.99

    def __post_init__(self) -> None:
        """Validate ``mode``."""

        if self.mode not in {"clip", "retained"}:
            raise ValueError(f"Spread mode must be 'clip' or 'retained', not {self.mode!r}.")
        if not 0.0 <= self.min_share < 1.0:
            raise ValueError(f"Spread min_share must be in [0, 1), not {self.min_share!r}.")


RowValue = str | int | float | list[str] | tuple[str, ...] | CellSurfaceOffset | Spread


def _required_value_fields(value: RowValue) -> set[str]:
    """Return GeoPackage fields needed by one value specification."""

    if isinstance(value, CellSurfaceOffset):
        return value.required_fields()
    if isinstance(value, str):
        return {value}
    return set(value) if isinstance(value, (list, tuple)) else set()


_SPEC_FACTORIES = {
    "chd": chd_spec,
    "drn": drn_spec,
    "evt": evt_spec,
    "ghb": ghb_spec,
    "rch": rch_spec,
    "riv": riv_spec,
    "wel": wel_spec,
}


def _metadata_value(value: RowValue):
    """Return a manifest-friendly representation of one value specification.

    Every dataclass in :data:`RowValue` needs a branch here. A run manifest is
    written with :func:`json.dumps` and no ``default=``, so a spec that reaches it
    unconverted fails at ``prepare_run`` -- long after the package built cleanly.
    ``tests/test_geopackage_spread.py`` asserts the whole union round-trips.
    """

    if isinstance(value, CellSurfaceOffset):
        return {
            "type": "CellSurfaceOffset",
            "reference": value.reference,
            "offset": value.offset,
            "minimum": value.minimum,
        }
    if isinstance(value, Spread):
        return {
            "type": "Spread",
            "field": value.field,
            "mode": value.mode,
            "min_share": value.min_share,
        }
    return value


@dataclass(slots=True)
class GeoPackageSource:
    """Map one GeoPackage layer onto a model grid and build reusable inputs.

    Layer values are assumed to be one-based by default because that is the
    usual convention in external GIS inputs. Stress-period values are assumed
    to be zero-based.
    """

    path: Path | str
    context: ModelContext
    nper: int
    layer: str | None = None
    name_field: str | None = "name"
    layer_field: str | None = "layer"
    period_field: str | None = None
    layer_base: int = 1
    period_base: int = 0
    _gdf: gpd.GeoDataFrame | None = field(default=None, init=False, repr=False)

    def __post_init__(self) -> None:
        """Coerce ``path`` and validate the file exists, a grid is present, and ``nper >= 1``."""

        self.path = Path(self.path)
        if not self.path.exists():
            raise FileNotFoundError(f"GeoPackage not found: {self.path}")
        if self.context.grid is None:
            raise ValueError("ModelContext.grid is required for GeoPackage mapping.")
        if self.nper < 1:
            raise ValueError("nper must be at least 1.")

    @property
    def grid(self):
        """Grid helper carried by the model context."""

        return self.context.grid

    @property
    def gdf(self) -> gpd.GeoDataFrame:
        """Read and align the GeoPackage layer with the model grid CRS."""

        if self._gdf is None:
            kwargs = {} if self.layer is None else {"layer": self.layer}
            gdf = gpd.read_file(self.path, **kwargs)
            grid_crs = getattr(self.grid, "crs", None)
            if grid_crs is not None and gdf.crs is not None and gdf.crs != grid_crs:
                gdf = gdf.to_crs(grid_crs)
            self._gdf = gdf
        return self._gdf

    def _active(self, layer: int, cell: int) -> bool:
        """Whether ``(layer, cell)`` can carry a boundary, per the context's idomain.

        MF6 idomain is three-valued, not a flag: ``> 0`` active, ``0`` inactive, and
        ``< 0`` **vertical passthrough** -- the cell is removed and flow passes
        straight through it. A passthrough cell holds no boundary; MF6 refuses one
        with ``Cell is outside active grid domain``.

        So the test is ``> 0``, not truthiness. ``bool(-1)`` is ``True``, which put
        347 DRN records into passthrough cells on a real model and failed the run at
        read time -- a `LayerStack` using ``pinch="passthrough"`` produces these in
        quantity (13,797 on the model in question).
        """

        domain = self.context.domain
        if domain is None:
            return True
        values = np.asarray(domain)
        if values.ndim == 1:
            return bool(values[cell] > 0)
        if layer >= values.shape[0]:
            raise ValueError(
                f"GeoPackage layer {layer} is outside domain with {values.shape[0]} layers."
            )
        return bool(values[layer, cell] > 0)

    def _layer(self, row) -> int:
        """The zero-based grid layer for a feature (from ``layer_field`` minus ``layer_base``)."""

        if self.layer_field is None:
            return 0
        layer = int(row[self.layer_field]) - self.layer_base
        if layer < 0:
            raise ValueError(f"GeoPackage layer resolves to negative index: {layer}")
        return layer

    def _periods(self, row) -> range | tuple[int]:
        """The stress periods a feature applies to: all periods, or its single ``period_field``."""

        if self.period_field is None:
            return range(self.nper)
        period = int(row[self.period_field]) - self.period_base
        if period not in range(self.nper):
            raise ValueError(f"GeoPackage period is outside the simulation: {period}")
        return (period,)

    def _cells(self, geometry, *, layer: int, edges_only: bool = False) -> list[int]:
        """The active grid cells a feature's ``geometry`` intersects (optionally edge cells only)."""

        return [cell for cell, _ in self._cell_shares(geometry, layer=layer, edges_only=edges_only)]

    def _cell_shares(
        self, geometry, *, layer: int, edges_only: bool = False
    ) -> list[tuple[int, float]]:
        """Active cells a feature covers, each with its fraction of the whole feature.

        The fraction is by length for a line, by area for a polygon, and 1.0 for a
        point or anything with no extent -- there is nothing to divide. Fractions are
        of the WHOLE feature, so they sum to less than one where part of it lies
        outside the active grid; :class:`Spread` decides what that means.

        Candidates come from the spatial index rather than an elementwise
        ``intersects`` over every cell, which was O(features x ncpl).
        """

        polys = self.grid.gdf_vorPolys
        candidates = polys.index[polys.sindex.query(geometry, predicate="intersects")]
        if edges_only:
            edge_cells = set(self.grid.get_grid_edge())
            candidates = [cell for cell in candidates if cell in edge_cells]

        measure = _extent(geometry)
        shares: list[tuple[int, float]] = []
        for cell in candidates:
            cell = int(cell)
            if not self._active(layer, cell):
                continue
            if measure <= 0:
                shares.append((cell, 1.0))
                continue
            piece = _extent(geometry.intersection(polys.geometry.loc[cell]))
            if piece <= 0:
                continue
            shares.append((cell, piece / measure))
        return shares

    def surface_value(self, reference: SurfaceReference, *, layer: int, cell: int) -> float:
        """Return a model/cell surface value for a mapped boundary cell.

        Delegates to the shared :class:`SurfaceResolver` so the ``gdf_topbtm``
        column layout lives in exactly one place (also used by ``mf.evt(context=...)``).
        """

        return SurfaceResolver(self.context).surface_value(reference, layer=layer, cell=cell)

    def _value(self, row, value: RowValue, *, period: int, layer: int, cell: int, share: float = 1.0):
        """Resolve one value spec for a cell.

        ``share`` is the cell's fraction of the feature and is used only by
        :class:`Spread`; everything else resolves the same value for every cell,
        which is what an elevation or a head should do.
        """

        if isinstance(value, CellSurfaceOffset):
            return value.resolve(self, row, layer=layer, cell=cell)
        if isinstance(value, Spread):
            return float(_period_value(row, value.field, period)) * share
        if isinstance(value, (int, float)):
            return value
        return _period_value(row, value, period)

    def _boundary_data(
        self,
        *fields: RowValue,
        edges_only: bool = False,
        boundnames: bool = False,
    ) -> dict[int, list[list[Any]]]:
        """Build MF6 stress-period data by mapping every feature to its cells and value fields.

        Validates the required columns exist, then for each feature emits one
        ``[(layer, cell), *values(, boundname)]`` record per applicable period and
        intersected cell, keyed by stress period.
        """

        spreads = tuple(f for f in fields if isinstance(f, Spread))
        required = {
            field_name
            for field in fields
            for field_name in _required_value_fields(field)
        }
        if self.layer_field is not None:
            required.add(self.layer_field)
        if self.period_field is not None:
            required.add(self.period_field)
        if boundnames and self.name_field is not None:
            required.add(self.name_field)
        missing = sorted(required - set(self.gdf.columns))
        if missing:
            raise ValueError(
                f"GeoPackage '{self.path.name}' is missing fields: {', '.join(missing)}"
            )

        data = {period: [] for period in range(self.nper)}
        for index, row in self.gdf.iterrows():
            layer = self._layer(row)
            shares = self._cell_shares(row.geometry, layer=layer, edges_only=edges_only)
            if spreads:
                shares = self._resolved_shares(shares, spreads, index)
            name = (
                str(row[self.name_field])
                if boundnames and self.name_field is not None
                else str(index)
            )
            for period in self._periods(row):
                for cell, share in shares:
                    values = [
                        self._value(row, field, period=period, layer=layer, cell=cell, share=share)
                        for field in fields
                    ]
                    record = [(layer, cell), *values]
                    if boundnames:
                        record.append(name)
                    data[period].append(record)
        return data

    def _resolved_shares(self, shares, spreads: tuple[Spread, ...], index: Any):
        """Apply the :class:`Spread` mode to one feature's cell shares, and say what was lost."""

        floor = max(spread.min_share for spread in spreads)
        if floor > 0:
            kept = [(cell, share) for cell, share in shares if share >= floor]
            if len(kept) < len(shares):
                logger.debug(
                    "feature %s: dropped %d cell(s) below min_share=%.3g",
                    index, len(shares) - len(kept), floor,
                )
            shares = kept or shares
        covered = sum(share for _, share in shares)
        mode = spreads[0].mode
        warn_below = min(spread.warn_below for spread in spreads)
        if covered <= 0:
            logger.warning(
                "feature %s covers no active cell; its spread fields contribute nothing", index
            )
            return shares
        if covered < warn_below:
            logger.warning(
                "feature %s lies %.1f%% inside the active grid; with mode=%r its spread "
                "fields %s",
                index,
                covered * 100.0,
                mode,
                "lose the remainder" if mode == "clip" else "keep their total anyway",
            )
        else:
            logger.debug("feature %s covered %.4f of itself over %d cell(s)", index, covered, len(shares))
        if mode == "retained":
            return [(cell, share / covered) for cell, share in shares]
        return shares

    def _metadata(self, method: str, **fields: Any) -> dict[str, Any]:
        """Return serializable provenance for a generated model input."""

        return {
            "source_type": "geopackage",
            "source_path": str(self.path),
            "source_layer": self.layer,
            "builder": method,
            "fields": {name: _metadata_value(value) for name, value in fields.items()},
            "name_field": self.name_field,
            "layer_field": self.layer_field,
            "period_field": self.period_field,
            "layer_base": self.layer_base,
            "period_base": self.period_base,
        }

    def _bc_from_features(
        self,
        package: str,
        *,
        name: str,
        boundnames: bool,
        edges_only: bool = False,
        options: dict[str, Any] | None = None,
        **fields: RowValue,
    ) -> PackageSpec:
        """Map features to cells and build one list-BC spec, driven by the registry.

        The seven public resolvers are thin wrappers over this. They keep their
        own signatures -- each names its package's real columns with real
        defaults, which is what makes ``mf.riv.gpkg(stage=..., rbot=...)``
        discoverable -- while the mapping, the spec call and the provenance
        metadata happen here once.

        ``fields`` must be passed in MF6 record order; the assertion below
        checks that against the descriptor rather than trusting the call site,
        because silently reordering a record is the one mistake in this file
        that would produce a model that runs and is wrong.
        """

        descriptor = _PACKAGE_EXPLORER_SPECS[package]
        expected = tuple(descriptor.gpkg_defaults)
        assert tuple(fields) == expected, (
            f"{package}: fields {tuple(fields)} are not in the descriptor's "
            f"record order {expected}"
        )

        data_options: dict[str, Any] = {"boundnames": boundnames}
        if descriptor.capabilities.edges_only:
            data_options["edges_only"] = edges_only
        elif edges_only:
            raise ValueError(f"{package} features cannot be restricted to grid edges.")

        return _SPEC_FACTORIES[package](
            self._boundary_data(*fields.values(), **data_options),
            name=name,
            boundnames=boundnames,
            **(options or {}),
        ).with_metadata(**self._metadata(package, **fields))

    def chd(
        self,
        *,
        head: RowValue = "head",
        name: str = "chd",
        edges_only: bool = False,
        boundnames: bool = True,
        **options,
    ) -> PackageSpec:
        """Return a CHD package spec from GeoPackage features."""

        return self._bc_from_features(
            "chd", name=name, boundnames=boundnames, edges_only=edges_only,
            options=options, head=head,
        )

    def ghb(
        self,
        *,
        head: RowValue = "head",
        conductance: RowValue = "conductance",
        name: str = "ghb",
        edges_only: bool = False,
        boundnames: bool = True,
        **options,
    ) -> PackageSpec:
        """Return a GHB package spec from GeoPackage features."""

        return self._bc_from_features(
            "ghb", name=name, boundnames=boundnames, edges_only=edges_only,
            options=options, head=head, conductance=conductance,
        )

    def drn(
        self,
        *,
        elevation: RowValue = "elevation",
        conductance: RowValue = "conductance",
        name: str = "drn",
        edges_only: bool = False,
        boundnames: bool = True,
        **options,
    ) -> PackageSpec:
        """Return a DRN package spec from GeoPackage features."""

        return self._bc_from_features(
            "drn", name=name, boundnames=boundnames, edges_only=edges_only,
            options=options, elevation=elevation, conductance=conductance,
        )

    def riv(
        self,
        *,
        stage: RowValue = "stage",
        conductance: RowValue = "conductance",
        rbot: RowValue = "rbot",
        name: str = "riv",
        edges_only: bool = False,
        boundnames: bool = True,
        **options,
    ) -> PackageSpec:
        """Return a RIV package spec from GeoPackage features."""

        return self._bc_from_features(
            "riv", name=name, boundnames=boundnames, edges_only=edges_only,
            options=options, stage=stage, conductance=conductance, rbot=rbot,
        )

    def wel(
        self,
        *,
        rate: RowValue = "rate",
        name: str = "wel",
        boundnames: bool = True,
        **options,
    ) -> PackageSpec:
        """Return a WEL package spec from GeoPackage features."""

        return self._bc_from_features(
            "wel", name=name, boundnames=boundnames, options=options, rate=rate,
        )

    def rch(
        self,
        *,
        recharge: RowValue = "recharge",
        name: str = "rch",
        boundnames: bool = True,
        **options,
    ) -> PackageSpec:
        """Return a list-based RCH package spec from GeoPackage features."""

        return self._bc_from_features(
            "rch", name=name, boundnames=boundnames, options=options,
            recharge=recharge,
        )

    def evt(
        self,
        *,
        surface: RowValue = "surface",
        rate: RowValue = "rate",
        depth: RowValue = "depth",
        name: str = "evt",
        boundnames: bool = True,
        **options,
    ) -> PackageSpec:
        """Return a list-based EVT package spec from GeoPackage features.

        Single-segment ET only: feature mapping emits fixed
        ``(cellid, surface, rate, depth)`` records, which is exactly the
        ``nseg=1`` record shape. Segmented ET needs ``nseg - 1`` extra
        ``pxdp``/``petm`` values per record that no feature mapping supplies,
        so it is rejected here rather than failing later inside FloPy.
        """

        if int(options.get("nseg", 1)) != 1:
            raise ValueError(
                "GeoPackage-driven EVT supports nseg=1 only: mapped features "
                "yield (cellid, surface, rate, depth) records with no pxdp/petm "
                "values. For segmented ET, assemble the records yourself and "
                "use mf.evt(...) / mf.evt.flopy(...) with nseg=."
            )

        return self._bc_from_features(
            "evt", name=name, boundnames=boundnames, options=options,
            surface=surface, rate=rate, depth=depth,
        )

    def k_array(
        self,
        *,
        value: str = "k",
        nlay: int,
        defaults: float | list[float],
    ) -> np.ndarray:
        """Return an ``(nlay, ncpl)`` K array from GeoPackage polygons."""

        fallback = [defaults] * nlay if np.isscalar(defaults) else list(defaults)
        if len(fallback) != nlay:
            raise ValueError("defaults must be a scalar or contain one value per layer.")
        result = np.array(
            [[float(fallback[layer])] * self.grid.ncpl for layer in range(nlay)],
            dtype=float,
        )
        for _, row in self.gdf.iterrows():
            layer = self._layer(row)
            if layer >= nlay:
                raise ValueError(f"GeoPackage layer {layer} is outside nlay={nlay}.")
            for cell in self._cells(row.geometry, layer=layer):
                result[layer, cell] = float(row[value])
        return result


__all__ = ["CellSurfaceOffset", "GeoPackageSource"]

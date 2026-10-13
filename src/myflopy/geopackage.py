"""GeoPackage-first construction of reusable MODFLOW model inputs."""

from __future__ import annotations

import warnings
from dataclasses import MISSING, dataclass, field
from dataclasses import fields as dataclass_fields
from pathlib import Path
from typing import Any, Protocol

import geopandas as gpd
import numpy as np

from myflopy._logging import get_logger
from myflopy._warnings import stacklevel_outside_package
from myflopy.advanced import (
    chd_spec,
    drn_spec,
    evt_spec,
    ghb_spec,
    rch_spec,
    riv_spec,
    wel_spec,
)
from myflopy.areas import _same_crs
from myflopy.modflow.mf6.boundary_support import bottom_clearance, river_pair_lift
from myflopy.modflow.mf6.grid.geometry import (
    _extent,  # noqa: F401 -- re-exported: callers and tests import it from here
    cell_shares,
    domain_active,
    miss_reason,
    polygon_cell_values,
    where_they_are,
)
from myflopy.modflow.mf6.grid.surfaces import (
    _fill_nearest,
    check_gaps,
    get_raster_vals_at_centroids,
    raster_extents_text,
)
from myflopy.modflow.mf6.package_registry import _PACKAGE_EXPLORER_SPECS
from myflopy.sources import ArrayProduct, BoundarySource, Polygons, RasterValues
from myflopy.specs import (
    GRID_PIN,
    GRID_PIN_OPTIONS,
    MOVE_RECIPE,
    SOURCE_SHA,
    ModelContext,
    PackageSpec,
    _json_value,
    _layer_elevations,
    grid_pin,
    source_fingerprint,
)

logger = get_logger(__name__)

SurfaceReference = str

#: The three surfaces an elevation can be measured from, in the order they are taught.
SURFACE_REFERENCES = ("model_top", "cell_top", "cell_bottom")


def surface_reference(reference: Any) -> str:
    """``reference`` lowercased, or a ``ValueError`` naming the three accepted spellings.

    The ONE rule for a surface reference: :class:`CellSurfaceOffset` and both
    :class:`SurfaceResolver` verbs (so :meth:`GeoPackageSource.surface_value` too) call
    it. The resolvers used to take anything -- the table decoder's last branch was a
    bare ``else`` -- so ``"Model_Top"`` or a typo read the CELL BOTTOM: measured 0.0
    where ``"model_top"`` read 10.0 on the same one-cell table, with nothing said
    (ledger 392).
    """

    spelled = reference.lower() if isinstance(reference, str) else None
    if spelled not in SURFACE_REFERENCES:
        raise ValueError(
            f"A surface reference is one of {', '.join(map(repr, SURFACE_REFERENCES))} "
            f"(any case), not {reference!r}: 'model_top' is the ground surface, "
            "'cell_top'/'cell_bottom' the top/bottom of the record's own cell."
        )
    return spelled


def _spec_payload(value: Any) -> dict[str, Any]:
    """Return the JSON-ready declaration of a flat ``RowValue`` dataclass.

    DEFAULTS ARE OMITTED, so adding a defaulted field later cannot move an existing
    package's digest -- since ledger 273 ``_canonical`` fingerprints this payload
    rather than the object.

    Every value goes through ``_json_value`` FIRST, which is where this deviates from
    ``Surface.to_dict``. ``type(np.float64(0.01)) is not float``, so an un-normalized
    numpy scalar fails the default test, gets written, and comes back a plain float
    that IS omitted -- two payloads for one object, and a digest that moves across a
    save. Normalizing here also makes a bare ``to_dict()`` safe for
    ``yaml.safe_dump``, which dispatches on the EXACT type and refuses a numpy scalar
    ``json.dumps`` takes happily (ledger 265).
    """

    payload: dict[str, Any] = {}
    for item in dataclass_fields(value):
        resolved = _json_value(getattr(value, item.name))
        default = (
            item.default_factory()
            if item.default is MISSING and item.default_factory is not MISSING
            else item.default
        )
        if type(resolved) is type(default) and resolved == default:
            continue
        payload[item.name] = resolved
    return payload


def _from_spec_payload(cls, data: Any):
    """Rebuild a flat ``RowValue`` dataclass from :func:`_spec_payload`'s output.

    Goes through ``__init__``, so ``__post_init__`` re-runs: a hand-edited payload is
    refused where the reader can act on it, rather than reviving an object that dies
    later inside a resolver. Serializing these is precisely what puts them in front of
    a text editor, so the validation is the feature.

    ``TypeError`` is translated too -- a missing required field otherwise reads
    ``Spread.__init__() missing 1 required positional argument: 'field'``, which names
    a constructor the reader never called.
    """

    if not isinstance(data, dict):
        raise ValueError(
            f"A {cls.__name__} payload must be a mapping of field names, "
            f"not a {type(data).__name__}."
        )
    try:
        return cls(**data)
    except TypeError as error:
        raise ValueError(f"{cls.__name__} payload is not usable: {error}") from error


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
        Which surface to measure from. **There is no default -- name it.** Any case
        is accepted and stored lowercased; another word raises ``ValueError`` naming
        the three.

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

    #: The key `myflopy.specs._json_value` writes this object under and `_spec_value`
    #: revives it from. A `$`-prefixed token, NOT a `kind` tag: the source dispatcher
    #: keys on `kind`, so one tagged as the other comes back a plain dict, silently
    #: (ledger 269). `tests/test_spec_token_registry.py` reads this by AST, so it must
    #: stay a string literal in the class body.
    _SPEC_TOKEN = "$cell_surface_offset"

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-ready declaration of this offset."""

        return _spec_payload(self)

    @classmethod
    def from_dict(cls, data: Any) -> CellSurfaceOffset:
        """Rebuild an offset from :meth:`to_dict`'s payload.

        Parameters
        ----------
        data : dict
            Field name to value: ``reference`` (required), and ``offset`` /
            ``minimum``, each of which takes the class default when omitted --
            :meth:`to_dict` leaves defaults out. It goes through the constructor, so
            ``reference`` is lowercased and validated exactly as it is there.

        Returns
        -------
        CellSurfaceOffset

        Raises
        ------
        ValueError
            ``data`` is not a ``dict``; it names a key that is not a field or lacks
            ``reference`` (the constructor's ``TypeError``, reworded); or
            ``reference`` is not ``model_top``, ``cell_top`` or ``cell_bottom``.
        """

        return _from_spec_payload(cls, data)

    def __post_init__(self) -> None:
        """Lowercase and validate ``reference`` (one of model_top/cell_top/cell_bottom)."""

        object.__setattr__(self, "reference", surface_reference(self.reference))

    @staticmethod
    def _row_or_value(
        row, value: int | float | str | None, *, role: str = "offset", source: str = ""
    ) -> float | None:
        """A numeric constant or a GeoPackage column, as a float (``None`` -> ``None``).

        **A NULL in the column is REFUSED, not passed on.** A nullable REAL column reads
        back as ``nan``, and ``max(value, nan)`` returns ``value`` in Python -- so a NULL
        in the ``minimum`` column SILENTLY DISABLED the floor, and the result was finite
        so ledger 247's own check let it through. Measured: a drain declared with
        ``minimum="floor"`` over a NULL column wrote output BYTE-IDENTICAL to no floor at
        all, 11 of 11 records below their own cell bottom, zero warnings (ledger 289).
        The other direction was no better -- a TEXT NULL raised a bare ``TypeError``
        naming neither the feature nor the field, where the same NULL in a plain value
        column gives ledger 247's full named error.
        """

        if value is None:
            return None
        if isinstance(value, str):
            return _finite(
                row[value],
                field=value,
                label=f"{row.get('name', getattr(row, 'name', '?'))!r} ({role}=)",
                source=source or "",
            )
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

        Parameters
        ----------
        source : GeoPackageSource or SurfaceResolver
            Any object with a ``surface_value(reference, *, layer, cell)`` method; it
            answers the surface this offset is measured from. When it has a ``path``
            (a :class:`GeoPackageSource`), that file's name appears in errors.
        row : pandas.Series or dict
            The feature's attribute row (the file-less
            :class:`~myflopy.modflow.mf6.evapotranspiration.EVTBuilder` passes an
            empty ``dict``). Read only where ``offset`` or ``minimum``
            names a column (``row[column]``); its ``"name"`` entry, else a Series's
            index label, identifies the feature in an error.
        layer : int
            Zero-based model layer of the record. It picks the cell for
            ``"cell_top"`` and ``"cell_bottom"``; ``"model_top"`` ignores it.
        cell : int
            Zero-based cell (CELL2D) number within the layer.

        Returns
        -------
        float
            ``surface + offset``, raised to ``minimum`` when that is larger.

        Raises
        ------
        ValueError
            A column named by ``offset`` or ``minimum`` holds a NULL or a non-number
            in ``row``; or ``source`` cannot resolve the surface -- it has no layer
            elevations, or none for ``reference`` at ``layer``.
        KeyError
            ``offset`` or ``minimum`` names a column ``row`` does not have.
        """

        surface = source.surface_value(self.reference, layer=layer, cell=cell)
        where = getattr(source, "path", None)
        where = f"'{where.name}'" if where is not None else ""
        value = surface + self._row_or_value(row, self.offset, role="offset", source=where)
        minimum = self._row_or_value(row, self.minimum, role="minimum", source=where)
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

        surfaces = _layer_elevations(self.context)
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
        """Return a model/cell surface value for a mapped boundary cell.

        Parameters
        ----------
        reference : {"model_top", "cell_top", "cell_bottom"}
            Which surface to read: the ground surface, or the top or bottom of the
            cell in ``layer``. Any case; another word raises ``ValueError`` naming
            the three, as :class:`CellSurfaceOffset` does -- one rule,
            :func:`surface_reference`.
        layer : int
            Zero-based model layer. ``"cell_top"`` of layer ``n > 0`` is the bottom of
            layer ``n - 1``; ``"model_top"`` ignores it.
        cell : int
            Zero-based cell (CELL2D) number.

        Returns
        -------
        float
            The elevation, read from ``context.surfaces`` or, when that is ``None``,
            ``grid.gdf_topbtm``.

        Raises
        ------
        ValueError
            ``reference`` is not one of the three; neither table is available; or the
            table has no column for ``reference`` at ``layer``.
        """

        reference = surface_reference(reference)
        return self._decode(self._surfaces(), reference, layer=layer, cell=cell)

    def surface_column(self, reference: SurfaceReference, *, layer: int, ncpl: int):
        """The whole per-column surface -- the ARRAY form's counterpart to :meth:`surface_value`.

        One value per column is what ``READASARRAYS`` wants, and resolving it a cell at
        a time through :meth:`surface_value` would re-read the elevation table ncpl
        times: ``_layer_elevations`` compares the context's table against the grid's
        live one on every call, so that is an O(ncpl x table) walk, not a loop of
        lookups. The table is fetched ONCE here and the same decoder runs per cell.

        Returns a plain ``float`` when every column agrees, so a flat land surface still
        writes MODFLOW's ``CONSTANT`` form -- and stays genuinely mesh-independent, which
        is what keeps it off the per-column pin (ledger 282).

        Parameters
        ----------
        reference : {"model_top", "cell_top", "cell_bottom"}
            As :meth:`surface_value`: any case, another word raises ``ValueError``.
        layer : int
            Zero-based model layer, as :meth:`surface_value`.
        ncpl : int
            Number of cells per layer: one value is resolved for each of ``0 .. ncpl - 1``.

        Returns
        -------
        float or numpy.ndarray
            A ``float`` when every column agrees, else a float array of length ``ncpl``.
        """

        reference = surface_reference(reference)
        surfaces = self._surfaces()
        values = np.array(
            [self._decode(surfaces, reference, layer=layer, cell=cell) for cell in range(int(ncpl))],
            dtype=float,
        )
        return float(values[0]) if values.size and np.all(values == values[0]) else values

    def _decode(self, surfaces, reference: SurfaceReference, *, layer: int, cell: int) -> float:
        """The ONE place the elevation table's layout is read. Both verbs go through it."""

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

    Parameters
    ----------
    field : str
        Attribute column holding each feature's TOTAL -- a conductance or a flux for
        the whole feature, in the package's own units. Each covered cell is written
        ``value * share``. A NULL or a non-number raises ``ValueError`` naming the
        feature; a column the file lacks raises ``ValueError`` naming the file and the
        column, before any record is built, as an unwrapped column does. Refused with
        ``TypeError`` on an elevation argument (``chd``/``ghb`` ``head``, ``drn``
        ``elevation``, ``riv`` ``stage``/``rbot``, ``evt`` ``surface``), where there is
        no total to divide.
    mode : {"clip", "retained"}, default "clip"
        What becomes of the part of a feature outside the active grid (above):
        ``"clip"`` keeps shares as fractions of the whole feature, so that part is
        not written; ``"retained"`` divides each kept share by the fraction covered,
        so the written cells still sum to the total. Anything else raises
        ``ValueError``. When one call wraps several fields they must all use the
        same mode -- a feature has one set of cell shares -- and two that differ
        raise ``ValueError`` naming both, at the call.
    min_share : float, default 0.01
        Cells holding less than this fraction of the feature are dropped before
        ``mode`` is applied -- dropped from the feature's records altogether, every
        field included, not only this one. When every cell falls below it none is
        dropped. Must be in ``[0, 1)`` (``ValueError`` otherwise); ``0`` keeps every
        cell. When one call wraps several fields, the largest applies.
    warn_below : float, default 0.99
        A feature whose kept cells hold less than this fraction of it is logged at
        WARNING through the ``myflopy`` logger, as a log line only; otherwise its
        coverage is logged at DEBUG. Separately, whatever this is, a feature lying
        only partly on usable cells is reported under the source's ``gaps`` (a
        ``UserWarning`` by default). Not validated. When one call wraps several
        fields, the smallest applies.
    """

    field: str
    mode: str = "clip"
    min_share: float = 0.01
    warn_below: float = 0.99

    #: The key `myflopy.specs._json_value` writes this object under and `_spec_value`
    #: revives it from. A `$`-prefixed token, NOT a `kind` tag: the source dispatcher
    #: keys on `kind`, so one tagged as the other comes back a plain dict, silently
    #: (ledger 269). `tests/test_spec_token_registry.py` reads this by AST, so it must
    #: stay a string literal in the class body.
    _SPEC_TOKEN = "$spread"

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-ready declaration of this split."""

        return _spec_payload(self)

    @classmethod
    def from_dict(cls, data: Any) -> Spread:
        """Rebuild a split from :meth:`to_dict`'s payload.

        Parameters
        ----------
        data : dict
            Field name to value: ``field`` (required), and ``mode`` / ``min_share``
            / ``warn_below``, each of which takes the class default when omitted --
            :meth:`to_dict` leaves defaults out. It goes through the constructor, so
            ``mode`` and ``min_share`` are validated exactly as they are there.

        Returns
        -------
        Spread

        Raises
        ------
        ValueError
            ``data`` is not a ``dict``; it names a key that is not a field or lacks
            ``field`` (the constructor's ``TypeError``, reworded); or ``mode`` or
            ``min_share`` is out of range.
        """

        return _from_spec_payload(cls, data)

    def __post_init__(self) -> None:
        """Validate ``mode``."""

        if self.mode not in {"clip", "retained"}:
            raise ValueError(f"Spread mode must be 'clip' or 'retained', not {self.mode!r}.")
        if not 0.0 <= self.min_share < 1.0:
            raise ValueError(f"Spread min_share must be in [0, 1), not {self.min_share!r}.")


RowValue = str | int | float | list[str] | tuple[str, ...] | CellSurfaceOffset | Spread


def _required_value_fields(value: RowValue) -> set[str]:
    """Return GeoPackage fields needed by one value specification.

    A :class:`Spread` names a column too: left out, a typo in one was missed by the
    up-front check and surfaced mid-build as a bare ``KeyError('typo')`` naming neither
    the file nor the feature.
    """

    if isinstance(value, CellSurfaceOffset):
        return value.required_fields()
    if isinstance(value, Spread):
        return _required_value_fields(value.field)
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


def _finite(value, *, field: RowValue, label: str, source: str):
    """Return ``value`` if it is a real number, else say which feature is unusable.

    A GeoPackage NULL reads back as ``None``, and neither path this guards used to
    look: ``k_array`` did ``float(row[value])``, which makes a NaN, and the boundary
    record kept the ``None``, which numpy casts to NaN when FloPy builds the
    recarray. Either way a non-finite number reached a written MF6 package under
    normal termination -- NaN K in NPF, or a NaN head in a CHD record -- and MF6
    answers a NaN budget while reporting success. The legacy ``KFromVector`` filled
    a null with the layer default (``fillna``, kflow.py); substituting a number for
    a value the user did not give is its own guess, so this names the row instead.
    """

    if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)):
        found = "empty (a GeoPackage NULL)" if value is None else f"{value!r}"
        raise ValueError(
            f"GeoPackage {source}: feature {label} has no usable number for "
            f"{_field_name(field)} -- it is {found}. MF6 would be written a "
            "non-finite value and report success. Fill the attribute, or drop the "
            "feature from the layer."
        )
    if not np.isfinite(value):
        raise ValueError(
            f"GeoPackage {source}: feature {label} has {value} for "
            f"{_field_name(field)}. MF6 would be written a non-finite value and "
            "report success. Fix the attribute, or drop the feature from the layer."
        )
    return value


def _numbered_index(
    value: Any, *, base: int, count: int | None, noun: str, owner: str, base_keyword: str
) -> tuple[int | None, str]:
    """``(zero-based index, "")`` for one numbered attribute value, or ``(None, why)``.

    The ONE reading of a GeoPackage column that NUMBERS something -- a layer
    (:func:`_layer_index`) or a stress period (:func:`_period_index`). A value is a
    WHOLE number at or above ``base`` -- and, where ``count`` is known, below
    ``base + count``. An integral float (``2.0``, what a column holding a NULL reads
    back as) and a numpy integer are that number; text that parses as one is too. A
    fraction is refused rather than truncated by ``int()`` into a number nobody wrote,
    and a NULL rather than raising numpy's bare "cannot convert float NaN". ``noun``
    names what is numbered, ``owner`` what has them, ``base_keyword`` the argument
    that says how the column counts.
    """

    if value is None:
        return None, "it is empty (a GeoPackage NULL)"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None, f"{value!r} is not a {noun} number"
    if not np.isfinite(number):
        return None, "it is empty (a GeoPackage NULL)" if np.isnan(number) else f"it is {number}"
    if not number.is_integer():
        return None, f"{value!r} is not a whole {noun} number"
    index = int(number) - int(base)
    top = None if count is None else int(count) + int(base) - 1
    if index < 0 or (count is not None and index >= int(count)):
        span = f"{int(base)} to {top}" if top is not None else f"{int(base)} up"
        # A zero base is also how a caller with no `<noun>_base=` asks (SFR's
        # reach_layer), so it is not named there.
        numbered = f"with {base_keyword}={int(base)} they run" if int(base) else "they run"
        return None, f"{int(number)}, a {noun} {owner} does not have -- {numbered} {span}"
    return index, ""


def _layer_index(value: Any, *, base: int, nlay: int | None) -> tuple[int | None, str]:
    """``(zero-based model layer, "")`` for one ``layer_field`` value, or ``(None, why)``.

    The reading of a feature's layer column for :func:`array_from_polygons`, every
    :class:`GeoPackageSource` builder and SFR's ``reach_layer``, by
    :func:`_numbered_index`. The polygon map checked only the upper bound, so a
    feature numbered below ``base`` matched no layer and painted NOTHING in silence:
    measured, two zones at ``lay`` 0 and 1 with ``layer_base=1`` gave the first zone's
    cells the default and no warning (ledger 393). A fractional value was truncated by
    ``int()`` into a layer nobody wrote, and a NULL raised numpy's bare "cannot convert
    float NaN".
    """

    return _numbered_index(
        value, base=base, count=nlay, noun="layer", owner="this model",
        base_keyword="layer_base",
    )


def _period_index(value: Any, *, base: int, nper: int) -> tuple[int | None, str]:
    """``(zero-based stress period, "")`` for one ``period_field`` value, or ``(None, why)``.

    The layer column's rule (:func:`_numbered_index`), with the simulation's ``nper``
    as the upper bound. ``GeoPackageSource`` read ``int(value)``, so ``1.5`` with
    ``period_base=1`` put the feature in period 0 beside the one declared for it, and a
    NULL raised numpy's bare "cannot convert float NaN to integer" (ledger 426).
    """

    return _numbered_index(
        value, base=base, count=nper, noun="stress period", owner="this simulation",
        base_keyword="period_base",
    )


def _why_it_paints_nothing(geometry, cells, recipe: Polygons) -> str:
    """Why a polygon that reaches the grid sets no cell of the array, in a few words.

    For :func:`array_from_polygons`' warning. A polygon meeting the cells only along an
    edge or at a point covers none of any cell's area, whatever the rule; otherwise the
    rule says it: ``centroid`` -- it holds no cell's centre; ``majority`` -- its value is
    outcovered in every cell it touches (by another value, or with ``default=`` by the
    part no polygon covers); ``mean`` -- only a zero-area touch is possible.
    """

    near = cells.geometry.iloc[cells.sindex.query(geometry, predicate="intersects")]
    if not near.intersection(geometry).area.gt(0.0).any():
        return "it only touches the grid's cells along an edge or at a point"
    if recipe.rule == "centroid":
        return "it holds no cell's centre"
    if recipe.default is None:
        return "its value is outcovered by another value in every cell it touches"
    return (
        "its value is outcovered in every cell it touches, by another value or by the "
        "default's vote -- the part no polygon covers, with any polygon carrying the "
        "default's value"
    )


def array_from_polygons(recipe: Polygons, grid, *, nlay: int, parameter: str) -> np.ndarray:
    """Resolve a :class:`~myflopy.sources.Polygons` recipe to an ``(nlay, ncpl)`` array.

    The build-time half of ``mf.npf(k=mf.Polygons(...))``. It reads the file, refuses a
    NULL naming the feature (or paints the recipe's one ``value`` into every polygon),
    paints each cell by the recipe's own ``rule``, and fills what no polygon reached
    with ``default`` -- or refuses when there is none. ``majority`` votes by VALUE: the
    features carrying one value are united before the values compete (ledger 553).
    With ``default`` the part of a cell no polygon covers COMPETES as the default
    (``majority``: a value takes a cell only when it covers more of it than that
    remainder plus any polygon carrying the default's own value; ``mean``: the
    remainder is averaged in at the default). A feature that reaches the grid but
    paints no cell is warned by name, with why, under every rule and with or without
    ``default`` (ledgers 551 and 585). Products and per-layer lists
    reach it through :func:`resolve_array_recipe`.

    **Refusing an unfilled cell is the point, not a convenience.** The array goes
    straight into MF6, where a NaN is answered with a NaN budget under "Normal
    termination" -- the failure ledger 247 exists to close. ``k_array`` cannot make
    that mistake because its ``defaults`` argument is required; here ``default=`` is
    optional, so the absence has to be caught.

    Layer handling mirrors ``k_array``'s except in the default: with no
    ``layer_field`` the map is AREAL and is broadcast to every layer, which is what a
    geology map digitized in plan view means. ``k_array`` requires the column.

    A feature that lands on no cell follows ledger 211 -- all of them missing raises
    with both extents and both CRSs, some of them warns by name -- through the same
    :func:`~myflopy.modflow.mf6.grid.geometry.miss_reason` and
    :func:`~myflopy.modflow.mf6.grid.geometry.where_they_are` every other miss in this
    library is worded by.
    """

    path = Path(recipe.path)
    read_kwargs = {} if recipe.layer is None else {"layer": recipe.layer}
    gdf = gpd.read_file(path, **read_kwargs)
    source_crs = gdf.crs
    grid_crs = getattr(grid, "crs", None)
    if grid_crs is not None and gdf.crs is not None and gdf.crs != grid_crs:
        gdf = gdf.to_crs(grid_crs)

    source = f"'{path.name}'" + (f" (layer {recipe.layer!r})" if recipe.layer else "")
    if recipe.column is not None and recipe.column not in gdf.columns:
        import difflib

        close = difflib.get_close_matches(recipe.column, list(gdf.columns), n=1, cutoff=0.6)
        raise ValueError(
            f"{source} has no column {recipe.column!r} to read {parameter} from"
            + (f" (did you mean {close[0]!r}?)" if close else "")
            + f". Its columns are: {', '.join(map(str, gdf.columns))}."
        )
    if len(gdf) == 0:
        raise ValueError(
            f"{source} has no features, so {parameter} would be its default everywhere."
        )

    cells = grid.gdf_vorPolys
    labels = [
        repr(str(row["name"])) if "name" in gdf.columns else f"#{index}"
        for index, row in gdf.iterrows()
    ]
    if recipe.column is None:
        # `value=`: one number for every polygon, finite by construction.
        numbers = np.full(len(gdf), float(recipe.value), dtype=float)
    else:
        numbers = np.array(
            [
                float(_finite(value, field=recipe.column, label=label, source=source))
                for value, label in zip(gdf[recipe.column], labels, strict=True)
            ],
            dtype=float,
        )

    painted = np.full((int(nlay), len(cells)), np.nan, dtype=float)
    # Whether each feature decides (under `majority`: carries the value that wins; under
    # `mean`: contributes to) any cell. With `default=` the part of a cell no polygon
    # covers competes as the default's vote (ledger 551), so a polygon smaller than
    # half of every cell it touches can be on the grid and paint nothing.
    painting = np.zeros(len(gdf), dtype=bool)
    if recipe.layer_field is None:
        values, _, painting = polygon_cell_values(
            cells, gdf, numbers, rule=recipe.rule, default=recipe.default, painters=True
        )
        painted[:] = values
    else:
        if recipe.layer_field not in gdf.columns:
            raise ValueError(
                f"{source} has no column {recipe.layer_field!r} to read the model "
                f"layer from. Its columns are: {', '.join(map(str, gdf.columns))}."
            )
        resolved = [
            _layer_index(one, base=recipe.layer_base, nlay=nlay)
            for one in gdf[recipe.layer_field]
        ]
        unusable = [
            f"{label} ({why})"
            for label, (layer, why) in zip(labels, resolved, strict=True)
            if layer is None
        ]
        if unusable:
            raise ValueError(
                f"{source}: {len(unusable)} of {len(gdf)} features name no model layer "
                f"in {recipe.layer_field!r}, so they would paint nothing into "
                f"{parameter}: {_listed(unusable)}. Correct the column, or set "
                "`layer_base=` to how its layers are numbered."
            )
        assigned = np.array([layer for layer, _ in resolved], dtype=int)
        for index in range(int(nlay)):
            chosen = assigned == index
            if not chosen.any():
                continue
            values, _, painters = polygon_cell_values(
                cells, gdf[chosen], numbers[chosen], rule=recipe.rule,
                default=recipe.default, painters=True,
            )
            painted[index] = values
            painting[np.flatnonzero(chosen)[painters]] = True

    # Per-FEATURE reach, which the per-CELL answer above cannot report: a feature can
    # lose every cell it touches to a larger neighbour under `majority` and still be on
    # the grid, which is not a miss. Intersection is the honest test of "reaches the
    # grid at all" (ledger 211).
    hit, _ = cells.sindex.query(gdf.geometry, predicate="intersects")
    landed = np.zeros(len(gdf), dtype=bool)
    landed[np.unique(hit)] = True
    where = where_they_are(
        gdf, cells, source_crs=source_crs, grid_crs=grid_crs,
        no_crs_subject="the file has NO CRS, so its coordinates",
    )
    if not landed.any():
        raise ValueError(
            f"No feature in {source} lands on the grid, so {parameter} would be its "
            f"default everywhere. {where}"
        )
    if not landed.all():
        missed = [
            f"{labels[position]} ({miss_reason(cells, gdf.geometry.iloc[position], target='cell')})"
            for position in np.flatnonzero(~landed)
        ]
        listing = ", ".join(missed[:10])
        if len(missed) > 10:
            listing += f", and {len(missed) - 10} more"
        logger.warning(
            "%d of %d features in %s land on no cell and paint nothing into %s: %s",
            len(missed), len(gdf), source, parameter, listing,
        )
        warnings.warn(
            f"{len(missed)} of {len(gdf)} features in {source} land on no cell and "
            f"paint nothing into {parameter}: {listing}.",
            UserWarning,
            stacklevel=stacklevel_outside_package(),
        )
    # A feature on the grid that sets no cell is warned by name under EVERY rule, with
    # or without `default=` (ledger 585): a unit outvoted in every cell it
    # touches by a larger neighbour vanished from the array in silence without one.
    silent = np.flatnonzero(landed & ~painting)
    if silent.size:
        reasons = [
            f"{labels[position]} ({_why_it_paints_nothing(gdf.geometry.iloc[position], cells, recipe)})"
            for position in silent[:10]
        ]
        listing = ", ".join(reasons)
        if silent.size > 10:
            listing += f", and {silent.size - 10} more"
        given = "" if recipe.default is None else f" with default={float(recipe.default)!r}"
        message = (
            f"{silent.size} of {len(gdf)} features in {source} reach the grid but "
            f"set no cell of {parameter} under rule={recipe.rule!r}{given}: {listing}. "
            "A finer mesh there, or rule='mean' to blend by area, keeps them."
        )
        logger.warning("%s", message)
        warnings.warn(message, UserWarning, stacklevel=stacklevel_outside_package())

    unfilled = ~np.isfinite(painted)
    if unfilled.any():
        if recipe.default is None:
            layers_short = sorted({int(one) for one in np.flatnonzero(unfilled.any(axis=1))})
            raise ValueError(
                f"{int(unfilled.sum())} of {painted.size} cells are reached by no "
                f"polygon in {source}, so {parameter} has no value for them "
                f"(layer(s) {', '.join(map(str, layers_short))}). MF6 would be written "
                "a non-finite value and report success. Pass `default=` to say what an "
                "uncovered cell means, or extend the polygons to cover the grid. "
                f"{where}"
            )
        painted[unfilled] = float(recipe.default)
    return painted


def array_from_raster(recipe: RasterValues, grid, *, parameter: str) -> np.ndarray:
    """Resolve a :class:`~myflopy.sources.RasterValues` recipe to an ``(ncpl,)`` array.

    The build-time half of ``mf.npf(k=mf.RasterValues(...))``, through the ONE raster
    sampler, :func:`~myflopy.modflow.mf6.grid.surfaces.get_raster_vals_at_centroids` --
    with the default ``mean="arithmetic"`` it returns exactly what
    ``Surface.raster(path).values(grid)`` does. A cell left without a value is never
    handed on: ``gaps="nearest"`` fills it from its nearest valued cell, ``default=``
    with that number, and otherwise it is REFUSED naming the count (ledger 247: MF6
    answers a NaN with a NaN budget under Normal termination). A raster giving NO cell
    a value is refused whatever the recipe says, naming both extents and both CRSs --
    ``default=`` everywhere would paint a misplaced raster in silence.

    Parameters
    ----------
    recipe : RasterValues
        What to sample; its ``path`` must already be anchored.
    grid : VoronoiGridPlus
        The grid to sample onto.
    parameter : str
        How the parameter is named in an error, e.g. ``"npf.k layer 2"``.

    Returns
    -------
    numpy.ndarray of float, shape (ncpl,)
        Every value finite.

    Raises
    ------
    ValueError
        A band the file lacks; a zero or negative pixel under a geometric or harmonic
        mean; cells left without a value as described above.
    """

    path = Path(recipe.path)
    sampled = get_raster_vals_at_centroids(
        grid, [path], ["value"], method=recipe.method, nodata=recipe.nodata,
        gaps="ignore", band=recipe.band, mean=recipe.mean,
    )
    values = sampled["value"].to_numpy(dtype=float)
    missing = np.flatnonzero(~np.isfinite(values))
    if missing.size == 0:
        return values
    if missing.size == values.size:
        raise ValueError(
            f"Raster '{path.name}' gave no value for ANY of the {values.size} cells, so "
            f"{parameter} would be entirely {'its default' if recipe.default is not None else 'no-data'}. "
            f"{raster_extents_text(grid, path)}"
        )
    if recipe.gaps == "nearest":
        xs, ys = np.asarray(grid.centroids[0]), np.asarray(grid.centroids[1])
        logger.info(
            "%s: filled %d of %d cells off raster '%s' from the nearest valued cell (%s)",
            parameter, missing.size, values.size, path.name, _listed([str(cell) for cell in missing]),
        )
        return _fill_nearest(values, xs, ys, missing)
    if recipe.default is not None:
        logger.info(
            "%s: %d of %d cells off raster '%s' take default=%g (%s)",
            parameter, missing.size, values.size, path.name, recipe.default,
            _listed([str(cell) for cell in missing]),
        )
        filled = values.copy()
        filled[missing] = recipe.default
        return filled
    raise ValueError(
        f"Raster '{path.name}' gives no value for {missing.size} of {values.size} cells, "
        f"so {parameter} has none there: cells {_listed([str(cell) for cell in missing])}. MF6 would be "
        "written a non-finite value and report success. Pass `gaps='nearest'` to take "
        "each from its nearest valued cell, or `default=` to say what an uncovered cell "
        "means."
    )


def _recipe_key(recipe: Any) -> str:
    """A hashable key for one recipe's declaration, for the per-build sample memo."""

    import json

    return json.dumps(_json_value(recipe), sort_keys=True)


def _resolve_factor(
    factor: Any, grid, *, nlay: int, areal: bool, parameter: str, cache: dict | None
) -> np.ndarray | float:
    """One recipe (or a number) resolved to ``(ncpl,)``, ``(nlay, ncpl)`` or a float."""

    if not isinstance(factor, (RasterValues, Polygons)):
        return float(factor)
    key = None
    if cache is not None:
        # Only a polygon map that paints model layers depends on the slot's shape; a
        # raster or an areal map is one value per cell wherever it is named.
        layered = isinstance(factor, Polygons) and factor.layer_field is not None
        key = (_recipe_key(factor), id(grid)) + ((nlay, areal) if layered else ())
        if key in cache:
            return cache[key]
    if isinstance(factor, RasterValues):
        resolved = array_from_raster(factor, grid, parameter=parameter)
    elif areal and factor.layer_field is not None:
        raise ValueError(
            f"{parameter} takes ONE value per column, and Polygons(layer_field="
            f"{factor.layer_field!r}) paints model layers. Drop layer_field= here."
        )
    elif areal or factor.layer_field is None:
        # An areal map paints every layer alike, so one layer is the whole answer.
        resolved = array_from_polygons(factor, grid, nlay=1, parameter=parameter)[0]
    else:
        resolved = array_from_polygons(factor, grid, nlay=nlay, parameter=parameter)
    if cache is not None:
        cache[key] = resolved
    return resolved


def resolve_array_recipe(
    value: Any,
    grid,
    *,
    nlay: int,
    parameter: str,
    areal: bool = False,
    cache: dict | None = None,
) -> np.ndarray:
    """Resolve an array recipe -- one, a product, or a per-layer list -- on ``grid``.

    The ONE resolver behind every array slot that takes :class:`~myflopy.sources.
    RasterValues`, :class:`~myflopy.sources.Polygons` or an
    :class:`~myflopy.sources.ArrayProduct` (ledger 541). A single recipe
    in a LAYERED slot is broadcast to every layer, as an areal polygon map always was;
    a per-layer LIST fills one layer per element, each a number or a recipe; a product
    multiplies its factors cell by cell after resolving each.

    Parameters
    ----------
    value : RasterValues, Polygons, ArrayProduct, or list of (those or float)
        The declaration. A list is accepted only when ``areal`` is ``False``.
    grid : VoronoiGridPlus
        The grid to resolve on.
    nlay : int
        The model's layer count. A list of another length raises ``ValueError``.
    parameter : str
        How the parameter is named in an error, e.g. ``"npf.k"``.
    areal : bool, default False
        ``True`` for a slot holding ONE value per column (``mf.rch.array``'s
        ``recharge``, ``mf.evt.array``'s fields): the result is ``(ncpl,)``.
    cache : dict, optional
        A memo shared across one build: each distinct recipe is sampled once however
        many periods or slots name it. ``None`` samples every time.

    Returns
    -------
    numpy.ndarray of float
        ``(ncpl,)`` when ``areal``, else ``(nlay, ncpl)``.

    Raises
    ------
    ValueError
        A per-layer list of the wrong length, and whatever resolving a recipe raises
        (:func:`array_from_raster`, :func:`array_from_polygons`).
    """

    nlay = int(nlay)
    ncpl = len(grid.gdf_vorPolys)

    def one(item: Any, label: str, *, item_areal: bool) -> np.ndarray:
        factors = item.factors if isinstance(item, ArrayProduct) else (item,)
        result: np.ndarray | float = 1.0
        for factor in factors:
            resolved = _resolve_factor(
                factor, grid, nlay=nlay, areal=item_areal, parameter=label, cache=cache
            )
            result = result * np.asarray(resolved, dtype=float)
        shape = (ncpl,) if item_areal else (nlay, ncpl)
        return np.broadcast_to(np.asarray(result, dtype=float), shape).copy()

    if isinstance(value, (list, tuple)):
        if areal:
            raise ValueError(
                f"{parameter} takes ONE value per column; a per-layer list means nothing "
                "there."
            )
        if len(value) != nlay:
            raise ValueError(
                f"{parameter} is a per-layer list of {len(value)} entries, and the model "
                f"has {nlay} layer(s). Give one entry -- a number or a recipe -- per "
                "layer."
            )
        return np.vstack(
            [
                one(item, f"{parameter} layer {index}", item_areal=True)
                for index, item in enumerate(value)
            ]
        )
    return one(value, parameter, item_areal=areal)


def _field_name(field: RowValue) -> str:
    """How to name the value spec that produced a bad number, in an error."""

    if isinstance(field, Spread):
        return f"'{field.field}' (spread)"
    if isinstance(field, CellSurfaceOffset):
        # A surface reference resolves to a number the file did not hold, so name the
        # columns it DID read -- a null offset is what makes the sum non-finite.
        columns = [c for c in (field.offset, field.minimum) if isinstance(c, str)]
        named = f"the {field.reference} surface reference"
        if columns:
            named += " (fields " + ", ".join(repr(c) for c in columns) + ")"
        return named
    if isinstance(field, (list, tuple)):
        return "its per-period fields " + ", ".join(repr(str(f)) for f in field)
    return repr(str(field))


def _metadata_value(value: RowValue):
    """Return a manifest-friendly representation of one value specification.

    Every dataclass in :data:`RowValue` needs a branch here. A run manifest is
    written with :func:`json.dumps` and no ``default=``, so a spec that reaches it
    unconverted fails at ``prepare_run`` -- long after the package built cleanly.
    ``tests/test_geopackage_spread.py`` asserts the whole union round-trips.

    **This is NOT ``to_dict``, deliberately** (ledger 274). Two converters, two jobs:
    ``to_dict`` round-trips, omits defaults and is revived into a live object;
    this one is human-readable provenance, spells every field out, and stays a plain
    dict. Delegating here to ``to_dict`` was measured wrong twice over --
    ``PackageSpec.from_dict`` runs ``metadata`` through ``_spec_value``, so a
    ``$``-tagged value would come back a LIVE OBJECT and change what
    ``metadata['fields']`` means for every consumer; and omitting defaults would make
    a manifest a reader consults sparser. They are tied together by a field-coverage
    test instead of by shared code.
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
            # `warn_below` was the one field of four this dropped, while the test
            # pinning the shape claimed it "keeps enough to reproduce the split". It
            # did not: whether a feature losing 40% of itself warns depended on a
            # number no persisted artifact recorded. `metadata` is excluded from the
            # digest, so adding it drifts nothing (ledger 274).
            "warn_below": value.warn_below,
        }
    return value


#: Coverage at or above this counts as whole -- a feature is not "clipped" by float noise.
_FULLY_COVERED = 1.0 - 1e-9

#: The ``(package, record field)`` pairs the library knows are EXTENSIVE, so a total of
#: what was NOT written means something. Deliberately one entry: a well's rate is
#: unambiguously a flux. A conductance is not -- ledger 169 measured
#: ``corr(cond, cell area) ~ 0`` over 205 real DRN records, so it is not a leakance and
#: summing it says nothing. Everything else reports per feature, with no total, unless
#: the caller wrapped it in ``mf.Spread``, which IS a declaration that it is extensive.
#: Deliberately NOT ``_fields_of_kind("extensive")``, which is EIGHT fields. Two
#: different questions share the word: the registry's ``kind`` asks "do several records
#: on one cell ADD?" (a conductance does), while this asks "does the DECLARED column
#: hold a total we can report a shortfall against?" -- and an unwrapped conductance does
#: not, because it is copied WHOLE to every covered cell, so declared-minus-written
#: calls a clean build a -300% loss (ledger 169/288). Collapsing them would put a false
#: loss figure on every drn/ghb/riv warning (ledger 297).
_EXTENSIVE_FIELDS = {("wel", "rate")}

#: The ``(package, GPKG-ARGUMENT field)`` pairs that are ELEVATIONS -- ``drn.elevation``
#: and ``ghb.head``, which are the names the HELPERS take, NOT the record names
#: ``drn.elev`` and ``ghb.bhead``. Both sets here and both their call sites use the
#: helper namespace; three of the seven packages differ between the two, so a rewrite
#: keyed on ``FieldSpec`` names would silently match NOTHING for those three rather
#: than failing (ledger 297). Derived from the registry below so the list cannot drift
#: from the descriptor that declares each field's kind. ``Spread`` divides a
#: feature's total across the cells it covers, which is meaningful for a flux or a
#: conductance and MEANINGLESS for an elevation: measured, ``elevation=mf.Spread("elev")``
#: on a declared 75.0 wrote 2.06 to 11.87 -- the elevation multiplied by each cell's area
#: share -- with all 11 records below their own cell bottom and no warning (ledger 289).
#: Ledger 274 refused exactly this on the ARRAY form and never reached the list form.
#: ``rch.recharge`` is deliberately absent: a rate column that really holds a total volume
#: IS spreadable, which is what ``Spread``'s own docstring teaches.
def _fields_of_kind(kind: str) -> frozenset[tuple[str, str]]:
    """``(package, gpkg-argument field)`` for every field the registry calls ``kind``.

    The registry declares kinds against RECORD names; these sets are keyed by the
    helper's argument names. The two are positionally aligned -- ``_bc_from_features``
    asserts exactly that before building a record -- so the translation is a zip, and
    doing it here once is what keeps the namespaces from being confused by hand.
    """

    found = set()
    for package, descriptor in _PACKAGE_EXPLORER_SPECS.items():
        gpkg_names = tuple(descriptor.gpkg_defaults or ())
        for position, spec in enumerate((descriptor.inputs or {}).values()):
            if spec.kind == kind and position < len(gpkg_names):
                found.add((package, gpkg_names[position]))
    return frozenset(found)


_ELEVATION_FIELDS = _fields_of_kind("elevation")


def refuse_spread_on_elevation(package: str, values, *, verb: str) -> None:
    """Refuse ``Spread`` on a field that is an ELEVATION, at the CALL.

    Needs no grid, so it belongs at declaration -- the rule `_lazy_bc_spec` already
    follows for a missing file and ``nper >= 1``. Called from BOTH paths: the eager one
    reaches `_bc_from_features`, and a lazy declaration never does, so a single site
    would have covered only half the API (ledger 289). ``verb`` is the call the message
    names -- the route actually taken, ``mf.<pkg>.gpkg`` or
    ``GeoPackageSource.<pkg>``.
    """

    # `values` is keyed by the HELPER's argument names, which is the namespace
    # `_ELEVATION_FIELDS` uses. Named for it, because calling it `record_field` is what
    # made the other namespace look authoritative here (ledger 297).
    for gpkg_field, spec in (values or {}).items():
        if isinstance(spec, Spread) and (package, gpkg_field) in _ELEVATION_FIELDS:
            raise TypeError(
                f"{verb}: {gpkg_field}=mf.Spread(...) divides a feature's "
                f"total across the cells it covers, and {gpkg_field} is an ELEVATION -- "
                "there is no total to divide, so each cell would get a FRACTION of the "
                "elevation you declared (measured: a declared 75.0 written as 2.06 to "
                f"11.87). Pass the column directly ({gpkg_field}={spec.field!r}), or "
                "mf.CellSurfaceOffset(...) to place it relative to a surface."
            )


def refuse_mixed_spread_modes(package: str, values, *, verb: str) -> None:
    """Refuse two :class:`Spread` fields of one call that disagree on ``mode``, at the CALL.

    A feature has ONE set of cell shares, so a call has one mode; the share step read
    the FIRST wrapped field's and dropped the rest in silence (``mf.evt.gpkg`` with
    ``rate`` clipped and ``depth`` retained wrote the depth clipped). Refused rather
    than resolved per field -- the two disagreeing is a declaration nobody meant
    (ledger 396). Called from the eager and the lazy path, as
    :func:`refuse_spread_on_elevation` is, and names ``verb`` as it does.
    ``min_share`` and ``warn_below`` are not refused: their combination (largest,
    smallest) is documented and loses nothing.
    """

    spreads = {
        name: spec for name, spec in (values or {}).items() if isinstance(spec, Spread)
    }
    modes = {spec.mode for spec in spreads.values()}
    if len(modes) > 1:
        named = ", ".join(f"{name}=Spread({spec.field!r}, mode={spec.mode!r})"
                          for name, spec in spreads.items())
        raise ValueError(
            f"{verb}: {named} -- one feature's cells take ONE set of shares, "
            "so every Spread in a call must use the same mode. Pick 'clip' (the part "
            "off the grid is not written) or 'retained' (the total is kept on the part "
            "that remains) for all of them."
        )


def check_placement_options(
    *, min_above_bottom: float | None, snap_within: float | None, verb: str
) -> None:
    """Refuse a negative ``min_above_bottom`` or a non-positive ``snap_within``, at the CALL.

    Neither needs a grid, so both belong at declaration. ONE function for the eager
    ``GeoPackageSource`` and the lazy ``mf.<pkg>.gpkg`` declaration: the lazy branch
    never reaches ``__post_init__``, and the copy it kept checked ``snap_within`` but not
    the clearance, so ``min_above_bottom=-1.0`` was accepted with no ``context=`` and
    failed at build.
    """

    if min_above_bottom is not None and float(min_above_bottom) < 0:
        raise ValueError(
            f"{verb}: min_above_bottom is a CLEARANCE above the cell bottom and cannot be "
            f"negative, not {min_above_bottom!r}. 0.0 means as close to the bottom as the "
            "written file can express."
        )
    if snap_within is not None and not (float(snap_within) > 0):
        raise ValueError(
            f"{verb}: snap_within is a distance in the grid's own coordinate units and "
            f"must be positive, not {snap_within!r}. Leave it None for one cell width, "
            "measured on the cell a feature would snap to."
        )


def _listed(items, limit: int = 10) -> str:
    """The first ``limit`` entries, then a count -- one wording for every listing here."""

    listing = ", ".join(items[:limit])
    return listing + (f", and {len(items) - limit} more" if len(items) > limit else "")


@dataclass
class _Coverage:
    """What became of every feature of one source: the message AND the record.

    A warning is best-effort -- `warnings.warn` dedupes, and on the lazy path it fires
    deep inside `prepare_run` where no stacklevel points at the declaration -- so the
    same numbers are stamped into `PackageSpec.metadata`, which is excluded from
    `package_digest` and therefore cannot report a stored run as drifted (ledger 288).
    """

    features: int = 0
    landed: int = 0
    snapped: list = field(default_factory=list)      # (label, distance)
    missed: list = field(default_factory=list)       # (index, row)
    clipped: list = field(default_factory=list)      # (label, covered fraction)
    below_bottom: list = field(default_factory=list)  # (label, record field, depth below)
    floored: list = field(default_factory=list)       # (label, record field, lift applied)
    lost: dict = field(default_factory=dict)         # extensive field -> quantity not written

    def quiet(self) -> bool:
        """Whether nothing worth saying happened."""

        return not (self.snapped or self.missed or self.clipped or self.below_bottom or self.floored)

    def summary(self) -> dict:
        """The durable record: counts, and the snap distances that reveal a systematic error.

        The MAX snap distance is the one that matters -- a uniform datum shift and
        ordinary digitizing slop are indistinguishable by the DIRECTION of their snap
        vectors (measured: both concentrate at R = 1.000), and separable only by
        magnitude.
        """

        record = {
            "features": self.features,
            "landed": self.landed,
            "snapped": len(self.snapped),
            "missed": len(self.missed),
            "clipped": len(self.clipped),
        }
        if self.below_bottom:
            record["below_bottom"] = len(self.below_bottom)
            record["worst_below_bottom"] = float(min(d for _, _, d in self.below_bottom))
        if self.floored:
            record["floored"] = len(self.floored)
            record["max_lift"] = float(max(lift for _, _, lift in self.floored))
        if self.snapped:
            distances = sorted(distance for _, distance in self.snapped)
            record["snap_max"] = float(distances[-1])
            record["snap_median"] = float(distances[len(distances) // 2])
        if self.lost:
            record["lost"] = {key: float(value) for key, value in self.lost.items()}
        return record


@dataclass(slots=True)
class GeoPackageSource:
    """Map one GeoPackage layer onto a model grid and build reusable inputs.

    Layer values are assumed to be one-based by default because that is the
    usual convention in external GIS inputs. Stress-period values are assumed
    to be zero-based.

    Parameters
    ----------
    path : Path or str
        The GeoPackage, or any file ``geopandas.read_file`` opens. Must exist when the
        source is created (``FileNotFoundError``). It is read on first use and cached
        on this object. Its features are reprojected to ``grid.crs`` when the file and
        the grid both declare a CRS and the two differ, and used as they are otherwise.
    context : ModelContext
        What the features are mapped onto. ``context.grid`` is required (``ValueError``
        without it) and supplies the cell polygons. ``context.domain`` decides which
        cells can take a record -- ``idomain > 0``; every cell when it is ``None`` --
        and where the active edge lies for ``edges_only``. The layer elevations
        (``context.surfaces``, else ``grid.gdf_topbtm``) serve
        :class:`CellSurfaceOffset` and ``min_above_bottom``.
    nper : int
        Number of stress periods. The boundary builders return one record list for
        each period ``0 .. nper - 1``, and with ``period_field`` unset every feature
        appears in every one. Must be at least 1 (``ValueError``).
    layer : str, optional
        The layer (table) to read, passed to ``geopandas.read_file(layer=)``. ``None``
        (the default) reads the file's default layer.
    name_field : str or None, default "name"
        Attribute column naming each feature. Its value identifies the feature in
        warnings and errors (the row index does when the column is absent or this is
        ``None``), and it becomes the record's boundname under ``boundnames=True``,
        which then requires the column (``ValueError`` naming it). ``None`` with
        ``boundnames=True`` writes the row index as the boundname.
    layer_field : str or None, default "layer"
        Attribute column giving each feature's model layer, read as
        ``value - layer_base``. A value that is not a whole number at or above
        ``layer_base`` -- a NULL, a fraction, one below the base -- raises
        ``ValueError`` naming the feature (:meth:`k_array` also refuses one past its
        ``nlay``). The column is required when this is set. ``None`` puts every
        feature in the top layer.
    period_field : str, optional
        Attribute column giving the ONE stress period a feature applies to, read as
        ``value - period_base``. A value that is not a whole number naming a period
        ``0 .. nper - 1`` -- a NULL, a fraction, one below ``period_base`` or past the
        last period -- raises ``ValueError`` naming the feature; ``2.0``, a numpy
        integer and the text ``"2"`` all read as 2, as in ``layer_field``. ``None``
        (the default) applies every feature to every period.
        :meth:`k_array` does not read it.
    layer_base : int, default 1
        Subtracted from each ``layer_field`` value: ``1`` reads one-based (GIS-style)
        layers, ``0`` zero-based ones.
    period_base : int, default 0
        Subtracted from each ``period_field`` value: ``0`` reads zero-based periods.
    gaps : {"warn", "nearest", "raise", "ignore"}, default "warn"
        What a feature reaching no usable cell means to the boundary builders.
        ``"warn"`` leaves it out and reports it -- a ``UserWarning`` and a log
        WARNING, in one message that also reports features clipped at the grid's
        edge, elevations below their own cell bottom and records lifted by
        ``min_above_bottom``. ``"nearest"`` moves a missed POINT onto the nearest
        usable cell within ``snap_within`` and reports the move; a line or polygon is
        never moved, and what cannot be moved is reported as under ``"warn"``.
        ``"raise"`` raises ``ValueError`` wherever ``"warn"`` would
        report. ``"ignore"`` says nothing. In every mode a source with no feature on
        a usable cell raises ``ValueError`` naming both extents and both CRSs. An
        unknown word raises ``ValueError`` here. :meth:`k_array` reads it the same
        way, a missed feature leaving its cells at ``defaults``.
    snap_within : float, optional
        ``gaps="nearest"`` only: the farthest a missed point may be moved, in the
        grid's coordinate units, measured to the cell POLYGON rather than its
        centroid. ``None`` (the default) is ``sqrt(area)`` of the cell it would move
        to. Must be positive (``ValueError``).
    min_above_bottom : float, optional
        Lift each elevation of a record to at least this far above its own cell
        bottom -- ``chd`` and ``ghb`` ``head``, ``drn`` ``elevation``, ``evt``
        ``surface``, and ``riv`` ``rbot``, whose ``stage`` is raised by the same
        amount only where it is itself below the new floor. The clearance is never
        less than ``max(|bottom| * 1e-6, 1e-6)``, so ``0.0`` means as close to the
        bottom as the written file can hold. ``None`` (the default) lifts nothing,
        but a record below its cell bottom is still reported (see ``gaps``). Both the
        report and the lift need layer elevations -- ``context.surfaces``, else
        ``grid.gdf_topbtm``. With neither, the report is skipped, and a value given
        here makes every builder with an elevation field raise ``ValueError`` saying
        where elevations come from, rather than write the records unlifted. Negative
        raises ``ValueError``.
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
    gaps: str = "warn"
    snap_within: float | None = None
    min_above_bottom: float | None = None
    _gdf: gpd.GeoDataFrame | None = field(default=None, init=False, repr=False)
    #: The CRS ``_gdf`` was aligned to, and the layer as read from the file. A grid
    #: relabelled after the first read is followed from the file's own frame rather
    #: than left on the old label (ledger 466).
    _gdf_crs: Any = field(default=None, init=False, repr=False)
    _raw: gpd.GeoDataFrame | None = field(default=None, init=False, repr=False)
    _edge_cells: dict[int, set[int]] = field(default_factory=dict, init=False, repr=False)
    _source_crs: Any = field(default=None, init=False, repr=False)
    #: What became of the features on the last build, for the metadata stamp.
    last_coverage: Any = field(default=None, init=False, repr=False, compare=False)
    #: The call a declaration-time refusal names, formatted with ``package=``. The
    #: ``mf.<pkg>.gpkg`` helpers set ``"mf.{package}.gpkg"`` on the source they build,
    #: so each route is named as the caller spelled it.
    _verb: str = field(default="GeoPackageSource.{package}", init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        """Coerce ``path`` and validate the file exists, a grid is present, and ``nper >= 1``."""

        self.path = Path(self.path)
        if not self.path.exists():
            raise FileNotFoundError(f"GeoPackage not found: {self.path}")
        if self.context.grid is None:
            raise ValueError("ModelContext.grid is required for GeoPackage mapping.")
        if self.nper < 1:
            raise ValueError("nper must be at least 1.")
        check_gaps(self.gaps, verb="GeoPackageSource")
        check_placement_options(
            min_above_bottom=self.min_above_bottom,
            snap_within=self.snap_within,
            verb="GeoPackageSource",
        )

    @property
    def grid(self):
        """Grid helper carried by the model context."""

        return self.context.grid

    @property
    def gdf(self) -> gpd.GeoDataFrame:
        """The GeoPackage layer, reprojected to the model grid's CRS as it is NOW.

        Read once and cached; the grid's CRS is read on every access, so a grid
        relabelled after the first read gets the layer reprojected to its new label
        from the file's own frame -- it stayed on the old one (ledger 466). The
        same object comes back while the grid's CRS is unchanged.
        """

        grid_crs = getattr(self.grid, "crs", None)
        if self._gdf is not None and self._gdf_crs is None:
            # A frame set from outside (`source._gdf = frame`) is taken to be in the
            # grid's CRS when it is first read, as `Boundaries` takes one (ledger 431).
            self._gdf_crs = grid_crs
        if self._gdf is not None and not _same_crs(self._gdf_crs, grid_crs):
            logger.debug("grid relabelled %s -> %s; reprojecting the GeoPackage layer",
                         self._gdf_crs, grid_crs)
            if self._raw is None:  # set from outside: it is all there is to move
                self._raw = self._gdf
            self._gdf = None
        if self._gdf is None:
            if self._raw is None:
                kwargs = {} if self.layer is None else {"layer": self.layer}
                self._raw = gpd.read_file(self.path, **kwargs)
                self._source_crs = self._raw.crs
            gdf = self._raw
            if grid_crs is not None and gdf.crs is not None and not _same_crs(gdf.crs, grid_crs):
                gdf = gdf.to_crs(grid_crs)
            self._gdf = gdf
            self._gdf_crs = grid_crs
        return self._gdf

    def _active(self, layer: int, cell: int) -> bool:
        """Whether ``(layer, cell)`` can carry a boundary, per the context's idomain.

        The three-valued idomain rule itself is
        :func:`myflopy.modflow.mf6.grid.geometry.domain_active`.
        """

        return domain_active(self.context.domain, layer, cell)

    def _layer(self, index: Any, row, *, nlay: int | None = None) -> int:
        """The zero-based grid layer for a feature (from ``layer_field`` minus ``layer_base``).

        Read by :func:`_layer_index`, the rule the polygon map shares; a value naming
        no layer raises naming the feature. ``nlay`` is checked where the caller knows
        it (:meth:`k_array`); the boundary builders leave the upper bound to the
        build's cellid range check.
        """

        if self.layer_field is None:
            return 0
        layer, why = _layer_index(row[self.layer_field], base=self.layer_base, nlay=nlay)
        if layer is None:
            raise ValueError(
                f"GeoPackage '{self.path.name}': feature {self._feature_label(index, row)} "
                f"names no model layer in {self.layer_field!r}: {why}. Correct the "
                "column, or set `layer_base=` to how its layers are numbered."
            )
        return layer

    def _periods(self, index: Any, row) -> range | tuple[int]:
        """The stress periods a feature applies to: all periods, or its single ``period_field``.

        Read by :func:`_period_index`, the layer column's rule: a value naming no stress
        period of the simulation -- a NULL, a fraction, one outside ``period_base`` to
        ``period_base + nper - 1`` -- raises naming the feature (ledger 426).
        """

        if self.period_field is None:
            return range(self.nper)
        period, why = _period_index(
            row[self.period_field], base=self.period_base, nper=self.nper
        )
        if period is None:
            raise ValueError(
                f"GeoPackage '{self.path.name}': feature {self._feature_label(index, row)} "
                f"names no stress period in {self.period_field!r}: {why}. Correct the "
                "column, or set `period_base=` to how its periods are numbered."
            )
        return (period,)

    def _cells(self, geometry, *, layer: int, edges_only: bool = False) -> list[int]:
        """The active grid cells a feature's ``geometry`` covers (optionally edge cells only)."""

        return [cell for cell, _ in self._cell_shares(geometry, layer=layer, edges_only=edges_only)]

    def _cell_shares(
        self, geometry, *, layer: int, edges_only: bool = False
    ) -> list[tuple[int, float]]:
        """Active cells a feature covers, each with its fraction of the whole feature.

        The measure, the tie rule and the sliver rule are
        :func:`myflopy.modflow.mf6.grid.geometry.cell_shares`, so that every
        caller asking which cells a geometry reaches answers alike. What belongs
        to this source is WHICH cells are eligible -- one layer's active cells,
        and only the grid edge under ``edges_only`` -- and both are resolved once
        per feature here rather than per candidate.
        """

        if geometry is None or geometry.is_empty:
            # Before the edge set, which walks every cell and reads the grid's own
            # idomain_path if it has one: a row with no geometry reaches nothing.
            return []
        edge_cells = self._edge_set(layer) if edges_only else None
        return cell_shares(
            self.grid.gdf_vorPolys,
            geometry,
            usable=lambda cells: self._usable_cells(cells, layer=layer, edge_cells=edge_cells),
        )

    def _snapped_shares(self, geometry, *, layer: int, edges_only: bool):
        """A zero-cell POINT moved onto its nearest usable cell: ``(shares, distance)``.

        ``(None, distance)`` when it cannot be moved -- too far, or not a point.

        **Points only, and that is measured rather than conservative.** Translating a
        feature by the minimal vector onto the grid recovers a POINT exactly (share
        1.0000 at 3 ft off and at 30 ft off), and recovers a LINE only when the line
        happens to lie PARALLEL to the boundary it missed: a line of the same length
        POINTING at the grid lands with its tip touching and resolves to **zero cells**,
        and a polygon translated by exactly its gap only touches too (ledger 288). That
        is orientation luck, not a rule. A point has no extent, so moving it is the one
        case where "the same feature, slightly over" is well defined; anything with a
        footprint would have that footprint changed or destroyed by the move.

        The tolerance is compared against the distance to the cell POLYGON, not its
        centroid -- "just off the edge" is the natural reading, and on one real miss the
        two differ by 3.07x (3.00 ft to the polygon, 9.21 ft to the centroid). The
        raster sibling uses centroid distance; the divergence is deliberate.
        """

        from shapely.affinity import translate
        from shapely.ops import nearest_points

        if geometry is None or geometry.is_empty or _extent(geometry) > 0:
            return None, None

        polys = self.grid.gdf_vorPolys
        edge_cells = self._edge_set(layer) if edges_only else None
        usable = self._usable_cells(list(polys.index), layer=layer, edge_cells=edge_cells)
        if not usable:
            return None, None
        candidates = polys.geometry.loc[usable]
        distances = candidates.distance(geometry)
        cell = distances.idxmin()
        distance = float(distances.loc[cell])

        # `None` means one cell width, measured on the cell it would snap TO. A grid-wide
        # statistic is the wrong instrument on a ragged mesh: cell areas span 3.9x on a
        # plain Voronoi grid and 18.5x on a refined one, so one number is ~2x too
        # generous beside the smallest cells and ~0.3x too tight beside the largest, in
        # the same model.
        tolerance = (
            float(self.snap_within)
            if self.snap_within is not None
            else float(np.sqrt(candidates.loc[cell].area))
        )
        if distance > tolerance:
            return None, distance

        start, end = nearest_points(geometry, candidates.loc[cell])
        moved = translate(geometry, end.x - start.x, end.y - start.y)
        shares = self._cell_shares(moved, layer=layer, edges_only=edges_only)
        return (shares or None), distance

    def _edge_set(self, layer: int) -> set[int]:
        """One layer's ACTIVE-domain edge cells, found once per layer.

        ``edges_only`` means the edge of the ACTIVE domain -- a cell ringing an
        inactive hole is on it (ledger 236) -- so ``get_grid_edge`` has to be told
        which cells are inactive. It was called with no argument, which left it
        the grid's own ``idomain_path`` and nothing else: a ``ModelContext``
        domain, the documented package-first route, never reached the edge test.
        On a 5x5 grid with one inactive cell that wrote 16 cells where the active
        edge is 21, and with the outer ring inactive it refused the source
        outright rather than using the ring inside it.

        The edge is per LAYER because activity is (``_active`` already is), and
        the answer is cached per layer because ``get_grid_edge`` walks every cell:
        a ten-feature source used to pay for all of that ten times.

        The context's idomain is the model's declared domain, so it WINS over the
        grid's ``idomain_path``; a grid built from such a file still works, because
        with no context domain nothing is passed and that fallback stands.
        """

        if layer not in self._edge_cells:
            domain = self.context.domain
            inactive = (
                None
                if domain is None
                else [
                    cell
                    for cell in range(self.grid.ncpl)
                    if not domain_active(domain, layer, cell)
                ]
            )
            self._edge_cells[layer] = set(self.grid.get_grid_edge(idomain=inactive))
        return self._edge_cells[layer]

    def _usable_cells(self, cells, *, layer: int, edge_cells: set[int] | None) -> list[int]:
        """The active cells among ``cells``, lowest-numbered first."""

        return sorted(
            cell
            for cell in cells
            if (edge_cells is None or cell in edge_cells) and self._active(layer, cell)
        )

    def surface_value(self, reference: SurfaceReference, *, layer: int, cell: int) -> float:
        """Return a model/cell surface value for a mapped boundary cell.

        Delegates to the shared :class:`SurfaceResolver` so the ``gdf_topbtm``
        column layout lives in exactly one place (also used by ``mf.evt(context=...)``).

        Parameters
        ----------
        reference : {"model_top", "cell_top", "cell_bottom"}
            Which surface to read: the ground surface, or the top or bottom of the
            cell in ``layer``. Any case; another word raises ``ValueError`` naming the
            three -- the rule :class:`CellSurfaceOffset` applies (ledger 392).
        layer : int
            Zero-based model layer. ``"cell_top"`` of layer ``n > 0`` is the bottom of
            layer ``n - 1``; ``"model_top"`` ignores it.
        cell : int
            Zero-based cell (CELL2D) number.

        Returns
        -------
        float
            The elevation, read from ``context.surfaces`` or, when that is ``None``,
            ``grid.gdf_topbtm``.

        Raises
        ------
        ValueError
            ``reference`` is not one of the three; neither table is available; or the
            table has no column for ``reference`` at ``layer``. A ``cell`` outside the table (or, in a dict table, a
            ``layer`` past its bottoms) is not checked: the lookup's own
            ``KeyError``/``IndexError`` propagates.
        """

        return SurfaceResolver(self.context).surface_value(reference, layer=layer, cell=cell)

    def _feature_label(self, index: Any, row) -> str:
        """A feature as a user would recognise it: its name if it has one, else its row."""

        if self.name_field is not None and self.name_field in row.index:
            return repr(str(row[self.name_field]))
        return f"#{index}"

    def _where(self) -> str:
        """Where the features are and where the grid is -- what an empty result is acted on with."""

        return where_they_are(
            self.gdf,
            self.grid.gdf_vorPolys,
            source_crs=self._source_crs,
            grid_crs=getattr(self.grid, "crs", None),
            no_crs_subject="the file has NO CRS, so its coordinates",
        )

    def _report(self, coverage, *, package: str, name: str, edges_only: bool, empty: str) -> None:
        """Say what happened to every feature -- refuse, warn, or stay quiet, per ``gaps``.

        **The all-missed refusal is ORTHOGONAL to the mode**, ``"ignore"`` included. The
        raster sibling returns early on ``"ignore"`` and so lets a wholly-missing input
        through, because ``Surface.sample`` has a documented next stage that fills it.
        A boundary package has no such caller: an empty one is rejected by MF6 naming
        the file, or runs without the boundary under ``maxbound=``. Ledger 211 decided
        that and this does not reopen it in one word (ledger 288).

        **Both channels, deliberately.** ``logger.warning`` stays for a configured
        application, and ``warnings.warn`` is added because it is the only one a default
        caller hears: ``import myflopy`` installs a ``NullHandler`` on the ``myflopy``
        logger, root has none, and ``logging.callHandlers`` therefore counts a handler
        as found and never reaches ``logging.lastResort`` -- measured, the partial-miss
        line produced zero bytes on stderr in a bare session and in a real ipykernel.

        The message carries the PACKAGE NAME and the quantity lost. That is the
        accounting, and it is also what stops `warnings.warn`'s per-(text, category,
        module, lineno) dedup collapsing one source driving four packages from four
        warnings to one -- the four lines used to be byte-identical.
        """

        source = f"'{self.path.name}'" + (f" (layer {self.layer!r})" if self.layer else "")
        target = "active grid-edge cell" if edges_only else "active cell"
        subject = package if name in ("", package) else f"{package} {name!r}"
        if coverage.features == 0:
            raise ValueError(f"GeoPackage {source} has no features; {empty}.")

        if coverage.missed and len(coverage.missed) == coverage.features:
            polys = self.grid.gdf_vorPolys
            reasons = [
                f"{self._feature_label(index, row)} "
                f"({miss_reason(polys, row.geometry, target=target)})"
                for index, row in coverage.missed
            ]
            raise ValueError(
                f"No feature in GeoPackage {source} lands on an {target}; {empty}. "
                f"{_listed(reasons)}. {self._where()}"
                + ("" if self.gaps != "ignore" else
                   " gaps='ignore' does not reach this: it quiets a feature that "
                   "misses, not a source none of whose features lands.")
            )

        if self.gaps == "ignore" or coverage.quiet():
            return

        sentence = self._sentence(coverage, target=target)
        if self.gaps == "raise":
            raise ValueError(
                f"GeoPackage {source} -> {subject}: {sentence} {self._where()} "
                "Pass gaps='warn' to be told instead, or gaps='nearest' to move a POINT "
                "that only just missed."
            )
        message = f"GeoPackage {source} -> {subject}: {sentence}"
        logger.warning("%s", message)
        warnings.warn(message, UserWarning, stacklevel=stacklevel_outside_package())

    def _sentence(self, coverage, *, target: str) -> str:
        """What happened, in one line that is DISTINCT per package.

        The package name and the lost quantity are in the message on purpose. Four
        packages driven from one GeoPackage used to emit four byte-identical lines, and
        `warnings.warn` dedupes on (text, category, module, lineno) -- so a default
        caller saw one of four. Making each line say what it is about fixes the dedup
        and IS the accounting, in one change.
        """

        polys = self.grid.gdf_vorPolys
        parts: list[str] = []
        if coverage.missed:
            reasons = [
                f"{self._feature_label(index, row)} "
                f"({miss_reason(polys, row.geometry, target=target)})"
                for index, row in coverage.missed
            ]
            parts.append(
                f"{len(coverage.missed)} of {coverage.features} features land on no "
                f"{target} and are left out: {_listed(reasons)}"
            )
        if coverage.snapped:
            moved = [f"{label} ({distance:,.3g})" for label, distance in coverage.snapped]
            parts.append(
                f"{len(coverage.snapped)} were MOVED onto the grid by gaps='nearest' "
                f"(distance in grid units): {_listed(moved)}"
            )
        if coverage.clipped:
            clipped = [f"{label} ({covered:.1%} inside)" for label, covered in coverage.clipped]
            parts.append(
                f"{len(coverage.clipped)} lie only partly inside the active grid and are "
                f"CLIPPED to the part that is: {_listed(clipped)}"
            )
        if coverage.below_bottom:
            worst = min(depth for _, _, depth in coverage.below_bottom)
            # Grouped by FEATURE: one drain line crossing six offending cells is one
            # problem the user can act on, not six identical listings.
            grouped: dict[tuple, list] = {}
            for label, name, depth in coverage.below_bottom:
                grouped.setdefault((label, name), []).append(depth)
            named = [
                f"{label} ({name}, {len(depths)} cell{'s' if len(depths) > 1 else ''}, "
                f"worst {min(depths):,.3g})"
                for (label, name), depths in grouped.items()
            ]
            parts.append(
                f"{len(coverage.below_bottom)} records place an ELEVATION BELOW their own "
                f"cell bottom (worst {worst:,.3g}): {_listed(named)} -- MF6 refuses these "
                "at icelltype != 0, and at icelltype 0 accepts them and over-drains; pass "
                "min_above_bottom= to lift them"
            )
        if coverage.floored:
            by_feature: dict[tuple, list] = {}
            for label, name, lift in coverage.floored:
                by_feature.setdefault((label, name), []).append(lift)
            lifted = [
                f"{label} ({name}, {len(lifts)} cell{'s' if len(lifts) > 1 else ''}, "
                f"up to +{max(lifts):,.3g})"
                for (label, name), lifts in by_feature.items()
            ]
            parts.append(
                f"{len(coverage.floored)} records were LIFTED to min_above_bottom: "
                f"{_listed(lifted)}"
            )
        if coverage.lost:
            lost = ", ".join(f"{key} {value:,.4g}" for key, value in sorted(coverage.lost.items()))
            parts.append(f"not written: {lost}")
        return "; ".join(parts) + "."

    def _value(
        self,
        row,
        value: RowValue,
        *,
        period: int,
        layer: int,
        cell: int,
        share: float = 1.0,
        label: str = "",
    ):
        """Resolve one value spec for a cell, refusing a value MF6 cannot use.

        ``share`` is the cell's fraction of the feature and is used only by
        :class:`Spread`; everything else resolves the same value for every cell,
        which is what an elevation or a head should do.

        The result is checked by :func:`_finite`, because a GeoPackage NULL used to
        travel all the way into the written package as a NaN (ledger 247).
        """

        if isinstance(value, CellSurfaceOffset):
            resolved = value.resolve(self, row, layer=layer, cell=cell)
        elif isinstance(value, Spread):
            raw = _finite(
                _period_value(row, value.field, period),
                field=value, label=label, source=f"'{self.path.name}'",
            )
            resolved = float(raw) * share
        elif isinstance(value, (int, float)):
            resolved = value
        else:
            resolved = _period_value(row, value, period)
        return _finite(resolved, field=value, label=label, source=f"'{self.path.name}'")

    def _boundary_data(
        self,
        *fields: RowValue,
        edges_only: bool = False,
        boundnames: bool = False,
        package: str = "",
        name: str = "",
    ) -> dict[int, list[list[Any]]]:
        """Build MF6 stress-period data by mapping every feature to its cells and value fields.

        Validates the required columns exist, then for each feature emits one
        ``[(layer, cell), *values(, boundname)]`` record per applicable period and
        covered cell, keyed by stress period.
        """

        coverage = _Coverage()
        # Detection needs the layer elevations. Where the model carries none there is
        # nothing to compare against, so it is skipped rather than raising -- a package
        # whose model declares no surfaces is not an error. A floor the caller ASKED for
        # is another matter: skipping it wrote every record exactly as declared, below
        # its cell bottom, with nothing said (ledger 395).
        elevation_positions = self._elevation_positions(package)
        elevations_available = bool(package) and _layer_elevations(self.context) is not None
        if elevation_positions and self.min_above_bottom is not None and not elevations_available:
            raise ValueError(
                f"GeoPackage '{self.path.name}' -> {package}: min_above_bottom="
                f"{self.min_above_bottom!r} lifts each elevation above its own cell "
                "bottom, and this model carries no layer elevations to read the bottom "
                "from. They come from ModelContext(surfaces=...) -- a LayerBuildResult or "
                "a {'top': ..., 'botm': ...} table -- or from the grid's gdf_topbtm; a "
                "package declared WITHOUT context= also reads the built discretization's "
                "top and botm. Supply one, or drop min_above_bottom."
            )
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
        missed: list[tuple[Any, Any]] = []
        for index, row in self.gdf.iterrows():
            layer = self._layer(index, row)
            label = self._feature_label(index, row)
            coverage.features += 1
            shares = self._cell_shares(row.geometry, layer=layer, edges_only=edges_only)
            if not shares and self.gaps == "nearest":
                shares, distance = self._snapped_shares(
                    row.geometry, layer=layer, edges_only=edges_only
                )
                if shares:
                    coverage.snapped.append((label, distance))
            if not shares:
                missed.append((index, row))
                self._account_for_missing(coverage, row, fields, package=package, covered=0.0)
                continue
            coverage.landed += 1
            covered = sum(share for _, share in shares)
            if covered < _FULLY_COVERED:
                coverage.clipped.append((label, covered))
                self._account_for_missing(
                    coverage, row, fields, package=package, covered=covered
                )
            if spreads:
                shares = self._resolved_shares(shares, spreads, index)
            boundname = (
                str(row[self.name_field])
                if boundnames and self.name_field is not None
                else str(index)
            )
            for period in self._periods(index, row):
                for cell, share in shares:
                    values = [
                        self._value(
                            row, field, period=period, layer=layer, cell=cell,
                            share=share, label=self._feature_label(index, row),
                        )
                        for field in fields
                    ]
                    if package == "riv":
                        # Needs no elevations, so it runs whether or not the model has
                        # any; it sat inside the floor and was skipped with it.
                        self._refuse_stage_below_rbot(
                            values, positions=elevation_positions, cell=cell, label=label
                        )
                    if elevations_available:
                        # AFTER the `gaps="nearest"` snap, never before: a snapped point
                        # lands in a DIFFERENT cell with a different bottom, so a floor
                        # applied earlier would use the wrong one (ledger 290).
                        values = self._floor_elevations(
                            values, package=package, layer=layer, cell=cell,
                            coverage=coverage, label=label,
                        )
                    record = [(layer, cell), *values]
                    if boundnames:
                        record.append(boundname)
                    data[period].append(record)
        coverage.missed = missed
        self._report(
            coverage,
            package=package,
            name=name,
            edges_only=edges_only,
            empty="the package would have no records -- MF6 rejects it, or runs "
            "without it under maxbound=",
        )
        self.last_coverage = coverage
        return data

    def _floor_elevations(self, values, *, package, layer, cell, coverage, label):
        """Keep an elevation inside its cell, and SAY when one is not -- per cell.

        **What MF6 does with a below-bottom elevation depends on `icelltype`, and both
        branches are live in this library** (ledger 290). At `icelltype = 0` -- which is
        what `mf.npf(k=...)` writes, because it names no icelltype at all -- all four
        packages ACCEPT one, terminate normally, and over-drain: measured 2.17x on a
        3x3 bench and 1.09x on a real 64-cell DISV model. At `icelltype != 0` -- which
        `CanonicalModelConfig` and 17 taught sites use -- MF6 REFUSES by name, NEWTON does
        not change it, and **the refusal is a per-period data read**, so an offender in
        period 3 dies after MF6 has already solved periods 1 and 2. On a 72-period model
        that is a late failure.

        So the DETECTION is default-on and the CORRECTION is opt-in. Nothing reported this
        before -- no warning, no metadata key, no validator -- and the silence is the part
        that is constant across models; the magnitude is not.

        **A floor exactly AT the bottom is not safe**, which is why the clearance is
        nudged. FloPy writes `botm` fixed-point and a list-BC value in E format, so a value
        set exactly equal to its cell bottom can round BELOW it in the written file and
        MF6's strict comparison trips: measured, 13 of 24 real runs FAILED at bottoms drawn
        from [40, 600], and which way the discarded digits fall is a coin flip. A CHD set
        exactly at its bottom is worse -- it aborts the SOLVER ("CONSTANT-HEAD CELL WENT
        DRY"), at round magnitudes too. The relative nudge is the one already at
        `myflopy.modflow.mf6.boundary_support.bottom_clearance`, which the MODFLOW-USG
        importer floors by as well -- one rule, since ledger 291 found the two copies
        had drifted apart.
        """

        positions = self._elevation_positions(package)
        if not positions:
            return values
        bottom = self.surface_value("cell_bottom", layer=layer, cell=cell)
        nudge = bottom_clearance(bottom)

        # RIV is a PAIR, not two independent elevations. MF6 checks `rbot` against the cell
        # bottom -- and names `stage` too when the stage is itself below it -- and enforces
        # `stage >= rbot` separately, so raising `rbot` alone turns a model that ran into a
        # hard failure whenever the stage was below the bottom as well. Measured, one
        # record: as declared q=-868.57; rbot alone REFUSED at both icelltypes; both
        # clamped q=-380.00 with the river depth collapsed to 0.000 from a declared 1.0;
        # carrying the stage up q=-108.57. Three readings, 3.5x apart.
        # A stage that is ALREADY inside its cell is left alone -- see `river_pair_lift`,
        # where lifting it unconditionally measured a reversed sign (ledger 291). A stage
        # below its rbot never reaches here: `_refuse_stage_below_rbot` ran first.
        pair = ("rbot", "stage") if package == "riv" else ()
        checked = ["rbot"] if pair else list(positions)

        for name in checked:
            index = positions[name]
            value = float(values[index])
            if value < bottom:
                coverage.below_bottom.append((label, name, value - bottom))
            if self.min_above_bottom is None:
                continue
            floor = bottom + max(float(self.min_above_bottom), nudge)
            if value >= floor:
                continue
            coverage.floored.append((label, name, floor - value))
            if pair:
                # Carry `stage` up by the same delta when it is itself below the floor,
                # so the declared river depth survives; leave a legal stage exactly where
                # it was drawn. One helper, because the MODFLOW-USG importer floors RIV
                # records by the same rule and used to carry its own copy of it -- which
                # clamped `stage` and left `rbot` below the bottom (ledger 291).
                stage = positions["stage"]
                values[stage], values[index] = river_pair_lift(
                    float(values[stage]), value, floor
                )
            else:
                values[index] = floor
        return values

    @staticmethod
    def _elevation_positions(package: str) -> dict[str, int]:
        """``{helper field: record position}`` for each ELEVATION field of ``package``."""

        if not package:
            return {}
        fields = _PACKAGE_EXPLORER_SPECS[package].gpkg_defaults
        return {name: i for i, name in enumerate(fields) if (package, name) in _ELEVATION_FIELDS}

    def _refuse_stage_below_rbot(self, values, *, positions, cell, label) -> None:
        """Refuse a RIV record whose stage lies below its own rbot.

        MF6 refuses this at BOTH icelltypes ("RIV BOUNDARY (n) RIVER STAGE ... IS LESS
        THAN RIVER BOTTOM"), and ``mf.riv.gpkg``'s own docstring promises the invariant
        (ledger 290). It compares two numbers of one record and needs no layer
        elevations, so it runs on every RIV record -- it used to sit inside the floor
        and was skipped with it on a model carrying none (ledger 395).
        """

        stage_now = float(values[positions["stage"]])
        rbot_now = float(values[positions["rbot"]])
        if stage_now < rbot_now:
            raise ValueError(
                f"GeoPackage '{self.path.name}': feature {label} gives cell {cell} a "
                f"river stage ({stage_now:,.4g}) BELOW its river bottom "
                f"({rbot_now:,.4g}). MF6 refuses that whatever icelltype is, and "
                "min_above_bottom cannot fix it -- a river's stage is its water "
                "surface and its rbot the bed beneath. Swap the two columns, or "
                "correct the attribute."
            )

    def _account_for_missing(self, coverage, row, fields, *, package: str, covered: float) -> None:
        """Add a feature's unwritten quantity to the running total, where that MEANS something.

        **Only where the library knows the field is extensive**, which is two cases: a
        field wrapped in :class:`Spread` -- the caller said so, and ``mode="clip"`` is
        the mode that actually discards the outside part -- and a well's rate, which is
        unambiguously a flux. Everything else is left to the per-feature listing with no
        total.

        A declared-minus-written identity over every field is measurably nonsense: an
        unwrapped conductance is copied WHOLE to each of the twelve cells a polygon
        covers, so a clean build reports -300% "loss" (ledger 169 measured
        ``corr(cond, cell area) ~ 0`` over 205 real records -- a conductance is not a
        leakance and does not divide by area).
        """

        # The HELPER's argument names (`drn.elevation`), not the record names
        # (`drn.elev`). This local was called `record_fields`, which is the other
        # namespace and reads as authoritative -- it is what sent a reader to conclude
        # the branch below was dead code, since `wel`'s record field is `q` and the set
        # says `rate` (ledger 297). Both are the helper namespace and they match.
        gpkg_fields = tuple(_PACKAGE_EXPLORER_SPECS[package].gpkg_defaults) if package else ()
        for position, spec in enumerate(fields):
            gpkg_field = gpkg_fields[position] if position < len(gpkg_fields) else None
            if isinstance(spec, Spread):
                column = spec.field if spec.mode == "clip" else None
            elif (package, gpkg_field) in _EXTENSIVE_FIELDS and isinstance(spec, str):
                column = spec
            else:
                column = None
            if column is None or column not in row:
                continue
            value = row[column]
            try:
                quantity = float(value)
            except (TypeError, ValueError):
                continue
            if not np.isfinite(quantity):
                continue
            key = f"{gpkg_field}={column}" if gpkg_field else column
            coverage.lost[key] = coverage.lost.get(key, 0.0) + quantity * (1.0 - covered)

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
            # The mesh these records' cell numbers refer to. Every list BC reaches
            # this one method, so pinning here pins all seven.
            GRID_PIN: grid_pin(self.grid),
            # The option this pin speaks for: an edit replacing the records voids it,
            # while renaming the file or toggling `save_flows` does not (ledger 275).
            # Measured across 98 (package, override) pairs on the seven list BCs: 91
            # left the built package's cellids byte-identical, and the 7 that moved
            # them were all `stress_period_data`.
            GRID_PIN_OPTIONS: ["stress_period_data"],
            "source_type": "geopackage",
            "source_path": str(self.path),
            # What the file CONTAINED when this package was declared. `source_path`
            # alone cannot answer "have my inputs changed since?" -- measured, a
            # library package built from a GeoPackage that was later edited reports
            # `unchanged` forever, and `drift()`'s own docstring promises that means
            # reproducible (ledger 276).
            SOURCE_SHA: source_fingerprint(self.path),
            "source_layer": self.layer,
            "builder": method,
            "fields": {name: _metadata_value(value) for name, value in fields.items()},
            "name_field": self.name_field,
            "layer_field": self.layer_field,
            "period_field": self.period_field,
            "layer_base": self.layer_base,
            "period_base": self.period_base,
            "gaps": self.gaps,
            "snap_within": self.snap_within,
            "min_above_bottom": self.min_above_bottom,
            # What became of the features -- counts, and the snap distances. A warning
            # is best-effort: `warnings.warn` dedupes, and a lazy declaration resolves
            # inside `prepare_run` where no stacklevel points at the caller. This
            # survives. `metadata` is excluded from `package_digest`, so recording it
            # cannot report a stored run as drifted (ledger 288).
            "coverage": dict(self.last_coverage.summary()) if self.last_coverage else None,
        }

    def _recipe(
        self, package: str, *, boundnames: bool, edges_only: bool, **fields: RowValue
    ) -> BoundarySource:
        """The lazy recipe for this EAGER declaration, recorded beside its records.

        An eager `.gpkg` resolves at declaration and holds finished records, so a grid
        move classifies it ``stale`` -- and that is the LARGEST stale class there is,
        because declaring with a ``context=`` is what every pre-B2 doc teaches. It does
        not need an overlay engine to move: everything needed to resolve it again is
        right here at the call, so it is written down and ``with_grid`` re-lazifies from
        it, EXACTLY, where an overlay would only approximate (ledger 296).

        **Recorded, not re-derived.** ``metadata['fields']`` looks like it already holds
        this, and reconstructing from it was the obvious plan -- but ``_metadata_value``
        is deliberately NOT a round-trip format (ledger 274, which measured two reasons
        to keep it apart from ``to_dict``), and it omits ``nper``, ``edges_only`` and
        ``boundnames`` entirely. ``edges_only`` is the one that matters: it decides
        WHICH CELLS are kept (ledger 246), so guessing it wrong writes a different
        package in silence.

        It rides in ``metadata``, which ``package_digest`` excludes, so recording it
        moves no stored digest -- and metadata is run through ``_spec_value`` on the way
        back, so the source revives as a live object rather than a plain dict.
        """

        return BoundarySource(
            path=self.path,
            package=package,
            nper=self.nper,
            values=dict(fields),
            layer=self.layer,
            name_field=self.name_field,
            layer_field=self.layer_field,
            period_field=self.period_field,
            layer_base=self.layer_base,
            period_base=self.period_base,
            gaps=self.gaps,
            snap_within=self.snap_within,
            min_above_bottom=self.min_above_bottom,
            edges_only=edges_only,
            boundnames=boundnames,
            source_sha=source_fingerprint(self.path),
        )

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

        verb = self._verb.format(package=package)
        refuse_spread_on_elevation(package, fields, verb=verb)
        refuse_mixed_spread_modes(package, fields, verb=verb)

        data_options: dict[str, Any] = {"boundnames": boundnames}
        if descriptor.capabilities.edges_only:
            data_options["edges_only"] = edges_only
        elif edges_only:
            raise ValueError(f"{package} features cannot be restricted to grid edges.")

        return _SPEC_FACTORIES[package](
            self._boundary_data(*fields.values(), package=package, name=name, **data_options),
            name=name,
            boundnames=boundnames,
            **(options or {}),
        ).with_metadata(
            **self._metadata(package, **fields),
            **{MOVE_RECIPE: self._recipe(package, boundnames=boundnames,
                                         edges_only=edges_only, **fields)},
        )

    def chd(
        self,
        *,
        head: RowValue = "head",
        name: str = "chd",
        edges_only: bool = False,
        boundnames: bool = True,
        **options,
    ) -> PackageSpec:
        """Return a CHD package spec from GeoPackage features.

        Each feature becomes one ``[(layer, cell), head(, boundname)]`` record per
        stress period it applies to, in every usable cell it covers.

        Parameters
        ----------
        head : str, float, list of str, or CellSurfaceOffset, default "head"
            The specified head, in model length units: a column name; a number used
            for every feature; a list of column names indexed by stress period
            (periods past its end read the last one); or a :class:`CellSurfaceOffset`,
            resolved per cell. Copied to every cell a feature covers, and lifted by
            ``min_above_bottom`` when that is set. A NULL or a non-number raises
            ``ValueError`` naming the feature.
        name : str, default "chd"
            Package name, passed to :func:`~myflopy.advanced.chd_spec`; it also
            names the package in this source's warnings.
        edges_only : bool, default False
            Keep only cells on the edge of the ACTIVE domain in each feature's layer
            -- those bordering the outside of the mesh or an inactive cell.
        boundnames : bool, default True
            End each record with a boundname (see ``name_field``) and write MF6's
            ``BOUNDNAMES`` option. ``False`` writes neither.
        **options
            Forwarded to :func:`~myflopy.advanced.chd_spec` and on to
            :class:`flopy.mf6.ModflowGwfchd`; an unknown name raises ``TypeError``.

        Returns
        -------
        PackageSpec
            Records keyed by period ``0 .. nper - 1``. Its ``metadata`` carries the
            grid pin, the file's path, content hash and this source's settings, a
            coverage summary, and the recipe ``with_grid`` re-declares it from.

        Raises
        ------
        ValueError
            A column this call reads is missing; a value is NULL or non-numeric;
            ``min_above_bottom`` is set and the model carries no layer elevations; no
            feature reaches a usable cell; or ``gaps="raise"`` found something to
            report.
        TypeError
            ``head`` is a :class:`Spread` -- a head has no total to divide.
        """

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
        """Return a GHB package spec from GeoPackage features.

        Each feature becomes one ``[(layer, cell), head, conductance(, boundname)]``
        record per stress period it applies to, in every usable cell it covers.

        Parameters
        ----------
        head : str, float, list of str, or CellSurfaceOffset, default "head"
            The boundary head (MF6's ``bhead``), in model length units: a column
            name; a number used for every feature; a list of column names indexed by
            stress period (periods past its end read the last one); or a
            :class:`CellSurfaceOffset`, resolved per cell. Copied to every cell a
            feature covers, and lifted by ``min_above_bottom`` when that is set.
        conductance : str, float, list of str, or Spread, default "conductance"
            The conductance (MF6's ``cond``), in length^2/time, as the same value
            forms. A plain value is copied WHOLE to every cell a feature covers;
            wrap the column in :class:`Spread` to divide a feature's total among its
            cells instead.
        name : str, default "ghb"
            Package name, passed to :func:`~myflopy.advanced.ghb_spec`; it also
            names the package in this source's warnings.
        edges_only : bool, default False
            Keep only cells on the edge of the ACTIVE domain in each feature's layer
            -- those bordering the outside of the mesh or an inactive cell.
        boundnames : bool, default True
            End each record with a boundname (see ``name_field``) and write MF6's
            ``BOUNDNAMES`` option. ``False`` writes neither.
        **options
            Forwarded to :func:`~myflopy.advanced.ghb_spec` and on to
            :class:`flopy.mf6.ModflowGwfghb`; an unknown name raises ``TypeError``.

        Returns
        -------
        PackageSpec
            Records keyed by period ``0 .. nper - 1``. Its ``metadata`` carries the
            grid pin, the file's path, content hash and this source's settings, a
            coverage summary, and the recipe ``with_grid`` re-declares it from.

        Raises
        ------
        ValueError
            A column this call reads is missing; a value is NULL or non-numeric;
            ``min_above_bottom`` is set and the model carries no layer elevations; no
            feature reaches a usable cell; or ``gaps="raise"`` found something to
            report.
        TypeError
            ``head`` is a :class:`Spread` -- a head has no total to divide.
        """

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
        """Return a DRN package spec from GeoPackage features.

        Each feature becomes one ``[(layer, cell), elevation, conductance(,
        boundname)]`` record per stress period it applies to, in every usable cell it
        covers.

        Parameters
        ----------
        elevation : str, float, list of str, or CellSurfaceOffset, default "elevation"
            The drain elevation (MF6's ``elev``), in model length units: a column
            name; a number used for every feature; a list of column names indexed by
            stress period (periods past its end read the last one); or a
            :class:`CellSurfaceOffset`, resolved per cell. Copied to every cell a
            feature covers, and lifted by ``min_above_bottom`` when that is set.
        conductance : str, float, list of str, or Spread, default "conductance"
            The conductance (MF6's ``cond``), in length^2/time, as the same value
            forms. A plain value is copied WHOLE to every cell a feature covers;
            wrap the column in :class:`Spread` to divide a feature's total among its
            cells instead.
        name : str, default "drn"
            Package name, passed to :func:`~myflopy.advanced.drn_spec`; it also
            names the package in this source's warnings.
        edges_only : bool, default False
            Keep only cells on the edge of the ACTIVE domain in each feature's layer
            -- those bordering the outside of the mesh or an inactive cell.
        boundnames : bool, default True
            End each record with a boundname (see ``name_field``) and write MF6's
            ``BOUNDNAMES`` option. ``False`` writes neither.
        **options
            Forwarded to :func:`~myflopy.advanced.drn_spec` and on to
            :class:`flopy.mf6.ModflowGwfdrn`; an unknown name raises ``TypeError``.

        Returns
        -------
        PackageSpec
            Records keyed by period ``0 .. nper - 1``. Its ``metadata`` carries the
            grid pin, the file's path, content hash and this source's settings, a
            coverage summary, and the recipe ``with_grid`` re-declares it from.

        Raises
        ------
        ValueError
            A column this call reads is missing; a value is NULL or non-numeric;
            ``min_above_bottom`` is set and the model carries no layer elevations; no
            feature reaches a usable cell; or ``gaps="raise"`` found something to
            report.
        TypeError
            ``elevation`` is a :class:`Spread` -- an elevation has no total to divide.
        """

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
        """Return a RIV package spec from GeoPackage features.

        Each feature becomes one ``[(layer, cell), stage, conductance, rbot(,
        boundname)]`` record per stress period it applies to, in every usable cell it
        covers.

        Parameters
        ----------
        stage : str, float, list of str, or CellSurfaceOffset, default "stage"
            The river stage, in model length units: a column name; a number used for
            every feature; a list of column names indexed by stress period (periods
            past its end read the last one); or a :class:`CellSurfaceOffset`,
            resolved per cell. Copied to every cell a feature covers. A stage below
            its record's ``rbot`` raises ``ValueError`` naming the feature, whether or
            not the model carries layer elevations.
        conductance : str, float, list of str, or Spread, default "conductance"
            The riverbed conductance (MF6's ``cond``), in length^2/time, as the same
            value forms. A plain value is copied WHOLE to every cell a feature
            covers; wrap the column in :class:`Spread` to divide a feature's total
            among its cells instead.
        rbot : str, float, list of str, or CellSurfaceOffset, default "rbot"
            The river bottom elevation, in model length units, as the same value
            forms as ``stage``. ``min_above_bottom`` lifts it, carrying ``stage`` up
            by the same amount only where the stage is itself below the new floor.
        name : str, default "riv"
            Package name, passed to :func:`~myflopy.advanced.riv_spec`; it also
            names the package in this source's warnings.
        edges_only : bool, default False
            Keep only cells on the edge of the ACTIVE domain in each feature's layer
            -- those bordering the outside of the mesh or an inactive cell.
        boundnames : bool, default True
            End each record with a boundname (see ``name_field``) and write MF6's
            ``BOUNDNAMES`` option. ``False`` writes neither.
        **options
            Forwarded to :func:`~myflopy.advanced.riv_spec` and on to
            :class:`flopy.mf6.ModflowGwfriv`; an unknown name raises ``TypeError``.

        Returns
        -------
        PackageSpec
            Records keyed by period ``0 .. nper - 1``. Its ``metadata`` carries the
            grid pin, the file's path, content hash and this source's settings, a
            coverage summary, and the recipe ``with_grid`` re-declares it from.

        Raises
        ------
        ValueError
            A column this call reads is missing; a value is NULL or non-numeric; a
            stage lies below its rbot; ``min_above_bottom`` is set and the model
            carries no layer elevations; no feature reaches a usable cell; or
            ``gaps="raise"`` found something to report.
        TypeError
            ``stage`` or ``rbot`` is a :class:`Spread` -- an elevation has no total
            to divide.
        """

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
        """Return a WEL package spec from GeoPackage features.

        Each feature becomes one ``[(layer, cell), rate(, boundname)]`` record per
        stress period it applies to, in every usable cell it covers.

        Parameters
        ----------
        rate : str, float, list of str, or Spread, default "rate"
            The well rate (MF6's ``q``), in length^3/time, negative for extraction:
            a column name; a number used for every feature; a list of column names
            indexed by stress period (periods past its end read the last one); or a
            :class:`Spread`, which divides a feature's rate among its cells. A plain
            value is copied WHOLE to every cell a feature covers, so a line or
            polygon feature pumps its rate once per cell. For a single column (plain,
            or a ``Spread`` in ``"clip"`` mode), a feature lying wholly or partly off
            the usable cells reports ``rate * (1 - fraction on them)`` as not
            written.
        name : str, default "wel"
            Package name, passed to :func:`~myflopy.advanced.wel_spec`; it also
            names the package in this source's warnings.
        boundnames : bool, default True
            End each record with a boundname (see ``name_field``) and write MF6's
            ``BOUNDNAMES`` option. ``False`` writes neither.
        **options
            Forwarded to :func:`~myflopy.advanced.wel_spec` and on to
            :class:`flopy.mf6.ModflowGwfwel`; an unknown name raises ``TypeError``.

        Returns
        -------
        PackageSpec
            Records keyed by period ``0 .. nper - 1``. Its ``metadata`` carries the
            grid pin, the file's path, content hash and this source's settings, a
            coverage summary, and the recipe ``with_grid`` re-declares it from.

        Raises
        ------
        ValueError
            A column this call reads is missing; a value is NULL or non-numeric; no
            feature reaches a usable cell; or ``gaps="raise"`` found something to
            report.
        """

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
        """Return a list-based RCH package spec from GeoPackage features.

        Each feature becomes one ``[(layer, cell), recharge(, boundname)]`` record per
        stress period it applies to, in every usable cell it covers.

        Parameters
        ----------
        recharge : str, float, list of str, or Spread, default "recharge"
            The recharge flux, in length/time: a column name; a number used for
            every feature; a list of column names indexed by stress period (periods
            past its end read the last one); or a :class:`Spread`, for a column that
            holds a feature's total rather than a rate. A plain value is copied to
            every cell a feature covers.
        name : str, default "rch"
            Package name, passed to :func:`~myflopy.advanced.rch_spec`; it also
            names the package in this source's warnings.
        boundnames : bool, default True
            End each record with a boundname (see ``name_field``) and write MF6's
            ``BOUNDNAMES`` option. ``False`` writes neither.
        **options
            Forwarded to :func:`~myflopy.advanced.rch_spec` and on to
            :class:`flopy.mf6.ModflowGwfrch`; an unknown name raises ``TypeError``.

        Returns
        -------
        PackageSpec
            Records keyed by period ``0 .. nper - 1``. Its ``metadata`` carries the
            grid pin, the file's path, content hash and this source's settings, a
            coverage summary, and the recipe ``with_grid`` re-declares it from.

        Raises
        ------
        ValueError
            A column this call reads is missing; a value is NULL or non-numeric;
            ``min_above_bottom`` is set and the model carries no layer elevations; no
            feature reaches a usable cell; or ``gaps="raise"`` found something to
            report.
        """

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

        Parameters
        ----------
        surface : str, float, list of str, or CellSurfaceOffset, default "surface"
            The ET surface elevation, in model length units: a column name; a number
            used for every feature; a list of column names indexed by stress period
            (periods past its end read the last one); or a
            :class:`CellSurfaceOffset`, resolved per cell (``"model_top"`` for land
            surface). Copied to every cell a feature covers, and lifted by
            ``min_above_bottom`` when that is set.
        rate : str, float, list of str, or Spread, default "rate"
            The maximum ET flux, in length/time, as the same value forms. A plain
            value is copied to every cell a feature covers.
        depth : str, float, or list of str, default "depth"
            The extinction depth below ``surface``, in model length units: a column
            name, a number used for every feature, or per-period column names.
            Copied to every cell a feature covers.
        name : str, default "evt"
            Package name, passed to :func:`~myflopy.advanced.evt_spec`; it also
            names the package in this source's warnings.
        boundnames : bool, default True
            End each record with a boundname (see ``name_field``) and write MF6's
            ``BOUNDNAMES`` option. ``False`` writes neither.
        **options
            Forwarded to :func:`~myflopy.advanced.evt_spec` and on to
            :class:`flopy.mf6.ModflowGwfevt`; an unknown name raises ``TypeError``.
            ``nseg`` other than 1 raises ``ValueError``, as above.

        Returns
        -------
        PackageSpec
            ``[(layer, cell), surface, rate, depth(, boundname)]`` records keyed by
            period ``0 .. nper - 1``. Its ``metadata`` carries the grid pin, the
            file's path, content hash and this source's settings, a coverage
            summary, and the recipe ``with_grid`` re-declares it from.

        Raises
        ------
        ValueError
            ``nseg`` is not 1; a column this call reads is missing; a value is NULL or
            non-numeric; no feature reaches a usable cell; or ``gaps="raise"`` found
            something to report.
        TypeError
            ``surface`` is a :class:`Spread` -- an elevation has no total to divide.
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
        """Return an ``(nlay, ncpl)`` K array from GeoPackage polygons.

        Each feature paints its value into every ACTIVE cell of its own layer
        (``layer_field``) that it covers, by the same cell rule the boundary builders
        use; where features overlap, the later row wins. Resolved now, against this
        source's context -- for a K map that resolves at build and survives
        ``Project.save()``, see :class:`~myflopy.sources.Polygons`.

        Parameters
        ----------
        value : str, default "k"
            Attribute column holding each feature's K, in length/time. A file without
            it raises ``ValueError`` naming the file and the column. A NULL or a
            non-number raises ``ValueError`` naming the feature; it is read only for
            features that reach a cell.
        nlay : int
            Number of model layers: the array's first dimension. A feature whose
            ``layer_field`` value is not a whole number from ``layer_base`` to
            ``layer_base + nlay - 1`` raises ``ValueError`` naming the feature.
        defaults : float or list of float
            K for every cell no feature paints: one value for every layer, or one per
            layer (any other length raises ``ValueError``).

        Returns
        -------
        numpy.ndarray
            Float array of shape ``(nlay, ncpl)``.

        Raises
        ------
        ValueError
            As above; the file lacks the ``layer_field`` column; no feature reaches an
            active cell (naming both extents and CRSs, whatever ``gaps`` says); or
            ``gaps="raise"`` and some feature misses.

        Notes
        -----
        A feature reaching no active cell is handled by this source's ``gaps``, as in
        the boundary builders: ``"warn"`` (the default) leaves its cells at
        ``defaults`` and names it in a ``UserWarning`` as well as the log;
        ``"nearest"`` first moves a missed POINT onto the nearest active cell within
        ``snap_within`` and reports the move; ``"raise"`` refuses; ``"ignore"`` says
        nothing. A polygon running past the grid is not reported -- a geology map
        routinely extends beyond the model.
        """

        fallback = [defaults] * nlay if np.isscalar(defaults) else list(defaults)
        if len(fallback) != nlay:
            raise ValueError("defaults must be a scalar or contain one value per layer.")
        required = {value} | ({self.layer_field} if self.layer_field is not None else set())
        missing = sorted(required - set(self.gdf.columns))
        if missing:
            raise ValueError(
                f"GeoPackage '{self.path.name}' is missing fields: {', '.join(missing)}"
            )
        result = np.array(
            [[float(fallback[layer])] * self.grid.ncpl for layer in range(nlay)],
            dtype=float,
        )
        # The same accounting and the same `_report` as the boundary builders, so `gaps=`
        # means one thing on this source whichever method reads it. The partial miss
        # went to `logger.warning` alone, which reaches nobody by default (ledger 288),
        # and `gaps=` was not read (ledger 394). CLIPPING is deliberately not
        # recorded: a geology map routinely runs past the model, and that is not news.
        coverage = _Coverage()
        for index, row in self.gdf.iterrows():
            layer = self._layer(index, row, nlay=nlay)
            label = self._feature_label(index, row)
            coverage.features += 1
            cells = self._cells(row.geometry, layer=layer)
            if not cells and self.gaps == "nearest":
                shares, distance = self._snapped_shares(row.geometry, layer=layer, edges_only=False)
                if shares:
                    coverage.snapped.append((label, distance))
                    cells = [cell for cell, _ in shares]
            if not cells:
                coverage.missed.append((index, row))
                continue
            coverage.landed += 1
            painted = _finite(
                row[value],
                field=value,
                label=label,
                source=f"'{self.path.name}'",
            )
            for cell in cells:
                result[layer, cell] = float(painted)
        self._report(
            coverage, package="k_array", name=value, edges_only=False,
            empty="the K array would be its defaults everywhere",
        )
        return result


__all__ = ["CellSurfaceOffset", "GeoPackageSource"]

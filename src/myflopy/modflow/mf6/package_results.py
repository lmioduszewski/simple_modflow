"""Generic result explorer classes behind model.packages."""

from __future__ import annotations

import re
from collections.abc import Iterable
from typing import TYPE_CHECKING

import pandas as pd

from myflopy._logging import get_logger

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase
from myflopy.modflow.mf6.package_budget import (
    budget_value_units,
    build_budget_result_table,
)
from myflopy.modflow.mf6.package_explorer_utils import (
    _filter_normalized_table,
    _normalize_iterable_filter,
    _normalize_term_filter,
)
from myflopy.modflow.mf6.package_plotting import (
    FieldMappable,
    SpatialView,
    _apply_backend,
    _symmetric_color_limit,
    build_cell_input_map_payload,
)
from myflopy.modflow.mf6.package_registry import (
    get_default_package_colorscale,
    get_package_explorer_spec,
    get_package_result_spec,
)
from myflopy.modflow.mf6.package_tables import (
    summarize_input_table,
)
from myflopy.modflow.utils.datatypes.hover import result_hover

logger = get_logger(__name__)


class CellBudgetResultsExplorer(SpatialView):
    """Normalized explorer for one cell-based package result term."""

    def __init__(
        self,
        model: SimulationBase,
        package_name: str,
        budget_text: str,
        value_name: str,
        result_name: str | None = None,
        label: str | None = None,
    ):
        """Bind one cell-based result term: its MF6 ``budget_text`` and output column.

        ``result_name`` is the registry/accessor identity (e.g. ``"q"``, the name
        behind ``results.q``); ``value_name`` is the emitted DataFrame COLUMN
        (e.g. ``"q_gwf"``), which names its reference frame. The two differ for
        the signed exchange results and are equal for everything else.
        ``package_name`` is lowercased.

        ``label`` overrides how :meth:`summary` names this view. It defaults to
        the package-grammar spelling, which is right for ``results.<field>`` but
        wrong for a model-budget term reached through ``model.budget.<term>`` --
        that is not a package result and should not claim to be one.
        """

        self.model = model
        self.package_name = str(package_name).lower()
        self.budget_text = str(budget_text)
        self.value_name = str(value_name)
        self.result_name = str(result_name) if result_name is not None else str(value_name)
        self.label = (
            str(label)
            if label is not None
            else f"{self.package_name}.results.{self.value_name}"
        )

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return the normalized result table for this budget term."""

        frame = build_budget_result_table(
            self.model,
            budget_text=self.budget_text,
            package_name=self.package_name,
            value_name=self.value_name,
        )
        return _filter_normalized_table(frame, per=per, layer=layer, cells=cells)

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of the available result rows."""

        frame = self.get()
        return summarize_input_table(
            frame,
            label=self.label,
            value_columns=[self.value_name],
        )

    def wide(
        self,
        *,
        per: int | Iterable[int] | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
        index: list[str] | tuple[str, ...] = ("layer", "cell"),
        values: str | None = None,
        agg: str = "sum",
    ) -> pd.DataFrame:
        """Pivot this result term to one column per stress period."""

        value_column = self.value_name if values is None else str(values)
        frame = self.get(layer=layer, cells=cells)
        per_values = _normalize_iterable_filter(per)
        if per_values is not None and "per" in frame.columns:
            frame = frame[frame["per"].isin(per_values)]
        if frame.empty:
            return pd.DataFrame(columns=[*index])
        if value_column not in frame.columns:
            raise KeyError(f"Result value column {value_column!r} was not found.")
        missing_index = [column for column in index if column not in frame.columns]
        if missing_index:
            raise KeyError(f"Wide result index columns were not found: {missing_index}")

        wide = frame.pivot_table(
            index=list(index),
            columns="per",
            values=value_column,
            aggfunc=agg,
        )
        wide.columns = [f"per_{int(column)}" for column in wide.columns]
        return wide.reset_index()

    def long(
        self,
        *,
        per: int | Iterable[int] | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
        values: str | None = None,
        agg: str = "sum",
    ) -> pd.Series:
        """Return a long result series indexed by ``kstpkper/layer/cell``."""

        value_column = self.value_name if values is None else str(values)
        frame = self.get(layer=layer, cells=cells)
        per_values = _normalize_iterable_filter(per)
        if per_values is not None and "per" in frame.columns:
            frame = frame[frame["per"].isin(per_values)]
        index_columns = ["kstpkper", "layer", "cell"]
        if frame.empty:
            empty_index = pd.MultiIndex.from_arrays(
                [[] for _ in index_columns],
                names=index_columns,
            )
            return pd.Series([], index=empty_index, dtype=float, name=value_column)
        if value_column not in frame.columns:
            raise KeyError(f"Result value column {value_column!r} was not found.")
        missing_index = [
            column for column in index_columns if column not in frame.columns
        ]
        if missing_index:
            raise KeyError(f"Long result index columns were not found: {missing_index}")

        series = (
            frame.groupby(index_columns, dropna=False)[value_column]
            .agg(agg)
            .sort_index()
        )
        series.name = value_column
        return series

    def stack(self, **kwargs) -> pd.Series:
        """Alias for :meth:`long`."""

        return self.long(**kwargs)

    # NOTE: the series view of this result is the unified grammar's ``plot()``
    # (SpatialView) -- ``results.q.plot(cells=[...])`` replaced plot_timeseries.

    def map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        backend: str = "plotly",
        **kwargs,
    ):
        """Build a choropleth for this cell-based result field."""

        selected = self.get(per=per, layer=layer)
        values, hover = build_cell_input_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            value_column=self.value_name,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        if self.result_name == "q":
            absmax = _symmetric_color_limit(values)
            kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
            kwargs.setdefault("zmax", absmax if absmax > 0 else None)
            kwargs.setdefault("zmid", 0.0)
        result_spec = get_package_result_spec(self.package_name, self.result_name)
        kwargs.setdefault(
            "hover_spec",
            result_hover(
                self.value_name,
                title=f"{self.package_name.upper()} {self.value_name}",
                units=(
                    {self.value_name: budget_value_units(self.model)}
                    if self.result_name == "q"
                    else None
                ),
            ),
        )
        choro = self.model.plot.map(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=(
                colorscale
                or (result_spec.colorscale if result_spec is not None else None)
                or ("RdBu" if self.value_name == "q" else None)
                or get_default_package_colorscale(self.package_name)
                or "earth"
            ),
            **kwargs,
        )
        return _apply_backend(choro, backend)


class StageResultsExplorer(SpatialView):
    """Normalized explorer for cell-mapped stage results such as LAK and SFR."""

    #: value label + series column for the unified grammar
    value_name = "stage"

    def __init__(self, model: SimulationBase, package_name: str, builder):
        """Bind a stage result explorer; ``builder(model)`` yields its normalized table."""

        self.model = model
        self.package_name = str(package_name).lower()
        self._builder = builder

    def _series_default_agg(self) -> str:
        """Collapse cells within a plotted line by mean (stage repeats per connected cell)."""

        return "mean"  # stage repeats per connected cell; summing is meaningless

    def get(
        self,
        *,
        per: int | None = None,
        layer: int | Iterable[int] | None = None,
        cells: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """Return the normalized stage table for the selected rows."""

        frame = self._builder(self.model)
        return _filter_normalized_table(frame, per=per, layer=layer, cells=cells)

    def summary(self) -> pd.DataFrame:
        """Return a compact summary of available stage results."""

        frame = self.get()
        return summarize_input_table(
            frame,
            label=f"{self.package_name}.results.stage",
            value_columns=["stage"],
        )

    def map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "first",
        colorscale: str | None = None,
        backend: str = "plotly",
        **kwargs,
    ):
        """Build a stage choropleth mapped to cells."""

        selected = self.get(per=per, layer=layer)
        values, hover = build_cell_input_map_payload(
            selected,
            ncpl=self.model.vor.ncpl,
            value_column="stage",
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        kwargs.setdefault(
            "hover_spec",
            result_hover(
                "stage",
                title=f"{self.package_name.upper()} stage",
                units={"stage": "ft"},
            ),
        )
        choro = self.model.plot.map(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or "earth",
            **kwargs,
        )
        return _apply_backend(choro, backend)


class CellPackageResultsNamespace(FieldMappable):
    """Namespace for cell-based package results represented by one budget term.

    Simple BC packages expose a single field ``q``; ``results.map()`` maps it and
    ``results.map(field="q")`` is explicit.
    """

    _default_field = "q"

    def _field_names(self):
        """Registry-declared result fields, always including ``"q"`` (the default term)."""

        names = self.fields["field"].tolist()
        return names if "q" in names else ["q", *names]  # .q always resolves

    def __init__(self, model: SimulationBase, package_name: str):
        """Bind a cell-based results namespace to ``model`` for one package."""

        self.model = model
        self.package_name = str(package_name).lower()

    @property
    def fields(self) -> pd.DataFrame:
        """Return registry metadata for supported package result fields."""

        spec = get_package_explorer_spec(self.package_name)
        if spec is None:
            return pd.DataFrame(
                columns=[
                    "field",
                    "budget_text",
                    "value_name",
                    "label",
                    "colorscale",
                    "diverging",
                ]
            )
        rows = [
            {
                "field": result_spec.name,
                "budget_text": result_spec.budget_text,
                "value_name": result_spec.value_name,
                "label": result_spec.label,
                "colorscale": result_spec.colorscale,
                "diverging": result_spec.diverging,
            }
            for result_spec in spec.results.values()
        ]
        return pd.DataFrame(rows)

    def summary(self) -> pd.DataFrame:
        """Return one compact summary row per supported package result field."""

        frames = []
        for field_name in self.fields["field"].tolist():
            frames.append(getattr(self, field_name).summary())
        if not frames:
            return pd.DataFrame(
                columns=[
                    "label",
                    "records",
                    "periods",
                    "layers",
                    "cells",
                    "value_columns",
                ]
            )
        summary = pd.concat(frames, ignore_index=True)
        summary["field_name"] = summary["value_columns"].apply(
            lambda values: values[0] if values else None
        )
        field_metadata = self.fields.rename(columns={"field": "field_name"})
        return summary.merge(field_metadata, on="field_name", how="left")

    def __getattr__(self, result_name: str) -> CellBudgetResultsExplorer:
        """Return a registry-backed named result explorer."""

        result_spec = get_package_result_spec(self.package_name, result_name)
        if result_spec is None:
            raise AttributeError(
                f"{type(self).__name__!s} has no result {result_name!r}"
            )
        return CellBudgetResultsExplorer(
            self.model,
            self.package_name,
            result_spec.budget_text,
            result_spec.value_name,
            result_name=result_name,
        )

    @property
    def q(self) -> CellBudgetResultsExplorer:
        """Return the primary package-exchange result explorer."""

        result_spec = get_package_result_spec(self.package_name, "q")
        if result_spec is None:
            return CellBudgetResultsExplorer(
                self.model, self.package_name, self.package_name.upper(), "q",
                result_name="q",
            )
        return CellBudgetResultsExplorer(
            self.model,
            self.package_name,
            result_spec.budget_text,
            result_spec.value_name,
            result_name="q",
        )


class UzfResultsNamespace(FieldMappable):
    """Namespace for UZF result explorers.

    Fields: ``gwrch`` (groundwater recharge, the default) and ``sat``
    (unsaturated-zone saturation). ``results.map(field="sat")`` or
    ``results.sat.map()``.
    """

    _default_field = "gwrch"

    def _field_names(self):
        """Registry-declared UZF result fields, defaulting to ``["gwrch", "sat"]``."""

        names = self.fields["field"].tolist()
        return names or ["gwrch", "sat"]

    def __init__(self, model: SimulationBase):
        """Bind the UZF results namespace to ``model``."""

        self.model = model

    @property
    def fields(self) -> pd.DataFrame:
        """Return registry metadata for supported UZF result fields."""

        spec = get_package_explorer_spec("uzf")
        if spec is None:
            return pd.DataFrame(
                columns=[
                    "field",
                    "budget_text",
                    "value_name",
                    "label",
                    "colorscale",
                    "diverging",
                ]
            )
        rows = [
            {
                "field": result_spec.name,
                "budget_text": result_spec.budget_text,
                "value_name": result_spec.value_name,
                "label": result_spec.label,
                "colorscale": result_spec.colorscale,
                "diverging": result_spec.diverging,
            }
            for result_spec in spec.results.values()
        ]
        return pd.DataFrame(rows)

    def summary(self) -> pd.DataFrame:
        """Return one compact summary row per supported UZF result field."""

        frames = []
        for field_name in self.fields["field"].tolist():
            frames.append(self._field(field_name).summary())
        if not frames:
            return pd.DataFrame(
                columns=[
                    "label",
                    "records",
                    "periods",
                    "layers",
                    "cells",
                    "value_columns",
                ]
            )
        summary = pd.concat(frames, ignore_index=True)
        summary["field_name"] = summary["value_columns"].apply(
            lambda values: values[0] if values else None
        )
        field_metadata = self.fields.rename(columns={"field": "field_name"})
        return summary.merge(field_metadata, on="field_name", how="left")

    def _field(self, field_name: str) -> CellBudgetResultsExplorer:
        """Return one registry-backed UZF result explorer."""

        result_spec = get_package_result_spec("uzf", field_name)
        if result_spec is None:
            raise AttributeError(
                f"{type(self).__name__!s} has no UZF result field {field_name!r}"
            )
        return CellBudgetResultsExplorer(
            self.model, "uzf", result_spec.budget_text, result_spec.value_name,
            result_name=field_name,
        )

    def __getattr__(self, field_name: str) -> CellBudgetResultsExplorer:
        """Return a registry-backed UZF result explorer."""

        return self._field(field_name)

    @property
    def gwrch(self) -> CellBudgetResultsExplorer:
        """Return groundwater recharge from the UZF package."""

        return self._field("gwrch")

    @property
    def sat(self) -> CellBudgetResultsExplorer:
        """Return normalized unsaturated-zone saturation results."""

        return self._field("sat")


class PackageBudgetTermExplorer:
    """Filtered helper for one package budget term or a small term family."""

    def __init__(
        self,
        namespace,
        *,
        term: str | Iterable[str],
        label: str,
    ):
        """Pin a parent budget namespace to one term (or family) under a display ``label``."""

        self._namespace = namespace
        self.term = term
        self.label = label

    @property
    def types(self) -> list[str]:
        """Return the normalized MF6 LAK term names covered by this helper."""

        return _normalize_term_filter(self.term) or []

    def get(
        self,
        *,
        per: int | Iterable[int] | None = None,
        **filters,
    ) -> pd.DataFrame:
        """Return the filtered package budget-term dataframe."""

        return self._namespace.get(term=self.term, per=per, **filters)

    def summary(
        self,
        *,
        per: int | Iterable[int] | None = None,
        by: list[str] | tuple[str, ...] | None = None,
        **filters,
    ) -> pd.DataFrame:
        """Summarize the filtered package budget terms."""

        return self._namespace.summary(term=self.term, per=per, by=by, **filters)

    def wide(
        self,
        *,
        per: int | Iterable[int] | None = None,
        index: list[str] | tuple[str, ...] = ("per", "lake"),
        values: str = "q",
        **filters,
    ) -> pd.DataFrame:
        """Pivot the filtered package budget terms to a wide dataframe."""

        return self._namespace.wide(
            term=self.term, per=per, index=index, values=values, **filters
        )


def budget_term_attribute(term: str) -> str:
    """Return the Python attribute name that reaches one MF6 budget record.

    ``"SOURCE-SINK MIX"`` -> ``source_sink_mix``, ``"FLOW-JA-FACE"`` ->
    ``flow_ja_face``, ``"STO-SS"`` -> ``sto_ss``, ``"UZF-GWRCH"`` ->
    ``uzf_gwrch``. MF6 separates words with hyphens and, in one case, a space;
    both collapse to the underscore.

    This is the FIRST normalizer in the codebase to run this direction. The
    package-level ``budget.<term>`` namespaces hand-write both spellings as
    literals, and their only normalizer (``_normalize_term_filter``) upcases
    without converting ``_`` back to ``-`` -- so ``budget.get(term="ext_inflow")``
    silently returns an empty frame there. Deriving the mapping in one place is
    what lets :meth:`ModelBudgetNamespace.__getitem__` accept either spelling.
    """

    return re.sub(r"[^0-9a-z]+", "_", str(term).strip().lower()).strip("_")


class ModelBudgetNamespace:
    """Every term in this model's OWN budget file, as a noun.

    ``model.budget.<term>`` returns a :class:`CellBudgetResultsExplorer`, so each
    term answers the full spatial verb set (``get``/``summary``/``plot``/``map``/
    ``xs``/``mosaic``/``animate``) rather than the reduced set the package-level
    ``<pkg>.budget.<term>`` namespaces offer.

    Terms are **discovered at runtime**, not declared. That is a deliberate
    departure from the hand-written package namespaces, because the model
    budget's term set genuinely varies with kind and configuration: GWT's storage
    term is ``STORAGE-AQUEOUS`` and GWE's is ``STORAGE-CELLBLK``, ``DECAY``
    appears only when MST declares decay or sorption, and a GWF model's terms are
    whichever boundary packages it happens to carry. Declaring them would mean
    hand-maintaining a list that is wrong for most models.

    Available on EVERY model kind. The plumbing underneath
    (``_get_budget_reader``, and the record converter) is kind-neutral, so gating
    this to transport would have been artificial -- and ``STO-SS``, ``STO-SY``,
    ``DATA-SAT`` and ``DATA-SPDIS`` have no package accessor at all, making this
    their only route that is not the legacy ``model.bud(...)`` wrapper.

    Not to be confused with its three neighbours:

    * ``model.bud(pkg)`` -- the legacy compatibility wrapper, raw frames.
    * ``model.budget_cumulative`` / ``model.budget_incremental`` -- whole-model
      totals from the LISTING file, not per-cell.
    * ``model.packages.<pkg>.budget.<term>`` -- a genuinely DIFFERENT file, the
      package-output budget, whose node layout is feature-first. Conflating the
      two is what produced the off-by-one in ledger 92.
    """

    def __init__(self, model: SimulationBase):
        """Bind the model-budget term namespace to ``model``."""

        self.model = model

    @property
    def types(self) -> list[str]:
        """Return the MF6 record names present in this model's budget file."""

        reader = self.model._get_budget_reader()
        return [
            str(name).strip()
            for name in reader.get_unique_record_names(decode=True)
        ]

    def _terms(self) -> dict[str, str]:
        """Map each attribute name to the MF6 record name it reaches."""

        return {budget_term_attribute(name): name for name in self.types}

    def __dir__(self) -> list[str]:
        """List the real terms alongside the normal members, for autocomplete."""

        try:
            terms = self._terms()
        except Exception:  # noqa: BLE001 - dir() must never raise
            # 7.3 left this broad on purpose. `_terms` opens the budget file,
            # and that chain was measured to raise types no builtin tuple
            # covers -- flopy's MFDataException, and NotImplementedError out of
            # the base Grid.shape that CellBudgetFile touches unconditionally
            # when a model has no discretization. An autocomplete that raises
            # is worse than one that comes back short.
            logger.debug("no budget terms for __dir__", exc_info=True)
            terms = {}
        return sorted({*super().__dir__(), *terms})

    def __getattr__(self, name: str):
        """Resolve one budget term to its explorer."""

        # `model` would recurse (it is looked up by _terms below) and dunder /
        # private probes must fail fast rather than open the budget file.
        if name.startswith("_") or name == "model":
            raise AttributeError(name)

        terms = self._terms()
        record = terms.get(name)
        if record is None:
            raise AttributeError(
                f"model {self.model.name!r} has no budget term {name!r}. "
                f"Available terms: {sorted(terms)}"
            )
        return CellBudgetResultsExplorer(
            self.model,
            # The MF6 record name, lowercased -- NOT the attribute spelling. For
            # a term that is also a package ("DRN") this keeps the registry
            # lookup working, so the term inherits that package's colorscale and
            # result spec; for the rest it makes the hover read "SOURCE-SINK MIX"
            # the way MF6 spells it, rather than "SOURCE_SINK_MIX".
            package_name=record.lower(),
            budget_text=record,
            value_name="q",
            result_name="q",
            label=f"budget.{name}",
        )

    def __getitem__(self, term: str):
        """Reach a term by EITHER its attribute name or its MF6 record name.

        ``model.budget["SOURCE-SINK MIX"]`` and ``model.budget["source_sink_mix"]``
        are the same object, so a term name copied straight out of ``types`` (or
        out of an MF6 listing file) always works.
        """

        return getattr(self, budget_term_attribute(term))

    def __repr__(self) -> str:
        """Show which terms this model's budget file actually carries."""

        try:
            terms = ", ".join(sorted(self._terms())) or "no terms"
        except Exception:  # noqa: BLE001 - repr must never raise
            # Broad on purpose, same reasoning as __dir__ above. A __repr__ that
            # raises breaks the debugger and the traceback you were reading when
            # you needed it -- so this one reports the failure in its own text.
            logger.debug("no budget terms for __repr__", exc_info=True)
            terms = "budget file unavailable"
        return f"<{type(self).__name__} {getattr(self.model, 'name', '?')!r}: {terms}>"


__all__ = [
    "CellBudgetResultsExplorer",
    "StageResultsExplorer",
    "CellPackageResultsNamespace",
    "UzfResultsNamespace",
    "PackageBudgetTermExplorer",
    "ModelBudgetNamespace",
    "budget_term_attribute",
]

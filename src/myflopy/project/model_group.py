"""Group-oriented helpers for comparing multiple models with one lazy API."""

from __future__ import annotations
from myflopy.viz import mpl_axes

from collections.abc import Mapping, Sequence
import hashlib
import warnings
from pathlib import Path
from typing import Generic, TypeVar

import numpy as np
import pandas as pd
from myflopy.modflow.mf6.package_explorer import (
    SpatialView,
    FieldMappable,
    LeafFieldSugar,
    get_package_input_field_names,
    _blue_white_red_diverging_colorscale,
    _normalize_connection_type_filter,
    build_surface_water_exchange_cell_table,
    build_surface_water_q_map_payload,
    build_cell_package_input_table,
    build_cell_input_map_payload,
    build_budget_result_table,
    build_lak_connection_table,
    build_lak_budget_result_table,
    build_lak_q_map_payload,
    LakStageResultsExplorer,
    SfrStageResultsExplorer,
    build_sfr_budget_result_table,
    build_sfr_q_map_payload,
    build_group_input_compare_map_payload,
    build_uzf_field_input_table,
    get_default_budget_term,
    get_default_group_compare_colorscale,
    get_default_package_colorscale,
    get_package_input_field_spec,
    get_default_package_value_column,
    _symmetric_color_limit,
)
from myflopy.modflow.mf6.package_surface_water import join_lak_stage, join_sfr_stage
from myflopy.modflow.utils.datatypes.hover import (
    cell_input_hover,
    compare_hover,
    lak_hover,
    result_hover,
    sfr_hover,
    surface_water_hover,
)

TResultsNamespace = TypeVar("TResultsNamespace")


def _warn_deprecated(old: str, new: str) -> None:
    """Emit a ``DeprecationWarning`` steering callers from ``old`` to ``new``."""

    warnings.warn(
        f"{old} is deprecated and will be removed in a future release; use {new}.",
        DeprecationWarning,
        stacklevel=3,
    )


def _coerce_kstpkper(model, per: int | None = None, kstpkper: tuple[int, int] | None = None):
    """Normalize stress-period selectors to a concrete ``(kstp, kper)`` tuple."""

    if kstpkper is not None:
        return tuple(int(value) for value in kstpkper)
    if per is not None:
        return tuple(int(value) for value in model.kstpkper[per])
    return None


def _normalize_iterable_filter(values) -> list[int] | None:
    """Normalize optional scalar-or-iterable selectors to integer lists."""

    if values is None:
        return None
    if isinstance(values, Sequence) and not isinstance(values, (str, bytes)):
        return [int(value) for value in values]
    return [int(values)]


def _filter_group_input_table(
    frame: pd.DataFrame,
    *,
    model_name: str | None = None,
    per: int | None = None,
    layer: int | Sequence[int] | None = None,
    cells: Sequence[int] | None = None,
) -> pd.DataFrame:
    """Apply standard grouped package filters to a normalized table."""

    selected = frame.copy()
    if model_name is not None and "model" in selected.columns:
        selected = selected[selected["model"] == str(model_name)]
    if per is not None and "per" in selected.columns:
        selected = selected[selected["per"] == int(per)]
    layer_values = _normalize_iterable_filter(layer)
    if layer_values is not None and "layer" in selected.columns:
        selected = selected[selected["layer"].isin(layer_values)]
    cell_values = _normalize_iterable_filter(cells)
    if cell_values is not None and "cell" in selected.columns:
        selected = selected[selected["cell"].isin(cell_values)]
    return selected.reset_index(drop=True)


def _stable_compare_keys(
    frame: pd.DataFrame,
    value_columns: Sequence[str],
    *,
    carried: Sequence[str] = ("kstpkper",),
) -> list[str]:
    """Return the stable integer-identity columns to merge two models on.

    A cell-budget table carries float geometry columns (``rlen``,
    ``distance_start/mid/end``) that differ between models by sub-unit
    floating-point noise -- e.g. a shared reach with ``rlen`` 16.000493 in one
    model and 16.000000 in another. Including any of those in an exact-match
    inner-join key drops the row entirely, silently emptying the comparison.

    So the join key is the integer identity only (``per``/``layer``/``cell``/
    ``node``/``node2``/``reach``/``package``): every float column, the value
    columns being diffed, ``model``, and time metadata such as ``kstpkper``
    (which varies with a model's time discretization) are excluded and instead
    carried through from the compared model.
    """

    skip = {"model", *value_columns, *carried}
    return [
        column
        for column in frame.columns
        if column not in skip and not pd.api.types.is_float_dtype(frame[column])
    ]


def _reduce_to_period_end(frame: pd.DataFrame) -> pd.DataFrame:
    """Collapse a per-timestep result table to one period-end row per identity.

    Models with the same physics but different ``nstp`` save their period-end
    output at a different ``kstp`` (e.g. kstp=13 vs kstp=9), so two runs cannot
    be aligned on the full ``(kstp, kper)`` tuple. Reducing each model to the
    largest-``kstp`` record within a stress period lets the comparison align on
    the period alone -- and, for a model that saves several steps per period,
    prevents a same-period self-join from exploding into a cartesian product.

    A ``per`` column (zero-based stress period) is derived from ``kstpkper``.
    This is a no-op for a model that already saves a single step per period.
    """

    if frame.empty or "kstpkper" not in frame.columns:
        return frame
    out = frame.copy()
    out["per"] = out["kstpkper"].map(lambda kk: int(kk[1]))
    out["_kstp"] = out["kstpkper"].map(lambda kk: int(kk[0]))
    identity = [
        column
        for column in ("model", "package", "layer", "cell", "node", "node2", "reach")
        if column in out.columns
    ]
    return (
        out.sort_values(["per", "_kstp"])
        .drop_duplicates(subset=[*identity, "per"], keep="last")
        .drop(columns="_kstp")
        .reset_index(drop=True)
    )


def _resolve_group_compare_target(group: "ModelGroup", model_name: str | None) -> str:
    """Resolve which non-reference model to use for a group difference map."""

    if model_name is not None:
        if model_name == group.reference:
            raise ValueError(
                f"Cannot map {model_name!r} against itself -- it is the group's "
                "reference. Choose a non-reference model."
            )
        if model_name not in group.models:
            raise KeyError(f"Model {model_name!r} is not in the group.")
        return str(model_name)

    non_reference = [name for name in group.run_ids if name != group.reference]
    if len(non_reference) == 1:
        return non_reference[0]
    raise ValueError(
        "A model_name is required to map a difference when the group has more "
        f"than one non-reference model. Pass one of: {non_reference}."
    )


def _ensure_group_map_compatible(group: "ModelGroup", model_name: str):
    """Confirm a selected model can be mapped against the reference grid."""

    reference_model = group.models[group.reference]
    selected_model = group.models[model_name]
    reference_signature = _grid_signature_for_model(reference_model)
    selected_signature = _grid_signature_for_model(selected_model)
    if reference_signature and selected_signature and reference_signature != selected_signature:
        raise ValueError(
            f"Model {model_name!r} does not match reference model {group.reference!r}; "
            "package comparison maps require the same grid."
        )
    if reference_model.grid_type != selected_model.grid_type:
        raise ValueError(
            f"Model {model_name!r} does not match reference model {group.reference!r}; "
            "package comparison maps require the same grid type."
        )
    if reference_model.vor.ncpl != selected_model.vor.ncpl:
        raise ValueError(
            f"Model {model_name!r} does not match reference model {group.reference!r}; "
            "package comparison maps require the same ncpl."
        )


def _default_show_layer_elevs(model) -> bool:
    """Return whether group choropleths should include layer elevation hover."""

    return getattr(model.vor, "gdf_topbtm", None) is not None


def _load_group_model(model_or_path, *, crs: str, verbosity_level: int):
    """Normalize one model/group entry into a model-like object.

    Parameters
    ----------
    model_or_path
        Either an existing model object or a path to an MF6 workspace.
    crs, verbosity_level
        Passed to :func:`myflopy.project.load_mf6_run` when a workspace
        path is provided.
    """

    if isinstance(model_or_path, (str, Path)):
        from myflopy.project.run_model import load_mf6_run

        return load_mf6_run(Path(model_or_path), crs=crs, verbosity_level=verbosity_level)
    return model_or_path


def _grid_signature_for_model(model) -> str | None:
    """Return a stable signature describing one model's grid/discretization."""

    if hasattr(model, "grid_signature"):
        signature = model.grid_signature()
        if signature is not None:
            return signature

    workspace = Path(model.workspace)
    grid_type = str(model.grid_type).lower()
    if grid_type in {"dis", "disu", "disv"}:
        path = workspace / f"{model.name}.{grid_type}"
        if path.exists():
            digest = hashlib.sha256()
            with path.open("r", encoding="utf-8", errors="ignore") as handle:
                for line in handle:
                    stripped = line.strip()
                    if not stripped or stripped.startswith("#"):
                        continue
                    digest.update(stripped.encode("utf-8"))
                    digest.update(b"\n")
            return digest.hexdigest()
    return None


def _coerce_models(models, *, crs: str, verbosity_level: int) -> dict[str, object]:
    """Normalize models or workspace paths into an ordered name->model mapping."""

    if isinstance(models, Mapping):
        coerced = {}
        for name, model in models.items():
            coerced[str(name)] = _load_group_model(model, crs=crs, verbosity_level=verbosity_level)
        return coerced

    if isinstance(models, Sequence) and not isinstance(models, (str, bytes)):
        coerced = {}
        for model in models:
            loaded = _load_group_model(model, crs=crs, verbosity_level=verbosity_level)
            coerced[str(loaded.name)] = loaded
        return coerced

    raise TypeError(
        "models must be a mapping of names to models/workspaces or a sequence of model objects/workspaces."
    )


class _GroupSpatialView(SpatialView):
    """:class:`SpatialView` wired for a :class:`ModelGroup`.

    Faceting defaults to one panel per model (``by="model"``) and the model
    axis is the group's members; the host implements ``_spatial_map(per, layer,
    model)`` to render one model's field as a ``Choro``. ``model=None`` means the
    reference model.
    """

    def _spatial_models(self):
        return list(self.group.models)

    def _spatial_reference_model(self):
        return self.group.models[self.group.reference]

    def _spatial_default_facet(self):
        return "model"

    def _spatial_periods(self):
        frame = self.get()
        columns = getattr(frame, "columns", [])
        if "per" in columns and not frame.empty:
            return sorted({int(value) for value in frame["per"].dropna().tolist()})
        return [0]

    def _group_target(self, model):
        """Resolve a model selector to a concrete model name (default reference)."""

        if model is None:
            return self.group.reference
        name = str(model)
        if name not in self.group.models:
            raise KeyError(f"Model {name!r} is not in the group.")
        return name


class GroupHeads(_GroupSpatialView):
    """Heads accessor for :class:`ModelGroup`.

    Mirrors the single-model ``model.hds`` leaf: aligned multi-model tables via
    :meth:`get`/:meth:`compare`, plus the unified grammar -- ``map``/``plot``/
    ``xs`` panels and ``mosaic``/``animate`` composers, faceting over the
    group's members (reference by default).
    """

    def __init__(self, group: "ModelGroup"):
        self.group = group

    # -- unified grammar hooks ----------------------------------------------
    def _spatial_map(self, *, per: int = 0, layer: int = 0, model=None, **kwargs):
        """Render one member's heads choropleth (raw heads, not deltas)."""

        target = self.group.models[self._group_target(model)]
        return target.hds.map(per=int(per), layer=int(layer), **kwargs)

    def _spatial_periods(self) -> list[int]:
        reference = self.group.models[self.group.reference]
        return sorted({int(key[1]) for key in reference.kstpkper})

    def _spatial_value_label(self) -> str:
        return "head"

    def _series_table(self) -> pd.DataFrame:
        # one head per (model, per, layer, cell): reduce to period-end saves
        return _reduce_to_period_end(self.get())

    def _series_value_column(self, frame) -> str:
        return "elev"

    def _series_default_agg(self) -> str:
        return "mean"

    def _sections(
        self,
        model=None,
        *,
        line=None,
        cells: int | list[int] | None = None,
        per: int | None = None,
        layer: int | list[int] = 0,
        **kwargs,
    ):
        """Return member :class:`XSection` objects for the ``xs`` verbs.

        ``model=None`` overlays every member (labeled by model name);
        ``model="F9b"`` sections just that member.
        """

        from myflopy.modflow.utils.datatypes.xsections import XSection

        return {
            name: XSection(
                model=self.group.models[name],
                section_name=name,
                line=line,
                cells=cells,
                per=per,
                layer=layer,
                **kwargs,
            )
            for name in self._resolve_models(model)
        }

    def get(
        self,
        *,
        per: int | None = None,
        kstpkper: tuple[int, int] | None = None,
        layer: int | list[int] | None = None,
        cells: list[int] | None = None,
    ) -> pd.DataFrame:
        """Return aligned heads for all models in the group.

        Parameters
        ----------
        per, kstpkper
            Optional stress-period selector. ``per`` is zero-based.
        layer
            Optional zero-based layer or layers to keep.
        cells
            Optional zero-based cell ids to keep.
        """

        frames = []
        layer_values = None
        if isinstance(layer, list):
            layer_values = [int(value) for value in layer]
        elif layer is not None:
            layer_values = [int(layer)]

        for model_name, model in self.group.models.items():
            frame = model.all_heads.reset_index().copy()
            frame["model"] = model_name
            selected = _coerce_kstpkper(model, per=per, kstpkper=kstpkper)
            if selected is not None:
                frame = frame[frame["kstpkper"] == selected]
            if layer_values is not None:
                frame = frame[frame["layer"].isin(layer_values)]
            if cells is not None:
                frame = frame[frame["cell"].isin([int(value) for value in cells])]
            frames.append(frame)

        if not frames:
            return pd.DataFrame(columns=["kstpkper", "layer", "cell", "elev", "model"])

        combined = pd.concat(frames, ignore_index=True)
        return combined[["model", "kstpkper", "layer", "cell", "elev"]]

    def compare(
        self,
        *,
        per: int | None = None,
        kstpkper: tuple[int, int] | None = None,
        layer: int | list[int] | None = None,
        cells: list[int] | None = None,
    ) -> pd.DataFrame:
        """Compare heads for all models against the group's reference model."""

        data = self.get(per=per, kstpkper=kstpkper, layer=layer, cells=cells)
        columns = [
            "model", "reference_model", "kstpkper", "per",
            "layer", "cell", "elev", "reference_elev", "diff",
        ]
        if data.empty:
            return pd.DataFrame(columns=columns)

        # Align on the stress PERIOD, not the full ``(kstp, kper)`` tuple: models
        # with the same physics but different time discretization save the
        # period-end head at a different ``kstp`` (e.g. kstp=13 vs kstp=9), so an
        # exact tuple match finds nothing. Reduce each model to its period-end
        # head and merge periods.
        data = _reduce_to_period_end(data)

        reference = self.group.reference
        ref = (
            data[data["model"] == reference]
            .rename(columns={"elev": "reference_elev"})
            .drop(columns=["model", "kstpkper"])
        )
        comp = data[data["model"] != reference].merge(
            ref,
            on=["per", "layer", "cell"],
            how="inner",
        )
        comp["reference_model"] = reference
        comp["diff"] = comp["elev"].astype(float) - comp["reference_elev"].astype(float)
        return comp[columns]

    def compare_map(
        self,
        *,
        model_name: str | None = None,
        per: int = 0,
        layer: int = 0,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        colorscale: str | None = None,
        **kwargs,
    ):
        """Diverging choropleth of head differences (model - reference).

        Colors every Voronoi cell by ``diff`` (Δhead) for one stress period and
        layer, using the same diff-map payload + rendering as the package diff
        maps. ``model_name`` selects the compared model (inferred when the group
        has a single non-reference model). ``per`` / ``layer`` are zero-based.
        """

        target_name = _resolve_group_compare_target(self.group, model_name)
        _ensure_group_map_compatible(self.group, target_name)
        reference_model = self.group.models[self.group.reference]
        selected = self.compare(per=per, layer=layer)
        if not selected.empty:
            selected = selected[selected["model"] == target_name]
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(reference_model))
        values, hover, absmax = build_group_input_compare_map_payload(
            selected,
            ncpl=reference_model.vor.ncpl,
            value_column="elev",
            diff_column="diff",
            model_name=target_name,
            reference_model=self.group.reference,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg="first",  # one head per cell per (per, layer) -- do not sum
        )
        kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
        kwargs.setdefault("zmax", absmax if absmax > 0 else None)
        kwargs.setdefault("zmid", 0.0)
        kwargs.setdefault(
            "hover_spec",
            compare_hover(
                "elev",
                "diff",
                title="Δ head vs reference",
                units={"elev": "ft", "reference_elev": "ft", "diff": "ft"},
                labels={"diff": "Δ head"},
            ),
        )
        return reference_model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or get_default_group_compare_colorscale(),
            **kwargs,
        )


class GroupBudget:
    """Budget accessor for :class:`ModelGroup`."""

    def __init__(self, group: "ModelGroup", package: str | None = None):
        self.group = group
        self.package = None if package is None else str(package).lower()

    @staticmethod
    def _normalize_nodes(frame: pd.DataFrame) -> pd.DataFrame:
        """Normalize budget node columns to zero-based indexing when needed."""

        normalized = frame.copy()
        for column in ("node", "node2"):
            if column in normalized.columns:
                numeric = pd.to_numeric(normalized[column], errors="coerce")
                mask = numeric.notna()
                if mask.any():
                    ints = numeric.loc[mask].astype(int)
                    if int(ints.min()) >= 1:
                        normalized[column] = numeric
                        normalized.loc[mask, column] = (ints - 1).astype(float)
                    else:
                        normalized[column] = numeric
                        normalized.loc[mask, column] = ints.astype(float)
        return normalized

    def get(
        self,
        *,
        package: str | None = None,
    ) -> pd.DataFrame:
        """Return aligned budget data for all models in the group.

        Parameters
        ----------
        package
            Optional package filter such as ``"rch"`` or ``"drn"``. If omitted,
            the package passed to :meth:`ModelGroup.bud` is used.
        """

        package_name = self.package if package is None else str(package).lower()
        if package_name is None:
            raise ValueError("A budget package name is required, for example group.bud('rch').get().")

        frames = []
        for model_name, model in self.group.models.items():
            budget = model.bud(package_name)
            frame = budget.df.reset_index().copy()
            renamed = {}
            if "level_0" in frame.columns:
                renamed["level_0"] = "node"
            if "level_1" in frame.columns:
                renamed["level_1"] = "kstpkper"
            if renamed:
                frame = frame.rename(columns=renamed)
            frame = self._normalize_nodes(frame)
            frame["model"] = model_name
            frames.append(frame)

        if not frames:
            return pd.DataFrame()
        return pd.concat(frames, ignore_index=True)

    def compare(
        self,
        *,
        package: str | None = None,
        value_column: str = "q",
    ) -> pd.DataFrame:
        """Compare budget values for all models against the reference model.

        Parameters
        ----------
        package
            Optional package override.
        value_column
            Budget value column to compare. Defaults to ``"q"``.
        """

        data = self.get(package=package)
        if data.empty:
            return pd.DataFrame()

        if value_column not in data.columns:
            raise KeyError(f"Budget value column {value_column!r} was not found.")

        reference = self.group.reference
        key_columns = [column for column in data.columns if column not in {"model", value_column}]
        ref = (
            data[data["model"] == reference]
            .rename(columns={value_column: f"reference_{value_column}"})
            .drop(columns=["model"])
        )
        comp = data[data["model"] != reference].merge(ref, on=key_columns, how="inner")
        comp["reference_model"] = reference
        comp["diff"] = comp[value_column].astype(float) - comp[f"reference_{value_column}"].astype(float)
        ordered = ["model", "reference_model", *key_columns, value_column, f"reference_{value_column}", "diff"]
        return comp[ordered]


class GroupPackageInputField(_GroupSpatialView):
    """One input field of a grouped package, pinned for the unified grammar.

    ``group.packages.ghb.inputs.cond`` -- same verbs as the parent accessor
    but every panel/series draws this field.
    """

    def __init__(self, parent: "GroupPackageInputs", field_name: str):
        self._parent = parent
        self.group = parent.group
        self.package_name = parent.package_name
        self.field_name = str(field_name).lower()

    def get(self, **kwargs) -> pd.DataFrame:
        """Return the aligned input rows narrowed to this field."""

        frame = self._parent.get(**kwargs)
        keep = [
            column
            for column in ("model", "package", "per", "layer", "cell", self.field_name)
            if column in getattr(frame, "columns", [])
        ]
        return frame[keep] if keep else frame

    def summary(self, **kwargs) -> pd.DataFrame:
        """Return the parent's per-model coverage summary."""

        return self._parent.summary(**kwargs)

    def _spatial_map(self, **kwargs):
        kwargs.setdefault("value_column", self.field_name)
        return self._parent._spatial_map(**kwargs)

    def _spatial_value_label(self) -> str:
        return self.field_name

    def _series_value_column(self, frame) -> str:
        return self.field_name


class GroupPackageInputs(LeafFieldSugar, _GroupSpatialView):
    """Grouped input accessor for simple cell-based MF6 stress-period packages.

    Inherits the unified grammar (``map``/``plot``/``mosaic``/``animate``;
    ``mosaic()`` defaults to one panel per model). Registry-backed fields are
    first-class nodes (``inputs.cond``) and every verb takes ``field=`` as
    sugar over them; without ``field=`` the package default field is drawn.
    """

    def __init__(self, group: "ModelGroup", package_name: str):
        self.group = group
        self.package_name = str(package_name).lower()

    def _field_names(self) -> list[str]:
        return get_package_input_field_names(self.package_name)

    def _field_node(self, name: str) -> GroupPackageInputField:
        return GroupPackageInputField(self, name)

    def get(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """Return aligned stress-period input data for all models."""

        frames = []
        for current_model_name, model in self.group.models.items():
            frame = build_cell_package_input_table(model, self.package_name)
            frame["model"] = current_model_name
            frames.append(frame)

        if not frames:
            return pd.DataFrame()
        combined = pd.concat(frames, ignore_index=True)
        return _filter_group_input_table(
            combined,
            model_name=model_name,
            per=per,
            layer=layer,
            cells=cells,
        )

    def summary(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """Return a per-model summary of this package's input coverage.

        Mirrors the single-model ``model.packages.<pkg>.inputs.summary()`` verb:
        one row per model with the number of applied cells and stress periods.
        """

        data = self.get(model_name=model_name, per=per, layer=layer)
        columns = ["model", "package", "cells", "periods"]
        if data.empty:
            return pd.DataFrame(columns=columns)
        rows = []
        for current_model_name, sub in data.groupby("model"):
            rows.append(
                {
                    "model": current_model_name,
                    "package": self.package_name,
                    "cells": int(sub["cell"].nunique()) if "cell" in sub.columns else 0,
                    "periods": int(sub["per"].nunique()) if "per" in sub.columns else 0,
                }
            )
        return pd.DataFrame(rows, columns=columns)

    def compare(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """Compare input values for all models against the reference model."""

        data = self.get(per=per, layer=layer, cells=cells)
        if data.empty:
            return pd.DataFrame()

        reference = self.group.reference
        value_columns = [
            column
            for column in data.columns
            if column not in {"model", "per", "layer", "cell"}
            and pd.api.types.is_numeric_dtype(data[column])
        ]
        key_columns = [column for column in data.columns if column not in {"model", *value_columns}]
        ref = data[data["model"] == reference].drop(columns=["model"]).copy()
        rename_map = {column: f"reference_{column}" for column in value_columns}
        ref = ref.rename(columns=rename_map)
        comp = data[data["model"] != reference].merge(ref, on=key_columns, how="inner")
        comp["reference_model"] = reference
        for column in value_columns:
            comp[f"{column}_diff"] = comp[column].astype(float) - comp[f"reference_{column}"].astype(float)
        ordered = ["model", "reference_model", *key_columns]
        for column in value_columns:
            ordered.extend([column, f"reference_{column}", f"{column}_diff"])
        result = comp[ordered]
        return _filter_group_input_table(
            result,
            model_name=model_name,
            per=per,
            layer=layer,
            cells=cells,
        )

    def _spatial_map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        model=None,
        value_column: str | None = None,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Render one model's package-input field as a raw-value ``Choro``."""

        target_name = self._group_target(model)
        target_model = self.group.models[target_name]
        selected = self.get(model_name=target_name, per=per, layer=layer)
        fallback = get_default_package_value_column(self.package_name)
        source = selected if not selected.empty else self.get(model_name=target_name)
        numeric_columns = [
            column
            for column in source.columns
            if column not in {"model", "package", "per", "layer", "cell"}
            and pd.api.types.is_numeric_dtype(source[column])
        ]
        chosen_value_column = (
            value_column
            or (fallback if fallback in source.columns else None)
            or (numeric_columns[0] if numeric_columns else None)
        )
        if chosen_value_column is None:
            raise ValueError(f"Could not infer a numeric value column for package {self.package_name!r}.")
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(target_model))
        values, hover = build_cell_input_map_payload(
            selected,
            ncpl=target_model.vor.ncpl,
            value_column=chosen_value_column,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        kwargs.setdefault("hover_spec", cell_input_hover(chosen_value_column))
        return target_model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or get_default_package_colorscale(self.package_name),
            **kwargs,
        )

    def compare_map(
        self,
        *,
        model_name: str | None = None,
        per: int = 0,
        layer: int = 0,
        value_column: str | None = None,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Build a diff choropleth against the group's reference model.

        Parameters
        ----------
        model_name
            Non-reference model to compare against the reference. When the
            group contains exactly one non-reference model, it is inferred.
        per, layer
            Zero-based stress period and layer to map.
        value_column
            Numeric input field whose difference should be mapped.
        multiplier, fill_value, agg
            Passed through to the shared diff-map payload builder.
        colorscale
            Optional diverging colorscale override.
        kwargs
            Forwarded to ``reference_model.cor(...)``.
        """

        target_name = _resolve_group_compare_target(self.group, model_name)
        _ensure_group_map_compatible(self.group, target_name)
        reference_model = self.group.models[self.group.reference]
        selected = self.compare(model_name=target_name, per=per, layer=layer)
        fallback = get_default_package_value_column(self.package_name)
        source = selected if not selected.empty else self.compare(model_name=target_name)
        numeric_columns = [
            column
            for column in source.columns
            if column not in {"per", "layer", "cell", "model", "reference_model", "package"}
            and pd.api.types.is_numeric_dtype(source[column])
            and not column.startswith("reference_")
            and not column.endswith("_diff")
        ]
        chosen_value_column = (
            value_column
            or (fallback if fallback in source.columns else None)
            or (numeric_columns[0] if numeric_columns else None)
        )
        if chosen_value_column is None:
            raise ValueError(f"Could not infer a numeric value column for package {self.package_name!r}.")
        diff_column = f"{chosen_value_column}_diff"
        if diff_column not in source.columns:
            raise KeyError(f"Difference column {diff_column!r} was not found.")
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(reference_model))
        values, hover, absmax = build_group_input_compare_map_payload(
            selected,
            ncpl=reference_model.vor.ncpl,
            value_column=chosen_value_column,
            diff_column=diff_column,
            model_name=target_name,
            reference_model=self.group.reference,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
        kwargs.setdefault("zmax", absmax if absmax > 0 else None)
        kwargs.setdefault("hover_spec", compare_hover(chosen_value_column, diff_column))
        return reference_model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or get_default_group_compare_colorscale(),
            **kwargs,
        )


class GroupCellPackageResults(_GroupSpatialView):
    """Grouped accessor for cell-based package result tables.

    Inherits the unified ``map`` / ``mosaic`` / ``animate`` grammar (``mosaic()``
    defaults to one panel per model).
    """

    def __init__(self, group: "ModelGroup", package_name: str, *, budget_text: str, value_name: str):
        self.group = group
        self.package_name = str(package_name).lower()
        self.budget_text = str(budget_text)
        self.value_name = str(value_name)

    def get(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """Return aligned package-result rows for all models in the group."""

        frames = []
        for current_model_name, model in self.group.models.items():
            frame = build_budget_result_table(
                model,
                budget_text=self.budget_text,
                package_name=self.package_name,
                value_name=self.value_name,
            )
            frame["model"] = current_model_name
            frames.append(frame)

        if not frames:
            return pd.DataFrame()
        combined = pd.concat(frames, ignore_index=True)
        return _filter_group_input_table(
            combined,
            model_name=model_name,
            per=per,
            layer=layer,
            cells=cells,
        )

    def compare(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """Compare grouped package results against the reference model."""

        data = self.get(per=per, layer=layer, cells=cells)
        if data.empty:
            return pd.DataFrame()
        data = _reduce_to_period_end(data)

        reference = self.group.reference
        value_columns = [self.value_name]
        key_columns = _stable_compare_keys(data, value_columns)
        ref = data[data["model"] == reference][key_columns + value_columns].copy()
        rename_map = {column: f"reference_{column}" for column in value_columns}
        ref = ref.rename(columns=rename_map)
        comp = data[data["model"] != reference].merge(ref, on=key_columns, how="inner")
        comp["reference_model"] = reference
        comp[f"{self.value_name}_diff"] = (
            comp[self.value_name].astype(float) - comp[f"reference_{self.value_name}"].astype(float)
        )
        ordered = [
            "model",
            "reference_model",
            *key_columns,
            self.value_name,
            f"reference_{self.value_name}",
            f"{self.value_name}_diff",
        ]
        result = comp[ordered]
        return _filter_group_input_table(
            result,
            model_name=model_name,
            per=per,
            layer=layer,
            cells=cells,
        )

    def _spatial_map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        model=None,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Render one model's package result field as a ``Choro``."""

        target_name = self._group_target(model)
        target_model = self.group.models[target_name]
        selected = self.get(model_name=target_name, per=per, layer=layer)
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(target_model))
        values, hover = build_cell_input_map_payload(
            selected,
            ncpl=target_model.vor.ncpl,
            value_column=self.value_name,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        if self.value_name == "q":
            absmax = _symmetric_color_limit(values)
            kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
            kwargs.setdefault("zmax", absmax if absmax > 0 else None)
            kwargs.setdefault("zmid", 0.0)
        kwargs.setdefault(
            "hover_spec",
            result_hover(
                self.value_name,
                title=f"{self.package_name.upper()} {self.value_name}",
                units={"q": "ft³/d"} if self.value_name == "q" else None,
            ),
        )
        return target_model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=(
                colorscale
                or ("RdBu" if self.value_name == "q" else None)
                or get_default_package_colorscale(self.package_name)
                or "earth"
            ),
            **kwargs,
        )

    def compare_map(
        self,
        *,
        model_name: str | None = None,
        per: int = 0,
        layer: int = 0,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Build a diff choropleth for one grouped package result field."""

        target_name = _resolve_group_compare_target(self.group, model_name)
        _ensure_group_map_compatible(self.group, target_name)
        reference_model = self.group.models[self.group.reference]
        selected = self.compare(model_name=target_name, per=per, layer=layer)
        diff_column = f"{self.value_name}_diff"
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(reference_model))
        values, hover, absmax = build_group_input_compare_map_payload(
            selected,
            ncpl=reference_model.vor.ncpl,
            value_column=self.value_name,
            diff_column=diff_column,
            model_name=target_name,
            reference_model=self.group.reference,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
        kwargs.setdefault("zmax", absmax if absmax > 0 else None)
        kwargs.setdefault("zmid", 0.0)
        kwargs.setdefault(
            "hover_spec",
            compare_hover(
                self.value_name,
                diff_column,
                title=f"Δ {self.package_name.upper()} {self.value_name} vs reference",
            ),
        )
        return reference_model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or get_default_group_compare_colorscale(),
            **kwargs,
        )

    # NOTE: the series view is the unified grammar's ``plot()`` (SpatialView) --
    # one line per model (and per cell with ``cells=[...]``).


class GroupSfrBudgetResults(GroupCellPackageResults):
    """Grouped SFR exchange accessor with reach-length-normalized maps."""

    def __init__(self, group: "ModelGroup", *, budget_text: str = "SFR", value_name: str = "q"):
        super().__init__(group, "sfr", budget_text=budget_text, value_name=value_name)

    def get(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """Return aligned SFR exchange rows, including ``q_per_length``."""

        frames = []
        for current_model_name, model in self.group.models.items():
            frame = build_sfr_budget_result_table(
                model,
                budget_text=self.budget_text,
                value_name=self.value_name,
            )
            frame["model"] = current_model_name
            frames.append(frame)

        if not frames:
            return pd.DataFrame()
        combined = pd.concat(frames, ignore_index=True)
        return _filter_group_input_table(
            combined,
            model_name=model_name,
            per=per,
            layer=layer,
            cells=cells,
        )

    def compare(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """Compare grouped SFR exchange rows against the reference model."""

        data = self.get(per=per, layer=layer, cells=cells)
        if data.empty:
            return pd.DataFrame()
        data = _reduce_to_period_end(data)

        reference = self.group.reference
        value_columns = ["q", "q_per_length"] if "q_per_length" in data.columns else ["q"]
        key_columns = _stable_compare_keys(data, value_columns)
        ref = data[data["model"] == reference][key_columns + value_columns].copy()
        rename_map = {column: f"reference_{column}" for column in value_columns}
        ref = ref.rename(columns=rename_map)
        comp = data[data["model"] != reference].merge(ref, on=key_columns, how="inner")
        comp["reference_model"] = reference
        comp["q_diff"] = comp["q"].astype(float) - comp["reference_q"].astype(float)
        if "q_per_length" in value_columns:
            comp["q_per_length_diff"] = (
                comp["q_per_length"].astype(float) - comp["reference_q_per_length"].astype(float)
            )
        ordered = [
            "model",
            "reference_model",
            *key_columns,
            "q",
            "reference_q",
            "q_diff",
        ]
        if "q_per_length" in value_columns:
            ordered.extend(["q_per_length", "reference_q_per_length", "q_per_length_diff"])
        result = comp[ordered]
        return _filter_group_input_table(
            result,
            model_name=model_name,
            per=per,
            layer=layer,
            cells=cells,
        )

    def _spatial_map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        model=None,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Render one model's SFR exchange (reach-length-normalized) as a ``Choro``."""

        del agg
        target_name = self._group_target(model)
        target_model = self.group.models[target_name]
        selected = join_sfr_stage(
            target_model, self.get(model_name=target_name, per=per, layer=layer), per=per
        )
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(target_model))
        values, hover = build_sfr_q_map_payload(
            selected,
            ncpl=target_model.vor.ncpl,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
        )
        absmax = _symmetric_color_limit(values)
        kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
        kwargs.setdefault("zmax", absmax if absmax > 0 else None)
        kwargs.setdefault("zmid", 0.0)
        kwargs.setdefault("hover_spec", sfr_hover())
        return target_model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or _blue_white_red_diverging_colorscale(),
            **kwargs,
        )

    def compare_map(
        self,
        *,
        model_name: str | None = None,
        per: int = 0,
        layer: int = 0,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Build a grouped SFR diff map using normalized exchange per unit length."""

        del agg
        target_name = _resolve_group_compare_target(self.group, model_name)
        _ensure_group_map_compatible(self.group, target_name)
        reference_model = self.group.models[self.group.reference]
        data = self.get(per=per, layer=layer)

        def _normalized_by_cell(frame: pd.DataFrame, current_model_name: str) -> pd.DataFrame:
            selected = frame[frame["model"] == current_model_name].copy()
            if selected.empty:
                return pd.DataFrame(columns=["cell", "q_per_length"])
            selected["q"] = pd.to_numeric(selected["q"], errors="coerce")
            selected["rlen"] = pd.to_numeric(selected["rlen"], errors="coerce")
            grouped = selected.groupby("cell", as_index=False).agg({"q": "sum", "rlen": "sum"})
            grouped["q_per_length"] = np.where(grouped["rlen"] > 0.0, grouped["q"] / grouped["rlen"], np.nan)
            return grouped[["cell", "q_per_length"]]

        reference = _normalized_by_cell(data, self.group.reference).rename(
            columns={"q_per_length": "reference_q_per_length"}
        )
        target = _normalized_by_cell(data, target_name)
        comparison = target.merge(reference, on="cell", how="inner")
        comparison["model"] = target_name
        comparison["reference_model"] = self.group.reference
        comparison["per"] = int(per)
        comparison["layer"] = int(layer)
        comparison["q_per_length_diff"] = (
            comparison["q_per_length"].astype(float) - comparison["reference_q_per_length"].astype(float)
        )
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(reference_model))
        values, hover, absmax = build_group_input_compare_map_payload(
            comparison,
            ncpl=reference_model.vor.ncpl,
            value_column="q_per_length",
            diff_column="q_per_length_diff",
            model_name=target_name,
            reference_model=self.group.reference,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg="sum",
        )
        kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
        kwargs.setdefault("zmax", absmax if absmax > 0 else None)
        kwargs.setdefault("zmid", 0.0)
        kwargs.setdefault(
            "hover_spec",
            compare_hover(
                "q_per_length",
                "q_per_length_diff",
                title="Δ stream exchange vs reference",
                units={
                    "q_per_length": "ft²/d",
                    "reference_q_per_length": "ft²/d",
                    "q_per_length_diff": "ft²/d",
                },
                labels={"q_per_length_diff": "Δ q / length"},
            ),
        )
        return reference_model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or get_default_group_compare_colorscale(),
            **kwargs,
        )

class GroupLakBudgetResults(GroupCellPackageResults):
    """Grouped LAK exchange accessor with area-normalized maps."""

    def __init__(self, group: "ModelGroup", *, budget_text: str = "GWF", value_name: str = "q"):
        super().__init__(group, "lak", budget_text=budget_text, value_name=value_name)

    def get(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
        connection_type: str | Sequence[str] | None = None,
    ) -> pd.DataFrame:
        """Return aligned LAK exchange rows, including ``q_per_area``."""

        frames = []
        for current_model_name, model in self.group.models.items():
            frame = build_lak_budget_result_table(
                model,
                budget_text=self.budget_text,
                value_name=self.value_name,
            )
            frame["model"] = current_model_name
            frames.append(frame)

        if not frames:
            return pd.DataFrame()
        combined = pd.concat(frames, ignore_index=True)
        combined = _filter_group_input_table(
            combined,
            model_name=model_name,
            per=per,
            layer=layer,
            cells=cells,
        )
        connection_types = _normalize_connection_type_filter(connection_type)
        if connection_types is not None and "claktype" in combined.columns:
            combined = combined.loc[combined["claktype"].astype("string").str.upper().isin(connection_types)].copy()
        return combined.reset_index(drop=True)

    def compare(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
        connection_type: str | Sequence[str] | None = None,
    ) -> pd.DataFrame:
        """Compare grouped LAK exchange rows against the reference model."""

        data = self.get(per=per, layer=layer, cells=cells, connection_type=connection_type)
        if data.empty:
            return pd.DataFrame()

        reference = self.group.reference
        value_columns = ["q", "q_per_area"] if "q_per_area" in data.columns else ["q"]
        key_columns = [column for column in data.columns if column not in {"model", *value_columns}]
        ref = data[data["model"] == reference].drop(columns=["model"]).copy()
        rename_map = {column: f"reference_{column}" for column in value_columns}
        ref = ref.rename(columns=rename_map)
        comp = data[data["model"] != reference].merge(ref, on=key_columns, how="inner")
        comp["reference_model"] = reference
        comp["q_diff"] = comp["q"].astype(float) - comp["reference_q"].astype(float)
        if "q_per_area" in value_columns:
            comp["q_per_area_diff"] = comp["q_per_area"].astype(float) - comp["reference_q_per_area"].astype(float)
        ordered = ["model", "reference_model", *key_columns, "q", "reference_q", "q_diff"]
        if "q_per_area" in value_columns:
            ordered.extend(["q_per_area", "reference_q_per_area", "q_per_area_diff"])
        result = comp[ordered]
        return _filter_group_input_table(
            result,
            model_name=model_name,
            per=per,
            layer=layer,
            cells=cells,
        )

    def _spatial_map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        model=None,
        connection_type: str | Sequence[str] | None = None,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Render one model's LAK exchange (area-normalized) as a ``Choro``."""

        del agg
        target_name = self._group_target(model)
        target_model = self.group.models[target_name]
        selected = join_lak_stage(
            target_model,
            self.get(model_name=target_name, per=per, layer=layer, connection_type=connection_type),
            per=per,
        )
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(target_model))
        values, hover = build_lak_q_map_payload(
            selected,
            ncpl=target_model.vor.ncpl,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
        )
        absmax = _symmetric_color_limit(values)
        kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
        kwargs.setdefault("zmax", absmax if absmax > 0 else None)
        kwargs.setdefault("zmid", 0.0)
        kwargs.setdefault("hover_spec", lak_hover())
        return target_model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            # match the SFR convention: gaining (negative q) blue, losing red
            colorscale=colorscale or _blue_white_red_diverging_colorscale(),
            **kwargs,
        )

    def compare_map(
        self,
        *,
        model_name: str | None = None,
        per: int = 0,
        layer: int = 0,
        connection_type: str | Sequence[str] | None = None,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Build a grouped LAK diff map using normalized exchange per area."""

        del agg
        target_name = _resolve_group_compare_target(self.group, model_name)
        _ensure_group_map_compatible(self.group, target_name)
        reference_model = self.group.models[self.group.reference]
        data = self.get(per=per, layer=layer, connection_type=connection_type)

        def _normalized_by_cell(frame: pd.DataFrame, current_model_name: str) -> pd.DataFrame:
            selected = frame[frame["model"] == current_model_name].copy()
            if selected.empty:
                return pd.DataFrame(columns=["cell", "q_per_area"])
            selected["q"] = pd.to_numeric(selected["q"], errors="coerce")
            selected["flow_area"] = pd.to_numeric(selected["flow_area"], errors="coerce")
            grouped = selected.groupby("cell", as_index=False).agg({"q": "sum", "flow_area": "sum"})
            grouped["q_per_area"] = np.where(grouped["flow_area"] > 0.0, grouped["q"] / grouped["flow_area"], np.nan)
            return grouped[["cell", "q_per_area"]]

        reference = _normalized_by_cell(data, self.group.reference).rename(
            columns={"q_per_area": "reference_q_per_area"}
        )
        target = _normalized_by_cell(data, target_name)
        comparison = target.merge(reference, on="cell", how="inner")
        comparison["model"] = target_name
        comparison["reference_model"] = self.group.reference
        comparison["per"] = int(per)
        comparison["layer"] = int(layer)
        comparison["q_per_area_diff"] = (
            comparison["q_per_area"].astype(float) - comparison["reference_q_per_area"].astype(float)
        )
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(reference_model))
        values, hover, absmax = build_group_input_compare_map_payload(
            comparison,
            ncpl=reference_model.vor.ncpl,
            value_column="q_per_area",
            diff_column="q_per_area_diff",
            model_name=target_name,
            reference_model=self.group.reference,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg="sum",
        )
        kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
        kwargs.setdefault("zmax", absmax if absmax > 0 else None)
        kwargs.setdefault("zmid", 0.0)
        kwargs.setdefault(
            "hover_spec",
            compare_hover(
                "q_per_area",
                "q_per_area_diff",
                title="Δ lake exchange vs reference",
                units={
                    "q_per_area": "ft/d",
                    "reference_q_per_area": "ft/d",
                    "q_per_area_diff": "ft/d",
                },
                labels={"q_per_area_diff": "Δ q / area"},
            ),
        )
        return reference_model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or get_default_group_compare_colorscale(),
            **kwargs,
        )


class GroupUzfFieldAccessor(_GroupSpatialView):
    """Grouped accessor for one UZF perioddata field such as ``finf``."""

    def __init__(self, group: "ModelGroup", field_name: str):
        self.group = group
        self.field_name = str(field_name)

    def get(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """Return one aligned UZF perioddata field for all models."""

        frames = []
        for current_model_name, model in self.group.models.items():
            frame = build_uzf_field_input_table(model, self.field_name)
            frame["model"] = current_model_name
            frames.append(frame)

        if not frames:
            return pd.DataFrame(columns=["model", "per", "ifno", "layer", "cell", self.field_name])

        combined = pd.concat(frames, ignore_index=True)
        combined["ifno"] = combined["ifno"].astype(int)
        combined["layer"] = combined["layer"].astype(int)
        combined["cell"] = combined["cell"].astype(int)
        return _filter_group_input_table(
            combined,
            model_name=model_name,
            per=per,
            layer=layer,
            cells=cells,
        )

    def compare(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """Compare the grouped UZF field against the reference model."""

        data = self.get(per=per, layer=layer, cells=cells)
        if data.empty:
            return pd.DataFrame()

        reference = self.group.reference
        ref = (
            data[data["model"] == reference]
            .rename(columns={self.field_name: f"reference_{self.field_name}"})
            .drop(columns=["model"])
        )
        comp = data[data["model"] != reference].merge(
            ref,
            on=["per", "ifno", "layer", "cell"],
            how="inner",
        )
        comp["reference_model"] = reference
        comp[f"{self.field_name}_diff"] = (
            comp[self.field_name].astype(float) - comp[f"reference_{self.field_name}"].astype(float)
        )
        result = comp[
            [
                "model",
                "reference_model",
                "per",
                "ifno",
                "layer",
                "cell",
                self.field_name,
                f"reference_{self.field_name}",
                f"{self.field_name}_diff",
            ]
        ]
        return _filter_group_input_table(
            result,
            model_name=model_name,
            per=per,
            layer=layer,
            cells=cells,
        )

    def _spatial_map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        model=None,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Render one model's UZF field as a raw-value ``Choro``."""

        target_name = self._group_target(model)
        target_model = self.group.models[target_name]
        selected = self.get(model_name=target_name, per=per, layer=layer)
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(target_model))
        values, hover = build_cell_input_map_payload(
            selected,
            ncpl=target_model.vor.ncpl,
            value_column=self.field_name,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        kwargs.setdefault("hover_spec", cell_input_hover(self.field_name))
        return target_model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or get_default_package_colorscale(f"uzf_{self.field_name}") or "earth",
            **kwargs,
        )

    def compare_map(
        self,
        *,
        model_name: str | None = None,
        per: int = 0,
        layer: int = 0,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale: str | None = None,
        **kwargs,
    ):
        """Build a diff choropleth for one grouped UZF field."""

        target_name = _resolve_group_compare_target(self.group, model_name)
        _ensure_group_map_compatible(self.group, target_name)
        reference_model = self.group.models[self.group.reference]
        selected = self.compare(model_name=target_name, per=per, layer=layer)
        diff_column = f"{self.field_name}_diff"
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(reference_model))
        values, hover, absmax = build_group_input_compare_map_payload(
            selected,
            ncpl=reference_model.vor.ncpl,
            value_column=self.field_name,
            diff_column=diff_column,
            model_name=target_name,
            reference_model=self.group.reference,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
        kwargs.setdefault("zmax", absmax if absmax > 0 else None)
        kwargs.setdefault("hover_spec", compare_hover(self.field_name, diff_column))
        return reference_model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or get_default_group_compare_colorscale(),
            **kwargs,
        )


class GroupUzfInputs:
    """Namespace for grouped UZF input fields."""

    def __init__(self, group: "ModelGroup"):
        self.group = group

    def _field(self, field_name: str) -> GroupUzfFieldAccessor:
        """Return one registry-backed grouped UZF field accessor."""

        field_spec = get_package_input_field_spec("uzf", field_name)
        if field_spec is None:
            raise AttributeError(f"{type(self).__name__!s} has no UZF input field {field_name!r}")
        return GroupUzfFieldAccessor(self.group, field_spec.name)

    def __getattr__(self, field_name: str) -> GroupUzfFieldAccessor:
        """Return a registry-backed grouped UZF field accessor."""

        return self._field(field_name)

    @property
    def finf(self) -> GroupUzfFieldAccessor:
        """Return grouped infiltration accessors."""

        return self._field("finf")

    @property
    def pet(self) -> GroupUzfFieldAccessor:
        """Return grouped potential evapotranspiration accessors."""

        return self._field("pet")

    @property
    def extdp(self) -> GroupUzfFieldAccessor:
        """Return grouped ET extinction-depth accessors."""

        return self._field("extdp")

    @property
    def extwc(self) -> GroupUzfFieldAccessor:
        """Return grouped ET extinction-water-content accessors."""

        return self._field("extwc")

    @property
    def ha(self) -> GroupUzfFieldAccessor:
        """Return grouped surface-depression-storage-depth accessors."""

        return self._field("ha")

    @property
    def hroot(self) -> GroupUzfFieldAccessor:
        """Return grouped root-zone-thickness accessors."""

        return self._field("hroot")

    @property
    def rootact(self) -> GroupUzfFieldAccessor:
        """Return grouped root-activity accessors."""

        return self._field("rootact")


class GroupLakOutputs:
    """Lake-output accessor for :class:`ModelGroup`."""

    def __init__(self, group: "ModelGroup"):
        self.group = group

    def stage(self) -> pd.DataFrame:
        """Return lake stages for all models as one aligned long-format table."""

        rows: list[pd.DataFrame] = []
        for model_name, model in self.group.models.items():
            stage_data = np.asarray(model.outputs.lak.stage.get(), dtype=float)
            if stage_data.ndim == 1:
                stage_data = stage_data.reshape(-1, 1)
            periods = list(model.kstpkper)
            frame = pd.DataFrame(stage_data, columns=list(range(stage_data.shape[1])))
            frame["kstpkper"] = periods[: len(frame)]
            frame = frame.melt(id_vars="kstpkper", var_name="lake", value_name="stage")
            frame["model"] = model_name
            rows.append(frame)

        if not rows:
            return pd.DataFrame(columns=["model", "kstpkper", "lake", "stage"])
        return pd.concat(rows, ignore_index=True)[["model", "kstpkper", "lake", "stage"]]


class GroupLakStageResults(_GroupSpatialView):
    """Grouped accessor for lake stages and stage comparisons.

    Inherits the :class:`SpatialView` grammar (``map``/``mosaic``/``animate``)
    with the model axis being the group's members; each panel delegates to the
    single-model LAK stage explorer so grouped stage maps match the single-model
    ``lak.results.stage.map`` exactly.
    """

    def __init__(self, group: "ModelGroup"):
        self.group = group

    def _spatial_value_label(self) -> str:
        return "stage"

    def _spatial_map(self, *, per: int = 0, layer: int = 0, model=None, **kwargs):
        """Render one model's lake stage broadcast to connected cells."""

        target_model = self.group.models[self._group_target(model)]
        return LakStageResultsExplorer(target_model).map(
            per=int(per), layer=int(layer), backend="plotly", **kwargs
        )

    def get(
        self,
        *,
        model_name: str | None = None,
        lake: int | None = None,
        per: int | None = None,
    ) -> pd.DataFrame:
        """Return aligned lake stages for all models."""

        rows: list[pd.DataFrame] = []
        for current_model_name, model in self.group.models.items():
            stage_data = np.asarray(model.outputs.lak.stage.get(), dtype=float)
            if stage_data.ndim == 1:
                stage_data = stage_data.reshape(-1, 1)
            periods = pd.DataFrame(stage_data, columns=list(range(stage_data.shape[1])))
            periods["per"] = periods.index.astype(int)
            frame = periods.melt(id_vars="per", var_name="lake", value_name="stage")
            frame["lake"] = frame["lake"].astype(int)
            frame["model"] = current_model_name
            rows.append(frame[["model", "per", "lake", "stage"]])

        if not rows:
            return pd.DataFrame(columns=["model", "per", "lake", "stage"])
        combined = pd.concat(rows, ignore_index=True)
        if model_name is not None:
            combined = combined.loc[combined["model"] == str(model_name)].copy()
        if lake is not None:
            combined = combined.loc[combined["lake"] == int(lake)].copy()
        if per is not None:
            combined = combined.loc[combined["per"] == int(per)].copy()
        return combined.reset_index(drop=True)

    def compare(
        self,
        *,
        model_name: str | None = None,
        lake: int | None = None,
        per: int | None = None,
    ) -> pd.DataFrame:
        """Compare lake stages against the reference model."""

        data = self.get(lake=lake, per=per)
        if data.empty:
            return pd.DataFrame(columns=["model", "reference_model", "per", "lake", "stage", "reference_stage", "stage_diff"])

        reference = self.group.reference
        ref = (
            data.loc[data["model"] == reference, ["per", "lake", "stage"]]
            .rename(columns={"stage": "reference_stage"})
            .copy()
        )
        comp = data.loc[data["model"] != reference].merge(ref, on=["per", "lake"], how="inner")
        comp["reference_model"] = reference
        comp["stage_diff"] = comp["stage"].astype(float) - comp["reference_stage"].astype(float)
        if model_name is not None:
            comp = comp.loc[comp["model"] == str(model_name)].copy()
        return comp[["model", "reference_model", "per", "lake", "stage", "reference_stage", "stage_diff"]]

    # NOTE: the series view is the unified grammar's ``plot()`` (SpatialView) --
    # one line per model and lake.


class GroupSfrStageResults(_GroupSpatialView):
    """Grouped accessor for SFR reach stages and stage comparisons.

    Mirrors :class:`GroupLakStageResults` for streams (keyed by ``reach``),
    reading ``model.outputs.sfr.stage``. Inherits the :class:`SpatialView`
    grammar; each panel delegates to the single-model SFR stage explorer.
    """

    def __init__(self, group: "ModelGroup"):
        self.group = group

    def _spatial_value_label(self) -> str:
        return "stage"

    def _spatial_map(self, *, per: int = 0, layer: int = 0, model=None, **kwargs):
        """Render one model's reach stage broadcast to reach cells."""

        target_model = self.group.models[self._group_target(model)]
        return SfrStageResultsExplorer(target_model).map(
            per=int(per), layer=int(layer), backend="plotly", **kwargs
        )

    def get(
        self,
        *,
        model_name: str | None = None,
        reach: int | None = None,
        per: int | None = None,
    ) -> pd.DataFrame:
        """Return aligned SFR reach stages for all models."""

        rows: list[pd.DataFrame] = []
        for current_model_name, model in self.group.models.items():
            stage_data = np.asarray(model.outputs.sfr.stage.get(), dtype=float)
            if stage_data.ndim == 1:
                stage_data = stage_data.reshape(-1, 1)
            periods = pd.DataFrame(stage_data, columns=list(range(stage_data.shape[1])))
            periods["per"] = periods.index.astype(int)
            frame = periods.melt(id_vars="per", var_name="reach", value_name="stage")
            frame["reach"] = frame["reach"].astype(int)
            frame["model"] = current_model_name
            rows.append(frame[["model", "per", "reach", "stage"]])

        if not rows:
            return pd.DataFrame(columns=["model", "per", "reach", "stage"])
        combined = pd.concat(rows, ignore_index=True)
        if model_name is not None:
            combined = combined.loc[combined["model"] == str(model_name)].copy()
        if reach is not None:
            combined = combined.loc[combined["reach"] == int(reach)].copy()
        if per is not None:
            combined = combined.loc[combined["per"] == int(per)].copy()
        return combined.reset_index(drop=True)

    def compare(
        self,
        *,
        model_name: str | None = None,
        reach: int | None = None,
        per: int | None = None,
    ) -> pd.DataFrame:
        """Compare SFR reach stages against the reference model."""

        columns = [
            "model", "reference_model", "per", "reach",
            "stage", "reference_stage", "stage_diff",
        ]
        data = self.get(reach=reach, per=per)
        if data.empty:
            return pd.DataFrame(columns=columns)

        reference = self.group.reference
        ref = (
            data.loc[data["model"] == reference, ["per", "reach", "stage"]]
            .rename(columns={"stage": "reference_stage"})
            .copy()
        )
        comp = data.loc[data["model"] != reference].merge(ref, on=["per", "reach"], how="inner")
        comp["reference_model"] = reference
        comp["stage_diff"] = comp["stage"].astype(float) - comp["reference_stage"].astype(float)
        if model_name is not None:
            comp = comp.loc[comp["model"] == str(model_name)].copy()
        return comp[columns]


class GroupLakConnections:
    """Grouped accessor for lake-connection geometry."""

    def __init__(self, group: "ModelGroup"):
        self.group = group

    def get(
        self,
        *,
        model_name: str | None = None,
        lake: int | None = None,
        layer: int | Sequence[int] | None = None,
        cells: Sequence[int] | None = None,
    ) -> pd.DataFrame:
        """Return aligned LAK connection rows for all models."""

        frames = []
        for current_model_name, model in self.group.models.items():
            frame = build_lak_connection_table(model)
            frame["model"] = current_model_name
            frames.append(frame)

        if not frames:
            return pd.DataFrame()
        combined = pd.concat(frames, ignore_index=True)
        combined = _filter_group_input_table(
            combined,
            model_name=model_name,
            per=None,
            layer=layer,
            cells=cells,
        )
        if lake is not None and "lake" in combined.columns:
            combined = combined.loc[combined["lake"] == int(lake)].copy()
        return combined.reset_index(drop=True)

    def map(
        self,
        *,
        model_name: str | None = None,
        lake: int | None = None,
        layer: int = 0,
        value_column: str = "connection_area",
        agg: str = "sum",
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        colorscale: str | None = None,
        **kwargs,
    ):
        """Build a grouped LAK connection-geometry map for one selected model."""

        target_name = self.group.reference if model_name is None else str(model_name)
        if target_name not in self.group.models:
            raise KeyError(f"Model {target_name!r} is not in the group.")
        target_model = self.group.models[target_name]
        selected = self.get(model_name=target_name, lake=lake, layer=layer)
        if value_column not in selected.columns:
            raise KeyError(f"LAK connection column {value_column!r} was not found.")
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(target_model))
        values, hover = build_cell_input_map_payload(
            selected,
            ncpl=target_model.vor.ncpl,
            value_column=value_column,
            per=None,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
            agg=agg,
        )
        kwargs.setdefault("hover_spec", cell_input_hover(value_column))
        return target_model.cor(
            per=0,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or "earth",
            **kwargs,
        )


class GroupOutputs:
    """Namespace for grouped package-specific output accessors."""

    def __init__(self, group: "ModelGroup"):
        self.group = group

    @property
    def lak(self) -> GroupLakOutputs:
        """Return grouped LAK output helpers."""

        return GroupLakOutputs(self.group)


class GroupPackageAccessor:
    """Namespace for one grouped package's preferred exploration helpers."""

    def __init__(self, accessor):
        self.inputs = accessor

    @property
    def results(self):
        """Return grouped result helpers for this package."""

        budget_text, value_name = get_default_budget_term(self.inputs.package_name) or (
            self.inputs.package_name.upper(),
            "q",
        )
        return GroupCellPackageResultsNamespace(
            GroupCellPackageResults(
                self.inputs.group,
                self.inputs.package_name,
                budget_text=budget_text,
                value_name=value_name,
            )
        )


class GroupResultsOnlyPackageAccessor(Generic[TResultsNamespace]):
    """Namespace for grouped packages that currently expose results only."""

    def __init__(self, results_namespace: TResultsNamespace):
        self._results_namespace = results_namespace

    @property
    def results(self) -> TResultsNamespace:
        """Return grouped result helpers for this package."""

        return self._results_namespace


class GroupUzfPackageAccessor:
    """Namespace for grouped UZF exploration helpers."""

    def __init__(self, accessor: GroupUzfInputs):
        self.inputs = accessor

    @property
    def results(self):
        """Return grouped UZF result helpers."""

        budget_text, value_name = get_default_budget_term("uzf_gwrch") or ("UZF-GWRCH", "gwrch")
        return GroupUzfResultsNamespace(
            GroupCellPackageResults(
                self.inputs.group,
                "uzf",
                budget_text=budget_text,
                value_name=value_name,
            )
        )


class GroupUzfResultsNamespace(FieldMappable):
    """Namespace for grouped UZF result accessors.

    Fields: ``gwrch`` (groundwater recharge, the default) and ``sat``.
    """

    _default_field = "gwrch"

    def _field_names(self):
        return ["gwrch", "sat"]

    def __init__(self, gwrch_accessor: GroupCellPackageResults, sat_accessor: GroupCellPackageResults | None = None):
        self._gwrch = gwrch_accessor
        self._sat = sat_accessor

    @property
    def gwrch(self) -> GroupCellPackageResults:
        """Return grouped groundwater-recharge results from UZF."""

        return self._gwrch

    @property
    def sat(self) -> GroupCellPackageResults:
        """Return grouped UZF saturation results."""

        if self._sat is None:
            budget_text, value_name = get_default_budget_term("uzf_sat") or ("DATA-SAT", "sat")
            self._sat = GroupCellPackageResults(
                self._gwrch.group,
                "uzf",
                budget_text=budget_text,
                value_name=value_name,
            )
        return self._sat


class GroupCellPackageResultsNamespace(FieldMappable):
    """Namespace for grouped cell-based package result accessors.

    ``results.map()`` maps the exchange field ``q`` for one model;
    ``results.mosaic(by="model")`` and ``results.animate(...)`` inherit the
    unified grammar. (SFR/LAK add a ``stage`` field via their subclasses.)
    """

    _default_field = "q"

    def _field_names(self):
        return ["q"]

    def __init__(self, result_accessor: GroupCellPackageResults):
        self._result_accessor = result_accessor

    @property
    def q(self) -> GroupCellPackageResults:
        """Return the primary grouped package-exchange result accessor."""

        return self._result_accessor


class GroupSfrResultsNamespace(GroupCellPackageResultsNamespace):
    """Namespace for grouped SFR result accessors."""

    def _field_names(self):
        return ["q", "stage"]

    @property
    def q(self) -> GroupSfrBudgetResults:
        """Return grouped SFR exchange helpers with normalized map behavior."""

        return self._result_accessor

    @property
    def stage(self) -> "GroupSfrStageResults":
        """Return grouped SFR reach-stage helpers."""

        return GroupSfrStageResults(self._result_accessor.group)


class GroupLakResultsNamespace(GroupCellPackageResultsNamespace):
    """Namespace for grouped LAK result accessors."""

    def _field_names(self):
        return ["q", "stage"]

    @property
    def stage(self) -> GroupLakStageResults:
        """Return grouped LAK stage helpers."""

        return GroupLakStageResults(self._result_accessor.group)

    @property
    def q(self) -> GroupLakBudgetResults:
        """Return grouped LAK exchange helpers with area-normalized map behavior."""

        return self._result_accessor


class GroupLakPackageAccessor:
    """Namespace for grouped LAK geometry and result helpers."""

    def __init__(self, group: "ModelGroup", results_namespace: GroupLakResultsNamespace):
        self.group = group
        self._results_namespace = results_namespace

    @property
    def connections(self) -> GroupLakConnections:
        """Return grouped LAK connection-geometry helpers."""

        return GroupLakConnections(self.group)

    @property
    def results(self) -> GroupLakResultsNamespace:
        """Return grouped LAK result helpers."""

        return self._results_namespace


class GroupSurfaceWaterExchangeResults(_GroupSpatialView):
    """Grouped combined SFR/LAK exchange accessor (unified map/mosaic/animate)."""

    def __init__(self, group: "ModelGroup"):
        self.group = group

    def get(
        self,
        *,
        model_name: str | None = None,
        per: int | None = None,
        layer: int | Sequence[int] | None = None,
        include: str | Sequence[str] | None = None,
        lak_connection_type: str | Sequence[str] | None = None,
    ) -> pd.DataFrame:
        """Return combined surface-water exchange rows for all selected models."""

        frames = []
        for current_model_name, model in self.group.models.items():
            frame = build_surface_water_exchange_cell_table(
                model,
                per=per,
                layer=layer,
                include=include,
                lak_connection_type=lak_connection_type,
            )
            if frame.empty:
                continue
            frame["model"] = current_model_name
            frames.append(frame)
        if not frames:
            return pd.DataFrame()
        combined = pd.concat(frames, ignore_index=True)
        return _filter_group_input_table(combined, model_name=model_name, per=per, layer=layer, cells=None)

    def _spatial_map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        model=None,
        include: str | Sequence[str] | None = None,
        lak_connection_type: str | Sequence[str] | None = None,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        colorscale=None,
        **kwargs,
    ):
        """Render one model's combined SFR/LAK exchange as a ``Choro``."""

        target_name = self._group_target(model)
        target_model = self.group.models[target_name]
        selected = self.get(
            model_name=target_name,
            per=per,
            layer=layer,
            include=include,
            lak_connection_type=lak_connection_type,
        )
        kwargs.setdefault("show_layer_elevs", _default_show_layer_elevs(target_model))
        values, hover = build_surface_water_q_map_payload(
            selected,
            ncpl=target_model.vor.ncpl,
            per=per,
            layer=layer,
            multiplier=multiplier,
            fill_value=fill_value,
        )
        absmax = _symmetric_color_limit(values)
        kwargs.setdefault("zmin", -absmax if absmax > 0 else None)
        kwargs.setdefault("zmax", absmax if absmax > 0 else None)
        kwargs.setdefault("zmid", 0.0)
        kwargs.setdefault("hover_spec", surface_water_hover())
        return target_model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            # signed exchange: gaining (negative) blue, losing red, like SFR/LAK
            colorscale=colorscale or _blue_white_red_diverging_colorscale(),
            **kwargs,
        )


class GroupSurfaceWaterResultsNamespace(FieldMappable):
    """Namespace for grouped combined surface-water result helpers (field ``q``)."""

    _default_field = "q"

    def _field_names(self):
        return ["q"]

    def __init__(self, group: "ModelGroup"):
        self.group = group

    @property
    def q(self) -> GroupSurfaceWaterExchangeResults:
        """Return grouped combined SFR/LAK exchange helpers."""

        return GroupSurfaceWaterExchangeResults(self.group)


class GroupPackages:
    """Preferred grouped package exploration namespace.

    This mirrors the single-model ``model.packages`` surface where practical
    while reusing the existing grouped ``get()/compare()`` accessors.
    """

    def __init__(self, group: "ModelGroup"):
        self.group = group

    @property
    def rch(self) -> GroupPackageAccessor:
        """Grouped recharge input helpers."""

        return GroupPackageAccessor(self.group._rch)

    @property
    def chd(self) -> GroupPackageAccessor:
        """Grouped constant-head input helpers."""

        return GroupPackageAccessor(self.group._chd)

    @property
    def drn(self) -> GroupPackageAccessor:
        """Grouped drain input helpers."""

        return GroupPackageAccessor(self.group._drn)

    @property
    def ghb(self) -> GroupPackageAccessor:
        """Grouped general-head-boundary input helpers."""

        return GroupPackageAccessor(self.group._ghb)

    @property
    def wel(self) -> GroupPackageAccessor:
        """Grouped well package input helpers."""

        return GroupPackageAccessor(self.group._wel)

    @property
    def uzf(self) -> GroupUzfPackageAccessor:
        """Grouped UZF input helpers."""

        return GroupUzfPackageAccessor(self.group._uzf)

    @property
    def sfr(self) -> GroupResultsOnlyPackageAccessor[GroupSfrResultsNamespace]:
        """Grouped SFR result helpers."""

        budget_text, value_name = get_default_budget_term("sfr") or ("SFR", "q")
        return GroupResultsOnlyPackageAccessor(
            GroupSfrResultsNamespace(
                GroupSfrBudgetResults(
                    self.group,
                    budget_text=budget_text,
                    value_name=value_name,
                )
            )
        )

    @property
    def lak(self) -> GroupLakPackageAccessor:
        """Grouped LAK geometry and result helpers."""

        budget_text, value_name = get_default_budget_term("lak") or ("GWF", "q")
        results_namespace = GroupLakResultsNamespace(
            GroupLakBudgetResults(
                self.group,
                budget_text=budget_text,
                value_name=value_name,
            )
        )
        return GroupLakPackageAccessor(self.group, results_namespace)

    @property
    def surface_water(self) -> GroupResultsOnlyPackageAccessor[GroupSurfaceWaterResultsNamespace]:
        """Grouped combined surface-water result helpers."""

        return GroupResultsOnlyPackageAccessor(GroupSurfaceWaterResultsNamespace(self.group))


class ModelGroup:
    """Collection of models with group-wise get/compare helpers.

    Parameters
    ----------
    models
        Mapping of names to models/workspace directories, or a sequence of
        model objects/workspace directories. When workspace directories are
        provided, they are opened lazily as :class:`LoadedMf6Run` objects.
        When a sequence is provided, each loaded model's ``name`` attribute is
        used as the group key.
    reference
        Optional reference model name. If omitted, the first model is used as
        the comparison baseline for ``compare()`` methods.
    crs
        CRS passed through when opening file-backed workspaces.
    verbosity_level
        Verbosity passed through when opening file-backed workspaces.
    shared_grid
        If ``True``, confirm all models have the same discretization and reuse
        the reference model's lazily built Voronoi/grid view for the whole
        group. This avoids rebuilding identical grids repeatedly.
    """

    def __init__(
        self,
        models,
        *,
        reference: str | None = None,
        crs: str = "EPSG:2927",
        verbosity_level: int = 0,
        shared_grid: bool = False,
    ):
        self.models = _coerce_models(models, crs=crs, verbosity_level=verbosity_level)
        if not self.models:
            raise ValueError("ModelGroup requires at least one model.")

        self.reference = reference or next(iter(self.models))
        if self.reference not in self.models:
            raise KeyError(f"Reference model {self.reference!r} is not in the group.")
        self.shared_grid = bool(shared_grid)

        if self.shared_grid:
            self._configure_shared_grid()

        self.hds = GroupHeads(self)
        self.outputs = GroupOutputs(self)
        self._rch = GroupPackageInputs(self, "rch")
        self._chd = GroupPackageInputs(self, "chd")
        self._drn = GroupPackageInputs(self, "drn")
        self._ghb = GroupPackageInputs(self, "ghb")
        self._wel = GroupPackageInputs(self, "wel")
        self._uzf = GroupUzfInputs(self)
        self.packages = GroupPackages(self)

    # -- deprecated flat input shortcuts --------------------------------------
    # These duplicate ``group.packages.<pkg>.inputs`` (which mirrors the
    # single-model ``model.packages.<pkg>.inputs``) and have no single-model
    # equivalent, so they are deprecated in favor of the mirrored path.
    _DEPRECATED_ATTRS = ("rch", "chd", "drn", "ghb", "wel", "uzf")

    def __dir__(self):
        """Hide the deprecated flat shortcuts from tab-completion / dir().

        They still work (with a DeprecationWarning) for back-compat, but should
        not be advertised -- use ``group.packages.<pkg>.inputs`` instead.
        """

        return [name for name in super().__dir__() if name not in self._DEPRECATED_ATTRS]

    @property
    def rch(self) -> "GroupPackageInputs":
        """Deprecated. Use ``group.packages.rch.inputs``."""

        _warn_deprecated("ModelGroup.rch", "group.packages.rch.inputs")
        return self._rch

    @property
    def chd(self) -> "GroupPackageInputs":
        """Deprecated. Use ``group.packages.chd.inputs``."""

        _warn_deprecated("ModelGroup.chd", "group.packages.chd.inputs")
        return self._chd

    @property
    def drn(self) -> "GroupPackageInputs":
        """Deprecated. Use ``group.packages.drn.inputs``."""

        _warn_deprecated("ModelGroup.drn", "group.packages.drn.inputs")
        return self._drn

    @property
    def ghb(self) -> "GroupPackageInputs":
        """Deprecated. Use ``group.packages.ghb.inputs``."""

        _warn_deprecated("ModelGroup.ghb", "group.packages.ghb.inputs")
        return self._ghb

    @property
    def wel(self) -> "GroupPackageInputs":
        """Deprecated. Use ``group.packages.wel.inputs``."""

        _warn_deprecated("ModelGroup.wel", "group.packages.wel.inputs")
        return self._wel

    @property
    def uzf(self) -> "GroupUzfInputs":
        """Deprecated. Use ``group.packages.uzf.inputs``."""

        _warn_deprecated("ModelGroup.uzf", "group.packages.uzf.inputs")
        return self._uzf

    @property
    def run_ids(self) -> list[str]:
        """Return the ordered model names in the group."""

        return list(self.models)

    def summary(self) -> pd.DataFrame:
        """Return per-model summaries stacked into one table."""

        frames = []
        for model_name, model in self.models.items():
            frame = model.summary().copy()
            frame["group_model"] = model_name
            frame["is_reference"] = model_name == self.reference
            frames.append(frame)
        return pd.concat(frames, ignore_index=True)

    def bud(self, package: str | None = None) -> GroupBudget:
        """Return a grouped budget accessor.

        Parameters
        ----------
        package
            Optional default package name such as ``"rch"`` or ``"drn"``.
        """

        return GroupBudget(self, package=package)

    def diff(self):
        """Return a :class:`~myflopy.project.model_diff.ModelDiff` for this group.

        Reference-star: every other model is compared against the group's
        reference. See :class:`~myflopy.project.model_diff.ModelDiff` for the
        structural + value difference tiers, ``summary()``, and ``report()``.
        """

        from myflopy.project.model_diff import ModelDiff

        return ModelDiff(self)

    def _configure_shared_grid(self):
        """Enable shared lazy grid reuse across models with identical grids."""

        anchor_name = self.reference
        anchor_model = self.models[anchor_name]
        anchor_signature = _grid_signature_for_model(anchor_model)
        if anchor_signature is None:
            raise ValueError(
                "shared_grid=True requires discretization files or grid signatures "
                f"for reference model {anchor_name!r}."
            )

        for model_name, model in self.models.items():
            signature = _grid_signature_for_model(model)
            if signature is None:
                raise ValueError(
                    "shared_grid=True requires discretization files or grid signatures "
                    f"for model {model_name!r}."
                )
            if signature != anchor_signature:
                raise ValueError(
                    "shared_grid=True requires all models to use the same discretization. "
                    f"Model {model_name!r} does not match reference model {anchor_name!r}."
                )

        for model_name, model in self.models.items():
            if model_name == anchor_name:
                continue
            setattr(model, "_shared_vor_source", anchor_model)

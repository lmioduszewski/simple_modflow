"""Group-oriented helpers for comparing multiple models with one lazy API."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import hashlib
from pathlib import Path
from typing import Generic, TypeVar

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors as mcolors
import numpy as np
import pandas as pd
from simple_modflow.modflow.mf6.package_explorer import (
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

TResultsNamespace = TypeVar("TResultsNamespace")


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


def _coerce_panel_model_names(group: "ModelGroup", model_names: Sequence[str] | None = None) -> list[str]:
    """Normalize optional subplot model ordering to a validated name list."""

    if model_names is None:
        return list(group.models.keys())
    normalized = [str(name) for name in model_names]
    missing = [name for name in normalized if name not in group.models]
    if missing:
        raise KeyError(f"Models not found in group: {missing!r}")
    return normalized


def _coerce_matplotlib_colormap(colorscale) -> mcolors.Colormap:
    """Convert a package-explorer colorscale into a Matplotlib colormap."""

    if colorscale is None:
        return plt.get_cmap("RdBu")
    if isinstance(colorscale, str):
        try:
            return plt.get_cmap(colorscale)
        except ValueError:
            return plt.get_cmap(colorscale.lower())
    if isinstance(colorscale, Sequence):
        color_values = []
        for entry in colorscale:
            if isinstance(entry, Sequence) and len(entry) >= 2:
                color_values.append(entry[1])
            else:
                color_values.append(entry)
        return mcolors.LinearSegmentedColormap.from_list("simple_modflow_surface_water", color_values)
    raise TypeError("colorscale must be None, a Matplotlib colormap name, or a Plotly-style colorscale list.")


def _plot_group_choropleth_subplots(
    group: "ModelGroup",
    panel_values: dict[str, np.ndarray],
    *,
    colorbar_label: str,
    title_prefix: str,
    colorscale,
    model_names: Sequence[str] | None = None,
    ncols: int | None = None,
    figsize: tuple[float, float] | None = None,
    symmetric: bool = True,
):
    """Plot one shared-scale choropleth panel per model and return the figure."""

    ordered_names = _coerce_panel_model_names(group, model_names)
    if not ordered_names:
        raise ValueError("At least one model must be selected for subplot_map().")

    ncols = int(ncols) if ncols is not None else min(3, max(1, len(ordered_names)))
    nrows = int(np.ceil(len(ordered_names) / ncols))
    if figsize is None:
        figsize = (5.0 * ncols, 4.75 * nrows)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
    axes_array = np.atleast_1d(axes).ravel()

    all_values = np.concatenate([np.asarray(panel_values[name], dtype=float) for name in ordered_names])
    finite = all_values[np.isfinite(all_values)]
    if symmetric:
        absmax = float(np.max(np.abs(finite))) if finite.size else 0.0
        if absmax <= 0.0:
            absmax = 1.0
        vmin, vmax = -absmax, absmax
    elif finite.size:
        vmin = float(np.min(finite))
        vmax = float(np.max(finite))
        if vmin == vmax:
            padding = abs(vmin) * 0.05 if vmin != 0.0 else 1.0
            vmin -= padding
            vmax += padding
    else:
        vmin, vmax = 0.0, 1.0
    cmap = _coerce_matplotlib_colormap(colorscale)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    for axis, model_name in zip(axes_array, ordered_names, strict=False):
        model = group.models[model_name]
        gdf = model.vor.gdf_vorPolys.copy()
        gdf["value"] = np.asarray(panel_values[model_name], dtype=float)
        gdf.plot(
            column="value",
            ax=axis,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            linewidth=0.3,
            edgecolor="#666666",
        )
        axis.set_title(str(model_name))
        axis.set_axis_off()
        axis.set_aspect("equal")

    for axis in axes_array[len(ordered_names):]:
        axis.set_visible(False)

    scalar_mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar_mappable.set_array([])
    colorbar = fig.colorbar(scalar_mappable, ax=axes_array[: len(ordered_names)], shrink=0.9)
    colorbar.set_label(colorbar_label)
    fig.suptitle(title_prefix, y=0.98)
    fig.subplots_adjust(top=0.9, wspace=0.08, hspace=0.12)
    fig._simple_modflow_panel_values = {name: np.asarray(panel_values[name], dtype=float) for name in ordered_names}
    return fig


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


def _resolve_group_compare_target(group: "ModelGroup", model_name: str | None) -> str:
    """Resolve which non-reference model to use for a group difference map."""

    if model_name is not None:
        if model_name == group.reference:
            raise ValueError("compare_map requires a non-reference model name.")
        if model_name not in group.models:
            raise KeyError(f"Model {model_name!r} is not in the group.")
        return str(model_name)

    non_reference = [name for name in group.run_ids if name != group.reference]
    if len(non_reference) == 1:
        return non_reference[0]
    raise ValueError(
        "compare_map requires model_name when the group has more than one non-reference model."
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
        Passed to :func:`simple_modflow.project.load_mf6_run` when a workspace
        path is provided.
    """

    if isinstance(model_or_path, (str, Path)):
        from simple_modflow.project.run_model import load_mf6_run

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


class GroupHeads:
    """Heads accessor for :class:`ModelGroup`.

    This mirrors the single-model ``model.hds`` surface conceptually while
    returning aligned multi-model data and derived comparisons.
    """

    def __init__(self, group: "ModelGroup"):
        self.group = group

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
        if data.empty:
            return pd.DataFrame(
                columns=[
                    "model",
                    "reference_model",
                    "kstpkper",
                    "layer",
                    "cell",
                    "elev",
                    "reference_elev",
                    "diff",
                ]
            )

        reference = self.group.reference
        ref = (
            data[data["model"] == reference]
            .rename(columns={"elev": "reference_elev"})
            .drop(columns=["model"])
        )
        comp = data[data["model"] != reference].merge(
            ref,
            on=["kstpkper", "layer", "cell"],
            how="inner",
        )
        comp["reference_model"] = reference
        comp["diff"] = comp["elev"].astype(float) - comp["reference_elev"].astype(float)
        return comp[
            ["model", "reference_model", "kstpkper", "layer", "cell", "elev", "reference_elev", "diff"]
        ]


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


class GroupPackageInputs:
    """Grouped input accessor for simple cell-based MF6 stress-period packages."""

    def __init__(self, group: "ModelGroup", package_name: str):
        self.group = group
        self.package_name = str(package_name).lower()

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

    def map(
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
        """Build a grouped raw-value choropleth for one selected model.

        Parameters
        ----------
        model_name
            Model to render. Defaults to the group's reference model.
        per, layer
            Zero-based stress period and layer to map.
        value_column
            Numeric package-input field to display. If omitted, a package
            default is used when available.
        multiplier, fill_value, agg
            Passed through to the shared cell-map payload builder.
        colorscale
            Optional colorscale override.
        kwargs
            Forwarded to ``model.cor(...)`` on the selected model.
        """

        target_name = self.group.reference if model_name is None else str(model_name)
        if target_name not in self.group.models:
            raise KeyError(f"Model {target_name!r} is not in the group.")
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


class GroupCellPackageResults:
    """Grouped accessor for cell-based package result tables."""

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

        reference = self.group.reference
        value_columns = [self.value_name]
        key_columns = [column for column in data.columns if column not in {"model", *value_columns}]
        ref = data[data["model"] == reference].drop(columns=["model"]).copy()
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

    def map(
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
        """Build a grouped raw-value choropleth for one selected result field."""

        target_name = self.group.reference if model_name is None else str(model_name)
        if target_name not in self.group.models:
            raise KeyError(f"Model {target_name!r} is not in the group.")
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
                or "Viridis"
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

    def subplot_map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        model_names: Sequence[str] | None = None,
        ncols: int | None = None,
        figsize: tuple[float, float] | None = None,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        agg: str = "sum",
        colorscale=None,
        symmetric: bool | None = None,
    ):
        """Plot one grouped package-result choropleth panel per model."""

        ordered_names = _coerce_panel_model_names(self.group, model_names)
        panel_values: dict[str, np.ndarray] = {}
        for current_model_name in ordered_names:
            model = self.group.models[current_model_name]
            selected = self.get(model_name=current_model_name, per=per, layer=layer)
            values, _hover = build_cell_input_map_payload(
                selected,
                ncpl=model.vor.ncpl,
                value_column=self.value_name,
                per=per,
                layer=layer,
                multiplier=multiplier,
                fill_value=fill_value,
                agg=agg,
            )
            panel_values[current_model_name] = np.asarray(values, dtype=float)

        use_symmetric = self.value_name == "q" if symmetric is None else bool(symmetric)
        return _plot_group_choropleth_subplots(
            self.group,
            panel_values,
            colorbar_label=f"{self.package_name.upper()} {self.value_name}",
            title_prefix=f"Grouped {self.package_name.upper()} {self.value_name} (per={per}, layer={layer})",
            colorscale=(
                colorscale
                or ("RdBu" if use_symmetric else None)
                or get_default_package_colorscale(self.package_name)
                or "Viridis"
            ),
            model_names=ordered_names,
            ncols=ncols,
            figsize=figsize,
            symmetric=use_symmetric,
        )

    def plot_timeseries(
        self,
        *,
        cells: int | Sequence[int] | None = None,
        layer: int | Sequence[int] | None = None,
        model_names: Sequence[str] | None = None,
        agg: str = "sum",
        ax=None,
        return_fig: bool = True,
    ):
        """Plot grouped package results by stress period for selected cells."""

        selected_cells = [int(cells)] if isinstance(cells, (int, np.integer)) else cells
        ordered_names = _coerce_panel_model_names(self.group, model_names)
        frame = self.get(layer=layer, cells=selected_cells)
        if ordered_names:
            frame = frame.loc[frame["model"].isin(ordered_names)].copy()
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 4))
        else:
            fig = ax.figure
        if frame.empty:
            ax.set_title(f"Grouped {self.package_name.upper()} {self.value_name} by stress period")
            ax.set_xlabel("Stress Period")
            ax.set_ylabel(self.value_name)
            if return_fig:
                return fig
            return None

        grouped_keys = [column for column in ("model", "layer", "cell") if column in frame.columns]
        for key, group in frame.groupby(grouped_keys, dropna=False):
            if not isinstance(key, tuple):
                key = (key,)
            key_map = dict(zip(grouped_keys, key, strict=False))
            series = group.groupby("per", as_index=False)[self.value_name].agg(agg).sort_values("per")
            model_label = str(key_map.get("model", "model"))
            layer_label = f"L{int(key_map['layer'])} " if "layer" in key_map and pd.notna(key_map["layer"]) else ""
            cell_label = f"C{int(key_map['cell'])}" if "cell" in key_map and pd.notna(key_map["cell"]) else "All cells"
            ax.plot(
                series["per"].astype(int).to_numpy(),
                series[self.value_name].astype(float).to_numpy(),
                marker="o",
                linewidth=2.0,
                label=f"{model_label} / {layer_label}{cell_label}",
            )

        ax.set_title(f"Grouped {self.package_name.upper()} {self.value_name} by stress period")
        ax.set_xlabel("Stress Period")
        ax.set_ylabel(self.value_name)
        ax.legend()
        fig.tight_layout()
        if return_fig:
            return fig
        return None


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

        reference = self.group.reference
        value_columns = ["q", "q_per_length"] if "q_per_length" in data.columns else ["q"]
        key_columns = [column for column in data.columns if column not in {"model", *value_columns}]
        ref = data[data["model"] == reference].drop(columns=["model"]).copy()
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

    def map(
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
        """Build a grouped SFR map normalized by total reach length per cell.

        Positive MF6 ``SFR``/``GWF`` exchange means flow from the stream to
        groundwater, so the default diverging colorscale is defined explicitly
        to render gaining reaches blue and losing reaches red.
        """

        del agg
        target_name = self.group.reference if model_name is None else str(model_name)
        if target_name not in self.group.models:
            raise KeyError(f"Model {target_name!r} is not in the group.")
        target_model = self.group.models[target_name]
        selected = self.get(model_name=target_name, per=per, layer=layer)
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

    def subplot_map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        model_names: Sequence[str] | None = None,
        ncols: int | None = None,
        figsize: tuple[float, float] | None = None,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        colorscale=None,
    ):
        """Plot one SFR exchange choropleth panel per model with a shared scale."""

        ordered_names = _coerce_panel_model_names(self.group, model_names)
        panel_values: dict[str, np.ndarray] = {}
        for current_model_name in ordered_names:
            model = self.group.models[current_model_name]
            selected = self.get(model_name=current_model_name, per=per, layer=layer)
            values, _hover = build_sfr_q_map_payload(
                selected,
                ncpl=model.vor.ncpl,
                per=per,
                layer=layer,
                multiplier=multiplier,
                fill_value=fill_value,
            )
            panel_values[current_model_name] = np.asarray(values, dtype=float)
        return _plot_group_choropleth_subplots(
            self.group,
            panel_values,
            colorbar_label="SFR exchange per reach length",
            title_prefix=f"Grouped SFR exchange (per={per}, layer={layer})",
            colorscale=colorscale or _blue_white_red_diverging_colorscale(),
            model_names=ordered_names,
            ncols=ncols,
            figsize=figsize,
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

    def map(
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
        """Build a grouped LAK map normalized by total exchange area per cell."""

        del agg
        target_name = self.group.reference if model_name is None else str(model_name)
        if target_name not in self.group.models:
            raise KeyError(f"Model {target_name!r} is not in the group.")
        target_model = self.group.models[target_name]
        selected = self.get(model_name=target_name, per=per, layer=layer, connection_type=connection_type)
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
        return target_model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or "RdBu",
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

    def subplot_map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        connection_type: str | Sequence[str] | None = None,
        model_names: Sequence[str] | None = None,
        ncols: int | None = None,
        figsize: tuple[float, float] | None = None,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        colorscale=None,
    ):
        """Plot one LAK exchange choropleth panel per model with a shared scale."""

        ordered_names = _coerce_panel_model_names(self.group, model_names)
        panel_values: dict[str, np.ndarray] = {}
        for current_model_name in ordered_names:
            model = self.group.models[current_model_name]
            selected = self.get(
                model_name=current_model_name,
                per=per,
                layer=layer,
                connection_type=connection_type,
            )
            values, _hover = build_lak_q_map_payload(
                selected,
                ncpl=model.vor.ncpl,
                per=per,
                layer=layer,
                multiplier=multiplier,
                fill_value=fill_value,
            )
            panel_values[current_model_name] = np.asarray(values, dtype=float)
        return _plot_group_choropleth_subplots(
            self.group,
            panel_values,
            colorbar_label="LAK exchange per flow area",
            title_prefix=f"Grouped LAK exchange (per={per}, layer={layer})",
            colorscale=colorscale or "RdBu",
            model_names=ordered_names,
            ncols=ncols,
            figsize=figsize,
        )


class GroupUzfFieldAccessor:
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

    def map(
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
        """Build a grouped raw-value choropleth for one selected UZF field."""

        target_name = self.group.reference if model_name is None else str(model_name)
        if target_name not in self.group.models:
            raise KeyError(f"Model {target_name!r} is not in the group.")
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
        return target_model.cor(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or get_default_package_colorscale(f"uzf_{self.field_name}") or "Viridis",
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


class GroupLakStageResults:
    """Grouped accessor for lake stages and stage comparisons."""

    def __init__(self, group: "ModelGroup"):
        self.group = group

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

    def plot_timeseries(
        self,
        *,
        lake: int | None = None,
        ax=None,
        return_fig: bool = True,
    ):
        """Plot lake stage over stress periods for every model in the group."""

        frame = self.get(lake=lake)
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 4))
        else:
            fig = ax.figure
        if frame.empty:
            ax.set_title("Grouped LAK stage by stress period")
            ax.set_xlabel("Stress Period")
            ax.set_ylabel("Stage")
            if return_fig:
                return fig
            return None

        for (model_name, lake_id), group in frame.groupby(["model", "lake"], dropna=False):
            group = group.sort_values("per")
            ax.plot(
                group["per"].astype(int).to_numpy(),
                group["stage"].astype(float).to_numpy(),
                marker="o",
                linewidth=2.0,
                label=f"{model_name} / Lake {int(lake_id)}",
            )
        ax.set_title("Grouped LAK stage by stress period")
        ax.set_xlabel("Stress Period")
        ax.set_ylabel("Stage")
        ax.legend()
        fig.tight_layout()
        if return_fig:
            return fig
        return None


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
        return target_model.cor(
            per=0,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            colorscale=colorscale or "Blues",
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


class GroupUzfResultsNamespace:
    """Namespace for grouped UZF result accessors."""

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


class GroupCellPackageResultsNamespace:
    """Namespace for grouped cell-based package result accessors."""

    def __init__(self, result_accessor: GroupCellPackageResults):
        self._result_accessor = result_accessor

    @property
    def q(self) -> GroupCellPackageResults:
        """Return the primary grouped package-exchange result accessor."""

        return self._result_accessor


class GroupSfrResultsNamespace(GroupCellPackageResultsNamespace):
    """Namespace for grouped SFR result accessors."""

    @property
    def q(self) -> GroupSfrBudgetResults:
        """Return grouped SFR exchange helpers with normalized map behavior."""

        return self._result_accessor


class GroupLakResultsNamespace(GroupCellPackageResultsNamespace):
    """Namespace for grouped LAK result accessors."""

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


class GroupSurfaceWaterExchangeResults:
    """Grouped combined SFR/LAK exchange accessor with shared L/T subplots."""

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

    def subplot_map(
        self,
        *,
        per: int = 0,
        layer: int = 0,
        include: str | Sequence[str] | None = None,
        lak_connection_type: str | Sequence[str] | None = None,
        model_names: Sequence[str] | None = None,
        ncols: int | None = None,
        figsize: tuple[float, float] | None = None,
        multiplier: float = 1.0,
        fill_value: float = 0.0,
        colorscale=None,
    ):
        """Plot one combined SFR/LAK exchange panel per model with a shared scale."""

        ordered_names = _coerce_panel_model_names(self.group, model_names)
        panel_values: dict[str, np.ndarray] = {}
        for current_model_name in ordered_names:
            model = self.group.models[current_model_name]
            selected = self.get(
                model_name=current_model_name,
                per=per,
                layer=layer,
                include=include,
                lak_connection_type=lak_connection_type,
            )
            values, _hover = build_surface_water_q_map_payload(
                selected,
                ncpl=model.vor.ncpl,
                per=per,
                layer=layer,
                multiplier=multiplier,
                fill_value=fill_value,
            )
            panel_values[current_model_name] = np.asarray(values, dtype=float)
        return _plot_group_choropleth_subplots(
            self.group,
            panel_values,
            colorbar_label="Surface-water exchange intensity",
            title_prefix=f"Grouped surface-water exchange (per={per}, layer={layer})",
            colorscale=colorscale or "RdBu",
            model_names=ordered_names,
            ncols=ncols,
            figsize=figsize,
        )


class GroupSurfaceWaterResultsNamespace:
    """Namespace for grouped combined surface-water result helpers."""

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

        return GroupPackageAccessor(self.group.rch)

    @property
    def chd(self) -> GroupPackageAccessor:
        """Grouped constant-head input helpers."""

        return GroupPackageAccessor(self.group.chd)

    @property
    def drn(self) -> GroupPackageAccessor:
        """Grouped drain input helpers."""

        return GroupPackageAccessor(self.group.drn)

    @property
    def ghb(self) -> GroupPackageAccessor:
        """Grouped general-head-boundary input helpers."""

        return GroupPackageAccessor(self.group.ghb)

    @property
    def wel(self) -> GroupPackageAccessor:
        """Grouped well package input helpers."""

        return GroupPackageAccessor(self.group.wel)

    @property
    def uzf(self) -> GroupUzfPackageAccessor:
        """Grouped UZF input helpers."""

        return GroupUzfPackageAccessor(self.group.uzf)

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
        self.rch = GroupPackageInputs(self, "rch")
        self.chd = GroupPackageInputs(self, "chd")
        self.drn = GroupPackageInputs(self, "drn")
        self.ghb = GroupPackageInputs(self, "ghb")
        self.wel = GroupPackageInputs(self, "wel")
        self.uzf = GroupUzfInputs(self)
        self.packages = GroupPackages(self)

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

"""Shared coercion/filtering helpers for the group package."""

from __future__ import annotations

import hashlib
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, TypeVar

import pandas as pd

if TYPE_CHECKING:
    from myflopy.project.group.core import ModelGroup


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


def _resolve_group_compare_target(group: ModelGroup, model_name: str | None) -> str:
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


def _ensure_group_map_compatible(group: ModelGroup, model_name: str):
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



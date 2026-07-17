"""Grouped spatial views: `_GroupSpatialView` base + `GroupHeads`."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pandas as pd

from myflopy.modflow.mf6.package_explorer import (
    SpatialView,
    build_group_input_compare_map_payload,
    get_default_group_compare_colorscale,
)
from myflopy.modflow.utils.datatypes.hover import compare_hover
from myflopy.project.group._shared import (
    _coerce_kstpkper,
    _default_show_layer_elevs,
    _ensure_group_map_compatible,
    _reduce_to_period_end,
    _resolve_group_compare_target,
)

if TYPE_CHECKING:
    from myflopy.project.group.core import ModelGroup


class _GroupSpatialView(SpatialView):
    """:class:`SpatialView` wired for a :class:`ModelGroup`.

    Faceting defaults to one panel per model (``by="model"``) and the model
    axis is the group's members; the host implements ``_spatial_map(per, layer,
    model)`` to render one model's field as a ``Choro``. ``model=None`` means the
    reference model.
    """

    def _spatial_models(self):
        """The group's member model names -- the model axis for faceting."""

        return list(self.group.models)

    def _spatial_reference_model(self):
        """The reference model whose grid/layers frame the panels."""

        return self.group.models[self.group.reference]

    def _spatial_default_facet(self):
        """Group mosaics default to one panel per model."""

        return "model"

    def _spatial_periods(self):
        """Stress periods present in the data (``[0]`` when there is no period axis)."""

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

    def __init__(self, group: ModelGroup):
        """Bind the group heads accessor to its :class:`ModelGroup`."""

        self.group = group

    # -- unified grammar hooks ----------------------------------------------
    def _spatial_map(self, *, per: int = 0, layer: int = 0, model=None, **kwargs):
        """Render one member's heads choropleth (raw heads, not deltas)."""

        target = self.group.models[self._group_target(model)]
        return target.hds.map(per=int(per), layer=int(layer), **kwargs)

    def _spatial_periods(self) -> list[int]:
        """Stress periods present in the reference model's saved heads."""

        reference = self.group.models[self.group.reference]
        return sorted({int(key[1]) for key in reference.kstpkper})

    def _spatial_value_label(self) -> str:
        """The mapped quantity's label -- head."""

        return "head"

    def _series_table(self) -> pd.DataFrame:
        """The period-end head table backing ``plot()`` (one head per model/per/layer/cell)."""

        # one head per (model, per, layer, cell): reduce to period-end saves
        return _reduce_to_period_end(self.get())

    def _series_value_column(self, frame) -> str:
        """The column ``plot()`` draws -- head elevation (``elev``)."""

        return "elev"

    def _series_default_agg(self) -> str:
        """Collapse cells within a plotted line by mean (averaging heads)."""

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



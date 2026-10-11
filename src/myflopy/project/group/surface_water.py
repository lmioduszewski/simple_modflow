"""Grouped combined SFR/LAK surface-water exchange results."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

import pandas as pd

from myflopy.modflow.mf6.package_explorer import (
    FieldMappable,
    _exchange_colorscale,
    _symmetric_color_limit,
    build_surface_water_exchange_cell_table,
    build_surface_water_q_map_payload,
)
from myflopy.modflow.utils.datatypes.hover import surface_water_hover
from myflopy.project.group._shared import _filter_group_input_table
from myflopy.project.group.spatial import _GroupSpatialView

if TYPE_CHECKING:
    from myflopy.project.group.core import ModelGroup


class GroupSurfaceWaterExchangeResults(_GroupSpatialView):
    """Grouped combined SFR/LAK exchange accessor (unified map/mosaic/animate)."""

    def __init__(self, group: ModelGroup):
        """Bind the grouped combined surface-water exchange accessor to ``group``."""

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
        """Return combined surface-water exchange rows for all selected models.

        Each member's ``build_surface_water_exchange_cell_table``: one row per
        stress period, cell and package TYPE (``"sfr"``/``"lak"``) -- the
        member's SFR packages combined by cell, and its LAK packages likewise
        (ledger 352) -- each period read at its LAST saved
        step on that member's own budget-file clock (ledgers 337, 348), plus a
        ``model`` column. A member with no SFR or LAK record adds no rows.

        Parameters
        ----------
        model_name : str or None, optional
            Keep one member's rows, by its name in the group. Default ``None``:
            every member. A name that is not in the group matches nothing and
            returns an empty frame.
        per : int or None, optional
            Keep one stress period, zero-based -- its last saved step on each
            member's own clock. Default ``None``: every period, each at its last
            saved step. A period no member saved returns an empty frame.
        layer : int or list of int or None, optional
            Keep the rows in these zero-based layers (an int, list or tuple; a
            NumPy array raises ``TypeError``). Default ``None``: every layer.
        include : {"sfr", "lak"} or list of them, optional
            Which packages' exchange to read, any case. Default ``None``: both.
            Any other name raises ``ValueError`` naming the two allowed.
        lak_connection_type : str or list of str, optional
            Keep only LAK connections of these MF6 types (``"VERTICAL"``,
            ``"HORIZONTAL"``, ``"EMBEDDEDH"``, ``"EMBEDDEDV"``), compared upper-cased
            with each connection's ``claktype``. Acts on the LAK rows only; SFR rows
            are unaffected. A spelling that is no MF6 type matches no connection
            and leaves no LAK rows, without an error. Default ``None``: every
            connection.

        Returns
        -------
        pandas.DataFrame
            ``model``, ``package``, ``per``, ``layer``, ``cell``, ``source``
            (``"sfr"``/``"lak"``), the cell's summed ``q`` (model length^3/time, in
            MF6's sign for that package) with ``rlen``/``q_per_length`` (SFR) or
            ``flow_area``/``q_per_area`` (LAK), and ``exchange_intensity`` --
            myflopy's one sign across both: positive where the feature GAINS from
            the aquifer (SFR's per-length flow negated, LAK's per-area flow as
            written). Empty when no member has a row.
        """

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
        return target_model.plot.map(
            per=per,
            layer=layer,
            type="custom",
            custom_zs=values,
            custom_hover=hover,
            hover_heads=False,
            hover_ks=False,
            # The combined map draws the normalized exchange_intensity field
            # (positive = the feature gains), so blue sits at the positive end --
            # the "feature" orientation, matching the single-model surface_water
            # map. The old blue-at-negative scale drew gaining cells red here.
            colorscale=colorscale or _exchange_colorscale("feature"),
            **kwargs,
        )


class GroupSurfaceWaterResultsNamespace(FieldMappable):
    """Namespace for grouped combined surface-water result helpers (field ``q``)."""

    _default_field = "q"

    def _field_names(self):
        """The single mappable grouped combined-exchange field: ``q``."""

        return ["q"]

    def __init__(self, group: ModelGroup):
        """Bind the grouped combined surface-water results namespace to ``group``."""

        self.group = group

    @property
    def q(self) -> GroupSurfaceWaterExchangeResults:
        """Return grouped combined SFR/LAK exchange helpers."""

        return GroupSurfaceWaterExchangeResults(self.group)



"""Grouped GWE temperature accessor (``group.temp``)."""

from __future__ import annotations

from myflopy.project.group.spatial import _GroupFieldView


class GroupTemp(_GroupFieldView):
    """Temperature accessor for :class:`ModelGroup` (``group.temp``).

    The energy-transport twin of
    :class:`~myflopy.project.group.spatial.GroupHeads`, configuring the shared
    :class:`~myflopy.project.group.spatial._GroupFieldView` exactly as
    :class:`~myflopy.project.group.conc.GroupConc` does: aligned multi-model
    tables via ``get``/``compare``, plus the unified grammar (``map``/``plot``/
    ``xs``/``mosaic``/``animate``) faceting over the group's members.

    Every member model must be a GWE model; the kind-gated ``temp`` reader
    supplies the error if not.

    Δtemp maps are reached through the one public diff verb --
    ``group.diff().temp.map("variant")``.
    """

    reader_attribute = "temp"
    table_attribute = "all_temp"
    value_column = "temp"
    value_label = "temp"
    #: matches ``temp_hover``'s default; a convention, not a physical law
    value_unit = "°C"


__all__ = ["GroupTemp"]

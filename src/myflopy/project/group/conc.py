"""Grouped GWT concentration accessor (``group.conc``)."""

from __future__ import annotations

from myflopy.project.group.spatial import _GroupFieldView


class GroupConc(_GroupFieldView):
    """Concentration accessor for :class:`ModelGroup` (``group.conc``).

    The transport twin of :class:`~myflopy.project.group.spatial.GroupHeads`,
    and deliberately just a configuration of the shared
    :class:`~myflopy.project.group.spatial._GroupFieldView`: aligned multi-model
    tables via ``get``/``compare``, plus the unified grammar (``map``/``plot``/
    ``xs``/``mosaic``/``animate``) faceting over the group's members.

    Every member model must be a GWT model. ``reader_attribute = "conc"`` is the
    kind-GATED reader, so a group holding a flow model fails with that reader's
    own clear "is a GWF model" error rather than silently mapping heads under a
    concentration label.

    Δconc maps are reached through the one public diff verb --
    ``group.diff().conc.map("variant")`` -- not through ``compare_map``, which is
    internal plumbing.
    """

    reader_attribute = "conc"
    table_attribute = "all_conc"
    value_column = "conc"
    value_label = "conc"
    #: matches ``conc_hover``'s default; model-dependent (mass/volume), a
    #: groundwater convention rather than a physical law
    value_unit = "mg/L"


__all__ = ["GroupConc"]

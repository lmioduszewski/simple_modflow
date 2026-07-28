"""Concentration observation targets for GWT calibration (``ConcTargets``)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

from myflopy.modflow.mf6.observations.heads import HeadTargets


@dataclass
class ConcTargets(HeadTargets):
    """Measured vs simulated CONCENTRATION at monitoring points, for PEST.

    The transport twin of :class:`~myflopy.modflow.mf6.observations.heads.\
HeadTargets`, and deliberately a *configuration* of it rather than a copy:
    concentration is sampled exactly the way head is -- a point value at a
    ``(layer, cell)`` per stress period -- so the only differences are which
    model table is read and what the field is called.

    That is why this subclasses ``HeadTargets`` instead of the named-series
    targets (``DrnFlowTargets`` and friends): those aggregate a FLUX over a zone
    of cells, which is the wrong semantics for a state variable, and they carry
    no layer in their identity while concentration is inherently three
    dimensional.

    Concentration is the observation that makes transport calibration worth
    doing: heads are insensitive to porosity and only weakly constrain flow
    paths, whereas concentration responds strongly to both. Measured on the
    canonical model, a 2.5x porosity change moves concentration by 54% and a
    2.5x K change by 38%.

    Examples
    --------
    ::

        wells = ConcTargets.from_cells(
            model=transport_view, cells=[120, 344], values=measured
        )
        cal.observe(wells)

    Notes
    -----
    The normalized frames keep ``head_target``/``sim_head`` as their internal
    value columns for every field kind. They are frame-internal names meaning
    "the target value" and "the simulated value"; sharing them is what lets the
    calibration-plot and PEST observation machinery serve concentration without
    change (ledger 102).
    """

    _TABLE_ATTRIBUTE: ClassVar[str] = "all_conc"
    _STORE_COLUMN: ClassVar[str] = "conc"
    _FIELD_LABEL: ClassVar[str] = "concentration"

    #: the long-form input column name callers supply measured values under
    value_column: str = "conc"

    def simulated_conc(self, model):
        """Return simulated concentration by stress period, wide by name.

        The concentration spelling of :meth:`HeadTargets.simulated_heads`; both
        read the same table, and either may be called.
        """

        return self.simulated_heads(model)


__all__ = ["ConcTargets"]

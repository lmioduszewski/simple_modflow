from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from myflopy.modflow.mf6.simulation.base import SimulationBase

import numpy as np
import pandas as pd

from myflopy.modflow.mf6.budget import LakBudget, SFRBudget
from myflopy.modflow.mf6.headsplus import _is_sentinel


class LakStageOutput:
    """Lake stage output accessor hosted under ``model.outputs.lak``."""

    def __init__(self, model: SimulationBase):
        """Bind the lake stage-output accessor to ``model``."""

        self.model = model

    def get(self, *, masked: bool = True) -> np.ndarray:
        """Return the LAK stage file as one row per saved record and one column per lake.

        MF6 writes a lake with no water surface as a number, not NaN: -1e30 for a
        DRY lake (its stage at or below the lake bottom) and +1e30 for an INACTIVE
        one (``STATUS INACTIVE``). By default both become ``NaN``, by the same rule
        as the head file (``headsplus._is_sentinel``, ``|value| >= 1e29``), so
        every stage view built on this -- ``packages.lak.results.stage``, its
        maps, series and stage changes, the group and diff tiers -- leaves a gap
        where there is no stage rather than drawing one 1e30 ft away (ledger 329).

        Parameters
        ----------
        masked : bool, default True
            Replace MF6's two no-stage sentinels (-1e30 dry, +1e30 inactive) with
            ``NaN``. ``False`` returns the file exactly as MF6 wrote it, which is
            the only way to tell a dry lake from an inactive one: the tidy tables
            and every picture show both as a gap, as a head view does for a dry
            and an inactive cell.

        Returns
        -------
        numpy.ndarray
            A float array of shape ``(n_records, nlakes)``, lakes in LAK order
            (zero-based column ``i`` is lake ``i + 1`` in the MF6 file). One row
            per SAVED stage record, which is one per stress period only when the
            ``OC`` saves heads once per period: LAK writes its stage whenever MF6
            saves HEAD. Nothing here labels the rows with a time (ledger 301's
            stage-clock family).

        Raises
        ------
        AttributeError
            If the model has no LAK package, or the package declares no
            ``stage_filerecord`` (FloPy's ``output.stage()`` is then ``None``).

        Examples
        --------
        >>> stage = model.outputs.lak.stage.get()          # NaN where a lake is dry
        >>> raw = model.outputs.lak.stage.get(masked=False)
        >>> dry = raw == -1e30                              # dry, not inactive
        """

        nlakes = self.model.fl.lak.nlakes.data
        stages = self.model.fl.lak.output.stage()
        values = np.asarray(stages.get_alldata(), dtype=float).reshape(-1, nlakes)
        if masked:
            values[_is_sentinel(values)] = np.nan
        return values


class SFRStageOutput:
    """Stream stage output accessor hosted under ``model.outputs.sfr``."""

    def __init__(self, model: SimulationBase):
        """Bind the stream stage-output accessor to ``model``."""

        self.model = model

    def get(self):
        """Return a reach-by-saved-timestep dataframe of stream stage output."""

        nreaches = self.model.fl.sfr.nreaches.data
        reader = self.model.fl.sfr.output.stage()
        kstpkpers = [tuple(value) for value in reader.get_kstpkper()]
        stages = reader.get_alldata()
        stage_data = stages.reshape(len(kstpkpers), nreaches).transpose()
        frame = pd.DataFrame(stage_data, columns=kstpkpers)
        frame.index.name = "reaches"
        frame.columns.name = "kstpkper"
        return frame


class LakOutputData:
    """Namespace for LAK-specific output helpers."""

    def __init__(self, model: SimulationBase):
        """Bind the LAK output namespace (``.bud`` / ``.stage``) to ``model``."""

        self.model = model

    @property
    def bud(self):
        """Return the LAK package-output budget accessor."""

        return LakBudget(self.model)

    @property
    def stage(self):
        """Return the LAK stage output accessor.

        The returned helper exposes ``get()`` so existing code can follow the
        same pattern as other output helpers while keeping stage logic under the
        ``outputs`` namespace.
        """

        return LakStageOutput(self.model)


class SFROutputData:
    """Namespace for SFR-specific output helpers."""

    def __init__(self, model: SimulationBase):
        """Bind the SFR output namespace (``.bud`` / ``.stage``) to ``model``."""

        self.model = model

    @property
    def bud(self):
        """Return the SFR package-output budget accessor."""

        return SFRBudget(self.model)

    @property
    def stage(self):
        """Return the SFR stage output accessor.

        The returned helper exposes ``get()`` so existing code can follow the
        same pattern as other output helpers while keeping stage logic under the
        ``outputs`` namespace.
        """

        return SFRStageOutput(self.model)

class UzfOutputData:
    """Namespace for UZF-specific output helpers."""

    def __init__(self, model: SimulationBase):
        """Bind the UZF output namespace to ``model``."""

        self.model = model
        self._ifno_to_cellid = None

    @property
    def ifno_to_cellid(self):
        """Map each UZF ``ifno`` to its ``(layer, cell)`` (vertex grids only; cached)."""

        if self._ifno_to_cellid is None:
            if self.model.modelgrid.grid_type == 'vertex':
                ifno_to_cellid = pd.DataFrame(self.model.fl.uzf.packagedata.get_data().cellid.tolist())
                ifno_to_cellid.index.name = 'ifno'
                ifno_to_cellid.columns = ['layer', 'cellid']
                self._ifno_to_cellid = ifno_to_cellid
            else:
                raise ValueError("Model grid type is not vertex. You'll have to figure this out manually...")
        return self._ifno_to_cellid

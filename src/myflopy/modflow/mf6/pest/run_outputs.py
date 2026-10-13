"""A completed run read through its OUTPUT files alone, for PEST's forward run (ledger 582).

The named-series forward run (``forward_run.write_named_series_targets``) opened
the just-run model with ``mf.load_mf6_run`` -- one full FloPy load per
realization, 8.3-9.0 s on a copy of Ten Trails ``ss`` beside a 17-20 s MF6 run
(ledger 540), 24.7 of 26.2 s of it in FloPy's external-list reader. Nothing
a DRN or RIV zone set reads there is an INPUT a parameter can change: the
packages' names and OBS declarations, the budget file the OC names, the TDIS
lengths and the grid's cell count per layer are fixed when the template is built.
So the build records them
(:func:`myflopy.modflow.mf6.pest.observations.snapshot_run_outputs`, into
``named_series_outputs.json``), and every realization reads its outputs through a
:class:`RunOutputs` built from that record -- the zone family's own
``simulated_series``, unchanged, on an object that answers the five things it
asks of a model and opens nothing but output files.

Light on purpose (numpy and the standard library at import; FloPy's
``CellBudgetFile`` when the budget file is first read), because PEST's forward run
imports it on every model run.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

__all__ = ["OUTPUTS_SNAPSHOT", "RUN_OUTPUT_KINDS", "RunOutputs", "node_to_lni"]

#: The file in the template that records what :class:`RunOutputs` needs.
OUTPUTS_SNAPSHOT = "named_series_outputs.json"

#: The named-series kinds read through :class:`RunOutputs`, and the package type each
#: reads. A lake or SFR set still loads the model: its readers go through FloPy's
#: package output objects and the LAK connection geometry (ledger 582).
RUN_OUTPUT_KINDS = {"drn_flow": "drn", "riv_flow": "riv"}


def node_to_lni(ncpl) -> dict[int, tuple[int, int]]:
    """Flat node -> ``(layer, cell)`` for ``ncpl[k]`` cells in layer ``k``, layer after layer.

    How MF6 numbers a DISV (or DIS) model's cells, and so a budget record's node:
    the build checks the model it records numbers them so
    (``pest.observations.snapshot_run_outputs``).

    Parameters
    ----------
    ncpl : sequence of int
        Cells in each layer, top down.

    Returns
    -------
    dict of int to tuple of (int, int)
        Zero-based node to zero-based ``(layer, cell)``.
    """

    out: dict[int, tuple[int, int]] = {}
    start = 0
    for layer, count in enumerate(ncpl):
        for cell in range(int(count)):
            out[start + cell] = (layer, cell)
        start += int(count)
    return out


class _Continuous:
    """An OBS package's continuous block, as FloPy's ``get_data()`` answers it."""

    def __init__(self, records: dict):
        self._data = {
            csv: [(name, kind, tuple(cellid)) for name, kind, cellid in rows]
            for csv, rows in records.items()
        }

    def get_data(self) -> dict:
        return self._data


class RunOutputs:
    """A completed run read through its output files and a build-time record of its inputs.

    What a DRN or RIV zone set's ``simulated_series`` asks of a model, answered
    without loading one: ``name``, ``workspace``, ``sim.tdis.perioddata``,
    ``node_to_lni``, ``_get_budget_reader()`` (the cell-budget file the OC
    declares) and ``_packages_of_type(type)`` (each package's name and its OBS
    declarations). Built by the forward run from :data:`OUTPUTS_SNAPSHOT`
    (:meth:`from_snapshot`); see the module docstring (ledger 582).

    Parameters
    ----------
    workspace : str or pathlib.Path
        The run directory holding the outputs.
    snapshot : dict
        ``pest.observations.snapshot_run_outputs``' record of the model.
    """

    def __init__(self, workspace, snapshot: dict):
        self.workspace = Path(workspace)
        self.name = str(snapshot["model_name"])
        self._snapshot = snapshot
        records = [tuple(row) for row in snapshot["periods"]]
        data = np.rec.fromrecords(records, names="perlen,nstp,tsmult") if records else None
        self.sim = SimpleNamespace(
            tdis=SimpleNamespace(perioddata=SimpleNamespace(get_data=lambda: data))
        )
        self.node_to_lni = node_to_lni(snapshot["ncpl"])
        self._budget = None

    @classmethod
    def from_snapshot(cls, workspace, path=OUTPUTS_SNAPSHOT) -> RunOutputs:
        """The run in ``workspace`` described by the snapshot file ``path``.

        Parameters
        ----------
        workspace : str or pathlib.Path
            The run directory (a PEST worker's, in the forward run).
        path : str or pathlib.Path, optional
            The snapshot file, relative to the current directory or absolute.
            Default :data:`OUTPUTS_SNAPSHOT`.

        Returns
        -------
        RunOutputs
            The run, read through its outputs.
        """

        return cls(workspace, json.loads(Path(path).read_text(encoding="utf-8")))

    def _get_budget_reader(self):
        """FloPy's reader of the cell-budget file the OC declared, opened once, double precision."""

        if self._budget is None:
            import flopy

            self._budget = flopy.utils.CellBudgetFile(
                str(self.workspace / self._snapshot["budget_file"]), precision="double"
            )
        return self._budget

    def _packages_of_type(self, package_type: str) -> list:
        """Each recorded package of the type: its ``package_name``, type and OBS declarations."""

        kind = str(package_type).lower()
        return [
            SimpleNamespace(
                package_name=entry["name"],
                package_type=kind,
                obs=[SimpleNamespace(continuous=_Continuous(records)) for records in entry["obs"]],
            )
            for entry in self._snapshot["packages"].get(kind, [])
        ]

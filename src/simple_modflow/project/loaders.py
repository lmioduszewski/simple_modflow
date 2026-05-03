"""Lazy loaders for MF6 outputs, pickled model objects, and file-backed runs."""

from __future__ import annotations

from pathlib import Path
import pickle

import flopy
import numpy as np
import pandas as pd

from simple_modflow.project.run_model import LoadedMf6Run, load_mf6_run
from simple_modflow.project.specs import RunRecord


class RunLoader:
    """Load run artifacts and file-backed models from a catalog or explorer."""

    def __init__(self, catalog):
        """Create a loader bound to a catalog-like object with ``load_run`` support."""

        self.catalog = catalog

    def _coerce_run(self, run: str | RunRecord) -> RunRecord:
        """Normalize a run id or record into a concrete :class:`RunRecord`."""

        if isinstance(run, RunRecord):
            return run
        return self.catalog.load_run(run)

    def get_path(self, run: str | RunRecord, key: str) -> Path:
        """Resolve a named path from a run record to an absolute path."""

        record = self._coerce_run(run)
        return record.get_path(key)

    def load_heads_file(self, run: str | RunRecord) -> flopy.utils.HeadFile:
        """Open the MF6 binary heads file for a run."""

        return flopy.utils.HeadFile(self.get_path(run, "heads_file"))

    def load_budget_file(self, run: str | RunRecord) -> flopy.utils.CellBudgetFile:
        """Open the MF6 binary cell-budget file for a run."""

        return flopy.utils.CellBudgetFile(self.get_path(run, "budget_file"))

    def load_model_object(self, run: str | RunRecord):
        """Load the legacy pickled ``.model`` object for a run, if present."""

        with self.get_path(run, "model_object_file").open("rb") as file:
            return pickle.load(file)

    def load_run_model(
        self,
        run: str | RunRecord,
        *,
        crs: str = "EPSG:2927",
        verbosity_level: int = 0,
    ) -> LoadedMf6Run:
        """Open a run as a file-backed :class:`LoadedMf6Run`."""

        record = self._coerce_run(run)
        if verbosity_level > 0:
            print(f"[RunLoader] Opening run '{record.run_id}' from {record.workspace}")
        return LoadedMf6Run.from_run_record(
            record,
            crs=crs,
            verbosity_level=verbosity_level,
        )

    def load_workspace_model(
        self,
        workspace: Path,
        *,
        crs: str = "EPSG:2927",
        verbosity_level: int = 0,
    ) -> LoadedMf6Run:
        """Open an MF6 workspace directly, without an existing run record."""

        if verbosity_level > 0:
            print(f"[RunLoader] Opening workspace model from {Path(workspace)}")
        return load_mf6_run(
            workspace,
            crs=crs,
            verbosity_level=verbosity_level,
        )

    def load_heads_array(
        self,
        run: str | RunRecord,
        *,
        kstpkper: tuple[int, int] | None = None,
    ) -> np.ndarray:
        """Load heads as a NumPy array, defaulting to the final timestep."""

        hds = self.load_heads_file(run)
        if kstpkper is None:
            kstpkper = hds.get_kstpkper()[-1]
        return np.asarray(hds.get_data(kstpkper=kstpkper)).squeeze()

    def load_heads_frame(
        self,
        run: str | RunRecord,
        *,
        kstpkper: tuple[int, int] | None = None,
        layer: int = 0,
        cells: list[int] | None = None,
    ) -> pd.DataFrame:
        """Load heads into a simple tabular form for quick analysis."""

        data = self.load_heads_array(run, kstpkper=kstpkper)
        flat = np.asarray(data).reshape(-1)
        frame = pd.DataFrame(
            {
                "cell": np.arange(flat.size, dtype=int),
                "layer": layer,
                "head": flat.astype(float),
            }
        )
        if cells is not None:
            frame = frame[frame["cell"].isin(cells)].reset_index(drop=True)
        return frame

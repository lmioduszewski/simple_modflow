"""High-level browsing helpers for directories full of MF6 runs."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from simple_modflow.project.compare import RunComparison
from simple_modflow.project.discovery import discover_existing_runs
from simple_modflow.project.loaders import RunLoader
from simple_modflow.project.model_group import ModelGroup
from simple_modflow.project.specs import RunRecord


@dataclass(slots=True)
class RunExplorer:
    """Browse discovered MF6 runs with a simple, notebook-friendly API."""

    root: Path
    records: list[RunRecord]
    name: str | None = None
    loader: RunLoader = field(init=False, repr=False)
    compare: RunComparison = field(init=False, repr=False)

    def __post_init__(self):
        """Initialize helper objects for loading and comparing discovered runs."""

        self.root = Path(self.root)
        self.name = self.name or self.root.name
        self.loader = RunLoader(self)
        self.compare = RunComparison(self)

    @classmethod
    def from_directory(
        cls,
        root: Path,
        *,
        model_spec: str = "explored_runs",
        status: str = "discovered",
        name: str | None = None,
    ):
        """Discover runs below ``root`` and wrap them in a :class:`RunExplorer`."""

        records = [
            discovered.to_run_record(model_spec=model_spec, status=status)
            for discovered in discover_existing_runs(root)
        ]
        return cls(root=Path(root), records=records, name=name)

    def __len__(self) -> int:
        return len(self.records)

    def __iter__(self):
        return iter(self.records)

    def __getitem__(self, run_id: str) -> RunRecord:
        return self.load_run(run_id)

    @property
    def run_ids(self) -> list[str]:
        """Return all discovered run ids."""

        return [record.run_id for record in self.records]

    def load_run(self, run_id: str) -> RunRecord:
        """Return the lightweight run record for one discovered run id."""

        for record in self.records:
            if record.run_id == run_id:
                return record
        raise FileNotFoundError(f"Run not found: {run_id}")

    def open(self, run: str | RunRecord, *, crs: str = "EPSG:2927", verbosity_level: int = 0):
        """Open a discovered run as a file-backed ``LoadedMf6Run`` object."""

        return self.loader.load_run_model(run, crs=crs, verbosity_level=verbosity_level)

    def group(
        self,
        runs: list[str | RunRecord],
        *,
        reference: str | None = None,
        crs: str = "EPSG:2927",
        verbosity_level: int = 0,
    ) -> ModelGroup:
        """Open several runs and return a :class:`ModelGroup`.

        Parameters
        ----------
        runs
            Run ids or records to include in the group, in the desired order.
        reference
            Optional reference model name. Defaults to the first requested run.
        crs, verbosity_level
            Passed through to :meth:`open` for each loaded run.
        """

        loaded = {}
        for run in runs:
            record = self.load_run(run) if not isinstance(run, RunRecord) else run
            loaded[record.run_id] = self.open(record, crs=crs, verbosity_level=verbosity_level)
        return ModelGroup(loaded, reference=reference)

    def summary(self) -> pd.DataFrame:
        """Return a lightweight table summarizing discovered runs."""

        rows = []
        for record in self.records:
            rows.append(
                {
                    "run_id": record.run_id,
                    "model_spec": record.model_spec,
                    "status": record.status,
                    "workspace": str(record.workspace),
                    "family": record.metadata.get("family"),
                    "grid_type": record.metadata.get("grid_type"),
                    "has_model_object": "model_object_file" in record.paths,
                    "has_heads": record.metadata.get("has_heads"),
                    "has_budget": record.metadata.get("has_budget"),
                    "tags": record.tags,
                }
            )
        return pd.DataFrame(rows)

    def families(self) -> pd.DataFrame:
        """Summarize discovered runs by inferred family label."""

        summary = self.summary()
        if summary.empty:
            return pd.DataFrame(columns=["family", "count"])
        family_counts = (
            summary.fillna({"family": "unknown"})
            .groupby("family", dropna=False)
            .size()
            .reset_index(name="count")
            .sort_values(["count", "family"], ascending=[False, True])
            .reset_index(drop=True)
        )
        return family_counts

    def filter(
        self,
        *,
        family: str | None = None,
        grid_type: str | None = None,
        tag: str | None = None,
        text: str | None = None,
        has_heads: bool | None = None,
        has_budget: bool | None = None,
    ):
        """Return a new explorer view filtered by simple metadata criteria."""

        records = list(self.records)
        if family is not None:
            records = [record for record in records if record.metadata.get("family") == family]
        if grid_type is not None:
            records = [record for record in records if record.metadata.get("grid_type") == grid_type]
        if tag is not None:
            records = [record for record in records if tag in record.tags]
        if has_heads is not None:
            records = [record for record in records if bool(record.metadata.get("has_heads")) is has_heads]
        if has_budget is not None:
            records = [record for record in records if bool(record.metadata.get("has_budget")) is has_budget]
        if text:
            lowered = text.lower()
            records = [
                record
                for record in records
                if lowered in record.run_id.lower()
                or lowered in str(record.workspace).lower()
                or any(lowered in tag_value.lower() for tag_value in record.tags)
            ]
        return RunExplorer(root=self.root, records=records, name=self.name)

    def head(self, n: int = 10) -> pd.DataFrame:
        """Return the first ``n`` rows of :meth:`summary`."""

        return self.summary().head(n)


def explore_runs(
    root: Path,
    *,
    model_spec: str = "explored_runs",
    status: str = "discovered",
    name: str | None = None,
) -> RunExplorer:
    """Discover runs below ``root`` and return a :class:`RunExplorer`."""

    return RunExplorer.from_directory(
        root,
        model_spec=model_spec,
        status=status,
        name=name,
    )

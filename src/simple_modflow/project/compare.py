from __future__ import annotations

import pandas as pd

from simple_modflow.project.loaders import RunLoader
from simple_modflow.project.specs import RunRecord


class RunComparison:

    def __init__(self, catalog):
        self.catalog = catalog
        self.loader = RunLoader(catalog)

    def _coerce_run(self, run: str | RunRecord) -> RunRecord:
        if isinstance(run, RunRecord):
            return run
        return self.catalog.load_run(run)

    def compare_package_versions(self, run_a: str | RunRecord, run_b: str | RunRecord) -> pd.DataFrame:
        record_a = self._coerce_run(run_a)
        record_b = self._coerce_run(run_b)
        packages = sorted(set(record_a.package_versions) | set(record_b.package_versions))
        rows = []
        for package in packages:
            version_a = record_a.package_versions.get(package)
            version_b = record_b.package_versions.get(package)
            rows.append(
                {
                    "package": package,
                    "run_a": version_a,
                    "run_b": version_b,
                    "changed": version_a != version_b,
                }
            )
        return pd.DataFrame(rows)

    def compare_heads(
        self,
        run_a: str | RunRecord,
        run_b: str | RunRecord,
        *,
        kstpkper: tuple[int, int] | None = None,
        cells: list[int] | None = None,
    ) -> pd.DataFrame:
        record_a = self._coerce_run(run_a)
        record_b = self._coerce_run(run_b)
        frame_a = self.loader.load_heads_frame(record_a, kstpkper=kstpkper, cells=cells).rename(
            columns={"head": record_a.run_id}
        )
        frame_b = self.loader.load_heads_frame(record_b, kstpkper=kstpkper, cells=cells).rename(
            columns={"head": record_b.run_id}
        )
        merged = frame_a.merge(frame_b, on=["cell", "layer"], how="inner")
        merged["head_diff"] = merged[record_b.run_id] - merged[record_a.run_id]
        return merged

    def compare_head_stats(
        self,
        run_a: str | RunRecord,
        run_b: str | RunRecord,
        *,
        kstpkper: tuple[int, int] | None = None,
        cells: list[int] | None = None,
    ) -> pd.DataFrame:
        head_diff = self.compare_heads(run_a, run_b, kstpkper=kstpkper, cells=cells)
        record_a = self._coerce_run(run_a)
        record_b = self._coerce_run(run_b)
        stats = [
            ("count", float(len(head_diff))),
            ("mean_a", float(head_diff[record_a.run_id].mean())),
            ("mean_b", float(head_diff[record_b.run_id].mean())),
            ("mean_diff", float(head_diff["head_diff"].mean())),
            ("max_abs_diff", float(head_diff["head_diff"].abs().max())),
        ]
        return pd.DataFrame(stats, columns=["metric", "value"])

    def compare_region_heads(
        self,
        run_a: str | RunRecord,
        run_b: str | RunRecord,
        *,
        region_cells: list[int],
        region_name: str = "region",
        kstpkper: tuple[int, int] | None = None,
    ) -> pd.DataFrame:
        comparison = self.compare_heads(run_a, run_b, kstpkper=kstpkper, cells=region_cells)
        comparison["region"] = region_name
        return comparison

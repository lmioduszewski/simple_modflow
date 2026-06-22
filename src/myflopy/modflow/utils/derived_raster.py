"""Provenance-aware caching for *derived* rasters (interpolated contours, etc.).

A derived raster is produced from one or more source files plus some parameters
(resolution, method, ...). We cache the produced raster on disk next to a small
``.provenance.json`` sidecar recording a content hash of the sources and the
parameters. On resolve:

- **missing** cache  -> produce it now,
- **fresh**   cache  -> reuse it,
- **stale**   cache  -> warn and reuse it (the source or params changed); the
  caller must explicitly ``refresh`` to rebuild.

This keeps the "auto-detect change, rebuild only when I tell you" workflow, and
is independent of how the raster is produced (so it is testable without GRASS).
"""

from __future__ import annotations

import hashlib
import json
import warnings
from pathlib import Path
from typing import Callable, Sequence


def _hash_file(path: Path) -> str | None:
    if not path.exists():
        return None
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


class DerivedRaster:
    """Cache + provenance for a single produced raster file."""

    def __init__(
        self,
        out: Path | str,
        sources: Sequence[Path | str],
        params: dict,
        produce: Callable[[], None],
    ):
        self.out = Path(out)
        self.sources = [Path(s) for s in sources]
        self.params = dict(params)
        self.produce = produce

    @property
    def sidecar(self) -> Path:
        return self.out.with_suffix(self.out.suffix + ".provenance.json")

    def _current_key(self) -> dict:
        key = {
            "sources": {str(s): _hash_file(s) for s in self.sources},
            "params": self.params,
        }
        # Normalize through JSON so it compares equal to the loaded sidecar.
        return json.loads(json.dumps(key, default=str, sort_keys=True))

    def _cached_key(self):
        if not self.sidecar.exists():
            return None
        try:
            return json.loads(self.sidecar.read_text())
        except (OSError, json.JSONDecodeError):
            return None

    def status(self) -> str:
        """``"missing"``, ``"fresh"``, or ``"stale"``."""
        if not self.out.exists():
            return "missing"
        return "fresh" if self._cached_key() == self._current_key() else "stale"

    def _write(self) -> None:
        self.produce()
        self.sidecar.write_text(
            json.dumps(self._current_key(), indent=2, sort_keys=True)
        )

    def ensure(self, *, refresh: bool = False, warn: bool = True) -> Path:
        """Return the cached raster path, producing/refreshing as needed.

        ``refresh=True`` rebuilds unconditionally. A *stale* cache is reused with
        a warning (unless ``warn=False``) -- it is only rebuilt on ``refresh``.
        """
        state = self.status()
        if refresh or state == "missing":
            self._write()
        elif state == "stale" and warn:
            warnings.warn(
                f"Derived raster '{self.out.name}' is stale (its source or "
                "parameters changed); reusing the cached file. Rebuild with "
                "refresh=True.",
                stacklevel=2,
            )
        return self.out

"""Provenance-aware caching for *derived* rasters (interpolated contours, etc.).

A derived raster is produced from one or more source files plus some parameters
(resolution, method, ...). We cache the produced raster on disk next to a small
``.provenance.json`` sidecar recording a content hash of the sources and the
parameters. On resolve:

- **missing** cache  -> produce it now,
- **fresh**   cache  -> reuse it,
- **stale**   cache  -> warn and reuse it (the source or params changed); the
  caller must explicitly ``refresh`` to rebuild -- unless its sidecar lacks a
  parameter the caller declared it ``requires``: such a cache was made by a rule
  that is now known to be WRONG, not an older version of the right answer, so it
  is refused rather than reused (ledger 602).

This keeps the "auto-detect change, rebuild only when I tell you" workflow, and
is independent of how the raster is produced (so it is testable without GRASS).
"""

from __future__ import annotations

import hashlib
import json
import warnings
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path


def _hash_file(path: Path) -> str | None:
    """The SHA-256 hex digest of a file's contents (``None`` if the file is missing)."""

    if not path.exists():
        return None
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


class DerivedRaster:
    """Cache + provenance for a single produced raster file.

    Parameters
    ----------
    out : str or pathlib.Path
        The produced raster; its sidecar is ``<out>.provenance.json``.
    sources : sequence of str or pathlib.Path
        Files whose CONTENT keys the cache (hashed; a missing one hashes ``None``).
    params : dict
        Values that key the cache, compared after a JSON round trip.
    produce : callable
        Called with no arguments to write ``out``.
    requires : mapping of str to str, optional
        ``{params key: message}``. Each key must be a key of ``params``. A STALE
        cache whose sidecar records no such parameter is refused by :meth:`ensure`
        with ``ValueError(message)`` instead of being reused -- for a parameter that
        marks a cache made by a rule since found wrong (ledger 602). A
        cache with no readable sidecar keeps the stale rule. Default None: the
        stale-reuse rule is unchanged.
    """

    def __init__(
        self,
        out: Path | str,
        sources: Sequence[Path | str],
        params: dict,
        produce: Callable[[], None],
        requires: Mapping[str, str] | None = None,
    ):
        """Bind an output raster path to its ``sources``, ``params``, and a ``produce()`` builder."""

        self.out = Path(out)
        self.sources = [Path(s) for s in sources]
        self.params = dict(params)
        self.produce = produce
        self.requires = dict(requires or {})
        unknown = sorted(set(self.requires) - set(self.params))
        if unknown:
            raise ValueError(
                f"requires= names {unknown}, which params does not hold: a requirement "
                "is a parameter the CURRENT key records."
            )

    @property
    def sidecar(self) -> Path:
        """The path of the ``.provenance.json`` sidecar next to the output raster."""

        return self.out.with_suffix(self.out.suffix + ".provenance.json")

    def _current_key(self) -> dict:
        """The current provenance key: source content hashes plus the build parameters."""

        key = {
            "sources": {str(s): _hash_file(s) for s in self.sources},
            "params": self.params,
        }
        # Normalize through JSON so it compares equal to the loaded sidecar.
        return json.loads(json.dumps(key, default=str, sort_keys=True))

    def _cached_key(self):
        """The provenance key stored in the sidecar, or ``None`` if absent/unreadable."""

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
        """Produce the raster and write its current provenance key to the sidecar."""

        self.produce()
        self.sidecar.write_text(
            json.dumps(self._current_key(), indent=2, sort_keys=True)
        )

    def ensure(self, *, refresh: bool = False, warn: bool = True) -> Path:
        """Return the cached raster path, producing/refreshing as needed.

        ``refresh=True`` rebuilds unconditionally. A *stale* cache is reused with
        a warning (unless ``warn=False``) -- it is only rebuilt on ``refresh``.

        Parameters
        ----------
        refresh : bool, default False
            Rebuild even a fresh cache. A missing cache is always built.
        warn : bool, default True
            Warn when a stale cache is reused; ``False`` reuses it silently. It does
            not silence a ``requires=`` refusal.

        Returns
        -------
        pathlib.Path
            ``out``.

        Raises
        ------
        ValueError
            A stale cache whose sidecar lacks a parameter named in ``requires=``,
            with that parameter's message, unless ``refresh=True``.
        """
        state = self.status()
        if not refresh and state == "stale" and self.requires:
            cached = self._cached_key()
            recorded = cached.get("params") if isinstance(cached, dict) else None
            if isinstance(recorded, dict):
                for key, message in self.requires.items():
                    if key not in recorded:
                        raise ValueError(message)
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

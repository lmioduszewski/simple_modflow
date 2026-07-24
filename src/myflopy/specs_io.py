"""YAML (de)serialization for :class:`~myflopy.specs.SimulationSpec` (plan §5.6).

A thin file-format wrapper over ``SimulationSpec.to_dict()`` / ``from_dict()``:
that dict is already fully JSON-safe (builders serialized as importable
references -- including the ``functools.partial`` list-BC builders -- sources
tagged, paths POSIX-encoded), so YAML support is just PyYAML ``safe_dump`` /
``safe_load`` around it. PyYAML is used in safe mode only and imported lazily, so
``import myflopy`` never pulls it in.

The same limits as the dict round-trip apply: a spec whose values are not
JSON-representable (e.g. a computed numpy array baked into a package's options)
cannot be serialized -- that raises in ``to_dict`` before YAML is involved.

TOML is intentionally not supported yet -- TOML has no null type (a ``None`` in
the dict would break the writer) and stdlib ``tomllib`` is 3.11+ while the
package targets ``>=3.10``. Recorded in ``docs/compromises_and_deferrals.md``.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from myflopy.specs import SimulationSpec


def simulation_to_yaml(spec: SimulationSpec) -> str:
    """Serialize a :class:`SimulationSpec` to a YAML string via its JSON-safe dict.

    ``sort_keys=False`` keeps the dict's own field order (``kind``/``name`` first),
    which reads far better than an alphabetized dump.
    """

    import yaml

    return yaml.safe_dump(spec.to_dict(), sort_keys=False, default_flow_style=False)


def simulation_from_yaml(source: str | Path) -> SimulationSpec:
    """Rebuild a :class:`SimulationSpec` from YAML text or a ``.yaml`` file path.

    A :class:`~pathlib.Path` (or a ``str`` naming an existing file) is read from
    disk; any other ``str`` is parsed as YAML text.
    """

    import yaml

    from myflopy.specs import SimulationSpec

    data = yaml.safe_load(_read_source(source))
    if not isinstance(data, dict):
        raise ValueError(
            "YAML does not describe a SimulationSpec (expected a mapping at the top level)."
        )
    return SimulationSpec.from_dict(data)


def _read_source(source: str | Path) -> str:
    """Return YAML text from a Path, an existing-file path string, or raw YAML text."""

    if isinstance(source, Path):
        return source.read_text()
    if isinstance(source, str):
        try:
            candidate = Path(source)
            if candidate.is_file():
                return candidate.read_text()
        except OSError:
            # e.g. an over-long YAML string that cannot be a filename.
            pass
        return source
    raise TypeError(
        f"from_yaml source must be a str or Path, not {type(source).__name__}."
    )

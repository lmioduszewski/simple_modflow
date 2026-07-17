"""myflopy's single deprecation mechanism (implementation plan 3.1).

Policy: ``docs/deprecation_policy.md``. The short version:

- Every deprecated name keeps working but emits a ``DeprecationWarning``
  naming its replacement and the release that deprecated it (``since``).
- Deprecated names are runtime-only (decision D12): they resolve exclusively
  through ``__getattr__`` — module-level (PEP 562) or class-level — and stay
  out of ``__all__``, ``dir()``, ``TYPE_CHECKING`` import blocks, and stubs,
  so IDEs and tab-completion never suggest them.
- A deprecated name lives at least two tagged releases after ``since``
  before it may be removed.
- ``myflopy.__compatibility__`` is the authoritative list of deprecated
  names; ``tests/test_deprecation.py`` cross-checks it against the registry
  the helpers below populate.

Warning messages start with the fully-qualified old name (``myflopy...``),
which lets the test suite escalate internally-triggered deprecations to
errors with a message-prefix filter (see ``[tool.pytest.ini_options]``
``filterwarnings`` in ``pyproject.toml``).
"""

from __future__ import annotations

import sys
import warnings
from collections.abc import Callable, Mapping
from importlib import import_module
from typing import Any

__all__ = [
    "deprecated_instance_getattr",
    "deprecated_module_getattr",
    "registered_deprecations",
    "warn_deprecated",
]

# Fully-qualified old name -> (replacement, since). Populated as deprecating
# modules import; compared against ``myflopy.__compatibility__`` by tests.
_REGISTRY: dict[str, tuple[str, str]] = {}


def registered_deprecations() -> dict[str, tuple[str, str]]:
    """Return a copy of the ``{old-name: (replacement, since)}`` registry."""

    return dict(_REGISTRY)


def warn_deprecated(old: str, new: str, *, since: str, stacklevel: int = 3) -> None:
    """Emit the standard myflopy ``DeprecationWarning`` for ``old``.

    ``old`` must be fully qualified and start with ``myflopy`` (the test
    suite's error filter keys on that message prefix). The default
    ``stacklevel`` of 3 attributes the warning to the caller's caller —
    correct for the one-intermediate-frame patterns here (a ``__getattr__``
    or a resolver function invoking this helper).
    """

    warnings.warn(
        f"{old} is deprecated since myflopy {since} and will be removed in a "
        f"later release (see docs/deprecation_policy.md); use {new} instead.",
        DeprecationWarning,
        stacklevel=stacklevel,
    )


def deprecated_module_getattr(
    mapping: Mapping[str, tuple[str, str, str]],
    module_name: str,
    *,
    fallback_getattr: Callable[[str], Any] | None = None,
    fallback_dir: Callable[[], list[str]] | None = None,
) -> tuple[Callable[[str], Any], Callable[[], list[str]]]:
    """Build a module ``(__getattr__, __dir__)`` pair serving warned aliases.

    ``mapping`` is ``{old_name: (target, replacement, since)}`` where
    ``target`` is ``"package.module:attr"`` for an attribute alias or
    ``"package.module"`` for a module alias. Assign both results at the
    bottom of the deprecating module::

        __getattr__, __dir__ = deprecated_module_getattr(_DEPRECATED, __name__)

    Facades that already have lazy-export machinery pass their existing
    functions via ``fallback_getattr`` / ``fallback_dir``; either way the
    returned ``__dir__`` excludes the deprecated names (D12), so hiding
    comes free at every facade.
    """

    for old, (_, replacement, since) in mapping.items():
        _REGISTRY[f"{module_name}.{old}"] = (replacement, since)

    def __getattr__(name: str) -> Any:
        entry = mapping.get(name)
        if entry is not None:
            target, replacement, since = entry
            warn_deprecated(f"{module_name}.{name}", replacement, since=since)
            target_module, _, target_attr = target.partition(":")
            module = import_module(target_module)
            return getattr(module, target_attr) if target_attr else module
        if fallback_getattr is not None:
            return fallback_getattr(name)
        raise AttributeError(f"module {module_name!r} has no attribute {name!r}")

    def __dir__() -> list[str]:
        if fallback_dir is not None:
            names = set(fallback_dir())
        else:
            names = set(vars(sys.modules[module_name]))
        return sorted(names - set(mapping))

    return __getattr__, __dir__


def deprecated_instance_getattr(
    mapping: Mapping[str, tuple[str, str, str]],
    owner: str,
) -> Callable[[Any, str], Any]:
    """Build a class-level ``__getattr__`` serving warned instance aliases.

    ``mapping`` is ``{old_name: (instance_attr, replacement, since)}``; a
    warned lookup returns ``getattr(self, instance_attr)``. ``owner`` is the
    fully-qualified class name used in warning text and the registry. Assign
    inside the class body::

        __getattr__ = deprecated_instance_getattr({...}, "myflopy.pkg.Cls")

    The deprecated names must not also exist as real class attributes:
    ``__getattr__`` only fires when normal lookup misses, and that miss is
    exactly what keeps them out of ``dir()`` and completion (D12).
    """

    for old, (_, replacement, since) in mapping.items():
        _REGISTRY[f"{owner}.{old}"] = (replacement, since)

    def __getattr__(self: Any, name: str) -> Any:
        entry = mapping.get(name)
        if entry is None:
            raise AttributeError(
                f"{type(self).__name__!r} object has no attribute {name!r}"
            )
        instance_attr, replacement, since = entry
        warn_deprecated(f"{owner}.{name}", replacement, since=since)
        return getattr(self, instance_attr)

    return __getattr__

"""Tests for the optional-dependency import helper and the declared extras."""

from __future__ import annotations

from pathlib import Path

import pytest

tomllib = pytest.importorskip("tomllib", reason="tomllib is stdlib only on Python 3.11+")

from myflopy import _optional

_REPO_ROOT = Path(__file__).resolve().parents[1]


def test_require_returns_installed_module():
    """``require`` imports and returns a present module."""

    import json

    assert _optional.require("json") is json


def test_require_missing_mentions_module_and_feature():
    """A missing module without a mapped extra still names the module + feature."""

    with pytest.raises(ModuleNotFoundError) as excinfo:
        _optional.require("myflopy_no_such_module_zzz", feature="unit testing")

    message = str(excinfo.value)
    assert "myflopy_no_such_module_zzz" in message
    assert "unit testing" in message
    # No extra maps to this fake module, so no pip hint is offered.
    assert "pip install" not in message


def test_require_missing_with_extra_gives_pip_hint(monkeypatch):
    """A missing module that maps to an extra yields a ``myflopy[extra]`` hint."""

    monkeypatch.setitem(
        _optional._EXTRA_FOR_MODULE, "myflopy_no_such_module_zzz", "pest"
    )
    with pytest.raises(ModuleNotFoundError) as excinfo:
        _optional.require("myflopy_no_such_module_zzz")

    assert 'pip install "myflopy[pest]"' in str(excinfo.value)


def test_known_optional_modules_map_to_their_extras():
    """pyemu and xugrid resolve to the extras that provide them."""

    assert _optional._EXTRA_FOR_MODULE["pyemu"] == "pest"
    assert _optional._EXTRA_FOR_MODULE["xugrid"] == "xugrid"


def _extras() -> dict[str, list[str]]:
    data = tomllib.loads((_REPO_ROOT / "pyproject.toml").read_text())
    return data["project"]["optional-dependencies"]


def test_pest_extra_declares_pyemu():
    extras = _extras()
    assert "pest" in extras
    assert any("pyemu" in requirement for requirement in extras["pest"])


def test_xugrid_extra_declares_xugrid():
    extras = _extras()
    assert "xugrid" in extras
    assert any("xugrid" in requirement for requirement in extras["xugrid"])

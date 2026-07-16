"""Executable resolution order for MP3DU (implementation plan 2.4).

Order: explicit arg > MYFLOPY_MP3DU_DIR env var > repo tools/mp3du/ >
deprecated in-package fallback (which warns).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from myflopy.modflow.mp3du import particles
from myflopy.modflow.mp3du.particles import resolve_mp3du_executable


def test_explicit_argument_wins(monkeypatch, tmp_path):
    monkeypatch.setenv("MYFLOPY_MP3DU_DIR", str(tmp_path / "envdir"))
    explicit = tmp_path / "elsewhere" / "mp3du.exe"
    assert resolve_mp3du_executable("mp3du.exe", explicit) == explicit


def test_env_var_beats_repo_dir(monkeypatch, tmp_path):
    env_dir = tmp_path / "envdir"
    monkeypatch.setenv("MYFLOPY_MP3DU_DIR", str(env_dir))
    assert resolve_mp3du_executable("mp3du.exe") == env_dir / "mp3du.exe"


def test_repo_tools_dir_when_present(monkeypatch, tmp_path):
    monkeypatch.delenv("MYFLOPY_MP3DU_DIR", raising=False)
    tools = tmp_path / "tools" / "mp3du"
    tools.mkdir(parents=True)
    (tools / "mp3du.exe").write_bytes(b"")
    monkeypatch.setattr(particles, "_REPO_TOOLS_DIR", tools)
    assert resolve_mp3du_executable("mp3du.exe") == tools / "mp3du.exe"


def test_legacy_package_fallback_warns(monkeypatch, tmp_path):
    monkeypatch.delenv("MYFLOPY_MP3DU_DIR", raising=False)
    monkeypatch.setattr(particles, "_REPO_TOOLS_DIR", tmp_path / "missing")
    legacy_dir = tmp_path / "pkg"
    legacy_dir.mkdir()
    (legacy_dir / "mp3du.exe").write_bytes(b"")
    monkeypatch.setattr(particles, "_MODULE_DIR", legacy_dir)
    with pytest.warns(DeprecationWarning, match="tools/mp3du"):
        resolved = resolve_mp3du_executable("mp3du.exe")
    assert resolved == legacy_dir / "mp3du.exe"


def test_unresolved_returns_repo_candidate(monkeypatch, tmp_path):
    monkeypatch.delenv("MYFLOPY_MP3DU_DIR", raising=False)
    monkeypatch.setattr(particles, "_REPO_TOOLS_DIR", tmp_path / "missing")
    monkeypatch.setattr(particles, "_MODULE_DIR", tmp_path / "also_missing")
    resolved = resolve_mp3du_executable("mp3du.exe")
    assert resolved == tmp_path / "missing" / "mp3du.exe"


def test_repo_tools_dir_points_at_repo_root():
    assert particles._REPO_TOOLS_DIR == Path(__file__).resolve().parents[1] / "tools" / "mp3du"

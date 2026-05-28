from __future__ import annotations

import os
import shutil
from pathlib import Path
from uuid import uuid4

import pytest


_REPO_ROOT = Path(__file__).resolve().parents[1]


def _pytest_temp_root() -> Path:
    configured = os.environ.get("SIMPLE_MODFLOW_PYTEST_TMP_ROOT")
    if configured:
        return Path(configured)
    return _REPO_ROOT / ".pytest-work" / "custom_tmp"


@pytest.fixture
def tmp_path():
    """Repo-local replacement for pytest's builtin tmp_path fixture.

    On this Windows setup, pytest's temp-path cleanup can intermittently leave
    behind folders with broken permissions. This fixture uses a normal temp
    directory we control and ignores cleanup failures instead of failing the
    whole test session.
    """

    root = _pytest_temp_root()
    root.mkdir(parents=True, exist_ok=True)
    path = root / f"case_{uuid4().hex[:10]}"
    path.mkdir(parents=True, exist_ok=False)
    yield path
    try:
        shutil.rmtree(path)
    except OSError:
        pass

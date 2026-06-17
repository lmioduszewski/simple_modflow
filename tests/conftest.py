from __future__ import annotations

import os
import shutil
from pathlib import Path
from uuid import uuid4

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")

from myflopy.modflow.mf6.canonical import CANONICAL_MODEL_CONTRACT
from myflopy.modflow.mf6.canonical_example import (
    CanonicalModelConfig,
    build_canonical_model,
)


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


@pytest.fixture
def canonical_config():
    """Return the scaled validation profile of the authoritative model."""

    return CanonicalModelConfig.validation()


@pytest.fixture
def canonical_model(tmp_path, canonical_config):
    """Build and contract-check the canonical model for integration tests."""

    model = build_canonical_model(tmp_path / "canonical", config=canonical_config)
    CANONICAL_MODEL_CONTRACT.validate(model)
    return model


@pytest.fixture
def canonical_run(canonical_model):
    """Run and return the contract-checked canonical integration model."""

    success, report = canonical_model.run_simulation()
    assert success, "\n".join(report[-30:])
    CANONICAL_MODEL_CONTRACT.validate(canonical_model)
    return canonical_model

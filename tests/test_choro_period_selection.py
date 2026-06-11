from __future__ import annotations

from types import SimpleNamespace

import pytest

from simple_modflow.modflow.utils.datatypes.choros import Choro


def _selector(kstpkper, per_timestep="last"):
    choro = Choro.__new__(Choro)
    values = list(kstpkper)
    choro.model = SimpleNamespace(kstpkper=values, hds=SimpleNamespace(kstpkper=values))
    choro.kstpkper = None
    choro._per = None
    choro.per_timestep = per_timestep
    return choro


def test_choro_period_defaults_to_final_saved_timestep():
    choro = _selector([(0, 0), (1, 0), (0, 1), (1, 1)])

    choro.per = 0

    assert choro.kstpkper == (1, 0)
    assert choro.per == 0


def test_choro_period_supports_first_and_explicit_timestep_selection():
    first = _selector([(0, 0), (2, 0), (4, 0)], per_timestep="first")
    first.per = 0
    assert first.kstpkper == (0, 0)

    exact = _selector([(0, 0), (2, 0), (4, 0)], per_timestep=2)
    exact.per = 0
    assert exact.kstpkper == (2, 0)

    indexed = _selector([(2, 0), (4, 0), (6, 0)], per_timestep=-1)
    indexed.per = 0
    assert indexed.kstpkper == (6, 0)


def test_choro_period_selection_reports_invalid_inputs():
    unavailable = _selector([(0, 0), (1, 0)])
    with pytest.raises(ValueError, match="Available periods"):
        unavailable.per = 4

    bad_selector = _selector([(0, 0), (1, 0)], per_timestep="middle")
    with pytest.raises(ValueError, match="per_timestep"):
        bad_selector.per = 0

"""Names are text, and an empty answer to a wrong question is refused (ledger 476).

``AoiLibrary.setdefault(1, box)`` stored ``'1'`` and then raised; a target set
registered as ``registry[1]`` was ``in`` the registry while ``registry[1]`` raised
``KeyError``; ``ModelConfig.from_model(<FloPy model>)`` and ``config.section(<typo>)``
came back empty.
"""

from __future__ import annotations

import pandas as pd
import pytest

import myflopy as mf
from myflopy.library import AoiLibrary
from myflopy.modflow.mf6.observations.registry import TargetRegistry
from myflopy.project.model_config import ModelConfig


@pytest.mark.parametrize("store", [
    lambda library: library.setdefault(1, (0, 0, 1, 1)),
    lambda library: library.__setitem__(1, (0, 0, 1, 1)),
    lambda library: library.update({1: (0, 0, 1, 1)}),
    lambda library: library._define(1, (0, 0, 1, 1)),
])
def test_an_area_is_named_by_text(store):
    library = AoiLibrary()
    with pytest.raises(TypeError, match="An area is named by text, not int"):
        store(library)
    assert list(library) == []
    assert library.setdefault("1", (0, 0, 1, 1)).geometry.bounds == (0, 0, 1, 1)


def _conc():
    return mf.ConcTargets(
        locations=pd.DataFrame({"name": ["w"], "cell": [3], "layer": [0]}),
        values=pd.DataFrame({"name": ["w"], "time": [0], "conc": [1.0]}),
    )


def test_a_target_set_is_named_by_text():
    registry = TargetRegistry(None)
    with pytest.raises(TypeError, match="A target set is named by text, not int"):
        registry[1] = _conc()
    registry["1"] = _conc()
    assert 1 not in registry and "1" in registry
    with pytest.raises(TypeError, match="named by text"):
        registry[1]
    with pytest.raises(TypeError, match="named by text"):
        del registry[1]
    assert registry["1"] is not None


def test_model_config_refuses_a_flopy_model_and_a_missing_section(canonical_model):
    with pytest.raises(TypeError, match="pass the myflopy model"):
        ModelConfig.from_model(canonical_model.gwf)
    config = canonical_model.config
    assert config.section("tdis")
    with pytest.raises(KeyError, match="no section 'zzz'"):
        config.section("zzz")
    with pytest.raises(KeyError, match="no section 'zzz'"):
        config.settings(section="zzz")
    assert config.section("TDIS") == config.tdis

"""Tests for model.config extraction and the ModelDiff config tier (Phase 2).

The extractor runs on tiny real flopy simulations (fast, no Voronoi). The
config-diff tier is tested on hand-built ``ModelConfig`` tables for precise
control over value differences and present-in-one-only settings.
"""

from __future__ import annotations

import tempfile

import flopy
import pandas as pd
import pytest

from myflopy.project.model_config import ModelConfig
from myflopy.project.model_group import ModelGroup


class _FlopyModelLike:
    """Minimal model exposing the flopy handles the extractor reads."""

    def __init__(self, sim, gwf, ims, name):
        self.sim = sim
        self.gwf = gwf
        self.ims = ims
        self.name = name


#: Where `_tiny_config` puts its simulations: one pytest-managed directory for the
#: module, rather than a bare `tempfile.mkdtemp()` per call that nothing removed.
_WORKSPACE = None


@pytest.fixture(scope="module", autouse=True)
def _workspace(tmp_path_factory):
    global _WORKSPACE
    _WORKSPACE = tmp_path_factory.mktemp("model_config")
    yield
    _WORKSPACE = None


def _tiny_config(
    *, kappa=0.1, outer=1e-2, tsmult=1.4, nstp=5, spec_discharge=True, newtonoptions=None
) -> ModelConfig:
    ws = tempfile.mkdtemp(dir=_WORKSPACE)
    sim = flopy.mf6.MFSimulation(sim_name="t", sim_ws=ws)
    flopy.mf6.ModflowTdis(
        sim, time_units="days", nper=2, perioddata=[(10.0, 1, 1.1), (100.0, nstp, tsmult)]
    )
    ims = flopy.mf6.ModflowIms(
        sim,
        complexity="MODERATE",
        outer_dvclose=outer,
        inner_dvclose=1e-4,
        linear_acceleration="BICGSTAB",
        under_relaxation="DBD",
        under_relaxation_kappa=kappa,
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname="t", newtonoptions=newtonoptions)
    flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=2, ncol=2)
    flopy.mf6.ModflowGwfnpf(
        gwf, save_flows=True, save_specific_discharge=spec_discharge, icelltype=1
    )
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord="t.hds", saverecord=[("HEAD", "ALL")])
    return ModelConfig.from_model(_FlopyModelLike(sim, gwf, ims, "t"))


# --- extractor ---------------------------------------------------------------
def test_extractor_captures_expected_sections_and_values():
    config = _tiny_config()
    assert set(["tdis", "ims", "npf", "oc"]).issubset(set(config.sections))
    assert config.tdis["time_units"] == "days"
    assert config.tdis["nper"] == 2
    # per-period timing is expanded to rows
    assert config.tdis["period[1].nstp"] == 5
    assert config.tdis["period[1].tsmult"] == pytest.approx(1.4)
    # solver block (across options/nonlinear/linear)
    assert config.ims["under_relaxation_kappa"] == pytest.approx(0.1)
    assert config.ims["linear_acceleration"] == "bicgstab"
    assert config.ims["outer_dvclose"] == pytest.approx(1e-2)


def test_the_config_sees_the_models_own_options_block():
    """`NEWTON` lives in the MODEL's name file, which is not in `gwf.packagelist`.

    The extractor walked that list, so the model-level options block was structurally
    invisible and two runs differing ONLY in `NEWTON` compared EQUAL through the config
    tier. That is the option whose absence converted 90,408 cells to DRY on a real
    model and diverged it -- the tool that should have caught it could not see it
    (ledger 286).
    """

    assert _tiny_config(newtonoptions="NEWTON").model
    assert _tiny_config(newtonoptions=None).model == {}


def test_the_config_diff_reports_a_newton_difference():
    """The refusing direction: the blind spot is only closed if a DIFF shows it."""

    with_newton = _tiny_config(newtonoptions="NEWTON").settings(section="model")
    without = _tiny_config(newtonoptions=None).settings(section="model")

    assert len(with_newton) == 1
    assert len(without) == 0
    assert with_newton["setting"].tolist() == ["newtonoptions"]


def test_extractor_excludes_output_filerecords():
    config = _tiny_config()
    oc = config.oc
    assert not any("filerecord" in key for key in oc)
    # but the save record (real config) is captured
    assert "saverecord" in oc


def test_extractor_reads_package_options():
    config = _tiny_config(spec_discharge=True)
    npf = config.section("npf")
    assert npf["save_specific_discharge"] is True


# --- config-diff tier (hand-built ModelConfig for precise control) -----------
def _cfg(rows) -> ModelConfig:
    return ModelConfig(pd.DataFrame(rows, columns=["section", "setting", "value"]))


class _FakeModel:
    def __init__(self, name, config, package_names=()):
        self.name = name
        self.config = config
        self.package_names = list(package_names)


def _group(configs, *, reference):
    models = {name: _FakeModel(name, cfg) for name, cfg in configs.items()}
    return ModelGroup(models, reference=reference)


def test_identical_configs_report_no_config_difference():
    rows = [("ims", "outer_dvclose", 0.01), ("tdis", "nper", 6)]
    diff = _group({"ref": _cfg(rows), "twin": _cfg(list(rows))}, reference="ref").diff()
    assert diff.config.settings().empty
    summary = diff.config.summary()
    assert bool(summary.iloc[0]["identical"])
    assert summary.iloc[0]["settings_differing"] == 0


def test_config_value_and_absent_differences_detected():
    ref = _cfg(
        [
            ("ims", "under_relaxation_kappa", 0.1),
            ("tdis", "period[1].nstp", 5),
            ("ims", "outer_dvclose", 0.01),  # absent in variant
        ]
    )
    variant = _cfg(
        [
            ("ims", "under_relaxation_kappa", 0.3),  # value differs
            ("tdis", "period[1].nstp", 2),  # value differs
        ]
    )
    diff = _group({"ref": ref, "variant": variant}, reference="ref").diff()

    settings = diff.config.settings().set_index("setting")
    assert settings.loc["under_relaxation_kappa", "reference_value"] == 0.1
    assert settings.loc["under_relaxation_kappa", "model_value"] == 0.3
    assert settings.loc["period[1].nstp", "model_value"] == 2
    # present only in the reference -> '<absent>' on the model side
    assert settings.loc["outer_dvclose", "model_value"] == "<absent>"

    assert diff.config.summary().iloc[0]["settings_differing"] == 3
    assert not diff.config.summary().iloc[0]["identical"]


def test_config_section_filter():
    ref = _cfg([("ims", "outer_dvclose", 0.01), ("tdis", "nper", 6)])
    variant = _cfg([("ims", "outer_dvclose", 0.5), ("tdis", "nper", 3)])
    diff = _group({"ref": ref, "variant": variant}, reference="ref").diff()

    ims_only = diff.config.settings(section="ims")
    assert set(ims_only["section"]) == {"ims"}
    assert set(ims_only["setting"]) == {"outer_dvclose"}


def test_config_differences_appear_in_report():
    ref = _cfg([("ims", "under_relaxation_kappa", 0.1)])
    variant = _cfg([("ims", "under_relaxation_kappa", 0.3)])
    report = _group({"ref": ref, "variant": variant}, reference="ref").diff().report()
    assert "Configuration differences" in report
    assert "under_relaxation_kappa" in report
    assert "differs from reference" in report


def test_config_self_reference_diff_errors():
    ref = _cfg([("ims", "outer_dvclose", 0.01)])
    variant = _cfg([("ims", "outer_dvclose", 0.5)])
    diff = _group({"ref": ref, "variant": variant}, reference="ref").diff()
    with pytest.raises(ValueError):
        diff.config.settings(model_name="ref")
    with pytest.raises(KeyError):
        diff.config.settings(model_name="ghost")

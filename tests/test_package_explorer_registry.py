import os

os.environ.setdefault("MPLBACKEND", "Agg")

from myflopy.modflow.mf6.package_explorer import (
    ModelPackages,
    get_default_budget_term,
    get_default_package_colorscale,
    get_default_package_value_column,
    get_package_explorer_spec,
)
from myflopy.modflow.mf6.package_registry import get_package_explorer_spec as registry_spec
from myflopy.project.model_group import GroupPackages


class _DummyModel:
    pass


def test_package_explorer_registry_preserves_existing_defaults():
    assert get_default_package_value_column("rch") == "recharge"
    assert get_default_package_value_column("drn") == "elev"
    assert get_default_package_colorscale("ghb") == "Portland"
    assert get_default_budget_term("uzf_gwrch") == ("UZF-GWRCH", "gwrch")
    assert get_package_explorer_spec("lak") is registry_spec("lak")


def test_package_explorer_registry_adds_wel_semantics():
    spec = get_package_explorer_spec("wel")

    assert spec is not None
    assert spec.kind == "cell_stress"
    assert spec.default_input == "q"
    assert spec.inputs["q"].label == "Well flow"
    assert get_default_budget_term("wel") == ("WEL", "q")


def test_package_explorer_registry_exposes_uzf_period_fields():
    spec = get_package_explorer_spec("uzf")

    assert spec is not None
    assert {"finf", "pet", "extdp", "extwc", "ha", "hroot", "rootact"}.issubset(spec.inputs)
    assert spec.inputs["pet"].label == "Potential evapotranspiration"
    assert get_default_package_colorscale("uzf_pet") == "YlOrRd"


def test_model_packages_exposes_registry_backed_wel_accessors():
    packages = ModelPackages(_DummyModel())

    assert packages.wel.package_name == "wel"
    assert packages.wel.inputs.q.field_name == "q"
    assert packages.wel.results.q.package_name == "wel"
    assert packages.wel.results.q.budget_text == "WEL"
    assert packages.rch.results.fields["field"].tolist() == ["q"]
    assert packages.rch.results.q.value_name == "q"
    assert packages.uzf.inputs.pet.field_name == "pet"
    assert packages.uzf.inputs.rootact.field_name == "rootact"
    assert packages.uzf.inputs.fields["field"].tolist() == [
        "finf",
        "pet",
        "extdp",
        "extwc",
        "ha",
        "hroot",
        "rootact",
    ]
    assert packages.uzf.results.fields["field"].tolist() == ["gwrch", "sat"]
    assert packages.uzf.results.gwrch.value_name == "gwrch"
    assert packages.uzf.results.sat.value_name == "sat"


def test_model_group_exposes_wel_package_accessor_without_full_initialization():
    group = _DummyModel()
    group.wel = "wel-inputs"
    packages = GroupPackages(group)

    assert packages.wel.inputs == "wel-inputs"

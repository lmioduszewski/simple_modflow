import os

os.environ.setdefault("MPLBACKEND", "Agg")

import myflopy.modflow.mf6.package_explorer as package_explorer
import myflopy.modflow.mf6.package_budget as package_budget
import myflopy.modflow.mf6.package_explorer_utils as package_explorer_utils
import myflopy.modflow.mf6.package_inputs as package_inputs
import myflopy.modflow.mf6.package_model as package_model
import myflopy.modflow.mf6.package_plotting as package_plotting
import myflopy.modflow.mf6.package_registry as package_registry
import myflopy.modflow.mf6.package_results as package_results
import myflopy.modflow.mf6.package_surface_water as package_surface_water
import myflopy.modflow.mf6.package_tables as package_tables
from myflopy.modflow.mf6.package_explorer import (
    ModelPackages,
    get_default_budget_term,
    get_default_package_colorscale,
    get_default_package_value_column,
    get_package_explorer_spec,
)
from myflopy.modflow.mf6.package_registry import (
    get_package_explorer_spec as registry_spec,
)
from myflopy.project.model_group import GroupPackages


class _DummyModel:
    pass


def test_package_explorer_facade_preserves_downstream_import_surface():
    names = [
        "_blue_white_red_diverging_colorscale",
        "_normalize_connection_type_filter",
        "build_surface_water_exchange_cell_table",
        "build_surface_water_q_map_payload",
        "build_cell_package_input_table",
        "build_cell_input_map_payload",
        "build_budget_result_table",
        "build_lak_connection_table",
        "build_lak_budget_result_table",
        "build_lak_q_map_payload",
        "build_sfr_budget_result_table",
        "build_sfr_q_map_payload",
        "build_group_input_compare_map_payload",
        "build_uzf_field_input_table",
        "get_default_budget_term",
        "get_default_group_compare_colorscale",
        "get_default_package_colorscale",
        "get_package_input_field_spec",
        "get_default_package_value_column",
        "_symmetric_color_limit",
        "ModelPackages",
    ]

    missing = [name for name in names if not hasattr(package_explorer, name)]

    assert missing == []


def test_package_explorer_facade_reexports_split_module_surfaces():
    modules = [
        package_registry,
        package_explorer_utils,
        package_tables,
        package_budget,
        package_plotting,
        package_inputs,
        package_results,
        package_surface_water,
        package_model,
    ]
    expected = {name for module in modules for name in module.__all__}

    assert expected.issubset(set(package_explorer.__all__))


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
    assert {"finf", "pet", "extdp", "extwc", "ha", "hroot", "rootact"}.issubset(
        spec.inputs
    )
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
    # GroupPackages reads the private input accessor; the public ``group.wel``
    # shortcut is deprecated in favor of ``group.packages.wel.inputs``.
    group._wel = "wel-inputs"
    packages = GroupPackages(group)

    assert packages.wel.inputs == "wel-inputs"

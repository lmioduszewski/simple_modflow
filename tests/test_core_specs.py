from __future__ import annotations

from types import SimpleNamespace

import flopy
import pytest

from myflopy import (
    ExchangeSpec,
    ModelSpec,
    ModelType,
    PackageBuilder,
    PackageSpec,
    SimulationSpec,
    build_gwf_gwt_exchange,
    build_ims,
)


def test_package_builder_contract_builds_from_stored_configuration_only():
    class ExampleBuilder:
        def __init__(self, value):
            self.value = value

        def build(self):
            return PackageSpec("example", dict, {"value": self.value})

    builder: PackageBuilder = ExampleBuilder(2.0)

    assert builder.build().options == {"value": 2.0}
    with pytest.raises(TypeError):
        builder.build(value=3.0)


def test_package_specs_are_replaced_by_name_without_mutating_the_base_model():
    base_npf = PackageSpec("npf", dict, {"k": 1.0})
    high_k_npf = base_npf.with_options(k=10.0)
    model = ModelSpec("flow", ModelType.GWF, packages=(base_npf,))

    variant = model.with_package(high_k_npf)

    assert model.packages[0].options["k"] == 1.0
    assert variant.packages[0].options["k"] == 10.0
    assert variant.package("npf").options["k"] == 10.0
    assert len(variant.packages) == 1

    with pytest.raises(KeyError, match="no package named 'missing'"):
        model.package("missing")


def test_spec_representations_are_not_large_data_dumps():
    package = PackageSpec(
        "drn",
        dict,
        {
            "stress_period_data": {
                period: [[(0, cell), 10.0, 1.0] for cell in range(50)]
                for period in range(3)
            },
            "boundnames": True,
        },
        metadata={"source_type": "geopackage"},
    )
    model = ModelSpec("flow", "gwf", packages=(package,))
    simulation = SimulationSpec("sim", models=(model,), packages=(PackageSpec("tdis", dict),))

    package_text = repr(package)
    model_text = repr(model)
    simulation_text = repr(simulation)

    assert "period_data(periods=3, records=150, max_per_period=50)" in package_text
    assert "(0, 49)" not in package_text
    assert "packages=('drn',)" in model_text
    assert "models=('flow',)" in simulation_text
    assert "PackageSpec: <code>drn</code>" in package._repr_html_()
    assert "period_data(periods=3, records=150, max_per_period=50)" in package._repr_html_()
    assert "ModelSpec: <code>flow</code>" in model._repr_html_()
    assert "SimulationSpec: <code>sim</code>" in simulation._repr_html_()


def test_model_spec_rejects_duplicate_package_names():
    package = PackageSpec("npf", dict)

    with pytest.raises(ValueError, match="Duplicate package names: npf"):
        ModelSpec("flow", "gwf", packages=(package, package))


def test_simulation_spec_builds_models_packages_and_exchanges_in_order():
    events: list[str] = []

    def build_simulation(*, sim_name, **options):
        events.append(f"simulation:{sim_name}")
        return SimpleNamespace(name=sim_name, options=options)

    def build_model(simulation, *, modelname, **options):
        events.append(f"model:{modelname}")
        return SimpleNamespace(name=modelname, simulation=simulation, options=options)

    def build_package(model, *, value):
        events.append(f"package:{model.name}:{value}")
        return SimpleNamespace(model=model, value=value)

    def build_exchange(simulation, models, *, factor):
        events.append(f"exchange:{models[0].name}:{models[1].name}")
        return SimpleNamespace(simulation=simulation, models=models, factor=factor)

    flow = ModelSpec(
        "flow",
        "gwf",
        builder=build_model,
        packages=(PackageSpec("npf", build_package, {"value": 1.0}),),
    )
    transport = ModelSpec("transport", "gwt", builder=build_model)
    exchange = ExchangeSpec(
        "flow_transport",
        build_exchange,
        models=("flow", "transport"),
        options={"factor": 2.0},
    )

    simulation = SimulationSpec(
        "coupled",
        models=(flow, transport),
        exchanges=(exchange,),
        builder=build_simulation,
    )
    built = simulation.build_flopy()

    assert simulation.model("flow") is flow
    assert list(built.models) == ["flow", "transport"]
    assert built.models["flow"].packages["npf"].value == 1.0
    assert built.exchanges["flow_transport"].factor == 2.0
    assert events == [
        "simulation:coupled",
        "model:flow",
        "package:flow:1.0",
        "model:transport",
        "exchange:flow:transport",
    ]


def test_simulation_spec_builds_real_coupled_gwf_and_gwt_models(tmp_path):
    flow = ModelSpec(
        "flow",
        "gwf",
        packages=(
            PackageSpec(
                "dis",
                flopy.mf6.ModflowGwfdis,
                {"nlay": 1, "nrow": 1, "ncol": 2, "delr": 1.0, "delc": 1.0, "top": 10.0, "botm": 0.0},
            ),
            PackageSpec("ic", flopy.mf6.ModflowGwfic, {"strt": 9.0}),
            PackageSpec("npf", flopy.mf6.ModflowGwfnpf, {"k": 1.0}),
        ),
    )
    transport = ModelSpec(
        "transport",
        "gwt",
        packages=(
            PackageSpec(
                "dis",
                flopy.mf6.ModflowGwtdis,
                {"nlay": 1, "nrow": 1, "ncol": 2, "delr": 1.0, "delc": 1.0, "top": 10.0, "botm": 0.0},
            ),
            PackageSpec("ic", flopy.mf6.ModflowGwtic, {"strt": 0.0}),
            PackageSpec("mst", flopy.mf6.ModflowGwtmst, {"porosity": 0.25}),
            PackageSpec("adv", flopy.mf6.ModflowGwtadv),
        ),
    )

    built = SimulationSpec(
        "coupled",
        models=(flow, transport),
        packages=(
            PackageSpec(
                "tdis",
                flopy.mf6.ModflowTdis,
                {"nper": 1, "perioddata": [(1.0, 1, 1.0)]},
            ),
            PackageSpec(
                "flow_solver",
                build_ims,
                {"models": ("flow",), "complexity": "SIMPLE"},
            ),
            PackageSpec(
                "transport_solver",
                build_ims,
                {"models": ("transport",), "complexity": "SIMPLE"},
            ),
        ),
        exchanges=(
            ExchangeSpec(
                "gwf_gwt",
                build_gwf_gwt_exchange,
                models=("flow", "transport"),
            ),
        ),
    ).build_flopy(tmp_path)

    assert isinstance(built.models["flow"].model, flopy.mf6.ModflowGwf)
    assert isinstance(built.models["transport"].model, flopy.mf6.ModflowGwt)
    assert built.sim is built.simulation
    assert built.gwf("flow") is built.models["flow"].gwf
    assert built.gwt("transport") is built.models["transport"].gwt
    assert built.model("flow") is built.gwf("flow")
    assert built.built_model("flow").sim is built.simulation
    assert built.package("tdis").parent is built.simulation
    assert built.models["flow"].package("npf").parent.name == "flow"
    assert built.models["flow"].packages["npf"].parent.name == "flow"
    assert built.exchanges["gwf_gwt"].exgtype == "GWF6-GWT6"
    with pytest.raises(TypeError, match="not a GWT model"):
        built.gwt("flow")

    built.simulation.write_simulation(silent=True)
    assert (tmp_path / "mfsim.nam").exists()
    assert (tmp_path / "flow.nam").exists()
    assert (tmp_path / "transport.nam").exists()


def test_exchange_rejects_unknown_model_names():
    exchange = ExchangeSpec("missing", lambda simulation, models: None, ("flow", "transport"))
    simulation = SimulationSpec("bad", models=(ModelSpec("flow", "gwf"),), exchanges=(exchange,))

    with pytest.raises(KeyError, match="unknown models: transport"):
        simulation.build_flopy()

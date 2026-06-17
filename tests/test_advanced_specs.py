from __future__ import annotations

from types import SimpleNamespace

import flopy
import pytest

from myflopy import (
    ModelContext,
    ModelSpec,
    PackageSpec,
    PostBuildHook,
    SimulationSpec,
    UZFBuilder,
    build_ims,
    ghb_spec,
    lak_spec,
    mvr_spec,
    sfr_spec,
    uzf_spec,
    wel_spec,
)


def test_model_context_and_post_build_hooks_receive_completed_model():
    events = []
    context = ModelContext(grid="grid", dates=["2026-01-01"], metadata={"scenario": "base"})

    def build_model(simulation, *, modelname):
        return SimpleNamespace(name=modelname, simulation=simulation)

    def build_package(model):
        events.append("package")
        return SimpleNamespace(parent=model)

    def inspect_model(model, packages, received_context):
        events.append("hook")
        assert model.myflopy_context is context
        assert list(packages) == ["npf"]
        return received_context.metadata["scenario"]

    model = ModelSpec(
        "flow",
        "gwf",
        builder=build_model,
        packages=(PackageSpec("npf", build_package),),
        context=context,
        hooks=(PostBuildHook("inspect", inspect_model),),
    )
    built = SimulationSpec(
        "context",
        models=(model,),
        builder=lambda *, sim_name: SimpleNamespace(name=sim_name),
    ).build_flopy()

    assert events == ["package", "hook"]
    assert built.models["flow"].context is context
    assert built.models["flow"].hook_results == {"inspect": "base"}


def test_mvr_spec_validates_declared_and_available_packages():
    with pytest.raises(ValueError, match="undeclared packages: lak"):
        mvr_spec([["sfr"]], {0: [["sfr", 0, "lak", 0, "FACTOR", 1.0]]})
    with pytest.raises(ValueError, match="maxmvr must be at least 1"):
        mvr_spec(
            [["sfr"], ["lak"]],
            {0: [["sfr", 0, "lak", 0, "FACTOR", 1.0]]},
            maxmvr=0,
        )

    mover = mvr_spec(
        [["sfr"], ["lak"]],
        {0: [["sfr", 0, "lak", 0, "FACTOR", 1.0]]},
    )
    with pytest.raises(ValueError, match="requires missing packages: lak"):
        ModelSpec("flow", "gwf", packages=(sfr_spec([], [], {}), mover))
    with pytest.raises(ValueError, match="must be declared after required packages: lak, sfr"):
        ModelSpec(
            "flow",
            "gwf",
            packages=(mover, sfr_spec([], [], {}), lak_spec([], [], {})),
        )


def test_uzf_builder_can_prepare_a_spec_without_live_model():
    data = UZFBuilder(
        context=ModelContext(domain=[[1, 1]]),
        nper=2,
        cells=[(0, 0), (0, 1)],
        vks=0.1,
        thtr=0.05,
        thts=0.3,
        thti=0.1,
        finf={0: [0.001, 0.002], 1: [0.003, 0.004]},
    )

    package = data.build()

    assert package.name == "uzf"
    assert package.options["nuzfcells"] == 2
    assert sorted(package.options["perioddata"]) == [0, 1]


def test_advanced_specs_build_and_write_with_flopy_310(tmp_path):
    flow = ModelSpec(
        "advanced",
        "gwf",
        packages=(
            PackageSpec(
                "dis",
                flopy.mf6.ModflowGwfdis,
                {
                    "nlay": 1,
                    "nrow": 1,
                    "ncol": 3,
                    "delr": 1.0,
                    "delc": 1.0,
                    "top": 10.0,
                    "botm": 0.0,
                },
            ),
            PackageSpec("ic", flopy.mf6.ModflowGwfic, {"strt": 9.0}),
            PackageSpec("npf", flopy.mf6.ModflowGwfnpf, {"k": 1.0}),
            PackageSpec("sto", flopy.mf6.ModflowGwfsto, {"steady_state": {0: True}}),
            wel_spec({0: [[(0, 0, 1), -0.1, "supply"]]}),
            ghb_spec({0: [[(0, 0, 2), 9.0, 1.0]]}),
            uzf_spec(
                [[0, (0, 0, 1), 1, -1, 0.0, 0.1, 0.05, 0.3, 0.1, 3.5]],
                {0: [[0, 0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]},
            ),
            lak_spec(
                [[0, 9.0, 1]],
                [[0, 0, (0, 0, 0), "VERTICAL", 1.0e-5, 0.0, 0.0, 0.0, 0.0]],
                {0: [[0, "STATUS", "ACTIVE"]]},
                mover=True,
            ),
            sfr_spec(
                [[0, (0, 0, 2), 1.0, 1.0, 0.001, 9.0, 1.0, 1.0, 0.03, 0, 1.0, 0]],
                [[0]],
                {0: [[0, "STATUS", "ACTIVE"], [0, "INFLOW", 0.1]]},
                mover=True,
            ),
            mvr_spec(
                [["sfr"], ["lak"]],
                {0: [["sfr", 0, "lak", 0, "FACTOR", 1.0]]},
            ),
            PackageSpec(
                "oc",
                flopy.mf6.ModflowGwfoc,
                {"saverecord": [("HEAD", "LAST"), ("BUDGET", "LAST")]},
            ),
        ),
    )
    simulation = SimulationSpec(
        "advanced",
        models=(flow,),
        packages=(
            PackageSpec(
                "tdis",
                flopy.mf6.ModflowTdis,
                {"nper": 1, "perioddata": [(1.0, 1, 1.0)]},
            ),
            PackageSpec("ims", build_ims, {"models": ("advanced",), "complexity": "SIMPLE"}),
        ),
    )

    built = simulation.build_flopy(tmp_path)
    built.simulation.write_simulation(silent=True)

    assert list(built.models["advanced"].packages)[4:10] == [
        "wel",
        "ghb",
        "uzf",
        "lak",
        "sfr",
        "mvr",
    ]
    for suffix in ("wel", "ghb", "uzf", "lak", "sfr", "mvr"):
        assert (tmp_path / f"advanced.{suffix}").exists()

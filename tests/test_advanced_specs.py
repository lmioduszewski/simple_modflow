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
    evt_spec,
    ghb_spec,
    lak_spec,
    mvr_spec,
    riv_spec,
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


def test_advanced_output_filerecords_carry_the_name_and_never_collide():
    """Every uzf/lak/sfr output filerecord includes ``{name}``, so two packages
    of the same type never share an output file.

    Regression: uzf ``budget_filerecord`` used to ignore ``name`` and two UZF
    packages both wrote ``{model}_budget.uzf``; the SFR artifact restore wrote
    ``sfr_budget.sfr`` with no model name, so two restored models collided.
    """

    from myflopy.modflow.mf6.package_registry import advanced_output_filerecords

    for pkg in ("uzf", "lak", "sfr"):
        one = advanced_output_filerecords(pkg, "one", "m")
        two = advanced_output_filerecords(pkg, "two", "m")
        assert not (set(one.values()) & set(two.values())), pkg
        assert all(rec.startswith("m_one_") for rec in one.values()), one
    # the SFR-restore bug: the model name must be present
    assert (
        advanced_output_filerecords("sfr", "sfr", "mdl")["budget_filerecord"]
        == "mdl_sfr_budget.sfr"
    )


def test_advanced_specs_route_filerecords_through_the_shared_helper():
    """Build and artifact-restore share ONE filerecord source, so a restored
    package's output files match the ones the original build wrote."""

    from myflopy.modflow.mf6.package_registry import advanced_output_filerecords

    for pkg, spec in (
        ("uzf", uzf_spec([], {})),
        ("lak", lak_spec([], [], {})),
        ("sfr", sfr_spec([], [], {})),
    ):
        for key, value in advanced_output_filerecords(pkg, pkg).items():
            assert spec.options[key] == value, f"{pkg}.{key}"


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
            riv_spec({0: [[(0, 0, 0), 9.5, 1.0, 8.5]]}),
            evt_spec({0: [[(0, 0, 0), 10.0, 0.001, 2.5]]}),
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

    assert list(built.models["advanced"].packages)[4:12] == [
        "wel",
        "ghb",
        "riv",
        "evt",
        "uzf",
        "lak",
        "sfr",
        "mvr",
    ]
    for suffix in ("wel", "ghb", "riv", "evt", "uzf", "lak", "sfr", "mvr"):
        assert (tmp_path / f"advanced.{suffix}").exists()


@pytest.mark.parametrize("packages", [["sfr", "lak"], ["sfr"], "sfr", [["sfr"], "lak"]])
def test_mvr_bare_package_names_are_refused_naming_the_record_spelling(packages):
    """A bare name was its own name for `requires=` and its first LETTER for the index
    marker (`["l", "s"]`), so a grid move missed the mover's reach numbers (ledger 303);
    FloPy reads the flat list as ONE record and refuses it at build (ledger 381)."""

    perioddata = {0: [["sfr", 0, "sfr", 1, "FACTOR", 1.0]]}
    with pytest.raises(ValueError, match=r"single-item list or tuple per package.*\[\['sfr'\]"):
        mvr_spec(packages, perioddata)


def test_mvr_index_marker_holds_whole_package_names():
    from myflopy.grid_identity import INDEXED_PACKAGES

    for packages in ([["sfr"], ["lak"]], [("sfr",), ("lak",)]):
        spec = mvr_spec(packages, {0: [["sfr", 0, "lak", 0, "FACTOR", 1.0]]})
        assert spec.metadata[INDEXED_PACKAGES] == ["lak", "sfr"]
        assert spec.requires == ("sfr", "lak")


#: A model's mover in the three record layouts ``modelnames=True`` came with: what
#: ``mf.mvr(modelnames=True)`` wrote (package-only records), and MF6's model-name
#: layout ``(mname, pname)`` with two packages and with one.
MODELNAMES_LAYOUTS = {
    "package-only records": ([["sfr"], ["lak"]], {0: [["sfr", 0, "lak", 0, "FACTOR", 1.0]]}),
    "model-name layout": ([["m", "sfr"], ["m", "lak"]], {0: [["m", "sfr", 0, "m", "lak", 0, "FACTOR", 1.0]]}),
    "model-name layout, one package": ([["m", "sfr"]], {0: [["m", "sfr", 0, "m", "sfr", 1, "FACTOR", 1.0]]}),
}


@pytest.mark.parametrize("layout", list(MODELNAMES_LAYOUTS))
def test_mvr_modelnames_is_refused_naming_mf6s_rule(layout):
    """The model-name layout was read as the single-model one -- its MODEL name became
    the package ``requires=`` and the index marker named, and the refusal was
    'unique package names' or 'undeclared packages: 0' -- and package-only records
    built a file MF6 refuses (ledger 413)."""

    packages, perioddata = MODELNAMES_LAYOUTS[layout]
    with pytest.raises(ValueError, match=r"modelnames=True: MF6 takes MODELNAMES only on the mover of a GWF-GWF exchange"):
        mvr_spec(packages, perioddata, modelnames=True)


def test_mf_mvr_takes_no_modelnames_and_refuses_one_given():
    import myflopy as mf
    from myflopy.grid_identity import INDEXED_PACKAGES

    moves = (mf.Move(mf.MoverConnection("sfr", 3), mf.MoverConnection("lak", 0)),)
    with pytest.raises(ValueError, match=r"modelnames=True: MF6 takes MODELNAMES only"):
        mf.mvr(nper=1, moves=moves, modelnames=True)
    with pytest.raises(TypeError, match="modelnames"):
        mf.MVRBuilder(nper=1, moves=moves, modelnames=True)
    spec = mf.mvr(nper=1, moves=moves, modelnames=False)  # MF6's default, passed on
    assert spec.requires == ("sfr", "lak") and spec.metadata[INDEXED_PACKAGES] == ["lak", "sfr"]


def test_mf6_refuses_modelnames_on_a_models_mover_in_either_layout(tmp_path):
    """What the refusal describes: MF6 stops reading a model's MVR file that declares
    MODELNAMES, whatever its records' layout (``docs/reproducers/mvr_modelnames.py``)."""

    import subprocess
    import warnings
    from pathlib import Path

    import numpy as np

    import myflopy as mf
    from tests.test_sfr_stream_inputs import _sfr_builder

    builder = _sfr_builder(mover=True)
    grid = builder.grid
    props = grid.get_disv_gridprops()
    xy = np.array([(c[1], c[2]) for c in props["cell2d"]])
    edge = np.flatnonzero((xy[:, 0] < 100.0) | (xy[:, 0] > 500.0))
    model = mf.gwf("m", grid=grid, packages=[
        mf.disv(nlay=1, ncpl=grid.ncpl, nvert=len(props["vertices"]), vertices=props["vertices"],
                cell2d=props["cell2d"], top=110.0, botm=[0.0]),
        mf.npf(k=10.0, icelltype=0), mf.ic(strt=90.0),
        mf.chd(stress_period_data={0: [[(0, int(c)), 90.0] for c in edge]}),
        builder.build(),
        mf.mvr(nper=1, moves=(mf.Move(mf.MoverConnection("sfr", 3), mf.MoverConnection("sfr", 4), value=0.5),)),
        mf.oc(head_filerecord="m.hds", saverecord=[("HEAD", "ALL")]),
    ])
    simulation = mf.simulation(model, name="s", packages=[mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)])])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = mf.Project(tmp_path / "p", name="p").run("r", simulation)
    assert run.success
    workspace = Path(run.model().workspace)
    written = (workspace / "m.mvr").read_text()
    declared = written.replace("BEGIN options\n", "BEGIN options\n  MODELNAMES\n")
    model_names = declared.replace("\n  sfr\n", "\n  m  sfr\n").replace("  sfr  4  sfr  5", "  m  sfr  4  m  sfr  5")
    assert model_names.count("  m  sfr") == 3, model_names
    for text in (declared, model_names):
        (workspace / "m.mvr").write_text(text)
        result = subprocess.run(["mf6"], cwd=workspace, capture_output=True, text=True, check=False)
        assert result.returncode != 0
        listing = " ".join((workspace / "mfsim.lst").read_text().split())
        assert "MODELNAMES cannot be specified unless the mover package is for an exchange" in listing

from __future__ import annotations

import flopy
import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon

import myflopy as mf
from myflopy.modflow.mf6.grid.voronoi import VoronoiGridPlus


def _two_cell_grid() -> VoronoiGridPlus:
    verts = np.array(
        [
            [0, 0],
            [1, 0],
            [2, 0],
            [0, 1],
            [1, 1],
            [2, 1],
        ],
        dtype=float,
    )
    iverts = [[0, 1, 4, 3], [1, 2, 5, 4]]
    xcyc = np.array([[0.5, 0.5], [1.5, 0.5]], dtype=float)
    grid = VoronoiGridPlus(verts=verts, iverts=iverts, xcyc=xcyc)
    grid.gdf_topbtm = gpd.GeoDataFrame(
        {0: [10.0, 9.0], 1: [0.0, 0.0]},
        geometry=grid.gdf_vorPolys.geometry,
        crs=grid.crs,
    )
    return grid


def test_package_api_builds_core_flopy_packages(tmp_path):
    flow = mf.gwf(
        "flow",
        packages=(
            mf.disv(
                nlay=1,
                ncpl=2,
                nvert=6,
                vertices=[
                    [0, 0.0, 0.0],
                    [1, 1.0, 0.0],
                    [2, 1.0, 1.0],
                    [3, 0.0, 1.0],
                    [4, 2.0, 0.0],
                    [5, 2.0, 1.0],
                ],
                cell2d=[
                    [0, 0.5, 0.5, 4, 0, 1, 2, 3],
                    [1, 1.5, 0.5, 4, 1, 4, 5, 2],
                ],
                top=[10.0, 9.0],
                botm=[[0.0, 0.0]],
            ),
            mf.ic(strt=[9.0, 8.0]),
            mf.npf(k=1.0),
            mf.oc(saverecord=[("HEAD", "ALL")]),
        ),
        model_nam_file="flow.nam",
        newtonoptions="under_relaxation",
        save_flows=True,
    )
    simulation = mf.SimulationSpec(
        "package_api",
        models=(flow,),
        packages=(
            mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
            mf.ims(models=("flow",), print_option="SUMMARY"),
        ),
    )

    built = simulation.build_flopy(tmp_path)
    built.simulation.write_simulation(silent=True)

    assert flow.options["newtonoptions"] == "under_relaxation"
    assert flow.options["save_flows"] is True
    assert isinstance(built.gwf("flow"), flopy.mf6.ModflowGwf)
    assert (tmp_path / "flow.disv").exists()
    assert "ims" in built.packages


def test_package_api_ims_can_be_shared_or_split_by_model():
    shared = mf.ims(
        models=("flow", "transport"),
        print_option="SUMMARY",
        outer_maximum=100,
        inner_maximum=50,
        linear_acceleration="BICGSTAB",
    )
    flow_solver = mf.ims(name="flow_solver", models=("flow",), complexity="SIMPLE")
    transport_solver = mf.ims(name="transport_solver", models=("transport",), complexity="COMPLEX")

    assert shared.options["models"] == ("flow", "transport")
    assert shared.options["outer_maximum"] == 100
    assert shared.options["linear_acceleration"] == "BICGSTAB"
    assert flow_solver.name == "flow_solver"
    assert flow_solver.options["pname"] == "flow_solver"
    assert transport_solver.options["complexity"] == "COMPLEX"


def test_package_api_builds_typed_model_specs_for_all_mf6_model_types(tmp_path):
    flow = mf.gwf("flow", newtonoptions="under_relaxation", save_flows=True)
    transport = mf.gwt(
        "transport",
        dependent_variable_scaling=True,
        model_nam_file="transport.nam",
    )
    energy = mf.gwe(
        "energy",
        dependent_variable_scaling=True,
        print_flows=True,
    )
    particles = mf.prt("particles", model_rel_path="prt", print_input=True)

    built = mf.SimulationSpec(
        "models",
        models=(flow, transport, energy, particles),
    ).build_flopy(tmp_path)

    assert flow.model_type == mf.ModelType.GWF
    assert transport.model_type == mf.ModelType.GWT
    assert energy.model_type == mf.ModelType.GWE
    assert particles.model_type == mf.ModelType.PRT
    assert flow.options["newtonoptions"] == "under_relaxation"
    assert transport.options["dependent_variable_scaling"] is True
    assert energy.options["dependent_variable_scaling"] is True
    assert particles.options["model_rel_path"] == "prt"
    assert isinstance(built.gwf("flow"), flopy.mf6.ModflowGwf)
    assert isinstance(built.gwt("transport"), flopy.mf6.ModflowGwt)
    assert isinstance(built.gwe("energy"), flopy.mf6.ModflowGwe)
    assert isinstance(built.prt("particles"), flopy.mf6.ModflowPrt)


def test_package_api_exposes_direct_and_geopackage_boundary_paths(tmp_path):
    grid = _two_cell_grid()
    context = mf.ModelContext(grid=grid, domain=np.array([[1, 1]]))
    gpkg = tmp_path / "drn.gpkg"
    gpd.GeoDataFrame(
        {
            "name": ["drain"],
            "layer": [1],
            "elevation": [8.5],
            "conductance": [25.0],
        },
        geometry=[Polygon([(0, 0), (2, 0), (2, 1), (0, 1)])],
        crs=grid.crs,
    ).to_file(gpkg, driver="GPKG")

    direct = mf.drn(stress_period_data={0: [[(0, 0), 8.5, 25.0]]})
    from_gpkg = mf.drn.gpkg(gpkg, context=context, nper=1)

    assert direct.name == "drn"
    assert from_gpkg.metadata["source_type"] == "geopackage"
    assert len(from_gpkg.options["stress_period_data"][0]) == 2


def test_package_api_rch_uses_builder_by_default_and_flopy_for_direct_data():
    context = mf.ModelContext(domain=np.array([[1, 1]]))

    built = mf.rch(context=context, nper=1, recharge=1.0e-4)
    direct = mf.rch.flopy(stress_period_data={0: [[(0, 0), 1.0e-4]]})

    assert built.metadata["builder"] == "RCHBuilder"
    assert built.options["stress_period_data"][0] == [[(0, 0), 1.0e-4], [(0, 1), 1.0e-4]]
    assert direct.metadata == {}


def test_package_api_advanced_helpers_and_flopy_escape_hatches():
    context = mf.ModelContext(domain=np.array([[1, 1]]))

    uzf = mf.uzf(
        context=context,
        nper=1,
        vks=0.1,
        thtr=0.05,
        thts=0.30,
        thti=0.15,
        finf=0.001,
    )
    uzf_direct = mf.uzf.flopy(packagedata=[[0, (0, 0), 1, -1, 0.001, 0.1, 0.05, 0.30, 0.15, 4.0]], perioddata={0: [[0, 0.001]]})
    mvr = mf.mvr(
        nper=1,
        moves=(mf.Move(mf.MoverConnection("sfr", 0), mf.MoverConnection("lak", 0)),),
    )
    mvr_direct = mf.mvr.flopy(
        packages=[["sfr"], ["lak"]],
        perioddata={0: [["sfr", 0, "lak", 0, "FACTOR", 1.0]]},
    )
    sfr_direct = mf.sfr.flopy(packagedata=[], connectiondata=[], perioddata={0: []})
    lak_direct = mf.lak.flopy(packagedata=[], connectiondata=[], perioddata={0: []})

    assert uzf.metadata["builder"] == "UZFBuilder"
    assert uzf_direct.name == "uzf"
    assert mvr.metadata["builder"] == "MVRBuilder"
    assert mvr_direct.requires == ("sfr", "lak")
    assert sfr_direct.name == "sfr"
    assert lak_direct.name == "lak"


def test_simulation_single_model_defaults():
    flow = mf.gwf("flow", packages=[])
    sim = mf.simulation(flow)

    assert [m.name for m in sim.models] == ["flow"]
    names = [p.name for p in sim.packages]
    assert names == ["tdis", "flow_ims"]            # steady tdis + one IMS for flow
    ims = next(p for p in sim.packages if p.name == "flow_ims")
    assert tuple(ims.options["models"]) == ("flow",)
    assert ims.options["complexity"] == "SIMPLE"


def test_simulation_one_ims_per_model_for_coupled():
    flow = mf.gwf("flow", packages=[])
    transport = mf.gwt("transport", packages=[])
    sim = mf.simulation(flow, transport)

    names = [p.name for p in sim.packages]
    assert names == ["tdis", "flow_ims", "transport_ims"]
    by_name = {p.name: p for p in sim.packages}
    assert tuple(by_name["flow_ims"].options["models"]) == ("flow",)
    assert tuple(by_name["transport_ims"].options["models"]) == ("transport",)


def test_simulation_overrides_tdis_and_solver():
    flow = mf.gwf("flow", packages=[])
    my_tdis = mf.tdis(nper=3, perioddata=[(1.0, 1, 1.0)] * 3)
    my_ims = mf.ims(models=["flow"], complexity="COMPLEX", name="custom_ims")
    sim = mf.simulation(flow, tdis=my_tdis, solver=my_ims)

    tdis_pkg = next(p for p in sim.packages if p.name == "tdis")
    assert tdis_pkg.options["nper"] == 3
    assert [p.name for p in sim.packages] == ["tdis", "custom_ims"]


def test_simulation_requires_a_model():
    import pytest

    with pytest.raises(ValueError):
        mf.simulation()


def test_drn_helper_runs_grid_only():
    # F3: a boundary helper given only a grid (no model) must keep the grid and
    # not null it out when there is no model.nper.
    from myflopy.modflow.mf6.drn import DRNFromVector

    vor = _two_cell_grid()
    helper = DRNFromVector(vor=vor)
    assert helper.vor is vor
    assert helper.nper is None

    data = helper.get_drn_stress_period_data(
        cells=[0, 1], conductance=10, bottom_addition=5, layer=0
    )
    assert len(data) == 2
    assert data[0] == [(0, 0), 5, 10]


def test_rch_accepts_stress_period_data():
    # F4: mf.rch should take stress_period_data= like mf.drn / mf.wel.
    spec = mf.rch(stress_period_data={0: [[(0, 0), 1.0e-3]]})
    assert spec.name == "rch"

    # The builder form still requires its arguments.
    import pytest

    with pytest.raises(TypeError):
        mf.rch()


def test_layerstack_build_attach_publishes_gdf_topbtm():
    from myflopy.layers import Array, LayerStack

    vor = _two_cell_grid()
    vor.gdf_topbtm = None  # build(attach=True) should repopulate it
    result = (
        LayerStack(vor, top=Array([10.0, 9.0]))
        .add("upper", thickness=4.0)
        .add("lower", thickness=6.0)
        .build(attach=True)
    )
    # Integer columns 0=top, 1..nlay=layer bottoms, matching the builder format.
    assert vor.gdf_topbtm is not None
    for col in (0, 1, 2):
        assert col in vor.gdf_topbtm.columns
    assert np.allclose(vor.gdf_topbtm[0].to_numpy(), result.top)
    assert np.allclose(vor.gdf_topbtm[2].to_numpy(), result.botm[1])


def test_semantic_mover_connections_resolve_reaches_and_lakes():
    import pytest

    from myflopy.package_api import lak_connection, sfr_connection
    from myflopy.specs import PackageSpec

    sfr_spec = PackageSpec(
        "sfr", flopy.mf6.ModflowGwfsfr,
        metadata={"sfr_index": {
            "package": "sfr",
            "outlets": {"main": 7}, "heads": {"main": 2},
            "by_stream": {"main": [2, 3, 7]},
            "centroids": {2: [0.0, 0.0], 3: [5.0, 0.0], 7: [10.0, 0.0]},
        }},
    )
    assert sfr_connection(sfr_spec, "main").index == 7              # outlet (default)
    assert sfr_connection(sfr_spec, "main", at="upstream").index == 2
    assert sfr_connection(sfr_spec, "main", at=(4.0, 0.0)).index == 3   # nearest reach
    assert sfr_connection(sfr_spec, "main").package == "sfr"

    lak_spec = PackageSpec(
        "lak", flopy.mf6.ModflowGwflak,
        metadata={"lak_index": {"package": "lak", "lakes": {"deep": 0, "shallow": 1}}},
    )
    assert lak_connection(lak_spec, "shallow").index == 1

    with pytest.raises(KeyError):
        sfr_connection(sfr_spec, "nope")
    with pytest.raises(ValueError):
        sfr_connection(PackageSpec("sfr", flopy.mf6.ModflowGwfsfr), "main")
    with pytest.raises(ValueError):
        lak_connection(PackageSpec("lak", flopy.mf6.ModflowGwflak), "deep")

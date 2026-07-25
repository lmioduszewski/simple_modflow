from __future__ import annotations

import flopy
import geopandas as gpd
import numpy as np
import pytest
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


def _disv_values() -> dict:
    """Minimal two-cell DISV grid dict, model-kind-agnostic.

    Vertices are listed CLOCKWISE per cell -- MF6 rejects counter-clockwise
    winding with "Calculated CELL2D area less than zero", so this grid both
    builds AND runs.
    """

    return dict(
        nlay=1,
        ncpl=2,
        nvert=6,
        vertices=[[0, 0.0, 0.0], [1, 1.0, 0.0], [2, 1.0, 1.0], [3, 0.0, 1.0], [4, 2.0, 0.0], [5, 2.0, 1.0]],
        cell2d=[[0, 0.5, 0.5, 4, 0, 3, 2, 1], [1, 1.5, 0.5, 4, 1, 2, 5, 4]],
        top=10.0,
        botm=0.0,
    )


def _dis_values() -> dict:
    """Minimal structured DIS grid dict, model-kind-agnostic (5.5)."""

    return dict(nlay=1, nrow=2, ncol=2, delr=100.0, delc=100.0, top=10.0, botm=0.0)


def test_core_helpers_dispatch_on_model_type(tmp_path):
    """mf.disv / mf.ic / mf.oc resolve the FloPy class from the model kind (5.3A).

    They stored ``ModflowGwf*`` before; now a module-level dispatch builder picks
    the gwf/gwt/gwe class off the built model's ``model_type``.
    """

    expected = {
        mf.gwf: ("gwf", flopy.mf6.ModflowGwfic, flopy.mf6.ModflowGwfdisv, flopy.mf6.ModflowGwfoc),
        mf.gwt: ("gwt", flopy.mf6.ModflowGwtic, flopy.mf6.ModflowGwtdisv, flopy.mf6.ModflowGwtoc),
        mf.gwe: ("gwe", flopy.mf6.ModflowGweic, flopy.mf6.ModflowGwedisv, flopy.mf6.ModflowGweoc),
    }
    for helper, (kind, ic_cls, disv_cls, oc_cls) in expected.items():
        model = helper(kind, packages=[
            mf.disv(**_disv_values()),
            mf.ic(strt=0.0),
            mf.oc(budget_filerecord=f"{kind}.cbc", saverecord=[("BUDGET", "ALL")]),
        ])
        built = mf.SimulationSpec(
            f"s_{kind}",
            models=(model,),
            packages=(mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]), mf.ims(models=(kind,))),
        ).build_flopy(tmp_path / kind)
        flopy_model = built.simulation.get_model(kind)
        assert isinstance(flopy_model.get_package("ic"), ic_cls)
        assert isinstance(flopy_model.get_package("disv"), disv_cls)
        assert isinstance(flopy_model.get_package("oc"), oc_cls)


def test_core_helper_dispatch_builders_round_trip_and_reject_prt():
    """The dispatch builders serialize by importable ref, and gate PRT correctly."""

    from types import SimpleNamespace

    from myflopy.builders import build_disu, build_ic

    assert mf.ic(strt=1.0).to_dict()["builder"] == "myflopy.builders:build_ic"
    assert mf.oc(budget_filerecord="m.cbc").to_dict()["builder"] == "myflopy.builders:build_oc"
    assert mf.disv(**_disv_values()).to_dict()["builder"] == "myflopy.builders:build_disv"
    assert mf.dis(**_dis_values()).to_dict()["builder"] == "myflopy.builders:build_dis"
    assert (
        mf.disu(nodes=2, nja=6, top=[10, 10], bot=[0, 0], area=[1, 1], iac=[2, 2], ja=[0, 1, 1, 0])
        .to_dict()["builder"]
        == "myflopy.builders:build_disu"
    )
    # round-trip resolves back to the same module-level function
    assert mf.PackageSpec.from_dict(mf.ic(strt=1.0).to_dict()).builder is build_ic

    # PRT is dispatched for dis/disv/oc (6.3A) but MF6 has no ModflowPrtic (PRT
    # needs no initial condition) and no ModflowPrtdisu; those must fail loudly,
    # not KeyError.
    prt_like = SimpleNamespace(model_type="prt6")
    for build in (build_ic, build_disu):
        with pytest.raises(ValueError, match="not available for model type 'prt6'"):
            build(prt_like)
    from myflopy.builders import _DIS_CLASSES, _DISV_CLASSES, _OC_CLASSES

    assert _OC_CLASSES["prt6"] is flopy.mf6.ModflowPrtoc
    assert _DISV_CLASSES["prt6"] is flopy.mf6.ModflowPrtdisv
    assert _DIS_CLASSES["prt6"] is flopy.mf6.ModflowPrtdis


def test_dis_dispatch_builds_structured_class_per_kind(tmp_path):
    """mf.dis resolves the GWF/GWT/GWE structured DIS class off the model kind (5.5)."""

    expected = {
        mf.gwf: ("gwf", flopy.mf6.ModflowGwfdis),
        mf.gwt: ("gwt", flopy.mf6.ModflowGwtdis),
        mf.gwe: ("gwe", flopy.mf6.ModflowGwedis),
    }
    for helper, (kind, dis_cls) in expected.items():
        model = helper(kind, packages=[
            mf.dis(**_dis_values()),
            mf.ic(strt=5.0),
            mf.oc(budget_filerecord=f"{kind}.cbc", saverecord=[("BUDGET", "ALL")]),
        ])
        built = mf.SimulationSpec(
            f"d_{kind}",
            models=(model,),
            packages=(mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]), mf.ims(models=(kind,))),
        ).build_flopy(tmp_path / kind)
        flopy_model = built.simulation.get_model(kind)
        assert isinstance(flopy_model.get_package("dis"), dis_cls)


def test_gwt_gwe_package_factories_build_real_models(tmp_path):
    """The 11 GWT/GWE package factories (5.3B) build the right FloPy packages."""

    disv = _disv_values()
    transport = mf.gwt("trans", packages=[
        mf.disv(**disv), mf.ic(strt=0.0),
        mf.adv(scheme="TVD"),
        mf.dsp(alh=1.0, ath1=0.1),
        mf.mst(porosity=0.25),
        mf.ist(porosity=0.05, volfrac=0.2, zetaim=1.0e-3),
        mf.cnc(stress_period_data={0: [[(0, 0), 100.0]]}),
        mf.src(stress_period_data={0: [[(0, 1), 0.5]]}),
        mf.oc(concentration_filerecord="trans.ucn", saverecord=[("CONCENTRATION", "ALL")]),
    ])
    energy = mf.gwe("energy", packages=[
        mf.disv(**disv), mf.ic(strt=10.0),
        mf.est(porosity=0.25, heat_capacity_solid=800.0, density_solid=2650.0),
        mf.cnd(ktw=0.6, kts=3.0),
        mf.ctp(stress_period_data={0: [[(0, 0), 25.0]]}),
        mf.esl(stress_period_data={0: [[(0, 1), 100.0]]}),
        mf.oc(temperature_filerecord="energy.ucn", saverecord=[("TEMPERATURE", "ALL")]),
    ])
    built = mf.SimulationSpec(
        "coupled_transport",
        models=(transport, energy),
        packages=(mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]), mf.ims(models=("trans", "energy"))),
    ).build_flopy(tmp_path)

    tm = built.simulation.get_model("trans")
    em = built.simulation.get_model("energy")
    assert isinstance(tm.get_package("adv"), flopy.mf6.ModflowGwtadv)
    assert isinstance(tm.get_package("dsp"), flopy.mf6.ModflowGwtdsp)
    assert isinstance(tm.get_package("mst"), flopy.mf6.ModflowGwtmst)
    assert isinstance(tm.get_package("ist"), flopy.mf6.ModflowGwtist)
    assert isinstance(tm.get_package("cnc"), flopy.mf6.ModflowGwtcnc)
    assert isinstance(tm.get_package("src"), flopy.mf6.ModflowGwtsrc)
    assert isinstance(em.get_package("est"), flopy.mf6.ModflowGweest)
    assert isinstance(em.get_package("cnd"), flopy.mf6.ModflowGwecnd)
    assert isinstance(em.get_package("ctp"), flopy.mf6.ModflowGwectp)
    assert isinstance(em.get_package("esl"), flopy.mf6.ModflowGweesl)
    # ssm references flow-model source packages by name; check the spec + its builder ref.
    ssm = mf.ssm(sources=[["chd", "AUX", "concentration"]])
    assert ssm.to_dict()["builder"] == "flopy.mf6.modflow.mfgwtssm:ModflowGwtssm"


def test_prt_model_fully_declarable_and_runs(tmp_path):
    """A coupled GWF+PRT simulation is declarable spec-first and runs MF6 (6.3A).

    mf.mip/mf.prp + prt6-dispatched mf.disv/mf.oc + mf.ems, coupled through
    build_gwf_prt_exchange with NO FMI package (flows pass through the exchange).
    mf.simulation()'s solver default gives the PRT model an EMS -- MF6 rejects
    explicit models under IMS6, so this is load-bearing, not cosmetic. Release
    points carry boundnames, which MF6 echoes (uppercased) into the track CSV's
    ``name`` column -- the group key the capture map reads.
    """

    import pandas as pd

    disv = _disv_values()
    flow = mf.gwf("flow", packages=[
        mf.disv(**disv), mf.ic(strt=9.0), mf.npf(k=1.0, save_specific_discharge=True),
        mf.chd(stress_period_data={0: [[(0, 0), 10.0], [(0, 1), 8.0]]}),
        mf.oc(head_filerecord="flow.hds", budget_filerecord="flow.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])
    particles = mf.prt("particles", packages=[
        mf.disv(**disv),
        mf.mip(porosity=0.25),
        mf.prp(packagedata=[(0, (0, 0), 0.35, 0.5, 9.0, "west_wells"),
                            (1, (0, 0), 0.45, 0.5, 9.0, "east_wells")]),
        mf.oc(trackcsv_filerecord="particles.trk.csv",
              budget_filerecord="particles.bud", saverecord=[("BUDGET", "ALL")]),
    ])
    sim = mf.simulation(
        flow, particles,
        exchanges=[mf.ExchangeSpec("gwfprt", mf.build_gwf_prt_exchange,
                                   models=("flow", "particles"))],
    )

    # the solver default is kind-aware: IMS for the GWF model, EMS for PRT
    solver_names = [p.name for p in sim.packages if p.name.endswith(("_ims", "_ems"))]
    assert solver_names == ["flow_ims", "particles_ems"]

    built = sim.build_flopy(tmp_path)
    pm = built.simulation.get_model("particles")
    assert isinstance(pm.get_package("mip"), flopy.mf6.ModflowPrtmip)
    assert isinstance(pm.get_package("prp"), flopy.mf6.ModflowPrtprp)
    assert isinstance(pm.get_package("disv"), flopy.mf6.ModflowPrtdisv)
    assert isinstance(pm.get_package("oc"), flopy.mf6.ModflowPrtoc)

    built.simulation.write_simulation()
    success, report = built.simulation.run_simulation(silent=True)
    assert success, "MF6 GWF+PRT did not converge:\n" + "\n".join(report[-20:])

    (track_csv,) = tmp_path.rglob("*.trk.csv")
    track = pd.read_csv(track_csv)
    assert len(track) > 0
    assert sorted(track["name"].dropna().unique()) == ["EAST_WELLS", "WEST_WELLS"]
    assert (track["ireason"] == 3).any()  # every particle terminates


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
            "stage": [9.5],
            "rbot": [7.5],
            "surface": [10.0],
            "et_rate": [0.002],
            "depth": [2.5],
        },
        geometry=[Polygon([(0, 0), (2, 0), (2, 1), (0, 1)])],
        crs=grid.crs,
    ).to_file(gpkg, driver="GPKG")

    direct = mf.drn(stress_period_data={0: [[(0, 0), 8.5, 25.0]]})
    from_gpkg = mf.drn.gpkg(gpkg, context=context, nper=1)

    assert direct.name == "drn"
    assert from_gpkg.metadata["source_type"] == "geopackage"
    assert len(from_gpkg.options["stress_period_data"][0]) == 2

    riv_direct = mf.riv(stress_period_data={0: [[(0, 0), 9.5, 25.0, 7.5]]})
    riv_gpkg = mf.riv.gpkg(gpkg, context=context, nper=1)

    assert riv_direct.name == "riv"
    assert riv_gpkg.metadata["source_type"] == "geopackage"
    assert riv_gpkg.options["stress_period_data"][0][0] == [(0, 0), 9.5, 25.0, 7.5, "drain"]

    # mf.evt.gpkg: the helper's ~10 keyword forwards are exercised here with
    # DISTINCT values per field, so a swapped/dropped forward (rate=depth,
    # missing surface=) changes the asserted record instead of passing silently.
    evt_direct = mf.evt(stress_period_data={0: [[(0, 0), 10.0, 0.002, 2.5]]})
    evt_gpkg = mf.evt.gpkg(gpkg, context=context, nper=1, rate="et_rate")

    assert evt_direct.name == "evt"
    assert evt_gpkg.metadata["source_type"] == "geopackage"
    assert evt_gpkg.options["stress_period_data"][0][0] == [(0, 0), 10.0, 0.002, 2.5, "drain"]

    # segmented ET cannot come from feature mapping (no pxdp/petm per record):
    # rejected up front rather than failing deep inside FloPy at build time
    with pytest.raises(ValueError, match="nseg=1 only"):
        mf.evt.gpkg(gpkg, context=context, nper=1, rate="et_rate", nseg=2)


def test_list_bc_specs_accept_every_native_flopy_input_shape():
    # The *_spec factories promise native FloPy inputs are taken "verbatim".
    # rch_spec and evt_spec used to pre-compute maxbound with
    # `max(len(r) for r in spd.values())`, which assumes a dict of sized lists --
    # so a bare list raised AttributeError and a None period raised TypeError,
    # while every sibling accepted both. FloPy computes MAXBOUND at write time,
    # so the inference was unnecessary as well as wrong. Cover all list BCs so
    # neither reintroduces it.
    shapes = {
        "bare list": [[(0, 0), 1.0]],
        "none period": {0: [[(0, 0), 1.0]], 1: None},
        "empty dict": {},
    }
    for name in ("chd", "ghb", "drn", "riv", "wel", "rch", "evt"):
        helper = getattr(mf, name)
        for label, spd in shapes.items():
            spec = helper.flopy(stress_period_data=spd)
            assert spec.name == name, f"{name} / {label}"
            # maxbound is FloPy's job; pre-computing it is what broke these
            assert "maxbound" not in spec.options, f"{name} / {label}"
            # "verbatim" means verbatim: the payload reaches the spec untouched,
            # including a period explicitly set to None
            assert spec.options["stress_period_data"] == spd, f"{name} / {label}"
    assert (
        mf.rch.flopy(stress_period_data=shapes["none period"])
        .options["stress_period_data"][1]
        is None
    )


def test_riv_evt_forward_optional_package_arguments():
    # nseg and auxiliary are forwarded, not silently defaulted: deleting either
    # forward in _EVTPackage/_RIVPackage would leave the spec at its default
    # and this is the only assertion that would notice.
    segmented = mf.evt(
        stress_period_data={0: [[(0, 3), 100.0, 2.0e-3, 2.5, 0.5, 0.3]]}, nseg=2
    )
    assert segmented.options["nseg"] == 2
    assert mf.evt.flopy(
        stress_period_data={0: [[(0, 3), 100.0, 2.0e-3, 2.5, 0.5, 0.3]]}, nseg=2
    ).options["nseg"] == 2

    riv_aux = mf.riv(
        stress_period_data={0: [[(0, 7), 98.0, 40.0, 96.5, 1.0]]}, auxiliary=["temp"]
    )
    assert riv_aux.options["auxiliary"] == ["temp"]
    assert mf.riv.flopy(
        stress_period_data={0: [[(0, 7), 98.0, 40.0, 96.5, 1.0]]}, auxiliary=["temp"]
    ).options["auxiliary"] == ["temp"]
    assert mf.evt(
        stress_period_data={0: [[(0, 3), 100.0, 2.0e-3, 2.5, 1.0]]}, auxiliary=["temp"]
    ).options["auxiliary"] == ["temp"]


def test_package_api_rch_uses_builder_by_default_and_flopy_for_direct_data():
    context = mf.ModelContext(domain=np.array([[1, 1]]))

    built = mf.rch(context=context, nper=1, recharge=1.0e-4)
    direct = mf.rch.flopy(stress_period_data={0: [[(0, 0), 1.0e-4]]})

    assert built.metadata["builder"] == "RCHBuilder"
    assert built.options["stress_period_data"][0] == [[(0, 0), 1.0e-4], [(0, 1), 1.0e-4]]
    assert direct.metadata == {}


def test_simple_list_bcs_have_a_real_flopy_escape_hatch():
    # D8: chd/ghb/drn/wel each carry .flopy(...) (thin raw-FloPy passthrough,
    # same spec as the direct () form) so all list BCs share the three entry
    # points their docstrings advertise: () / .gpkg / .flopy.
    cases = {
        "chd": {0: [[(0, 0), 10.0]]},
        "ghb": {0: [[(0, 5), 86.0, 50.0]]},
        "drn": {0: [[(0, 12), 95.0, 30.0]]},
        "riv": {0: [[(0, 7), 98.0, 40.0, 96.5]]},
        "evt": {0: [[(0, 3), 100.0, 2.0e-3, 2.5]]},
        "wel": {0: [[(0, 42), -500.0]]},
    }
    for name, spd in cases.items():
        helper = getattr(mf, name)
        direct = helper(stress_period_data=spd)
        hatch = helper.flopy(stress_period_data=spd)
        assert hatch.name == name
        assert hatch.options == direct.options
        # builders are functools.partial instances (never compare equal
        # directly); same build function + same FloPy package class
        assert hatch.builder.func is direct.builder.func
        assert hatch.builder.args == direct.builder.args


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

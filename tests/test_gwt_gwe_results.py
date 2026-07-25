"""GWT/GWE results tier -- ``model.conc`` / ``model.temp`` (plan §6.0/6.1/6.2).

The 6.0 keystone factored the heads reader into a generic dependent-variable
surface and made the model view kind-aware; 6.1/6.2 instantiate it for GWT
concentration and GWE temperature. These end-to-end tests build a small coupled
GWF+GWT and GWF+GWE model on a Voronoi grid, run MF6, and assert the transport
field is readable through the SAME grammar heads use (get/summary/array/map),
with a field-appropriate hover and colorscale -- and that the kind guard keeps
``.hds`` off a transport model.
"""

from __future__ import annotations

import numpy as np
import pytest

import myflopy as mf
from myflopy.modflow.mf6.canonical import irregular_voronoi_grid
from myflopy.specs import ModelContext


def _coupled_run(tmp_path, kind: str):
    """Build + run a small coupled GWF+(GWT|GWE) model; return the built ``Run``."""

    vor = irregular_voronoi_grid(nrow=8, ncol=8, cell_size=100.0)
    gp = vor.get_disv_gridprops()
    ctx = ModelContext(grid=vor, domain=np.ones((1, vor.ncpl), dtype=int))

    cx = np.asarray(vor.centroids_x, dtype=float)
    span = cx.max() - cx.min()
    west = [int(c) for c in np.where(cx <= cx.min() + 0.15 * span)[0]]
    east = [int(c) for c in np.where(cx >= cx.max() - 0.15 * span)[0]]

    def disv():
        return mf.disv(nlay=1, ncpl=gp["ncpl"], nvert=len(gp["vertices"]),
                       vertices=gp["vertices"], cell2d=gp["cell2d"], top=10.0, botm=0.0)

    aux = "concentration" if kind == "gwt" else "temperature"
    chd_spd = [[(0, c), 10.0, 1.0] for c in west] + [[(0, c), 8.0, 0.0] for c in east]
    flow = mf.gwf("flow", context=ctx, packages=[
        disv(), mf.ic(strt=9.0), mf.npf(k=1.0, save_specific_discharge=True),
        mf.chd(stress_period_data={0: chd_spd}, auxiliary=aux),
        mf.oc(head_filerecord="flow.hds", budget_filerecord="flow.cbc",
              saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]),
    ])

    if kind == "gwt":
        transport = mf.gwt("trans", context=ctx, packages=[
            disv(), mf.ic(strt=0.0), mf.adv(scheme="UPSTREAM"), mf.mst(porosity=0.2),
            mf.ssm(sources=[["chd", "AUX", "concentration"]]),
            mf.oc(concentration_filerecord="trans.ucn", saverecord=[("CONCENTRATION", "ALL")]),
        ])
        exchange = mf.ExchangeSpec("gwfgwt", mf.build_gwf_gwt_exchange, models=("flow", "trans"))
    else:
        transport = mf.gwe("trans", context=ctx, packages=[
            disv(), mf.ic(strt=0.0), mf.adv(scheme="UPSTREAM"),
            mf.est(porosity=0.2, heat_capacity_solid=800.0, density_solid=2650.0,
                   heat_capacity_water=4184.0, density_water=1000.0),
            mf.cnd(ktw=0.6, kts=0.5, alh=1.0, ath1=0.1),
            mf.ssm(sources=[["chd", "AUX", "temperature"]]),
            mf.oc(temperature_filerecord="trans.ucn", saverecord=[("TEMPERATURE", "ALL")]),
        ])
        exchange = mf.ExchangeSpec("gwfgwe", mf.build_gwf_gwe_exchange, models=("flow", "trans"))

    sim = mf.SimulationSpec(f"{kind}cpl", models=[flow, transport], packages=[
        mf.tdis(nper=1, perioddata=[(1.0, 1, 1.0)]),
        mf.ims(name="ims_flow", models=["flow"], complexity="SIMPLE"),
        mf.ims(name="ims_trans", models=["trans"], complexity="SIMPLE",
               linear_acceleration="BICGSTAB"),
    ], exchanges=[exchange])

    project = mf.Project(tmp_path / f"{kind}proj", name="p")
    run = project.prepare_run("r", sim, overwrite=True)
    success, report = run.execute()
    assert success, f"{kind} MF6 did not converge:\n" + "\n".join(report[-20:])
    return run


@pytest.fixture(scope="module")
def gwt_run(tmp_path_factory):
    return _coupled_run(tmp_path_factory.mktemp("gwt"), "gwt")


@pytest.fixture(scope="module")
def gwe_run(tmp_path_factory):
    return _coupled_run(tmp_path_factory.mktemp("gwe"), "gwe")


@pytest.mark.slow
def test_gwt_conc_reader_grammar(gwt_run):
    """model.conc reads concentration through the shared grammar (get/summary/array)."""

    model = gwt_run.model("trans")
    assert model.model_type == "gwt6"

    frame = model.conc.get()
    assert list(frame.columns) == ["per", "layer", "cell", "conc"]
    assert frame["conc"].notna().any()

    summary = model.conc.summary()
    assert summary.loc[0, "label"] == "conc"
    assert summary.loc[0, "max"] > 0  # the plume carries mass off the west source

    arr = model.conc.array(layer=0)
    assert arr.shape == (model.vor.ncpl,)
    assert np.isfinite(arr).any()


@pytest.mark.slow
def test_gwt_conc_map_hook_and_hover(gwt_run):
    """model.conc.map() wires the value-kind hook: conc reader, conc hover, 'earth'."""

    choro = gwt_run.model("trans").conc.map()
    # value-kind hook: the map reads the conc reader and labels the layer table 'conc'.
    from myflopy.modflow.mf6.headsplus import ConcResults

    assert choro._value_name == "conc"
    assert choro._value_column == "conc"
    assert isinstance(choro._depvar_reader, ConcResults)
    assert choro.colorscale == "earth"  # plan §6.1.2: concentration default
    # the default hover is the concentration spec (not heads).
    assert choro.hover_spec.primary == "conc"
    assert choro.hover_spec.title == "Concentration"
    # the per-cell z-values the choropleth colors are the concentration field.
    zs = np.asarray(choro.zs, dtype=float)
    assert np.nanmax(zs) > 0


@pytest.mark.slow
def test_hds_guard_on_transport_model(gwt_run):
    """A GWT model exposes .conc, not .hds -- the guard raises a clear error."""

    model = gwt_run.model("trans")
    with pytest.raises(AttributeError, match="is a GWT model.*heads"):
        model.hds


@pytest.mark.slow
def test_gwe_temp_reader_and_map(gwe_run):
    """model.temp reads temperature through the same grammar, with a temp hover."""

    model = gwe_run.model("trans")
    assert model.model_type == "gwe6"

    frame = model.temp.get()
    assert list(frame.columns) == ["per", "layer", "cell", "temp"]
    assert model.temp.summary().loc[0, "label"] == "temp"

    choro = model.temp.map()
    assert choro._value_name == "temp"
    assert choro._value_column == "temp"
    assert choro.colorscale == "earth"
    assert choro.hover_spec.title == "Temperature"

    with pytest.raises(AttributeError, match="is a GWE model.*heads"):
        model.hds

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


def _coupled_run(tmp_path, kind: str, *, porosity: float = 0.2):
    """Build + run a small coupled GWF+(GWT|GWE) model; return the built ``Run``.

    ``porosity`` is a knob so callers can build two models that differ in one
    physical parameter -- what a grouped comparison needs to be meaningful.
    """

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
            disv(), mf.ic(strt=0.0), mf.adv(scheme="UPSTREAM"), mf.mst(porosity=porosity),
            mf.ssm(sources=[["chd", "AUX", "concentration"]]),
            mf.oc(concentration_filerecord="trans.ucn", budget_filerecord="trans.cbc",
                  saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")]),
        ])
        exchange = mf.ExchangeSpec("gwfgwt", mf.build_gwf_gwt_exchange, models=("flow", "trans"))
    else:
        transport = mf.gwe("trans", context=ctx, packages=[
            disv(), mf.ic(strt=0.0), mf.adv(scheme="UPSTREAM"),
            mf.est(porosity=porosity, heat_capacity_solid=800.0, density_solid=2650.0,
                   heat_capacity_water=4184.0, density_water=1000.0),
            mf.cnd(ktw=0.6, kts=0.5, alh=1.0, ath1=0.1),
            mf.ssm(sources=[["chd", "AUX", "temperature"]]),
            mf.oc(temperature_filerecord="trans.ucn", budget_filerecord="trans.cbc",
                  saverecord=[("TEMPERATURE", "ALL"), ("BUDGET", "ALL")]),
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


def test_a_multi_model_workspace_is_discoverable_without_loading_flopy(tmp_path):
    """Reopening ANY model of a coupled run used to raise RecursionError.

    A simulation with several models gives each its own subdirectory, so nothing
    but `mfsim.nam` sits at the top. Package discovery globbed the workspace
    ROOT, found nothing, left the grid type "unknown", and the `grid_type`
    property then fell back to inspecting `self.gwf` -- whose loader consults
    `grid_type`. Every reopened coupled model, transport or flow, recursed to
    death.

    Deliberately a file-tree fixture, not a real run: the failure was in cheap
    filename discovery, long before FloPy was asked for anything, and pinning it
    that way keeps it fast and independent of MF6.
    """

    from myflopy.project.run_model import load_mf6_run

    (tmp_path / "mfsim.nam").write_text(
        "BEGIN models\n"
        "  gwf6  flow/flow.nam  flow\n"
        "  gwt6  trans/trans.nam  trans\n"
        "END models\n",
        encoding="utf-8",
    )
    for name, suffix in (("flow", "disv"), ("trans", "disv")):
        (tmp_path / name).mkdir()
        (tmp_path / name / f"{name}.nam").touch()
        (tmp_path / name / f"{name}.{suffix}").touch()

    flow = load_mf6_run(tmp_path, model_name="flow")
    trans = load_mf6_run(tmp_path, model_name="trans")

    assert flow.grid_type == "disv" and trans.grid_type == "disv"
    assert "DISV" in flow.package_names

    # ...and the reopened transport model must know it is transport, or it
    # advertises `.hds` and hands back the wrong physics under a familiar name.
    assert trans.model_type == "gwt6"
    assert flow.model_type == "gwf6"

    # a single-model workspace stays flat and must keep working
    flat = tmp_path / "solo"
    flat.mkdir()
    (flat / "mfsim.nam").write_text(
        "BEGIN models\n  gwf6  solo.nam  solo\nEND models\n", encoding="utf-8"
    )
    (flat / "solo.nam").touch()
    (flat / "solo.disv").touch()
    assert load_mf6_run(flat).name == "solo"
    assert load_mf6_run(flat).grid_type == "disv"


def test_the_lazy_loader_does_not_recurse_when_no_grid_type_is_discovered(tmp_path):
    """The recursion guard, pinned independently of the discovery fix.

    Fixing discovery hides this: with a grid type discovered from filenames,
    `grid_type` answers from the override and never reaches `self.gwf`. But the
    cycle is still there for any workspace whose grid package is not
    recognizable -- `grid_type` -> `gwf` -> `_ensure_core_loaded` -> `grid_type`.
    So assert the loader uses the cheap file-derived answer, not the property.
    """

    from myflopy.project.run_model import load_mf6_run

    (tmp_path / "mfsim.nam").write_text(
        "BEGIN models\n  gwf6  solo.nam  solo\nEND models\n", encoding="utf-8"
    )
    (tmp_path / "solo.nam").touch()

    model = load_mf6_run(tmp_path, model_name="solo")
    assert model._grid_type_override is None, "no grid package -> nothing discovered"

    try:
        model._ensure_core_loaded()
    except RecursionError:
        pytest.fail("the loader re-entered itself through the grid_type property")
    except Exception:
        # FloPy cannot load this stub workspace, which is fine and beside the
        # point: what matters is that it got as far as trying.
        pass


@pytest.mark.slow
def test_a_transport_model_actually_writes_its_budget(gwt_run, gwe_run):
    """MF6 writes a ZERO-BYTE .cbc for a transport model unless SAVE_FLOWS is set
    on the model itself -- and FloPy then raises "file is empty" rather than
    saying what is missing. `mf.gwt`/`mf.gwe` default `save_flows=True` so that
    asking OC for a budget actually yields one.

    Asserted against real coupled runs of BOTH kinds, because the term names
    differ between them and are not what the plan assumed: there is no term
    called "SSM" (the SSM package's record is "SOURCE-SINK MIX"), and GWE's
    storage term is STORAGE-CELLBLK, not STORAGE-AQUEOUS.
    """

    for run, storage in ((gwt_run, "STORAGE-AQUEOUS"), (gwe_run, "STORAGE-CELLBLK")):
        model = run.model("trans")
        reader = model._get_budget_reader()
        assert reader is not None, "no budget file was written at all"

        terms = {name.strip() for name in reader.get_unique_record_names(decode=True)}
        assert terms == {storage, "FLOW-JA-FACE", "SOURCE-SINK MIX"}


@pytest.mark.slow
def test_a_transport_storage_term_builds_a_real_table_not_an_empty_one(gwt_run, gwe_run):
    """Two of a transport model's three budget terms silently returned NOTHING.

    MF6 writes storage (imeth=1) as a plain ``(nlay, 1, ncpl)`` float array with
    no ``node`` column. ``pd.DataFrame.from_records`` turns that into a nonsense
    ``(1, 1)`` frame with an integer column name, which the table builder's
    ``"node" not in frame.columns`` guard then skipped -- so
    ``STORAGE-AQUEOUS``/``STORAGE-CELLBLK`` came back as a ``(0, 7)`` table with
    no error at all, and a caller could not tell "no flow" from "not read".

    The strong assertion here is not the shape but the MASS BALANCE: over one
    steady period the storage term and the source-sink term must cancel. A frame
    of the right shape carrying the wrong values would fail that, so it pins the
    values and not merely the plumbing.
    """

    from myflopy.modflow.mf6.package_budget import build_budget_result_table

    for run, storage in ((gwt_run, "STORAGE-AQUEOUS"), (gwe_run, "STORAGE-CELLBLK")):
        model = run.model("trans")
        ncpl = int(model.vor.ncpl)

        storage_table = build_budget_result_table(
            model, budget_text=storage, package_name="sto", value_name="q"
        )
        assert len(storage_table) == ncpl, f"{storage}: one row per cell"
        assert sorted(storage_table["cell"].unique()) == list(range(ncpl))

        ssm_table = build_budget_result_table(
            model, budget_text="SOURCE-SINK MIX", package_name="ssm", value_name="q"
        )
        assert not ssm_table.empty

        storage_total = float(storage_table["q"].sum())
        ssm_total = float(ssm_table["q"].sum())
        scale = max(abs(storage_total), abs(ssm_total))
        assert scale > 0, f"{storage}: the run moved no mass, so this proves nothing"
        assert abs(storage_total + ssm_total) / scale < 1e-6, (
            f"{storage}: budget does not balance against SOURCE-SINK MIX "
            f"({storage_total} vs {ssm_total}) -- the table's VALUES are wrong"
        )


@pytest.mark.slow
def test_a_connection_indexed_transport_term_refuses_instead_of_going_quiet(gwt_run):
    """FLOW-JA-FACE is indexed by cell CONNECTION, so it is not a cell table.

    On this 64-cell grid it is a full array of 388 values. There is no honest
    mapping onto cells, and the bug being fixed was *silence* -- so the
    replacement has to be an error naming the mismatch, not an empty frame.
    """

    from myflopy.modflow.mf6.package_budget import build_budget_result_table

    with pytest.raises(ValueError, match="indexed by cell CONNECTION"):
        build_budget_result_table(
            gwt_run.model("trans"),
            budget_text="FLOW-JA-FACE",
            package_name="npf",
            value_name="q",
        )


@pytest.mark.slow
def test_transport_budget_nodes_agree_with_the_boundary_that_made_them(gwt_run, gwe_run):
    """Both budget paths must land on the same zero-based cells as the CHD source.

    ``package_budget`` zero-bases with its own rule while ``model.bud()`` used a
    registry lookup keyed on package name -- which no transport record matches,
    so the legacy path handed back MF6's raw 1-based ids while the modern path
    handed back zero-based ones. Ground truth is the CHD package that injects
    the mass: the SSM record's cells are exactly the cells CHD occupies.
    """

    from myflopy.modflow.mf6.package_budget import build_budget_result_table

    for run in (gwt_run, gwe_run):
        model = run.model("trans")
        chd = run.model("flow").gwf.get_package("chd").stress_period_data.get_data(key=0)
        truth = {int(row[0][1]) for row in chd}

        modern = build_budget_result_table(
            model, budget_text="SOURCE-SINK MIX", package_name="ssm", value_name="q"
        )
        modern_cells = {int(value) for value in modern["cell"].unique()}

        legacy = model.bud("SOURCE-SINK MIX").df.reset_index()
        legacy_cells = {int(value) for value in legacy["node"].dropna().unique()}

        assert modern_cells == truth, "the modern explorer path drifted off the CHD cells"
        assert legacy_cells == truth, "model.bud() returned raw 1-based node ids"


@pytest.mark.slow
def test_the_ssm_budget_is_reachable_by_its_package_name(gwt_run):
    """``model.bud("ssm")`` was unreachable on every transport model.

    MF6 names the record for the PROCESS ("SOURCE-SINK MIX"), not for the
    package that wrote it, and the lookup substring-matched the package name
    against record names -- so the obvious call raised "not included in budget
    file" while the non-obvious literal worked.
    """

    model = gwt_run.model("trans")

    by_alias = model.bud("ssm").df
    by_record = model.bud("SOURCE-SINK MIX").df
    assert by_alias.equals(by_record)

    # a genuinely absent package must still fail, and now says what IS available
    with pytest.raises(ValueError, match="Available records"):
        model.bud("drn")


@pytest.mark.slow
def test_the_xs_verb_works_on_the_transport_readers(gwt_run, gwe_run):
    """``model.conc.section()`` / ``.temp.section()`` raised from the day they shipped.

    ``XSection`` cached its value table from ``model.hds``, which the §6.0 kind
    guard refuses on a transport model -- so ``xs`` raised "is a GWT model;
    '.hds' is only available on GWF models" while the docs advertised the full
    grammar. It now reads the model's OWN dependent variable. The column name
    differs per kind (elev/conc/temp) but never mattered: the one consumer reads
    values positionally (ledger 99).
    """

    from myflopy.viz import Fig

    for run, reader_name in ((gwt_run, "conc"), (gwe_run, "temp")):
        model = run.model("trans")
        reader = getattr(model, reader_name)
        assert isinstance(reader.section(cells=[0, 5, 10]), Fig)

        # the section really carries THIS model's field, not heads
        section = reader._sections(cells=[0, 5, 10])[model.name]
        assert reader_name in section.all_heads.columns


@pytest.mark.slow
def test_model_budget_terms_follow_the_model_kind(gwt_run, gwe_run):
    """``model.budget.<term>`` is discovered, so it names each kind's real terms.

    This is why the namespace cannot be hand-declared the way the package-level
    ``<pkg>.budget.<term>`` namespaces are: GWT's storage term is
    ``STORAGE-AQUEOUS`` and GWE's is ``STORAGE-CELLBLK``, so a single declared
    list would be wrong for one of them. Closes plan §6.1/6.2 item 3.
    """

    from myflopy.viz import Fig

    for run, storage_attr in (
        (gwt_run, "storage_aqueous"),
        (gwe_run, "storage_cellblk"),
    ):
        model = run.model("trans")
        exposed = {name for name in dir(model.budget) if not name.startswith("_")}
        assert {storage_attr, "source_sink_mix", "flow_ja_face"} <= exposed
        # the OTHER kind's storage term must not be there
        assert "storage_cellblk" in exposed or "storage_aqueous" in exposed
        assert not {"storage_aqueous", "storage_cellblk"} <= exposed

        source = model.budget.source_sink_mix
        assert source.summary().loc[0, "label"] == "budget.source_sink_mix"
        assert isinstance(source.map(), object) and source.map().colorscale == "RdBu"
        cells = sorted(source.get()["cell"].unique())[:2]
        assert isinstance(source.plot(cells=cells), Fig)

        # the budget still balances when read through the new noun, which pins
        # the VALUES rather than only the plumbing
        storage_total = float(getattr(model.budget, storage_attr).get()["q"].sum())
        source_total = float(source.get()["q"].sum())
        scale = max(abs(storage_total), abs(source_total))
        assert scale > 0
        assert abs(storage_total + source_total) / scale < 1e-6


@pytest.mark.slow
def test_budget_hover_units_follow_the_model_kind(gwt_run, gwe_run, canonical_run):
    """A transport budget is not measured in cubic feet.

    Every hover hardcoded ``"ft³/d"``, which was wrong twice over: the wrong
    DIMENSION on a transport model (GWT carries mass per time, GWE energy per
    time) and the wrong UNITS on any GWF model not declared in feet and days.
    The canonical model does declare feet and days, so it must still read
    ``ft³/d`` -- the fix derives that label rather than hardcoding it.
    """

    from myflopy.modflow.mf6.package_budget import budget_value_units

    assert budget_value_units(gwt_run.model("trans")) == "M/T"
    assert budget_value_units(gwe_run.model("trans")) == "E/T"
    assert budget_value_units(canonical_run) == "ft³/d"
    # ...and a model that declares no units says so dimensionally instead of
    # inventing feet and days.
    assert budget_value_units(gwt_run.model("flow")) == "L³/T"

import os

os.environ.setdefault("MPLBACKEND", "Agg")

import pytest

import myflopy.modflow.mf6.package_budget as package_budget
import myflopy.modflow.mf6.package_explorer as package_explorer
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


def test_every_cell_bc_default_input_is_pinned():
    """``default_input`` drives every bare ``inputs.map()`` -- pin all of them.

    Table-driven on purpose (same shape as ``test_colorscale_policy``): rch/drn/wel
    were pinned individually, so chd/ghb/riv/evt drifted in unpinned. Swapping
    riv "stage" -> "cond" or evt "rate" -> "depth" passed the whole suite before
    this existed, silently changing what every bare riv/evt map draws.

    The rule the values follow: head-like BCs default to their driving head,
    flux-like BCs default to their flux.
    """

    expected = {
        "chd": "head",       # driving head
        "drn": "elev",       # driving head (drain elevation)
        "ghb": "bhead",      # driving head
        "riv": "stage",      # driving head
        "rch": "recharge",   # flux
        "wel": "q",          # flux
        "evt": "rate",       # flux
    }
    from myflopy.modflow.mf6.package_registry import _PACKAGE_EXPLORER_SPECS

    cell_bcs = {
        name for name, spec in _PACKAGE_EXPLORER_SPECS.items()
        if spec.kind == "cell_stress"
    }
    assert cell_bcs == set(expected), (
        "a cell BC was added or removed without pinning its default_input: "
        f"{sorted(cell_bcs ^ set(expected))}"
    )
    for package, field in expected.items():
        assert get_default_package_value_column(package) == field, package
        # the default must be a real declared field, not a typo
        assert field in _PACKAGE_EXPLORER_SPECS[package].inputs, package


def test_package_explorer_registry_preserves_existing_defaults():
    assert get_default_package_value_column("rch") == "recharge"
    assert get_default_package_value_column("drn") == "elev"
    # colorscale policy: non-signed data uses the house brown-to-blue scale
    assert get_default_package_colorscale("ghb") == "earth"
    assert get_default_budget_term("uzf_gwrch") == ("UZF-GWRCH", "gwrch")
    assert get_package_explorer_spec("lak") is registry_spec("lak")


def test_package_explorer_registry_adds_wel_semantics():
    spec = get_package_explorer_spec("wel")

    assert spec is not None
    assert spec.kind == "cell_stress"
    assert spec.default_input == "q"
    assert spec.inputs["q"].label == "Well flow"
    assert get_default_budget_term("wel") == ("WEL", "q_gwf")


def test_package_explorer_registry_exposes_uzf_period_fields():
    spec = get_package_explorer_spec("uzf")

    assert spec is not None
    assert {"finf", "pet", "extdp", "extwc", "ha", "hroot", "rootact"}.issubset(
        spec.inputs
    )
    assert spec.inputs["pet"].label == "Potential evapotranspiration"
    assert get_default_package_colorscale("uzf_pet") == "earth"


def test_model_packages_exposes_registry_backed_wel_accessors():
    packages = ModelPackages(_DummyModel())

    assert packages.wel.package_name == "wel"
    assert packages.wel.inputs.q.field_name == "q"
    assert packages.wel.results.q.package_name == "wel"
    assert packages.wel.results.q.budget_text == "WEL"
    assert packages.rch.results.fields["field"].tolist() == ["q"]
    # accessor stays results.q; the emitted COLUMN names its frame (q_gwf)
    assert packages.rch.results.q.result_name == "q"
    assert packages.rch.results.q.value_name == "q_gwf"
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


# --- model.packages.summary() / .mosaic() (2026-08-30) ---------------------- #

def test_summary_names_every_package_and_marks_what_can_be_drawn(canonical_model):
    """The map from "what is in this model" to "what can I look at".

    Without it you have to know that record packages answer
    `.inputs.<field>.map()` while array packages answer `.<field>.map()`, and
    that several packages have no explorer at all.
    """

    frame = canonical_model.packages.summary()

    assert list(frame.columns) == [
        "package", "kind", "mappable", "fields", "results",
    ]
    assert list(frame["package"]) == list(canonical_model.package_names), \
        "every package in the model gets a row, in model order"

    by_name = frame.set_index("package")
    assert by_name.loc["CHD", "mappable"]
    assert by_name.loc["CHD", "fields"] == "head"
    assert by_name.loc["NPF", "kind"] == "static_array"
    assert by_name.loc["NPF", "fields"] == "k, k22, k33"
    # OC has no explorer and nothing to draw; saying so is the point.
    assert not by_name.loc["OC", "mappable"]


def test_summary_is_cheap_by_default(canonical_model, monkeypatch):
    """`LoadedMf6Run` overrides package discovery precisely to avoid a full
    load, so the default must not read package records. `detail='data'` is where
    that cost is opted into."""

    calls = []
    namespace = type(canonical_model.packages)
    original = namespace._package_data_counts
    monkeypatch.setattr(
        namespace, "_package_data_counts",
        lambda self, row: calls.append(row["package"]) or original(self, row),
    )

    canonical_model.packages.summary()
    assert calls == [], "the cheap path read package data"

    canonical_model.packages.summary(detail="data")
    assert calls, "detail='data' did not read package data"


def test_summary_data_reports_the_layers_a_package_actually_occupies(canonical_model):
    """The column that explains an empty map before you draw one.

    `wel.inputs.map()` defaults to layer 0, and the canonical wells are in
    layers 1 and 3 -- which used to surface as a `KeyError` naming a column the
    package certainly has.
    """

    frame = canonical_model.packages.summary(detail="data").set_index("package")

    assert frame.loc["WEL", "layers"] == "1, 3"
    assert int(frame.loc["WEL", "records"]) == 12
    # UZF keeps its records per FIELD, with no whole-package get(); counting the
    # first declared field is still the right answer for a summary.
    assert int(frame.loc["UZF", "records"]) > 0
    # A package with no readable records stays blank rather than turning the
    # whole column into floats.
    assert str(frame["records"].dtype) == "Int64"


def test_summary_rejects_an_unknown_detail(canonical_model):
    with pytest.raises(ValueError, match="detail must be 'fields' or 'data'"):
        canonical_model.packages.summary(detail="everything")


def test_mosaic_draws_every_mappable_package(canonical_model):
    """The combinator over `summary()`: it draws what that table says is
    mappable, rather than making you name them."""

    table = canonical_model.packages.summary()
    expected = int(table["mappable"].sum())

    fig = canonical_model.packages.mosaic()
    assert len(fig.data) == expected

    subset = canonical_model.packages.mosaic(packages=["chd", "drn", "rch"])
    assert len(subset.data) == 3


def test_package_summary_is_gone(canonical_model):
    """Deleted rather than deprecated: zero callers, and not in the api
    snapshot, `__all__` or `__compatibility__`. `model.packages.summary()`
    replaces it, in the grammar the rest of the namespace already uses."""

    assert not hasattr(canonical_model, "package_summary")
    assert not hasattr(type(canonical_model), "package_summary")


# --- input maps must not need results (2026-08-30) -------------------------- #

def test_an_input_map_does_not_need_the_model_to_have_run(tmp_path):
    """Drawing INPUTS needs no results, and used to demand them anyway.

    `_grid_of` probes `hasattr(source, "hds")` to decide whether a source can
    serve a results field. `hasattr` swallows only `AttributeError`, so on an
    unrun model FloPy's `FileNotFoundError` escaped the CAPABILITY PROBE and
    killed every `map()` -- including the input maps that need no results file
    at all.
    """

    from myflopy.modflow.mf6.canonical_example import (
        CanonicalModelConfig,
        build_canonical_model,
    )

    model = build_canonical_model(tmp_path / "unrun", config=CanonicalModelConfig.testing())
    assert not (tmp_path / "unrun" / model.name / f"{model.name}.hds").exists(), \
        "fixture ran the model; the regression cannot reproduce"

    for package in ("chd", "drn", "rch"):
        picture = getattr(model.packages, package).inputs.map()
        assert picture.__class__.__name__ == "Choro"

    # And the whole namespace still answers, which is what makes summary/mosaic
    # usable while building a model rather than only after running one.
    assert len(model.packages.summary()) == len(model.package_names)
    assert model.packages.mosaic() is not None


def test_a_period_or_layer_with_no_records_draws_an_empty_map(canonical_model):
    """`wel` sits in layers 1 and 3, so the default `layer=0` selects nothing.

    An empty selection comes back WITHOUT its value column, so guarding the
    column before the emptiness turned "no records here" into
    `KeyError: Value column 'q' was not found` -- naming a column the package
    certainly has, and making the all-fill branch below it unreachable.
    """

    empty = canonical_model.packages.wel.inputs.get(per=0, layer=0)
    assert empty.empty and "q" not in empty.columns, "the trap stopped reproducing"

    picture = canonical_model.packages.wel.inputs.map()          # layer 0
    assert picture.__class__.__name__ == "Choro"
    assert canonical_model.packages.wel.inputs.map(layer=3) is not None


def test_a_genuinely_missing_value_column_still_raises(canonical_model):
    """Reordering the guards must not lose the guard: a NON-empty frame missing
    the column is still a caller error, and now says what the frame does have."""

    from myflopy.modflow.mf6.package_plotting import build_cell_input_map_payload

    frame = canonical_model.packages.drn.inputs.get().drop(columns=["elev"])
    with pytest.raises(KeyError, match="the selection has"):
        build_cell_input_map_payload(
            frame, ncpl=canonical_model.vor.ncpl, value_column="elev",
            per=0, layer=0,
        )

def test_mappable_is_probed_not_inferred_from_the_registry():
    """HFB is deliberately NOT in `package_registry` -- it is face-indexed, so it
    has no cellid, and MF6 writes it no budget record (ledger 162). Its explorer
    is hand-written and it DRAWS: `model.packages.hfb.inputs.map()`.

    Inferring `mappable` from registry fields reported exactly that package --
    the one whose entry is entirely hand-written -- as undrawable. So the column
    probes for a real `map()` instead.
    """

    is_mappable = package_model.ModelPackages._is_mappable

    class HfbShaped:
        """The shape HFB actually has: `.inputs` is self, and self maps."""

        @property
        def inputs(self):
            return self

        def map(self):
            return object()

    # No registry fields, no results, no static arrays -- and still mappable.
    assert is_mappable(HfbShaped(), "", [])

    class ResultsOnly:
        """LAK/SFR: nothing to pick on the inputs side, but results draw."""

        class _R:
            def map(self):
                return object()

        results = _R()

    assert is_mappable(ResultsOnly(), "surface_water", [])

    class Inert:
        """An explorer that cannot draw anything."""

    assert not is_mappable(Inert(), "", [])
    assert not is_mappable(None, "cell_stress", ["head"])
    # Static arrays draw through `<field>.map()`, one level shallower.
    assert is_mappable(Inert(), "static_array", ["k"])
    assert not is_mappable(Inert(), "static_array", [])


def test_summary_uses_the_probe_for_a_registry_absent_package(canonical_model, monkeypatch):
    """The predicate above, reached through real `summary()` output.

    `ModelPackages.hfb` binds an explorer to any model, so naming HFB in the
    package list is enough to get its row -- no barrier needs to exist. That
    keeps this pinned to the summary path rather than to `_is_mappable` alone.
    """

    monkeypatch.setattr(
        type(canonical_model), "package_names",
        property(lambda self: ["NPF", "HFB"]),
    )
    row = canonical_model.packages.summary().set_index("package").loc["HFB"]

    assert row["mappable"], "the one hand-written package was reported undrawable"
    assert row["fields"] == "", "HFB maps as a whole; there is no field to pick"
    assert row["results"] == "", "MF6 writes HFB no budget record at all"


# --- the *_summary() family, retired onto nouns (2026-08-30) ---------------- #

@pytest.mark.parametrize("retired", ["package_summary", "output_summary", "result_summary"])
def test_the_zero_caller_summary_methods_are_gone(canonical_model, retired):
    """Deleted outright rather than warned: zero callers anywhere, absent from
    `tests/api_snapshot.json`, `__all__` and `__compatibility__`.

    `output_summary` and `result_summary` also hardcoded `f"{self.name}.hds"`,
    so on a GWT or GWE model they reported `has_heads=False` and null statistics
    silently -- where `model._field_reader` already dispatches
    HeadsPlus/ConcResults/TempResults correctly.
    """

    assert not hasattr(canonical_model, retired)
    assert not hasattr(type(canonical_model), retired)


def test_the_noun_that_replaces_result_summary_is_strictly_better(canonical_run):
    """`result_summary` is a pure delete, not a migration: `model.hds.summary()`
    already existed, is already tested, is kind-dispatched, and reports MORE.

    It also disagreed on the numbers -- `result_summary` read `get_data()` with
    no `kstpkper`, i.e. the LAST time step only -- for reasons no column
    explained.
    """

    columns = set(canonical_run.hds.summary().columns)
    assert {"records", "periods", "layers", "cells", "min", "max", "mean"} <= columns


# --- an input map's hover shows the whole record (2026-08-30) --------------- #

def _hover_text(picture) -> str:
    """The hovertemplate with its markup stripped, for readable assertions."""

    import re

    template = picture.fig.data[0].hovertemplate
    return re.sub(r"<span[^>]*>|</span>|</?b>|<extra></extra>", "", template)


@pytest.mark.parametrize(
    ("package", "coloured", "siblings"),
    [
        ("ghb", "bhead", ["cond"]),
        ("drn", "elev", ["cond"]),
        ("riv", "stage", ["cond", "rbot"]),
        ("evt", "rate", ["surface", "depth"]),
    ],
)
def test_an_input_hover_carries_the_packages_other_fields(
    canonical_model, package, coloured, siblings
):
    """Reading a GHB map means asking "what head, against what conductance".

    The payload already aggregated every field; only the hover spec was leaving
    them on the floor, so the hover showed the coloured field and nothing else.
    """

    text = _hover_text(getattr(canonical_model.packages, package).inputs.map())

    for sibling in siblings:
        assert sibling in text, f"{package} hover dropped {sibling}"
    # Twice and no more: once as the hover title, once as the primary value.
    # A third would mean the coloured field was also listed among its siblings.
    assert text.count(coloured) == 2, f"{package}: {coloured!r} is duplicated"


def test_the_hover_says_which_layer_you_are_looking_at(canonical_model):
    """The commonest reason a boundary map comes back empty is that the package
    has no records in the default layer 0 -- so the map should say which layer it
    is drawing."""

    assert "layer" in _hover_text(canonical_model.packages.ghb.inputs.map())


def test_field_selects_what_is_coloured_and_the_rest_ride_along(canonical_model):
    """`field=` drives the colour; the siblings stay in the hover either way."""

    default = _hover_text(canonical_model.packages.riv.inputs.map())
    by_cond = _hover_text(canonical_model.packages.riv.inputs.map(field="cond"))

    # The title is the coloured field -- `stage` by default (the registry's
    # default_input), `cond` when asked for.
    assert default.strip().startswith("stage")
    assert by_cond.strip().startswith("cond")
    for field in ("stage", "cond", "rbot"):
        assert field in default and field in by_cond


def test_record_count_appears_only_when_records_actually_merged():
    """`records 1` on every cell of every map is noise; the count earns its line
    only when aggregation happened, which is the one thing that explains a
    summed value."""

    from myflopy.modflow.mf6.package_inputs import _input_context_fields

    assert _input_context_fields({"Record Count": [1, 1, 1]}) == ("Layer",)
    assert _input_context_fields({"Record Count": [1, 3, 1]}) == ("Layer", "Record Count")
    assert _input_context_fields({}) == ("Layer",)


def test_a_single_field_package_gets_no_empty_sibling_block(canonical_model):
    """RCH has one field, and UZF/NPF payloads carry one array -- `Fields` skips
    names the payload lacks, so nothing renders an empty block."""

    for picture in (
        canonical_model.packages.rch.inputs.map(),
        canonical_model.packages.uzf.inputs.finf.map(),
        canonical_model.packages.npf.k.map(),
    ):
        text = _hover_text(picture)
        assert "layer" in text
        assert "%{customdata[2]}" in text     # layer is the 3rd column, not the 4th

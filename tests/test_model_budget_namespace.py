"""``model.budget.<term>`` -- every term in a model's own budget file, as a noun.

Plan §6.1/6.2 item 3. The terms are DISCOVERED from the file rather than declared,
because the set genuinely varies with model kind and configuration (GWT's storage
term is ``STORAGE-AQUEOUS`` and GWE's is ``STORAGE-CELLBLK``; a flow model carries
whichever boundary packages it has). These tests pin the namespace mechanics on a
GWF model; the kind-specific term names are pinned in ``test_gwt_gwe_results.py``.
"""

from __future__ import annotations

import pandas as pd
import pytest

from myflopy.modflow.mf6.package_results import (
    CellBudgetResultsExplorer,
    budget_term_attribute,
)


def test_the_record_name_normalizer_covers_mf6s_real_vocabulary():
    """MF6 separates words with hyphens, and in one case a space."""

    assert budget_term_attribute("SOURCE-SINK MIX") == "source_sink_mix"
    assert budget_term_attribute("STORAGE-AQUEOUS") == "storage_aqueous"
    assert budget_term_attribute("STORAGE-CELLBLK") == "storage_cellblk"
    assert budget_term_attribute("FLOW-JA-FACE") == "flow_ja_face"
    assert budget_term_attribute("UZF-GWRCH") == "uzf_gwrch"
    assert budget_term_attribute("STO-SS") == "sto_ss"
    assert budget_term_attribute("DATA-SPDIS") == "data_spdis"
    assert budget_term_attribute("CHD") == "chd"
    # already-normalized input is a fixed point, so __getitem__ accepts either
    assert budget_term_attribute("source_sink_mix") == "source_sink_mix"


def test_the_namespace_discovers_the_models_real_terms(canonical_run):
    """``dir()`` and ``types`` describe the same file, unlike the package namespaces.

    The package-level ``<pkg>.budget.<term>`` namespaces hand-write their terms,
    so ``dir()`` (declared) and ``types`` (present in the file) can silently
    disagree -- a term in the file but not declared is unreachable by attribute.
    Deriving one from the other makes that impossible here.
    """

    namespace = canonical_run.budget
    reader_terms = {
        str(name).strip()
        for name in canonical_run._get_budget_reader().get_unique_record_names(
            decode=True
        )
    }

    assert set(namespace.types) == reader_terms
    exposed = {name for name in dir(namespace) if not name.startswith("_")}
    assert {budget_term_attribute(term) for term in reader_terms} <= exposed


def test_a_term_answers_the_full_spatial_verb_set(canonical_run):
    """Each term is a spatial noun: get/summary/plot/map, per view_layer_conventions.

    This is why the namespace hands back ``CellBudgetResultsExplorer`` (a
    ``SpatialView``) rather than the ``PackageBudgetTermExplorer`` the package
    namespaces use -- that one has no ``plot()`` at all.
    """

    from myflopy.viz import Fig

    term = canonical_run.budget.drn
    assert isinstance(term, CellBudgetResultsExplorer)

    frame = term.get()
    assert isinstance(frame, pd.DataFrame) and not frame.empty
    assert {"per", "layer", "cell", "q"} <= set(frame.columns)

    # summary names the noun by where the user reached it, NOT as a package result
    assert term.summary().loc[0, "label"] == "budget.drn"

    assert isinstance(term.plot(cells=sorted(frame["cell"].unique())[:2]), Fig)
    assert term.map().colorscale == "RdBu"  # a signed flow gets the diverging scale


def test_a_term_agrees_with_the_package_results_path(canonical_run):
    """The new noun must not become a second, disagreeing source of truth.

    ``model.budget.drn`` and ``model.packages.drn.results.q`` read the same
    record. They differ ONLY in column name -- the package path renames ``q`` to
    ``q_gwf`` so the column carries its reference frame, while this namespace is
    explicitly the raw model-budget view and keeps MF6's own ``q``.
    """

    from_budget = canonical_run.budget.drn.get()
    from_package = canonical_run.packages.drn.results.q.get()

    assert "q" in from_budget.columns and "q_gwf" in from_package.columns
    assert set(from_budget["cell"]) == set(from_package["cell"])
    assert float(from_budget["q"].sum()) == pytest.approx(
        float(from_package["q_gwf"].sum())
    )


def test_terms_with_no_package_accessor_are_reachable(canonical_run):
    """The reason this is not gated to transport models.

    ``STO-SS``/``STO-SY``/``DATA-SPDIS``/``DATA-SAT`` have no package accessor at
    all, so before this noun the only route to them was the legacy
    ``model.bud("STO-SS").df`` wrapper.
    """

    storage = canonical_run.budget.sto_ss.get(per=0)
    assert not storage.empty
    # a full-array record covers every node of that period exactly once
    assert len(storage) == len(canonical_run.node_to_lni)
    assert sorted(storage["cell"].unique()) == sorted(
        {cell for _, cell in canonical_run.node_to_lni.values()}
    )
    # ...and get() spans every period, one full grid each
    assert len(canonical_run.budget.sto_ss.get()) == len(
        canonical_run.node_to_lni
    ) * len(canonical_run.kstpkper)

    assert not canonical_run.budget.data_spdis.get().empty

    # the hover names the term the way MF6 spells it, not the attribute spelling
    assert canonical_run.budget.sto_ss.map().hover_spec.title == "STO-SS q"


def test_either_spelling_reaches_the_same_term(canonical_run):
    """A name copied straight out of ``types`` must work, hyphens and all.

    The package namespaces' term filter upcases without converting ``_`` back to
    ``-``, so ``budget.get(term="ext_inflow")`` silently returns an empty frame
    there. Indexing normalizes both spellings to one answer instead.
    """

    by_record = canonical_run.budget["STO-SS"]
    by_attribute = canonical_run.budget["sto_ss"]

    assert by_record.budget_text == by_attribute.budget_text == "STO-SS"
    assert by_record.budget_text == canonical_run.budget.sto_ss.budget_text


def test_an_unknown_term_says_what_the_model_actually_has(canonical_run):
    """A typo must name the alternatives, not just fail."""

    with pytest.raises(AttributeError, match="has no budget term 'nope'"):
        canonical_run.budget.nope

    with pytest.raises(AttributeError, match="Available terms:.*drn"):
        canonical_run.budget.nope


def test_flow_ja_face_is_visible_but_refuses_to_pretend(canonical_run):
    """FLOW-JA-FACE stays in the namespace and explains itself (ledger 91/96).

    Hiding it would put ``dir()`` and ``types`` back in disagreement -- the term
    IS in the file. It is connection-indexed, so it cannot be a cell table, and
    the error says exactly that.
    """

    assert "flow_ja_face" in dir(canonical_run.budget)
    assert "FLOW-JA-FACE" in canonical_run.budget.types

    with pytest.raises(ValueError, match="indexed by cell CONNECTION"):
        canonical_run.budget.flow_ja_face.get()


def test_the_namespace_works_on_a_reopened_run(canonical_run):
    """A noun that works live but not on reopen is a known past bug class.

    ``ModelView`` and ``LoadedMf6Run`` both subclass ``SimulationBase``, so the
    property lands on both -- but the plumbing underneath differs
    (``LoadedMf6Run`` overrides ``_get_budget_reader`` to read the ``.cbc``
    directly). Assert the reopened run gives the same answer, not merely that the
    attribute exists.
    """

    from myflopy.project.run_model import load_mf6_run

    reopened = load_mf6_run(canonical_run.workspace, verbosity_level=0)

    assert set(reopened.budget.types) == set(canonical_run.budget.types)
    assert float(reopened.budget.drn.get()["q"].sum()) == pytest.approx(
        float(canonical_run.budget.drn.get()["q"].sum())
    )

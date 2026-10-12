"""The LAK/SFR package-budget views name their terms, and refuse what cannot act (ledger 473).

* ``model.sfr.budget.<term>.wide()`` defaulted to the LAK index ``("per", "lake")`` and
  raised ``KeyError: 'lake'`` on every SFR term;
* ``budget.get(term="ext_inflow")`` matched nothing and returned an empty table in
  silence; a term MODFLOW 6 never writes for the package now raises, naming the terms;
* ``model.sfr.budget.rain`` was pinned to ``"RAIN"``, which no SFR budget file holds
  (MF6 writes ``RAINFALL``), so it read nothing on every model;
* ``model.lak.connections.map(per=...)`` raised ``TypeError`` naming ``ModelPlots.map``
  while its spliced docstring documented ``per``.
"""

from __future__ import annotations

import inspect

import pytest

pytestmark = pytest.mark.slow


def test_an_sfr_terms_wide_indexes_by_reach(canonical_run):
    wide = canonical_run.sfr.budget.gwf.wide()
    assert {"per", "reach"} <= set(wide.columns)
    assert not wide.empty
    assert {"per", "lake"} <= set(canonical_run.lak.budget.gwf.wide().columns)


@pytest.mark.parametrize("package", ["lak", "sfr"])
def test_a_misspelled_term_is_refused_naming_the_mf6_spelling(canonical_run, package):
    budget = getattr(canonical_run, package).budget
    with pytest.raises(ValueError, match="'ext_inflow' is spelled 'EXT-INFLOW'"):
        budget.get(term="ext_inflow")
    with pytest.raises(ValueError, match="the file holds"):
        budget.summary(term="NOT-A-TERM")
    assert not budget.get(term="ext-inflow").empty


def test_a_term_the_file_lacks_is_empty_not_refused(canonical_run):
    budget = canonical_run.lak.budget
    absent = sorted({"WITHDRAWAL", "CONSTANT", "AUXILIARY"} - set(budget.types))
    assert absent, "the fixture's LAK file must lack one of MF6's optional terms"
    assert budget.get(term=absent[0]).empty
    assert canonical_run.lak.budget.auxiliary.get().empty


def test_the_sfr_rainfall_helper_reads_mf6s_record(canonical_run):
    budget = canonical_run.sfr.budget
    assert "RAINFALL" in budget.types
    rows = budget.rainfall.get()
    assert not rows.empty
    assert set(rows["term"]) == {"RAINFALL"}
    assert "RAINFALL" in set(budget.stream_fluxes.get()["term"])
    assert not hasattr(budget, "rain")


def test_lake_connections_take_no_period(canonical_run):
    connections = canonical_run.lak.connections
    with pytest.raises(TypeError, match=r"results\.q\.map\(per="):
        connections.map(per=1)
    doc = inspect.getdoc(type(connections).map)
    assert "\nper :" not in doc

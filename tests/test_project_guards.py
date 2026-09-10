"""``Project`` accepted specs into the wrong library, and went quiet when it lost work.

Two failures with the same shape -- every signal said success and the loss showed
up in a later session:

``project.add_package("flow_hi", mf.gwf(...))`` was accepted. ``validate()``
returned clean, ``save()`` reported success while silently writing no
simulations, and the NEXT session's ``Project.load`` raised ``KeyError:
'builder'`` -- the whole project unopenable, from one wrong verb.

A simulation carrying an inline array cannot be serialized, so ``save()`` skips
it and ``Project.load`` never sees it. That was announced as
``SaveReport(simulations=0, skipped=1)`` -- which is truthy, raises nothing, and
is read by nobody -- and at DEBUG on the run path. It is what left a real
project's ``specs/`` directory empty.
"""

from __future__ import annotations

import logging

import numpy as np
import pytest

import myflopy as mf

_DIS = dict(nlay=1, nrow=1, ncol=1, delr=1.0, delc=1.0, top=10.0, botm=0.0)


@pytest.fixture
def project(tmp_path):
    return mf.Project(tmp_path / "proj", name="p")


def _model_with_inline_array():
    return mf.gwf(
        "flow",
        packages=[
            mf.dis(**_DIS),
            mf.ic(strt=5.0),
            mf.npf(k=np.full((1, 1, 1), 2.5)),
            mf.oc(head_filerecord="flow.hds", saverecord=[("HEAD", "ALL")]),
        ],
    )


# --- the wrong-library guard ------------------------------------------------ #


def test_a_model_handed_to_add_package_is_refused(project):
    with pytest.raises(TypeError, match="was given a model"):
        project.add_package("flow_hi", mf.gwf("flow"))


def test_a_model_handed_to_add_grid_is_refused(project):
    with pytest.raises(TypeError, match="was given a model"):
        project.add_grid("mesh", mf.gwf("flow"))


def test_a_package_handed_to_add_simulation_is_refused(project):
    with pytest.raises(TypeError, match="was given a package"):
        project.add_simulation(mf.npf(k=1.0))


def test_the_refusal_names_the_verb_that_was_wanted(project):
    with pytest.raises(TypeError) as excinfo:
        project.add_package("flow_hi", mf.gwf("flow"))
    assert "project.add_simulation()" in str(excinfo.value)


def test_the_project_is_left_unchanged_by_a_refusal(project):
    with pytest.raises(TypeError):
        project.add_package("flow_hi", mf.gwf("flow"))
    assert project.packages.names() == []


def test_the_right_things_are_still_accepted(project):
    project.add_package("npf/base", mf.npf(k=1.0))
    project.add_grid("mesh", mf.GridSpec.python("g.py", function="build"))
    project.add_simulation(mf.simulation(mf.gwf("flow"), name="baseline"))

    assert project.packages.names() == ["npf/base"]
    assert project.grids.names() == ["mesh"]
    assert project.simulations.names() == ["baseline"]


def test_a_built_grid_object_is_not_mistaken_for_a_spec(project):
    """``add_grid`` legitimately takes a mesh, a GridSpec, or a GridRef."""

    from types import SimpleNamespace

    project.add_grid("mesh", SimpleNamespace(ncpl=3))
    assert project.grids.names() == ["mesh"]


# --- saying so at the moment of loss ---------------------------------------- #


def test_registering_an_unserializable_simulation_warns(project, caplog):
    spec = mf.simulation(_model_with_inline_array(), name="baseline")
    with caplog.at_level(logging.WARNING, logger="myflopy.workspace"):
        project.add_simulation(spec)

    assert "will NOT persist" in caplog.text
    assert "baseline" in caplog.text


def test_the_warning_names_the_fix(project, caplog):
    with caplog.at_level(logging.WARNING, logger="myflopy.workspace"):
        project.add_simulation(
            mf.simulation(_model_with_inline_array(), name="baseline")
        )
    assert "mf.ref(key)" in caplog.text


def test_a_serializable_simulation_warns_about_nothing(project, caplog):
    spec = mf.simulation(
        mf.gwf(
            "flow",
            packages=[mf.dis(**_DIS), mf.ic(strt=5.0), mf.ref("npf/base")],
        ),
        name="baseline",
    )
    with caplog.at_level(logging.WARNING, logger="myflopy.workspace"):
        project.add_simulation(spec)
    assert caplog.text == ""


def test_a_run_whose_spec_cannot_be_stored_warns_once(project, caplog):
    """DEBUG before, so the run silently stopped being re-runnable.

    De-duplicated because ``_record_status`` re-summarizes on every status
    change, which would otherwise repeat it three or four times per run.
    """

    project.add_simulation(mf.simulation(_model_with_inline_array(), name="baseline"))
    caplog.clear()

    with caplog.at_level(logging.WARNING, logger="myflopy.workspace"):
        run = project.prepare_run("baseline", "baseline")
        run.save_manifest()
        run.save_manifest()

    assert caplog.text.count("is not re-runnable") == 1
    assert "reopen_run('baseline').spec will be None" in caplog.text

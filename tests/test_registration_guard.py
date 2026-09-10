"""Registering a key that exists used to destroy the old entry, silently.

`packages.fork` refused the collision; `add_package` replaced without a word. So
whether re-using a key raised or destroyed something depended on which verb you
reached for -- and the guarded one was the deliberate operation while the
unguarded one was the accident: a second script run, a copied line, a notebook
cell executed twice.

`packages.edit(key, ...)` stays deliberately unguarded. Changing what a key MEANS
is the supported operation, and `project.runs.drift()` is what makes its
consequences visible.
"""

from __future__ import annotations

import pytest

import myflopy as mf


@pytest.fixture
def project(tmp_path):
    return mf.Project(tmp_path / "proj", name="p")


def _sim(name="baseline"):
    return mf.simulation(mf.gwf("flow", packages=[]), name=name)


# --- the guard -------------------------------------------------------------- #


def test_re_registering_a_package_key_is_refused(project):
    project.add_package("npf/base", mf.npf(k=1.0))
    with pytest.raises(KeyError, match="already has a package 'npf/base'"):
        project.add_package("npf/base", mf.npf(k=99.0))


def test_re_registering_a_grid_key_is_refused(project):
    project.add_grid("mesh", mf.GridSpec.python("g.py", function="build"))
    with pytest.raises(KeyError, match="already has a grid 'mesh'"):
        project.add_grid("mesh", mf.GridSpec.python("other.py", function="build"))


def test_re_registering_a_simulation_name_is_refused(project):
    project.add_simulation(_sim())
    with pytest.raises(KeyError, match="already has a simulation 'baseline'"):
        project.add_simulation(_sim())


def test_the_refused_entry_is_left_untouched(project):
    project.add_package("npf/base", mf.npf(k=1.0))
    with pytest.raises(KeyError):
        project.add_package("npf/base", mf.npf(k=99.0))
    assert project.package("npf/base").options["k"] == 1.0


# --- and the ways to mean it ------------------------------------------------ #


def test_overwrite_replaces(project):
    project.add_package("npf/base", mf.npf(k=1.0))
    project.add_package("npf/base", mf.npf(k=99.0), overwrite=True)
    assert project.package("npf/base").options["k"] == 99.0


def test_overwrite_works_for_grids_and_simulations(project):
    project.add_grid("mesh", mf.GridSpec.python("g.py", function="build"))
    project.add_grid("mesh", mf.GridSpec.python("g2.py", function="build"), overwrite=True)
    project.add_simulation(_sim())
    project.add_simulation(_sim(), overwrite=True)
    assert project.grids["mesh"].script.endswith("g2.py")


def test_edit_is_never_refused(project):
    """Changing what a key means is the supported operation, not the accident."""

    project.add_package("npf/base", mf.npf(k=1.0))
    project.packages.edit("npf/base", k=50.0)
    assert project.package("npf/base").options["k"] == 50.0


def test_the_message_points_at_edit_for_a_deliberate_change(project):
    project.add_package("npf/base", mf.npf(k=1.0))
    with pytest.raises(KeyError) as excinfo:
        project.add_package("npf/base", mf.npf(k=99.0))
    assert "packages.edit(...)" in str(excinfo.value)


def test_a_fresh_key_is_unaffected(project):
    project.add_package("npf/base", mf.npf(k=1.0))
    project.add_package("npf/high", mf.npf(k=99.0))
    assert project.packages.names() == ["npf/base", "npf/high"]


def test_loading_a_saved_project_still_works(project, tmp_path):
    """`Project.load` registers into a FRESH project, so nothing collides."""

    project.add_package("npf/base", mf.npf(k=1.0))
    project.add_simulation(_sim())
    project.save()

    reloaded = mf.Project.load(project.root)
    assert reloaded.packages.names() == ["npf/base"]
    assert reloaded.simulations.names() == ["baseline"]

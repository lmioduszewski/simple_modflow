"""GRASS progress output, and the cache invariant that makes the flag safe.

`mf.Contours(..., progress=True)` exists because GRASS modules are subprocesses
that write their percentages to the inherited file descriptors. In a terminal
that already reaches the screen; in Jupyter it reaches the KERNEL's console, not
the cell, and `redirect_stdout` cannot help because the writing happens below
Python. The flag swaps the fds for a pipe and re-emits through `sys.stdout`.

The tests that matter here are the two invariants, not the plumbing: the flag
must not change what gets computed, and must not invalidate the cached raster --
`grass_kwargs` feeds the derived-raster signature, so routing a display concern
through it would silently re-interpolate the first time anyone asked to watch.
"""

from __future__ import annotations

import io
import sys

import pytest

from myflopy.modflow.utils.contour_interp import _teed_console
from myflopy.surfaces import Surface


def test_progress_is_off_by_default():
    """Library code logs; it does not print. The flag is the human asking."""

    assert Surface.from_contours("c.gpkg").progress is False


def test_progress_is_not_part_of_the_cache_signature():
    """Two surfaces differing ONLY in `progress` must share a cache entry.

    Otherwise turning the flag on silently triggers a full re-interpolation --
    the opposite of what someone asking to watch the progress wants.
    """

    quiet = Surface.from_contours("c.gpkg", z="Elev", resolution=25)
    loud = Surface.from_contours("c.gpkg", z="Elev", resolution=25, progress=True)

    assert loud.progress and not quiet.progress
    assert quiet.grass_kwargs == loud.grass_kwargs == {}
    assert quiet._derived_raster().params == loud._derived_raster().params


def test_progress_does_not_leak_into_grass_kwargs():
    """It is a named parameter precisely so it cannot ride the passthrough."""

    surface = Surface.from_contours("c.gpkg", progress=True, region_vector="d.gpkg")
    assert "progress" not in surface.grass_kwargs
    assert surface.grass_kwargs == {"region_vector": "d.gpkg"}


def test_the_tee_is_a_no_op_when_disabled():
    """Off must not capture anything.

    Asserting `sys.stdout is unchanged` looked like the check and was not: the
    tee swaps FILE DESCRIPTORS and never touches `sys.stdout`, so that assertion
    held even with the flag ignored. Writing to fd 1 and finding nothing captured
    is the property that actually distinguishes on from off.
    """

    import os

    captured, real = io.StringIO(), sys.stdout
    sys.stdout = captured
    try:
        with _teed_console(False):
            os.write(1, b"this belongs on the terminal\n")
    finally:
        sys.stdout = real

    assert captured.getvalue() == ""


def test_the_tee_routes_low_level_writes_into_python_stdout():
    """The whole point: output written to fd 1 by a SUBPROCESS-style write must
    come back through `sys.stdout`, which is what a Jupyter cell displays."""

    import os

    captured, real = io.StringIO(), sys.stdout
    sys.stdout = captured                       # no fileno -> the Jupyter shape
    try:
        with _teed_console(True):
            os.write(1, b"37%  74%  100%\n")    # bypasses sys.stdout entirely
    finally:
        sys.stdout = real

    assert "100%" in captured.getvalue()


def test_the_tee_restores_the_file_descriptors():
    """A leaked dup would silently swallow every later write in the process."""

    import os

    saved = os.dup(1)
    try:
        with _teed_console(True):
            os.write(1, b"x")
        # fd 1 must once more be the same file as the duplicate taken before.
        assert os.fstat(1).st_ino == os.fstat(saved).st_ino
    finally:
        os.close(saved)


@pytest.mark.parametrize("enabled", [True, False])
def test_the_tee_restores_even_when_the_body_raises(enabled):
    import os

    saved = os.dup(1)
    try:
        with pytest.raises(RuntimeError), _teed_console(enabled):
            raise RuntimeError("boom")
        assert os.fstat(1).st_ino == os.fstat(saved).st_ino
    finally:
        os.close(saved)

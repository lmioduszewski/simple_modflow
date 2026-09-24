"""A cell with no head is blank in every view of the head file (ledger 315).

MF6 writes a DRY cell's head as -1e30 and an INACTIVE cell's as +1e30. The readers
share one rule for them, ``headsplus._is_sentinel`` (``|value| >= 1e29``, either
sign). The map and the section each tested a sign of their own, or none. On a
non-Newton run a dry cell mapped and hovered as -1e30 and dragged the colour range
of the map, the static map, ``mosaic`` and ``animate`` down to -1e30. It also made
the contour levels impossible, so no contour drew at all. An inactive cell hovered
as 1e30 of mounding, and a section drew both kinds as points 1e30 ft away.

Every test reads ``sentinel_run`` (tests/conftest.py): one real MF6 run with a dry
patch, an idomain-0 block and heads on both sides of 10,000 ft. What it holds is
asserted here, not assumed.
"""

from __future__ import annotations

import numpy as np
import pytest

#: Crosses the idomain-0 block in the middle and the dry north-east corner.
DIAGONAL = [(10.0, 10.0), (590.0, 590.0)]


def _heads(run):
    """``(raw, blank)``: the head file as MF6 wrote it, and where it holds no head."""

    raw = np.asarray(run.hds.array(masked=False), dtype=float)
    blank = np.isnan(np.asarray(run.hds.array(), dtype=float))
    assert (raw == -1e30).any(), "the fixture must hold a DRY cell"
    assert (raw == 1e30).any(), "the fixture must hold an INACTIVE cell"
    assert np.array_equal(blank, np.abs(raw) >= 1e29)
    return raw, blank


def _numbers(customdata):
    """Every numeric value in a trace's ``customdata``, however it is nested."""

    for row in customdata:
        for value in row if isinstance(row, (list, tuple, np.ndarray)) else [row]:
            try:
                yield float(value)
            except (TypeError, ValueError):
                continue


def test_a_head_map_is_blank_exactly_where_the_head_file_has_no_head(sentinel_run):
    """The mapped values, both hovers and the drawn trace agree with ``hds.array()``.

    The map blanked only ``== 1e30``, so a dry cell mapped and hovered as -1e30.
    The sectioned hover blanked the same one sign. The legacy hover blanked neither.
    """

    _, blank = _heads(sentinel_run)
    expected = sentinel_run.hds.array()
    choro = sentinel_run.hds.map()

    np.testing.assert_array_equal(np.asarray(choro.linear_zs, dtype=float), expected)
    context = choro._build_hover_context()
    np.testing.assert_array_equal(np.asarray(context.layer_fields["head"][0], float), expected)
    np.testing.assert_array_equal(np.asarray(choro.hover_dict["Layer 1 Head"], float), expected)

    trace = choro.fig.data[0]
    assert np.array_equal(np.isnan(np.asarray(trace.z, dtype=float)), blank)
    assert not [v for v in _numbers(trace.customdata) if abs(v) >= 1e29], (
        "the rendered hover still carries a sentinel"
    )


def test_a_mounding_map_reads_zero_at_a_dry_cell_and_blank_at_an_inactive_one(sentinel_run):
    """The user's call (ledger 315). A DRY cell has no water above the datum, so it
    reads 0, like every cell whose head sits below the datum. An INACTIVE cell has
    no head, so it is blank.

    The map, the hover and ``_mounding_series`` must give one number. The map used
    to compute its own and blanked +1e30, while the hover said 1e30.
    """

    raw, _ = _heads(sentinel_run)
    bottom = np.asarray(sentinel_run.gwf.modelgrid.botm[0], dtype=float)
    wet = np.abs(raw) < 1e29
    expected = np.full(raw.shape, np.nan)
    expected[wet] = np.clip(raw[wet] - bottom[wet], 0.0, None)
    expected[raw == -1e30] = 0.0
    choro = sentinel_run.hds.map(show_mounding=True)

    np.testing.assert_allclose(np.asarray(choro.linear_zs, dtype=float), expected)
    np.testing.assert_allclose(np.asarray(choro._mounding_series(), dtype=float), expected)
    payload = choro._build_hover_context().payload[choro._mounding_label]
    np.testing.assert_allclose(np.asarray(payload, dtype=float), np.round(expected, 2))


def test_a_depth_below_ground_map_is_blank_wherever_the_layer_has_no_head(sentinel_run):
    """The user's call (ledger 315). With no head in the layer there is no depth to
    water to show, for dry and inactive cells alike, as on the plain head map. The
    map drew -1e30 and +1e30, so every real cell was one colour."""

    raw, blank = _heads(sentinel_run)
    top = np.asarray(sentinel_run.gwf.modelgrid.top, dtype=float)
    expected = np.where(blank, np.nan, raw - top)
    choro = sentinel_run.hds.map(bgs=True)

    np.testing.assert_allclose(np.asarray(choro.linear_zs, dtype=float), expected)
    np.testing.assert_allclose(
        np.asarray(choro.hover_dict["Layer 1 Below Ground"], dtype=float), expected
    )


@pytest.mark.parametrize("contours", [True, "heads"])
def test_contours_still_draw_on_a_run_with_dry_cells(sentinel_run, contours):
    """A dry cell put -1e30 in the contoured values, so the levels ran from -1e30 in
    steps of 10. numpy refused the array ("Maximum allowed size exceeded"), the
    error was logged at DEBUG, and the map drew no contour at all."""

    _, blank = _heads(sentinel_run)
    choro = sentinel_run.hds.map(contours=contours)

    values = np.asarray(choro._contour_vector(), dtype=float)
    assert np.array_equal(np.isnan(values), blank)
    assert [trace for trace in choro.fig.data[1:] if trace.type == "scattermap"], (
        "no contour was drawn"
    )


def test_every_head_map_colours_the_real_heads_only(sentinel_run):
    """The static map's norm and the shared colour axis of ``mosaic`` and
    ``animate`` all read the mapped values, and all ran from -1e30."""

    raw, blank = _heads(sentinel_run)
    lo, hi = raw[~blank].min(), raw[~blank].max()

    static = sentinel_run.hds.map(backend="mpl")
    norms = [
        collection.norm
        for axes in static.axes
        for collection in axes.collections
        if collection.get_array() is not None
    ]
    assert norms, "the static map drew no coloured cells"
    for norm in norms:
        assert lo - 1e-6 <= norm.vmin and norm.vmax <= hi + 1e-6, (norm.vmin, norm.vmax)
    for figure in (sentinel_run.hds.mosaic(per=[0]), sentinel_run.hds.animate()):
        axis = figure.layout.coloraxis
        assert lo - 1e-6 <= axis.cmin and axis.cmax <= hi + 1e-6, (axis.cmin, axis.cmax)


def test_a_head_section_leaves_a_gap_where_there_is_no_head(sentinel_run):
    """``section()`` without ``interpolator=`` drew each crossed cell's raw head:
    a dry cell as a point at -1e30 and an inactive one at +1e30. With
    ``interpolator=`` it samples the head surface, which was already right."""

    from myflopy.modflow.utils.datatypes.xsections import XSection

    raw, blank = _heads(sentinel_run)
    section = XSection(model=sentinel_run, line=DIAGONAL, kstpkper=(0, 0))
    crossed = section.xs.index.to_numpy()
    assert (raw[crossed] == -1e30).any(), "the line must cross a DRY cell"
    assert (raw[crossed] == 1e30).any(), "the line must cross an INACTIVE cell"

    _, profile = section.xsect
    assert np.array_equal(np.isnan(np.asarray(profile[0], dtype=float)), blank[crossed])
    drawn = np.asarray(sentinel_run.hds.section(line=DIAGONAL).data[0].y, dtype=float)
    np.testing.assert_array_equal(drawn, sentinel_run.hds.array()[crossed])

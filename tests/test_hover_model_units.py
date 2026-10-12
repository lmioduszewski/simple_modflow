"""Every hover label reads the model's own length and time units (ledger 434).

``lak_hover``/``sfr_hover`` hardcoded ``ft``/``ft²``/``ft³/d``/``ft²/d`` on their
fields, the group diff maps ``ft/d``/``ft²/d``/``ft``, and the stage map ``ft`` --
so a model in metres was labelled in feet, beside a footer written in ``m²``
(ledger 401). The integration tests draw the canonical run REDECLARED in metres:
MF6 computes nothing from the unit, so its outputs are the canonical run's
(``docs/reproducers/mf6_length_units.py``) and only the labels may move.
"""

from __future__ import annotations

import re
import shutil
import subprocess

import pytest

import myflopy as mf
from myflopy.modflow.mf6.hydraulics import unit_labels
from myflopy.modflow.utils.datatypes.hover import (
    LayerTable,
    dimension_units,
    lak_hover,
    sfr_hover,
    trans_hover,
)

# --- the one rule ------------------------------------------------------------ #


def test_dimension_units_writes_each_field_from_the_two_units():
    dims = {"stage": (1, 0), "area": (2, 0), "q": (3, -1), "per_len": (2, -1), "flux": (1, -1)}
    assert dimension_units("m", "d", dims) == {
        "stage": "m", "area": "m²", "q": "m³/d", "per_len": "m²/d", "flux": "m/d",
    }
    # A unit the model does not declare is never assumed: those fields carry none.
    assert dimension_units("m", None, dims) == {"stage": "m", "area": "m²"}
    assert dimension_units(None, "d", dims) == {}
    with pytest.raises(ValueError, match="length\\^4"):
        dimension_units("m", "d", {"volume_squared": (4, 0)})


def test_lak_and_sfr_hovers_write_the_units_they_are_given():
    assert lak_hover(length="m", time="d").units == {
        "q_per_area": "m/d", "q": "m³/d", "flow_area": "m²", "stage": "m",
    }
    assert sfr_hover(length="cm", time="s").units == {
        "q_per_length": "cm²/s", "q": "cm³/s", "rlen": "cm", "stage": "cm", "depth": "cm",
    }
    assert lak_hover().units == {} and sfr_hover().units == {}
    # The transmissivity hover reads through the same rule, unchanged.
    assert trans_hover(length="ft", time="d").units == {
        "T": "ft²/d", "K": "ft/d", "thickness": "ft",
    }
    assert trans_hover(length="ft").units == {"thickness": "ft"}


def test_a_layer_table_has_no_layers_field_that_nothing_reads():
    """``LayerTable.layers`` was accepted and drew the same table whatever it said."""

    with pytest.raises(TypeError, match="layers"):
        LayerTable(fields=("head",), layers="active")


# --- drawn maps of a model declared in metres ------------------------------- #


@pytest.fixture(scope="module")
def metres_canonical(canonical_run, tmp_path_factory):
    """The canonical run's files, its DISV redeclared ``METERS`` and run again."""

    if shutil.which("mf6") is None:
        pytest.skip("needs the mf6 executable")
    source = canonical_run.sim.sim_path
    workspace = tmp_path_factory.mktemp("metres") / "canonical"
    shutil.copytree(source, workspace)
    for disv in workspace.glob("*.disv"):
        disv.write_text(re.sub(r"LENGTH_UNITS\s+\S+", "LENGTH_UNITS  meters", disv.read_text()))
    run = subprocess.run(["mf6"], cwd=workspace, capture_output=True, text=True, check=False)
    assert "Normal termination" in run.stdout, run.stdout[-1500:]
    return workspace


def _units_written(trace) -> str:
    return trace.hovertemplate + "|" + "|".join(str(v) for row in trace.customdata for v in row)


@pytest.mark.slow
def test_the_lak_and_sfr_exchange_maps_are_labelled_in_metres(metres_canonical):
    model = mf.load_mf6_run(str(metres_canonical))
    assert unit_labels(model) == ("m", "d")
    for picture, expected in (
        (model.lak.results.q.map(per=0), ("m/d", "m³/d", "m²")),
        (model.sfr.results.q.map(per=0), ("m²/d", "m³/d")),
        (model.lak.results.stage.map(per=0), ("\xa0m",)),
    ):
        written = _units_written(picture.get_choropleth())
        assert "ft" not in written
        for unit in expected:
            assert unit in written, unit


@pytest.mark.slow
def test_the_group_maps_are_labelled_in_metres(metres_canonical):
    first = mf.load_mf6_run(str(metres_canonical))
    second = mf.load_mf6_run(str(metres_canonical))
    group = mf.ModelGroup({"a": first, "b": second}, reference="a")
    pictures = [
        group.packages.lak.results.q.map(per=0),
        group.packages.sfr.results.q.map(per=0),
        group.diff().packages.lak.results.q.map("b", per=0),
        group.diff().packages.sfr.results.q.map("b", per=0),
        group.diff().hds.map("b", per=0),
    ]
    for picture in pictures:
        written = _units_written(picture.get_choropleth())
        assert "ft" not in written
        # a value is written "<number>\xa0<unit>"
        assert re.search(r"\xa0m[²³]?(/d)?(?![A-Za-z])", written), written[:400]


@pytest.mark.slow
def test_every_cross_section_axis_is_labelled_in_metres(metres_canonical):
    """``ModelCrossSectionStyle`` and ``XSection.plot_mpl`` wrote ``(ft)`` whatever
    the model declared (ledger 464)."""

    import matplotlib.pyplot as plt
    from shapely.geometry import LineString

    from myflopy.modflow.mf6.cross_section_plotting import plot_model_cross_section

    model = mf.load_mf6_run(str(metres_canonical))
    xmin, ymin, xmax, ymax = model.vor.gdf_vorPolys.total_bounds
    line = LineString([(xmin, (ymin + ymax) / 2), (xmax, (ymin + ymax) / 2)])
    try:
        axes = [
            model.plot.section(line=line).axes[0],
            plot_model_cross_section(model, line)[1],
            model.hds._sections(line=line)[model.name].plot_mpl()[1],
        ]
        for ax in axes:
            assert (ax.get_xlabel(), ax.get_ylabel()) == ("Distance (m)", "Elevation (m)")
        layout = model.plot.section(line=line, backend="plotly").fig.layout
        assert (layout.xaxis.title.text, layout.yaxis.title.text) == (
            "Distance (m)", "Elevation (m)",
        )
    finally:
        plt.close("all")

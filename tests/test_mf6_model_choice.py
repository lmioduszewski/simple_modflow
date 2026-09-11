"""``load_mf6_run`` opens the FLOW model, not the alphabetically first.

The model name was the stem of the first ``*.nam`` in sorted order, so in a coupled
GWF+GWT/GWE/PRT run whose transport model sorted first, the loader opened that
model as if it were flow. ``mfsim.nam``'s models block names each model's type.
"""

from __future__ import annotations

from myflopy.project.run_model import _infer_model_name


def _workspace(tmp_path, models_block: str, stems: tuple[str, ...]):
    (tmp_path / "mfsim.nam").write_text(f"BEGIN models\n{models_block}END models\n")
    for stem in stems:
        (tmp_path / f"{stem}.nam").write_text("BEGIN packages\nEND packages\n")
    return tmp_path


def test_the_flow_model_is_chosen_even_when_transport_sorts_first(tmp_path):
    ws = _workspace(tmp_path, "  gwt6  atrans.nam  atrans\n  gwf6  zflow.nam  zflow\n", ("atrans", "zflow"))
    assert _infer_model_name(ws) == "zflow"


def test_a_model_declared_without_a_name_takes_its_name_file_stem(tmp_path):
    ws = _workspace(tmp_path, "  gwt6  atrans.nam  atrans\n  gwf6  zflow.nam\n", ("atrans", "zflow"))
    assert _infer_model_name(ws) == "zflow"


def test_without_mfsim_the_name_file_still_decides(tmp_path):
    (tmp_path / "only.nam").write_text("BEGIN packages\nEND packages\n")
    assert _infer_model_name(tmp_path) == "only"

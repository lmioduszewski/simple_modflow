from __future__ import annotations

from simple_modflow.project.manifest_io import read_toml, write_toml


def test_toml_roundtrip_preserves_nested_none_values(tmp_path):
    path = tmp_path / "nested_none.toml"
    payload = {
        "package_data": {
            "rows": [
                [0, 1.0, None],
                [1, None, ["nested", None]],
            ]
        }
    }

    write_toml(payload, path)
    restored = read_toml(path)

    assert restored == payload

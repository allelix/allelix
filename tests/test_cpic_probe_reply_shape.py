import io
import json

import pytest

from allelix.databases import cpic_loader


@pytest.mark.parametrize("payload", [{"message": "upstream error"}, 1, True])
def test_cpic_probe_bad_json_shape_is_unavailable(monkeypatch, payload):
    monkeypatch.setattr(
        cpic_loader.urllib.request,
        "urlopen",
        lambda *a, **k: io.BytesIO(json.dumps(payload).encode()),
    )
    assert cpic_loader.fetch_cpic_remote_signal("https://synthetic.invalid") is None

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""A usable padded placement beats a no-call at the same normalized site."""

import pytest
from click.testing import CliRunner
from test_extract_reference_calls import _input, _row
from validated_call_helpers import BUILDS, RSID, setup

from allelix.cli import main, utility


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("indexed", [False, True], ids=["sequential", "indexed"])
@pytest.mark.parametrize("named", [False, True], ids=["unnamed", "named"])
def test_usable_padded_call_wins_over_mapped_no_call(tmp_path, monkeypatch, build, indexed, named):
    data = tmp_path / "data"
    pos = setup(data, build)
    found = []
    monkeypatch.setattr(
        utility, "_render_extract_table", lambda path, wanted, rows: found.append(rows)
    )
    path = _input(
        tmp_path,
        build,
        [
            _row(pos - 1, "0/1", name=RSID if named else ".", ref="CG", alt="CA"),
            _row(pos, "./.", name=RSID),
        ],
        indexed=indexed,
    )
    result = CliRunner().invoke(
        main,
        ["extract", str(path), "--snps", RSID, "--build", build.lower(), "--data-dir", str(data)],
    )
    assert result.exit_code == 0, (result.output, result.exception)
    assert RSID in found[0]
    chosen = found[0][RSID]
    assert not chosen.is_no_call
    assert (chosen.position, chosen.ref, chosen.allele1, chosen.allele2) == (
        pos - 1,
        "CG",
        "CG",
        "CA",
    )

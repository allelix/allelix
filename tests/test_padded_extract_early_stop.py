# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Extract's early stop must agree with its normalized placement choice."""

import re

import pytest
from click.testing import CliRunner
from test_extract_reference_calls import _input, _row
from validated_call_helpers import BUILDS, RSID, setup

from allelix.cli import main, utility


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("named", [False, True], ids=["unnamed", "named"])
@pytest.mark.parametrize("lengths", [False, True], ids=["contig-order-only", "contig-lengths"])
@pytest.mark.parametrize("padded_first", [False, True], ids=["plain-first", "padded-first"])
def test_early_stop_agrees_with_normalized_placement_and_file_order(
    tmp_path, monkeypatch, build, named, lengths, padded_first
):
    data = tmp_path / "data"
    pos = setup(data, build)
    plain = _row(pos, "0/1", name=RSID)
    padded = _row(pos - 1, "1/1", name=RSID if named else ".", ref="CG", alt="CA")
    rows = [padded if padded_first else plain]
    # Go beyond resolver read-ahead. Both calls stand at the same placement;
    # their file order breaks the tie, as if both were written trimmed.
    rows.extend(_row(pos + 10 + i, "0/1", name="rs900260382") for i in range(1100))
    rows.append(plain if padded_first else padded)
    _input(tmp_path, build, rows, indexed=False)
    path = tmp_path / "input.vcf"
    if not lengths:
        path.write_text(re.sub(r",length=\d+", "", path.read_text()))
    found = []
    monkeypatch.setattr(
        utility, "_render_extract_table", lambda path, wanted, rows: found.append(rows)
    )
    result = CliRunner().invoke(
        main,
        ["extract", str(path), "--snps", RSID, "--build", build.lower(), "--data-dir", str(data)],
    )
    assert result.exit_code == 0, (result.output, result.exception)
    chosen = found[0][RSID]
    expected = (pos - 1, "CG", "CA", "CA") if padded_first else (pos, "G", "G", "A")
    assert (chosen.position, chosen.ref, chosen.allele1, chosen.allele2) == expected

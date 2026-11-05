# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Padding and missing IDs do not hide extract's reference or rejected calls."""

import pytest
from click.testing import CliRunner
from test_extract_reference_calls import _input, _row
from validated_call_helpers import BUILDS, RSID, setup

from allelix.cli import main, utility


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("indexed", [False, True], ids=["sequential", "indexed"])
@pytest.mark.parametrize("shape", ["prefix", "suffix", "both"])
@pytest.mark.parametrize(
    "gt,filt,gq",
    [
        ("./.", "PASS", 50),
        (".", "PASS", 50),
        ("0/.", "PASS", 50),
        ("0/0", "PASS", 50),
        ("0", "PASS", 50),
        ("0/1", "LowQual", 50),
        ("0/1", "PASS", 1),
    ],
    ids=["no-call", "haploid-no-call", "half-ref", "hom-ref", "haploid-ref", "filtered", "low-gq"],
)
def test_extract_keeps_padded_reference_and_rejected_calls(
    tmp_path, monkeypatch, build, indexed, shape, gt, filt, gq
):
    data = tmp_path / "data"
    pos = setup(data, build)
    ref, alt = {"prefix": ("CG", "CA"), "suffix": ("GC", "AC"), "both": ("CGC", "CAC")}[shape]
    start = pos - (shape in ("prefix", "both"))
    captured = []
    monkeypatch.setattr(
        utility, "_render_extract_table", lambda path, wanted, found: captured.append(found)
    )
    for named in (True, False):
        root = tmp_path / ("named" if named else "unnamed")
        root.mkdir()
        path = _input(
            root,
            build,
            [_row(start, gt, name=RSID if named else ".", ref=ref, alt=alt, filt=filt, gq=gq)],
            indexed=indexed,
        )
        result = CliRunner().invoke(
            main,
            [
                "extract",
                str(path),
                "--snps",
                RSID,
                "--build",
                build.lower(),
                "--data-dir",
                str(data),
            ],
        )
        assert result.exit_code == 0, (result.output, result.exception)
        assert RSID in captured[-1], (named, result.output)
    named, unnamed = (found[RSID] for found in captured)
    assert (unnamed.chromosome, unnamed.position, unnamed.ref) == (
        named.chromosome,
        named.position,
        named.ref,
    )
    assert (unnamed.allele1, unnamed.allele2, unnamed.haploid) == (
        named.allele1,
        named.allele2,
        named.haploid,
    )

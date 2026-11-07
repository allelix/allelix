# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""A gVCF's reference-block ALT does not change a record's identity in compare.

A gVCF writes ``<NON_REF>`` (or ``<*>``) as an extra ALT on every variant
record. It names no allele the record carries, so ``rs1 A G,<NON_REF>`` in a
gVCF and ``rs1 A G`` in a plain VCF are the same record and are compared.
"""

from __future__ import annotations

import re

import pytest
from click.testing import CliRunner

from allelix.cli import main
from allelix.compare import compare_variants
from allelix.parsers.vcf import VcfParser


def _vcf(tmp_path, name, rows, gvcf=False):
    path = tmp_path / name
    path.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=1,length=248956422>\n"
        + ('##ALT=<ID=NON_REF,Description="Any other allele">\n' if gvcf else "")
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        + "".join(
            f"1\t{pos}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n"
            for pos, rsid, ref, alt, gt in rows
        )
    )
    return path, list(VcfParser().parse(path))


@pytest.mark.parametrize("symbol", ["<NON_REF>", "<*>"])
@pytest.mark.parametrize(
    "ref,alt,gt",
    [
        ("A", "G", "0/1"),
        ("A", "G", "1/1"),
        ("AT", "A", "0/1"),
        ("A", "AT", "1/1"),
        ("A", "G,C", "1/2"),
    ],
)
@pytest.mark.parametrize("named", [True, False])
@pytest.mark.parametrize("swap", [False, True])
def test_gvcf_record_pairs_with_the_plain_record(tmp_path, symbol, ref, alt, gt, named, swap):
    rsid = "rs3999980401" if named else "."
    _, gvcf = _vcf(tmp_path, "g.vcf", [(100, rsid, ref, f"{alt},{symbol}", gt)], gvcf=True)
    _, plain = _vcf(tmp_path, "p.vcf", [(100, rsid, ref, alt, gt)])
    left, right = (plain, gvcf) if swap else (gvcf, plain)
    result = compare_variants(left, right, build1="GRCh38", build2="GRCh38")
    assert (result.shared, result.concordant) == (1, 1)
    assert (result.file1_only, result.file2_only) == (0, 0)
    assert (result.shared_named, result.shared_unnamed) == ((1, 0) if named else (0, 1))


@pytest.mark.parametrize("named", [True, False])
def test_gvcf_record_with_another_real_alt_stays_distinct(tmp_path, named):
    # Control: only the reference-block symbol is ignored; a different real ALT
    # is still a different record.
    rsid = "rs3999980402" if named else "."
    _, gvcf = _vcf(tmp_path, "g.vcf", [(100, rsid, "A", "G,<NON_REF>", "0/1")], gvcf=True)
    _, plain = _vcf(tmp_path, "p.vcf", [(100, rsid, "A", "C", "0/1")])
    result = compare_variants(gvcf, plain, build1="GRCh38", build2="GRCh38")
    assert result.shared == 0
    assert (result.file1_only, result.file2_only) == (1, 1)


def test_gvcf_record_and_its_discordant_plain_record_are_compared(tmp_path):
    _, gvcf = _vcf(
        tmp_path, "g.vcf", [(100, "rs3999980403", "A", "G,<NON_REF>", "0/1")], gvcf=True
    )
    _, plain = _vcf(tmp_path, "p.vcf", [(100, "rs3999980403", "A", "G", "1/1")])
    result = compare_variants(gvcf, plain, build1="GRCh38", build2="GRCh38")
    assert (result.shared, result.discordant) == (1, 1)


def test_cli_compares_a_gvcf_with_a_plain_vcf(tmp_path):
    g, _ = _vcf(tmp_path, "g.vcf", [(100, "rs3999980404", "A", "G,<NON_REF>", "0/1")], gvcf=True)
    p, _ = _vcf(tmp_path, "p.vcf", [(100, "rs3999980404", "A", "G", "0/1")])
    result = CliRunner().invoke(main, ["compare", str(g), str(p)])
    assert result.exit_code == 0, result.output
    assert re.search(r"Concordant\D*1\b", result.output), result.output
    assert not re.search(r"only in file [12]\D*[1-9]", result.output, re.I), result.output

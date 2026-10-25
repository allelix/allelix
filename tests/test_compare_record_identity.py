# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Literal input identity and truthful comparison coverage."""

from __future__ import annotations

import gzip
import json
from dataclasses import replace

import pytest
from click.testing import CliRunner

from allelix.cli import main
from allelix.compare import compare_variants
from allelix.models import Variant
from allelix.parsers.vcf import VcfParser, _parse_data_line
from allelix.reports.json_report import SCHEMA_VERSION, render_json
from tests.test_original_rsid_never_serializes import _result_with_sentinel


def _header(build="GRCh38"):
    length = 249250621 if build == "GRCh37" else 248956422
    return (
        f"##fileformat=VCFv4.2\n##contig=<ID=1,length={length}>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    )


def _row(pos=100, rsid=".", ref="A", alt="G", gt="0/1", chrom="1", filt="PASS"):
    return f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t.\t{filt}\t.\tGT:GQ\t{gt}\n"


def _parse(tmp_path, name, rows, build="GRCh38", compressed=False):
    path = tmp_path / name
    text = _header(build) + "".join(rows)
    if compressed:
        path = path.with_suffix(".vcf.gz")
        path.write_bytes(gzip.compress(text.encode(), mtime=0))
    else:
        path.write_text(text)
    return path, list(VcfParser().parse(path))


@pytest.mark.parametrize("compressed", [False, True])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("shared", [False, True])
def test_distinct_unnamed_records_do_not_collapse(tmp_path, compressed, reverse, shared):
    left = [_row(100), _row(200, alt="C")]
    right = left[:] if shared else [_row(300), _row(400)]
    _, a = _parse(tmp_path, "a.vcf", left[::-1] if reverse else left, compressed=compressed)
    _, b = _parse(tmp_path, "b.vcf", right, compressed=compressed)
    result = compare_variants(a, b, build1="GRCh38", build2="GRCh38")
    assert (result.file1_total, result.file2_total) == (2, 2)
    assert (result.shared, result.concordant) == ((2, 2) if shared else (0, 0))
    assert (result.file1_only, result.file2_only) == ((0, 0) if shared else (2, 2))
    assert result.shared_named == 0
    assert result.shared_unnamed == (2 if shared else 0)
    assert (result.file1_input_unnamed, result.file1_unnamed) == (2, 2)


@pytest.mark.parametrize(
    "left,right,shared,no_call",
    [
        (_row(ref="A"), _row(ref="C"), 0, 0),
        (_row(alt="G"), _row(alt="C"), 0, 0),
        (_row(alt="G,C"), _row(alt="G,T"), 0, 0),
        (_row(alt="G,C"), _row(alt="C,G", gt="0/2"), 0, 0),
        (_row(alt="G", gt="./."), _row(alt="C", gt="./."), 0, 0),
        (_row(alt="G", gt="0/1:3"), _row(alt="C", gt="0/1:3"), 0, 0),
        (_row(gt="./."), _row(gt="0/1"), 1, 1),
        (_row(gt="0/1:3"), _row(gt="0/1:99"), 1, 1),
        (_row(ref="a", alt="g,c"), _row(ref="A", alt="G,C"), 1, 0),
        (_row(chrom="chr1"), _row(chrom="1"), 1, 0),
        (_row(rsid="1:100:A:G"), _row(), 1, 0),
    ],
)
def test_full_literal_allele_identity_is_independent_of_called_genotype(
    tmp_path, left, right, shared, no_call
):
    _, a = _parse(tmp_path, "a.vcf", [left])
    _, b = _parse(tmp_path, "b.vcf", [right])
    result = compare_variants(a, b, build1="GRCh38", build2="GRCh38")
    assert result.shared == shared
    assert result.no_call == no_call
    assert result.file1_only == result.file2_only == 1 - shared


def test_same_literal_record_with_different_genotype_is_compared(tmp_path):
    _, a = _parse(tmp_path, "a.vcf", [_row(gt="0/1")])
    _, b = _parse(tmp_path, "b.vcf", [_row(gt="1/1")])
    result = compare_variants(a, b, build1="GRCh38", build2="GRCh38")
    assert (result.shared, result.discordant) == (1, 1)


@pytest.mark.parametrize(
    "builds", [("GRCh37", "GRCh38"), ("unknown", "unknown"), ("GRCh38", "unknown")]
)
def test_unnamed_omissions_preserve_named_identifier_matching(tmp_path, builds):
    _, a = _parse(tmp_path, "a.vcf", [_row(), _row(200, "rs900400001")])
    _, b = _parse(tmp_path, "b.vcf", [_row(), _row(300, "rs900400001")])
    result = compare_variants(a, b, build1=builds[0], build2=builds[1])
    assert (result.file1_total, result.file2_total, result.shared) == (1, 1, 1)
    assert result.shared_named == 1 and result.shared_unnamed == 0
    assert result.file1_unnamed_omitted == result.file2_unnamed_omitted == 1
    assert result.file1_input_named == result.file1_input_unnamed == 1
    assert (result.file1_named, result.file1_unnamed) == (1, 0)


def test_named_and_unnamed_namespaces_remain_separate(tmp_path):
    _, a = _parse(tmp_path, "a.vcf", [_row(rsid="rs900400001")])
    _, b = _parse(tmp_path, "b.vcf", [_row()])
    assert compare_variants(a, b).shared == 0
    vendor = Variant("i5000123", "1", 100, "A", "G")
    result = compare_variants([vendor], [replace(vendor, position=200)])
    assert result.shared_named == result.concordant == 1


@pytest.mark.parametrize("missing", ["alt", "ref", "input-ref"])
def test_insufficient_unnamed_identity_is_omitted(tmp_path, missing):
    _, calls = _parse(tmp_path, "a.vcf", [_row()])
    variant = calls[0]
    if missing == "alt":
        variant.input_alt = None
    elif missing == "ref":
        variant.ref = None
    else:
        variant.ref_from_input = False
    result = compare_variants([variant], [variant], build1="GRCh38", build2="GRCh38")
    assert result.shared == result.file1_total == result.file2_total == 0
    assert result.file1_unnamed_omitted == result.file2_unnamed_omitted == 1


def test_duplicate_keys_keep_existing_last_record_convention(tmp_path):
    _, calls = _parse(tmp_path, "a.vcf", [_row(gt="0/1"), _row(gt="1/1")])
    first, second = calls
    result = compare_variants([first, second], [first])
    assert (result.file1_total, result.file1_input_unnamed, result.discordant) == (1, 2, 1)
    reversed_result = compare_variants([second, first], [first])
    assert reversed_result.file1_total == reversed_result.concordant == 1


def test_cli_reports_input_eligible_and_positional_counts(tmp_path):
    a, _ = _parse(tmp_path, "a.vcf", [_row(), _row(200), _row(300, "rs900400001")])
    b, _ = _parse(tmp_path, "b.vcf", [_row(), _row(400), _row(500, "rs900400001")])
    result = CliRunner(env={"COLUMNS": "140"}).invoke(main, ["compare", str(a), str(b)])
    assert result.exit_code == 0, result.output
    text = " ".join(result.output.lower().split())
    assert "input named / unnamed (file 1)" in text
    assert "eligible named / unnamed (file 1)" in text
    assert "unnamed variants compared by position" in text
    assert "shared named identifiers" in text
    assert "shared eligible variants" in text
    assert "100.00%" in text


def test_cli_all_unnamed_different_builds_has_no_false_concordance(tmp_path):
    a, _ = _parse(tmp_path, "a.vcf", [_row()], build="GRCh37")
    b, _ = _parse(tmp_path, "b.vcf", [_row()], build="GRCh38")
    result = CliRunner(env={"COLUMNS": "140"}).invoke(main, ["compare", str(a), str(b)])
    assert result.exit_code == 0, result.output
    text = " ".join(result.output.lower().split())
    assert "omitted unnamed records (file 1)" in text
    assert "omitted unnamed records (file 2)" in text
    assert "100.00%" not in text and "100.0%" not in text and "100%" not in text
    assert "position-based comparisons may be unreliable" not in text


def test_input_alt_plumbing_preserves_model_equality_and_public_serialization(tmp_path):
    variant = _parse_data_line(_row(alt="g,c").rstrip(), 0)
    assert variant is not None
    assert getattr(variant, "input_alt", None) == "G,C"
    assert variant == replace(variant, input_alt="G,T")
    legacy = Variant("rs1", "1", 100, "A", "G", "GRCh38", "A", True, None, False)
    assert legacy.input_alt is None and legacy.haploid is False
    result = _result_with_sentinel(tmp_path / "source.vcf")
    result.hv_variants[0].input_alt = "G,C"
    out = tmp_path / "report.json"
    render_json(result, output_path=out)
    payload = json.loads(out.read_text())
    assert payload["schema_version"] == SCHEMA_VERSION
    assert '"input_alt"' not in out.read_text()

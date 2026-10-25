# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Source identity, record type and ploidy of reconstructed reference calls."""

from __future__ import annotations

import contextlib
import gzip
import json
import sqlite3
import sys
import zipfile
from collections.abc import Mapping
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner

from allelix.annotators.clinvar import ClinVarAnnotator, clinvar_db_filename
from allelix.annotators.rsid_identity import RsidIdentityAnnotator
from allelix.cli import main
from allelix.databases.rsid_identity_loader import RSID_IDENTITY_DB_FILENAME, write_identity_db
from allelix.parsers.vcf import VcfParser, _LocusWatcher
from allelix.reports import vcf_quality_lines
from allelix.reports._pipeline import run_analysis
from allelix.reports.diff import load_previous_report
from allelix.reports.json_report import SCHEMA_VERSION, render_json
from tests.test_clinvar_matching_batch import _write_clinvar

if TYPE_CHECKING:
    from pathlib import Path

_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=2,length=243199373>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
)


@pytest.fixture
def data(tmp_path: Path) -> Path:
    directory = tmp_path / "data"
    directory.mkdir()
    _write_clinvar(directory)
    with contextlib.closing(sqlite3.connect(directory / clinvar_db_filename("GRCh37"))) as conn:
        conn.executemany(
            "INSERT INTO clinvar_variants (rsid, chromosome, position, ref, alt, "
            "clinical_significance, condition, gene, review_status, allele_id, aggregate, "
            "submitters_same, submitters_total) VALUES (?, ?, 1000, 'G', 'G', 'Pathogenic', "
            "'Synthetic condition', 'GENE', 'criteria_provided,_single_submitter', ?, "
            "'Pathogenic', 1, 1)",
            [("rs900230002", "X", 20), ("rs900230003", "MT", 21)],
        )
        conn.commit()
    return directory


def _write(path: Path, rows: list[str], container: str = "plain", extra: str = "") -> Path:
    text = _HEADER.replace("#CHROM", extra + "#CHROM") + "\n".join(rows) + "\n"
    if container == "gz":
        path = path.with_suffix(".vcf.gz")
        path.write_bytes(gzip.compress(text.encode(), mtime=0))
    elif container == "zip":
        path = path.with_suffix(".zip")
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr("synthetic.vcf", text)
    else:
        path.write_text(text)
    return path


def _row(
    pos: int = 2000,
    ref: str = "T",
    alt: str = "C",
    gt: str = "0/0",
    info: str = ".",
    chrom: str = "2",
) -> str:
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\tGT:GQ\t{gt}:99"


@pytest.mark.parametrize("container", ["plain", "gz", "zip"])
@pytest.mark.parametrize("ref", ["A", "C", "t", "T"])
def test_reference_call_uses_the_input_ref(
    data: Path, tmp_path: Path, container: str, ref: str
) -> None:
    path = _write(tmp_path / "ref.vcf", [_row(ref=ref)], container)
    result = run_analysis(path, VcfParser(), [ClinVarAnnotator(data)])
    assert result.vcf_quality is not None
    if ref.upper() == "T":
        assert [(a.rsid, a.genotype_match) for a in result.annotations] == [("rs900230001", "TT")]
        assert result.vcf_quality.reference_carried_calls == 1
        assert result.vcf_quality.clinvar_ref_mismatches == 0
    else:
        assert result.annotations == []
        assert result.vcf_quality.reference_carried_calls == 0
        assert result.vcf_quality.clinvar_ref_mismatches == 1
        assert "1 record(s) at a ClinVar position have a REF" in " ".join(
            vcf_quality_lines(result.vcf_quality)
        )


@pytest.mark.parametrize(
    "alt",
    ["<DEL>", "<DUP>", "<DUP:TANDEM>", "<INV>", "<INS>", "<CNV>", "G]3:1000]", "<DEL>,<NON_REF>"],
)
@pytest.mark.parametrize("start", [1990, 2000])
def test_symbolic_spans_never_reconstruct_reference_calls(
    data: Path, tmp_path: Path, alt: str, start: int
) -> None:
    path = _write(
        tmp_path / "sv.vcf", [_row(1000, "G", "A", "0/1"), _row(start, "T", alt, info="END=2010")]
    )
    parser = VcfParser()
    result = run_analysis(path, parser, [ClinVarAnnotator(data)])
    assert not any(a.rsid == "rs900230001" for a in result.annotations)
    assert result.vcf_quality is not None
    assert result.vcf_quality.is_gvcf is False
    assert result.vcf_quality.symbolic_rows == 1
    assert parser.vcf_type(path) == "plain"


@pytest.mark.parametrize("alt", ["<NON_REF>", "<*>", "."])
@pytest.mark.parametrize("gt", ["0", "0/0"])
def test_real_reference_blocks_keep_their_source_ploidy(
    data: Path, tmp_path: Path, alt: str, gt: str
) -> None:
    # The block anchor is a different base from the reference-carried site.
    path = _write(tmp_path / "block.vcf", [_row(1990, "A", alt, gt, "END=2010")])
    result = run_analysis(path, VcfParser(), [ClinVarAnnotator(data)])
    assert result.vcf_quality is not None and result.vcf_quality.is_gvcf
    assert len(result.annotations) == 1
    row = result.annotations[0]
    assert row.genotype_match == ("T" if gt == "0" else "TT")
    assert row.zygosity_code == ("hemizygous" if gt == "0" else "homozygous")


@pytest.mark.parametrize(
    "name,length,skipped",
    [("chrM", 16571, 2), ("M", None, 2), ("chrM", 16569, 0), ("MT", None, 0), ("chrMT", None, 0)],
)
def test_reference_calls_follow_mitochondrial_exclusions(
    data: Path, tmp_path: Path, name: str, length: int | None, skipped: int
) -> None:
    extra = f"##contig=<ID={name},length={length}>\n" if length else ""
    path = _write(
        tmp_path / "mt.vcf",
        [_row(2000, "G", "A", "1/1", chrom=name), _row(1000, "G", "A", chrom=name)],
        extra=extra,
    )
    result = run_analysis(path, VcfParser(), [ClinVarAnnotator(data)])
    assert result.build_diagnostics is not None
    assert result.build_diagnostics.evidence()["mt_records_skipped"] == skipped
    assert any(a.rsid == "rs900230003" for a in result.annotations) is (skipped == 0)
    assert result.vcf_quality is not None
    assert result.vcf_quality.reference_carried_calls == (0 if skipped else 1)


def test_haploid_explicit_reference_call(data: Path, tmp_path: Path) -> None:
    path = _write(tmp_path / "x.vcf", [_row(1000, "G", "A", "0", chrom="X")])
    result = run_analysis(path, VcfParser(), [ClinVarAnnotator(data)])
    (row,) = result.annotations
    assert (row.genotype_match, row.zygosity_code) == ("G", "hemizygous")


@pytest.mark.parametrize("container", ["plain", "gz", "zip"])
def test_reference_only_cli_and_stats(data: Path, tmp_path: Path, container: str) -> None:
    path = _write(tmp_path / "only-ref.vcf", [_row()], container)
    out = tmp_path / "report.json"
    result = CliRunner().invoke(
        main, ["analyze", str(path), "--data-dir", str(data), "--no-update", "--output", str(out)]
    )
    assert result.exit_code == 0, result.output
    assert "This file has no non-reference calls." in result.output
    payload = json.loads(out.read_text())
    assert [a["genotype_match"] for a in payload["annotations"]] == ["TT"]
    quality = payload["input"]["vcf_quality"]
    assert (
        quality["genotype_records"],
        quality["hom_ref_calls"],
        quality["non_reference_calls"],
    ) == (1, 1, 0)
    assert payload["schema_version"] == SCHEMA_VERSION == "10"
    html = tmp_path / "report.html"
    rendered = CliRunner().invoke(
        main, ["analyze", str(path), "--data-dir", str(data), "--no-update", "--output", str(html)]
    )
    assert rendered.exit_code == 0, rendered.output
    assert "This file has no non-reference calls." in html.read_text()
    stats = CliRunner().invoke(main, ["stats", str(path), "--data-dir", str(data)])
    assert stats.exit_code == 0, stats.output
    assert "Genotype records" in stats.output
    assert "Homozygous-reference" in stats.output
    assert "This file has no non-reference calls." in stats.output


def test_header_only_vcf_still_refuses(data: Path, tmp_path: Path) -> None:
    path = _write(tmp_path / "empty.vcf", [])
    for command in ("analyze", "stats"):
        result = CliRunner().invoke(
            main,
            [
                command,
                str(path),
                "--data-dir",
                str(data),
                *(["--no-update"] if command == "analyze" else []),
            ],
        )
        assert result.exit_code == 1
        assert "contains no genotype records" in result.output


def _identity_rows(
    data: Path,
    refs: list[str],
    *,
    same_id: bool = False,
    chrom: str = "2",
    pos37: int = 2000,
    pos38: int = 2100,
    palindromic: bool = False,
) -> None:
    write_identity_db(
        [
            (
                "rs900230001" if same_id and i == 0 else f"rs90023009{i}",
                chrom,
                pos37,
                ref,
                "A",
                pos38,
                "C",
                "G",
                int(palindromic),
            )
            for i, ref in enumerate(refs)
        ],
        data / RSID_IDENTITY_DB_FILENAME,
        version="synthetic",
        source_url="test://",
        downloaded_at="2001-01-01",
    )


@pytest.mark.parametrize("same_id", [False, True])
@pytest.mark.parametrize(
    "refs,input_ref,kept,mismatches,conflicts",
    [
        (["T"], "T", True, 0, 0),
        (["T"], "A", False, 1, 0),
        (["A"], "T", False, 0, 1),
        (["A"], "A", False, 0, 1),
        (["T", "T"], "T", True, 0, 0),
        (["T", "T"], "A", False, 1, 0),
        (["T", "A"], "T", False, 0, 1),
        (["A", "T"], "T", False, 0, 1),
        (["T", "T", "A"], "T", False, 0, 1),
        (["A", "A", "T"], "A", False, 0, 1),
        (["N"], "T", True, 0, 0),
        (["TT"], "T", True, 0, 0),
        (["-"], "T", True, 0, 0),
        ([""], "T", True, 0, 0),
        (["<DEL>"], "T", True, 0, 0),
        (["T", "N"], "A", False, 1, 0),
    ],
)
def test_coordinate_reference_provider_consistency(
    data: Path,
    tmp_path: Path,
    same_id: bool,
    refs: list[str],
    input_ref: str,
    kept: bool,
    mismatches: int,
    conflicts: int,
) -> None:
    _identity_rows(data, refs, same_id=same_id)
    path = _write(tmp_path / "identity.vcf", [_row(ref=input_ref)])
    with RsidIdentityAnnotator(data) as identity:
        result = run_analysis(path, VcfParser(), [ClinVarAnnotator(data)], identity=identity)
    assert bool(result.annotations) is kept
    if kept:
        assert result.annotations[0].rsid == "rs900230001"
    assert result.vcf_quality is not None
    assert result.vcf_quality.clinvar_ref_mismatches == mismatches
    assert getattr(result.vcf_quality, "reference_source_conflicts", 0) == conflicts
    notes = " ".join(vcf_quality_lines(result.vcf_quality))
    if conflicts:
        assert "1 reference call(s) were not assessed because reference sources disagree" in notes
        assert "REF that is not ClinVar's" not in notes


@pytest.mark.parametrize("chrom,pos37,pos38", [("3", 2000, 2100), ("2", 2100, 2000)])
def test_identity_reference_uses_only_the_input_build_and_coordinate(
    data: Path,
    tmp_path: Path,
    chrom: str,
    pos37: int,
    pos38: int,
) -> None:
    _identity_rows(data, ["A"], chrom=chrom, pos37=pos37, pos38=pos38)
    path = _write(tmp_path / "identity.vcf", [_row()])
    with RsidIdentityAnnotator(data) as identity:
        result = run_analysis(path, VcfParser(), [ClinVarAnnotator(data)], identity=identity)
    assert [a.rsid for a in result.annotations] == ["rs900230001"]
    assert result.vcf_quality is not None
    assert result.vcf_quality.clinvar_ref_mismatches == 0
    assert getattr(result.vcf_quality, "reference_source_conflicts", 0) == 0


@pytest.mark.parametrize("ref,kept", [("T", True), ("A", False)])
def test_same_build_palindromic_reference_is_usable(
    data: Path, tmp_path: Path, ref: str, kept: bool
) -> None:
    _identity_rows(data, [ref], palindromic=True)
    path = _write(tmp_path / "palindromic.vcf", [_row()])
    with RsidIdentityAnnotator(data) as identity:
        result = run_analysis(path, VcfParser(), [ClinVarAnnotator(data)], identity=identity)
    assert bool(result.annotations) is kept
    assert result.vcf_quality is not None
    assert result.vcf_quality.clinvar_ref_mismatches == 0
    assert getattr(result.vcf_quality, "reference_source_conflicts", 0) == int(not kept)


@pytest.mark.parametrize("ref,kept", [("T", True), ("A", False)])
def test_block_interior_uses_its_own_reference_sources(
    data: Path, tmp_path: Path, ref: str, kept: bool
) -> None:
    _identity_rows(data, [ref])
    path = _write(tmp_path / "interior.vcf", [_row(1990, "A", "<NON_REF>", info="END=2010")])
    with RsidIdentityAnnotator(data) as identity:
        result = run_analysis(path, VcfParser(), [ClinVarAnnotator(data)], identity=identity)
    assert bool(result.annotations) is kept
    assert result.vcf_quality is not None
    assert result.vcf_quality.clinvar_ref_mismatches == 0
    assert getattr(result.vcf_quality, "reference_source_conflicts", 0) == int(not kept)


def test_coordinate_rejection_check_does_not_grant_array_reference(data: Path) -> None:
    from allelix.annotators._resolution import _populate_ref_from_identity
    from allelix.models import Variant

    _identity_rows(data, ["T"])
    variants = [
        Variant(rsid, "2", 2000, "T", "T", build="GRCh37")
        for rsid in ("rs900230090", "rs900230001", "i5000123")
    ]
    with RsidIdentityAnnotator(data) as identity:
        _populate_ref_from_identity(identity, variants)
    assert [v.ref for v in variants] == ["T", None, None]
    assert [v.rsid for v in variants] == ["rs900230090", "rs900230001", "i5000123"]


@pytest.mark.parametrize(
    "refs,ref,mismatches,conflicts", [(["T", "T"], "A", 2, 0), (["A", "T", "T"], "T", 0, 2)]
)
def test_rejections_count_input_observations_not_alias_rows(
    data: Path, tmp_path: Path, refs: list[str], ref: str, mismatches: int, conflicts: int
) -> None:
    _identity_rows(data, refs)
    path = _write(tmp_path / "duplicates.vcf", [_row(ref=ref), _row(ref=ref)])
    with RsidIdentityAnnotator(data) as identity:
        result = run_analysis(path, VcfParser(), [ClinVarAnnotator(data)], identity=identity)
    assert result.annotations == []
    assert result.vcf_quality is not None
    assert result.vcf_quality.clinvar_ref_mismatches == mismatches
    assert result.vcf_quality.reference_source_conflicts == conflicts


def test_reference_source_conflict_is_reported_truthfully(data: Path, tmp_path: Path) -> None:
    _identity_rows(data, ["A", "A"])
    path = _write(tmp_path / "source-conflict.vcf", [_row()])
    for suffix in ("json", "html"):
        out = tmp_path / f"conflict.{suffix}"
        result = CliRunner().invoke(
            main,
            ["analyze", str(path), "--data-dir", str(data), "--no-update", "--output", str(out)],
        )
        assert result.exit_code == 0, result.output
        assert (
            "1 reference call(s) were not assessed because reference sources disagree"
            in " ".join(result.output.split())
        )
        assert "REF that is not ClinVar's" not in result.output
        if suffix == "json":
            payload = json.loads(out.read_text())
            assert payload["annotations"] == []
            assert payload["input"]["vcf_quality"]["reference_source_conflicts"] == 1
            assert payload["input"]["vcf_quality"]["clinvar_ref_mismatches"] == 0
        else:
            assert (
                "1 reference call(s) were not assessed because reference sources disagree"
                in out.read_text()
            )


def test_rendered_reference_report_is_a_supported_diff_baseline(
    data: Path, tmp_path: Path
) -> None:
    path = _write(tmp_path / "reference.vcf", [_row()])
    with ClinVarAnnotator(data) as clinvar:
        analysis = run_analysis(path, VcfParser(), [clinvar])
    baseline = tmp_path / "baseline.json"
    render_json(analysis, output_path=baseline)
    payload = json.loads(baseline.read_text())
    assert payload["schema_version"] == SCHEMA_VERSION == "10"
    assert load_previous_report(baseline) == payload


def test_cli_current_report_round_trips_through_diff(data: Path, tmp_path: Path) -> None:
    path = _write(tmp_path / "reference.vcf", [_row()])
    baseline = tmp_path / "baseline.json"
    output = tmp_path / "diff.json"
    runner = CliRunner()
    args = ["analyze", str(path), "--data-dir", str(data), "--no-update"]
    first = runner.invoke(main, [*args, "--output", str(baseline)])
    assert first.exit_code == 0, first.output
    second = runner.invoke(main, [*args, "--diff", str(baseline), "--output", str(output)])
    assert second.exit_code == 0, second.output
    payload = json.loads(output.read_text())
    assert payload["schema_version"] == SCHEMA_VERSION == "10"
    assert payload["diff"]["new"] == []
    assert payload["diff"]["removed"] == []
    assert payload["diff"]["changed"] == []


@pytest.mark.parametrize("version", [str(n) for n in range(1, 11)])
def test_all_previous_report_schemas_remain_supported(tmp_path: Path, version: str) -> None:
    baseline = tmp_path / "legacy.json"
    payload = {"schema_version": version, "annotations": []}
    baseline.write_text(json.dumps(payload))
    assert load_previous_report(baseline) == payload


@pytest.mark.parametrize("version", ["11", "12", "99"])
def test_future_report_schema_is_not_assumed_compatible(tmp_path: Path, version: str) -> None:
    baseline = tmp_path / "future.json"
    baseline.write_text(json.dumps({"schema_version": version, "annotations": []}))
    with pytest.raises(ValueError, match="schema version"):
        load_previous_report(baseline)


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("block", ["explicit", "anchor", "interior"])
@pytest.mark.parametrize("aliases", [False, True])
@pytest.mark.parametrize("repeats", [1, 120])
def test_contradictory_reference_ploidy_abstains_independent_of_order(
    data: Path, tmp_path: Path, reverse: bool, block: str, aliases: bool, repeats: int
) -> None:
    if aliases:
        _identity_rows(data, ["T", "T", "T"])
    rows = [
        (
            _row(gt="0")
            if block == "explicit"
            else _row(
                1990 if block == "interior" else 2000,
                "A" if block == "interior" else "T",
                "<NON_REF>",
                "0",
                "END=2010",
            )
        ),
        _row(gt="0/0"),
        _row(gt="0/0"),
        _row(ref="A", gt="0"),
    ]
    rows *= repeats
    path = _write(tmp_path / "ploidy.vcf", rows[::-1] if reverse else rows)
    with ClinVarAnnotator(data) as clinvar, RsidIdentityAnnotator(data) as identity:
        result = run_analysis(path, VcfParser(), [clinvar], identity=identity)
    assert result.annotations == []
    assert result.vcf_quality is not None
    assert result.vcf_quality.reference_ploidy_conflicts == 3 * repeats
    assert result.vcf_quality.clinvar_ref_mismatches == repeats
    assert result.vcf_quality.reference_source_conflicts == 0
    assert result.vcf_quality.reference_carried_calls == 0
    notes = " ".join(vcf_quality_lines(result.vcf_quality)).lower()
    assert f"{3 * repeats} reference call(s)" in notes and "ploidy" in notes
    output = tmp_path / "ploidy.json"
    render_json(result, output_path=output)
    payload = json.loads(output.read_text())
    assert payload["schema_version"] == SCHEMA_VERSION == "10"
    assert payload["input"]["vcf_quality"]["reference_ploidy_conflicts"] == 3 * repeats


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("gt", ["0", "0/0"])
def test_matching_reference_ploidy_ignores_ref_mismatches(
    data: Path, tmp_path: Path, reverse: bool, gt: str
) -> None:
    rows = [_row(gt=gt), _row(ref="A", gt="0/0" if gt == "0" else "0")]
    path = _write(tmp_path / "mismatch.vcf", rows[::-1] if reverse else rows)
    with ClinVarAnnotator(data) as clinvar:
        result = run_analysis(path, VcfParser(), [clinvar])
    assert [(a.genotype_match, a.zygosity_code) for a in result.annotations] == [
        ("T", "hemizygous") if gt == "0" else ("TT", "homozygous")
    ]
    assert result.vcf_quality is not None
    assert result.vcf_quality.clinvar_ref_mismatches == 1
    assert getattr(result.vcf_quality, "reference_ploidy_conflicts", -1) == 0


@pytest.mark.parametrize("gt", ["0", "0/0"])
def test_agreeing_duplicate_and_overlapping_reference_ploidy(
    data: Path, tmp_path: Path, gt: str
) -> None:
    rows = [_row(gt=gt)] * 120 + [_row(1990, "A", "<NON_REF>", gt, "END=2010")] * 120
    path = _write(tmp_path / "agreement.vcf", rows)
    with ClinVarAnnotator(data) as clinvar:
        result = run_analysis(path, VcfParser(), [clinvar])
    assert [(a.genotype_match, a.zygosity_code) for a in result.annotations] == [
        ("T", "hemizygous") if gt == "0" else ("TT", "homozygous")
    ]
    assert result.vcf_quality is not None
    assert result.vcf_quality.hom_ref_calls == 240
    assert result.vcf_quality.reference_carried_calls == 1
    assert getattr(result.vcf_quality, "reference_ploidy_conflicts", -1) == 0


def test_reference_provider_conflict_takes_precedence_over_ploidy(
    data: Path, tmp_path: Path
) -> None:
    _identity_rows(data, ["T", "A", "T"])
    rows = [_row(gt="0"), _row(gt="0/0"), _row(ref="A"), _row(gt="0")]
    path = _write(tmp_path / "providers.vcf", rows)
    with ClinVarAnnotator(data) as clinvar, RsidIdentityAnnotator(data) as identity:
        result = run_analysis(path, VcfParser(), [clinvar], identity=identity)
    assert result.annotations == []
    assert result.vcf_quality is not None
    assert result.vcf_quality.reference_source_conflicts == len(rows)
    assert result.vcf_quality.clinvar_ref_mismatches == 0
    assert getattr(result.vcf_quality, "reference_ploidy_conflicts", -1) == 0


@pytest.mark.parametrize("ref", ["A", "N", ".", "-", ""])
def test_repeated_explicit_ref_mismatches_keep_exact_observation_count(
    data: Path, tmp_path: Path, ref: str
) -> None:
    path = _write(tmp_path / "invalid-ref.vcf", [_row(ref=ref)] * 240)
    with ClinVarAnnotator(data) as clinvar:
        result = run_analysis(path, VcfParser(), [clinvar])
    assert result.annotations == []
    assert result.vcf_quality is not None
    assert result.vcf_quality.clinvar_ref_mismatches == 240
    assert result.vcf_quality.reference_source_conflicts == 0


def _retained_size(value: object, seen: set[int] | None = None) -> int:
    """Measure retained state without depending on its container representation."""
    if seen is None:
        seen = set()
    if id(value) in seen:
        return 0
    seen.add(id(value))
    size = sys.getsizeof(value)
    if isinstance(value, Mapping):
        size += sum(_retained_size(k, seen) + _retained_size(v, seen) for k, v in value.items())
    elif isinstance(value, (list, tuple, set, frozenset)):
        size += sum(_retained_size(item, seen) for item in value)
    elif hasattr(value, "__dict__"):
        size += _retained_size(vars(value), seen)
    else:
        for slot in getattr(type(value), "__slots__", ()):
            if hasattr(value, slot):
                size += _retained_size(getattr(value, slot), seen)
    return size


@pytest.mark.parametrize("overlap", [False, True])
def test_reference_observation_retention_is_bounded_per_watched_locus(
    tmp_path: Path, overlap: bool
) -> None:
    sizes = []
    for count in (200, 2000, 20000):
        rows = [_row(ref="A", gt="0"), _row(gt="0/0")]
        if overlap:
            rows.append(_row(1990, "A", "<NON_REF>", "0", "END=2010"))
        # Exercise the supported unsorted path, not just a same-start buffer.
        rows.append(_row(pos=2001))
        path = _write(tmp_path / "repeated.vcf", rows * count)
        parser = VcfParser()
        parser.watch({("2", 2000)}, set())
        assert list(parser.parse(path)) == []
        stats = parser.quality_stats
        assert stats.hom_ref_calls == count * len(rows)
        assert len(stats.locus_states) == 1
        assert ("2", 2000) in stats.reference_calls
        sizes.append(_retained_size(stats.reference_calls))
    assert max(sizes) < 8192, sizes
    assert max(sizes) - min(sizes) < 256, sizes


@pytest.mark.parametrize("vcf_out", [False, True])
def test_reference_record_provenance_is_kept_only_for_vcf_output(
    data: Path, tmp_path: Path, vcf_out: bool
) -> None:
    """Records differing only in their sample field keep reference state bounded.

    The records behind each reference call are kept for ``--vcf-out`` alone;
    without it the total reference state does not grow with the file.
    """
    sizes = []
    for count in (200, 2000, 20000):
        rows = []
        for n in range(count):
            # GQ varies, so every record is distinct; 2:2001 takes the unsorted path.
            rows += [_row(gt=f"0/0:{20 + n % 50_000}"), _row(pos=2001)]
        path = _write(tmp_path / "varying.vcf", rows)
        parser = VcfParser()
        with ClinVarAnnotator(data) as clinvar:
            result = run_analysis(path, parser, [clinvar], collect_vcf_write_inputs=vcf_out)
        assert [a.rsid for a in result.annotations] == ["rs900230001"]
        stats = parser.quality_stats
        assert len(stats.reference_rows) == (count if vcf_out else 0)
        sizes.append(_retained_size(stats.reference_calls) + _retained_size(stats.reference_rows))
    if not vcf_out:
        assert max(sizes) < 8192, sizes
        assert max(sizes) - min(sizes) < 256, sizes


@pytest.mark.parametrize("end", [1999, 2000, 2010])
@pytest.mark.parametrize("reverse", [False, True])
def test_watcher_retains_only_reference_blocks_covering_a_later_start(
    end: int, reverse: bool
) -> None:
    watcher = _LocusWatcher({"2": [2000, 2020]})
    block = _row(1990, "A", "<NON_REF>", "0", f"END={end}").split("\t")
    explicit = _row(gt="0/0").split("\t")
    rows = [block, explicit]
    covering = []
    for row in rows[::-1] if reverse else rows:
        covering.extend(watcher.feed(row))
    covering.extend(watcher.feed(_row(pos=2011).split("\t")))
    if not reverse:
        assert all(start == 2011 for _, start in watcher._held)
    covering.extend(watcher.finish())
    assert [(row, loci) for row, loci in covering if row == explicit] == [
        (explicit, [("2", 2000)])
    ]
    assert [(row, loci) for row, loci in covering if row == block] == (
        [(block, [("2", 2000)])] if end >= 2000 else []
    )

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Synthetic pipeline and decoded-output controls for annotated VCF refresh."""

from __future__ import annotations

import contextlib
import json
import sqlite3
from dataclasses import replace
from datetime import UTC, datetime
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner

from allelix.annotators.clinvar import ClinVarAnnotator, clinvar_db_filename
from allelix.cli import main
from allelix.parsers.vcf import VcfParser
from allelix.reports._pipeline import run_analysis
from allelix.reports.vcf import render_vcf
from tests.test_clinvar_matching_batch import _write_clinvar

if TYPE_CHECKING:
    from pathlib import Path

META = """##fileformat=VCFv4.2
##source=synthetic ALLELIX_ text remains
##contig=<ID=2,length=243199373>
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=FLAG,Number=0,Type=Flag,Description="Flag">
##INFO=<ID=OTHER,Number=1,Type=String,Description="Other">
##ALT=<ID=NON_REF,Description="Other allele">
##ALT=<ID=DEL,Description="Deletion">
##FILTER=<ID=q10,Description="Filtered">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
"""
OWNED = """##INFO=<ID=ALLELIX_CLINVAR,Number=A,Type=String,Description="Old">
##INFO=<ID=ALLELIX_RETIRED,Number=0,Type=Flag,Description="Old">
##ALLELIX_Version=old
##ALLELIX_Retired=old
"""
PLAIN = "DP=30;FLAG;OTHER=x%3By;DP=31"
STALE = "DP=30;ALLELIX_CLINVAR=OLD_PATHOGENIC;FLAG;ALLELIX_RETIRED;OTHER=x%3By;DP=31"
FIXED = datetime(2001, 1, 1, tzinfo=UTC)


def _row(
    pos: int = 1000,
    ref: str = "G",
    alt: str = "A",
    gt: str = "0/1",
    rsid: str = "rs900227001",
    info: str = PLAIN,
    filt: str = "PASS",
    gq: int = 99,
    dp: int = 20,
    second: str | None = None,
) -> str:
    sample = f"{gt}:{gq}:{dp}" + (f"\t{second}" if second else "")
    return f"2\t{pos}\t{rsid}\t{ref}\t{alt}\t50\t{filt}\t{info}\tGT:GQ:DP\t{sample}\n"


def _input(tmp_path: Path, rows: list[str], *, samples: str = "S", owned: bool = False) -> Path:
    path = tmp_path / "input.vcf"
    path.write_text(
        META
        + (OWNED if owned else "")
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + samples
        + "\n"
        + "".join(rows),
        encoding="utf-8",
    )
    return path


@pytest.fixture
def data(tmp_path: Path) -> Path:
    path = tmp_path / "data"
    path.mkdir()
    _write_clinvar(path)
    return path


def _analysis(
    src: Path,
    data: Path,
    *,
    sample: str | None = None,
    min_gq: int = 20,
    keep_filtered: bool = False,
    build: str = "GRCh37",
):
    return run_analysis(
        src,
        VcfParser(sample=sample, min_gq=min_gq, keep_filtered=keep_filtered),
        [ClinVarAnnotator(data)],
        build_override=build,
        collect_vcf_write_inputs=True,
    )


def _write(src: Path, out: Path, result) -> int:
    wi = result.vcf_write_inputs
    assert wi is not None
    return render_vcf(
        input_path=src,
        output_path=out,
        build=result.build,
        annotators_used=result.annotators_used,
        annotations_by_key=wi.annotations_by_key,
        recovered_rsids=wi.recovered_rsids,
        reference_rows=wi.reference_rows,
        run_date=FIXED,
    )


def _body(path: Path) -> list[list[str]]:
    return [
        line.split("\t")
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]


@pytest.mark.parametrize("kind", ["ordinary", "filtered", "sv", "block", "dot-block"])
def test_owned_cleanup_precedes_all_early_returns(tmp_path: Path, data: Path, kind: str) -> None:
    changes = {
        "ordinary": {},
        "filtered": {"filt": "q10"},
        "sv": {"alt": "<DEL>", "info": "END=2010;" + STALE},
        "block": {
            "pos": 1990,
            "ref": "A",
            "alt": "<NON_REF>",
            "gt": "0/0",
            "info": "END=2010;" + STALE,
        },
        "dot-block": {
            "pos": 1990,
            "ref": "A",
            "alt": ".",
            "gt": "0/0",
            "info": "END=2010;" + STALE,
        },
    }[kind]
    src = _input(
        tmp_path,
        [_row(**({"pos": 1100, "rsid": ".", "info": STALE} | changes)), _row(pos=5000, rsid=".")],
        owned=True,
    )
    before = src.read_bytes()
    result = _analysis(src, data)
    assert len(result.annotations) == int(kind in {"block", "dot-block"})
    out = tmp_path / "out.vcf"
    assert _write(src, out, result) == 0
    old, new = _body(src)[0], _body(out)[0]
    assert new[:7] == old[:7] and new[8:] == old[8:]
    assert new[7] == ("END=2010;" if kind in {"sv", "block", "dot-block"} else "") + PLAIN
    text = out.read_text(encoding="utf-8")
    assert "OLD_PATHOGENIC" not in text and "ALLELIX_RETIRED" not in text
    assert "##ALLELIX_Retired" not in text and text.count("##ALLELIX_Version=") == 1
    assert "##source=synthetic ALLELIX_ text remains\n" in text
    assert src.read_bytes() == before
    pysam = pytest.importorskip("pysam")
    with pysam.VariantFile(str(out)) as stream:
        records = list(stream)
        assert len(records) == 2
        assert not any(key.startswith("ALLELIX_") for record in records for key in record.info)


def test_stale_original_id_is_removed_without_restoring_id(tmp_path: Path) -> None:
    src = _input(tmp_path, [_row(info="ALLELIX_ORIGINAL_ID=older;ALLELIX_RETIRED")])
    out = tmp_path / "out.vcf"
    render_vcf(
        input_path=src, output_path=out, build="GRCh37", annotators_used=[], annotations_by_key={}
    )
    assert _body(out)[0][2] == "rs900227001"
    assert _body(out)[0][7] == "."


def test_three_analysis_passes_stabilize_after_recovery_history_cleanup(
    tmp_path: Path,
    data: Path,
) -> None:
    pysam = pytest.importorskip("pysam")
    src = _input(tmp_path, [_row(rsid=".")])
    before = src.read_bytes()
    outputs = [tmp_path / f"pass{number}.vcf" for number in range(1, 4)]
    for number, output in enumerate(outputs):
        result = _analysis(src if number == 0 else outputs[number - 1], data)
        assert len(result.annotations) == 1
        assert result.vcf_write_inputs is not None
        assert result.vcf_write_inputs.recovered_rsids == (
            {("2", 1000, "G", "A"): "rs900227001"} if number == 0 else {}
        )
        assert _write(src if number == 0 else outputs[number - 1], output, result) == 1
        record = _body(output)[0]
        assert record[2] == "rs900227001" and record[3:7] == _body(src)[0][3:7]
        assert record[8:] == _body(src)[0][8:]
        assert ("ALLELIX_ORIGINAL_ID=." in record[7]) == (number == 0)
        assert ("##INFO=<ID=ALLELIX_ORIGINAL_ID," in output.read_text(encoding="utf-8")) == (
            number == 0
        )
        with pysam.VariantFile(str(output)) as stream:
            decoded = next(stream)
            assert decoded.id == "rs900227001"
            assert decoded.info["ALLELIX_CLINVAR"] == ("clinvar_pathogenic|rs900227001|GENE",)
            assert ("ALLELIX_ORIGINAL_ID" in decoded.info) == (number == 0)
    assert outputs[0].read_bytes() != outputs[1].read_bytes()
    assert outputs[1].read_bytes() == outputs[2].read_bytes()
    assert src.read_bytes() == before


@pytest.mark.parametrize(
    "ref,alt,gt,pos,expected",
    [
        ("G", "A", "0/1", 1000, ("rs900227001",)),
        ("g", "a", "0/1", 1000, ("rs900227001",)),
        ("A", "G,T", "0/2", 3000, (None, "rs900244002")),
        ("a", "g,t", "0/2", 3000, (None, "rs900244002")),
        ("A", "G,T", "2/1", 3000, ("rs900244001", "rs900244002")),
        ("a", "g,t", "2/1", 3000, ("rs900244001", "rs900244002")),
        ("A", "G,<NON_REF>,T", "0/3", 3000, (None, None, "rs900244002")),
        ("a", "g,<NON_REF>,t", "0/3", 3000, (None, None, "rs900244002")),
        ("TAC", "TA", "1/1", 3999, ("rs900244003",)),
        ("tac", "ta", "1/1", 3999, ("rs900244003",)),
    ],
)
def test_pipeline_alt_slots_and_literal_alleles(
    tmp_path: Path, data: Path, ref, alt, gt, pos, expected
) -> None:
    pysam = pytest.importorskip("pysam")
    src = _input(tmp_path, [_row(pos=pos, ref=ref, alt=alt, gt=gt, rsid=".")])
    before = src.read_bytes()
    result = _analysis(src, data)
    assert len(result.annotations) == sum(value is not None for value in expected)
    if pos == 3999:
        assert result.annotations[0].record_key == ("2", pos, ref.upper(), alt.upper())
        assert result.annotations[0].locus == ("2", 4000)
    out = tmp_path / "out.vcf"
    assert _write(src, out, result) == 1
    old, new = _body(src)[0], _body(out)[0]
    assert new[:2] == old[:2] and new[3:7] == old[3:7] and new[8:] == old[8:]
    assert src.read_bytes() == before
    with pysam.VariantFile(str(out)) as stream:
        assert stream.header.info["ALLELIX_CLINVAR"].number == "A"
        record = next(stream)
        values = record.info["ALLELIX_CLINVAR"]
        assert len(values) == len(expected)
        assert (
            tuple(value.split("|")[1] if value not in (None, ".") else None for value in values)
            == expected
        )
        assert record.ref == ref and record.alts == tuple(alt.split(","))


@pytest.mark.parametrize(
    "gt,ref,alt,rsid",
    [
        ("0/0", "T", "C", "rs900230001"),
        ("0/0", "T", "C", "."),
        ("0|0", "T", "C", "."),
        ("0", "T", "C", "."),
        ("0/0", "t", "c", "."),
        ("0/0", "T", "C,G", "."),
        ("0/0", "T", "C,<NON_REF>,G", "."),
    ],
)
def test_pipeline_explicit_reference_decodes_in_each_concrete_slot(
    tmp_path: Path,
    data: Path,
    gt: str,
    ref: str,
    alt: str,
    rsid: str,
) -> None:
    pysam = pytest.importorskip("pysam")
    src = _input(tmp_path, [_row(pos=2000, ref=ref, alt=alt, gt=gt, rsid=rsid)])
    result = _analysis(src, data)
    assert len(result.annotations) == 1 and result.annotations[0].reference_carries_alt
    assert result.annotations[0].haploid == (gt == "0")
    out = tmp_path / "out.vcf"
    assert _write(src, out, result) == 1
    assert _body(src)[0][:7] == _body(out)[0][:7]
    assert _body(src)[0][8:] == _body(out)[0][8:]
    with pysam.VariantFile(str(out)) as stream:
        declaration = stream.header.info["ALLELIX_CLINVAR"]
        assert declaration.number == "A"
        assert declaration.description == (
            "ClinVar significance | rsID | gene. On a homozygous-reference or "
            "haploid reference record, the finding is for the allele the "
            "reference itself carries."
        )
        record = next(stream)
        assert record.samples["S"]["GT"] == ((0,) if gt == "0" else (0, 0))
        assert record.info["ALLELIX_CLINVAR"] == tuple(
            "." if allele == "<NON_REF>" else "clinvar_pathogenic|rs900230001|GENE"
            for allele in alt.split(",")
        )


@pytest.mark.parametrize("gt,gq,min_gq", [("0/0", 99, 20), ("0/0", 3, 0), ("0", 99, 20)])
def test_cli_passes_selected_sample_and_actual_quality_policy(
    tmp_path: Path,
    data: Path,
    gt: str,
    gq: int,
    min_gq: int,
) -> None:
    src = _input(
        tmp_path,
        [_row(pos=2000, ref="T", alt="C", gt="1/1", rsid=".", second=f"{gt}:{gq}:20")],
        samples="FIRST\tSECOND",
    )
    out = tmp_path / "out.vcf"
    result = CliRunner().invoke(
        main,
        [
            "analyze",
            str(src),
            "--no-update",
            "--data-dir",
            str(data),
            "--build",
            "grch37",
            "--sample",
            "SECOND",
            "--min-gq",
            str(min_gq),
            "--vcf-out",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    assert "ALLELIX_CLINVAR=clinvar_pathogenic|rs900230001|GENE" in _body(out)[0][7]
    assert _body(src)[0][8:] == _body(out)[0][8:]
    assert "1 stamped row(s)" in result.output


def test_reference_enrichment_is_not_an_alt_annotation(tmp_path: Path, data: Path) -> None:
    src = _input(tmp_path, [_row(pos=2000, ref="T", alt="C,G", gt="0/0", rsid=".")])
    result = _analysis(src, data)
    assert len(result.annotations) == 1
    annotation = result.annotations[0]
    assert annotation.reference_carries_alt
    annotation.allele_frequency = 0.25
    annotation.am_pathogenicity = 0.99
    annotation.am_class = "likely_pathogenic"
    annotation.cadd_phred = 30.0
    out = tmp_path / "out.vcf"
    assert _write(src, out, result) == 1
    info = _body(out)[0][7]
    assert "ALLELIX_CLINVAR=" in info
    assert "ALLELIX_GNOMAD_AF=" not in info
    assert "ALLELIX_ALPHAMISSENSE=" not in info
    assert "ALLELIX_CADD=" not in info
    assert (
        annotation.allele_frequency,
        annotation.am_pathogenicity,
        annotation.am_class,
        annotation.cadd_phred,
    ) == (0.25, 0.99, "likely_pathogenic", 30.0)


@pytest.mark.parametrize("gt,gq,min_gq", [("0/0", 99, 20), ("0/0", 3, 0), ("0", 99, 20)])
def test_selected_sample_output_decodes(
    tmp_path: Path,
    data: Path,
    gt: str,
    gq: int,
    min_gq: int,
) -> None:
    pysam = pytest.importorskip("pysam")
    test_cli_passes_selected_sample_and_actual_quality_policy(tmp_path, data, gt, gq, min_gq)
    with pysam.VariantFile(str(tmp_path / "out.vcf")) as stream:
        assert stream.header.info["ALLELIX_CLINVAR"].number == "A"
        record = next(stream)
        assert record.samples["FIRST"]["GT"] == (1, 1)
        assert record.samples["SECOND"]["GT"] == ((0,) if gt == "0" else (0, 0))
        assert record.info["ALLELIX_CLINVAR"] == ("clinvar_pathogenic|rs900230001|GENE",)


@pytest.mark.parametrize(
    "changes",
    [
        {"source": "gwas"},
        {
            "reference_carries_alt": False,
            "alt": "C",
            "genotype_match": "CT",
            "carried_alleles": ("C", "T"),
        },
        {"locus": ("2", 2001)},
        {"locus_ref": "A"},
    ],
)
def test_reference_payload_requires_its_exact_clinvar_identity(
    tmp_path: Path, data: Path, changes
) -> None:
    src = _input(tmp_path, [_row(pos=2000, ref="T", alt="C", gt="0/0", rsid=".")])
    result = _analysis(src, data)
    assert result.vcf_write_inputs is not None
    rows = result.vcf_write_inputs.reference_rows
    out = tmp_path / "out.vcf"
    for annotation, expected in (
        (result.annotations[0], 1),
        (replace(result.annotations[0], **changes), 0),
    ):
        assert (
            render_vcf(
                input_path=src,
                output_path=out,
                build="GRCh37",
                annotators_used=result.annotators_used,
                annotations_by_key={("2", 2000, "T", "T"): [annotation]},
                reference_rows=rows,
            )
            == expected
        )
        out.unlink()


@pytest.mark.parametrize("finding", [False, True])
def test_owned_declaration_replaced_and_lookalike_kept(
    tmp_path: Path,
    data: Path,
    finding: bool,
) -> None:
    pysam = pytest.importorskip("pysam")
    declaration = (
        '##INFO=<ID=ALLELIX_CLINVAR,Number=A,Type=String,Description="Retired description">'
    )
    lookalike = (
        '##INFO=<Number=1,Type=String,Description="Keep comma,ID=ALLELIX_CLINVAR '
        'and \\"quoted\\" text",ID=KEEP>'
    )
    src = _input(tmp_path, [_row(info="DP=30;ALLELIX_CLINVAR=OLD_PATHOGENIC;KEEP=value")])
    src.write_text(
        src.read_text(encoding="utf-8").replace(
            "#CHROM", declaration + "\n" + lookalike + "\n#CHROM"
        ),
        encoding="utf-8",
    )
    with pysam.VariantFile(str(src)) as stream:
        assert stream.header.info["ALLELIX_CLINVAR"].description == "Retired description"
        assert next(stream).info["ALLELIX_CLINVAR"] == ("OLD_PATHOGENIC",)
    result = _analysis(src, data)
    wi = result.vcf_write_inputs
    assert wi is not None
    first, second = tmp_path / "first.vcf", tmp_path / "second.vcf"
    for inp, output in [(src, first), (first, second)]:
        assert render_vcf(
            input_path=inp,
            output_path=output,
            build="GRCh37",
            annotators_used=result.annotators_used,
            annotations_by_key=wi.annotations_by_key if finding else {},
            run_date=FIXED,
        ) == int(finding)
    assert first.read_bytes() == second.read_bytes()
    text = first.read_text(encoding="utf-8")
    assert declaration not in text
    assert lookalike + "\n" in text
    with pysam.VariantFile(str(first)) as stream:
        assert "Keep comma,ID=ALLELIX_CLINVAR" in stream.header.info["KEEP"].description
        record = next(stream)
        assert record.info["KEEP"] == "value"
        if finding:
            assert record.info["ALLELIX_CLINVAR"] == ("clinvar_pathogenic|rs900227001|GENE",)
            assert "homozygous-reference" in stream.header.info["ALLELIX_CLINVAR"].description
        else:
            assert "ALLELIX_CLINVAR" not in stream.header.info
            assert "ALLELIX_CLINVAR" not in record.info


@pytest.mark.parametrize(
    "changes",
    [
        {"gq": 3},
        {"dp": 0},
        {"filt": "q10"},
        {"gt": "0"},
        {"gt": "./0"},
        {"gt": "0/0/0"},
        {"gt": "1/1"},
        {"alt": "<DEL>"},
        {"alt": "<NON_REF>", "info": "END=2010"},
    ],
)
def test_reference_payload_does_not_override_row_eligibility(
    tmp_path: Path, data: Path, changes
) -> None:
    src = _input(tmp_path, [_row(pos=2000, ref="T", alt="C", gt="0/0", rsid=".")])
    result = _analysis(src, data)
    src.write_text(
        src.read_text(encoding="utf-8").split("#CHROM")[0]
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        + _row(**({"pos": 2000, "ref": "T", "alt": "C", "gt": "0/0", "rsid": "."} | changes)),
        encoding="utf-8",
    )
    out = tmp_path / "out.vcf"
    assert _write(src, out, result) == 0
    assert _body(src) == _body(out)


def test_wrong_build_does_not_create_reference_payload(tmp_path: Path, data: Path) -> None:
    src = _input(tmp_path, [_row(pos=2000, ref="T", alt="C", gt="0/0", rsid=".")])
    src.write_text(
        src.read_text(encoding="utf-8").replace(",length=243199373", ""), encoding="utf-8"
    )
    result = _analysis(src, data, build="GRCh38")
    assert result.annotations == []
    out = tmp_path / "out.vcf"
    assert _write(src, out, result) == 0
    assert _body(src) == _body(out)


@pytest.mark.parametrize("finding", [True, False])
def test_pipeline_reannotation_replaces_stale_values_and_is_idempotent(
    tmp_path: Path,
    data: Path,
    finding: bool,
) -> None:
    pysam = pytest.importorskip("pysam")
    src = _input(tmp_path, [_row(pos=1000 if finding else 1100, info=STALE)], owned=True)
    first, second = tmp_path / "first.vcf", tmp_path / "second.vcf"
    assert _write(src, first, _analysis(src, data)) == int(finding)
    assert _write(first, second, _analysis(first, data)) == int(finding)
    assert first.read_bytes() == second.read_bytes()
    text = first.read_text(encoding="utf-8")
    assert text.count("##INFO=<ID=ALLELIX_CLINVAR,") == int(finding)
    assert "OLD_PATHOGENIC" not in text
    with pysam.VariantFile(str(second)) as stream:
        record = next(stream)
        if finding:
            assert record.info["ALLELIX_CLINVAR"] == ("clinvar_pathogenic|rs900227001|GENE",)
        else:
            assert "ALLELIX_CLINVAR" not in record.info


@pytest.mark.parametrize(
    "alt,pos", [("<NON_REF>", 1990), ("<*>", 1990), (".", 1990), ("<NON_REF>", 2000), (".", 2000)]
)
def test_reference_blocks_and_zero_alt_slots_remain_unstamped(
    tmp_path: Path,
    data: Path,
    alt: str,
    pos: int,
) -> None:
    info = PLAIN if alt == "." and pos == 2000 else "END=2010;" + PLAIN
    src = _input(
        tmp_path,
        [_row(pos=pos, ref="T" if pos == 2000 else "A", alt=alt, gt="0/0", rsid=".", info=info)],
    )
    result = _analysis(src, data)
    assert len(result.annotations) == 1 and result.annotations[0].reference_carries_alt
    out = tmp_path / "out.vcf"
    assert _write(src, out, result) == 0
    assert _body(src) == _body(out)


@pytest.mark.parametrize(
    "changes", [{"gq": 3}, {"dp": 0}, {"filt": "q10"}, {"ref": "A"}, {"alt": "<DEL>"}]
)
def test_pipeline_reference_exclusions_are_not_reconstructed(
    tmp_path: Path, data: Path, changes
) -> None:
    src = _input(
        tmp_path,
        [
            _row(**({"pos": 2000, "ref": "T", "alt": "C", "gt": "0/0", "rsid": "."} | changes)),
            _row(),
        ],
    )
    result = _analysis(src, data)
    assert not any(a.reference_carries_alt for a in result.annotations)
    out = tmp_path / "out.vcf"
    assert _write(src, out, result) == 1
    assert "ALLELIX_" not in _body(out)[0][7]


@pytest.mark.parametrize("second", ["0/0:3:20", "0:99:20"])
def test_conflicting_reference_observations_do_not_gain_stamps(
    tmp_path: Path, data: Path, second: str
) -> None:
    src = _input(
        tmp_path,
        [
            _row(pos=2000, ref="T", alt="C", gt="0/0", rsid="."),
            _row(
                pos=2000,
                ref="T",
                alt="C",
                gt=second.split(":")[0],
                gq=int(second.split(":")[1]),
                rsid=".",
            ),
        ],
    )
    result = _analysis(src, data)
    assert not result.annotations
    out = tmp_path / "out.vcf"
    assert _write(src, out, result) == 0
    assert _body(src) == _body(out)


def _cli(src: Path, data: Path, tmp_path: Path) -> tuple[dict, list[list[str]]]:
    """Run analyze with --output and --vcf-out; return the JSON and VCF records."""
    report, out = tmp_path / "report.json", tmp_path / "out.vcf"
    args = ["analyze", str(src), "--no-update", "--data-dir", str(data), "--build", "grch37"]
    result = CliRunner().invoke(main, [*args, "--output", str(report), "--vcf-out", str(out)])
    assert result.exit_code == 0, result.output
    return json.loads(report.read_text(encoding="utf-8")), _body(out)


def test_excluded_ucsc_mito_record_does_not_borrow_the_mt_stamp(
    tmp_path: Path, data: Path
) -> None:
    """MT and chrM records at one position: only the record the analysis read is stamped."""
    with contextlib.closing(sqlite3.connect(data / clinvar_db_filename("GRCh37"))) as conn:
        conn.execute(
            "INSERT INTO clinvar_variants (rsid, chromosome, position, ref, alt, "
            "clinical_significance, condition, gene, review_status, allele_id, aggregate, "
            "submitters_same, submitters_total) VALUES ('rs900250001', 'MT', 1000, 'G', 'G', "
            "'Pathogenic', 'Synthetic condition', 'GENE', "
            "'criteria_provided,_single_submitter', 6, 'Pathogenic', 1, 1)"
        )
        conn.commit()
    src = tmp_path / "input.vcf"
    src.write_text(
        META
        + "##contig=<ID=MT,length=16569>\n"
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        + "MT\t1000\t.\tG\tA\t50\tPASS\t.\tGT:GQ\t0/0:99\n"
        + "chrM\t1000\t.\tG\tA\t50\tPASS\t.\tGT:GQ\t0/0:99\n",
        encoding="utf-8",
    )
    report, records = _cli(src, data, tmp_path)
    assert [a["rsid"] for a in report["annotations"]] == ["rs900250001"]
    assert report["input"]["build_evidence"]["mt_records_skipped"] == 1
    assert [(r[0], r[7]) for r in records] == [
        ("MT", "ALLELIX_CLINVAR=clinvar_pathogenic|rs900250001|GENE"),
        ("chrM", "."),
    ]


@pytest.mark.parametrize("pos,ref,alt", [(2000, "TC", "T"), (1999, "AT", "A")])
def test_multi_base_ref_record_carries_its_reference_finding(
    tmp_path: Path, data: Path, pos: int, ref: str, alt: str
) -> None:
    """The finding at 2:2000 goes on a 0/0 record whose REF starts there or covers it."""
    src = _input(tmp_path, [_row(pos=pos, ref=ref, alt=alt, gt="0/0", rsid=".", info=".")])
    report, records = _cli(src, data, tmp_path)
    assert [(a["rsid"], a["genotype_match"]) for a in report["annotations"]] == [
        ("rs900230001", "TT")
    ]
    assert records[0][3:5] == [ref, alt]
    assert records[0][7] == "ALLELIX_CLINVAR=clinvar_pathogenic|rs900230001|GENE"

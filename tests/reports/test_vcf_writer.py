# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for the annotated VCF output writer (#150 PR 1).

Writer is invoked directly with hand-built annotation dicts — no
pipeline coupling in PR 1. Fixture files under ``tests/fixtures/``
provide the shapes we care about: single-sample, multi-sample,
multi-allelic (``mock_vcf.vcf`` line 8: ``A → G,C``), and gVCF
``<NON_REF>`` + reference blocks.
"""

from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path
from typing import TYPE_CHECKING

from allelix import __version__
from allelix.models import Annotation
from allelix.reports.vcf import VCF_LICENSE, VCF_SCHEMA_VERSION, render_vcf

if TYPE_CHECKING:
    import pytest

FIXTURES = Path(__file__).resolve().parent.parent / "fixtures"


def _ann(**overrides: object) -> Annotation:
    """Annotation factory with sane defaults for a ClinVar call."""
    defaults: dict[str, object] = {
        "source": "clinvar",
        "rsid": "rs1801133",
        "significance": "clinvar_pathogenic",
        "category": "clinical",
        "magnitude": 9.0,
        "description": "MTHFR C677T",
        "attribution": "ClinVar",
        "genotype_match": "AG",
        "gene": "MTHFR",
        "condition": "MTHFR deficiency",
        "alt": "A",
    }
    defaults.update(overrides)
    return Annotation(**defaults)  # type: ignore[arg-type]


def _fixed_date() -> datetime:
    return datetime(2026, 7, 2, 14, 30, 0, tzinfo=UTC)


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _split_header_and_body(text: str) -> tuple[list[str], list[str]]:
    header: list[str] = []
    body: list[str] = []
    for line in text.splitlines():
        if line.startswith("#"):
            header.append(line)
        else:
            body.append(line)
    return header, body


def _declared_contigs(text: str) -> list[str]:
    """Extract raw IDs from VCF ``##contig`` declarations."""
    prefix = "##contig=<ID="
    return [
        line[len(prefix) :].split(",", 1)[0].split(">", 1)[0]
        for line in text.splitlines()
        if line.startswith(prefix)
    ]


def _emitted_contigs(text: str) -> list[str]:
    """Return unique raw CHROM values in first-data-line order."""
    result: list[str] = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        contig = line.split("\t", 1)[0]
        if contig not in result:
            result.append(contig)
    return result


class TestHeaderPreservation:
    def test_mock_fixture_itself_contains_mixed_undeclared_contigs(self) -> None:
        """The undeclared-contig case originates in the fixture, not writer normalization."""
        text = _read(FIXTURES / "mock_vcf.vcf")
        assert _emitted_contigs(text) == ["1", "22", "19", "17", "chr1", "chrX", "chrM"]
        assert _declared_contigs(text) == ["1", "22"]

    def test_input_meta_lines_preserved_verbatim(self, tmp_path: Path) -> None:
        out = tmp_path / "out.vcf"
        render_vcf(
            input_path=FIXTURES / "mock_vcf.vcf",
            output_path=out,
            build="GRCh38",
            annotators_used=[("clinvar", "20260101")],
            annotations_by_key={},
            run_date=_fixed_date(),
        )
        text = _read(out)
        assert "##fileformat=VCFv4.2" in text
        assert '##FILTER=<ID=PASS,Description="All filters passed">' in text
        assert "##contig=<ID=1,length=249250621,assembly=GRCh38>" in text
        assert '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">' in text

    def test_missing_emitted_contigs_receive_minimal_declarations(self, tmp_path: Path) -> None:
        out = tmp_path / "out.vcf"
        render_vcf(
            input_path=FIXTURES / "mock_vcf.vcf",
            output_path=out,
            build="GRCh38",
            annotators_used=[("clinvar", "20260101")],
            annotations_by_key={},
            run_date=_fixed_date(),
        )
        text = _read(out)
        declared = _declared_contigs(text)
        emitted = _emitted_contigs(text)

        assert set(emitted) <= set(declared)
        assert declared == ["1", "22", "19", "17", "chr1", "chrX", "chrM"]
        assert declared.count("1") == 1
        assert declared.count("22") == 1
        for contig in ("17", "19", "chr1", "chrM", "chrX"):
            assert f"##contig=<ID={contig}>" in text

    def test_chrom_line_preserved_and_appears_after_injected_block(self, tmp_path: Path) -> None:
        out = tmp_path / "out.vcf"
        render_vcf(
            input_path=FIXTURES / "mock_vcf.vcf",
            output_path=out,
            build="GRCh38",
            annotators_used=[("clinvar", "20260101")],
            annotations_by_key={
                ("1", 11796321, "G", "A"): [_ann()],
            },
            run_date=_fixed_date(),
        )
        text = _read(out)
        chrom_idx = text.index("#CHROM\tPOS\tID\tREF\tALT")
        clinvar_info_idx = text.index("##INFO=<ID=ALLELIX_CLINVAR")
        version_idx = text.index("##ALLELIX_Version=")
        # Injected block appears before #CHROM.
        assert clinvar_info_idx < chrom_idx
        assert version_idx < chrom_idx


class TestProvenanceBlock:
    def _run(self, tmp_path: Path) -> str:
        out = tmp_path / "out.vcf"
        render_vcf(
            input_path=FIXTURES / "mock_vcf.vcf",
            output_path=out,
            build="GRCh38",
            annotators_used=[("clinvar", "20260101"), ("gnomad", "4.1")],
            annotations_by_key={},
            run_date=_fixed_date(),
        )
        return _read(out)

    def test_version_line(self, tmp_path: Path) -> None:
        text = self._run(tmp_path)
        assert f"##ALLELIX_Version={__version__}" in text

    def test_schema_version_line(self, tmp_path: Path) -> None:
        text = self._run(tmp_path)
        assert f"##ALLELIX_VCF_SchemaVersion={VCF_SCHEMA_VERSION}" in text

    def test_license_line(self, tmp_path: Path) -> None:
        text = self._run(tmp_path)
        assert f"##ALLELIX_License={VCF_LICENSE}" in text

    def test_run_date_line(self, tmp_path: Path) -> None:
        text = self._run(tmp_path)
        assert "##ALLELIX_RunDate=2026-07-02T14:30:00+00:00" in text

    def test_build_line(self, tmp_path: Path) -> None:
        text = self._run(tmp_path)
        assert "##ALLELIX_Build=GRCh38" in text

    def test_database_lines_per_annotator(self, tmp_path: Path) -> None:
        text = self._run(tmp_path)
        assert "##ALLELIX_Database=<Name=clinvar,Version=20260101>" in text
        assert "##ALLELIX_Database=<Name=gnomad,Version=4.1>" in text

    def test_attribution_lines_carry_spdx_and_url(self, tmp_path: Path) -> None:
        text = self._run(tmp_path)
        # gnomAD carries ODbL-1.0; ClinVar has a custom license identifier.
        # Full display names and URLs come from the annotator classes.
        assert "##ALLELIX_Attribution=<Name=gnomAD," in text
        assert "License=ODbL-1.0" in text
        assert "##ALLELIX_Attribution=<Name=ClinVar," in text

    def test_missing_annotator_version_renders_unknown(self, tmp_path: Path) -> None:
        out = tmp_path / "out.vcf"
        render_vcf(
            input_path=FIXTURES / "mock_vcf.vcf",
            output_path=out,
            build="GRCh37",
            annotators_used=[("clinvar", None)],
            annotations_by_key={},
            run_date=_fixed_date(),
        )
        text = _read(out)
        assert "##ALLELIX_Database=<Name=clinvar,Version=unknown>" in text


class TestInfoDeclarations:
    def test_declaration_only_emitted_for_used_sources(self, tmp_path: Path) -> None:
        out = tmp_path / "out.vcf"
        render_vcf(
            input_path=FIXTURES / "mock_vcf.vcf",
            output_path=out,
            build="GRCh38",
            annotators_used=[
                ("clinvar", "20260101"),
                ("pharmgkb", "2026-05-11"),
                ("gnomad", "4.1"),
            ],
            annotations_by_key={
                ("1", 11796321, "G", "A"): [_ann()],
            },
            run_date=_fixed_date(),
        )
        text = _read(out)
        # ClinVar appears in an annotation, so its INFO line is declared.
        assert "##INFO=<ID=ALLELIX_CLINVAR,Number=A,Type=String" in text
        # No pharmgkb / gnomad annotations were produced this run, so
        # their INFO declarations are omitted from the header.
        assert "##INFO=<ID=ALLELIX_CLINPGX" not in text
        assert "##INFO=<ID=ALLELIX_GNOMAD_AF" not in text

    def test_all_used_source_declarations_present(self, tmp_path: Path) -> None:
        out = tmp_path / "out.vcf"
        render_vcf(
            input_path=FIXTURES / "mock_vcf.vcf",
            output_path=out,
            build="GRCh38",
            annotators_used=[
                ("clinvar", "20260101"),
                ("gnomad", "4.1"),
                ("cadd", "1.7"),
            ],
            annotations_by_key={
                # The pipeline mutates enrichment fields onto the primary
                # (ClinVar) Annotation object — it does NOT create new
                # source="gnomad" / source="cadd" annotations. This test
                # matches that real shape so the ##INFO declaration gate
                # exercises the enrichment detection path.
                ("1", 11796321, "G", "A"): [
                    _ann(allele_frequency=0.001, cadd_phred=28.3),
                ],
            },
            run_date=_fixed_date(),
        )
        text = _read(out)
        assert 'Type=String,Description="ClinVar significance | rsID | gene"' in text
        assert 'Type=Float,Description="gnomAD allele frequency"' in text
        assert 'Type=Float,Description="CADD PHRED score"' in text

    def test_original_id_declaration_only_when_recovery_supplied(self, tmp_path: Path) -> None:
        out_without = tmp_path / "without.vcf"
        render_vcf(
            input_path=FIXTURES / "mock_vcf.vcf",
            output_path=out_without,
            build="GRCh38",
            annotators_used=[("clinvar", "20260101")],
            annotations_by_key={("1", 11796321, "G", "A"): [_ann()]},
            run_date=_fixed_date(),
        )
        assert "ALLELIX_ORIGINAL_ID" not in _read(out_without)

        out_with = tmp_path / "with.vcf"
        render_vcf(
            input_path=FIXTURES / "mock_vcf_rsidless.vcf",
            output_path=out_with,
            build="GRCh38",
            annotators_used=[("clinvar", "20260101")],
            annotations_by_key={},
            recovered_rsids={("1", 11856378, "G", "A"): "rs1801133"},
            run_date=_fixed_date(),
        )
        assert "##INFO=<ID=ALLELIX_ORIGINAL_ID,Number=1" in _read(out_with)


class TestDataLineStamping:
    def test_single_sample_annotation_stamped(self, tmp_path: Path) -> None:
        out = tmp_path / "out.vcf"
        count = render_vcf(
            input_path=FIXTURES / "mock_vcf.vcf",
            output_path=out,
            build="GRCh38",
            annotators_used=[("clinvar", "20260101")],
            annotations_by_key={
                ("1", 11796321, "G", "A"): [_ann()],
            },
            run_date=_fixed_date(),
        )
        assert count == 1
        text = _read(out)
        # Find the rs1801133 data line and confirm the ALLELIX_CLINVAR
        # field was appended to its INFO column.
        target = [
            line for line in text.splitlines() if "rs1801133" in line and "\t11796321\t" in line
        ]
        assert len(target) == 1
        assert "ALLELIX_CLINVAR=clinvar_pathogenic|rs1801133|MTHFR" in target[0]
        # Original INFO content (DP=30) preserved.
        assert "DP=30" in target[0]

    def test_variant_with_no_hits_passes_through_unchanged(self, tmp_path: Path) -> None:
        out = tmp_path / "out.vcf"
        render_vcf(
            input_path=FIXTURES / "mock_vcf.vcf",
            output_path=out,
            build="GRCh38",
            annotators_used=[("clinvar", "20260101")],
            annotations_by_key={},  # No hits anywhere.
            run_date=_fixed_date(),
        )
        text = _read(out)
        # rs4680 data line has no ALLELIX_ prefix in its INFO column.
        rs4680_line = next(line for line in text.splitlines() if "rs4680" in line)
        info_col = rs4680_line.split("\t")[7]
        assert "ALLELIX_" not in info_col

    def test_gnomad_af_stamped_from_enriched_primary_annotation(self, tmp_path: Path) -> None:
        """gnomAD AF is an enrichment field the pipeline mutates onto a
        primary Annotation (typically ClinVar). The writer must extract
        it via ``_extract_enrichment`` and emit ``ALLELIX_GNOMAD_AF``
        even though no ``source="gnomad"`` Annotation exists."""
        out = tmp_path / "out.vcf"
        render_vcf(
            input_path=FIXTURES / "mock_vcf.vcf",
            output_path=out,
            build="GRCh38",
            annotators_used=[("clinvar", "20260101"), ("gnomad", "4.1")],
            annotations_by_key={
                ("1", 11796321, "G", "A"): [_ann(allele_frequency=0.0012)],
            },
            run_date=_fixed_date(),
        )
        target = next(line for line in _read(out).splitlines() if "\t11796321\t" in line)
        assert "ALLELIX_GNOMAD_AF=0.0012" in target

    def test_all_enrichment_fields_stamped_together(self, tmp_path: Path) -> None:
        """One primary ClinVar annotation with all three enrichment
        fields populated emits ALLELIX_CLINVAR, ALLELIX_GNOMAD_AF,
        ALLELIX_ALPHAMISSENSE, and ALLELIX_CADD on the same row."""
        out = tmp_path / "out.vcf"
        render_vcf(
            input_path=FIXTURES / "mock_vcf.vcf",
            output_path=out,
            build="GRCh38",
            annotators_used=[
                ("clinvar", "20260101"),
                ("gnomad", "4.1"),
                ("alphamissense", "1.0"),
                ("cadd", "1.7"),
            ],
            annotations_by_key={
                ("1", 11796321, "G", "A"): [
                    _ann(
                        allele_frequency=0.0012,
                        am_class="likely_pathogenic",
                        am_pathogenicity=0.92,
                        cadd_phred=28.3,
                    ),
                ],
            },
            run_date=_fixed_date(),
        )
        info = next(line for line in _read(out).splitlines() if "\t11796321\t" in line).split(
            "\t"
        )[7]
        assert "ALLELIX_CLINVAR=clinvar_pathogenic|rs1801133|MTHFR" in info
        assert "ALLELIX_GNOMAD_AF=0.0012" in info
        assert "ALLELIX_ALPHAMISSENSE=likely_pathogenic|0.92" in info
        assert "ALLELIX_CADD=28.3" in info

    def test_stamped_count_reflects_stamped_rows(self, tmp_path: Path) -> None:
        out = tmp_path / "out.vcf"
        count = render_vcf(
            input_path=FIXTURES / "mock_vcf.vcf",
            output_path=out,
            build="GRCh38",
            annotators_used=[("clinvar", "20260101"), ("gnomad", "4.1")],
            annotations_by_key={
                ("1", 11796321, "G", "A"): [_ann()],
                ("22", 19963748, "G", "A"): [
                    # gnomAD AF as an enrichment on a primary ClinVar
                    # annotation (matches how the pipeline actually
                    # produces this shape).
                    _ann(allele_frequency=0.05, rsid="rs4680", alt="A"),
                ],
            },
            run_date=_fixed_date(),
        )
        assert count == 2

    def test_newline_in_annotation_value_does_not_corrupt_row(self, tmp_path: Path) -> None:
        """Annotation carrying an embedded ``\\n`` (adversarial input from
        a source with multi-line text) must be sanitized before it hits
        the file — otherwise the newline splits the VCF data line into
        two, producing an invalid record. Validity, not fidelity."""
        out = tmp_path / "out.vcf"
        render_vcf(
            input_path=FIXTURES / "mock_vcf.vcf",
            output_path=out,
            build="GRCh38",
            annotators_used=[("clinvar", "20260101")],
            annotations_by_key={
                ("1", 11796321, "G", "A"): [
                    _ann(gene="MTHFR\nInjectedLine\rMore"),
                ],
            },
            run_date=_fixed_date(),
        )
        text = _read(out)
        # There must be exactly one output data line for position 11796321.
        data_lines_for_pos = [
            line
            for line in text.splitlines()
            if line and not line.startswith("#") and "\t11796321\t" in line
        ]
        assert len(data_lines_for_pos) == 1
        # And no orphan garbage lines survived the split.
        assert "InjectedLine" in data_lines_for_pos[0]
        assert not any(line == "InjectedLine" for line in text.splitlines())
        # The stamped INFO value must have neutralized the newline / CR.
        info_col = data_lines_for_pos[0].split("\t")[7]
        assert "\n" not in info_col
        assert "\r" not in info_col
        assert "InjectedLine" in info_col

    def test_dot_info_column_is_replaced_not_prefixed(self, tmp_path: Path) -> None:
        out = tmp_path / "out.vcf"
        render_vcf(
            input_path=FIXTURES / "mock_multisample.vcf",
            output_path=out,
            build="GRCh38",
            annotators_used=[("clinvar", "20260101")],
            annotations_by_key={
                ("1", 11796321, "G", "A"): [_ann()],
            },
            run_date=_fixed_date(),
        )
        target = next(line for line in _read(out).splitlines() if "\t11796321\t" in line)
        info_col = target.split("\t")[7]
        # Input was `.` — output should be pure Allelix INFO, no leading dot.
        assert not info_col.startswith(".")
        assert info_col.startswith("ALLELIX_CLINVAR=")


class TestMultiAllelic:
    def test_per_allele_semantics_alt_order(self, tmp_path: Path) -> None:
        """mock_vcf line 8 is ``1 100000 rs900000001 A G,C ...`` — two ALTs.

        Annotate only ALT=C (the second slot). Expect
        ``ALLELIX_CLINVAR=.,<value>`` — empty first slot, value second.
        """
        out = tmp_path / "out.vcf"
        render_vcf(
            input_path=FIXTURES / "mock_vcf.vcf",
            output_path=out,
            build="GRCh38",
            annotators_used=[("clinvar", "20260101")],
            annotations_by_key={
                ("1", 100000, "A", "C"): [
                    _ann(rsid="rs900000001", alt="C"),
                ],
            },
            run_date=_fixed_date(),
        )
        target = next(line for line in _read(out).splitlines() if "\t100000\t" in line)
        info_col = target.split("\t")[7]
        # ALT order is G,C → first slot empty (.), second slot the value.
        assert "ALLELIX_CLINVAR=.,clinvar_pathogenic|rs900000001|MTHFR" in info_col

    def test_both_alts_annotated(self, tmp_path: Path) -> None:
        """Multi-allelic ``A → G,C`` with gnomAD AF enrichment on each
        ALT's primary Annotation produces per-allele
        ``ALLELIX_GNOMAD_AF=0.001,0.0005`` in ALT order."""
        out = tmp_path / "out.vcf"
        render_vcf(
            input_path=FIXTURES / "mock_vcf.vcf",
            output_path=out,
            build="GRCh38",
            annotators_used=[("clinvar", "20260101"), ("gnomad", "4.1")],
            annotations_by_key={
                ("1", 100000, "A", "G"): [
                    _ann(allele_frequency=0.001, alt="G"),
                ],
                ("1", 100000, "A", "C"): [
                    _ann(allele_frequency=0.0005, alt="C"),
                ],
            },
            run_date=_fixed_date(),
        )
        target = next(line for line in _read(out).splitlines() if "\t100000\t" in line)
        info_col = target.split("\t")[7]
        assert "ALLELIX_GNOMAD_AF=0.001,0.0005" in info_col


class TestChrPrefixNormalization:
    """Every real GRCh38 VCF from modern callers (GATK / DeepVariant /
    DRAGEN) uses ``chr``-prefixed contigs (``chr22`` etc.). The parser
    normalizes to bare form via ``normalize_chromosome`` at ingest, so
    the pipeline's ``annotations_by_key`` is keyed on ``"22"``, not
    ``"chr22"``. The writer reads the raw input line's column 0 and
    would look up the un-normalized ``"chr22"`` — a silent miss on the
    entire flagship use case if the normalization boundary drifts.

    This test locks in the invariant against a chr-prefixed fixture:
    the writer must find the annotation via a normalized lookup and
    stamp the row.
    """

    def test_chr_prefixed_input_still_stamps_annotations(self, tmp_path: Path) -> None:
        out = tmp_path / "out.vcf"
        # The chr-prefix fixture carries rs1801133-shaped data at chr1
        # 11796321 (matches our _ann() default). Key uses BARE form
        # ("1") because that's what the pipeline builds after parser
        # normalization.
        render_vcf(
            input_path=FIXTURES / "mock_vcf_chr_prefix_grch38.vcf",
            output_path=out,
            build="GRCh38",
            annotators_used=[("clinvar", "20260101")],
            annotations_by_key={
                ("1", 11796321, "G", "A"): [_ann()],
            },
            run_date=_fixed_date(),
        )
        text = _read(out)
        # The chr1 data line must appear stamped, with the output CHROM
        # preserved as ``chr1`` (writer keys via bare, writes verbatim).
        stamped_line = next(
            line for line in text.splitlines() if line.startswith("chr1\t11796321\t")
        )
        assert "ALLELIX_CLINVAR=clinvar_pathogenic|rs1801133|MTHFR" in stamped_line

    def test_chrm_normalizes_to_mt_for_lookup(self, tmp_path: Path) -> None:
        """Mitochondrial nomenclature drift: ``chrM`` (GRCh38 convention)
        must key to ``MT`` (Allelix canonical) or the annotation misses."""
        input_vcf = tmp_path / "input_with_chrM.vcf"
        input_vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chrM,length=16569>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chrM\t3243\trs199474657\tA\tG\t100\tPASS\t.\n",
            encoding="utf-8",
        )
        out = tmp_path / "out.vcf"
        render_vcf(
            input_path=input_vcf,
            output_path=out,
            build="GRCh38",
            annotators_used=[("clinvar", "20260101")],
            annotations_by_key={
                ("MT", 3243, "A", "G"): [
                    _ann(rsid="rs199474657", alt="G", gene="MT-TL1"),
                ],
            },
            run_date=_fixed_date(),
        )
        text = _read(out)
        stamped = next(line for line in text.splitlines() if line.startswith("chrM\t"))
        # Output CHROM preserved as chrM; lookup succeeded via MT
        # normalization.
        assert "ALLELIX_CLINVAR=clinvar_pathogenic|rs199474657|MT-TL1" in stamped


class TestBomAndEncoding:
    def test_bom_input_produces_bom_free_output(self, tmp_path: Path) -> None:
        """GH #126/#133: real Sequencing.com and DRAGEN VCFs ship with a
        UTF-8 BOM (``\\xef\\xbb\\xbf``). The parser strips it via
        ``utf-8-sig``; the writer must do the same so the output starts
        with a clean ``##fileformat`` and downstream tools that check
        for the exact header don't reject it."""
        out = tmp_path / "annotated.vcf"
        render_vcf(
            input_path=FIXTURES / "mock_vcf_bom.vcf",
            output_path=out,
            build="GRCh38",
            annotators_used=[("clinvar", "20260101")],
            annotations_by_key={},
            run_date=_fixed_date(),
        )
        raw_bytes = out.read_bytes()
        # BOM must NOT appear anywhere in the output — not at line 1, not
        # embedded in a header line, not smuggled anywhere.
        assert b"\xef\xbb\xbf" not in raw_bytes
        # Output starts with the standard fileformat directive, no BOM
        # or other prefix.
        assert raw_bytes.startswith(b"##fileformat=VCFv4.2")

    def test_non_utf8_byte_in_header_does_not_raise(self, tmp_path: Path) -> None:
        """GH #121: real VCFs commonly carry non-UTF-8 bytes in tool
        headers (a latin-1 ``\\xa9`` copyright in ``##CL=`` from GATK /
        DeepVariant / bcftools). The parser uses ``errors="replace"``
        to degrade gracefully; the writer must too, otherwise the
        run dies mid-stamp with UnicodeDecodeError."""
        input_vcf = tmp_path / "with_bad_byte.vcf"
        input_vcf.write_bytes(
            b"##fileformat=VCFv4.2\n"
            b"##CL=gatk-HC \xa9 2024\n"
            b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            b"1\t100\trs1\tA\tG\t.\t.\t.\n"
        )
        out = tmp_path / "annotated.vcf"
        # Must not raise. The bad byte should render as U+FFFD in the
        # output header line, matching parser behavior.
        render_vcf(
            input_path=input_vcf,
            output_path=out,
            build="GRCh38",
            annotators_used=[("clinvar", "20260101")],
            annotations_by_key={},
            run_date=_fixed_date(),
        )
        text = _read(out)
        assert "##fileformat=VCFv4.2" in text
        assert "##CL=" in text


class TestGzipInput:
    def test_gzipped_input_produces_same_output_as_plain(self, tmp_path: Path) -> None:
        """Gzip-compress the plain fixture on the fly and confirm the
        writer produces byte-identical output from both forms."""
        import gzip
        import shutil

        gz_input = tmp_path / "mock_vcf.vcf.gz"
        with (
            (FIXTURES / "mock_vcf.vcf").open("rb") as src,
            gzip.open(gz_input, "wb") as dst,
        ):
            shutil.copyfileobj(src, dst)

        out_plain = tmp_path / "plain.vcf"
        out_gz_in = tmp_path / "from_gz.vcf"
        annotations_by_key = {("1", 11796321, "G", "A"): [_ann()]}
        render_vcf(
            input_path=FIXTURES / "mock_vcf.vcf",
            output_path=out_plain,
            build="GRCh38",
            annotators_used=[("clinvar", "20260101")],
            annotations_by_key=annotations_by_key,
            run_date=_fixed_date(),
        )
        render_vcf(
            input_path=gz_input,
            output_path=out_gz_in,
            build="GRCh38",
            annotators_used=[("clinvar", "20260101")],
            annotations_by_key=annotations_by_key,
            run_date=_fixed_date(),
        )
        assert _read(out_plain) == _read(out_gz_in)


class TestSymbolicAlleles:
    def test_non_ref_symbolic_alt_is_not_annotated(self, tmp_path: Path) -> None:
        """gVCF ``<NON_REF>`` slots skip annotation lookup entirely."""
        out = tmp_path / "out.vcf"
        render_vcf(
            input_path=FIXTURES / "mock_gvcf.g.vcf",
            output_path=out,
            build="GRCh38",
            annotators_used=[("clinvar", "20260101")],
            annotations_by_key={
                # Legitimate hit on the first ALT of a mixed row.
                ("1", 11796321, "G", "A"): [_ann()],
            },
            run_date=_fixed_date(),
        )
        text = _read(out)
        # The mixed row `A,<NON_REF>` should stamp the A slot but leave
        # the symbolic slot as `.` — two alts, so Number=A produces two
        # comma-joined values.
        target = next(line for line in text.splitlines() if "\t11796321\t" in line)
        info_col = target.split("\t")[7]
        assert "ALLELIX_CLINVAR=clinvar_pathogenic|rs1801133|MTHFR,." in info_col

    def test_reference_block_line_passes_through(self, tmp_path: Path) -> None:
        """gVCF reference blocks (``.  A  <NON_REF>  ...  END=...``) carry
        no ALT annotation — they should pass through untouched."""
        out = tmp_path / "out.vcf"
        render_vcf(
            input_path=FIXTURES / "mock_gvcf.g.vcf",
            output_path=out,
            build="GRCh38",
            annotators_used=[("clinvar", "20260101")],
            annotations_by_key={},
            run_date=_fixed_date(),
        )
        text = _read(out)
        ref_block = next(line for line in text.splitlines() if "\t10001\t" in line)
        info_col = ref_block.split("\t")[7]
        assert info_col == "END=11796320"


class TestMultiSample:
    def test_all_sample_columns_preserved_verbatim(self, tmp_path: Path) -> None:
        out = tmp_path / "out.vcf"
        render_vcf(
            input_path=FIXTURES / "mock_multisample.vcf",
            output_path=out,
            build="GRCh38",
            annotators_used=[("clinvar", "20260101")],
            annotations_by_key={
                ("1", 11796321, "G", "A"): [_ann()],
            },
            run_date=_fixed_date(),
        )
        target = next(line for line in _read(out).splitlines() if "\t11796321\t" in line)
        # mock_multisample columns: CHROM POS ID REF ALT QUAL FILTER INFO
        # FORMAT SAMPLE_A SAMPLE_B SAMPLE_C
        cols = target.split("\t")
        assert cols[8] == "GT"
        assert cols[9] == "0/1"
        assert cols[10] == "1/1"
        assert cols[11] == "0/0"


class TestRecoveredRsid:
    def test_id_column_receives_recovered_rsid_and_original_stashed(self, tmp_path: Path) -> None:
        out = tmp_path / "out.vcf"
        render_vcf(
            input_path=FIXTURES / "mock_vcf_rsidless.vcf",
            output_path=out,
            build="GRCh38",
            annotators_used=[("clinvar", "20260101")],
            annotations_by_key={
                ("1", 11856378, "G", "A"): [_ann()],
            },
            recovered_rsids={
                ("1", 11856378, "G", "A"): "rs1801133",
            },
            run_date=_fixed_date(),
        )
        target = next(line for line in _read(out).splitlines() if "\t11856378\t" in line)
        cols = target.split("\t")
        # ID column stamped with recovered rsID.
        assert cols[2] == "rs1801133"
        info_col = cols[7]
        assert "ALLELIX_ORIGINAL_ID=." in info_col
        assert "ALLELIX_CLINVAR=" in info_col

    def test_stamped_count_includes_recovery_only_rows(self, tmp_path: Path) -> None:
        """A row with only a recovered rsID (no annotation) still counts."""
        out = tmp_path / "out.vcf"
        count = render_vcf(
            input_path=FIXTURES / "mock_vcf_rsidless.vcf",
            output_path=out,
            build="GRCh38",
            annotators_used=[],
            annotations_by_key={},
            recovered_rsids={
                ("1", 11856378, "G", "A"): "rs1801133",
            },
            run_date=_fixed_date(),
        )
        assert count == 1


class TestFormattingHelpers:
    def test_primary_clinvar(self) -> None:
        from allelix.reports.vcf import _format_primary

        result = _format_primary(_ann())
        assert result == ("clinvar", "clinvar_pathogenic|rs1801133|MTHFR")

    def test_primary_gwas_pipe_delimited(self) -> None:
        from allelix.reports.vcf import _format_primary

        result = _format_primary(
            _ann(
                source="gwas",
                trait="LDL cholesterol",
                p_value=1.2e-9,
                attribution="Kessler et al 2021",
            )
        )
        # Sanitization replaces the space in trait/attribution with underscores.
        assert result == ("gwas", "LDL_cholesterol|1.2e-09|Kessler_et_al_2021")

    def test_primary_pharmgkb(self) -> None:
        from allelix.reports.vcf import _format_primary

        result = _format_primary(
            _ann(source="pharmgkb", gene="CYP2D6", significance="Poor_Metabolizer")
        )
        assert result == ("pharmgkb", "CYP2D6|Poor_Metabolizer")

    def test_primary_snpedia(self) -> None:
        from allelix.reports.vcf import _format_primary

        result = _format_primary(
            _ann(source="snpedia", rsid="rs1234", significance="risk_variant")
        )
        assert result == ("snpedia", "rs1234|risk_variant")

    def test_primary_enrichment_sources_return_none(self) -> None:
        """Enrichment sources are never source-dispatched — they go
        through ``_extract_enrichment`` instead."""
        from allelix.reports.vcf import _format_primary

        assert _format_primary(_ann(source="gnomad", allele_frequency=0.01)) is None
        assert _format_primary(_ann(source="cadd", cadd_phred=28.3)) is None
        assert _format_primary(_ann(source="alphamissense", am_pathogenicity=0.9)) is None

    def test_primary_unknown_source_returns_none(self) -> None:
        from allelix.reports.vcf import _format_primary

        assert _format_primary(_ann(source="future_source")) is None

    def test_extract_enrichment_from_primary_annotation(self) -> None:
        """The pipeline mutates gnomAD AF / CADD / AlphaMissense fields
        onto existing PRIMARY Annotation objects — not new ones. This
        test mirrors that real shape: a ClinVar annotation with
        enrichment fields populated. All three must extract."""
        from allelix.reports.vcf import _extract_enrichment

        ann = _ann(
            source="clinvar",
            allele_frequency=0.0012,
            am_class="likely_pathogenic",
            am_pathogenicity=0.92,
            cadd_phred=28.3,
        )
        result = _extract_enrichment([ann])
        assert result == {
            "gnomad": "0.0012",
            "alphamissense": "likely_pathogenic|0.92",
            "cadd": "28.3",
        }

    def test_extract_enrichment_takes_first_non_none(self) -> None:
        """Multiple annotations at the same variant share the same
        enrichment; ``_extract_enrichment`` picks the first non-None
        value it finds."""
        from allelix.reports.vcf import _extract_enrichment

        anns = [
            _ann(source="clinvar", allele_frequency=None, cadd_phred=None),
            _ann(source="pharmgkb", allele_frequency=0.05, cadd_phred=15.5),
            _ann(source="gwas", allele_frequency=0.99, cadd_phred=99.9),
        ]
        result = _extract_enrichment(anns)
        assert result["gnomad"] == "0.05"
        assert result["cadd"] == "15.5"

    def test_extract_enrichment_empty_when_no_fields_populated(self) -> None:
        from allelix.reports.vcf import _extract_enrichment

        assert _extract_enrichment([_ann()]) == {}
        assert _extract_enrichment([]) == {}

    def test_sanitize_replaces_reserved_chars(self) -> None:
        from allelix.reports.vcf import _sanitize_info_value

        assert _sanitize_info_value("a;b,c=d e\tf") == "a_b_c_d_e_f"

    def test_sanitize_neutralizes_record_structural_chars(self) -> None:
        """``\\n`` and ``\\r`` are record-structural — an embedded newline
        in an INFO value would split the data line into two and corrupt
        the file. Both must be substituted."""
        from allelix.reports.vcf import _sanitize_info_value

        assert _sanitize_info_value("line1\nline2") == "line1_line2"
        assert _sanitize_info_value("cr\rlf\r\nend") == "cr_lf__end"

    def test_build_info_field_drops_all_empty_source(self) -> None:
        from allelix.reports.vcf import _build_info_field

        # ClinVar has values, gnomad is all-empty → gnomad dropped.
        result = _build_info_field(
            {
                "clinvar": ["foo", ""],
                "gnomad": ["", ""],
            }
        )
        assert result == "ALLELIX_CLINVAR=foo,."


class TestMalformedInput:
    def test_short_data_line_passes_through(self, tmp_path: Path) -> None:
        """A row with fewer than 8 columns is malformed; the writer
        passes it through rather than crashing the whole file."""
        input_vcf = tmp_path / "in.vcf"
        input_vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t100\trs1\tA\tG\t.\t.\t.\n"
            "1\t200\tshort_row\n",  # only 3 columns
            encoding="utf-8",
        )
        out = tmp_path / "out.vcf"
        count = render_vcf(
            input_path=input_vcf,
            output_path=out,
            build="GRCh38",
            annotators_used=[],
            annotations_by_key={},
            run_date=_fixed_date(),
        )
        assert count == 0
        text = _read(out)
        assert "1\t200\tshort_row\n" in text

    def test_non_integer_pos_passes_through(self, tmp_path: Path) -> None:
        input_vcf = tmp_path / "in.vcf"
        input_vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\tNOT_A_NUMBER\trs1\tA\tG\t.\t.\t.\n",
            encoding="utf-8",
        )
        out = tmp_path / "out.vcf"
        count = render_vcf(
            input_path=input_vcf,
            output_path=out,
            build="GRCh38",
            annotators_used=[],
            annotations_by_key={},
            run_date=_fixed_date(),
        )
        assert count == 0
        assert "NOT_A_NUMBER" in _read(out)


class TestUnknownAnnotator:
    def test_unknown_annotator_name_logs_and_skips_attribution(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        out = tmp_path / "out.vcf"
        with caplog.at_level("WARNING"):
            render_vcf(
                input_path=FIXTURES / "mock_vcf.vcf",
                output_path=out,
                build="GRCh38",
                annotators_used=[("does_not_exist", "1.0"), ("clinvar", "20260101")],
                annotations_by_key={},
                run_date=_fixed_date(),
            )
        assert any("does_not_exist" in r.message for r in caplog.records)
        text = _read(out)
        # Unknown annotator gets a Database line (best-effort) but no
        # Attribution line (no LicenseDescriptor to source from).
        assert "##ALLELIX_Database=<Name=does_not_exist,Version=1.0>" in text
        assert "##ALLELIX_Attribution=<Name=does_not_exist" not in text

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for the ClinVar annotator."""

from __future__ import annotations

import contextlib
import sqlite3
from typing import TYPE_CHECKING

import pytest

from allelix.annotators.clinvar import ClinVarAnnotator, clinvar_db_filename, clinvar_record_name
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
from allelix.databases.manager import load_clinvar_vcf
from allelix.models import Variant

if TYPE_CHECKING:
    from pathlib import Path


@pytest.fixture
def annotator(clinvar_data_dir: Path):
    """Yield an annotator and ensure its SQLite connection is closed.

    N-1: without explicit teardown, _connection() opens a sqlite3.Connection
    that's only reaped by GC. Yield + close() pins the contract every test
    relies on.
    """
    ann = ClinVarAnnotator(clinvar_data_dir)
    try:
        yield ann
    finally:
        ann.close()


@pytest.fixture
def annotator_with_benign(clinvar_data_dir: Path):
    """Annotator that includes benign/likely_benign annotations."""
    ann = ClinVarAnnotator(clinvar_data_dir, include_benign=True)
    try:
        yield ann
    finally:
        ann.close()


class TestSetupAndStatus:
    def test_unconfigured_is_not_ready(self, tmp_path: Path):
        ann = ClinVarAnnotator(tmp_path)
        assert ann.is_ready() is False
        assert ann.version() is None

    def test_configured_is_ready(self, annotator: ClinVarAnnotator):
        assert annotator.is_ready() is True
        assert annotator.version() is not None


class TestSignalGuard:
    def test_setup_aborts_when_signal_fetch_fails(self, tmp_path: Path, monkeypatch):
        """setup() raises RuntimeError when remote signal is None."""
        ann = ClinVarAnnotator(tmp_path)
        monkeypatch.setattr(
            ClinVarAnnotator, "_fetch_remote_signal_for", staticmethod(lambda _build: None)
        )
        with pytest.raises(RuntimeError, match="cannot verify remote freshness signal"):
            ann.setup()


class TestInterpreterVersionStamp:
    """CLINVAR_INTERPRETER_VERSION stamp in cache's local_version_tag."""

    def test_is_ready_accepts_matching_iv_stamp(self, annotator: ClinVarAnnotator):
        """Freshly loaded cache has the current iv stamp — is_ready returns True."""
        assert annotator.is_ready() is True

    def test_is_ready_rejects_cache_without_tag(
        self, tmp_path: Path, mock_clinvar_grch37_vcf: Path
    ):
        """Cache with no local_version_tag is self-healed by one-shot migration."""
        build = "GRCh37"
        db_path = tmp_path / clinvar_db_filename(build)
        load_clinvar_vcf(
            mock_clinvar_grch37_vcf,
            db_path,
            source_url="test://mock",
            record_name=clinvar_record_name(build),
        )
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            conn.execute(
                "UPDATE database_versions SET remote_signal = 'md5:abc', "
                "local_version_tag = NULL WHERE name = ?",
                (clinvar_record_name(build),),
            )
            conn.commit()
        ann = ClinVarAnnotator(tmp_path, builds=(build,))
        assert ann.is_ready() is True
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            tag = conn.execute(
                "SELECT local_version_tag FROM database_versions WHERE name = ?",
                (clinvar_record_name(build),),
            ).fetchone()[0]
            assert tag == f"iv:{CLINVAR_INTERPRETER_VERSION}"

    def test_is_ready_rejects_old_iv_stamp(self, tmp_path: Path, mock_clinvar_grch37_vcf: Path):
        """Cache stamped with an older iv version is rejected."""
        build = "GRCh37"
        db_path = tmp_path / clinvar_db_filename(build)
        load_clinvar_vcf(
            mock_clinvar_grch37_vcf,
            db_path,
            source_url="test://mock",
            record_name=clinvar_record_name(build),
        )
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            conn.execute(
                "UPDATE database_versions SET local_version_tag = 'iv:0' WHERE name = ?",
                (clinvar_record_name(build),),
            )
            conn.commit()
        ann = ClinVarAnnotator(tmp_path, builds=(build,))
        assert ann.is_ready() is False


class TestGenotypeMatching:
    """ADR-0007: ClinVar entries trigger only when the user carries ALT."""

    def test_heterozygous_carrier_triggers(self, annotator: ClinVarAnnotator):
        # mock ClinVar: rs1801133, REF=G, ALT=A, Pathogenic
        v = Variant("rs1801133", "1", 11796321, "G", "A")
        results = annotator.annotate(v)
        assert len(results) == 1
        a = results[0]
        assert a.significance == "clinvar_pathogenic"
        assert a.attribution == "ClinVar"
        assert a.source == "clinvar"
        assert a.category == "clinical"
        assert a.gene == "MTHFR"
        # ADR-0023: genotype_match shows the user's diploid (sorted), not
        # the matched ALT base. G/A → "AG".
        assert a.genotype_match == "AG"
        assert a.magnitude == 9.0

    def test_homozygous_alt_triggers(self, annotator: ClinVarAnnotator):
        # mock ClinVar: rs4680, REF=G, ALT=A, Drug_response
        v = Variant("rs4680", "22", 19963748, "A", "A")
        results = annotator.annotate(v)
        assert len(results) == 1
        assert results[0].significance == "clinvar_drug_response"
        assert results[0].magnitude == 6.5

    def test_homozygous_reference_does_not_trigger(self, annotator: ClinVarAnnotator):
        # mock ClinVar: rs121918506, REF=G, ALT=T, Pathogenic; mock has G/G
        v = Variant("rs121918506", "17", 7577538, "G", "G")
        assert annotator.annotate(v) == []

    def test_no_call_does_not_trigger(self, annotator: ClinVarAnnotator):
        v = Variant("rs1801133", "1", 11796321, "-", "-")
        assert annotator.annotate(v) == []

    def test_asymmetric_no_call_does_not_trigger(self, annotator: ClinVarAnnotator):
        """r-2: one good allele + one no-call must short-circuit before lookup.

        Catches mutations like `if variant.is_no_call` → `if variant.allele1 == "-"`
        that pass the both-no-call test but leak through here.
        """
        v_left = Variant("rs1801133", "1", 11796321, "-", "A")
        v_right = Variant("rs1801133", "1", 11796321, "A", "-")
        assert annotator.annotate(v_left) == []
        assert annotator.annotate(v_right) == []

    def test_unknown_rsid_does_not_trigger(self, annotator: ClinVarAnnotator):
        v = Variant("rs999000111", "1", 1000, "A", "T")
        assert annotator.annotate(v) == []


class TestAttribution:
    """ADR-0003: Significance and attribution must be source-prefixed."""

    def test_all_annotations_attribute_to_clinvar(self, annotator: ClinVarAnnotator):
        v = Variant("rs1801133", "1", 11796321, "G", "A")
        results = annotator.annotate(v)
        for a in results:
            assert a.attribution == "ClinVar"
            assert a.significance.startswith("clinvar_")
            assert a.category == "clinical"

    def test_description_attributes_to_clinvar(self, annotator: ClinVarAnnotator):
        v = Variant("rs1801133", "1", 11796321, "G", "A")
        results = annotator.annotate(v)
        assert results[0].description.startswith("ClinVar classifies")

    def test_review_status_populated(self, annotator: ClinVarAnnotator):
        """CLNREVSTAT is surfaced on the Annotation."""
        v = Variant("rs1801133", "1", 11796321, "G", "A")
        results = annotator.annotate(v)
        assert len(results) >= 1
        assert results[0].review_status == "criteria_provided,_single_submitter"


class TestRegistryMetadata:
    def test_class_attributes(self):
        assert ClinVarAnnotator.name == "clinvar"
        assert ClinVarAnnotator.display_name == "ClinVar"
        assert ClinVarAnnotator.attribution == "ClinVar"
        assert ClinVarAnnotator.requires_download is True


class TestIndelMatching:
    """M-4: ClinVar contains pathogenic indels (e.g., CFTR ΔF508). Must match."""

    def test_indel_carrier_triggers(self, annotator: ClinVarAnnotator):
        # mock ClinVar: rs113993960 REF=CTT ALT=C, Pathogenic CFTR
        v = Variant("rs113993960", "7", 117199644, "CTT", "C")
        results = annotator.annotate(v)
        assert len(results) == 1
        assert results[0].gene == "CFTR"
        # ADR-0023: indel diploid passes through as `"CTT/C"` to keep
        # multi-base alleles readable rather than concatenating them.
        assert results[0].genotype_match == "CTT/C"
        assert results[0].significance == "clinvar_pathogenic"

    def test_indel_homozygous_reference_does_not_trigger(self, annotator: ClinVarAnnotator):
        v = Variant("rs113993960", "7", 117199644, "CTT", "CTT")
        assert annotator.annotate(v) == []


class TestIndelAnchorProtection:
    """ADR-0011: indel rows must NOT fire on single-base array readouts.

    ClinVar encodes indels with anchor-base notation (REF=CTT ALT=C). Array
    parsers report single bases at probe positions. Pre-v0.4.2 the carrier
    rule's `alt in {allele1, allele2}` matched ClinVar's single-character
    anchor against an array's single-character readout, producing categorical
    false-positive "Pathogenic" calls in cancer-predisposition genes for
    users who carried only the wild-type sequence.
    """

    def test_array_single_base_does_not_fire_on_indel_row(self, annotator: ClinVarAnnotator):
        # mock fixture: rs113993960 REF=CTT ALT=C Pathogenic CFTR.
        # Array reads a single C at the probe position; the user does NOT
        # carry the deletion. Pre-v0.4.2 incorrectly fired.
        v = Variant("rs113993960", "7", 117199644, "C", "C")
        assert annotator.annotate(v) == []

    def test_indel_calling_parser_still_fires(self, annotator: ClinVarAnnotator):
        # A multi-base genotype like CTT/C indicates a parser that actually
        # calls indels (future VCF parser). Indel matching must still work.
        v = Variant("rs113993960", "7", 117199644, "CTT", "C")
        results = annotator.annotate(v)
        assert len(results) == 1
        assert results[0].significance == "clinvar_pathogenic"

    def test_homozygous_alt_indel_still_fires_for_multibase_parser(
        self, annotator: ClinVarAnnotator
    ):
        # Hypothetical homozygous deletion. The user's genotype carries the
        # multi-base form on at least one side, so the indel filter doesn't
        # short-circuit; the carrier rule still applies.
        v = Variant("rs113993960", "7", 117199644, "C", "CTT")
        results = annotator.annotate(v)
        assert len(results) == 1


class TestMultiAllelicMatching:
    """C-2: Multi-allelic ClinVar rows must match per-ALT, not as the joined string."""

    def test_carrier_of_pathogenic_alt_triggers(self, annotator: ClinVarAnnotator):
        # mock ClinVar: rs1065852 ALT=A,C with CLNSIG=Drug_response|Benign.
        # MHG fixture has G/A — carries A only.
        v = Variant("rs1065852", "22", 42526694, "G", "A")
        results = annotator.annotate(v)
        # Should match exactly the A-allele record (Drug_response), not the C one.
        sigs = {r.significance for r in results}
        assert "clinvar_drug_response" in sigs
        assert "clinvar_benign" not in sigs

    def test_carrier_of_benign_alt_only(self, annotator_with_benign: ClinVarAnnotator):
        # User carries G/C — only the C-allele record should fire (Benign).
        v = Variant("rs1065852", "22", 42526694, "G", "C")
        results = annotator_with_benign.annotate(v)
        sigs = {r.significance for r in results}
        assert "clinvar_benign" in sigs
        assert "clinvar_drug_response" not in sigs


class TestRemoteSignal:
    """Freshness signal: ClinVar uses the .md5 sidecar file (ADR-0012)."""

    def test_fetch_returns_md5_prefixed_signal(self, annotator: ClinVarAnnotator, monkeypatch):
        """ADR-0021: signal is composite across managed builds."""
        from allelix.annotators import clinvar as clinvar_module

        valid_md5 = "abcdef0123456789abcdef0123456789"  # 32 hex digits
        monkeypatch.setattr(
            clinvar_module,
            "fetch_remote_text",
            lambda url: f"{valid_md5}  clinvar.vcf.gz\n",
        )
        # Default annotator manages both builds → composite signal.
        assert annotator.fetch_remote_signal() == (
            f"GRCh37:md5:{valid_md5}|GRCh38:md5:{valid_md5}"
        )

    def test_fetch_returns_none_on_network_error(self, annotator: ClinVarAnnotator, monkeypatch):
        from allelix.annotators import clinvar as clinvar_module

        monkeypatch.setattr(clinvar_module, "fetch_remote_text", lambda url: None)
        assert annotator.fetch_remote_signal() is None

    def test_fetch_returns_none_on_empty_md5_body(self, annotator: ClinVarAnnotator, monkeypatch):
        from allelix.annotators import clinvar as clinvar_module

        monkeypatch.setattr(clinvar_module, "fetch_remote_text", lambda url: "   \n")
        assert annotator.fetch_remote_signal() is None

    def test_fetch_returns_none_on_html_error_page(self, annotator: ClinVarAnnotator, monkeypatch):
        """GH #21: a CDN 503 returning HTML must not be accepted as a hash.

        Before the fix the first whitespace-separated token of the body
        was treated as the md5, so ``<!DOCTYPE`` would become the
        "signal" and propagate to ``verify_file_hash``, which would then
        delete the freshly downloaded VCF on the resulting mismatch.
        """
        from allelix.annotators import clinvar as clinvar_module

        html = (
            "<!DOCTYPE html>\n<html><head><title>503 Service "
            "Unavailable</title></head><body>...</body></html>\n"
        )
        monkeypatch.setattr(clinvar_module, "fetch_remote_text", lambda url: html)
        assert annotator.fetch_remote_signal() is None

    def test_fetch_returns_none_on_short_hex(self, annotator: ClinVarAnnotator, monkeypatch):
        """GH #21: tokens that are hex but the wrong length are rejected."""
        from allelix.annotators import clinvar as clinvar_module

        monkeypatch.setattr(
            clinvar_module,
            "fetch_remote_text",
            lambda url: "abcdef01  clinvar.vcf.gz\n",  # only 8 hex chars
        )
        assert annotator.fetch_remote_signal() is None

    def test_fetch_returns_none_on_non_hex_token(self, annotator: ClinVarAnnotator, monkeypatch):
        """GH #21: 32-char-long but non-hex tokens are rejected."""
        from allelix.annotators import clinvar as clinvar_module

        monkeypatch.setattr(
            clinvar_module,
            "fetch_remote_text",
            lambda url: "zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz  clinvar.vcf.gz\n",
        )
        assert annotator.fetch_remote_signal() is None

    def test_cached_returns_none_for_unconfigured(self, tmp_path: Path):
        ann = ClinVarAnnotator(tmp_path)
        assert ann.cached_remote_signal() is None

    def test_cached_returns_none_for_v041_cache(self, annotator: ClinVarAnnotator):
        """v0.4.1 caches were populated without a remote_signal column."""
        # The clinvar_data_dir fixture writes via load_clinvar_vcf without
        # passing remote_signal, so the column exists (new schema) but the
        # value is NULL — cached_remote_signal should return None.
        assert annotator.cached_remote_signal() is None

    def test_cached_round_trip_after_setup(self, tmp_path: Path, mock_clinvar_vcf: Path):
        """ADR-0021: composite cached signal is `GRCh37:<sig>|GRCh38:<sig>`.

        For a single-build annotator the composite collapses to one part.
        """
        from allelix.annotators.clinvar import clinvar_db_filename, clinvar_record_name
        from allelix.databases.manager import load_clinvar_vcf

        load_clinvar_vcf(
            mock_clinvar_vcf,
            tmp_path / clinvar_db_filename("GRCh37"),
            source_url="test",
            remote_signal="md5:deadbeef",
            record_name=clinvar_record_name("GRCh37"),
        )
        ann = ClinVarAnnotator(tmp_path, builds=("GRCh37",))
        try:
            assert ann.cached_remote_signal() == "GRCh37:md5:deadbeef"
        finally:
            ann.close()


class TestConstructorValidation:
    def test_unsupported_build_raises(self, tmp_path: Path):
        with pytest.raises(ValueError, match="Unsupported"):
            ClinVarAnnotator(tmp_path, builds=("GRCh99",))


class TestCloseable:
    """C-1: ClinVarAnnotator must release its SQLite connections deterministically."""

    def test_close_releases_connection(self, annotator: ClinVarAnnotator):
        # Touch the connection
        annotator.annotate(Variant("rs1801133", "1", 11796321, "G", "A"))
        assert annotator._conns, "expected at least one open per-build connection"
        annotator.close()
        assert annotator._conns == {}

    def test_close_is_idempotent(self, annotator: ClinVarAnnotator):
        annotator.close()
        annotator.close()  # must not raise

    def test_context_manager_closes_on_exit(self, clinvar_data_dir: Path):
        ann = ClinVarAnnotator(clinvar_data_dir)
        with ann as bound:
            assert bound is ann
            ann.annotate(Variant("rs1801133", "1", 11796321, "G", "A"))
            assert ann._conns
        assert ann._conns == {}


class TestBatchAnnotateParity:
    """batch_annotate(vs) must return identical results to flatmap(annotate, vs).

    Pins the contract from base.py: a chunked SQL override must produce
    the same annotations in the same order as the per-variant path.
    """

    def test_parity_mixed_carriers_and_non_carriers(self, annotator: ClinVarAnnotator):
        variants = [
            Variant("rs1801133", "1", 11796321, "G", "A"),  # carrier (het, pathogenic)
            Variant("rs4680", "22", 19963748, "A", "A"),  # carrier (homo alt)
            Variant("rs121918506", "17", 7577538, "G", "G"),  # homo ref, no trigger
            Variant("rs999000111", "1", 1000, "A", "T"),  # unknown rsid
            Variant("rs1801133", "1", 11796321, "-", "-"),  # no-call
        ]
        per_variant = [a for v in variants for a in annotator.annotate(v)]
        batched = list(annotator.batch_annotate(variants))
        assert per_variant == batched

    def test_parity_empty_input(self, annotator: ClinVarAnnotator):
        assert list(annotator.batch_annotate([])) == []

    def test_parity_all_no_calls(self, annotator: ClinVarAnnotator):
        variants = [
            Variant("rs1801133", "1", 11796321, "-", "-"),
            Variant("rs4680", "22", 19963748, "-", "-"),
        ]
        per_variant = [a for v in variants for a in annotator.annotate(v)]
        batched = list(annotator.batch_annotate(variants))
        assert per_variant == batched == []

    def test_parity_duplicate_rsids(self, annotator: ClinVarAnnotator):
        """Same rsid appearing twice in input yields the right per-occurrence results."""
        variants = [
            Variant("rs1801133", "1", 11796321, "G", "A"),
            Variant("rs1801133", "1", 11796321, "G", "A"),
            Variant("rs1801133", "1", 11796321, "G", "G"),  # different genotype, no trigger
        ]
        per_variant = [a for v in variants for a in annotator.annotate(v)]
        batched = list(annotator.batch_annotate(variants))
        assert per_variant == batched

    def test_parity_mixed_builds(self, annotator: ClinVarAnnotator):
        """Mixed GRCh37 + GRCh38 variants in one batch hit their respective caches.

        ClinVar holds one SQLite cache per build (ADR-0021). The batch
        path groups variants by build internally so each group queries
        its dedicated connection.
        """
        variants = [
            Variant("rs1801133", "1", 11796321, "G", "A", build="GRCh37"),
            Variant("rs4680", "22", 19963748, "A", "A", build="GRCh37"),
            Variant("rs1801133", "1", 11796321, "G", "A", build="GRCh38"),
            Variant("rs4680", "22", 19963748, "A", "A", build="GRCh38"),
        ]
        per_variant = [a for v in variants for a in annotator.annotate(v)]
        batched = list(annotator.batch_annotate(variants))
        assert per_variant == batched

    def test_chunk_boundary_500_and_501_resolve(self, annotator: ClinVarAnnotator):
        """Position chunking at the _BATCH_CHUNK boundary stays correct.

        bulk_resolve_rsids chunks positions per chromosome at 500/query.
        500 → one query; 501 → two. Both must surface the real matches.
        """
        # 498 filler positions that won't match + 2 real
        boundary_500 = [Variant("", "1", 50_000_000 + i, "A", "T") for i in range(498)] + [
            Variant("", "1", 11856378, "G", "A"),
            Variant("", "22", 19951271, "G", "A"),
        ]
        assert len(boundary_500) == 500
        resolved_500 = annotator.bulk_resolve_rsids(boundary_500)
        # Two real positions resolve; chr22 above lives in a separate per-chrom
        # query so chunk arithmetic on chr1 alone is what's being exercised.
        assert ("1", 11856378, "G", "A") in resolved_500
        assert ("22", 19951271, "G", "A") in resolved_500

        boundary_501 = [*boundary_500, Variant("", "1", 50_000_999, "C", "G")]
        assert len(boundary_501) == 501
        resolved_501 = annotator.bulk_resolve_rsids(boundary_501)
        assert ("1", 11856378, "G", "A") in resolved_501
        assert ("22", 19951271, "G", "A") in resolved_501

    def test_chunk_boundary_500_and_501(self, annotator: ClinVarAnnotator):
        """SQL chunking at the _BATCH_CHUNK boundary (500 rsIDs) is correct.

        500 rsIDs → one SQL query; 501 → two. Both must yield the same
        result set as the per-variant path. Pins the chunking arithmetic
        and ensures the second-chunk path isn't broken.
        """
        # Most are filler rsIDs that won't match; a couple of real ones
        # carry the actual annotations.
        boundary_variants_500 = [
            Variant(f"rs90{i:07d}", "1", 1000 + i, "A", "T") for i in range(498)
        ] + [
            Variant("rs1801133", "1", 11796321, "G", "A"),
            Variant("rs4680", "22", 19963748, "A", "A"),
        ]
        assert len(boundary_variants_500) == 500
        boundary_variants_501 = [
            *boundary_variants_500,
            Variant("rs121918506", "17", 7577538, "G", "G"),
        ]
        assert len(boundary_variants_501) == 501

        per_variant_500 = [a for v in boundary_variants_500 for a in annotator.annotate(v)]
        batched_500 = list(annotator.batch_annotate(boundary_variants_500))
        assert per_variant_500 == batched_500

        per_variant_501 = [a for v in boundary_variants_501 for a in annotator.annotate(v)]
        batched_501 = list(annotator.batch_annotate(boundary_variants_501))
        assert per_variant_501 == batched_501


class TestBulkResolveRsids:
    """Resolves rsIDs by (chrom, pos) for VCFs from callers that emit ID=. (GH #8).

    Variant callers like GATK HaplotypeCaller and DeepVariant write `.` to the
    ID column. The annotation pipeline is rsID-keyed, so without resolution
    every rsID-keyed annotator returns zero hits on these files. This class
    pins the resolver's contract: position lookup, carrier-allele match,
    in-place rsid mutation, multi-build dispatch.
    """

    def test_single_variant_resolves(self, annotator: ClinVarAnnotator):
        """rsID-less variant at a known ClinVar position gets its rsid stamped."""
        v = Variant("", "1", 11856378, "G", "A", build="GRCh37")
        resolved = annotator.bulk_resolve_rsids([v])
        assert resolved == {("1", 11856378, "G", "A"): "rs1801133"}
        assert v.rsid == "rs1801133"

    def test_unknown_position_unchanged(self, annotator: ClinVarAnnotator):
        """Position not in ClinVar leaves the variant's rsid empty."""
        v = Variant("", "1", 99_999_999, "A", "T", build="GRCh37")
        resolved = annotator.bulk_resolve_rsids([v])
        assert resolved == {}
        assert v.rsid == ""

    def test_no_call_skipped(self, annotator: ClinVarAnnotator):
        """No-call variants don't participate in resolution."""
        v = Variant("", "1", 11856378, "-", "-", build="GRCh37")
        resolved = annotator.bulk_resolve_rsids([v])
        assert resolved == {}
        assert v.rsid == ""

    def test_empty_input(self, annotator: ClinVarAnnotator):
        """Empty input → empty dict, no SQL queries fired."""
        assert annotator.bulk_resolve_rsids([]) == {}

    def test_carrier_allele_must_match_clinvar_alt(self, annotator: ClinVarAnnotator):
        """Variant at the right position but carrying the wrong ALT doesn't resolve.

        ClinVar at 1:11856378 has REF=G, ALT=A. A user carrying G/C at that
        position is not the same variant and must NOT pick up rs1801133.
        """
        v = Variant("", "1", 11856378, "G", "C", build="GRCh37")
        resolved = annotator.bulk_resolve_rsids([v])
        assert resolved == {}
        assert v.rsid == ""

    def test_multi_allelic_disambiguation(self, annotator: ClinVarAnnotator):
        """Multi-allelic positions resolve to the row matching the carrier allele.

        Mock ClinVar at 22:42526694 has two ALT alleles: A (drug response, rs1065852
        from the synthetic row) and C (benign). A user carrying G/A must pick
        the A-row, not the C-row. Both rows share rsid in the synthetic data,
        so the test asserts on the (ref, alt) key shape.
        """
        v = Variant("", "22", 42526694, "G", "A", build="GRCh37")
        resolved = annotator.bulk_resolve_rsids([v])
        assert ("22", 42526694, "G", "A") in resolved
        assert v.rsid != ""

    def test_multi_build_dispatch(self, annotator: ClinVarAnnotator):
        """Each variant queries its own build's cache (ADR-0021).

        GRCh37 and GRCh38 caches hold the same variants at different
        coordinates (lift-over). Resolution must dispatch by build so each
        position lookup hits the right cache.
        """
        variants = [
            # GRCh37 coordinates for rs1801133 and rs4680
            Variant("", "1", 11856378, "G", "A", build="GRCh37"),
            # GRCh38 coordinates for the same two rsIDs (lifted over)
            Variant("", "1", 11796321, "G", "A", build="GRCh38"),
            Variant("", "22", 19963748, "G", "A", build="GRCh38"),
        ]
        resolved = annotator.bulk_resolve_rsids(variants)
        assert variants[0].rsid == "rs1801133"
        assert variants[1].rsid == "rs1801133"
        assert variants[2].rsid == "rs4680"
        assert ("1", 11856378, "G", "A") in resolved
        assert ("1", 11796321, "G", "A") in resolved
        assert ("22", 19963748, "G", "A") in resolved

    def test_empty_chromosome_skipped(self, annotator: ClinVarAnnotator):
        """Variants with empty chrom or non-positive position are skipped."""
        variants = [
            Variant("", "", 11856378, "G", "A", build="GRCh37"),
            Variant("", "1", 0, "G", "A", build="GRCh37"),
        ]
        resolved = annotator.bulk_resolve_rsids(variants)
        assert resolved == {}
        assert variants[0].rsid == ""
        assert variants[1].rsid == ""

    def test_post_resolution_batch_annotate_yields_hits(self, annotator: ClinVarAnnotator):
        """Once rsIDs are stamped, the existing batch_annotate path produces hits.

        End-to-end pin: this is the contract the pipeline depends on for
        rsID-less VCFs. Resolution must mutate variants such that the
        unchanged downstream annotators see real rsIDs.
        """
        v = Variant("", "1", 11856378, "G", "A", build="GRCh37")
        # Before resolution: no annotations
        assert list(annotator.batch_annotate([v])) == []
        annotator.bulk_resolve_rsids([v])
        # After resolution: annotations appear via the resolved rsid
        annotations = list(annotator.batch_annotate([v]))
        assert annotations
        assert all(a.rsid == "rs1801133" for a in annotations)

    def test_position_index_self_heals_existing_caches(self, annotator: ClinVarAnnotator):
        """First `_connection()` call creates idx_clinvar_position if missing.

        Older caches predate this index. The lazy CREATE INDEX IF NOT EXISTS
        in `_connection()` migrates them on first use without a full rebuild.
        """
        # Force connection open + index creation
        annotator._connection("GRCh37")
        conn = annotator._conns["GRCh37"]
        cursor = conn.execute(
            "SELECT name FROM sqlite_master WHERE type='index' AND name='idx_clinvar_position'"
        )
        assert cursor.fetchone() is not None

    def test_multi_allelic_resolution_is_deterministic(self, tmp_path):
        """Multi-allelic positions with different rsIDs resolve deterministically.

        In production ClinVar data rsID is typically shared across ALTs at
        a position, but the resolver's output must not depend on SQLite's
        physical row order. Sort by (ref, alt, rsid) before carrier match.
        """
        # Build a synthetic single-build ClinVar with two co-located indel
        # rows sharing an ALT but differing in REF anchor length — the
        # exact shape that the ALT-only resolver mis-routed before the
        # subset fix.
        from allelix.databases.schema import CLINVAR_SCHEMA

        db_path = tmp_path / clinvar_db_filename("GRCh37")
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            for stmt in CLINVAR_SCHEMA.split(";"):
                if stmt.strip():
                    conn.execute(stmt)
            # Two anchor-base deletions at the same position. With ALT-only
            # matching, both rows would match a user carrying ALT=C; sort
            # picks the wrong one. The subset check requires REF and ALT
            # to both be consistent with the user's diploid call.
            conn.executemany(
                "INSERT INTO clinvar_variants (rsid, chromosome, position, ref, alt, "
                "clinical_significance, condition, gene, review_status, allele_id) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                [
                    ("rs_DEL_A", "1", 5000, "CA", "C", "Benign", "X", "X", "ok", 1),
                    ("rs_DEL_AG", "1", 5000, "CAG", "C", "Pathogenic", "X", "X", "ok", 1),
                ],
            )
            conn.execute(
                "INSERT INTO database_versions (name, source_url, version, "
                "downloaded_at, record_count, local_version_tag) "
                "VALUES (?, ?, ?, ?, ?, ?)",
                (
                    "clinvar.GRCh37",
                    "test://mock",
                    "test",
                    "2026-01-01",
                    2,
                    f"iv:{CLINVAR_INTERPRETER_VERSION}",
                ),
            )
            conn.commit()

        ann = ClinVarAnnotator(tmp_path, builds=("GRCh37",))
        try:
            # User carries the longer deletion as 0/1 het: VCF REF=CAG,
            # ALT=C, GT=0/1 → parser yields allele1="CAG", allele2="C".
            # The correct rsid is rs_DEL_AG (Pathogenic), NOT rs_DEL_A
            # (Benign) — even though rs_DEL_A sorts first lexicographically
            # and shares the ALT "C".
            v = Variant("", "1", 5000, "CAG", "C", build="GRCh37")
            ann.bulk_resolve_rsids([v])
            assert v.rsid == "rs_DEL_AG", (
                f"expected rs_DEL_AG (REF=CAG matches user), got {v.rsid!r}"
            )

            # Re-run with a fresh open to confirm idempotence and that the
            # REF-aware match doesn't depend on SQLite row order.
            ann.close()
            v2 = Variant("", "1", 5000, "CAG", "C", build="GRCh37")
            ann.bulk_resolve_rsids([v2])
            assert v2.rsid == "rs_DEL_AG"
        finally:
            ann.close()

    def test_indel_anchor_homalt_abstains_when_ambiguous(self, tmp_path):
        """Hom-alt at a shared-anchor position abstains rather than coin-flip.

        Pairs the het-anchor test above. For a hom-alt user (allele1 ==
        allele2), the Variant model loses which REF the user carried —
        the subset check passes for ALL co-located rows because their REF
        is allowed to be anything when ``{ALT} ⊆ {REF, ALT}``. Without
        Variant.ref to break the tie, the resolver MUST abstain rather
        than pick by sort order. ``len(matches) > 1`` → no rsid stamped,
        variant flows through without rsID-keyed annotations. A future
        Variant.ref field (v2.1+) would close this residual.
        """
        from allelix.databases.schema import CLINVAR_SCHEMA

        db_path = tmp_path / clinvar_db_filename("GRCh37")
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            for stmt in CLINVAR_SCHEMA.split(";"):
                if stmt.strip():
                    conn.execute(stmt)
            conn.executemany(
                "INSERT INTO clinvar_variants (rsid, chromosome, position, ref, alt, "
                "clinical_significance, condition, gene, review_status, allele_id) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                [
                    ("rs_DEL_A", "1", 6000, "CA", "C", "Benign", "X", "X", "ok", 1),
                    ("rs_DEL_AG", "1", 6000, "CAG", "C", "Pathogenic", "X", "X", "ok", 1),
                ],
            )
            conn.execute(
                "INSERT INTO database_versions (name, source_url, version, "
                "downloaded_at, record_count, local_version_tag) "
                "VALUES (?, ?, ?, ?, ?, ?)",
                (
                    "clinvar.GRCh37",
                    "test://mock",
                    "test",
                    "2026-01-01",
                    2,
                    f"iv:{CLINVAR_INTERPRETER_VERSION}",
                ),
            )
            conn.commit()

        ann = ClinVarAnnotator(tmp_path, builds=("GRCh37",))
        try:
            # User hom-alt C/C at this position. Parser yields (C, C).
            # Both ClinVar rows pass the subset check ({C} ⊆ {CA, C} and
            # {C} ⊆ {CAG, C}) because the user's alleles are a subset of
            # either row's {REF, ALT}. With two ambiguous matches and no
            # REF info on Variant to break the tie, the resolver
            # abstains — better an unannotated variant than the wrong
            # rsid stamped onto it.
            v = Variant("", "1", 6000, "C", "C", build="GRCh37")
            ann.bulk_resolve_rsids([v])
            assert v.rsid == "", (
                f"hom-alt at shared-anchor must abstain (two ambiguous "
                f"subset matches), got {v.rsid!r}"
            )
        finally:
            ann.close()

    def test_multi_allelic_1_2_genotype_no_resolution(self, tmp_path):
        """1/2 multi-allelic skip: user carries two ALTs at a position where
        each ClinVar row has the position's REF. ``user_alleles ⊆ {ref, alt}``
        rejects all rows (REF isn't in the user's alleles), so no rsid is
        stamped. Conservative behavior — beats silently picking one of the
        rsids the user partially carries.
        """
        from allelix.databases.schema import CLINVAR_SCHEMA

        db_path = tmp_path / clinvar_db_filename("GRCh37")
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            for stmt in CLINVAR_SCHEMA.split(";"):
                if stmt.strip():
                    conn.execute(stmt)
            conn.executemany(
                "INSERT INTO clinvar_variants (rsid, chromosome, position, ref, alt, "
                "clinical_significance, condition, gene, review_status, allele_id) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                [
                    ("rs_AG", "1", 7000, "A", "G", "Pathogenic", "X", "X", "ok", 1),
                    ("rs_AT", "1", 7000, "A", "T", "Pathogenic", "X", "X", "ok", 1),
                ],
            )
            conn.execute(
                "INSERT INTO database_versions (name, source_url, version, "
                "downloaded_at, record_count, local_version_tag) "
                "VALUES (?, ?, ?, ?, ?, ?)",
                (
                    "clinvar.GRCh37",
                    "test://mock",
                    "test",
                    "2026-01-01",
                    2,
                    f"iv:{CLINVAR_INTERPRETER_VERSION}",
                ),
            )
            conn.commit()

        ann = ClinVarAnnotator(tmp_path, builds=("GRCh37",))
        try:
            # VCF: REF=A, ALT=G,T, GT=1/2 → parser yields (G, T).
            v = Variant("", "1", 7000, "G", "T", build="GRCh37")
            ann.bulk_resolve_rsids([v])
            assert v.rsid == "", (
                f"1/2 multi-allelic at multi-rsid position should not resolve "
                f"(REF info would be needed), got {v.rsid!r}"
            )
        finally:
            ann.close()

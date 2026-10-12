# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Cross-build rsID identity table (ADR-0040).

Covers the builder logic (RefSeq accession mapping, dbSNP line parsing,
cross-build allele pairing), the annotator's locus-keyed queries, and the
pipeline uses: GRCh37 enrichment through translation, ``Variant.ref`` on
GRCh37 arrays, and rsID recovery for rsID-less GRCh37 VCFs. All sites are
synthetic.
"""

from __future__ import annotations

import contextlib
import gzip
import shutil
import sqlite3
from typing import TYPE_CHECKING

import pytest

from allelix.annotators import get_annotator_class, get_annotators
from allelix.annotators.alphamissense import AlphaMissenseAnnotator
from allelix.annotators.cadd import CaddAnnotator
from allelix.annotators.clinvar import ClinVarAnnotator, clinvar_db_filename, clinvar_record_name
from allelix.annotators.gnomad import GnomadAnnotator
from allelix.annotators.rsid_identity import RsidIdentityAnnotator
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
from allelix.databases.rsid_identity_loader import (
    RSID_IDENTITY_DB_FILENAME,
    IdentityRow,
    build_identity_rows,
    install_prebuilt_cache,
    pair_alleles,
    parse_dbsnp_line,
    refseq_to_chrom,
    write_identity_db,
)
from allelix.databases.schema import (
    ALPHAMISSENSE_SCHEMA,
    CADD_SCHEMA,
    CLINVAR_SCHEMA,
    GNOMAD_SCHEMA,
)
from allelix.parsers.myhappygenes import MyHappyGenesParser
from allelix.parsers.vcf import VcfParser
from allelix.reports._pipeline import run_analysis

if TYPE_CHECKING:
    from pathlib import Path

    from allelix.models import Annotation

_SINGLE = "criteria_provided,_single_submitter"


# ── Builder logic ──────────────────────────────────────────────


class TestRefseqToChrom:
    @pytest.mark.parametrize(
        ("accession", "chrom"),
        [
            ("NC_000001.10", "1"),
            ("NC_000001.11", "1"),
            ("NC_000022.11", "22"),
            ("NC_000023.10", "X"),
            ("NC_000024.10", "Y"),
            ("NC_012920.1", "MT"),
        ],
    )
    def test_primary_accessions(self, accession: str, chrom: str) -> None:
        assert refseq_to_chrom(accession) == chrom

    @pytest.mark.parametrize("accession", ["NT_187361.1", "NW_003315905.1", "NC_000025.1", "X"])
    def test_non_primary_rejected(self, accession: str) -> None:
        assert refseq_to_chrom(accession) is None


class TestParseDbsnpLine:
    wanted = frozenset({"rs1", "rs2"})

    def test_wanted_primary_placement(self) -> None:
        line = "NC_000001.11\t100\trs1\tc\tA,G\t.\t.\tVC=SNV\n"
        assert parse_dbsnp_line(line, self.wanted) == ("rs1", "1", 100, "C", ("A", "G"))

    @pytest.mark.parametrize(
        "line",
        [
            "##fileformat=VCFv4.2\n",
            "NC_000001.11\t100\trs3\tC\tA\t.\t.\t.\n",  # not wanted
            "NT_187361.1\t100\trs1\tC\tA\t.\t.\t.\n",  # not primary
            "NC_000001.11\tx\trs1\tC\tA\t.\t.\t.\n",  # bad position
            "NC_000001.11\t100\trs1\tC\t.\t.\t.\t.\n",  # no alt
            "NC_000001.11\t100\trs1\n",  # truncated
        ],
    )
    def test_rejected_lines(self, line: str) -> None:
        assert parse_dbsnp_line(line, self.wanted) is None


class TestPairAlleles:
    def test_same_strand(self) -> None:
        rows, palindromic = pair_alleles("C", ("T",), "C", ("T",))  # type: ignore[misc]
        assert rows == [("C", "T", "C", "T")]
        assert not palindromic

    def test_strand_inverted_snv(self) -> None:
        """GRCh37 C>G is GRCh38 G>C on the inverted strand: G pairs with C."""
        rows, palindromic = pair_alleles("C", ("A",), "G", ("T",))  # type: ignore[misc]
        assert rows == [("C", "A", "G", "T")]
        assert not palindromic

    def test_reference_swap_keeps_both_alleles(self) -> None:
        """The GRCh38 reference adopted the GRCh37 alt: both alleles get a row."""
        rows, _ = pair_alleles("A", ("G",), "G", ("A",))  # type: ignore[misc]
        assert rows == [("A", "A", "G", "A"), ("A", "G", "G", "G")]

    def test_palindromic_is_flagged(self) -> None:
        rows, palindromic = pair_alleles("A", ("T",), "A", ("T",))  # type: ignore[misc]
        assert rows == [("A", "T", "A", "T")]
        assert palindromic

    def test_indel_same_strand(self) -> None:
        rows, palindromic = pair_alleles("CA", ("C",), "CA", ("C",))  # type: ignore[misc]
        assert rows == [("CA", "C", "CA", "C")]
        assert not palindromic

    def test_indel_on_inverted_strand_is_unpairable(self) -> None:
        assert pair_alleles("CA", ("C",), "TG", ("G",)) is None

    def test_different_allele_sets_are_unpairable(self) -> None:
        assert pair_alleles("C", ("T",), "C", ("A",)) is None


def test_build_identity_rows_counts_exclusions() -> None:
    p37 = {
        "rs1": [("rs1", "1", 100, "C", ("T",))],
        "rs2": [("rs2", "1", 200, "C", ("T",)), ("rs2", "1", 300, "C", ("T",))],
        "rs3": [("rs3", "1", 400, "C", ("T",))],
        "rs4": [("rs4", "1", 500, "C", ("T",))],
        "rs5": [("rs5", "1", 600, "C", ("T",))],
    }
    p38 = {
        "rs1": [("rs1", "1", 1100, "C", ("T",))],
        "rs2": [("rs2", "1", 1200, "C", ("T",))],
        "rs3": [("rs3", "2", 1400, "C", ("T",))],
        "rs4": [("rs4", "1", 1500, "C", ("A",))],
        "rs6": [("rs6", "1", 1600, "C", ("T",))],
    }
    stats: dict[str, int] = {}
    rows = list(build_identity_rows(p37, p38, stats))
    assert rows == [("rs1", "1", 100, "C", "T", 1100, "C", "T", 0)]
    assert stats == {
        "rsids": 1,
        "rows": 1,
        "multi_placement": 1,
        "chrom_changed": 1,
        "unpairable": 1,
        "only_grch37": 1,
        "only_grch38": 1,
    }


# ── Synthetic identity cache ───────────────────────────────────

# (rsid, chrom, pos37, ref37, alt37, pos38, ref38, alt38, palindromic)
_ROWS: tuple[IdentityRow, ...] = (
    # Same strand.
    ("rs900173001", "1", 1_000_100, "G", "A", 1_100_100, "G", "A", 0),
    # Strand inverted between builds: GRCh37 C>A is GRCh38 G>T.
    ("rs900173002", "1", 1_000_200, "C", "A", 1_100_200, "G", "T", 0),
    # Palindromic: translation must abstain.
    ("rs900173003", "1", 1_000_300, "A", "T", 1_100_300, "A", "T", 1),
    # Reference swap: GRCh37 alt G is the GRCh38 reference.
    ("rs900173004", "1", 1_000_400, "A", "A", 1_100_400, "G", "A", 0),
    ("rs900173004", "1", 1_000_400, "A", "G", 1_100_400, "G", "G", 0),
    # Multi-allelic same strand (rsID-less resolution of a 1/2 call).
    ("rs900173005", "2", 2_000_500, "C", "A", 2_100_500, "C", "A", 0),
    ("rs900173005", "2", 2_000_500, "C", "T", 2_100_500, "C", "T", 0),
)


def _write_identity(data_dir: Path) -> None:
    write_identity_db(
        _ROWS,
        data_dir / RSID_IDENTITY_DB_FILENAME,
        version="dbSNP test",
        source_url="test://identity",
        downloaded_at="2026-10-10T00:00:00Z",
    )


@pytest.fixture
def identity(tmp_path: Path, cm_stack: contextlib.ExitStack) -> RsidIdentityAnnotator:
    _write_identity(tmp_path)
    return cm_stack.enter_context(RsidIdentityAnnotator(tmp_path))


class TestAnnotator:
    def test_registered_and_ready(self, identity: RsidIdentityAnnotator) -> None:
        assert get_annotator_class("rsid_identity") is RsidIdentityAnnotator
        assert identity.is_ready()
        assert identity.version() == "dbSNP test"
        assert identity.record_count() == len(_ROWS)
        assert identity.annotate(None) == []  # type: ignore[arg-type]
        assert identity.fetch_remote_signal() is None
        assert identity.cached_remote_signal() is None

    def test_in_default_annotator_set(self, tmp_path: Path) -> None:
        names = [a.name for a in get_annotators(tmp_path)]
        assert "rsid_identity" in names
        for a in get_annotators(tmp_path):
            a.close()

    def test_palindromic_provider_wired_to_pharmgkb(self, tmp_path: Path) -> None:
        """get_annotators wires rsid_identity.palindromic_rsids
        into PharmGKBAnnotator so palindromic sites are guarded."""
        annotators = get_annotators(tmp_path)
        try:
            pharmgkb = next(a for a in annotators if a.name == "pharmgkb")
            assert pharmgkb._palindromic_provider is not None
            # Must be bound to the rsid_identity instance, not some other callable
            rsid_identity = next(a for a in annotators if a.name == "rsid_identity")
            assert pharmgkb._palindromic_provider == rsid_identity.palindromic_rsids
        finally:
            for a in annotators:
                a.close()

    def test_missing_cache_is_not_ready(self, tmp_path: Path) -> None:
        with RsidIdentityAnnotator(tmp_path) as missing:
            assert not missing.is_ready()
            assert missing.version() is None
            assert missing.record_count() is None
            with pytest.raises(FileNotFoundError):
                missing.refs_at_loci("GRCh37", {("1", 1)})

    def test_setup_refuses_without_a_pinned_cache(self, tmp_path: Path, monkeypatch) -> None:
        monkeypatch.setattr("allelix.annotators.rsid_identity.RSID_IDENTITY_EXPECTED_SHA256", "")
        with (
            RsidIdentityAnnotator(tmp_path) as ann,
            pytest.raises(RuntimeError, match="no published cache"),
        ):
            ann.setup()

    def test_setup_downloads_verifies_and_installs(self, tmp_path: Path, monkeypatch) -> None:
        build_dir = tmp_path / "build"
        build_dir.mkdir()
        _write_identity(build_dir)
        gz = build_dir / "rsid_identity.sqlite.gz"
        with (build_dir / RSID_IDENTITY_DB_FILENAME).open("rb") as f_in, gzip.open(gz, "wb") as f:
            shutil.copyfileobj(f_in, f)
        monkeypatch.setattr(
            "allelix.annotators.rsid_identity.download",
            lambda _url, dest: shutil.copyfile(gz, dest),
        )
        verified: list[str] = []
        monkeypatch.setattr(
            "allelix.annotators.rsid_identity.verify_file_hash",
            lambda _p, algo, digest: verified.append(f"{algo}:{digest}"),
        )
        monkeypatch.setattr(
            "allelix.annotators.rsid_identity.RSID_IDENTITY_EXPECTED_SHA256", "ab" * 32
        )
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        with RsidIdentityAnnotator(data_dir) as ann:
            ann.setup()
            assert ann.is_ready()
        assert verified == [f"sha256:{'ab' * 32}"]
        assert not (data_dir / "rsid_identity.sqlite.gz").exists()

    def test_install_prebuilt_cache_stamps_schema(self, tmp_path: Path) -> None:
        _write_identity(tmp_path)
        gz = tmp_path / "x.gz"
        with (tmp_path / RSID_IDENTITY_DB_FILENAME).open("rb") as f_in, gzip.open(gz, "wb") as f:
            shutil.copyfileobj(f_in, f)
        dest = tmp_path / "out"
        dest.mkdir()
        install_prebuilt_cache(gz, dest / RSID_IDENTITY_DB_FILENAME, source_url="test://x")
        with RsidIdentityAnnotator(dest) as ann:
            assert ann.is_ready()

    def test_refs_at_loci_per_build(self, identity: RsidIdentityAnnotator) -> None:
        assert identity.refs_at_loci("GRCh37", {("1", 1_000_200)}) == {("1", 1_000_200): {"C"}}
        assert identity.refs_at_loci("GRCh38", {("1", 1_100_200)}) == {("1", 1_100_200): {"G"}}
        assert identity.refs_at_loci("GRCh36", {("1", 1_000_200)}) == {}
        assert identity.refs_at_loci("GRCh37", set()) == {}

    def test_alleles_at_loci(self, identity: RsidIdentityAnnotator) -> None:
        got = identity.alleles_at_loci("GRCh37", {("2", 2_000_500)})
        assert sorted(got[("2", 2_000_500)]) == [
            ("C", "A", "rs900173005"),
            ("C", "T", "rs900173005"),
        ]
        assert identity.alleles_at_loci("GRCh36", {("2", 2_000_500)}) == {}

    def test_translate_to_grch38(self, identity: RsidIdentityAnnotator) -> None:
        keys = {
            ("1", 1_000_100, "G", "A"),  # same strand
            ("1", 1_000_200, "C", "A"),  # inverted strand
            ("1", 1_000_300, "A", "T"),  # palindromic: abstain
            ("1", 1_000_400, "A", "G"),  # GRCh38 reference: abstain
            ("1", 1_000_400, "A", "A"),  # GRCh38 alt A after the swap
            ("1", 1_000_100, "G", "C"),  # allele not in the table
        }
        assert identity.translate_to_grch38(keys) == {
            ("1", 1_000_100, "G", "A"): ("1", 1_100_100, "G", "A"),
            ("1", 1_000_200, "C", "A"): ("1", 1_100_200, "G", "T"),
            ("1", 1_000_400, "A", "A"): ("1", 1_100_400, "G", "A"),
        }

    def test_many_loci_are_chunked(self, identity: RsidIdentityAnnotator) -> None:
        loci = {("1", 5_000_000 + i) for i in range(1000)} | {("1", 1_000_100)}
        assert identity.refs_at_loci("GRCh37", loci) == {("1", 1_000_100): {"G"}}


# ── Pipeline: GRCh37 input with the identity table ─────────────


def _clinvar37(data_dir: Path) -> None:
    rows = [
        ("rs900173001", "1", 1_000_100, "G", "A", "Pathogenic", "Ident_A", "IDA", _SINGLE, 1),
        ("rs900173002", "1", 1_000_200, "C", "A", "Pathogenic", "Ident_B", "IDB", _SINGLE, 2),
        ("rs900173003", "1", 1_000_300, "A", "T", "Pathogenic", "Ident_C", "IDC", _SINGLE, 3),
        ("rs900173005", "2", 2_000_500, "C", "A", "Pathogenic", "Ident_E", "IDE", _SINGLE, 5),
        ("rs900173005", "2", 2_000_500, "C", "T", "Pathogenic", "Ident_E2", "IDE", _SINGLE, 6),
    ]
    with contextlib.closing(sqlite3.connect(data_dir / clinvar_db_filename("GRCh37"))) as conn:
        conn.executescript(CLINVAR_SCHEMA)
        conn.executemany(
            "INSERT INTO clinvar_variants (rsid, chromosome, position, ref, alt, "
            "clinical_significance, condition, gene, review_status, allele_id) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            rows,
        )
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
            "record_count, remote_signal, local_version_tag) "
            "VALUES (?, 'test://id', '1', '2026-10-10', ?, NULL, ?)",
            (clinvar_record_name("GRCh37"), len(rows), f"iv:{CLINVAR_INTERPRETER_VERSION}"),
        )
        conn.commit()


def _enrichment38(data_dir: Path) -> None:
    """GRCh38-keyed gnomAD / AlphaMissense / CADD rows for the translated sites."""
    with contextlib.closing(sqlite3.connect(data_dir / "gnomad.sqlite")) as conn:
        conn.executescript(GNOMAD_SCHEMA)
        conn.executemany(
            "INSERT INTO gnomad_frequencies (chrom, pos, ref, alt, rsid, af) VALUES (?,?,?,?,?,?)",
            [
                ("1", 1_100_100, "G", "A", "rs900173001", 0.11),
                ("1", 1_100_200, "G", "T", "rs900173002", 0.22),
                ("1", 1_100_300, "A", "T", "rs900173003", 0.33),
                # Decoy: the GRCh37 coordinate numerically, in GRCh38 space.
                ("1", 1_000_100, "G", "A", "rs_decoy", 0.99),
            ],
        )
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
            "record_count, local_version_tag) VALUES ('gnomad', 't', '4.1', 'x', 4, 'sv:1')"
        )
        conn.commit()
    with contextlib.closing(sqlite3.connect(data_dir / "alphamissense.sqlite")) as conn:
        conn.executescript(ALPHAMISSENSE_SCHEMA)
        conn.execute(
            "INSERT INTO alphamissense_scores VALUES "
            "('1', 1100200, 'G', 'T', 'rs900173002', 'P1', 'T1', 'A1V', 0.91, 'likely_pathogenic')"
        )
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
            "record_count, local_version_tag) "
            "VALUES ('alphamissense', 't', '2023', 'x', 1, 'sv:1')"
        )
        conn.commit()
    with contextlib.closing(sqlite3.connect(data_dir / "cadd.sqlite")) as conn:
        conn.executescript(CADD_SCHEMA)
        conn.execute(
            "INSERT INTO cadd_scores (chrom, pos, ref, alt, phred) VALUES "
            "('1', 1100100, 'G', 'A', 23.5)"
        )
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
            "record_count, local_version_tag) VALUES ('cadd', 't', 'v1.7', 'x', 1, 'sv:1')"
        )
        conn.commit()


@pytest.fixture
def grch37_data_dir(tmp_path: Path) -> Path:
    data = tmp_path / "data"
    data.mkdir()
    _clinvar37(data)
    _enrichment38(data)
    _write_identity(data)
    return data


_CHIP37 = """\
# MyHappyGenes [TEMPUS]
# Sample ID\tID000001
SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward
rs900173001\t1\t1000100\tG\tA
rs900173002\t1\t1000200\tC\tA
rs900173003\t1\t1000300\tA\tT
"""


def _run_chip37(
    tmp_path: Path, data: Path, cm_stack: contextlib.ExitStack, *, with_identity: bool
) -> tuple[dict[str, Annotation], list[tuple[str, str | None]]]:
    chip = tmp_path / "chip37.txt"
    chip.write_text(_CHIP37)
    result = run_analysis(
        chip,
        MyHappyGenesParser(),
        [cm_stack.enter_context(ClinVarAnnotator(data, builds=("GRCh37",)))],
        build_override="GRCh37",
        gnomad=cm_stack.enter_context(GnomadAnnotator(data)),
        alphamissense=cm_stack.enter_context(AlphaMissenseAnnotator(data)),
        cadd=cm_stack.enter_context(CaddAnnotator(data)),
        identity=cm_stack.enter_context(RsidIdentityAnnotator(data)) if with_identity else None,
    )
    return {a.rsid: a for a in result.annotations}, result.annotators_used


class TestGrch37Enrichment:
    def test_translated_keys_enrich_grch37_annotations(
        self, tmp_path: Path, grch37_data_dir: Path, cm_stack: contextlib.ExitStack
    ) -> None:
        by_rsid, used = _run_chip37(tmp_path, grch37_data_dir, cm_stack, with_identity=True)
        same = by_rsid["rs900173001"]
        assert same.allele_frequency == pytest.approx(0.11)  # not the 0.99 decoy
        assert same.cadd_phred == pytest.approx(23.5)
        inverted = by_rsid["rs900173002"]
        assert inverted.allele_frequency == pytest.approx(0.22)
        assert inverted.am_pathogenicity == pytest.approx(0.91)
        assert inverted.am_class == "likely_pathogenic"
        names = {name for name, _version in used}
        assert {"rsid_identity", "gnomad", "alphamissense", "cadd"} <= names

    def test_palindromic_site_abstains(
        self, tmp_path: Path, grch37_data_dir: Path, cm_stack: contextlib.ExitStack
    ) -> None:
        by_rsid, _ = _run_chip37(tmp_path, grch37_data_dir, cm_stack, with_identity=True)
        pal = by_rsid["rs900173003"]
        assert pal.allele_frequency is None
        assert pal.cadd_phred is None

    def test_without_identity_grch37_stays_unenriched(
        self, tmp_path: Path, grch37_data_dir: Path, cm_stack: contextlib.ExitStack
    ) -> None:
        by_rsid, used = _run_chip37(tmp_path, grch37_data_dir, cm_stack, with_identity=False)
        assert by_rsid
        assert all(a.allele_frequency is None for a in by_rsid.values())
        assert all(a.cadd_phred is None for a in by_rsid.values())
        names = {name for name, _version in used}
        assert names.isdisjoint({"rsid_identity", "gnomad", "alphamissense", "cadd"})

    def test_grch37_chip_gets_variant_ref(
        self, tmp_path: Path, grch37_data_dir: Path, cm_stack: contextlib.ExitStack
    ) -> None:
        chip = tmp_path / "chip37.txt"
        chip.write_text(_CHIP37)
        identity = cm_stack.enter_context(RsidIdentityAnnotator(grch37_data_dir))
        result = run_analysis(
            chip,
            MyHappyGenesParser(),
            [],
            build_override="GRCh37",
            identity=identity,
            high_value_rsids={"rs900173001", "rs900173002"},
        )
        assert {v.rsid: v.ref for v in result.hv_variants} == {
            "rs900173001": "G",
            "rs900173002": "C",
        }
        assert result.no_ref_variant_count == 0


_VCF37 = """\
##fileformat=VCFv4.2
##contig=<ID=1,length=249250621,assembly=GRCh37>
##contig=<ID=2,length=243199373,assembly=GRCh37>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
1\t1000100\t.\tG\tA\t50\tPASS\t.\tGT\t0/1
2\t2000500\t.\tC\tA,T\t50\tPASS\t.\tGT\t1/2
1\t1000400\t.\tA\tG\t50\tPASS\t.\tGT\t0/1
"""


def test_rsidless_grch37_vcf_resolves_through_identity(
    tmp_path: Path, cm_stack: contextlib.ExitStack
) -> None:
    """No ClinVar, no gnomAD: rsIDs come from the identity table alone, on GRCh37."""
    data = tmp_path / "data"
    data.mkdir()
    _write_identity(data)
    vcf = tmp_path / "rsidless37.vcf"
    vcf.write_text(_VCF37)
    result = run_analysis(
        vcf,
        VcfParser(),
        [],
        build_override="GRCh37",
        identity=cm_stack.enter_context(RsidIdentityAnnotator(data)),
        panel_rsids=frozenset({"rs900173001", "rs900173005", "rs900173004"}),
    )
    coverage = result.panel_coverage()
    assert coverage is not None
    assert coverage["missing"] == []
    assert result.rsidless_variant_count == 0


# ── Low-memory build pipeline ──────────────────────────────────


class TestRsidBloom:
    def test_no_false_negatives(self) -> None:
        from allelix.databases.rsid_identity_loader import RsidBloom

        bloom = RsidBloom(bits=20)
        members = [f"rs{n}" for n in range(1, 200_000, 7)]
        for rsid in members:
            bloom.add(rsid)
        assert all(rsid in bloom for rsid in members)

    def test_false_positive_rate_is_small(self) -> None:
        from allelix.databases.rsid_identity_loader import RsidBloom

        bloom = RsidBloom(bits=20)
        for n in range(0, 100_000, 2):
            bloom.add(f"rs{n}")
        outsiders = [f"rs{n}" for n in range(1_000_001, 1_100_001)]
        assert sum(r in bloom for r in outsiders) / len(outsiders) < 0.01

    @pytest.mark.parametrize("value", ["", "i5000123", "rs", "rs12a", 12, None])
    def test_non_rsids_are_never_members(self, value: object) -> None:
        from allelix.databases.rsid_identity_loader import RsidBloom

        assert value not in RsidBloom(bits=10)


_DBSNP37 = """\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
NC_000001.10\t100\trs1\tC\tT\t.\t.\tVC=SNV
NC_000001.10\t200\trs2\tC\tT\t.\t.\tVC=SNV
NT_187361.1\t300\trs1\tC\tT\t.\t.\tVC=SNV
NC_000001.10\t400\trs9\tC\tT\t.\t.\tVC=SNV
NC_000002.11\t500\trs3\tC\tA\t.\t.\tVC=SNV
NC_000001.10\t600\trs4\tA\tT\t.\t.\tVC=SNV
"""
_DBSNP38 = """\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
NC_000001.11\t1100\trs1\tC\tT\t.\t.\tVC=SNV
NC_000002.12\t2500\trs3\tG\tT\t.\t.\tVC=SNV
NC_000001.11\t1600\trs4\tA\tT\t.\t.\tVC=SNV
NC_000001.11\t1700\trs5\tG\tA\t.\t.\tVC=SNV
"""


def _sqlite_pipeline(universe: list[str], bloom_bits: int = 20) -> tuple[list, dict]:
    import io

    from allelix.databases.rsid_identity_loader import (
        add_universe,
        bloom_from_universe,
        create_work_db,
        filter_dbsnp_stream,
        grouped_placements,
        load_placements,
        pair_grouped_placements,
    )

    with contextlib.closing(sqlite3.connect(":memory:")) as work:
        create_work_db(work)
        add_universe(work, universe)
        bloom = bloom_from_universe(work, bits=bloom_bits)
        for build, text in (("GRCh37", _DBSNP37), ("GRCh38", _DBSNP38)):
            out = io.StringIO()
            filter_dbsnp_stream(text.splitlines(keepends=True), bloom, out)
            load_placements(work, build, out.getvalue().splitlines(keepends=True))
        stats: dict[str, int] = {}
        return list(pair_grouped_placements(grouped_placements(work), stats)), stats


class TestSqlitePipeline:
    def test_matches_in_memory_reference(self) -> None:
        universe = ["rs1", "rs2", "rs3", "rs4", "rs5", "i123", "rs1"]
        rows, stats = _sqlite_pipeline(universe)
        wanted = {"rs1", "rs2", "rs3", "rs4", "rs5"}
        p37: dict = {}
        p38: dict = {}
        for text, target in ((_DBSNP37, p37), (_DBSNP38, p38)):
            for line in text.splitlines(keepends=True):
                parsed = parse_dbsnp_line(line, wanted)
                if parsed:
                    target.setdefault(parsed[0], []).append(parsed)
        ref_stats: dict[str, int] = {}
        assert rows == list(build_identity_rows(p37, p38, ref_stats))
        assert stats == ref_stats
        assert rows == [
            ("rs1", "1", 100, "C", "T", 1100, "C", "T", 0),
            ("rs3", "2", 500, "C", "A", 2500, "G", "T", 0),
            ("rs4", "1", 600, "A", "T", 1600, "A", "T", 1),
        ]
        assert stats["only_grch37"] == 1  # rs2
        assert stats["only_grch38"] == 1  # rs5

    def test_bloom_false_positives_are_removed_exactly(self) -> None:
        """A 1-bit Bloom filter admits everything; the universe join must not."""
        rows, _ = _sqlite_pipeline(["rs1"], bloom_bits=3)
        assert [r[0] for r in rows] == ["rs1"]


def test_build_script_end_to_end(tmp_path: Path) -> None:
    """The script builds a working cache from tiny caches and dbSNP files."""
    import subprocess
    import sys
    from pathlib import Path as _Path

    data = tmp_path / "data"
    data.mkdir()
    for filename, table in (
        (clinvar_db_filename("GRCh37"), "clinvar_variants"),
        (clinvar_db_filename("GRCh38"), "clinvar_variants"),
        ("pharmgkb.sqlite", "pharmgkb_annotations"),
        ("gwas.sqlite", "gwas_associations"),
        ("snpedia.sqlite", "snpedia_genotypes"),
    ):
        with contextlib.closing(sqlite3.connect(data / filename)) as conn:
            conn.execute(f"CREATE TABLE {table} (rsid TEXT)")
            conn.executemany(f"INSERT INTO {table} VALUES (?)", [("rs1",), ("rs3",), ("rs4",)])
            conn.commit()
    sources = {}
    for name, text in (("d37.vcf.gz", _DBSNP37), ("d38.vcf.gz", _DBSNP38)):
        path = tmp_path / name
        with gzip.open(path, "wt") as fh:
            fh.write(text)
        sources[name] = path
    out = tmp_path / "out" / "rsid_identity.sqlite"
    out.parent.mkdir()
    script = _Path(__file__).resolve().parents[1] / "scripts" / "build_rsid_identity_cache.py"
    proc = subprocess.run(
        [
            sys.executable,
            str(script),
            "--data-dir",
            str(data),
            "--grch37-vcf",
            str(sources["d37.vcf.gz"]),
            "--grch38-vcf",
            str(sources["d38.vcf.gz"]),
            "--output",
            str(out),
            "--gzip",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    assert "sha256=" in proc.stderr
    assert (out.parent / "rsid_identity.sqlite.gz").exists()
    assert not (out.parent / "rsid_identity_work").exists()
    with RsidIdentityAnnotator(out.parent) as ann:
        assert ann.is_ready()
        assert ann.record_count() == 3
        assert ann.translate_to_grch38({("2", 500, "C", "A")}) == {
            ("2", 500, "C", "A"): ("2", 2500, "G", "T")
        }

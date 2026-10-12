# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Allele identity: Annotation invariants and locus-keyed ClinVar selection.

1. An Annotation never names an ``alt`` the user carries zero copies of
   (every source).
2. ClinVar, ClinPGx and SNPedia never emit ``alt`` equal to the reference
   allele. GWAS is exempt: its ``alt`` is the study's effect allele.

Plus the field defect behind invariant 2: ClinVar publishes reference-
allele records (HGVS ``=``, e.g. ``c.1601=``) whose VCF columns carry
REF == ALT. Matched on rsID against a homozygous-reference chip call they
produced "Pathogenic" rows with ``alt`` equal to the reference base. The
sites and rsIDs below are synthetic; only the row shape mirrors ClinVar.

ClinVar rows are selected by ``(chrom, pos, ref, alt)`` on every path,
never by rsID alone, and every position resolver applies one rule
(``resolve_rsid_at_locus``).
"""

from __future__ import annotations

import contextlib
import sqlite3
from typing import TYPE_CHECKING

import pytest

from allelix.annotators.clinvar import (
    ClinVarAnnotator,
    clinvar_db_filename,
    clinvar_record_name,
)
from allelix.annotators.gwas import GWASCatalogAnnotator
from allelix.annotators.pharmgkb import PharmGKBAnnotator
from allelix.annotators.snpedia import SNPediaAnnotator
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
from allelix.databases.schema import CLINVAR_SCHEMA
from allelix.models import (
    ALT_IS_NON_REFERENCE_SOURCES,
    Annotation,
    AnnotationInvariantError,
    Variant,
    annotation_to_public_dict,
)
from allelix.parsers.myhappygenes import MyHappyGenesParser
from allelix.reports._pipeline import run_analysis
from allelix.utils.allele import carried_in_source_orientation, resolve_rsid_at_locus

if TYPE_CHECKING:
    from pathlib import Path


def _annotation(**overrides: object) -> Annotation:
    fields: dict[str, object] = {
        "source": "clinvar",
        "rsid": "rs1",
        "significance": "clinvar_pathogenic",
        "category": "clinical",
        "magnitude": 9.0,
        "description": "d",
        "attribution": "ClinVar",
        "genotype_match": "AG",
        "alt": "G",
        "locus_ref": "A",
        "carried_alleles": ("A", "G"),
    }
    fields.update(overrides)
    return Annotation(**fields)  # type: ignore[arg-type]


class TestConstructorInvariants:
    def test_valid_annotation_constructs(self) -> None:
        assert _annotation().alt == "G"

    @pytest.mark.parametrize("source", ["clinvar", "pharmgkb", "snpedia", "gwas"])
    def test_zero_carried_copies_rejected_for_every_source(self, source: str) -> None:
        with pytest.raises(AnnotationInvariantError, match="not carried"):
            _annotation(source=source, alt="T", carried_alleles=("C", "C"), locus_ref="C")

    @pytest.mark.parametrize("source", sorted(ALT_IS_NON_REFERENCE_SOURCES))
    def test_alt_equal_ref_rejected(self, source: str) -> None:
        with pytest.raises(AnnotationInvariantError, match="reference allele"):
            _annotation(source=source, alt="C", locus_ref="C", carried_alleles=("C", "C"))

    def test_gwas_effect_allele_may_be_reference(self) -> None:
        a = _annotation(
            source="gwas",
            significance="gwas_association",
            alt="C",
            locus_ref="C",
            carried_alleles=("C", "T"),
        )
        assert a.alt == "C"

    def test_gwas_is_the_only_exempt_source(self) -> None:
        assert frozenset({"clinvar", "pharmgkb", "snpedia"}) == ALT_IS_NON_REFERENCE_SOURCES

    def test_unknown_context_is_not_checked(self) -> None:
        """Rows without context (diff baselines, hand-built test rows) construct."""
        a = _annotation(alt="T", locus_ref=None, carried_alleles=())
        assert a.alt == "T"

    def test_empty_alt_is_not_checked(self) -> None:
        a = _annotation(alt="", locus_ref="A", carried_alleles=("A", "A"))
        assert a.alt == ""

    def test_context_fields_never_serialize(self) -> None:
        d = annotation_to_public_dict(_annotation())
        assert "locus_ref" not in d
        assert "carried_alleles" not in d

    def test_context_fields_do_not_affect_equality(self) -> None:
        assert _annotation() == _annotation(locus_ref=None, carried_alleles=())


class TestCarriedInSourceOrientation:
    def test_direct_match_returns_alleles_unchanged(self) -> None:
        assert carried_in_source_orientation("A", "G", "G") == ("A", "G")

    def test_flipped_match_returns_complement(self) -> None:
        assert carried_in_source_orientation("T", "C", "G") == ("A", "G")


# ── ClinVar reference-allele records ───────────────────────────

_SINGLE = "criteria_provided,_single_submitter"

# (rsid, chrom, pos, ref, alt, clnsig, condition, gene, review, allele_id)
_REFERENCE_RECORD_ROWS: tuple[tuple[object, ...], ...] = (
    # Multi-allelic pathogenic site with a reference-allele record ("=").
    ("rs900170001", "1", 1_700_100, "C", "T", "Pathogenic", "Ref_rec_A", "REFA", _SINGLE, 1),
    ("rs900170001", "1", 1_700_100, "C", "A", "Pathogenic", "Ref_rec_A", "REFA", _SINGLE, 2),
    ("rs900170001", "1", 1_700_100, "C", "C", "Pathogenic", "Ref_rec_A", "REFA", _SINGLE, 3),
    # Heterozygous user at a site whose "=" record is pathogenic.
    ("rs900170002", "2", 2_700_200, "A", "C", "Pathogenic", "Ref_rec_B", "REFB", _SINGLE, 4),
    ("rs900170002", "2", 2_700_200, "A", "A", "Pathogenic", "Ref_rec_B2", "REFB", _SINGLE, 5),
)


@pytest.fixture
def reference_record_data_dir(tmp_path: Path) -> Path:
    for build in ("GRCh37", "GRCh38"):
        with contextlib.closing(sqlite3.connect(tmp_path / clinvar_db_filename(build))) as conn:
            conn.executescript(CLINVAR_SCHEMA)
            conn.executemany(
                "INSERT INTO clinvar_variants (rsid, chromosome, position, ref, alt, "
                "clinical_significance, condition, gene, review_status, allele_id) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                _REFERENCE_RECORD_ROWS,
            )
            conn.execute(
                "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
                "record_count, remote_signal, local_version_tag) "
                "VALUES (?, 'test://ref-record', '20261010', '2026-10-10', ?, NULL, ?)",
                (
                    clinvar_record_name(build),
                    len(_REFERENCE_RECORD_ROWS),
                    f"iv:{CLINVAR_INTERPRETER_VERSION}",
                ),
            )
            conn.commit()
    return tmp_path


class TestClinVarReferenceAlleleRecords:
    @pytest.mark.parametrize("batch", [False, True])
    def test_hom_ref_chip_call_is_not_annotated(
        self, reference_record_data_dir: Path, cm_stack: contextlib.ExitStack, batch: bool
    ) -> None:
        clinvar = cm_stack.enter_context(ClinVarAnnotator(reference_record_data_dir))
        v = Variant("rs900170001", "1", 1_700_100, "C", "C", build="GRCh38")
        got = list(clinvar.batch_annotate([v])) if batch else clinvar.annotate(v)
        assert got == []

    @pytest.mark.parametrize("batch", [False, True])
    def test_het_call_annotates_only_the_carried_alt(
        self, reference_record_data_dir: Path, cm_stack: contextlib.ExitStack, batch: bool
    ) -> None:
        clinvar = cm_stack.enter_context(ClinVarAnnotator(reference_record_data_dir))
        v = Variant("rs900170002", "2", 2_700_200, "A", "C", build="GRCh38")
        got = list(clinvar.batch_annotate([v])) if batch else clinvar.annotate(v)
        assert [(a.alt, a.condition) for a in got] == [("C", "Ref_rec_B")]
        assert got[0].locus_ref == "A"
        assert got[0].carried_alleles == ("A", "C")


# ── Every annotator supplies the invariant context ─────────────


def test_pipeline_annotations_carry_invariant_context(
    mock_mhg_path: Path, all_annotators_data_dir: Path, cm_stack: contextlib.ExitStack
) -> None:
    """ClinVar, ClinPGx and GWAS rows from a real pipeline run all carry
    ``carried_alleles``, so the zero-carried-copies check is live in
    production rather than only for hand-built rows."""
    annotators = [
        cm_stack.enter_context(ClinVarAnnotator(all_annotators_data_dir)),
        cm_stack.enter_context(PharmGKBAnnotator(all_annotators_data_dir)),
        cm_stack.enter_context(GWASCatalogAnnotator(all_annotators_data_dir)),
    ]
    result = run_analysis(mock_mhg_path, MyHappyGenesParser(), annotators)
    sources = {a.source for a in result.annotations}
    assert {"clinvar", "pharmgkb", "gwas"} <= sources
    for a in result.annotations:
        assert len(a.carried_alleles) == 2, a
        if a.alt:
            assert a.alt in a.carried_alleles, a


def test_snpedia_annotations_carry_invariant_context(
    tmp_path: Path, cm_stack: contextlib.ExitStack
) -> None:
    with contextlib.closing(sqlite3.connect(tmp_path / "snpedia.sqlite")) as conn:
        conn.executescript(
            "CREATE TABLE snpedia_genotypes (rsid TEXT NOT NULL, allele1 TEXT NOT NULL, "
            "allele2 TEXT NOT NULL, magnitude REAL, repute TEXT, summary TEXT, gene TEXT, "
            "orientation TEXT, scraped_at TEXT);"
            "CREATE TABLE database_versions (name TEXT PRIMARY KEY, source_url TEXT NOT NULL, "
            "version TEXT, downloaded_at TEXT NOT NULL, record_count INTEGER NOT NULL, "
            "remote_signal TEXT, local_version_tag TEXT);"
        )
        conn.execute(
            "INSERT INTO snpedia_genotypes VALUES ('rs900170003', 'C', 'T', 3.0, 'Bad', "
            "'synthetic', 'SNPA', 'plus', '2026-10-10T00:00:00')"
        )
        conn.commit()
    snpedia = cm_stack.enter_context(SNPediaAnnotator(tmp_path))
    v = Variant("rs900170003", "1", 1_700_300, "C", "T", build="GRCh38", ref="C")
    got = list(snpedia.batch_annotate([v])) + snpedia.annotate(v)
    assert len(got) == 2
    for a in got:
        assert a.alt == "T"
        assert a.locus_ref == "C"
        assert a.carried_alleles == ("C", "T")


# ── Shared position-resolver rule ──────────────────────────────


class TestResolveRsidAtLocus:
    def test_duplicate_per_scv_rows_are_one_rsid(self) -> None:
        rows = [("T", "G", "rs1"), ("T", "G", "rs1"), ("T", "C", "rs1")]
        assert resolve_rsid_at_locus("T", "T", "G", rows) == ("rs1", ("T", "G"))

    def test_two_alts_of_one_rsid_resolve_without_a_pair(self) -> None:
        rows = [("A", "G", "rs1"), ("A", "T", "rs1")]
        assert resolve_rsid_at_locus("A", "G", "T", rows) == ("rs1", None)

    def test_two_alts_of_different_rsids_abstain(self) -> None:
        rows = [("A", "G", "rs1"), ("A", "T", "rs2")]
        assert resolve_rsid_at_locus("A", "G", "T", rows) is None

    def test_carried_alt_with_two_rsids_abstains(self) -> None:
        rows = [("A", "G", "rs1"), ("A", "G", "rs2")]
        assert resolve_rsid_at_locus("A", "A", "G", rows) is None

    def test_reference_allele_records_never_resolve(self) -> None:
        assert resolve_rsid_at_locus("C", "C", "C", [("C", "C", "rs1")]) is None

    def test_ref_mismatch_abstains(self) -> None:
        assert resolve_rsid_at_locus("G", "G", "A", [("C", "A", "rs1")]) is None

    def test_unknown_ref_picks_the_anchor_that_explains_the_call(self) -> None:
        rows = [("CA", "C", "rs_short"), ("CAG", "C", "rs_long")]
        assert resolve_rsid_at_locus(None, "CAG", "C", rows) == ("rs_long", ("CAG", "C"))

    def test_uncarried_alt_does_not_resolve(self) -> None:
        assert resolve_rsid_at_locus("A", "A", "C", [("A", "G", "rs1")]) is None


# ── ClinVar rows are selected by locus, never by rsID alone ────


_LOCUS_ROWS: tuple[tuple[object, ...], ...] = (
    ("rs900170010", "3", 3_700_100, "G", "A", "Pathogenic", "Locus_A", "LOCA", _SINGLE, 10),
    ("rs900170011", "X", 3_700_200, "T", "C", "Pathogenic", "Locus_X", "LOCX", _SINGLE, 11),
    ("rs900170012", "MT", 3_700, "A", "G", "Pathogenic", "Locus_MT", "LOCMT", _SINGLE, 12),
    # A different variant, recorded by ClinVar at another position.
    ("rs900170013", "3", 3_700_500, "C", "T", "Pathogenic", "Locus_Other", "LOCO", _SINGLE, 13),
)


@pytest.fixture
def locus_clinvar(tmp_path: Path, cm_stack: contextlib.ExitStack) -> ClinVarAnnotator:
    with contextlib.closing(sqlite3.connect(tmp_path / clinvar_db_filename("GRCh38"))) as conn:
        conn.executescript(CLINVAR_SCHEMA)
        conn.executemany(
            "INSERT INTO clinvar_variants (rsid, chromosome, position, ref, alt, "
            "clinical_significance, condition, gene, review_status, allele_id) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            _LOCUS_ROWS,
        )
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
            "record_count, remote_signal, local_version_tag) "
            "VALUES (?, 'test://locus', '20261010', '2026-10-10', ?, NULL, ?)",
            (clinvar_record_name("GRCh38"), len(_LOCUS_ROWS), f"iv:{CLINVAR_INTERPRETER_VERSION}"),
        )
        conn.commit()
    return cm_stack.enter_context(ClinVarAnnotator(tmp_path, builds=("GRCh38",)))


def _both_paths(clinvar: ClinVarAnnotator, v: Variant) -> list[Annotation]:
    single = clinvar.annotate(v)
    assert single == list(clinvar.batch_annotate([v]))
    return single


class TestLocusKeyedClinVarSelection:
    def test_right_rsid_at_wrong_position_is_not_annotated(
        self, locus_clinvar: ClinVarAnnotator
    ) -> None:
        v = Variant("rs900170010", "3", 9_999_999, "G", "A", build="GRCh38")
        assert _both_paths(locus_clinvar, v) == []

    def test_locus_match_annotates_under_the_inputs_rsid(
        self, locus_clinvar: ClinVarAnnotator
    ) -> None:
        """A merged / alias rsID at the same locus is the same variant."""
        v = Variant("rs111111111", "3", 3_700_100, "G", "A", build="GRCh38")
        got = _both_paths(locus_clinvar, v)
        assert [(a.rsid, a.alt) for a in got] == [("rs111111111", "A")]

    def test_rsidless_call_takes_the_rows_rsid(self, locus_clinvar: ClinVarAnnotator) -> None:
        v = Variant("", "3", 3_700_100, "G", "A", build="GRCh38", ref="G")
        assert [a.rsid for a in _both_paths(locus_clinvar, v)] == ["rs900170010"]

    def test_ref_disagreement_is_not_annotated(self, locus_clinvar: ClinVarAnnotator) -> None:
        v = Variant("rs900170010", "3", 3_700_100, "T", "A", build="GRCh38", ref="T")
        assert _both_paths(locus_clinvar, v) == []

    def test_reverse_strand_ref_is_accepted(self, locus_clinvar: ClinVarAnnotator) -> None:
        """REF given as the complement (reverse-stranded input) still identifies the locus."""
        v = Variant("rs900170010", "3", 3_700_100, "C", "T", build="GRCh38", ref="C")
        got = _both_paths(locus_clinvar, v)
        assert [(a.alt, a.carried_alleles) for a in got] == [("A", ("G", "A"))]

    @pytest.mark.parametrize(
        ("chrom", "pos", "a1", "a2"),
        [("chrX", 3_700_200, "T", "C"), ("x", 3_700_200, "C", "C"), ("M", 3_700, "G", "G")],
    )
    def test_chromosome_names_are_normalized(
        self, locus_clinvar: ClinVarAnnotator, chrom: str, pos: int, a1: str, a2: str
    ) -> None:
        v = Variant("rs1", chrom, pos, a1, a2, build="GRCh38")
        assert len(_both_paths(locus_clinvar, v)) == 1


class TestArrayProbeAndForeignLocusGuards:
    """Edge cases for locus-keyed selection.

    A: array vendor probes without an rsID (23andMe ``i5000123``)
       produced pathogenic findings from whatever ClinVar row sat at their
       position, first with no reference allele, then with one inferred
       from the identity table (strand unknown, so a homozygous-reference
       call read on the other strand matched as homozygous pathogenic).
    B: a chip rsID positioned on a different ClinVar variant reported that
       variant's record under the chip rsID (field case: rs369737664 shown
       with the record for rs876660301).
    """

    def test_vendor_probe_without_ref_is_not_locus_matched(
        self, locus_clinvar: ClinVarAnnotator
    ) -> None:
        v = Variant("i5000123", "3", 3_700_100, "G", "A", build="GRCh38")
        assert _both_paths(locus_clinvar, v) == []

    def test_vendor_probe_without_ref_is_not_resolved(
        self, locus_clinvar: ClinVarAnnotator
    ) -> None:
        v = Variant("i5000123", "3", 3_700_100, "A", "A", build="GRCh38")
        assert locus_clinvar.bulk_resolve_rsids([v]) == {}
        assert v.rsid == "i5000123"

    def test_rsidless_call_with_input_ref_is_matched(
        self, locus_clinvar: ClinVarAnnotator
    ) -> None:
        """A no-rsID call with a REF is a VCF record (ID ``.`` or a non-rs ID)
        whose REF came from the file: the REF check applies and it matches.
        The pipeline never gives a vendor probe a REF (see the end-to-end
        test below)."""
        v = Variant("chr3_3700100", "3", 3_700_100, "G", "A", build="GRCh38", ref="G")
        assert [(a.rsid, a.alt) for a in _both_paths(locus_clinvar, v)] == [("rs900170010", "A")]

    def test_chip_rsid_known_elsewhere_does_not_take_the_locus_record(
        self, locus_clinvar: ClinVarAnnotator
    ) -> None:
        """rs900170013 is ClinVar's variant at 3:3700500; at 3:3700100 the row
        belongs to rs900170010 and must not be reported as rs900170013."""
        v = Variant("rs900170013", "3", 3_700_100, "G", "A", build="GRCh38")
        assert _both_paths(locus_clinvar, v) == []

    def test_chip_rsid_at_its_own_clinvar_locus_still_matches(
        self, locus_clinvar: ClinVarAnnotator
    ) -> None:
        v = Variant("rs900170013", "3", 3_700_500, "C", "T", build="GRCh38")
        assert [(a.rsid, a.alt) for a in _both_paths(locus_clinvar, v)] == [("rs900170013", "T")]

    def test_retired_rsid_still_recovers_by_locus(self, locus_clinvar: ClinVarAnnotator) -> None:
        """An rsID ClinVar does not hold (merged / retired, like rs10735810 ->
        rs2228570) still matches the ClinVar variant at its locus."""
        v = Variant("rs222222222", "3", 3_700_100, "G", "A", build="GRCh38")
        got = _both_paths(locus_clinvar, v)
        assert [(a.rsid, a.alt, a.condition) for a in got] == [("rs222222222", "A", "Locus_A")]


_PROBE_CHIP = """\
# MyHappyGenes [TEMPUS]
# Sample ID\tPROBE000001
SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward
i5000123\t3\t3700100\tG\tA
rs900170013\t3\t3700100\tG\tA
i5000124\t3\t3700300\tC\tC
rs900170020\t3\t3700300\tC\tC
"""

# A G>C pathogenic site the identity table covers. A probe reported on the
# minus strand as C/C is homozygous reference (G/G); with a REF of G
# inferred from position it would read as homozygous pathogenic.
_STRAND_ROW = (
    "rs900170020",
    "3",
    3_700_300,
    "G",
    "C",
    "Pathogenic",
    "Locus_S",
    "LOCS",
    _SINGLE,
    20,
)


@pytest.mark.parametrize("with_identity", [False, True])
def test_pipeline_vendor_probe_and_foreign_locus(
    tmp_path: Path, cm_stack: contextlib.ExitStack, with_identity: bool
) -> None:
    """End to end through run_analysis: vendor probes stay silent whether or
    not the identity table covers their position, including a probe whose
    opposite-strand call would match a pathogenic row once given the table's
    REF. The same genotype under its rsID is matched (the rsID confirms the
    variant). The foreign-locus chip rsID is never reported either way."""
    from allelix.annotators.rsid_identity import RsidIdentityAnnotator
    from allelix.databases.rsid_identity_loader import (
        RSID_IDENTITY_DB_FILENAME,
        write_identity_db,
    )

    data = tmp_path / "data"
    data.mkdir()
    with contextlib.closing(sqlite3.connect(data / clinvar_db_filename("GRCh38"))) as conn:
        conn.executescript(CLINVAR_SCHEMA)
        conn.executemany(
            "INSERT INTO clinvar_variants (rsid, chromosome, position, ref, alt, "
            "clinical_significance, condition, gene, review_status, allele_id) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (*_LOCUS_ROWS, _STRAND_ROW),
        )
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
            "record_count, remote_signal, local_version_tag) "
            "VALUES (?, 'test://probe', '1', '2026-10-10', ?, NULL, ?)",
            (
                clinvar_record_name("GRCh38"),
                len(_LOCUS_ROWS) + 1,
                f"iv:{CLINVAR_INTERPRETER_VERSION}",
            ),
        )
        conn.commit()
    write_identity_db(
        [
            ("rs900170010", "3", 3_600_100, "G", "A", 3_700_100, "G", "A", 0),
            ("rs900170020", "3", 3_600_300, "G", "C", 3_700_300, "G", "C", 1),
        ],
        data / RSID_IDENTITY_DB_FILENAME,
        version="dbSNP test",
        source_url="test://probe",
        downloaded_at="2026-10-10T00:00:00Z",
    )
    chip = tmp_path / "probe_chip.txt"
    chip.write_text(_PROBE_CHIP)
    result = run_analysis(
        chip,
        MyHappyGenesParser(),
        [cm_stack.enter_context(ClinVarAnnotator(data, builds=("GRCh38",)))],
        build_override="GRCh38",
        identity=cm_stack.enter_context(RsidIdentityAnnotator(data)) if with_identity else None,
    )
    got = sorted((a.rsid, a.alt, a.zygosity) for a in result.annotations)
    assert got == [("rs900170020", "C", "Homozygous")]
    assert not any(a.rsid.startswith("i") for a in result.annotations)

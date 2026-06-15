# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Mock-data-as-spec invariants (ADR-0015).

The mock data generators define what real source data looks like. Code that
operates on parser output must work against the GENERATED fixtures. Fixtures
that don't represent the real-world format silently hide categorical bugs
(see the v0.4.2 indel-anchor incident, where the MHG mock generator put
"CTT"/"C" at rs113993960 — a multi-base genotype that no real array
produces — and the ClinVar carrier rule appeared correct in tests while
emitting hundreds of false positives in production).

These tests pin invariants the generator output MUST satisfy. If a check
here fails, the generator is wrong and needs fixing (or, rarely, the
invariant itself is wrong and ADR-0015 needs amending).
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

from allelix.parsers.myhappygenes import MyHappyGenesParser

if TYPE_CHECKING:
    from pathlib import Path

_VALID_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
_SINGLE_BASE = re.compile(r"^[ACGT]$")
_NO_CALL = "-"


class TestMhgGeneratorInvariants:
    """The MHG mock generator must produce data shaped like real MHG output.

    Per the MyHappyGenes format spec: tab-delimited, 5 columns,
    `Allele1 - Forward` and `Allele2 - Forward` columns hold single
    nucleotide bases (A/T/G/C) or the no-call marker `-`. Arrays cannot
    report multi-base genotypes — they read one base per probe.
    """

    def test_every_genotype_is_single_base_or_no_call(self, mock_mhg_path: Path):
        """The bug that hid the v0.4.2 indel-anchor incident.

        If this fails: the MHG mock generator is producing genotypes a
        real array can't produce, and any test using that fixture is
        masking real-world behavior. Fix the generator, don't relax this.
        """
        parser = MyHappyGenesParser()
        violations: list[tuple[str, str, str]] = []
        for variant in parser.parse(mock_mhg_path):
            for label, allele in (("allele1", variant.allele1), ("allele2", variant.allele2)):
                if allele == _NO_CALL:
                    continue
                if not _SINGLE_BASE.match(allele):
                    violations.append((variant.rsid, label, allele))
        assert not violations, (
            f"MHG mock generator produced {len(violations)} non-single-base "
            f"alleles. Real MyHappyGenes (Tempus) arrays cannot call indels "
            f"or multi-base alleles — every genotype is a single A/T/G/C or "
            f"a no-call ('-'). Examples: {violations[:5]}. Fix the generator "
            f"at tests/generate_mock_data.py."
        )

    def test_every_chromosome_is_valid(self, mock_mhg_path: Path):
        parser = MyHappyGenesParser()
        invalid = [
            v.chromosome for v in parser.parse(mock_mhg_path) if v.chromosome not in _VALID_CHROMS
        ]
        assert not invalid, f"Invalid chromosomes in mock: {set(invalid)!r}"

    def test_positions_are_positive(self, mock_mhg_path: Path):
        parser = MyHappyGenesParser()
        bad = [(v.rsid, v.position) for v in parser.parse(mock_mhg_path) if v.position <= 0]
        assert not bad, f"Non-positive positions: {bad[:5]}"

    def test_rsids_are_well_formed(self, mock_mhg_path: Path):
        parser = MyHappyGenesParser()
        pattern = re.compile(r"^rs\d+$")
        bad = [v.rsid for v in parser.parse(mock_mhg_path) if not pattern.match(v.rsid)]
        assert not bad, f"Malformed rsIDs: {bad[:5]}"


class TestClinvarMockInvariants:
    """The ClinVar mock fixture must look like real ClinVar VCF output.

    Real ClinVar VCFs contain both SNVs (single-base REF and ALT) and
    indels (multi-base REF or ALT, anchor-base encoded). The mock MUST
    include at least one of each so the indel-anchor protection
    (ADR-0011) has something to fire against in integration tests.
    """

    def test_includes_at_least_one_snv(self, mock_clinvar_vcf: Path):
        from allelix.databases.manager import iter_clinvar_records

        snvs = [
            r
            for r in iter_clinvar_records(mock_clinvar_vcf)
            if len(r["ref"]) == 1 and len(r["alt"]) == 1
        ]
        assert snvs, "Mock ClinVar lacks SNVs — integration tests would be incomplete."

    def test_includes_at_least_one_indel(self, mock_clinvar_vcf: Path):
        """Without an indel row, ADR-0011's protection has nothing to fire on."""
        from allelix.databases.manager import iter_clinvar_records

        indels = [
            r
            for r in iter_clinvar_records(mock_clinvar_vcf)
            if len(r["ref"]) > 1 or len(r["alt"]) > 1
        ]
        assert indels, (
            "Mock ClinVar lacks indel (anchor-base) rows. The indel-anchor "
            "protection (ADR-0011) only has meaningful integration coverage "
            "if at least one ClinVar indel row exists. Add one to "
            "tests/generate_clinvar_fixture.py."
        )

    def test_includes_at_least_one_multi_allelic(self, mock_clinvar_vcf: Path):
        """Without a multi-allelic row, ADR's C-2 split logic has no fixture coverage."""
        with mock_clinvar_vcf.open() as fh:
            multi = [
                line for line in fh if not line.startswith("#") and "," in line.split("\t")[4]
            ]
        assert multi, "Mock ClinVar lacks multi-allelic (comma-ALT) rows."


class TestPharmgkbMockInvariants:
    """The ClinPGx mock fixtures must include rows representative of each
    structured `Allele Function` value the loader's classifier handles.

    Without examples of each function class, ADR-0016's structured
    classifier (function_class column, is_nonfinding flag) can't be
    exercised in integration tests — regressions could slip through.
    """

    def test_includes_at_least_one_carrier_finding(
        self, mock_pharmgkb_dir: Path, mock_cpic_lookup: dict[tuple[str, str], str]
    ):
        from allelix.databases.pharmgkb_loader import iter_pharmgkb_records

        findings = [
            r
            for r in iter_pharmgkb_records(mock_pharmgkb_dir, mock_cpic_lookup)
            if not r["is_nonfinding"]
        ]
        assert findings, (
            "Mock ClinPGx lacks carrier-finding rows. Need at least one "
            "(rsid, genotype) where the user carries a non-Normal allele "
            "per the CPIC lookup."
        )

    def test_includes_at_least_one_nonfinding(
        self, mock_pharmgkb_dir: Path, mock_cpic_lookup: dict[tuple[str, str], str]
    ):
        from allelix.databases.pharmgkb_loader import iter_pharmgkb_records

        nonfindings = [
            r
            for r in iter_pharmgkb_records(mock_pharmgkb_dir, mock_cpic_lookup)
            if r["is_nonfinding"]
        ]
        assert nonfindings, (
            "Mock ClinPGx lacks non-finding rows. ADR-0020 structured "
            "filter has no integration coverage. Add a (rsid, genotype) "
            "where both bases map to Normal function in the CPIC lookup."
        )

    def test_snv_rows_have_empty_allele_function(self, mock_pharmgkb_dir: Path):
        """ADR-0017: real ClinPGx SNV rows have empty `Allele Function`.

        The mock fixture must mirror this — every 2-letter A/C/G/T genotype
        row in `clinical_ann_alleles.tsv` has Allele Function = "". An
        earlier fixture revision violated this and concealed the v0.6.0
        regression (the structured classifier appeared to work in tests
        because the fixture populated Allele Function on every row, the
        inverse of real ClinPGx). Pin the real-data shape here.
        """
        import csv

        snv_re = re.compile(r"^[ACGT]{2}$")
        snv_with_allele_function: list[tuple[str, str, str]] = []
        with (mock_pharmgkb_dir / "clinical_ann_alleles.tsv").open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                genotype = row.get("Genotype/Allele", "").strip()
                allele_function = row.get("Allele Function", "").strip()
                if snv_re.match(genotype) and allele_function:
                    snv_with_allele_function.append(
                        (row.get("Clinical Annotation ID", ""), genotype, allele_function)
                    )
        assert not snv_with_allele_function, (
            "Mock ClinPGx has SNV genotype rows with populated Allele "
            "Function. Real ClinPGx leaves this field empty on every SNV "
            "row (it's populated only on haplotype rows like *1, *2). The "
            "earlier inverted-shape fixture is the proximate cause of the "
            "v0.6.0 production regression — see ADR-0015 + ADR-0017. "
            f"Offending rows: {snv_with_allele_function[:3]}"
        )

    def test_haplotype_rows_may_have_populated_allele_function(self, mock_pharmgkb_dir: Path):
        """The fixture should retain at least one haplotype row with
        populated Allele Function so the structured path is exercised by
        unit tests (even though the loader rejects haplotype genotypes
        at `_normalize_genotype`).
        """
        import csv

        snv_re = re.compile(r"^[ACGT]{2}$")
        populated_haplotypes = []
        with (mock_pharmgkb_dir / "clinical_ann_alleles.tsv").open() as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                genotype = row.get("Genotype/Allele", "").strip()
                allele_function = row.get("Allele Function", "").strip()
                if not snv_re.match(genotype) and allele_function:
                    populated_haplotypes.append((genotype, allele_function))
        assert populated_haplotypes, (
            "Mock ClinPGx lacks haplotype rows with populated Allele "
            "Function. Need at least one (e.g. *1/*2 with Normal function) "
            "to model the real-data shape where structured signal IS "
            "available — just for rows the loader rejects."
        )


class TestGwasMockInvariants:
    """The GWAS Catalog mock fixture must include rows representative of each
    code path the annotator exercises: known risk allele, unknown risk allele,
    multiple p-value tiers, and at least one non-rsID SNP (haplotype skip).

    ADR-0015: without examples of each path in the fixture, integration tests
    can't verify that the annotator handles real-world data correctly.
    """

    def test_includes_single_base_risk_allele(self, mock_gwas_tsv: Path) -> None:
        """Positive carrier-rule path requires a row with a known risk allele."""
        from allelix.databases.gwas_loader import iter_gwas_records

        with_allele = [r for r in iter_gwas_records(mock_gwas_tsv) if r["risk_allele"] is not None]
        assert with_allele, (
            "Mock GWAS fixture lacks rows with a single-base risk allele. "
            "The carrier rule (ADR-0007) can't be exercised without one."
        )

    def test_includes_unknown_risk_allele(self, mock_gwas_tsv: Path) -> None:
        """Unknown-risk path fires on rsID match with magnitude cap (ADR-0024)."""
        from allelix.databases.gwas_loader import iter_gwas_records

        unknown = [r for r in iter_gwas_records(mock_gwas_tsv) if r["risk_allele"] is None]
        assert unknown, (
            "Mock GWAS fixture lacks rows with unknown risk allele (? in "
            "STRONGEST SNP-RISK ALLELE). The unknown-risk-allele magnitude "
            "cap path (ADR-0024) has no fixture coverage."
        )

    def test_spans_pvalue_tiers(self, mock_gwas_tsv: Path) -> None:
        """The fixture must exercise at least three p-value magnitude tiers."""
        from allelix.databases.gwas_loader import iter_gwas_records

        p_values = [
            r["p_value"] for r in iter_gwas_records(mock_gwas_tsv) if r["p_value"] is not None
        ]
        tiers_hit = set()
        for p in p_values:
            if p < 5e-20:
                tiers_hit.add("strong")
            elif p < 5e-8:
                tiers_hit.add("genome_wide")
            elif p < 5e-6:
                tiers_hit.add("suggestive")
            else:
                tiers_hit.add("nominal")
        assert len(tiers_hit) >= 3, (
            f"Mock GWAS fixture only covers {len(tiers_hit)} p-value tier(s): "
            f"{tiers_hit}. Need at least 3 to exercise magnitude scoring. "
            f"See _magnitude() thresholds in allelix/annotators/gwas.py."
        )

    def test_raw_tsv_includes_non_rsid_snp(self, mock_gwas_tsv: Path) -> None:
        """The raw TSV must include at least one row where SNPS is not rs-format.

        iter_gwas_records() skips these (haplotypes, multi-SNP). Without a
        non-rsID row in the fixture, the skip path has no coverage.
        """
        import csv

        rs_re = re.compile(r"^rs\d+$")
        non_rsid_count = 0
        with mock_gwas_tsv.open("r", encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                snp_field = (row.get("SNPS") or "").strip()
                if snp_field and not rs_re.match(snp_field):
                    non_rsid_count += 1
        assert non_rsid_count > 0, (
            "Mock GWAS TSV lacks non-rsID SNP rows (haplotypes, interactions). "
            "The loader's skip path has no fixture coverage. Add a row with "
            "a SNPS value like 'rs123 x rs456' or a haplotype identifier."
        )


class TestClinvarFixturePositionInvariants:
    """ADR-0021 + ADR-0015: each per-build ClinVar fixture must use
    build-correct positions for every rsID present in the build-detect
    table. The fixture itself MUST satisfy the same invariant the
    auto-detector enforces on real input files — otherwise the mock
    has the very bug it's supposed to test against.
    """

    def _read_vcf_positions(self, path: Path) -> dict[str, tuple[str, int]]:
        import re

        rs_re = re.compile(r"RS=(\d+)")
        positions: dict[str, tuple[str, int]] = {}
        for line in path.read_text().splitlines():
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 8:
                continue
            chrom, pos = cols[0], int(cols[1])
            info = cols[7]
            match = rs_re.search(info)
            if match:
                positions[f"rs{match.group(1)}"] = (chrom, pos)
        return positions

    def test_grch37_fixture_positions_match_table(self, mock_clinvar_grch37_vcf: Path):
        from allelix.utils.build_detect import BUILD_GRCH37, KNOWN_SNP_POSITIONS

        actual = self._read_vcf_positions(mock_clinvar_grch37_vcf)
        mismatches: list[str] = []
        for rsid, (chrom, pos) in actual.items():
            if rsid not in KNOWN_SNP_POSITIONS:
                continue  # not a detector-table SNP; positions are free
            expected_chrom, expected_pos = KNOWN_SNP_POSITIONS[rsid][BUILD_GRCH37]
            if (chrom, pos) != (expected_chrom, expected_pos):
                mismatches.append(
                    f"{rsid}: fixture has {chrom}:{pos}, build_detect "
                    f"says {expected_chrom}:{expected_pos}"
                )
        assert not mismatches, (
            "Mock GRCh37 ClinVar fixture has positions that don't match the "
            "build_detect table. This is exactly the bug ADR-0021 was written "
            "to catch in production — hardcoded into our test fixture.\n  "
            + "\n  ".join(mismatches)
        )

    def test_grch38_fixture_positions_match_table(self, mock_clinvar_grch38_vcf: Path):
        from allelix.utils.build_detect import BUILD_GRCH38, KNOWN_SNP_POSITIONS

        actual = self._read_vcf_positions(mock_clinvar_grch38_vcf)
        mismatches: list[str] = []
        for rsid, (chrom, pos) in actual.items():
            if rsid not in KNOWN_SNP_POSITIONS:
                continue
            expected_chrom, expected_pos = KNOWN_SNP_POSITIONS[rsid][BUILD_GRCH38]
            if (chrom, pos) != (expected_chrom, expected_pos):
                mismatches.append(
                    f"{rsid}: fixture has {chrom}:{pos}, build_detect "
                    f"says {expected_chrom}:{expected_pos}"
                )
        assert not mismatches, (
            "Mock GRCh38 ClinVar fixture has positions that don't match the "
            "build_detect table.\n  " + "\n  ".join(mismatches)
        )

    def test_nipa1_has_strand_inverted_ref_alt_across_builds(
        self, mock_clinvar_grch37_vcf: Path, mock_clinvar_grch38_vcf: Path
    ):
        """ADR-0021 motivation pinned: rs104894490's REF/ALT inverts
        between builds. The fixture MUST preserve this — otherwise the
        per-build dispatch test below is meaningless.
        """
        import re

        rs_re = re.compile(r"RS=104894490")

        def find_ref_alt(path: Path) -> tuple[str, str] | None:
            for line in path.read_text().splitlines():
                if not line or line.startswith("#"):
                    continue
                cols = line.split("\t")
                if len(cols) < 8 or not rs_re.search(cols[7]):
                    continue
                return cols[3], cols[4]
            return None

        g37 = find_ref_alt(mock_clinvar_grch37_vcf)
        g38 = find_ref_alt(mock_clinvar_grch38_vcf)
        assert g37 == ("C", "G"), f"GRCh37 NIPA1 expected REF=C ALT=G, got {g37}"
        assert g38 == ("G", "A"), f"GRCh38 NIPA1 expected REF=G ALT=A, got {g38}"


class TestCpicLookupMockInvariants:
    """ADR-0015 + m-6: MOCK_CPIC_LOOKUP must mirror real CPIC's shape.

    A fixture that's structurally different from real CPIC output
    silently hides categorical bugs (see the v0.6.0 inverted-shape
    incident). Pin the key invariants here so a future refactor of
    either the mock or the loader can't drift apart without a test
    failure.
    """

    def test_keys_are_rsid_and_single_base_tuples(
        self, mock_cpic_lookup: dict[tuple[str, str], str]
    ):
        rsid_re = re.compile(r"^rs\d+$")
        for key in mock_cpic_lookup:
            assert isinstance(key, tuple) and len(key) == 2, (
                f"MOCK_CPIC_LOOKUP key must be (rsid, base) tuple, got {key!r}. "
                "Real CPIC fetch returns the same shape — see "
                "fetch_cpic_allele_functions()."
            )
            rsid, base = key
            assert rsid_re.match(rsid), (
                f"MOCK_CPIC_LOOKUP rsid {rsid!r} doesn't match real CPIC's "
                f"sequence_location.dbsnpid format (rs followed by digits)."
            )
            assert base in {"A", "C", "G", "T"}, (
                f"MOCK_CPIC_LOOKUP base {base!r} for {rsid} isn't a single "
                f"A/C/G/T nucleotide. fetch_cpic_allele_functions() filters "
                f"multi-base alleles — the mock must mirror that."
            )

    def test_values_are_known_function_classes(self, mock_cpic_lookup: dict[tuple[str, str], str]):
        from allelix.databases.cpic_loader import (
            FUNCTION_CLASS_DECREASED,
            FUNCTION_CLASS_INCREASED,
            FUNCTION_CLASS_NO_FUNCTION,
            FUNCTION_CLASS_NORMAL,
            FUNCTION_CLASS_UNCERTAIN,
        )

        valid = {
            FUNCTION_CLASS_NORMAL,
            FUNCTION_CLASS_DECREASED,
            FUNCTION_CLASS_NO_FUNCTION,
            FUNCTION_CLASS_INCREASED,
            FUNCTION_CLASS_UNCERTAIN,
        }
        for key, value in mock_cpic_lookup.items():
            assert value in valid, (
                f"MOCK_CPIC_LOOKUP[{key!r}] = {value!r} isn't a recognized "
                f"function_class. Real CPIC's clinicalfunctionalstatus maps "
                f"to one of {sorted(valid)}."
            )

    def test_includes_both_normal_and_non_normal(
        self, mock_cpic_lookup: dict[tuple[str, str], str]
    ):
        """The fixture must exercise both branches of the filter — pure
        Normal homozygotes (suppression) and at least one non-Normal
        allele (carrier emission). A mock missing either branch would
        give false confidence in the filter's correctness.
        """
        from allelix.databases.cpic_loader import FUNCTION_CLASS_NORMAL

        values = set(mock_cpic_lookup.values())
        assert FUNCTION_CLASS_NORMAL in values, (
            "MOCK_CPIC_LOOKUP has no Normal-function alleles. Non-finding "
            "suppression can't be exercised."
        )
        non_normal = values - {FUNCTION_CLASS_NORMAL}
        assert non_normal, (
            "MOCK_CPIC_LOOKUP has no non-Normal alleles. Carrier emission can't be exercised."
        )

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for core data models."""

from __future__ import annotations

from allelix.models import Annotation, Variant


class TestVariant:
    def test_homozygous(self):
        v = Variant("rs1", "1", 100, "A", "A")
        assert not v.is_heterozygous
        assert not v.is_no_call
        assert v.genotype == "A/A"

    def test_heterozygous(self):
        v = Variant("rs1", "1", 100, "A", "G")
        assert v.is_heterozygous
        assert not v.is_no_call
        assert v.genotype == "A/G"

    def test_no_call_both_alleles(self):
        v = Variant("rs1", "1", 100, "-", "-")
        assert v.is_no_call
        assert not v.is_heterozygous

    def test_no_call_one_allele(self):
        v = Variant("rs1", "1", 100, "A", "-")
        assert v.is_no_call
        assert not v.is_heterozygous

    def test_default_build_is_grch37(self):
        v = Variant("rs1", "1", 100, "A", "A")
        assert v.build == "GRCh37"

    def test_explicit_build(self):
        v = Variant("rs1", "1", 100, "A", "A", build="GRCh38")
        assert v.build == "GRCh38"

    def test_indel_genotype(self):
        v = Variant("rs1", "7", 117199644, "CTT", "C")
        assert v.is_heterozygous
        assert v.genotype == "CTT/C"

    def test_lowercase_alleles_normalized_to_uppercase(self):
        """GH #14 (case half): a lowercase user allele used to silently
        fail the raw-set-membership carrier match in clinical annotators
        (which check `alt in {variant.allele1, variant.allele2}` against
        uppercase REF/ALT from the source databases). Normalize at the
        ``Variant`` model boundary so the invariant is impossible to
        violate downstream — even if a custom panel / filter file leaks
        lowercase through.
        """
        v = Variant("rs1801133", "1", 11796321, "a", "g")
        assert v.allele1 == "A"
        assert v.allele2 == "G"

    def test_mixed_case_alleles_normalized(self):
        v = Variant("rs1", "1", 100, "Aa", "gT")  # synthetic mixed-case indel
        assert v.allele1 == "AA"
        assert v.allele2 == "GT"

    def test_no_call_marker_preserved(self):
        """``-`` must not be uppercased (it would still equal itself,
        but be explicit about the carve-out)."""
        v = Variant("rs1", "1", 100, "-", "-")
        assert v.allele1 == "-"
        assert v.allele2 == "-"
        assert v.is_no_call

    def test_empty_allele_preserved(self):
        """Empty string passes through unchanged (defensive — should not
        happen, but ``"".upper()`` would still be ``""``)."""
        v = Variant("rs1", "1", 100, "", "A")
        assert v.allele1 == ""
        assert v.allele2 == "A"


class TestAnnotation:
    def _minimal(self, **overrides):
        defaults = {
            "source": "clinvar",
            "rsid": "rs1",
            "significance": "clinvar_pathogenic",
            "category": "clinical",
            "magnitude": 5.0,
            "description": "Test annotation",
            "attribution": "ClinVar",
            "genotype_match": "A/A",
        }
        defaults.update(overrides)
        return Annotation(**defaults)

    def test_required_fields(self):
        a = self._minimal()
        assert a.source == "clinvar"
        assert a.attribution == "ClinVar"

    def test_default_optional_fields(self):
        a = self._minimal()
        assert a.references == []
        assert a.condition == ""
        assert a.gene == ""

    def test_default_references_independent_per_instance(self):
        a1 = self._minimal()
        a2 = self._minimal()
        a1.references.append("PMID:1")
        assert a2.references == []

    def test_optional_fields_set(self):
        a = self._minimal(references=["PMID:1234"], condition="Sickle cell", gene="HBB")
        assert a.references == ["PMID:1234"]
        assert a.condition == "Sickle cell"
        assert a.gene == "HBB"

    def test_zygosity_homozygous(self):
        a = self._minimal(genotype_match="A/A")
        assert a.zygosity == "Homozygous"

    def test_zygosity_heterozygous(self):
        a = self._minimal(genotype_match="A/G")
        assert a.zygosity == "Heterozygous"

    def test_zygosity_no_call(self):
        a = self._minimal(genotype_match="A/-")
        assert a.zygosity == "No Call"

    def test_zygosity_no_call_both(self):
        a = self._minimal(genotype_match="-/-")
        assert a.zygosity == "No Call"

    def test_zygosity_single_allele(self):
        a = self._minimal(genotype_match="A")
        assert a.zygosity == "Homozygous"

    def test_zygosity_single_allele_heterozygous(self):
        a = self._minimal(genotype_match="AG")
        assert a.zygosity == "Heterozygous"

    def test_zygosity_concat_no_call(self):
        a = self._minimal(genotype_match="A-")
        assert a.zygosity == "No Call"

    def test_zygosity_concat_no_call_both(self):
        a = self._minimal(genotype_match="--")
        assert a.zygosity == "No Call"

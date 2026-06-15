# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Core data models for genotype variants and reference annotations.

Trust boundary: parsers are responsible for validating raw input. Model
constructors do not enforce chromosome names, position bounds, or allele
encodings — they trust their caller. If a Variant or Annotation is
constructed by code outside the `allelix.parsers` package, the caller owns
the validation.
"""

from __future__ import annotations

from dataclasses import dataclass, field

NO_CALL_MARKER = "-"
DEFAULT_BUILD = "GRCh37"


@dataclass
class Variant:
    """A single genotype call: which alleles a sample carries at one position.

    All parsers normalize to this representation. Downstream code (annotators,
    reports) only sees Variants, never raw file formats.

    Attributes:
        rsid: dbSNP reference identifier (e.g., "rs1801133").
        chromosome: Chromosome name. "1"-"22", "X", "Y", or "MT".
        position: 1-based genomic coordinate in the given build.
        allele1: First observed allele. A/T/G/C, multi-base for indels, or "-" for no-call.
        allele2: Second observed allele. Same encoding as allele1.
        build: Reference genome build. "GRCh37" (hg19) or "GRCh38" (hg38).
    """

    rsid: str
    chromosome: str
    position: int
    allele1: str
    allele2: str
    build: str = DEFAULT_BUILD

    @property
    def is_heterozygous(self) -> bool:
        """True if the two alleles differ (and neither is a no-call)."""
        if self.is_no_call:
            return False
        return self.allele1 != self.allele2

    @property
    def is_no_call(self) -> bool:
        """True if either allele is the no-call marker.

        Typically indicates assay failure at this position, but the precise
        meaning is format-dependent (some VCFs use `-` for indel deletions).
        """
        return self.allele1 == NO_CALL_MARKER or self.allele2 == NO_CALL_MARKER

    @property
    def genotype(self) -> str:
        """Human-readable genotype string (e.g., "C/T")."""
        return f"{self.allele1}/{self.allele2}"


@dataclass
class Annotation:
    """A claim about a variant sourced from a specific reference database.

    Allelix never asserts variant significance directly — every Annotation is
    attributed to its source database. See README § Regulatory Posture.

    Attributes:
        source: Lowercase database identifier (e.g., "clinvar", "pharmgkb").
        rsid: The variant this annotation applies to.
        significance: Source-prefixed classification (e.g., "clinvar_pathogenic").
        category: Coarse filter bucket. Use non-diagnostic labels: "clinical",
            "pharma", "carrier", "trait", "methylation". Never bare medical terms
            like "pathogenic" — those would read as Allelix's own classification.
        magnitude: 0-10 importance score (SNPedia-style).
        description: Human-readable explanation.
        attribution: Display name of the source ("ClinVar", "ClinPGx", ...).
        genotype_match: Which genotype triggers this annotation. For SNVs this
            is a concatenated, sorted allele pair (e.g., "AG", "TT"); the slash
            form (e.g., "AT/A") appears only for indels.
        references: PubMed IDs or URLs supporting the claim.
        condition: Disease or condition name, if applicable.
        gene: Gene symbol, if known.
        review_status: ClinVar review status (CLNREVSTAT), empty for non-ClinVar.
        is_must_include: Internal flag for GWAS rollup; excluded from public output.
    """

    source: str
    rsid: str
    significance: str
    category: str
    magnitude: float
    description: str
    attribution: str
    genotype_match: str
    references: list[str] = field(default_factory=list)
    condition: str = ""
    gene: str = ""
    review_status: str = ""
    alt: str = ""
    is_must_include: bool = False
    allele_frequency: float | None = None
    am_pathogenicity: float | None = None
    am_class: str = ""
    cadd_phred: float | None = None

    @property
    def zygosity(self) -> str:
        """Classify the genotype call as Heterozygous, Homozygous, or No Call."""
        if NO_CALL_MARKER in self.genotype_match:
            return "No Call"
        parts = self.genotype_match.split("/")
        if len(parts) != 2:
            return "Homozygous" if len(set(self.genotype_match)) == 1 else "Heterozygous"
        return "Heterozygous" if parts[0] != parts[1] else "Homozygous"

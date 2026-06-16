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

from dataclasses import asdict, dataclass, field

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
        ref: Reference allele on the forward strand at this position, when
            known. VCF parsers populate from the REF column; array parsers
            currently leave it None. Required by strand-aware carrier matching
            (ADR-0035) and downstream consumers that need to identify the risk
            allele within the user's pair. None means "cannot resolve" —
            consumers degrade gracefully (return ambiguous, skip strand check,
            etc.) rather than guess.
    """

    rsid: str
    chromosome: str
    position: int
    allele1: str
    allele2: str
    build: str = DEFAULT_BUILD
    ref: str | None = None

    def __post_init__(self) -> None:
        """Normalize allele case at construction (GH #14, ADR-0035).

        Reference databases (ClinVar, gnomAD, ClinPGx, etc.) all ship
        uppercase alleles, and carrier matching is raw set membership
        against ``{allele1, allele2}`` — a lowercase user allele would
        silently fail to match and zero annotations would be produced
        for a real carrier. Production parsers all emit uppercase
        today, but a user-supplied filter file (custom panel) or a
        future format variant could leak lowercase through. Normalize
        at the model boundary so the invariant is impossible to
        violate downstream. The no-call marker is left as-is;
        multi-base alleles (indels) are uppercased in place.

        ``ref`` (ADR-0035) is a sibling allele field feeding the same
        downstream matching path (PR 4 strand-aware carrier matching
        compares ``ref`` against ``{allele1, allele2}`` and against
        gnomAD's uppercased REF). VCFs derived from soft-masked
        references emit lowercase bases in the REF column, so the
        same normalization applies here for the same reason.
        """
        if self.allele1 and self.allele1 != NO_CALL_MARKER:
            self.allele1 = self.allele1.upper()
        if self.allele2 and self.allele2 != NO_CALL_MARKER:
            self.allele2 = self.allele2.upper()
        if self.ref is not None and self.ref != NO_CALL_MARKER:
            self.ref = self.ref.upper()

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
    # Internal flag — see ``_INTERNAL_ANNOTATION_FIELDS`` below. Excluded
    # from every public serialization path (JSON report, diff entries,
    # diff change records). If you add another internal-only field,
    # extend the frozenset; the three call sites all read from it.
    is_must_include: bool = False
    allele_frequency: float | None = None
    am_pathogenicity: float | None = None
    am_class: str = ""
    cadd_phred: float | None = None
    # ADR-0035 cluster manifest (PR 3): structured GWAS fields promoted out of
    # the rendered ``description`` prose. `description` still carries the same
    # rendered text for HTML / terminal display, but downstream code (rollup,
    # diff, future consumers) reads the structured fields directly instead of
    # regex-parsing prose. Populated by the GWAS annotator; empty / None on
    # non-GWAS rows.
    trait: str = ""
    p_value: float | None = None
    phecode: str = ""

    @property
    def zygosity(self) -> str:
        """Classify the genotype call as Heterozygous, Homozygous, or No Call."""
        if NO_CALL_MARKER in self.genotype_match:
            return "No Call"
        parts = self.genotype_match.split("/")
        if len(parts) != 2:
            return "Homozygous" if len(set(self.genotype_match)) == 1 else "Heterozygous"
        return "Heterozygous" if parts[0] != parts[1] else "Homozygous"


# Annotation fields that are internal-only and MUST be stripped from every
# public serialization path (JSON report annotations, JSON diff entries,
# diff change records). Adding a new internal field? Extend this frozenset
# and every existing serializer picks up the exclusion automatically — no
# new "remember to strip the field" obligation across N call sites.
#
# GH allelix-dev #3 (the v2.1.1 "is_must_include leak by omission risk"
# cleanup) eliminated the three duplicated inline filters that previously
# carried this rule. If you find yourself writing
# ``{k: v for k, v in asdict(a).items() if k != "<some_field>"}`` in a
# new public-output path, you're re-introducing the bug — add the field
# here and use ``annotation_to_public_dict`` instead.
_INTERNAL_ANNOTATION_FIELDS: frozenset[str] = frozenset({"is_must_include"})


def annotation_to_public_dict(a: Annotation) -> dict:
    """Serialize an Annotation for public output, stripping internal fields.

    Canonical helper for every public Annotation serialization path. The
    JSON report's main annotation list, JSON diff entries, and diff change
    records all route through this — so any future internal-only field
    just needs adding to ``_INTERNAL_ANNOTATION_FIELDS`` and every consumer
    inherits the exclusion.

    Does NOT add display-derived fields (``zygosity``, ``am_caveat``,
    diff-specific ``previous_*``) — those are layered by the calling
    serializer after this helper has stripped internal state.
    """
    return {k: v for k, v in asdict(a).items() if k not in _INTERNAL_ANNOTATION_FIELDS}

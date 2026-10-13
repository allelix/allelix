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
    # pre-resolution rsID stashed when a resolver
    # (`ClinVarAnnotator.bulk_resolve_rsids` or
    # `GnomadAnnotator.bulk_resolve_rsids_from_positions`, both invoked
    # in ``_pipeline._flush()``) recovers the real ``rs...`` ID for a
    # variant that arrived without one — e.g. a DeepVariant gVCF row
    # whose ID column is ``.``, parsed as a positional pseudo-ID like
    # ``chr1:11856378:G:A``. The stamp rewrites ``rsid`` to the
    # recovered value so every downstream call site that compares
    # ``v.rsid`` against a requested set (panel coverage, high-value
    # match, extract, PLINK ref/alt, GWAS must-include) works
    # transparently; ``original_rsid`` holds the pre-stamp value
    # purely for audit / debug.
    #
    # **DEBUG-ONLY — do not branch production logic on this field.**
    # Reading ``original_rsid`` outside of logging, debugging, or
    # ``__repr__``-style introspection couples production behavior to
    # provenance metadata and turns this into an unintended API
    # surface. The serialization-guard test in
    # ``tests/test_original_rsid_never_serializes.py`` asserts this
    # field never reaches JSON / HTML / ``--diff`` output.
    original_rsid: str | None = None

    def __post_init__(self) -> None:
        """Normalize allele case at construction (ADR-0035).

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
        downstream matching path (strand-aware carrier matching
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


class AnnotationInvariantError(ValueError):
    """An Annotation was constructed in violation of an allele-identity invariant."""


# Sources whose ``Annotation.alt`` names a non-reference allele the user
# carries. GWAS is deliberately absent: its ``alt`` is the effect (risk)
# allele reported by the study, which may be the reference allele.
ALT_IS_NON_REFERENCE_SOURCES: frozenset[str] = frozenset({"clinvar", "pharmgkb", "snpedia"})


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
        alt: The allele the annotation is about. For ClinVar, ClinPGx and
            SNPedia this is a non-reference allele the user carries. For the
            GWAS Catalog it is the study's effect (risk) allele, which may be
            the reference allele; empty when the study does not report one.
        is_must_include: Internal flag for GWAS rollup; excluded from public output.
        plausibility: Quality-control flag from the call's zygosity and the
            allele's gnomAD frequency (ADR-0041): "plausible", "suspicious",
            "implausible", or empty when not assessed. Never suppresses the
            annotation.
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
    plausibility: str = ""
    # ADR-0035 cluster manifest: structured GWAS fields promoted out of
    # the rendered ``description`` prose. `description` still carries the same
    # rendered text for HTML / terminal display, but downstream code (rollup,
    # diff, future consumers) reads the structured fields directly instead of
    # regex-parsing prose. Populated by the GWAS annotator; empty / None on
    # non-GWAS rows.
    trait: str = ""
    p_value: float | None = None
    phecode: str = ""
    # Internal allele-identity context for the construction invariants
    # below. Excluded from every public serialization via
    # ``_INTERNAL_ANNOTATION_FIELDS`` and from equality / repr, so two
    # Annotations with the same public content still compare equal.
    #
    # ``locus_ref``: forward-strand reference allele at the locus, when the
    #   annotator knows it (the ClinVar row's REF, or ``Variant.ref``).
    # ``carried_alleles``: the user's two alleles in the same orientation
    #   as ``alt`` (complemented when a strand-flipped match fired).
    locus_ref: str | None = field(default=None, compare=False, repr=False)
    carried_alleles: tuple[str, ...] = field(default=(), compare=False, repr=False)
    # ``locus``: the input variant's ``(chromosome, position)`` on the
    #   analysis build, chromosome in Allelix's bare form. Lets enrichment
    #   key every annotation by position rather than by rsID (ADR-0040).
    locus: tuple[str, int] | None = field(default=None, compare=False, repr=False)
    # ``reference_carries_alt``: ``alt`` is the reference base because this
    #   build's reference sequence carries the variant allele (F5 c.1601G>A
    #   on GRCh37): a ClinVar row written REF == ALT, or a ClinPGx / SNPedia
    #   genotype whose alt is derived against the normal allele from the
    #   other build. Exempts the annotation from invariant 2 only.
    reference_carries_alt: bool = field(default=False, compare=False, repr=False)

    def __post_init__(self) -> None:
        """Enforce the allele-identity invariants.

        1. An Annotation whose ``alt`` is set never names an allele the user
           carries zero copies of. Applies to every source. Checked when
           the annotator supplied ``carried_alleles``.
        2. For sources whose ``alt`` means a non-reference allele
           (``ALT_IS_NON_REFERENCE_SOURCES``), ``alt`` is never the
           reference base. GWAS is exempt: its ``alt`` is the study's
           effect (risk) allele, which is often the reference allele.
           Checked when ``locus_ref`` is known. An annotation whose alt is
           the reference base because the build's reference carries the
           variant allele is exempt (``reference_carries_alt``).

        Annotators filter offending source rows before construction; a
        violation here is a programming error, not a data condition.
        """
        if self.alt and self.carried_alleles and self.alt not in self.carried_alleles:
            raise AnnotationInvariantError(
                f"{self.source} annotation for {self.rsid}: alt {self.alt!r} is not "
                f"carried (user alleles {self.carried_alleles!r})"
            )
        if (
            self.source in ALT_IS_NON_REFERENCE_SOURCES
            and self.alt
            and self.locus_ref is not None
            and self.alt == self.locus_ref
            and not self.reference_carries_alt
        ):
            raise AnnotationInvariantError(
                f"{self.source} annotation for {self.rsid}: alt {self.alt!r} is the "
                "reference allele"
            )

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
# (the v2.1.1 "is_must_include leak by omission risk"
# cleanup) eliminated the three duplicated inline filters that previously
# carried this rule. If you find yourself writing
# ``{k: v for k, v in asdict(a).items() if k != "<some_field>"}`` in a
# new public-output path, you're re-introducing the bug — add the field
# here and use ``annotation_to_public_dict`` instead.
_INTERNAL_ANNOTATION_FIELDS: frozenset[str] = frozenset(
    {"is_must_include", "locus_ref", "carried_alleles", "locus", "reference_carries_alt"}
)


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

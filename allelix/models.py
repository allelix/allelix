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
# Array codes for a deletion and an insertion allele (23andMe ``DI``), not bases.
DI_ALLELES = frozenset({"D", "I"})
DEFAULT_BUILD = "GRCh37"
# A half call's missing allele as reports show it (``T/.``), and its zygosity.
MISSING_ALLELE_DISPLAY = "."
HALF_CALL = "Half call"

# Zygosity: a machine code in JSON, a display string in HTML and the
# terminal. ``hemizygous`` is a call the file declares haploid (a one-index
# VCF GT, a one-letter array call): one copy of the allele.
HOMOZYGOUS = "homozygous"
HETEROZYGOUS = "heterozygous"
HEMIZYGOUS = "hemizygous"
HALF_CALL_CODE = "half_call"
NO_CALL_CODE = "no_call"
ZYGOSITY_DISPLAY: dict[str, str] = {
    HOMOZYGOUS: "Homozygous",
    HETEROZYGOUS: "Heterozygous",
    HEMIZYGOUS: "Hemizygous",
    HALF_CALL_CODE: HALF_CALL,
    NO_CALL_CODE: "No call",
}


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
        ref_from_input: True when ``ref`` came from the input file (a VCF
            REF column), False when it is unknown or was filled in from a
            reference table. ClinVar rows without an rsID match only calls
            whose REF came from the input (ADR-0046).
    """

    rsid: str
    chromosome: str
    position: int
    allele1: str
    allele2: str
    build: str = DEFAULT_BUILD
    ref: str | None = None
    ref_from_input: bool = False
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
    # The file declares one copy: a one-index VCF GT, or a one-letter array
    # call. The parsers still give both alleles (the called one twice), so
    # matching is unchanged; reports show the call as hemizygous.
    haploid: bool = False

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
    def is_half_call(self) -> bool:
        """True for a VCF call with one allele called and the other missing.

        GT ``./1`` or ``1/.`` (ADR-0047). ``is_no_call`` is also true for
        it: sources that need both alleles (ClinPGx, SNPedia) skip it, and
        those that need only a carried allele (ClinVar, GWAS) assess the
        called one. An array call never counts: its REF does not come from
        the input.
        """
        return self.ref_from_input and (self.allele1 == NO_CALL_MARKER) != (
            self.allele2 == NO_CALL_MARKER
        )

    @property
    def carried_pair(self) -> tuple[str, str]:
        """The two alleles for allele-membership checks.

        A half call's called allele stands for both, so a check of "which
        alleles does the user carry" sees only real alleles. Never use it
        for zygosity or copy counts.
        """
        if self.is_half_call:
            called = self.allele2 if self.allele1 == NO_CALL_MARKER else self.allele1
            return called, called
        return self.allele1, self.allele2

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
            allele's gnomAD frequency (ADR-0041, ADR-0050): "plausible",
            "suspicious", "implausible" or "not_assessed". Never suppresses
            the annotation.
        plausibility_reason: Why plausibility was not assessed (ADR-0050):
            "no_frequency", "gnomad_not_run", "chromosome", "half_call",
            "no_allele"; empty when it was assessed.
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
    plausibility_reason: str = ""
    # ADR-0035 cluster manifest: structured GWAS fields promoted out of
    # the rendered ``description`` prose. `description` still carries the same
    # rendered text for HTML / terminal display, but downstream code (rollup,
    # diff, future consumers) reads the structured fields directly instead of
    # regex-parsing prose. Populated by the GWAS annotator; empty / None on
    # non-GWAS rows.
    trait: str = ""
    p_value: float | None = None
    phecode: str = ""
    # GWAS rows only: what the study's effect allele does,
    # ``increases``, ``decreases`` or ``unknown``; empty on other sources.
    effect_direction: str = ""
    # ADR-0045, ClinVar rows only (empty / None / False elsewhere):
    # ``clinvar_aggregate``: ClinVar's own classification of the variant.
    # ``submitters_same`` of ``submitters_total``: how many of the
    #   submitters it counts gave this row's classification.
    # ``aggregate_capped``: the magnitude was lowered to the aggregate's.
    # ``outside_aggregate``: an outranked P/LP submission, reported because
    #   the aggregate is not a disease classification (drug response, ...).
    # ``common_allele_capped`` (ADR-0052): a zero-star P/LP submission on an
    #   allele whose gnomAD frequency is at least 5%, scored below the
    #   default floor.
    # ``clinvar_allele_af`` (ADR-0052): gnomAD's frequency of ClinVar's own
    #   allele, set on a zero-star P/LP row when it could be looked up.
    clinvar_aggregate: str = ""
    submitters_same: int | None = None
    submitters_total: int | None = None
    aggregate_capped: bool = False
    outside_aggregate: bool = False
    common_allele_capped: bool = False
    clinvar_allele_af: float | None = None
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
    # ``record_key``: ``(chromosome, position, REF, ALT)`` of the input VCF
    #   record and allele as written, when the annotation matched only after
    #   trimming that allele (ADR-0047), so ``locus`` and ``alt`` are the
    #   trimmed form. ``--vcf-out`` stamps the record by it.
    record_key: tuple[str, int, str, str] | None = field(default=None, compare=False, repr=False)
    # ``haploid``: the call is one the file declares haploid (Variant.haploid);
    #   ``genotype_match`` then holds the one allele.
    haploid: bool = field(default=False, compare=False, repr=False)
    # ``effect_measure``: ``odds_ratio`` or ``beta`` when the GWAS row's
    #   effect direction was read from one, for the report's plain words.
    effect_measure: str = field(default="", compare=False, repr=False)

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
    def zygosity_code(self) -> str:
        """The call's zygosity as a JSON code (homozygous, hemizygous, half_call, ...)."""
        if NO_CALL_MARKER in self.genotype_match:
            return NO_CALL_CODE
        if self.haploid:
            return HEMIZYGOUS
        parts = self.genotype_match.split("/")
        if MISSING_ALLELE_DISPLAY in parts:
            return HALF_CALL_CODE
        if len(parts) != 2:
            return HOMOZYGOUS if len(set(self.genotype_match)) == 1 else HETEROZYGOUS
        return HETEROZYGOUS if parts[0] != parts[1] else HOMOZYGOUS

    @property
    def zygosity(self) -> str:
        """The call's zygosity for display (Homozygous, Hemizygous, Half call, ...)."""
        return ZYGOSITY_DISPLAY[self.zygosity_code]


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
    {
        "is_must_include",
        "locus_ref",
        "carried_alleles",
        "locus",
        "reference_carries_alt",
        "record_key",
        "haploid",
        "effect_measure",
    }
)


def report_rank_key(annotation: Annotation) -> tuple[float, str]:
    """The order reports list rows in: highest magnitude first, then rsID.

    Sorts using it are stable, so rows that tie keep the pipeline's order.
    """
    return (-annotation.magnitude, annotation.rsid)


def half_call_genotype(variant: Variant) -> str | None:
    """A half call's genotype as reports show it (``T/.``), or None for any other call."""
    if not variant.is_half_call:
        return None
    return f"{variant.carried_pair[0]}/{MISSING_ALLELE_DISPLAY}"


def diploid_label(variant: Variant) -> str:
    """The user's call as reports show it (ADR-0023), the same for every source.

    Two one-base alleles of an SNV are sorted and joined (``("G", "A")``
    gives ``AG``). An indel's alleles keep a slash (``AC/A``), including a
    homozygous one-base deletion from a VCF record whose REF is longer
    (``AC>A`` 1/1 is ``A/A``, not an SNV-like ``AA``). A half call is
    ``T/.`` (ADR-0047).
    """
    half = half_call_genotype(variant)
    if half is not None:
        return half
    a1, a2 = variant.allele1, variant.allele2
    record_is_indel = variant.ref_from_input and len(variant.ref or "") > 1
    if len(a1) == 1 and len(a2) == 1 and not record_is_indel:
        return "".join(sorted((a1, a2)))
    return f"{a1}/{a2}"


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
    d = {k: v for k, v in asdict(a).items() if k not in _INTERNAL_ANNOTATION_FIELDS}
    # ADR-0046: where the call sits, so a finding without an rsID still
    # names its variant. ``ref`` is the forward-strand reference allele
    # when the annotator knows it.
    d["chromosome"], d["position"] = a.locus if a.locus is not None else ("", None)
    d["ref"] = a.locus_ref or ""
    return d


def variant_label(rsid: str, chromosome: str, position: int | None, ref: str, alt: str) -> str:
    """The rsID, or ``chrom:pos ref>alt`` for a variant without one (ADR-0046)."""
    if rsid or not chromosome or position is None:
        return rsid
    if ref and alt:
        return f"{chromosome}:{position} {ref}>{alt}"
    return f"{chromosome}:{position}"


def annotation_label(a: Annotation) -> str:
    """How reports name the variant of ``a``: its rsID, or its locus without one."""
    chromosome, position = a.locus if a.locus is not None else ("", None)
    return variant_label(a.rsid, chromosome, position, a.locus_ref or "", a.alt)


def public_dict_label(d: dict) -> str:
    """``annotation_label`` for an annotation read back from a JSON report."""
    return variant_label(
        d.get("rsid", ""),
        d.get("chromosome", ""),
        d.get("position"),
        d.get("ref", ""),
        d.get("alt", ""),
    )

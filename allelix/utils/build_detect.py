# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Genome build detection from position data.

ADR-0021: Allelix detects the build of an input genotype file from a
handful of well-known SNP positions rather than trusting the file header.
A real-world MyHappyGenes/Tempus file was confirmed to label its build
as "37.1" while shipping GRCh38 coordinates; cross-build REF/ALT
comparison produced a false-positive pathogenic call on NIPA1.

The detection table holds authoritative (chromosome, 1-based position)
pairs for all three builds (GRCh36, GRCh37, GRCh38) across ~10 SNPs
spread over chromosomes 1, 10, 11, 12, 17, 19, and 22. Each entry's positions
differ by tens of thousands to millions of bases — there is no
coordinate ambiguity. A single matched rsID identifies a candidate build,
but ADR-0037 requires three concordant seeds before that candidate may
override other build signals.

Position data is normative; headers are not.
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING, NamedTuple

if TYPE_CHECKING:
    from collections.abc import Iterable

    from allelix.models import Variant

BUILD_GRCH36 = "GRCh36"
BUILD_GRCH37 = "GRCh37"
BUILD_GRCH38 = "GRCh38"

# Authoritative 1-based positions (ADR-0044). GRCh37 and GRCh38 come from
# the dbSNP b157 VCF releases for GRCh37.p13 and GRCh38.p14, cross-checked
# with the NCBI Variation API (refsnp SPDI placements); GRCh36 comes from
# dbSNP build 130 as published for hg18 (UCSC snp130), the last dbSNP build
# with NCBI36 placements. tests/fixtures/dbsnp_seed_placements.json records
# those sources and a test checks every entry against them. SNPs were
# chosen for:
#   - presence on virtually every consumer array
#   - clinical or pharmacogenomic relevance (so coverage is high)
#   - distribution across chromosomes so partial-coverage files still
#     hit at least one entry
#
# If the API ever returns inverted labels (mine did for chr11/12/19 due
# to NC accession version quirks), THIS table is the source of truth.
# Verify against dbSNP's web view before editing.
KNOWN_SNP_POSITIONS: dict[str, dict[str, tuple[str, int]]] = {
    # MTHFR — methylation pathway, chromosome 1 short arm
    "rs1801133": {
        BUILD_GRCH36: ("1", 11778965),
        BUILD_GRCH37: ("1", 11856378),
        BUILD_GRCH38: ("1", 11796321),
    },
    "rs1801131": {
        BUILD_GRCH36: ("1", 11777063),
        BUILD_GRCH37: ("1", 11854476),
        BUILD_GRCH38: ("1", 11794419),
    },
    # CYP2C9 / CYP2C19 cluster — chromosome 10 long arm
    "rs1799853": {
        BUILD_GRCH36: ("10", 96692037),
        BUILD_GRCH37: ("10", 96702047),
        BUILD_GRCH38: ("10", 94942290),
    },
    "rs1057910": {
        BUILD_GRCH36: ("10", 96731043),
        BUILD_GRCH37: ("10", 96741053),
        BUILD_GRCH38: ("10", 94981296),
    },
    "rs4244285": {
        BUILD_GRCH36: ("10", 96531606),
        BUILD_GRCH37: ("10", 96541616),
        BUILD_GRCH38: ("10", 94781859),
    },
    # SLCO1B1 — statin myopathy, chromosome 12
    "rs4149056": {
        BUILD_GRCH36: ("12", 21222816),
        BUILD_GRCH37: ("12", 21331549),
        BUILD_GRCH38: ("12", 21178615),
    },
    # DRD2/ANKK1 — chromosome 11
    "rs1800497": {
        BUILD_GRCH36: ("11", 112776038),
        BUILD_GRCH37: ("11", 113270828),
        BUILD_GRCH38: ("11", 113400106),
    },
    # BRCA1 — hereditary cancer, chromosome 17. An insertion: the position
    # is the VCF anchor base (dbSNP SPDI start), as ClinVar and arrays report
    # it. dbSNP build 130 (the last with NCBI36 placements) has no record for
    # it, so it has no GRCh36 entry and votes only for GRCh37 or GRCh38.
    "rs80357906": {
        BUILD_GRCH37: ("17", 41209079),
        BUILD_GRCH38: ("17", 43057062),
    },
    # APOE — chromosome 19, near telomere
    "rs429358": {
        BUILD_GRCH36: ("19", 50103781),
        BUILD_GRCH37: ("19", 45411941),
        BUILD_GRCH38: ("19", 44908684),
    },
    "rs7412": {
        BUILD_GRCH36: ("19", 50103919),
        BUILD_GRCH37: ("19", 45412079),
        BUILD_GRCH38: ("19", 44908822),
    },
    # COMT — chromosome 22
    "rs4680": {
        BUILD_GRCH36: ("22", 18331271),
        BUILD_GRCH37: ("22", 19951271),
        BUILD_GRCH38: ("22", 19963748),
    },
}


_MIN_CONFIDENT_MATCHES = 3
_MIN_CONFIDENT_RATIO = 0.8


class BuildDetectionResult(NamedTuple):
    """Outcome of build detection on an input file.

    `build` is `"GRCh36"`, `"GRCh37"`, `"GRCh38"`, or None if no known SNPs were
    found in the input. `matched` counts how many table entries matched
    the winning build; `inspected` counts how many table entries were
    found in the input (regardless of which build their positions
    matched). When `matched < inspected` the file is internally
    inconsistent (e.g., one rsID matches GRCh37, another matches
    GRCh38). The majority is authoritative only when ``is_confident``;
    otherwise callers surface diagnostics and retain their fallback build.

    `tied` names the builds that share the highest non-zero vote count when
    no single build wins (``build`` is then None and ``matched`` 0, as
    before); `tied_votes` is that shared count. Both are empty / 0 when
    there is a winner or no vote.
    """

    build: str | None
    matched: int
    inspected: int
    tied: tuple[str, ...] = ()
    tied_votes: int = 0

    @property
    def is_confident(self) -> bool:
        """True iff enough rsIDs matched and the winning ratio is at least 80%.

        Requires at least ``_MIN_CONFIDENT_MATCHES`` (3) concordant
        positions before declaring confident. ``_MIN_CONFIDENT_RATIO``
        tolerates one stale or repositioned seed once the majority is
        strong enough, while a 3/4 vote remains below threshold.
        """
        return (
            self.inspected > 0
            and self.matched >= _MIN_CONFIDENT_MATCHES
            and self.matched / self.inspected >= _MIN_CONFIDENT_RATIO
        )


# (chromosome, position) -> (seed rsID, build) for every seed on
# every build. Lets an rsID-less record (VCF ``ID=.``) vote by position.
# No coordinate repeats across seeds or builds (pinned by a test).
_SEED_AT_LOCUS: dict[tuple[str, int], tuple[str, str]] = {
    locus: (rsid, build)
    for rsid, entry in KNOWN_SNP_POSITIONS.items()
    for build, locus in entry.items()
}


def known_seed_loci() -> list[tuple[str, int]]:
    """Every seed's ``(chromosome, position)`` on every build, for indexed lookups."""
    return sorted(_SEED_AT_LOCUS)


def known_snp_seed(variant: Variant) -> str | None:
    """The build-detection seed this variant can vote for, or None.

    A variant carrying a seed rsID votes through its rsID. An rsID-less
    variant (no ``rs`` ID) at a seed's coordinate on any build votes by
    position. rsID-less resolution cannot run first because it needs the
    build being detected; the seed table already knows every build's
    coordinate, so no lookup is needed.
    """
    if variant.rsid in KNOWN_SNP_POSITIONS:
        return variant.rsid
    if variant.rsid.startswith("rs"):
        return None
    hit = _SEED_AT_LOCUS.get((variant.chromosome, variant.position))
    return hit[0] if hit is not None else None


def detect_build(variants: Iterable[Variant]) -> BuildDetectionResult:
    """Detect the genome build of an iterable of `Variant` records.

    Iterates the input, looking for any seed in `KNOWN_SNP_POSITIONS`
    (by rsID, or by position for rsID-less records; see
    `known_snp_seed`), and tallies which build's (chromosome, position)
    each match votes for. Returns when every entry in the table has been
    seen OR the input is exhausted. Streaming-friendly — does not
    materialize the full variant list.
    """
    votes: dict[str, int] = {BUILD_GRCH36: 0, BUILD_GRCH37: 0, BUILD_GRCH38: 0}
    inspected = 0
    remaining = set(KNOWN_SNP_POSITIONS)
    for variant in variants:
        seed = known_snp_seed(variant)
        if seed not in remaining:
            continue
        assert seed is not None
        entry = KNOWN_SNP_POSITIONS[seed]
        remaining.discard(seed)
        inspected += 1
        for build, (chrom, pos) in entry.items():
            if variant.chromosome == chrom and variant.position == pos:
                votes[build] += 1
                break
        if not remaining:
            break

    if inspected == 0:
        return BuildDetectionResult(build=None, matched=0, inspected=0)

    winner = max(votes, key=votes.__getitem__)
    if votes[winner] == 0:
        return BuildDetectionResult(build=None, matched=0, inspected=inspected)
    # Tie between builds with equal non-zero votes: no winner, but say which.
    top = votes[winner]
    tied = tuple(sorted(b for b, n in votes.items() if n == top))
    if len(tied) > 1:
        return BuildDetectionResult(
            build=None, matched=0, inspected=inspected, tied=tied, tied_votes=top
        )
    return BuildDetectionResult(build=winner, matched=votes[winner], inspected=inspected)


_HG_TO_BUILD = {"18": BUILD_GRCH36, "19": BUILD_GRCH37, "38": BUILD_GRCH38}
_NUM_TO_BUILD = {"36": BUILD_GRCH36, "37": BUILD_GRCH37, "38": BUILD_GRCH38}

# build-label patterns are anchored on word boundaries so a date
# (`2038-01-01`) or version string (`v37`) doesn't false-match by bare
# substring containment. Each pattern captures one of the canonical
# numeric tokens; a label that surfaces *multiple distinct builds* (e.g.
# `hg19/hg38`, `38 (37 liftover)`) is ambiguous and returns None — the
# safer answer than picking one substring-wins.
_BUILD_TOKEN_PATTERNS: tuple[tuple[re.Pattern[str], dict[str, str]], ...] = (
    (re.compile(r"\bgrch(36|37|38)\b", re.IGNORECASE), _NUM_TO_BUILD),
    (re.compile(r"\bhg(18|19|38)\b", re.IGNORECASE), _HG_TO_BUILD),
    (re.compile(r"\bncbi\s*0*(36|37|38)\b", re.IGNORECASE), _NUM_TO_BUILD),
    (re.compile(r"\bbuild\s*0*(36|37|38)\b", re.IGNORECASE), _NUM_TO_BUILD),
    (re.compile(r"\b(36|37|38)\b"), _NUM_TO_BUILD),
)


def normalize_build_label(label: str | None) -> str | None:
    """Map a human-written build label to canonical `GRCh36`, `GRCh37`, or `GRCh38`.

    Examples that map to GRCh36: `"GRCh36"`, `"hg18"`, `"build 36"`,
    `"NCBI 36"`. Examples for GRCh37: `"GRCh37"`, `"grch37"`, `"hg19"`,
    `"37.1"`, `"build 37.1"`, `"NCBI 37"`. Examples for GRCh38: `"GRCh38"`,
    `"hg38"`, `"38"`. Unrecognized labels return None.

    Used to compare a file's header-claimed build against the detected
    build. The label space is informal and provider-specific; this
    function only recognizes well-known aliases.

    Tokens are matched on word boundaries (so `2038-01-01` does not
    register as GRCh38, and `v37` does not register as GRCh37). If a
    label mentions more than one build (`hg19/hg38`, `38 (37 liftover)`),
    that is treated as ambiguous and returns None — picking either one
    risks shipping a mislabeled-build comparison further down the chain.
    """
    if not label:
        return None
    found: set[str] = set()
    for regex, mapping in _BUILD_TOKEN_PATTERNS:
        for m in regex.finditer(label):
            found.add(mapping[m.group(1)])
    if len(found) == 1:
        return next(iter(found))
    return None


# ── Contig lengths (ADR-0044) ─────────────────────────────────────────
# Nuclear contig lengths of the two supported assemblies, as declared in the
# ``##contig`` lines of the GIAB v4.2.1 HG002 benchmark VCFs (GRCh37: hs37d5;
# GRCh38: no-alt analysis set). Every nuclear contig differs between the two
# builds, so any one of them identifies the build. MT is 16,569 bp in both
# and is not evidence.
_CONTIG_LENGTHS: dict[str, dict[str, int]] = {
    BUILD_GRCH37: {
        "1": 249250621,
        "2": 243199373,
        "3": 198022430,
        "4": 191154276,
        "5": 180915260,
        "6": 171115067,
        "7": 159138663,
        "8": 146364022,
        "9": 141213431,
        "10": 135534747,
        "11": 135006516,
        "12": 133851895,
        "13": 115169878,
        "14": 107349540,
        "15": 102531392,
        "16": 90354753,
        "17": 81195210,
        "18": 78077248,
        "19": 59128983,
        "20": 63025520,
        "21": 48129895,
        "22": 51304566,
        "X": 155270560,
        "Y": 59373566,
    },
    BUILD_GRCH38: {
        "1": 248956422,
        "2": 242193529,
        "3": 198295559,
        "4": 190214555,
        "5": 181538259,
        "6": 170805979,
        "7": 159345973,
        "8": 145138636,
        "9": 138394717,
        "10": 133797422,
        "11": 135086622,
        "12": 133275309,
        "13": 114364328,
        "14": 107043718,
        "15": 101991189,
        "16": 90338345,
        "17": 83257441,
        "18": 80373285,
        "19": 58617616,
        "20": 64444167,
        "21": 46709983,
        "22": 50818468,
        "X": 156040895,
        "Y": 57227415,
    },
}
_NUCLEAR = frozenset(_CONTIG_LENGTHS[BUILD_GRCH37])
# UCSC hg19's chrM is the 16,571 bp Yoruba sequence, not the rCRS that
# GRCh37 MT and ClinVar positions use.
HG19_CHRM_LENGTH = 16571


class LengthVerdict(NamedTuple):
    """What a VCF's ``##contig`` lengths say about its build.

    ``status`` is ``"none"`` (no nuclear contig declares a length),
    ``"match"`` (every declared nuclear length is one build's; ``build``
    names it), ``"mixed"`` (some lengths are GRCh37's, others GRCh38's) or
    ``"unsupported"`` (a nuclear length matches neither build). ``detail``
    is a short human-readable account; ``skip_mt`` is True when chrM/MT is
    declared with the hg19 length, whose coordinates are not rCRS.
    """

    status: str
    build: str | None
    detail: str
    skip_mt: bool = False


def _bare_contig(name: str) -> str:
    """``chr1`` -> ``1``, ``chrX`` -> ``X``, ``chrM`` / ``M`` -> ``MT``."""
    bare = (name[3:] if name.lower().startswith("chr") else name).upper()
    return "MT" if bare in ("M", "MT") else bare


def contig_length_verdict(lengths: dict[str, int]) -> LengthVerdict:
    """Classify declared contig lengths (``{contig ID: length}``)."""
    seen: dict[str, set[str]] = {}
    unsupported: list[str] = []
    skip_mt = False
    for name, length in lengths.items():
        bare = _bare_contig(name)
        if bare == "MT":
            skip_mt = skip_mt or length == HG19_CHRM_LENGTH
            continue
        if bare not in _NUCLEAR:
            continue  # alt, decoy, random and unplaced contigs are ignored
        builds = {b for b, table in _CONTIG_LENGTHS.items() if table[bare] == length}
        if not builds:
            unsupported.append(f"{name}={length:,}")
            continue
        for build in builds:
            seen.setdefault(build, set()).add(name)
    if unsupported:
        return LengthVerdict(
            "unsupported",
            None,
            "contig lengths match neither GRCh37 nor GRCh38 (" + ", ".join(unsupported[:3]) + ")",
            skip_mt,
        )
    if not seen:
        return LengthVerdict("none", None, "no nuclear contig declares a length", skip_mt)
    if len(seen) > 1:
        parts = [
            f"{', '.join(sorted(names)[:3])} = {build}" for build, names in sorted(seen.items())
        ]
        return LengthVerdict(
            "mixed", None, "contig lengths disagree (" + "; ".join(parts) + ")", skip_mt
        )
    build, names = next(iter(seen.items()))
    return LengthVerdict("match", build, f"{len(names)} contig length(s) match {build}", skip_mt)


# ── The build decision (ADR-0044) ─────────────────────────────────────
# Fixed vocabulary for where the build came from (JSON ``input.build_source``).
SOURCE_OVERRIDE = "override"
SOURCE_CONTIG_LENGTHS = "contig_lengths"
SOURCE_SEEDS = "seeds"
SOURCE_ASSEMBLY_HEADER = "assembly_header"
SOURCE_FILE_LABEL = "file_label"
SOURCE_DEFAULT = "default"
BUILD_SOURCES = (
    SOURCE_OVERRIDE,
    SOURCE_CONTIG_LENGTHS,
    SOURCE_SEEDS,
    SOURCE_ASSEMBLY_HEADER,
    SOURCE_FILE_LABEL,
    SOURCE_DEFAULT,
)
# Short reader-facing wording for each source (HTML report, ``stats``).
BUILD_SOURCE_LABELS = {
    SOURCE_OVERRIDE: "forced with --build",
    SOURCE_CONTIG_LENGTHS: "from the contig lengths",
    SOURCE_SEEDS: "from known-SNP positions",
    SOURCE_ASSEMBLY_HEADER: "from the VCF assembly header",
    SOURCE_FILE_LABEL: "from the file's build line",
    SOURCE_DEFAULT: "format default; the file states no build",
}

KIND_VCF = "vcf"
KIND_ARRAY = "array"


# The builds --build (and the Python API's override) may force. GRCh36 runs
# with position lookups off (ADR-0025).
SUPPORTED_OVERRIDES = (BUILD_GRCH36, BUILD_GRCH37, BUILD_GRCH38)


class BuildDecisionError(ValueError):
    """The build cannot be decided safely; the run must stop.

    The message says what is wrong, what evidence the file carries and
    what to supply. Nothing has been written when it is raised.
    """


class BuildDecision(NamedTuple):
    """The decided build, where it came from, and the evidence behind it."""

    build: str
    source: str
    file_label: str | None
    seeds: BuildDetectionResult
    lengths: LengthVerdict | None
    label_overridden: bool = False
    tentative_grch36: bool = False
    skip_mt: bool = False


NO_SEEDS = BuildDetectionResult(build=None, matched=0, inspected=0)


def _seed_text(seeds: BuildDetectionResult) -> str:
    if seeds.inspected == 0:
        return "no known-SNP positions in the file"
    if seeds.tied:
        return (
            f"{seeds.inspected} known-SNP positions found, tied between "
            f"{' and '.join(seeds.tied)} ({seeds.tied_votes} each)"
        )
    if seeds.build is None:
        return (
            f"{seeds.inspected} known-SNP position(s) found, none at any supported "
            "build's coordinate"
        )
    strength = "confident" if seeds.is_confident else "tentative"
    return (
        f"{seeds.matched} of {seeds.inspected} known-SNP positions match {seeds.build} "
        f"({strength})"
    )


def _stop(file_name: str, problem: str, evidence: list[str]) -> BuildDecisionError:
    return BuildDecisionError(
        f"Cannot use {file_name}: {problem}. Evidence: {'; '.join(evidence)}. "
        "Tell Allelix which genome build the file uses (command line: --build grch37 "
        "or --build grch38; for an old GRCh36 / hg18 file, --build grch36, see "
        "docs/grch36-liftover.md), and run it again."
    )


def decide_build(
    *,
    file_name: str,
    kind: str,
    override: str | None,
    file_label: str | None,
    default: str | None,
    lengths: LengthVerdict | None,
    seeds: BuildDetectionResult,
) -> BuildDecision:
    """Decide the build from the evidence the file carries (ADR-0044).

    ``kind`` is ``"vcf"`` or ``"array"``. ``file_label`` is the build the
    file states (a VCF's ``assembly=``, an array's header line), or None.
    ``default`` is the parser's documented default (arrays only).
    ``lengths`` is the VCF contig-length verdict (None for arrays).
    ``seeds`` is the known-SNP evidence from the whole file (pre-scan), or
    ``NO_SEEDS`` when contig lengths decide and the seeds are only checked
    during the annotation pass (``check_late_seeds``).

    Raises ``BuildDecisionError`` when the run must stop.
    """
    if override and override not in SUPPORTED_OVERRIDES:
        raise ValueError(
            f"Unsupported build override {override!r}: use one of {', '.join(SUPPORTED_OVERRIDES)}"
        )
    if override:
        return BuildDecision(override, SOURCE_OVERRIDE, file_label, seeds, lengths)
    label_text = f"the file says {file_label}" if file_label else "the file states no build"
    if lengths is not None and lengths.status != "none":
        if lengths.status == "unsupported":
            raise _stop(file_name, "unsupported assembly", [lengths.detail, label_text])
        if lengths.status == "mixed":
            raise _stop(file_name, "its contig lengths contradict each other", [lengths.detail])
        assert lengths.build is not None
        if file_label and file_label != lengths.build:
            raise _stop(
                file_name,
                f"its header says {file_label} but its contig lengths are {lengths.build}",
                [lengths.detail, label_text],
            )
        if seeds.build is not None and seeds.build != lengths.build:
            raise _stop(
                file_name,
                f"its contig lengths say {lengths.build} "
                f"but its known-SNP positions say {seeds.build}",
                [lengths.detail, _seed_text(seeds)],
            )
        if seeds.tied and lengths.build not in seeds.tied:
            raise _stop(
                file_name,
                f"its contig lengths say {lengths.build} but its known-SNP positions "
                f"are tied between {' and '.join(seeds.tied)}",
                [lengths.detail, _seed_text(seeds)],
            )
        return BuildDecision(
            lengths.build,
            SOURCE_CONTIG_LENGTHS,
            file_label,
            seeds,
            lengths,
            skip_mt=lengths.skip_mt,
        )
    if seeds.inspected > 0 and seeds.build is None and not seeds.tied:
        # Seeds found but placed on no build are evidence against
        # every build, not an absence of evidence; the label or the default
        # cannot be trusted either. (When contig lengths decided, they are
        # decisive and this point is not reached.)
        raise _stop(
            file_name,
            "its known-SNP positions match none of the supported builds",
            [_seed_text(seeds), label_text],
        )
    if seeds.is_confident:
        assert seeds.build is not None
        return BuildDecision(
            seeds.build,
            SOURCE_SEEDS,
            file_label,
            seeds,
            lengths,
            label_overridden=bool(file_label and file_label != seeds.build),
        )
    if BUILD_GRCH36 in seeds.tied:
        # When contig lengths did not decide, a tie that includes
        # GRCh36 takes the GRCh36 fail-safe (ADR-0025) whatever the label
        # says. A wrong GRCh36 only turns position lookups off, and says so;
        # a wrong modern build on a real GRCh36 file attaches wrong records.
        return BuildDecision(
            BUILD_GRCH36,
            SOURCE_SEEDS,
            file_label,
            seeds,
            lengths,
            label_overridden=bool(file_label and file_label != BUILD_GRCH36),
            tentative_grch36=True,
        )
    tentative = seeds.build
    if tentative == BUILD_GRCH36:
        # Fail-safe (ADR-0025): GRCh36 turns position lookups off, so a
        # tentative GRCh36 winner cannot attach another build's data.
        return BuildDecision(
            BUILD_GRCH36,
            SOURCE_SEEDS,
            file_label,
            seeds,
            lengths,
            label_overridden=bool(file_label and file_label != BUILD_GRCH36),
            tentative_grch36=True,
        )
    if file_label:
        if tentative is not None and tentative != file_label:
            raise _stop(
                file_name,
                f"its known-SNP positions point to {tentative} but {label_text}",
                [_seed_text(seeds), label_text],
            )
        if seeds.tied and file_label not in seeds.tied:
            raise _stop(
                file_name,
                f"its known-SNP positions are tied between {' and '.join(seeds.tied)} "
                f"but {label_text}",
                [_seed_text(seeds), label_text],
            )
        source = SOURCE_ASSEMBLY_HEADER if kind == KIND_VCF else SOURCE_FILE_LABEL
        return BuildDecision(file_label, source, file_label, seeds, lengths)
    if kind == KIND_VCF or default is None:
        evidence = [_seed_text(seeds), label_text]
        if lengths is not None:
            evidence.insert(0, lengths.detail)
        raise _stop(
            file_name,
            "the file carries no contig lengths, no assembly and not enough known-SNP positions",
            evidence,
        )
    if tentative is not None and tentative != default:
        raise _stop(
            file_name,
            f"its known-SNP positions point to {tentative}, "
            f"not the {default} default for this format",
            [_seed_text(seeds), label_text],
        )
    # A tie that excludes the default needs a third build besides the default;
    # with GRCh37 as every format's default that tie includes GRCh36 and took
    # the GRCh36 fail-safe above, so it never reaches this point.
    if seeds.inspected == 0:
        # No label and no evidence at all (for example an X-only
        # export, which can carry none of the autosomal seeds). The default
        # needs at least one seed position that does not contradict it.
        raise _stop(
            file_name,
            "the file states no genome build and carries no known-SNP positions to check one",
            [_seed_text(seeds), label_text],
        )
    return BuildDecision(default, SOURCE_DEFAULT, None, seeds, lengths)


def check_late_seeds(decision: BuildDecision, seeds: BuildDetectionResult, file_name: str) -> None:
    """Stop when seeds seen during the annotation pass contradict contig lengths.

    Contig lengths decide a VCF's build without a pre-scan; the seeds the
    annotation pass meets are checked afterwards. Any seed majority for
    another build stops the run before anything is written, and so does a
    tie that does not include the contig-length build (ADR-0044). A tie
    that includes it lets the run proceed; the caller shows the tie as an
    inconclusive warning.
    """
    if decision.source != SOURCE_CONTIG_LENGTHS:
        return
    if seeds.tied and decision.build not in seeds.tied:
        assert decision.lengths is not None
        raise _stop(
            file_name,
            f"its contig lengths say {decision.build} but its known-SNP positions are "
            f"tied between {' and '.join(seeds.tied)} (found during the annotation pass, "
            "which ran in full before this check; nothing was written)",
            [decision.lengths.detail, _seed_text(seeds)],
        )
    if seeds.build is None:
        return
    if seeds.build != decision.build:
        assert decision.lengths is not None
        raise _stop(
            file_name,
            f"its contig lengths say {decision.build} but its known-SNP positions say "
            f"{seeds.build} (found during the annotation pass, which ran in full before "
            "this check; nothing was written)",
            [decision.lengths.detail, _seed_text(seeds)],
        )

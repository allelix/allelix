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
ambiguity. A single matched rsID identifies the build; multiple are
confirmatory.

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

# Authoritative 1-based positions per NCBI dbSNP / Variation API. Each
# entry was cross-checked against the API's SPDI (0-based) + 1 and the
# correct NC accession version for each build. SNPs were chosen for:
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
        BUILD_GRCH36: ("10", 96692448),
        BUILD_GRCH37: ("10", 96702047),
        BUILD_GRCH38: ("10", 94942290),
    },
    "rs1057910": {
        BUILD_GRCH36: ("10", 96731043),
        BUILD_GRCH37: ("10", 96741053),
        BUILD_GRCH38: ("10", 94981296),
    },
    "rs4244285": {
        BUILD_GRCH36: ("10", 96532017),
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
    # BRCA1 — hereditary cancer, chromosome 17
    "rs80357906": {
        BUILD_GRCH36: ("17", 38449327),
        BUILD_GRCH37: ("17", 41209080),
        BUILD_GRCH38: ("17", 43057063),
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


class BuildDetectionResult(NamedTuple):
    """Outcome of build detection on an input file.

    `build` is `"GRCh36"`, `"GRCh37"`, `"GRCh38"`, or None if no known SNPs were
    found in the input. `matched` counts how many table entries matched
    the winning build; `inspected` counts how many table entries were
    found in the input (regardless of which build their positions
    matched). When `matched < inspected` the file is internally
    inconsistent (e.g., one rsID matches GRCh37, another matches
    GRCh38) — surface a warning but pick the majority.
    """

    build: str | None
    matched: int
    inspected: int

    @property
    def is_confident(self) -> bool:
        """True iff enough rsIDs matched and all matches agreed.

        Requires at least ``_MIN_CONFIDENT_MATCHES`` (3) concordant
        positions before declaring confident. A single-SNP match
        could be a table error; three concordant matches across
        different chromosomes eliminates that risk.
        """
        return self.matched >= _MIN_CONFIDENT_MATCHES and self.matched == self.inspected


def detect_build(variants: Iterable[Variant]) -> BuildDetectionResult:
    """Detect the genome build of an iterable of `Variant` records.

    Iterates the input, looking for any rsID in `KNOWN_SNP_POSITIONS`,
    and tallies which build's (chromosome, position) each match votes
    for. Returns when every entry in the table has been seen OR the
    input is exhausted. Streaming-friendly — does not materialize the
    full variant list.
    """
    votes: dict[str, int] = {BUILD_GRCH36: 0, BUILD_GRCH37: 0, BUILD_GRCH38: 0}
    inspected = 0
    remaining = set(KNOWN_SNP_POSITIONS)
    for variant in variants:
        if variant.rsid not in remaining:
            continue
        entry = KNOWN_SNP_POSITIONS[variant.rsid]
        remaining.discard(variant.rsid)
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
    # Tie between two builds with equal non-zero votes — don't pick.
    top_counts = sorted(votes.values(), reverse=True)
    if top_counts[0] == top_counts[1]:
        return BuildDetectionResult(build=None, matched=0, inspected=inspected)
    return BuildDetectionResult(build=winner, matched=votes[winner], inspected=inspected)


_HG_TO_BUILD = {"18": BUILD_GRCH36, "19": BUILD_GRCH37, "38": BUILD_GRCH38}
_NUM_TO_BUILD = {"36": BUILD_GRCH36, "37": BUILD_GRCH37, "38": BUILD_GRCH38}

# GH #16: build-label patterns are anchored on word boundaries so a date
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

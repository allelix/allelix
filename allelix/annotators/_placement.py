# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Place the rows of files that carry no positions, from the rsID identity table.

An rsID and genotype list gives only an ID and a call. Each row whose
rsID the identity table places at one locus on GRCh38 (ADR-0040) is
placed there when exactly one strand reading of its call fits the site's
forward alleles, REF plus every ALT (ADR-0066, through
``forward_strand_reading``; ADR-0069):

    - only the call as written fits: placed as given;
    - only its complement fits: placed, the call complemented to the
      forward strand and counted (``complemented``);
    - both fit (``palindromic`` at an A/T or C/G site, ``both_fit`` at a
      site with more alleles) or neither fits (``no_fit``): unplaced.

Each row is read on its own: the file's other rows never decide its
strand. Rows are also left unplaced, and counted by reason, when the ID is
not an rsID, when the table holds no single placement for it, or when the
call is a D/I call. An unplaced row keeps chromosome "" and position 0:
position-keyed sources (ClinVar, gnomAD, AlphaMissense, CADD) skip it, and
rsID-keyed sources still look it up, under their own strand checks. Its
findings publish a null chromosome and position. A no-call is placed
without a check, as it carries no alleles to orient.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING

from allelix.models import NO_CALL_MARKER
from allelix.utils.allele import (
    STRAND_BOTH_FIT,
    STRAND_COMPLEMENT,
    STRAND_DIRECT,
    STRAND_NO_FIT,
    STRAND_PALINDROMIC,
    complement,
    forward_strand_reading,
)

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator

    from allelix.annotators.rsid_identity import ForwardSite, RsidIdentityAnnotator
    from allelix.models import Variant

_RSID = re.compile(r"rs\d+")
_CHUNK = 5000  # rows per identity lookup; the pipeline's batch size

# Reasons a row stays unplaced, in the order reports list them.
NOT_AN_RSID = "not_an_rsid"
NOT_IN_TABLE = "not_in_table"
DI_CALL = "di_call"
PALINDROMIC = STRAND_PALINDROMIC
BOTH_FIT = STRAND_BOTH_FIT
NO_FIT = STRAND_NO_FIT
UNPLACED_REASONS = (NOT_AN_RSID, NOT_IN_TABLE, PALINDROMIC, BOTH_FIT, NO_FIT, DI_CALL)

# Why a row of a position-less file stayed unplaced, in report words (ADR-0069).
UNPLACED_WORDS = {
    "not_an_rsid": "no rsID",
    "not_in_table": "rsID not placed once on GRCh38 in the identity table",
    "palindromic": "A/T or C/G site (strand cannot be told)",
    "both_fit": "call fits both strands at the site",
    "no_fit": "call fits neither strand at the site",
    "di_call": "D/I call",
}


@dataclass
class PlacementStats:
    """What placement did with a position-less file's rows."""

    build: str
    rows: int = 0
    placed: int = 0
    complemented: int = 0
    unplaced: dict[str, int] = field(default_factory=lambda: dict.fromkeys(UNPLACED_REASONS, 0))
    # Each unplaced row with an rsID: rsID -> (call as written, reason), so
    # the report can name the ones at high-value or P/LP rsIDs. Not in
    # ``to_dict``; one entry per rsID.
    unplaced_rows: dict[str, tuple[str, str]] = field(default_factory=dict)

    @property
    def unplaced_total(self) -> int:
        """Rows left unplaced, for every reason."""
        return sum(self.unplaced.values())

    def to_dict(self) -> dict[str, object]:
        """JSON-ready record (``input.placement``)."""
        return {
            "build": self.build,
            "rows": self.rows,
            "placed": self.placed,
            "complemented": self.complemented,
            "unplaced": dict(self.unplaced),
        }


def place_stream(
    variants: Iterable[Variant],
    identity: RsidIdentityAnnotator,
    build: str,
    stats: PlacementStats,
) -> Iterator[Variant]:
    """Yield ``variants`` placed on ``build`` where one reading fits, else unplaced.

    Streams in chunks of ``_CHUNK`` rows, one identity query per chunk.
    """
    chunk: list[Variant] = []
    for variant in variants:
        chunk.append(variant)
        if len(chunk) >= _CHUNK:
            yield from _place_chunk(chunk, identity, build, stats)
            chunk = []
    if chunk:
        yield from _place_chunk(chunk, identity, build, stats)


def _place_chunk(
    chunk: list[Variant],
    identity: RsidIdentityAnnotator,
    build: str,
    stats: PlacementStats,
) -> Iterator[Variant]:
    sites = identity.forward_sites({v.rsid for v in chunk if _RSID.fullmatch(v.rsid)}, build)
    for variant in chunk:
        stats.rows += 1
        placed, reason = _place(variant, sites.get(variant.rsid))
        if placed is None:
            stats.unplaced[reason] += 1
            if reason != NOT_AN_RSID:
                stats.unplaced_rows.setdefault(variant.rsid, (variant.genotype, reason))
            yield variant
            continue
        stats.placed += 1
        if (placed.allele1, placed.allele2) != (variant.allele1, variant.allele2):
            stats.complemented += 1
        yield replace(placed, build=build)


def _place(variant: Variant, site: ForwardSite | None) -> tuple[Variant | None, str]:
    """``(placed variant, "")``, or ``(None, reason)`` when the row stays unplaced."""
    if not _RSID.fullmatch(variant.rsid):
        return None, NOT_AN_RSID
    if site is None:
        return None, NOT_IN_TABLE
    chrom, pos, forward = site
    at_locus = replace(variant, chromosome=chrom, position=pos)
    if variant.allele1 == NO_CALL_MARKER:
        return at_locus, ""
    if {variant.allele1, variant.allele2} & {"D", "I"}:
        return None, DI_CALL
    reading = forward_strand_reading((variant.allele1, variant.allele2), forward)
    if reading == STRAND_DIRECT:
        return at_locus, ""
    if reading == STRAND_COMPLEMENT:
        return (
            replace(
                at_locus,
                allele1=complement(variant.allele1),
                allele2=complement(variant.allele2),
            ),
            "",
        )
    return None, reading

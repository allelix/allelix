# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Which requested rsID an input record speaks for (2.5.2).

One predicate for panel coverage and the high-value scan. When the rsID
identity table places the rsID on the build, a record speaks for it when
its chromosome and position are one of those placements and its REF equals
that placement's REF (case-insensitive). The record's ID and ALT do not
decide it: a record with another ID, or none, at the placement is that
rsID's record, and one carrying the rsID's ID elsewhere is not. With no
placement for the rsID, the record's ID decides, as before. An input that
supplies no REF (array and text-chip files) is matched by ID too
(supervisor ruling): annotation of such files matches by rsID, and
coverage must agree with it.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping

# (chromosome, position, REF in upper case) of one identity placement.
Placement = tuple[str, int, str]
Placements = dict[str, frozenset[Placement]]

_BUILDS = frozenset({"GRCh37", "GRCh38"})


def build_placements(
    rows: Mapping[str, Iterable[tuple[str, int, str, str, int, str, str]]], build: str
) -> Placements:
    """Each rsID's placements on ``build`` from identity rows.

    ``rows`` is ``RsidIdentityAnnotator.placements``: ``(chrom, pos37,
    ref37, alt37, pos38, ref38, alt38)`` per row. Several rows at one
    placement (a multi-allelic site) give one placement.
    """
    if build not in _BUILDS:
        return {}
    grch37 = build == "GRCh37"
    out: Placements = {}
    for rsid, rsid_rows in rows.items():
        out[rsid] = frozenset(
            (chrom, pos37, ref37.upper()) if grch37 else (chrom, pos38, ref38.upper())
            for chrom, pos37, ref37, _alt37, pos38, ref38, _alt38 in rsid_rows
        )
    return out


def placement_index(placements: Placements) -> dict[tuple[str, int], frozenset[str]]:
    """The rsIDs placed at each ``(chromosome, position)``."""
    index: dict[tuple[str, int], set[str]] = {}
    for rsid, places in placements.items():
        for chrom, pos, _ref in places:
            index.setdefault((chrom, pos), set()).add(rsid)
    return {locus: frozenset(rsids) for locus, rsids in index.items()}


def attributed(
    rsid: str,
    chrom: str,
    pos: int,
    ref: str | None,
    record_id: str,
    placements: Placements,
) -> bool:
    """Whether a record at ``chrom:pos`` with ``ref`` and ID ``record_id`` speaks for ``rsid``.

    ``ref`` is None when the input supplied no REF; such a record is
    matched by ID. ``chrom`` is the normalized chromosome name.
    """
    places = placements.get(rsid)
    if not places or ref is None:
        return record_id == rsid
    upper = ref.upper()
    return any((chrom, pos) == (c, p) and upper == r for c, p, r in places)


def attributed_rsids(
    wanted: Iterable[str] | frozenset[str] | set[str],
    chrom: str,
    pos: int,
    ref: str | None,
    record_id: str,
    placements: Placements,
    index: Mapping[tuple[str, int], frozenset[str]],
) -> set[str]:
    """The requested rsIDs a record speaks for: its own ID, and rsIDs placed at its locus."""
    wanted_set = wanted if isinstance(wanted, (set, frozenset)) else set(wanted)
    candidates = set(index.get((chrom, pos), ()))
    if record_id in wanted_set:
        candidates.add(record_id)
    return {
        rsid
        for rsid in candidates
        if rsid in wanted_set and attributed(rsid, chrom, pos, ref, record_id, placements)
    }

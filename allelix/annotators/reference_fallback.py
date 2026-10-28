# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""The normal allele at a reference-carried site without the input build's ClinVar.

ClinPGx and SNPedia treat a homozygous-reference call as a non-finding, and
derive their ``alt`` against the normal allele (ADR-0023). Where the build's
reference itself carries a variant allele (T at F5 rs6025 on GRCh37), the
normal allele is not the REF, and ClinVar's cache for that build says so.
Without that cache, this module answers for those sites only: rsIDs whose
GRCh37 and GRCh38 single-base REFs differ in the rsID identity table. Every
other site is left as it was (no suppression without ClinVar).

The normal allele comes from, in order:

1. the other build's ClinVar, when it classifies the input build's REF as
   the variant allele (a row REF = other build's REF, ALT = this build's
   REF): the other build's REF is normal;
2. gnomAD's frequencies at the site: the major allele, among the REF and
   every ALT gnomAD lists, is normal (a tie decides nothing);
3. otherwise nothing: the site abstains, nothing is suppressed, and the
   abstention is counted.

Which caches the run reads is decided once, at the first call: a chip
makes this call for every rsID, and a readiness check opens the cache.
``prefetch`` decides a batch's sites with one chunked query per source.

GRCh38's REF is never taken as normal on its own evidence: hiding a real
pharmacogenomic finding is worse than showing a normal genotype's row.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from allelix.utils.build_detect import BUILD_GRCH37, BUILD_GRCH38

if TYPE_CHECKING:
    from collections.abc import Iterable

    from allelix.annotators.clinvar import ClinVarAnnotator
    from allelix.annotators.gnomad import GnomadAnnotator
    from allelix.annotators.rsid_identity import RsidIdentityAnnotator
    from allelix.models import Variant

# One identity row's placement and REF on each build.
_Site = tuple[str, int, str, int, str]


class ReferenceCarriedFallback:
    """The ``reference_at`` / ``normal_allele_for`` pair ClinPGx and SNPedia use.

    ClinVar answers whenever it has a readable cache for the call's build;
    otherwise the fallback above does, at reference-carried sites only.
    """

    def __init__(
        self,
        clinvar: ClinVarAnnotator,
        identity: RsidIdentityAnnotator,
        gnomad: GnomadAnnotator,
    ) -> None:
        """Bind the sources; nothing is read until a call needs it."""
        self._clinvar = clinvar
        self._identity = identity
        self._gnomad = gnomad
        self._sites: dict[str, _Site] | None = None
        self._reads: dict[str, bool] = {}
        self._gnomad_ready: bool | None = None
        self._normal: dict[tuple[str, str], str | None] = {}
        self._abstained: set[tuple[str, str]] = set()

    def reference_at(self, rsid: str, build: str, locus: tuple[str, int]) -> str | None:
        """The allele a homozygous non-finding call carries at ``locus``, or None."""
        if self._reads_build(build):
            return self._clinvar.reference_at(rsid, build, locus)
        if not self._placed(rsid, build, locus):
            return None
        return self._normal_for(rsid, build)

    def normal_allele_for(
        self, rsid: str, build: str, locus: tuple[str, int] | None = None
    ) -> str | None:
        """The normal allele where ``build``'s reference carries a variant, or None.

        Without the build's ClinVar, only a call at the reference-differing
        placement itself gets one: the rsID's other placements, and a call
        with no position, keep deriving their ``alt`` against their REF.
        """
        if self._reads_build(build):
            return self._clinvar.normal_allele_for(rsid, build, locus)
        if locus is None or not self._placed(rsid, build, locus):
            return None
        return self._normal_for(rsid, build)

    def _placed(self, rsid: str, build: str, locus: tuple[str, int]) -> bool:
        """True when ``locus`` is the rsID's reference-differing placement on ``build``."""
        site = self._site(rsid, build)
        return site is not None and locus == (
            site[0],
            site[1] if build == BUILD_GRCH37 else site[3],
        )

    def prefetch(self, variants: Iterable[Variant]) -> None:
        """Decide every reference-carried site among ``variants`` in bulk.

        ``reference_at`` and ``normal_allele_for`` then answer those sites
        from memory; a site not prefetched is still decided on its own.
        """
        wanted: dict[str, set[str]] = {}
        for variant in variants:
            build = variant.build
            key = (variant.rsid, build)
            if key in self._normal or self._reads_build(build):
                continue
            if self._site(variant.rsid, build) is not None:
                wanted.setdefault(build, set()).add(variant.rsid)
        for build, rsids in wanted.items():
            self._decide_many(sorted(rsids), build)

    def abstentions(self) -> int:
        """Reference-carried sites whose normal allele could not be told, per build."""
        return len(self._abstained)

    def _reads_build(self, build: str) -> bool:
        if build not in self._reads:
            self._reads[build] = self._clinvar.reads_build(build)
        return self._reads[build]

    def _site(self, rsid: str, build: str) -> _Site | None:
        if build not in (BUILD_GRCH37, BUILD_GRCH38):
            return None
        if self._sites is None:
            ready = self._identity.is_ready()
            self._sites = self._identity.reference_differing_sites() if ready else {}
        return self._sites.get(rsid)

    def _normal_for(self, rsid: str, build: str) -> str | None:
        key = (rsid, build)
        if key not in self._normal:
            self._decide_many([rsid], build)
        return self._normal[key]

    def _gnomad_reads(self) -> bool:
        if self._gnomad_ready is None:
            self._gnomad_ready = self._gnomad.is_ready()
        return self._gnomad_ready

    def _decide_many(self, rsids: list[str], build: str) -> None:
        """Decide ``rsids`` (reference-carried sites on ``build``) with one query per source."""
        assert self._sites is not None
        other_build = BUILD_GRCH38 if build == BUILD_GRCH37 else BUILD_GRCH37
        rows: dict[str, set[tuple[str, str]]] = {}
        if self._reads_build(other_build):
            rows = self._clinvar.snv_alleles_many(rsids, other_build)
        undecided: list[str] = []
        for rsid in rsids:
            _chrom, _pos37, ref37, _pos38, ref38 = self._sites[rsid]
            here, other = (ref37, ref38) if build == BUILD_GRCH37 else (ref38, ref37)
            alleles = rows.get(rsid, set())
            # ClinVar's variant allele is this build's REF: the other REF is normal.
            if (other, here) in alleles:
                self._record(rsid, build, other)
            # The other build's reference carries ClinVar's variant (REF == ALT).
            elif (other, other) in alleles:
                self._record(rsid, build, here)
            else:
                undecided.append(rsid)
        frequencies: dict[tuple[str, int, str], dict[str, float | None]] = {}
        if undecided and self._gnomad_reads():
            frequencies = self._gnomad.sites_frequencies(
                {(s[0], s[3], s[4]) for s in (self._sites[r] for r in undecided)}
            )
        for rsid in undecided:
            chrom, _pos37, _ref37, pos38, ref38 = self._sites[rsid]
            alts = frequencies.get((chrom, pos38, ref38), {})
            self._record(rsid, build, _major_allele(ref38, alts))

    def _record(self, rsid: str, build: str, normal: str | None) -> None:
        self._normal[(rsid, build)] = normal
        if normal is None:
            self._abstained.add((rsid, build))


# Frequencies closer than this are a tie: gnomAD's AFs are decimals, and the
# REF's, computed as what the ALTs leave, carries floating-point error.
_TIE = 1e-9


def _major_allele(ref: str, alts: dict[str, float | None]) -> str | None:
    """The allele gnomAD shows most often at a site, or None when that is not known.

    ``alts`` holds every ALT gnomAD lists at the site (keyed on GRCh38) with
    its frequency; the REF's is what they leave. Nothing is decided when the
    site is not listed, when any ALT's frequency is unknown (it is not zero),
    when the frequencies sum past one, or when the two highest are within
    ``_TIE``. One ALT's low frequency does not make the REF major: another
    ALT can be.
    """
    if not alts or any(af is None for af in alts.values()):
        return None
    known = {alt: float(af) for alt, af in alts.items() if af is not None}
    total = sum(known.values())
    if total > 1.0 + _TIE:
        return None
    frequencies = {**known, ref: max(0.0, 1.0 - total)}
    ranked = sorted(frequencies.items(), key=lambda item: item[1], reverse=True)
    if len(ranked) > 1 and ranked[0][1] - ranked[1][1] <= _TIE:
        return None
    return ranked[0][0]

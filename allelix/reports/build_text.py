# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Wording for the build decision, shared by every command and report (ADR-0044).

``source_text`` says where the build came from; ``build_notices`` lists the
notices the decision carries (label mismatch, ``--build`` against tentative
evidence, format default, tie, hg19 chrM, GRCh36). The CLI prints them for
``analyze``, ``stats``, ``extract``, ``export plink`` and ``compare``; the HTML
report shows the mismatch and tie notices. One place, so a tie is never
described as "decided by known-SNP positions" on one surface and correctly on
another.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, NamedTuple

from allelix.utils.build_detect import (
    BUILD_GRCH36,
    BUILD_SOURCE_LABELS,
    SOURCE_ASSEMBLY_HEADER,
    SOURCE_CONTIG_LENGTHS,
    SOURCE_DEFAULT,
    SOURCE_FILE_LABEL,
    SOURCE_SEEDS,
)

if TYPE_CHECKING:
    from allelix.reports._pipeline import BuildDiagnostics

# Notice kinds, so a surface can pick the ones it shows.
MISMATCH = "mismatch"
DEFAULT = "default"
INCONCLUSIVE = "inconclusive"
CHRM = "chrm"
GRCH36 = "grch36"


class Notice(NamedTuple):
    """One notice: its kind (see the constants above) and its text."""

    kind: str
    text: str


def _join(builds: tuple[str, ...]) -> str:
    """``A and B`` / ``A, B and C``."""
    return builds[0] if len(builds) == 1 else ", ".join(builds[:-1]) + " and " + builds[-1]


def _tie(diag: BuildDiagnostics) -> str:
    return (
        f"tied between {_join(diag.tied_builds)} "
        f"({diag.tied_votes} each of {diag.inspected_count})"
    )


def _is_grch36_fail_safe_tie(diag: BuildDiagnostics) -> bool:
    return (
        diag.source == SOURCE_SEEDS
        and bool(diag.tied_builds)
        and diag.effective_build == BUILD_GRCH36
    )


def _seed_count(diag: BuildDiagnostics) -> str:
    if diag.tied_builds:
        return f"known-SNP positions {_tie(diag)}"
    if diag.detected_build:
        strength = "" if diag.position_confident else ", tentative"
        return (
            f"{diag.matched_count}/{diag.inspected_count} known-SNP positions match "
            f"{diag.detected_build}{strength}"
        )
    return f"{diag.inspected_count} known-SNP position(s) found"


def _default_evidence(diag: BuildDiagnostics) -> str:
    """The seed counts behind a format-default decision, in words."""
    if diag.tied_builds:
        return f"its known-SNP positions are {_tie(diag)}"
    return (
        f"{diag.matched_count} of {diag.inspected_count} known-SNP positions match "
        f"{diag.detected_build}"
    )


def short_source(diag: BuildDiagnostics) -> str:
    """Where the build came from, in a few words (HTML card, ``stats`` row)."""
    if diag.override:
        return BUILD_SOURCE_LABELS[diag.source]
    if _is_grch36_fail_safe_tie(diag):
        return "GRCh36 safety choice; known-SNP positions tied"
    if diag.tied_builds:
        return f"{BUILD_SOURCE_LABELS[diag.source]}; known-SNP positions tied"
    return BUILD_SOURCE_LABELS[diag.source]


def source_text(diag: BuildDiagnostics) -> str:
    """Where the build came from, with the evidence (the ``Build:`` banner)."""
    if diag.override:
        text = "forced with --build"
        if diag.file_label and diag.file_label != diag.effective_build:
            text += f"; the file says {diag.file_label}"
        return text
    if diag.source == SOURCE_CONTIG_LENGTHS:
        text = f"from the contig lengths: {diag.lengths_detail}"
        return text + (f"; known-SNP positions {_tie(diag)}" if diag.tied_builds else "")
    if _is_grch36_fail_safe_tie(diag):
        return f"GRCh36 safety choice: known-SNP positions {_tie(diag)}"
    if diag.source == SOURCE_SEEDS:
        strength = "" if diag.position_confident else ", tentative"
        return (
            f"from known-SNP positions: {diag.matched_count}/{diag.inspected_count} match "
            f"{diag.detected_build}{strength}"
        )
    if diag.source == SOURCE_ASSEMBLY_HEADER:
        return f"from the VCF assembly header; {_seed_count(diag)}"
    if diag.source == SOURCE_FILE_LABEL:
        return f"from the file's build line; {_seed_count(diag)}"
    return f"this format's default; the file states no build; {_seed_count(diag)}"


_TIE_REASON = {
    SOURCE_CONTIG_LENGTHS: "which the contig lengths decide",
    SOURCE_ASSEMBLY_HEADER: "which the VCF assembly header states",
    SOURCE_FILE_LABEL: "which the file's build line states",
    SOURCE_DEFAULT: "this format's default, as the file states no build",
}


# The builds a user can name, in the same words as the stop message.
_BUILD_CHOICES = (
    "--build grch37 or --build grch38 (for an old GRCh36 / hg18 file, --build grch36, "
    "see docs/grch36-liftover.md)"
)


def build_notices(diag: BuildDiagnostics) -> list[Notice]:
    """The notices this decision carries, in the order they are shown."""
    notices: list[Notice] = []
    if diag.mismatch and diag.tied_builds:
        notices.append(
            Notice(
                MISMATCH,
                f"Build mismatch: the file says {diag.file_label}, but its known-SNP "
                f"positions are {_tie(diag)}. Using {diag.effective_build} as the safe "
                "choice: it turns position lookups off rather than risk the wrong build's "
                "records. Pass --build to choose the build yourself.",
            )
        )
    elif diag.mismatch:
        strength = "" if diag.position_confident else "tentatively "
        notices.append(
            Notice(
                MISMATCH,
                f"Build mismatch: the file says {diag.file_label}, but its known-SNP "
                f"positions {strength}say {diag.detected_build} ({diag.matched_count}/"
                f"{diag.inspected_count}). Using {diag.effective_build}. The provider may "
                "have mislabeled the build; pass --build to choose it yourself.",
            )
        )
    elif diag.source == SOURCE_DEFAULT and not diag.override:
        notices.append(
            Notice(
                DEFAULT,
                f"Note: the file states no genome build, and {_default_evidence(diag)}: not "
                f"enough to be certain, so this format's default ({diag.effective_build}) "
                f"is used. If the file uses another build, pass {_BUILD_CHOICES}.",
            )
        )
    if diag.override and diag.effective_build != BUILD_GRCH36:
        # Confident contrary evidence stopped the run (ADR-0044);
        # tentative evidence against --build is shown, not obeyed.
        if diag.tied_builds and diag.effective_build not in diag.tied_builds:
            against = f"the file's known-SNP positions are {_tie(diag)}"
        elif diag.detected_build and diag.detected_build != diag.effective_build:
            against = (
                f"the file's known-SNP positions tentatively say {diag.detected_build} "
                f"({diag.matched_count}/{diag.inspected_count})"
            )
        else:
            against = ""
        if against:
            notices.append(
                Notice(
                    MISMATCH,
                    f"Build warning: --build {diag.effective_build.lower()} was given, "
                    f"but {against}. Using {diag.effective_build} as forced. If the file "
                    f"uses another build, pass {_BUILD_CHOICES}.",
                )
            )
    if diag.inconclusive:
        if _is_grch36_fail_safe_tie(diag):
            using = "Using GRCh36 as the safe choice: position lookups are off."
        else:
            reason = _TIE_REASON.get(diag.source, "")
            using = f"Using {diag.effective_build}, {reason}; it is one of the tied builds."
        notices.append(
            Notice(
                INCONCLUSIVE,
                f"Build evidence inconclusive: the known-SNP positions are {_tie(diag)}. "
                f"{using} If the file uses a different build, pass {_BUILD_CHOICES}.",
            )
        )
    if diag.mt_skipped and diag.mt_ucsc_name:
        notices.append(
            Notice(
                CHRM,
                f"Note: {diag.mt_skipped:,} mitochondrial record(s) named chrM or M were not "
                "annotated. On GRCh37 that name is UCSC hg19's chrM (16,571 bp), whose "
                "coordinates are not the reference mitochondrial sequence, unless the file "
                "declares it with the reference length (16,569 bp).",
            )
        )
    elif diag.mt_skipped:
        notices.append(
            Notice(
                CHRM,
                "Note: chrM is declared with the hg19 length (16,571 bp), whose coordinates "
                f"are not the reference mitochondrial sequence; {diag.mt_skipped:,} "
                "mitochondrial record(s) were not annotated.",
            )
        )
    if diag.effective_build == BUILD_GRCH36:
        if diag.override:
            how = "GRCh36 (hg18) was forced with --build."
        elif _is_grch36_fail_safe_tie(diag):
            how = (
                "GRCh36 (hg18) was chosen as the safe option because the known-SNP "
                "positions are tied."
            )
        elif diag.source == SOURCE_SEEDS and not diag.position_confident:
            how = (
                "GRCh36 (hg18) was chosen as the safe option: the known-SNP positions "
                f"lean to it ({diag.matched_count}/{diag.inspected_count}) but not "
                "enough to be certain."
            )
        elif diag.source == SOURCE_SEEDS:
            how = "GRCh36 (hg18) detected from the known-SNP positions."
        else:
            how = "GRCh36 (hg18): the file states this build."
        notices.append(
            Notice(
                GRCH36,
                f"Warning: {how} ClinPGx matches by rsID and still runs; GWAS Catalog risk "
                "alleles and single-base SNPedia genotypes are not assessed on GRCh36 "
                "(counted below). ClinVar is skipped (no GRCh36 cache, see ADR-0025), and "
                "gnomAD, AlphaMissense and CADD give no values (they need GRCh38 "
                "coordinates). For full coverage, lift the file over to GRCh38 first: "
                "docs/grch36-liftover.md",
            )
        )
    return notices

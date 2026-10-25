# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Report rendering: terminal, JSON, and HTML."""

from __future__ import annotations

import contextlib
import os
import stat
import sys
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterator
    from typing import BinaryIO, TextIO

    from allelix.parsers.base import ArrayInputStats
    from allelix.parsers.vcf import VcfQualityStats

# Single source of truth for the project's compliance / regulatory contract.
# Surfaced verbatim in JSON `regulatory_notice` and in the HTML's Technical
# details section; a two-sentence version is at the top of the HTML. See ADR-0003.
REGULATORY_NOTICE = (
    "This report is informational research output. It surfaces classifications "
    "made by external databases (ClinVar, ClinPGx, …) for variants present in "
    "the input genotype file. It is not medical advice and not a diagnosis. "
    "Every classification is attributed to its source database; Allelix does "
    "not independently classify variants."
)


def validate_output_path(input_path: Path, output_path: Path, *, legacy_temp: bool = True) -> None:
    """Refuse input aliases, including the former fixed temp when applicable.

    Missing paths can be compared by resolved identity. Real resolution/stat
    errors propagate rather than certifying an unknown path safe.
    """
    source = input_path.resolve()
    try:
        input_path.stat()
    except FileNotFoundError:
        source_exists = False
    else:
        source_exists = True
    candidates = [output_path]
    if legacy_temp:
        candidates.append(output_path.with_name(output_path.name + ".tmp"))
    for candidate in candidates:
        same = candidate.resolve() == source
        if not same:
            try:
                candidate.stat()
            except FileNotFoundError:
                pass
            else:
                if source_exists:
                    same = candidate.samefile(input_path)
        if same:
            raise ValueError(
                f"Output {output_path} aliases input {input_path} through {candidate}."
            )


class StagedOutput:
    """One invocation's allocated file and descriptor, published only after close."""

    def __init__(self, destination: Path) -> None:
        """Allocate a private same-directory temp and remember ordinary destination mode."""
        self.destination = destination
        try:
            existing = destination.stat()
        except FileNotFoundError:
            self.mode = None
        else:
            self.mode = existing.st_mode & 0o777 if stat.S_ISREG(existing.st_mode) else None
        descriptor, name = tempfile.mkstemp(
            dir=destination.parent, prefix="." + destination.name + ".", suffix=".tmp"
        )
        self.temporary = Path(name)
        self.published = False
        self.descriptor: int | None = descriptor
        self.stream: TextIO | BinaryIO | None = None

    def open_text(self, encoding: str | None = "utf-8", newline: str | None = None) -> TextIO:
        """Transfer the allocated descriptor to a text stream exactly once."""
        if self.descriptor is None:
            raise RuntimeError("Staged descriptor already consumed")
        stream = os.fdopen(self.descriptor, "w", encoding=encoding, newline=newline)
        self.stream = stream
        self.descriptor = None
        return stream

    def open_binary(self) -> BinaryIO:
        """Transfer the allocated descriptor to a binary stream exactly once."""
        if self.descriptor is None:
            raise RuntimeError("Staged descriptor already consumed")
        stream = os.fdopen(self.descriptor, "wb")
        self.stream = stream
        self.descriptor = None
        return stream

    def publish(self) -> None:
        """Replace the destination with complete, closed content, retaining its mode."""
        if self.published:
            raise RuntimeError("Staged output already published")
        if self.stream is None or not self.stream.closed:
            raise RuntimeError("Staged output must close successfully before publication")
        if self.mode is not None:
            self.temporary.chmod(self.mode)
        os.replace(self.temporary, self.destination)
        self.published = True

    def cleanup(self) -> None:
        """Close and remove only resources allocated by this invocation."""
        error: BaseException | None = None
        try:
            if self.stream is not None:
                self.stream.close()
            elif self.descriptor is not None:
                descriptor, self.descriptor = self.descriptor, None
                os.close(descriptor)
        except BaseException as exc:
            error = exc
        try:
            if not self.published:
                try:
                    self.temporary.unlink(missing_ok=True)
                except PermissionError:
                    if os.name != "nt":
                        raise
                    # Only this still-owned unpublished temp may need its
                    # Windows readonly attribute cleared after failed replace.
                    self.temporary.chmod(stat.S_IREAD | stat.S_IWRITE)
                    self.temporary.unlink(missing_ok=True)
        except BaseException as exc:
            if error is None:
                error = exc
        if error is not None:
            raise error


@contextlib.contextmanager
def staged_output(destination: Path) -> Iterator[StagedOutput]:
    """Own one staged output through failures and interrupts; never delete final files."""
    staged = StagedOutput(destination)
    try:
        yield staged
    finally:
        primary_error = sys.exc_info()[0] is not None
        try:
            staged.cleanup()
        except BaseException:
            if not primary_error:
                raise


def atomic_write_text(
    path: Path, content: str, encoding: str = "utf-8", *, input_path: Path | None = None
) -> None:
    """Publish complete text through a private owned temp in the target directory.

    Preserve an allowed ordinary destination's mode. Failures before publication
    keep prior output; cleanup never touches unrelated old fixed temp paths.
    Ownership, ACLs, xattrs, crash durability and concurrent path swaps are not
    guaranteed.
    """
    if input_path is not None:
        validate_output_path(input_path, path)
    with staged_output(path) as staged:
        with staged.open_text(encoding=encoding) as stream:
            stream.write(content)
        staged.publish()


def vcf_quality_lines(stats: VcfQualityStats | None) -> list[str]:
    """Plain-language lines describing what the VCF quality gate did.

    Shared by the terminal, HTML and CLI output so every surface says the
    same thing. Empty when the input was not a VCF or nothing was gated.
    """
    if stats is None:
        return []
    lines: list[str] = []
    if stats.hom_ref_calls and not stats.non_reference_calls:
        lines.append("This file has no non-reference calls.")
    if stats.filtered_rows:
        labels = ", ".join(
            f"{label} {count:,}" for label, count in sorted(stats.filtered_by_label.items())
        )
        if stats.kept_filtered:
            lines.append(
                f"{stats.filtered_rows:,} row(s) failed the variant caller's filters "
                f"({labels}) and were kept because of --keep-filtered."
            )
        else:
            lines.append(
                f"{stats.filtered_rows:,} row(s) excluded by the variant caller's filters "
                f"({labels}). Use --keep-filtered to include them."
            )
    if stats.low_quality_calls:
        lines.append(
            f"{stats.low_quality_calls:,} call(s) reported as no-calls: genotype quality "
            f"below {stats.min_gq} or no reads."
        )
    if stats.spanning_deletion_calls:
        lines.append(
            f"{stats.spanning_deletion_calls:,} call(s) with a spanning deletion (*) were "
            "not assessed; the JSON report lists their positions."
        )
    if stats.half_calls:
        lines.append(
            f"{stats.half_calls:,} half call(s) with one allele not called (such as ./1): "
            "ClinVar and the GWAS Catalog assessed the called allele; ClinPGx and "
            "SNPedia, which need both alleles, did not."
        )
    if stats.trimmed_matches:
        lines.append(
            f"{stats.trimmed_matches:,} call(s) written with untrimmed alleles matched "
            "ClinVar after trimming."
        )
    if stats.clinvar_ref_mismatches:
        lines.append(
            f"{stats.clinvar_ref_mismatches:,} record(s) at a ClinVar position have a "
            "REF that is not ClinVar's there (not matched; a VCF is forward strand, so a "
            "complemented REF is a mismatch). Check the file's build and coordinates."
        )
    if stats.reference_carried_calls:
        lines.append(
            f"{stats.reference_carried_calls:,} homozygous-reference call(s) at sites where "
            "a ClinVar record describes the allele the reference genome itself carries "
            "were assessed with the ploidy declared by the file."
        )
    if stats.reference_source_conflicts:
        lines.append(
            f"{stats.reference_source_conflicts:,} reference call(s) were not assessed "
            "because reference sources disagree at the same build and position. "
            "The reference databases do not establish one reference base; "
            "this is not an established input REF mismatch."
        )
    if stats.reference_ploidy_conflicts:
        lines.append(
            f"{stats.reference_ploidy_conflicts:,} reference call(s) were not assessed "
            "because matching input observations disagree about ploidy "
            "(haploid versus diploid) at the same position."
        )
    if stats.reference_carried_not_in_file:
        lines.append(
            f"{stats.reference_carried_not_in_file:,} site(s) where the reference genome "
            "carries a ClinVar variant allele have no record in this VCF. A person "
            "homozygous for that allele has no record there either, so these sites were "
            "not assessed; a gVCF or an all-sites VCF says which."
        )
    if stats.malformed_rows:
        lines.append(f"{stats.malformed_rows:,} malformed row(s) were skipped.")
    if stats.no_genotype_rows:
        lines.append(
            f"{stats.no_genotype_rows:,} row(s) had no genotype for this sample and were "
            "not assessed."
        )
    if stats.unusable_genotype_calls:
        lines.append(
            f"{stats.unusable_genotype_calls:,} call(s) had a genotype that names no allele "
            "of its record and were not assessed."
        )
    if stats.symbolic_allele_calls:
        lines.append(
            f"{stats.symbolic_allele_calls:,} call(s) named a symbolic allele (such as "
            "<NON_REF>) and were not assessed."
        )
    if stats.symbolic_rows:
        lines.append(
            f"{stats.symbolic_rows:,} structural-variant / copy-number record(s) were "
            "not assessed: Allelix annotates SNVs and indels."
        )
    if not stats.is_gvcf:
        lines.append(
            "This VCF has no reference blocks, so a position without a record may be "
            "reference or may not have been covered. Only explicit no-calls are listed."
        )
    return lines


def array_quality_lines(
    stats: ArrayInputStats | None, probe_classes: dict[str, int] | None
) -> list[str]:
    """Plain-language lines on the array calls ClinVar cannot match (ADR-0051).

    Shared by the terminal and HTML output. Empty for VCF input, or when
    there is nothing to say.
    """
    lines: list[str] = []
    i_probes = probe_classes.get("i_probes_at_plp_loci", 0) if probe_classes else 0
    if i_probes:
        lines.append(
            f"{i_probes:,} vendor probe(s) without an rsID (such as i5000123) sit where "
            "ClinVar classifies the variant as Pathogenic or Likely pathogenic. They are "
            "not matched to ClinVar: nothing confirms that a probe tests that variant."
        )
    if stats is not None and stats.di_calls:
        at_plp = probe_classes.get("di_calls_at_plp_rsids", 0) if probe_classes else 0
        lines.append(
            f"{stats.di_calls:,} call(s) are written as D or I (a deletion or an "
            "insertion), which Allelix cannot match to ClinVar"
            + (
                f"; {at_plp:,} of them at rsIDs ClinVar classifies as Pathogenic or "
                "Likely pathogenic."
                if at_plp
                else "."
            )
        )
    return lines


def input_left_unchecked(
    stats: VcfQualityStats | None,
    array_stats: ArrayInputStats | None,
    probe_classes: dict[str, int] | None,
) -> bool:
    """True when the input-quality lines report calls that were not checked.

    The counterpart of ``vcf_quality_lines`` and ``array_quality_lines``
    for the HTML summary line: rows excluded, no-calls, and calls not
    assessed by every source count. The notes that only describe the
    file (no reference blocks, kept filtered rows, trimmed matches,
    reference-carried calls assessed) do not.
    """
    if stats is not None and (
        (stats.filtered_rows and not stats.kept_filtered)
        or stats.low_quality_calls
        or stats.spanning_deletion_calls
        or stats.half_calls
        or stats.clinvar_ref_mismatches
        or stats.reference_source_conflicts
        or stats.reference_ploidy_conflicts
        or stats.reference_carried_not_in_file
        or stats.malformed_rows
        or stats.no_genotype_rows
        or stats.unusable_genotype_calls
        or stats.symbolic_allele_calls
        or stats.symbolic_rows
    ):
        return True
    if probe_classes and probe_classes.get("i_probes_at_plp_loci", 0):
        return True
    return bool(array_stats is not None and array_stats.di_calls)


# How each declared input method is described to the user.
_METHOD_LABELS: dict[str, str] = {
    "low-pass-wgs-imputed": (
        "low-pass whole genome sequencing (about 2x coverage); genotypes are imputed"
    ),
    "array-imputed": (
        "a genotyping array extended by imputation; most genotypes are imputed, not measured"
    ),
    "vcf-imputed": (
        "imputed: the VCF header declares imputation fields (Minimac, Beagle or GLIMPSE), so most "
        "genotypes are imputed, not measured"
    ),
}


def input_method_label(method: str) -> str:
    """Plain-language provenance for a declared input method, or "" if none."""
    if not method:
        return ""
    return _METHOD_LABELS.get(method, method)


def is_imputed(method: str) -> bool:
    """True when the input's genotypes are imputed rather than directly called."""
    return method.endswith("-imputed")

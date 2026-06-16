# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Utility subcommands: stats, extract, compare, export."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

import click
from rich.table import Table

from allelix.cli import _helpers, main
from allelix.cli._helpers import console
from allelix.cli._options import (
    _BUILD_OPT,
    _DATA_DIR_OPT,
    _FILE_ARG,
    _FORMAT_OPT,
    _SAMPLE_OPT,
)
from allelix.databases import resolve_data_dir

if TYPE_CHECKING:
    from allelix.models import Variant
    from allelix.parsers.base import GenotypeParser


@main.command()
@_FILE_ARG
@_FORMAT_OPT
@_SAMPLE_OPT
def stats(file_path: Path, fmt: str | None, sample: str | None) -> None:
    """Show summary statistics for a genotype file."""
    from allelix.reports.high_value import format_warnings, load_high_value_snps, scan_no_calls

    parser = _helpers._resolve_parser(file_path, fmt, sample=sample)
    counter, stderr_handler, snapshot = _helpers._wire_parser_logging()

    high_value = load_high_value_snps()
    hv_rsids = set(high_value)
    hv_variants: list[Variant] = []

    total = 0
    no_calls = 0
    het = 0
    hom = 0
    chrom_counts: dict[str, int] = {}
    try:
        metadata = parser.get_metadata(file_path)
        for variant in parser.parse(file_path):
            total += 1
            if variant.rsid in hv_rsids:
                hv_variants.append(variant)
            if variant.is_no_call:
                no_calls += 1
            elif variant.is_heterozygous:
                het += 1
            else:
                hom += 1
            chrom_counts[variant.chromosome] = chrom_counts.get(variant.chromosome, 0) + 1
    finally:
        _helpers._unwire_parser_logging(counter, stderr_handler, snapshot)

    summary = Table(title=f"Genotype File Stats: {file_path.name}")
    summary.add_column("Metric", style="cyan", no_wrap=True)
    summary.add_column("Value", style="green")
    summary.add_row("Format", parser.display_name)
    summary.add_row("Sample ID", metadata["sample_id"] or "(unknown)")
    summary.add_row("Build", metadata["build"])
    summary.add_row("Total SNPs", f"{total:,}")
    summary.add_row("No-calls", f"{no_calls:,} ({_helpers._percent(no_calls, total)})")
    summary.add_row("Heterozygous", f"{het:,} ({_helpers._percent(het, total)})")
    summary.add_row("Homozygous", f"{hom:,} ({_helpers._percent(hom, total)})")
    if counter.count:
        summary.add_row(
            "Skipped (malformed)",
            f"[yellow]{counter.count:,}[/yellow] (see warnings on stderr)",
        )

    hv_warnings = scan_no_calls(hv_variants, high_value)
    if hv_warnings:
        summary.add_row(
            "High-value no-calls",
            f"[red]{len(hv_warnings)}[/red]",
        )
    console.print(summary)

    if hv_warnings:
        for line in format_warnings(hv_warnings):
            console.print(f"  [red]⚠[/red] {line}")

    chrom_table = Table(title="Variants per Chromosome")
    chrom_table.add_column("Chromosome", style="cyan", no_wrap=True)
    chrom_table.add_column("Count", style="green", justify="right")
    for chrom in sorted(chrom_counts, key=_helpers._chrom_sort_key):
        chrom_table.add_row(chrom, f"{chrom_counts[chrom]:,}")
    console.print(chrom_table)


@main.command()
@_FILE_ARG
@_FORMAT_OPT
@click.option(
    "--snps",
    required=True,
    help="Comma-separated rsIDs to extract (e.g., rs1801133,rs4680).",
)
@_DATA_DIR_OPT
@_SAMPLE_OPT
def extract(
    file_path: Path,
    fmt: str | None,
    snps: str,
    data_dir: Path | None,
    sample: str | None,
) -> None:
    """Print diploid genotypes for specific rsIDs — spot-check carrier status.

    Useful for verifying ClinVar / ClinPGx hits against the actual file
    before trusting them. The "Genotype" column shows the diploid call as
    the array (or VCF) reported it; "Het?" and "No-call?" answer the
    questions the carrier rule (ADR-0007) actually checks.

    For VCF input, an indexed fast-path is used automatically when the
    file is gzipped, a ``.tbi`` index is present alongside, ``pysam`` is
    installed (``pip install allelix[vcf-index]``), and gnomAD is set up
    to resolve rsIDs to coordinates. Falls back to a sequential scan
    otherwise — correct for any VCF, just O(file size).
    """
    parser = _helpers._resolve_parser(file_path, fmt, sample=sample)
    wanted = {s.strip() for s in snps.split(",") if s.strip()}
    if not wanted:
        raise click.ClickException("--snps cannot be empty.")

    found = _try_tabix_extract(file_path, wanted, sample, data_dir)
    if found is None:
        # Sequential fallback — works for any parser.
        found = _sequential_extract(parser, file_path, wanted)
        _maybe_print_vcf_index_hint(file_path)

    _render_extract_table(file_path, wanted, found)


def _sequential_extract(
    parser: GenotypeParser, file_path: Path, wanted: set[str]
) -> dict[str, Variant]:
    """Stream the file once; collect Variants whose rsID is in ``wanted``.

    Streaming early-exit once all wanted rsIDs are found — for arrays
    that yield <1 M variants this is fast enough on any input.
    """
    counter, stderr_handler, snapshot = _helpers._wire_parser_logging()
    found: dict[str, Variant] = {}
    try:
        for variant in parser.parse(file_path):
            if variant.rsid in wanted:
                found[variant.rsid] = variant
                if len(found) == len(wanted):
                    break
    finally:
        _helpers._unwire_parser_logging(counter, stderr_handler, snapshot)
    return found


def _try_tabix_extract(
    file_path: Path,
    wanted: set[str],
    sample: str | None,
    data_dir: Path | None,
) -> dict[str, Variant] | None:
    """Fast-path VCF lookup via pysam tabix. Returns None if not viable.

    Required: ``.vcf.gz`` input, a ``.tbi`` index alongside, pysam
    importable, and gnomAD ready (used to resolve rsID → genomic
    coordinate for the tabix query).

    The function never raises on a check failure — every "not viable"
    path returns None so the caller falls back to sequential scan.
    """
    # Check 1: gzipped VCF with a tabix index alongside.
    if "".join(file_path.suffixes[-2:]).lower() != ".vcf.gz":
        return None
    tbi_path = file_path.with_suffix(file_path.suffix + ".tbi")
    if not tbi_path.exists():
        return None

    # Check 2: pysam importable.
    try:
        import pysam  # pragma: no cover
    except ImportError:
        return None

    # The remaining body requires pysam + gnomAD ready, neither of which
    # the test environment has. Exercised in integration / manual
    # validation on the droplet; unit tests cover the not-viable
    # return-None paths above.
    return _execute_tabix_extract(  # pragma: no cover
        file_path, wanted, sample, data_dir, pysam
    )


def _execute_tabix_extract(  # pragma: no cover
    file_path: Path,
    wanted: set[str],
    sample: str | None,
    data_dir: Path | None,
    pysam: Any,  # noqa: ANN401 — pysam ships no type stubs
) -> dict[str, Variant] | None:
    """Real tabix execution.

    Separated from ``_try_tabix_extract`` so the pre-flight checks can
    be unit-tested without pysam installed. Returns None on any
    in-execution failure (gnomAD not ready, invalid sample binding,
    etc.) so the caller falls back to sequential scan.
    """
    from allelix.annotators.gnomad import GnomadAnnotator
    from allelix.parsers.vcf import _open_vcf, _parse_data_line, _read_header

    resolved_data_dir = resolve_data_dir(data_dir)
    gnomad = GnomadAnnotator(resolved_data_dir)
    try:
        if not gnomad.is_ready():
            return None
        coord_map = gnomad.bulk_resolve_coordinates(wanted)
    finally:
        gnomad.close()
    if not coord_map:
        return None

    found: dict[str, Variant] = {}
    with _open_vcf(file_path) as handle:
        header = _read_header(handle)
    if not header.samples:
        sample_idx = -1
    elif sample is not None:
        if sample not in header.samples:
            return None
        sample_idx = header.samples.index(sample)
    elif len(header.samples) > 1:
        return None
    else:
        sample_idx = 0

    tabix_file = pysam.TabixFile(str(file_path))
    try:
        for rsid, coords in coord_map.items():
            if rsid in found or not coords:
                continue
            for chrom, pos, _ref, _alt in coords:
                # gnomAD chromosomes are bare ('1', 'X'); the VCF might
                # use 'chr1'. Try both.
                for chrom_candidate in (chrom, f"chr{chrom}"):
                    try:
                        rows = list(tabix_file.fetch(chrom_candidate, pos - 1, pos))
                    except (ValueError, OSError):
                        continue
                    for row in rows:
                        variant = _parse_data_line(row, sample_idx)
                        if variant is None:
                            continue
                        if variant.rsid == rsid:
                            # _parse_data_line already normalizes chromosome;
                            # no re-normalization needed.
                            found[rsid] = variant
                            break
                    if rsid in found:
                        break
                if rsid in found:
                    break
    finally:
        tabix_file.close()
    return found


def _maybe_print_vcf_index_hint(file_path: Path) -> None:
    """One-line stderr hint when a VCF fell back to sequential scan.

    Only printed when:
      - input file is gzipped VCF (likely a large clinical export)
      - pysam is NOT importable (so installing the extras is actionable)

    Skips arrays, .vcf (often small test files), and the case where
    pysam IS importable (the user already has the dep — index missing
    or gnomAD absent is on them).
    """
    if "".join(file_path.suffixes[-2:]).lower() != ".vcf.gz":
        return
    try:
        import pysam  # noqa: F401
    except ImportError:
        console.print(
            "[dim]Tip: pip install allelix[vcf-index] for fast indexed "
            "extraction on large VCFs.[/dim]"
        )


def _render_extract_table(file_path: Path, wanted: set[str], found: dict[str, Variant]) -> None:
    """Render the extract result as a Rich table."""
    table = Table(title=f"Genotypes from {file_path.name}")
    table.add_column("rsID", style="cyan", no_wrap=True)
    table.add_column("Chr", no_wrap=True)
    table.add_column("Position", justify="right")
    table.add_column("Genotype", style="yellow", no_wrap=True)
    table.add_column("Het?", justify="center")
    table.add_column("No-call?", justify="center")
    for rsid in sorted(wanted):
        variant = found.get(rsid)
        if variant is None:
            table.add_row(rsid, "—", "—", "[red]not in file[/red]", "—", "—")
            continue
        table.add_row(
            variant.rsid,
            variant.chromosome,
            f"{variant.position:,}",
            variant.genotype,
            "yes" if variant.is_heterozygous else "no",
            "[red]yes[/red]" if variant.is_no_call else "no",
        )
    console.print(table)


@main.command()
@click.argument("file1", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.argument("file2", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option("--format1", "fmt1", default=None, help="Force parser for file 1.")
@click.option("--format2", "fmt2", default=None, help="Force parser for file 2.")
@click.option("--sample1", default=None, help="VCF sample column for file 1 (multi-sample VCFs).")
@click.option("--sample2", default=None, help="VCF sample column for file 2 (multi-sample VCFs).")
def compare(
    file1: Path,
    file2: Path,
    fmt1: str | None,
    fmt2: str | None,
    sample1: str | None,
    sample2: str | None,
) -> None:
    """Compare two genotype files — coverage overlap and concordance.

    Reports shared rsIDs, file-specific rsIDs, genotype agreement,
    strand-flip matches (complementary alleles on opposite strands),
    discordant calls, and strand-ambiguous positions.
    """
    from allelix.compare import compare_variants
    from allelix.utils.build_detect import detect_build

    parser1 = _helpers._resolve_parser(file1, fmt1, sample=sample1)
    parser2 = _helpers._resolve_parser(file2, fmt2, sample=sample2)
    variants1 = list(parser1.parse(file1))
    variants2 = list(parser2.parse(file2))

    det1 = detect_build(variants1)
    det2 = detect_build(variants2)
    build1 = det1.build or parser1.get_metadata(file1).get("build", "unknown")
    build2 = det2.build or parser2.get_metadata(file2).get("build", "unknown")

    result = compare_variants(variants1, variants2, build1=build1, build2=build2)

    if result.build1 != result.build2:
        console.print(
            f"[yellow]Warning: builds differ ({result.build1} vs {result.build2}). "
            "Position-based comparisons may be unreliable.[/yellow]"
        )

    table = Table(title="Coverage Summary")
    table.add_column("Metric", style="bold")
    table.add_column("Value", justify="right")
    table.add_row("File 1", f"{file1.name} ({result.file1_total:,} variants)")
    table.add_row("File 2", f"{file2.name} ({result.file2_total:,} variants)")
    table.add_row("Build (file 1)", result.build1)
    table.add_row("Build (file 2)", result.build2)
    table.add_row("Shared rsIDs", f"{result.shared:,}")
    table.add_row("File 1 only", f"{result.file1_only:,}")
    table.add_row("File 2 only", f"{result.file2_only:,}")
    console.print(table)

    conc_table = Table(title="Genotype Concordance")
    conc_table.add_column("Category", style="bold")
    conc_table.add_column("Count", justify="right")
    conc_table.add_column("%", justify="right")
    for label, count in [
        ("Concordant", result.concordant),
        ("Strand-flip match", result.strand_flip_match),
        ("Discordant", result.discordant),
        ("Strand-ambiguous", result.strand_ambiguous),
        ("No-call (either file)", result.no_call),
    ]:
        pct = _helpers._percent(count, result.shared) if result.shared else "—"
        conc_table.add_row(label, f"{count:,}", pct)
    console.print(conc_table)

    if result.chromosome_counts:
        chrom_table = Table(title="Per-Chromosome Breakdown")
        chrom_table.add_column("Chr", style="cyan", no_wrap=True)
        chrom_table.add_column("Concordant", justify="right")
        chrom_table.add_column("Flip", justify="right")
        chrom_table.add_column("Discordant", justify="right")
        chrom_table.add_column("Ambiguous", justify="right")
        chrom_table.add_column("No-call", justify="right")
        for chrom in sorted(result.chromosome_counts, key=_helpers._chrom_sort_key):
            c = result.chromosome_counts[chrom]
            chrom_table.add_row(
                chrom,
                str(c.get("concordant", 0)),
                str(c.get("strand_flip_match", 0)),
                str(c.get("discordant", 0)),
                str(c.get("strand_ambiguous", 0)),
                str(c.get("no_call", 0)),
            )
        console.print(chrom_table)


@main.group()
def export() -> None:
    """Export parsed genotype data to other formats."""


@export.command("plink")
@_FILE_ARG
@click.option(
    "--output-prefix",
    "-o",
    type=click.Path(path_type=Path),
    default=None,
    help="Base path for .bed/.bim/.fam (default: input stem).",
)
@_FORMAT_OPT
@_BUILD_OPT
@_DATA_DIR_OPT
@_SAMPLE_OPT
def export_plink_cmd(
    file_path: Path,
    output_prefix: Path | None,
    fmt: str | None,
    build: str,
    data_dir: Path | None,
    sample: str | None,
) -> None:
    """Convert to PLINK1 binary format (.bed/.bim/.fam).

    Produces a single-sample, SNP-major .bed file suitable for downstream
    tools (plink2 PCA, ADMIXTURE, PRSice). Uses gnomAD ref/alt for allele
    coding when available; falls back to monomorphic (A2=0) for positions
    without gnomAD coverage.
    """
    from allelix.exporters.plink import export_plink, resolve_ref_alt_via_gnomad

    parser = _helpers._resolve_parser(file_path, fmt, sample=sample)
    prefix = output_prefix if output_prefix else file_path.with_suffix("")
    build_override = _helpers._normalize_cli_build(build)
    metadata = parser.get_metadata(file_path)
    effective_build = build_override or metadata.get("build", "GRCh37")
    resolved = resolve_data_dir(data_dir)

    variants = list(parser.parse(file_path))

    # Sort by chromosome then position so the .bim has contiguous
    # chromosome blocks — PLINK1.9 rejects split chromosomes.
    chrom_order = {str(i): i for i in range(1, 23)}
    chrom_order.update({"X": 23, "Y": 24, "XY": 25, "MT": 26})
    variants.sort(
        key=lambda v: (chrom_order.get(v.chromosome, 99), v.chromosome, v.position),
    )

    variant_by_rsid: dict[str, Variant] = {v.rsid: v for v in variants if not v.is_no_call}

    try:
        ref_alt_map = resolve_ref_alt_via_gnomad(variant_by_rsid, resolved)
    except Exception:
        console.print(
            "[yellow]gnomAD coordinate resolution failed; using fallback allele coding.[/yellow]"
        )
        ref_alt_map = {}

    written, skipped, indel_skip, mono = export_plink(
        iter(variants), prefix, effective_build, ref_alt_map or None
    )
    skip_parts = []
    if skipped:
        skip_parts.append(f"{skipped:,} no-calls")
    if indel_skip:
        skip_parts.append(f"{indel_skip:,} indels")
    skip_msg = f" ({', '.join(skip_parts)} skipped)" if skip_parts else ""
    console.print(f"Wrote {written:,} variants to {prefix}.bed/.bim/.fam{skip_msg}")
    if mono > 0:
        pct = mono / written * 100 if written else 0
        console.print(
            f"[dim]{mono:,} markers ({pct:.0f}%) exported as monomorphic "
            f"(A2=0, ref/alt unknown or ambiguous).[/dim]"
        )
    if not ref_alt_map:
        console.print(
            "[yellow]gnomAD not available — all homozygous markers exported "
            "as monomorphic.[/yellow]"
        )
        console.print("[yellow]Run `allelix db update` first for proper allele coding.[/yellow]")
    console.print(
        "[dim]Single-sample export. Merging with other samples requires "
        "allele harmonization (--merge-mode or set-all-var-ids).[/dim]"
    )

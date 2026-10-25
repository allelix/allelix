# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Utility subcommands: stats, extract, compare, export."""

from __future__ import annotations

from dataclasses import replace
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
from allelix.parsers._helpers import normalize_chromosome
from allelix.reports._pipeline import BuildDiagnostics, UnusableInputError
from allelix.reports.build_text import short_source
from allelix.utils.build_detect import (
    BUILD_GRCH38,
    SOURCE_CONTIG_LENGTHS,
    SOURCE_OVERRIDE,
    BuildDecisionError,
    check_late_seeds,
    detect_build,
    needs_late_seeds,
)

if TYPE_CHECKING:
    from collections.abc import Callable

    from allelix.exporters.plink import RecordKey
    from allelix.models import Variant
    from allelix.parsers.base import GenotypeParser
    from allelix.utils.build_detect import BuildDecision


@main.command()
@_FILE_ARG
@_FORMAT_OPT
@_SAMPLE_OPT
@_DATA_DIR_OPT
@_BUILD_OPT
def stats(
    file_path: Path, fmt: str | None, sample: str | None, data_dir: Path | None, build: str
) -> None:
    """Show summary statistics for a genotype file.

    High-value no-call warnings cover rsID-less VCF records too: calls are
    resolved to rsIDs with the installed databases first.

    VCF genotype counts include reference-only records, before filtering.
    """
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
        with _helpers._resolved_variants(
            parser,
            file_path,
            data_dir=data_dir,
            build_override=_helpers._normalize_cli_build(build),
            high_value=hv_rsids,
        ) as resolved:
            for variant in resolved.variants:
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
    decision = resolved.decision
    diag = resolved.build_diagnostics()
    # The same decision as ``analyze`` (ADR-0044), with where it came from.
    summary.add_row("Build", f"{diag.effective_build} ({short_source(diag)})")
    summary.add_row("Total SNPs", f"{total:,}")
    summary.add_row("No-calls", f"{no_calls:,} ({_helpers._percent(no_calls, total)})")
    summary.add_row("Heterozygous", f"{het:,} ({_helpers._percent(het, total)})")
    summary.add_row("Homozygous", f"{hom:,} ({_helpers._percent(hom, total)})")
    if counter.count:
        summary.add_row(
            "Skipped (malformed)",
            f"[yellow]{counter.count:,}[/yellow] (see warnings on stderr)",
        )

    hv_warnings = scan_no_calls(
        hv_variants, high_value, resolved.high_value_states(decision.build)
    )
    if hv_warnings:
        summary.add_row(
            "High-value no-calls",
            f"[red]{len(hv_warnings)}[/red]",
        )
    console.print(summary)
    quality = getattr(parser, "quality_stats", None)
    if quality is not None:
        counts = Table(title="VCF Genotype Records")
        counts.add_column("Metric", style="cyan")
        counts.add_column("Count", style="green", justify="right")
        counts.add_row("Genotype records", f"{quality.genotype_records:,}")
        counts.add_row("Homozygous-reference", f"{quality.hom_ref_calls:,}")
        counts.add_row("Non-reference", f"{quality.non_reference_calls:,}")
        console.print(counts)
        if quality.hom_ref_calls and not quality.non_reference_calls:
            console.print("[yellow]Note:[/yellow] This file has no non-reference calls.")
    _helpers._emit_build_notices(diag, banner=False)

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
@_BUILD_OPT
def extract(
    file_path: Path,
    fmt: str | None,
    snps: str,
    data_dir: Path | None,
    sample: str | None,
    build: str,
) -> None:
    """Print diploid genotypes for specific rsIDs — spot-check carrier status.

    Useful for verifying ClinVar / ClinPGx hits against the actual file
    before trusting them. The "Genotype" column shows the diploid call as
    the array (or VCF) reported it; "Het?" and "No-call?" answer the
    questions the carrier rule (ADR-0007) actually checks.

    For VCF input, an indexed fast-path is used automatically when the
    file is gzipped, a ``.tbi`` index is present alongside, ``pysam`` is
    installed (``pip install allelix[vcf-index]``), and gnomAD is set up
    to resolve rsIDs to coordinates. Requests still unfound after indexed
    lookup share one sequential pass, including mapped coordinate misses.
    Unavailable indexes fall back to the normal sequential scan.
    """
    parser = _helpers._resolve_parser(file_path, fmt, sample=sample)
    wanted = {s.strip() for s in snps.split(",") if s.strip()}
    if not wanted:
        raise click.ClickException("--snps cannot be empty.")

    build_override = _helpers._normalize_cli_build(build)
    try:
        found = _try_tabix_extract(file_path, wanted, sample, data_dir, build_override)
    except BuildDecisionError as exc:
        raise click.ClickException(str(exc)) from exc
    if found is None:
        # Sequential fallback — works for any parser.
        found = _sequential_extract(parser, file_path, wanted, data_dir, build_override)
        _maybe_print_vcf_index_hint(file_path)

    _render_extract_table(file_path, wanted, found)


def _sequential_extract(
    parser: GenotypeParser,
    file_path: Path,
    wanted: set[str],
    data_dir: Path | None = None,
    build_override: str | None = None,
) -> dict[str, Variant]:
    """Stream the file once; collect Variants whose rsID is in ``wanted``.

    Streaming early-exit once all wanted rsIDs are found — for arrays
    that yield <1 M variants this is fast enough on any input. No early
    exit when contig lengths decided the build: the known-SNP check at
    the end of the stream must run before anything is printed (ADR-0044).

    Variant callers emit ``ID=.`` or positional pseudo-IDs
    (DeepVariant, DRAGEN, the GIAB GRCh38 benchmark). The stream goes
    through the shared rsID-less resolution step (``_resolved_variants``:
    build detection, then ClinVar / identity table / gnomAD), so a wanted
    rsID matches whether the file named it or a resolver recovered it.
    """
    counter, stderr_handler, snapshot = _helpers._wire_parser_logging()
    found: dict[str, Variant] = {}
    try:
        with _helpers._resolved_variants(
            parser, file_path, data_dir=data_dir, build_override=build_override
        ) as resolved:
            for variant in resolved.variants:
                if variant.rsid in wanted and variant.rsid not in found:
                    found[variant.rsid] = variant
                    if len(found) == len(wanted) and not resolved.must_exhaust:
                        break
            _helpers._emit_build_notices(resolved.build_diagnostics())
    finally:
        _helpers._unwire_parser_logging(counter, stderr_handler, snapshot)
    return found


def _try_tabix_extract(
    file_path: Path,
    wanted: set[str],
    sample: str | None,
    data_dir: Path | None,
    build_override: str | None = None,
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
    # the test environment has. The tabix fetch was validated manually
    # with pysam installed and a gnomAD cache present; unit tests cover
    # the not-viable return-None paths above.
    return _execute_tabix_extract(  # pragma: no cover
        file_path, wanted, sample, data_dir, pysam, build_override
    )


def _execute_tabix_extract(
    file_path: Path,
    wanted: set[str],
    sample: str | None,
    data_dir: Path | None,
    pysam: Any,  # noqa: ANN401 — pysam ships no type stubs
    build_override: str | None = None,
) -> dict[str, Variant] | None:
    """Real tabix execution.

    Separated from ``_try_tabix_extract`` so the pre-flight checks can
    be unit-tested without pysam installed. Returns None on any
    in-execution failure (gnomAD not ready, invalid sample binding,
    etc.) so the caller falls back to sequential scan.

    gnomAD gives each wanted rsID's GRCh38 coordinates, which is where
    the rows are fetched, so the fast path runs only on a VCF whose build
    is GRCh38: forced with ``--build``, or decided by its contig lengths
    (``_tabix_fast_path_allowed``, ADR-0044). Any other VCF, including one
    labelled ``assembly=GRCh38`` without contig lengths and one whose header
    evidence stops the run, falls back to the sequential scan, which applies
    the full decision. When the lengths decided, the rows at the known-SNP
    positions are fetched and checked against them; a contradiction
    raises ``BuildDecisionError``. The fetched rows then go through the shared
    rsID-less resolution step: a row matches a wanted rsID
    when the file names it or a resolver recovers it at that locus. The
    same rule applies here as in every other command, rather than a
    separate (ref, alt) comparison.
    """
    from allelix.annotators._resolution import resolve_batch
    from allelix.annotators.clinvar import ClinVarAnnotator
    from allelix.annotators.gnomad import GnomadAnnotator
    from allelix.annotators.rsid_identity import RsidIdentityAnnotator
    from allelix.parsers.vcf import VcfParser, _open_vcf, _parse_data_line, _read_header
    from allelix.utils.build_detect import (
        KIND_VCF,
        NO_SEEDS,
        check_late_seeds,
        contig_length_verdict,
        decide_build,
        known_seed_loci,
    )

    with _open_vcf(file_path) as handle:
        header = _read_header(handle)
    try:
        decision = decide_build(
            file_name=file_path.name,
            kind=KIND_VCF,
            override=build_override,
            file_label=header.build,
            default=None,
            lengths=contig_length_verdict(header.contig_lengths),
            seeds=NO_SEEDS,
        )
    except BuildDecisionError:
        return None
    if not _tabix_fast_path_allowed(decision):
        return None
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

    resolved_data_dir = resolve_data_dir(data_dir)
    clinvar = ClinVarAnnotator(resolved_data_dir)
    gnomad = GnomadAnnotator(resolved_data_dir)
    identity = RsidIdentityAnnotator(resolved_data_dir)
    found: dict[str, Variant] = {}
    try:
        if not gnomad.is_ready():
            return None
        coord_map = gnomad.bulk_resolve_coordinates(wanted)
        if not coord_map:
            return None
        tabix_file = pysam.TabixFile(str(file_path))
        try:
            aliases = _tabix_contig_aliases(tabix_file)
            row_cache: dict[tuple[str, int], list[Variant]] = {}
            seed_rows: list[Variant] = []
            for chrom, pos in known_seed_loci():
                seed_rows.extend(
                    _fetch_rows(
                        tabix_file,
                        chrom,
                        pos,
                        sample_idx,
                        _parse_data_line,
                        contig_aliases=aliases,
                        row_cache=row_cache,
                    )
                )
            late = detect_build(seed_rows)
            check_late_seeds(decision, late, file_path.name)
            for rsid, coords in coord_map.items():
                for chrom, pos, _ref, _alt in coords:
                    hit = _fetch_resolved(
                        tabix_file,
                        chrom,
                        pos,
                        rsid,
                        sample_idx,
                        lambda batch: resolve_batch(
                            batch, clinvar=clinvar, gnomad=gnomad, identity=identity
                        ),
                        parse_row=_parse_data_line,
                        build=BUILD_GRCH38,
                        contig_aliases=aliases,
                        row_cache=row_cache,
                    )
                    if hit is not None:
                        found[rsid] = hit
                        break
        finally:
            tabix_file.close()
    finally:
        for source in (clinvar, gnomad, identity):
            source.close()
    pending = wanted - found.keys()
    if pending:
        # One pass for every unresolved request, including mapped fetch misses.
        # The pending set excludes indexed no-calls and all other accepted hits.
        found.update(
            _sequential_extract(
                VcfParser(sample=sample), file_path, pending, resolved_data_dir, build_override
            )
        )
    else:
        _helpers._emit_build_notices(BuildDiagnostics.from_decision(decision, late))
    return found


def _tabix_fast_path_allowed(decision: BuildDecision) -> bool:
    """Whether the indexed extract may fetch rows at GRCh38 coordinates.

    Only for a GRCh38 build that came from ``--build`` or from the VCF's
    contig lengths. The fast path reads no records before fetching, so a
    build that rests on the ``assembly=`` label alone (or on seeds, which
    it never collects) is not enough: such a file goes through the
    sequential scan, which applies the full decision (ADR-0044). So does a
    file whose chrM has the hg19 length: the scan drops its MT records and
    reports them with the chrM note.
    """
    return (
        decision.build == BUILD_GRCH38
        and decision.source in (SOURCE_OVERRIDE, SOURCE_CONTIG_LENGTHS)
        and not decision.skip_mt
    )


def _fetch_resolved(
    tabix_file: Any,  # noqa: ANN401 — pysam.TabixFile, no type stubs
    chrom: str,
    pos: int,
    rsid: str,
    sample_idx: int,
    resolve: Callable[[list[Variant]], None],
    *,
    parse_row: Callable[[str, int], Variant | None],
    build: str,
    contig_aliases: dict[str, tuple[str, ...]] | None = None,
    row_cache: dict[tuple[str, int], list[Variant]] | None = None,
) -> Variant | None:
    """The row at ``chrom:pos`` that is ``rsid`` after resolution, or None.

    Contig aliases come from the index in its own order, using existing
    chromosome normalization. Fetched rows take ``build`` (the
    coordinates' build) before resolution.
    """
    variants = [
        replace(v, build=build)
        for v in _fetch_rows(
            tabix_file,
            chrom,
            pos,
            sample_idx,
            parse_row,
            contig_aliases=contig_aliases,
            row_cache=row_cache,
        )
    ]
    resolve(variants)
    for variant in variants:
        if variant.rsid == rsid:
            return variant
    return None


def _tabix_contig_aliases(tabix_file: Any) -> dict[str, tuple[str, ...]]:  # noqa: ANN401
    """Index-order spellings grouped by the existing chromosome normalization."""
    names: dict[str, list[str]] = {}
    for name in dict.fromkeys(tabix_file.contigs):
        names.setdefault(normalize_chromosome(name), []).append(name)
    return {chrom: tuple(actual) for chrom, actual in names.items()}


def _fetch_rows(
    tabix_file: Any,  # noqa: ANN401 — pysam.TabixFile, no type stubs
    chrom: str,
    pos: int,
    sample_idx: int,
    parse_row: Callable[[str, int], Variant | None],
    *,
    contig_aliases: dict[str, tuple[str, ...]] | None = None,
    row_cache: dict[tuple[str, int], list[Variant]] | None = None,
) -> list[Variant]:
    """Parsed rows at an exact position, using actual index contigs in file order.

    The optional cache belongs to one extraction and one sample binding.
    """
    canonical = normalize_chromosome(chrom)
    key = (canonical, pos)
    if row_cache is not None and key in row_cache:
        return row_cache[key]
    if contig_aliases is None:
        contig_aliases = _tabix_contig_aliases(tabix_file)
    found: list[Variant] = []
    for chrom_candidate in contig_aliases.get(canonical, ()):
        try:
            rows = list(tabix_file.fetch(chrom_candidate, pos - 1, pos))
        except (ValueError, OSError):
            continue
        found.extend(
            v
            for v in (parse_row(row, sample_idx) for row in rows)
            if v is not None and v.position == pos
        )
    if row_cache is not None:
        row_cache[key] = found
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
@click.option(
    "--format1", "fmt1", default=None, help="Force the reader for file 1, at your own risk."
)
@click.option(
    "--format2", "fmt2", default=None, help="Force the reader for file 2, at your own risk."
)
@click.option("--sample1", default=None, help="VCF sample column for file 1 (multi-sample VCFs).")
@click.option("--sample2", default=None, help="VCF sample column for file 2 (multi-sample VCFs).")
@click.option(
    "--build1",
    type=click.Choice(["grch37", "grch38", "grch36", "auto"], case_sensitive=False),
    default="auto",
    help=(
        "Genome build of file 1: grch37, grch38 or grch36 ('auto' decides from "
        "the file's evidence)."
    ),
)
@click.option(
    "--build2",
    type=click.Choice(["grch37", "grch38", "grch36", "auto"], case_sensitive=False),
    default="auto",
    help=(
        "Genome build of file 2: grch37, grch38 or grch36 ('auto' decides from "
        "the file's evidence)."
    ),
)
def compare(
    file1: Path,
    file2: Path,
    fmt1: str | None,
    fmt2: str | None,
    sample1: str | None,
    sample2: str | None,
    build1: str,
    build2: str,
) -> None:
    """Compare two genotype files — coverage overlap and concordance.

    Reports shared named identifiers and same-build unnamed records with
    literal input REF/ALT identity, file-specific records, genotype agreement,
    strand-flip matches (complementary alleles on opposite strands),
    discordant calls, and strand-ambiguous positions. Coverage separates
    parsed named/unnamed counts, eligible unique keys and omitted records.
    """
    from allelix.compare import compare_variants
    from allelix.reports._pipeline import decide_file_build

    parser1 = _helpers._resolve_parser(file1, fmt1, sample=sample1)
    parser2 = _helpers._resolve_parser(file2, fmt2, sample=sample2)
    variants1 = list(parser1.parse(file1))
    variants2 = list(parser2.parse(file2))

    # The same decision as every other command (ADR-0044); a file that
    # must stop stops here rather than taking a tentative winner.
    decided: list[str] = []
    for number, flag, parser, path, variants, choice in (
        (1, "--build1", parser1, file1, variants1, build1),
        (2, "--build2", parser2, file2, variants2, build2),
    ):
        # A structural-only VCF or a file with no genotype records gets its own
        # error whatever decides the build.
        _helpers._reject_unusable(parser, path, len(variants))
        try:
            decision = decide_file_build(
                parser,
                path,
                parser.get_metadata(path),
                override=_helpers._normalize_cli_build(choice),
                variants=variants,
            )
            late = detect_build(variants)
            check_late_seeds(decision, late, path.name)
        except UnusableInputError as exc:
            raise click.ClickException(str(exc)) from exc
        except BuildDecisionError as exc:
            raise click.ClickException(
                f"{exc} For compare, give the build of this file with {flag}."
            ) from exc
        late_record = late if needs_late_seeds(decision) else None
        _helpers._emit_build_notices(
            BuildDiagnostics.from_decision(decision, late_record),
            heading=f"File {number} ({path.name})",
        )
        decided.append(decision.build)

    result = compare_variants(variants1, variants2, build1=decided[0], build2=decided[1])

    if result.build1 != result.build2:
        console.print(
            f"[yellow]Warning: builds differ ({result.build1} vs {result.build2}). "
            "Unnamed records are omitted from comparison.[/yellow]"
        )

    table = Table(title="Coverage Summary")
    table.add_column("Metric", style="bold")
    table.add_column("Value", justify="right")
    table.add_row("File 1", f"{file1.name} ({result.file1_total:,} eligible variants)")
    table.add_row("File 2", f"{file2.name} ({result.file2_total:,} eligible variants)")
    for number in (1, 2):
        table.add_row(
            f"Input named / unnamed (file {number})",
            f"{getattr(result, f'file{number}_input_named'):,} / "
            f"{getattr(result, f'file{number}_input_unnamed'):,}",
        )
        table.add_row(
            f"Eligible named / unnamed (file {number})",
            f"{getattr(result, f'file{number}_named'):,} / "
            f"{getattr(result, f'file{number}_unnamed'):,}",
        )
        table.add_row(
            f"Omitted unnamed records (file {number})",
            f"{getattr(result, f'file{number}_unnamed_omitted'):,}",
        )
    table.add_row("Build (file 1)", result.build1)
    table.add_row("Build (file 2)", result.build2)
    table.add_row("Shared named identifiers", f"{result.shared_named:,}")
    table.add_row("Unnamed variants compared by position", f"{result.shared_unnamed:,}")
    table.add_row("Shared eligible variants", f"{result.shared:,}")
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
    help=(
        "Base path for .bed/.bim/.fam (default: input stem). Refuses input aliases; "
        "stages all three before atomic publication per file."
    ),
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
    tools (plink2 PCA, ADMIXTURE, PRSice). REF/ALT for allele coding come
    from the export's own build: the rsID identity table on GRCh37, gnomAD
    on GRCh38. A site without a known pair falls back to monomorphic (A2=0)
    for homozygous calls. A call whose alleles are not the site's REF and
    ALT as written is exported as missing and counted. A variant without an
    rsID gets the ID chr:pos:ref:alt.
    """
    from allelix.exporters.plink import (
        EmptyPlinkSetError,
        coded_against,
        default_plink_prefix,
        export_plink,
        plink_paths,
        plink_record_key,
        plink_writable,
        resolve_ref_alt_via_gnomad,
        resolve_ref_alt_via_gnomad_grch37,
        resolve_ref_alt_via_identity,
    )

    parser = _helpers._resolve_parser(file_path, fmt, sample=sample)
    prefix = output_prefix if output_prefix else default_plink_prefix(file_path)
    from allelix.reports import validate_output_path

    try:
        for destination in plink_paths(prefix).values():
            validate_output_path(file_path, destination)
    except (OSError, RuntimeError, ValueError) as exc:
        raise click.ClickException(str(exc)) from exc
    build_override = _helpers._normalize_cli_build(build)
    resolved = resolve_data_dir(data_dir)

    # rsID-less records (VCF ID=.) are resolved to rsIDs before
    # the rsID-keyed gnomAD ref/alt lookup below, so they get proper A1/A2
    # coding instead of exporting as monomorphic.
    with _helpers._resolved_variants(
        parser, file_path, data_dir=resolved, build_override=build_override
    ) as resolved_stream:
        variants = list(resolved_stream.variants)
        effective_build = resolved_stream.decision.build
        diagnostics = resolved_stream.build_diagnostics()
        _helpers._emit_build_notices(diagnostics)

    # Sort by chromosome then position so the .bim has contiguous
    # chromosome blocks — PLINK1.9 rejects split chromosomes.
    chrom_order = {str(i): i for i in range(1, 23)}
    chrom_order.update({"X": 23, "Y": 24, "XY": 25, "MT": 26})
    variants.sort(
        key=lambda v: (chrom_order.get(v.chromosome, 99), v.chromosome, v.position),
    )

    # Coordinates, REF/ALT and genotype coding all come from one build.
    grch38 = effective_build == "GRCh38"
    identity_missing = False
    ref_alt_map: dict[str, tuple[str, str]] = {}
    # GRCh37 (and other builds): chosen per call, never shared by rsID.
    record_alleles: dict[RecordKey, tuple[str, str]] | None = None if grch38 else {}
    identity_map: dict[RecordKey, tuple[str, str]] | None = {}
    gnomad_map: dict[RecordKey, tuple[str, str]] = {}
    if grch38:
        variant_by_rsid: dict[str, Variant] = {v.rsid: v for v in variants if not v.is_no_call}
        try:
            ref_alt_map = resolve_ref_alt_via_gnomad(variant_by_rsid, resolved)
        except Exception:
            console.print(
                "[yellow]gnomAD coordinate resolution failed; "
                "using fallback allele coding.[/yellow]"
            )
    elif effective_build == "GRCh37":
        try:
            identity_map = resolve_ref_alt_via_identity(variants, resolved)
        except Exception:
            console.print(
                "[yellow]rsID identity table lookup failed; using fallback allele coding.[/yellow]"
            )
            identity_map = {}
        identity_missing = identity_map is None
        record_alleles = dict(identity_map or {})
        # Then gnomAD's pair where it is not palindromic and the call fits it
        # as written (a strand-inverted site cannot); else the fallback.
        try:
            gnomad_map = resolve_ref_alt_via_gnomad_grch37(variants, resolved, record_alleles)
        except Exception:
            console.print("[yellow]gnomAD lookup failed; using fallback allele coding.[/yellow]")
            gnomad_map = {}
        record_alleles.update(gnomad_map)

    try:
        counts = export_plink(
            iter(variants),
            prefix,
            effective_build,
            ref_alt_map or None,
            record_alleles=record_alleles,
            # Only an ID with a colon can collide with a generated one.
            supplied_ids={v.rsid for v in variants if ":" in v.rsid},
        )
    except EmptyPlinkSetError as exc:
        reasons = [f"{exc.no_calls:,} no-calls", f"{exc.indels:,} indels"]
        if diagnostics.mt_skipped:
            reasons.append(
                f"{diagnostics.mt_skipped:,} mitochondrial records not exported "
                "(chrM has the hg19 length)"
            )
        raise click.ClickException(
            f"Nothing to export from {file_path.name}: no variant can be written to a "
            f"PLINK set (skipped: {'; '.join(reasons)}). No files were written."
        ) from exc
    written, skipped, indel_skip, mono = counts[:4]
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
    source = "gnomAD" if grch38 else "rsID identity table"
    if counts.allele_mismatches:
        console.print(
            f"[yellow]{counts.allele_mismatches:,} call(s) fit neither {source} allele at "
            "their site and were written as missing.[/yellow]"
        )
    unnamed = sum(1 for v in variants if not v.rsid and plink_writable(v))
    if unnamed:
        console.print(f"{unnamed:,} variant(s) without an rsID were given the ID chr:pos:ref:alt.")
    if grch38 and not ref_alt_map:
        console.print(
            "[yellow]gnomAD not available — all homozygous markers exported "
            "as monomorphic.[/yellow]"
        )
        console.print("[yellow]Run `allelix db update` first for proper allele coding.[/yellow]")
    elif not grch38:
        writable = [v for v in variants if plink_writable(v)]
        # As the writer decides: a call that fits a palindromic identity pair
        # gets the fallback; one that does not fit is written missing against it.
        identity_pairs = identity_map or {}
        from_identity = sum(
            1
            for v in writable
            if plink_record_key(v) in identity_pairs
            and coded_against(v, identity_pairs[plink_record_key(v)])
        )
        from_gnomad = sum(1 for v in writable if plink_record_key(v) in gnomad_map)
        fallback = len(writable) - from_identity - from_gnomad
        if effective_build == "GRCh37":
            console.print(
                f"Allele coding: {from_identity:,} site(s) from the rsID identity table, "
                f"{from_gnomad:,} from gnomAD, {fallback:,} with fallback coding."
            )
        elif fallback:
            console.print(
                f"[yellow]{fallback:,} site(s) used fallback coding on this build.[/yellow]"
            )
        if identity_missing:
            console.print(
                "[yellow]Run `allelix db update` first for proper allele coding.[/yellow]"
            )
    console.print(
        "[dim]Single-sample export. Merging with other samples requires "
        "allele harmonization (--merge-mode or set-all-var-ids).[/dim]"
    )

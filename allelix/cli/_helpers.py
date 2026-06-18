# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""CLI orchestration helpers and the shared Rich console.

This module is the home for everything subcommands reach for that isn't
itself a click command: parser/annotator wiring, the analysis pipeline
orchestrator, the build diagnostics emitter, filter-file parsing, etc.

The Rich ``console`` instance is also defined here. All subcommand modules
import it from here for consistent output styling.
"""

from __future__ import annotations

import logging
import re
import sys
import time
from typing import TYPE_CHECKING, NamedTuple, cast

import click
from rich.console import Console

from allelix.annotators import get_annotators
from allelix.databases import resolve_data_dir
from allelix.parsers import ParserNotFoundError, detect_parser, get_parser_by_name
from allelix.reports._pipeline import rollup_gwas_duplicates, run_analysis
from allelix.reports.diff import compute_diff, load_previous_report
from allelix.reports.high_value import format_warnings, load_high_value_snps, scan_no_calls
from allelix.reports.html import render_html
from allelix.reports.json_report import render_json
from allelix.reports.terminal import render_terminal, render_terminal_diff

if TYPE_CHECKING:
    from pathlib import Path

    from allelix.annotators.alphamissense import AlphaMissenseAnnotator
    from allelix.annotators.base import Annotator
    from allelix.annotators.cadd import CaddAnnotator
    from allelix.annotators.gnomad import GnomadAnnotator
    from allelix.parsers.base import GenotypeParser

console = Console()

# Sort 1-22 numerically, then X, Y, MT, then anything else alphabetically.
_NAMED_CHROM_ORDER = {"X": 0, "Y": 1, "MT": 2}


def _chrom_sort_key(chrom: str) -> tuple[int, int, str]:
    """Sort key: autosomes (1-22), then X/Y/MT, then unknowns alphabetically."""
    if chrom.isdigit():
        return (0, int(chrom), "")
    if chrom in _NAMED_CHROM_ORDER:
        return (1, _NAMED_CHROM_ORDER[chrom], "")
    return (2, 0, chrom)


def _percent(part: int, total: int) -> str:
    if total == 0:
        return "0.00%"
    return f"{part / total * 100:.2f}%"


class _WarningCounter(logging.Handler):
    """Count warning records emitted by the parser pipeline."""

    def __init__(self) -> None:
        super().__init__(level=logging.WARNING)
        self.count = 0

    def emit(self, record: logging.LogRecord) -> None:
        self.count += 1


class _LoggerSnapshot(NamedTuple):
    """Captured state of a Python logger for restoration after CLI mutates it."""

    level: int
    propagate: bool


def _wire_parser_logging() -> tuple[_WarningCounter, logging.Handler, _LoggerSnapshot]:
    """Attach warning capture + stderr surfacing to the parsers logger."""
    parser_logger = logging.getLogger("allelix.parsers")
    counter = _WarningCounter()
    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setLevel(logging.WARNING)
    stderr_handler.setFormatter(logging.Formatter("warning: %(message)s"))
    snapshot = _LoggerSnapshot(level=parser_logger.level, propagate=parser_logger.propagate)
    parser_logger.addHandler(counter)
    parser_logger.addHandler(stderr_handler)
    parser_logger.setLevel(logging.WARNING)
    parser_logger.propagate = False
    return counter, stderr_handler, snapshot


def _unwire_parser_logging(
    counter: _WarningCounter,
    stderr_handler: logging.Handler,
    snapshot: _LoggerSnapshot,
) -> None:
    parser_logger = logging.getLogger("allelix.parsers")
    parser_logger.removeHandler(counter)
    parser_logger.removeHandler(stderr_handler)
    parser_logger.setLevel(snapshot.level)
    parser_logger.propagate = snapshot.propagate


def _resolve_parser(file_path: Path, fmt: str | None, sample: str | None = None) -> GenotypeParser:
    """Resolve the parser for ``file_path``, with VCF sample binding.

    For VCF inputs, pre-flight checks the sample column header:

    - Multi-sample VCF without ``--sample`` → fail loudly with the list
      of available samples (instead of letting the parser raise
      mid-stream from the analysis pipeline).
    - ``--sample`` provided but not in the file → fail with the available
      list.
    - Single-sample VCF → ``--sample`` is silently ignored if also given
      (it's not wrong, just redundant).
    - ``--sample`` on a non-VCF format → silently ignored (no error;
      the option is VCF-specific by design).
    """
    try:
        parser = get_parser_by_name(fmt) if fmt else detect_parser(file_path)
    except ParserNotFoundError as exc:
        raise click.ClickException(str(exc)) from exc

    from allelix.parsers.vcf import (
        MultiSampleError,
        SampleNotFoundError,
        VcfParser,
    )

    if isinstance(parser, VcfParser):
        bound = VcfParser(sample=sample) if sample is not None else parser
        try:
            bound.validate_sample(file_path)
        except (MultiSampleError, SampleNotFoundError) as exc:
            raise click.ClickException(str(exc)) from exc
        return bound

    return parser


def _ready_annotators(
    data_dir: Path | None,
    *,
    include_benign: bool = False,
    gwas_filter_traits: bool = True,
    cadd_full: bool = False,
) -> tuple[Path, list[Annotator], list[Annotator]]:
    resolved = resolve_data_dir(data_dir)
    annotators = get_annotators(
        resolved,
        include_benign=include_benign,
        gwas_filter_traits=gwas_filter_traits,
        cadd_full=cadd_full,
    )
    ready: list[Annotator] = []
    not_ready: list[Annotator] = []
    for a in annotators:
        if a.is_ready():
            ready.append(a)
        else:
            not_ready.append(a)
    if not ready:
        names = ", ".join(a.name for a in annotators)
        raise click.ClickException(
            f"No annotators are ready. Run `allelix db update` first. Registered: {names}"
        )
    return resolved, ready, not_ready


_STALENESS_SECONDS = 7 * 24 * 60 * 60  # 7 days


def _run_setup(annotator: Annotator) -> bool:
    """Invoke annotator.setup(). Returns True on success, False on failure."""
    try:
        annotator.setup()
    except Exception as exc:
        if hasattr(exc, "close"):
            exc.close()
        console.print(f"  [red]{annotator.name}: {exc}[/red]")
        return False
    sig = getattr(annotator, "cached_remote_signal", lambda: None)()
    if sig and "cpic:unavailable" in sig:
        console.print(
            f"  [yellow]{annotator.name}: updated (CPIC unavailable — "
            "non-finding filter degraded, retry later)[/yellow]"
        )
    return True


def _maybe_refresh_databases(data_dir: Path) -> None:
    """Check database mtimes; refresh any that are stale and have a changed remote signal.

    Only runs for annotators that download data (SNPedia excluded).
    If the network is unreachable, warns and continues with stale caches.
    """
    now = time.time()
    annotators = get_annotators(data_dir)
    for annotator in annotators:
        with annotator:
            if not annotator.requires_download or not annotator.is_ready():
                continue
            # Code-driven sources (commit-pinned HF caches) never change
            # at a fixed URL — skip the HEAD request. See ADR-0030.
            if not annotator.server_driven_freshness:
                continue
            db_files = list(data_dir.glob(f"{annotator.name}*sqlite*"))
            if not db_files:
                continue
            newest_mtime = max(f.stat().st_mtime for f in db_files)
            age = now - newest_mtime
            if age <= _STALENESS_SECONDS:
                continue

            remote = annotator.fetch_remote_signal()
            if remote is None:
                age_days = int(age / 86400)
                console.print(
                    f"[yellow]{annotator.display_name} database is {age_days} days old. "
                    "Run `allelix db update` when online.[/yellow]"
                )
                continue

            cached = annotator.cached_remote_signal()
            if cached == remote:
                continue

            console.print(f"[bold]Updating {annotator.display_name}…[/bold]")
            if _run_setup(annotator):
                console.print(
                    f"[green]✓ {annotator.display_name} updated[/green] "
                    f"(version {annotator.version() or '(unknown)'})"
                )


def _format_from_path(output: Path, override: str | None) -> str:
    if override:
        return override.lower()
    suffix = output.suffix.lower()
    if suffix == ".html":
        return "html"
    if suffix == ".json":
        return "json"
    raise click.ClickException(
        f"Cannot infer report format from {output.name!r}. "
        "Pass --report-format html|json explicitly."
    )


_RSID_PATTERN = re.compile(r"^rs\d+$", re.IGNORECASE)


def _parse_filter_file(path: Path) -> tuple[frozenset[str], frozenset[str]]:
    r"""Parse a filter file into ``(gene_names, rsids)``.

    Lines matching ``^rs\d+$`` (case-insensitive) are rsIDs. Everything
    else is a gene name. Lines starting with ``#`` and blank lines are
    ignored. Gene names starting with ``RS`` (e.g., RSPO1, RSF1) are
    correctly classified as gene names, not rsIDs.

    Input is case-tolerant; output is canonical: rsIDs are normalized to
    lowercase (``rs1801133``), gene names to uppercase (``MTHFR``). The
    filter recorded in JSON output therefore looks identical regardless
    of how the user typed the entries in the filter file.
    """
    genes: set[str] = set()
    rsids: set[str] = set()
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if _RSID_PATTERN.match(line):
            rsids.add(line.lower())
        else:
            genes.add(line.upper())
    return frozenset(genes), frozenset(rsids)


def _resolve_clinvar_builds(value: str) -> tuple[str, ...]:
    """Map a `db update --build` value to a tuple of build identifiers."""
    v = (value or "both").strip().lower()
    if v == "both":
        return ("GRCh37", "GRCh38")
    if v == "grch37":
        return ("GRCh37",)
    if v == "grch38":
        return ("GRCh38",)
    raise click.ClickException(f"Unknown --build value {value!r}")


def _normalize_cli_build(value: str | None) -> str | None:
    """Map a --build CLI value to a canonical build identifier or None for auto."""
    if value is None:
        return None
    v = value.strip().lower()
    if v in ("", "auto"):
        return None
    if v == "grch37":
        return "GRCh37"
    if v == "grch38":
        return "GRCh38"
    raise click.ClickException(f"Unknown --build value {value!r}")


def _pop_enrichment_annotator(
    ready: list[Annotator], name: str, *, skip: bool
) -> tuple[Annotator | None, list[Annotator]]:
    """Pop a named enrichment annotator off the ready list.

    Used for gnomAD, AlphaMissense, and CADD — annotators that the pipeline
    routes via dedicated parameters rather than the generic annotator loop.
    Returns ``(annotator_or_None, ready_without_it)``. When ``skip`` is True,
    the annotator is still removed from ``ready`` but the returned reference
    is None — the enrichment is being explicitly disabled by the caller.
    """
    found = None
    if not skip:
        for a in ready:
            if a.name == name:
                found = a
                break
    return found, [a for a in ready if a.name != name]


def _emit_build_diagnostics(result: object) -> None:
    """Print a one-line build banner and a warning on header/data mismatch."""
    diag = getattr(result, "build_diagnostics", None)
    if diag is None:
        return
    matched = f"{diag.matched_count}/{diag.inspected_count}" if diag.inspected_count else "0/0"
    if diag.override:
        source = "override"
    elif diag.detected_build:
        source = "detected"
    elif diag.header_build:
        source = "header (no position confirmation)"
    elif diag.chr_prefix_inferred:
        # GH #38: chr-prefixed contig names ("chr1", "chrX", ...) reliably
        # indicate GRCh38 in modern caller output. We DID detect a build;
        # the banner and the warning should say so instead of reading as
        # a blind default.
        source = "inferred from chr-prefixed contig names"
    else:
        source = "fallback (no known SNPs matched)"
    console.print(
        f"[dim]Build: {diag.effective_build} ({source}; "
        f"{matched} known-SNP positions matched)[/dim]"
    )
    if diag.mismatch:
        console.print(
            f"[yellow]Build mismatch: file header claims {diag.header_build} but "
            f"position data is {diag.detected_build}. Using {diag.detected_build}. "
            f"This is a real-world data-quality issue — your provider may have "
            f"mislabeled the build (see ADR-0021).[/yellow]"
        )
    elif diag.chr_prefix_inferred:
        # GH #38: positive, accurate message — the inference path
        # actually fired. Still recommend `--build` for users who
        # want to lock in the answer; chr-prefix is a strong signal
        # but UCSC hg19 also uses `chr` prefixes, so the heuristic
        # isn't guaranteed against a hg19-converted file.
        console.print(
            f"[dim]Inferred {diag.effective_build} from chr-prefixed contig "
            f"names (GRCh38 convention). Pass --build grch37 if this file is "
            f"UCSC hg19 with chr-prefixed contigs instead.[/dim]"
        )
    elif not diag.override and diag.detected_build is None and diag.header_build is None:
        # Common shape: VCF from a variant caller where the ID column is `.`
        # and the header has no ##contig assembly tag, AND no chr-prefix
        # signal was observed. All three auto-detect paths failed.
        # Loudly recommend an explicit --build because picking the wrong one
        # silently means every annotation lookup uses wrong coordinates.
        console.print(
            f"[yellow]Could not auto-detect genome build (no rsIDs in input, "
            f"no ##contig assembly tag, no chr-prefixed contigs). Defaulted to "
            f"{diag.effective_build}. If the file is the other build, pass "
            f"--build grch37 or --build grch38 explicitly — annotation "
            f"coordinates differ between builds and silently using the wrong "
            f"one will miss every hit.[/yellow]"
        )
    elif (
        not diag.override
        and diag.detected_build is None
        and diag.header_build is not None
        and diag.inspected_count > 0
    ):
        # Position-detection inspected known-rsID rows but couldn't pick a
        # build — either votes tied across builds or no row matched any
        # build's reference position. Without this warning, the pipeline
        # silently falls through to header_build, and a GRCh36 file with a
        # GRCh37-mislabeled header gets the GRCh37 ClinVar cache (the
        # silent-coords trap #15). The dim "header (no position
        # confirmation)" status line shows the same facts but reads as
        # routine — yellow is what the situation deserves.
        console.print(
            f"[yellow]Build detection inconclusive: "
            f"{diag.inspected_count} known-rsID position checks ran but "
            f"did not converge on a build. Using the file's header-claimed "
            f"build ({diag.header_build}), which has not been confirmed "
            f"against your position data. If the file is actually a "
            f"different build, pass --build grch37 or --build grch38 to "
            f"force — wrong coordinates will silently mis-annotate every "
            f"variant.[/yellow]"
        )
    if diag.effective_build == "GRCh36":
        console.print(
            "[yellow]Warning: GRCh36 (hg18) detected. rsID-based annotations "
            "(ClinPGx, GWAS Catalog, SNPedia, gnomAD) are complete. ClinVar "
            "position-matching is skipped (no GRCh36 cache — see ADR-0025). "
            "For full ClinVar coverage, liftOver to GRCh38 first: "
            "docs/grch36-liftover.md[/yellow]"
        )


def _emit_runtime_nudges(result: object) -> None:
    """Surface runtime-degradation notes the analyze run can detect.

    Both are once-per-run, post-pipeline. Neither blocks the report; the
    intent is transparency for cases where the silent default behavior
    used to leave the user uncertain about WHY their count looked low.

    - **GH #90: strand-aware matching inactive.** When carrier/genotype
      matching ran with ``Variant.ref`` unavailable across the input —
      typically array data with ``--no-gnomad`` or no gnomAD cache —
      the strand-flip path in ``utils/allele.py`` falls back to the
      v2.0.x direct-match-only behavior. Behavior is correct (forward-
      normalized array data matches the direct path regardless), just
      degraded relative to v2.1+'s contract. VCF inputs always carry
      REF from the file, so this never fires for them.
    - **GH #91: GRCh38 rsID-less undercount.** When the effective
      build is GRCh38 AND most variants entered the pipeline without
      rsIDs (the variant-caller ``ID=.`` case), annotation coverage
      is materially lower than the GRCh37 equivalent because rsID
      resolution falls back to ClinVar position lookup — a smaller
      surface than dbSNP. Real fix is #62; this is the interim
      honesty signal until then.
    """
    no_ref = getattr(result, "no_ref_variant_count", 0)
    rsidless = getattr(result, "rsidless_variant_count", 0)
    total = getattr(result, "total_variants", 0)
    parser_name = getattr(result, "parser_name", "")
    diag = getattr(result, "build_diagnostics", None)
    effective_build = diag.effective_build if diag is not None else None

    # GH #90: array inputs (everything except the VCF parser) running
    # without a usable gnomAD cache leave every variant ref=None. Gate
    # on "the whole input lacked ref" rather than "any variant lacked
    # ref" so a hybrid file (e.g. a few ref-less rows in an otherwise
    # populated VCF) doesn't flap the message on.
    if parser_name != "vcf" and total > 0 and no_ref == total:
        console.print(
            "[dim]strand-aware matching inactive — no reference context "
            "(run `allelix db update` for gnomAD-backed ref resolution, "
            "or this is expected with `--no-gnomad`)[/dim]"
        )

    # GH #91: GRCh38 + a high rsID-less fraction. The threshold (>50%)
    # distinguishes a variant-caller VCF (mostly ID=. — fires) from a
    # 23andMe-style array (rsID-bearing — no warning) or a mostly-
    # rsID-bearing VCF where a handful of novel calls came in with
    # ID=. (no warning, ID=. for novel sites is expected).
    if effective_build == "GRCh38" and total > 0 and rsidless / total > 0.5:
        console.print(
            "[yellow]GRCh38 input without rsIDs: annotation coverage is "
            "currently lower than GRCh37 because rsID resolution falls "
            "back to ClinVar positions. Full dbSNP resolution is planned "
            "(#62). If you have a GRCh37 build available, it currently "
            "surfaces more annotations.[/yellow]"
        )


def _run_analysis_command(
    file_path: Path,
    fmt: str | None,
    data_dir: Path | None,
    output: Path | None,
    report_format: str | None,
    min_magnitude: float,
    category: str | None,
    genes: frozenset[str] | None,
    rsids: frozenset[str] | None = None,
    build: str | None = None,
    include_benign: bool = False,
    gwas_min_magnitude: float | None = None,
    snpedia_min_magnitude: float | None = None,
    exclude_sources: frozenset[str] | None = None,
    gwas_all: bool = False,
    diff_path: Path | None = None,
    no_update: bool = False,
    no_gnomad: bool = False,
    no_alphamissense: bool = False,
    no_cadd: bool = False,
    sample: str | None = None,
) -> None:
    resolved = resolve_data_dir(data_dir)
    if not no_update:
        _maybe_refresh_databases(resolved)
    parser = _resolve_parser(file_path, fmt, sample=sample)

    from allelix.config import load_config

    cfg = load_config(resolved)
    _, ready, not_ready = _ready_annotators(
        data_dir,
        include_benign=include_benign,
        gwas_filter_traits=not gwas_all,
        cadd_full=cfg.cadd_full,
    )
    annotator_classes = {type(a).name: type(a) for a in ready}
    ready = [a for a in ready if cfg.is_enabled(a.name, annotator_classes)]

    if exclude_sources:
        ready = [a for a in ready if a.name not in exclude_sources]

    gnomad_annotator, ready = _pop_enrichment_annotator(ready, "gnomad", skip=no_gnomad)
    am_annotator, ready = _pop_enrichment_annotator(ready, "alphamissense", skip=no_alphamissense)
    cadd_annotator, ready = _pop_enrichment_annotator(ready, "cadd", skip=no_cadd)

    if not_ready:
        names = [a.name for a in not_ready]
        console.print(
            f"[yellow]Skipping unready annotators: {', '.join(names)}[/yellow] "
            "(run `allelix db update` to populate)"
        )

    all_active: list[Annotator] = list(ready)
    if gnomad_annotator is not None and gnomad_annotator.is_ready():
        all_active.append(gnomad_annotator)
    if am_annotator is not None and am_annotator.is_ready():
        all_active.append(am_annotator)
    if cadd_annotator is not None and cadd_annotator.is_ready():
        all_active.append(cadd_annotator)
    versions = ", ".join(f"{a.display_name} ({a.version() or 'unknown'})" for a in all_active)
    console.print(f"[dim]Analyzing against: {versions}[/dim]")

    high_value = load_high_value_snps()
    hv_rsids = set(high_value)

    counter, stderr_handler, snapshot = _wire_parser_logging()
    try:
        result = run_analysis(
            file_path,
            parser,
            ready,
            skipped_count_provider=lambda: counter.count,
            build_override=build,
            gnomad=cast("GnomadAnnotator | None", gnomad_annotator),
            alphamissense=cast("AlphaMissenseAnnotator | None", am_annotator),
            cadd=cast("CaddAnnotator | None", cadd_annotator),
            high_value_rsids=hv_rsids,
            panel_rsids=rsids,
        )
    finally:
        _unwire_parser_logging(counter, stderr_handler, snapshot)

    _emit_build_diagnostics(result)
    _emit_runtime_nudges(result)

    hv_warnings = scan_no_calls(result.hv_variants, high_value)
    if hv_warnings:
        console.print(
            f"[bold red]Warning:[/bold red] {len(hv_warnings)} high-value SNP(s) returned no-call:"
        )
        for line in format_warnings(hv_warnings):
            console.print(f"  [red]⚠[/red] {line}")

    if counter.count:
        console.print(
            f"[yellow]Note:[/yellow] {counter.count:,} malformed line(s) skipped "
            "(see warnings on stderr)."
        )

    source_floors: dict[str, float] | None = None
    if gwas_min_magnitude is not None or snpedia_min_magnitude is not None:
        source_floors = {}
        if gwas_min_magnitude is not None:
            source_floors["gwas"] = gwas_min_magnitude
        if snpedia_min_magnitude is not None:
            source_floors["snpedia"] = snpedia_min_magnitude

    diff_result = None
    if diff_path is not None:
        try:
            prev = load_previous_report(diff_path)
        except ValueError as exc:
            raise click.ClickException(str(exc)) from exc
        filtered_for_diff = result.filter(
            min_magnitude=min_magnitude,
            category=category,
            genes=genes,
            rsids=rsids,
            source_min_magnitudes=source_floors,
        )
        filtered_for_diff = rollup_gwas_duplicates(filtered_for_diff)
        diff_result = compute_diff(
            filtered_for_diff,
            prev["annotations"],
            prev.get("generated_at", ""),
        )

    if output is None:
        if diff_result is not None:
            rendered = render_terminal_diff(diff_result, console)
        else:
            rendered = render_terminal(
                result,
                console=console,
                min_magnitude=min_magnitude,
                category=category,
                genes=genes,
                rsids=rsids,
                source_min_magnitudes=source_floors,
            )
    else:
        chosen = _format_from_path(output, report_format)
        hv_warning_lines = format_warnings(hv_warnings) if hv_warnings else None
        if chosen == "json":
            hv_dicts = (
                [{"rsid": w.snp.rsid, "gene": w.snp.gene, "note": w.snp.note} for w in hv_warnings]
                if hv_warnings
                else None
            )
            rendered = render_json(
                result,
                output_path=output,
                min_magnitude=min_magnitude,
                category=category,
                genes=genes,
                rsids=rsids,
                source_min_magnitudes=source_floors,
                diff=diff_result,
                high_value_no_calls=hv_dicts,
            )
        else:
            rendered = render_html(
                result,
                output_path=output,
                min_magnitude=min_magnitude,
                category=category,
                genes=genes,
                rsids=rsids,
                source_min_magnitudes=source_floors,
                diff=diff_result,
                high_value_no_calls=hv_warning_lines,
            )
        console.print(f"[green]Wrote {rendered:,} annotation(s) to {output}[/green]")

    console.print(
        f"[dim]{len(result.annotations):,} total annotation(s) from {len(ready)} "
        f"database(s) across {result.total_variants:,} variant(s).[/dim]"
    )

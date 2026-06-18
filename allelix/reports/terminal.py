# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Terminal report rendering for `allelix analyze`.

GH #9: the terminal output is a quick-eyeball view (custom panels,
sanity checks, extracts). Deep-dive enrichment data (allele frequency,
AlphaMissense, CADD, ClinVar review status, derived zygosity) belongs
in the HTML/JSON report — not in a 12-column-wide terminal table that
Rich squeezes to hairline-zero-width on typical 100-120 col terminals.

The bare-minimum column set kept here is:

    rsID | Gene? | Source | Significance | Mag | GT | Condition?

Gene and Condition are conditional on at least one filtered row having
a value. Source uses a short attribution (``GWAS Catalog`` → ``GWAS``)
and Significance drops the redundant ``source_`` prefix that the
Source column already conveys.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from rich.table import Table

from allelix.reports._pipeline import rollup_gwas_duplicates

if TYPE_CHECKING:
    from collections.abc import Iterable

    from rich.console import Console

    from allelix.models import Annotation
    from allelix.reports._pipeline import AnalysisResult
    from allelix.reports.diff import DiffResult


def render_terminal(
    result: AnalysisResult,
    console: Console,
    *,
    min_magnitude: float = 0.0,
    category: str | None = None,
    genes: Iterable[str] | None = None,
    rsids: Iterable[str] | None = None,
    source_min_magnitudes: dict[str, float] | None = None,
) -> int:
    """Render an AnalysisResult as a Rich table. Returns annotation count.

    Per ADR-0003 (regulatory posture), every row shows the source attribution
    in its own column — no rendered claim is unattributed.
    """
    filtered = result.filter(
        min_magnitude=min_magnitude,
        category=category,
        genes=genes,
        rsids=rsids,
        source_min_magnitudes=source_min_magnitudes,
    )
    filtered = rollup_gwas_duplicates(filtered)
    _print_table(filtered, console)
    _print_panel_coverage_warning(result, console, filtered)
    _print_regulatory_notice(console)
    return len(filtered)


def _print_regulatory_notice(console: Console) -> None:
    """Emit the ADR-0003 regulatory notice in terminal output.

    Evaluator defect 1: HTML and JSON already carry the same notice.
    The terminal is the default analyze output; without this line the
    disclaimer that lives in the other surfaces is missing from the
    most common one.
    """
    from allelix.reports import REGULATORY_NOTICE

    console.print(f"\n[dim italic]{REGULATORY_NOTICE}[/dim italic]")


def _print_panel_coverage_warning(
    result: AnalysisResult,
    console: Console,
    filtered_annotations: list[Annotation],
) -> None:
    """GH #75: surface panel rsIDs that weren't in the user's input file.

    "Not on your chip" is critically different from "homozygous reference"
    — making this invisible was the original audit complaint. Quiet when
    no panel was supplied, or when every panel rsID was genotyped.

    GH #106: the post-filter annotation list is threaded through so
    `panel_coverage()` derives "annotated" from what was actually rendered,
    not from the unfiltered set. Without this thread, panel rsids whose
    only annotations were below the magnitude floor showed in "found" but
    not in any rendered surface (the 9-of-20 limbo bug).
    """
    coverage = result.panel_coverage(filtered_annotations)
    if coverage is None:
        return
    missing = coverage["missing"]
    requested = coverage["requested"]
    if not missing:
        return
    sample = ", ".join(missing[:3])
    if len(missing) > 3:
        sample = f"{sample}, …"
    console.print(
        f"[yellow]⚠[/yellow] {len(missing)}/{requested} panel variants not "
        f"found in input: {sample}"
    )


def render_terminal_diff(
    diff: DiffResult,
    console: Console,
) -> int:
    """Render a diff summary and tables for new/changed/removed annotations.

    Uses the same bare-min column set as the main annotation table.
    """
    from allelix.reports.diff import summarize_diff

    summary = summarize_diff(diff)
    if not diff.has_changes:
        console.print(f"[green]{summary}[/green]")
        return 0

    console.print(f"[bold]{summary}[/bold]")
    total = 0

    if diff.new:
        has_gene = any(a.gene for a in diff.new)
        has_condition = any(a.condition for a in diff.new)
        table = Table(title=f"New Annotations ({len(diff.new)})")
        table.add_column("rsID", style="cyan", no_wrap=True, min_width=11)
        if has_gene:
            table.add_column("Gene", style="magenta", no_wrap=True)
        table.add_column("Source", style="blue", no_wrap=True)
        table.add_column("Significance", style="yellow")
        table.add_column("Mag", justify="right", min_width=4)
        table.add_column("GT", no_wrap=True)
        if has_condition:
            table.add_column("Condition", overflow="fold")
        for a in diff.new:
            row = [a.rsid]
            if has_gene:
                row.append(a.gene or "—")
            row.extend(
                [
                    _compact_source(a.attribution),
                    _compact_significance(a.significance, a.source),
                    f"{a.magnitude:.1f}",
                    a.genotype_match,
                ]
            )
            if has_condition:
                row.append(a.condition or "—")
            table.add_row(*row)
        console.print(table)
        total += len(diff.new)

    if diff.changed:
        has_gene = any(c.current.gene for c in diff.changed)
        has_condition = any(c.current.condition for c in diff.changed)
        table = Table(title=f"Changed Annotations ({len(diff.changed)})")
        table.add_column("rsID", style="cyan", no_wrap=True, min_width=11)
        if has_gene:
            table.add_column("Gene", style="magenta", no_wrap=True)
        table.add_column("Source", style="blue", no_wrap=True)
        table.add_column("Old Sig", style="dim")
        table.add_column("New Sig", style="yellow")
        table.add_column("Old Mag", justify="right", style="dim", min_width=4)
        table.add_column("New Mag", justify="right", min_width=4)
        if has_condition:
            table.add_column("Condition", overflow="fold")
        for c in diff.changed:
            prev_mag_str = "—" if c.previous_magnitude is None else f"{c.previous_magnitude:.1f}"
            row = [c.current.rsid]
            if has_gene:
                row.append(c.current.gene or "—")
            row.extend(
                [
                    _compact_source(c.current.attribution),
                    _compact_significance(c.previous_significance, c.current.source),
                    _compact_significance(c.current.significance, c.current.source),
                    prev_mag_str,
                    f"{c.current.magnitude:.1f}",
                ]
            )
            if has_condition:
                row.append(c.current.condition or "—")
            table.add_row(*row)
        console.print(table)
        total += len(diff.changed)

    if diff.removed:
        has_gene = any(d.get("gene") for d in diff.removed)
        has_condition = any(d.get("condition") for d in diff.removed)
        table = Table(title=f"Removed Annotations ({len(diff.removed)})")
        table.add_column("rsID", style="dim cyan", no_wrap=True, min_width=11)
        if has_gene:
            table.add_column("Gene", style="dim magenta", no_wrap=True)
        table.add_column("Source", style="dim blue", no_wrap=True)
        table.add_column("Significance", style="dim")
        table.add_column("Mag", justify="right", style="dim", min_width=4)
        if has_condition:
            table.add_column("Condition", overflow="fold", style="dim")
        for d in diff.removed:
            source = d.get("source", "")
            row = [d.get("rsid", "")]
            if has_gene:
                row.append(d.get("gene", "") or "—")
            row.extend(
                [
                    _compact_source(d.get("attribution", "")),
                    _compact_significance(d.get("significance", ""), source),
                    f"{d.get('magnitude', 0.0):.1f}",
                ]
            )
            if has_condition:
                row.append(d.get("condition", "") or "—")
            table.add_row(*row)
        console.print(table)
        total += len(diff.removed)

    return total


def _compact_significance(significance: str, source: str) -> str:
    """Strip the redundant ``source_`` prefix from significance.

    The Source column already shows the database, so
    ``clinvar_pathogenic`` shown next to ``ClinVar`` is wasteful — show
    just ``pathogenic``. Falls back to the raw value if the prefix
    doesn't match.
    """
    prefix = f"{source}_"
    if source and significance.startswith(prefix):
        return significance[len(prefix) :]
    return significance


def _compact_source(attribution: str) -> str:
    """Shorten multi-word source names that reliably get truncated.

    ``GWAS Catalog`` renders as ``GWAS Ca…`` in narrow terminals; show
    ``GWAS`` instead. Other attributions pass through unchanged.
    """
    if attribution == "GWAS Catalog":
        return "GWAS"
    return attribution


def _print_table(filtered: list[Annotation], console: Console) -> None:
    if not filtered:
        console.print("[yellow]No annotations matched the current filters.[/yellow]")
        return

    has_gene = any(a.gene for a in filtered)
    has_condition = any(a.condition for a in filtered)

    table = Table(title=f"Annotations ({len(filtered)})")
    table.add_column("rsID", style="cyan", no_wrap=True, min_width=11)
    if has_gene:
        table.add_column("Gene", style="magenta", no_wrap=True)
    table.add_column("Source", style="blue", no_wrap=True)
    table.add_column("Significance", style="yellow")
    table.add_column("Mag", justify="right", min_width=4)
    table.add_column("GT", no_wrap=True)
    if has_condition:
        table.add_column("Condition", overflow="fold")

    for a in filtered:
        row = [a.rsid]
        if has_gene:
            row.append(a.gene or "—")
        row.extend(
            [
                _compact_source(a.attribution),
                _compact_significance(a.significance, a.source),
                f"{a.magnitude:.1f}",
                a.genotype_match,
            ]
        )
        if has_condition:
            row.append(a.condition or "—")
        table.add_row(*row)
    console.print(table)

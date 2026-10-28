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

import contextlib
import logging
import re
import sys
import time
from typing import TYPE_CHECKING, NamedTuple, cast

import click
from rich.console import Console
from rich.markup import escape

from allelix.annotators import get_annotators
from allelix.annotators.clinvar import CLINVAR_SUPPORTED_BUILDS, ClinVarAnnotator
from allelix.databases import resolve_data_dir
from allelix.parsers import ParserNotFoundError, detect_parser, get_parser_by_name
from allelix.parsers._helpers import (
    INPUT_NOTICE_CODES,
    INPUT_NOTICE_LOGGER,
    quiet_input_notices,
)
from allelix.reports import (
    array_quality_lines,
    input_method_label,
    is_imputed,
    vcf_quality_lines,
)
from allelix.reports._pipeline import (
    _IDENTITY_BUILDS,
    UnusableInputError,
    rollup_gwas_duplicates,
    run_analysis,
)
from allelix.reports._sources import ReportWarning, build_sources, clinvar_skipped_warning
from allelix.reports.association_text import strand_abstention_lines
from allelix.reports.clinvar_text import (
    common_allele_line,
    scored_down_line,
    scored_down_variant_text,
)
from allelix.reports.diff import compute_diff, load_previous_report
from allelix.reports.high_value import (
    DI_CALLS_NOT_CHECKED,
    UNCHECKED_DI_CALL,
    format_warnings,
    load_high_value_snps,
    scan_no_calls,
    unchecked_sites,
)
from allelix.reports.html import render_html
from allelix.reports.json_report import render_json
from allelix.reports.terminal import render_terminal, render_terminal_diff
from allelix.utils.build_detect import BuildDecisionError, needs_late_seeds

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

    from allelix.annotators.alphamissense import AlphaMissenseAnnotator
    from allelix.annotators.base import Annotator
    from allelix.annotators.cadd import CaddAnnotator
    from allelix.annotators.gnomad import GnomadAnnotator
    from allelix.annotators.rsid_identity import RsidIdentityAnnotator
    from allelix.attribution import Placements
    from allelix.config import AllelixConfig
    from allelix.databases.cache_state import CacheStatus
    from allelix.models import Variant
    from allelix.parsers.base import GenotypeParser
    from allelix.reports._pipeline import BuildDiagnostics, StreamBuild
    from allelix.utils.build_detect import BuildDecision

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


# Report-warning codes the parser logger prints itself, on stderr.
_LOGGED_INPUT_NOTICES = frozenset(INPUT_NOTICE_CODES.values())
# report_warnings code: --build forced a build no known-SNP position matches.
BUILD_OVERRIDE_UNMATCHED = "build_override_unmatched"


class _WarningCounter(logging.Handler):
    """Count the parser pipeline's per-row warnings (skipped rows).

    Notices about the whole input file (ignored trailing bytes, a possible
    cut) are not skipped rows; they are kept in ``input_notices`` as
    ``(code, message)`` so the report can list them (``report_warnings``).
    """

    def __init__(self) -> None:
        super().__init__(level=logging.WARNING)
        self.count = 0
        self.input_notices: list[tuple[str, str]] = []

    def emit(self, record: logging.LogRecord) -> None:
        if record.name != INPUT_NOTICE_LOGGER:
            self.count += 1
            return
        code = getattr(record, "input_notice", None)
        if code is not None:
            self.input_notices.append((code, record.getMessage()))


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


def _resolve_parser(
    file_path: Path,
    fmt: str | None,
    sample: str | None = None,
    *,
    keep_filtered: bool = False,
    min_gq: int | None = None,
) -> GenotypeParser:
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

    ``keep_filtered`` / ``min_gq`` set the VCF quality gate;
    ``min_gq=None`` keeps the parser default. Both are ignored for
    non-VCF formats.

    ``--format`` names a reader whose own detection may not recognise the
    file. The run goes ahead, at the user's risk: the warning is printed
    now, before anything else, and kept on the parser (``forced_warning``)
    so every report carries it (``forced_format`` in ``report_warnings``).
    """
    try:
        parser = get_parser_by_name(fmt) if fmt else detect_parser(file_path)
    except ParserNotFoundError as exc:
        raise click.ClickException(str(exc)) from exc
    warning = _forced_format_warning(file_path, parser) if fmt else None

    from allelix.parsers.vcf import (
        DEFAULT_MIN_GQ,
        MultiSampleError,
        SampleNotFoundError,
        VcfParser,
    )

    # Every call gets its own parser instance. The registry holds one shared
    # instance per format, and a parser keeps per-parse state (the VCF
    # quality stats, watched loci); two files read by one command (compare)
    # must not see each other's.
    if isinstance(parser, VcfParser):
        bound = VcfParser(
            sample=sample,
            keep_filtered=keep_filtered,
            min_gq=DEFAULT_MIN_GQ if min_gq is None else min_gq,
        )
        try:
            bound.validate_sample(file_path)
        except (MultiSampleError, SampleNotFoundError) as exc:
            raise click.ClickException(str(exc)) from exc
        return _mark_forced(bound, warning)

    return _mark_forced(type(parser)(), warning)


FORCED_FORMAT_CODE = "forced_format"


def _forced_format_warning(file_path: Path, parser: GenotypeParser) -> str | None:
    """The warning for ``--format`` on a file the reader does not recognise, else None."""
    with quiet_input_notices():
        recognised = parser.can_parse(file_path)
    if recognised:
        return None
    name = parser.display_name
    return (
        f"{file_path.name!r} doesn't look like a {name} file. Reading it as {name} because "
        f"you chose --format {parser.name}. Allelix can't check its layout, strand or build "
        "the usual way, so the results may be wrong. Use the original file from your "
        "testing company if you can."
    )


def _mark_forced(parser: GenotypeParser, warning: str | None) -> GenotypeParser:
    """Print the forced-format warning and keep it on the parser for the reports."""
    if warning is not None:
        console.print(f"[bold yellow]Warning:[/bold yellow] {escape(warning)}")
        parser.forced_warning = warning
    return parser


def forced_warning(parser: GenotypeParser) -> ReportWarning | None:
    """The ``forced_format`` report warning a forced, unrecognised file carries."""
    message = parser.forced_warning
    return ReportWarning(FORCED_FORMAT_CODE, message) if message else None


def _ready_annotators(
    data_dir: Path | None,
    *,
    include_benign: bool = False,
    gwas_filter_traits: bool = True,
    cadd_full: bool = False,
) -> tuple[Path, list[Annotator], list[Annotator], list[Annotator]]:
    """The data directory, the ready and not-ready annotators, and every annotator in order."""
    resolved = resolve_data_dir(data_dir)
    # ClinVar is ready when one build's cache is: a user who installed
    # only GRCh37 can analyze GRCh37 input. Input on the other build
    # finds ClinVar skipped once the build is decided.
    ready_builds = tuple(
        build
        for build in CLINVAR_SUPPORTED_BUILDS
        if ClinVarAnnotator(resolved, builds=(build,)).cache_status().ready
    )
    annotators = get_annotators(
        resolved,
        include_benign=include_benign,
        gwas_filter_traits=gwas_filter_traits,
        cadd_full=cadd_full,
        clinvar_builds=ready_builds or CLINVAR_SUPPORTED_BUILDS,
    )
    ready: list[Annotator] = []
    not_ready: list[Annotator] = []
    for a in annotators:
        if a.is_ready():
            ready.append(a)
        else:
            not_ready.append(a)
    if not ready:
        reasons = "\n".join(f"  {unready_line(a, a.cache_status())}" for a in annotators)
        raise click.ClickException(
            f"No annotators are ready. Run `allelix db update` first.\n{reasons}"
        )
    return resolved, ready, not_ready, annotators


def unready_line(annotator: Annotator, status: CacheStatus) -> str:
    """One console line: why a source's cache cannot be used, and the fix (ADR-0048)."""
    from allelix.databases.cache_state import CacheState

    name = annotator.display_name
    where = status.path if status.path is not None else annotator.data_dir
    detail = f" ({status.detail})" if status.detail else ""
    if status.state is CacheState.OUTDATED:
        return f"{name}: out of date{detail}. Run `allelix db update`."
    if status.state is CacheState.UNREADABLE:
        return f"{name}: cannot read {where}{detail}. Check the permissions on {where}."
    if status.state is CacheState.CORRUPT:
        return f"{name}: {where} is damaged{detail}. Run `allelix db update`."
    if status.state is CacheState.WAL:
        return (
            f"{name}: {where} is in WAL mode from an older Allelix. Run `allelix db update` "
            "once to convert it (converting it downloads nothing)."
        )
    if status.state is CacheState.READONLY_WAL:
        return (
            f"{name}: {where} is in WAL mode, which cannot be read in a data directory "
            "you cannot write. Run `allelix db update` once with write access to "
            f"{where.parent}."
        )
    return f"{name}: not installed. Run `allelix db update`."


_STALENESS_SECONDS = 7 * 24 * 60 * 60  # 7 days


class _ResolvedStream:
    """Variants from ``_resolved_variants``, the build decision, and the high-value watch results.

    When contig lengths decided the build or stand beside ``--build``, the
    known-SNP positions are checked after the last variant, so a consumer must exhaust
    ``variants`` (``must_exhaust``) before it writes or prints a result.
    """

    def __init__(
        self,
        variants: Iterator[Variant],
        parser: GenotypeParser,
        placements: dict[str, list[tuple[str, int, str, str, int, str, str]]] | None,
        high_value: set[str],
        decision: BuildDecision,
        record: StreamBuild,
    ) -> None:
        self.variants = variants
        self.decision = decision
        self._record = record
        self._parser = parser
        self._placements = placements
        self._high_value = high_value
        self._attribution: dict[str, tuple[Placements, dict[tuple[str, int], frozenset[str]]]] = {}

    def build_diagnostics(self) -> BuildDiagnostics:
        """The build decision with what the stream learnt; call after exhausting it."""
        return self._record.diagnostics()

    @property
    def must_exhaust(self) -> bool:
        """True when stopping early would skip the late known-SNP check.

        Contig lengths deciding the build, or confirming a ``--build``, both
        leave the known-SNP check to the end of the pass (``needs_late_seeds``).
        """
        return needs_late_seeds(self.decision)

    def high_value_targets(self, variant: Variant) -> set[str]:
        """The high-value rsIDs a streamed call stands for (``allelix.attribution``).

        VCF input (REF from the file): a call at one of an rsID's placements on
        its build with that placement's REF, whatever its ID; by ID when the
        rsID has no placement. Other input: by ID, as before.
        """
        from allelix.attribution import attributed_rsids, build_placements, placement_index

        if not variant.ref_from_input or not self._placements:
            return {variant.rsid} if variant.rsid in self._high_value else set()
        cached = self._attribution.get(variant.build)
        if cached is None:
            places = build_placements(self._placements, variant.build)
            cached = (places, placement_index(places))
            self._attribution[variant.build] = cached
        places, index = cached
        return attributed_rsids(
            self._high_value,
            variant.chromosome,
            variant.position,
            variant.ref,
            variant.rsid,
            places,
            index,
        )

    def high_value_states(self, build: str) -> dict[str, tuple[str, str]] | None:
        """VCF input: what the file says at each high-value position.

        Call after ``variants`` is exhausted, with the effective build.
        None for non-VCF input or when no high-value set was watched.
        """
        from allelix.reports._pipeline import high_value_states

        return high_value_states(
            getattr(self._parser, "quality_stats", None),
            self._placements,
            self._high_value,
            build,
        )


@contextlib.contextmanager
def _resolved_variants(
    parser: GenotypeParser,
    file_path: Path,
    *,
    data_dir: Path | None,
    build_override: str | None = None,
    high_value: set[str] | None = None,
) -> Iterator[_ResolvedStream]:
    """Parse ``file_path`` with the decided build set and rsID-less calls resolved.

    This is the command-level form of ``_flush()``'s resolution step, for
    ``stats``, ``extract`` and ``export plink``. The build is decided first
    with the same rule as ``analyze`` (``decide_file_build``, ADR-0044;
    resolution needs the build); a run that must stop raises
    ``click.ClickException`` before any variant is yielded, or, for the
    late known-SNP check, when the stream ends. Then ``resolve_rsidless`` applies the
    shared ClinVar / identity-table / gnomAD chain. Sources whose caches
    are not installed are skipped, so without databases the variants pass
    through unchanged.

    ``high_value``: rsIDs whose positions a VCF parser should watch, as
    ``analyze`` does, so rsID-less no-call records and gVCF blocks at those
    positions are classified (``_ResolvedStream.high_value_states``).

    This context owns only the parser/build/resolution/output iterators it
    creates. They close before the providers even if the returned stream is
    retained or a consumer exits early; remaining rows are not consumed.
    """
    from allelix.annotators._resolution import resolve_rsidless
    from allelix.annotators.clinvar import ClinVarAnnotator
    from allelix.annotators.gnomad import GnomadAnnotator
    from allelix.annotators.rsid_identity import RsidIdentityAnnotator
    from allelix.parsers.vcf import VcfParser
    from allelix.reports._pipeline import (
        StreamBuild,
        _cleanup_on_exit,
        _owned_annotators,
        decide_file_build,
        watch_high_value,
        with_effective_build,
    )

    try:
        decision = decide_file_build(
            parser, file_path, parser.get_metadata(file_path), override=build_override
        )
    except (BuildDecisionError, UnusableInputError) as exc:
        raise click.ClickException(str(exc)) from exc
    resolved = resolve_data_dir(data_dir)
    clinvar = ClinVarAnnotator(resolved)
    gnomad = GnomadAnnotator(resolved)
    identity = RsidIdentityAnnotator(resolved)
    hv = high_value or set()
    with _owned_annotators((clinvar, gnomad, identity)), _cleanup_on_exit() as cleanup:
        if hv and isinstance(parser, VcfParser):
            cleanup.callback(parser.watch, set(), set())
        placements = watch_high_value(parser, hv, identity)
        record = StreamBuild(decision, parser=parser)

        def own_iterator(iterator: Iterator[Variant]) -> Iterator[Variant]:
            close = getattr(iterator, "close", None)
            if callable(close):
                cleanup.callback(close)
            return iterator

        parsed = own_iterator(parser.parse(file_path))
        stream = own_iterator(
            with_effective_build(
                parsed, decision=decision, file_name=file_path.name, record=record
            )
        )
        variants = own_iterator(
            resolve_rsidless(stream, clinvar=clinvar, gnomad=gnomad, identity=identity)
        )
        output = own_iterator(_stop_as_click(variants, parser, file_path, record))
        yield _ResolvedStream(
            output,
            parser,
            placements,
            hv,
            decision,
            record,
        )


def _stop_as_click(
    variants: Iterator[Variant],
    parser: GenotypeParser | None = None,
    file_path: Path | None = None,
    record: StreamBuild | None = None,
) -> Iterator[Variant]:
    """Re-raise a late ``BuildDecisionError`` as the CLI's error type.

    When ``parser`` and ``file_path`` are given and the stream yields
    nothing, a structural-only VCF or a file with no genotype records is
    reported with its own error after the stream ends (the same check
    ``analyze`` makes), so ``stats``, ``extract`` and ``export plink`` never
    show an empty table or write an empty PLINK set, whatever decided the
    build. MT records dropped for an hg19-length chrM (``record``) are
    records: such a file is not empty, and gets the chrM note instead.
    """
    count = 0
    try:
        for variant in variants:
            count += 1
            yield variant
    except BuildDecisionError as exc:
        raise click.ClickException(str(exc)) from exc
    if parser is not None and file_path is not None:
        dropped = record.mt_skipped + record.ucsc_mito_dropped if record is not None else 0
        _reject_unusable(parser, file_path, count + dropped)


def _reject_unusable(parser: GenotypeParser, file_path: Path, records: int) -> None:
    """Report an unusable file (``reject_unusable_input``) as the CLI's error.

    A structural-only VCF or a file with no genotype records, whatever
    decided the build.
    """
    from allelix.reports._pipeline import reject_unusable_input

    try:
        reject_unusable_input(parser, file_path, records)
    except UnusableInputError as exc:
        raise click.ClickException(str(exc)) from exc


def _run_setup(annotator: Annotator, *, force: bool = False) -> bool:
    """Invoke annotator.setup(). Returns True on success, False on failure.

    ``force`` reaches the sources whose setup can refresh without a freshness
    signal (GWAS Catalog, ClinPGx: no checksum depends on it). ClinVar's
    signal is the checksum that verifies its download, so it never can.
    """
    from allelix.annotators.gwas import GWASCatalogAnnotator
    from allelix.annotators.pharmgkb import PharmGKBAnnotator

    without_signal = False
    try:
        if force and isinstance(annotator, (GWASCatalogAnnotator, PharmGKBAnnotator)):
            without_signal = True
            annotator.setup(force=True)
        else:
            annotator.setup()
    except Exception as exc:
        if hasattr(exc, "close"):
            exc.close()
        console.print(f"  [red]{annotator.name}: {exc}[/red]")
        return False
    sig = getattr(annotator, "cached_remote_signal", lambda: None)()
    # No freshness signal stored: none at all, or (ClinPGx) only the CPIC mark.
    no_pgkb_part = isinstance(annotator, PharmGKBAnnotator) and not (sig or "").startswith("pgkb:")
    if without_signal and (sig is None or no_pgkb_part):
        console.print(
            f"  [yellow]{annotator.name}: refreshed without a freshness signal; none is "
            "stored, so the next `db update` refreshes again once it is available.[/yellow]"
        )
    if sig and "cpic:unavailable" in sig:
        console.print(
            f"  [yellow]{annotator.name}: updated (CPIC unavailable — "
            "non-finding filter degraded, retry later)[/yellow]"
        )
    return True


def _prepare_caches(annotators: list[Annotator], data_dir: Path) -> None:
    """The update step's local work: no network, may write (ADR-0048).

    Converts every WAL cache in the data directory to rollback-journal
    mode, so a later read-only analysis can read it, then lets each
    outdated cache bring itself up to date from local data (a legacy
    stamp, a GWAS re-ingest from the cached TSV, a SNPedia parse).
    """
    from allelix.databases.cache_state import CacheState
    from allelix.databases.loader_utils import convert_to_rollback_journal

    for db_path in sorted(data_dir.glob("*.sqlite")):
        if convert_to_rollback_journal(db_path):
            console.print(f"  [dim]{db_path.name}: converted to rollback-journal mode[/dim]")
    for annotator in annotators:
        with annotator:
            if annotator.cache_status().state is CacheState.OUTDATED:
                annotator.prepare_cache()


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
    else is a gene entry (comma/semicolon lists are split). A ``#`` starts
    a comment, on its own line or after
    an entry (``rs1801133  # MTHFR``), and the rest of the line is ignored;
    rsIDs and gene names never contain ``#``. Blank lines are ignored.
    Gene names starting with ``RS`` (e.g., RSPO1, RSF1) are correctly
    classified as gene names, not rsIDs.

    Input is case-tolerant; output is canonical: rsIDs are normalized to
    lowercase (``rs1801133``), gene names to uppercase (``MTHFR``). The
    filter recorded in JSON output therefore looks identical regardless
    of how the user typed the entries in the filter file.
    """
    from allelix.reports._pipeline import gene_tokens

    genes: set[str] = set()
    rsids: set[str] = set()
    # utf-8-sig: a list saved with a byte order mark keeps its first entry.
    for raw in path.read_text(encoding="utf-8-sig").splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if _RSID_PATTERN.match(line):
            rsids.add(line.lower())
        else:
            genes.update(gene_tokens(line))
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
    if v == "grch36":
        return "GRCh36"
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


def _source_choices(
    registry: list[Annotator],
    cfg: AllelixConfig,
    *,
    exclude_sources: frozenset[str] | None,
    flags: dict[str, bool],
) -> dict[str, str]:
    """Each source the user chose not to run, with the reason (in the order applied)."""
    from allelix.annotators.base import Permission

    classes = {type(a).name: type(a) for a in registry}
    choices: dict[str, str] = {}
    for annotator in registry:
        name = annotator.name
        if not cfg.is_enabled(name, classes):
            permission = cfg.permission_for(name, classes)
            denied = permission is not None and permission is not Permission.ALLOW
            choices[name] = "license" if denied else "disabled"
        elif exclude_sources and name in exclude_sources:
            choices[name] = "excluded"
        elif flags.get(name):
            choices[name] = "flag"
    return choices


def _emit_build_diagnostics(result: object) -> None:
    """``analyze``: the build banner and notices for the run's decision (ADR-0044)."""
    diag = getattr(result, "build_diagnostics", None)
    if diag is not None:
        _emit_build_notices(diag)


def _emit_build_notices(
    diag: BuildDiagnostics, *, heading: str | None = None, banner: bool = True
) -> None:
    """Every command's build banner and notices, from one place (``build_text``).

    ``heading`` names the file when a command reads more than one
    (``compare``); ``banner=False`` when the command shows the build itself
    (the ``stats`` table).
    """
    from allelix.reports.build_text import build_notices, source_text

    prefix = f"{heading}: " if heading else ""
    if banner:
        console.print(f"[dim]{prefix}Build: {diag.effective_build} ({source_text(diag)})[/dim]")
    for notice in build_notices(diag):
        console.print(f"[yellow]{prefix}{notice.text}[/yellow]")


def _emit_runtime_nudges(result: object, missing_resolvers: tuple[str, ...] = ()) -> None:
    """Surface runtime-degradation notes the analyze run can detect.

    Both are once-per-run, post-pipeline. Neither blocks the report; the
    intent is transparency for cases where the silent default behavior
    used to leave the user uncertain about WHY their count looked low.

    - **strand-aware matching inactive or limited.** When
      carrier/genotype matching ran with ``Variant.ref`` unavailable
      across the input, the strand-flip path in ``utils/allele.py``
      falls back to the v2.0.x direct-match-only behavior. On GRCh37
      and GRCh38 arrays this means the rsID identity table (ADR-0040)
      and gnomAD are both missing; on other builds no reference source
      exists. VCF inputs always carry REF from the file, so this never
      fires for them.
    - **ClinVar rows scored down or hidden** (ADR-0045). How many
      Pathogenic or Likely pathogenic rows left the default report because
      ClinVar's own classification lowered them, and how many rows carry
      a classification Allelix does not score, and how many zero-star
      Pathogenic or Likely pathogenic rows on a common allele left it
      (ADR-0052).
    - **rsID-less GRCh38 input without the resolvers.** When the
      effective build is GRCh38, most variants entered the pipeline
      without rsIDs (the variant-caller ``ID=.`` case), and a resolver
      named in ``missing_resolvers`` (the rsID identity table, gnomAD) is
      not ready, rsIDs are recovered from ClinVar positions alone and
      fewer findings appear. With both installed, coverage matches the
      GRCh37 equivalent, so the note is shown only when one is missing.
    """
    scored_down = scored_down_line(getattr(result, "clinvar_scored_down", 0))
    if scored_down:
        console.print(f"[dim]{scored_down}[/dim]")
        for ann in getattr(result, "clinvar_scored_down_variants", []):
            console.print(f"[dim]  {escape(scored_down_variant_text(ann))}[/dim]")
    common_allele = common_allele_line(getattr(result, "clinvar_common_allele", 0))
    if common_allele:
        console.print(f"[dim]{common_allele}[/dim]")
    unrecognised = getattr(result, "clinvar_unrecognised", 0)
    if unrecognised:
        console.print(
            f"[yellow]{unrecognised} ClinVar row(s) not shown: their classification "
            "is not one Allelix scores yet.[/yellow]"
        )
    for line in strand_abstention_lines(getattr(result, "strand_abstentions", {})):
        console.print(f"[dim]{escape(line)}[/dim]")

    no_ref = getattr(result, "no_ref_variant_count", 0)
    rsidless = getattr(result, "rsidless_variant_count", 0)
    total = getattr(result, "total_variants", 0)
    parser_name = getattr(result, "parser_name", "")
    diag = getattr(result, "build_diagnostics", None)
    effective_build = diag.effective_build if diag is not None else None

    # array inputs (everything except the VCF parser) running
    # without compatible gnomAD reference context leave every variant
    # ref=None. Gate
    # on "the whole input lacked ref" rather than "any variant lacked
    # ref" so a hybrid file (e.g. a few ref-less rows in an otherwise
    # populated VCF) doesn't flap the message on.
    if parser_name != "vcf" and total > 0 and no_ref == total:
        if effective_build in ("GRCh37", "GRCh38"):
            console.print(
                "[dim]strand-aware matching inactive — no reference context "
                "(run `allelix db update` for the rsID identity table and gnomAD "
                "reference resolution)[/dim]"
            )
        else:
            console.print(
                "[dim]strand-aware matching inactive — no reference context is "
                f"available for {effective_build or 'this'} input[/dim]"
            )

    # GRCh38 + a high rsID-less fraction. The threshold (>50%)
    # distinguishes a variant-caller VCF (mostly ID=. — fires) from a
    # 23andMe-style array (rsID-bearing — no warning) or a mostly-
    # rsID-bearing VCF where a handful of novel calls came in with
    # ID=. (no warning, ID=. for novel sites is expected).
    if missing_resolvers and effective_build == "GRCh38" and total > 0 and rsidless / total > 0.5:
        console.print(
            "[yellow]GRCh38 input without rsIDs: "
            f"{' and '.join(missing_resolvers)} not ready, so rsIDs are recovered "
            "from fewer sources and fewer findings are shown. Run `allelix db update` "
            "to install them.[/yellow]"
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
    vcf_out: Path | None = None,
    keep_filtered: bool = False,
    min_gq: int | None = None,
    filter_path: Path | None = None,
) -> None:
    from allelix.reports import check_writable, same_destination, validate_output_path

    # Before any work (2.5.2): no output may be a file this run reads (the
    # input, the --diff baseline, the --filter-file), and each must be writable.
    inputs = [path for path in (file_path, diff_path, filter_path) if path is not None]
    try:
        for destination, legacy_temp in ((output, True), (vcf_out, False)):
            if destination is None:
                continue
            for source in inputs:
                validate_output_path(source, destination, legacy_temp=legacy_temp)
        if output is not None and vcf_out is not None and same_destination(output, vcf_out):
            msg = (
                f"--output {output} and --vcf-out {vcf_out} name the same file; "
                "choose two different paths."
            )
            raise ValueError(msg)
        for destination in (output, vcf_out):
            if destination is not None:
                check_writable(destination)
    except (OSError, RuntimeError, ValueError) as exc:
        raise click.ClickException(str(exc)) from exc

    # The parser first: a --format warning is printed before anything else.
    parser = _resolve_parser(
        file_path, fmt, sample=sample, keep_filtered=keep_filtered, min_gq=min_gq
    )
    # --vcf-out preflight. Bail before any expensive setup so a
    # user who mistyped the flag or aimed it at a 23andMe file finds
    # out immediately, not after the analyze pipeline runs.
    if vcf_out is not None:
        if parser.name != "vcf":
            raise click.ClickException(
                f"--vcf-out requires VCF or gVCF input; got {parser.display_name} "
                f"({parser.name!r}). Convert to VCF first, or drop --vcf-out to "
                "produce the default report."
            )
        if vcf_out.exists():
            raise click.ClickException(
                f"--vcf-out refuses to overwrite an existing file: {vcf_out}. "
                "Move or delete the existing file, or choose a different path."
            )

    resolved = resolve_data_dir(data_dir)
    from allelix.config import load_config

    # An unreadable config.toml stops the run before the update step writes.
    cfg = load_config(resolved, create=not no_update)
    if not no_update:
        _prepare_caches(get_annotators(resolved), resolved)
        _maybe_refresh_databases(resolved)
    _, ready, not_ready, registry = _ready_annotators(
        data_dir,
        include_benign=include_benign,
        gwas_filter_traits=not gwas_all,
        cadd_full=cfg.cadd_full,
    )
    choices = _source_choices(
        registry,
        cfg,
        exclude_sources=exclude_sources,
        flags={"gnomad": no_gnomad, "alphamissense": no_alphamissense, "cadd": no_cadd},
    )
    ready = [a for a in ready if a.name not in choices]
    # gnomAD, AlphaMissense, CADD and the identity table run through
    # their own parameters, not the generic annotator loop.
    gnomad_annotator, ready = _pop_enrichment_annotator(ready, "gnomad", skip=False)
    am_annotator, ready = _pop_enrichment_annotator(ready, "alphamissense", skip=False)
    cadd_annotator, ready = _pop_enrichment_annotator(ready, "cadd", skip=False)
    identity_annotator, ready = _pop_enrichment_annotator(ready, "rsid_identity", skip=False)

    # A source the user turned off is not reported as unusable.
    unready = [(a, a.cache_status()) for a in not_ready if a.name not in choices]
    for annotator, status in unready:
        console.print(f"[yellow]Skipping {escape(unready_line(annotator, status))}[/yellow]")

    all_active: list[Annotator] = list(ready)
    if gnomad_annotator is not None and gnomad_annotator.is_ready():
        all_active.append(gnomad_annotator)
    if am_annotator is not None and am_annotator.is_ready():
        all_active.append(am_annotator)
    if cadd_annotator is not None and cadd_annotator.is_ready():
        all_active.append(cadd_annotator)
    if identity_annotator is not None and identity_annotator.is_ready():
        all_active.append(identity_annotator)
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
            identity=cast("RsidIdentityAnnotator | None", identity_annotator),
            plausibility=cfg.plausibility,
            high_value_rsids=hv_rsids,
            panel_rsids=rsids,
            collect_vcf_write_inputs=vcf_out is not None,
            resource_owners=registry,
        )
    except (BuildDecisionError, UnusableInputError) as exc:
        # ADR-0044: stop before any report, annotated VCF or diff is written.
        raise click.ClickException(str(exc)) from exc
    finally:
        _unwire_parser_logging(counter, stderr_handler, snapshot)
    failures = {a.name: (status.state.value, unready_line(a, status)) for a, status in unready}
    not_applicable = (
        [identity_annotator.name]
        if identity_annotator is not None and result.build not in _IDENTITY_BUILDS
        else []
    )
    for annotator in ready:
        if annotator.serves_build(result.build):
            continue
        # ClinVar is ready on the other build only.
        if result.build in CLINVAR_SUPPORTED_BUILDS:
            other = ClinVarAnnotator(resolved, builds=(result.build,))
            status = other.cache_status()
            failures[annotator.name] = (status.state.value, unready_line(other, status))
            console.print(
                f"[yellow]Skipping {escape(unready_line(other, status))}"
                f" (input is {escape(result.build)})[/yellow]"
            )
        else:
            not_applicable.append(annotator.name)
    result.sources = build_sources(
        [(a.name, a.display_name) for a in registry],
        choices=choices,
        failures=failures,
        used=dict(result.annotators_used),
        not_applicable=not_applicable,
    )
    # Notices about the input file; already printed as "warning:" lines.
    for code, message in counter.input_notices:
        result.report_warnings.append(ReportWarning(code, message))
    warning = clinvar_skipped_warning(result.sources)
    if warning is not None:
        result.report_warnings.insert(0, warning)
    forced = forced_warning(parser)
    if forced is not None:
        result.report_warnings.insert(0, forced)
    # --build forced a build the file's known-SNP positions match none of
    # (2.5.2): the build notice is already on the console; the report keeps it.
    from allelix.reports.build_text import forced_unmatched_text

    unmatched = (
        forced_unmatched_text(result.build_diagnostics)
        if result.build_diagnostics is not None
        else None
    )
    if unmatched is not None:
        result.report_warnings.insert(0, ReportWarning(BUILD_OVERRIDE_UNMATCHED, unmatched))
    unchecked = unchecked_sites(
        high_value,
        result.hv_variants,
        result.di_calls_at_plp,
        result.vcf_quality.reference_carried_not_in_file_rsids
        if result.vcf_quality is not None
        else [],
    )
    di_unchecked = [site for site in unchecked if site.reason == UNCHECKED_DI_CALL]
    if di_unchecked:
        result.report_warnings.append(
            ReportWarning(
                DI_CALLS_NOT_CHECKED,
                f"{len(di_unchecked)} call(s) at high-value or Pathogenic sites are written "
                "as D or I and could not be checked: "
                + "; ".join(site.text for site in di_unchecked),
                len(di_unchecked),
            )
        )
    for report_warning in result.report_warnings:
        if report_warning.code in _LOGGED_INPUT_NOTICES or report_warning.code == (
            BUILD_OVERRIDE_UNMATCHED
        ):
            continue  # already printed: a "warning:" line, or the build notice
        if report_warning.code == DI_CALLS_NOT_CHECKED:
            continue  # named in the high-value block below
        console.print(f"[bold yellow]Warning:[/bold yellow] {escape(report_warning.message)}")

    missing_resolvers = tuple(
        label
        for label, source in (
            ("the rsID identity table", identity_annotator),
            ("gnomAD", gnomad_annotator),
        )
        if source is None or not source.is_ready()
    )
    _reject_unusable(parser, file_path, result.total_variants)
    _emit_build_diagnostics(result)
    _emit_runtime_nudges(result, missing_resolvers)

    hv_warnings = scan_no_calls(result.hv_variants, high_value, result.high_value_vcf_states)
    if hv_warnings:
        console.print(
            f"[bold red]Warning:[/bold red] {len(hv_warnings)} high-value SNP(s) returned no-call:"
        )
        for line in format_warnings(hv_warnings):
            console.print(f"  [red]⚠[/red] {line}")
    if unchecked:
        console.print(
            f"[bold red]Warning:[/bold red] {len(unchecked)} high-value or Pathogenic "
            "site(s) could not be checked:"
        )
        for site in unchecked:
            console.print(f"  [red]⚠[/red] {escape(site.text)}")

    if counter.count:
        console.print(
            f"[yellow]Note:[/yellow] {counter.count:,} malformed line(s) skipped "
            "(see warnings on stderr)."
        )
    for line in vcf_quality_lines(result.vcf_quality) + array_quality_lines(
        result.array_quality, result.probe_classes
    ):
        console.print(f"[yellow]Note:[/yellow] {line}")
    if result.input_method:
        console.print(
            f"[yellow]Note:[/yellow] input method: {input_method_label(result.input_method)}."
            + (
                " Rare variants are the least reliable imputed calls; confirm any finding"
                " by direct genotyping."
                if is_imputed(result.input_method)
                else ""
            )
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
        # Below the display floor means the floors exclude the finding,
        # before the GWAS rollup hides a PheCode sibling (2.5.2).
        above_floor = {id(a) for a in filtered_for_diff}
        filtered_for_diff = rollup_gwas_duplicates(filtered_for_diff)
        # Removal means absent from the current findings in scope, before any
        # display floor: a finding scored below the floor is changed. The pool
        # is the displayed findings, rolled up as the report shows them, plus
        # the in-scope findings the report does not display; a hidden GWAS
        # PheCode sibling never replaces a displayed parent.
        displayed = {id(a) for a in filtered_for_diff}
        undisplayed = [
            a
            for a in result.filter(min_magnitude=0.0, category=category, genes=genes, rsids=rsids)
            if id(a) not in displayed
        ]
        diff_result = compute_diff(
            filtered_for_diff,
            prev["annotations"],
            prev.get("generated_at", ""),
            unfiltered=filtered_for_diff + undisplayed,
            below_floor={id(a) for a in undisplayed if id(a) not in above_floor},
        )
        previous_input = (prev.get("input") or {}).get("file")
        if previous_input and previous_input != file_path.name:
            console.print(
                f"[bold yellow]Warning:[/bold yellow] the --diff baseline was made from "
                f"{escape(previous_input)}, not {escape(file_path.name)}; the comparison "
                "may be between two different people or files."
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
                [
                    {
                        "rsid": w.snp.rsid,
                        "gene": w.snp.gene,
                        "note": w.snp.note,
                        "reason": w.reason,
                    }
                    for w in hv_warnings
                ]
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
                high_value_unchecked=[site.to_dict() for site in unchecked] or None,
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
                high_value_unchecked=[site.text for site in unchecked] or None,
            )
        console.print(f"[green]Wrote {rendered:,} annotation(s) to {output}[/green]")

    if vcf_out is not None:
        from allelix.reports.vcf import render_vcf

        # Pipeline populates result.vcf_write_inputs whenever the
        # ``collect_vcf_write_inputs`` flag is set — which the preflight
        # above guarantees for this branch — so the None case here is
        # a defensive assert, not a runtime path.
        assert result.vcf_write_inputs is not None
        stamped = render_vcf(
            input_path=file_path,
            output_path=vcf_out,
            build=result.build,
            build_source=(
                result.build_diagnostics.source if result.build_diagnostics is not None else None
            ),
            annotators_used=result.annotators_used,
            skipped_sources=[
                (s.name, s.reason or "") for s in result.sources if s.state == "skipped"
            ],
            annotations_by_key=result.vcf_write_inputs.annotations_by_key,
            recovered_rsids=result.vcf_write_inputs.recovered_rsids,
            keep_filtered=keep_filtered,
            reference_rows=result.vcf_write_inputs.reference_rows,
            record_rows=result.vcf_write_inputs.record_rows,
            record_recovered=result.vcf_write_inputs.record_recovered,
        )
        console.print(
            f"[green]Wrote annotated VCF with {stamped:,} stamped row(s) to {vcf_out}[/green]"
        )

    console.print(
        f"[dim]{len(result.annotations):,} total annotation(s) from {len(ready)} "
        f"database(s) across {result.total_variants:,} variant(s).[/dim]"
    )

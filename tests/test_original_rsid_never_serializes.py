# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Serialization-guard test for ``Variant.original_rsid`` (GH #128).

``original_rsid`` is a debug-only audit field stamped by the rsID-less
resolver paths in ``_pipeline._flush()`` and ``cli.utility``. It must
never reach JSON, HTML, or ``--diff`` output — leaking it would couple
production behavior to provenance metadata and turn an internal field
into an unintended API surface, exactly the failure mode the field's
docstring warns against.

This test exercises every renderer path that could leak the field:

  - ``render_json`` (JSON report)
  - ``render_html`` (HTML report)
  - ``render_diff`` (``--diff`` JSON output)
  - ``render_terminal`` (Rich terminal report; reads ``hv_variants``
    for the high-value no-call banner)
  - ``_render_extract_table`` (``allelix extract`` output; this PR
    stamps ``original_rsid`` directly on Variants returned from
    ``_sequential_extract``, so the extract-side renderer is the most
    immediately exposed surface)

For each path it builds a Variant with ``original_rsid`` set to a
sentinel string, runs the renderer, and asserts the sentinel does not
appear in the rendered output.

If a renderer change in a future PR accidentally includes the field —
e.g. by adding ``asdict(variant)`` somewhere that serializes — this
test fails immediately at PR time. Adding new renderer paths to this
file is the only ongoing maintenance burden ``original_rsid`` carries.
"""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from pathlib import Path

from allelix.models import Annotation, Variant
from allelix.reports._pipeline import AnalysisResult, BuildDiagnostics
from allelix.reports.html import render_html
from allelix.reports.json_report import render_json

_SENTINEL = "chr1:11796321:G:A__ORIGINAL_RSID_LEAK_CANARY__"


def _variant_with_sentinel() -> Variant:
    """Construct a Variant carrying the sentinel in original_rsid."""
    v = Variant(
        rsid="rs1801133",
        chromosome="1",
        position=11796321,
        allele1="G",
        allele2="A",
        build="GRCh38",
        ref="G",
    )
    v.original_rsid = _SENTINEL
    return v


def _result_with_sentinel(file_path: Path) -> AnalysisResult:
    """Build an AnalysisResult whose hv_variants carry the sentinel."""
    v = _variant_with_sentinel()
    annotation = Annotation(
        source="clinvar",
        rsid="rs1801133",
        significance="clinvar_drug_response",
        category="clinical",
        magnitude=6.5,
        description="Test annotation for the GH #128 serialization guard.",
        attribution="ClinVar",
        genotype_match="AG",
        condition="MTHFR deficiency",
        gene="MTHFR",
    )
    return AnalysisResult(
        file_path=file_path,
        parser_name="vcf",
        parser_display_name="VCF / gVCF",
        sample_id="HG002",
        build="GRCh38",
        total_variants=1,
        skipped_count=0,
        annotators_used=[("clinvar", "2026-06-18")],
        annotations=[annotation],
        build_diagnostics=BuildDiagnostics(
            header_build="GRCh38",
            detected_build="GRCh38",
            effective_build="GRCh38",
            override=False,
            matched_count=1,
            inspected_count=1,
            chr_prefix_inferred=False,
        ),
        hv_variants=[v],
        panel_rsids=frozenset({"rs1801133"}),
        genotyped_panel_rsids=frozenset({"rs1801133"}),
        panel_annotated_rsids=frozenset({"rs1801133"}),
        panel_genotypes={"rs1801133": "G/A"},
    )


@pytest.fixture
def file_path(tmp_path: Path) -> Path:
    """A throwaway path label for the AnalysisResult (no file I/O needed)."""
    return tmp_path / "test_input.vcf.gz"


class TestOriginalRsidDoesNotLeak:
    """Every renderer must scrub ``original_rsid`` from its output."""

    def test_field_string_never_renders_in_json(self, tmp_path: Path, file_path: Path) -> None:
        """JSON report must not contain ``original_rsid`` key or sentinel value."""
        result = _result_with_sentinel(file_path)
        out = tmp_path / "report.json"
        render_json(result, output_path=out)
        as_text = out.read_text()
        # Defensive sanity — the JSON must be parseable.
        json.loads(as_text)
        assert "original_rsid" not in as_text, (
            "JSON report exposes original_rsid — production logic could now "
            "branch on a debug-only audit field (GH #128)."
        )
        assert _SENTINEL not in as_text, (
            "JSON report leaks the sentinel originally stashed in "
            "Variant.original_rsid; some serializer is including the field."
        )

    def test_field_string_never_renders_in_html(self, tmp_path: Path, file_path: Path) -> None:
        """HTML report must not contain the sentinel anywhere in markup."""
        result = _result_with_sentinel(file_path)
        out = tmp_path / "report.html"
        render_html(result, output_path=out)
        as_text = out.read_text()
        assert "original_rsid" not in as_text, (
            "HTML report exposes original_rsid in markup — debug-only field "
            "leaked into user-visible output (GH #128)."
        )
        assert _SENTINEL not in as_text, (
            "HTML report leaks the sentinel originally stashed in "
            "Variant.original_rsid; some renderer interpolated the field."
        )

    def test_field_string_never_renders_in_diff(self, tmp_path: Path, file_path: Path) -> None:
        """JSON ``--diff`` output must not surface original_rsid changes.

        Builds a baseline JSON, mutates the in-memory result, renders the
        current JSON with ``diff=`` pointed at the baseline, and asserts
        the resulting diff payload contains neither the field name nor
        the sentinel.
        """
        from allelix.reports.diff import compute_diff, load_previous_report

        prior_result = _result_with_sentinel(file_path)
        baseline_path = tmp_path / "baseline.json"
        render_json(prior_result, output_path=baseline_path)
        baseline_payload = load_previous_report(baseline_path)

        current_result = _result_with_sentinel(file_path)
        # Touch the sentinel so the diff has a reason to flag the field
        # if it leaks. Vary the annotation description so the diff has
        # at least one real change to compute against.
        current_result.hv_variants[0].original_rsid = _SENTINEL + "_CURRENT"
        current_result.annotations[0] = Annotation(
            source=current_result.annotations[0].source,
            rsid=current_result.annotations[0].rsid,
            significance=current_result.annotations[0].significance,
            category=current_result.annotations[0].category,
            magnitude=current_result.annotations[0].magnitude,
            description="Test annotation -- modified for diff (GH #128 guard).",
            attribution=current_result.annotations[0].attribution,
            genotype_match=current_result.annotations[0].genotype_match,
            condition=current_result.annotations[0].condition,
            gene=current_result.annotations[0].gene,
        )
        current_path = tmp_path / "current.json"
        diff = compute_diff(
            current_result.annotations,
            baseline_payload.get("annotations", []),
            baseline_payload.get("generated_at", ""),
        )
        render_json(current_result, output_path=current_path, diff=diff)
        as_text = current_path.read_text()
        assert "original_rsid" not in as_text, (
            "Diff output exposes original_rsid — a debug-only field is now "
            "visible as a tracked change between report runs (GH #128)."
        )
        assert _SENTINEL not in as_text, (
            "Diff output leaks the sentinel stashed in Variant.original_rsid."
        )

    def test_field_string_never_renders_in_terminal(self, file_path: Path) -> None:
        """The Rich terminal renderer reads ``result.hv_variants`` (for
        the high-value no-call banner). That path could ``str(variant)``
        a Variant in a future change — if it ever does, the field would
        leak into terminal output and this test fails. Capture via
        Rich's ``record=True`` mode and ``export_text()``.
        """
        from io import StringIO

        from rich.console import Console

        from allelix.reports.terminal import render_terminal

        result = _result_with_sentinel(file_path)
        # ``record=True`` captures every printable into an internal
        # buffer that ``export_text()`` returns verbatim — including
        # everything Rich would have written to the terminal.
        console = Console(file=StringIO(), force_terminal=False, record=True, width=120)
        render_terminal(result, console)
        as_text = console.export_text()
        assert "original_rsid" not in as_text, (
            "Terminal report exposes original_rsid — debug-only audit "
            "field reached user-visible CLI output (GH #128)."
        )
        assert _SENTINEL not in as_text, (
            "Terminal report leaks the sentinel originally stashed in "
            "Variant.original_rsid; some path is calling str() / repr() "
            "on a Variant in a way that surfaces the field."
        )

    def test_field_string_never_renders_in_extract_table(self) -> None:
        """``allelix extract`` is the most immediately exposed surface
        for this field — ``_sequential_extract`` and
        ``_execute_tabix_extract`` BOTH stamp ``original_rsid`` directly
        on the Variants they return, and ``_render_extract_table``
        renders those Variants via Rich. A future "show provenance in
        a new column" change here would silently turn the audit field
        into an API surface; this test fails the moment that happens.
        """
        from io import StringIO
        from pathlib import Path as _RuntimePath

        from rich.console import Console

        from allelix.cli import utility
        from allelix.cli.utility import _render_extract_table

        # Stamp the sentinel onto a Variant exactly the way the extract
        # paths do — the rsID was recovered from gnomAD, original_rsid
        # holds the pre-stamp pseudo-ID.
        variant = Variant(
            rsid="rs1801133",
            chromosome="1",
            position=11796321,
            allele1="G",
            allele2="A",
            build="GRCh38",
            ref="G",
        )
        variant.original_rsid = _SENTINEL

        recording = Console(file=StringIO(), force_terminal=False, record=True, width=120)
        # ``_render_extract_table`` uses the module-level ``console``
        # singleton; swap in the recording console for the duration of
        # the call so we capture its output, then restore.
        original_console = utility.console
        utility.console = recording
        try:
            _render_extract_table(
                _RuntimePath("HG002.child.g.vcf.gz"),
                {"rs1801133"},
                {"rs1801133": variant},
            )
        finally:
            utility.console = original_console
        as_text = recording.export_text()
        assert "original_rsid" not in as_text, (
            "Extract table renderer surfaces original_rsid — the debug-"
            "only field stamped by _sequential_extract / "
            "_execute_tabix_extract is now visible in CLI output (GH #128)."
        )
        assert _SENTINEL not in as_text, (
            "Extract table leaks the sentinel stashed in "
            "Variant.original_rsid; a column or formatter is reading the "
            "field directly."
        )

    def test_variant_repr_may_include_field_but_renderers_must_not(self) -> None:
        """``repr(Variant)`` is allowed to show the field — that's the
        debug surface the field exists for. This test pins the asymmetry
        so a future refactor doesn't conclude ``original_rsid`` shouldn't
        appear in ``repr`` either (which would defeat the audit purpose)."""
        v = _variant_with_sentinel()
        rendered = repr(v)
        assert _SENTINEL in rendered, (
            "Variant.__repr__ no longer shows original_rsid — if this was "
            "intentional, the field has lost its only justified read site. "
            "Reconsider whether the field still earns its keep."
        )

"""The HTML report's variant data sits in a ``<script type="application/json">`` block.

No database text may change how the browser reads that block: the text between the
block's open tag and the first ``</script`` holds no ``<`` at all, it is the whole
JSON and parses back to the exact text, the page has exactly its two script
elements (the data block and the controller), and the footer with the licence
attributions is still there. Checked per payload in each text field the block
carries (condition, description, gene, references), and once through the CLI with
a ClinVar condition ``Foo <!--<script> bar``.
"""

from __future__ import annotations

import contextlib
import json
import re
import sqlite3
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner
from test_annotation_parity_and_named_disclosures import _clinvar, _identity

from allelix.annotators.clinvar import ClinVarAnnotator, clinvar_db_filename
from allelix.cli import main
from allelix.models import Annotation
from allelix.reports._pipeline import AnalysisResult
from allelix.reports.html import render_html

if TYPE_CHECKING:
    from pathlib import Path

OPEN_TAG = '<script id="variant-data" type="application/json">'

PAYLOADS = {
    "ordinary": "Hereditary synthetic condition",
    "close-script": "a</script>b",
    "close-script-upper-space": "a</SCRIPT b",
    "comment-then-script": "Foo <!--<script> bar",
    "comment-open": "a<!--b",
    "comment-close": "a-->b",
    "open-script": "a<script>b",
    "backslash-then-lt": "a\\<b",
    "non-bmp-and-u2028": "a" + chr(0x1F9EC) + chr(0x2028) + "b",
}
FIELDS = ("condition", "description", "gene", "references")


def _annotation(field: str, payload: str) -> Annotation:
    values = {
        "condition": "Synthetic condition",
        "description": "Synthetic description",
        "gene": "SYNTH",
        "references": ["PMID:1"],
    }
    values[field] = [payload] if field == "references" else payload
    return Annotation(
        source="clinvar",
        rsid="rs900443001",
        significance="Pathogenic",
        category="disease",
        magnitude=8.0,
        attribution="ClinVar",
        genotype_match="AG",
        alt="A",
        **values,
    )


def _render(tmp_path: Path, annotation: Annotation) -> str:
    out = tmp_path / "report.html"
    result = AnalysisResult(
        file_path=tmp_path / "chip.txt",
        parser_name="23andme",
        parser_display_name="23andMe",
        sample_id="SYNTH443",
        build="GRCh37",
        total_variants=1,
        skipped_count=0,
        annotators_used=[("clinvar", "1")],
        annotations=[annotation],
    )
    render_html(result, output_path=out)
    return out.read_text(encoding="utf-8")


def _data_block(page: str) -> str:
    start = page.index(OPEN_TAG) + len(OPEN_TAG)
    end = page.lower().index("</script", start)
    return page[start:end]


def _assert_page_intact(page: str) -> list:
    block = _data_block(page)
    assert "<" not in block, block
    data = json.loads(block)  # the whole of the block, up to the first </script
    # exactly two script elements: the data block and the controller
    opens = [m.start() for m in re.finditer(r"<script\b", page, flags=re.IGNORECASE)]
    assert len(opens) == 2, len(opens)
    assert page.startswith(OPEN_TAG, opens[0]) or page.startswith(OPEN_TAG, opens[1])
    assert sum(page.startswith("<script>", at) for at in opens) == 1
    assert len(re.findall(r"</script\s*>", page, flags=re.IGNORECASE)) == 2
    # the footer and its licence attributions follow both scripts
    footer = page.rindex("<footer>")
    assert footer > max(opens)
    assert page.rstrip().endswith("</footer></body></html>")
    assert ClinVarAnnotator.license.attribution_text in page[footer:]
    return data


@pytest.mark.parametrize("field", FIELDS)
@pytest.mark.parametrize("payload", list(PAYLOADS.values()), ids=list(PAYLOADS))
def test_variant_data_holds_no_lt_and_round_trips_the_text(tmp_path, field, payload):
    page = _render(tmp_path, _annotation(field, payload))
    data = _assert_page_intact(page)
    (entry,) = data
    if field == "gene":
        assert entry["gene"] == payload
    else:
        (ann,) = entry["annotations"]
        assert ann[field] == ([payload] if field == "references" else payload)


# ---------------------------------------------------------------------------
# Through the CLI (outside reviewer A's case): a ClinVar condition
# "Foo <!--<script> bar" left one script element, no footer and no panel.
# ---------------------------------------------------------------------------

CONFIG = (
    "[sources]\nclinvar = true\npharmgkb = false\nsnpedia = false\ngwas = false\n"
    "gnomad = false\nalphamissense = false\ncadd = false\nrsid_identity = true\n"
)
P = 3_000_000


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
def test_cli_html_with_a_comment_and_script_in_a_clinvar_condition(tmp_path, build):
    data = tmp_path / "data"
    data.mkdir()
    (data / "config.toml").write_text(CONFIG)
    _clinvar(
        data,
        build,
        [("rs900443101", "1", P, "G", "A", "SYNTH", "Pathogenic", "Pathogenic", 1, 1)],
    )
    with contextlib.closing(sqlite3.connect(data / clinvar_db_filename(build))) as conn, conn:
        conn.execute("UPDATE clinvar_variants SET condition = ?", ("Foo <!--<script> bar",))
    _identity(data, [("rs900443199", "5", 1000, "G", "A", 1000, "G", "A", 0)])
    source = tmp_path / "chip.txt"
    source.write_text(
        "# This data file generated by 23andMe at: synthetic\n"
        f"# Reference human assembly {build}\n"
        "# rsid\tchromosome\tposition\tgenotype\n"
        f"rs900443101\t1\t{P}\tAG\n"
    )
    out = tmp_path / "report.html"
    result = CliRunner().invoke(
        main,
        [
            "analyze",
            str(source),
            "--data-dir",
            str(data),
            "--no-update",
            "--build",
            build.lower(),
            "--output",
            str(out),
        ],
    )
    assert result.exit_code == 0, (result.output, result.exception)
    page = out.read_text(encoding="utf-8")
    entries = _assert_page_intact(page)
    conditions = [a.get("condition") for e in entries for a in e["annotations"]]
    assert "Foo <!--<script> bar" in conditions, conditions

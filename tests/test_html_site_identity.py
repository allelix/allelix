"""HTML groups findings by their originating site, never by rsID alone.

Record zero is a real record; a record takes precedence over a trimmed locus.
Positionless lists keep the legacy grouping. All data are synthetic.
"""

import json
import re
from dataclasses import replace
from html.parser import HTMLParser
from pathlib import Path

import pytest
from test_html_variant_data_script import OPEN_TAG

from allelix.models import Annotation
from allelix.reports._pipeline import AnalysisResult
from allelix.reports.html import render_html

RS = "rs900529001"


def ann(**kw):
    values = dict(
        source="clinvar",
        rsid=RS,
        significance="clinvar_pathogenic",
        category="clinical",
        magnitude=8.0,
        description="Synthetic evidence",
        attribution="ClinVar",
        genotype_match="AG",
        gene="FIRST",
        condition="First condition",
        alt="A",
        locus=("1", 1000),
        record_index=0,
        allele_frequency=0.4,
        cadd_phred=1.0,
        am_pathogenicity=0.1,
        am_class="likely_benign",
        plausibility="plausible",
    )
    values.update(kw)
    return Annotation(**values)


def result(annotations, build="GRCh38"):
    return AnalysisResult(
        file_path=Path("synthetic.vcf"),
        parser_name="vcf",
        parser_display_name="VCF",
        sample_id="SYNTH",
        build=build,
        total_variants=2,
        skipped_count=0,
        annotators_used=[(s, "synthetic") for s in sorted({a.source for a in annotations})],
        annotations=annotations,
    )


class Rows(HTMLParser):
    def __init__(self):
        super().__init__()
        self.rows = []
        self.current = None

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "tr" and "data-row-id" in attrs:
            self.current = {"attrs": attrs, "text": ""}
            self.rows.append(self.current)

    def handle_endtag(self, tag):
        if tag == "tr":
            self.current = None

    def handle_data(self, text):
        if self.current is not None:
            self.current["text"] += text


def page(tmp_path, annotations, **kw):
    output = tmp_path / "report.html"
    render_html(result(annotations), output_path=output, **kw)
    text = output.read_text()
    start = text.index(OPEN_TAG) + len(OPEN_TAG)
    block = text[start : text.lower().index("</script", start)]
    assert "<" not in block
    data = json.loads(block)
    assert len(re.findall(r"<script\b", text, flags=re.I)) == 2
    assert len(re.findall(r"</script\s*>", text, flags=re.I)) == 2
    assert "<script>" in text
    assert "<footer>" in text and text.rstrip().endswith("</footer></body></html>")
    rows = Rows()
    rows.feed(text)
    assert len(rows.rows) == len(data)
    return text, data, rows.rows


def pair(mode, *, genotype="AG", named=True):
    first = ann(genotype_match=genotype, haploid=genotype == "A", rsid=RS if named else "")
    second = replace(
        first,
        locus=("1", 5000),
        record_index=1,
        gene="SECOND",
        condition="Second condition",
        allele_frequency=0.000001,
        cadd_phred=35.0,
        am_pathogenicity=0.99,
        am_class="likely_pathogenic",
        plausibility="suspicious",
    )
    if mode == "chip":
        first, second = replace(first, record_index=None), replace(second, record_index=None)
    if mode == "same-position-records":
        second = replace(second, locus=first.locus)
    if mode == "other-chromosome":
        second = replace(second, locus=("2", 1000))
    return first, second


@pytest.mark.parametrize("mode", ["vcf", "chip", "same-position-records", "other-chromosome"])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("genotype", ["AG", "AA", "A", "GA/A"])
def test_each_site_retains_its_metrics_sources_and_warning(tmp_path, mode, reverse, genotype):
    first, second = pair(mode, genotype=genotype)
    annotations = []
    for a in (first, second):
        # Two sources at each site: no summary or search text can leak across sites.
        annotations.extend(
            [
                a,
                replace(
                    a,
                    source="snpedia",
                    attribution="SNPedia",
                    significance="snpedia_bad",
                    magnitude=5.0,
                ),
            ]
        )
    if reverse:
        annotations.reverse()
    _, data, rows = page(tmp_path, annotations)
    assert len(data) == 2, "distinct sites were merged into one HTML group"
    for entry, row in zip(data, rows, strict=True):
        expected = first if entry["gene"] == "FIRST" else second
        assert entry["gene"] == expected.gene
        assert entry["allele_frequency"] == expected.allele_frequency
        assert entry["cadd_phred"] == expected.cadd_phred
        assert entry["am_pathogenicity"] == expected.am_pathogenicity
        assert entry["am_class"] == expected.am_class
        assert entry["genotype"] == genotype
        assert entry["zygosity"] == expected.zygosity
        assert entry["plausibility"] == expected.plausibility
        assert {a["source"] for a in entry["annotations"]} == {"clinvar", "snpedia"}
        assert {a["condition"] for a in entry["annotations"]} == {expected.condition}
        assert entry["allele_metrics"] == []
        assert "+1]" in row["text"] and "+3]" not in row["text"]
        other = second if expected is first else first
        assert other.condition not in row["text"]
        assert other.gene.lower() not in row["attrs"]["data-search-text"]
        assert row["attrs"]["data-plausibility"] == expected.plausibility
        assert ("⚠" in row["text"]) == (expected.plausibility == "suspicious")


@pytest.mark.parametrize("named", [True, False])
@pytest.mark.parametrize("reverse", [False, True])
def test_two_records_at_one_position_do_not_pool_same_alt_metrics(tmp_path, named, reverse):
    first, second = pair("same-position-records", named=named)
    # Both carry A, but originated in different records/REFs. ALT alone is no key.
    first = replace(first, locus_ref="G", record_key=("1", 1000, "G", "A"))
    second = replace(second, locus_ref="C", record_key=("1", 1000, "C", "A"))
    _, data, _ = page(tmp_path, [second, first] if reverse else [first, second])
    assert len(data) == 2, "record identity was lost at a shared position/ALT"
    assert sorted(e["cadd_phred"] for e in data) == [1.0, 35.0]


@pytest.mark.parametrize("reverse", [False, True])
def test_multi_allelic_metrics_are_local_to_the_record(tmp_path, reverse):
    annotations = []
    for i, pos in enumerate((1000, 5000)):
        for alt, score in (("A", 1.0 + 30 * i), ("T", 2.0 + 30 * i)):
            annotations.append(
                ann(
                    locus=("1", pos),
                    record_index=i,
                    genotype_match="AT",
                    alt=alt,
                    cadd_phred=score,
                    condition=f"site {pos}",
                )
            )
    if reverse:
        annotations.reverse()
    _, data, _ = page(tmp_path, annotations)
    assert len(data) == 2, "multi-allelic groups crossed record boundaries"
    got = sorted(sorted((m["allele"], m["cadd_phred"]) for m in e["allele_metrics"]) for e in data)
    assert got == [[("A", 1.0), ("T", 2.0)], [("A", 31.0), ("T", 32.0)]]


@pytest.mark.parametrize("repute,significance", [("bad", "snpedia_bad"), ("good", "snpedia_good")])
def test_cards_count_distinct_sites(tmp_path, repute, significance):
    first, second = pair("vcf")
    text, _, _ = page(
        tmp_path,
        [replace(a, source="snpedia", significance=significance) for a in (first, second)],
    )
    card = re.search(
        rf'<div class="card card-{repute}">.*?<div class="value">(.*?)</div>', text, re.S
    )
    assert card is not None
    assert card[1] == "2", "summary card counted merged loci once"


@pytest.mark.parametrize("kind", ["distinct-locus", "same-locus", "different-genotypes"])
def test_shared_label_has_position_in_panel_data_and_controller(tmp_path, kind):
    first, second = pair("same-position-records" if kind == "same-locus" else "vcf")
    if kind == "different-genotypes":
        second = replace(second, genotype_match="AA")
    text, data, _ = page(tmp_path, [first, second])
    assert [e.get("locus") for e in data] == [
        "1:1000",
        "1:1000" if kind == "same-locus" else "1:5000",
    ]
    # Data alone is not a user-visible Position field. The existing renderer uses fieldRow.
    assert 'fieldRow("Position",' in text or "fieldRow('Position'," in text


@pytest.mark.parametrize("mode", ["one-record", "trimmed-record", "chip", "positionless"])
def test_one_site_keeps_legacy_grouping_and_no_position_row(tmp_path, mode):
    first = ann()
    second = replace(first, source="snpedia", attribution="SNPedia", magnitude=5.0)
    if mode == "trimmed-record":
        first = replace(first, locus=("1", 1001), record_key=("1", 1000, "GA", "G"))
    if mode in {"chip", "positionless"}:
        first, second = replace(first, record_index=None), replace(second, record_index=None)
    if mode == "positionless":
        first, second = replace(first, locus=None), replace(second, locus=None)
    _, data, _ = page(tmp_path, [first, second])
    assert len(data) == 1
    assert len(data[0]["annotations"]) == 2
    assert not data[0].get("locus"), "ordinary panels must not gain a Position row"


def test_positionless_and_positioned_findings_are_not_one_site(tmp_path):
    first = ann(record_index=None, locus=None)
    second = ann(record_index=None, condition="Known locus")
    _, data, _ = page(tmp_path, [first, second])
    assert len(data) == 2, "unknown location is not evidence of the known location"


def test_filtering_a_locus_restores_the_ordinary_panel(tmp_path):
    first, second = pair("vcf")
    second = replace(second, magnitude=1.0)
    _, data, _ = page(tmp_path, [first, second], min_magnitude=4.0)
    assert len(data) == 1 and not data[0].get("locus")


def test_script_escape_survives_site_separation(tmp_path):
    first, second = pair("vcf")
    payload = "Foo <!--<script> bar \\< 🧬"
    text, data, _ = page(
        tmp_path, [replace(a, condition=payload, references=[payload]) for a in (first, second)]
    )
    assert len(data) == 2
    assert all(e["annotations"][0]["condition"] == payload for e in data)
    assert all(e["annotations"][0]["references"] == [payload] for e in data)
    assert "<footer>" in text


def test_panel_coverage_remains_rsid_based(tmp_path):
    first, second = pair("vcf")
    report = result([first, second])
    report.panel_rsids = frozenset({RS})
    report.genotyped_panel_rsids = frozenset({RS})
    report.panel_annotated_rsids = frozenset({RS})
    report.panel_genotypes = {RS: "A/G"}
    output = tmp_path / "panel.html"
    render_html(report, output_path=output)
    text = output.read_text()
    assert "1 requested: 1 annotated, 0 typed (no annotation), 0 no call, 0 not in file." in text
    assert re.search(r"<td>annotated</td><td>[^<]*</td><td>2</td>", text)

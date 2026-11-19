"""The position-first split-site pass gives exactly what the per-record pass gives.

2.7.1: a record is held back and the split-site component checks run
only when the next record has the same CHROM and POS. The per-record pass
of 2.7.0 stays as ``VcfParser._POSITION_FIRST = False``, the reference here:
on every committed VCF fixture and on the split shapes, both passes give the
same calls in the same order, with the same record numbers and counts.
"""

import gzip
from pathlib import Path

import pytest
from split_site_contract_helpers import BUILDS, OTHER, POS, RSID, record, write_vcf

from allelix.parsers.vcf import VcfParser

FIXTURES = Path(__file__).parent / "fixtures"
VCFS = sorted(p for p in FIXTURES.rglob("*") if p.name.endswith((".vcf", ".vcf.gz", ".vcf.bgz")))


def _samples(path):
    opener = gzip.open if path.suffix in (".gz", ".bgz") else open
    with opener(path, "rt", encoding="utf-8-sig", errors="replace") as stream:
        for line in stream:
            if line.startswith("#CHROM\t"):
                return line.rstrip("\n").split("\t")[9:]
    return []


def _read(path, position_first, monkeypatch, *, sample=None, drop_mito=False):
    monkeypatch.setattr(VcfParser, "_POSITION_FIRST", position_first)
    parser = VcfParser(sample=sample)
    parser.combine_split_sites = True
    if drop_mito:
        parser.drop_ucsc_mito = True
    calls = [
        (
            v.rsid,
            v.chromosome,
            v.position,
            v.allele1,
            v.allele2,
            v.ref,
            v.input_alt,
            v.haploid,
            v.record_index,
            v.allele_records,
        )
        for v in parser.parse(path)
    ]
    stats = parser.quality_stats
    counts = (
        parser.split_site_calls,
        stats.half_calls,
        stats.split_site_half_calls,
        stats.malformed_rows,
        stats.filtered_rows,
        stats.low_quality_calls,
        parser.ucsc_mito_dropped,
        tuple(parser.split_site_withheld),
        tuple(sorted(parser.split_pair_records)),
    )
    return calls, counts


@pytest.mark.parametrize("fixture", VCFS, ids=lambda p: str(p.relative_to(FIXTURES)))
@pytest.mark.parametrize("drop_mito", [False, True])
def test_every_vcf_fixture_reads_the_same(tmp_path, monkeypatch, fixture, drop_mito):
    samples = _samples(fixture)
    sample = samples[0] if len(samples) > 1 else None
    try:
        before = _read(fixture, False, monkeypatch, sample=sample, drop_mito=drop_mito)
    except Exception as exc:  # a fixture the parser refuses must refuse the same way
        with pytest.raises(type(exc)):
            _read(fixture, True, monkeypatch, sample=sample, drop_mito=drop_mito)
        return
    assert _read(fixture, True, monkeypatch, sample=sample, drop_mito=drop_mito) == before


def _indel(pos, gt="0/1"):
    """An insertion (G>GA): at the same position, but never a split-site component."""
    return record("GA", gt, pos=pos)


SHAPES = {
    "exact": [record("A", "1/0"), record("C", "0/1")],
    "vt": [record("A", "1/."), record("C", "./1")],
    "phased_vt": [record("A", "1|.", ps="7"), record("C", ".|1", ps="7")],
    "agreeing": [record("A", "0/1"), record("C", "0/1")],
    "rejected_sibling": [record("A", "0/1"), record("C", "0/1", filt="q10")],
    "low_gq_sibling": [record("A", "0/1"), record("C", "0/1", gq=3)],
    "missing_sibling": [record("A", "0/1"), record("C", "./.")],
    "contradicting": [record("A", "1/1"), record("C", "0/1")],
    "third_alt": [record("A", "1/0"), record("C", "0/1"), record("T", "0/1")],
    "haploid": [record("A", "1"), record("C", "0/1")],
    "two_rsids": [record("A", "1/0", rsid=RSID), record("C", "0/1", rsid=OTHER)],
    "same_alt_repeat": [record("A", "0/1"), record("A", "1/1")],
    "lone": [record("A", "0/1")],
    "lone_rejected": [record("A", "0/1", filt="q10")],
    "half_alone": [record("A", "1/.")],
    "reference": [record("A", "0/0"), record("C", "0/1")],
}


def _layout(name, build):
    """The shape, surrounded by neighbours, an indel at the same position and a later site."""
    rows = SHAPES[name]
    here = POS[build]
    before = [record("A", "0/1", rsid=".", pos=here - 50)]
    middle = [_indel(here)]  # same position, not a component: it keeps its own call
    after = [record("A", "1/0", pos=here + 50), record("C", "0/1", pos=here + 50)]
    return before + rows[:1] + middle + rows[1:] + after


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("name", SHAPES)
@pytest.mark.parametrize("layout", ["alone", "surrounded"])
def test_split_shapes_read_the_same(tmp_path, monkeypatch, build, name, layout):
    rows = SHAPES[name] if layout == "alone" else _layout(name, build)
    source = write_vcf(tmp_path / "input.vcf", build, rows)
    assert _read(source, True, monkeypatch) == _read(source, False, monkeypatch)


@pytest.mark.parametrize("build", BUILDS)
def test_malformed_and_unsorted_rows_end_a_site_the_same(tmp_path, monkeypatch, build):
    source = write_vcf(
        tmp_path / "input.vcf",
        build,
        [record("A", "1/0"), record("C", "0/1"), record("A", "0/1", pos=POS[build] + 9)],
    )
    lines = source.read_text().splitlines(keepends=True)
    # A malformed row between two records of one position, then the same
    # position again after another position (an unsorted file).
    broken = [*lines[:-2], "1\tbroken\n", lines[-2], lines[-1], lines[-3]]
    source.write_text("".join(broken))
    assert _read(source, True, monkeypatch) == _read(source, False, monkeypatch)

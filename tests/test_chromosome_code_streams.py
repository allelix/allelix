"""Reference watch, tabix and raw GenomeStudio envelopes share locus semantics."""

import sqlite3
from contextlib import closing

import pytest
from test_annotation_parity_and_named_disclosures import _vcf
from validated_call_helpers import RSID, analyze, reference_effect, signatures

from allelix.annotators.base import variant_locus
from allelix.cli.utility import _fetch_rows
from allelix.parsers import detect_parser
from allelix.parsers.vcf import _parse_data_line

ALIASES = [("23", "X"), ("24", "Y"), ("25", "X"), ("26", "MT"), ("XY", "X")]


@pytest.mark.parametrize(("code", "canonical"), ALIASES)
@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("gt", ["0/0", "0"])
def test_reference_watch_coded_contig_matches_canonical_target(
    tmp_path, code, canonical, build, gt
):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    for filename, table, column in [
        (f"clinvar.{build}.sqlite", "clinvar_variants", "chromosome"),
        ("rsid_identity.sqlite", "rsid_identity", "chrom"),
        ("gnomad.sqlite", "gnomad_frequencies", "chrom"),
    ]:
        with closing(sqlite3.connect(data / filename)) as conn:
            conn.execute(f"UPDATE {table} SET {column}=?", (canonical,))
            conn.commit()
    answers = []
    for name, chrom in [("control", canonical), ("coded", code)]:
        path = _vcf(tmp_path / f"{name}.vcf", [(RSID, chrom, pos, "G", "A", gt, ".")])
        answers.append(analyze(data, path, build, ("gwas", "pharmgkb", "snpedia")))
    assert {a.source for a in answers[0].annotations} == {"gwas", "pharmgkb", "snpedia"}
    assert signatures(answers[1]) == signatures(answers[0])
    assert all(a.haploid == (gt == "0") for a in answers[1].annotations)


@pytest.mark.parametrize(("code", "canonical"), ALIASES)
def test_indexed_extract_finds_raw_numeric_contig_by_canonical_coordinate(
    tmp_path, code, canonical
):
    import pysam

    raw = _vcf(tmp_path / "source.vcf", [(RSID, code, 65001, "G", "A", "0/1", ".")])
    compressed = tmp_path / "source.vcf.gz"
    pysam.tabix_compress(str(raw), str(compressed))
    pysam.tabix_index(str(compressed), preset="vcf")
    with pysam.TabixFile(str(compressed)) as tabix:
        rows = _fetch_rows(tabix, canonical, 65001, 0, _parse_data_line)
    assert len(rows) == 1
    assert rows[0].rsid == RSID
    assert variant_locus(rows[0]) == (canonical, 65001)


@pytest.mark.parametrize(("code", "canonical"), ALIASES)
def test_genomestudio_envelope_reaches_canonical_locus(tmp_path, code, canonical):
    path = tmp_path / "synthetic.txt"
    path.write_text(
        "# MyHappyGenes\n"
        "SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward\n"
        f"{RSID}\t{code}\t65001\tG\tA\n"
    )
    parser = detect_parser(path)
    assert parser.name == "myhappygenes"
    rows = list(parser.parse(path))
    assert len(rows) == 1
    assert variant_locus(rows[0]) == (canonical, 65001)


def test_genomestudio_zero_remains_unknown(tmp_path):
    path = tmp_path / "synthetic.txt"
    path.write_text(
        "# MyHappyGenes\n"
        "SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward\n"
        f"{RSID}\t0\t65001\tG\tA\n"
    )
    rows = list(detect_parser(path).parse(path))
    assert len(rows) == 1
    assert rows[0].chromosome == "0"
    assert variant_locus(rows[0]) == ("0", 65001)

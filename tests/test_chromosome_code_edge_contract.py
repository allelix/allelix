"""Haploid guards, chromosome summaries and explicit PAR survive aliases."""

from collections import Counter

import pytest
from click.testing import CliRunner
from test_annotation_parity_and_named_disclosures import _gnomad_site
from test_chromosome_code_commands import CONFIG, POS, RSID, _data

from allelix.cli import main
from allelix.compare import compare_variants
from allelix.parsers import detect_parser
from allelix.parsers._helpers import normalize_chromosome
from allelix.parsers.ftdna_famfinder import FTDNAFamFinderParser

ALIASES = [("23", "X"), ("24", "Y"), ("25", "X"), ("26", "MT"), ("XY", "X")]


def _famfinder(path, chrom, second, allele="G"):
    path.write_text(
        "# Family Tree DNA - FamFinder\n"
        "RSID\tCHROMOSOME\tPOSITION\tALLELE1\tALLELE2\n"
        f"{RSID}\t{chrom}\t{POS}\t{allele}{second}\n"
    )
    parser = FTDNAFamFinderParser()
    assert parser.can_parse(path)
    return list(parser.parse(path))


@pytest.mark.parametrize(("code", "canonical"), [("24", "Y"), ("26", "MT")])
@pytest.mark.parametrize("prefix", ["", "chr", "CHR"])
@pytest.mark.parametrize("second", ["", "\t"], ids=["absent-allele2", "empty-allele2"])
@pytest.mark.parametrize("allele", ["G", "-"], ids=["called", "no-call"])
def test_famfinder_haploid_guard_maps_codes(tmp_path, code, canonical, prefix, second, allele):
    control = _famfinder(tmp_path / "control.txt", canonical, second, allele)
    actual = _famfinder(tmp_path / "coded.txt", prefix + code, second, allele)
    assert len(control) == 1
    assert control[0].haploid is (allele == "G")
    assert actual == control
    assert len(actual) == 1
    assert actual[0].haploid == control[0].haploid
    assert actual[0].is_no_call == control[0].is_no_call


@pytest.mark.parametrize("chrom", ["1", "X", "23", "25", "XY", "0", "GL00023.1"])
@pytest.mark.parametrize("second", ["", "\t"])
def test_famfinder_truncated_non_y_mt_is_still_rejected(tmp_path, chrom, second):
    assert _famfinder(tmp_path / "truncated.txt", chrom, second) == []


def _chip(path, chroms, reader="myhappygenes"):
    headers = {
        "myhappygenes": "# MyHappyGenes\n"
        "SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward\n",
        "ancestrydna": "#AncestryDNA raw data download\n"
        "rsid\tchromosome\tposition\tallele1\tallele2\n",
        "famfinder": "# Family Tree DNA - FamFinder\n"
        "RSID\tCHROMOSOME\tPOSITION\tALLELE1\tALLELE2\n",
        "decodeme": "Name,Variation,Chromosome,Position,Strand,YourCode\n",
    }
    text = headers[reader]
    for i, chrom in enumerate(chroms):
        rsid = RSID if i == 0 else "rs9000000092"
        text += (
            f"{rsid},G/A,{chrom},{POS + i},+,GA\n"
            if reader == "decodeme"
            else f"{rsid}\t{chrom}\t{POS + i}\tG\tA\n"
        )
    path.write_text(text)
    return path


@pytest.mark.parametrize(("code", "canonical"), ALIASES)
@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
def test_compare_breakdown_coalesces_aliases_in_either_order(tmp_path, code, canonical, build):
    paths = [
        _chip(tmp_path / "mixed.txt", [code, canonical]),
        _chip(tmp_path / "letters.txt", [canonical, canonical]),
    ]
    calls = [list(detect_parser(p).parse(p)) for p in paths]
    assert [len(c) for c in calls] == [2, 2]
    control = compare_variants(calls[1], calls[1], build1=build, build2=build)
    expected = {canonical: Counter(concordant=2)}
    assert control.chromosome_counts == expected
    for first, second in (calls, calls[::-1]):
        result = compare_variants(first, second, build1=build, build2=build)
        assert result.shared == result.concordant == 2
        assert result.chromosome_counts == expected


@pytest.mark.parametrize(("code", "canonical"), [*ALIASES, ("0", "0"), ("GL00023.1", "GL00023.1")])
@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
def test_stats_breakdown_coalesces_aliases_but_keeps_unknowns(tmp_path, code, canonical, build):
    data = tmp_path / "data"
    data.mkdir()
    (data / "config.toml").write_text(CONFIG.replace("= true", "= false"))
    for label, chroms in [("letters", [canonical, canonical]), ("mixed", [code, canonical])]:
        source = _chip(tmp_path / f"{label}.txt", chroms)
        result = CliRunner().invoke(
            main, ["stats", str(source), "--build", build.lower(), "--data-dir", str(data)]
        )
        assert result.exit_code == 0, (result.output, result.exception)
        rows = [" ".join(line.replace("│", " ").split()) for line in result.output.splitlines()]
        assert f"{canonical} 2" in rows, result.output
        if code != canonical:
            assert f"{code} 1" not in rows, result.output


@pytest.mark.parametrize("reader", ["myhappygenes", "ancestrydna", "famfinder", "decodeme"])
@pytest.mark.parametrize("code", ["25", "XY"])
@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
def test_plink_explicit_par_is_kept_from_chip_readers(tmp_path, reader, code, build):
    data = _data(tmp_path, "X")
    _gnomad_site(data, RSID, "X", POS, "G", "A", 0.1)
    (data / "config.toml").write_text(CONFIG.replace("gnomad = false", "gnomad = true"))
    source = _chip(tmp_path / "source.txt", [code], reader)
    prefix = tmp_path / "out"
    result = CliRunner().invoke(
        main,
        [
            "export",
            "plink",
            str(source),
            "--build",
            build.lower(),
            "--data-dir",
            str(data),
            "--output-prefix",
            str(prefix),
        ],
    )
    assert result.exit_code == 0, (result.output, result.exception)
    rows = [line.split() for line in prefix.with_suffix(".bim").read_text().splitlines()]
    assert len(rows) == 1
    assert rows[0][0] == "25", rows
    assert rows[0][1] == RSID


@pytest.mark.parametrize("code", ["xy", "xY", "Xy", "chrxy", "CHRxY", "cHrXY"])
def test_par_alias_is_case_insensitive(code):
    assert normalize_chromosome(code) == "X"

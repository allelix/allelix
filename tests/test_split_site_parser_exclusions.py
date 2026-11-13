"""Split-site buffering preserves input exclusions and malformed-row handling."""

import json

import pytest
from click.testing import CliRunner
from test_split_person_genotype import _data

from allelix.cli import main
from allelix.parsers.vcf import VcfParser

_HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"


def _rows(chrom, pair):
    rows = [f"{chrom}\t100\trs1801133\tG\tA\t.\tPASS\t.\tGT\t1/0"]
    if pair:
        rows.append(f"{chrom}\t100\trs1801133\tG\tT\t.\tPASS\t.\tGT\t0/1")
    return rows


@pytest.mark.parametrize("chrom", ("M", "chrM"))
@pytest.mark.parametrize("pair", (False, True))
@pytest.mark.parametrize("rcrs", (False, True))
def test_split_buffer_keeps_the_parser_mito_exclusion(tmp_path, chrom, pair, rcrs):
    source = tmp_path / "calls.vcf"
    source.write_text(_HEADER + "\n".join(_rows(chrom, pair)) + "\n")
    parsed = []
    for combine in (False, True):
        parser = VcfParser()
        parser.combine_split_sites = combine
        parser.drop_ucsc_mito = True
        parser.rcrs_mito_names = frozenset({chrom}) if rcrs else frozenset()
        calls = list(parser.parse(source))
        parsed.append((calls, parser.ucsc_mito_dropped, parser.split_site_calls))
    assert bool(parsed[0][0]) is rcrs, "Establish the existing exclusion decision"
    calls, dropped, split_count = parsed[1]
    assert bool(calls) is rcrs, "Grouping must not re-admit an excluded mitochondrial call"
    assert dropped == (0 if rcrs else len(_rows(chrom, pair)))
    assert split_count == 0, "Build exclusion is not a split-genotype abstention"
    if rcrs:
        assert len(calls) == 1
        assert set(calls[0].carried_pair) == ({"A", "T"} if pair else {"G", "A"})


@pytest.mark.parametrize("chrom", ("MT", "chrMT"))
@pytest.mark.parametrize("pair", (False, True))
def test_grch37_mito_names_still_give_calls(tmp_path, chrom, pair):
    source = tmp_path / "calls.vcf"
    source.write_text(_HEADER + "\n".join(_rows(chrom, pair)) + "\n")
    parser = VcfParser()
    parser.combine_split_sites = True
    parser.drop_ucsc_mito = True
    calls = list(parser.parse(source))
    assert len(calls) == 1
    assert parser.ucsc_mito_dropped == 0
    assert parser.split_site_calls == 0


@pytest.mark.parametrize("chrom", ("M", "chrM"))
@pytest.mark.parametrize("pair", (False, True))
@pytest.mark.parametrize("explicit", (False, True))
def test_analysis_does_not_report_excluded_mito_calls(tmp_path, chrom, pair, explicit):
    data = tmp_path / "data"
    _data(data, "GRCh37")
    (data / "config.toml").write_text(
        "[sources]\nclinvar=false\nsnpedia=false\npharmgkb=true\ngwas=false\n"
        "gnomad=false\nalphamissense=false\ncadd=false\nrsid_identity=false\n"
    )
    source = tmp_path / "calls.vcf"
    source.write_text(
        _HEADER.replace(
            "##fileformat=VCFv4.2\n", "##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n"
        )
        + "1\t200\trs999999991\tG\tA\t.\tPASS\t.\tGT\t0/1\n"
        + "\n".join(_rows(chrom, pair))
        + "\n"
    )
    report = tmp_path / "report.json"
    args = [
        "analyze",
        str(source),
        "--data-dir",
        str(data),
        "--no-update",
        "--output",
        str(report),
    ]
    if explicit:
        args += ["--build", "grch37"]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, (result.output, result.exception)
    payload = json.loads(report.read_text())
    assert not payload["annotations"], "The GRCh37 excluded call must not make a ClinPGx finding"
    assert "mitochondrial" in result.output.lower(), "Keep the existing exclusion notice"


@pytest.mark.parametrize("columns", (1, 2, 3, 4, 5, 7))
@pytest.mark.parametrize("before", (False, True))
def test_short_malformed_record_is_still_skipped_and_counted(tmp_path, columns, before):
    malformed = "\t".join(["1", "300", "rs999999992", "G", "A", ".", "PASS"][:columns])
    good = _rows("1", False)[0]
    records = [malformed, good] if before else [good, malformed]
    source = tmp_path / "calls.vcf"
    source.write_text(_HEADER + "\n".join(records) + "\n")
    outcomes = []
    for combine in (False, True):
        parser = VcfParser()
        parser.combine_split_sites = combine
        try:
            calls = list(parser.parse(source))
        except (IndexError, ValueError) as exc:
            pytest.fail(
                f"A row already rejected as malformed must not crash split buffering: {exc}"
            )
        outcomes.append(
            (
                [(v.rsid, v.position, v.carried_pair) for v in calls],
                parser.quality_stats.malformed_rows,
            )
        )
    assert outcomes[0][1] == 1
    assert len(outcomes[0][0]) == 1
    assert outcomes[1] == outcomes[0], "Grouping must preserve the ordinary malformed-row policy"

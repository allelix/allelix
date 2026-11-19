"""The split-site warning names a site split because its ALTs belong to two rsIDs (2.7.1, B-07).

Two records that are one exact complementary pair (G>A 1/0 + G>C 0/1) give
the person's genotype, AC. When their ALTs belong to two different rsIDs,
written on the records or found by position, the alleles are assessed apart,
as half calls, but the genotype is known. RC2's one warning said for every
counted site that "the person's genotype there cannot be told from them".
SUP's ruling: such a site gets its own sentence; the "cannot be told"
sentence stays for the sites that really are not one genotype.
"""

import html
import re

import pytest
from click.testing import CliRunner
from split_site_contract_helpers import BUILDS, OTHER, RSID, data_dir, record, run, write_vcf
from test_split_pair_two_position_rsids import _pair, _two_rsid_data

from allelix.cli import main

UNKNOWABLE = "cannot be told from them"
KNOWN = re.compile(r"genotype there is known", re.IGNORECASE)
SEPARATE = re.compile(r"belong to two different rsIDs", re.IGNORECASE)


def _flat(text):
    return " ".join(html.unescape(re.sub(r"<[^>]+>", " ", text)).split())


def _split_message(report):
    (warning,) = [w for w in report["report_warnings"] if w["code"] == "split_site_calls"]
    return warning


def _surfaces(tmp_path, data, rows, build, name):
    source = write_vcf(tmp_path / f"{name}.vcf", build, rows)
    report, console = run(data, source, tmp_path / f"{name}.json", build)
    page = tmp_path / f"{name}.html"
    result = CliRunner(env={"COLUMNS": "240"}).invoke(
        main,
        [
            "analyze",
            str(source),
            "--data-dir",
            str(data),
            "--no-update",
            "--output",
            str(page),
            "--build",
            build.lower(),
            "--min-magnitude",
            "0",
        ],
    )
    assert result.exit_code == 0, (result.output, repr(result.exception))
    return report, _flat(console), _flat(page.read_text())


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("ids", ["position", "written"])
def test_a_pair_split_by_two_rsids_is_a_known_genotype(tmp_path, build, ids):
    data = _two_rsid_data(tmp_path / "data")
    rows = _pair((".", ".") if ids == "position" else (RSID, OTHER))
    report, console, page = _surfaces(tmp_path, data, rows, build, "pair")
    warning = _split_message(report)
    assert warning["count"] == 2
    for surface, text in (("json", warning["message"]), ("terminal", console), ("html", page)):
        assert SEPARATE.search(text) and KNOWN.search(text), (surface, text)
        assert UNKNOWABLE not in text, (surface, "The two-ALT genotype is known here")


@pytest.mark.parametrize("build", BUILDS)
def test_records_that_only_agree_keep_the_cannot_be_told_sentence(tmp_path, build):
    data = data_dir(tmp_path / "data", identity=True, alts="AC")
    rows = [record("A", "0/1"), record("C", "0/1")]
    report, console, page = _surfaces(tmp_path, data, rows, build, "agree")
    warning = _split_message(report)
    assert warning["count"] == 2
    for surface, text in (("json", warning["message"]), ("terminal", console), ("html", page)):
        assert UNKNOWABLE in text, (surface, text)
        assert not KNOWN.search(text), (surface, "Two 0/1 records are not one genotype")


@pytest.mark.parametrize("build", BUILDS)
def test_both_kinds_in_one_file_give_both_sentences_with_their_counts(tmp_path, build):
    data = _two_rsid_data(tmp_path / "data")
    pair = _pair((RSID, OTHER))
    probe = write_vcf(tmp_path / "probe.vcf", build, pair)
    elsewhere = int(probe.read_text().splitlines()[-1].split("\t")[1]) + 10
    agree = [record(alt, "0/1", rsid=".", pos=elsewhere) for alt in "AC"]
    report, _console, _page = _surfaces(tmp_path, data, pair + agree, build, "both")
    message = _split_message(report)["message"]
    assert _split_message(report)["count"] == 4
    assert SEPARATE.search(message) and KNOWN.search(message) and UNKNOWABLE in message
    assert re.search(r"\b2\b[^.]*belong to two different rsIDs", message), message

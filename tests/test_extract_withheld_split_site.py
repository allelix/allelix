"""extract says a contradicting split site is present but withheld (2.7.1).

Reviewer A's A-04a: analyze withholds a site whose records contradict each
other and names it; extract printed "not in file" for the rsID those records
name. It now shows the site as a no-call, with a note that its records
contradict each other, on the sequential and the indexed path. It is never
turned into a genotype.
"""

import pytest
from click.testing import CliRunner
from test_extract_duplicate_placements import _cache
from test_extract_split_site_contract import POS, RS, _input
from test_homref_ref_by_locus import _write_identity

from allelix.cli import main, utility


def _run(tmp_path, monkeypatch, build, gts, *, indexed):
    data = _cache(tmp_path, POS)
    _write_identity(data, [(RS, "1", POS, "G", a, POS, "G", a, 0) for a in "AT"])
    pairs = [(a, gt, "42", RS, "PASS") for a, gt in zip("AT", gts, strict=True)]
    source = _input(tmp_path, build, pairs, indexed=indexed)
    captured = []
    monkeypatch.setattr(utility, "_render_extract_table", lambda _p, _w, f: captured.append(f))
    args = ["extract", str(source), "--data-dir", str(data), "--snps", RS, "--build", build]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, (result.output, result.exception)
    (found,) = captured
    return found, " ".join(result.output.split())


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("indexed", [False, True])
@pytest.mark.parametrize("gts", [("1/1", "0/1"), ("1", "0/1"), ("0|1", "0|1")])
def test_contradicting_site_is_a_present_no_call_with_a_note(
    tmp_path, monkeypatch, build, indexed, gts
):
    found, text = _run(tmp_path, monkeypatch, build, gts, indexed=indexed)
    assert RS in found, "The records are in the file: not 'not in file'"
    assert found[RS].is_no_call, "A contradicting site is never shown as a genotype"
    assert (found[RS].chromosome, found[RS].position) == ("1", POS)
    assert RS in text and "contradict" in text


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("indexed", [False, True])
def test_exact_pair_has_no_withheld_note(tmp_path, monkeypatch, build, indexed):
    found, text = _run(tmp_path, monkeypatch, build, ("1/0", "0/1"), indexed=indexed)
    assert sorted((found[RS].allele1, found[RS].allele2)) == ["A", "T"]
    assert "contradict" not in text

"""Extract must not print the first component as a person's split-site genotype."""

import pysam
import pytest
from click.testing import CliRunner
from test_extract_duplicate_placements import _cache
from test_homref_ref_by_locus import _vcf, _write_identity

from allelix.cli import main, utility

RS = "rs900252201"
POS = 200400


def _input(tmp_path, build, pairs, *, indexed=False, reverse=False):
    rows = [
        ("1", POS, rs, "G", alt, 99, filt, ".", "GT:GQ:PS", f"{gt}:99:{ps}")
        for alt, gt, ps, rs, filt in pairs
    ]
    if reverse:
        rows.reverse()
    source = _vcf(tmp_path / "input.vcf", rows, build=build)
    if indexed:
        gz = source.with_suffix(".vcf.gz")
        pysam.tabix_compress(str(source), str(gz), force=True)
        pysam.tabix_index(str(gz), preset="vcf", force=True)
        return gz
    return source


def _extract(tmp_path, monkeypatch, build, pairs, *, indexed, auto, reverse=False):
    data = _cache(tmp_path, POS)
    _write_identity(data, [(RS, "1", POS, "G", a, POS, "G", a, 0) for a in "AT"])
    source = _input(tmp_path, build, pairs, indexed=indexed, reverse=reverse)
    captures = []
    fast_results = []
    original_fast = utility._try_tabix_extract

    def fast(*args, **kwargs):
        answer = original_fast(*args, **kwargs)
        fast_results.append(answer)
        return answer

    monkeypatch.setattr(utility, "_try_tabix_extract", fast)
    monkeypatch.setattr(utility, "_render_extract_table", lambda _p, _w, f: captures.append(f))
    args = [
        "extract",
        str(source),
        "--data-dir",
        str(data),
        "--snps",
        RS,
        "--build",
        "auto" if auto else build,
    ]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, (result.output, result.exception)
    assert len(captures) == 1
    if indexed and build == "GRCh38":
        assert fast_results and fast_results[0] is not None, "Exercise the real indexed path"
    return captures[0]


def _pair(gts, ids=(RS, RS), phase=("42", "42")):
    return [(a, gt, ps, rs, "PASS") for a, gt, ps, rs in zip("AT", gts, phase, ids, strict=True)]


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("indexed", [False, True])
@pytest.mark.parametrize("auto", [False, True])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("gts", [("1/0", "0/1"), ("1/.", "./1"), ("1|0", "0|1"), ("1|.", ".|1")])
def test_exact_split_pair_extracts_the_full_person_call(
    tmp_path, monkeypatch, build, indexed, auto, reverse, gts
):
    found = _extract(
        tmp_path, monkeypatch, build, _pair(gts), indexed=indexed, auto=auto, reverse=reverse
    )
    assert RS in found
    call = found[RS]
    assert sorted((call.allele1, call.allele2)) == ["A", "T"], call
    assert not call.is_no_call and not call.is_half_call


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("indexed", [False, True])
@pytest.mark.parametrize("ids", [(RS, "."), (".", RS), (RS, "22:200400:G:T")])
def test_compatible_id_labels_do_not_lose_the_second_allele(
    tmp_path, monkeypatch, build, indexed, ids
):
    found = _extract(
        tmp_path, monkeypatch, build, _pair(("1/0", "0/1"), ids), indexed=indexed, auto=True
    )
    assert RS in found
    assert sorted((found[RS].allele1, found[RS].allele2)) == ["A", "T"]


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("indexed", [False, True])
@pytest.mark.parametrize("gts", [("1/1", "0/1"), ("1", "0/1"), ("0|1", "0|1")])
def test_contradicting_records_cannot_print_a_component_as_the_person(
    tmp_path, monkeypatch, build, indexed, gts
):
    found = _extract(tmp_path, monkeypatch, build, _pair(gts), indexed=indexed, auto=False)
    assert RS not in found or found[RS].is_no_call, found


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("indexed", [False, True])
def test_uncertain_agreeing_site_does_not_claim_the_first_full_genotype(
    tmp_path, monkeypatch, build, indexed
):
    found = _extract(
        tmp_path, monkeypatch, build, _pair(("0/1", "0/1")), indexed=indexed, auto=True
    )
    # The exact UI representation is not prescribed. A fully called AG or GT
    # is false: analyze supplies passing alleles only, not those genotypes.
    assert RS not in found or found[RS].is_no_call or found[RS].is_half_call, found


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("indexed", [False, True])
def test_merged_record_control_is_unchanged(tmp_path, monkeypatch, build, indexed):
    found = _extract(
        tmp_path, monkeypatch, build, [("A,T", "1/2", ".", RS, "PASS")], indexed=indexed, auto=True
    )
    assert sorted((found[RS].allele1, found[RS].allele2)) == ["A", "T"]


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("indexed", [False, True])
def test_distinct_rsids_do_not_make_a_full_pair(tmp_path, monkeypatch, build, indexed):
    found = _extract(
        tmp_path,
        monkeypatch,
        build,
        _pair(("1/0", "0/1"), (RS, "rs900252202")),
        indexed=indexed,
        auto=True,
    )
    assert RS not in found or found[RS].is_no_call or found[RS].is_half_call, found

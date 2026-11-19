"""A row with a malformed position still ends a pending split site (2.7.1).

An out-of-range POS is skipped as a malformed row, as a non-numeric POS is.
Either way, as in 2.7.0, the row is a boundary in the file's record order:
records of one position on either side of it are not one site, so only the
row's own call is dropped.
"""

import pytest

from allelix.parsers.vcf import VcfParser

HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
INVALID = ["0", "-1", "-5", str(2**63), str(2**63 + 1), str(10**80), "not-a-position"]


def _parse(tmp_path, rows):
    path = tmp_path / "input.vcf"
    path.write_text(HEADER + "".join("\t".join(r) + "\n" for r in rows))
    parser = VcfParser()
    parser.combine_split_sites = True
    return parser, list(parser.parse(path))


def _row(pos, alt, gt, rsid="rs999000555", chrom="1"):
    return (chrom, pos, rsid, "G", alt, "50", "PASS", ".", "GT", gt)


@pytest.mark.parametrize("invalid", INVALID)
@pytest.mark.parametrize("bad_shape", ["component", "reference", "indel"])
def test_malformed_position_row_ends_the_pending_site(tmp_path, invalid, bad_shape):
    bad = {
        "component": _row(invalid, "C", "0/1", "rs999000556"),
        "reference": _row(invalid, "C", "0/0", "rs999000556"),
        "indel": ("1", invalid, "rs999000556", "G", "GA", "50", "PASS", ".", "GT", "0/1"),
    }[bad_shape]
    parser, calls = _parse(tmp_path, [_row("200500", "A", "1/0"), bad, _row("200500", "T", "0/1")])
    assert parser.quality_stats.malformed_rows == 1
    assert [sorted((v.allele1, v.allele2)) for v in calls] == [["A", "G"], ["G", "T"]], (
        "Records on either side of a malformed row must not become one call"
    )
    assert not any(v.allele_records for v in calls)
    assert parser.split_site_calls == 0


@pytest.mark.parametrize("invalid", INVALID)
def test_a_site_after_the_malformed_row_still_pairs(tmp_path, invalid):
    parser, calls = _parse(
        tmp_path,
        [
            _row("200400", "A", "1/0"),
            _row(invalid, "C", "0/1", "rs999000556"),
            _row("200500", "A", "1/0"),
            _row("200500", "T", "0/1"),
        ],
    )
    assert parser.quality_stats.malformed_rows == 1
    assert [sorted((v.allele1, v.allele2)) for v in calls] == [["A", "G"], ["A", "T"]]
    assert calls[1].allele_records is not None


def test_without_a_malformed_row_the_pair_is_one_call(tmp_path):
    parser, calls = _parse(tmp_path, [_row("200500", "A", "1/0"), _row("200500", "T", "0/1")])
    assert [sorted((v.allele1, v.allele2)) for v in calls] == [["A", "T"]]
    assert parser.quality_stats.malformed_rows == 0

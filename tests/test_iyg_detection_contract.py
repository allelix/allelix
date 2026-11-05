"""Strict two-column rsID/genotype envelopes, with invented rows only."""

import gzip

import pytest

from allelix.parsers import ParserNotFoundError, detect_parser

ROWS = [
    "rs3999980001\tAG",
    "rs3999980002\tCC",
    "rs3999980003\tT",
    "rs3999980004\t--",
    "rs3999980005\tDI",
]


def _file(tmp_path, rows, ending="\n", compressed=False):
    body = (ending.join(["# Synthetic two-column export", *rows]) + ending).encode()
    path = tmp_path / ("synthetic.txt.gz" if compressed else "synthetic.txt")
    path.write_bytes(gzip.compress(body) if compressed else body)
    return path


@pytest.mark.parametrize("ending", ("\n", "\r\n", "\r"), ids=("lf", "crlf", "cr"))
@pytest.mark.parametrize("compressed", (False, True), ids=("plain", "gzip"))
@pytest.mark.parametrize("quoted", (False, True))
def test_strict_two_column_input_is_detected_and_preserves_calls(
    tmp_path, ending, compressed, quoted
):
    rows = ['"' + row.replace("\t", '"\t"') + '"' for row in ROWS] if quoted else ROWS
    path = _file(tmp_path, rows, ending, compressed)
    parser = detect_parser(path)
    assert parser.name == "rsid-genotype"
    values = list(parser.parse(path))
    assert [(v.rsid, v.allele1, v.allele2, v.haploid) for v in values] == [
        ("rs3999980001", "A", "G", False),
        ("rs3999980002", "C", "C", False),
        ("rs3999980003", "T", "T", True),
        ("rs3999980004", "-", "-", False),
        ("rs3999980005", "D", "I", False),
    ]
    assert all(not v.chromosome and not v.position for v in values)
    assert parser.get_metadata(path)["build"] is None


@pytest.mark.parametrize("rs_count", (3, 4, 5))
def test_four_in_five_identifiers_must_be_rsids(tmp_path, rs_count):
    rows = ROWS[:rs_count] + [f"i90000000{i}\tAG" for i in range(5 - rs_count)]
    path = _file(tmp_path, rows)
    if rs_count >= 4:
        parser = detect_parser(path)
        assert parser.name == "rsid-genotype"
        assert len(list(parser.parse(path))) == 5
    else:
        with pytest.raises(ParserNotFoundError):
            detect_parser(path)


@pytest.mark.parametrize(
    "stray",
    (
        "rs3999980090 AG",
        "rs3999980090,AG",
        "rs3999980090\t1\tAG",
        "rs3999980090\t1\t100\tAG",
        "rs3999980090\tA G",
        "rs3999980090\tAGT",
        "rs3999980090\tNN",
        "rsid\tgenotype",
    ),
)
@pytest.mark.parametrize("at_start", (False, True))
def test_other_layouts_are_not_claimed_by_the_two_column_reader(tmp_path, stray, at_start):
    rows = [stray, *ROWS] if at_start else [*ROWS, stray]
    with pytest.raises(ParserNotFoundError):
        detect_parser(_file(tmp_path, rows))


@pytest.mark.parametrize("ending", ("\n", "\r\n", "\r"), ids=("lf", "crlf", "cr"))
@pytest.mark.parametrize("layout", ("chromosome-columns", "allele-columns"))
def test_unknown_source_column_layouts_remain_rejected(tmp_path, ending, layout):
    if layout == "chromosome-columns":
        rows = ["SNP,Chr,Allele1,Allele2"] + [f"rs399998000{i},1,A,G" for i in range(1, 6)]
    else:
        rows = ["SNP_ID\tAllele1\tAllele2"] + [f"rs399998000{i}\tA\tG" for i in range(1, 6)]
    with pytest.raises(ParserNotFoundError):
        detect_parser(_file(tmp_path, rows, ending))

"""An array row at position 0 is an unplaced rsID call, as in 2.7.0; a VCF POS stays 1 or more."""

import logging

import pytest
from test_position_integer_bounds import READERS, _source

from allelix.parsers import detect_parser

RS = "rs900618101"
ARRAYS = tuple(r for r in READERS if r != "vcf")


@pytest.mark.parametrize("reader", ARRAYS)
@pytest.mark.parametrize("chromosome", ["1", "0"])
def test_array_position_zero_is_an_unplaced_rsid_call(tmp_path, caplog, reader, chromosome):
    source = _source(tmp_path / "input.txt", reader, [("rs900618102", 12345), (RS, 0)])
    if chromosome == "0":
        text = source.read_text()
        for sep in ("\t", ","):
            text = text.replace(f"{RS}{sep}1{sep}0{sep}", f"{RS}{sep}0{sep}0{sep}")
            text = text.replace(f'"{RS}"{sep}"1"{sep}"0"', f'"{RS}"{sep}"0"{sep}"0"')
        text = text.replace(f"{RS},A/G,1,0,", f"{RS},A/G,0,0,")  # deCODEme
        source.write_text(text)
    parser = detect_parser(source)
    assert parser.name == reader
    with caplog.at_level(logging.WARNING):
        calls = list(parser.parse(source))
    assert [(v.rsid, v.position) for v in calls] == [("rs900618102", 12345), (RS, 0)], (
        "2.7.0 read an array rsID at position 0 as an unplaced call; the bounds "
        "check must not drop it as malformed"
    )
    assert calls[1].chromosome == chromosome
    assert calls[1].allele1 and calls[1].allele1 != "-"
    assert not [r for r in caplog.records if "position" in r.message.lower()], (
        "An unplaced array row is not malformed and gets no warning"
    )


@pytest.mark.parametrize("reader", ARRAYS)
@pytest.mark.parametrize("position", ["-1", "-0005", str(2**63)])
def test_array_negative_and_out_of_range_positions_stay_malformed(
    tmp_path, caplog, reader, position
):
    source = _source(tmp_path / "input.txt", reader, [("rs900618102", 12345), (RS, position)])
    with caplog.at_level(logging.WARNING):
        calls = list(detect_parser(source).parse(source))
    assert [v.rsid for v in calls] == ["rs900618102"]


def test_vcf_pos_zero_stays_malformed(tmp_path):
    source = _source(tmp_path / "input.vcf", "vcf", [("rs900618102", 12345), (RS, 0)])
    parser = detect_parser(source)
    assert [v.rsid for v in parser.parse(source)] == ["rs900618102"]
    assert parser.quality_stats.malformed_rows == 1

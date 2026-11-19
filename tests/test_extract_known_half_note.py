"""A selected known allele must not borrow a withheld group's note at its locus."""

from click.testing import CliRunner
from test_extract_duplicate_placements import RSID, _cache, _input, _row

from allelix.cli import main, utility


def test_known_half_at_revisited_locus_does_not_borrow_withheld_note(tmp_path, monkeypatch):
    data = _cache(tmp_path, 100)
    _input(tmp_path, [_row(100, "0/1")])
    source = tmp_path / "input.vcf"
    header = source.read_text().split("#CHROM")[0]
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    source.write_text(
        header
        + _row(100, "1/1")
        + _row(100, "0/1").replace("\tG\tA\t", "\tG\tT\t")
        # A different position ends the first group, including for a 0/0 row.
        + _row(200, "0/0", rsid=".")
        + _row(100, "0/1").replace("\tPASS\t", "\tLowQual\t")
        + _row(100, "0/1").replace("\tG\tA\t", "\tG\tT\t")
    )
    captured = []
    monkeypatch.setattr(
        utility, "_render_extract_table", lambda _p, _w, found: captured.append(found)
    )
    result = CliRunner().invoke(
        main,
        ["extract", str(source), "--snps", RSID, "--data-dir", str(data), "--build", "grch38"],
    )
    assert result.exit_code == 0, (result.output, result.exception)
    chosen = captured[0][RSID]
    assert chosen.position == 100
    alleles = (chosen.allele1, chosen.allele2)
    assert alleles in (("T", "-"), ("-", "-")), alleles
    text = " ".join(result.output.split())
    # The baseline selects the full no-call; the new ranking selects T/-.
    # In either case the explanation must describe the call actually shown.
    if alleles == ("T", "-"):
        assert "contradict" not in text, (
            "The selected T/- has a known allele. A separate withheld group at "
            "this same locus must not lend it its no-genotype explanation."
        )
    else:
        assert "contradict" in text

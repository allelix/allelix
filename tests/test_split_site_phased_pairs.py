"""Phase-confirmed complementary split calls follow the owner's exact-pair rule."""

import pytest
from test_annotation_parity_and_named_disclosures import _vcf
from test_split_person_genotype import SOURCES, _calls, _data, _row
from validated_call_helpers import BUILDS, analyze


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("swap_haplotypes", (False, True))
def test_complementary_phased_calls_in_one_phase_set_combine(
    tmp_path, build, reverse, swap_haplotypes
):
    gts = ("0|1", "1|0") if swap_haplotypes else ("1|0", "0|1")
    result = _run(tmp_path, build, gts, ("300", "300"), reverse)
    assert _calls(result) == [("pharmgkb", "AT"), ("snpedia", "AT")]
    assert not any(w.code == "split_site_calls" for w in result.report_warnings)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize(
    "gts,phase_sets",
    (
        (("0|1", "0|1"), ("300", "300")),
        (("1|0", "1|0"), ("300", "300")),
        (("0|1", "1|0"), ("300", "301")),
        (("1|0", "0|1"), ("300", "301")),
        (("0|1", "1/0"), ("300", ".")),
        (("1/0", "0|1"), (".", "300")),
    ),
)
def test_unrelated_or_contradictory_phase_abstains_and_is_counted(
    tmp_path, build, reverse, gts, phase_sets
):
    result = _run(tmp_path, build, gts, phase_sets, reverse)
    assert not result.annotations
    assert any(w.code == "split_site_calls" and w.count for w in result.report_warnings)


def _run(tmp_path, build, gts, phase_sets, reverse):
    data = tmp_path / "data"
    pos = _data(data, build)
    rows = [_row(pos, "A", gts[0]), _row(pos, "T", gts[1])]
    path = _vcf(tmp_path / "input.vcf", rows)
    lines = path.read_text().splitlines()
    for n, ps in enumerate(phase_sets):
        cols = lines[2 + n].split("\t")
        cols[8] = "GT:PS"
        cols[9] += f":{ps}"
        lines[2 + n] = "\t".join(cols)
    if reverse:
        lines[2:] = reversed(lines[2:])
    path.write_text("\n".join(lines) + "\n")
    return analyze(data, path, build, SOURCES)

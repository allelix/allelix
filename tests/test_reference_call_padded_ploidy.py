"""Reference records use the same ploidy locus in padded and minimal spelling."""

import pytest
from test_annotation_parity_and_named_disclosures import _vcf
from validated_call_helpers import BUILDS, analyze, reference_effect, sql, warning_count


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("sources", ((), ("pharmgkb",)), ids=("no_sources", "pharmgkb"))
@pytest.mark.parametrize("shape", ("prefix", "suffix", "both"))
@pytest.mark.parametrize("gts", (("0", "0/0"), ("0/0", "1")))
def test_padded_reference_ploidy_uses_the_same_site(tmp_path, build, reverse, sources, shape, gts):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    padded = {
        "prefix": (pos - 1, "CG", "CA"),
        "suffix": (pos, "GC", "AC"),
        "both": (pos - 1, "CGC", "CAC"),
    }[shape]
    counts = []
    for label, (at, ref, alt) in (("minimal", (pos, "G", "A")), ("padded", padded)):
        rows = [(".", "1", at, ref, alt, gts[0], "."), (".", "1", pos, "G", "A", gts[1], ".")]
        if reverse:
            rows.reverse()
        result = analyze(data, _vcf(tmp_path / f"{label}.vcf", rows), build, sources)
        counts.append(warning_count(result, "contradictory_calls"))
    assert counts[0] == (2 if sources else 0)
    assert counts[1] == counts[0]


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("sources", ((), ("pharmgkb",)), ids=("no_sources", "pharmgkb"))
@pytest.mark.parametrize("ref,alt", (("G", "GA"), ("GC", "G")))
@pytest.mark.parametrize("gts", (("0", "1/1"), ("0/0", "1")))
def test_reference_and_alt_indel_records_keep_their_written_ploidy_site(
    tmp_path, build, reverse, sources, ref, alt, gts
):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    ref_column, alt_column = ("ref37", "alt37") if build == "GRCh37" else ("ref38", "alt38")
    sql(
        data,
        "rsid_identity.sqlite",
        f"UPDATE rsid_identity SET {ref_column}=?, {alt_column}=?",
        (ref, alt),
    )
    rows = [(".", "1", pos, ref, alt, gt, ".") for gt in gts]
    if reverse:
        rows.reverse()
    result = analyze(data, _vcf(tmp_path / "input.vcf", rows), build, sources)
    assert warning_count(result, "contradictory_calls") == (2 if sources else 0)

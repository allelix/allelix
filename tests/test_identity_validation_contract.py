"""Known contradictory placement is rejected before any source matches the call."""

import pytest
from test_haploid_and_homozygous_warning import _array
from validated_call_helpers import BUILDS, RSID, SOURCES, analyze, call, setup, warning_count

from allelix.parsers.twentythreeandme import TwentyThreeAndMeParser


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("source", SOURCES)
@pytest.mark.parametrize("conflict", ["position", "chromosome"])
@pytest.mark.parametrize("format_name", ["vcf", "chip"])
def test_identity_contradiction_rejected_and_counted(
    tmp_path, build, source, conflict, format_name
):
    data = tmp_path / "data"
    pos = setup(data, build)
    control = analyze(data, call(tmp_path / "control.vcf", pos), build, (source,))
    assert control.annotations
    chrom = "2" if conflict == "chromosome" else "1"
    pos += 100 if conflict == "position" else 0
    if format_name == "chip":
        path = _array(tmp_path / "input.txt", [(RSID, chrom, pos, "AG")])
        parser = TwentyThreeAndMeParser()
    else:
        path = call(tmp_path / "input.vcf", pos, chrom=chrom)
        parser = None
    result = analyze(data, path, build, (source,), parser=parser)
    assert not result.annotations
    assert warning_count(result, "identity_locus_conflict") == 1


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("state", ["absent", "no_rsid"])
def test_unavailable_identity_is_not_a_contradiction(tmp_path, build, state):
    data = tmp_path / "data"
    pos = setup(data, build)
    if state == "absent":
        (data / "rsid_identity.sqlite").unlink()
    else:
        from validated_call_helpers import sql

        sql(data, "rsid_identity.sqlite", "DELETE FROM rsid_identity")
    result = analyze(data, call(tmp_path / "input.vcf", pos), build, ("pharmgkb",))
    assert len(result.annotations) == 1
    assert warning_count(result, "identity_locus_conflict") == 0


def test_all_sources_share_one_record_conflict_count(tmp_path):
    data = tmp_path / "data"
    pos = setup(data, "GRCh38")
    result = analyze(data, call(tmp_path / "input.vcf", pos + 100), "GRCh38")
    assert not result.annotations
    assert warning_count(result, "identity_locus_conflict") == 1


@pytest.mark.parametrize("reverse", [False, True])
def test_any_identity_placement_is_valid(tmp_path, reverse):
    from test_annotation_parity_and_named_disclosures import _identity

    data = tmp_path / "data"
    pos = setup(data, "GRCh38")
    (data / "rsid_identity.sqlite").unlink()
    rows = [(RSID, "1", 11856378, "G", "A", p, "G", "A", 0) for p in (pos, pos + 100)]
    _identity(data, list(reversed(rows)) if reverse else rows)
    result = analyze(data, call(tmp_path / "input.vcf", pos + 100), "GRCh38", ("pharmgkb",))
    assert len(result.annotations) == 1
    assert warning_count(result, "identity_locus_conflict") == 0


def test_matching_clinvar_row_cannot_override_contradictory_identity(tmp_path):
    from validated_call_helpers import sql

    data = tmp_path / "data"
    pos = setup(data, "GRCh38")
    sql(data, "clinvar.GRCh38.sqlite", "UPDATE clinvar_variants SET position=position+100")
    result = analyze(data, call(tmp_path / "input.vcf", pos + 100), "GRCh38", ("clinvar",))
    assert not result.annotations
    assert warning_count(result, "identity_locus_conflict") == 1

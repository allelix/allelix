"""Named reference calls retain chip parity when identity is absent or contradicts."""

import pytest
from test_annotation_parity_and_named_disclosures import _vcf
from test_haploid_and_homozygous_warning import _array
from validated_call_helpers import (
    BUILDS,
    RSID,
    analyze,
    reference_effect,
    signatures,
    sql,
    warning_count,
)

from allelix.parsers.twentythreeandme import TwentyThreeAndMeParser


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("gt", ("0/0", "0|0", "0"))
@pytest.mark.parametrize(
    "identity_state", ("present", "absent", "missing_rsid", "unplaced", "conflict")
)
def test_named_reference_call_matches_chip_identity_decision(tmp_path, build, gt, identity_state):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    if identity_state == "absent":
        (data / "rsid_identity.sqlite").unlink()
    elif identity_state == "missing_rsid":
        sql(data, "rsid_identity.sqlite", "DELETE FROM rsid_identity")
    elif identity_state == "unplaced":
        column = "pos37" if build == "GRCh37" else "pos38"
        sql(data, "rsid_identity.sqlite", f"UPDATE rsid_identity SET {column}=0")
    elif identity_state == "conflict":
        pos += 100

    chip = _array(tmp_path / "chip.txt", [(RSID, "1", pos, "G" if gt == "0" else "GG")])
    control = analyze(data, chip, build, ("pharmgkb",), parser=TwentyThreeAndMeParser())
    path = _vcf(tmp_path / "input.vcf", [(RSID, "1", pos, "G", "A", gt, ".")])
    result = analyze(data, path, build, ("pharmgkb",))
    expected_count = int(identity_state == "conflict")
    assert warning_count(control, "identity_locus_conflict") == expected_count
    assert bool(control.annotations) == (not expected_count)
    assert signatures(result) == signatures(control)
    assert warning_count(result, "identity_locus_conflict") == expected_count

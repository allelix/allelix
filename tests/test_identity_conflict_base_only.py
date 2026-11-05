# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Only a base call against a base placement can contradict identity (ADR-0063)."""

import pytest
from test_annotation_parity_and_named_disclosures import _identity, _vcf
from test_haploid_and_homozygous_warning import _array
from validated_call_helpers import BUILDS, RSID, analyze, call, setup, warning_count

from allelix.parsers.twentythreeandme import TwentyThreeAndMeParser


def _placements(data, rows):
    (data / "rsid_identity.sqlite").unlink()
    _identity(data, [(RSID, "1", pos, ref, alt, pos, ref, alt, 0) for pos, ref, alt in rows])


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("identity_row", ["base", "indel"])
@pytest.mark.parametrize("offset", [1, 3, 40])
def test_a_vcf_indel_call_is_never_a_contradiction(tmp_path, build, identity_row, offset):
    data = tmp_path / "data"
    pos = setup(data, build)
    _placements(data, [(pos, "G", "A") if identity_row == "base" else (pos, "GT", "G")])
    path = _vcf(tmp_path / "input.vcf", [(RSID, "1", pos + offset, "GT", "G", "0/1", ".")])
    result = analyze(data, path, build, ("pharmgkb",))
    assert warning_count(result, "identity_locus_conflict") == 0


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("genotype", ["DI", "DD", "II"])
@pytest.mark.parametrize("identity_row", ["base", "indel"])
def test_a_chip_di_call_is_never_a_contradiction(tmp_path, build, genotype, identity_row):
    data = tmp_path / "data"
    pos = setup(data, build)
    _placements(data, [(pos, "G", "A") if identity_row == "base" else (pos, "GT", "G")])
    path = _array(tmp_path / "input.txt", [(RSID, "1", pos + 1, genotype)])
    result = analyze(data, path, build, ("pharmgkb",), parser=TwentyThreeAndMeParser())
    assert warning_count(result, "identity_locus_conflict") == 0


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("format_name", ["vcf", "chip"])
def test_identity_indel_rows_alone_are_identity_unavailable(tmp_path, build, format_name):
    data = tmp_path / "data"
    pos = setup(data, build)
    _placements(data, [(pos + 500, "GT", "G"), (pos + 500, "G", "GA")])
    if format_name == "chip":
        path = _array(tmp_path / "input.txt", [(RSID, "1", pos, "AG")])
        result = analyze(data, path, build, ("pharmgkb",), parser=TwentyThreeAndMeParser())
    else:
        result = analyze(data, call(tmp_path / "input.vcf", pos), build, ("pharmgkb",))
    assert warning_count(result, "identity_locus_conflict") == 0
    assert len(result.annotations) == 1


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("format_name", ["vcf", "chip"])
def test_a_base_call_off_its_base_placement_stays_a_contradiction(tmp_path, build, format_name):
    data = tmp_path / "data"
    pos = setup(data, build)
    # An indel row elsewhere does not make the base placement unavailable.
    _placements(data, [(pos, "G", "A"), (pos + 500, "GT", "G")])
    if format_name == "chip":
        path = _array(tmp_path / "input.txt", [(RSID, "1", pos + 1, "AG")])
        result = analyze(data, path, build, ("pharmgkb",), parser=TwentyThreeAndMeParser())
    else:
        result = analyze(data, call(tmp_path / "input.vcf", pos + 1), build, ("pharmgkb",))
    assert warning_count(result, "identity_locus_conflict") == 1
    assert not result.annotations


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("shape", ["prefix", "suffix", "both"])
def test_a_padded_base_call_stands_at_its_trimmed_site(tmp_path, build, shape):
    data = tmp_path / "data"
    pos = setup(data, build)
    result = analyze(data, call(tmp_path / "input.vcf", pos, shape=shape), build, ("pharmgkb",))
    assert warning_count(result, "identity_locus_conflict") == 0
    assert len(result.annotations) == 1


@pytest.mark.parametrize("build", BUILDS)
def test_a_padded_base_placement_is_read_at_its_trimmed_site(tmp_path, build):
    data = tmp_path / "data"
    pos = setup(data, build)
    _placements(data, [(pos - 1, "CG", "CA")])
    result = analyze(data, call(tmp_path / "input.vcf", pos), build, ("pharmgkb",))
    assert warning_count(result, "identity_locus_conflict") == 0
    assert len(result.annotations) == 1


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("shape", ["prefix", "suffix", "both"])
def test_a_padded_reference_only_call_has_no_base_form(tmp_path, build, shape):
    data = tmp_path / "data"
    pos = setup(data, build)
    path = call(tmp_path / "input.vcf", pos + 100, shape=shape, gt="0/0")
    result = analyze(data, path, build, ("pharmgkb",))
    assert warning_count(result, "identity_locus_conflict") == 0

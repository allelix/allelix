"""Explicit reference calls reach source matching; absent sites are never invented."""

import pytest
from test_annotation_parity_and_named_disclosures import _vcf
from test_haploid_and_homozygous_warning import _array
from validated_call_helpers import BUILDS, RSID, analyze, reference_effect, signatures

from allelix.parsers.twentythreeandme import TwentyThreeAndMeParser
from allelix.parsers.vcf import VcfParser

SOURCES = ("gwas", "pharmgkb", "snpedia")


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("gt", ["0/0", "0|0", "0"])
@pytest.mark.parametrize("named", [False, True])
@pytest.mark.parametrize("alt,padding", [("A", False), (".", False), ("A", True)])
def test_explicit_reference_calls_equal_chip_calls(tmp_path, build, gt, named, alt, padding):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    array = _array(tmp_path / "chip.txt", [(RSID, "1", pos, "G" if gt == "0" else "GG")])
    control = analyze(data, array, build, SOURCES, parser=TwentyThreeAndMeParser())
    assert {a.source for a in control.annotations} == set(SOURCES)
    path = _vcf(
        tmp_path / "input.vcf",
        [
            (
                RSID if named else ".",
                "1",
                pos,
                "GC" if padding else "G",
                "AC" if padding else alt,
                gt,
                ".",
            )
        ],
    )
    result = analyze(data, path, build, SOURCES)
    assert signatures(result) == signatures(control)
    assert all(a.haploid == (gt == "0") for a in result.annotations)


@pytest.mark.parametrize("kind", ["absent", "block", "no_call", "low_gq", "filtered"])
@pytest.mark.parametrize("build", BUILDS)
def test_absence_blocks_and_rejected_reference_calls_are_not_fabricated(tmp_path, kind, build):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    if kind == "absent":
        row = ("rs900260001", "2", 100, "G", "A", "0/0", ".")
    else:
        row = (
            RSID,
            "1",
            pos,
            "G",
            "<NON_REF>" if kind == "block" else "A",
            "./." if kind == "no_call" else "0/0",
            f"END={pos + 10}" if kind == "block" else ".",
        )
    path = _vcf(tmp_path / "input.vcf", [row], block=kind == "block")
    if kind == "low_gq":
        path.write_text(path.read_text().replace("\tGT\t0/0", "\tGT:GQ\t0/0:5"))
    if kind == "filtered":
        path.write_text(path.read_text().replace("\tPASS\t", "\tLowQual\t"))
    result = analyze(data, path, build, SOURCES, parser=VcfParser())
    assert not result.annotations

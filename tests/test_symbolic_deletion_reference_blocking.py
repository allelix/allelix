"""A carried <DEL>'s END span cannot simultaneously support a reference call."""

import json

import pytest
from click.testing import CliRunner
from test_annotation_parity_and_named_disclosures import _vcf
from validated_call_helpers import BUILDS, RSID, analyze, reference_effect

from allelix.cli import main
from allelix.parsers.vcf import VcfParser

SOURCES = ("gwas", "snpedia", "pharmgkb")


def _reference(pos, gt="0/0"):
    return RSID, "1", pos, "G", "A", gt, "."


def _deletion(pos, gt="0/1", *, end=None, alt="<DEL>", chrom="1"):
    return ".", chrom, pos - 10, "G", alt, gt, f"END={pos + 10 if end is None else end}"


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("gt", ("0/1", "1/1", "1", "0|1", "1|0"))
@pytest.mark.parametrize("at_end", (False, True))
def test_carried_symbolic_deletion_blocks_reference_through_its_end(
    tmp_path, build, reverse, gt, at_end
):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    alone = analyze(data, _vcf(tmp_path / "alone.vcf", [_reference(pos)]), build, SOURCES)
    assert {a.source for a in alone.annotations} == set(SOURCES)
    rows = [_deletion(pos, gt, end=pos if at_end else pos + 10), _reference(pos)]
    result = analyze(
        data, _vcf(tmp_path / "input.vcf", rows[::-1] if reverse else rows), build, SOURCES
    )
    assert not result.annotations, (
        "A called END deletion contradicts the reference genotype inside it"
    )


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize(
    "alt,gt",
    (
        ("<DEL>,A", "1/2"),
        ("A,<DEL>", "0/2"),
        ("<DEL>,<NON_REF>", "0/1"),
        ("<NON_REF>,<DEL>", "0/2"),
    ),
)
def test_the_selected_alt_controls_whether_a_symbolic_deletion_blocks(tmp_path, build, alt, gt):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    source = _vcf(tmp_path / "input.vcf", [_deletion(pos, gt, alt=alt), _reference(pos)])
    result = analyze(data, source, build, SOURCES)
    assert not result.annotations


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize(
    "control",
    (
        "no_deletion",
        "other_chromosome",
        "before",
        "after_end",
        "uncarried_del",
        "other_alt_carried",
        "dup",
        "inv",
        "cnv",
        "non_ref",
        "missing_end",
        "bad_end",
        "backwards_end",
        "literal_deletion",
    ),
)
def test_scope_controls_preserve_nonoverlaps_and_existing_literal_behavior(
    tmp_path, build, control
):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    blocker = _deletion(pos)
    if control == "no_deletion":
        blocker = None
    elif control == "other_chromosome":
        blocker = _deletion(pos, chrom="2")
    elif control == "before":
        blocker = (".", "1", pos + 1, "G", "<DEL>", "0/1", f"END={pos + 10}")
    elif control == "after_end":
        blocker = _deletion(pos, end=pos - 1)
    elif control == "uncarried_del":
        blocker = _deletion(pos, "0/0")
    elif control == "other_alt_carried":
        blocker = _deletion(pos, "0/2", alt="<DEL>,A")
    elif control in ("dup", "inv", "cnv", "non_ref"):
        blocker = _deletion(pos, alt=f"<{control.upper()}>")
    elif control == "missing_end":
        blocker = (*blocker[:6], ".")
    elif control == "bad_end":
        blocker = (*blocker[:6], "END=not-a-position")
    elif control == "backwards_end":
        blocker = _deletion(pos, end=pos - 20)
    else:
        blocker = (".", "1", pos - 10, "G" * 21, "G", "0/1", ".")
    rows = ([blocker] if blocker is not None else []) + [_reference(pos)]
    result = analyze(data, _vcf(tmp_path / "input.vcf", rows), build, SOURCES)
    expected = set() if control == "literal_deletion" else set(SOURCES)
    assert {a.source for a in result.annotations} == expected


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("sample,blocked", (("DEL_SAMPLE", True), ("REF_SAMPLE", False)))
def test_another_samples_deletion_cannot_block_the_selected_sample(
    tmp_path, build, sample, blocked
):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    source = tmp_path / "samples.vcf"
    source.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDEL_SAMPLE\tREF_SAMPLE\n"
        f"1\t{pos - 10}\t.\tG\t<DEL>\t.\tPASS\tEND={pos + 10}\tGT\t0/1\t0/0\n"
        f"1\t{pos}\t{RSID}\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/0\n"
    )
    result = analyze(data, source, build, SOURCES, parser=VcfParser(sample=sample))
    assert {a.source for a in result.annotations} == (set() if blocked else set(SOURCES))


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("explicit", (False, True))
@pytest.mark.parametrize("reference_gt", ("0/0", "0"))
def test_cli_and_annotated_vcf_do_not_publish_the_conflicting_reference_finding(
    tmp_path, build, explicit, reference_gt
):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    source = _vcf(tmp_path / "input.vcf", [_deletion(pos), _reference(pos, reference_gt)])
    length = 249250621 if build == "GRCh37" else 248956422
    source.write_text(
        source.read_text().replace(
            "##fileformat=VCFv4.2\n", f"##fileformat=VCFv4.2\n##contig=<ID=1,length={length}>\n"
        )
    )
    # Isolate the three reference-genotype sources while using the actual CLI.
    (data / "config.toml").write_text(
        "[sources]\nclinvar=false\npharmgkb=true\nsnpedia=true\ngwas=true\ngnomad=false\nalphamissense=false\ncadd=false\nrsid_identity=true\n"
    )
    output = tmp_path / "report.json"
    annotated = tmp_path / "annotated.vcf"
    args = [
        "analyze",
        str(source),
        "--data-dir",
        str(data),
        "--no-update",
        "--min-magnitude",
        "0",
        "--gwas-min-magnitude",
        "0",
        "--output",
        str(output),
        "--vcf-out",
        str(annotated),
    ]
    if explicit:
        args += ["--build", build.lower()]
    run = CliRunner().invoke(main, args)
    assert run.exit_code == 0, (run.output, run.exception)
    report = json.loads(output.read_text())
    assert report["input"]["build"] == build
    assert not [r for r in report["annotations"] if r["source"] in SOURCES]
    snv = next(
        line.split("\t")
        for line in annotated.read_text().splitlines()
        if not line.startswith("#") and line.split("\t")[1] == str(pos)
    )
    assert not any(f"ALLELIX_{source}=" in snv[7] for source in ("GWAS", "SNPEDIA", "CLINPGX"))


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("reverse", (False, True))
@pytest.mark.parametrize("reference_gt", ("0/0", "0"))
@pytest.mark.parametrize("kind", ("symbolic", "literal"))
def test_clinvar_reference_carried_finding_obeys_the_same_block(
    tmp_path, build, reverse, reference_gt, kind
):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    reference = _reference(pos, reference_gt)
    alone = analyze(data, _vcf(tmp_path / "alone.vcf", [reference]), build, ("clinvar",))
    assert [a.source for a in alone.annotations] == ["clinvar"]
    deletion = (
        _deletion(pos) if kind == "symbolic" else (".", "1", pos - 10, "G" * 21, "G", "0/1", ".")
    )
    rows = [deletion, reference]
    result = analyze(
        data, _vcf(tmp_path / "input.vcf", rows[::-1] if reverse else rows), build, ("clinvar",)
    )
    assert not result.annotations

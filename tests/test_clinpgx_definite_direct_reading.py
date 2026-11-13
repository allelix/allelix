"""A definite direct site cannot admit a second reading through the legacy REF rule."""

import contextlib
import json

import pytest
from click.testing import CliRunner
from test_clinpgx_ambiguous_reading_abstention import (
    BUILDS,
    COUNT,
    LENGTH,
    POS,
    RS,
    _call,
    _data,
    _labels,
)

from allelix.annotators import get_annotators
from allelix.cli import main
from allelix.utils.allele import STRAND_DIRECT, forward_strand_reading

CASES = (
    ("GAC", "AC", "GT", "G"),
    ("CAG", "AG", "CT", "C"),
    ("ACT", "CT", "AG", "A"),
    ("CGT", "GT", "AC", "C"),
)


def _run(data, calls, path, *, switched_off=()):
    # The CLI supplies configuration switches to get_annotators; unit callers
    # must supply them too. Merely writing config.toml does not disable a source.
    with contextlib.ExitStack() as stack:
        registry = get_annotators(data, switched_off=switched_off)
        for annotator in registry:
            stack.enter_context(annotator)
        pgx = next(a for a in registry if a.name == "pharmgkb")
        assert pgx.is_ready()
        assert pgx.run_counts().get(COUNT, 0) == 0
        findings = (
            list(pgx.batch_annotate(calls))
            if path == "batch"
            else [a for call in calls for a in pgx.annotate(call)]
        )
        return findings, pgx.run_counts().get(COUNT, 0)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("path", ("single", "batch"))
@pytest.mark.parametrize("site,source,opposite,ref", CASES)
@pytest.mark.parametrize("evidence", ("direct", "missing", "disabled", "elsewhere"))
def test_only_available_definite_direct_evidence_overrides_the_legacy_ref_reading(
    tmp_path, build, path, site, source, opposite, ref, evidence
):
    assert forward_strand_reading(set(source), frozenset(site)) == STRAND_DIRECT
    # The helper's no_fit mode installs the given site; the supplied alleles
    # here deliberately make that site DIRECT rather than the helper default.
    data = _data(
        tmp_path / "data",
        rows=(source,),
        mode="no_fit" if evidence == "direct" else evidence,
        site=site,
    )
    findings, count = _run(
        data,
        [_call(build, opposite, ref)],
        path,
        switched_off=frozenset({"rsid_identity"}) if evidence == "disabled" else (),
    )
    assert _labels(findings) == ([] if evidence == "direct" else [f"row {source}"]), (
        "A decided direct reading must not fall back to complement matching by REF"
    )
    assert count == 0, "A definite reading (or a single legacy match) is not ambiguous"


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("path", ("single", "batch"))
@pytest.mark.parametrize("site,source,opposite,ref", CASES)
def test_definite_direct_site_keeps_the_written_carrier(
    tmp_path, build, path, site, source, opposite, ref
):
    data = _data(tmp_path / "data", rows=(source,), mode="no_fit", site=site)
    findings, count = _run(data, [_call(build, source, ref)], path)
    assert _labels(findings) == [f"row {source}"]
    assert count == 0


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("auto", (False, True))
def test_vcf_forward_call_does_not_get_a_complement_row_at_a_direct_site(tmp_path, build, auto):
    data = _data(tmp_path / "data", rows=("AC",), mode="no_fit", site="GAC")
    source = tmp_path / "input.vcf"
    source.write_text(
        "##fileformat=VCFv4.2\n"
        f"##contig=<ID=1,length={LENGTH[build]}>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        f"1\t{POS[build]}\t{RS}\tG\tT\t.\tPASS\t.\tGT\t0/1\n"
    )
    output = tmp_path / "report.json"
    run = CliRunner().invoke(
        main,
        [
            "analyze",
            str(source),
            "--data-dir",
            str(data),
            "--no-update",
            "--build",
            "auto" if auto else build.lower(),
            "--min-magnitude",
            "0",
            "--output",
            str(output),
        ],
    )
    assert run.exit_code == 0, run.output
    report = json.loads(output.read_text())
    assert report["input"]["build"] == build
    assert not [a for a in report["annotations"] if a["source"] == "pharmgkb"], (
        "The VCF GT call is not the direct-strand ClinPGx AC genotype"
    )

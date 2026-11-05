"""Selected sample FT and site FILTER share one quality gate."""

import pytest
from validated_call_helpers import BUILDS, RSID, analyze, setup

from allelix.parsers.vcf import VcfParser


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("gt", ["0/1", "0|1", "1", "./1", "0/0"])
@pytest.mark.parametrize(
    "ft,site,labels",
    [
        ("LowAD", "PASS", {"LowAD"}),
        ("LowAD;LowDP", "LowQual", {"LowAD", "LowDP", "LowQual"}),
        ("PASS", "PASS", set()),
        (".", ".", set()),
    ],
)
@pytest.mark.parametrize("keep", [False, True])
def test_sample_and_site_filter_count_once(tmp_path, build, gt, ft, site, labels, keep):
    data = tmp_path / "data"
    pos = setup(data, build)
    path = tmp_path / "input.vcf"
    path.write_text(
        "##fileformat=VCFv4.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        f"1\t{pos}\t{RSID}\tG\tA\t.\t{site}\t.\tGT:FT:GQ:DP\t{gt}:{ft}:99:30\n"
    )
    result = analyze(data, path, build, parser=VcfParser(keep_filtered=keep))
    quality = result.vcf_quality
    assert quality.filtered_rows == bool(labels)
    assert dict(quality.filtered_by_label) == {label: 1 for label in labels}
    if labels and not keep:
        assert not result.annotations
    elif gt != "0/0":
        assert result.annotations


@pytest.mark.parametrize("selected", ["GOOD", "BAD"])
@pytest.mark.parametrize("order", [False, True])
def test_only_selected_samples_ft_is_consulted(tmp_path, selected, order):
    data = tmp_path / "data"
    pos = setup(data, "GRCh38")
    samples = [("GOOD", "0/1:PASS"), ("BAD", "0/1:LowAD")]
    if order:
        samples.reverse()
    path = tmp_path / "multi.vcf"
    path.write_text(
        "##fileformat=VCFv4.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(n for n, _ in samples)
        + "\n"
        + f"1\t{pos}\t{RSID}\tG\tA\t.\tPASS\t.\tGT:FT\t"
        + "\t".join(v for _, v in samples)
        + "\n"
    )
    result = analyze(data, path, "GRCh38", parser=VcfParser(sample=selected))
    assert bool(result.annotations) == (selected == "GOOD")
    assert result.vcf_quality.filtered_rows == (selected == "BAD")


@pytest.mark.parametrize("named", [False, True])
@pytest.mark.parametrize("shape", ["plain", "prefix"])
@pytest.mark.parametrize("keep", [False, True])
def test_filter_overrides_survive_normalization_and_resolution(tmp_path, named, shape, keep):
    from test_annotation_parity_and_named_disclosures import _run
    from validated_call_helpers import call

    data = tmp_path / "data"
    pos = setup(data, "GRCh38")
    path = call(tmp_path / "input.vcf", pos, rsid=RSID if named else ".", shape=shape)
    path.write_text(path.read_text().replace("\tGT\t0/1", "\tGT:FT\t0/1:LowAD"))
    output = tmp_path / "annotated.vcf"
    extra = ("--min-magnitude", "0", "--vcf-out", str(output))
    if keep:
        extra += ("--keep-filtered",)
    report, _ = _run(tmp_path, data, path, "GRCh38", "json", extra=extra)
    assert bool(report["annotations"]) == keep
    assert report["input"]["vcf_quality"]["filtered_rows"] == 1
    row = next(r for r in output.read_text().splitlines() if not r.startswith("#"))
    assert ("ALLELIX_CLINVAR=" in row) == keep


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("indexed", [False, True])
@pytest.mark.parametrize("gt", ["0/1", "1", "0/0", "./1"])
def test_extract_observes_sample_filter_on_both_paths(tmp_path, monkeypatch, build, indexed, gt):
    from click.testing import CliRunner
    from test_extract_reference_calls import _input, _row

    from allelix.cli import main, utility

    data = tmp_path / "data"
    pos = setup(data, build)
    row = _row(pos, gt, name=RSID).replace("\tGT:GQ:DP\t", "\tGT:GQ:DP:FT\t")
    row = row.rstrip("\n") + ":LowAD\n"
    path = _input(tmp_path, build, [row], indexed=indexed)
    found = []
    monkeypatch.setattr(utility, "_render_extract_table", lambda p, w, f: found.append(f))
    result = CliRunner().invoke(
        main,
        ["extract", str(path), "--snps", RSID, "--build", build.lower(), "--data-dir", str(data)],
    )
    assert result.exit_code == 0, (result.output, result.exception)
    assert len(found) == 1 and RSID in found[0]
    record = found[0][RSID]
    assert record.is_no_call and not record.is_half_call

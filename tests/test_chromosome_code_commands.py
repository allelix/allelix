"""Chromosome aliases must retain findings and raw VCF record/header spelling."""

import json

import pytest
from click.testing import CliRunner
from test_annotation_parity_and_named_disclosures import _clinvar, _gnomad_site, _identity, _vcf

from allelix.cli import main

ALIASES = [("23", "X"), ("24", "Y"), ("25", "X"), ("26", "MT"), ("XY", "X")]
RSID = "rs9000000091"
POS = 65001
CONFIG = "[sources]\n" + "".join(
    f"{name} = {str(name in ('clinvar', 'rsid_identity')).lower()}\n"
    for name in (
        "clinvar",
        "pharmgkb",
        "snpedia",
        "gwas",
        "gnomad",
        "cadd",
        "alphamissense",
        "rsid_identity",
    )
)


def _data(tmp_path, canonical):
    data = tmp_path / "data"
    data.mkdir()
    (data / "config.toml").write_text(CONFIG)
    for build in ("GRCh37", "GRCh38"):
        _clinvar(
            data,
            build,
            [(RSID, canonical, POS, "G", "A", "SYN", "Pathogenic", "Pathogenic", 1, 1)],
        )
    _identity(data, [(RSID, canonical, POS, "G", "A", POS, "G", "A", 0)])
    return data


def _input(path, chromosome, gt, form):
    if form == "vcf":
        _vcf(path, [(RSID, chromosome, POS, "G", "A", gt, ".")])
        path.write_text(path.read_text().replace("#CHROM", f"##contig=<ID={chromosome}>\n#CHROM"))
    else:
        call = {"0/1": "GA", "1/1": "AA", "1": "A"}[gt]
        path.write_text(
            "Name,Variation,Chromosome,Position,Strand,YourCode\n"
            + f"{RSID},G/A,{chromosome},{POS},+,{call}\n"
        )
    return path


def _analyze(tmp_path, data, source, build, label, vcf_out=False):
    out = tmp_path / f"{label}.json"
    args = [
        "analyze",
        str(source),
        "--build",
        build.lower(),
        "--data-dir",
        str(data),
        "--no-update",
        "--output",
        str(out),
    ]
    annotated = tmp_path / f"{label}.annotated.vcf"
    if vcf_out:
        args += ["--vcf-out", str(annotated)]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, (result.output, result.exception)
    return json.loads(out.read_text()), annotated


@pytest.mark.parametrize(("code", "canonical"), ALIASES)
@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("gt", ["0/1", "1/1", "1"])
@pytest.mark.parametrize("form", ["vcf", "decodeme"])
def test_alias_has_same_clinvar_findings_as_letters(tmp_path, code, canonical, build, gt, form):
    data = _data(tmp_path, canonical)
    reports = []
    for label, chrom in [("control", canonical), ("coded", code)]:
        source = _input(tmp_path / f"{label}.{form}", chrom, gt, form)
        report, _ = _analyze(tmp_path, data, source, build, label)
        reports.append(report)
    control, actual = reports
    assert len(control["annotations"]) == 1
    assert actual["annotations"] == control["annotations"]
    assert actual["report_warnings"] == control["report_warnings"]


@pytest.mark.parametrize(("code", "canonical"), ALIASES)
@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
def test_vcf_out_uses_canonical_lookup_and_preserves_raw_chrom(tmp_path, code, canonical, build):
    data = _data(tmp_path, canonical)
    source = _input(tmp_path / "source.vcf", code, "1", "vcf")
    report, annotated = _analyze(tmp_path, data, source, build, "report", True)
    lines = annotated.read_text().splitlines()
    assert f"##contig=<ID={code}>" in lines
    records = [line.split("\t") for line in lines if not line.startswith("#")]
    assert len(records) == 1 and records[0][0] == code
    assert "ALLELIX_CLINVAR=" in records[0][7]
    assert len(report["annotations"]) == 1
    assert report["annotations"][0]["chromosome"] == canonical


@pytest.mark.parametrize(("code", "canonical"), ALIASES)
@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
def test_plink_codes_and_alleles_match_canonical_input(tmp_path, code, canonical, build):
    data = _data(tmp_path, canonical)
    _gnomad_site(data, RSID, canonical, POS, "G", "A", 0.1)
    (data / "config.toml").write_text(CONFIG.replace("gnomad = false", "gnomad = true"))
    outputs = []
    for label, chrom in [("control", canonical), ("coded", code)]:
        source = _input(tmp_path / f"{label}.vcf", chrom, "0/1", "vcf")
        prefix = tmp_path / label
        result = CliRunner().invoke(
            main,
            [
                "export",
                "plink",
                str(source),
                "--build",
                build.lower(),
                "--data-dir",
                str(data),
                "--output-prefix",
                str(prefix),
            ],
        )
        assert result.exit_code == 0, (result.output, result.exception)
        outputs.append(
            (prefix.with_suffix(".bim").read_bytes(), prefix.with_suffix(".bed").read_bytes())
        )
    expected_bim, expected_bed = outputs[0]
    if code in {"25", "XY"}:
        # Lookup uses X; PLINK must retain an explicitly supplied PAR designation.
        chrom, separator, rest = expected_bim.partition(b"\t")
        assert chrom == b"23" and separator
        expected_bim = b"25" + separator + rest
    assert outputs[1] == (expected_bim, expected_bed)


@pytest.mark.parametrize("code", ["23", "25", "XY"])
def test_x_one_copy_plausibility_uses_x_rule(tmp_path, code):
    data = _data(tmp_path, "X")
    _gnomad_site(data, RSID, "X", POS, "G", "A", 0.0)
    (data / "config.toml").write_text(CONFIG.replace("gnomad = false", "gnomad = true"))
    reports = []
    for label, chrom in [("control", "X"), ("coded", code)]:
        source = _input(tmp_path / f"{label}.vcf", chrom, "0/1", "vcf")
        report, _ = _analyze(tmp_path, data, source, "GRCh38", label)
        reports.append(report)
    assert reports[0]["annotations"][0]["plausibility"] == "suspicious"
    assert reports[1]["annotations"] == reports[0]["annotations"]

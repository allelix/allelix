# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""An unusable tabix sidecar cannot prevent a readable VCF being extracted."""

import gzip

import pysam
import pytest


@pytest.mark.parametrize("build", ["grch37", "grch38"])
@pytest.mark.parametrize("compression", ["gzip", "bgzf"])
@pytest.mark.parametrize("index_state", ["corrupt", "unreadable"])
@pytest.mark.parametrize("call,genotype", [("0/1", "A/G"), ("0/0", "A/A")])
def test_bad_index_warns_and_returns_the_sequential_selected_sample(
    rc_case,
    build,
    compression,
    index_state,
    call,
    genotype,
):
    case = rc_case
    header = case.header.replace("\tSAMPLE\n", "\tFIRST\tSECOND\n")
    if build == "grch38":
        header = header.replace("249250621", "248956422")
    case.source.write_text(header + case.row.replace("GT\t0/1\n", f"GT\t./.\t{call}\n"))
    source = case.root / "input.vcf.gz"
    if compression == "bgzf":
        pysam.tabix_compress(str(case.source), str(source))
    else:
        with gzip.open(source, "wb") as stream:
            stream.write(case.source.read_bytes())
    args = [
        "extract",
        source,
        "--data-dir",
        case.data,
        "--build",
        build,
        "--sample",
        "SECOND",
        "--snps",
        case.rsid,
    ]
    sequential = case.invoke(args)
    assert sequential.exit_code == 0, (sequential.output, repr(sequential.exception))
    assert genotype in sequential.output, sequential.output
    expected = next(line.strip() for line in sequential.output.splitlines() if case.rsid in line)
    index = source.with_suffix(".gz.tbi")
    if compression == "bgzf":
        pysam.tabix_index(str(source), preset="vcf")
    else:
        index.write_bytes(b"unusable index\n")
    if index_state == "corrupt":
        index.write_bytes(b"not a tabix index\n")
    else:
        index.chmod(0)
    before = source.read_bytes()
    result = case.invoke(args)
    assert result.exit_code == 0, (result.output, repr(result.exception))
    assert expected in result.output, result.output
    assert "warn" in result.output.lower(), result.output
    assert index.name in result.output, result.output
    assert "Traceback" not in result.output
    assert source.read_bytes() == before


@pytest.mark.parametrize("build", ["grch37", "grch38"])
def test_valid_index_stays_quiet_and_preserves_the_call(rc_case, build):
    case = rc_case
    if build == "grch38":
        case.source.write_text(case.source.read_text().replace("249250621", "248956422"))
    source = case.root / "valid.vcf.gz"
    pysam.tabix_compress(str(case.source), str(source))
    pysam.tabix_index(str(source), preset="vcf")
    result = case.invoke(
        ["extract", source, "--data-dir", case.data, "--build", build, "--snps", case.rsid]
    )
    assert result.exit_code == 0, (result.output, repr(result.exception))
    assert case.rsid in result.output and "A/G" in result.output
    assert "warn" not in result.output.lower(), result.output


@pytest.mark.parametrize("compression", ["gzip", "bgzf"])
def test_unreadable_vcf_is_a_clean_path_error(rc_case, compression):
    case = rc_case
    source = case.root / "unreadable.vcf.gz"
    if compression == "bgzf":
        pysam.tabix_compress(str(case.source), str(source))
    else:
        with gzip.open(source, "wb") as stream:
            stream.write(case.source.read_bytes())
    source.with_suffix(".gz.tbi").write_bytes(b"bad index")
    source.chmod(0)
    result = case.invoke(
        ["extract", source, "--data-dir", case.data, "--build", "grch37", "--snps", case.rsid]
    )
    assert result.exit_code in (1, 2), result.output
    assert isinstance(result.exception, SystemExit), repr(result.exception)
    assert "error" in result.output.lower() and source.name in result.output
    assert "Traceback" not in result.output

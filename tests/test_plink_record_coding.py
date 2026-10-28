# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""PLINK exports preserve each input record's observed alleles."""

from __future__ import annotations

import contextlib
import os
import re
import shutil
import sqlite3
import subprocess
from collections import Counter
from pathlib import Path

import pytest
from click.testing import CliRunner

from allelix.annotators.gnomad import GnomadAnnotator
from allelix.cli import main
from allelix.databases._versions import GNOMAD_SCHEMA_VERSION
from allelix.databases.schema import GNOMAD_SCHEMA

ROWS = [
    ("2", 1000, "G", "A", "rs900961001", 0.001),
    ("2", 1000, "G", "T", "rs900961001", 0.001),
]


def create_cache(directory: Path) -> Path:
    """Create a new fixture cache, refusing to overwrite an existing cache."""
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "gnomad.sqlite"
    if path.exists():
        raise FileExistsError(path)
    with contextlib.closing(sqlite3.connect(path)) as conn:
        conn.executescript(GNOMAD_SCHEMA)
        conn.executemany(
            "INSERT INTO gnomad_frequencies (chrom,pos,ref,alt,rsid,af) VALUES (?,?,?,?,?,?)",
            ROWS,
        )
        conn.execute(
            "INSERT INTO database_versions "
            "(name,source_url,version,downloaded_at,record_count,remote_signal,local_version_tag) "
            "VALUES ('gnomad','test://synthetic','synthetic',"
            "'2026-10-25',?,NULL,?)",
            (len(ROWS), f"sv:{GNOMAD_SCHEMA_VERSION}"),
        )
        conn.commit()
    with GnomadAnnotator(directory) as annotator:
        assert annotator.is_ready(), "The synthetic gnomAD cache must be usable"
    return path


FIXTURES = Path(__file__).parent / "fixtures"
EXPECTED = {
    "joint": Counter({("A", "T"): 1}),
    "split": Counter({("A", "G"): 1, ("G", "T"): 1}),
}


def _export(tmp_path: Path, shape: str, with_gnomad: bool) -> Path:
    data = tmp_path / "data"
    data.mkdir()
    if with_gnomad:
        create_cache(data)
    prefix = tmp_path / "exported"
    result = CliRunner().invoke(
        main,
        [
            "export",
            "plink",
            str(FIXTURES / f"plink_{shape}.vcf"),
            "--data-dir",
            str(data),
            "--output-prefix",
            str(prefix),
        ],
        env={"TERM": "xterm", "COLUMNS": "130"},
    )
    assert result.exit_code == 0, result.output
    return prefix


def _decode_one_sample_bed(prefix: Path) -> Counter:
    """Decode PLINK1 SNP-major data independently of Allelix's writer.

    For one sample each variant occupies one byte. Low bits 00/10/11
    encode A1/A1, A1/A2, A2/A2; 01 is missing. Compare bases, not A1/A2
    orientation or record order, and preserve repeated IDs as separate rows.
    """
    fam = Path(str(prefix) + ".fam").read_text().splitlines()
    assert len(fam) == 1
    bim = [line.split() for line in Path(str(prefix) + ".bim").read_text().splitlines()]
    bed = Path(str(prefix) + ".bed").read_bytes()
    assert bed[:3] == b"\x6c\x1b\x01", "Expected SNP-major PLINK1 BED"
    assert len(bed) == 3 + len(bim)
    calls = []
    for fields, byte in zip(bim, bed[3:], strict=True):
        assert len(fields) == 6
        a1, a2 = fields[4:6]
        call = {0: (a1, a1), 1: (".", "."), 2: (a1, a2), 3: (a2, a2)}[byte & 3]
        calls.append(tuple(sorted(call)))
    return Counter(calls)


@pytest.mark.parametrize("shape", ["joint", "split"])
@pytest.mark.parametrize("with_gnomad", [False, True], ids=["without_gnomad", "with_gnomad"])
def test_plink_bed_preserves_carried_alleles(tmp_path: Path, shape: str, with_gnomad: bool):
    prefix = _export(tmp_path, shape, with_gnomad)
    assert _decode_one_sample_bed(prefix) == EXPECTED[shape]


def _plink2_executable() -> str:
    executable = os.environ.get("PLINK2_BIN") or shutil.which("plink2")
    if executable is None:
        pytest.fail("PLINK2 is required: set PLINK2_BIN or put plink2 on PATH; no auto-download")
    return executable


def _decode_roundtrip_vcf(path: Path) -> Counter:
    calls = []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        assert len(fields) == 10, "Expected one sample in PLINK2's exported VCF"
        alleles = [fields[3], *fields[4].split(",")]
        gt_column = fields[8].split(":").index("GT")
        gt = fields[9].split(":")[gt_column]
        indices = re.split(r"[/|]", gt)
        assert len(indices) == 2, "These fixtures contain diploid calls"
        call = tuple(sorted("." if index == "." else alleles[int(index)] for index in indices))
        calls.append(call)
    return Counter(calls)


@pytest.mark.integration
@pytest.mark.parametrize("shape", ["joint", "split"])
@pytest.mark.parametrize("with_gnomad", [False, True], ids=["without_gnomad", "with_gnomad"])
def test_plink2_roundtrip_preserves_carried_alleles(tmp_path: Path, shape: str, with_gnomad: bool):
    prefix = _export(tmp_path, shape, with_gnomad)
    out = tmp_path / "roundtrip"
    result = subprocess.run(
        [
            _plink2_executable(),
            "--bfile",
            str(prefix),
            "--export",
            "vcf",
            "--out",
            str(out),
            "--memory",
            "640",
            "--threads",
            "1",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert _decode_roundtrip_vcf(Path(str(out) + ".vcf")) == EXPECTED[shape]


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("shape", ["joint", "split", "placements"])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("cache_order", [False, True])
def test_coding_does_not_depend_on_record_or_cache_order(
    tmp_path, build, shape, reverse, cache_order
):
    """Each literal record survives, even when no reference pair fits its call."""
    data = tmp_path / "data"
    create_cache(data)
    coords = [
        ("2", 1000, "G", "A", "rs900961001", 0.001),
        ("2", 1000, "G", "T", "rs900961001", 0.001),
    ]
    if shape == "placements":
        coords[1] = ("2", 2000, "C", "T", "rs900961001", 0.001)
    with contextlib.closing(sqlite3.connect(data / "gnomad.sqlite")) as conn:
        conn.execute("DELETE FROM gnomad_frequencies")
        conn.executemany(
            "INSERT INTO gnomad_frequencies (chrom,pos,ref,alt,rsid,af) VALUES (?,?,?,?,?,?)",
            list(reversed(coords)) if cache_order else coords,
        )
        conn.commit()
    header = (FIXTURES / "plink_joint.vcf").read_text().split("#CHROM")[0]
    header = header.replace("242193529", "243199373") if build == "GRCh37" else header
    rows = {
        "joint": ["2\t1000\trs900961001\tG\tA,T\t99\tPASS\t.\tGT:GQ\t2|1:99"],
        "split": [
            "2\t1000\trs900961001\tG\tA\t99\tPASS\t.\tGT:GQ\t0/1:99",
            "2\t1000\trs900961001\tG\tT\t99\tPASS\t.\tGT:GQ\t1|0:99",
        ],
        "placements": [
            "2\t1000\trs900961001\tG\tA\t99\tPASS\t.\tGT:GQ\t0/1:99",
            "2\t2000\trs900961001\tC\tT\t99\tPASS\t.\tGT:GQ\t0/1:99",
        ],
    }[shape]
    if reverse:
        rows.reverse()
    src = tmp_path / "input.vcf"
    src.write_text(
        header
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        + "\n".join(rows)
        + "\n"
    )
    prefix = tmp_path / "out"
    result = CliRunner().invoke(
        main,
        [
            "export",
            "plink",
            str(src),
            "--data-dir",
            str(data),
            "--build",
            build.lower(),
            "--output-prefix",
            str(prefix),
        ],
    )
    assert result.exit_code == 0, result.output
    expected = (
        Counter({("A", "G"): 1, ("C", "T"): 1}) if shape == "placements" else EXPECTED[shape]
    )
    assert _decode_one_sample_bed(prefix) == expected
    assert len(Path(str(prefix) + ".bim").read_text().splitlines()) == len(rows)


@pytest.mark.parametrize("with_gnomad", [False, True])
def test_two_alt_call_at_named_grch38_locus(tmp_path, with_gnomad):
    """A synthetic call at the published multiallelic placement keeps A/T."""
    data = tmp_path / "data"
    data.mkdir()
    if with_gnomad:
        create_cache(data)
        with contextlib.closing(sqlite3.connect(data / "gnomad.sqlite")) as conn:
            conn.execute(
                "UPDATE gnomad_frequencies SET chrom='22',pos=20010383,rsid='rs933702824'"
            )
            conn.commit()
    src = tmp_path / "input.vcf"
    src.write_text(
        (FIXTURES / "plink_joint.vcf")
        .read_text()
        .replace("ID=2,length=242193529", "ID=22,length=50818468")
        .replace("2\t1000\trs900961001", "22\t20010383\trs933702824")
    )
    prefix = tmp_path / "out"
    result = CliRunner().invoke(
        main,
        ["export", "plink", str(src), "--data-dir", str(data), "--output-prefix", str(prefix)],
    )
    assert result.exit_code == 0, result.output
    assert _decode_one_sample_bed(prefix) == Counter({("A", "T"): 1})

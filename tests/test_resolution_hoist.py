# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""rsID-less resolution as one shared step for every command.

``stats``, ``extract`` and ``export plink`` used to compare ``v.rsid``
against rsID-keyed sets without ever reaching the resolver that lived
inside ``_pipeline._flush()``; build detection only matched seeds by
rsID. These tests run each path on rsID-less VCF input (``ID=.``).

All caches and genotypes are synthetic. The rsIDs and coordinates are the
public dbSNP placements of well-known SNPs; no genotype data from a real
file is used.
"""

from __future__ import annotations

import contextlib
import sqlite3
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner

from allelix.annotators._resolution import resolve_batch, resolve_rsidless
from allelix.annotators.gnomad import GnomadAnnotator
from allelix.cli import main
from allelix.databases._versions import GNOMAD_SCHEMA_VERSION
from allelix.databases.gnomad_loader import GNOMAD_DB_FILENAME
from allelix.databases.rsid_identity_loader import RSID_IDENTITY_DB_FILENAME, write_identity_db
from allelix.databases.schema import GNOMAD_SCHEMA
from allelix.models import Variant
from allelix.reports._pipeline import _BuildDetectionState
from allelix.utils.build_detect import (
    BUILD_GRCH37,
    BUILD_GRCH38,
    KNOWN_SNP_POSITIONS,
    detect_build,
    known_snp_seed,
)

if TYPE_CHECKING:
    from pathlib import Path

# APOE seeds: rs7412 C>T at GRCh38 19:44908822, rs429358 T>C at 19:44908684.
_RS7412_38 = ("19", 44_908_822, "C", "T")
_RS429358_38 = ("19", 44_908_684, "T", "C")


def _gnomad(data: Path, rows: list[tuple[str, int, str, str, str, float]]) -> None:
    with contextlib.closing(sqlite3.connect(data / GNOMAD_DB_FILENAME)) as conn:
        conn.executescript(GNOMAD_SCHEMA)
        conn.executemany(
            "INSERT INTO gnomad_frequencies (chrom, pos, ref, alt, rsid, af) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            rows,
        )
        conn.execute(
            "INSERT INTO database_versions (name, source_url, version, downloaded_at, "
            "record_count, local_version_tag) VALUES ('gnomad', 'test://hoist', '4.1', "
            "'2026-10-10', ?, ?)",
            (len(rows), f"sv:{GNOMAD_SCHEMA_VERSION}"),
        )
        conn.commit()


def _identity(data: Path, rows: list[tuple]) -> None:
    write_identity_db(
        rows,
        data / RSID_IDENTITY_DB_FILENAME,
        version="dbSNP test",
        source_url="test://hoist",
        downloaded_at="2026-10-10T00:00:00Z",
    )


def _vcf(path: Path, rows: list[str], assembly: str = "GRCh38", chrom_prefix: str = "") -> Path:
    header = (
        "##fileformat=VCFv4.2\n"
        f"##contig=<ID={chrom_prefix}19,length=58617616,assembly={assembly}>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    )
    path.write_text(header + "".join(r + "\n" for r in rows))
    return path


# ── The shared step ───────────────────────────────────────────────


class TestSharedResolutionStep:
    def test_stream_resolves_in_order_and_keeps_unresolvable(self, tmp_path: Path) -> None:
        data = tmp_path / "data"
        data.mkdir()
        _gnomad(data, [(*_RS7412_38, "rs7412", 0.08)])
        variants = [
            Variant("", "19", 44_908_822, "C", "T", build=BUILD_GRCH38, ref="C"),
            Variant("", "19", 100, "A", "G", build=BUILD_GRCH38, ref="A"),
            Variant("rs1", "19", 200, "A", "G", build=BUILD_GRCH38, ref="A"),
        ]
        with contextlib.closing(GnomadAnnotator(data)) as gnomad:
            out = list(resolve_rsidless(iter(variants), clinvar=None, gnomad=gnomad, batch_size=2))
        assert [v.rsid for v in out] == ["rs7412", "", "rs1"]
        assert out[0].original_rsid == ""

    def test_no_sources_is_a_no_op(self) -> None:
        v = Variant("", "19", 44_908_822, "C", "T", build=BUILD_GRCH38, ref="C")
        resolve_batch([v], clinvar=None, gnomad=None, identity=None)
        assert v.rsid == ""

    def test_vendor_probe_is_never_resolved(self, tmp_path: Path) -> None:
        """A call with no REF and no rsID stays unresolved (vendor-probe guard)."""
        data = tmp_path / "data"
        data.mkdir()
        _gnomad(data, [(*_RS7412_38, "rs7412", 0.08)])
        v = Variant("i5000123", "19", 44_908_822, "C", "T", build=BUILD_GRCH38)
        with contextlib.closing(GnomadAnnotator(data)) as gnomad:
            resolve_batch([v], clinvar=None, gnomad=gnomad)
        assert v.rsid == "i5000123"


# ── stats: high-value warnings on rsID-less VCF ───────────────────


class TestStatsHighValue:
    def test_rsidless_no_call_and_resolved_partner(self, tmp_path: Path) -> None:
        """rs429358 is a ./. record with ID=. and rs7412 a called record with
        ID=. The no-call is placed by the identity table; rs7412 is only
        known through gnomAD resolution, which makes the APOE cluster
        incomplete rather than silently absent."""
        data = tmp_path / "data"
        data.mkdir()
        _gnomad(data, [(*_RS7412_38, "rs7412", 0.08)])
        _identity(data, [("rs429358", "19", 45_411_941, "T", "C", 44_908_684, "T", "C", 0)])
        vcf = _vcf(
            tmp_path / "apoe.vcf",
            [
                "19\t44908684\t.\tT\tC\t50\tPASS\t.\tGT\t./.",
                "19\t44908822\t.\tC\tT\t50\tPASS\t.\tGT\t0/1",
            ],
        )
        result = CliRunner(env={"COLUMNS": "200"}).invoke(
            main, ["stats", str(vcf), "--data-dir", str(data)]
        )
        assert result.exit_code == 0, result.output
        assert "High-value no-calls" in result.output
        line = next(ln for ln in result.output.splitlines() if "rs429358" in ln)
        assert "APOE genotype cannot be determined" in line

    def test_without_databases_nothing_is_invented(self, tmp_path: Path) -> None:
        vcf = _vcf(tmp_path / "apoe.vcf", ["19\t44908684\t.\tT\tC\t50\tPASS\t.\tGT\t./."])
        result = CliRunner(env={"COLUMNS": "200"}).invoke(
            main, ["stats", str(vcf), "--data-dir", str(tmp_path / "empty")]
        )
        assert result.exit_code == 0, result.output
        assert "rs429358" not in result.output


# ── export plink: ref/alt for gnomAD-resolvable rsID-less calls ───


class TestPlinkRefAlt:
    def test_rsidless_hom_alt_gets_ref_alt_not_monomorphic(self, tmp_path: Path) -> None:
        data = tmp_path / "data"
        data.mkdir()
        _gnomad(data, [(*_RS7412_38, "rs7412", 0.08)])
        vcf = _vcf(tmp_path / "plink.vcf", ["19\t44908822\t.\tC\tT\t50\tPASS\t.\tGT\t1/1"])
        prefix = tmp_path / "out"
        result = CliRunner(env={"COLUMNS": "200"}).invoke(
            main,
            ["export", "plink", str(vcf), "-o", str(prefix), "--data-dir", str(data)],
        )
        assert result.exit_code == 0, result.output
        fields = prefix.with_suffix(".bim").read_text().strip().split("\t")
        assert fields[1] == "rs7412"
        assert (fields[4], fields[5]) == ("C", "T")
        assert prefix.with_suffix(".bed").read_bytes()[3] == 0b11
        assert "monomorphic" not in result.output

    def test_unresolvable_call_still_falls_back(self, tmp_path: Path) -> None:
        data = tmp_path / "data"
        data.mkdir()
        _gnomad(data, [(*_RS7412_38, "rs7412", 0.08)])
        vcf = _vcf(tmp_path / "plink.vcf", ["19\t100\t.\tA\tG\t50\tPASS\t.\tGT\t1/1"])
        prefix = tmp_path / "out"
        result = CliRunner(env={"COLUMNS": "200"}).invoke(
            main,
            ["export", "plink", str(vcf), "-o", str(prefix), "--data-dir", str(data)],
        )
        assert result.exit_code == 0, result.output
        assert prefix.with_suffix(".bim").read_text().split("\t")[5].strip() == "0"


# ── Build detection on rsID-less input ────────────────────────────


def _seed_records(build: str, n: int = 6) -> list[Variant]:
    """rsID-less records at the first ``n`` seeds' coordinates on ``build``."""
    return [
        Variant("", entry[build][0], entry[build][1], "A", "G", ref="A")
        for entry in list(KNOWN_SNP_POSITIONS.values())[:n]
    ]


class TestRsidlessBuildDetection:
    def test_seed_coordinates_are_unique_across_builds(self) -> None:
        loci = [locus for entry in KNOWN_SNP_POSITIONS.values() for locus in entry.values()]
        assert len(loci) == len(set(loci))

    @pytest.mark.parametrize("build", [BUILD_GRCH37, BUILD_GRCH38])
    def test_detect_build_votes_by_position(self, build: str) -> None:
        result = detect_build(_seed_records(build))
        assert result.build == build
        assert result.is_confident

    def test_rsid_bearing_record_does_not_vote_by_position(self) -> None:
        """A record that names a different rsID at a seed coordinate is a
        different variant; only rsID-less records vote by position."""
        rsid, entry = next(iter(KNOWN_SNP_POSITIONS.items()))
        chrom, pos = entry[BUILD_GRCH37]
        assert known_snp_seed(Variant("rs999", chrom, pos, "A", "G")) is None
        assert known_snp_seed(Variant("", chrom, pos, "A", "G")) == rsid
        assert known_snp_seed(Variant("chr1:5:A:G", chrom, pos, "A", "G")) == rsid

    def test_detection_state_converges_on_rsidless_grch37(self) -> None:
        """No header build and no chr prefix: the seeds alone decide."""
        diag = _BuildDetectionState(override=None, header_build=None)
        out: list[Variant] = []
        for v in _seed_records(BUILD_GRCH37):
            _ready, batch = diag.feed(v)
            out.extend(batch)
        out.extend(diag.flush())
        assert diag.effective_build == BUILD_GRCH37
        assert diag.position_confident
        assert {v.build for v in out} == {BUILD_GRCH37}

    def test_positions_beat_a_chr_prefix_on_ucsc_hg19(self, tmp_path: Path) -> None:
        """A GRCh37 gVCF with UCSC chr-prefixed contigs used to be called
        GRCh38 by the chr-prefix heuristic with no dissent."""
        from allelix.parsers.vcf import VcfParser
        from allelix.reports._pipeline import run_analysis

        rows = [
            f"chr{v.chromosome}\t{v.position}\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"
            for v in _seed_records(BUILD_GRCH37)
        ]
        header = (
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=249250621>\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        )
        vcf = tmp_path / "hg19.vcf"
        vcf.write_text(header + "\n".join(rows) + "\n")
        result = run_analysis(vcf, VcfParser(), [])
        assert result.build == BUILD_GRCH37
        assert result.build_diagnostics is not None
        assert result.build_diagnostics.position_confident


# ── extract via tabix: fetched rows go through the shared step ────


class _FakeTabix:
    """Stands in for pysam.TabixFile: serves fixed rows by contig."""

    def __init__(self, rows: dict[str, list[str]]) -> None:
        self._rows = rows

    def fetch(self, chrom: str, start: int, end: int) -> list[str]:
        if chrom not in self._rows:
            raise ValueError(chrom)
        return [r for r in self._rows[chrom] if start < int(r.split("\t")[1]) <= end]

    def close(self) -> None:
        pass


class TestTabixExtract:
    def test_rsidless_row_resolves_to_the_wanted_rsid(self, tmp_path: Path) -> None:
        from types import SimpleNamespace

        from allelix.cli.utility import _execute_tabix_extract

        data = tmp_path / "data"
        data.mkdir()
        _gnomad(data, [(*_RS7412_38, "rs7412", 0.08)])
        row = "chr19\t44908822\t.\tC\tT\t50\tPASS\t.\tGT\t0/1"
        vcf = _vcf(tmp_path / "x.vcf", [row], chrom_prefix="chr")
        pysam = SimpleNamespace(TabixFile=lambda _path: _FakeTabix({"chr19": [row]}))
        found = _execute_tabix_extract(vcf, {"rs7412"}, None, data, pysam)
        assert found is not None
        assert found["rs7412"].genotype in {"C/T", "CT"}
        assert found["rs7412"].original_rsid == ""

    def test_a_different_variant_at_the_coordinate_is_not_taken(self, tmp_path: Path) -> None:
        from types import SimpleNamespace

        from allelix.cli.utility import _execute_tabix_extract

        data = tmp_path / "data"
        data.mkdir()
        _gnomad(data, [(*_RS7412_38, "rs7412", 0.08)])
        row = "19\t44908822\t.\tC\tA\t50\tPASS\t.\tGT\t0/1"  # C>A, not rs7412's C>T
        vcf = _vcf(tmp_path / "x.vcf", [row])
        pysam = SimpleNamespace(TabixFile=lambda _path: _FakeTabix({"19": [row]}))
        assert _execute_tabix_extract(vcf, {"rs7412"}, None, data, pysam) == {}

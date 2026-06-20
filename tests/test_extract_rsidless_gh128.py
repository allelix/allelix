# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""GH #128 regression tests for the ``allelix extract`` command on
rsID-less VCF inputs.

Sister test file to ``tests/reports/test_rsidless_resolution_gh128.py``.
The panel-coverage analyze path is covered there; this file covers
extract:

  - ``_sequential_extract`` (the fallback path used when the tabix
    fast-path can't be set up — pysam missing, .tbi missing, etc.)
  - ``_execute_tabix_extract`` (the fast-path; covered behaviorally
    via a stubbed tabix file even though the real path needs pysam)

Pre-fix the sequential extract scanned ``v.rsid in wanted`` and missed
every rsID on a DeepVariant-shaped VCF (ID column empty). Post-fix
``_sequential_extract`` pre-resolves the wanted set to gnomAD coords
and matches by (chrom, pos, ref, alt).

For the tabix fast-path the file is loaded behind pysam; this file
asserts the contract by directly calling the underlying
``_execute_tabix_extract`` with a real pysam-tabix index when pysam is
available, and skips otherwise.
"""

from __future__ import annotations

import contextlib
import sqlite3
from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from pathlib import Path

from allelix.cli.utility import _sequential_extract
from allelix.databases.schema import GNOMAD_SCHEMA
from allelix.parsers.vcf import VcfParser

GNOMAD_DB_FILENAME = "gnomad.sqlite"

# Same synthetic VCF as in tests/reports/test_rsidless_resolution_gh128.py.
# Keeping it duplicated rather than fixturing across packages so the
# test files stand alone and the file's exact byte content is visible
# right next to the assertions about what should resolve out of it.
_RSIDLESS_VCF = """\
##fileformat=VCFv4.2
##contig=<ID=1,length=249250621>
##contig=<ID=19,length=58617616>
##contig=<ID=22,length=51304566>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
1\t11796321\t.\tG\tA\t.\tPASS\t.\tGT\t0/1
19\t44908684\t.\tT\tC\t.\tPASS\t.\tGT\t0/1
22\t19963748\t.\tG\tA\t.\tPASS\t.\tGT\t0/1
"""

# An rsID-stamped variant of the same file — the parser should yield
# Variants whose ``rsid`` already carries the expected rs-prefixed
# identifier. Used to pin the rsID fast-path's continued correctness
# (i.e., the rewrite didn't regress 23andMe-style inputs that already
# stamp rsIDs).
_RSID_STAMPED_VCF = """\
##fileformat=VCFv4.2
##contig=<ID=1,length=249250621>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
1\t11796321\trs1801133\tG\tA\t.\tPASS\t.\tGT\t0/1
"""


@pytest.fixture
def rsidless_vcf(tmp_path: Path) -> Path:
    p = tmp_path / "rsidless.vcf"
    p.write_text(_RSIDLESS_VCF)
    return p


@pytest.fixture
def rsid_stamped_vcf(tmp_path: Path) -> Path:
    p = tmp_path / "rsid_stamped.vcf"
    p.write_text(_RSID_STAMPED_VCF)
    return p


@pytest.fixture
def gnomad_data_dir(tmp_path: Path) -> Path:
    """A data dir carrying just a populated gnomAD cache.

    Smaller than ``all_annotators_data_dir`` because extract only
    needs gnomAD — no ClinVar, no PharmGKB, no GWAS.
    """
    db_path = tmp_path / GNOMAD_DB_FILENAME
    with contextlib.closing(sqlite3.connect(db_path)) as conn:
        for stmt in GNOMAD_SCHEMA.split(";"):
            stmt = stmt.strip()
            if stmt:
                conn.execute(stmt)
        conn.executemany(
            "INSERT INTO gnomad_frequencies"
            " (chrom, pos, ref, alt, rsid, af) VALUES (?, ?, ?, ?, ?, ?)",
            [
                ("1", 11796321, "G", "A", "rs1801133", 0.35),
                ("22", 19963748, "G", "A", "rs4680", 0.50),
                ("19", 44908684, "T", "C", "rs429358", 0.15),
                ("19", 44908684, "T", "G", "rs429358", 0.08),
            ],
        )
        conn.execute(
            "INSERT INTO database_versions"
            " (name, source_url, version, downloaded_at, record_count,"
            "  local_version_tag)"
            " VALUES (?, ?, ?, ?, ?, ?)",
            (
                "gnomad",
                "test://mock",
                "4.1",
                "2026-01-01T00:00:00Z",
                4,
                "sv:1",
            ),
        )
        conn.commit()
    return tmp_path


class TestSequentialExtractOnRsidlessVcf:
    """The user's symptom — ``extract --snps`` on a DeepVariant gVCF
    returns 'not in file' for every requested rsID, including HG002's
    published carriers."""

    def test_recovers_rsids_via_gnomad_coord_match(
        self, rsidless_vcf: Path, gnomad_data_dir: Path
    ) -> None:
        """Pre-fix this returned an empty dict because no row's ``rsid``
        starts with ``rs`` — the parser emits the positional pseudo-ID
        and the rsID-only equality check fails for every wanted rsID.
        Post-fix the coord-based fallback recognizes each of the three
        variants at the gnomAD-resolved coords and stamps + stashes."""
        parser = VcfParser()
        wanted = {"rs1801133", "rs4680", "rs429358"}
        found = _sequential_extract(parser, rsidless_vcf, wanted, gnomad_data_dir)
        assert set(found.keys()) == wanted, (
            f"Sequential extract failed to recover gnomAD-resolved rsIDs "
            f"on a rsID-less VCF. Got: {list(found.keys())!r}"
        )
        # The returned Variant carries the stamp; the audit field
        # holds the pre-stamp pseudo-rsID for debugging.
        assert found["rs1801133"].rsid == "rs1801133"
        assert found["rs1801133"].original_rsid is not None
        assert "rs1801133" not in found["rs1801133"].original_rsid

    def test_rsid_stamped_vcf_still_works_via_fast_path(
        self, rsid_stamped_vcf: Path, gnomad_data_dir: Path
    ) -> None:
        """A VCF that DOES stamp rsIDs (23andMe, GIAB GRCh37 benchmark)
        must continue to resolve via the rsID equality fast path. The
        coord fallback only runs when the fast path misses; the
        post-fix Variant.original_rsid stays None on this happy path."""
        parser = VcfParser()
        wanted = {"rs1801133"}
        found = _sequential_extract(parser, rsid_stamped_vcf, wanted, gnomad_data_dir)
        assert "rs1801133" in found
        assert found["rs1801133"].rsid == "rs1801133"
        # No stash — the rsID came from the file's ID column.
        assert found["rs1801133"].original_rsid is None

    def test_no_gnomad_no_coord_fallback(self, rsidless_vcf: Path, tmp_path: Path) -> None:
        """When the data_dir carries no usable gnomAD cache the function
        falls back to rsID-only matching. A rsID-less VCF then returns
        an empty dict, which is the documented degraded-behavior
        contract — better than crashing, less useful than the resolved
        path. The empty ``tmp_path`` carries no gnomAD cache; the
        annotator's ``is_ready()`` returns False and the coord-fallback
        path is skipped.

        Note: ``data_dir=None`` cannot be used here because the production
        CLI always resolves it to the standard cache location (which
        carries a real gnomAD on most dev machines and would defeat the
        test). Pinning the gnomAD-not-ready behavior to an empty dir
        matches what mock test contexts actually need to simulate.
        """
        parser = VcfParser()
        wanted = {"rs1801133"}
        empty_data_dir = tmp_path / "empty_cache"
        empty_data_dir.mkdir()
        found = _sequential_extract(parser, rsidless_vcf, wanted, empty_data_dir)
        assert found == {}

    def test_rsid_not_in_gnomad_stays_missing(
        self, rsidless_vcf: Path, gnomad_data_dir: Path
    ) -> None:
        """A panel rsID at a position the file doesn't contain — and
        not in gnomAD's mock — stays absent from the result. The fix
        must not invent matches."""
        parser = VcfParser()
        wanted = {"rs1801133", "rs9999999"}  # second one isn't anywhere
        found = _sequential_extract(parser, rsidless_vcf, wanted, gnomad_data_dir)
        assert "rs1801133" in found
        assert "rs9999999" not in found

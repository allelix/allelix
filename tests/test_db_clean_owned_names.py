# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""`db clean` removes only the names allelix creates in the data directory.

Caches, their SQLite side files, the files kept beside them, staging entries
(current and earlier releases' names) and earlier releases' downloads are
removed. Everything else (CADD full-mode files, a genome file, a reports
folder, config.toml) is kept.
"""

from __future__ import annotations

import contextlib
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner

from allelix.cli import main
from allelix.databases.staging import owned_directory, owned_staging

if TYPE_CHECKING:
    from pathlib import Path

_CADD_FULL_MODE = (
    "whole_genome_SNVs.tsv.gz",
    "whole_genome_SNVs.tsv.gz.tbi",
    "gnomad.genomes.r4.0.indel.tsv.gz",
    "gnomad.genomes.r4.0.indel.tsv.gz.tbi",
)

_CACHES = (
    "clinvar.GRCh37.sqlite",
    "clinvar.GRCh38.sqlite",
    "gnomad.sqlite",
    "alphamissense.sqlite",
    "cadd.sqlite",
    "snpedia.sqlite",
    "gwas.sqlite",
    "pharmgkb.sqlite",
    "rsid_identity.sqlite",
)


def _clean(data: Path, *args: str):
    return CliRunner().invoke(main, ["db", "clean", "--data-dir", str(data), *args])


def _flat(text: str) -> str:
    return " ".join(text.split())


def test_keep_cadd_keeps_cadd_full_mode_files(tmp_path):
    """README: full mode needs whole_genome_SNVs.tsv.gz (~81 GB), the indel file
    and their .tbi in the data directory; `db update` cannot fetch them."""
    for name in (*_CADD_FULL_MODE, "cadd.sqlite", "clinvar.GRCh37.sqlite"):
        (tmp_path / name).write_bytes(b"x")

    result = _clean(tmp_path, "--keep-cadd", "--yes")

    assert result.exit_code == 0, result.output
    assert not (tmp_path / "clinvar.GRCh37.sqlite").exists()
    assert (tmp_path / "cadd.sqlite").exists()
    lost = [name for name in _CADD_FULL_MODE if not (tmp_path / name).exists()]
    assert not lost, f"--keep-cadd deleted CADD full mode's files: {lost}"


def test_plain_clean_keeps_cadd_full_mode_files(tmp_path):
    for name in (*_CADD_FULL_MODE, "cadd.sqlite"):
        (tmp_path / name).write_bytes(b"x")

    result = _clean(tmp_path, "--yes")

    assert result.exit_code == 0, result.output
    assert not (tmp_path / "cadd.sqlite").exists()
    assert all((tmp_path / name).exists() for name in _CADD_FULL_MODE)
    assert "whole_genome_SNVs.tsv.gz" in _flat(result.output)


def test_files_allelix_did_not_create_are_kept(tmp_path):
    (tmp_path / "gnomad.sqlite").write_bytes(b"x")
    (tmp_path / "genome.txt").write_text("rs1\t1\t1\tAA\n")
    (tmp_path / "gnomad_notes.txt").write_text("mine")
    (tmp_path / "clinvar.GRCh38.sqlite.bak").write_bytes(b"x")
    reports = tmp_path / "reports"
    reports.mkdir()
    (reports / "report.html").write_text("<html></html>")
    hidden = tmp_path / ".gnomad-mine"
    hidden.mkdir()
    (tmp_path / "config.toml").write_text("[sources]\n")

    result = _clean(tmp_path, "--yes")

    assert result.exit_code == 0, result.output
    assert not (tmp_path / "gnomad.sqlite").exists()
    for kept in (
        "genome.txt",
        "gnomad_notes.txt",
        "clinvar.GRCh38.sqlite.bak",
        "reports/report.html",
        ".gnomad-mine",
        "config.toml",
    ):
        assert (tmp_path / kept).exists(), kept
    assert "Kept 5 item(s)" in _flat(result.output)


def test_unrelated_directory_is_left_alone(tmp_path):
    (tmp_path / "resume.pdf").write_bytes(b"x")
    (tmp_path / "config.toml").write_text("[package]\n")

    for args in (("--yes",), ("--dry-run",), ("--force", "--yes")):
        result = _clean(tmp_path, *args)
        assert result.exit_code == 0, result.output
        assert "Nothing to clean" in _flat(result.output)
    assert (tmp_path / "resume.pdf").exists()
    assert (tmp_path / "config.toml").exists()


@pytest.mark.parametrize("cache", _CACHES)
def test_every_cache_and_its_side_files_are_removed(tmp_path, cache):
    names = [cache, *(f"{cache}{side}" for side in ("-wal", "-shm", "-journal"))]
    for name in names:
        (tmp_path / name).write_bytes(b"x")

    result = _clean(tmp_path, "--yes")

    assert result.exit_code == 0, result.output
    assert [name for name in names if (tmp_path / name).exists()] == []


def test_files_kept_beside_caches_are_removed(tmp_path):
    for name in ("clinvar-build.log", "gwas_catalog_associations.tsv"):
        (tmp_path / name).write_text("x")

    result = _clean(tmp_path, "--yes")

    assert result.exit_code == 0, result.output
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize(
    "name",
    [
        # Downloads earlier releases wrote into the data directory.
        "gnomad.sqlite.gz",
        "alphamissense.sqlite.gz",
        "cadd.sqlite.gz",
        "snpedia.sqlite.gz",
        "rsid_identity.sqlite.gz",
        "gwas_catalog_associations.zip",
        "clinicalAnnotations.zip",
        "clinvar.vcf.gz",
        "clinvar.vcf",
        "gnomad.sqlite.gz.part",
        # Staging names shared by every writer, before 2.6.0.
        "gwas.sqlite.tmp",
        "gwas.sqlite.tmp-journal",
        "clinvar.GRCh38.sqlite.partial",
    ],
)
def test_earlier_releases_names_are_removed(tmp_path, name):
    (tmp_path / name).write_bytes(b"x")

    result = _clean(tmp_path, "--yes")

    assert result.exit_code == 0, result.output
    assert not (tmp_path / name).exists()


def test_staging_entries_are_removed_without_force(tmp_path):
    """Entries a killed update leaves (hidden staging files and directories) are
    allelix's, so a directory holding only those is cleaned without --force."""
    stack = contextlib.ExitStack()
    made: list[Path] = []
    with stack:
        for cache in _CACHES:
            for suffix in (".tmp", ".part", ".partial"):
                made.append(stack.enter_context(owned_staging(tmp_path / cache, suffix)))
        for prefix in (
            "gnomad",
            "alphamissense",
            "cadd",
            "snpedia",
            "gwas",
            "pharmgkb",
            "rsid_identity",
            "clinvar",
        ):
            directory = stack.enter_context(owned_directory(tmp_path, prefix))
            # What a killed download leaves there (2.6.1: only downloads qualify).
            (directory / "variant_summary.txt.gz").write_bytes(b"x")
            made.append(directory)
        journal = made[0].with_name(made[0].name + "-journal")
        journal.write_bytes(b"x")
        made.append(journal)

        result = _clean(tmp_path, "--yes")

        assert result.exit_code == 0, result.output
        assert [path.name for path in made if path.exists()] == []


def test_keep_cadd_keeps_every_cadd_cache_name(tmp_path):
    with owned_staging(tmp_path / "cadd.sqlite") as staging:
        for name in ("cadd.sqlite", "cadd.sqlite-wal", "cadd.sqlite.gz"):
            (tmp_path / name).write_bytes(b"x")
        (tmp_path / "gnomad.sqlite").write_bytes(b"x")

        result = _clean(tmp_path, "--keep-cadd", "--yes")

        assert result.exit_code == 0, result.output
        assert not (tmp_path / "gnomad.sqlite").exists()
        assert staging.exists()
        for name in ("cadd.sqlite", "cadd.sqlite-wal", "cadd.sqlite.gz"):
            assert (tmp_path / name).exists(), name


def test_help_names_what_is_kept_and_what_cannot_be_fetched_again():
    result = CliRunner().invoke(main, ["db", "clean", "--help"])
    text = _flat(result.output)
    assert "every cache removed here can be re-fetched" not in text
    assert "allelix did not create" in text
    assert "CADD full-mode" in text


def test_completion_message_does_not_promise_everything_comes_back(tmp_path):
    (tmp_path / "gnomad.sqlite").write_bytes(b"x")
    (tmp_path / "clinvar-build.log").write_text("x")

    result = _clean(tmp_path, "--yes")

    text = _flat(result.output)
    assert "Re-populate with allelix db update" not in text
    assert "build log" in text


def test_a_staging_shaped_directory_holding_a_built_cache_is_kept(tmp_path):
    """Only a staging directory holding names allelix writes there (downloads
    and partial files) is removed. One holding a cache under its own name,
    built there or copied there as a backup, or an extraction, is kept and
    listed (2.6.1)."""
    (tmp_path / "gnomad.sqlite").write_bytes(b"x")
    kept = tmp_path / ".clinvar-abc12345"
    kept.mkdir()
    (kept / "variant_summary.txt.gz").write_bytes(b"x")
    (kept / "clinvar.GRCh38.sqlite").write_bytes(b"a partly built cache")
    extraction = tmp_path / ".gwas-abc98765"
    (extraction / "extracted").mkdir(parents=True)
    (extraction / "extracted" / "gwas-catalog-associations.tsv").write_bytes(b"x")
    removed = tmp_path / ".gwas-xyz98765"
    removed.mkdir()
    (removed / "gwas_catalog_associations.zip").write_bytes(b"x")
    (removed / ".gwas_catalog_associations.zip.0123456789abcdef.part").write_bytes(b"x")

    result = _clean(tmp_path, "--yes")

    text = _flat(result.output)
    assert result.exit_code == 0, result.output
    assert (kept / "clinvar.GRCh38.sqlite").exists()
    # An extraction's name comes from EBI's zip and cannot be proven allelix's.
    assert (extraction / "extracted" / "gwas-catalog-associations.tsv").exists()
    for name in (".clinvar-abc12345", ".gwas-abc98765"):
        assert f"Kept {name}: may be an interrupted update's leftovers" in text
    assert not removed.exists()


def test_a_staging_directory_holding_a_cache_partial_or_an_empty_extraction_is_removed(tmp_path):
    """A cache's partial file (with its side file) and an emptied ``extracted/``
    are names allelix writes in a staging directory: a build or a load killed
    there leaves only those, and ``db clean`` removes them. A partial file of a
    name allelix does not write keeps its directory (2.6.1)."""
    building = tmp_path / ".clinvar-abc12345"
    building.mkdir()
    (building / "variant_summary.txt.gz").write_bytes(b"x")
    (building / ".clinvar.GRCh38.sqlite.0123456789abcdef.partial").write_bytes(b"x")
    (building / ".clinvar.GRCh38.sqlite.0123456789abcdef.partial-journal").write_bytes(b"x")
    loading = tmp_path / ".gwas-abc98765"
    (loading / "extracted").mkdir(parents=True)
    (loading / "gwas_catalog_associations.zip").write_bytes(b"x")
    (loading / "gwas_catalog_associations.tsv").write_bytes(b"x")
    other = tmp_path / ".clinvar-xyz98765"
    other.mkdir()
    (other / ".family.sqlite.0123456789abcdef.partial").write_bytes(b"x")

    result = _clean(tmp_path, "--yes")

    assert result.exit_code == 0, result.output
    assert not building.exists()
    assert not loading.exists()
    assert (other / ".family.sqlite.0123456789abcdef.partial").exists()
    assert "Kept .clinvar-xyz98765: may be an interrupted update's leftovers" in _flat(
        result.output
    )

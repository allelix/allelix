# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Caches are read-only at analyze time (ADR-0048).

``analyze --no-update`` writes nothing to the data directory; each
source it cannot use is named with its state and the fix; the installer
leaves caches in rollback-journal mode. All rows are synthetic.
"""

from __future__ import annotations

import base64
import contextlib
import gzip
import hashlib
import os
import shutil
import sqlite3
import stat
import subprocess
from pathlib import Path

import pytest
from click.testing import CliRunner

from allelix.annotators.clinvar import ClinVarAnnotator, clinvar_db_filename
from allelix.cli import main
from allelix.databases._versions import (
    ALPHAMISSENSE_SCHEMA_VERSION,
    CADD_SCHEMA_VERSION,
    GNOMAD_SCHEMA_VERSION,
    PHARMGKB_INTERPRETER_VERSION,
    RSID_IDENTITY_SCHEMA_VERSION,
)
from allelix.databases.cache_state import CacheState, is_wal
from allelix.databases.loader_utils import install_prebuilt_gz_cache
from allelix.databases.schema import (
    ALPHAMISSENSE_SCHEMA,
    CADD_SCHEMA,
    GNOMAD_SCHEMA,
    RSID_IDENTITY_SCHEMA,
)
from allelix.databases.snpedia_parser import _PARSER_VERSION

pytestmark = pytest.mark.skipif(
    hasattr(os, "geteuid") and os.geteuid() == 0, reason="root ignores file permissions"
)

_VERSIONS = (
    "CREATE TABLE IF NOT EXISTS database_versions (name TEXT PRIMARY KEY, "
    "source_url TEXT NOT NULL, version TEXT, downloaded_at TEXT NOT NULL, "
    "record_count INTEGER NOT NULL, remote_signal TEXT, local_version_tag TEXT);"
)
_SNPEDIA_SCHEMA = (
    "CREATE TABLE snpedia_genotypes (rsid TEXT NOT NULL, allele1 TEXT NOT NULL, "
    "allele2 TEXT NOT NULL, magnitude REAL, repute TEXT, summary TEXT, gene TEXT, "
    "orientation TEXT, scraped_at TEXT);"
)


def _cache(
    path: Path,
    schema: str,
    name: str,
    tag: str,
    insert: tuple[str, tuple[object, ...]],
    *,
    wal: bool = False,
) -> Path:
    """One-row cache with its version row; WAL journal mode when ``wal``."""
    with contextlib.closing(sqlite3.connect(path)) as conn:
        if wal:
            conn.execute("PRAGMA journal_mode=WAL")
        conn.executescript(schema + _VERSIONS)
        conn.execute(insert[0], insert[1])
        conn.execute(
            "INSERT OR REPLACE INTO database_versions VALUES (?, 'test://', 'v1', "
            "'2026-10-10', 1, NULL, ?)",
            (name, tag),
        )
        conn.commit()
    assert is_wal(path) is wal
    return path


def _gnomad(path: Path, *, wal: bool = False) -> Path:
    return _cache(
        path,
        GNOMAD_SCHEMA,
        "gnomad",
        f"sv:{GNOMAD_SCHEMA_VERSION}",
        (
            "INSERT INTO gnomad_frequencies (chrom, pos, ref, alt, rsid, af) "
            "VALUES ('1', 100, 'A', 'G', 'rs900780001', 0.1)",
            (),
        ),
        wal=wal,
    )


def _alphamissense(path: Path, *, wal: bool = False) -> Path:
    return _cache(
        path,
        ALPHAMISSENSE_SCHEMA,
        "alphamissense",
        f"sv:{ALPHAMISSENSE_SCHEMA_VERSION}",
        (
            "INSERT INTO alphamissense_scores (chrom, pos, ref, alt, rsid, uniprot_id, "
            "transcript_id, protein_variant, am_pathogenicity, am_class) "
            "VALUES ('1', 100, 'A', 'G', 'rs900780001', 'P0', 'T0', 'A1G', 0.9, "
            "'likely_pathogenic')",
            (),
        ),
        wal=wal,
    )


def _data_dir(root: Path, build_clinvar, mock_pharmgkb_dir: Path, mock_gwas_tsv: Path) -> Path:
    """Every cache type: ClinVar, ClinPGx, GWAS, SNPedia, gnomAD installed through
    the installer from a WAL file, AlphaMissense placed directly in WAL mode,
    CADD and the rsID identity table."""
    from allelix.databases.gwas_loader import load_gwas_tsv
    from allelix.databases.pharmgkb_loader import load_pharmgkb_tsv

    data = root / "data"
    data.mkdir()
    for build in ("GRCh37", "GRCh38"):
        build_clinvar(data / clinvar_db_filename(build), build)
    load_pharmgkb_tsv(mock_pharmgkb_dir, data / "pharmgkb.sqlite", source_url="test://")
    load_gwas_tsv(mock_gwas_tsv, data / "gwas.sqlite", source_url="test://")
    _cache(
        data / "snpedia.sqlite",
        _SNPEDIA_SCHEMA,
        "snpedia",
        f"pv:{_PARSER_VERSION}",
        (
            "INSERT INTO snpedia_genotypes (rsid, allele1, allele2, magnitude) "
            "VALUES ('rs900780009', 'A', 'G', 2.0)",
            (),
        ),
    )
    staged = _gnomad(root / "staged_gnomad.sqlite", wal=True)
    gz = root / "gnomad.sqlite.gz"
    with staged.open("rb") as f_in, gzip.open(gz, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    install_prebuilt_gz_cache(
        gz, data / "gnomad.sqlite", "gnomad", schema_version_tag=f"sv:{GNOMAD_SCHEMA_VERSION}"
    )
    _alphamissense(data / "alphamissense.sqlite", wal=True)
    _cache(
        data / "cadd.sqlite",
        CADD_SCHEMA,
        "cadd",
        f"sv:{CADD_SCHEMA_VERSION}",
        (
            "INSERT INTO cadd_scores (chrom, pos, ref, alt, phred) "
            "VALUES ('1', 100, 'A', 'G', 20.0)",
            (),
        ),
    )
    _cache(
        data / "rsid_identity.sqlite",
        RSID_IDENTITY_SCHEMA,
        "rsid_identity",
        f"sv:{RSID_IDENTITY_SCHEMA_VERSION}",
        (
            "INSERT INTO rsid_identity (rsid, chrom, pos37, ref37, alt37, pos38, ref38, alt38, "
            "palindromic) VALUES ('rs900780001', '1', 100, 'A', 'G', 200, 'A', 'G', 0)",
            (),
        ),
    )
    return data


def _snapshot(directory: Path) -> dict[str, tuple[int, int, str]]:
    return {
        p.name: (
            p.stat().st_size,
            p.stat().st_mtime_ns,
            hashlib.sha256(p.read_bytes()).hexdigest(),
        )
        for p in sorted(directory.iterdir())
    }


def _windows_acl(path: Path, restore: str | None = None) -> str:
    """Read or restore the fixture's DACL, retaining owner/group and inheritance."""
    shell = shutil.which("powershell") or shutil.which("pwsh")
    assert shell is not None, "Windows ACL fixtures require the installed PowerShell"
    script = (
        "$ErrorActionPreference='Stop'; $ProgressPreference='SilentlyContinue';"
        "$acl=Get-Acl -LiteralPath $env:ALLELIX_TEST_ACL_PATH;"
    )
    if restore is not None:
        script += (
            "$acl.SetSecurityDescriptorSddlForm($env:ALLELIX_TEST_ACL_SDDL,"
            "[System.Security.AccessControl.AccessControlSections]::Access);"
            "Set-Acl -LiteralPath $env:ALLELIX_TEST_ACL_PATH -AclObject $acl;"
            "$acl=Get-Acl -LiteralPath $env:ALLELIX_TEST_ACL_PATH;"
        )
    script += "[Console]::Write($acl.Sddl)"
    environment = dict(
        os.environ, ALLELIX_TEST_ACL_PATH=str(path), ALLELIX_TEST_ACL_SDDL=restore or ""
    )
    result = subprocess.run(
        [
            shell,
            "-NoLogo",
            "-NoProfile",
            "-NonInteractive",
            "-EncodedCommand",
            base64.b64encode(script.encode("utf-16-le")).decode("ascii"),
        ],
        env=environment,
        capture_output=True,
        text=True,
        encoding="utf-8",
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.strip(), "ACL snapshot is empty"
    return result.stdout.strip()


def _restored_acl_matches(original: str, restored: str) -> bool:
    """Require identical SDDL except an OS-added DACL AUTO_INHERITED marker.

    SetSecurityInfo can add SE_DACL_AUTO_INHERITED while converting the
    descriptor to the current inheritance model without changing semantics:
    https://learn.microsoft.com/en-us/windows/win32/secauthz/automatic-propagation-of-inheritable-aces
    Preserve every ACE (including order/flags), owner, group and other control
    marker. Only adding DACL AI is allowed; removing it or changing AR/P is not.
    """
    if restored == original:
        return True
    prefix, dacl_marker, dacl = restored.partition("D:")
    flags, ace_marker, entries = dacl.partition("(")
    if not dacl_marker or "(" in prefix or not ace_marker or flags.count("AI") != 1:
        return False
    previous_flags = flags.replace("AI", "", 1)
    if previous_flags not in ("", "P", "AR", "PAR"):
        return False
    return original == f"{prefix}D:{previous_flags}({entries}"


@contextlib.contextmanager
def _denied_permissions(path: Path, *, read: bool = False):
    """Prove denial and restore original access rules, ownership and inheritance."""
    if hasattr(os, "geteuid") and os.geteuid() == 0:
        pytest.skip("root ignores file permissions")
    probe = path / ".permission-denial-probe" if not read else None

    def allowed():
        if probe is not None:
            probe.write_bytes(b"synthetic permission probe")
            probe.unlink()
        else:
            with path.open("rb") as fh:
                fh.read(1)

    allowed()
    mode = path.stat().st_mode
    acl = _windows_acl(path) if os.name == "nt" else None
    try:
        if acl is not None:
            rights = "(RD)" if read else "(WD,AD,WEA,WA,DE,DC)"
            result = subprocess.run(
                ["icacls", str(path), "/deny", f"*S-1-1-0:{rights}"],
                capture_output=True,
                text=True,
                check=False,
            )
            assert result.returncode == 0, result.stdout + result.stderr
        else:
            path.chmod(mode & ~(0o444 if read else 0o222))
        if probe is None:
            with pytest.raises(PermissionError):
                allowed()
        else:
            with pytest.raises(PermissionError):
                probe.write_bytes(b"synthetic permission probe")
            assert not probe.exists(), "Directory denial must refuse file creation"
        yield
    finally:
        if acl is not None:
            restored = _windows_acl(path, acl)
            assert _restored_acl_matches(acl, restored), "Original ACL was not restored"
        else:
            path.chmod(mode)
        if probe is not None:
            probe.unlink(missing_ok=True)
        allowed()


def _read_only(directory: Path):
    """A real directory-create denial, on POSIX and Windows alike."""
    return _denied_permissions(directory)


def _flat(text: str) -> str:
    return "".join(text.split())


@pytest.fixture
def connects(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """Every ``sqlite3.connect`` target from here on (clear it after setup)."""
    calls: list[str] = []
    real = sqlite3.connect

    def _record(database, *args, **kwargs):  # type: ignore[no-untyped-def]
        calls.append(str(database))
        return real(database, *args, **kwargs)

    monkeypatch.setattr(sqlite3, "connect", _record)
    return calls


def _analyze(data: Path, out: Path, *extra: str) -> tuple[int, str]:
    result = CliRunner().invoke(
        main,
        [
            "analyze",
            str(Path(__file__).parent / "fixtures" / "mock_23andme.txt"),
            "--build",
            "grch37",
            "--no-update",
            "--data-dir",
            str(data),
            "--output",
            str(out),
            *extra,
        ],
    )
    if result.exception is not None and not isinstance(result.exception, SystemExit):
        raise result.exception
    # Rich wraps long paths anywhere, so compare without whitespace.
    return result.exit_code, "".join(result.output.split())


class TestNoUpdateWritesNothing:
    def test_read_only_data_directory(
        self,
        tmp_path: Path,
        build_synthetic_clinvar_cache,
        mock_pharmgkb_dir: Path,
        mock_gwas_tsv: Path,
    ) -> None:
        data = _data_dir(tmp_path, build_synthetic_clinvar_cache, mock_pharmgkb_dir, mock_gwas_tsv)
        assert not is_wal(data / "gnomad.sqlite")
        before = _snapshot(data)
        with _read_only(data):
            code, output = _analyze(data, tmp_path / "report.json")
        assert code == 0, output
        assert _snapshot(data) == before
        assert not [
            p.name for p in data.iterdir() if p.name.endswith(("-wal", "-shm", "-journal"))
        ]
        path = data / "alphamissense.sqlite"
        assert (
            _flat(
                f"Skipping AlphaMissense: {path} is in WAL mode, which cannot be read in a data "
                f"directory you cannot write. Run `allelix db update` once with write access "
                f"to {data}."
            )
            in output
        )
        assert _flat("Skipping gnomAD") not in output

    def test_writable_data_directory_skips_a_wal_cache(
        self,
        tmp_path: Path,
        build_synthetic_clinvar_cache,
        mock_pharmgkb_dir: Path,
        mock_gwas_tsv: Path,
    ) -> None:
        data = _data_dir(tmp_path, build_synthetic_clinvar_cache, mock_pharmgkb_dir, mock_gwas_tsv)
        before = _snapshot(data)
        code, output = _analyze(data, tmp_path / "report.json")
        assert code == 0, output
        assert _snapshot(data) == before
        path = data / "alphamissense.sqlite"
        assert (
            _flat(
                f"Skipping AlphaMissense: {path} is in WAL mode from an older Allelix. "
                "Run `allelix db update` once to convert it (converting it downloads nothing)."
            )
            in output
        )

    def test_version_reads_never_open_a_wal_cache(
        self, tmp_path: Path, connects: list[str]
    ) -> None:
        """AlphaMissense's check of gnomAD's version treats a WAL gnomAD as unavailable."""
        from allelix.annotators.alphamissense import AlphaMissenseAnnotator
        from allelix.databases.manager import get_database_info

        _gnomad(tmp_path / "gnomad.sqlite", wal=True)
        path = _alphamissense(tmp_path / "alphamissense.sqlite")
        with contextlib.closing(sqlite3.connect(path)) as conn:
            conn.execute(
                "INSERT INTO database_versions VALUES ('alphamissense_gnomad_source', "
                "'test://', 'v0', '2026-10-10', 1, NULL, NULL)"
            )
            conn.commit()
        before = _snapshot(tmp_path)
        connects.clear()
        assert get_database_info(tmp_path / "gnomad.sqlite", "gnomad") is None
        with AlphaMissenseAnnotator(tmp_path) as annotator:
            annotator._connection()
        assert not [c for c in connects if "gnomad" in c]
        assert _snapshot(tmp_path) == before

    def test_symlinked_wal_cache_reports_wal_and_opens_nothing(
        self, tmp_path: Path, connects: list[str]
    ) -> None:
        """A WAL cache reached through a symlink, with 99 of its 100 rows still in
        the -wal file beside the real file: the version read reports WAL and
        opens nothing, so it cannot read the one checkpointed row as the cache."""
        from allelix.databases.manager import get_database_info

        build = tmp_path / "build"
        real = tmp_path / "real"
        data = tmp_path / "data"
        for d in (build, real, data):
            d.mkdir()
        src = _gnomad(build / "gnomad.sqlite", wal=True)
        with contextlib.closing(sqlite3.connect(src)) as conn:
            conn.execute("PRAGMA wal_autocheckpoint=0")
            conn.executemany(
                "INSERT INTO gnomad_frequencies (chrom, pos, ref, alt, rsid, af) "
                "VALUES ('1', ?, 'A', 'G', NULL, 0.1)",
                [(1000 + i,) for i in range(99)],
            )
            conn.commit()
            for name in ("gnomad.sqlite", "gnomad.sqlite-wal"):
                shutil.copy(build / name, real / name)
        (data / "gnomad.sqlite").symlink_to(real / "gnomad.sqlite")
        assert is_wal(data / "gnomad.sqlite")
        before = _snapshot(real)
        connects.clear()

        assert get_database_info(data / "gnomad.sqlite", "gnomad") is None
        # Wide enough for one table row per line.
        result = CliRunner(env={"COLUMNS": "300"}).invoke(
            main, ["db", "status", "--data-dir", str(data)]
        )
        assert result.exit_code == 0, result.output
        assert _flat(
            "WAL mode from an older Allelix; run `allelix db update` to convert it"
        ) in _flat(result.output)
        assert not [c for c in connects if "gnomad" in c]
        assert _snapshot(real) == before
        # `db status` creates config.toml as before; nothing else appears.
        assert sorted(p.name for p in data.iterdir()) == ["config.toml", "gnomad.sqlite"]

    def test_analyze_update_step_converts_without_downloading(
        self,
        tmp_path: Path,
        build_synthetic_clinvar_cache,
        mock_pharmgkb_dir: Path,
        mock_gwas_tsv: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        """Without --no-update, a WAL cache is converted even when nothing downloads."""
        data = _data_dir(tmp_path, build_synthetic_clinvar_cache, mock_pharmgkb_dir, mock_gwas_tsv)
        refreshed: list[Path] = []
        monkeypatch.setattr("allelix.cli._helpers._maybe_refresh_databases", refreshed.append)
        result = CliRunner().invoke(
            main,
            [
                "analyze",
                str(Path(__file__).parent / "fixtures" / "mock_23andme.txt"),
                "--build",
                "grch37",
                "--data-dir",
                str(data),
                "--output",
                str(tmp_path / "report.json"),
            ],
        )
        assert result.exit_code == 0, result.output
        assert refreshed == [data]
        assert not is_wal(data / "alphamissense.sqlite")
        assert "alphamissense.sqlite: converted to rollback-journal mode" in result.output
        assert "Skipping AlphaMissense" not in result.output


class TestClinVarConnection:
    def test_connection_is_read_only_and_builds_no_index(
        self, tmp_path: Path, build_synthetic_clinvar_cache
    ) -> None:
        db = tmp_path / clinvar_db_filename("GRCh37")
        build_synthetic_clinvar_cache(db, "GRCh37")
        with contextlib.closing(sqlite3.connect(db)) as conn:
            indexes = {r[0] for r in conn.execute("SELECT name FROM sqlite_master")}
        assert "idx_clinvar_position" in indexes
        before = _snapshot(tmp_path)
        with ClinVarAnnotator(tmp_path, builds=("GRCh37",)) as clinvar:
            conn = clinvar._connection("GRCh37")
            assert conn is not None
            with pytest.raises(sqlite3.OperationalError, match="readonly"):
                conn.execute("CREATE TABLE probe (x)")
        assert _snapshot(tmp_path) == before

    def test_read_only_file(self, tmp_path: Path, build_synthetic_clinvar_cache) -> None:
        db = tmp_path / clinvar_db_filename("GRCh37")
        build_synthetic_clinvar_cache(db, "GRCh37")
        db.chmod(stat.S_IRUSR)
        before = _snapshot(tmp_path)
        with ClinVarAnnotator(tmp_path, builds=("GRCh37",)) as clinvar:
            assert clinvar.is_ready()
            conn = clinvar._connection("GRCh37")
            assert conn is not None
            assert conn.execute("SELECT COUNT(*) FROM clinvar_variants").fetchone()[0] > 0
        assert _snapshot(tmp_path) == before


class TestSkippedSourceLines:
    """One console line per state, with the fix (ADR-0048)."""

    @pytest.fixture
    def data(self, tmp_path: Path, build_synthetic_clinvar_cache, mock_pharmgkb_dir: Path) -> Path:
        from allelix.databases.pharmgkb_loader import load_pharmgkb_tsv

        data = tmp_path / "data"
        data.mkdir()
        build_synthetic_clinvar_cache(data / clinvar_db_filename("GRCh37"), "GRCh37")
        build_synthetic_clinvar_cache(data / clinvar_db_filename("GRCh38"), "GRCh38")
        load_pharmgkb_tsv(mock_pharmgkb_dir, data / "pharmgkb.sqlite", source_url="test://")
        return data

    def test_missing(self, data: Path, tmp_path: Path) -> None:
        (data / "pharmgkb.sqlite").unlink()
        code, output = _analyze(data, tmp_path / "r.json")
        assert code == 0, output
        assert _flat("Skipping ClinPGx: not installed. Run `allelix db update`.") in output

    def test_outdated(self, data: Path, tmp_path: Path) -> None:
        with contextlib.closing(sqlite3.connect(data / "pharmgkb.sqlite")) as conn:
            conn.execute("UPDATE database_versions SET local_version_tag = 'iv:0'")
            conn.commit()
        code, output = _analyze(data, tmp_path / "r.json")
        assert code == 0, output
        assert (
            _flat(
                "Skipping ClinPGx: out of date (built for iv:0, needs "
                f"iv:{PHARMGKB_INTERPRETER_VERSION}). Run `allelix db update`."
            )
            in output
        )

    @pytest.mark.parametrize("metadata_allows", [False, True], ids=["metadata", "open-denied"])
    def test_unreadable(
        self,
        data: Path,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
        metadata_allows: bool,
    ) -> None:
        path = data / "pharmgkb.sqlite"
        if metadata_allows:
            access = os.access

            def metadata_access(candidate, mode, *args, **kwargs):
                if candidate == path and mode == os.R_OK:
                    return True
                return access(candidate, mode, *args, **kwargs)

            monkeypatch.setattr(os, "access", metadata_access)
        with _denied_permissions(path, read=True):
            if metadata_allows:
                assert os.access(path, os.R_OK)
            with pytest.raises(PermissionError) as denied:
                path.read_bytes()
            code, output = _analyze(data, tmp_path / "r.json")
        assert code == 0, output
        if metadata_allows:
            assert _flat(f"({denied.value})") in output
        # Metadata can reject the file before open, or the actual denied open
        # supplies useful OS detail. Keep the entire diagnostic together.
        assert any(
            _flat(
                f"Skipping ClinPGx: cannot read {path}{detail}. Check the permissions on {path}."
            )
            in output
            for detail in ("", f" ({denied.value})")
        ), output

    def test_corrupt(self, data: Path, tmp_path: Path) -> None:
        path = data / "pharmgkb.sqlite"
        path.write_bytes(bytes(range(256)) * 32)
        code, output = _analyze(data, tmp_path / "r.json")
        assert code == 0, output
        assert (
            _flat(
                f"Skipping ClinPGx: {path} is damaged (not an SQLite database). "
                "Run `allelix db update`."
            )
            in output
        )

    def test_corrupt_clinvar_names_the_source(self, data: Path, tmp_path: Path) -> None:
        path = data / clinvar_db_filename("GRCh37")
        with contextlib.closing(sqlite3.connect(path)) as conn:
            conn.execute("DROP TABLE clinvar_variants")
            conn.commit()
        code, output = _analyze(data, tmp_path / "r.json")
        assert code == 0, output
        assert (
            _flat(
                f"Skipping ClinVar: {path} is damaged (missing table clinvar_variants). "
                "Run `allelix db update`."
            )
            in output
        )

    def test_readonly_wal(self, data: Path, tmp_path: Path) -> None:
        path = data / "pharmgkb.sqlite"
        with contextlib.closing(sqlite3.connect(path)) as conn:
            conn.execute("PRAGMA journal_mode=WAL")
        with _read_only(data):
            code, output = _analyze(data, tmp_path / "r.json")
        assert code == 0, output
        assert (
            _flat(
                f"Skipping ClinPGx: {path} is in WAL mode, which cannot be read in a data "
                f"directory you cannot write. Run `allelix db update` once with write access "
                f"to {data}."
            )
            in output
        )

    def test_nothing_ready_lists_each_reason(self, tmp_path: Path) -> None:
        data = tmp_path / "data"
        data.mkdir()
        (data / "pharmgkb.sqlite").write_bytes(b"x" * 100)
        code, output = _analyze(data, tmp_path / "r.json")
        assert code == 1
        assert _flat("No annotators are ready. Run `allelix db update` first.") in output
        assert _flat("ClinVar: not installed. Run `allelix db update`.") in output
        assert _flat(f"ClinPGx: {data / 'pharmgkb.sqlite'} is damaged") in output


class TestInstaller:
    def test_wal_in_rollback_out_and_stale_side_files_removed(self, tmp_path: Path) -> None:
        staged = _gnomad(tmp_path / "staged.sqlite", wal=True)
        gz = tmp_path / "gnomad.sqlite.gz"
        with staged.open("rb") as f_in, gzip.open(gz, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        data = tmp_path / "data"
        data.mkdir()
        for suffix in ("-wal", "-shm"):
            (data / f"gnomad.sqlite{suffix}").write_bytes(b"left from the old cache")
        install_prebuilt_gz_cache(gz, data / "gnomad.sqlite", "gnomad")
        assert not is_wal(data / "gnomad.sqlite")
        assert sorted(p.name for p in data.iterdir()) == ["gnomad.sqlite"]

    def test_update_step_converts_an_installed_wal_cache(self, tmp_path: Path) -> None:
        from allelix.cli._helpers import _prepare_caches

        path = _alphamissense(tmp_path / "alphamissense.sqlite", wal=True)
        _prepare_caches([], tmp_path)
        assert not is_wal(path)
        assert sorted(p.name for p in tmp_path.iterdir()) == ["alphamissense.sqlite"]

    def test_update_step_leaves_a_read_only_directory_alone(self, tmp_path: Path) -> None:
        from allelix.cli._helpers import _prepare_caches

        data = tmp_path / "data"
        data.mkdir()
        path = _alphamissense(data / "alphamissense.sqlite", wal=True)
        before = _snapshot(data)
        with _read_only(data):
            _prepare_caches([], data)
        assert is_wal(path)
        assert _snapshot(data) == before


class TestSchemaChecksReadOnly:
    """The loaders' schema checks open a converted cache read-only and never open a WAL one."""

    def _check(self, path: Path, schema_is_current, connects: list[str]) -> None:  # type: ignore[no-untyped-def]
        connects.clear()
        assert schema_is_current(path) is True
        assert connects
        assert all(c.startswith("file:") and c.endswith("?mode=ro") for c in connects), connects
        with contextlib.closing(sqlite3.connect(path)) as conn:
            conn.execute("PRAGMA journal_mode=WAL")
        assert is_wal(path)
        before = _snapshot(path.parent)
        connects.clear()
        assert schema_is_current(path) is False
        assert connects == []
        assert _snapshot(path.parent) == before

    def test_gwas(self, tmp_path: Path, mock_gwas_tsv: Path, connects: list[str]) -> None:
        from allelix.databases.gwas_loader import load_gwas_tsv, schema_is_current

        path = tmp_path / "gwas.sqlite"
        load_gwas_tsv(mock_gwas_tsv, path, source_url="test://")
        self._check(path, schema_is_current, connects)

    def test_pharmgkb(self, tmp_path: Path, mock_pharmgkb_dir: Path, connects: list[str]) -> None:
        from allelix.databases.pharmgkb_loader import load_pharmgkb_tsv, schema_is_current

        path = tmp_path / "pharmgkb.sqlite"
        load_pharmgkb_tsv(mock_pharmgkb_dir, path, source_url="test://")
        self._check(path, schema_is_current, connects)


def test_base_defaults_for_an_annotator_without_a_cache() -> None:
    """``cache_status()`` follows ``is_ready()``; ``prepare_cache()`` changes nothing."""
    from allelix.annotators.base import Annotator
    from allelix.annotators.gnomad import GnomadAnnotator

    class _NoCache(Annotator):
        name = "nocache"
        display_name = "No cache"
        license = GnomadAnnotator.license
        ready = True

        def setup(self) -> None: ...

        def annotate(self, variant):  # type: ignore[no-untyped-def]
            return []

        def is_ready(self) -> bool:
            return self.ready

        def version(self) -> str | None:
            return None

        def close(self) -> None: ...

        def fetch_remote_signal(self) -> str | None:
            return None

        def cached_remote_signal(self) -> str | None:
            return None

    annotator = _NoCache(Path("/nonexistent"))
    assert annotator.cache_status().state is CacheState.READY
    annotator.ready = False
    assert annotator.cache_status().state is CacheState.MISSING
    assert annotator.prepare_cache() is False


def test_sources_on_the_analysis_result(
    tmp_path: Path,
    build_synthetic_clinvar_cache,
    mock_pharmgkb_dir: Path,
    mock_gwas_tsv: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The result keeps each source's part in the run for the reports."""
    import dataclasses

    from allelix.cli import _helpers
    from allelix.reports._pipeline import AnalysisResult

    fields = {f.name: f for f in dataclasses.fields(AnalysisResult)}
    assert fields["sources"].default_factory() == []  # type: ignore[misc]
    assert fields["report_warnings"].default_factory() == []  # type: ignore[misc]

    data = _data_dir(tmp_path, build_synthetic_clinvar_cache, mock_pharmgkb_dir, mock_gwas_tsv)
    seen: list[AnalysisResult] = []
    real = _helpers.run_analysis

    def _capture(*args, **kwargs):  # type: ignore[no-untyped-def]
        result = real(*args, **kwargs)
        seen.append(result)
        return result

    monkeypatch.setattr(_helpers, "run_analysis", _capture)
    code, output = _analyze(data, tmp_path / "r.json")
    assert code == 0, output
    skipped = [(s.name, s.reason) for s in seen[0].sources if s.state == "skipped"]
    assert skipped == [("alphamissense", "wal"), ("cadd", "disabled")]


def test_connect_readonly_reads_every_committed_row_through_a_symlink(tmp_path: Path) -> None:
    """A plain ``mode=ro`` open sees rows still in the -wal file beside the real
    file; ``immutable=1`` would read only the checkpointed ones (PR review)."""
    from allelix.databases.cache_state import connect_readonly

    build = tmp_path / "build"
    real = tmp_path / "real"
    data = tmp_path / "data"
    for d in (build, real, data):
        d.mkdir()
    src = _gnomad(build / "gnomad.sqlite", wal=True)
    with contextlib.closing(sqlite3.connect(src)) as conn:
        conn.execute("PRAGMA wal_autocheckpoint=0")
        conn.executemany(
            "INSERT INTO gnomad_frequencies (chrom, pos, ref, alt, rsid, af) "
            "VALUES ('1', ?, 'A', 'G', NULL, 0.1)",
            [(1000 + i,) for i in range(99)],
        )
        conn.commit()
        for name in ("gnomad.sqlite", "gnomad.sqlite-wal"):
            shutil.copy(build / name, real / name)
    (data / "gnomad.sqlite").symlink_to(real / "gnomad.sqlite")
    with contextlib.closing(connect_readonly(data / "gnomad.sqlite")) as conn:
        assert conn.execute("SELECT COUNT(*) FROM gnomad_frequencies").fetchone()[0] == 100


def test_cache_state_values_are_stable() -> None:
    assert [s.value for s in CacheState] == [
        "ready",
        "missing",
        "outdated",
        "unreadable",
        "corrupt",
        "wal",
        "readonly_wal",
    ]


_ANNOTATORS = [
    "clinvar",
    "pharmgkb",
    "gwas",
    "snpedia",
    "gnomad",
    "alphamissense",
    "cadd",
    "rsid_identity",
]


def _annotator(name: str, data: Path):
    from allelix.annotators.alphamissense import AlphaMissenseAnnotator
    from allelix.annotators.cadd import CaddAnnotator
    from allelix.annotators.gnomad import GnomadAnnotator
    from allelix.annotators.gwas import GWASCatalogAnnotator
    from allelix.annotators.pharmgkb import PharmGKBAnnotator
    from allelix.annotators.rsid_identity import RsidIdentityAnnotator
    from allelix.annotators.snpedia import SNPediaAnnotator

    return {
        "clinvar": lambda: ClinVarAnnotator(data, builds=("GRCh37",)),
        "pharmgkb": lambda: PharmGKBAnnotator(data),
        "gwas": lambda: GWASCatalogAnnotator(data),
        "snpedia": lambda: SNPediaAnnotator(data),
        "gnomad": lambda: GnomadAnnotator(data),
        "alphamissense": lambda: AlphaMissenseAnnotator(data),
        "cadd": lambda: CaddAnnotator(data),
        "rsid_identity": lambda: RsidIdentityAnnotator(data),
    }[name]()


class TestEveryAnnotator:
    @pytest.fixture
    def data(
        self,
        tmp_path: Path,
        build_synthetic_clinvar_cache,
        mock_pharmgkb_dir: Path,
        mock_gwas_tsv: Path,
    ) -> Path:
        from allelix.databases.loader_utils import convert_to_rollback_journal

        data = _data_dir(tmp_path, build_synthetic_clinvar_cache, mock_pharmgkb_dir, mock_gwas_tsv)
        assert convert_to_rollback_journal(data / "alphamissense.sqlite")
        return data

    @pytest.mark.parametrize("name", _ANNOTATORS)
    def test_connection_is_read_only(self, data: Path, name: str) -> None:
        before = _snapshot(data)
        with _annotator(name, data) as annotator:
            assert annotator.is_ready()
            conn = (
                annotator._connection("GRCh37") if name == "clinvar" else annotator._connection()
            )
            with pytest.raises(sqlite3.OperationalError, match="readonly"):
                conn.execute("CREATE TABLE probe (x)")
        assert _snapshot(data) == before

    @pytest.mark.parametrize("name", _ANNOTATORS)
    def test_missing_and_corrupt(self, tmp_path: Path, name: str) -> None:
        with _annotator(name, tmp_path) as annotator:
            assert annotator.cache_status().state is CacheState.MISSING
        filename = "clinvar.GRCh37.sqlite" if name == "clinvar" else f"{name}.sqlite"
        (tmp_path / filename).write_bytes(b"\x00" * 4096)
        with _annotator(name, tmp_path) as annotator:
            assert annotator.cache_status().state is CacheState.CORRUPT
            assert not annotator.is_ready()

    @pytest.mark.parametrize("name", _ANNOTATORS)
    def test_empty_main_table_is_corrupt(self, data: Path, name: str) -> None:
        file, table = _MAIN_TABLE[name]
        with contextlib.closing(sqlite3.connect(data / file)) as conn:
            conn.execute(f"DELETE FROM {table}")
            conn.commit()
        with _annotator(name, data) as annotator:
            assert annotator.cache_status().state is CacheState.CORRUPT
            assert not annotator.is_ready()


_MAIN_TABLE = {
    "clinvar": ("clinvar.GRCh37.sqlite", "clinvar_variants"),
    "pharmgkb": ("pharmgkb.sqlite", "pharmgkb_annotations"),
    "gwas": ("gwas.sqlite", "gwas_associations"),
    "snpedia": ("snpedia.sqlite", "snpedia_genotypes"),
    "gnomad": ("gnomad.sqlite", "gnomad_frequencies"),
    "alphamissense": ("alphamissense.sqlite", "alphamissense_scores"),
    "cadd": ("cadd.sqlite", "cadd_scores"),
    "rsid_identity": ("rsid_identity.sqlite", "rsid_identity"),
}


class TestDbUpdateDownloadsOnlyWhatIsNeeded:
    """Plain ``db update`` replaces a damaged cache only, and never re-downloads a WAL one."""

    @pytest.fixture
    def setups(self, monkeypatch: pytest.MonkeyPatch) -> list[str]:
        """Record each source ``db update`` fetches; no network, no remote change."""
        import socket

        from allelix.cli import db as db_module

        calls: list[str] = []
        monkeypatch.setattr(
            "allelix.cli._helpers._run_setup",
            lambda annotator, **kwargs: calls.append(annotator.name) or True,
        )
        real = db_module.get_annotators

        def _annotators(*args, **kwargs):  # type: ignore[no-untyped-def]
            annotators = real(*args, **kwargs)
            for annotator in annotators:
                annotator.fetch_remote_signal = annotator.cached_remote_signal  # type: ignore[method-assign]
            return annotators

        monkeypatch.setattr(db_module, "get_annotators", _annotators)

        def _no_network(*_a, **_k):  # type: ignore[no-untyped-def]
            raise AssertionError("network used")

        monkeypatch.setattr(socket, "socket", _no_network)
        return calls

    def test_damaged_cache_is_the_only_source_fetched(
        self,
        tmp_path: Path,
        build_synthetic_clinvar_cache,
        mock_pharmgkb_dir: Path,
        mock_gwas_tsv: Path,
        setups: list[str],
    ) -> None:
        from allelix.databases.loader_utils import convert_to_rollback_journal

        data = _data_dir(tmp_path, build_synthetic_clinvar_cache, mock_pharmgkb_dir, mock_gwas_tsv)
        assert convert_to_rollback_journal(data / "alphamissense.sqlite")
        (data / "pharmgkb.sqlite").write_bytes(b"x" * 4096)
        others = {k: v for k, v in _snapshot(data).items() if k != "pharmgkb.sqlite"}
        result = CliRunner().invoke(main, ["db", "update", "--data-dir", str(data)])
        assert result.exit_code == 0, result.output
        assert setups == ["pharmgkb"]
        after = _snapshot(data)
        assert {k: after[k] for k in others} == others
        assert sorted(after) == sorted([*others, "pharmgkb.sqlite", "config.toml"])

    def test_wal_cache_that_cannot_be_converted_is_not_downloaded(
        self,
        tmp_path: Path,
        build_synthetic_clinvar_cache,
        mock_pharmgkb_dir: Path,
        mock_gwas_tsv: Path,
        setups: list[str],
    ) -> None:
        data = _data_dir(tmp_path, build_synthetic_clinvar_cache, mock_pharmgkb_dir, mock_gwas_tsv)
        before = _snapshot(data)
        with _read_only(data):
            result = CliRunner().invoke(main, ["db", "update", "--data-dir", str(data)])
        assert result.exit_code == 1, result.output
        assert setups == []
        # Rich wraps at 80 columns and the message holds the tmp path.
        output = " ".join(result.output.split())
        assert "alphamissense:" in output
        assert "Nothing was downloaded" in output
        assert _snapshot(data) == before


class TestUpdateStepOnly:
    def test_snpedia_raw_pages_are_not_parsed_by_is_ready(self, tmp_path: Path) -> None:
        from allelix.annotators.snpedia import SNPediaAnnotator

        path = tmp_path / "snpedia.sqlite"
        with contextlib.closing(sqlite3.connect(path)) as conn:
            conn.execute("CREATE TABLE _raw_pages (title TEXT, category TEXT, content TEXT)")
            conn.execute("INSERT INTO _raw_pages VALUES ('Rs900780009', 'snp', '')")
            conn.commit()
        before = _snapshot(tmp_path)
        with SNPediaAnnotator(tmp_path) as snpedia:
            status = snpedia.cache_status()
            assert snpedia.is_ready() is False
        assert (status.state, status.detail) == (CacheState.OUTDATED, "pages not parsed yet")
        assert _snapshot(tmp_path) == before

    def test_db_update_converts_a_wal_cache_without_downloading(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        path = _alphamissense(tmp_path / "alphamissense.sqlite", wal=True)
        monkeypatch.setattr("allelix.cli.db.get_annotators", lambda *_a, **_k: [])
        result = CliRunner().invoke(main, ["db", "update", "--data-dir", str(tmp_path)])
        assert result.exit_code == 0, result.output
        assert "alphamissense.sqlite: converted to rollback-journal mode" in result.output
        assert not is_wal(path)

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Small offline CLI fixtures for destination and diagnostic contracts."""

from __future__ import annotations

import contextlib
import os
import socket
import sqlite3
from types import SimpleNamespace

import pytest
from click.testing import CliRunner

from allelix.annotators.clinvar import clinvar_db_filename, clinvar_record_name
from allelix.cli import _helpers, main
from allelix.config import AllelixConfig, save_config
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION, GNOMAD_SCHEMA_VERSION
from allelix.databases.schema import CLINVAR_SCHEMA, GNOMAD_SCHEMA


@pytest.fixture
def rc_case(tmp_path, monkeypatch):
    """Use real parsers, annotators and renderers with one synthetic site."""
    for name in ("data", "out", "home", "xdg", "cache"):
        (tmp_path / name).mkdir()
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg"))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    monkeypatch.setenv("ALLELIX_DATA_DIR", str(tmp_path / "data"))
    monkeypatch.setenv("COLUMNS", "240")
    monkeypatch.setenv("NO_COLOR", "1")

    def no_network(*args, **kwargs):
        raise AssertionError("synthetic CLI checks must not use the network")

    monkeypatch.setattr(socket.socket, "connect", no_network)
    monkeypatch.setattr(socket, "create_connection", no_network)
    data = tmp_path / "data"
    cfg = AllelixConfig()
    cfg.sources = {name: name == "clinvar" for name in cfg.sources}
    save_config(data, cfg)
    rsid = "rs900252501"
    for build in ("GRCh37", "GRCh38"):
        with contextlib.closing(sqlite3.connect(data / clinvar_db_filename(build))) as conn:
            conn.executescript(CLINVAR_SCHEMA)
            conn.execute(
                "INSERT INTO clinvar_variants (rsid,chromosome,position,ref,alt,"
                "clinical_significance,condition,gene,allele_id) VALUES (?,?,?,?,?,?,?,?,?)",
                (
                    rsid,
                    "1",
                    1000,
                    "A",
                    "G",
                    "Pathogenic",
                    "Synthetic condition",
                    "MTHFR",
                    900252501,
                ),
            )
            conn.execute(
                "INSERT INTO database_versions (name,source_url,version,downloaded_at,"
                "record_count,local_version_tag) VALUES (?,'test://synthetic','synthetic',"
                "'2001-01-01',1,?)",
                (clinvar_record_name(build), f"iv:{CLINVAR_INTERPRETER_VERSION}"),
            )
            conn.commit()
    with contextlib.closing(sqlite3.connect(data / "gnomad.sqlite")) as conn:
        conn.executescript(GNOMAD_SCHEMA)
        conn.execute(
            "INSERT INTO gnomad_frequencies (chrom,pos,ref,alt,rsid,af) "
            "VALUES ('1',1000,'A','G',?,0.01)",
            (rsid,),
        )
        conn.execute(
            "INSERT INTO database_versions VALUES "
            "('gnomad','test://synthetic','1','2001-01-01',1,NULL,?)",
            (f"sv:{GNOMAD_SCHEMA_VERSION}",),
        )
        conn.commit()
    source = tmp_path / "input.vcf"
    header = (
        "##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    )
    row = f"1\t1000\t{rsid}\tA\tG\t.\tPASS\t.\tGT\t0/1\n"
    source.write_text(header + row)
    runner = CliRunner()

    def invoke(args):
        return runner.invoke(main, [str(arg) for arg in args], terminal_width=240)

    def report_args(command):
        return [command, source, "--data-dir", data, "--build", "grch37", "--no-update"]

    yield SimpleNamespace(
        root=tmp_path,
        data=data,
        out=tmp_path / "out",
        source=source,
        rsid=rsid,
        header=header,
        row=row,
        invoke=invoke,
        report_args=report_args,
    )
    # Restore permissions so failed assertions never strand pytest's own files.
    for path in tmp_path.rglob("*"):
        if not path.is_symlink():
            with contextlib.suppress(FileNotFoundError):
                path.chmod(0o700 if path.is_dir() else 0o600)


@pytest.fixture
def rc_snapshot():
    def snapshot(root):
        result = {}
        for path in root.rglob("*"):
            key = str(path.relative_to(root))
            if path.is_symlink():
                result[key] = ("symlink", os.readlink(path))
            elif path.is_file():
                info = path.stat()
                result[key] = (path.read_bytes(), info.st_ino, info.st_mtime_ns, info.st_mode)
            else:
                result[key] = ("directory",)
        return result

    return snapshot


@pytest.fixture
def rc_clean_error():
    def check(result, path, code=1):
        assert result.exit_code == code, (result.output, repr(result.exception))
        assert isinstance(result.exception, SystemExit), repr(result.exception)
        assert "error" in result.output.lower(), result.output
        assert str(path) in result.output, result.output
        assert "Traceback" not in result.output, result.output

    return check


@pytest.fixture
def rc_no_analysis(monkeypatch):
    calls = []

    def forbidden(*args, **kwargs):
        calls.append("analysis started")
        raise AssertionError("refusal must precede analysis")

    monkeypatch.setattr(_helpers, "run_analysis", forbidden)
    monkeypatch.setattr(_helpers, "_resolved_variants", forbidden)
    return calls

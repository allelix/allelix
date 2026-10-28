# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Permission errors cover new data directories and non-ClinVar cache files."""

import importlib

import pytest

from allelix.annotators.gnomad import GnomadAnnotator
from allelix.config import AllelixConfig, save_config

COMMANDS = [
    "analyze",
    "methylation",
    "pharmacogenomics",
    "stats",
    "extract",
    "export plink",
    "db update",
    "db status",
    "db clean",
    "db path",
    "db path --check",
    "config get",
    "config set",
    "config show",
]


def _args(case, command):
    args = command.split()
    if command in ("analyze", "methylation", "pharmacogenomics"):
        return [
            *args,
            case.source,
            "--build",
            "grch37",
            "--no-update",
            "--output",
            case.out / "report.json",
        ]
    if command == "stats":
        return [*args, case.source, "--build", "grch37"]
    if command == "extract":
        return [*args, case.source, "--build", "grch37", "--snps", case.rsid]
    if command == "export plink":
        return [*args, case.source, "--build", "grch37", "--output-prefix", case.out / "export"]
    if command == "db update":
        return [*args, "--only", "clinvar", "--force"]
    if command == "db clean":
        return [*args, "--yes"]
    if command == "config get":
        return [*args, "sources.clinvar"]
    if command == "config set":
        return [*args, "sources.clinvar", "false"]
    return args


@pytest.mark.parametrize("command", COMMANDS)
@pytest.mark.parametrize("path_source", ["option", "environment"])
def test_new_data_directory_below_readonly_parent_is_a_clean_error(
    rc_case, rc_snapshot, rc_clean_error, monkeypatch, command, path_source
):
    case = rc_case
    parent = case.root / "readonly-parent"
    parent.mkdir()
    parent.chmod(0o555)
    destination = parent / "new-data"
    args = _args(case, command)
    if path_source == "option":
        args += ["--data-dir", destination]
    else:
        monkeypatch.setenv("ALLELIX_DATA_DIR", str(destination))
    before = rc_snapshot(case.root)
    result = case.invoke(args)
    rc_clean_error(result, destination)
    assert rc_snapshot(case.root) == before


@pytest.mark.parametrize("path_source", ["option", "environment"])
def test_db_path_check_keeps_its_existing_readonly_directory_diagnostic(
    rc_case, rc_snapshot, rc_clean_error, monkeypatch, path_source
):
    case = rc_case
    destination = case.root / "existing-readonly"
    destination.mkdir()
    destination.chmod(0o555)
    args = ["db", "path", "--check"]
    if path_source == "option":
        args += ["--data-dir", destination]
    else:
        monkeypatch.setenv("ALLELIX_DATA_DIR", str(destination))
    before = rc_snapshot(case.root)
    result = case.invoke(args)
    rc_clean_error(result, destination)
    assert "path is not writable" in result.output.lower()
    assert rc_snapshot(case.root) == before


def test_nonclinvar_readonly_cache_stops_before_download(
    rc_case, rc_snapshot, rc_clean_error, monkeypatch
):
    case = rc_case
    cfg = AllelixConfig()
    cfg.sources = {name: name == "gnomad" for name in cfg.sources}
    save_config(case.data, cfg)
    annotator = GnomadAnnotator(case.data)
    db = importlib.import_module("allelix.cli.db")
    monkeypatch.setattr(db, "get_annotators", lambda *args, **kwargs: [annotator])
    calls = []

    def forbidden(*args, **kwargs):
        calls.append("download")
        raise AssertionError("readonly non-ClinVar cache reached source setup")

    monkeypatch.setattr(annotator, "setup", forbidden)
    cache = case.data / "gnomad.sqlite"
    cache.chmod(0o444)
    before = rc_snapshot(case.root)
    result = case.invoke(["db", "update", "--data-dir", case.data, "--force"])
    assert not calls, result.output
    rc_clean_error(result, cache)
    assert rc_snapshot(case.root) == before

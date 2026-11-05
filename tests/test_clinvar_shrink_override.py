# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""An explicit one-run shrink permission does not weaken other build safeguards."""

from __future__ import annotations

import re

import pytest
from click.testing import CliRunner

from allelix.annotators import clinvar
from allelix.cli import main
from allelix.databases import manager
from tests.test_clinvar_integrity import BUILDS, SS, VS, _assert_old, _network, _old_caches


def _update(data, build, *, allow=False, force=False):
    args = ["db", "update", "--only", "clinvar", "--data-dir", str(data), "--build", build]
    if allow:
        args.append("--allow-clinvar-shrink")
    if force:
        args.append("--force")
    return CliRunner().invoke(main, args)


def _counts(data, expected):
    for build, count in expected.items():
        info = manager.get_database_info(
            data / clinvar.clinvar_db_filename(build), clinvar.clinvar_record_name(build)
        )
        assert info is not None and info["record_count"] == count
        assert info["local_version_tag"] == f"iv:{clinvar.CLINVAR_INTERPRETER_VERSION}"


def _count_evidence(text, build, old, new):
    text = " ".join(text.replace(",", "").split())
    # Counts must be attributable to a build; exact prose and punctuation are free.
    start = text.find(build)
    assert start >= 0, text
    window = text[start:]
    assert re.search(
        rf"\b(?:old|previous)(?: count| rows| row count)?\s*[:=]?\s*{old}\b", window, re.I
    ), text
    assert re.search(
        rf"\b(?:new|candidate|final)(?: has| count| rows| row count)?\s*[:=]?\s*{new}\b",
        window,
        re.I,
    ), text


@pytest.mark.parametrize("build", ["grch37", "grch38", "both"])
@pytest.mark.parametrize("force", [False, True])
def test_relative_floor_can_be_lifted_for_this_update(tmp_path, monkeypatch, build, force):
    data = tmp_path / "data"
    _old_caches(data, dict.fromkeys(BUILDS, 10))
    _network(monkeypatch, dict.fromkeys(BUILDS, 8))
    result = _update(data, build, allow=True, force=force)
    assert result.exit_code == 0, result.output
    selected = BUILDS if build == "both" else ("GRCh" + build[-2:],)
    _counts(data, {b: 8 if b in selected else 10 for b in BUILDS})
    log = (data / manager.CLINVAR_BUILD_LOG).read_text()
    for b in selected:
        _count_evidence(result.output, b, 10, 8)
        _count_evidence(log, b, 10, 8)


@pytest.mark.parametrize("build", ["grch37", "grch38", "both"])
@pytest.mark.parametrize("force", [False, True])
def test_relative_refusal_preserves_cache_and_names_the_override(
    tmp_path, monkeypatch, build, force
):
    data = tmp_path / "data"
    before = _old_caches(data, dict.fromkeys(BUILDS, 10))
    _network(monkeypatch, dict.fromkeys(BUILDS, 8))
    result = _update(data, build, force=force)
    assert result.exit_code == 1, result.output
    _assert_old(data, before)
    assert "--allow-clinvar-shrink" in result.output
    selected = "GRCh37" if build == "both" else "GRCh" + build[-2:]
    _count_evidence(result.output, selected, 10, 8)
    _count_evidence((data / manager.CLINVAR_BUILD_LOG).read_text(), selected, 10, 8)


@pytest.mark.parametrize("build", ["grch37", "grch38", "both"])
@pytest.mark.parametrize("allow", [False, True])
@pytest.mark.parametrize("new,floor", [(2, 3), (0, 0)], ids=["fixed-floor", "zero"])
def test_fixed_floor_and_zero_are_never_overridden(
    tmp_path, monkeypatch, build, allow, new, floor
):
    data = tmp_path / "data"
    before = _old_caches(data, dict.fromkeys(BUILDS, 10))
    _network(monkeypatch, dict.fromkeys(BUILDS, new))
    monkeypatch.setattr(manager, "CLINVAR_MIN_ROWS", dict.fromkeys(BUILDS, floor))
    result = _update(data, build, allow=allow, force=True)
    assert result.exit_code == 1, result.output
    _assert_old(data, before)
    assert "no such option" not in result.output.lower()
    assert "existing caches kept" in result.output.lower()


@pytest.mark.parametrize("build", ["grch37", "grch38", "both"])
def test_permission_is_not_saved_or_used_by_the_next_run(tmp_path, monkeypatch, build):
    data = tmp_path / "data"
    _old_caches(data, dict.fromkeys(BUILDS, 10))
    config = data / "config.toml"
    config.write_text("[sources]\nclinvar = true\n")
    original_config = config.read_bytes()
    _network(monkeypatch, dict.fromkeys(BUILDS, 8))
    first = _update(data, build, allow=True, force=True)
    assert first.exit_code == 0, first.output
    snapshots = {b: (data / clinvar.clinvar_db_filename(b)).read_bytes() for b in BUILDS}
    _network(monkeypatch, dict.fromkeys(BUILDS, 6))
    second = _update(data, build, force=True)
    assert second.exit_code == 1, second.output
    _assert_old(data, snapshots)
    assert config.read_bytes() == original_config
    assert "--allow-clinvar-shrink" in second.output


@pytest.mark.parametrize("source", [VS, SS])
def test_shrink_permission_never_bypasses_download_verification(tmp_path, monkeypatch, source):
    data = tmp_path / "data"
    before = _old_caches(data, dict.fromkeys(BUILDS, 10))
    _, sidecars, _ = _network(monkeypatch, dict.fromkeys(BUILDS, 8))
    sidecars[source] = f"{'0' * 32}  {source}\n"
    result = _update(data, "both", allow=True, force=True)
    assert result.exit_code == 1, result.output
    assert "integrity check failed" in result.output.lower()
    _assert_old(data, before)


def test_both_builds_are_kept_when_the_second_candidate_is_empty(tmp_path, monkeypatch):
    data = tmp_path / "data"
    before = _old_caches(data, dict.fromkeys(BUILDS, 10))
    _network(monkeypatch, {"GRCh37": 8, "GRCh38": 0})
    result = _update(data, "both", allow=True, force=True)
    assert result.exit_code == 1, result.output
    assert "0 rows" in result.output
    _assert_old(data, before)


@pytest.mark.parametrize("new,success", [(10, True), (9, False)])
def test_default_relative_floor_keeps_its_rounding(tmp_path, monkeypatch, new, success):
    data = tmp_path / "data"
    before = _old_caches(data, dict.fromkeys(BUILDS, 11))
    _network(monkeypatch, dict.fromkeys(BUILDS, new))
    result = _update(data, "both", force=True)
    assert result.exit_code == (0 if success else 1), result.output
    if success:
        _counts(data, dict.fromkeys(BUILDS, new))
    else:
        _assert_old(data, before)


def test_permission_does_not_itself_force_a_download(tmp_path, monkeypatch):
    data = tmp_path / "data"
    before = _old_caches(data, dict.fromkeys(BUILDS, 10))
    _, _, downloads = _network(monkeypatch, dict.fromkeys(BUILDS, 10))
    monkeypatch.setattr(clinvar.ClinVarAnnotator, "fetch_remote_signal", lambda self: "md5:old")
    result = _update(data, "both", allow=True)
    assert result.exit_code == 0, result.output
    assert downloads == []
    _assert_old(data, before)


def test_help_states_the_narrow_one_run_scope():
    result = CliRunner().invoke(main, ["db", "update", "--help"])
    assert result.exit_code == 0
    text = " ".join(result.output.lower().split())
    assert "--allow-clinvar-shrink" in text
    assert "90%" in text and "fixed" in text and "zero" in text

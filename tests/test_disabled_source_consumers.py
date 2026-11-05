# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Source choices also isolate resolver, extraction and export dependencies."""

import pytest
from click.testing import CliRunner
from test_disabled_source_reports import _configure, _park, _prepare
from test_extract_reference_calls import _input, _row
from validated_call_helpers import BUILDS, RSID

from allelix.annotators import get_annotators
from allelix.cli import main, utility

RESOLVERS = ("clinvar", "gnomad", "rsid_identity")


def _isolate(root, build, source):
    data, path = _prepare(root, build)
    for other in RESOLVERS:
        if other != source:
            _park(data, other)
    _configure(data, source)
    return data, path


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("source", RESOLVERS)
@pytest.mark.parametrize("indexed", [False, True])
@pytest.mark.parametrize("gt", ["0/1", "./."])
def test_extract_source_off_equals_cache_absent(tmp_path, monkeypatch, build, source, indexed, gt):
    from validated_call_helpers import POSITIONS

    data, _ = _isolate(tmp_path, build, source)
    path = _input(tmp_path, build, [_row(POSITIONS[build], gt, name=".")], indexed=indexed)
    found = []
    monkeypatch.setattr(utility, "_render_extract_table", lambda p, w, f: found.append(f))
    for present in (True, False):
        if not present:
            _park(data, source)
        result = CliRunner().invoke(
            main,
            [
                "extract",
                str(path),
                "--snps",
                RSID,
                "--build",
                build.lower(),
                "--data-dir",
                str(data),
            ],
        )
        assert result.exit_code == 0, (result.output, result.exception)
    assert len(found) == 2 and found[0] == found[1]


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("source", RESOLVERS)
@pytest.mark.parametrize("named", [False, True])
def test_plink_source_off_equals_cache_absent(tmp_path, build, source, named):
    from validated_call_helpers import POSITIONS, call

    data, _ = _isolate(tmp_path, build, source)
    path = call(tmp_path / "input.vcf", POSITIONS[build], rsid=RSID if named else ".", gt="1/1")
    outputs = []
    for present in (True, False):
        if not present:
            _park(data, source)
        prefix = tmp_path / str(present)
        result = CliRunner().invoke(
            main,
            [
                "export",
                "plink",
                str(path),
                "--output-prefix",
                str(prefix),
                "--build",
                build.lower(),
                "--data-dir",
                str(data),
            ],
        )
        assert result.exit_code == 0, (result.output, result.exception)
        outputs.append(
            tuple(prefix.with_suffix(suffix).read_bytes() for suffix in (".bed", ".bim", ".fam"))
        )
    assert outputs[0] == outputs[1]


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("source", RESOLVERS)
@pytest.mark.parametrize("command", ["stats", "extract", "plink"])
def test_non_report_commands_never_query_disabled_sources(
    tmp_path, monkeypatch, build, source, command
):
    from validated_call_helpers import POSITIONS, call

    data, _ = _isolate(tmp_path, build, source)
    path = call(tmp_path / "input.vcf", POSITIONS[build], rsid=".")
    registry = get_annotators(data)
    cls = type(next(a for a in registry if a.name == source))
    queried = []
    original = cls._connection

    def track(self, *args, **kwargs):
        queried.append(source)
        return original(self, *args, **kwargs)

    monkeypatch.setattr(cls, "_connection", track)
    common = [str(path), "--build", build.lower(), "--data-dir", str(data)]
    if command == "plink":
        args = ["export", "plink", *common, "--output-prefix", str(tmp_path / "export")]
    else:
        args = [command, *common]
        if command == "extract":
            args += ["--snps", RSID]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, (result.output, result.exception)
    assert not queried

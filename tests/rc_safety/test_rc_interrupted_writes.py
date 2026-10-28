# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""A real SIGINT after the first write preserves destinations and removes stages."""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

import allelix

# Intercept filesystem streams, not a renderer or its output. Tracking os.open
# descriptors also covers staged writers without requiring Linux /proc paths.
INTERRUPT_CHILD = r"""
import builtins
import io
import json
import os
import signal
import socket
import sys
from pathlib import Path

parent = Path(sys.argv[1]).resolve()
args = json.loads(sys.argv[2])
fd_paths = {}
original_os_open = os.open
original_fdopen = os.fdopen
original_io_open = io.open
original_open = builtins.open

def no_network(*args, **kwargs):
    raise AssertionError("no network in a synthetic interruption check")
socket.socket.connect = no_network
socket.create_connection = no_network

class Interrupted:
    def __init__(self, stream):
        self.stream = stream
    def __getattr__(self, name):
        return getattr(self.stream, name)
    def __enter__(self):
        return self
    def __exit__(self, *args):
        return self.stream.__exit__(*args)
    def write(self, value):
        self.stream.write(value[:max(1, min(16, len(value)))])
        self.stream.flush()
        print("INTERRUPT_REACHED", file=sys.stderr, flush=True)
        os.kill(os.getpid(), signal.SIGINT)
        raise AssertionError("SIGINT did not interrupt")
    def writelines(self, values):
        for value in values:
            self.write(value)

def wrap(stream, file, mode):
    path = fd_paths.get(file) if isinstance(file, int) else Path(file).resolve()
    if path is not None and path.parent == parent and any(c in mode for c in "wax+"):
        return Interrupted(stream)
    return stream

def os_open(file, *args, **kwargs):
    fd = original_os_open(file, *args, **kwargs)
    fd_paths[fd] = Path(file).resolve()
    return fd

def fdopen(fd, mode="r", *args, **kwargs):
    return wrap(original_fdopen(fd, mode, *args, **kwargs), fd, mode)

def io_open(file, mode="r", *args, **kwargs):
    return wrap(original_io_open(file, mode, *args, **kwargs), file, mode)

def builtin_open(file, mode="r", *args, **kwargs):
    return wrap(original_open(file, mode, *args, **kwargs), file, mode)

os.open = os_open
os.fdopen = fdopen
io.open = io_open
builtins.open = builtin_open
from allelix.cli import main
main(args=args, prog_name="allelix")
"""

WRITERS = [
    (command, kind)
    for command in ("analyze", "methylation", "pharmacogenomics")
    for kind in ("json", "html")
] + [("export plink", "plink"), ("config set", "config")]


def _interrupted(case, args, directory):
    env = os.environ.copy()
    env["PYTHONPATH"] = str(Path(allelix.__file__).resolve().parent.parent)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    return subprocess.run(
        [
            sys.executable,
            "-c",
            INTERRUPT_CHILD,
            str(directory),
            json.dumps([str(arg) for arg in args]),
        ],
        cwd=case.root,
        env=env,
        capture_output=True,
        text=True,
        timeout=20,
    )


@pytest.mark.parametrize("command,kind", WRITERS)
@pytest.mark.parametrize("existing", [False, True])
def test_interrupted_write_preserves_the_previous_state(
    rc_case,
    rc_snapshot,
    command,
    kind,
    existing,
):
    case = rc_case
    directory = case.out
    if kind == "config":
        directory = case.data
        if not existing:
            (case.data / "config.toml").unlink()
        args = ["config", "set", "--data-dir", case.data, "sources.clinvar", "false"]
    elif kind == "plink":
        prefix = case.out / "export"
        if existing:
            for suffix in (".bed", ".bim", ".fam"):
                prefix.with_suffix(suffix).write_bytes(b"previous destination\n")
        args = [
            "export",
            "plink",
            case.source,
            "--data-dir",
            case.data,
            "--build",
            "grch37",
            "--output-prefix",
            prefix,
        ]
    else:
        output = case.out / f"report.{kind}"
        if existing:
            output.write_text("previous destination\n")
        args = [*case.report_args(command), "--output", output]
    before = rc_snapshot(directory)
    result = _interrupted(case, args, directory)
    assert "INTERRUPT_REACHED" in result.stderr, (result.stdout, result.stderr)
    assert result.returncode != 0
    assert rc_snapshot(directory) == before, (result.stdout, result.stderr)
    assert not list(directory.glob("*.tmp"))
    assert not list(directory.glob(".*.tmp"))


@pytest.mark.parametrize("suffix", [".vcf", ".vcf.gz"])
def test_interrupted_new_vcf_leaves_no_destination_or_stage(rc_case, rc_snapshot, suffix):
    case = rc_case
    output = case.out / ("annotated" + suffix)
    before = rc_snapshot(case.out)
    result = _interrupted(case, [*case.report_args("analyze"), "--vcf-out", output], case.out)
    assert "INTERRUPT_REACHED" in result.stderr, (result.stdout, result.stderr)
    assert result.returncode != 0
    assert rc_snapshot(case.out) == before
    assert not output.exists()

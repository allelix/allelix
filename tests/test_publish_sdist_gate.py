"""The uploaded sdist must be installed and checked independently of the wheel."""

import os
import re
import subprocess
import sys

import pytest
from test_publish_release_gates import _gate_steps, _workflow


def _sdist_step():
    steps = [
        s
        for s in _gate_steps()
        if "pip install" in s.get("run", "") and re.search(r"tar\.gz|SDIST", s.get("run", ""))
    ]
    # The wheel step mentions tar.gz only in its dist-shape precondition;
    # require the actual install command to name the sdist variable or file.
    steps = [s for s in steps if re.search(r"pip install[^\n]*(?:tar\.gz|SDIST)", s["run"])]
    assert len(steps) == 1, "Exactly one required sdist install/version gate must precede upload"
    return steps[0]


def test_sdist_uses_a_separate_clean_environment_and_both_version_checks():
    step = _sdist_step()
    body = step["run"]
    assert "-m venv" in body and "--system-site-packages" not in body
    assert "wheel-venv" not in body, "Do not reuse the installed wheel's environment"
    assert "--version" in body and "importlib.metadata" in body
    assert "VERSION" in body and "cd " in body
    assert not step.get("continue-on-error") and "always()" not in str(step.get("if", ""))


def test_publisher_still_uploads_the_verified_artifacts_without_rebuild_or_tests():
    _sdist_step()
    workflow = _workflow()
    upload = workflow["jobs"]["publish"]
    assert upload["needs"] == "verify"
    assert any(s.get("with", {}).get("name") == "dist" for s in upload["steps"])
    all_runs = "\n".join(s.get("run", "") for s in _gate_steps())
    assert len(re.findall(r"^\s*python(?:3)? -m build\s*$", all_runs, re.M)) == 1
    assert not re.search(r"\bpytest\b|\btox\b", all_runs)
    assert not any("-m build" in s.get("run", "") for s in upload["steps"])


@pytest.mark.parametrize("failure", ["none", "install", "cli", "metadata"])
def test_sdist_shell_gate_fails_closed_and_installs_the_actual_dist(tmp_path, failure):
    body = _sdist_step()["run"]
    checkout = tmp_path / "checkout"
    checkout.mkdir()
    (checkout / "dist").mkdir()
    sdist = checkout / "dist/allelix-2.7.0.tar.gz"
    sdist.write_bytes(b"synthetic artifact; installer is stubbed")
    runner = tmp_path / "runner"
    runner.mkdir()
    binary = tmp_path / "bin"
    binary.mkdir()
    # No network or actual installation: execute the workflow's real shell
    # with stand-ins that trace the artifact and expose each failure mode.
    stub = (
        f"#!{sys.executable}\n"
        + r"""
import json, os, pathlib, sys
args = sys.argv[1:]
log = pathlib.Path(os.environ["PROBE_LOG"])
with log.open("a") as out:
    entry = {"program": pathlib.Path(sys.argv[0]).name, "args": args, "cwd": os.getcwd()}
    out.write(json.dumps(entry) + "\n")
if args[:2] == ["-m", "venv"]:
    bindir = pathlib.Path(args[-1]) / "bin"
    bindir.mkdir(parents=True)
    for name in ("python", "python3", "allelix"):
        dest = bindir / name
        dest.write_text(pathlib.Path(sys.argv[0]).read_text())
        dest.chmod(0o755)
elif "install" in args:
    assert pathlib.Path(args[-1]).resolve() == pathlib.Path(os.environ["SDIST_EXPECTED"])
    if os.environ["FAILURE"] == "install": sys.exit(7)
elif args == ["--version"] or ("-c" in args):
    assert os.getcwd() != os.environ["SOURCE_ROOT"], "version checks must run outside checkout"
    kind = "cli" if args == ["--version"] else "metadata"
    version = "0.0.0" if os.environ["FAILURE"] == kind else os.environ["VERSION"]
    print("allelix, version " + version if args == ["--version"] else version)
else:
    raise AssertionError(args)
"""
    )
    for name in ("python", "python3"):
        path = binary / name
        path.write_text(stub)
        path.chmod(0o755)
    env = dict(
        os.environ,
        PATH=str(binary) + os.pathsep + os.environ["PATH"],
        RUNNER_TEMP=str(runner),
        VERSION="2.7.0",
        FAILURE=failure,
        PROBE_LOG=str(tmp_path / "trace.jsonl"),
        SOURCE_ROOT=str(checkout),
        SDIST_EXPECTED=str(sdist.resolve()),
    )
    run = subprocess.run(
        ["bash", "-euo", "pipefail", "-c", body],
        cwd=checkout,
        env=env,
        capture_output=True,
        text=True,
        timeout=15,
    )
    assert (run.returncode == 0) == (failure == "none"), (run.stdout, run.stderr)

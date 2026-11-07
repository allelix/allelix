# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""A db update stopped by SIGTERM or SIGHUP removes its staging entries.

SIGTERM and SIGHUP end the command through its normal exit, so its own
cleanup runs. A signal the caller ignores (SIGHUP under nohup) or handles
stays as it was, and only the handlers installed here are changed. SIGKILL
leaves the entry, which `db clean` removes (it deletes only allelix's names).

The update runs in a child process (`db update --only clinvar --force`) with
the network replaced and its first download blocked, so each signal arrives
while the update's staging directory exists. The last test is review's
(``test_termination_restores_a_caller_owned_other_signal_handler``), unchanged
except that it imports the helper inside the test, so the file collects on a
tree without it.
"""

from __future__ import annotations

import os
import signal
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest
from click.testing import CliRunner

from allelix.cli import main

pytestmark = pytest.mark.skipif(os.name != "posix", reason="POSIX signals and flock")

_REPO = Path(__file__).resolve().parents[1]

_CHILD = textwrap.dedent(
    """
    import sys, time
    from pathlib import Path
    from allelix.annotators import clinvar as cv
    from allelix.annotators.clinvar import ClinVarAnnotator

    def blocking(url, dest, *args, **kwargs):
        Path(dest).write_bytes(b"x" * 4096)
        print("READY", flush=True)
        time.sleep(120)

    cv.download = blocking
    cv._fetch_checksum = lambda url: "0" * 32
    ClinVarAnnotator.remote_signal = lambda self: "md5:" + "0" * 32
    from allelix.cli import main
    data = sys.argv[1]
    sys.argv = ["allelix", "db", "update", "--data-dir", data, "--only", "clinvar", "--force"]
    main()
    """
)


def _staging(data: Path) -> list[str]:
    return sorted(entry.name for entry in data.iterdir() if entry.name.startswith("."))


class _Update:
    """A blocked `db update` in a child process."""

    def __init__(self, tmp_path: Path, data: Path, *, ignore_hup: bool = False) -> None:
        script = tmp_path / "child.py"
        script.write_text(_CHILD)
        env = {**os.environ, "PYTHONPATH": str(_REPO)}

        def _preexec() -> None:
            if ignore_hup:
                signal.signal(signal.SIGHUP, signal.SIG_IGN)

        self.process = subprocess.Popen(
            [sys.executable, str(script), str(data)],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            env=env,
            preexec_fn=_preexec,
        )
        lines = []
        assert self.process.stdout is not None
        for line in self.process.stdout:
            lines.append(line)
            if "READY" in line:
                break
        else:
            self.process.wait(timeout=30)
            pytest.fail("the update never reached its download:\n" + "".join(lines))

    def kill(self, signum: int) -> int:
        self.process.send_signal(signum)
        return self.process.wait(timeout=30)

    def stop(self) -> None:
        if self.process.poll() is None:
            self.process.kill()
            self.process.wait(timeout=30)
        if self.process.stdout is not None:
            self.process.stdout.close()


@pytest.fixture
def update(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    started: list[_Update] = []

    def start(**kwargs) -> _Update:
        child = _Update(tmp_path, data, **kwargs)
        started.append(child)
        return child

    yield data, start
    for child in started:
        child.stop()


@pytest.mark.parametrize("name", ["SIGTERM", "SIGHUP"])
def test_sigterm_and_sighup_run_the_update_cleanup(update, name):
    data, start = update
    child = start()
    assert _staging(data), "the update made no staging entry"

    signum = getattr(signal, name)
    assert child.kill(signum) == 128 + signum
    assert _staging(data) == []


def test_ignored_sighup_stays_ignored(update):
    """Under nohup the update goes on after SIGHUP; SIGTERM still cleans up."""
    data, start = update
    child = start(ignore_hup=True)

    child.process.send_signal(signal.SIGHUP)
    with pytest.raises(subprocess.TimeoutExpired):
        child.process.wait(timeout=1)
    assert _staging(data)
    assert child.kill(signal.SIGTERM) == 128 + signal.SIGTERM
    assert _staging(data) == []


def test_db_clean_after_sigkill_empties_the_data_directory(update):
    data, start = update
    child = start()
    child.kill(signal.SIGKILL)

    result = CliRunner().invoke(main, ["db", "clean", "--data-dir", str(data), "--yes"])

    assert result.exit_code == 0, result.output
    assert sorted(entry.name for entry in data.iterdir()) == ["config.toml"]


def test_termination_restores_a_caller_owned_other_signal_handler():
    if not hasattr(signal, "SIGHUP"):
        return  # SIGHUP is not a Windows signal.
    old_hup = signal.getsignal(signal.SIGHUP)
    old_term = signal.getsignal(signal.SIGTERM)

    def caller_handler(signum, frame):
        pass

    try:
        signal.signal(signal.SIGHUP, caller_handler)
        signal.signal(signal.SIGTERM, signal.SIG_DFL)
        try:
            from allelix.cli import _signals_end_normally

            with _signals_end_normally():
                signal.raise_signal(signal.SIGTERM)
        except SystemExit as exc:
            assert exc.code == 128 + signal.SIGTERM
        else:
            raise AssertionError("SIGTERM did not end the command")
        assert signal.getsignal(signal.SIGHUP) is caller_handler, (
            "a handler the command did not install was permanently replaced with SIG_IGN"
        )
        assert signal.getsignal(signal.SIGTERM) is signal.SIG_DFL
    finally:
        signal.signal(signal.SIGHUP, old_hup)
        signal.signal(signal.SIGTERM, old_term)

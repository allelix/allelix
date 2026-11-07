# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Allelix command-line interface.

The CLI is organized as a package: `main` (the click group) is defined
here, and each subcommand module decorates against it. Importing this
package triggers all subcommand registrations as a side effect.

Public API: ``from allelix.cli import main`` — this is what the package
entry point ``allelix = "allelix.cli:main"`` calls.
"""

from __future__ import annotations

import codecs
import contextlib
import signal
import sys
import threading
from typing import TYPE_CHECKING

import click

from allelix import __version__
from allelix.config import ConfigFileError
from allelix.parsers import DamagedFileError, InputContainerError
from allelix.reports import UnwritableDestinationError

if TYPE_CHECKING:
    from collections.abc import Iterator
    from types import FrameType

# Signals that end a command through its normal exit, so every staging entry
# and partial output it owns is removed on the way out (2.6.1).
_ENDING_SIGNALS = tuple(
    getattr(signal, name) for name in ("SIGTERM", "SIGHUP") if hasattr(signal, name)
)


def _end_normally(signum: int, frame: FrameType | None) -> None:
    # A second signal must not cut the cleanup short; SIGKILL still can. Only
    # the handlers installed here are changed: a caller's own stays as it was.
    for ending in _ENDING_SIGNALS:
        if signal.getsignal(ending) is _end_normally:
            signal.signal(ending, signal.SIG_IGN)
    raise SystemExit(128 + signum)


@contextlib.contextmanager
def _signals_end_normally() -> Iterator[None]:
    """Turn SIGTERM and SIGHUP into SystemExit while a command runs.

    Python's default action ends the process without running ``finally``
    blocks. A signal the caller ignores (SIGHUP under ``nohup``) or handles
    stays as it was. The previous handlers are restored afterwards.
    """
    if threading.current_thread() is not threading.main_thread():
        yield
        return
    previous = {}
    for ending in _ENDING_SIGNALS:
        if signal.getsignal(ending) is signal.SIG_DFL:
            previous[ending] = signal.signal(ending, _end_normally)
    try:
        yield
    finally:
        for ending, handler in previous.items():
            signal.signal(ending, handler)


def _replace_unencodable_output(stream: object) -> None:
    """Make ``stream`` replace a character its encoding cannot hold, instead of raising (2.6.1).

    For a console, pipe or redirect that is not UTF-8 (Windows' cp1252, an
    ASCII locale): `analyze` crashed with UnicodeEncodeError on a character
    such as the plausibility marker "○". A UTF-8 stream, one with a
    non-strict error mode (``PYTHONIOENCODING=cp1252:ignore`` or
    ``:backslashreplace``), and one that cannot be reconfigured are left as
    they are; a strict one, given (``cp1252:strict``) or by default, is
    switched. Files are written as UTF-8.
    """
    encoding = getattr(stream, "encoding", None)
    reconfigure = getattr(stream, "reconfigure", None)
    if not encoding or reconfigure is None or getattr(stream, "errors", None) != "strict":
        return
    try:
        if codecs.lookup(encoding).name == "utf-8":
            return
    except LookupError:
        return
    reconfigure(errors="replace")


class _Group(click.Group):
    """Click group that reports a damaged or unreadable file as a plain error.

    A cut or corrupt gzip file is found while a command reads it, which can
    be deep in any command. Converting it here gives every command the same
    "Error: <file> looks truncated or damaged ..." and exit status 1. So do
    an unreadable config.toml and an output that cannot be written (2.5.2).
    SIGTERM and SIGHUP end the command as a normal exit would, with status
    128 + the signal number (2.6.1).
    """

    def invoke(self, ctx: click.Context) -> object:
        for stream in (sys.stdout, sys.stderr):
            _replace_unencodable_output(stream)
        try:
            with _signals_end_normally():
                return super().invoke(ctx)
        except (
            DamagedFileError,
            InputContainerError,
            ConfigFileError,
            UnwritableDestinationError,
        ) as exc:
            raise click.ClickException(str(exc)) from exc


@click.group(cls=_Group)
@click.version_option(version=__version__, prog_name="allelix")
def main() -> None:
    """Allelix: open-source genotype analysis toolkit."""


# Side-effect imports: each module decorates against `main` at import
# time. Order matters only for --help listing, not correctness.
from allelix.cli import (  # noqa: E402,F401
    analyze,
    config,
    db,
    focused,
    utility,
)

__all__ = ["main"]

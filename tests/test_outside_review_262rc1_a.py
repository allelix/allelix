# Outside review of Allelix 2.6.2 rc1 (focus A: integrity and operations).
#
# One test per confirmed finding. Each test FAILS on rc1 at the reported
# behaviour; the comment on each says how it behaves on base (2.6.1).
#
# Run (rc1):
#   PYTHONPATH=/root/outside-262/A/rc1:/root/outside-262/A/rc1/tests \
#     venv-rc1/bin/python -m pytest out/test_review_rc1.py -q -p no:cacheprovider
# Run (base), for the base column of the comments:
#   PYTHONPATH=/root/outside-262/A/base:/root/outside-262/A/base/tests \
#     venv-base/bin/python -m pytest out/test_review_rc1.py -q -p no:cacheprovider
#
# Synthetic data only; no network (subprocesses run with --no-update).

from __future__ import annotations

import contextlib
import os
import re
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest
from click.testing import CliRunner

import allelix
from allelix.annotators.clinvar import clinvar_db_filename
from allelix.cli import main
from test_annotation_parity_and_named_disclosures import _mthfr, _vcf

# The source tree of the allelix under test (rc1/ or base/).
REPO = Path(allelix.__file__).resolve().parents[1]

# A POSIX ASCII locale as the 2.6.2 CHANGELOG means it: stdout ascii/surrogateescape.
C_LOCALE = {"LC_ALL": "C", "LANG": "C", "PYTHONUTF8": "0", "PYTHONCOERCECLOCALE": "0"}
CLINVAR_ONLY = (
    "[sources]\nclinvar = true\npharmgkb = false\nsnpedia = false\ngwas = false\n"
    "gnomad = false\nalphamissense = false\ncadd = false\nrsid_identity = true\n"
)
CHILD = (
    "import sys; from allelix.cli import main; "
    "print(f'streams: {sys.stdout.encoding}/{sys.stdout.errors}', file=sys.stderr); "
    "sys.argv = ['allelix', *sys.argv[1:]]; main()"
)


def _child_env(**extra: str) -> dict[str, str]:
    env = {
        k: v
        for k, v in os.environ.items()
        if not k.startswith(("LC_", "LANG", "PYTHONIO", "PYTHONUTF8", "PYTHONCOERCE"))
    }
    env.update(extra, PYTHONPATH=str(REPO))
    return env


def _allelix(*args: str, env: dict[str, str]) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        [sys.executable, "-c", CHILD, *args],
        capture_output=True,
        env=env,
        timeout=180,
        check=False,
    )


def _data(tmp_path: Path, *, config: str | None = None) -> tuple[Path, int]:
    data = tmp_path / "data"
    pos = _mthfr(data, "GRCh37")
    if config is not None:
        (data / "config.toml").write_text(config)
    return data, pos


# ---------------------------------------------------------------------------
# A-01 (P3, new in 2.6.2): rich.markup.escape() doubles a trailing backslash,
# so database text that ends in "\" is shown with two backslashes wherever the
# escaped text is not followed by a closing tag (table cells, "(version ...)").
# base (2.6.1): PASSES (the text was not escaped and printed as written).
# rc1: FAILS ("GENE\\", "Condition ends \\", "Trait ends \\").
# ---------------------------------------------------------------------------
def test_a01_database_text_ending_in_a_backslash_is_shown_as_written(tmp_path):
    data, pos = _data(tmp_path)
    with contextlib.closing(sqlite3.connect(data / "gwas.sqlite")) as conn:
        conn.execute("UPDATE gwas_associations SET trait = 'Trait ends \\'")
        conn.commit()
    with contextlib.closing(sqlite3.connect(data / clinvar_db_filename("GRCh37"))) as conn:
        conn.execute("UPDATE clinvar_variants SET condition = 'Condition ends \\', gene = 'GENE\\'")
        conn.commit()
    vcf = _vcf(tmp_path / "input.vcf", [("rs1801133", "1", pos, "G", "A", "0/1", ".")])
    args = [
        "analyze", str(vcf), "--build", "grch37", "--data-dir", str(data), "--no-update",
        "--min-magnitude", "0", "--gwas-min-magnitude", "0",
    ]
    result = CliRunner(env={"COLUMNS": "300"}).invoke(main, args)
    assert result.exit_code == 0, result.output
    out = result.output
    assert "Condition ends \\" in out and "Trait ends \\" in out and "GENE\\" in out
    assert "\\\\" not in out, [line for line in out.splitlines() if "\\\\" in line]


# ---------------------------------------------------------------------------
# A-03 (P3, new in 2.6.2, documentation): the CHANGELOG gives the trigger as
# "a POSIX ASCII locale (LC_ALL=C, with PYTHONIOENCODING unset)". On the
# supported Pythons (>=3.12) LC_ALL=C alone turns on UTF-8 mode (PEP 540), so
# stdout is utf-8/surrogateescape and 2.6.1 did not crash; the crash needs
# PYTHONUTF8=0 as well (as rc1's own test sets). The entry does not say so.
# base: FAILS (no 2.6.2 section). rc1: FAILS on the CHANGELOG assertion.
# ---------------------------------------------------------------------------
def test_a03_changelog_names_the_setting_that_makes_the_console_ascii():
    probe = subprocess.run(
        [sys.executable, "-c", "import sys; print(sys.stdout.encoding, sys.stdout.errors)"],
        capture_output=True,
        env=_child_env(LC_ALL="C"),
        check=True,
    )
    assert probe.stdout.split() == [b"utf-8", b"surrogateescape"]  # LC_ALL=C alone is not ASCII
    text = (REPO / "CHANGELOG.md").read_text(encoding="utf-8")
    section = text.split("## [2.6.2]", 1)[1].split("\n## [", 1)[0]
    entry = next(p for p in section.split("\n- ") if "UnicodeEncodeError" in p)
    assert "PYTHONUTF8" in entry or "-X utf8" in entry, entry

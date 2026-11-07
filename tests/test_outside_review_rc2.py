"""Outside review A of allelix 2.6.1 RC2: re-check of the RC1 fixes and the RC1->RC2 delta.

Each test fails on RC2 for the reason in its docstring. Offline and synthetic:
downloads, checksums and freshness probes are replaced, and the "killed"
updates are real `db update` runs in a child process that SIGKILLs itself at
the named phase (no timing involved). Run with:

    venv-rc2/bin/python -m pytest -c /dev/null -p no:cacheprovider out/test_review_rc2.py
"""

from __future__ import annotations

import io
import json
import os
import subprocess
import sys
import textwrap
import urllib.request
from email.message import Message
from typing import TYPE_CHECKING

import pytest
from click.testing import CliRunner

from allelix.cli import main
from allelix.databases import manager

if TYPE_CHECKING:
    from pathlib import Path

pytestmark = pytest.mark.skipif(os.name != "posix", reason="POSIX SIGKILL")

# A child `db update` of one source that SIGKILLs itself at the named phase.
_CHILD = textwrap.dedent(
    r"""
    import gzip, hashlib, importlib, os, signal, sys, zipfile
    from pathlib import Path
    import allelix

    source, data = sys.argv[1], sys.argv[2]
    repo = Path(allelix.__file__).resolve().parents[1]

    def killed(*args, **kwargs):
        os.kill(os.getpid(), signal.SIGKILL)

    db_module = importlib.import_module("allelix.cli.db")
    args = []
    if source == "gwas":
        from allelix.annotators import gwas as module
        from allelix.annotators.gwas import GWASCatalogAnnotator

        table = repo / "tests" / "fixtures" / "mock_gwas_catalog.tsv"

        def download(url, dest):
            with zipfile.ZipFile(dest, "w") as archive:
                # EBI's member name, as in the real zip.
                archive.write(table, arcname="gwas-catalog-associations_ontology-annotated.tsv")

        module.download = download
        module.GWAS_MIN_ROWS = 0
        module.head_request_headers = lambda url, *a, **k: {"ETag": '"v1"'}
        module.load_gwas_tsv = killed  # killed while the table loads
        db_module.get_annotators = lambda d, **kw: [GWASCatalogAnnotator(d)]
    else:
        from allelix.annotators import clinvar as cv
        from allelix.databases import manager

        def tsv(header_map, width, first, rows):
            header = [header_map.get(i, f"unused{i}") for i in range(width)]
            header[0] = first
            lines = ["\t".join(header)]
            for values in rows:
                row = [""] * width
                for index, value in values.items():
                    row[index] = value
                lines.append("\t".join(row))
            return gzip.compress(("\n".join(lines) + "\n").encode(), mtime=0)

        vs = tsv(manager._VS_REQUIRED_COLUMNS, 43, "#AlleleID", [
            {0: "101", 1: "single nucleotide variant", 2: "NM_X:c.1A>G", 4: "GENE",
             6: "Pathogenic", 9: "900253001", 16: "GRCh37", 18: "1",
             24: "criteria provided, single submitter", 30: "101", 31: "1000",
             32: "A", 33: "G"},
        ])
        ss = tsv(manager._SS_REQUIRED_COLUMNS, 16, "#VariationID", [
            {0: "101", 1: "Pathogenic", 5: "na:Synthetic",
             6: "criteria provided, single submitter",
             10: "SCV900252901.1", 15: "yes"},
        ])
        payloads = {"variant_summary.txt.gz": vs, "submission_summary.txt.gz": ss}
        digests = {name: hashlib.md5(body).hexdigest() for name, body in payloads.items()}
        cv._fetch_checksum = lambda url: digests[url.rsplit("/", 1)[-1]]
        cv.download = lambda url, dest: Path(dest).write_bytes(payloads[url.rsplit("/", 1)[-1]])
        manager.CLINVAR_MIN_ROWS = {"GRCh37": 1, "GRCh38": 1}
        real_groups = manager._variation_groups

        def groups(path, counts=None):
            for item in real_groups(path, counts):
                killed()  # killed while the cache builds
                yield item

        manager._variation_groups = groups
        args = ["--only", "clinvar", "--build", "grch37"]

    from allelix.cli import main
    sys.argv = ["allelix", "db", "update", "--data-dir", data, *args]
    main()
    """
)


def _killed_update(tmp_path: Path, source: str) -> Path:
    data = tmp_path / "data"
    data.mkdir()
    script = tmp_path / "child.py"
    script.write_text(_CHILD)
    done = subprocess.run(
        [sys.executable, str(script), source, str(data)],
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )
    assert done.returncode == -9, done.stdout + done.stderr  # it was killed, as planned
    return data


def _staging_dirs(data: Path) -> list[Path]:
    return [p for p in data.iterdir() if p.is_dir() and p.name.startswith(".")]


def test_db_clean_removes_what_a_killed_gwas_update_leaves(tmp_path):
    """RC2 `_is_owned_staging_directory` (allelix/cli/db.py) keeps a staging dir
    unless every entry is a download name. Its own comment lists "the GWAS table
    once moved out of `extracted/`", but `GWASCatalogAnnotator.setup`
    (annotators/gwas.py:286-287) moves the table out with `rename` and leaves
    the empty `extracted/` directory, which is not a download name. So every
    GWAS update killed after extraction (the load is the long phase) leaves
    `.gwas-*/{gwas_catalog_associations.zip, gwas_catalog_associations.tsv,
    extracted/}`: nothing but allelix's own names, and `db clean` keeps it
    ("delete it by hand"). RC1 removed it. CHANGELOG: "One killed by SIGKILL or
    a crash still leaves them; `db clean` removes them."
    """
    data = _killed_update(tmp_path, "gwas")
    left = _staging_dirs(data)
    assert len(left) == 1
    assert sorted(p.name for p in left[0].iterdir()) == [
        "extracted",
        "gwas_catalog_associations.tsv",
        "gwas_catalog_associations.zip",
    ]
    assert list((left[0] / "extracted").iterdir()) == []

    result = CliRunner().invoke(main, ["db", "clean", "--data-dir", str(data), "--yes"])

    assert result.exit_code == 0, result.output
    assert _staging_dirs(data) == [], " ".join(result.output.split())


def test_db_clean_removes_what_a_killed_clinvar_build_leaves(tmp_path):
    """A ClinVar update killed while it builds (the long phase, after ~700 MB
    of downloads) leaves `.clinvar-*/{variant_summary.txt.gz,
    submission_summary.txt.gz, .clinvar.GRCh37.sqlite.<16 hex>.partial}`. The
    partial file has `owned_staging`'s random-token name, which no user backup
    has. RC2 keeps the whole directory and says "delete it by hand if so".
    RC1 removed it. ADR-0056 (RC2 amendment) now says such a directory is
    kept, but the CHANGELOG ("`db clean` removes them", and for Windows
    "`db clean` removes what they leave"), the README and `db clean --help`
    ("staging entries a killed update left") still say it is removed.
    """
    data = _killed_update(tmp_path, "clinvar")
    left = _staging_dirs(data)
    assert len(left) == 1
    names = sorted(p.name for p in left[0].iterdir())
    assert "variant_summary.txt.gz" in names and "submission_summary.txt.gz" in names
    assert any(n.startswith(".clinvar.GRCh37.sqlite.") and n.endswith(".partial") for n in names)

    result = CliRunner().invoke(main, ["db", "clean", "--data-dir", str(data), "--yes"])

    assert result.exit_code == 0, result.output
    assert _staging_dirs(data) == [], " ".join(result.output.split())


class _Reply(io.BytesIO):
    def __init__(self, body: object, content_range: str) -> None:
        super().__init__(json.dumps(body).encode())
        self.headers = Message()
        self.headers["Content-Range"] = content_range

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


@pytest.mark.parametrize(
    "row",
    [
        {"definitionid": 1, "clinicalfunctionalstatus": 3},  # .strip() on an int
        {"id": [1], "dbsnpid": "rs1"},  # an unhashable key
    ],
)
def test_cpic_row_with_a_wrong_typed_value_takes_the_cpic_unavailable_path(monkeypatch, row):
    """RC2 closes F6 for rows that are not objects. A row that is an object
    with a value of the wrong type still passes `_check_whole_table` and then
    raises AttributeError or TypeError in `fetch_cpic_allele_functions`
    (cpic_loader.py:255-275), outside PROBE_ERRORS, so the ClinPGx update
    fails instead of being built without CPIC rows (ADR-0020 amendment).
    Low likelihood: the same class as F6, not yet closed.
    """
    from allelix.databases import cpic_loader

    monkeypatch.setattr(cpic_loader, "CPIC_RETRY_BACKOFF_SECONDS", (0.0, 0.0, 0.0))
    monkeypatch.setattr(
        urllib.request, "urlopen", lambda request, timeout=None: _Reply([row], "0-0/1")
    )

    try:
        cpic_loader.fetch_cpic_allele_functions()
    except getattr(manager, "PROBE_ERRORS", (OSError, ValueError)):
        return  # the ClinPGx setup's fallback path catches it
    except Exception as exc:  # the point of the test
        pytest.fail(f"{type(exc).__name__} escapes the CPIC-unavailable path: {exc}")

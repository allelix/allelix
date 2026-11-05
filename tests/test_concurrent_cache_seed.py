"""Deterministic offline reproducer: two refreshes share gwas.sqlite.tmp.

Only the record provider and scheduling at os.replace are controlled. The
production cache writer, schema, commits, row floor and cleanup all run.
No network or real genomic data is used.
"""

from __future__ import annotations

import contextlib
import sqlite3
import threading
from typing import TYPE_CHECKING
from unittest.mock import patch

from allelix.databases import gwas_loader

if TYPE_CHECKING:
    from pathlib import Path


def reproduce(root: Path) -> dict:
    root.mkdir(parents=True, exist_ok=True)
    db = root / "gwas.sqlite"
    events = {name: threading.Event() for name in ("a_ready", "b_open", "a_done")}
    results: dict[str, object] = {}
    real_replace = gwas_loader.os.replace
    record = dict(
        rsid="rs900000001",
        risk_allele="A",
        trait="Synthetic trait",
        p_value=1e-9,
        or_beta=1.1,
        ci_text="",
        gene="TEST",
        study_accession="GCST000000",
        pubmed_id="1",
        risk_allele_frequency=0.1,
        context="",
        mapped_trait_uri="",
        trait_category="other",
    )

    def records(path: Path):
        if path.name == "b.tsv":
            # B has now unlinked A's committed staging file, created its own
            # database and schema, and is about to read its first record.
            events["b_open"].set()
            if not events["a_done"].wait(15):
                raise TimeoutError("A did not finish")
            # Model a failed parse/read after A has published B's staging DB.
            raise OSError("simulated B input read failure")
        yield record

    def scheduled_replace(src, dst):
        if threading.current_thread().name == "writer-a":
            events["a_ready"].set()
            if not events["b_open"].wait(15):
                raise TimeoutError("B did not open its staging database")
        return real_replace(src, dst)

    def writer_a():
        try:
            results["a_returned_count"] = gwas_loader.load_gwas_tsv(root / "a.tsv", db, min_rows=1)
        except Exception as exc:
            results["a_exception"] = f"{type(exc).__name__}: {exc}"
        finally:
            events["a_done"].set()

    def writer_b():
        try:
            results["b_returned_count"] = gwas_loader.load_gwas_tsv(root / "b.tsv", db, min_rows=1)
        except Exception as exc:
            results["b_exception"] = f"{type(exc).__name__}: {exc}"

    with patch.object(gwas_loader, "iter_gwas_records", records):
        # Start with a known-good cache; the failure must not destroy it.
        gwas_loader.load_gwas_tsv(root / "seed.tsv", db, min_rows=1)
        with patch.object(gwas_loader.os, "replace", scheduled_replace):
            a = threading.Thread(target=writer_a, name="writer-a")
            b = threading.Thread(target=writer_b, name="writer-b")
            a.start()
            if not events["a_ready"].wait(15):
                raise TimeoutError("A never reached its publish step")
            b.start()
            a.join(20)
            b.join(20)
            if a.is_alive() or b.is_alive():
                raise TimeoutError("A writer failed to finish")

    with contextlib.closing(sqlite3.connect(db)) as conn:
        results["installed_associations"] = conn.execute(
            "SELECT COUNT(*) FROM gwas_associations"
        ).fetchone()[0]
        results["installed_version_rows"] = conn.execute(
            "SELECT COUNT(*) FROM database_versions"
        ).fetchone()[0]
        results["integrity_check"] = conn.execute("PRAGMA integrity_check").fetchone()[0]
    return results


def test_successful_writer_never_publishes_the_failed_writers_empty_cache(tmp_path):
    result = reproduce(tmp_path / "race")
    assert result.get("a_returned_count") == 1, result
    assert result.get("b_exception", "").startswith("OSError: simulated B input read failure"), (
        result
    )
    assert result["installed_associations"] == 1, result
    assert result["installed_version_rows"] == 1, result
    assert result["integrity_check"] == "ok", result

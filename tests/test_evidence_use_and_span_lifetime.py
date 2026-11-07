# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Track the new evidence dependency and limit span caching to its connection.

No findings need be emitted for gnomAD to have changed the report. Closing
an annotator also ends its SQLite snapshot; a later cache must be measured
again. Synthetic, offline, with disabled/missing and unchanged-cache controls.
"""

import contextlib
import json
import sqlite3

import pytest
from click.testing import CliRunner
from test_allele_evidence_provider import _own_only
from test_outside_review_rc1_b import BUILDS, POS_A, RSA, RSB, _chip, _clinvar, _config

from allelix.annotators.clinvar import ClinVarAnnotator, clinvar_db_filename
from allelix.annotators.gnomad import GnomadAnnotator
from allelix.cli import main
from allelix.models import Variant


@pytest.mark.parametrize("kind", ["23andme", "ancestrydna"])
@pytest.mark.parametrize(
    "state", ["queried-hit", "queried-empty", "unqueried", "disabled", "missing", "flag"]
)
def test_lazy_evidence_query_is_reported_as_source_use(tmp_path, monkeypatch, kind, state):
    data = _own_only(tmp_path, "enabled")
    config = data / "config.toml"
    if state == "disabled":
        config.write_text(config.read_text().replace("gnomad = true", "gnomad = false"))
    elif state == "missing":
        (data / "gnomad.sqlite").unlink()
    elif state == "queried-empty":
        with contextlib.closing(sqlite3.connect(data / "gnomad.sqlite")) as conn:
            conn.execute("UPDATE gnomad_frequencies SET pos = pos + 100")
            conn.commit()
    asked = []
    original = GnomadAnnotator.bulk_resolve_coordinates

    def record_query(owner, rsids):
        asked.append(set(rsids))
        return original(owner, rsids)

    monkeypatch.setattr(GnomadAnnotator, "bulk_resolve_coordinates", record_query)
    call = "GG" if state == "unqueried" else "CC"
    source = _chip(tmp_path / "chip.txt", [(RSA, "1", POS_A, call)], build="GRCh38", kind=kind)
    out = tmp_path / "out.json"
    args = [
        "analyze",
        str(source),
        "--build",
        "grch38",
        "--data-dir",
        str(data),
        "--no-update",
        "--output",
        str(out),
    ]
    if state == "flag":
        args.append("--no-gnomad")
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, result.output
    report = json.loads(out.read_text())
    assert report["annotations"] == [], "this probes dependency-only use, not enrichment"
    assert report["input"]["probe_classes"]["foreign_allele_calls_at_plp_loci"] == int(
        state not in {"queried-hit", "unqueried"}
    )
    queried = state.startswith("queried-")
    assert asked == ([{RSA}] if queried else [])
    record = next(s for s in report["sources"] if s["name"] == "gnomad")
    expected_state = (
        "used" if queried else ("ready_not_used" if state == "unqueried" else "skipped")
    )
    assert record["state"] == expected_state, record
    if queried:
        assert record["version"] == "synthetic"
        assert [a["name"] for a in report["annotators"]].count("gnomad") == 1
    else:
        assert record["reason"] == (None if state == "unqueried" else state)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("path", ["single", "batch"])
@pytest.mark.parametrize("initial_ref", ["G", "GA"], ids=["no-indels", "short-indel"])
@pytest.mark.parametrize("refresh", [False, True], ids=["same-cache-control", "longer-new-cache"])
def test_span_cache_ends_when_its_connection_closes(tmp_path, build, path, initial_ref, refresh):
    data = tmp_path / "data"
    _config(data)
    _clinvar(
        data,
        build,
        [
            (
                RSB,
                "1",
                POS_A,
                initial_ref,
                "A" if initial_ref == "G" else "G",
                "GENEP",
                "Pathogenic",
                "Pathogenic",
            )
        ],
    )
    probe = Variant("i3999261199", "1", POS_A + 75, "D", "I", build=build)
    owner = ClinVarAnnotator(data, builds=(build,))

    def names():
        with owner:
            rows = owner.annotate(probe) if path == "single" else owner.batch_annotate([probe])
            assert list(rows) == []
            return set(owner.unchecked())

    assert names() == set()
    # A refresh occurs only AFTER the old connection has closed.
    if refresh:
        with contextlib.closing(sqlite3.connect(data / clinvar_db_filename(build))) as conn:
            conn.execute("UPDATE clinvar_variants SET ref = ?, alt = 'G'", ("G" + "A" * 79,))
            conn.commit()
    assert names() == ({(probe.rsid, "di_call")} if refresh else set())

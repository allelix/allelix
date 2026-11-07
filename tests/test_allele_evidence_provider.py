# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""gnomAD's own-rsID alleles count before any foreign or unknown decision.

ClinVar asks an enabled, ready gnomAD for the call's own rsID at its GRCh38
locus whenever the identity table's SNV rows and ClinVar's own rows do not
hold every called base, whether or not another rsID's row stands there, so
the counts and names agree with the full evidence. A gnomAD switched off,
or with no cache, is never consulted, and the call's own rsID's rows match
either way. Synthetic and offline.
"""

from __future__ import annotations

import contextlib
import sqlite3

import pytest
from test_outside_review_rc1_b import (
    POS_A,
    RSA,
    _analyze,
    _chip,
    _clinvar,
    _clinvar_findings,
    _config,
    _identity,
    _named,
)

from allelix.annotators.clinvar import ClinVarAnnotator
from allelix.databases._versions import GNOMAD_SCHEMA_VERSION
from allelix.databases.schema import GNOMAD_SCHEMA
from allelix.models import Variant


def _gnomad(data, rows):
    with contextlib.closing(sqlite3.connect(data / "gnomad.sqlite")) as conn:
        conn.executescript(GNOMAD_SCHEMA)
        conn.executemany(
            "INSERT INTO gnomad_frequencies (chrom,pos,ref,alt,rsid,af) VALUES (?,?,?,?,?,?)",
            rows,
        )
        conn.execute(
            "INSERT INTO database_versions (name,source_url,version,downloaded_at,"
            "record_count,local_version_tag) VALUES ('gnomad','test://synthetic',"
            "'synthetic','2001-01-01',?,?)",
            (len(rows), f"sv:{GNOMAD_SCHEMA_VERSION}"),
        )
        conn.commit()


def _own_only(tmp_path, gnomad_state):
    """Identity and ClinVar give rsA G>A Pathogenic only; gnomAD gives rsA G>C."""
    data = tmp_path / "data"
    _config(data)
    if gnomad_state == "enabled":
        config = data / "config.toml"
        config.write_text(config.read_text().replace("gnomad = false", "gnomad = true"))
    _clinvar(data, "GRCh38", [(RSA, "1", POS_A, "G", "A", "GENEA", "Pathogenic", "Pathogenic")])
    _identity(data, [(RSA, "1", POS_A, "G", "A", POS_A, "G", "A", 0)])
    if gnomad_state != "missing":
        _gnomad(data, [("1", POS_A, "G", "C", RSA, 0.01)])
    return data


@pytest.mark.parametrize("gnomad_state", ["enabled", "disabled", "missing"])
def test_no_other_rsid_row_still_takes_gnomads_evidence_before_counting(tmp_path, gnomad_state):
    # No other rsID's row: CC is outside {G, A} but inside {G, A, C} with an
    # enabled gnomAD, so it is neither counted nor named then. Off or missing,
    # gnomAD is not consulted and the call is counted and named as foreign.
    data = _own_only(tmp_path, gnomad_state)
    source = _chip(tmp_path / "chip.txt", [(RSA, "1", POS_A, "CC")], build="GRCh38")
    report, output = _analyze(tmp_path, data, source, "GRCh38")
    foreign = gnomad_state != "enabled"
    assert report["input"]["probe_classes"]["foreign_allele_calls_at_plp_loci"] == int(foreign)
    assert ((RSA, "foreign_alleles") in _named(report)) is foreign, (_named(report), output)
    assert _clinvar_findings(report) == [], output


@pytest.mark.parametrize("gnomad_state", ["enabled", "disabled", "missing"])
def test_own_rsid_finding_does_not_depend_on_the_provider(tmp_path, gnomad_state):
    data = _own_only(tmp_path, gnomad_state)
    source = _chip(tmp_path / "chip.txt", [(RSA, "1", POS_A, "AA")], build="GRCh38")
    report, output = _analyze(tmp_path, data, source, "GRCh38")
    assert _clinvar_findings(report) == [(RSA, "AA", "clinvar_pathogenic")], output


@pytest.mark.parametrize("path", ["single", "batch"])
@pytest.mark.parametrize("wired", [True, False])
def test_single_and_batch_ask_the_provider_alike(tmp_path, path, wired):
    data = _own_only(tmp_path, "missing")
    asked: list[str] = []

    def provider(v):
        asked.append(v.rsid)
        return frozenset("GC")

    call = Variant(RSA, "1", POS_A, "C", "C", build="GRCh38", ref="G")
    call.allele_evidence = frozenset("GA")
    with ClinVarAnnotator(data, builds=("GRCh38",)) as clinvar:
        if wired:
            clinvar.set_allele_evidence_provider(provider)
        list(clinvar.annotate(call) if path == "single" else clinvar.batch_annotate([call]))
        counts = clinvar.run_counts()
    assert asked == ([RSA] if wired else [])
    assert counts["foreign_allele_calls_at_plp_loci"] == int(not wired)


def test_grch37_never_asks_gnomad(tmp_path):
    # gnomAD's coordinates are GRCh38: a GRCh37 call keeps the identity and
    # ClinVar evidence only, and is counted as foreign.
    data = tmp_path / "data"
    _config(data)
    config = data / "config.toml"
    config.write_text(config.read_text().replace("gnomad = false", "gnomad = true"))
    _clinvar(data, "GRCh37", [(RSA, "1", POS_A, "G", "A", "GENEA", "Pathogenic", "Pathogenic")])
    _identity(data, [(RSA, "1", POS_A, "G", "A", POS_A, "G", "A", 0)])
    _gnomad(data, [("1", POS_A, "G", "C", RSA, 0.01)])
    source = _chip(tmp_path / "chip.txt", [(RSA, "1", POS_A, "CC")], build="GRCh37")
    report, output = _analyze(tmp_path, data, source, "GRCh37")
    assert report["input"]["probe_classes"]["foreign_allele_calls_at_plp_loci"] == 1, output

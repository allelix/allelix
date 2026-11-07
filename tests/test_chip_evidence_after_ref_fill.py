# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""An earlier REF fill must not hide enabled own-rsID gnomAD evidence.

The same chip call can match the other rsID's allele when its own rsID's
SNV evidence includes it. Adding an older identity row must not remove that
evidence. Disabled or missing gnomAD supplies no evidence, even if its
cache exists. Everything here is synthetic and offline.
"""

import contextlib
import sqlite3

import pytest
from test_outside_review_rc1_b import (
    POS_A,
    RSA,
    RSB,
    _analyze,
    _chip,
    _clinvar,
    _clinvar_findings,
    _config,
    _identity,
)

from allelix.databases._versions import GNOMAD_SCHEMA_VERSION
from allelix.databases.schema import GNOMAD_SCHEMA


@pytest.mark.parametrize("identity_case", ["absent", "snv", "indel", "disabled"])
@pytest.mark.parametrize("gnomad_state", ["enabled", "disabled", "missing"])
@pytest.mark.parametrize("kind", ["23andme", "ancestrydna"])
@pytest.mark.parametrize("call", ["CC", "GC"])
def test_enabled_own_rsid_evidence_survives_an_earlier_ref_fill(
    tmp_path, identity_case, gnomad_state, kind, call
):
    data = tmp_path / "data"
    _config(data, identity=identity_case != "disabled")
    config = data / "config.toml"
    if gnomad_state != "disabled":
        config.write_text(config.read_text().replace("gnomad = false", "gnomad = true"))
    _clinvar(
        data,
        "GRCh38",
        [(RSB, "1", POS_A, "G", "C", "OTHER", "Pathogenic", "Pathogenic")],
    )
    alt = "GA" if identity_case == "indel" else "A"
    rows = [] if identity_case == "absent" else [(RSA, "1", POS_A, "G", alt, POS_A, "G", alt, 0)]
    _identity(data, rows)
    if gnomad_state != "missing":
        with contextlib.closing(sqlite3.connect(data / "gnomad.sqlite")) as conn:
            conn.executescript(GNOMAD_SCHEMA)
            conn.execute(
                "INSERT INTO gnomad_frequencies (chrom,pos,ref,alt,rsid,af) VALUES (?,?,?,?,?,?)",
                ("1", POS_A, "G", "C", RSA, 0.001),
            )
            conn.execute(
                "INSERT INTO database_versions (name,source_url,version,downloaded_at,"
                "record_count,local_version_tag) VALUES ('gnomad','test://synthetic',"
                "'synthetic','2001-01-01',1,?)",
                (f"sv:{GNOMAD_SCHEMA_VERSION}",),
            )
            conn.commit()
    source = _chip(tmp_path / "chip.txt", [(RSA, "1", POS_A, call)], build="GRCh38", kind=kind)
    report, output = _analyze(tmp_path, data, source, "GRCh38")
    findings = _clinvar_findings(report)
    assert bool(findings) is (gnomad_state == "enabled"), (findings, output)
    if gnomad_state == "enabled":
        assert len(findings) == 1
        assert findings[0][2] == "clinvar_pathogenic"
        counts = report["input"]["probe_classes"]
        assert counts["foreign_allele_calls_at_plp_loci"] == 0
        assert counts["unknown_allele_calls_at_plp_loci"] == 0
        assert not [e for e in report.get("high_value_unchecked") or [] if e["rsid"] == RSA]

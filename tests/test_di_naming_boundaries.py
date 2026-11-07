# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Pin the REF-span boundary and the aggregate-P/LP exclusion.

Naming does not interpret D/I. Only a call already named by its own
aggregate-P/LP rsID bypasses positional disclosure; a Benign row is not
such a disclosure. Exercise the fixed window and the separate long rows.
"""

import pytest
from test_outside_review_rc1_b import (
    BUILDS,
    _analyze,
    _chip,
    _clinvar,
    _clinvar_findings,
    _config,
    _identity,
    _named,
)

from allelix.annotators.clinvar import ClinVarAnnotator
from allelix.models import Variant

POS = 2_100_000
TARGET = "rs3999485001"
CALL = "rs3999485002"
PROBE = "i3999485003"


def _observe(tmp_path, data, build, path, call_id, pos):
    if path == "cli":
        source = _chip(tmp_path / "chip.txt", [(call_id, "13", pos, "DI")], build=build)
        report, _ = _analyze(tmp_path, data, source, build)
        assert not _clinvar_findings(report), "disclosure must not become a finding"
        return _named(report)
    call = Variant(call_id, "13", pos, "D", "I", build=build)
    with ClinVarAnnotator(data, builds=(build,)) as clinvar:
        findings = clinvar.annotate(call) if path == "single" else clinvar.batch_annotate([call])
        assert not list(findings), "disclosure must not become a finding"
        names = list(clinvar.unchecked())
        # The existing by-rsID channel is separate from positional disclosure.
        names.extend((rsid, "di_call") for rsid in clinvar.di_calls_at_plp())
        return names


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("path", ["single", "batch", "cli"])
@pytest.mark.parametrize("call_id", [PROBE, CALL], ids=["vendor", "unknown-rsid"])
@pytest.mark.parametrize(
    "ref,alt",
    [
        pytest.param("G", "G" + "A" * 100, id="long-insertion"),
        pytest.param("GA", "G", id="short-deletion"),
        pytest.param("G" + "A" * 63, "G", id="64-base-deletion"),
        pytest.param("G" + "A" * 64, "G", id="65-base-deletion"),
    ],
)
@pytest.mark.parametrize(
    "distance,named", [(5, True), (6, False)], ids=["last-slack", "past-slack"]
)
def test_ref_end_has_five_bases_of_slack(
    tmp_path, build, path, call_id, ref, alt, distance, named
):
    data = tmp_path / "data"
    _config(data)
    _identity(data, [])
    _clinvar(data, build, [(TARGET, "13", POS, ref, alt, "GENEP", "Pathogenic", "Pathogenic")])
    # The last REF base is inclusive. Insertions have only their anchor.
    pos = POS + len(ref) - 1 + distance
    names = _observe(tmp_path, data, build, path, call_id, pos)
    assert names.count((call_id, "di_call")) == int(named), names


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("path", ["single", "batch", "cli"])
@pytest.mark.parametrize("aggregate", ["Benign", "Uncertain significance", "Pathogenic"])
@pytest.mark.parametrize("target_ref", ["GA", "G" + "A" * 79], ids=["short", "long"])
def test_only_own_plp_disclosure_suppresses_positional_naming(
    tmp_path, build, path, aggregate, target_ref
):
    data = tmp_path / "data"
    _config(data)
    _identity(data, [])
    # A cache row under the call's ID on another chromosome is not enough:
    # it must itself be aggregate-P/LP to provide the existing rsID notice.
    _clinvar(
        data,
        build,
        [
            (CALL, "2", 800_000, "G", "A", "GENEOWN", aggregate, aggregate),
            (TARGET, "13", POS, target_ref, "G", "GENEP", "Pathogenic", "Pathogenic"),
        ],
    )
    names = _observe(tmp_path, data, build, path, CALL, POS + len(target_ref) - 1)
    assert names.count((CALL, "di_call")) == 1, names

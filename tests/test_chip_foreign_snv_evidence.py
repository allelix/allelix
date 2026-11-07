# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Only complete own-rsID SNV pairs judge a chip's base calls.

An indel-only rsID remains unjudged, even if one side of its identity pair
is a base. An indel row cannot add a permitted base to an own SNV pair.
REF population itself keeps its existing rules; these tests concern only
the new ClinVar exclusion and its count.
"""

import pytest
from test_chip_foreign_alleles import BUILDS, POS, RSA, RSB, _analyze, _chip, _data

from allelix.annotators._resolution import resolve_batch
from allelix.annotators.rsid_identity import RsidIdentityAnnotator
from allelix.models import Variant

CASES = [
    pytest.param([("G", "GA")], "AA", False, id="insertion-only-AA"),
    pytest.param([("G", "GA")], "CC", False, id="insertion-only-CC"),
    pytest.param([("GC", "C")], "AA", False, id="deletion-only-AA"),
    pytest.param([("G", "A"), ("C", "CT")], "CC", True, id="indel-ref-not-evidence"),
    pytest.param([("G", "A"), ("TC", "C")], "CC", True, id="indel-alt-not-evidence"),
    pytest.param([("G", "A"), ("C", "CT")], "AA", False, id="own-snv-control"),
]


def _database(tmp_path, pairs):
    rows = [(RSA, "1", POS, ref, alt, POS, ref, alt, 0) for ref, alt in pairs]
    rows.append((RSB, "1", POS, "G", "C", POS, "G", "C", 1))
    return _data(tmp_path, identity_rows=rows)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("pairs,call,foreign", CASES)
def test_resolved_chip_uses_only_complete_snv_pairs(tmp_path, build, pairs, call, foreign):
    data = _database(tmp_path, pairs)
    variant = Variant(
        rsid=RSA,
        chromosome="1",
        position=POS,
        allele1=call[0],
        allele2=call[1],
        build=build,
    )
    with RsidIdentityAnnotator(data) as identity:
        resolve_batch([variant], clinvar=None, gnomad=None, identity=identity)
    assert variant.foreign_alleles is foreign


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("pairs,call,foreign", CASES)
def test_cli_counts_only_complete_snv_disagreements(tmp_path, build, pairs, call, foreign):
    # in the analysis, ClinVar's own G>A row of rsA (in ``_data``) is
    # evidence too, so a call outside {G, A} is counted even where the
    # identity table holds no SNV row of rsA (the insertion-only CC case).
    foreign = foreign or not set(call) <= {"G", "A"}
    data = _database(tmp_path, pairs)
    report, _ = _analyze(tmp_path, data, _chip(tmp_path, call, build=build), build)
    assert report["input"]["probe_classes"]["foreign_allele_calls_at_plp_loci"] == int(foreign)
    if foreign:
        assert not [row for row in report["annotations"] if row["source"] == "clinvar"]

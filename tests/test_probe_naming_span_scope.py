# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Extend D/I probe naming to covering spans without losing exact-locus names.

A long deletion can cover the probe despite starting more than 69 bases
before it. Existing exact-locus P/LP naming also survives the extension.
Neither is interpretation; no ClinVar finding may be emitted for D/I.
"""

import pytest
from test_outside_review_rc1_b import (
    BUILDS,
    POS_A,
    RSB,
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

PROBE = "i3999261099"


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("path", ["single", "batch", "cli"])
@pytest.mark.parametrize(
    "ref,alt,offset,named",
    [
        pytest.param("G", "A", 0, True, id="existing-exact-snv-locus"),
        pytest.param("G" + "A" * 79, "G", 75, True, id="inside-long-deletion"),
        pytest.param("GA", "G", 1, True, id="inside-short-deletion-control"),
        pytest.param("GA", "G", 20, False, id="outside-span-control"),
    ],
)
def test_probe_naming_covers_the_span_and_keeps_existing_names(
    tmp_path, build, path, ref, alt, offset, named
):
    data = tmp_path / "data"
    _config(data)
    _clinvar(data, build, [(RSB, "1", POS_A, ref, alt, "GENEP", "Pathogenic", "Pathogenic")])
    _identity(data, [])
    pos = POS_A + offset
    if path == "cli":
        source = _chip(tmp_path / "chip.txt", [(PROBE, "1", pos, "DI")], build=build)
        report, _ = _analyze(tmp_path, data, source, build)
        assert not _clinvar_findings(report), "a D/I probe must never become a ClinVar finding"
        names = _named(report)
    else:
        call = Variant(
            rsid=PROBE, chromosome="1", position=pos, allele1="D", allele2="I", build=build
        )
        with ClinVarAnnotator(data, builds=(build,)) as clinvar:
            findings = (
                clinvar.annotate(call) if path == "single" else clinvar.batch_annotate([call])
            )
            assert not list(findings), "a D/I probe must never become a ClinVar finding"
            # The old accessor permits a direct RC1 preservation check too.
            getter = getattr(clinvar, "unchecked", None) or clinvar.unchecked_at_plp_loci
            names = list(getter())
    assert ((PROBE, "di_call") in names) is named, names
    assert len([item for item in names if item[0] == PROBE]) == int(named)

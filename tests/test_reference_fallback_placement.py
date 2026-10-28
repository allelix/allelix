# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Normal-allele evidence belongs to the reference-differing placement only."""

from __future__ import annotations

from contextlib import ExitStack

import pytest

from allelix.annotators import get_annotators
from allelix.models import Variant
from tests.test_annotation_parity_and_named_disclosures import (
    _gnomad_site,
    _identity,
    _pharm,
    _run,
    _vcf,
)


@pytest.mark.parametrize("batch", [False, True], ids=["single", "batch"])
@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("placement", ["matched", "position", "chromosome", "other_identity"])
def test_normal_allele_fallback_does_not_recode_another_placement(
    tmp_path, batch, build, placement
):
    rsid = "rs900252901"
    identity = [(rsid, "1", 100, "T", "C", 300, "C", "T", 0)]
    if placement == "other_identity":
        # This second placement has the same reference on both builds.
        # The new reference-carried fallback must not change its finding.
        identity.append((rsid, "1", 200, "C", "T", 400, "C", "T", 0))
    _identity(tmp_path, identity)
    _gnomad_site(tmp_path, rsid, "1", 300, "C", "T", 0.9)
    _pharm(tmp_path, rsid, "SYN", ["TT"])
    chrom = "2" if placement == "chromosome" else "1"
    position = 100 if build == "GRCh37" else 300
    if placement in {"position", "other_identity"}:
        position = 200 if build == "GRCh37" else 400
    ref = "T" if placement == "matched" and build == "GRCh37" else "C"
    variant = Variant(rsid, chrom, position, "T", "T", build=build, ref=ref, ref_from_input=True)
    with ExitStack() as stack:
        sources = get_annotators(tmp_path)
        for source in sources:
            stack.enter_context(source)
        pharm = next(source for source in sources if source.name == "pharmgkb")
        rows = list(pharm.batch_annotate([variant])) if batch else pharm.annotate(variant)
    if placement == "matched":
        # gnomAD resolves T as normal at this placement, so TT is suppressed.
        assert rows == []
    else:
        # ClinPGx matched this genotype before the fallback was introduced.
        # Preserve both the finding and its carried ALT outside the target site.
        assert len(rows) == 1
        assert rows[0].alt == "T"


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
def test_cli_keeps_the_carried_alt_outside_the_reference_differing_placement(tmp_path, build):
    rsid = "rs900252902"
    data = tmp_path / "data"
    data.mkdir()
    _identity(data, [(rsid, "1", 100, "T", "C", 300, "C", "T", 0)])
    _gnomad_site(data, rsid, "1", 300, "C", "T", 0.9)
    _pharm(data, rsid, "SYN", ["TT"])
    source = _vcf(tmp_path / "calls.vcf", [(rsid, "1", 200, "C", "T", "1/1", ".")])
    report, _ = _run(tmp_path, data, source, build, "json")
    rows = [row for row in report["annotations"] if row["source"] == "pharmgkb"]
    assert len(rows) == 1
    assert rows[0]["alt"] == "T"

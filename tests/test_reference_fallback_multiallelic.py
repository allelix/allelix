# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""A minor allele cannot be called normal from one alternative's low frequency."""

import contextlib
import sqlite3

import pytest
from test_annotation_parity_and_named_disclosures import _gnomad_site, _identity, _pharm

from allelix.annotators import get_annotators
from allelix.models import Variant


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("batch", [False, True])
def test_minor_reference_is_not_suppressed_at_three_allele_site(tmp_path, build, batch):
    data = tmp_path / "data"
    data.mkdir()
    rsid = "rs900252843"
    _identity(data, [(rsid, "1", 200, "T", "C", 1000, "C", "T", 0)])
    _pharm(data, rsid, "SYNTHETIC", ("CC",), level="2B")
    _gnomad_site(data, rsid, "1", 1000, "C", "T", 0.1)
    with contextlib.closing(sqlite3.connect(data / "gnomad.sqlite")) as conn:
        conn.execute(
            "INSERT INTO gnomad_frequencies (chrom,pos,ref,alt,rsid,af) "
            "VALUES ('1',1000,'C','G',?,0.8)",
            (rsid,),
        )
        conn.commit()
    with contextlib.ExitStack() as stack:
        registry = get_annotators(data, clinvar_builds=("GRCh37", "GRCh38"))
        for annotator in registry:
            stack.enter_context(annotator)
        pharm = next(a for a in registry if a.name == "pharmgkb")
        variant = Variant(
            rsid,
            "1",
            200 if build == "GRCh37" else 1000,
            "C",
            "C",
            ref="T" if build == "GRCh37" else "C",
            build=build,
            ref_from_input=True,
        )
        rows = list(pharm.batch_annotate([variant])) if batch else pharm.annotate(variant)
    # G is the major allele (0.8); C and T are each 0.1. C/C is not normal.
    assert len(rows) == 1

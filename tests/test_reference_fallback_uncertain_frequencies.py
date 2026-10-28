# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Unknown ALT frequencies and tied major alleles do not suppress real rows."""

import contextlib
import sqlite3

import pytest
from test_annotation_parity_and_named_disclosures import _gnomad_site, _identity, _pharm

from allelix.annotators import get_annotators
from allelix.models import Variant


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("batch", [False, True])
@pytest.mark.parametrize(
    "af_t,af_g,genotype,shown",
    [
        (0.1, None, "CC", True),
        (0.4, 0.2, "TT", True),
        (None, None, "CC", True),
        (0.6, 0.6, "TT", True),
        (0.1, 0.3, "CC", False),
        (0.7, 0.1, "TT", False),
    ],
    ids=["unknown-alt", "ref-alt-tie", "all-unknown", "sum-over-one", "major-ref", "major-alt"],
)
def test_fallback_requires_a_known_unique_major_allele(
    tmp_path, build, batch, af_t, af_g, genotype, shown
):
    data = tmp_path / "data"
    data.mkdir()
    rsid = "rs900252844"
    _identity(data, [(rsid, "1", 200, "T", "C", 1000, "C", "T", 0)])
    _pharm(data, rsid, "SYNTHETIC", (genotype,), level="2B")
    _gnomad_site(data, rsid, "1", 1000, "C", "T", af_t)
    with contextlib.closing(sqlite3.connect(data / "gnomad.sqlite")) as conn:
        conn.execute(
            "INSERT INTO gnomad_frequencies (chrom,pos,ref,alt,rsid,af) "
            "VALUES ('1',1000,'C','G',?,?)",
            (rsid, af_g),
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
            *genotype,
            ref="T" if build == "GRCh37" else "C",
            build=build,
            ref_from_input=True,
        )
        rows = list(pharm.batch_annotate([variant])) if batch else pharm.annotate(variant)
    assert len(rows) == int(shown)

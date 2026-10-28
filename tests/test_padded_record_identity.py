# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Trimming preserves each originating call when records share a normalized locus."""

import contextlib

import pytest
from test_annotation_parity_and_named_disclosures import _mthfr

from allelix.annotators import get_annotators
from allelix.models import Variant
from allelix.reports._pipeline import _annotate_batch


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("source", ["pharmgkb", "snpedia", "gwas"])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("shape", ["different-padding", "different-call"])
def test_trimmed_findings_keep_their_own_record(tmp_path, build, source, reverse, shape):
    data = tmp_path / "data"
    pos = _mthfr(data, build)
    first = Variant("rs1801133", "1", pos, "GC", "AC", ref="GC", build=build, ref_from_input=True)
    second = Variant(
        "rs1801133",
        "1",
        pos,
        "GCC" if shape == "different-padding" else "AC",
        "ACC" if shape == "different-padding" else "AC",
        ref="GCC" if shape == "different-padding" else "GC",
        build=build,
        ref_from_input=True,
    )
    variants = [first, second][:: -1 if reverse else 1]
    with contextlib.ExitStack() as stack:
        registry = get_annotators(data, clinvar_builds=(build,))
        for annotator in registry:
            stack.enter_context(annotator)
        annotator = next(a for a in registry if a.name == source)
        expected = [a for v in variants for a in _annotate_batch(annotator, [v])]
        actual = _annotate_batch(annotator, variants)

    def signature(a):
        return a.genotype_match, a.record_key, a.alt, a.carried_alleles

    assert expected
    assert [signature(a) for a in actual] == [signature(a) for a in expected]

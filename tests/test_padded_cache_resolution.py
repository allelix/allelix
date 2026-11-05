# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Recover calls when a source cache preserves the same padded VCF spelling."""

import pytest
from test_annotation_parity_and_named_disclosures import _clinvar, _gnomad_site, _identity, _vcf

from allelix.annotators import get_annotators
from allelix.annotators._resolution import resolve_batch
from allelix.parsers.vcf import VcfParser


@pytest.mark.parametrize(
    "provider,build",
    [
        ("identity", "GRCh37"),
        ("identity", "GRCh38"),
        ("clinvar", "GRCh37"),
        ("clinvar", "GRCh38"),
        ("gnomad", "GRCh38"),
    ],
)
@pytest.mark.parametrize("ref,alt", [("GC", "AC"), ("CG", "CA"), ("CGC", "CAC")])
def test_padded_cache_call_keeps_its_previously_recoverable_rsid(
    tmp_path, provider, build, ref, alt
):
    data = tmp_path / "data"
    data.mkdir()
    rsid = "rs900260381"
    pos = 100 if build == "GRCh37" else 200
    if provider == "identity":
        _identity(data, [(rsid, "1", 100, ref, alt, 200, ref, alt, 0)])
    elif provider == "gnomad":
        _gnomad_site(data, rsid, "1", pos, ref, alt, 0.01)
    else:
        _clinvar(
            data,
            build,
            [(rsid, "1", pos, ref, alt, "SYNTHETIC", "Pathogenic", "Pathogenic", 1, 1)],
        )
    path = _vcf(tmp_path / "input.vcf", [(".", "1", pos, ref, alt, "0/1", ".")])
    variants = list(VcfParser().parse(path))
    variants[0].build = build
    original = variants[0].rsid
    registry = get_annotators(data, clinvar_builds=(build,))
    by = {a.name: a for a in registry}
    try:
        resolve_batch(
            variants,
            clinvar=by["clinvar"] if provider == "clinvar" else None,
            gnomad=by["gnomad"] if provider == "gnomad" else None,
            identity=by["rsid_identity"] if provider == "identity" else None,
        )
        assert variants[0].rsid == rsid
        assert variants[0].original_rsid == original
        assert (variants[0].position, variants[0].ref) == (pos, ref)
        assert (variants[0].allele1, variants[0].allele2) == (ref, alt)
    finally:
        for annotator in registry:
            annotator.close()

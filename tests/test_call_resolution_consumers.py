"""Shared resolver and both extraction paths resolve the same normalized call."""

import pytest
from click.testing import CliRunner
from test_extract_reference_calls import _input, _row
from validated_call_helpers import BUILDS, RSID, call, setup

from allelix.annotators import get_annotators
from allelix.annotators._resolution import resolve_batch, resolve_rsidless
from allelix.cli import main, utility
from allelix.parsers.vcf import VcfParser


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("shape", ["prefix", "suffix", "both"])
@pytest.mark.parametrize("stream", [False, True])
def test_shared_resolution_before_each_command_consumer(tmp_path, build, shape, stream):
    data = tmp_path / "data"
    pos = setup(data, build)
    path = call(tmp_path / "input.vcf", pos, shape=shape, rsid=".")
    variants = list(VcfParser().parse(path))
    assert len(variants) == 1
    variants[0].build = build
    registry = get_annotators(data, clinvar_builds=(build,))
    by = {a.name: a for a in registry}
    try:
        providers = {n: by[n] for n in ("clinvar", "gnomad")}
        providers["identity"] = by["rsid_identity"]
        if stream:
            variants = list(resolve_rsidless(variants, **providers))
        else:
            resolve_batch(variants, **providers)
        assert variants[0].rsid == RSID
    finally:
        for annotator in registry:
            annotator.close()


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("shape", ["prefix", "suffix"])
@pytest.mark.parametrize("indexed", [False, True])
@pytest.mark.parametrize("gt", ["0/1", "./1"])
def test_extract_resolves_padded_unnamed_calls(tmp_path, monkeypatch, build, shape, indexed, gt):
    data = tmp_path / "data"
    pos = setup(data, build)
    prefix = shape == "prefix"
    ref, alt = ("CG", "CA") if prefix else ("GC", "AC")
    path = _input(
        tmp_path, build, [_row(pos - prefix, gt, name=".", ref=ref, alt=alt)], indexed=indexed
    )
    found = []
    monkeypatch.setattr(utility, "_render_extract_table", lambda p, w, f: found.append(f))
    result = CliRunner().invoke(
        main,
        ["extract", str(path), "--snps", RSID, "--build", build.lower(), "--data-dir", str(data)],
    )
    assert result.exit_code == 0, (result.output, result.exception)
    assert len(found) == 1 and RSID in found[0]
    assert found[0][RSID].is_half_call == (gt == "./1")


@pytest.mark.parametrize("build", BUILDS)
def test_repeat_indel_is_not_left_aligned_during_resolution(tmp_path, build):
    from test_annotation_parity_and_named_disclosures import _identity, _vcf

    data = tmp_path / "data"
    data.mkdir()
    _identity(data, [(RSID, "1", 100, "TT", "T", 100, "TT", "T", 0)])
    path = _vcf(tmp_path / "input.vcf", [(".", "1", 101, "TT", "T", "0/1", ".")])
    variants = list(VcfParser().parse(path))
    variants[0].build = build
    registry = get_annotators(data, clinvar_builds=(build,))
    by = {a.name: a for a in registry}
    try:
        resolve_batch(
            variants, clinvar=by["clinvar"], gnomad=by["gnomad"], identity=by["rsid_identity"]
        )
        assert variants[0].rsid != RSID
    finally:
        for annotator in registry:
            annotator.close()

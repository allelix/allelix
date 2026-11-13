"""Disabled sources cannot influence another source through bound callbacks."""

import json
from collections import Counter

import pytest
from test_annotation_parity_and_named_disclosures import _f5, _run, _vcf
from validated_call_helpers import call, setup

from allelix.annotators.clinvar import ClinVarAnnotator
from allelix.annotators.gnomad import GnomadAnnotator
from allelix.annotators.rsid_identity import RsidIdentityAnnotator
from allelix.config import load_config, save_config

CLASSES = {
    "gnomad": GnomadAnnotator,
    "clinvar": ClinVarAnnotator,
    "rsid_identity": RsidIdentityAnnotator,
}


def prepare(tmp_path, source, method):
    data = tmp_path / "data"
    if source == "rsid_identity":
        pos = setup(data, "GRCh38")
        path, build = call(tmp_path / "input.vcf", pos), "GRCh38"
    else:
        _f5(data, ("GRCh37",) if source == "clinvar" else (), with_gnomad=source == "gnomad")
        path = _vcf(tmp_path / "input.vcf", [("rs6025", "1", 169519049, "T", "C", "1/1", ".")])
        build = "GRCh37"
    extra = ("--min-magnitude", "0", "--gwas-min-magnitude", "0", "--snpedia-min-magnitude", "0")
    if method == "flag":
        extra += ("--no-gnomad",)
    else:
        config = load_config(data, create=False)
        config.sources[source] = False
        save_config(data, config)
    return data, path, build, extra


@pytest.mark.parametrize(
    "source,method",
    [("gnomad", "flag"), ("gnomad", "config"), ("clinvar", "config"), ("rsid_identity", "config")],
)
@pytest.mark.parametrize("command", ["analyze", "pharmacogenomics", "methylation"])
def test_disabled_dependency_cache_presence_cannot_change_answers(
    tmp_path, source, method, command
):
    data, path, build, extra = prepare(tmp_path, source, method)
    if command == "pharmacogenomics":
        # Per-source GWAS/SNPedia floors cannot apply to this focused report.
        extra = ("--min-magnitude", "0") + (("--no-gnomad",) if method == "flag" else ())
    reports = []
    for present in (True, False):
        if not present:
            pattern = "clinvar.*.sqlite" if source == "clinvar" else source + ".sqlite"
            for cache in data.glob(pattern):
                cache.rename(cache.with_suffix(".parked"))
        report, _ = _run(
            tmp_path, data, path, build, "json", command=command, extra=extra, label=str(present)
        )
        state = next(s for s in report["sources"] if s["name"] == source)
        assert state["state"] == "skipped" and state["version"] is None
        assert not any(a["source"] == source for a in report["annotations"])
        reports.append(Counter(json.dumps(a, sort_keys=True) for a in report["annotations"]))
    assert reports[0] == reports[1]


@pytest.mark.parametrize("source", list(CLASSES))
def test_disabled_source_is_not_queried_as_dependency(tmp_path, monkeypatch, source):
    data, path, build, extra = prepare(tmp_path, source, "config")
    calls = []
    cls = CLASSES[source]
    original = cls._connection

    def tracked(self, *args, **kwargs):
        calls.append(source)
        return original(self, *args, **kwargs)

    monkeypatch.setattr(cls, "_connection", tracked)
    _run(tmp_path, data, path, build, "json", extra=extra)
    assert calls == []


@pytest.mark.parametrize("surface", ["json", "html", "terminal"])
def test_unavailable_dependency_is_disclosed(tmp_path, surface):
    data, path, build, extra = prepare(tmp_path, "gnomad", "flag")
    result, console = _run(tmp_path, data, path, build, surface, extra=extra)
    if surface == "json":
        warnings = [
            w for w in result["report_warnings"] if w["code"] == "source_dependency_unavailable"
        ]
        assert warnings
        text = " ".join(w["message"] for w in warnings)
    else:
        text = result
    assert "gnomad" in text.lower()
    # Must explain the dependent interpretation, not merely repeat 'gnomAD skipped'.
    assert "normal allele" in text.lower()
    assert "normal allele" in console.lower()

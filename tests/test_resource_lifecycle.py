# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Analysis owns every acquired handle through the last provider use."""

from __future__ import annotations

import contextlib
import dataclasses
import sqlite3
from collections.abc import Iterator

import pytest
from test_readonly_caches import _data_dir
from test_sources import _run

import allelix.reports._pipeline as pipeline
from allelix.annotators import get_annotators
from allelix.annotators.alphamissense import AlphaMissenseAnnotator
from allelix.annotators.cadd import CaddAnnotator
from allelix.annotators.clinvar import ClinVarAnnotator
from allelix.annotators.gnomad import GnomadAnnotator
from allelix.annotators.gwas import GWASCatalogAnnotator
from allelix.annotators.pharmgkb import PharmGKBAnnotator
from allelix.annotators.rsid_identity import RsidIdentityAnnotator
from allelix.annotators.snpedia import SNPediaAnnotator
from allelix.parsers.vcf import VcfParser

_CLASSES = (
    ClinVarAnnotator,
    PharmGKBAnnotator,
    GWASCatalogAnnotator,
    SNPediaAnnotator,
    GnomadAnnotator,
    AlphaMissenseAnnotator,
    CaddAnnotator,
    RsidIdentityAnnotator,
)
_ENRICHMENT = {"gnomad", "alphamissense", "cadd", "rsid_identity"}


class Handles:
    def __init__(self, monkeypatch):
        self.owners = {}
        self.connections = {}
        self.closes = []
        for cls in _CLASSES:
            connection, close = cls._connection, cls.close

            def tracked_connection(owner, *args, original=connection, **kwargs):
                self.owners[id(owner)] = owner
                conn = original(owner, *args, **kwargs)
                if conn is not None:
                    self.connections[id(conn)] = (owner, conn)
                return conn

            def tracked_close(owner, original=close):
                self.owners[id(owner)] = owner
                self.closes.append(id(owner))
                original(owner)

            monkeypatch.setattr(cls, "_connection", tracked_connection)
            monkeypatch.setattr(cls, "close", tracked_close)

    def assert_closed(self):
        assert self.connections, "the control must actually acquire a SQLite handle"
        for owner, conn in self.connections.values():
            with pytest.raises(sqlite3.ProgrammingError, match="closed"):
                conn.execute("SELECT 1")
            assert self.closes.count(id(owner)) >= 1

    def cleanup(self):
        # Only after assertions: retain and close even a handle discarded from its slot.
        for owner in self.owners.values():
            owner.close()
        for _, conn in self.connections.values():
            conn.close()


@pytest.fixture
def handles(monkeypatch):
    tracked = Handles(monkeypatch)
    try:
        yield tracked
    finally:
        tracked.cleanup()


@pytest.fixture
def enriched_data(tmp_path, build_synthetic_clinvar_cache, mock_pharmgkb_dir, mock_gwas_tsv):
    data = _data_dir(tmp_path, build_synthetic_clinvar_cache, mock_pharmgkb_dir, mock_gwas_tsv)
    with contextlib.closing(sqlite3.connect(data / "alphamissense.sqlite")) as conn:
        conn.execute("PRAGMA journal_mode=DELETE")
    inserts = {
        "rsid_identity": (
            "INSERT INTO rsid_identity VALUES "
            "('rs1801133', '1', 11856378, 'G', 'A', 11796321, 'G', 'A', 0)"
        ),
        "gnomad": (
            "INSERT INTO gnomad_frequencies (chrom,pos,ref,alt,rsid,af) "
            "VALUES ('1',11796321,'G','A','rs1801133',0.02)"
        ),
        "alphamissense": (
            "INSERT INTO alphamissense_scores "
            "(chrom,pos,ref,alt,rsid,am_pathogenicity,am_class) "
            "VALUES ('1',11796321,'G','A','rs1801133',0.8,'likely_pathogenic')"
        ),
        "cadd": "INSERT INTO cadd_scores VALUES ('1',11796321,'G','A',23.0)",
    }
    for name, sql in inserts.items():
        with contextlib.closing(sqlite3.connect(data / f"{name}.sqlite")) as conn:
            conn.execute(sql)
            conn.commit()
    (data / "config.toml").write_text("[sources]\ncadd = true\n", encoding="utf-8")
    return data


def _vcf(tmp_path):
    path = tmp_path / "synthetic.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        "1\t11856378\trs1801133\tG\tA\t50\tPASS\t.\tGT\t0/1\n",
        encoding="utf-8",
    )
    return path


def _direct(path, registry, **kwargs):
    by_name = {a.name: a for a in registry}
    return pipeline.run_analysis(
        path,
        kwargs.pop("parser", VcfParser()),
        [a for a in registry if a.name not in _ENRICHMENT],
        build_override="GRCh37",
        **{name: by_name[name] for name in ("gnomad", "alphamissense", "cadd")},
        identity=by_name["rsid_identity"],
        **kwargs,
    )


def _after(monkeypatch, target, name, error):
    original = getattr(target, name)

    def fail(*args, **kwargs):
        result = original(*args, **kwargs)
        if isinstance(result, Iterator):
            list(result)  # Execute generator-backed annotation before injecting failure.
        raise error

    monkeypatch.setattr(target, name, fail)


def test_enrichment_success_repeat_and_deduplicated_ownership(enriched_data, tmp_path, handles):
    registry = get_annotators(enriched_data)
    path = _vcf(tmp_path)
    snapshots = []
    for run in (1, 2):
        result = _direct(path, registry)
        handles.assert_closed()
        assert all(handles.closes.count(id(owner)) == run for owner in registry)
        assert {name for name, _ in result.annotators_used} >= _ENRICHMENT
        enriched = [a for a in result.annotations if a.rsid == "rs1801133"]
        assert enriched
        assert all(
            (a.allele_frequency, a.am_pathogenicity, a.cadd_phred) == (0.02, 0.8, 23.0)
            for a in enriched
        )
        snapshots.append([dataclasses.asdict(a) for a in result.annotations])
    assert snapshots[0] == snapshots[1]


@pytest.mark.parametrize(
    "name,method",
    [
        ("rsid_identity", "translate_to_grch38"),
        ("gnomad", "bulk_lookup_by_position"),
        ("alphamissense", "bulk_lookup_by_position"),
        ("cadd", "bulk_lookup"),
    ],
)
@pytest.mark.parametrize("error_type", [RuntimeError, KeyboardInterrupt])
def test_each_enrichment_post_query_failure_closes_actual_handles(
    enriched_data,
    tmp_path,
    handles,
    monkeypatch,
    name,
    method,
    error_type,
):
    registry = get_annotators(enriched_data)
    owner = next(a for a in registry if a.name == name)
    _after(monkeypatch, owner, method, error_type("post-query failure"))
    with pytest.raises(error_type, match="post-query failure"):
        _direct(_vcf(tmp_path), registry)
    assert any(a is owner for a, _ in handles.connections.values())
    handles.assert_closed()


@pytest.mark.parametrize(
    "phase",
    [
        "resolve",
        "annotate",
        "parse",
        "placements",
        "late_build",
        "plausibility",
        "scoring",
        "version",
        "readiness",
        "counts",
        "result",
    ],
)
def test_analysis_failure_boundaries_close_handles(
    enriched_data, tmp_path, handles, monkeypatch, phase
):
    registry = get_annotators(enriched_data)
    by_name = {a.name: a for a in registry}
    parser = VcfParser()
    kwargs = {}
    error = RuntimeError("analysis phase failure")
    if phase == "resolve":
        _after(monkeypatch, pipeline, "resolve_batch", error)
    elif phase == "annotate":
        _after(monkeypatch, by_name["clinvar"], "batch_annotate", error)
    elif phase == "parse":
        monkeypatch.setattr(pipeline, "_BATCH_SIZE", 1)
        original = parser.parse

        def parse(path):
            yield from original(path)
            raise error

        monkeypatch.setattr(parser, "parse", parse)
    elif phase == "placements":
        kwargs["high_value_rsids"] = {"rs1801133"}
        _after(monkeypatch, by_name["rsid_identity"], "loci_of_rsids", error)
    elif phase == "late_build":
        _after(monkeypatch, pipeline, "check_late_seeds", error)
    elif phase in {"plausibility", "scoring"}:
        _after(
            monkeypatch,
            pipeline,
            "_flag_plausibility" if phase == "plausibility" else "_score_down_common_zero_star",
            error,
        )
    elif phase in {"version", "counts"}:
        _after(
            monkeypatch,
            by_name["clinvar"],
            "version" if phase == "version" else "run_counts",
            error,
        )
    elif phase == "readiness":
        original = pipeline._score_down_common_zero_star

        def score(*args):
            original(*args)
            _after(monkeypatch, by_name["rsid_identity"], "is_ready", error)

        monkeypatch.setattr(pipeline, "_score_down_common_zero_star", score)
    else:

        def fail():
            raise error

        kwargs["skipped_count_provider"] = fail
    with pytest.raises(RuntimeError, match="analysis phase failure"):
        _direct(_vcf(tmp_path), registry, parser=parser, **kwargs)
    handles.assert_closed()
    assert parser._watch_loci == {}


@pytest.mark.parametrize("consumer", ["pharmgkb", "snpedia", "both"])
@pytest.mark.parametrize("fail", [False, True])
def test_cli_owns_outdated_callback_clinvar_without_activating_it(
    enriched_data,
    tmp_path,
    handles,
    monkeypatch,
    consumer,
    fail,
):
    for build in ("GRCh37", "GRCh38"):
        with contextlib.closing(
            sqlite3.connect(enriched_data / f"clinvar.{build}.sqlite")
        ) as conn:
            conn.execute("UPDATE database_versions SET local_version_tag='iv:4'")
            conn.commit()
    if fail:
        _after(monkeypatch, ClinVarAnnotator, "reference_at", RuntimeError("callback failure"))
    sources = "[sources]\ncadd = true\ngwas = false\n"
    if consumer != "both":
        sources += "snpedia = false\n" if consumer == "pharmgkb" else "pharmgkb = false\n"
    (enriched_data / "config.toml").write_text(sources, encoding="utf-8")
    extra = []
    out = tmp_path / "callback.json"
    if fail:
        with pytest.raises(RuntimeError, match="callback failure"):
            _run(enriched_data, out, *extra)
    else:
        import json

        code, output = _run(enriched_data, out, *extra)
        assert code == 0, output
        result = json.loads(out.read_text(encoding="utf-8"))
        assert next(a for a in result["sources"] if a["name"] == "clinvar")["reason"] == "outdated"
        assert "clinvar" not in {a["name"] for a in result["annotators"]}
        assert any(w["code"] == "clinvar_skipped" for w in result["report_warnings"])
    assert any(isinstance(a, ClinVarAnnotator) for a, _ in handles.connections.values())
    handles.assert_closed()


def test_explicit_callback_owners_identity_dedup_and_unused_laziness(
    enriched_data, tmp_path, handles
):
    registry = get_annotators(enriched_data)
    clinvar, pharm, _, _, gnomad, _, _, identity = registry
    unused = GnomadAnnotator(enriched_data)  # Same source name, distinct owner.
    result = pipeline.run_analysis(
        _vcf(tmp_path),
        VcfParser(),
        [pharm],
        build_override="GRCh37",
        gnomad=gnomad,
        identity=identity,
        resource_owners=[*registry, clinvar, identity, unused],
    )
    handles.assert_closed()
    assert all(handles.closes.count(id(owner)) == 1 for owner in [*registry, unused])
    assert not any(owner is unused for owner, _ in handles.connections.values())
    assert "clinvar" not in {name for name, _ in result.annotators_used}


def test_resolved_stream_placement_query_failure_is_owned(
    enriched_data, tmp_path, handles, monkeypatch
):
    from allelix.cli._helpers import _resolved_variants

    _after(monkeypatch, RsidIdentityAnnotator, "loci_of_rsids", RuntimeError("placement failure"))
    with (
        pytest.raises(RuntimeError, match="placement failure"),
        _resolved_variants(
            VcfParser(),
            _vcf(tmp_path),
            data_dir=enriched_data,
            build_override="GRCh37",
            high_value={"rs1801133"},
        ),
    ):
        pytest.fail("failed acquisition must not yield a stream")
    handles.assert_closed()


def test_panel_placement_failure_clears_registered_watches(
    enriched_data, tmp_path, handles, monkeypatch
):
    registry = get_annotators(enriched_data)
    identity = registry[-1]
    parser = VcfParser()
    original = identity.forward_sites

    def fail(*args):
        assert parser._watch_rsids == frozenset({"rs1801133"})
        original(*args)
        raise RuntimeError("panel placement failed after registration")

    monkeypatch.setattr(identity, "forward_sites", fail)
    with pytest.raises(RuntimeError, match="panel placement failed after registration"):
        _direct(
            _vcf(tmp_path),
            registry,
            parser=parser,
            high_value_rsids={"rs1801133"},
            panel_rsids=frozenset({"rs1801133"}),
        )
    handles.assert_closed()
    assert parser._watch_loci == {}
    assert parser._watch_rsids == frozenset()
    assert parser._panel_rsids == frozenset()


@pytest.mark.parametrize("failure", [None, "GRCh37", "GRCh38"])
def test_zero_star_rescoring_owns_both_clinvar_builds(tmp_path, handles, monkeypatch, failure):
    from test_clinvar_common_zero_star import (
        _write_clinvar,
        _write_gnomad,
        _write_identity,
        _write_vcf,
    )

    data = tmp_path / "common"
    data.mkdir()
    for build in ("GRCh37", "GRCh38"):
        _write_clinvar(data, build)
    _write_gnomad(data)
    _write_identity(data)
    clinvar, gnomad, identity = (
        ClinVarAnnotator(data),
        GnomadAnnotator(data),
        RsidIdentityAnnotator(data),
    )
    original = clinvar.records_of_rsids

    def records(rsids, build):
        result = original(rsids, build)
        if build == failure:
            raise RuntimeError("rescoring post-query failure")
        return result

    monkeypatch.setattr(clinvar, "records_of_rsids", records)

    def run():
        return pipeline.run_analysis(
            _write_vcf(tmp_path / "common.vcf", "GRCh37"),
            VcfParser(),
            [clinvar],
            build_override="GRCh37",
            gnomad=gnomad,
            identity=identity,
        )

    if failure:
        with pytest.raises(RuntimeError, match="rescoring post-query failure"):
            run()
    else:
        result = run()
        assert any(a.common_allele_capped for a in result.annotations)
    handles.assert_closed()
    if failure != "GRCh37":
        assert sum(owner is clinvar for owner, _ in handles.connections.values()) >= 2


@pytest.mark.parametrize("consumer", ["clinvar", "pharmgkb", "gwas", "snpedia"])
def test_identity_reached_only_through_callback_is_owned(
    enriched_data, tmp_path, handles, consumer
):
    if consumer == "clinvar":
        # Current reference-carried row, without the other build: ClinVar's
        # callback to identity supplies the normal allele to its consumer.
        with contextlib.closing(sqlite3.connect(enriched_data / "clinvar.GRCh37.sqlite")) as conn:
            conn.execute("UPDATE clinvar_variants SET ref='A',alt='A' WHERE rsid='rs1801133'")
            conn.commit()
        (enriched_data / "clinvar.GRCh38.sqlite").unlink()
    elif consumer == "snpedia":
        with contextlib.closing(sqlite3.connect(enriched_data / "snpedia.sqlite")) as conn:
            conn.execute(
                "INSERT INTO snpedia_genotypes (rsid,allele1,allele2,magnitude) "
                "VALUES ('rs1801133','G','A',3)"
            )
            conn.commit()
    registry = get_annotators(enriched_data)
    active = next(
        owner
        for owner in registry
        if owner.name == ("pharmgkb" if consumer == "clinvar" else consumer)
    )
    identity = registry[-1]
    result = pipeline.run_analysis(
        _vcf(tmp_path),
        VcfParser(),
        [active],
        build_override="GRCh37",
        resource_owners=registry,
    )
    assert any(owner is identity for owner, _ in handles.connections.values())
    assert "rsid_identity" not in {name for name, _ in result.annotators_used}
    handles.assert_closed()


@pytest.mark.parametrize("preopen", [False, True])
def test_missing_sources_and_unsupported_build_do_not_activate_owners(
    tmp_path,
    handles,
    preopen,
):
    # One valid but unsupported cache plus entirely absent providers.
    cache = tmp_path / "gnomad.sqlite"
    with contextlib.closing(sqlite3.connect(cache)) as conn:
        conn.execute("CREATE TABLE synthetic (n INTEGER)")
    registry = get_annotators(tmp_path)
    gnomad = next(owner for owner in registry if owner.name == "gnomad")
    if preopen:
        gnomad._connection().execute("SELECT 1")
    result = pipeline.run_analysis(
        _vcf(tmp_path),
        VcfParser(),
        [],
        build_override="GRCh36",
        resource_owners=registry,
    )
    assert result.annotators_used == []
    assert result.annotations == []
    assert all(handles.closes.count(id(owner)) == 1 for owner in registry)
    if preopen:
        handles.assert_closed()
    else:
        assert handles.connections == {}


def test_primary_computation_error_survives_owner_cleanup_error(
    enriched_data, tmp_path, handles, monkeypatch
):
    owner = GnomadAnnotator(enriched_data)
    owner._connection().execute("SELECT 1")
    close = owner.close
    parser = VcfParser()
    with monkeypatch.context() as patcher:

        def cleanup():
            close()
            raise RuntimeError("secondary provider cleanup")

        patcher.setattr(owner, "close", cleanup)
        _after(patcher, parser, "get_metadata", ValueError("primary computation refusal"))
        with pytest.raises(ValueError, match="primary computation refusal"):
            pipeline.run_analysis(_vcf(tmp_path), parser, [], gnomad=owner)
        handles.assert_closed()


@pytest.mark.parametrize("fail", [False, True])
def test_cadd_full_backend_handles_close_after_enrichment(tmp_path, monkeypatch, fail):
    from allelix.models import Annotation

    class Tabix:
        def __init__(self):
            self.closed = False

        def fetch(self, *args):
            assert not self.closed
            return iter(["1\t100\tA\tG\t0\t19.0", "1\t200\tA\tAG\t0\t21.0"])

        def close(self):
            self.closed = True

    provider = CaddAnnotator(tmp_path, full_mode=True)
    snv, indel = Tabix(), Tabix()

    def snv_open():
        provider._tabix = snv
        return snv

    def indel_open():
        provider._indel_tabix = indel
        return indel

    monkeypatch.setattr(provider, "_open_tabix", snv_open)
    monkeypatch.setattr(provider, "_open_indel_tabix", indel_open)
    monkeypatch.setattr(provider, "is_ready", lambda: True)
    annotations = [
        Annotation(
            source="synthetic",
            rsid="rs990001",
            significance="synthetic",
            category="clinical",
            magnitude=5,
            description="synthetic",
            attribution="synthetic",
            genotype_match="A/G",
            locus=("1", pos),
            locus_ref="A",
            alt=alt,
            carried_alleles=("A", alt),
        )
        for pos, alt in [(100, "G"), (200, "AG")]
    ]
    monkeypatch.setattr(
        pipeline,
        "_enrichment_keys",
        lambda *args: [(ann, (*ann.locus, "A", ann.alt)) for ann in annotations],
    )
    if fail:
        _after(monkeypatch, provider, "bulk_lookup", RuntimeError("post-tabix query"))
    try:
        if fail:
            with pytest.raises(RuntimeError, match="post-tabix query"):
                pipeline.run_analysis(
                    _vcf(tmp_path), VcfParser(), [], build_override="GRCh37", cadd=provider
                )
        else:
            pipeline.run_analysis(
                _vcf(tmp_path), VcfParser(), [], build_override="GRCh37", cadd=provider
            )
            assert [ann.cadd_phred for ann in annotations] == [19.0, 21.0]
        assert snv.closed and indel.closed
    finally:
        provider.close()

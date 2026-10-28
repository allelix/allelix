# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Reference extraction does not add records to annotation or export streams."""

from __future__ import annotations

import contextlib
import json
import sqlite3

import pytest
from click.testing import CliRunner
from test_extract_reference_calls import RSID, _cache, _input, _row

from allelix.annotators.clinvar import clinvar_db_filename, clinvar_record_name
from allelix.cli import _helpers, main, utility
from allelix.databases._versions import CLINVAR_INTERPRETER_VERSION
from allelix.databases.schema import CLINVAR_SCHEMA
from allelix.parsers.vcf import VcfParser


@pytest.fixture(params=["GRCh37", "GRCh38"])
def reference_input(tmp_path, request):
    build = request.param
    placements = [("1", p, p + 800, "G", "A") for p in (200, 300, 400, 500, 600)]
    data = _cache(tmp_path, placements)
    positions = [p[1 if build == "GRCh37" else 2] for p in placements]
    with contextlib.closing(sqlite3.connect(data / clinvar_db_filename(build))) as conn:
        conn.executescript(CLINVAR_SCHEMA)
        conn.executemany(
            "INSERT INTO clinvar_variants "
            "(rsid,chromosome,position,ref,alt,clinical_significance,condition,gene,"
            "review_status,allele_id) VALUES (?, '1', ?, 'G', ?, 'Pathogenic', "
            "'Synthetic condition', 'SYNTHETIC', 'criteria provided, single submitter', ?)",
            [
                (RSID, pos, "G" if index == 0 else "A", 900252400 + index)
                for index, pos in enumerate(positions)
            ],
        )
        conn.execute(
            "INSERT INTO database_versions "
            "(name,source_url,version,downloaded_at,record_count,local_version_tag) "
            "VALUES (?, 'test://synthetic', 'synthetic', '2001-01-01', 5, ?)",
            (clinvar_record_name(build), f"iv:{CLINVAR_INTERPRETER_VERSION}"),
        )
        conn.commit()
    source = _input(
        tmp_path,
        build,
        [
            _row(positions[0]),  # Legitimate reference-carried ClinVar finding.
            _row(positions[1]),  # Ordinary ALT allele is absent.
            _row(positions[2], "0/1"),
            _row(
                positions[3] - 10,
                name=".",
                ref="A",
                alt="<NON_REF>",
                info=f"END={positions[3] + 10}",
            ),
            _row(positions[4], filt="LowQual"),
        ],
        indexed=False,
    )
    return build, data, source, positions


def _calls(variants):
    return [(v.rsid, v.position, v.allele1, v.allele2, v.record_index) for v in variants]


def test_extract_watch_never_changes_the_parser_yield_stream(reference_input):
    _build, _data, source, positions = reference_input
    parser = VcfParser()
    baseline = _calls(parser.parse(source))
    assert [(v[1], v[2], v[3]) for v in baseline] == [(positions[2], "G", "A")]
    parser.watch_extract(frozenset({RSID}), {"1": [(pos, "G") for pos in positions]})
    assert _calls(parser.parse(source)) == baseline
    assert len(parser.extract_candidates) == 4
    parser.watch_extract(frozenset())
    assert _calls(parser.parse(source)) == baseline
    assert parser.extract_candidates == []


def test_reusing_a_parser_after_extract_keeps_the_normal_resolved_stream(reference_input):
    build, data, source, positions = reference_input
    parser = VcfParser()
    found = utility._sequential_extract(parser, source, {RSID}, data, build)
    assert (found[RSID].position, found[RSID].allele1, found[RSID].allele2) == (
        positions[0],
        "G",
        "G",
    )
    with _helpers._resolved_variants(
        parser, source, data_dir=data, build_override=build
    ) as resolved:
        calls = _calls(resolved.variants)
    assert [(v[1], v[2], v[3]) for v in calls] == [(positions[2], "G", "A")]
    assert parser.extract_candidates == []


def test_extract_does_not_change_analysis_vcf_or_plink_outputs(tmp_path, reference_input):
    build, data, source, positions = reference_input
    runner = CliRunner()
    outputs = []
    for phase in ("before", "after"):
        if phase == "after":
            extracted = runner.invoke(
                main,
                [
                    "extract",
                    str(source),
                    "--snps",
                    RSID,
                    "--data-dir",
                    str(data),
                    "--build",
                    build.lower(),
                ],
            )
            assert extracted.exit_code == 0, extracted.output
            assert "G/G" in extracted.output
        report = tmp_path / f"{phase}.json"
        vcf = tmp_path / f"{phase}.vcf"
        args = [str(source), "--data-dir", str(data), "--build", build.lower()]
        analyzed = runner.invoke(
            main,
            [
                "analyze",
                *args,
                "--no-update",
                "--min-magnitude",
                "0",
                "--output",
                str(report),
                "--vcf-out",
                str(vcf),
            ],
        )
        assert analyzed.exit_code == 0, analyzed.output
        payload = json.loads(report.read_text())
        assert sorted(a["position"] for a in payload["annotations"]) == [
            positions[0],
            positions[2],
        ]
        prefix = tmp_path / phase
        exported = runner.invoke(main, ["export", "plink", *args, "--output-prefix", str(prefix)])
        assert exported.exit_code == 0, exported.output
        bim = prefix.with_suffix(".bim").read_text().splitlines()
        assert len(bim) == 1 and int(bim[0].split()[3]) == positions[2]
        records = [line for line in vcf.read_text().splitlines() if not line.startswith("#")]
        assert "ALLELIX_CLINVAR_REF=" in records[0]
        assert all("ALLELIX_CLINVAR" not in records[i] for i in (1, 3, 4))
        outputs.append(
            (
                payload["annotations"],
                payload["input"]["vcf_quality"],
                records,
                [prefix.with_suffix(suffix).read_bytes() for suffix in (".bed", ".bim", ".fam")],
            )
        )
    assert outputs[0] == outputs[1]

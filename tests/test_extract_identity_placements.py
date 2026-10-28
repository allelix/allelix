# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Missing gnomAD placements must not hide calls placed by identity.

All rows, coordinates, identifiers and calls in this module are synthetic.
"""

from __future__ import annotations

import contextlib
import sqlite3
from types import SimpleNamespace

import pysam
import pytest
from click.testing import CliRunner

from allelix.annotators._resolution import resolve_batch
from allelix.annotators.gnomad import GnomadAnnotator
from allelix.annotators.rsid_identity import RsidIdentityAnnotator
from allelix.cli import main, utility
from allelix.databases._versions import GNOMAD_SCHEMA_VERSION
from allelix.databases.rsid_identity_loader import RSID_IDENTITY_DB_FILENAME, write_identity_db
from allelix.databases.schema import GNOMAD_SCHEMA
from allelix.exporters.plink import (
    export_plink,
    resolve_ref_alt_via_gnomad_grch37,
    resolve_ref_alt_via_gnomad_records,
    resolve_ref_alt_via_identity,
)
from allelix.models import Variant

TARGET = "rs900252901"
ANCHOR = "rs900252902"


def _data(
    tmp_path,
    chrom="1",
    *,
    gnomad=True,
    identity=True,
    target_mapped=False,
    target_position=900,
    multiallelic=False,
    extra=(),
):
    data = tmp_path / "data"
    data.mkdir()
    if multiallelic:
        extra = [*extra, (TARGET, chrom, 200, "G", "T", 1000, "G", "T", 0)]
    if identity:
        write_identity_db(
            [(TARGET, chrom, 200, "G", "A", 1000, "G", "A", 0), *extra],
            data / RSID_IDENTITY_DB_FILENAME,
            version="synthetic",
            source_url="test://synthetic",
            downloaded_at="2001-01-01",
        )
    if gnomad:
        rows = [(chrom, 100, "C", "T", ANCHOR, 0.1)]
        if target_mapped:
            rows.append((chrom, target_position, "G", "A", TARGET, 0.1))
            if multiallelic:
                rows.append((chrom, target_position, "G", "T", TARGET, 0.1))
        with contextlib.closing(sqlite3.connect(data / "gnomad.sqlite")) as conn:
            conn.executescript(GNOMAD_SCHEMA)
            conn.executemany(
                "INSERT INTO gnomad_frequencies (chrom,pos,ref,alt,rsid,af) VALUES (?,?,?,?,?,?)",
                rows,
            )
            conn.execute(
                "INSERT INTO database_versions "
                "(name,source_url,version,downloaded_at,record_count,local_version_tag) "
                "VALUES ('gnomad','test://synthetic','1','2001-01-01',?,?)",
                (len(rows), f"sv:{GNOMAD_SCHEMA_VERSION}"),
            )
            conn.commit()
    return data


def _row(
    pos,
    *,
    chrom="1",
    name=".",
    ref="G",
    alt="A",
    gt="0/0",
    filt="PASS",
    gq=50,
    info=".",
    second=None,
):
    sample = f"{gt}:{gq}:20"
    if second is not None:
        sample += f"\t{second}:50:20"
    return f"{chrom}\t{pos}\t{name}\t{ref}\t{alt}\t.\t{filt}\t{info}\tGT:GQ:DP\t{sample}\n"


def _vcf(tmp_path, rows, chrom="1", *, indexed=False, second=False, build="GRCh38"):
    lengths = (
        {"1": 248956422, "X": 156040895, "Y": 57227415}
        if build == "GRCh38"
        else {"1": 249250621, "X": 155270560, "Y": 59373566}
    )
    source = tmp_path / "calls.vcf"
    source.write_text(
        "##fileformat=VCFv4.2\n"
        f"##contig=<ID={chrom},length={lengths[chrom]}>\n"
        '##ALT=<ID=NON_REF,Description="Any other allele">\n'
        '##INFO=<ID=END,Number=1,Type=Integer,Description="Block end">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"
        + ("\tOTHER" if second else "")
        + "\n"
        + "".join(rows)
    )
    compressed = tmp_path / "calls.vcf.gz"
    pysam.tabix_compress(str(source), str(compressed))
    if indexed:
        pysam.tabix_index(str(compressed), preset="vcf")
    return compressed


def _extract(monkeypatch, path, data, *, sample=None, build="GRCh38"):
    captured = []
    render = utility._render_extract_table

    def capture(source, wanted, found):
        captured.append(found)
        render(source, wanted, found)

    monkeypatch.setattr(utility, "_render_extract_table", capture)
    args = [
        "extract",
        str(path),
        "--snps",
        TARGET,
        "--data-dir",
        str(data),
        "--build",
        build.lower(),
    ]
    if sample:
        args += ["--sample", sample]
    result = CliRunner(env={"COLUMNS": "180"}).invoke(main, args)
    assert result.exit_code == 0, (result.output, result.exception)
    assert len(captured) == 1
    return captured[0], result.output


@pytest.mark.parametrize("chrom", ["1", "Y"])
@pytest.mark.parametrize(
    "build,pos,provider",
    [
        ("GRCh37", 200, "identity-only"),
        ("GRCh37", 200, "gnomad-and-identity"),
        ("GRCh38", 1000, "identity-only"),
        ("GRCh38", 1000, "gnomad-and-identity"),
        ("GRCh38", 1000, "gnomad-only"),
    ],
)
@pytest.mark.parametrize("indexed", [False, True], ids=["sequential", "index-present"])
@pytest.mark.parametrize(
    "shape,gt,expected",
    [
        ("het", "0/1", ("G", "A")),
        ("haploid-alt", "1", ("A", "A")),
        ("reference", "0/0", ("G", "G")),
        ("phased-reference", "0|0", ("G", "G")),
        ("haploid-ref", "0", ("G", "G")),
        ("alt-dot", "0/0", ("G", "G")),
        ("half-reference", "0/.", "G"),
        ("half-reference-right", "./0", "G"),
        ("half-alt", "./1", "A"),
        ("half-alt-left", "1/.", "A"),
        ("phased-half-reference", "0|.", "G"),
        ("phased-half-alt", ".|1", "A"),
        ("alt-dot-half", "0/.", "G"),
        ("no-call", "./.", None),
        ("haploid-no-call", ".", None),
        ("filtered", "0/0", None),
        ("low-gq", "0/0", None),
        ("block", "0/0", ("G", "G")),
        ("multi-het", "1/2", ("A", "T")),
        ("multi-reverse", "2/1", ("T", "A")),
        ("multi-half", "./2", "T"),
        ("multi-half-left", "2/.", "T"),
        ("multi-haploid-alt", "2", ("T", "T")),
        ("multi-no-call", "./.", None),
    ],
)
def test_placed_record_keeps_its_call(
    tmp_path, monkeypatch, chrom, indexed, shape, gt, expected, provider, build, pos
):
    multi = shape.startswith("multi")
    data = _data(
        tmp_path,
        chrom,
        target_mapped=provider != "identity-only",
        identity=provider != "gnomad-only",
        target_position=1000,
        multiallelic=multi,
    )
    row = _row(
        pos,
        chrom=chrom,
        gt=gt,
        alt="." if shape.startswith("alt-dot") else "A,T" if multi else "A",
        filt="LowQual" if shape == "filtered" else "PASS",
        gq=5 if shape == "low-gq" else 50,
    )
    if shape == "block":
        # Use the placement's G, not the reference block anchor's C.
        row = _row(pos - 10, chrom=chrom, ref="C", alt="<NON_REF>", info=f"END={pos + 10}")
    path = _vcf(tmp_path, [row], chrom, indexed=indexed, build=build)
    found, output = _extract(monkeypatch, path, data, build=build)
    assert TARGET in found, output
    call = found[TARGET]
    assert call.genotype in output and "not in file" not in output
    assert (call.chromosome, call.position) == (chrom, pos)
    if expected is None:
        assert call.is_no_call and not call.is_half_call
    elif isinstance(expected, str):
        assert call.is_half_call and expected in call.carried_pair
    else:
        assert sorted(call.carried_pair) == sorted(expected)


@pytest.mark.parametrize("shape", ["reference", "alt-dot", "block", "haploid", "no-call"])
def test_index_fetches_identity_only_sites_in_mixed_request(tmp_path, monkeypatch, shape):
    data = _data(tmp_path, "Y")
    row = _row(
        1000,
        chrom="Y",
        gt={"haploid": "0", "no-call": "./."}.get(shape, "0/0"),
        alt="." if shape == "alt-dot" else "A",
    )
    if shape == "block":
        row = _row(990, chrom="Y", ref="C", alt="<NON_REF>", info="END=1010")
    path = _vcf(
        tmp_path,
        [_row(100, chrom="Y", name=ANCHOR, ref="C", alt="T", gt="0/1"), row],
        "Y",
        indexed=True,
    )
    fetches = []
    scans = []
    sequential = utility._sequential_extract

    class TracedIndex:
        def __init__(self, path):
            self.inner = pysam.TabixFile(path)
            self.contigs = self.inner.contigs

        def fetch(self, chrom, start, end):
            fetches.append((chrom, start, end))
            return self.inner.fetch(chrom, start, end)

        def close(self):
            self.inner.close()

    def scan(parser, source, wanted, *args):
        scans.append(wanted)
        return sequential(parser, source, wanted, *args)

    monkeypatch.setattr(utility, "_sequential_extract", scan)
    found = utility._execute_tabix_extract(
        path, {TARGET, ANCHOR}, None, data, SimpleNamespace(TabixFile=TracedIndex), "GRCh38"
    )
    assert found is not None and set(found) == {TARGET, ANCHOR}
    assert (found[TARGET].chromosome, found[TARGET].position) == ("Y", 1000)
    assert any(c == "Y" and start <= 999 < end for c, start, end in fetches)
    assert scans == [], "An indexed identity placement must not need a full-file scan"


@pytest.mark.parametrize("indexed", [False, True])
@pytest.mark.parametrize("usable_first", [False, True])
def test_identity_placements_keep_mapped_usable_priority(
    tmp_path, monkeypatch, indexed, usable_first
):
    extra = [(TARGET, "1", 300, "G", "A", 1100, "G", "A", 0)]
    data = _data(tmp_path, extra=extra)
    good, bad = (1000, 1100) if usable_first else (1100, 1000)
    rows = [
        _row(100, name=TARGET, gt="./."),
        _row(1000, gt="0/0" if good == 1000 else "./."),
        _row(1100, gt="0/0" if good == 1100 else "./."),
    ]
    assert good != bad
    found, output = _extract(monkeypatch, _vcf(tmp_path, rows, indexed=indexed), data)
    assert TARGET in found, output
    assert found[TARGET].position == good and found[TARGET].carried_pair == ("G", "G")


@pytest.mark.parametrize("indexed", [False, True])
def test_identity_placement_obeys_selected_sample(tmp_path, monkeypatch, indexed):
    data = _data(tmp_path)
    path = _vcf(tmp_path, [_row(1000, gt="0/0", second="./.")], indexed=indexed, second=True)
    found, output = _extract(monkeypatch, path, data, sample="OTHER")
    assert TARGET in found and found[TARGET].is_no_call, output


@pytest.mark.parametrize(
    "gnomad,identity,name,expected",
    [
        (False, True, ".", True),
        (True, False, ".", False),
        (True, False, TARGET, True),
        (True, True, TARGET, True),
    ],
)
@pytest.mark.parametrize("indexed", [False, True])
def test_missing_cache_and_named_call_controls(
    tmp_path, monkeypatch, gnomad, identity, name, expected, indexed
):
    data = _data(tmp_path, gnomad=gnomad, identity=identity)
    path = _vcf(tmp_path, [_row(1000, name=name)], indexed=indexed)
    found, output = _extract(monkeypatch, path, data)
    assert (TARGET in found) == expected, output
    if expected:
        assert found[TARGET].carried_pair == ("G", "G")


@pytest.mark.parametrize("indexed", [False, True])
def test_existing_gnomad_placement_keeps_precedence(tmp_path, monkeypatch, indexed):
    data = _data(tmp_path, target_mapped=True)
    path = _vcf(tmp_path, [_row(900), _row(1000)], indexed=indexed)
    found, output = _extract(monkeypatch, path, data)
    assert TARGET in found and found[TARGET].position == 900, output


@pytest.mark.parametrize("build,pos", [("GRCh37", 200), ("GRCh38", 1000)])
@pytest.mark.parametrize("named", [False, True])
def test_shared_resolution_already_uses_identity_without_gnomad_placement(
    tmp_path, build, pos, named
):
    data = _data(tmp_path)
    call = Variant(
        TARGET if named else ".",
        "1",
        pos,
        "G",
        "A",
        build=build,
        ref=None if named else "G",
        ref_from_input=not named,
    )
    with GnomadAnnotator(data) as gnomad, RsidIdentityAnnotator(data) as identity:
        resolve_batch([call], clinvar=None, gnomad=gnomad, identity=identity)
    assert (call.rsid, call.ref, call.position, call.carried_pair) == (
        TARGET,
        "G",
        pos,
        ("G", "A"),
    )


@pytest.mark.parametrize("build,pos", [("GRCh37", 200), ("GRCh38", 1000)])
@pytest.mark.parametrize("pair", [("G", "A"), ("A", "A")])
def test_plink_missing_gnomad_placement_retains_observed_call(tmp_path, build, pos, pair):
    data = _data(tmp_path)
    call = Variant(TARGET, "1", pos, *pair, build=build, ref="G")
    if build == "GRCh37":
        coding = resolve_ref_alt_via_identity([call], data)
        assert coding
        assert resolve_ref_alt_via_gnomad_grch37([call], data, coding) == {}
    else:
        coding = resolve_ref_alt_via_gnomad_records([call], data)
        assert coding == {}
    prefix = tmp_path / "plink"
    counts = export_plink(iter([call]), prefix, build, record_alleles=coding)
    assert counts.written == 1 and counts.allele_mismatches == 0
    fields = prefix.with_suffix(".bim").read_text().split()
    code = prefix.with_suffix(".bed").read_bytes()[3] & 3
    assert code != 1
    decoded = {0: (fields[4], fields[4]), 2: (fields[4], fields[5]), 3: (fields[5], fields[5])}[
        code
    ]
    assert sorted(decoded) == sorted(pair)
    assert fields[1] == TARGET and int(fields[3]) == pos

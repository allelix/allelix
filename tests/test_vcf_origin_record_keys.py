"""O1: the finding's writer key must belong to its originating record too.

A locus, rsID and carried ALT can agree while the records' REFs differ.
Filtering by record number after assigning the first record's key loses the
other record's real finding. Exercise actual SNPedia annotations and the VCF
writer, with identity present/absent, single-base and deletion records, and batching.
"""

import pytest
from test_vcf_origin_isolation import RS, cache

from allelix.annotators.rsid_identity import RsidIdentityAnnotator
from allelix.annotators.snpedia import SNPediaAnnotator
from allelix.databases.rsid_identity_loader import write_identity_db
from allelix.parsers.vcf import VcfParser
from allelix.reports import _pipeline
from allelix.reports.vcf import render_vcf


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("batch", [1, 5000])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("identity_mode", ["absent-named", "present-named", "present-unnamed"])
@pytest.mark.parametrize("second_ref", ["G", "C", "GA"])
def test_writer_key_follows_the_originating_record(
    tmp_path, monkeypatch, build, batch, reverse, identity_mode, second_ref
):
    monkeypatch.setattr(_pipeline, "_BATCH_SIZE", batch)
    data = tmp_path / "data"
    cache(data, second=True)
    present = identity_mode.startswith("present")
    if present:
        pairs = [("G", "A")]
        if second_ref != "G":
            pairs.append((second_ref, "A"))
        write_identity_db(
            [(RS, "1", 1000, ref, alt, 1000, ref, alt, 0) for ref, alt in pairs],
            data / "rsid_identity.sqlite",
            version="synthetic",
            source_url="test://synthetic",
            downloaded_at="2026-11-03",
        )
    rsid = "." if identity_mode.endswith("unnamed") else RS
    records = [("G", "0/1", "snpedia_bad"), (second_ref, "1/1", "snpedia_good")]
    if reverse:
        records.reverse()
    src = tmp_path / "input.vcf"
    src.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        + "".join(f"1\t1000\t{rsid}\t{ref}\tA\t.\tPASS\t.\tGT\t{gt}\n" for ref, gt, _ in records)
    )
    result = _pipeline.run_analysis(
        src,
        VcfParser(),
        [SNPediaAnnotator(data)],
        identity=RsidIdentityAnnotator(data) if present else None,
        build_override=build,
        collect_vcf_write_inputs=True,
    )
    # Establish actual source findings before checking the writer boundary.
    # O4 keeps the slash for a deletion record, even when its called ALT is one base.
    hom_label = "A/A" if second_ref == "GA" else "AA"
    assert sorted(a.genotype_match for a in result.annotations) == [hom_label, "AG"]
    assert sorted(a.record_index for a in result.annotations) == [0, 1]
    assert not any(w.code == "identity_locus_conflict" for w in result.report_warnings)
    wi = result.vcf_write_inputs
    assert wi is not None
    out = tmp_path / "output.vcf"
    render_vcf(
        input_path=src,
        output_path=out,
        build=result.build,
        annotators_used=result.annotators_used,
        annotations_by_key=wi.annotations_by_key,
        recovered_rsids=wi.recovered_rsids,
        reference_rows=wi.reference_rows,
        record_rows=wi.record_rows,
        record_recovered=wi.record_recovered,
    )
    lines = [line.split("\t") for line in out.read_text().splitlines() if not line.startswith("#")]
    assert len(lines) == 2
    for cols, (ref, gt, significance) in zip(lines, records, strict=True):
        assert (cols[3], cols[4], cols[9]) == (ref, "A", gt)
        assert cols[2] == RS
        info = dict(token.split("=", 1) for token in cols[7].split(";") if "=" in token)
        assert info.get("ALLELIX_SNPEDIA") == f"{RS}|{significance}", (
            "a finding must not be lost under another record's REF/ALT key"
        )

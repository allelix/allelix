"""Opt-in real-bundle gate (``-m bundle``); aggregate assertions only, never downloads data."""

import json
import os
from pathlib import Path

import pytest
from click.testing import CliRunner

from allelix.cli import main
from allelix.parsers.vcf import VcfParser, _canonical_rsid


@pytest.mark.slow
@pytest.mark.bundle
def test_hg00096_chr22_has_eight_pairs_and_58102_calls(
    tmp_path, all_annotators_data_dir, monkeypatch
):
    default = Path(__file__).parents[1] / "test_data/real/vcf/thousandG_chr22.vcf.gz"
    source = Path(os.environ.get("ALLELIX_THOUSANDG_VCF", str(default)))
    assert source.is_file(), (
        "Mount the real bundle or set ALLELIX_THOUSANDG_VCF; this gate cannot skip"
    )
    data = all_annotators_data_dir
    (data / "config.toml").write_text(
        "[sources]\nclinvar=true\npharmgkb=false\ngwas=false\nsnpedia=false\n"
        "gnomad=false\nalphamissense=false\ncadd=false\nrsid_identity=false\n"
    )
    original = VcfParser._split_site_calls
    counts = {"merged": 0, "positional_pairs": 0}
    isolated = []

    def observe(self, site, sample_idx):
        calls = list(original(self, site, sample_idx))
        counts["merged"] += sum(len(v.allele_records or ()) == 2 for v in calls)
        if len(site) == 2:
            c1, c2 = site[0][2], site[1][2]
            fields = [
                dict(zip(c[8].split(":"), c[9 + sample_idx].split(":"), strict=False))
                for c in (c1, c2)
            ]
            if (
                sorted(f.get("GT", "") for f in fields) == ["0|1", "1|0"]
                and c1[3] == c2[3]
                and c1[4] != c2[4]
                and fields[0].get("PS") == fields[1].get("PS")
                and c1[6] == c2[6] == "PASS"
                and c1[2] != c2[2]
                and not _canonical_rsid(c1[2])
                and not _canonical_rsid(c2[2])
            ):
                counts["positional_pairs"] += 1
                for c in (c1, c2):
                    # Only the selected public sample; retained temporarily for
                    # the same-record ID='.' control. Never written to reports.
                    isolated.append([*c[:2], ".", *c[3:9], c[9 + sample_idx]])
        yield from calls

    monkeypatch.setattr(VcfParser, "_split_site_calls", observe)
    output = tmp_path / "summary.json"
    result = CliRunner().invoke(
        main,
        [
            "analyze",
            str(source),
            "--sample",
            "HG00096",
            "--build",
            "grch38",
            "--data-dir",
            str(data),
            "--no-update",
            "--output",
            str(output),
        ],
    )
    assert result.exit_code == 0, (result.output, type(result.exception).__name__)
    report = json.loads(output.read_text())
    actual = (
        report["input"]["total_variants"],
        counts["merged"],
        sum(w["count"] for w in report["report_warnings"] if w["code"] == "split_site_calls"),
    )
    assert counts["positional_pairs"] == 8, counts
    # Verify Addendum section 1's dot-ID control independently of the failing
    # aggregate assertion. All other per-record/sample fields are preserved.
    monkeypatch.setattr(VcfParser, "_split_site_calls", original)
    control = tmp_path / "dot-control.vcf"
    control.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=22,length=50818468>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG00096\n"
        + "".join("\t".join(c) + "\n" for c in isolated)
    )
    parser = VcfParser(sample="HG00096")
    parser.combine_split_sites = True
    controls = list(parser.parse(control))
    control_counts = (
        len(controls),
        sum(len(v.allele_records or ()) == 2 for v in controls),
        parser.split_site_calls,
    )
    print(
        json.dumps(
            {
                "original": actual,
                "dot_id_control": control_counts,
                "positional_pairs": counts["positional_pairs"],
            }
        )
    )
    assert control_counts == (8, 8, 0)
    assert actual == (58102, 8, 0), "calls read / merged calls / records counted as non-pairs"

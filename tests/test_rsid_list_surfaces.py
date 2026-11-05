"""rsID and genotype lists (ADR-0069): the placement lines, the build text, other commands stop."""

import json

import pytest
from click.testing import CliRunner
from test_annotation_parity_and_named_disclosures import _identity, _pharm

from allelix.annotators._placement import PlacementStats
from allelix.cli import main
from allelix.reports import placement_lines
from allelix.reports._pipeline import BuildDiagnostics
from allelix.reports.build_text import short_source, source_text
from allelix.utils.build_detect import SOURCE_RSID_PLACEMENT, decide_placement_build

ROWS = "\n".join(f"rs39999801{i:02d}\tAG" for i in range(5)) + "\n"


def test_no_lines_for_a_file_with_positions():
    assert placement_lines(None) == []


def test_lines_count_placed_complemented_and_each_unplaced_reason():
    stats = PlacementStats(build="GRCh38", rows=10, placed=6, complemented=2)
    stats.unplaced.update(not_in_table=1, both_fit=2, no_fit=1)
    first, second = placement_lines(stats)
    assert first == (
        "The file has no positions: 6 of 10 rows were placed on GRCh38 from the rsID "
        "identity table, 2 of them complemented to the forward strand."
    )
    assert second.startswith("4 row(s) could not be placed (1 rsID not placed once")
    assert "2 call fits both strands at the site; 1 call fits neither strand" in second
    assert "palindromic" not in second and "D/I" not in second
    assert "SNPedia, ClinPGx and the GWAS Catalog still look them up by rsID" in second


def test_every_placed_row_gives_one_line():
    assert placement_lines(PlacementStats(build="GRCh38", rows=3, placed=3)) == [
        "The file has no positions: 3 of 3 rows were placed on GRCh38 from the rsID "
        "identity table."
    ]


def test_build_text_names_the_identity_table():
    decision = decide_placement_build(file_name="x.txt", override=None)
    assert decision.source == SOURCE_RSID_PLACEMENT
    diag = BuildDiagnostics.from_decision(decision)
    assert source_text(diag) == (
        "the file has no positions; its rsIDs are placed from the rsID identity table"
    )
    assert short_source(diag) == (
        "rsIDs placed from the rsID identity table; the file has no positions"
    )


@pytest.mark.parametrize(
    "command",
    (
        ["stats", "{path}"],
        ["stats", "{path}", "--build", "grch38"],
        ["export", "plink", "{path}", "-o", "{out}"],
        ["extract", "{path}", "--snps", "rs3999980100"],
    ),
    ids=("stats", "stats-grch38", "export-plink", "extract"),
)
def test_commands_without_placement_stop_on_a_position_free_file(tmp_path, command):
    path = tmp_path / "list.txt"
    path.write_text(ROWS)
    argv = [a.format(path=path, out=tmp_path / "out") for a in command]
    result = CliRunner().invoke(main, [*argv, "--data-dir", str(tmp_path / "data")])
    assert result.exit_code != 0
    assert "has no positions" in result.output
    assert "Only analyze, methylation and pharmacogenomics read it" in result.output
    assert not list(tmp_path.glob("out*"))


def test_compare_stops_on_a_position_free_file(tmp_path):
    path = tmp_path / "list.txt"
    path.write_text(ROWS)
    result = CliRunner().invoke(main, ["compare", str(path), str(path)])
    assert result.exit_code != 0
    assert "has no positions" in result.output
    assert "--build1" not in result.output


def test_analyze_reports_placement_on_every_surface(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    _identity(data, [("rs3999980100", "1", 900, "A", "G", 1000, "A", "G", 0)])
    _pharm(data, "rs3999980100", "MOCK", ("AG",))
    path = tmp_path / "list.txt"
    path.write_text(ROWS)
    runs = {}
    for suffix in ("html", "json"):
        out = tmp_path / f"report.{suffix}"
        result = CliRunner().invoke(
            main,
            ["analyze", str(path), "--data-dir", str(data), "--no-update", "--output", str(out)],
        )
        assert result.exit_code == 0, result.output
        runs[suffix] = (" ".join(result.output.split()), out)
    console = runs["html"][0]
    assert "The file has no positions: 1 of 5 rows were placed on GRCh38" in console
    assert "4 row(s) could not be placed (4 rsID not placed once on GRCh38" in console
    html = " ".join(runs["html"][1].read_text().split())
    assert "The file has no positions: 1 of 5 rows were placed on GRCh38" in html
    report = json.loads(runs["json"][1].read_text())
    assert report["input"]["build_source"] == SOURCE_RSID_PLACEMENT
    assert report["input"]["placement"]["unplaced"]["not_in_table"] == 4

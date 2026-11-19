"""A withheld-site reason belongs to the selected call, not every call of its rsID."""

import sqlite3
from contextlib import closing

import pytest
from click.testing import CliRunner
from test_extract_duplicate_placements import RSID, _cache, _input, _row

from allelix.cli import main, utility


@pytest.mark.parametrize("indexed", [False, True], ids=["sequential", "indexed"])
@pytest.mark.parametrize("selected", ["missing", "hom_ref", "heterozygote", "withheld"])
def test_withheld_note_is_only_for_the_selected_candidate(
    tmp_path, monkeypatch, indexed, selected
):
    data = _cache(tmp_path, 100)
    with closing(sqlite3.connect(data / "gnomad.sqlite")) as conn, conn:
        conn.execute(
            "INSERT INTO gnomad_frequencies (chrom,pos,ref,alt,rsid,af) VALUES (?,?,?,?,?,?)",
            ("1", 300, "G", "A", RSID, 0.01),
        )
    held_position = 100 if selected == "withheld" else 300
    ordinary_position = 300 if selected == "withheld" else 100
    gt = {"missing": "./.", "hom_ref": "0/0", "heterozygote": "0/1", "withheld": "./."}[selected]
    rows = [
        (ordinary_position, _row(ordinary_position, gt)),
        (held_position, _row(held_position, "1/1")),
        (held_position, _row(held_position, "0/1").replace("\tG\tA\t", "\tG\tT\t")),
    ]
    source = _input(tmp_path, [line for _pos, line in sorted(rows, key=lambda row: row[0])])
    if not indexed:
        source.with_suffix(source.suffix + ".tbi").unlink()
    captured = []
    monkeypatch.setattr(
        utility, "_render_extract_table", lambda _p, _w, found: captured.append(found)
    )
    result = CliRunner().invoke(
        main,
        ["extract", str(source), "--snps", RSID, "--data-dir", str(data), "--build", "grch38"],
    )
    assert result.exit_code == 0, (result.output, result.exception)
    assert captured[0][RSID].position == 100
    text = " ".join(result.output.split())
    if selected == "withheld":
        assert captured[0][RSID].is_no_call
        assert "contradict" in text
    else:
        assert "contradict" not in text, (
            "The chosen call at 100 is ordinary; the contradictory records at 300 "
            "must not lend it their reason"
        )

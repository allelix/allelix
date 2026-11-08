# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Pin the gene and the words of a D/I call named by position.

The line says how the naming row reaches the call: a variant at its
position, an indel whose REF span covers it, or one whose span ends a
stated number of bases before it (within the five-base allowance).

A D/I call under an rsID ClinVar holds without an aggregate-P/LP row keeps
that rsID's own ClinVar gene, and its line says which Pathogenic row named
it. A vendor probe's or an unheld rsID's call keeps the naming row's gene
and says the same. A P/LP rsID's own D/I call keeps its text. Which calls
are named, and every count, are unchanged.
"""

import json

import pytest
from click.testing import CliRunner
from test_outside_review_rc1_b import BUILDS, _analyze, _chip, _clinvar, _config, _identity

from allelix.cli import main

POS = 2_489_000
PLP = "rs3999489001"
HELD = "rs3999489002"
UNHELD = "rs3999489003"
PROBE = "i3999489004"
COVERS = "named because a Pathogenic or Likely pathogenic indel in GENEP covers its position"
AT = "named because a Pathogenic or Likely pathogenic variant in GENEP is at its position"
NOT_ITSELF = "ClinVar does not classify this rsID itself as Pathogenic or Likely pathogenic"


def _entries(tmp_path, build, rows, calls, *, show_unchecked=False):
    data = tmp_path / "data"
    _config(data)
    _identity(data, [])
    _clinvar(data, build, rows)
    source = _chip(tmp_path / "chip.txt", calls, build=build)
    if show_unchecked:
        out = tmp_path / "report.json"
        result = CliRunner(env={"COLUMNS": "500"}).invoke(
            main,
            [
                "analyze",
                str(source),
                "--data-dir",
                str(data),
                "--no-update",
                "--output",
                str(out),
                "--build",
                build.lower(),
                "--show-unchecked",
            ],
        )
        assert result.exit_code == 0, (result.output, result.exception)
        report, console = json.loads(out.read_text()), result.output
    else:
        report, console = _analyze(tmp_path, data, source, build)
    return {e["rsid"]: e for e in report["high_value_unchecked"]}, report, console


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("own_gene", ["GENEB", ""], ids=["own-gene", "no-own-gene"])
def test_held_rsid_keeps_its_own_gene_and_says_why(tmp_path, build, own_gene):
    rows = [
        (PLP, "13", POS, "GAT", "G", "GENEP", "Pathogenic", "Pathogenic"),
        (HELD, "13", POS + 1, "AT", "A", own_gene, "Benign", "Benign"),
    ]
    named, report, console = _entries(
        tmp_path, build, rows, [(HELD, "13", POS + 1, "II")], show_unchecked=True
    )
    entry = named[HELD]
    assert (entry["reason"], entry["gene"]) == ("di_call", own_gene)
    assert COVERS in entry["message"] and NOT_ITSELF in entry["message"]
    assert COVERS in console
    warning = {w["code"]: w for w in report["report_warnings"]}["di_calls_not_checked"]
    assert warning["count"] == 1 and COVERS in warning["message"]


@pytest.mark.parametrize("build", BUILDS)
def test_unheld_rsid_and_probe_keep_the_naming_rows_gene(tmp_path, build):
    rows = [
        (PLP, "13", POS, "GAT", "G", "GENEP", "Pathogenic", "Pathogenic"),
        ("rs3999489005", "13", POS + 40, "C", "T", "GENEP", "Pathogenic", "Pathogenic"),
    ]
    calls = [(UNHELD, "13", POS + 2, "DD"), (PROBE, "13", POS + 40, "DI")]
    named, report, _ = _entries(tmp_path, build, rows, calls)
    assert named[UNHELD]["gene"] == named[PROBE]["gene"] == "GENEP"
    # ClinVar does not hold UNHELD, so its line says nothing of ClinVar classifying it.
    assert COVERS in named[UNHELD]["message"] and NOT_ITSELF not in named[UNHELD]["message"]
    assert AT in named[PROBE]["message"] and "at a vendor probe" in named[PROBE]["message"]
    assert NOT_ITSELF not in named[PROBE]["message"]
    assert report["input"]["probe_classes"]["di_calls_at_plp_rsids"] == 0


@pytest.mark.parametrize("build", BUILDS)
def test_plp_rsid_d_i_call_keeps_todays_text(tmp_path, build):
    rows = [(PLP, "13", POS, "GAT", "G", "GENEP", "Pathogenic", "Pathogenic")]
    named, report, _ = _entries(tmp_path, build, rows, [(PLP, "13", POS, "DI")])
    entry = named[PLP]
    assert (entry["reason"], entry["gene"]) == ("di_call", "GENEP")
    assert entry["message"] == (
        f"GENEP ({PLP}): called D/I, a deletion or insertion code that Allelix cannot "
        "check against its sources. Not assessed: this is neither a finding nor a "
        "normal result."
    )
    assert report["input"]["probe_classes"]["di_calls_at_plp_rsids"] == 1


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize(
    "ref,alt,offset,words",
    [
        pytest.param(
            "GAT", "G", 1, "indel in GENEP covers its position", id="deletion-inside-span"
        ),
        pytest.param(
            "GAT", "G", 2, "indel in GENEP covers its position", id="deletion-last-ref-base"
        ),
        pytest.param(
            "GAT", "G", 3, "indel in GENEP ends 1 base before its position", id="deletion-plus-1"
        ),
        pytest.param(
            "GAT", "G", 7, "indel in GENEP ends 5 bases before its position", id="deletion-plus-5"
        ),
        pytest.param(
            "A", "AC", 4, "indel in GENEP ends 4 bases before its position", id="insertion-plus-4"
        ),
        pytest.param("GAT", "G", 0, "variant in GENEP is at its position", id="at-position"),
    ],
)
@pytest.mark.parametrize("call_id", [PROBE, UNHELD, HELD], ids=["probe", "unheld", "held"])
def test_line_says_how_the_naming_row_reaches_the_call(
    tmp_path, build, ref, alt, offset, words, call_id
):
    rows = [(PLP, "13", POS, ref, alt, "GENEP", "Pathogenic", "Pathogenic")]
    if call_id == HELD:
        rows.append((HELD, "2", 800_000, "G", "A", "GENEB", "Benign", "Benign"))
    named, report, _ = _entries(tmp_path, build, rows, [(call_id, "13", POS + offset, "DI")])
    message = named[call_id]["message"]
    assert f"named because a Pathogenic or Likely pathogenic {words}" in message, message
    if not words.endswith("covers its position"):
        assert "covers" not in message, message
    assert (NOT_ITSELF in message) == (call_id == HELD), message
    assert named[call_id]["gene"] == ("GENEB" if call_id == HELD else "GENEP")
    assert report["input"]["probe_classes"]["di_calls_at_plp_rsids"] == 0

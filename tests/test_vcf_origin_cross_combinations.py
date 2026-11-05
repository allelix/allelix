"""The annotated VCF's record ownership against its other per-record rules (ADR-0059).

A record takes only the findings its own call produced (``record_index``).
These cases cross that rule with per-ALT slots, explicit reference records,
recovered rsIDs, findings without a record number and the writer's mode
without record rows.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from allelix.models import Annotation
from allelix.parsers.vcf import reference_row_key
from allelix.reports.vcf import render_vcf

if TYPE_CHECKING:
    from pathlib import Path

_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=1>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
)


def _finding(alt: str, origin: int | None, label: str, *, reference: bool = False) -> Annotation:
    return Annotation(
        source="clinvar",
        rsid="rs900000001",
        significance=label,
        category="clinical",
        magnitude=6,
        description=label,
        attribution="Synthetic",
        genotype_match="AG",
        gene="GENE",
        alt=alt,
        locus=("1", 1000),
        locus_ref="G",
        record_index=origin,
        reference_carries_alt=reference,
    )


def _write(tmp_path: Path, rows: list[str], **kwargs: object) -> list[list[str]]:
    src = tmp_path / "input.vcf"
    src.write_text(_HEADER + "".join(f"{r}\n" for r in rows))
    out = tmp_path / "output.vcf"
    defaults: dict[str, object] = {
        "recovered_rsids": {},
        "reference_rows": {},
        "record_recovered": {},
    }
    render_vcf(
        input_path=src,
        output_path=out,
        build="GRCh38",
        annotators_used=[("clinvar", "synthetic")],
        **{**defaults, **kwargs},  # type: ignore[arg-type]
    )
    return [line.split("\t") for line in out.read_text().splitlines() if not line.startswith("#")]


def _clinvar(cols: list[str]) -> str | None:
    for token in cols[7].split(";"):
        if token.startswith("ALLELIX_CLINVAR="):
            return token.split("=", 1)[1]
    return None


def test_each_alt_slot_takes_its_own_records_finding(tmp_path: Path) -> None:
    """Two records share ALT A; the second also carries C. Each slot is per record."""
    a, c = ("1", 1000, "G", "A"), ("1", 1000, "G", "C")
    rows = [
        "1\t1000\trs900000001\tG\tA\t.\tPASS\t.\tGT\t0/1",
        "1\t1000\trs900000001\tG\tA,C\t.\tPASS\t.\tGT\t1/2",
    ]
    out = _write(
        tmp_path,
        rows,
        annotations_by_key={a: [_finding("A", 0, "a0")], c: [_finding("C", 1, "c1")]},
        record_rows={a: {0, 1}, c: {1}},
    )
    assert _clinvar(out[0]) == "a0|rs900000001|GENE"
    # Record 1 carries A but did not produce the A finding: its A slot is empty.
    assert _clinvar(out[1]) == ".,c1|rs900000001|GENE"


def test_explicit_reference_records_keep_their_rule(tmp_path: Path) -> None:
    """``reference_rows`` stamps every listed explicit record, whatever its number."""
    rows = [
        "1\t1000\t.\tG\t.\t.\tPASS\t.\tGT\t0/0",
        "1\t1000\t.\tG\t.\t.\tPASS\t.\tGT\t0/0",
    ]
    key = reference_row_key(rows[0].split("\t"))
    out = _write(
        tmp_path,
        rows,
        annotations_by_key={("1", 1000, "G", "G"): [_finding("G", 0, "ref0", reference=True)]},
        record_rows={},
        reference_rows={key: {("1", 1000)}},
    )
    assert all("ALLELIX_CLINVAR_REF=ref0|rs900000001|GENE" in cols[7] for cols in out)


def test_recovered_rsid_stays_with_its_record_without_a_finding(tmp_path: Path) -> None:
    """A record whose rsID was recovered keeps it when the finding is another record's."""
    a = ("1", 1000, "G", "A")
    rows = [
        "1\t1000\t.\tG\tA\t.\tPASS\t.\tGT\t0/1",
        "1\t1000\t.\tG\tA\t.\tPASS\t.\tGT\t1/1",
    ]
    out = _write(
        tmp_path,
        rows,
        annotations_by_key={a: [_finding("A", 0, "a0")]},
        record_rows={a: {0, 1}},
        record_recovered={0: "rs900000001", 1: "rs900000001"},
    )
    assert [cols[2] for cols in out] == ["rs900000001", "rs900000001"]
    assert _clinvar(out[0]) is not None
    assert _clinvar(out[1]) is None
    assert "ALLELIX_ORIGINAL_ID=." in out[1][7]


def test_finding_without_a_record_number_is_not_stamped(tmp_path: Path) -> None:
    """With record rows, a finding no record owns goes on no record."""
    a = ("1", 1000, "G", "A")
    out = _write(
        tmp_path,
        ["1\t1000\trs900000001\tG\tA\t.\tPASS\t.\tGT\t0/1"],
        annotations_by_key={a: [_finding("A", None, "orphan")]},
        record_rows={a: {0}},
    )
    assert _clinvar(out[0]) is None


def test_without_record_rows_the_writer_keeps_its_key_lookup(tmp_path: Path) -> None:
    """A direct caller passing no record rows gets the key lookup, as before."""
    a = ("1", 1000, "G", "A")
    out = _write(
        tmp_path,
        ["1\t1000\trs900000001\tG\tA\t.\tPASS\t.\tGT\t0/1"] * 2,
        annotations_by_key={a: [_finding("A", 0, "a0")]},
        record_rows=None,
    )
    assert [_clinvar(cols) for cols in out] == ["a0|rs900000001|GENE"] * 2

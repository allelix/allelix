# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""High-value SNP no-call detection.

Loads a list of clinically important SNPs from a YAML data file and
scans parsed variants for no-calls on those positions. When a
high-value SNP is a no-call, the scanner produces a warning that
surfaces in ``stats``, ``analyze``, and focused reports.

The loader supports merging multiple YAML files (built-in first, then
user-provided overrides keyed by rsid) so custom high-value lists can
be layered on without modifying the shipped data.
"""

from __future__ import annotations

from dataclasses import dataclass
from importlib import resources
from typing import TYPE_CHECKING

import yaml

if TYPE_CHECKING:
    from pathlib import Path

    from allelix.models import Variant


@dataclass(frozen=True)
class HighValueSNP:
    """A clinically important SNP that warrants explicit no-call warnings."""

    rsid: str
    gene: str
    cluster: str
    note: str


@dataclass(frozen=True)
class NoCallWarning:
    """A high-value SNP that returned a no-call in the input file.

    ``reason`` (VCF input) says why: ``no_call`` (the caller emitted
    a missing genotype), ``filtered``, ``low_gq_call``, ``low_gq_block``,
    ``spanning_deletion``, ``symbolic_allele`` or ``no_coverage`` (gVCF:
    no record covers the position). Empty for array input.
    """

    snp: HighValueSNP
    cluster_incomplete: bool
    reason: str = ""


# Plain-language text for each no-call reason.
_REASON_TEXT: dict[str, str] = {
    "no_call": "no genotype called",
    "filtered": "rejected by the variant caller's filters",
    "low_gq_call": "genotype quality below the floor",
    "low_gq_block": "reference call with low genotype quality or no reads",
    "spanning_deletion": "overlapped by a deletion on one copy",
    "symbolic_allele": "structural or symbolic allele",
    "no_coverage": "no record covers this position",
}


def load_high_value_snps(
    extra_paths: list[Path] | None = None,
) -> dict[str, HighValueSNP]:
    """Load the built-in high-value SNP list, optionally merging user files.

    Returns a dict keyed by rsid. User-provided files override built-in
    entries with the same rsid.
    """
    result: dict[str, HighValueSNP] = {}
    builtin = resources.files("allelix.data").joinpath("high_value_snps.yaml")
    with resources.as_file(builtin) as p:
        result.update(_load_yaml(p))
    for path in extra_paths or []:
        result.update(_load_yaml(path))
    return result


def _load_yaml(path: Path) -> dict[str, HighValueSNP]:
    """Parse a single YAML file into HighValueSNP entries."""
    try:
        with open(path, encoding="utf-8") as fh:
            entries = yaml.safe_load(fh) or []
    except (OSError, yaml.YAMLError) as exc:
        msg = f"Failed to read high-value SNP file {path}: {exc}"
        raise ValueError(msg) from exc
    if not isinstance(entries, list):
        msg = f"Expected a YAML list in {path}, got {type(entries).__name__}"
        raise ValueError(msg)
    out: dict[str, HighValueSNP] = {}
    for i, entry in enumerate(entries):
        if not isinstance(entry, dict) or "rsid" not in entry:
            msg = f"Entry {i} in {path} is missing required 'rsid' field"
            raise ValueError(msg)
        rsid = entry["rsid"]
        out[rsid] = HighValueSNP(
            rsid=rsid,
            gene=entry.get("gene", ""),
            cluster=entry.get("cluster", ""),
            note=entry.get("note", ""),
        )
    return out


def scan_no_calls(
    variants: list[Variant],
    high_value: dict[str, HighValueSNP] | None = None,
    vcf_states: dict[str, tuple[str, str]] | None = None,
) -> list[NoCallWarning]:
    """Scan variants for no-calls on high-value SNPs.

    Returns a list of ``NoCallWarning`` for each high-value SNP that is
    a no-call. Cluster-incomplete detection: if any member of a cluster
    is a no-call while another member was called, the warning notes that
    the cluster result is incomplete (e.g., "APOE genotype cannot be
    determined").

    ``vcf_states`` (VCF input) adds what the file says at each
    high-value position, including rsID-less reference blocks:
    ``(called | reference | no_call, reason)``. A confident reference
    call counts as called. A variant call always wins over a no-call.
    """
    if high_value is None:
        high_value = load_high_value_snps()

    no_call_reasons: dict[str, str] = {}
    called_rsids: set[str] = set()
    for v in variants:
        if v.rsid in high_value:
            if v.is_no_call:
                no_call_reasons.setdefault(v.rsid, "")
            else:
                called_rsids.add(v.rsid)
    for rsid, (state, reason) in (vcf_states or {}).items():
        if rsid not in high_value:
            continue
        if state == "no_call":
            no_call_reasons[rsid] = reason or no_call_reasons.get(rsid, "")
        else:
            called_rsids.add(rsid)
    no_call_rsids = set(no_call_reasons) - called_rsids

    warnings: list[NoCallWarning] = []
    for rsid in sorted(no_call_rsids):
        snp = high_value[rsid]
        cluster_incomplete = False
        if snp.cluster:
            cluster_members = {r for r, s in high_value.items() if s.cluster == snp.cluster}
            cluster_incomplete = bool(cluster_members & called_rsids)
        warnings.append(
            NoCallWarning(
                snp=snp, cluster_incomplete=cluster_incomplete, reason=no_call_reasons[rsid]
            )
        )
    return warnings


def format_warnings(warnings: list[NoCallWarning]) -> list[str]:
    """Format no-call warnings as human-readable strings."""
    lines: list[str] = []
    for w in warnings:
        line = f"{w.snp.rsid} ({w.snp.gene}): {w.snp.note}"
        if w.reason:
            line += f" [{_REASON_TEXT.get(w.reason, w.reason)}]"
        if w.cluster_incomplete:
            line += f" — {w.snp.cluster} genotype cannot be determined"
        lines.append(line)
    return lines

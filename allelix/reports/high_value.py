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
    from collections.abc import Iterable
    from pathlib import Path

    from allelix.models import Variant


@dataclass(frozen=True)
class HighValueSNP:
    """A clinically important SNP that warrants explicit no-call warnings."""

    rsid: str
    gene: str
    cluster: str
    note: str
    # The variant's common name (``6174delT``, ``F508del``), when it has one.
    name: str = ""


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
            name=entry.get("name", ""),
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


# The report warning that names D/I calls at high-value or Pathogenic sites.
DI_CALLS_NOT_CHECKED = "di_calls_not_checked"

# Why a high-value or Pathogenic site could not be checked.
UNCHECKED_DI_CALL = "di_call"
UNCHECKED_NOT_IN_FILE = "reference_carried_not_in_file"
# 2.6.1: calls the counts already covered, now also named.
UNCHECKED_UNPLACED = "unplaced"
UNCHECKED_IDENTITY_CONFLICT = "identity_locus_conflict"
UNCHECKED_FOREIGN_ALLELES = "foreign_alleles"
UNCHECKED_VENDOR_PROBE = "vendor_probe"
UNCHECKED_UNKNOWN_ALLELES = "unknown_alleles"
# ``detail`` of a foreign- or unknown-allele call when another rsID's ClinVar
# row stands at its position.
OTHER_RSID_AT_LOCUS = "other_rsid"


@dataclass(frozen=True)
class UncheckedCall:
    """A call the run could not check, before the high-value / P/LP filter.

    ``site`` is the call's rsID or vendor probe ID, ``call`` the call as the
    file wrote it, ``reason`` one of the ``UNCHECKED_*`` reasons, and
    ``detail`` the words that say more (an unplaced row's reason, a
    conflict's position, the row that names a D/I call by position).
    ``plp_gene`` is ClinVar's gene where ClinVar classifies the site
    Pathogenic or Likely pathogenic, None where it does not; a call is named
    when it is P/LP or its rsID is high-value. ``own_plp`` (2.6.2) says
    whether the call's own rsID has a P/LP row, not only a row at its locus.
    """

    site: str
    call: str
    reason: str
    detail: str = ""
    plp_gene: str | None = None
    own_plp: bool = False


@dataclass(frozen=True)
class UncheckedSite:
    """A high-value or ClinVar Pathogenic/Likely pathogenic site the run could not check.

    ``reason`` is ``di_call`` (an array call written as D or I, a deletion
    or insertion code no source can be matched to; ``call`` is the call as
    written), ``reference_carried_not_in_file`` (a plain VCF with no
    record where the reference genome carries the ClinVar variant allele,
    so a person homozygous for it has no record either), or since 2.6.1
    ``unplaced`` (a row of an rsID and genotype list that placement could
    not place), ``identity_locus_conflict`` (a call at a position the rsID
    identity table does not give for its rsID) or ``foreign_alleles`` (a
    chip call naming a base that is not one of its own rsID's alleles there)
    or ``vendor_probe`` (a base call at a vendor probe without an rsID
    where ClinVar classifies the variant Pathogenic or Likely pathogenic), or
    ``unknown_alleles`` (an rsID call whose alleles at its position no source
    knows).
    ``rsid`` is the probe ID (``i…``) for a call at a vendor probe.
    ``high_value`` (2.6.2) says the rsID is on the high-value list and
    ``own_plp`` that ClinVar classifies the rsID itself P/LP; they decide
    whether the reports name the site or count it (``named_inline``).
    """

    rsid: str
    gene: str
    name: str
    reason: str
    call: str = ""
    detail: str = ""
    high_value: bool = False
    own_plp: bool = False

    @property
    def named_inline(self) -> bool:
        """Whether the terminal and HTML reports name the site (2.6.2).

        A high-value site always is; any other site when ClinVar classifies
        its own rsID P/LP and the call is not a D or I code. The rest are
        counted by kind (``unchecked_rollup``).
        """
        return self.high_value or (self.own_plp and self.reason != UNCHECKED_DI_CALL)

    @property
    def rollup_group(self) -> str:
        """The count a site not named inline goes into (``_ROLLUP``)."""
        if self.reason == UNCHECKED_DI_CALL:
            if not self.rsid.startswith("rs"):
                return "probe_di"
            return "di_plp" if self.own_plp else "di_near"
        if self.reason == UNCHECKED_VENDOR_PROBE:
            return "probe_base"
        if self.reason == UNCHECKED_UNKNOWN_ALLELES:
            return "unknown"
        if self.reason == UNCHECKED_FOREIGN_ALLELES:
            return "foreign"
        return self.reason

    @property
    def text(self) -> str:
        """The site and why it was not checked, in one line."""
        label = " ".join(part for part in (self.gene, self.name) if part)
        site = f"{label} ({self.rsid})" if label else self.rsid
        not_assessed = "Not assessed: this is neither a finding nor a normal result."
        if self.reason == UNCHECKED_DI_CALL and not self.rsid.startswith("rs"):
            # A vendor probe: only ClinVar's side is said; SNPedia may read the
            # probe's own page.
            if self.detail:  # named by the P/LP row at or near its position
                return (
                    f"{site}: called {self.call}, a deletion or insertion code at a vendor "
                    f"probe, named because {self.detail}. ClinVar does not match a probe "
                    "to that variant, so ClinVar did not check it."
                )
            return (
                f"{site}: called {self.call}, a deletion or insertion code at a vendor "
                "probe. ClinVar does not match a probe to the variant at its position, "
                "so ClinVar did not check it."
            )
        if self.reason == UNCHECKED_VENDOR_PROBE:
            return (
                f"{site}: called {self.call} at a vendor probe without an rsID. ClinVar "
                "does not match a probe to the variant at its position (nothing confirms "
                "the probe tests it), so ClinVar did not check it."
            )
        if self.reason == UNCHECKED_DI_CALL and self.detail:
            # an rsID with no P/LP row, named by position. The detail says
            # which row named it, and for an rsID ClinVar holds, that ClinVar
            # does not classify the rsID itself as P/LP.
            return (
                f"{site}: called {self.call}, a deletion or insertion code that Allelix "
                f"cannot check against its sources, named because {self.detail}. "
                f"{not_assessed}"
            )
        if self.reason == UNCHECKED_DI_CALL:
            return (
                f"{site}: called {self.call}, a deletion or insertion code that Allelix "
                f"cannot check against its sources. {not_assessed}"
            )
        if self.reason == UNCHECKED_UNPLACED:
            return (
                f"{site}: called {self.call}, a row of the rsID and genotype list that "
                f"could not be placed ({self.detail}), so ClinVar could not check it. "
                "This does not affect a source that looks a call up by its rsID."
            )
        if self.reason == UNCHECKED_IDENTITY_CONFLICT:
            return (
                f"{site}: called {self.call} at {self.detail}, a position the rsID "
                "identity table does not give for this rsID, so no source read it. "
                f"{not_assessed}"
            )
        other = self.detail == OTHER_RSID_AT_LOCUS
        if self.reason == UNCHECKED_FOREIGN_ALLELES:
            return (
                f"{site}: called {self.call}, a base that is not one of its own rsID's "
                "alleles at this position (in the rsID identity table, ClinVar or gnomAD), "
                + (
                    "so ClinVar did not match it to another rsID's row there, which could "
                    "have been taken for it."
                    if other
                    else "so ClinVar compared it only with its own rsID's rows."
                )
            )
        if self.reason == UNCHECKED_UNKNOWN_ALLELES:
            return (
                f"{site}: called {self.call}. No source knows this rsID's alleles at its "
                "position, "
                + (
                    "so ClinVar did not match it to another rsID's row there, which could "
                    "have been taken for it."
                    if other
                    else "so ClinVar did not check it against the variant there."
                )
            )
        return (
            f"{site}: no record in this VCF. The reference genome carries the ClinVar "
            "variant allele here, so a person homozygous for it has no record either; "
            "not assessed. A gVCF or an all-sites VCF says which."
        )

    def to_dict(self) -> dict[str, str]:
        """The JSON report's ``high_value_unchecked`` entry."""
        return {
            "rsid": self.rsid,
            "gene": self.gene,
            "name": self.name,
            "reason": self.reason,
            "call": self.call,
            "message": self.text,
        }


# 2.6.2: the counts of the sites the reports do not name, in this order, each
# with why the run could not check them. A count is of sites (one per ID and
# reason, as the header counts them), not of calls: an ID on two rows is one.
_ROLLUP: dict[str, str] = {
    "di_plp": (
        "site(s) at rsIDs ClinVar classifies Pathogenic or Likely pathogenic with a call "
        "written as D or I, a deletion or insertion code Allelix cannot check against its "
        "sources. Not assessed: neither a finding nor a normal result."
    ),
    "di_near": (
        "site(s) at rsIDs with no Pathogenic or Likely pathogenic row of their own with a "
        "call written as D or I, at or near a Pathogenic or Likely pathogenic variant. "
        "Allelix cannot check a deletion or insertion code against its sources."
    ),
    "probe_di": (
        "site(s) at vendor probes without an rsID with a call written as D or I, at or "
        "near a Pathogenic or Likely pathogenic variant. ClinVar does not match a probe to "
        "a variant, so ClinVar did not check them."
    ),
    "probe_base": (
        "site(s) at vendor probes without an rsID with a base call, at a Pathogenic or "
        "Likely pathogenic position. ClinVar does not match a probe to the variant there, "
        "so ClinVar did not check them."
    ),
    "unknown": (
        "site(s) at a Pathogenic or Likely pathogenic position with a call whose alleles "
        "no source knows for its rsID, so ClinVar did not check them against the variant "
        "there."
    ),
    "foreign": (
        "site(s) at a Pathogenic or Likely pathogenic position with a call naming a base "
        "that is not one of its own rsID's alleles, so ClinVar did not check them against "
        "that variant."
    ),
}


def unchecked_rollup(sites: list[UncheckedSite]) -> tuple[list[UncheckedSite], list[str]]:
    """The sites the reports name, and one count line for each kind of the rest.

    2.6.2: a list of every unchecked call can run to thousands of lines ahead
    of the findings, so the terminal and HTML reports name only the sites
    ``named_inline`` keeps; the JSON report and ``--show-unchecked`` keep
    every site.
    """
    inline = [site for site in sites if site.named_inline]
    counts: dict[str, int] = {}
    for site in sites:
        if not site.named_inline:
            counts[site.rollup_group] = counts.get(site.rollup_group, 0) + 1
    order = list(_ROLLUP)
    lines = [
        f"{counts[key]:,} {_ROLLUP.get(key, f'call(s) not checked ({key}).')}"
        for key in sorted(counts, key=lambda k: (order.index(k) if k in order else len(order), k))
    ]
    return inline, lines


def unchecked_sites(
    high_value: dict[str, HighValueSNP],
    variants: list[Variant],
    di_calls_at_plp: dict[str, tuple[str, str]],
    not_in_file: list[str],
    other_calls: Iterable[UncheckedCall] = (),
) -> list[UncheckedSite]:
    """The high-value and ClinVar P/LP sites that could not be checked, by rsID.

    A D/I array call is named at a high-value rsID (``variants``) or at an
    rsID ClinVar classifies Pathogenic or Likely pathogenic
    (``di_calls_at_plp``: rsID -> call, gene). A reference-carried site with
    no record in a plain VCF (``not_in_file``) is named when it is high-value.
    Each of ``other_calls`` is named when ClinVar classifies its site P/LP
    or its rsID is high-value, once per site and reason.
    """
    from allelix.models import DI_ALLELES

    calls: dict[str, tuple[str, str]] = dict(di_calls_at_plp)
    for v in variants:
        if v.rsid in high_value and not {v.allele1, v.allele2}.isdisjoint(DI_ALLELES):
            calls.setdefault(v.rsid, (v.genotype, ""))
    sites: list[UncheckedSite] = []
    for rsid, (call, gene) in sorted(calls.items()):
        snp = high_value.get(rsid)
        sites.append(
            UncheckedSite(
                rsid=rsid,
                gene=snp.gene if snp is not None and snp.gene else gene,
                name=snp.name if snp is not None else "",
                reason=UNCHECKED_DI_CALL,
                call=call,
                high_value=snp is not None,
                own_plp=rsid in di_calls_at_plp,
            )
        )
    for rsid in sorted(set(not_in_file) & set(high_value)):
        snp = high_value[rsid]
        sites.append(
            UncheckedSite(
                rsid=rsid,
                gene=snp.gene,
                name=snp.name,
                reason=UNCHECKED_NOT_IN_FILE,
                high_value=True,
            )
        )
    named = {(site.rsid, site.reason) for site in sites}
    for c in sorted(other_calls, key=lambda c: (c.site, c.reason)):
        snp = high_value.get(c.site)
        if (c.site, c.reason) in named or (c.plp_gene is None and snp is None):
            continue
        named.add((c.site, c.reason))
        sites.append(
            UncheckedSite(
                rsid=c.site,
                gene=snp.gene if snp is not None and snp.gene else c.plp_gene or "",
                name=snp.name if snp is not None else "",
                reason=c.reason,
                call=c.call,
                detail=c.detail,
                high_value=snp is not None,
                own_plp=c.own_plp,
            )
        )
    return sites


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

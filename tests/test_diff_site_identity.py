"""A positioned baseline finding pairs only at its published locus.

A baseline with no position retains the legacy identity, independently of
whether it has an ALT. Record numbers and untrimmed writer keys are not the
public position. Synthetic fixtures; no cache or source access.
"""

from dataclasses import replace

import pytest

from allelix.models import Annotation, annotation_to_public_dict
from allelix.reports.diff import compute_diff, diff_annotation_to_dict

RS = "rs900527001"


def ann(locus=("1", 1000), *, good=False, **kwargs):
    values = dict(
        source="clinvar",
        rsid=RS,
        significance="benign" if good else "pathogenic",
        category="clinical",
        magnitude=2.0 if good else 8.0,
        description="Synthetic evidence",
        attribution="ClinVar",
        genotype_match="AG",
        condition="Synthetic condition",
        alt="A",
        locus=locus,
    )
    values.update(kwargs)
    return Annotation(**values)


def public(a):
    return annotation_to_public_dict(a)


def ordered(rows, reverse):
    return list(reversed(rows)) if reverse else rows


@pytest.mark.parametrize("source", ["clinvar", "snpedia", "pharmgkb", "gwas"])
@pytest.mark.parametrize(
    "reverse_old,reverse_new", [(False, False), (True, False), (False, True), (True, True)]
)
@pytest.mark.parametrize("genotype", ["AG", "AA", "A", "GA/A"])
def test_swapped_classifications_are_two_changes_at_their_own_loci(
    source, reverse_old, reverse_new, genotype
):
    first = ann(source=source, genotype_match=genotype, haploid=genotype == "A")
    second = ann(
        ("1", 5000), good=True, source=source, genotype_match=genotype, haploid=genotype == "A"
    )
    old = ordered([public(first), public(second)], reverse_old)
    current = ordered(
        [
            replace(first, significance=second.significance, magnitude=2.0),
            replace(second, significance=first.significance, magnitude=8.0),
        ],
        reverse_new,
    )
    diff = compute_diff(current, old, "2001-01-01")
    assert not diff.new and not diff.removed
    assert len(diff.changed) == 2, "matching identical content across loci hid both changes"
    by_locus = {c.current.locus: c for c in diff.changed}
    assert by_locus[("1", 1000)].previous_significance == "pathogenic"
    assert by_locus[("1", 5000)].previous_significance == "benign"
    assert {diff_annotation_to_dict(c)["position"] for c in diff.changed} == {1000, 5000}


@pytest.mark.parametrize("new_locus", [("1", 5000), ("2", 1000), None])
@pytest.mark.parametrize("same_content", [True, False])
@pytest.mark.parametrize("has_alt", [True, False])
def test_a_positioned_finding_moved_elsewhere_is_removed_and_new(new_locus, same_content, has_alt):
    before = public(ann())
    if not has_alt:
        del before["alt"]
    after = ann(new_locus, good=not same_content)
    diff = compute_diff([after], [before], "")
    assert len(diff.new) == len(diff.removed) == 1, "a moved site must not be paired"
    assert diff.new[0] is after and diff.removed[0] == before
    assert not diff.changed


@pytest.mark.parametrize("position_shape", ["missing", "null"])
@pytest.mark.parametrize("has_alt", [False, True])
@pytest.mark.parametrize("good", [False, True])
def test_legacy_baseline_without_position_keeps_the_old_pairing(position_shape, has_alt, good):
    before = public(ann())
    if position_shape == "missing":
        before.pop("position")
    else:
        before["position"] = None
    before["chromosome"] = None
    if not has_alt:
        before.pop("alt")
    after = ann(("2", 5000), good=good)
    diff = compute_diff([after], [before], "")
    assert not diff.new and not diff.removed
    assert len(diff.changed) == int(good)


@pytest.mark.parametrize(
    "reverse_old,reverse_new", [(False, False), (True, False), (False, True), (True, True)]
)
def test_positioned_match_has_priority_over_positionless_fallback(reverse_old, reverse_new):
    modern = public(ann())
    legacy = public(ann(None, good=True))
    current = [ann(good=True), ann(("1", 5000))]
    diff = compute_diff(ordered(current, reverse_new), ordered([modern, legacy], reverse_old), "")
    assert not diff.new and not diff.removed
    assert len(diff.changed) == 2, "a legacy row must not steal a positioned counterpart"
    by_locus = {c.current.locus: c.previous_significance for c in diff.changed}
    assert by_locus == {("1", 1000): "pathogenic", ("1", 5000): "benign"}


@pytest.mark.parametrize("reverse", [False, True])
def test_row_multiplicity_is_preserved_within_each_locus(reverse):
    before = [ann(), ann(), ann(("1", 5000), good=True)]
    after = [ann(good=True), ann(good=True), ann(("1", 5000))]
    diff = compute_diff(ordered(after, reverse), [public(a) for a in before], "")
    assert not diff.new and not diff.removed
    assert len(diff.changed) == 3
    assert [c.current.locus for c in diff.changed].count(("1", 1000)) == 2
    assert [c.current.locus for c in diff.changed].count(("1", 5000)) == 1


@pytest.mark.parametrize("moved", [False, True])
def test_below_floor_findings_do_not_mask_a_removed_other_site(moved):
    before = public(ann())
    after = ann(("1", 5000) if moved else ("1", 1000), good=True)
    diff = compute_diff([], [before], "", unfiltered=[after], below_floor={id(after)})
    assert not diff.new
    if moved:
        assert diff.removed == [before] and not diff.changed
    else:
        assert not diff.removed and len(diff.changed) == 1
        assert diff.changed[0].below_display_floor


@pytest.mark.parametrize("named", [False, True])
@pytest.mark.parametrize("same_locus", [False, True])
def test_trimmed_public_locus_not_writer_key_or_record_index_decides_pairing(named, same_locus):
    before = ann(
        ("1", 1001),
        rsid=RS if named else "",
        locus_ref="A",
        alt="T",
        genotype_match="AT",
        record_index=0,
        record_key=("1", 1000, "GA", "GT"),
    )
    after = replace(
        before,
        locus=("1", 1001 if same_locus else 5001),
        record_index=99,
        record_key=("1", 5000, "GA", "GT"),
    )
    diff = compute_diff([after], [public(before)], "")
    assert not diff.changed
    assert len(diff.new) == len(diff.removed) == int(not same_locus)


@pytest.mark.parametrize("alt_mode", ["both", "legacy"])
def test_same_site_alleles_and_source_condition_boundaries_stay_distinct(alt_mode):
    rows = [
        ann(alt="A"),
        ann(alt="T", genotype_match="GT"),
        ann(source="snpedia"),
        ann(condition="Other condition"),
    ]
    before = [public(a) for a in rows]
    if alt_mode == "legacy":
        for a in before:
            a.pop("alt")
    diff = compute_diff(list(reversed(rows)), before, "")
    assert not diff.has_changes


def test_same_position_alleles_cannot_swap_across_loci():
    first = ann(alt="A", genotype_match="AT")
    second = ann(("1", 5000), alt="A", genotype_match="AT", good=True)
    control = ann(alt="T", genotype_match="AT")
    diff = compute_diff(
        [
            replace(first, magnitude=2.0, significance="benign"),
            replace(second, magnitude=8.0, significance="pathogenic"),
            control,
        ],
        [public(a) for a in (first, second, control)],
        "",
    )
    assert len(diff.changed) == 2 and not diff.new and not diff.removed
    assert all(c.current.alt == "A" for c in diff.changed)

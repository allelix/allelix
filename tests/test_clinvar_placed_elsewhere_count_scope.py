"""A placement guard counts only an otherwise exact carried VCF allele."""

import pytest
from test_clinvar_own_placement_disclosure import BUILDS, ORDINARY, POS, TARGET, _data

from allelix.annotators.clinvar import ClinVarAnnotator
from allelix.models import NO_CALL_MARKER, Variant


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("batch", (False, True))
@pytest.mark.parametrize("padded", (False, True))
@pytest.mark.parametrize(
    "ref,alleles,exact",
    [
        ("G", ("G", "A"), True),
        ("G", (NO_CALL_MARKER, "A"), True),
        ("C", ("C", "A"), False),
        ("G", ("G", "T"), False),
        ("G", ("G", "G"), False),
        ("G", ("G", NO_CALL_MARKER), False),
    ],
    ids=["exact", "half-alt", "wrong-ref", "alt-not-carried", "hom-ref", "half-ref"],
)
def test_placed_elsewhere_count_requires_an_exact_carried_row(
    tmp_path, build, batch, padded, ref, alleles, exact
):
    data = tmp_path / "data"
    _data(data, build, ORDINARY)

    def call(rsid):
        prefix = "C" if padded else ""
        a1, a2 = (prefix + a if a != NO_CALL_MARKER else a for a in alleles)
        return Variant(
            rsid,
            "1",
            POS - int(padded),
            a1,
            a2,
            ref=prefix + ref,
            ref_from_input=True,
            build=build,
        )

    def run(annotator, variant):
        return list(annotator.batch_annotate([variant])) if batch else annotator.annotate(variant)

    # Establish that the same allele would match (or would not match) if
    # its ID were correctly placed. Wrong REF or an uncarried ALT is not
    # a finding lost specifically to the placed-elsewhere guard.
    with ClinVarAnnotator(data, builds=(build,)) as control:
        assert bool(run(control, call(TARGET))) is exact
    with ClinVarAnnotator(data, builds=(build,)) as annotator:
        assert not run(annotator, call(ORDINARY))
        assert annotator.run_counts()["placed_elsewhere_calls"] == int(exact)
        assert ((ORDINARY, "placed_elsewhere") in annotator.unchecked()) is exact

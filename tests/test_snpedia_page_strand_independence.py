"""Each SNPedia page settles its own strand, independently of sibling pages."""

import sqlite3
from contextlib import closing

import pytest
from validated_call_helpers import BUILDS, POSITIONS, RSID, identity, setup

from allelix.annotators.rsid_identity import RsidIdentityAnnotator
from allelix.annotators.snpedia import SNPediaAnnotator
from allelix.models import Variant

# Site alleles, call, target orientation, sibling genotype/orientation,
# matching page labels, and per-site abstentions. The target page is A/A.
CASES = [
    ("AGT", "AA", "", None, set(), {"both_fit": 1}),
    ("AGT", "AA", "", ("GG", ""), set(), {"both_fit": 1}),
    ("AGT", "TT", "", ("CC", ""), set(), {"both_fit": 1}),
    ("AGT", "AA", "", ("GG", "plus"), set(), {"both_fit": 1}),
    ("AGT", "TT", "", ("CC", "minus"), set(), {"both_fit": 1}),
    ("AG", "AA", "", ("CC", ""), {"target"}, {}),
    ("AG", "AA", "", ("AT", ""), {"target"}, {"no_fit": 1}),
    ("GT", "TT", "", ("GG", ""), {"target"}, {}),
    ("AG", "AA", "", ("TT", ""), {"target", "sibling"}, {}),
    ("AGT", "AA", "", ("TT", ""), set(), {"both_fit": 1}),
    ("CG", "GG", "", ("TT", ""), set(), {"no_fit": 1}),
    ("AGT", "AA", "plus", ("GG", ""), {"target"}, {}),
    ("AGT", "TT", "minus", ("CC", ""), {"target"}, {}),
]
IDS = [
    "ambiguous-alone",
    "ambiguous-beside-direct-page",
    "ambiguous-beside-complement-page",
    "ambiguous-beside-stated-plus",
    "ambiguous-beside-stated-minus",
    "direct-beside-complement-page",
    "direct-beside-invalid-page",
    "complement-beside-direct-page",
    "two-pages-read-opposite-ways",
    "two-ambiguous-pages-count-once",
    "two-invalid-pages-count-once",
    "stated-plus-control",
    "stated-minus-control",
]


def _pages(data, target_orientation, sibling, reverse):
    with closing(sqlite3.connect(data / "snpedia.sqlite")) as conn:
        conn.row_factory = sqlite3.Row
        seed = dict(conn.execute("SELECT * FROM snpedia_genotypes LIMIT 1").fetchone())
        conn.execute("DELETE FROM snpedia_genotypes")
        specs = [("AA", target_orientation, "target")]
        if sibling is not None:
            specs.append((*sibling, "sibling"))
        for genotype, orientation, label in reversed(specs) if reverse else specs:
            row = dict(seed)
            row.update(
                allele1=genotype[0],
                allele2=genotype[1],
                orientation=orientation,
                summary=f"Synthetic {label} page",
            )
            conn.execute(
                f"INSERT INTO snpedia_genotypes ({','.join(row)}) "
                f"VALUES ({','.join('?' for _ in row)})",
                tuple(row.values()),
            )
        conn.commit()


@pytest.mark.parametrize("site,genotype,orientation,sibling,expected,counts", CASES, ids=IDS)
@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("batch", [False, True], ids=["single", "batch"])
@pytest.mark.parametrize("reverse", [False, True], ids=["forward-order", "reverse-order"])
def test_each_page_strand_is_independent(
    tmp_path, site, genotype, orientation, sibling, expected, counts, build, batch, reverse
):
    data = tmp_path / "data"
    setup(data, build)
    identity(data, site)
    _pages(data, orientation, sibling, reverse)
    variant = Variant(
        RSID,
        "1",
        POSITIONS[build],
        genotype[0],
        genotype[1],
        build=build,
        ref="G",
        ref_from_input=True,
    )
    with (
        RsidIdentityAnnotator(data) as ident,
        SNPediaAnnotator(data, forward_site_provider=ident.forward_sites) as source,
    ):
        assert ident.is_ready() and source.is_ready()
        rows = list(source.batch_annotate([variant])) if batch else source.annotate(variant)
        assert {row.description for row in rows} == {
            f"SNPedia: Synthetic {label} page" for label in expected
        }
        assert {reason: count for reason, count in source.run_counts().items() if count} == counts

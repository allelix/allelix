# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""The long P/LP indels for D/I-call naming are read in one scan per build.

RC3 read them with one query per chromosome, each through the position
index and a sort: about 18 s of CPU on a full GRCh37 cache for a chip with
D/I calls on every chromosome, against about 1 s for one scan of the
table. The calls named, and their genes, are the same.
"""

import pytest
from test_outside_review_rc1_b import BUILDS, _clinvar, _config, _identity

from allelix.annotators.clinvar import ClinVarAnnotator
from allelix.models import Variant

POS = 3_100_000
LONG_REF = "G" + "A" * 79  # longer than the fixed window: read with the long indels
CHROMS = ("2", "13", "X")


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("path", ["single", "batch"])
def test_one_long_indel_scan_per_build(tmp_path, build, path):
    data = tmp_path / "data"
    _config(data)
    _identity(data, [])
    plp = ("Pathogenic", "Pathogenic")
    rows = [(f"rs39994910{i}", c, POS, LONG_REF, "G", f"GENE{c}", *plp)
            for i, c in enumerate(CHROMS)]  # fmt: skip
    # A long Benign indel is not read for naming.
    rows.append(("rs3999491099", "7", POS, LONG_REF, "G", "GENEB", "Benign", "Benign"))
    _clinvar(data, build, rows)
    end = POS + len(LONG_REF) - 1
    calls = [Variant(f"i39994911{i}", chrom, end, "D", "D", build=build)
             for i, chrom in enumerate((*CHROMS, "7", "5"))]  # fmt: skip
    scans: list[str] = []
    with ClinVarAnnotator(data, builds=(build,)) as clinvar:
        conn = clinvar._connection(build)
        assert conn is not None
        conn.set_trace_callback(lambda sql: scans.append(sql) if "length(ref)" in sql else None)
        if path == "single":
            for call in calls:
                assert not list(clinvar.annotate(call))
        else:
            assert not list(clinvar.batch_annotate(calls))
        named = {site: gene for (site, reason), (_, gene, _) in clinvar.unchecked().items()
                 if reason == "di_call"}  # fmt: skip
    assert len(scans) == 1, scans
    assert named == {f"i39994911{i}": f"GENE{chrom}" for i, chrom in enumerate(CHROMS)}

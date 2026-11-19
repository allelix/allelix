"""extract shows a split site's passing half call, not a rejected sibling's no-call."""

import pytest
from test_extract_duplicate_placements import RSID, _both, _cache, _input, _row

from allelix.models import NO_CALL_MARKER

REJECTED = {
    "filter": lambda line: line.replace("\tPASS\t", "\tLowQual\t"),
    "sample_ft": lambda line: line.replace("\tGT:GQ\t1/0:99", "\tGT:GQ:FT\t1/0:99:LowQual"),
}


@pytest.mark.parametrize("rejected", sorted(REJECTED))
@pytest.mark.parametrize("rejected_first", [True, False], ids=["rejected_first", "passing_first"])
def test_passing_half_call_beats_rejected_record(tmp_path, monkeypatch, rejected, rejected_first):
    bad = REJECTED[rejected](_row(100, "1/0"))
    good = _row(100, "0/1").replace("\tG\tA\t", "\tG\tC\t")
    source = _input(tmp_path, [bad, good] if rejected_first else [good, bad])
    # Indexed and sequential, each the same.
    for found in _both(tmp_path, monkeypatch, _cache(tmp_path, 100), source):
        chrom, pos, a1, a2, *_ = found[RSID]
        assert (chrom, pos) == ("1", 100)
        assert sorted((a1, a2)) == sorted(("C", NO_CALL_MARKER)), (
            "The G>C record passes and analyze assesses C/.; the rejected G>A record "
            "must not hide the known C behind a full no-call, whatever the record order"
        )


def test_full_no_call_at_placement_still_beats_usable_call_elsewhere(tmp_path, monkeypatch):
    source = _input(tmp_path, [_row(100, "./."), _row(300, "0/1")])
    for found in _both(tmp_path, monkeypatch, _cache(tmp_path, 100), source):
        assert found[RSID][:4] == ("1", 100, NO_CALL_MARKER, NO_CALL_MARKER)

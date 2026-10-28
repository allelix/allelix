# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""A readable tabix header does not make its chunk offsets usable."""

import gzip
import struct

import pysam
import pytest


def test_bad_chunk_offset_warns_and_falls_back_to_sequential(rc_case):
    case = rc_case
    case.source.write_text(case.source.read_text().replace("249250621", "248956422"))
    source = case.root / "bad-chunk.vcf.gz"
    pysam.tabix_compress(str(case.source), str(source))
    pysam.tabix_index(str(source), preset="vcf")
    index = source.with_suffix(".gz.tbi")
    raw = bytearray(gzip.decompress(index.read_bytes()))
    assert raw[:4] == b"TBI\x01"
    # The header's name length is followed by the first reference's bins.
    offset = 36 + struct.unpack_from("<i", raw, 32)[0]
    bins = struct.unpack_from("<i", raw, offset)[0]
    assert bins > 0
    offset += 4
    bin_id, chunks = struct.unpack_from("<Ii", raw, offset)
    assert bin_id < 37450 and chunks > 0
    # A valid index container points into the middle of a compressed block.
    invalid_offset = 1 << 16
    struct.pack_into("<QQ", raw, offset + 8, invalid_offset, invalid_offset + 100)
    with pysam.BGZFile(str(index), "wb") as stream:
        stream.write(bytes(raw))
    with pysam.TabixFile(str(source)) as reader, pytest.raises((OSError, ValueError)):
        list(reader.fetch("1", 999, 1000))
    result = case.invoke(
        ["extract", source, "--data-dir", case.data, "--build", "grch38", "--snps", case.rsid]
    )
    assert result.exit_code == 0, (result.output, repr(result.exception))
    assert case.rsid in result.output and "A/G" in result.output, result.output
    assert "warn" in result.output.lower() and index.name in result.output, result.output
    assert "Traceback" not in result.output

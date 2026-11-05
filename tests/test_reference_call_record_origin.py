# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""An explicit reference-carried finding keeps its input record identity."""

import pytest
from validated_call_helpers import BUILDS, analyze, call, reference_effect


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("gt", ["0", "0/0", "0/1", "0/."])
def test_reference_carried_finding_keeps_its_record_index(tmp_path, build, gt):
    data = tmp_path / "data"
    pos = reference_effect(data, build)
    result = analyze(data, call(tmp_path / "input.vcf", pos, gt=gt), build)
    found = [a for a in result.annotations if a.source == "clinvar"]
    assert found
    assert all(a.record_index == 0 for a in found)

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Pin the real ClinPGx zip layout against the loader's hardcoded filenames.

Captured manually from: curl -A "allelix" \
  https://api.pharmgkb.org/v1/download/file/data/clinicalAnnotations.zip | unzip -l
"""

from __future__ import annotations

from allelix.databases.pharmgkb_loader import (
    CLINICAL_ANN_ALLELES_FILENAME,
    CLINICAL_ANN_FILENAME,
)

_REAL_ZIP_MEMBERS = {
    "LICENSE.txt",
    "CREATED_2025-07-05.txt",
    "clinical_annotations.tsv",
    "clinical_ann_alleles.tsv",
    "clinical_ann_history.tsv",
    "clinical_ann_evidence.tsv",
    "README.pdf",
}


def test_loader_filenames_match_real_pharmgkb_zip():
    assert CLINICAL_ANN_FILENAME in _REAL_ZIP_MEMBERS
    assert CLINICAL_ANN_ALLELES_FILENAME in _REAL_ZIP_MEMBERS

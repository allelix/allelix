# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Methylation pathway gene panel.

The set is intentionally small and curated — covering the one-carbon /
folate / methylation cycle genes most often discussed in consumer methylation
reports. Not medical guidance; see ADR-0003 (regulatory posture).
"""

from __future__ import annotations

# Folate / one-carbon / methylation cycle genes. Add via PR + ADR if expanding.
METHYLATION_PANEL_GENES: frozenset[str] = frozenset(
    {
        "ACAT1",
        "AHCY",
        "BHMT",
        "BHMT2",
        "CBS",
        "COMT",
        "DHFR",
        "DNMT1",
        "DNMT3A",
        "DNMT3B",
        "FOLR1",
        "FOLR2",
        "FUT2",
        "GNMT",
        "GSTM1",
        "GSTP1",
        "MAOA",
        "MAT1A",
        "MAT2A",
        "MAT2B",
        "MTHFD1",
        "MTHFD1L",
        "MTHFR",
        "MTR",
        "MTRR",
        "NOS3",
        "PEMT",
        "SHMT1",
        "SHMT2",
        "SLC19A1",
        "SUOX",
        "TCN1",
        "TCN2",
        "VDR",
    }
)

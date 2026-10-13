# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Interpreter version stamps for annotator cache invalidation.

Increment the constant for an annotator when its emit/suppression logic
changes in a way that should invalidate prior reports built against
existing caches.  The stamp is stored in the ``local_version_tag``
column of ``database_versions`` (e.g. ``iv:1``) so ``is_ready()`` can
reject stale caches without forcing a full re-download.
"""

# v2.2: per-SCV TSV loader (variant_summary + submission_summary)
# v2.2 follow-up: iv:3→iv:4 invalidates
# caches built by the broken loader that left placeholder CLNSIG values
# ("-", "not specified", "no classification provided", etc.) in the cache.
# Without this bump, existing iv:3 caches stay poisoned — is_ready() would
# treat them as current and the fix would only reach fresh / --force builds.
# iv:4→iv:5: the loader classifies records by their HGVS name (2.4.2).
# Normal-allele records (c.N=) lose their REF == ALT rows and their
# Pathogenic / Likely pathogenic rows, and REF == ALT non-SNV rows are
# dropped, so a REF == ALT row means this build's reference carries the
# variant allele. An iv:4 cache still holds the old rows and must not be
# read with the new rule.
CLINVAR_INTERPRETER_VERSION = 5
# Also bump for a change to how CPIC allele-function data is loaded: a bump
# makes existing ClinPGx caches not ready, and `allelix db update` rebuilds
# them, fetching the CPIC table again, so the change reaches them.
PHARMGKB_INTERPRETER_VERSION = 1
GNOMAD_SCHEMA_VERSION = 1
ALPHAMISSENSE_SCHEMA_VERSION = 1
CADD_SCHEMA_VERSION = 1
RSID_IDENTITY_SCHEMA_VERSION = 1

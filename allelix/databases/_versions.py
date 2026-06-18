# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Interpreter version stamps for annotator cache invalidation.

Increment the constant for an annotator when its emit/suppression logic
changes in a way that should invalidate prior reports built against
existing caches.  The stamp is stored in the ``local_version_tag``
column of ``database_versions`` (e.g. ``iv:1``) so ``is_ready()`` can
reject stale caches without forcing a full re-download.
"""

# v2.2 #42 stage B: per-SCV TSV loader (variant_summary + submission_summary)
# v2.2 #42 follow-up (evaluator defect 5, PR #101): iv:3→iv:4 invalidates
# caches built by the broken loader that left placeholder CLNSIG values
# ("-", "not specified", "no classification provided", etc.) in the cache.
# Without this bump, existing iv:3 caches stay poisoned — is_ready() would
# treat them as current and the fix would only reach fresh / --force builds.
CLINVAR_INTERPRETER_VERSION = 4
PHARMGKB_INTERPRETER_VERSION = 1
GNOMAD_SCHEMA_VERSION = 1
ALPHAMISSENSE_SCHEMA_VERSION = 1
CADD_SCHEMA_VERSION = 1

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Interpreter version stamps for annotator cache invalidation.

Increment the constant for an annotator when its emit/suppression logic
changes in a way that should invalidate prior reports built against
existing caches.  The stamp is stored in the ``local_version_tag``
column of ``database_versions`` (e.g. ``iv:1``) so ``is_ready()`` can
reject stale caches without forcing a full re-download.
"""

CLINVAR_INTERPRETER_VERSION = 1
PHARMGKB_INTERPRETER_VERSION = 1
GNOMAD_SCHEMA_VERSION = 1
ALPHAMISSENSE_SCHEMA_VERSION = 1
CADD_SCHEMA_VERSION = 1

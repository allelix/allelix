# ADR-0031: Centralized License Descriptors on Annotator Base Class

**Status:** Accepted
**Date:** 2026-06-10

## Context

License metadata for each data source was duplicated across four
hand-maintained locations with no single source of truth:

1. `allelix/config.py` — `NON_COMMERCIAL_SOURCES` frozenset
2. `allelix/reports/html.py` — `_LICENSE_ATTRIBUTIONS` dict
3. `allelix/reports/json_report.py` — `_LICENSE_ATTRIBUTIONS` dict
4. Annotator module docstrings (informal, not machine-readable)

Adding or changing a source required editing all four, and nothing
enforced agreement. The GWAS Catalog was a live instance of this
failure mode: it was a registered annotator but absent from both
attribution maps. The problem grows linearly with each new source
(CADD in v1.6.0 would have been a five-location edit).

## Decision

Each annotator declares a `license` ClassVar of type
`LicenseDescriptor` — a frozen dataclass with these fields:

- `spdx`: SPDX license identifier (e.g. `"CC-BY-SA-4.0"`,
  `"CC-BY-NC-SA-3.0-US"`, `"custom-clinvar"`, `"LicenseRef-CADD"`)
- `license_url`: license deed or terms-of-use URL
- `attribution_text`: one-line human-readable attribution string
- `source_url`: optional source website URL (e.g. pharmgkb.org)
- `citation`: optional citation string (e.g. AlphaMissense paper)
- `commercial_ok`: optional bool — required when `spdx` starts with
  `LicenseRef-` or `custom-` (enforced by `__init_subclass__`).
  Overrides the SPDX allowlist when set.

The `Annotator` base class declares `license: ClassVar[LicenseDescriptor]`
with no default, so any subclass that omits it will raise
`AttributeError` on access.

Non-commercial gating is derived from `commercial_ok` first (if set),
then falls back to the SPDX identifier via
`is_non_commercial(spdx)`, which checks against a frozen set of
known NC identifiers. The hand-maintained `NON_COMMERCIAL_SOURCES`
frozenset in `config.py` is deleted.

Report renderers (`html.py`, `json_report.py`) generate attribution
text from the `LicenseDescriptor` fields at render time, replacing
the hardcoded `_LICENSE_ATTRIBUTIONS` dicts.

## Consequences

- Adding a new data source requires declaring its license on the
  annotator class. Omitting it is a construction-time error, not a
  silent gap.
- Non-commercial gating is a property of the license, not a
  separate manually maintained set.
- Report attributions are always in sync with the annotator's
  declared license.
- JSON schema bumped from 2 to 3. The `license_attributions` block
  now carries `source_url` (source website) and `license_url`
  (license deed) as separate keys. The `license` field uses SPDX
  identifiers. The schema version is a monotonic counter for breaking
  output changes; it does not correspond to the "JSON v3" roadmap
  feature set (structured caveats, conflict flags, provenance).

# ADR-0003: Source-attributed annotations (regulatory posture)

- **Date:** 2026-05-11
- **Status:** Accepted

## Context

Tools that classify genetic variants and produce "health reports" land in regulator-adjacent territory. In the US, the FDA has historically treated software that interprets clinical genetic data as a medical device. Promethease, an established commercial analog, is positioned explicitly as informational rather than diagnostic.

If Allelix labels variants as "pathogenic" in its own voice, every report becomes a medical claim Allelix made. That is the wrong posture for an open-source tool that runs on the user's laptop with no clinical oversight.

## Decision

Allelix never asserts variant significance directly. Every annotation is attributed to its source database. The tool says "ClinVar classifies this variant as pathogenic," never "this variant is pathogenic."

Concrete implications across the codebase:

- `Annotation.attribution` (e.g., `"ClinVar"`) is a required field — every annotation carries the source it came from in user-facing language.
- `Annotation.significance` is source-prefixed: `"clinvar_pathogenic"`, `"pharmgkb_interaction"`, never bare `"pathogenic"`.
- `Annotation.category` uses non-diagnostic filter buckets (`"clinical"`, `"pharma"`, `"carrier"`, `"trait"`, `"methylation"`) — never bare medical terms.
- Report rendering (HTML, terminal, JSON) must surface attribution prominently. A user reading a report sees what ClinVar/PharmGKB/etc. said, not what Allelix decided.
- CLI flags and category filters use the same vocabulary. `--category clinical` filters a bucket; it is not a clinical claim.

## Consequences

- This is not a disclaimer afterthought. It shapes data model field names, report templates, and CLI vocabulary throughout the codebase. PR review must catch language drift.
- Code that reads "pathogenic" as a bare label in a model docstring, test fixture, or report template fails this ADR.
- Allelix accepts that some users will mentally collapse "ClinVar says pathogenic" into "is pathogenic." That is unavoidable; the attribution is what protects the tool.
- A separate `--exclude-snpedia` flag exists for commercial users (license reason — see ADR-0004) but has no bearing on regulatory posture.

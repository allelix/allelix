# Security Policy

## Supported versions

The current minor release receives security fixes.

| Version | Supported |
|---------|-----------|
| 1.8.x   | ✓         |
| < 1.8   | ✗         |

## Reporting a vulnerability

Report security issues through
[GitHub's private vulnerability reporting](https://github.com/allelix/allelix/security/advisories/new).
Do not open a public issue for security vulnerabilities.

Best-effort response, typically within two weeks.

## Scope

**In scope:** allelix source code, download integrity (ADR-0029),
local file handling, CLI behavior.

**Out of scope:** upstream database content correctness. ClinVar
misclassifications, PharmGKB annotation errors, GWAS Catalog data
quality, and similar issues are third-party data — report them to
the source database. Allelix verifies download integrity but does
not and cannot verify the clinical accuracy of what those databases
contain.

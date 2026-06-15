#!/usr/bin/env bash
# Regenerate the sample reports from sample_input/demo_23andme.txt.
#
# Requires:
#   - allelix installed (pip install -e ".[dev]" or from PyPI)
#   - reference databases downloaded (allelix db update)
#   - SNPedia must NOT be downloaded if you intend to share these reports
#     (SNPedia is CC BY-NC-SA 3.0 — non-commercial only). The --exclude-snpedia
#     flag below guarantees no SNPedia content reaches the output regardless
#     of what's in the cache.
#
# Database versions will differ from when the committed reports were
# generated; the annotations may shift slightly between runs as upstream
# databases publish new releases. That's expected.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

INPUT="$SCRIPT_DIR/sample_input/demo_23andme.txt"
OUT="$SCRIPT_DIR/sample_reports"

if [ ! -f "$INPUT" ]; then
    echo "Error: $INPUT not found." >&2
    exit 1
fi

if ! command -v allelix >/dev/null 2>&1; then
    echo "Error: 'allelix' not on PATH. Install with: pip install allelix" >&2
    exit 1
fi

mkdir -p "$OUT"

echo "Generating HTML report..."
allelix analyze "$INPUT" --exclude-snpedia --output "$OUT/demo_report.html"

echo "Generating JSON report..."
allelix analyze "$INPUT" --exclude-snpedia --output "$OUT/demo_report.json"

echo "Generating stats..."
allelix stats "$INPUT" > "$OUT/demo_stats.txt"

echo ""
echo "Done. Open $OUT/demo_report.html in a browser to view."

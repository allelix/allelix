#!/usr/bin/env bash
set -euo pipefail

RELEASE_URL="https://github.com/dial481/allelix/releases/download/v1.1.1/test_data.tar.gz"
GWAS_URL="https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip"
DEST="test_data"

# Real genotype fixtures (from GitHub release asset)
if [ -d "$DEST/real" ] && [ -d "$DEST/transcoded" ]; then
    echo "Test data already present."
else
    echo "Downloading test data..."
    curl -L -o test_data.tar.gz "$RELEASE_URL"
    tar -xzf test_data.tar.gz
    rm test_data.tar.gz
    echo "Done. $(find test_data/real -type f | wc -l) files extracted."
fi

# GWAS Catalog (for @slow integration tests)
if [ -f "$DEST/gwas_catalog.zip" ]; then
    echo "GWAS Catalog already present."
else
    echo "Downloading GWAS Catalog (~66 MB)..."
    curl -L -o "$DEST/gwas_catalog.zip" "$GWAS_URL"
    echo "Done."
fi

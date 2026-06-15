#!/bin/bash
# scripts/tag-release.sh
set -e

REMOTE="${1:-private}"
VERSION=$(python3 -c "import tomllib, pathlib; print(tomllib.loads(pathlib.Path('pyproject.toml').read_text())['project']['version'])")
TAG="v${VERSION}"

if git rev-parse "$TAG" >/dev/null 2>&1; then
    echo "ABORT: tag $TAG already exists"
    exit 1
fi

echo "Tagging $TAG and pushing main + tag to $REMOTE"
git tag -a "$TAG" -m "$TAG"
git push "$REMOTE" main "$TAG"

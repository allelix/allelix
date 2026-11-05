# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""ADR-0029: prebuilt HF caches use immutable revisions, with no network or git.

Upstream source data and live API releases are not prebuilt HF caches. Their
URLs are inventoried separately; this contract does not freeze those releases.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path
from urllib.parse import urlsplit

import pytest

import allelix
from allelix.databases import (
    alphamissense_loader,
    cadd_loader,
    gnomad_loader,
    rsid_identity_loader,
    snpedia_loader,
)

HF_CACHE_URLS = {
    "alphamissense": alphamissense_loader.ALPHAMISSENSE_CACHE_URL,
    "cadd": cadd_loader.CADD_CACHE_URL,
    "gnomad": gnomad_loader.GNOMAD_CACHE_URL,
    "rsid_identity": rsid_identity_loader.RSID_IDENTITY_CACHE_URL,
    "snpedia": snpedia_loader.SNPEDIA_CACHE_URL,
}


@pytest.mark.parametrize("source", sorted(HF_CACHE_URLS))
def test_prebuilt_hf_cache_url_has_a_full_commit_pin(source: str) -> None:
    """Reject moving branches, tags and abbreviated commit names."""
    url = HF_CACHE_URLS[source]
    assert re.fullmatch(
        r"https://huggingface\.co/datasets/[^/?#]+/[^/?#]+/resolve/"
        r"[0-9a-fA-F]{40}/[^/?#]+",
        url,
    ), f"{source}: expected resolve/<40-hex commit>/<asset>, got {url}"


def test_all_hf_url_literals_are_in_the_prebuilt_cache_inventory() -> None:
    """A new HF endpoint cannot silently escape the per-source contract.

    Python's AST joins adjacent literals, including the loaders' split URLs.
    This scans package code, not tests or documentation; it needs no git tags
    and works in a shallow CI checkout. Runtime constants are tested above.
    """
    package = Path(allelix.__file__).resolve().parent
    found: set[str] = set()
    unexpected: list[str] = []
    for path in sorted(package.rglob("*.py")):
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if not isinstance(node, ast.Constant) or not isinstance(node.value, str):
                continue
            url = node.value
            if not url.startswith(("https://", "http://")):
                continue
            if urlsplit(url).hostname != "huggingface.co":
                continue
            found.add(url)
            if url not in HF_CACHE_URLS.values():
                unexpected.append(f"{path.relative_to(package)}:{node.lineno}: {url}")
    assert not unexpected, "HF URL missing from the pin inventory: " + "; ".join(unexpected)
    assert found == set(HF_CACHE_URLS.values()), "Recheck the download-URL class inventory"


def test_identity_url_pin_keeps_the_released_asset_digest() -> None:
    """R1 pins the existing identity bytes; it does not publish a new cache."""
    assert rsid_identity_loader.RSID_IDENTITY_EXPECTED_SHA256 == (
        "2a00787cd0be4e7572efd21686013c49f3d0180137cd4403730524d50a5db285"
    )
    assert urlsplit(rsid_identity_loader.RSID_IDENTITY_CACHE_URL).path.endswith(
        "/rsid_identity.sqlite.gz"
    )

# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for the Annotator abstract base class default behaviors."""

from __future__ import annotations

from typing import TYPE_CHECKING

from allelix.annotators.base import Annotator, LicenseDescriptor
from allelix.models import Annotation, Variant

if TYPE_CHECKING:
    from pathlib import Path


class _StubAnnotator(Annotator):
    """Minimal Annotator subclass for testing default base-class methods.

    Implements the abstract surface and returns one Annotation per
    variant from ``annotate()``. Does NOT override ``batch_annotate``,
    so it exercises the default loop fallback.
    """

    name = "stub"
    display_name = "Stub"
    attribution = "Stub"
    requires_download = False
    license = LicenseDescriptor(
        spdx="CC0-1.0",
        license_url="https://example.test/license",
        attribution_text="Stub",
        commercial_ok=True,
    )

    def setup(self) -> None:
        pass

    def annotate(self, variant: Variant) -> list[Annotation]:
        return [
            Annotation(
                source=self.name,
                rsid=variant.rsid,
                significance="stub_sig",
                category="clinical",
                magnitude=1.0,
                description=f"stub for {variant.rsid}",
                attribution=self.attribution,
                genotype_match=variant.genotype,
                gene=None,
            )
        ]

    def is_ready(self) -> bool:
        return True

    def version(self) -> str | None:
        return "stub-1"

    def close(self) -> None:
        pass

    def fetch_remote_signal(self) -> str | None:
        return None

    def cached_remote_signal(self) -> str | None:
        return None


def _variant(rsid: str) -> Variant:
    return Variant(
        rsid=rsid,
        chromosome="1",
        position=1000,
        allele1="A",
        allele2="G",
        build="GRCh37",
    )


class TestBatchAnnotateDefault:
    """Default `batch_annotate` falls back to looping over `annotate()`.

    Existing annotators don't override `batch_annotate`; the pipeline still
    needs to work with them. This pins the contract: in arrival order, every
    annotation from the per-variant path also appears via the batch path.
    """

    def test_default_yields_annotations_for_each_variant(self, tmp_path: Path):
        ann = _StubAnnotator(tmp_path)
        variants = [_variant(f"rs{i}") for i in (1, 2, 3, 4)]
        results = list(ann.batch_annotate(variants))
        assert len(results) == 4
        assert [a.rsid for a in results] == ["rs1", "rs2", "rs3", "rs4"]

    def test_default_preserves_order(self, tmp_path: Path):
        """Order matters for downstream rollup and rendering."""
        ann = _StubAnnotator(tmp_path)
        variants = [_variant(f"rs{i}") for i in (42, 7, 99, 1)]
        results = list(ann.batch_annotate(variants))
        assert [a.rsid for a in results] == ["rs42", "rs7", "rs99", "rs1"]

    def test_default_empty_input_yields_empty(self, tmp_path: Path):
        ann = _StubAnnotator(tmp_path)
        assert list(ann.batch_annotate([])) == []

    def test_default_matches_per_variant_path(self, tmp_path: Path):
        """`batch_annotate(vs)` == flattened `annotate(v)` for each v.

        This is the contract subclasses must preserve when they override
        with a batched SQL query.
        """
        ann = _StubAnnotator(tmp_path)
        variants = [_variant(f"rs{i}") for i in (1, 2, 3)]
        per_variant = [a for v in variants for a in ann.annotate(v)]
        batched = list(ann.batch_annotate(variants))
        assert per_variant == batched

    def test_default_handles_iterator_input(self, tmp_path: Path):
        """`batch_annotate` accepts any Iterable, not just lists."""
        ann = _StubAnnotator(tmp_path)
        variants_iter = iter([_variant("rs1"), _variant("rs2")])
        results = list(ann.batch_annotate(variants_iter))
        assert [a.rsid for a in results] == ["rs1", "rs2"]

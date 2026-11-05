"""Annotation.record_index is a declared, internal field (ADR-0062)."""

from dataclasses import fields

from allelix.models import Annotation, annotation_to_public_dict


def test_record_index_is_declared_internal_and_not_public():
    declared = {f.name: f for f in fields(Annotation)}
    assert "record_index" in declared
    assert declared["record_index"].default is None
    annotation = Annotation(
        source="gwas",
        rsid="rs1",
        significance="",
        category="trait",
        magnitude=0.0,
        description="",
        attribution="GWAS Catalog",
        genotype_match="AG",
    )
    annotation.record_index = 3
    assert "record_index" not in annotation_to_public_dict(annotation)

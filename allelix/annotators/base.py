# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Abstract base class for reference-database annotators."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum, auto
from typing import TYPE_CHECKING, ClassVar

from allelix.databases.cache_state import READY, CacheState, CacheStatus
from allelix.parsers._helpers import normalize_chromosome
from allelix.utils.allele import derive_alt_from_diploid

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Iterator
    from pathlib import Path
    from types import TracebackType

    from allelix.models import Annotation, Variant

    # ``(rsids, build) -> {rsid: (chrom, pos, forward alleles)}``; wired to
    # ``RsidIdentityAnnotator.forward_sites``.
    ForwardSiteProvider = Callable[[set[str], str], dict[str, tuple[str, int, frozenset[str]]]]


@dataclass(frozen=True)
class LicenseDescriptor:
    """Single source of truth for a data source's license terms."""

    spdx: str
    license_url: str
    attribution_text: str
    source_url: str | None = None
    citation: str | None = None
    commercial_ok: bool | None = None
    licensable: bool = False
    purchase_url: str | None = None


_NON_COMMERCIAL_SPDX: frozenset[str] = frozenset(
    {
        "CC-BY-NC-SA-3.0-US",
        "CC-BY-NC-SA-4.0",
        "CC-BY-NC-4.0",
    }
)


def is_non_commercial(descriptor: LicenseDescriptor) -> bool:
    """Return True if the license prohibits commercial use."""
    if descriptor.commercial_ok is not None:
        return not descriptor.commercial_ok
    return descriptor.spdx in _NON_COMMERCIAL_SPDX


class Permission(Enum):
    """Three-state permission result for a source in the current license context."""

    ALLOW = auto()
    BLOCK_FINAL = auto()
    BLOCK_PURCHASABLE = auto()


def permission(
    descriptor: LicenseDescriptor,
    *,
    commercial: bool,
    license_held: bool,
) -> Permission:
    """Determine whether a source is permitted under the current license context."""
    if not commercial:
        return Permission.ALLOW
    if not is_non_commercial(descriptor):
        return Permission.ALLOW
    if not descriptor.licensable:
        return Permission.BLOCK_FINAL
    if license_held:
        return Permission.ALLOW
    return Permission.BLOCK_PURCHASABLE


def variant_locus(variant: Variant) -> tuple[str, int]:
    """The variant's ``(chromosome, position)`` with the chromosome in bare form.

    Normalized here rather than trusted from each parser, so an input that
    writes ``chrX`` or ``M`` still meets the caches' ``X`` / ``MT`` keys.
    """
    return normalize_chromosome(variant.chromosome), variant.position


def forward_alleles_at_call(
    variant: Variant, sites: dict[str, tuple[str, int, frozenset[str]]]
) -> frozenset[str] | None:
    """The forward alleles of the call's rsID, or None when it is not placed there.

    None when the identity table does not hold the rsID, or places it at
    another locus than the call's.
    """
    site = sites.get(variant.rsid)
    if site is None:
        return None
    chrom, pos, alleles = site
    if variant.position and variant_locus(variant) != (chrom, pos):
        return None
    return alleles


def is_clinvar_homref(
    variant: Variant,
    clinvar_ref_provider: Callable[[str, str, tuple[str, int]], str | None] | None,
) -> bool:
    """Return True if the variant is homozygous reference per ClinVar (ADR-0023).

    ``clinvar_ref_provider`` is ``(rsid, build, locus) -> REF | None``
    (``ClinVarAnnotator.reference_at``): ClinVar's REF counts only where
    ClinVar places the rsID at the call's own locus. A call with no
    position is never suppressed.
    """
    if clinvar_ref_provider is None or not variant.position:
        return False
    ref = clinvar_ref_provider(variant.rsid, variant.build, variant_locus(variant))
    return ref is not None and len(ref) == 1 and variant.allele1 == ref and variant.allele2 == ref


def carried_alt(
    variant: Variant,
    allele1: str,
    allele2: str,
    normal_allele_provider: Callable[[str, str], str | None] | None,
) -> tuple[str, bool]:
    """The non-normal allele of a matched genotype, for per-genotype sources.

    Returns ``(alt, reference_carries_alt)``. The alt is derived against the
    normal allele, which is ``variant.ref`` except where the build's
    reference carries a ClinVar variant allele: there
    ``normal_allele_provider`` (``ClinVarAnnotator.normal_allele_for``) gives
    the normal allele from the other build (C at F5 rs6025 on GRCh37, whose
    reference is the Leiden T). The alt is then the same on both builds,
    and it can equal ``variant.ref``, which the Annotation is told through
    ``reference_carries_alt``.
    """
    normal = variant.ref
    if normal is not None and normal_allele_provider is not None and variant.rsid.startswith("rs"):
        other = normal_allele_provider(variant.rsid, variant.build)
        if other is not None and len(other) == 1 and other != normal:
            normal = other
    alt = derive_alt_from_diploid(normal, allele1, allele2)
    return alt, bool(alt) and alt == variant.ref


class Annotator(ABC):
    """Base class for all reference database annotators.

    Annotators bind to a `data_dir` at construction. `setup()` performs the
    one-time download/parse into the cache. `is_ready()` reports whether the
    cache exists and is queryable. `annotate(variant)` returns zero or more
    `Annotation` objects for the variant — checking both rsid AND genotype, per
    the regulatory posture (ADR-0003) and the genotype-matching rule (ADR-0007).

    Annotators hold resources (SQLite connections, file handles). Always close
    them via `close()` or the context manager protocol; the CLI uses
    `contextlib.ExitStack` to guarantee deterministic cleanup.

    Attributes:
        name: Lowercase identifier (e.g., "clinvar").
        display_name: Human-readable name ("ClinVar").
        attribution: Display label used in user-facing reports ("ClinVar").
            Equal to `display_name` for first-party single-source annotators.
        requires_download: Whether `setup()` needs network/disk space.
    """

    name: ClassVar[str]
    display_name: ClassVar[str]
    attribution: ClassVar[str]
    requires_download: ClassVar[bool] = True
    server_driven_freshness: ClassVar[bool] = True
    license: ClassVar[LicenseDescriptor]

    def __init_subclass__(cls, **kwargs: object) -> None:
        """Enforce required ClassVars at subclass definition time."""
        super().__init_subclass__(**kwargs)
        is_abstract = any(getattr(v, "__isabstractmethod__", False) for v in cls.__dict__.values())
        if not is_abstract and not hasattr(cls, "license"):
            msg = f"{cls.__name__} must declare a 'license' ClassVar of type LicenseDescriptor"
            raise TypeError(msg)
        if hasattr(cls, "license"):
            desc = cls.license
            if (
                desc.spdx.startswith("LicenseRef-") or desc.spdx.startswith("custom-")
            ) and desc.commercial_ok is None:
                msg = (
                    f"{cls.__name__} uses custom SPDX '{desc.spdx}' but "
                    f"does not declare commercial_ok (True or False)"
                )
                raise TypeError(msg)
            if desc.licensable and desc.purchase_url is None:
                msg = (
                    f"{cls.__name__} declares licensable=True but "
                    f"purchase_url is None — set it explicitly"
                )
                raise TypeError(msg)

    def __init__(self, data_dir: Path) -> None:
        """Bind the annotator to a data directory (created elsewhere)."""
        self.data_dir = data_dir

    def __enter__(self) -> Annotator:
        """Return self for `with` usage."""
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: TracebackType | None,
    ) -> None:
        """Release any open resources via `close()`."""
        self.close()

    @abstractmethod
    def setup(self) -> None:
        """Download and prepare the reference database. Idempotent."""
        ...

    @abstractmethod
    def annotate(self, variant: Variant) -> list[Annotation]:
        """Return all annotations for this variant.

        Implementations MUST verify both rsid AND genotype — presence in the
        database is not enough. The user must carry the flagged allele.
        """
        ...

    def batch_annotate(self, variants: Iterable[Variant]) -> Iterator[Annotation]:
        """Annotate a batch of variants. Yields annotations in arrival order.

        Default implementation loops over ``annotate(variant)`` so any
        existing annotator works unchanged. Subclasses with rsID-based
        SQLite lookups should override this with a chunked ``WHERE rsid
        IN (...)`` query to avoid per-variant round-trips at WGS scale
        (4-6M variants per VCF).

        The default keeps the pipeline single-path: callers always use
        ``batch_annotate``; the loop fallback gives backward compatibility
        for annotators that haven't grown a batched query path yet.
        """
        for variant in variants:
            yield from self.annotate(variant)

    @abstractmethod
    def is_ready(self) -> bool:
        """Whether the local cache exists and is queryable.

        Annotators with a cache return ``self.cache_status().ready``, so
        readiness has one truth (ADR-0048). Never writes.
        """
        ...

    def cache_status(self) -> CacheStatus:
        """The cache's state, and why it cannot be used (ADR-0048).

        Annotators with a cache override this; the default serves one
        without (ready or missing).
        """
        return READY if self.is_ready() else CacheStatus(CacheState.MISSING)

    def prepare_cache(self) -> bool:
        """Bring an outdated cache up to date from local data only (ADR-0048).

        Runs in the update step (``db update``, or ``analyze`` without
        ``--no-update``), never in ``is_ready()``. It may write; it never
        uses the network. Returns True when it changed the cache.
        """
        return False

    def serves_build(self, build: str) -> bool:
        """Whether this source can annotate input on ``build`` in this run."""
        return True

    def run_counts(self) -> dict[str, int]:
        """What this source counted during the run, by name (empty by default)."""
        return {}

    @abstractmethod
    def version(self) -> str | None:
        """Return the cached database version, or None if not set up."""
        ...

    @abstractmethod
    def close(self) -> None:
        """Release any open resources (database connections, file handles)."""
        ...

    @abstractmethod
    def fetch_remote_signal(self) -> str | None:
        """Fetch a small remote freshness signal (md5 hash, ETag, Last-Modified).

        Implementations MUST return a prefixed, opaque string (e.g.
        `"md5:abcdef…"`, `"lm:Wed, 21 Oct 2025 …"`, `"etag:…"`) so that a
        future server-side switch in signal type triggers a refresh
        rather than a silent miss.

        Returns None on any failure (network error, timeout, missing
        header, source doesn't expose a signal). Never raises — `db update`
        treats None as "can't verify freshness" and falls through to skip
        with a notice. See ADR-0012.
        """
        ...

    @abstractmethod
    def cached_remote_signal(self) -> str | None:
        """Return the remote signal stored at last successful download, or None.

        Returns None if the cache is missing entirely OR if the cache was
        written by a pre-v0.4.2 release that didn't capture signals.
        """
        ...

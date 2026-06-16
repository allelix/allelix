# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Abstract base class for reference-database annotators."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum, auto
from typing import TYPE_CHECKING, ClassVar

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Iterator
    from pathlib import Path
    from types import TracebackType

    from allelix.models import Annotation, Variant


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


def is_clinvar_homref(
    variant: Variant,
    clinvar_ref_provider: Callable[[str, str], str | None] | None,
) -> bool:
    """Return True if the variant is homozygous reference per ClinVar (ADR-0023)."""
    if clinvar_ref_provider is None:
        return False
    ref = clinvar_ref_provider(variant.rsid, variant.build)
    return ref is not None and len(ref) == 1 and variant.allele1 == ref and variant.allele2 == ref


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
        """Whether the local cache exists and is queryable."""
        ...

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

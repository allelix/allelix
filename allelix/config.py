# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""User configuration for source management and license mode."""

from __future__ import annotations

import logging
import math
import tomllib
from dataclasses import dataclass, field, fields, replace
from typing import TYPE_CHECKING

from allelix.annotators.base import Permission
from allelix.annotators.base import permission as check_permission
from allelix.utils.plausibility import PlausibilityThresholds

if TYPE_CHECKING:
    from pathlib import Path

    from allelix.annotators.base import Annotator

logger = logging.getLogger(__name__)

CONFIG_FILENAME = "config.toml"

_DEFAULT_SOURCES: dict[str, bool] = {
    "clinvar": True,
    "pharmgkb": True,
    "gwas": True,
    "gnomad": True,
    "alphamissense": True,
    "snpedia": True,
    "cadd": False,
    "rsid_identity": True,
}


# The source names a config may name, under [sources] and [license].
KNOWN_SOURCES: frozenset[str] = frozenset(_DEFAULT_SOURCES)


class ConfigFileError(ValueError):
    """config.toml cannot be read or parsed; the message names the file (2.5.2)."""


@dataclass
class AllelixConfig:
    """Persistent user configuration loaded from ``config.toml``."""

    sources: dict[str, bool] = field(default_factory=lambda: dict(_DEFAULT_SOURCES))
    commercial: bool = False
    cadd_full: bool = False
    license_overrides: dict[str, bool] = field(default_factory=dict)
    plausibility: PlausibilityThresholds = field(default_factory=PlausibilityThresholds)

    def license_held(self, source_name: str) -> bool:
        """Return whether the user asserts they hold a license for this source."""
        return self.license_overrides.get(source_name, False)

    def permission_for(
        self,
        source_name: str,
        annotator_classes: dict[str, type[Annotator]] | None = None,
    ) -> Permission | None:
        """Compute the three-state ``Permission`` for ``source_name``.

        Resolves the annotator class for ``source_name`` and applies the
        current license context (``commercial`` flag + per-source override).
        Returns ``None`` when there is no annotator class registered for the
        name. Callers that need a bool answer should use :meth:`is_enabled`
        instead. The CLI's ``config show`` calls this directly because it
        renders different table cells per permission state.
        """
        cls = None
        if annotator_classes:
            cls = annotator_classes.get(source_name)
        if cls is None:
            from allelix.annotators import get_annotator_class

            cls = get_annotator_class(source_name)

        if cls is None:
            return None
        return check_permission(
            cls.license,
            commercial=self.commercial,
            license_held=self.license_held(source_name),
        )

    def is_enabled(
        self,
        source_name: str,
        annotator_classes: dict[str, type[Annotator]] | None = None,
    ) -> bool:
        """Check if a source is enabled, respecting the permission ladder.

        Resolves the annotator class, computes the three-state
        ``Permission``, and returns ``False`` for any non-ALLOW result.
        Falls through to the source toggle for allowed sources.
        """
        perm = self.permission_for(source_name, annotator_classes)
        if perm is not None:
            if perm is not Permission.ALLOW:
                return False
        elif source_name not in _DEFAULT_SOURCES:
            return False

        return self.sources.get(source_name, True)


# Tunable plausibility cut-offs (ADR-0041), in config.toml's [plausibility].
PLAUSIBILITY_KEYS: tuple[str, ...] = tuple(f.name for f in fields(PlausibilityThresholds))


def parse_frequency(value: object) -> float | None:
    """An allele frequency in [0, 1] from a TOML value or CLI string, or None."""
    if isinstance(value, bool):
        return None
    try:
        number = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    if math.isnan(number) or not 0.0 <= number <= 1.0:
        return None
    return number


def _load_plausibility(raw: object) -> PlausibilityThresholds:
    thresholds = PlausibilityThresholds()
    if not isinstance(raw, dict):
        return thresholds
    for name in PLAUSIBILITY_KEYS:
        if name not in raw:
            continue
        number = parse_frequency(raw[name])
        if number is None:
            logger.warning(
                "config.toml: plausibility.%s = %r is not a frequency in [0, 1]; "
                "using the default",
                name,
                raw[name],
            )
            continue
        thresholds = replace(thresholds, **{name: number})
    return thresholds


def _config_path(data_dir: Path) -> Path:
    return data_dir / CONFIG_FILENAME


def _serialize(config: AllelixConfig) -> str:
    """Serialize config to TOML string. Manual formatting — no dependency."""
    lines = ["[sources]"]
    for name, enabled in sorted(config.sources.items()):
        lines.append(f"{name} = {str(enabled).lower()}")
    lines.append("")
    lines.append("[license]")
    lines.append(f"commercial = {str(config.commercial).lower()}")
    for name, held in sorted(config.license_overrides.items()):
        if held:
            lines.append(f"{name} = {str(held).lower()}")
    lines.append("")
    lines.append("[options]")
    lines.append(f"cadd_full = {str(config.cadd_full).lower()}")
    lines.append("")
    lines.append("[plausibility]")
    for name in PLAUSIBILITY_KEYS:
        lines.append(f"{name} = {getattr(config.plausibility, name)!r}")
    lines.append("")
    return "\n".join(lines)


def save_config(data_dir: Path, config: AllelixConfig) -> None:
    """Write config to disk.

    Through a staged temp, so an interrupted or failed write keeps the
    previous file (2.5.2). A destination that cannot be written raises
    ``UnwritableDestinationError`` (an ``OSError``) naming the path.
    """
    from allelix.reports import atomic_write_text, check_writable

    path = _config_path(data_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    check_writable(path)
    atomic_write_text(path, _serialize(config))


def load_config(data_dir: Path, *, create: bool = True) -> AllelixConfig:
    """Load config from disk. Creates default config if absent.

    With ``create`` False (``analyze --no-update``, ADR-0048) an absent
    config is not written; the defaults are used. A data directory that
    cannot be written gets the defaults too.
    """
    path = _config_path(data_dir)
    if not path.exists():
        config = AllelixConfig()
        if create:
            try:
                save_config(data_dir, config)
            except OSError as exc:
                logger.info("Could not create default config at %s: %s", path, exc)
            else:
                logger.info("Created default config at %s", path)
        return config

    try:
        with open(path, "rb") as fh:
            raw = tomllib.load(fh)
    except (OSError, UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        msg = (
            f"Cannot read the config file {path}: {exc}. Fix it, or delete it to "
            "start again from the defaults."
        )
        raise ConfigFileError(msg) from exc

    # Only known source names are read; any other key is ignored (2.5.2).
    sources = dict(_DEFAULT_SOURCES)
    if "sources" in raw and isinstance(raw["sources"], dict):
        for key, val in raw["sources"].items():
            if key in KNOWN_SOURCES and isinstance(val, bool):
                sources[key] = val

    commercial = False
    license_overrides: dict[str, bool] = {}
    if "license" in raw and isinstance(raw["license"], dict):
        commercial = bool(raw["license"].get("commercial", False))
        for key, val in raw["license"].items():
            if key in KNOWN_SOURCES and isinstance(val, bool) and val:
                license_overrides[key] = True

    cadd_full = False
    if "options" in raw and isinstance(raw["options"], dict):
        cadd_full = bool(raw["options"].get("cadd_full", False))

    return AllelixConfig(
        sources=sources,
        commercial=commercial,
        cadd_full=cadd_full,
        license_overrides=license_overrides,
        plausibility=_load_plausibility(raw.get("plausibility")),
    )

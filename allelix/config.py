# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""User configuration for source management and license mode."""

from __future__ import annotations

import logging
import tomllib
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from allelix.annotators.base import Permission
from allelix.annotators.base import permission as check_permission

if TYPE_CHECKING:
    from pathlib import Path

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
}


@dataclass
class AllelixConfig:
    """Persistent user configuration loaded from ``config.toml``."""

    sources: dict[str, bool] = field(default_factory=lambda: dict(_DEFAULT_SOURCES))
    commercial: bool = False
    cadd_full: bool = False
    license_overrides: dict[str, bool] = field(default_factory=dict)

    def license_held(self, source_name: str) -> bool:
        """Return whether the user asserts they hold a license for this source."""
        return self.license_overrides.get(source_name, False)

    def is_enabled(
        self,
        source_name: str,
        annotator_classes: dict[str, type] | None = None,
    ) -> bool:
        """Check if a source is enabled, respecting the permission ladder.

        Resolves the annotator class, computes the three-state
        ``Permission``, and returns ``False`` for any non-ALLOW result.
        Falls through to the source toggle for allowed sources.
        """
        cls = None
        if annotator_classes:
            cls = annotator_classes.get(source_name)
        if cls is None:
            from allelix.annotators import get_annotator_class

            cls = get_annotator_class(source_name)

        if cls is not None:
            perm = check_permission(
                cls.license,
                commercial=self.commercial,
                license_held=self.license_held(source_name),
            )
            if perm is not Permission.ALLOW:
                return False
        elif source_name not in _DEFAULT_SOURCES:
            return False

        return self.sources.get(source_name, True)


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
    return "\n".join(lines)


def save_config(data_dir: Path, config: AllelixConfig) -> None:
    """Write config to disk."""
    path = _config_path(data_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_serialize(config), encoding="utf-8")


def load_config(data_dir: Path) -> AllelixConfig:
    """Load config from disk. Creates default config if absent."""
    path = _config_path(data_dir)
    if not path.exists():
        config = AllelixConfig()
        save_config(data_dir, config)
        logger.info("Created default config at %s", path)
        return config

    with open(path, "rb") as fh:
        raw = tomllib.load(fh)

    sources = dict(_DEFAULT_SOURCES)
    if "sources" in raw and isinstance(raw["sources"], dict):
        for key, val in raw["sources"].items():
            if isinstance(val, bool):
                sources[key] = val

    commercial = False
    license_overrides: dict[str, bool] = {}
    if "license" in raw and isinstance(raw["license"], dict):
        commercial = bool(raw["license"].get("commercial", False))
        for key, val in raw["license"].items():
            if key != "commercial" and isinstance(val, bool) and val:
                license_overrides[key] = True

    cadd_full = False
    if "options" in raw and isinstance(raw["options"], dict):
        cadd_full = bool(raw["options"].get("cadd_full", False))

    return AllelixConfig(
        sources=sources,
        commercial=commercial,
        cadd_full=cadd_full,
        license_overrides=license_overrides,
    )

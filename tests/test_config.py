# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for the config file system."""

from __future__ import annotations

from typing import TYPE_CHECKING, ClassVar

from allelix.annotators.base import (
    Annotator,
    LicenseDescriptor,
    Permission,
    is_non_commercial,
    permission,
)
from allelix.config import AllelixConfig, load_config, save_config

if TYPE_CHECKING:
    from pathlib import Path


def _annotator_classes() -> dict[str, type]:
    """Build annotator class map for commercial-mode tests."""
    from allelix.annotators import (
        AlphaMissenseAnnotator,
        CaddAnnotator,
        ClinVarAnnotator,
        GnomadAnnotator,
        GWASCatalogAnnotator,
        PharmGKBAnnotator,
        SNPediaAnnotator,
    )

    return {
        cls.name: cls
        for cls in [
            ClinVarAnnotator,
            PharmGKBAnnotator,
            GWASCatalogAnnotator,
            SNPediaAnnotator,
            GnomadAnnotator,
            AlphaMissenseAnnotator,
            CaddAnnotator,
        ]
    }


class TestAllelixConfig:
    def test_default_sources_all_enabled(self):
        config = AllelixConfig()
        assert config.is_enabled("clinvar")
        assert config.is_enabled("pharmgkb")
        assert config.is_enabled("gwas")
        assert config.is_enabled("gnomad")
        assert config.is_enabled("alphamissense")
        assert config.is_enabled("snpedia")

    def test_cadd_disabled_by_default(self):
        config = AllelixConfig()
        assert not config.is_enabled("cadd")

    def test_cadd_full_defaults_false(self):
        config = AllelixConfig()
        assert not config.cadd_full

    def test_unknown_source_not_in_defaults_disabled(self):
        config = AllelixConfig()
        assert not config.is_enabled("future_source")

    def test_disabled_source(self):
        config = AllelixConfig(sources={"gnomad": False})
        assert not config.is_enabled("gnomad")

    def test_commercial_mode_disables_snpedia(self):
        config = AllelixConfig(commercial=True)
        classes = _annotator_classes()
        assert not config.is_enabled("snpedia", classes)

    def test_commercial_mode_keeps_open_sources(self):
        config = AllelixConfig(commercial=True)
        classes = _annotator_classes()
        assert config.is_enabled("clinvar", classes)
        assert config.is_enabled("pharmgkb", classes)
        assert config.is_enabled("gwas", classes)
        assert config.is_enabled("gnomad", classes)
        assert config.is_enabled("alphamissense", classes)

    def test_commercial_overrides_explicit_enable(self):
        config = AllelixConfig(
            sources={"snpedia": True},
            commercial=True,
        )
        classes = _annotator_classes()
        assert not config.is_enabled("snpedia", classes)

    def test_missing_license_raises_at_definition(self):
        """An annotator subclass without a license ClassVar fails at class definition."""
        import pytest

        with pytest.raises(TypeError, match="must declare a 'license' ClassVar"):

            class BrokenAnnotator(Annotator):
                name: ClassVar[str] = "broken"
                display_name: ClassVar[str] = "Broken"
                attribution: ClassVar[str] = "Broken"

                def setup(self) -> None: ...
                def annotate(self, variant):
                    return []

                def is_ready(self) -> bool:
                    return True

                def version(self) -> str | None:
                    return None

                def close(self) -> None: ...
                def fetch_remote_signal(self) -> str | None:
                    return None

                def cached_remote_signal(self) -> str | None:
                    return None

    def test_commercial_mode_disables_nc_without_explicit_classes(self):
        """Commercial mode blocks NC sources even without passing annotator_classes."""
        config = AllelixConfig(commercial=True)
        assert not config.is_enabled("snpedia")

    def test_nc_spdx_disabled_in_commercial_mode(self):
        """An annotator with a non-commercial SPDX is disabled in commercial mode."""

        class FakeAnnotator(Annotator):
            name: ClassVar[str] = "fake_nc"
            display_name: ClassVar[str] = "Fake NC"
            attribution: ClassVar[str] = "Fake NC"
            license: ClassVar[LicenseDescriptor] = LicenseDescriptor(
                spdx="CC-BY-NC-4.0",
                license_url="https://creativecommons.org/licenses/by-nc/4.0/",
                attribution_text="Fake.",
                source_url="https://example.com",
            )

            def setup(self) -> None: ...
            def annotate(self, variant):
                return []

            def is_ready(self) -> bool:
                return True

            def version(self) -> str | None:
                return None

            def close(self) -> None: ...
            def fetch_remote_signal(self) -> str | None:
                return None

            def cached_remote_signal(self) -> str | None:
                return None

        config = AllelixConfig(commercial=True)
        assert not config.is_enabled("fake_nc", {"fake_nc": FakeAnnotator})
        assert is_non_commercial(FakeAnnotator.license)

    def test_custom_spdx_without_commercial_ok_raises(self):
        """A custom SPDX annotator without commercial_ok fails at definition time."""
        import pytest

        with pytest.raises(TypeError, match="does not declare commercial_ok"):

            class BadCustomAnnotator(Annotator):
                name: ClassVar[str] = "bad_custom"
                display_name: ClassVar[str] = "Bad"
                attribution: ClassVar[str] = "Bad"
                license: ClassVar[LicenseDescriptor] = LicenseDescriptor(
                    spdx="LicenseRef-BadSource",
                    license_url="https://example.com/license",
                    attribution_text="Bad.",
                )

                def setup(self) -> None: ...
                def annotate(self, variant):
                    return []

                def is_ready(self) -> bool:
                    return True

                def version(self) -> str | None:
                    return None

                def close(self) -> None: ...
                def fetch_remote_signal(self) -> str | None:
                    return None

                def cached_remote_signal(self) -> str | None:
                    return None

    def test_commercial_ok_overrides_spdx_allowlist(self):
        """commercial_ok=False on a non-allowlisted SPDX is still non-commercial."""
        desc = LicenseDescriptor(
            spdx="LicenseRef-Custom",
            license_url="https://example.com",
            attribution_text="Custom.",
            commercial_ok=False,
        )
        assert is_non_commercial(desc)

    def test_commercial_ok_true_overrides_nc_spdx(self):
        """commercial_ok=True overrides even an NC SPDX (hypothetical)."""
        desc = LicenseDescriptor(
            spdx="CC-BY-NC-4.0",
            license_url="https://example.com",
            attribution_text="NC but ok.",
            commercial_ok=True,
        )
        assert not is_non_commercial(desc)

    def test_license_completeness_guard(self):
        """Every registered annotator with custom SPDX has explicit commercial_ok."""
        from allelix.annotators import _ANNOTATOR_CLASSES

        for name, cls in _ANNOTATOR_CLASSES.items():
            desc = cls.license
            if desc.spdx.startswith("LicenseRef-") or desc.spdx.startswith("custom-"):
                assert desc.commercial_ok is not None, (
                    f"{name} uses custom SPDX '{desc.spdx}' without commercial_ok"
                )

    def test_licensable_without_purchase_url_raises(self):
        """An annotator with licensable=True but no purchase_url fails at definition."""
        import pytest

        with pytest.raises(TypeError, match="licensable=True but purchase_url is None"):

            class BadLicensable(Annotator):
                name: ClassVar[str] = "bad_lic"
                display_name: ClassVar[str] = "Bad"
                attribution: ClassVar[str] = "Bad"
                license: ClassVar[LicenseDescriptor] = LicenseDescriptor(
                    spdx="LicenseRef-Bad",
                    license_url="https://example.com/license",
                    attribution_text="Bad.",
                    commercial_ok=False,
                    licensable=True,
                )

                def setup(self) -> None: ...
                def annotate(self, variant):
                    return []

                def is_ready(self) -> bool:
                    return True

                def version(self) -> str | None:
                    return None

                def close(self) -> None: ...
                def fetch_remote_signal(self) -> str | None:
                    return None

                def cached_remote_signal(self) -> str | None:
                    return None

    def test_permission_commercial_cadd_no_assertion(self):
        """Commercial mode + CADD + no license assertion → BLOCK_PURCHASABLE."""
        from allelix.annotators.cadd import CaddAnnotator

        perm = permission(CaddAnnotator.license, commercial=True, license_held=False)
        assert perm is Permission.BLOCK_PURCHASABLE

    def test_permission_commercial_cadd_with_assertion(self):
        """Commercial mode + CADD + license asserted → ALLOW."""
        from allelix.annotators.cadd import CaddAnnotator

        perm = permission(CaddAnnotator.license, commercial=True, license_held=True)
        assert perm is Permission.ALLOW

    def test_permission_commercial_snpedia_assertion_inert(self):
        """Commercial mode + SNPedia + assertion → still BLOCK_FINAL."""
        from allelix.annotators.snpedia import SNPediaAnnotator

        perm = permission(SNPediaAnnotator.license, commercial=True, license_held=True)
        assert perm is Permission.BLOCK_FINAL

    def test_permission_noncommercial_everything_allowed(self):
        """Non-commercial mode → ALLOW for all sources."""
        classes = _annotator_classes()
        for cls in classes.values():
            perm = permission(cls.license, commercial=False, license_held=False)
            assert perm is Permission.ALLOW, f"{cls.name} should be ALLOW in non-commercial mode"

    def test_is_enabled_distinguishes_toggle_from_license(self):
        """is_enabled returns False for both; permission() tells why."""
        from allelix.annotators.snpedia import SNPediaAnnotator

        classes = _annotator_classes()

        toggle_off = AllelixConfig(sources={"gnomad": False})
        assert not toggle_off.is_enabled("gnomad", classes)
        perm_gnomad = permission(classes["gnomad"].license, commercial=False, license_held=False)
        assert perm_gnomad is Permission.ALLOW

        license_block = AllelixConfig(commercial=True)
        assert not license_block.is_enabled("snpedia", classes)
        perm_snpedia = permission(SNPediaAnnotator.license, commercial=True, license_held=False)
        assert perm_snpedia is Permission.BLOCK_FINAL


class TestLoadSaveConfig:
    def test_creates_default_on_missing(self, tmp_path: Path):
        config = load_config(tmp_path)
        assert config.is_enabled("clinvar")
        assert not config.commercial
        assert (tmp_path / "config.toml").exists()

    def test_roundtrip(self, tmp_path: Path):
        original = AllelixConfig(
            sources={"clinvar": True, "gnomad": False, "snpedia": True},
            commercial=True,
        )
        save_config(tmp_path, original)
        loaded = load_config(tmp_path)
        assert not loaded.is_enabled("gnomad")
        assert loaded.commercial
        assert not loaded.is_enabled("snpedia", _annotator_classes())

    def test_cadd_full_roundtrip(self, tmp_path: Path):
        original = AllelixConfig(cadd_full=True)
        save_config(tmp_path, original)
        loaded = load_config(tmp_path)
        assert loaded.cadd_full

    def test_cadd_full_missing_defaults_false(self, tmp_path: Path):
        (tmp_path / "config.toml").write_text(
            "[sources]\nclinvar = true\n\n[license]\ncommercial = false\n"
        )
        config = load_config(tmp_path)
        assert not config.cadd_full

    def test_missing_keys_get_defaults(self, tmp_path: Path):
        (tmp_path / "config.toml").write_text(
            "[sources]\nclinvar = false\n\n[license]\ncommercial = false\n"
        )
        config = load_config(tmp_path)
        assert not config.is_enabled("clinvar")
        assert config.is_enabled("gnomad")

    def test_corrupt_values_ignored(self, tmp_path: Path):
        (tmp_path / "config.toml").write_text(
            '[sources]\nclinvar = "yes"\n\n[license]\ncommercial = false\n'
        )
        config = load_config(tmp_path)
        assert config.is_enabled("clinvar")

    def test_license_overrides_roundtrip(self, tmp_path: Path):
        original = AllelixConfig(
            commercial=True,
            license_overrides={"cadd": True},
        )
        save_config(tmp_path, original)
        loaded = load_config(tmp_path)
        assert loaded.license_held("cadd")
        assert not loaded.license_held("snpedia")

    def test_license_overrides_absent_is_false(self, tmp_path: Path):
        (tmp_path / "config.toml").write_text(
            "[sources]\nclinvar = true\n\n[license]\ncommercial = true\n"
        )
        config = load_config(tmp_path)
        assert not config.license_held("cadd")

    def test_serialized_format(self, tmp_path: Path):
        config = AllelixConfig(commercial=True)
        save_config(tmp_path, config)
        text = (tmp_path / "config.toml").read_text()
        assert "[sources]" in text
        assert "[license]" in text
        assert "commercial = true" in text

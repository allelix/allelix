# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for the parser registry and auto-detection."""

from __future__ import annotations

import pytest

from allelix.parsers import (
    PARSERS,
    ParserNotFoundError,
    detect_parser,
    get_parser_by_name,
)
from allelix.parsers.myhappygenes import MyHappyGenesParser


class TestRegistry:
    def test_mhg_parser_registered(self):
        names = [p.name for p in PARSERS]
        assert "myhappygenes" in names

    def test_get_parser_by_name(self):
        parser = get_parser_by_name("myhappygenes")
        assert isinstance(parser, MyHappyGenesParser)

    def test_get_parser_by_name_unknown_raises(self):
        with pytest.raises(ParserNotFoundError, match="Unknown parser"):
            get_parser_by_name("does_not_exist")


class TestAutoDetect:
    def test_detects_mhg(self, mock_mhg_path):
        parser = detect_parser(mock_mhg_path)
        assert parser.name == "myhappygenes"

    def test_no_match_raises(self, tmp_path):
        f = tmp_path / "garbage.txt"
        f.write_text("nothing recognizable here\n", encoding="utf-8")
        with pytest.raises(ParserNotFoundError, match="No parser recognized"):
            detect_parser(f)

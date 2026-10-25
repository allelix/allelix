# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Keep ACL restoration strict while accepting Windows inheritance bookkeeping."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest
import test_readonly_caches as readonly


@pytest.mark.parametrize("flags", ["", "P"])
def test_permission_fixture_accepts_only_added_auto_inherited_marker(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, flags: str
) -> None:
    """Exercise the real fixture's restoration path with the observed OS normalization."""
    original = f"O:SYG:BAD:{flags}(A;OICIID;FA;;;SY)(A;ID;FR;;;BA)"
    restored = original.replace(f"D:{flags}", f"D:{flags}AI", 1)
    denied = False
    calls = []
    write_bytes = Path.write_bytes

    def write(path, data):
        if denied and path.name == ".permission-denial-probe":
            raise PermissionError("synthetic ACL denial")
        return write_bytes(path, data)

    def acl(path, restore=None):
        nonlocal denied
        calls.append((path, restore))
        if restore is None:
            return original
        assert restore == original
        denied = False
        return restored

    def deny(command, **kwargs):
        nonlocal denied
        assert command == ["icacls", str(tmp_path), "/deny", "*S-1-1-0:(WD,AD,WEA,WA,DE,DC)"]
        denied = True
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(readonly, "os", SimpleNamespace(name="nt", geteuid=lambda: 1000))
    monkeypatch.setattr(readonly, "_windows_acl", acl)
    monkeypatch.setattr(readonly.subprocess, "run", deny)
    monkeypatch.setattr(Path, "write_bytes", write)
    with readonly._denied_permissions(tmp_path):
        assert denied
        assert not (tmp_path / ".permission-denial-probe").exists()
    assert not denied
    assert calls == [(tmp_path, None), (tmp_path, original)]
    assert not (tmp_path / ".permission-denial-probe").exists()
    assert os.access(tmp_path, os.W_OK)


_RULES = "(A;OICIID;FA;;;SY)(A;ID;FR;;;BA)"
_ORIGINAL = f"O:SYG:BAD:P{_RULES}"


@pytest.mark.parametrize("flags", ["", "P", "AR", "PAR", "AI", "PAI", "AIAR", "PAIAR"])
def test_identical_acl_stays_identical(flags: str) -> None:
    descriptor = f"O:SYG:BAD:{flags}{_RULES}"
    assert readonly._restored_acl_matches(descriptor, descriptor)


@pytest.mark.parametrize(
    "before,after", [("", "AI"), ("P", "PAI"), ("AR", "AIAR"), ("PAR", "PAIAR")]
)
def test_only_dacl_auto_inherited_may_be_added(before: str, after: str) -> None:
    assert readonly._restored_acl_matches(
        f"O:SYG:BAD:{before}{_RULES}", f"O:SYG:BAD:{after}{_RULES}"
    )


@pytest.mark.parametrize(
    "changed",
    [
        f"O:BAG:BAD:PAI{_RULES}",  # owner
        f"O:SYG:SYD:PAI{_RULES}",  # group
        f"O:SYG:BAD:AI{_RULES}",  # protection removed
        f"O:SYG:BAD:PAIAR{_RULES}",  # inheritance request added
        "O:SYG:BAD:PAI(A;OICIID;FR;;;SY)(A;ID;FR;;;BA)",  # rights
        "O:SYG:BAD:PAI(A;OICIID;FA;;;BA)(A;ID;FR;;;BA)",  # trustee
        "O:SYG:BAD:PAI(A;CIID;FA;;;SY)(A;ID;FR;;;BA)",  # propagation flags
        "O:SYG:BAD:PAI(A;OICI;FA;;;SY)(A;ID;FR;;;BA)",  # inherited ACE flag
        "O:SYG:BAD:PAI(A;ID;FR;;;BA)(A;OICIID;FA;;;SY)",  # ACE order
        f"O:SYG:BAD:PAI(D;;FW;;;WD){_RULES}",  # added deny
        "O:SYG:BAD:PAI(A;OICIID;FA;;;SY)",  # missing ACE
        f"O:SYG:BAD:PAIAI{_RULES}",  # duplicate marker
        f"O:SYG:BAD:PXAI{_RULES}",  # unknown control marker
        f"O:SYG:BAD:PAI{_RULES}S:(AU;SA;FR;;;WD)",  # added SACL
        "O:SYG:BA",  # missing DACL
        "O:SYG:BAD:PAI",  # empty DACL
        "O:SYG:BAD:NO_ACCESS_CONTROL",  # null DACL
    ],
)
def test_acl_normalization_rejects_security_changes(changed: str) -> None:
    assert not readonly._restored_acl_matches(_ORIGINAL, changed)


@pytest.mark.parametrize(
    "before,after",
    [
        (f"O:SYG:BAD:PAI{_RULES}", _ORIGINAL),  # AI removal
        (f"O:SYG:BAD:PAR{_RULES}", f"O:SYG:BAD:PAI{_RULES}"),  # AR removal
        (f"O:SYG:BAD:{_RULES}", f"O:SYG:BAD:PAI{_RULES}"),  # protection added
        (f"{_ORIGINAL}S:(AU;SA;FR;;;WD)", f"{_ORIGINAL}S:AI(AU;SA;FR;;;WD)"),
        ("O:SYG:BAD:", "O:SYG:BAD:AI"),  # do not normalize an empty DACL
        ("O:SYG:BAS:(AU;SA;FR;;;WD)", "O:SYG:BAS:AI(AU;SA;FR;;;WD)"),
    ],
)
def test_normalization_is_directional_and_dacl_only(before: str, after: str) -> None:
    assert not readonly._restored_acl_matches(before, after)


def test_unchanged_sacl_is_preserved_when_dacl_marker_is_added() -> None:
    sacl = "S:(AU;SA;FR;;;WD)"
    assert readonly._restored_acl_matches(_ORIGINAL + sacl, f"O:SYG:BAD:PAI{_RULES}" + sacl)

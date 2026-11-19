"""The start-up hard-link probe removes only the names it created (2.7.1).

``check_hard_links`` makes a temp in the destination's folder and links it to
a second name. If that link fails because the name already exists, the name
is another program's and must be left as it is; the probe's own temp is
removed in every case.
"""

import errno
from pathlib import Path

import pytest

from allelix import reports


def _leftovers(folder: Path, destination: Path) -> list[str]:
    return sorted(p.name for p in folder.iterdir() if p.name.startswith("." + destination.name))


@pytest.mark.parametrize("kind", ["file", "symlink", "dangling"])
def test_a_competing_name_is_kept_and_the_probe_temp_removed(tmp_path, monkeypatch, kind):
    destination = tmp_path / "out.vcf"
    target = tmp_path / "other-file"
    target.write_bytes(b"other program's content")
    claimed = []

    def competing_link(source, name, *args, **kwargs):
        name = Path(name)
        claimed.append(name)
        if kind == "file":
            name.write_bytes(b"other program's claim")
        else:
            name.symlink_to(target if kind == "symlink" else tmp_path / "absent")
        raise FileExistsError(errno.EEXIST, "File exists", str(name))

    monkeypatch.setattr(reports.os, "link", competing_link)
    with pytest.raises(OSError):
        reports.check_hard_links(destination)
    (name,) = claimed
    assert name.is_symlink() or name.exists(), "A name the probe never made is not removed"
    if kind == "file":
        assert name.read_bytes() == b"other program's claim"
    assert target.read_bytes() == b"other program's content"
    # Only the competing name is left: the probe's own temp is gone.
    assert _leftovers(tmp_path, destination) == [name.name]


def test_no_hard_links_leaves_nothing(tmp_path, monkeypatch):
    destination = tmp_path / "out.vcf"

    def unsupported(source, name, *args, **kwargs):
        raise OSError(errno.EOPNOTSUPP, "Operation not supported", str(name))

    monkeypatch.setattr(reports.os, "link", unsupported)
    with pytest.raises(reports.UnwritableDestinationError, match="hard links"):
        reports.check_hard_links(destination)
    assert _leftovers(tmp_path, destination) == []


def test_a_working_link_removes_both_names(tmp_path):
    destination = tmp_path / "out.vcf"
    reports.check_hard_links(destination)
    assert _leftovers(tmp_path, destination) == []

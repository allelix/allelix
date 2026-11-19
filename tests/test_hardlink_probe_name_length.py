"""The start-up hard-link probe accepts every name the writer can write (2.7.1).

Reviewer A's A-03: the probe linked ``.NAME.XXXXXXXX.tmp`` to a name 5 bytes
longer, so a 237-241 byte ``--vcf-out`` name the writer itself handles failed
with a raw ENAMETOOLONG naming the probe's paths. The probe's names must not
depend on the destination's length, and any other error must name the
destination, not a probe path.
"""

from pathlib import Path

import pytest

from allelix import reports


def _left(folder: Path) -> list[str]:
    return sorted(p.name for p in folder.iterdir())


@pytest.mark.parametrize("length", [200, 236, 237, 240, 241])
def test_probe_accepts_names_the_writer_can_write(tmp_path, length):
    destination = tmp_path / ("v" * (length - len(".vcf")) + ".vcf")
    assert len(destination.name) == length
    reports.check_hard_links(destination)
    assert _left(tmp_path) == []
    # The writer can write that name: stage, publish without replacing, done.
    staged = reports.StagedOutput(destination)
    with staged.open_text() as stream:
        stream.write("##fileformat=VCFv4.2\n")
    staged.publish(replace_new=False, replace_existing=False)
    assert destination.read_text() == "##fileformat=VCFv4.2\n"


def test_other_probe_errors_name_the_destination_not_a_probe_path(tmp_path, monkeypatch):
    destination = tmp_path / "out.vcf"

    def denied(source, name, *args, **kwargs):
        raise PermissionError(13, "Permission denied", str(name))

    # EACCES is not a "no hard links" errno: it must still read as a write error.
    monkeypatch.setattr(reports.os, "link", denied)
    with pytest.raises(reports.UnwritableDestinationError) as caught:
        reports.check_hard_links(destination)
    message = str(caught.value)
    assert message.startswith(f"Cannot write {destination}")
    assert ".tmp" not in message and ".link" not in message, message
    assert _left(tmp_path) == []

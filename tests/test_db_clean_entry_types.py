"""A filename owned by a cache writer is not ownership of a directory tree."""

import pytest
from click.testing import CliRunner

from allelix.cli import main


@pytest.mark.parametrize(
    "name",
    [
        "gnomad.sqlite",
        "gwas_catalog_associations.tsv",
        "clinvar-build.log",
        "gnomad.sqlite.tmp",
        ".gnomad.sqlite.0123456789abcdef.tmp",
    ],
)
@pytest.mark.parametrize("force", [False, True])
def test_cache_file_name_does_not_authorize_directory_deletion(tmp_path, name, force):
    # Allelix's writers create files at these names, never directory trees.
    # A name collision in a user-selected directory must not recurse into it.
    directory = tmp_path / name
    directory.mkdir()
    evidence = directory / "personal-report.txt"
    evidence.write_bytes(b"user data\n")
    (tmp_path / "snpedia.sqlite").write_bytes(b"disposable cache")
    args = ["db", "clean", "--data-dir", str(tmp_path), "--yes"]
    if force:
        args.append("--force")
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, result.output
    assert evidence.is_file(), "a file-name match authorized recursive deletion of user data"
    assert evidence.read_bytes() == b"user data\n"
    assert not (tmp_path / "snpedia.sqlite").exists()
    assert name in " ".join(result.output.split())


def test_cache_symlink_removes_only_the_link(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    evidence = outside / "personal-report.txt"
    evidence.write_bytes(b"user data\n")
    link = data / "gnomad.sqlite"
    link.symlink_to(outside, target_is_directory=True)
    result = CliRunner().invoke(main, ["db", "clean", "--data-dir", str(data), "--yes"])
    assert result.exit_code == 0, result.output
    assert evidence.read_bytes() == b"user data\n"
    assert not link.is_symlink()


def test_real_owned_staging_directory_is_removed(tmp_path):
    directory = tmp_path / ".gnomad-abc12345"
    directory.mkdir()
    (directory / "gnomad.sqlite.gz").write_bytes(b"cache staging")
    result = CliRunner().invoke(main, ["db", "clean", "--data-dir", str(tmp_path), "--yes"])
    assert result.exit_code == 0, result.output
    assert not directory.exists()

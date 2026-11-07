"""A user-supplied backend that cannot be read is not a successful update."""

import os

import pytest
from click.testing import CliRunner
from test_cadd_full_mode_files import _place_full_mode_files, full_mode  # noqa: F401

from allelix.annotators.cadd import CaddAnnotator
from allelix.cli import main
from allelix.databases.cache_state import CacheState


@pytest.mark.parametrize("name", ["whole_genome_SNVs.tsv.gz", "whole_genome_SNVs.tsv.gz.tbi"])
@pytest.mark.parametrize("force", [False, True])
def test_unreadable_full_mode_backend_fails_without_download(
    full_mode,  # noqa: F811
    monkeypatch,
    name,
    force,
):
    data, downloads = full_mode
    _place_full_mode_files(data)
    inaccessible = data / name
    real_access = os.access
    monkeypatch.setattr(
        os,
        "access",
        lambda path, mode, *args, **kwargs: (
            False
            if path == inaccessible and mode == os.R_OK
            else real_access(path, mode, *args, **kwargs)
        ),
    )
    source = CaddAnnotator(data, full_mode=True)
    assert source.cache_status().state is CacheState.UNREADABLE
    args = ["db", "update", "--cadd", "--data-dir", str(data)]
    if force:
        args.append("--force")
    result = CliRunner().invoke(main, args)
    assert downloads == []
    assert result.exit_code == 1, result.output
    assert name in " ".join(result.output.split()), result.output
    assert "cadd: ready" not in result.output

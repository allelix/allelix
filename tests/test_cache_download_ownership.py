"""O5: setup owns downloaded bytes until verification, ingestion and cleanup finish.

Fake transports supply tiny synthetic payloads. Setup and its filesystem paths
are production; lower-level database installation races are covered separately.
"""

import importlib
import zipfile
from pathlib import Path

import pytest


@pytest.mark.parametrize(
    "source,class_name",
    [
        ("gnomad", "GnomadAnnotator"),
        ("cadd", "CaddAnnotator"),
        ("alphamissense", "AlphaMissenseAnnotator"),
        ("snpedia", "SNPediaAnnotator"),
        ("rsid_identity", "RsidIdentityAnnotator"),
    ],
)
@pytest.mark.parametrize("inner_fails", [False, True])
def test_prebuilt_setup_keeps_its_download_until_installed(
    tmp_path, monkeypatch, source, class_name, inner_fails
):
    module = importlib.import_module("allelix.annotators." + source)
    cls = getattr(module, class_name)
    paths = []
    installed = []
    inside = False
    nested = False
    errors = []

    def download(url, path):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"synthetic verified gzip")
        paths.append(path)

    def install(path, destination, **kwargs):
        if inside and inner_fails:
            raise OSError("synthetic installer failure")
        assert path.read_bytes() == b"synthetic verified gzip"
        installed.append("inner" if inside else "outer")

    def verify(path, *args):
        nonlocal inside, nested
        if not nested:
            nested = True
            inside = True
            try:
                with cls(tmp_path) as other:
                    other.setup()
            except BaseException as exc:
                errors.append(exc)
            finally:
                inside = False
        assert path.read_bytes() == b"synthetic verified gzip"

    monkeypatch.setattr(module, "download", download)
    monkeypatch.setattr(module, "verify_file_hash", verify)
    monkeypatch.setattr(module, "install_prebuilt_cache", install)
    if source == "snpedia":
        monkeypatch.setattr(cls, "prepare_cache", lambda self: None)
        # 2.6.1: setup parses and publishes the staged file after the install;
        # the install is faked here, so there is nothing to parse.
        monkeypatch.setattr(cls, "_parse_and_publish", lambda self, staged: None)
    outer_error = None
    try:
        with cls(tmp_path) as outer:
            outer.setup()
    except BaseException as exc:
        outer_error = exc
    assert nested and len(paths) == 2
    assert outer_error is None, f"another setup deleted the verified download: {outer_error!r}"
    assert installed == (["outer"] if inner_fails else ["inner", "outer"])
    assert len(errors) == int(inner_fails)
    if inner_fails:
        assert isinstance(errors[0], OSError)


@pytest.mark.parametrize(
    "source,class_name", [("gwas", "GWASCatalogAnnotator"), ("pharmgkb", "PharmGKBAnnotator")]
)
@pytest.mark.parametrize("inner_fails", [False, True])
def test_mutable_source_ingests_the_bytes_its_own_download_delivered(
    tmp_path, monkeypatch, source, class_name, inner_fails
):
    module = importlib.import_module("allelix.annotators." + source)
    cls = getattr(module, class_name)
    inside = False
    nested = False
    installed = []
    errors = []
    paths = []

    def download(url, path):
        path.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr("rows.tsv", "B" if inside else "A")
        paths.append(path)

    def load(path, destination, **kwargs):
        nonlocal inside, nested
        my_label = "B" if inside else "A"
        if not nested:
            nested = True
            inside = True
            try:
                with cls(tmp_path) as other:
                    other.setup()
            except BaseException as exc:
                errors.append(exc)
            finally:
                inside = False
        if my_label == "B" and inner_fails:
            raise OSError("synthetic ingestion failure")
        if source == "gwas":
            content = Path(path).read_text()
        else:
            with zipfile.ZipFile(path) as archive:
                content = archive.read("rows.tsv").decode()
        assert content == my_label, "setup ingested another refresh's source bytes"
        installed.append(my_label)
        return 1

    monkeypatch.setattr(cls, "fetch_remote_signal", lambda self: 'etag:"synthetic"')
    monkeypatch.setattr(module, "download", download)
    monkeypatch.setattr(module, "load_gwas_tsv" if source == "gwas" else "load_pharmgkb_tsv", load)
    if source == "pharmgkb":
        monkeypatch.setattr(module, "fetch_cpic_allele_functions", lambda: {})
        monkeypatch.setattr(module, "pharmgkb_content_hash", lambda *args: None)
    outer_error = None
    try:
        with cls(tmp_path) as outer:
            outer.setup()
    except BaseException as exc:
        outer_error = exc
    assert nested and len(paths) == 2
    assert outer_error is None, f"own downloaded input was overwritten: {outer_error!r}"
    assert installed == (["A"] if inner_fails else ["B", "A"])
    assert len(errors) == int(inner_fails)

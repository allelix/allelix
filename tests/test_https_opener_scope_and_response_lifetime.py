"""A transport policy owns its opener and closes replies it refuses."""

import io
import socket
import urllib.request
import urllib.response
from email.message import Message

import pytest
from test_https_redirect_policy import BODY, PATHS, _invoke, _url

from allelix.databases import cpic_loader


@pytest.fixture
def replies(monkeypatch):
    route = {}
    seen = []
    opened = []

    def no_network(*args, **kwargs):
        pytest.fail("The synthetic redirect check tried to use a real socket")

    def transport(self, request, data=None):
        seen.append(request.full_url)
        location = route.get(request.full_url)
        status = 302 if location else 200
        headers = Message()
        headers["Content-Length"] = str(0 if location else len(BODY))
        headers["Content-Range"] = "0-0/1"
        if location:
            headers["Location"] = location
        response = urllib.response.addinfourl(
            io.BytesIO(b"" if location else BODY), headers, request.full_url, status
        )
        response.msg = "Redirect" if location else "OK"
        opened.append(response)
        return response

    monkeypatch.setattr(socket.socket, "connect", no_network)
    monkeypatch.setattr(urllib.request.OpenerDirector, "_open", transport)
    monkeypatch.setattr(cpic_loader.time, "sleep", lambda _: None)
    # Isolate urllib state even if the implementation changes it.
    monkeypatch.setattr(urllib.request, "_opener", None)
    try:
        yield route, seen, opened
    finally:
        for response in opened:
            response.close()


@pytest.mark.parametrize("kind", PATHS)
def test_a_source_request_does_not_replace_the_process_opener(
    tmp_path, monkeypatch, replies, kind
):
    previous = urllib.request.build_opener()
    monkeypatch.setattr(urllib.request, "_opener", previous)
    result = _invoke(kind, tmp_path / "download")
    _, seen, opened = replies
    assert seen == [_url(kind)], "The controlled successful request actually ran"
    if kind == "download":
        assert (tmp_path / "download").read_bytes() == BODY
    else:
        assert result is not None
    assert all(response.closed for response in opened)
    assert urllib.request._opener is previous, (
        "A source request must not replace another user's urllib opener; "
        "use the source policy's own opener directly"
    )


@pytest.mark.parametrize("kind", PATHS)
@pytest.mark.parametrize("scheme", ("http", "ftp"))
def test_a_refused_redirect_closes_its_response(tmp_path, replies, kind, scheme):
    route, seen, opened = replies
    route[_url(kind)] = f"{scheme}://insecure.invalid/payload"
    if kind in ("download", "cpic_table"):
        with pytest.raises(OSError):
            _invoke(kind, tmp_path / "download")
    else:
        assert _invoke(kind, tmp_path / "download") is None
    assert opened, "The controlled redirect response was received"
    assert set(seen) == {_url(kind)}, "No forbidden target was contacted"
    assert all(response.closed for response in opened), (
        "Refusal must close every redirect reply, including each CPIC retry; "
        "URLError does not own the response as HTTPError does"
    )
    assert not (tmp_path / "download").exists()

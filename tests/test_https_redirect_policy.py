"""All cache transport paths reject a downgrade before contacting its target."""

import io
import json
import socket
import urllib.request
import urllib.response
from email.message import Message
from urllib.parse import urlsplit

import pytest

from allelix.databases import cpic_loader, manager

PATHS = ("text", "head", "download", "cpic_signal", "cpic_table")
BODY = json.dumps([{"date": "2026-11-10", "definitionid": 1}]).encode()


@pytest.fixture
def transport(monkeypatch):
    """Exercise real urllib redirects; replace only the final network transport."""
    seen = []
    route = {}

    def no_network(*args, **kwargs):
        pytest.fail("The offline redirect test attempted a real socket connection")

    def open_reply(self, request, data=None):
        url = request.full_url
        seen.append((url, request.get_method(), dict(request.header_items()), request.timeout))
        status, location = route.get(url, (200, None))
        headers = Message()
        headers["Content-Length"] = str(len(BODY) if status == 200 else 0)
        headers["Content-Range"] = "0-0/1"
        headers["ETag"] = '"synthetic"'
        if location is not None:
            headers["Location"] = location
        reply = urllib.response.addinfourl(
            io.BytesIO(BODY if status == 200 else b""), headers, url, status
        )
        reply.msg = "OK" if status == 200 else "Redirect"
        return reply

    monkeypatch.setattr(socket.socket, "connect", no_network)
    monkeypatch.setattr(urllib.request, "_opener", None)
    monkeypatch.setattr(urllib.request.OpenerDirector, "_open", open_reply)
    monkeypatch.setattr(cpic_loader.time, "sleep", lambda _: None)
    return route, seen


def _url(kind):
    if kind == "cpic_signal":
        return "https://source.invalid/change_log?select=date&order=date.desc&limit=1"
    if kind == "cpic_table":
        return "https://source.invalid/allele"
    return "https://source.invalid/cache"


def _invoke(kind, dest):
    if kind == "text":
        return manager.fetch_remote_text(_url(kind))
    if kind == "head":
        return manager.head_request_headers(_url(kind))
    if kind == "download":
        return manager.download(_url(kind), dest)
    if kind == "cpic_signal":
        return cpic_loader.fetch_cpic_remote_signal("https://source.invalid")
    return cpic_loader._http_get_json(_url(kind))


@pytest.mark.parametrize("kind", PATHS)
@pytest.mark.parametrize("status", (301, 302, 303, 307, 308))
def test_http_redirect_is_refused_before_target_contact(tmp_path, transport, kind, status):
    route, seen = transport
    route[_url(kind)] = status, "http://insecure.invalid/payload"
    _assert_refusal(tmp_path, kind, seen)


@pytest.mark.parametrize("kind", PATHS)
def test_later_redirect_hop_cannot_downgrade(tmp_path, transport, kind):
    route, seen = transport
    route[_url(kind)] = 302, "https://cdn.invalid/first"
    route["https://cdn.invalid/first"] = 307, "http://insecure.invalid/payload"
    _assert_refusal(tmp_path, kind, seen)
    assert any(url == "https://cdn.invalid/first" for url, *_ in seen)


@pytest.mark.parametrize("kind", PATHS)
@pytest.mark.parametrize("scheme", ("ftp", "file", "data", "javascript"))
def test_other_non_https_redirect_schemes_are_also_refused(tmp_path, transport, kind, scheme):
    route, seen = transport
    route[_url(kind)] = 302, f"{scheme}://insecure.invalid/payload"
    _assert_refusal(tmp_path, kind, seen)


def _assert_refusal(tmp_path, kind, seen):
    dest = tmp_path / "retained.data"
    dest.write_bytes(b"prior verified data")
    dest.chmod(0o640)
    before = dest.stat()
    error = None
    result = None
    try:
        result = _invoke(kind, dest)
    except (OSError, ValueError) as exc:
        error = exc
    assert seen, "The HTTPS request must reach the controlled transport"
    assert all(urlsplit(url).scheme.lower() == "https" for url, *_ in seen), (
        "A redirect must be rejected before reading any bytes from its non-HTTPS target",
        seen,
    )
    if kind in ("text", "head", "cpic_signal"):
        assert error is None and result is None, "Freshness failure stays unavailable, not fatal"
    else:
        assert error is not None, "An installer must not report a refused download as success"
    assert dest.read_bytes() == b"prior verified data"
    assert (dest.stat().st_ino, dest.stat().st_mode) == (before.st_ino, before.st_mode)
    assert sorted(p.name for p in tmp_path.iterdir()) == [dest.name]


@pytest.mark.parametrize("kind", PATHS)
@pytest.mark.parametrize(
    "target", (None, "https://cdn.invalid/payload", "/relative", "//cdn.invalid/payload")
)
def test_https_and_relative_redirects_keep_requests_and_results(tmp_path, transport, kind, target):
    route, seen = transport
    if target is not None:
        route[_url(kind)] = 302, target
    dest = tmp_path / "output.data"
    result = _invoke(kind, dest)
    assert len(seen) == (1 if target is None else 2)
    for url, method, headers, timeout in seen:
        assert urlsplit(url).scheme == "https"
        assert method == ("HEAD" if kind == "head" else "GET")
        normalized = {key.lower(): value for key, value in headers.items()}
        assert normalized["user-agent"] == manager.USER_AGENT
        assert timeout > 0
        if kind.startswith("cpic"):
            assert normalized["accept"] == "application/json"
        if kind == "cpic_table":
            assert normalized["range"] == f"0-{cpic_loader.CPIC_MAX_ROWS}"
            assert normalized["prefer"] == "count=exact"
    if kind == "text":
        assert result == BODY.decode()
    elif kind == "head":
        assert result["etag"] == '"synthetic"'
    elif kind == "download":
        assert dest.read_bytes() == BODY
    elif kind == "cpic_signal":
        assert result == "lastchange:2026-11-10"
    else:
        assert result == json.loads(BODY)


@pytest.mark.parametrize("kind", PATHS)
@pytest.mark.parametrize("status", (301, 302, 303, 307, 308))
@pytest.mark.parametrize("legacy_parent", (False, True), ids=("native-urllib", "legacy-urllib"))
def test_allelix_preserves_method_on_every_https_hop(
    tmp_path, transport, monkeypatch, kind, status, legacy_parent
):
    """Allelix owns the method policy even if its urllib parent changes HEAD to GET.

    Python 3.12's parent discards HEAD; Python 3.14's parent preserves it.
    The legacy-parent control exercises the former behavior on both versions,
    so a newer interpreter cannot hide reliance on the parent's policy.
    """
    if legacy_parent:
        parent = urllib.request.HTTPRedirectHandler.redirect_request

        def legacy_redirect(self, req, fp, code, msg, headers, newurl):
            redirected = parent(self, req, fp, code, msg, headers, newurl)
            if redirected is not None and req.get_method() == "HEAD":
                redirected.method = "GET"
            return redirected

        monkeypatch.setattr(
            urllib.request.HTTPRedirectHandler, "redirect_request", legacy_redirect
        )

    route, seen = transport
    route[_url(kind)] = status, "https://cdn.invalid/first"
    route["https://cdn.invalid/first"] = 307, "/second"
    route["https://cdn.invalid/second"] = 308, "https://final.invalid/payload"
    dest = tmp_path / "output.data"
    result = _invoke(kind, dest)
    assert [url for url, *_ in seen] == [
        _url(kind),
        "https://cdn.invalid/first",
        "https://cdn.invalid/second",
        "https://final.invalid/payload",
    ]
    expected = "HEAD" if kind == "head" else "GET"
    assert [method for _, method, *_ in seen] == [expected] * 4, (
        "Allelix must preserve the original HEAD/GET method at every HTTPS redirect hop",
        seen,
    )
    if kind == "head":
        assert result["etag"] == '"synthetic"'
    elif kind == "download":
        assert dest.read_bytes() == BODY
    elif kind == "text":
        assert result == BODY.decode()
    elif kind == "cpic_signal":
        assert result == "lastchange:2026-11-10"
    else:
        assert result == json.loads(BODY)

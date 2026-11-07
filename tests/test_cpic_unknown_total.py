"""A server can cap a response below the requested window; a short body proves nothing."""

import pytest
from test_cpic_complete_response import _Reply, cpic  # noqa: F401

from allelix.databases import cpic_loader


@pytest.mark.parametrize(
    "endpoint", ["/sequence_location?", "/allele_location_value?", "/allele?"]
)
@pytest.mark.parametrize("status", [200, 206])
@pytest.mark.parametrize("content_range", ["unknown", None])
def test_short_unknown_total_is_not_a_verified_whole_table(
    cpic,  # noqa: F811
    endpoint,
    status,
    content_range,
):
    replies, _ = cpic

    def partial(rows):
        kept = rows[:-1]
        header = (f"0-{len(kept) - 1}/*" if kept else "*/*") if content_range else None
        return _Reply(kept, status=status, content_range=header)

    replies[endpoint] = partial
    with pytest.raises(ValueError):
        cpic_loader.fetch_cpic_allele_functions()

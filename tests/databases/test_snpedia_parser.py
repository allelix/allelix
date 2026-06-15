# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Tests for the SNPedia raw-to-structured parser."""

from __future__ import annotations

import contextlib
import sqlite3

from allelix.databases.snpedia_parser import (
    _PARSER_VERSION,
    _dedupe_existing,
    detect_raw_table,
    has_structured_table,
    parse_raw_pages,
    parser_is_current,
)


def _make_db(
    tmp_path,
    table_name="pages",
    snp_pages=None,
    genotype_pages=None,
):
    """Create a SQLite db with raw pages for parsing."""
    db = tmp_path / "snpedia.sqlite"
    with contextlib.closing(sqlite3.connect(db)) as conn:
        conn.execute(
            f"CREATE TABLE {table_name} "
            "(title TEXT PRIMARY KEY, category TEXT, content TEXT, scraped_at TEXT)"
        )
        for title, content in snp_pages or []:
            conn.execute(
                f"INSERT INTO {table_name} VALUES (?, 'snp', ?, '2026-01-01T00:00:00')",
                (title, content),
            )
        for title, content in genotype_pages or []:
            conn.execute(
                f"INSERT INTO {table_name} VALUES (?, 'genotype', ?, '2026-01-01T00:00:00')",
                (title, content),
            )
        conn.commit()
    return str(db)


class TestDetectRawTable:
    """Identify which raw-pages table exists in the SQLite file."""

    def test_raw_pages_table(self, tmp_path):
        db = tmp_path / "t.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.execute("CREATE TABLE _raw_pages (title TEXT)")
            assert detect_raw_table(conn) == "_raw_pages"

    def test_pages_table(self, tmp_path):
        db = tmp_path / "t.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.execute("CREATE TABLE pages (title TEXT)")
            assert detect_raw_table(conn) == "pages"

    def test_no_raw_table(self, tmp_path):
        db = tmp_path / "t.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.execute("CREATE TABLE other (id INTEGER)")
            assert detect_raw_table(conn) is None

    def test_prefers_raw_pages_over_pages(self, tmp_path):
        db = tmp_path / "t.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.execute("CREATE TABLE _raw_pages (title TEXT)")
            conn.execute("CREATE TABLE pages (title TEXT)")
            assert detect_raw_table(conn) == "_raw_pages"


class TestHasStructuredTable:
    """Check for presence and non-emptiness of snpedia_genotypes."""

    def test_has_rows(self, tmp_path):
        db = tmp_path / "t.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.execute(
                "CREATE TABLE snpedia_genotypes "
                "(rsid TEXT, allele1 TEXT, allele2 TEXT, magnitude REAL, "
                "repute TEXT, summary TEXT, gene TEXT, scraped_at TEXT)"
            )
            conn.execute(
                "INSERT INTO snpedia_genotypes VALUES ('rs1','A','A',0,NULL,'x',NULL,'2026-01-01')"
            )
            conn.commit()
            assert has_structured_table(conn) is True

    def test_empty_table(self, tmp_path):
        db = tmp_path / "t.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.execute(
                "CREATE TABLE snpedia_genotypes "
                "(rsid TEXT, allele1 TEXT, allele2 TEXT, magnitude REAL, "
                "repute TEXT, summary TEXT, gene TEXT, scraped_at TEXT)"
            )
            assert has_structured_table(conn) is False

    def test_no_table(self, tmp_path):
        db = tmp_path / "t.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            assert has_structured_table(conn) is False


class TestParseRawPages:
    """End-to-end parsing of raw wiki markup into structured rows."""

    def test_no_raw_table_returns_zero(self, tmp_path):
        db = tmp_path / "snpedia.sqlite"
        with contextlib.closing(sqlite3.connect(db)) as conn:
            conn.execute("CREATE TABLE other (id INTEGER)")
        assert parse_raw_pages(str(db)) == 0

    def test_parse_basic_genotype(self, tmp_path):
        snp_pages = [("Rs1801133", "{{rsnum\n|Gene=MTHFR\n}}")]
        genotype_pages = [
            (
                "Rs1801133(C;T)",
                "{{Genotype\n|allele1=C\n|allele2=T\n"
                "|magnitude=2.2\n|repute=Bad\n|summary=1 copy of C677T\n}}",
            ),
        ]
        db_path = _make_db(tmp_path, snp_pages=snp_pages, genotype_pages=genotype_pages)
        count = parse_raw_pages(db_path)
        assert count == 1

        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            row = conn.execute("SELECT * FROM snpedia_genotypes").fetchone()
            rsid, a1, a2, mag, repute, summary, gene, _ = row
            assert rsid == "rs1801133"
            assert a1 == "C"
            assert a2 == "T"
            assert mag == 2.2
            assert repute == "Bad"
            assert summary == "1 copy of C677T"
            assert gene == "MTHFR"

    def test_allele_sorting(self, tmp_path):
        genotype_pages = [
            ("Rs100(T;C)", "{{Genotype|allele1=T|allele2=C|magnitude=1.0|summary=test}}"),
        ]
        db_path = _make_db(tmp_path, genotype_pages=genotype_pages)
        parse_raw_pages(db_path)
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            row = conn.execute("SELECT allele1, allele2 FROM snpedia_genotypes").fetchone()
            assert row == ("C", "T")

    def test_missing_magnitude(self, tmp_path):
        genotype_pages = [
            ("Rs100(A;G)", "{{Genotype|allele1=A|allele2=G|summary=test}}"),
        ]
        db_path = _make_db(tmp_path, genotype_pages=genotype_pages)
        parse_raw_pages(db_path)
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            assert conn.execute("SELECT magnitude FROM snpedia_genotypes").fetchone()[0] is None

    def test_invalid_magnitude(self, tmp_path):
        genotype_pages = [
            ("Rs100(A;G)", "{{Genotype|allele1=A|allele2=G|magnitude=NaN|summary=test}}"),
        ]
        db_path = _make_db(tmp_path, genotype_pages=genotype_pages)
        parse_raw_pages(db_path)
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            assert conn.execute("SELECT magnitude FROM snpedia_genotypes").fetchone()[0] is None

    def test_parses_i_probe_genotype(self, tmp_path):
        """I-prefixed 23andMe probe IDs are parsed into structured rows."""
        genotype_pages = [
            ("I5006212(A;G)", "{{Genotype|allele1=A|allele2=G|summary=carrier}}"),
        ]
        db_path = _make_db(tmp_path, genotype_pages=genotype_pages)
        assert parse_raw_pages(db_path) == 1
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            row = conn.execute(
                "SELECT rsid, allele1, allele2, summary FROM snpedia_genotypes"
            ).fetchone()
            assert row == ("i5006212", "A", "G", "carrier")

    def test_i_probe_gene_from_23andme_snp_template(self, tmp_path):
        """Gene mapped from {{23andMe SNP}} template on I-probe SNP pages."""
        snp_pages = [("I3000001", "{{23andMe SNP\n|Gene_s=CFTR\n}}")]
        genotype_pages = [
            ("I3000001(A;G)", "{{Genotype|allele1=A|allele2=G|summary=carrier}}"),
        ]
        db_path = _make_db(tmp_path, snp_pages=snp_pages, genotype_pages=genotype_pages)
        parse_raw_pages(db_path)
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            assert conn.execute("SELECT gene FROM snpedia_genotypes").fetchone()[0] == "CFTR"

    def test_skips_non_rs_non_i_title(self, tmp_path):
        genotype_pages = [("X12345(A;G)", "{{Genotype|allele1=A|allele2=G|summary=x}}")]
        assert parse_raw_pages(_make_db(tmp_path, genotype_pages=genotype_pages)) == 0

    def test_skips_empty_content(self, tmp_path):
        genotype_pages = [("Rs100(A;G)", "")]
        assert parse_raw_pages(_make_db(tmp_path, genotype_pages=genotype_pages)) == 0

    def test_title_fallback_when_template_alleles_empty(self, tmp_path):
        """Alleles missing from template are extracted from the page title."""
        genotype_pages = [
            ("Rs100(A;G)", "{{Genotype|magnitude=7|repute=Bad|summary=test}}"),
        ]
        db_path = _make_db(tmp_path, genotype_pages=genotype_pages)
        assert parse_raw_pages(db_path) == 1
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            row = conn.execute(
                "SELECT allele1, allele2, summary FROM snpedia_genotypes"
            ).fetchone()
            assert row == ("A", "G", "test")

    def test_title_fallback_sorts_alleles(self, tmp_path):
        """Alleles extracted from title are still sorted."""
        genotype_pages = [
            ("Rs100(T;C)", "{{Genotype|magnitude=0|summary=test}}"),
        ]
        db_path = _make_db(tmp_path, genotype_pages=genotype_pages)
        parse_raw_pages(db_path)
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            row = conn.execute("SELECT allele1, allele2 FROM snpedia_genotypes").fetchone()
            assert row == ("C", "T")

    def test_skips_empty_alleles_and_unparseable_title(self, tmp_path):
        """No alleles in template AND title has empty alleles — skip."""
        genotype_pages = [
            ("Rs100(;)", "{{Genotype|magnitude=0|summary=test}}"),
        ]
        assert parse_raw_pages(_make_db(tmp_path, genotype_pages=genotype_pages)) == 0

    def test_skips_no_genotype_template(self, tmp_path):
        genotype_pages = [("Rs100(A;G)", "{{OtherTemplate|allele1=A|allele2=G}}")]
        assert parse_raw_pages(_make_db(tmp_path, genotype_pages=genotype_pages)) == 0

    def test_null_repute_and_summary(self, tmp_path):
        genotype_pages = [("Rs100(A;G)", "{{Genotype|allele1=A|allele2=G|magnitude=1.0}}")]
        db_path = _make_db(tmp_path, genotype_pages=genotype_pages)
        parse_raw_pages(db_path)
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            row = conn.execute("SELECT repute, summary FROM snpedia_genotypes").fetchone()
            assert row == (None, None)

    def test_database_versions_row(self, tmp_path):
        genotype_pages = [
            ("Rs100(A;G)", "{{Genotype|allele1=A|allele2=G|summary=test}}"),
        ]
        db_path = _make_db(tmp_path, genotype_pages=genotype_pages)
        parse_raw_pages(db_path)
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            row = conn.execute(
                "SELECT name, source_url, record_count FROM database_versions WHERE name='snpedia'"
            ).fetchone()
            assert row[0] == "snpedia"
            assert "snpedia.com" in row[1]
            assert row[2] == 1

    def test_verbose_mode(self, tmp_path):
        genotype_pages = [
            ("Rs100(A;G)", "{{Genotype|allele1=A|allele2=G|summary=test}}"),
        ]
        db_path = _make_db(tmp_path, genotype_pages=genotype_pages)
        assert parse_raw_pages(db_path, verbose=True) == 1

    def test_raw_pages_table_name(self, tmp_path):
        genotype_pages = [
            ("Rs100(A;G)", "{{Genotype|allele1=A|allele2=G|summary=test}}"),
        ]
        db_path = _make_db(tmp_path, table_name="_raw_pages", genotype_pages=genotype_pages)
        assert parse_raw_pages(db_path) == 1

    def test_gene_from_snp_template(self, tmp_path):
        snp_pages = [("Rs4680", "{{snp\n|Gene=COMT\n}}")]
        genotype_pages = [
            ("Rs4680(A;G)", "{{Genotype|allele1=A|allele2=G|summary=intermediate}}"),
        ]
        db_path = _make_db(tmp_path, snp_pages=snp_pages, genotype_pages=genotype_pages)
        parse_raw_pages(db_path)
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            assert conn.execute("SELECT gene FROM snpedia_genotypes").fetchone()[0] == "COMT"

    def test_snp_page_without_gene(self, tmp_path):
        snp_pages = [("Rs100", "{{rsnum}}")]
        genotype_pages = [
            ("Rs100(A;G)", "{{Genotype|allele1=A|allele2=G|summary=test}}"),
        ]
        db_path = _make_db(tmp_path, snp_pages=snp_pages, genotype_pages=genotype_pages)
        parse_raw_pages(db_path)
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            assert conn.execute("SELECT gene FROM snpedia_genotypes").fetchone()[0] is None

    def test_snp_page_non_rs_title_skipped(self, tmp_path):
        snp_pages = [("SomeOtherPage", "{{rsnum|Gene=FOO}}")]
        genotype_pages = [
            ("Rs100(A;G)", "{{Genotype|allele1=A|allele2=G|summary=test}}"),
        ]
        db_path = _make_db(tmp_path, snp_pages=snp_pages, genotype_pages=genotype_pages)
        parse_raw_pages(db_path)
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            assert conn.execute("SELECT gene FROM snpedia_genotypes").fetchone()[0] is None

    def test_reparse_clears_old_data(self, tmp_path):
        genotype_pages = [
            ("Rs100(A;G)", "{{Genotype|allele1=A|allele2=G|summary=test}}"),
        ]
        db_path = _make_db(tmp_path, genotype_pages=genotype_pages)
        parse_raw_pages(db_path)
        parse_raw_pages(db_path)
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            assert conn.execute("SELECT COUNT(*) FROM snpedia_genotypes").fetchone()[0] == 1

    def test_multiple_genotypes_with_gene_map(self, tmp_path):
        snp_pages = [("Rs100", "{{rsnum|Gene=MTHFR}}")]
        genotype_pages = [
            (
                "Rs100(A;A)",
                "{{Genotype|allele1=A|allele2=A|magnitude=0|repute=Good|summary=normal}}",
            ),
            (
                "Rs100(A;G)",
                "{{Genotype|allele1=A|allele2=G|magnitude=2.0|repute=Bad|summary=het}}",
            ),
            (
                "Rs100(G;G)",
                "{{Genotype|allele1=G|allele2=G|magnitude=3.0|repute=Bad|summary=hom}}",
            ),
        ]
        db_path = _make_db(tmp_path, snp_pages=snp_pages, genotype_pages=genotype_pages)
        assert parse_raw_pages(db_path) == 3

        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            rows = conn.execute("SELECT gene FROM snpedia_genotypes").fetchall()
            assert all(r[0] == "MTHFR" for r in rows)

    def test_snp_page_empty_content_skipped(self, tmp_path):
        snp_pages = [("Rs100", "")]
        genotype_pages = [
            ("Rs100(A;G)", "{{Genotype|allele1=A|allele2=G|summary=test}}"),
        ]
        db_path = _make_db(tmp_path, snp_pages=snp_pages, genotype_pages=genotype_pages)
        parse_raw_pages(db_path)
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            assert conn.execute("SELECT gene FROM snpedia_genotypes").fetchone()[0] is None

    def test_parser_dedupes_identical_genotype_rows(self, tmp_path):
        """Two SNPedia pages with whitespace-differing titles but same genotype collapse."""
        genotype_pages = [
            (
                "Rs4950928(C;G)",
                "{{Genotype|allele1=C|allele2=G|magnitude=2.5|summary=het variant}}",
            ),
            (
                "Rs4950928 (C;G)",
                "{{Genotype|allele1=C|allele2=G|magnitude=2.5|summary=het variant}}",
            ),
        ]
        db_path = _make_db(tmp_path, genotype_pages=genotype_pages)
        parse_raw_pages(db_path)
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            count = conn.execute("SELECT COUNT(*) FROM snpedia_genotypes").fetchone()[0]
            assert count == 1

    def test_backfill_dedupes_existing_rows(self, tmp_path):
        """Pre-existing duplicate rows are collapsed on next parser run."""
        db = tmp_path / "snpedia.sqlite"
        with contextlib.closing(sqlite3.connect(str(db))) as conn:
            conn.execute(
                "CREATE TABLE pages "
                "(title TEXT PRIMARY KEY, category TEXT, content TEXT, scraped_at TEXT)"
            )
            conn.execute(
                "CREATE TABLE snpedia_genotypes "
                "(rsid TEXT, allele1 TEXT, allele2 TEXT, magnitude REAL, "
                "repute TEXT, summary TEXT, gene TEXT, scraped_at TEXT)"
            )
            conn.execute(
                "INSERT INTO snpedia_genotypes VALUES "
                "('rs100','A','G',1.0,NULL,'test',NULL,'2026-01-01')"
            )
            conn.execute(
                "INSERT INTO snpedia_genotypes VALUES "
                "('rs100','A','G',1.0,NULL,'test',NULL,'2026-01-01')"
            )
            conn.commit()
            assert conn.execute("SELECT COUNT(*) FROM snpedia_genotypes").fetchone()[0] == 2
        with contextlib.closing(sqlite3.connect(str(db))) as conn:
            removed = _dedupe_existing(conn)
            assert removed == 1

    def test_backfill_dedupes_null_summary_rows(self, tmp_path):
        """NULL-summary duplicates collapse despite SQLite NULL≠NULL semantics."""
        db = tmp_path / "snpedia.sqlite"
        with contextlib.closing(sqlite3.connect(str(db))) as conn:
            conn.execute(
                "CREATE TABLE pages "
                "(title TEXT PRIMARY KEY, category TEXT, content TEXT, scraped_at TEXT)"
            )
            conn.execute(
                "CREATE TABLE snpedia_genotypes "
                "(rsid TEXT, allele1 TEXT, allele2 TEXT, magnitude REAL, "
                "repute TEXT, summary TEXT, gene TEXT, scraped_at TEXT)"
            )
            conn.execute(
                "INSERT INTO snpedia_genotypes VALUES "
                "('rs1050828','-','G',NULL,NULL,NULL,'G6PD','2026-01-01')"
            )
            conn.execute(
                "INSERT INTO snpedia_genotypes VALUES "
                "('rs1050828','-','G',NULL,NULL,NULL,'G6PD','2026-01-01')"
            )
            conn.commit()
            assert conn.execute("SELECT COUNT(*) FROM snpedia_genotypes").fetchone()[0] == 2
        with contextlib.closing(sqlite3.connect(str(db))) as conn:
            removed = _dedupe_existing(conn)
            assert removed == 1


class TestParserVersionStamp:
    """Parser version is stamped into database_versions.remote_signal."""

    def test_parse_stamps_version(self, tmp_path):
        genotype_pages = [
            ("Rs100(A;G)", "{{Genotype|allele1=A|allele2=G|summary=test}}"),
        ]
        db_path = _make_db(tmp_path, genotype_pages=genotype_pages)
        parse_raw_pages(db_path)
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            tag = conn.execute(
                "SELECT local_version_tag FROM database_versions WHERE name='snpedia'"
            ).fetchone()[0]
            assert tag == f"pv:{_PARSER_VERSION}"

    def test_parser_is_current_true(self, tmp_path):
        genotype_pages = [
            ("Rs100(A;G)", "{{Genotype|allele1=A|allele2=G|summary=test}}"),
        ]
        db_path = _make_db(tmp_path, genotype_pages=genotype_pages)
        parse_raw_pages(db_path)
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            assert parser_is_current(conn) is True

    def test_parser_is_current_false_old_version(self, tmp_path):
        genotype_pages = [
            ("Rs100(A;G)", "{{Genotype|allele1=A|allele2=G|summary=test}}"),
        ]
        db_path = _make_db(tmp_path, genotype_pages=genotype_pages)
        parse_raw_pages(db_path)
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            conn.execute(
                "UPDATE database_versions SET local_version_tag='pv:0' WHERE name='snpedia'"
            )
            conn.commit()
            assert parser_is_current(conn) is False

    def test_parser_is_current_false_no_tag(self, tmp_path):
        genotype_pages = [
            ("Rs100(A;G)", "{{Genotype|allele1=A|allele2=G|summary=test}}"),
        ]
        db_path = _make_db(tmp_path, genotype_pages=genotype_pages)
        parse_raw_pages(db_path)
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            conn.execute(
                "UPDATE database_versions SET local_version_tag=NULL WHERE name='snpedia'"
            )
            conn.commit()
            assert parser_is_current(conn) is False

    def test_parser_is_current_legacy_remote_signal_migrates(self, tmp_path):
        """Legacy cache with |pv:N in remote_signal self-heals without re-parse."""
        genotype_pages = [
            ("Rs100(A;G)", "{{Genotype|allele1=A|allele2=G|summary=test}}"),
        ]
        db_path = _make_db(tmp_path, genotype_pages=genotype_pages)
        parse_raw_pages(db_path)
        with contextlib.closing(sqlite3.connect(db_path)) as conn:
            conn.execute(
                "UPDATE database_versions "
                "SET remote_signal = ?, local_version_tag = NULL "
                "WHERE name = 'snpedia'",
                (f"etag:abc|pv:{_PARSER_VERSION}",),
            )
            conn.commit()
            assert parser_is_current(conn) is True
            row = conn.execute(
                "SELECT remote_signal, local_version_tag "
                "FROM database_versions WHERE name='snpedia'"
            ).fetchone()
            assert row[0] == "etag:abc"
            assert row[1] == f"pv:{_PARSER_VERSION}"

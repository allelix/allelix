# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
"""Shared helpers for parsers with CSV or concatenated-genotype formats.

Used by FTDNA, MyHeritage, and Living DNA parsers. Extracted here to avoid
duplicating the genotype-splitting and CSV-line-splitting logic across
structurally similar formats.
"""

from __future__ import annotations

import copy
import io
import logging
import zipfile
import zlib
from contextlib import contextmanager
from contextvars import ContextVar
from typing import IO, TYPE_CHECKING, TextIO

from allelix.models import NO_CALL_MARKER

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

logger = logging.getLogger(__name__)

# Notices about the input file as a whole (a possible cut, ignored trailing
# bytes), as opposed to per-line warnings about skipped rows. Callers that
# count skipped rows from the parser warnings leave this logger out.
INPUT_NOTICE_LOGGER = "allelix.parsers.input"
_notice_logger = logging.getLogger(INPUT_NOTICE_LOGGER)
# Each input notice's stable code, set on its log record as ``input_notice``
# so a report can list it among its warnings (ADR-0049 ``report_warnings``).
INPUT_NOTICE_CODES = {
    "trailing": "input_trailing_bytes",
    "bgzf-eof": "input_possibly_truncated",
}

GZIP_MAGIC = b"\x1f\x8b"
ZIP_MAGIC = b"PK\x03\x04"
# A zip with no entries at all is only its end-of-directory record.
_EMPTY_ZIP_MAGIC = b"PK\x05\x06"
_ZIP_MAGICS = (ZIP_MAGIC, _EMPTY_ZIP_MAGIC)
# A zip member is streamed, never extracted, so an oversized one costs time,
# not memory or disk. Past this many bytes read out of the zip the read stops
# (for a gzipped member, the gzip bytes); the size the zip directory declares
# is not trusted. The largest real genotype
# text seen is about 0.4 GB; this leaves room for a zipped whole-genome VCF.
ZIP_MAX_MEMBER_BYTES = 4 << 30
# Compression methods the standard library reads (stored, deflate, bzip2, lzma).
_ZIP_METHODS = {
    zipfile.ZIP_STORED,
    zipfile.ZIP_DEFLATED,
    zipfile.ZIP_BZIP2,
    zipfile.ZIP_LZMA,
}
_ZIP_METHOD_NAMES = {9: "Deflate64", 10: "PKWARE implode", 93: "Zstandard", 95: "XZ", 98: "PPMd"}
# The 28-byte empty block that ends every complete BGZF file (SAM/BAM
# specification, section 4.1.2). A file cut at a block boundary is still a
# valid gzip stream; only this block's absence shows the cut.
BGZF_EOF_BLOCK = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
_READ_CHUNK = 1 << 16

# Files already warned about in this process: a file is opened several times
# (detection, metadata, one or more parse passes) and each warning should
# appear once.
_warned: set[tuple[str, int, int, str]] = set()
# True while format detection sniffs a file: a small file can be read to its
# end there, before the caller has set up its warning output.
_quiet: ContextVar[bool] = ContextVar("allelix_parsers_quiet", default=False)


class DamagedFileError(ValueError):
    """A compressed genotype file ends early or holds corrupt data."""


class InputContainerError(ValueError):
    """A zip download Allelix cannot read as one genotype file.

    Several files in one zip, a password-protected member, a compression
    method Python cannot read, or a member past ``ZIP_MAX_MEMBER_BYTES``.
    The message says what to do instead.
    """


@contextmanager
def quiet_input_notices() -> Iterator[None]:
    """Hold back input notices (used by format detection)."""
    token = _quiet.set(True)
    try:
        yield
    finally:
        _quiet.reset(token)


def _warn_once(file_path: Path, kind: str, message: str) -> None:
    if _quiet.get():
        return
    try:
        stat = file_path.stat()
        key = (str(file_path.resolve()), stat.st_size, stat.st_mtime_ns, kind)
    except OSError:
        key = (str(file_path), -1, -1, kind)
    if key in _warned or not _notice_logger.isEnabledFor(logging.WARNING):
        return
    record = _notice_logger.makeRecord(
        _notice_logger.name,
        logging.WARNING,
        __file__,
        0,
        message,
        (),
        None,
        extra={"input_notice": INPUT_NOTICE_CODES[kind]},
    )
    # A read whose records a filter drops (the build pre-scan) does not use up
    # the one notice; the next read shows it.
    if not _notice_logger.filter(record):
        return
    _warned.add(key)
    _notice_logger.handle(record)


def _damaged(file_path: Path, detail: str) -> DamagedFileError:
    return DamagedFileError(
        f"{file_path.name} looks truncated or damaged ({detail}). "
        "Download or copy it again; the records after the damage cannot be read."
    )


def _is_bgzf_header(header: bytes) -> bool:
    """True when ``header`` starts a BGZF block (gzip with the ``BC`` extra subfield)."""
    return len(header) >= 14 and header[3] & 0x04 != 0 and header[12:14] == b"BC"


class _GzipMembers(io.RawIOBase):
    """Decompress every gzip member of a file, naming damage instead of raising zlib errors.

    ``gzip.GzipFile`` raises ``EOFError`` for a stream cut inside a member
    and ``BadGzipFile`` for bytes after the last member that are not gzip
    data. Here a cut or corrupt member raises ``DamagedFileError``, and
    non-gzip bytes after at least one complete member are ignored with a
    warning, as the ``gzip`` command does. zlib checks each member's CRC
    and length trailer. Reaching the end of a BGZF file runs the
    end-of-file block check; detection never reads that far, so the
    warning comes from a parse.
    """

    def __init__(self, file_path: Path, raw: IO[bytes] | None = None) -> None:
        super().__init__()
        self._path = file_path
        # ``raw`` is a stream inside a container (a zip member): it is read
        # from its start, and the BGZF end-block check, which seeks in the
        # file itself, is skipped.
        self._in_container = raw is not None
        self._raw = raw if raw is not None else open(file_path, "rb")  # noqa: SIM115
        self._in = b""  # compressed bytes read but not yet decompressed
        self._out = b""  # decompressed bytes; those before _pos are already returned
        self._pos = 0
        self._decomp: zlib._Decompress | None = None
        self._members = 0
        self._done = False

    def readable(self) -> bool:
        return True

    def close(self) -> None:
        if not self.closed:
            self._raw.close()
        super().close()

    def readinto(self, buffer: memoryview) -> int:  # type: ignore[override]
        while self._pos >= len(self._out) and not self._done:
            self._fill()
        size = min(len(buffer), len(self._out) - self._pos)
        buffer[:size] = self._out[self._pos : self._pos + size]
        self._pos += size
        return size

    def _fill(self) -> None:
        if self._decomp is None:
            self._start_member()
            return
        data = self._in or self._raw.read(_READ_CHUNK)
        self._in = b""
        if not data:
            raise _damaged(self._path, "the compressed data ends in the middle of a block")
        try:
            self._out = self._decomp.decompress(data)
            self._pos = 0
        except zlib.error as exc:
            raise _damaged(self._path, f"corrupt compressed data: {exc}") from exc
        if self._decomp.eof:
            self._in = self._decomp.unused_data
            self._decomp = None
            self._members += 1

    def _start_member(self) -> None:
        while len(self._in) < 2:
            more = self._raw.read(_READ_CHUNK)
            if not more:
                break
            self._in += more
        if not self._in:
            self._done = True
            if self._members and not self._in_container:
                _check_bgzf_end(self._path)
            return
        if self._in[:2] == GZIP_MAGIC:
            self._decomp = zlib.decompressobj(wbits=31)
            return
        if self._members == 0:
            raise _damaged(self._path, "not gzip data")
        trailing = len(self._in) + len(self._raw.read())
        self._in = b""
        self._done = True
        _warn_once(
            self._path,
            "trailing",
            f"{self._path.name}: ignored {trailing} bytes after the end of the "
            "compressed data (they are not gzip data); every complete record was read",
        )


def _check_bgzf_end(file_path: Path) -> None:
    """Warn when a BGZF file lacks its end-of-file block (a cut at a block boundary)."""
    with open(file_path, "rb") as raw:
        header = raw.read(14)
        if not _is_bgzf_header(header):
            return
        raw.seek(0, io.SEEK_END)
        size = raw.tell()
        if size >= len(BGZF_EOF_BLOCK):
            raw.seek(size - len(BGZF_EOF_BLOCK))
            if raw.read() == BGZF_EOF_BLOCK:
                return
    _warn_once(
        file_path,
        "bgzf-eof",
        f"{file_path.name} does not end with the BGZF end-of-file block, so it may be "
        "truncated: an interrupted download or copy loses the records after the cut "
        "without any other sign",
    )


def _is_debris(info: zipfile.ZipInfo) -> bool:
    """Directories and macOS archive debris are not members (``__MACOSX/``, ``.DS_Store``)."""
    name = info.filename
    return info.is_dir() or name.startswith("__MACOSX/") or name.rsplit("/", 1)[-1] == ".DS_Store"


def _zip_member(file_path: Path, archive: zipfile.ZipFile) -> zipfile.ZipInfo:
    """The one genotype member of ``archive``, or the error that says why not."""
    members = [info for info in archive.infolist() if not _is_debris(info)]
    if not members:
        raise InputContainerError(f"{file_path.name} is a zip file with no file in it.")
    if len(members) > 1:
        names = ", ".join(info.filename for info in members[:5])
        more = f" and {len(members) - 5} more" if len(members) > 5 else ""
        raise InputContainerError(
            f"{file_path.name} is a zip file holding {len(members)} files ({names}{more}). "
            "Allelix reads a zip with one genotype file in it; unzip this one and pass "
            "the genotype file itself."
        )
    (member,) = members
    if member.flag_bits & 0x1:
        raise InputContainerError(
            f"{file_path.name} is a password-protected zip file. Unzip it first and pass "
            "the genotype file itself."
        )
    if member.compress_type not in _ZIP_METHODS:
        method = _ZIP_METHOD_NAMES.get(member.compress_type, f"method {member.compress_type}")
        raise InputContainerError(
            f"{file_path.name} is compressed with {method}, which Allelix cannot read. "
            "Unzip it first and pass the genotype file itself."
        )
    if member.file_size > ZIP_MAX_MEMBER_BYTES:
        raise _too_big(file_path)
    return member


def _too_big(file_path: Path) -> InputContainerError:
    return InputContainerError(
        f"{file_path.name} holds a file over {ZIP_MAX_MEMBER_BYTES / (1 << 30):g} GiB "
        "unzipped, more than Allelix reads from a zip. Unzip it, or gzip it instead, "
        "and pass that file."
    )


def check_container(file_path: Path) -> None:
    """Raise now when ``file_path`` is a zip Allelix cannot read as one genotype file.

    Reads the zip directory only. Format detection calls it first, so a
    damaged zip or one with several files gets its own error instead of
    "No parser recognized" (each parser's ``can_parse`` treats an unreadable
    file as not its format). Anything that is not a zip passes.
    """
    try:
        with open(file_path, "rb") as raw:
            if raw.read(4) not in _ZIP_MAGICS:
                return
    except OSError:
        return
    try:
        with zipfile.ZipFile(file_path) as archive:
            _zip_member(file_path, archive)
    except zipfile.BadZipFile as exc:
        raise _damaged(file_path, f"the zip directory is missing or corrupt: {exc}") from exc


class _ZipMember(io.RawIOBase):
    """Stream the single member of a zip, naming damage and stopping at the size cap.

    zipfile stops reading a member at the size its entry declares and drops
    the rest without a word when the CRC was written for the shorter data.
    The member is opened with that size lifted, so every byte of the
    compressed data is read (and the CRC still checked), and a member whose
    data is longer or shorter than its entry declares is damage.
    """

    def __init__(self, file_path: Path) -> None:
        super().__init__()
        self._path = file_path
        try:
            self._archive = zipfile.ZipFile(file_path)
        except zipfile.BadZipFile as exc:
            raise _damaged(file_path, f"the zip directory is missing or corrupt: {exc}") from exc
        try:
            info = _zip_member(file_path, self._archive)
            self._declared = info.file_size
            unbounded = copy.copy(info)
            # Far past the cap, so the cap and size checks below run first.
            unbounded.file_size = 1 << 62
            self._member = self._archive.open(unbounded)
        except BaseException:
            self._archive.close()
            raise
        self._read = 0

    def readable(self) -> bool:
        return True

    def close(self) -> None:
        if not self.closed:
            self._member.close()
            self._archive.close()
        super().close()

    def readinto(self, buffer: memoryview) -> int:  # type: ignore[override]
        try:
            data = self._member.read(len(buffer))
        except (zipfile.BadZipFile, zlib.error, EOFError, ValueError) as exc:
            # A CRC failure, corrupt deflate data or a member cut short.
            raise _damaged(self._path, f"corrupt zip member: {exc}") from exc
        self._read += len(data)
        if self._read > self._declared:
            if self._read > ZIP_MAX_MEMBER_BYTES:
                raise _too_big(self._path)
            raise _damaged(
                self._path,
                f"the zip member holds more than the {self._declared:,} bytes its entry declares",
            )
        if not data and len(buffer) and self._read < self._declared:
            raise _damaged(
                self._path,
                f"the zip member ends after {self._read:,} of the {self._declared:,} bytes "
                "its entry declares",
            )
        buffer[: len(data)] = data
        return len(data)


def open_genotype_text(file_path: Path, errors: str = "strict") -> TextIO:
    """Open a genotype file in text mode, gzipped or not, without a byte order mark.

    Sniffs the gzip magic bytes (``1f 8b``) rather than trusting the
    extension: FTDNA delivers ``.csv.gz`` and users rename files. Decodes
    ``utf-8-sig``, which drops a leading UTF-8 byte order mark and reads
    unmarked input exactly as ``utf-8`` does, so a signature or header
    check on the first line sees the same text either way.

    Gzip input is read by ``_GzipMembers``: a cut or corrupt stream raises
    ``DamagedFileError`` while reading, and a BGZF file read to its end
    without finding its end-of-file block is warned about once.

    A zip download (``PK`` magic bytes) is read from its one member,
    streamed and never extracted; the member may itself be gzipped. A zip
    with several files, a password or an unreadable compression method
    raises ``InputContainerError``; a damaged one ``DamagedFileError``.
    """
    with open(file_path, "rb") as raw:
        magic = raw.read(4)
    if magic in _ZIP_MAGICS:
        member = io.BufferedReader(_ZipMember(file_path), _READ_CHUNK)
        try:
            head = member.peek(4)[:4]
            if head in _ZIP_MAGICS:
                raise InputContainerError(
                    f"{file_path.name} is a zip file holding another zip file. Unzip it "
                    "until you reach the genotype file, and pass that file."
                )
        except BaseException:
            member.close()
            raise
        gzipped = head[:2] == GZIP_MAGIC
        stream: IO[bytes] = (
            io.BufferedReader(_GzipMembers(file_path, raw=member), _READ_CHUNK)
            if gzipped
            else member
        )
        return io.TextIOWrapper(stream, encoding="utf-8-sig", errors=errors)
    if magic[:2] == GZIP_MAGIC:
        return io.TextIOWrapper(
            io.BufferedReader(_GzipMembers(file_path), _READ_CHUNK),
            encoding="utf-8-sig",
            errors=errors,
        )
    return open(file_path, encoding="utf-8-sig", errors=errors)


def first_nonblank_line(fh: TextIO, limit: int) -> str:
    """The first line of ``fh`` with any non-whitespace in it, or "" (within ``limit`` lines)."""
    for _ in range(limit):
        line = fh.readline()
        if not line or line.strip():
            return line
    return ""


def split_csv_line(line: str, columns: int | None = None) -> list[str]:
    """Split a comma-delimited line and strip surrounding quotes from each field.

    Implementation is ``line.split(",")`` followed by a per-field
    ``strip().strip('"')``. This is NOT a real CSV parser: a quoted field
    containing a literal comma yields the wrong column count and is
    silently dropped by callers' ``len(parts) != EXPECTED_COLUMNS``
    guard.

    Adequate for FTDNA / MyHeritage / Living DNA because every value in
    those exports is either an rsID, chromosome identifier, integer
    position, or concatenated genotype string — none of which contain
    commas. If a future format ever ships embedded commas in quoted
    fields, swap to ``csv.reader`` rather than relying on this helper.

    Strips both surrounding double quotes (``"rs1"``) and the
    double-double-quote variant some MyHeritage exports produce
    (``""rs1""``) — the latter via two iterations of the trailing
    ``strip('"')``.

    ``columns``: the format's column count. Trailing empty fields beyond
    it are dropped: a spreadsheet re-save pads every row to the widest
    column count (``"rs1","1","82154","AA",,,``), which would otherwise
    read as seven columns and be skipped. Fields up to ``columns`` are
    kept even when empty, so an empty genotype stays a no-call.
    """
    fields = [field.strip().strip('"') for field in line.split(",")]
    return fields if columns is None else strip_trailing_empty(fields, keep=columns)


def strip_trailing_empty(fields: list[str], keep: int = 0) -> list[str]:
    """Drop trailing empty fields (spreadsheet padding such as 'a,b,,,' -> ['a','b']).

    Never shortens the list below ``keep`` fields.
    """
    while len(fields) > keep and fields[-1] == "":
        fields.pop()
    return fields


def is_padding_only(line: str) -> bool:
    """True for a row of nothing but delimiters (``,,,``): an emptied spreadsheet row."""
    return not line.replace(",", "").replace("\t", "").replace('"', "").strip()


# Array files write a no-call as ``--``, ``00`` or ``NN`` (``-``, ``0`` or
# ``N`` per allele column). None of them is a base: each reads as a no-call.
NO_CALL_CODES = frozenset({NO_CALL_MARKER, "0", "N"})


def no_call_allele(allele: str) -> str:
    """``allele``, or ``NO_CALL_MARKER`` when it is a no-call code."""
    return NO_CALL_MARKER if allele.strip().upper() in NO_CALL_CODES else allele


def split_genotype(genotype: str) -> tuple[str, str]:
    """Split a concatenated genotype field into two alleles.

    ``"AG"`` -> ``("A", "G")``; ``"--"``, ``"00"`` and ``"NN"`` ->
    ``("-", "-")``; ``"A"`` -> ``("A", "A")``: a one-letter call is
    haploid, and ``is_haploid_call`` says so (the Variant carries
    ``haploid``). A one-letter ``-``, ``0`` or ``N`` is a no-call.
    """
    if len(genotype) == 2:
        return no_call_allele(genotype[0]), no_call_allele(genotype[1])
    if len(genotype) == 1:
        allele = no_call_allele(genotype)
        return allele, allele
    logger.warning("Unexpected genotype format %r — treating as no-call", genotype)
    return NO_CALL_MARKER, NO_CALL_MARKER


def is_haploid_call(genotype: str) -> bool:
    """True for a one-letter call: the file declares one copy (ADR-0050)."""
    return len(genotype) == 1 and genotype.upper() not in NO_CALL_CODES


def normalize_chromosome(chrom: str) -> str:
    """Normalize a chromosome identifier to Allelix's canonical bare form.

    VCF/UCSC sources commonly prefix chromosomes with ``chr`` (``chr1``,
    ``chrX``, ``chrM``). Allelix's internal representation uses the bare
    form (``1``, ``X``, ``MT``). Returns the bare form regardless of input.

    Mappings:
      - ``chr1`` / ``CHR1`` / ``1`` -> ``1``
      - ``chrX`` / ``chrx`` / ``X`` / ``x`` -> ``X`` (same for Y)
      - ``chrM`` / ``chrm`` / ``chrMT`` / ``chrmt`` / ``M`` / ``MT`` /
        ``m`` / ``mt`` -> ``MT``
      - Unknown contigs (``chrUn_*``, ``GL00*``, etc.) pass through unchanged
        after the optional ``chr`` strip — downstream code decides whether
        to skip them. Case is preserved for unknown contigs.

    Mitochondrial nomenclature varies: GRCh37 chrMT, GRCh38 chrM, some
    VCFs ship bare M. All collapse to MT for Allelix.

    Standard chromosome names (1-22, X, Y, M, MT) are normalized to
    upper case so a pipeline that lowercases the CHROM column doesn't
    silently miss every annotation database lookup. Unknown contigs
    preserve case because their identifiers are opaque and the
    canonical-vs-non-canonical distinction is what matters.
    """
    if not chrom:
        return chrom
    # Strip case-insensitive 'chr' prefix
    stripped = chrom[3:] if chrom[:3].lower() == "chr" else chrom
    upper = stripped.upper()
    # Collapse mitochondrial variants to canonical MT
    if upper in ("M", "MT"):
        return "MT"
    # Standard sex chromosomes normalize to uppercase
    if upper in ("X", "Y"):
        return upper
    # Autosomes (digits) are unchanged by upper(); return as-is
    if stripped.isdigit():
        return stripped
    # Unknown contigs (GL00*, chrUn_*, alts) preserve case
    return stripped

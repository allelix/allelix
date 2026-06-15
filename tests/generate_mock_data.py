# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Allelix
#!/usr/bin/env python3
"""Generate synthetic MyHappyGenes genotype files for Allelix testing.

This script produces deterministic mock genotype data that:
- Matches the real MyHappyGenes/Tempus file format exactly
- Contains known clinically significant SNPs with hardcoded genotypes
- Covers all chromosomes including X, Y, MT
- Includes deliberate no-calls, heterozygous, and homozygous variants
- Uses fake rsIDs (rs900XXXXXX range) for filler to avoid ClinVar collisions
- Is fully reproducible (deterministic, no randomness in committed fixture)

Usage:
    python generate_mock_data.py                    # 2000-line test fixture
    python generate_mock_data.py --large             # 50000-line perf fixture
    python generate_mock_data.py --output my_file.txt
"""

import argparse
import hashlib
import shutil
from pathlib import Path

# ============================================================================
# KNOWN SNPS: Hardcoded genotypes for deterministic test assertions
# ============================================================================
# ADR-0021: each SNP carries a (chrom, GRCh37 position, GRCh38 position)
# triple so the generator can emit either build cleanly. The build-
# detection feature uses the GRCh37/GRCh38 position differences for the
# SNPs that are also in `allelix.utils.build_detect.KNOWN_SNP_POSITIONS`.
# For SNPs outside that table, the positions don't strictly need to be
# real per-build — but we keep them honest where authoritative data is
# readily available.
#
# Each entry: (rsid, chrom, position_grch37, position_grch38, allele1,
# allele2, note).
#
# When `position_grch37 == position_grch38`, the position happens to be
# build-invariant for that SNP (rare for autosomes; common for tiny
# chromosomes or pre-2009 coordinates).

KNOWN_SNPS = [
    # ── Methylation panel ────────────────────────────────────────────
    # MTHFR C677T: heterozygous (one risk allele, reduced function)
    ("rs1801133", "1", 11856378, 11796321, "G", "A", "MTHFR C677T het"),
    # MTHFR A1298C: wild type (no risk)
    ("rs1801131", "1", 11854476, 11794419, "T", "T", "MTHFR A1298C wt"),
    # MTR A2756G: wild type
    ("rs1805087", "1", 236885200, 236721900, "A", "A", "MTR wt"),
    # MTRR A66G: homozygous variant (both risk alleles)
    ("rs1801394", "5", 7870860, 7870973, "G", "G", "MTRR hom variant"),
    # COMT V158M: homozygous slow (Met/Met)
    ("rs4680", "22", 19951271, 19963748, "A", "A", "COMT slow"),
    # CBS C699T: heterozygous
    ("rs234706", "21", 43065240, 41645140, "G", "A", "CBS het"),
    # ── Pharmacogenomics ─────────────────────────────────────────────
    # CYP2D6: heterozygous (intermediate metabolizer)
    ("rs1065852", "22", 42526694, 42130692, "G", "A", "CYP2D6 het"),
    # CYP2C19*2: wild type (normal metabolizer)
    ("rs4244285", "10", 96541616, 94781859, "G", "G", "CYP2C19 wt"),
    # CYP2C9*2: heterozygous (reduced warfarin metabolism)
    ("rs1799853", "10", 96702047, 94942290, "C", "T", "CYP2C9*2 het"),
    # SLCO1B1: heterozygous (statin myopathy risk)
    ("rs4149056", "12", 21331549, 21178615, "T", "C", "SLCO1B1 het"),
    # ── ClinVar pathogenic (should trigger alerts) ───────────────────
    # BRCA1: heterozygous carrier (pathogenic)
    ("rs80357906", "17", 41209080, 43057063, "G", "A", "BRCA1 carrier"),
    # TP53: wild type (should NOT be flagged)
    ("rs121918506", "17", 7577538, 7674222, "G", "G", "TP53 wt"),
    # ── Carrier status (recessive, one copy) ─────────────────────────
    # CFTR delta F508: arrays can't reliably call this indel — it's a
    # 3-base deletion (ΔF508), not a SNV. Real MyHappyGenes files report
    # a no-call here, NOT a multi-base genotype. Earlier versions of this
    # generator put "CTT"/"C" here, which violated the MHG format
    # invariant (single-base alleles only) and silently hid the
    # indel-anchor bug ClinVar surfaced in v0.4.2. See ADR-0015.
    ("rs113993960", "7", 117199644, 117559590, "-", "-", "CFTR delta F508 (no-call on array)"),
    # NIPA1 — ADR-0021 strand-inversion regression case. User has G/G.
    # GRCh37: REF=C ALT=G (Pathogenic). GRCh38: REF=G ALT=A (Pathogenic).
    # On GRCh38 the user is homozygous reference (NO annotation). On
    # GRCh37 the user's G matches ALT=G (incorrect pathogenic emission
    # if cross-build comparison happens). The auto-detector picks
    # GRCh38 from the file's positions; the per-build cache dispatch
    # ensures the right ClinVar entry is queried.
    ("rs104894490", "15", 23060816, 22812251, "G", "G", "NIPA1 strand-inversion regression"),
    # ── No-call examples ─────────────────────────────────────────────
    # Deliberate no-calls on known SNPs to test handling
    ("rs12248560", "10", 96521657, 94761900, "-", "-", "CYP2C19*17 no-call"),
    ("rs3918290", "1", 97915614, 97450058, "-", "-", "DPYD no-call"),
]


def _build_header_label(build_in_header: str) -> str:
    """Map a target build to the literal MHG header phrasing."""
    return "37.1" if build_in_header == "GRCh37" else "38"


# Chromosome sort order matching real MyHappyGenes files
CHROM_ORDER = {str(i): i for i in range(1, 23)}
CHROM_ORDER.update({"X": 23, "Y": 24, "MT": 25})

# Approximate chromosome sizes (hg19/GRCh37) for realistic position generation
CHROM_SIZES = {
    "1": 249250621,
    "2": 243199373,
    "3": 198022430,
    "4": 191154276,
    "5": 180915260,
    "6": 171115067,
    "7": 159138663,
    "8": 146364022,
    "9": 141213431,
    "10": 135534747,
    "11": 135006516,
    "12": 133851895,
    "13": 115169878,
    "14": 107349540,
    "15": 102531392,
    "16": 90354753,
    "17": 81195210,
    "18": 78077248,
    "19": 59128983,
    "20": 63025520,
    "21": 48129895,
    "22": 51304566,
    "X": 155270560,
    "Y": 59373566,
    "MT": 16569,
}

ALLELES = ["A", "T", "G", "C"]


def deterministic_genotype(seed_value: int, het_rate: float = 0.3) -> tuple[str, str]:
    """Generate a deterministic genotype from a seed value.

    Uses the seed to pick alleles consistently. No randomness.
    """
    a1_idx = seed_value % 4
    a1 = ALLELES[a1_idx]

    # Use a different bit of the seed for het/hom decision
    is_het = (seed_value % 100) < (het_rate * 100)
    if is_het:
        a2_idx = (seed_value // 4) % 3
        remaining = [a for a in ALLELES if a != a1]
        a2 = remaining[a2_idx]
    else:
        a2 = a1

    return a1, a2


SnpRow = tuple[str, str, int, str, str]


def generate_filler_snps(count: int, no_call_rate: float = 0.05) -> list[SnpRow]:
    """Generate deterministic filler SNPs across all chromosomes.

    Uses rs900XXXXXXX range to avoid collisions with real rsIDs.
    Real rsIDs go up to ~900M as of 2025 but the 900_000_000+ range
    is sparse enough to be safe for synthetic data.
    """
    filler = []
    seen_rsids = set()

    for i in range(count):
        # Deterministic chromosome distribution (weighted toward larger chroms)
        chrom_idx = i % 25
        if chrom_idx < 22:
            chrom = str(chrom_idx + 1)
        elif chrom_idx == 22:
            chrom = "X"
        elif chrom_idx == 23:
            chrom = "Y"
        else:
            chrom = "MT"

        # Deterministic position spread within realistic chromosome bounds
        pos_hash = int(hashlib.md5(f"pos_{i}".encode()).hexdigest()[:8], 16)
        max_pos = CHROM_SIZES.get(chrom, 248000000)
        position = 1000 + (pos_hash % (max_pos - 1000))

        # Fake rsID in the 900M+ range
        rsid = f"rs{900_000_001 + i}"
        if rsid in seen_rsids:
            continue
        seen_rsids.add(rsid)

        # Deterministic no-calls at fixed interval
        if i % int(1 / no_call_rate) == 7:
            a1, a2 = "-", "-"
        else:
            a1, a2 = deterministic_genotype(pos_hash)

        filler.append((rsid, chrom, position, a1, a2))

    return filler


def write_mhg_file(
    filepath: str,
    snps: list[SnpRow],
    *,
    header_build: str = "GRCh37",
) -> None:
    """Write SNPs in exact MyHappyGenes/Tempus format.

    `header_build` controls only the build-label comment in the header,
    NOT the positions in the rows. To produce the deliberate-mislabel
    fixture (header claims GRCh37, positions are GRCh38), pass
    `header_build="GRCh37"` and supply GRCh38 positions in `snps`.
    """
    label = _build_header_label(header_build)
    # Sort by chromosome then position
    sorted_snps = sorted(snps, key=lambda x: (CHROM_ORDER.get(x[1], 99), x[2]))

    with open(filepath, "w") as f:
        f.write("# MyHappyGenes [TEMPUS]\n")
        f.write("# This file was generated by MyHappyGenes, Inc.\n")
        f.write("# Below is a text version of your DNA from Tempus.\n")
        f.write("# THIS INFORMATION IS FOR YOUR PERSONAL USE AND IS INTENDED FOR RESEARCH ONLY\n")
        f.write(
            "# IT IS NOT INTENDED FOR MEDICAL, DIAGNOSTIC, OR HEALTH PURPOSES. THE EXPORTED DATA\n"
        )
        f.write(
            "# IS SUBJECT TO THE MyHappyGenes TERMS AND CONDITIONS, BUT PLEASE BE AWARE THAT\n"
        )
        f.write("# THE DOWNLOADED DATA WILL NO LONGER BE PROTECTED BY OUR SECURITY MEASURES.\n")
        f.write("# WHEN YOU DOWNLOAD YOUR RAW DATA, YOU ASSUME ALL RISK OF STORING,\n")
        f.write(
            "# SECURING AND PROTECTING YOUR DATA. FOR MORE INFORMATION, SEE MyHappyGenes FAQs.\n"
        )
        f.write("# Genetic data is provided below as five TAB delimited columns. Each line\n")
        f.write(
            "# corresponds to an SNP. Column one provides the SNP"
            " identifier (rsid where possible).\n"
        )
        f.write(
            "# Columns two and three contain the chromosome and basepair position of the SNP\n"
        )
        f.write(f"# using human reference build {label} coordinates. Columns four and five\n")
        f.write(
            "# contain the two alleles observed at this SNP (genotype). The genotype is reported\n"
        )
        f.write("# on the forward (+) strand with respect to the human reference.\n")
        f.write("# Sample ID\tMHG000001\n")
        f.write("SNP Name\tChr\tPosition\tAllele1 - Forward\tAllele2 - Forward\n")

        for rsid, chrom, pos, a1, a2 in sorted_snps:
            f.write(f"{rsid}\t{chrom}\t{pos}\t{a1}\t{a2}\n")

    total = len(sorted_snps)
    no_calls = sum(1 for _, _, _, a1, a2 in sorted_snps if a1 == "-" or a2 == "-")
    het = sum(1 for _, _, _, a1, a2 in sorted_snps if a1 != a2 and a1 != "-" and a2 != "-")
    known = len(KNOWN_SNPS)

    print(f"Generated {filepath}")
    print(f"  Total SNPs:    {total}")
    print(f"  Known SNPs:    {known}")
    print(f"  Filler SNPs:   {total - known}")
    print(f"  No-calls:      {no_calls} ({no_calls / total * 100:.1f}%)")
    print(f"  Heterozygous:  {het} ({het / total * 100:.1f}%)")


def generate_mock_gnomad(fixtures_dir: Path) -> None:
    """Build a tiny gzipped gnomAD SQLite for test_cli.py db update tests."""
    import contextlib
    import gzip
    import sqlite3

    from allelix.databases.schema import GNOMAD_SCHEMA

    raw_db = fixtures_dir / "mock_gnomad.sqlite"
    gz_path = fixtures_dir / "mock_gnomad.sqlite.gz"

    with contextlib.closing(sqlite3.connect(raw_db)) as conn:
        for stmt in GNOMAD_SCHEMA.split(";"):
            stmt = stmt.strip()
            if stmt:
                conn.execute(stmt)
        conn.executemany(
            "INSERT INTO gnomad_frequencies"
            " (chrom, pos, ref, alt, rsid, af) VALUES (?, ?, ?, ?, ?, ?)",
            [
                ("1", 11796321, "G", "A", "rs1801133", 0.35),
                ("22", 19963748, "G", "A", "rs4680", 0.50),
                ("19", 44908684, "T", "C", "rs429358", 0.15),
            ],
        )
        conn.execute(
            "INSERT INTO database_versions"
            " (name, source_url, version, downloaded_at, record_count)"
            " VALUES (?, ?, ?, ?, ?)",
            ("gnomad", "test://mock", "4.1", "2026-01-01T00:00:00Z", 3),
        )
        conn.commit()

    with raw_db.open("rb") as f_in, gzip.open(gz_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    raw_db.unlink()
    print(f"Generated {gz_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate synthetic MyHappyGenes genotype files for Allelix testing."
    )
    parser.add_argument(
        "--output",
        "-o",
        default="tests/fixtures/mock_myhappygenes.txt",
        help="Output file path (default: tests/fixtures/mock_myhappygenes.txt)",
    )
    parser.add_argument(
        "--large",
        action="store_true",
        help="Generate large fixture (50K lines) for performance testing",
    )
    parser.add_argument(
        "--build",
        choices=("grch37", "grch38"),
        default="grch38",
        help=(
            "Which build's positions to emit in the SNP rows. Default 'grch38' "
            "matches real-world MyHappyGenes files (which ship GRCh38 data "
            "despite their header — see ADR-0021)."
        ),
    )
    parser.add_argument(
        "--header-build",
        choices=("grch37", "grch38"),
        default=None,
        help=(
            "Build label to write into the file header. When omitted, matches "
            "--build. To replicate the real-world MyHappyGenes mislabel "
            "(header claims '37.1' but positions are GRCh38), pass "
            "--build grch38 --header-build grch37."
        ),
    )
    args = parser.parse_args()

    if args.large:
        output = args.output.replace(".txt", "_large.txt")
        filler_count = 50_000
    else:
        output = args.output
        filler_count = 2_000

    Path(output).parent.mkdir(parents=True, exist_ok=True)

    target_build = "GRCh37" if args.build == "grch37" else "GRCh38"
    header_build = (
        target_build
        if args.header_build is None
        else ("GRCh37" if args.header_build == "grch37" else "GRCh38")
    )

    # KNOWN_SNPS is (rsid, chrom, pos_grch37, pos_grch38, a1, a2, note).
    position_idx = 2 if target_build == "GRCh37" else 3
    known = [
        (rsid, chrom, entry[position_idx], a1, a2)
        for entry in KNOWN_SNPS
        for rsid, chrom, _p37, _p38, a1, a2, _note in (entry,)
    ]
    filler = generate_filler_snps(filler_count)
    all_snps = known + filler

    write_mhg_file(output, all_snps, header_build=header_build)

    generate_mock_gnomad(Path("tests/fixtures"))


if __name__ == "__main__":
    main()

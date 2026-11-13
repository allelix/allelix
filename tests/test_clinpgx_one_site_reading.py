"""A resolved reverse site supplies one genotype reading, never two."""

import contextlib
import itertools
import json
import sqlite3

import pytest
from click.testing import CliRunner
from test_annotation_parity_and_named_disclosures import _clinvar, _identity, _pharm

from allelix.annotators import get_annotators
from allelix.cli import main
from allelix.models import Variant
from allelix.utils.allele import complement

RS = "rs900111001"
BUILDS = ("GRCh37", "GRCh38")
POSITIONS = {"GRCh37": 1000000, "GRCh38": 5000000}
LENGTHS = {"GRCh37": 249250621, "GRCh38": 248956422}
SITES = ("AC", "AG", "CT", "GT", "ACG", "ACT", "AGT", "CGT")
CASES = [
    (site, "".join(pair))
    for site in SITES
    for pair in itertools.combinations_with_replacement(site, 2)
]


def _pair(text):
    return "".join(sorted(text))


def _data(path, site, *, clinvar=False, reading="complement", reverse=False, nonfinding=None):
    path.mkdir(exist_ok=True)
    genotypes = [
        _pair(complement("".join(g))) if reading == "complement" else "".join(g)
        for g in itertools.combinations_with_replacement(site, 2)
    ]
    _pharm(path, RS, "SYNPGX", genotypes[::-1] if reverse else genotypes, level="1A")
    with contextlib.closing(sqlite3.connect(path / "pharmgkb.sqlite")) as conn, conn:
        for genotype in genotypes:
            conn.execute(
                "UPDATE pharmgkb_annotations SET annotation_text=? WHERE genotype=?",
                (f"row {genotype}", genotype),
            )
        if nonfinding is not None:
            conn.execute(
                "UPDATE pharmgkb_annotations SET is_nonfinding=1 WHERE genotype=?", (nonfinding,)
            )
    _identity(
        path,
        [
            (
                RS,
                "1",
                POSITIONS["GRCh37"],
                site[0],
                a,
                POSITIONS["GRCh38"],
                site[0],
                a,
                int(site in {"AT", "CG"}),
            )
            for a in site[1:]
        ],
    )
    if clinvar:
        for build in BUILDS:
            _clinvar(
                path,
                build,
                [
                    (
                        RS,
                        "1",
                        POSITIONS[build],
                        site[0],
                        site[1],
                        "SYNPGX",
                        "Pathogenic",
                        "Pathogenic",
                        1,
                        1,
                    )
                ],
            )
    names = (
        "clinvar",
        "pharmgkb",
        "snpedia",
        "gwas",
        "gnomad",
        "cadd",
        "alphamissense",
        "rsid_identity",
    )
    enabled = {"pharmgkb", "rsid_identity"} | ({"clinvar"} if clinvar else set())
    (path / "config.toml").write_text(
        "[sources]\n" + "".join(f"{n}={str(n in enabled).lower()}\n" for n in names)
    )
    return path


@pytest.fixture(scope="module")
def cache(tmp_path_factory):
    """Reuse read-only synthetic databases across the genotype matrix."""
    root = tmp_path_factory.mktemp("site-reading")
    paths = {}

    def get(site, **options):
        key = (site, tuple(sorted(options.items())))
        if key not in paths:
            paths[key] = _data(root / f"cache-{len(paths)}", site, **options)
        return paths[key]

    return get


def _call(site, genotype, build):
    return Variant(RS, "1", POSITIONS[build], genotype[0], genotype[1], build=build, ref=site[0])


def _find(data, call, path):
    with contextlib.ExitStack() as stack:
        registry = get_annotators(data)
        for annotator in registry:
            stack.enter_context(annotator)
        pgx = next(a for a in registry if a.name == "pharmgkb")
        assert pgx.is_ready(), "Synthetic cache must be readable before the matching assertion"
        return pgx.annotate(call) if path == "single" else list(pgx.batch_annotate([call]))


def _labels(rows):
    return sorted(row.description.split(" — ")[-1] for row in rows)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("path", ["single", "batch"])
@pytest.mark.parametrize("clinvar", [False, True])
@pytest.mark.parametrize("site,genotype", CASES)
def test_every_genotype_has_only_its_reverse_site_row(cache, build, path, clinvar, site, genotype):
    data = cache(site, clinvar=clinvar)
    call = _call(site, genotype, build)
    rows = _find(data, call, path)
    expected = (
        [] if clinvar and genotype == site[0] * 2 else [f"row {_pair(complement(genotype))}"]
    )
    assert _labels(rows) == expected, (
        "Definitive reverse evidence cannot also admit the exact-spelling genotype"
    )
    assert all(
        r.genotype_match == genotype and r.carried_alleles == (genotype[0], genotype[1])
        for r in rows
    )


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("path", ["single", "batch"])
@pytest.mark.parametrize(
    "site,genotype",
    [(s, "".join(g)) for s in ("AT", "CG") for g in itertools.combinations_with_replacement(s, 2)],
)
@pytest.mark.parametrize("clinvar", [False, True])
def test_palindromic_site_keeps_the_exact_only_reading(
    cache, build, path, site, genotype, clinvar
):
    data = cache(site, clinvar=clinvar)
    rows = _find(data, _call(site, genotype, build), path)
    expected = [] if clinvar and genotype == site[0] * 2 else [f"row {genotype}"]
    assert _labels(rows) == expected


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("site,genotype", CASES)
def test_direct_site_controls_keep_the_written_genotype(cache, build, site, genotype):
    data = cache(site, reading="direct")
    assert _labels(_find(data, _call(site, genotype, build), "batch")) == [f"row {genotype}"]


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("path", ["single", "batch"])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("genotype", ["CC", "GG"])
def test_a_nonfinding_reverse_row_does_not_fall_back_to_the_wrong_exact_row(
    cache, build, path, reverse, genotype
):
    data = cache("ACG", reverse=reverse, nonfinding=_pair(complement(genotype)))
    assert not _find(data, _call("ACG", genotype, build), path)


@pytest.mark.parametrize("build", BUILDS)
@pytest.mark.parametrize("source_kind", ["vcf", "chip"])
@pytest.mark.parametrize("clinvar", [False, True])
@pytest.mark.parametrize("genotype", ["AA", "AC", "AG", "CC", "CG", "GG"])
def test_pipeline_rows_use_the_same_single_reading(
    tmp_path, build, source_kind, clinvar, genotype
):
    site = "GAC"
    data = _data(tmp_path / "data", site, clinvar=clinvar)
    source = tmp_path / ("input.vcf" if source_kind == "vcf" else "input.txt")
    if source_kind == "vcf":
        gt = "/".join(str(site.index(a)) for a in genotype)
        source.write_text(
            "##fileformat=VCFv4.2\n"
            + f"##contig=<ID=1,length={LENGTHS[build]}>\n"
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
            + f"1\t{POSITIONS[build]}\t{RS}\tG\tA,C\t.\tPASS\t.\tGT\t{gt}\n"
        )
    else:
        source.write_text(
            "# This data file generated by 23andMe\n"
            "# rsid\tchromosome\tposition\tgenotype\n"
            f"{RS}\t1\t{POSITIONS[build]}\t{genotype}\n"
        )
    out = tmp_path / "report.json"
    args = [
        "analyze",
        str(source),
        "--data-dir",
        str(data),
        "--no-update",
        "--build",
        build.lower(),
        "--min-magnitude",
        "0",
        "--output",
        str(out),
    ]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, (result.output, result.exception)
    rows = [a for a in json.loads(out.read_text())["annotations"] if a["source"] == "pharmgkb"]
    expected = [] if clinvar and genotype == "GG" else [f"row {_pair(complement(genotype))}"]
    assert sorted(a["description"].split(" — ")[-1] for a in rows) == expected
    assert all(a["genotype_match"] == genotype for a in rows)

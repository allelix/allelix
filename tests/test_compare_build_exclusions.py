"""O3: compare must apply each file's build decision before comparing calls."""

import pytest
from click.testing import CliRunner
from rich.console import Console

import allelix.compare as compare_module
from allelix.cli import main, utility


def vcf(path, kind, *, named, second, gvcf, build):
    length = 249250621 if build == "grch37" else 248956422
    chrom = "M" if kind == "name-M" else "chrM"
    declaration = {"hg19": 16571, "rcrs": 16569}.get(kind)
    meta = f"##contig=<ID=chr1,length={length}>\n"
    if declaration:
        meta += f"##contig=<ID={chrom},length={declaration}>\n"
    alt = "A,<NON_REF>" if gvcf else "A"
    rsid = "rs900000099" if named else "."
    path.write_text(
        "##fileformat=VCFv4.2\n"
        + meta
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        + f"chr1\t1000\trs900000001\tG\t{alt}\t.\tPASS\t.\tGT\t0/1\n"
        + f"{chrom}\t100\t{rsid}\tG\t{alt}\t.\tPASS\t.\tGT\t{'1/1' if second else '0/1'}\n"
    )
    return path


@pytest.mark.parametrize(
    "kinds",
    [
        ("hg19", "hg19"),
        ("hg19", "rcrs"),
        ("name-chrM", "rcrs"),
        ("name-M", "rcrs"),
        ("rcrs", "rcrs"),
    ],
)
@pytest.mark.parametrize("named", [False, True])
@pytest.mark.parametrize("explicit", [False, True])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("gvcf", [False, True])
def test_compare_obeys_each_input_mito_exclusion(
    tmp_path, monkeypatch, kinds, named, explicit, reverse, gvcf
):
    monkeypatch.setattr(utility, "console", Console(width=160, force_terminal=False))
    paths = [
        vcf(
            tmp_path / f"file{i}.vcf", kind, named=named, second=bool(i), gvcf=gvcf, build="grch37"
        )
        for i, kind in enumerate(kinds)
    ]
    if reverse:
        paths.reverse()
        kinds = kinds[::-1]
    observed = []
    original = compare_module.compare_variants

    def capture(first, second, **kwargs):
        first, second = list(first), list(second)
        observed.append((first, second))
        return original(first, second, **kwargs)

    monkeypatch.setattr(compare_module, "compare_variants", capture)
    args = ["compare", *(str(p) for p in paths)]
    if explicit:
        args += ["--build1", "grch37", "--build2", "grch37"]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, (result.output, result.exception)
    assert len(observed) == 1
    expected = [kind == "rcrs" for kind in kinds]
    for calls, eligible in zip(observed[0], expected, strict=True):
        assert sum(v.chromosome == "MT" for v in calls) == int(eligible)
        assert sum(v.chromosome == "1" for v in calls) == 1
    if all(expected):
        assert "50.00%" in result.output, result.output
    else:
        assert "100.00%" in result.output, result.output
        assert "│ MT " not in result.output
        assert "mitochondrial" in result.output.lower() or "hg19" in result.output.lower()


@pytest.mark.parametrize("named", [False, True])
@pytest.mark.parametrize("explicit", [False, True])
@pytest.mark.parametrize("kind", ["rcrs", "name-chrM", "name-M"])
def test_grch38_mito_call_is_not_discarded(tmp_path, monkeypatch, named, explicit, kind):
    monkeypatch.setattr(utility, "console", Console(width=160, force_terminal=False))
    paths = [
        vcf(
            tmp_path / f"file{i}.vcf",
            kind,
            named=named,
            second=bool(i),
            gvcf=False,
            build="grch38",
        )
        for i in range(2)
    ]
    args = ["compare", *(str(p) for p in paths)]
    if explicit:
        args += ["--build1", "grch38", "--build2", "grch38"]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, (result.output, result.exception)
    assert "50.00%" in result.output, result.output
    assert "│ MT " in result.output

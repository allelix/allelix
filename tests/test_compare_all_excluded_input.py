"""O3: build exclusions do not turn a nonempty file into an empty-input error."""

import pytest
from click.testing import CliRunner
from rich.console import Console
from test_compare_build_exclusions import vcf

import allelix.compare as compare_module
from allelix.cli import main, utility


@pytest.mark.parametrize("build", ["grch37", "grch38"])
@pytest.mark.parametrize("kind", ["hg19", "rcrs", "name-chrM", "name-M"])
@pytest.mark.parametrize("named", [False, True])
@pytest.mark.parametrize("explicit", [False, True])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("gvcf", [False, True])
def test_compare_counts_a_file_whose_records_are_all_excluded(
    tmp_path, monkeypatch, build, kind, named, explicit, reverse, gvcf
):
    monkeypatch.setattr(utility, "console", Console(width=180, force_terminal=False))
    candidate = vcf(
        tmp_path / "mitochondrial.vcf",
        kind,
        named=named,
        second=False,
        gvcf=gvcf,
        build=build,
    )
    candidate.write_text(
        "\n".join(
            line for line in candidate.read_text().splitlines() if not line.startswith("chr1\t")
        )
        + "\n"
    )
    control = vcf(
        tmp_path / "nuclear.vcf",
        "rcrs",
        named=True,
        second=False,
        gvcf=gvcf,
        build=build,
    )
    control.write_text(
        "\n".join(
            line for line in control.read_text().splitlines() if not line.startswith("chrM\t")
        )
        + "\n"
    )
    paths = [candidate, control]
    if reverse:
        paths.reverse()
    observed = []
    original = compare_module.compare_variants

    def capture(first, second, **kwargs):
        first, second = list(first), list(second)
        observed.append((first, second))
        return original(first, second, **kwargs)

    monkeypatch.setattr(compare_module, "compare_variants", capture)
    args = ["compare", *(str(path) for path in paths)]
    if explicit:
        args += ["--build1", build, "--build2", build]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, (
        "records excluded by the build decision are not an empty input",
        result.output,
        result.exception,
    )
    assert len(observed) == 1
    excluded = kind == "hg19" or (build == "grch37" and kind.startswith("name-"))
    candidate_calls = observed[0][int(reverse)]
    control_calls = observed[0][int(not reverse)]
    assert len(candidate_calls) == int(not excluded)
    assert len(control_calls) == 1 and control_calls[0].chromosome == "1"
    assert all(call.chromosome == "MT" for call in candidate_calls)
    assert "contains no genotype records" not in result.output
    if excluded:
        text = " ".join(result.output.split())
        heading = f"File {2 if reverse else 1} (mitochondrial.vcf):"
        notices = [part for part in text.split(heading)[1:] if "1 mitochondrial record(s)" in part]
        assert notices, "the correct file must disclose the excluded record count"
        assert "0 eligible variants" in text


@pytest.mark.parametrize("kind", ["empty", "structural"])
@pytest.mark.parametrize("reverse", [False, True])
def test_truly_unusable_compare_input_still_stops(tmp_path, monkeypatch, kind, reverse):
    monkeypatch.setattr(utility, "console", Console(width=180, force_terminal=False))
    bad = tmp_path / "unusable.vcf"
    bad.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=249250621>\n"
        '##ALT=<ID=DEL,Description="Deletion">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        + ("chr1\t1000\t.\tG\t<DEL>\t.\tPASS\tEND=2000\tGT\t0/1\n" if kind == "structural" else "")
    )
    good = vcf(tmp_path / "good.vcf", "rcrs", named=True, second=False, gvcf=False, build="grch37")
    paths = [bad, good]
    if reverse:
        paths.reverse()
    result = CliRunner().invoke(main, ["compare", *(str(path) for path in paths)])
    assert result.exit_code == 1
    assert ("no genotype records" if kind == "empty" else "structural-variant") in result.output

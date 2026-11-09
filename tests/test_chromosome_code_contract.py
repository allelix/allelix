"""Human chromosome codes have one meaning at shared lookup boundaries.

Synthetic only. Zero and unknown contigs are not sex-chromosome aliases.
"""

import pytest

from allelix.annotators.base import variant_locus
from allelix.compare import compare_variants
from allelix.models import Variant
from allelix.parsers._helpers import normalize_chromosome
from allelix.parsers.decodeme import DeCODEmeParser
from allelix.parsers.ftdna_famfinder import FTDNAFamFinderParser
from allelix.parsers.ftdna_illumina import FTDNAIlluminaParser
from allelix.parsers.vcf import VcfParser
from allelix.utils.build_detect import contig_length_verdict

ALIASES = [("23", "X"), ("24", "Y"), ("25", "X"), ("26", "MT"), ("XY", "X")]


@pytest.mark.parametrize(("code", "canonical"), ALIASES)
@pytest.mark.parametrize("prefix", ["", "chr", "CHR"])
def test_shared_chromosome_code_meaning(code, canonical, prefix):
    assert normalize_chromosome(prefix + code) == canonical


@pytest.mark.parametrize("chrom", ["0", "27", "100", "GL00023.1", "Un_XY", "", "1", "22"])
def test_unknown_contigs_and_autosomes_are_not_relabelled(chrom):
    assert normalize_chromosome(chrom) == chrom


@pytest.mark.parametrize(("code", "canonical"), ALIASES)
def test_variant_locus_meets_canonical_cache_keys(code, canonical):
    variant = Variant("rs900001", code, 65001, "G", "A", build="GRCh38")
    assert variant_locus(variant) == (canonical, 65001)


@pytest.mark.parametrize(("code", "canonical"), ALIASES)
@pytest.mark.parametrize("reader", ["vcf", "decodeme", "illumina", "famfinder"])
def test_shared_function_readers_yield_canonical_calls(tmp_path, code, canonical, reader):
    content = {
        "vcf": "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        f"{code}\t65001\trs900001\tG\tA\t.\tPASS\t.\tGT\t0/1\n",
        "decodeme": "Name,Variation,Chromosome,Position,Strand,YourCode\n"
        f"rs900001,G/A,{code},65001,+,GA\n",
        "illumina": f"RSID\tCHROMOSOME\tPOSITION\tRESULT\nrs900001\t{code}\t65001\tGA\n",
        "famfinder": "# Family Tree DNA - FamFinder\n"
        "RSID\tCHROMOSOME\tPOSITION\tALLELE1\tALLELE2\n"
        f"rs900001\t{code}\t65001\tG\tA\n",
    }[reader]
    parser = {
        "vcf": VcfParser,
        "decodeme": DeCODEmeParser,
        "illumina": FTDNAIlluminaParser,
        "famfinder": FTDNAFamFinderParser,
    }[reader]()
    path = tmp_path / "synthetic.txt"
    path.write_text(content)
    assert parser.can_parse(path)
    calls = list(parser.parse(path))
    assert len(calls) == 1
    assert calls[0].chromosome == canonical
    assert {calls[0].allele1, calls[0].allele2} == {"G", "A"}


@pytest.mark.parametrize(("code", "canonical"), ALIASES)
@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize("padded", [False, True])
def test_unnamed_compare_pairs_codes_with_letters(code, canonical, build, padded):
    kwargs = {
        "ref": "CG" if padded else "G",
        "input_alt": "CA" if padded else "A",
        "build": build,
        "ref_from_input": True,
    }
    a = Variant(
        "",
        code,
        65000 if padded else 65001,
        "CG" if padded else "G",
        "CA" if padded else "A",
        **kwargs,
    )
    b = Variant(
        "", canonical, 65001, "G", "A", ref="G", input_alt="A", build=build, ref_from_input=True
    )
    control = compare_variants([b], [b], build1=build, build2=build)
    assert control.shared == 1 and control.concordant == 1
    for first, second in [(a, b), (b, a)]:
        result = compare_variants([first], [second], build1=build, build2=build)
        assert result.shared == 1 and result.concordant == 1


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
@pytest.mark.parametrize(
    ("code", "canonical"), [("23", "X"), ("24", "Y"), ("25", "X"), ("XY", "X")]
)
def test_contig_length_build_vote_uses_same_codes(build, code, canonical):
    lengths = {
        "GRCh37": {"X": 155270560, "Y": 59373566},
        "GRCh38": {"X": 156040895, "Y": 57227415},
    }
    control = contig_length_verdict({canonical: lengths[build][canonical]})
    actual = contig_length_verdict({code: lengths[build][canonical]})
    assert control.build == build
    if code in {"25", "XY"}:
        # PAR labels name a region, not a full X chromosome: no length vote.
        assert actual.status == "none"
        assert actual.build is None
    else:
        assert actual.build == control.build


def test_numeric_mito_contig_keeps_hg19_exclusion():
    control = contig_length_verdict({"MT": 16571})
    actual = contig_length_verdict({"26": 16571})
    assert control.skip_mt
    assert actual.skip_mt

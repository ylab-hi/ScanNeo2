"""Tests for workflow/scripts/summarize.py (#227)."""

import csv
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "workflow/scripts"))

import summarize  # noqa: E402


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return str(path)


def table(path, header, rows):
    lines = ["\t".join(header)] + ["\t".join(r) for r in rows]
    return write(path, "\n".join(lines) + "\n")


def test_input_records_vcf_skips_header(tmp_path):
    p = write(tmp_path / "a.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\nchr1\t1\nchr1\t2\n")
    assert summarize.count_input_records(p) == 2


def test_input_records_fusion_tsv_hash_header(tmp_path):
    p = write(tmp_path / "f.tsv", "#gene1\tgene2\nA\tB\nC\tD\nE\tF\n")
    assert summarize.count_input_records(p) == 3


def test_input_records_protein_tsv_plain_header(tmp_path):
    p = write(tmp_path / "p.tsv", "id\twildtype_protein\tmutant_protein\nx\tMK\tMR\n")
    assert summarize.count_input_records(p) == 1


@pytest.mark.parametrize(
    "records, effects, neoepitopes, expected",
    [
        (0, 0, [0], "no_input"),
        (11, 0, [0], "no_effects"),
        (11, 4, [0], "no_neoepitopes"),
        (11, 4, [0, 3], "ok"),
        (11, 4, [7], "ok"),
    ],
)
def test_status(records, effects, neoepitopes, expected):
    assert summarize.status(records, effects, neoepitopes) == expected


def source_dir(tmp_path, source, effects_rows, epitopes):
    d = tmp_path / "prioritization" / source
    table(d / f"{source}_variant_effects.tsv", ["chrom"], [["chr1"]] * effects_rows)
    table(d / f"{source}_mhc-I_neoepitopes.txt", ["allele", "mt_epitope_seq"], epitopes)
    return str(d)


def test_write_summary(tmp_path):
    vcf = write(tmp_path / "exitrons.vcf", "#CHROM\nchr1\t1\nchr1\t2\n")
    entries = [
        {"sample": "S1", "source": "exitrons", "inputs": [vcf],
         "dir": source_dir(tmp_path, "exitrons", 2, [])},
        {"sample": "S1", "source": "somatic.snvs", "inputs": [vcf],
         "dir": source_dir(tmp_path, "somatic.snvs", 2,
                           [["A*01:01", "KLMNPQRST"], ["A*02:01", "KLMNPQRST"], ["A*01:01", "AAAAAAAAA"]])},
    ]
    out = tmp_path / "summary.tsv"
    summarize.write_summary(entries, ["I"], str(out))

    assert b"\r" not in out.read_bytes()  # plain \n lines for awk/cut
    rows = list(csv.DictReader(open(out), delimiter="\t"))
    assert list(rows[0]) == ["sample", "source", "input_records", "variant_effects",
                             "neoepitopes_mhc-I", "distinct_peptides_mhc-I", "status"]
    exitrons, snvs = rows
    assert (exitrons["input_records"], exitrons["variant_effects"],
            exitrons["neoepitopes_mhc-I"], exitrons["status"]) == ("2", "2", "0", "no_neoepitopes")
    assert (snvs["neoepitopes_mhc-I"], snvs["distinct_peptides_mhc-I"], snvs["status"]) == ("3", "2", "ok")


def test_parse_entry():
    e = summarize.parse_entry("S1|fusions|results/S1/prioritization/fusions/|a.tsv,b.tsv")
    assert e == {"sample": "S1", "source": "fusions",
                 "dir": "results/S1/prioritization/fusions/", "inputs": ["a.tsv", "b.tsv"]}

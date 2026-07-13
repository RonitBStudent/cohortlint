from __future__ import annotations

import csv
import gzip
from pathlib import Path

from .model import CohortLintError


def _write_fastq(path: Path, records: list[tuple[str, str, str]]) -> None:
    opener = gzip.open if path.suffix == ".gz" else path.open
    with opener(path, "wt", encoding="ascii", newline="") if path.suffix == ".gz" else opener("w", encoding="ascii", newline="") as handle:  # type: ignore[call-overload]
        for identifier, sequence, quality in records:
            handle.write(f"@{identifier}\n{sequence}\n+\n{quality}\n")


def write_demo(directory: Path) -> tuple[Path, Path]:
    if directory.exists() and any(directory.iterdir()):
        raise CohortLintError(f"demo directory is not empty: {directory}")
    directory.mkdir(parents=True, exist_ok=True)
    reads = directory / "reads"
    variants = directory / "variants"
    reference = directory / "reference"
    reads.mkdir()
    variants.mkdir()
    reference.mkdir()

    good_r1 = [
        ("INST:1:FC:1:1101:100:100 1:N:0:ACGT", "ACGTACGTACGT", "IIIIIIIIIIII"),
        ("INST:1:FC:1:1101:101:100 1:N:0:ACGT", "GCGTGCGTGCGT", "IIIIIIIIIIII"),
        ("INST:1:FC:1:1101:102:100 1:N:0:ACGT", "TTGCAATGCAAT", "IIIIIIIIIIII"),
    ]
    good_r2 = [(name.replace(" 1:", " 2:"), sequence[::-1], quality) for name, sequence, quality in good_r1]
    bad_r2 = list(good_r2)
    bad_r2[1] = ("INST:1:FC:1:1101:999:999 2:N:0:ACGT", bad_r2[1][1], bad_r2[1][2])
    _write_fastq(reads / "ALPHA_R1.fastq.gz", good_r1)
    _write_fastq(reads / "ALPHA_R2.fastq.gz", good_r2)
    _write_fastq(reads / "BETA_R1.fastq.gz", good_r1)
    _write_fastq(reads / "BETA_R2.fastq.gz", bad_r2)

    alpha_vcf = variants / "ALPHA.vcf"
    alpha_vcf.write_text(
        "##fileformat=VCFv4.3\n"
        "##contig=<ID=chr1,length=999>\n"
        "##contig=<ID=chr2,length=800>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tWRONG_ALPHA\n"
        "chr1\t101\t.\tA\tG\t60\tPASS\t.\tGT\t0/1\n",
        encoding="utf-8",
    )
    beta_vcf = variants / "BETA.vcf"
    beta_vcf.write_text(
        "##fileformat=VCFv4.3\n"
        "##contig=<ID=1,length=1000>\n"
        "##contig=<ID=2,length=800>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tBETA\n"
        "1\t205\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\n",
        encoding="utf-8",
    )

    fai = reference / "study.fa.fai"
    fai.write_text("chr1\t1000\t0\t60\t61\nchr2\t800\t1017\t60\t61\n", encoding="utf-8")

    manifest = directory / "cohort.csv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("sample_id", "site", "fastq_1", "fastq_2", "alignment", "variants"))
        writer.writerow(("ALPHA", "Site-East", "reads/ALPHA_R1.fastq.gz", "reads/ALPHA_R2.fastq.gz", "", "variants/ALPHA.vcf"))
        writer.writerow(("BETA", "Site-West", "reads/BETA_R1.fastq.gz", "reads/BETA_R2.fastq.gz", "", "variants/BETA.vcf"))
    return manifest, fai

from __future__ import annotations

import gzip
from pathlib import Path
import struct
import tempfile
import unittest
import zlib

from cohortlint.hts import inspect_bam, inspect_vcf
from cohortlint.model import CohortLintError, Severity
from cohortlint.reference import load_fai


BGZF_EOF = bytes.fromhex(
    "1f8b08040000000000ff0600424302001b0003000000000000000000"
)


def _bgzf_block(data: bytes) -> bytes:
    compressor = zlib.compressobj(level=6, wbits=-15)
    payload = compressor.compress(data) + compressor.flush()
    block_size = 18 + len(payload) + 8
    if block_size > 65_536:
        raise ValueError("test BGZF block is too large")
    header = (
        b"\x1f\x8b\x08\x04"
        + struct.pack("<I", 0)
        + b"\x00\xff"
        + struct.pack("<H", 6)
        + b"BC"
        + struct.pack("<HH", 2, block_size - 1)
    )
    trailer = struct.pack("<II", zlib.crc32(data), len(data) & 0xFFFFFFFF)
    return header + payload + trailer


def _bam_bytes(
    *,
    sample_names: tuple[str, ...] = ("sample-1",),
    contigs: tuple[tuple[str, int], ...] = (("chr1", 1000), ("chr2", 500)),
) -> bytes:
    header_lines = ["@HD\tVN:1.6\tSO:coordinate"]
    header_lines.extend(f"@SQ\tSN:{name}\tLN:{length}" for name, length in contigs)
    header_lines.extend(
        f"@RG\tID:rg{index}\tSM:{sample}"
        for index, sample in enumerate(sample_names, start=1)
    )
    header_text = ("\n".join(header_lines) + "\n").encode("utf-8")

    data = bytearray(b"BAM\x01")
    data.extend(struct.pack("<i", len(header_text)))
    data.extend(header_text)
    data.extend(struct.pack("<i", len(contigs)))
    for name, length in contigs:
        encoded_name = name.encode("utf-8") + b"\x00"
        data.extend(struct.pack("<i", len(encoded_name)))
        data.extend(encoded_name)
        data.extend(struct.pack("<i", length))
    return bytes(data)


class HtsInspectionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.fai = self.root / "reference.fa.fai"
        self.fai.write_text(
            "chr1\t1000\t0\t80\t81\nchr2\t500\t1013\t80\t81\n",
            encoding="utf-8",
        )
        self.reference = load_fai(self.fai)

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    @staticmethod
    def codes(inspection) -> set[str]:
        return {finding.code for finding in inspection.findings}

    def write_vcf(self, name: str, records: str, sample: str = "sample-1") -> Path:
        path = self.root / name
        content = (
            "##fileformat=VCFv4.3\n"
            "##contig=<ID=chr1,length=1000>\n"
            "##contig=<ID=chr2,length=500>\n"
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
            + records
        )
        if name.endswith(".gz"):
            with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
                handle.write(content)
        else:
            path.write_text(content, encoding="utf-8")
        return path

    def write_bam(
        self,
        name: str = "sample-1.bam",
        *,
        sample_names: tuple[str, ...] = ("sample-1",),
        contigs: tuple[tuple[str, int], ...] = (("chr1", 1000), ("chr2", 500)),
        eof: bool = True,
        index: bool = True,
    ) -> Path:
        path = self.root / name
        raw_bam = _bam_bytes(sample_names=sample_names, contigs=contigs)
        path.write_bytes(_bgzf_block(raw_bam) + (BGZF_EOF if eof else b""))
        if index:
            Path(f"{path}.bai").touch()
        return path

    def test_load_fai_fingerprints_only_name_length_dictionary(self) -> None:
        alternate = self.root / "alternate.fai"
        alternate.write_text(
            "chr1\t1000\t999\t20\t21\nchr2\t500\t3000\t20\t21\n",
            encoding="utf-8",
        )

        result = load_fai(alternate)

        self.assertEqual(result.contigs, (("chr1", 1000), ("chr2", 500)))
        self.assertEqual(result.fingerprint, self.reference.fingerprint)
        self.assertEqual(len(result.fingerprint), 64)

    def test_load_fai_rejects_duplicate_and_invalid_contigs(self) -> None:
        invalid = self.root / "invalid.fai"
        invalid.write_text("chr1\t1000\nchr1\t20\n", encoding="utf-8")

        with self.assertRaisesRegex(CohortLintError, "duplicate contig"):
            load_fai(invalid)

    def test_valid_vcf_reports_samples_contigs_and_scan_metrics(self) -> None:
        path = self.write_vcf(
            "sample.vcf",
            "chr1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n"
            "chr1\t20\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\n"
            "chr2\t5\t.\tN\tA\t.\tPASS\t.\tGT\t0/1\n",
        )

        result = inspect_vcf(
            path,
            sample_id="sample-1",
            reference=self.reference,
            full=True,
        )

        self.assertEqual(result.findings, ())
        self.assertEqual(result.sample_names, ("sample-1",))
        self.assertEqual(result.contigs, (("chr1", 1000), ("chr2", 500)))
        self.assertEqual(result.metrics["records_inspected"], 3)
        self.assertTrue(result.metrics["scan_complete"])

    def test_gzip_vcf_needs_index_and_obeys_record_budget(self) -> None:
        path = self.write_vcf(
            "sample.vcf.gz",
            "".join(
                f"chr1\t{position}\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n"
                for position in range(1, 6)
            ),
        )

        result = inspect_vcf(path, sample_id="sample-1", max_records=2)

        self.assertIn("HTS_VCF_INDEX_MISSING", self.codes(result))
        self.assertEqual(result.metrics["records_inspected"], 2)
        self.assertFalse(result.metrics["scan_complete"])
        self.assertTrue(result.metrics["compressed"])

    def test_vcf_detects_identity_reference_ref_and_sort_failures(self) -> None:
        path = self.write_vcf(
            "bad.vcf",
            "chr1\t20\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n"
            "chr1\t10\t.\tR\tT\t.\tPASS\t.\tGT\t0/1\n",
            sample="someone-else",
        )
        mismatched_fai = self.root / "mismatch.fai"
        mismatched_fai.write_text(
            "chr1\t999\nchr2\t500\n",
            encoding="utf-8",
        )

        result = inspect_vcf(
            path,
            sample_id="sample-1",
            reference=load_fai(mismatched_fai),
            full=True,
        )

        self.assertTrue(
            {
                "HTS_VCF_SAMPLE_MISMATCH",
                "HTS_VCF_REF_INVALID",
                "HTS_VCF_UNSORTED",
                "REF_LENGTH_MISMATCH",
            }.issubset(self.codes(result))
        )
        self.assertTrue(
            all(
                finding.sample_id == "sample-1"
                for finding in result.findings
            )
        )

    def test_corrupt_vcf_returns_findings_instead_of_raising(self) -> None:
        path = self.root / "broken.vcf.gz"
        path.write_bytes(b"\x1f\x8b\x08broken")

        result = inspect_vcf(path, sample_id="sample-1")

        self.assertIn("HTS_VCF_READ_ERROR", self.codes(result))
        self.assertIn("HTS_VCF_HEADER_MISSING", self.codes(result))
        self.assertFalse(result.metrics["scan_complete"])

    def test_vcf_requires_fileformat_declaration(self) -> None:
        path = self.root / "missing-fileformat.vcf"
        path.write_text(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample-1\n"
            "chr1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
            encoding="utf-8",
        )

        result = inspect_vcf(path, sample_id="sample-1", full=True)

        self.assertIn("HTS_VCF_FILEFORMAT_MISSING", self.codes(result))

    def test_vcf_requires_exact_mandatory_header_columns(self) -> None:
        path = self.root / "wrong-columns.vcf"
        path.write_text(
            "##fileformat=VCFv4.3\n"
            "#CHROM\tPOS\tID\tALT\tREF\tQUAL\tFILTER\tINFO\tFORMAT\tsample-1\n"
            "chr1\t1\t.\tG\tA\t.\tPASS\t.\tGT\t0/1\n",
            encoding="utf-8",
        )

        result = inspect_vcf(path, sample_id="sample-1", full=True)

        finding = next(
            finding
            for finding in result.findings
            if finding.code == "HTS_VCF_HEADER_INVALID"
        )
        self.assertEqual(finding.severity, Severity.ERROR)
        self.assertIn("expected", finding.detail)

    def test_plain_text_named_vcf_gz_reports_compression_mismatch(self) -> None:
        path = self.root / "mislabelled.vcf.gz"
        path.write_text(
            "##fileformat=VCFv4.3\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample-1\n"
            "chr1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
            encoding="utf-8",
        )

        result = inspect_vcf(path, sample_id="sample-1", full=True)

        self.assertIn("HTS_VCF_COMPRESSION_MISMATCH", self.codes(result))
        self.assertFalse(result.metrics["compressed"])
        self.assertTrue(result.metrics["scan_complete"])

    def test_observed_contigs_are_not_reported_as_a_sequence_dictionary(self) -> None:
        path = self.root / "no-contig-dictionary.vcf"
        path.write_text(
            "##fileformat=VCFv4.3\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample-1\n"
            "chr1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
            encoding="utf-8",
        )

        result = inspect_vcf(path, sample_id="sample-1", full=True)

        self.assertEqual(result.contigs, ())
        self.assertEqual(result.metrics["observed_contigs"], ("chr1",))

    def test_valid_bgzf_bam_extracts_sample_and_reference_dictionary(self) -> None:
        path = self.write_bam()

        result = inspect_bam(path, sample_id="sample-1", reference=self.reference)

        self.assertEqual(result.findings, ())
        self.assertEqual(result.sample_names, ("sample-1",))
        self.assertEqual(result.contigs, (("chr1", 1000), ("chr2", 500)))
        self.assertTrue(result.metrics["bgzf"])
        self.assertTrue(result.metrics["bgzf_eof"])
        self.assertTrue(result.metrics["indexed"])

    def test_bam_detects_sample_reference_eof_and_index_failures(self) -> None:
        path = self.write_bam(
            "bad.bam",
            sample_names=("someone-else",),
            contigs=(("chr1", 999), ("chr3", 100)),
            eof=False,
            index=False,
        )

        result = inspect_bam(path, sample_id="sample-1", reference=self.reference)

        expected = {
            "HTS_BAM_SAMPLE_MISMATCH",
            "HTS_BAM_EOF_MISSING",
            "HTS_BAM_INDEX_MISSING",
            "REF_CONTIG_UNKNOWN",
            "REF_CONTIG_MISSING",
            "REF_LENGTH_MISMATCH",
        }
        self.assertTrue(expected.issubset(self.codes(result)))
        eof_finding = next(
            finding
            for finding in result.findings
            if finding.code == "HTS_BAM_EOF_MISSING"
        )
        self.assertEqual(eof_finding.severity, Severity.WARNING)

    def test_corrupt_bam_returns_findings_instead_of_raising(self) -> None:
        path = self.root / "broken.bam"
        path.write_bytes(_bgzf_block(b"not a BAM") + BGZF_EOF)

        result = inspect_bam(path, sample_id="sample-1")

        self.assertIn("HTS_BAM_READ_ERROR", self.codes(result))
        self.assertIn("HTS_BAM_SAMPLE_MISSING", self.codes(result))
        self.assertIn("HTS_BAM_INDEX_MISSING", self.codes(result))

    def test_missing_hts_files_return_unreadable_findings(self) -> None:
        vcf = inspect_vcf(self.root / "missing.vcf", sample_id="sample-1")
        bam = inspect_bam(self.root / "missing.bam", sample_id="sample-1")

        self.assertIn("HTS_VCF_UNREADABLE", self.codes(vcf))
        self.assertIn("HTS_BAM_UNREADABLE", self.codes(bam))


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from cohortlint.external import inspect_cram
from cohortlint.reference import load_fai


class InspectCramTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.cram = self.root / "S1.cram"
        self.cram.touch()
        Path(f"{self.cram}.crai").touch()

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    @staticmethod
    def codes(inspection) -> set[str]:
        return {finding.code for finding in inspection.findings}

    def inspect_header(self, header_text: str, *, sample_id: str = "S1", reference=None):
        quickcheck = subprocess.CompletedProcess(
            args=["samtools", "quickcheck"],
            returncode=0,
            stdout="",
            stderr="",
        )
        header = subprocess.CompletedProcess(
            args=["samtools", "view"],
            returncode=0,
            stdout=header_text,
            stderr="",
        )
        with (
            patch("cohortlint.external.tool_path", return_value="/mock/samtools"),
            patch("cohortlint.external.tool_version", return_value="samtools 1.20"),
            patch("cohortlint.external.subprocess.run", side_effect=(quickcheck, header)),
        ):
            return inspect_cram(
                self.cram,
                sample_id=sample_id,
                reference=reference,
            )

    def test_valid_header_matches_sample_and_reference(self) -> None:
        fai = self.root / "reference.fa.fai"
        fai.write_text("chr1\t100\nchr2\t200\n", encoding="utf-8")

        result = self.inspect_header(
            "@SQ\tSN:chr1\tLN:100\n"
            "@SQ\tSN:chr2\tLN:200\n"
            "@RG\tID:rg1\tSM:S1\n",
            reference=load_fai(fai),
        )

        self.assertEqual(result.findings, ())
        self.assertEqual(result.sample_names, ("S1",))
        self.assertEqual(result.contigs, (("chr1", 100), ("chr2", 200)))
        self.assertEqual(result.metrics["samtools"], "samtools 1.20")

    def test_header_process_errors_become_findings(self) -> None:
        quickcheck = subprocess.CompletedProcess(
            args=["samtools", "quickcheck"],
            returncode=0,
            stdout="",
            stderr="",
        )
        errors = (
            OSError("samtools disappeared"),
            subprocess.TimeoutExpired("samtools view -H", 30),
        )

        for error in errors:
            with self.subTest(error=type(error).__name__):
                with (
                    patch(
                        "cohortlint.external.tool_path",
                        return_value="/mock/samtools",
                    ),
                    patch(
                        "cohortlint.external.tool_version",
                        return_value="samtools 1.20",
                    ),
                    patch(
                        "cohortlint.external.subprocess.run",
                        side_effect=(quickcheck, error),
                    ),
                ):
                    result = inspect_cram(self.cram, sample_id="S1")

                self.assertIn("CRAM_HEADER_UNREADABLE", self.codes(result))
                finding = next(
                    item
                    for item in result.findings
                    if item.code == "CRAM_HEADER_UNREADABLE"
                )
                self.assertIn(type(error).__name__, finding.detail)

    def test_missing_multiple_and_mismatched_samples_are_reported(self) -> None:
        multiple = self.inspect_header(
            "@SQ\tSN:chr1\tLN:100\n"
            "@RG\tID:rg1\tSM:S1\n"
            "@RG\tID:rg2\tSM:S2\n"
            "@RG\tID:rg3\n"
        )
        self.assertIn("CRAM_MULTIPLE_SAMPLES", self.codes(multiple))
        self.assertIn("CRAM_READ_GROUP_SAMPLE_MISSING", self.codes(multiple))
        self.assertNotIn("CRAM_SAMPLE_MISMATCH", self.codes(multiple))

        mismatched = self.inspect_header(
            "@SQ\tSN:chr1\tLN:100\n@RG\tID:rg1\tSM:other\n"
        )
        self.assertIn("CRAM_SAMPLE_MISMATCH", self.codes(mismatched))

        missing = self.inspect_header(
            "@SQ\tSN:chr1\tLN:100\n@RG\tID:rg1\n"
        )
        self.assertIn("CRAM_READ_GROUP_SAMPLE_MISSING", self.codes(missing))
        self.assertIn("CRAM_SAMPLE_MISSING", self.codes(missing))

    def test_reference_name_length_and_order_mismatches_are_reported(self) -> None:
        fai = self.root / "reference.fa.fai"
        fai.write_text(
            "chr1\t100\nchr2\t200\nchr3\t300\n",
            encoding="utf-8",
        )

        result = self.inspect_header(
            "@SQ\tSN:chr2\tLN:201\n"
            "@SQ\tSN:chr1\tLN:100\n"
            "@SQ\tSN:chr4\tLN:400\n"
            "@RG\tID:rg1\tSM:S1\n",
            reference=load_fai(fai),
        )

        expected = {
            "REF_CONTIG_UNKNOWN",
            "REF_CONTIG_MISSING",
            "REF_LENGTH_MISMATCH",
            "REF_CONTIG_ORDER_MISMATCH",
        }
        self.assertTrue(expected.issubset(self.codes(result)))

    def test_malformed_and_duplicate_contigs_do_not_raise(self) -> None:
        result = self.inspect_header(
            "@SQ\tLN:100\n"
            "@SQ\tSN:invalid-length\tLN:nope\n"
            "@SQ\tSN:chr1\tLN:100\n"
            "@SQ\tSN:chr1\tLN:100\n"
            "@SQ\tSN:chr2\tSN:chr2\tLN:200\n"
            "@RG\tID:rg1\tSM:S1\n"
        )

        self.assertIn("CRAM_CONTIG_HEADER_INVALID", self.codes(result))
        self.assertIn("CRAM_CONTIG_DUPLICATE", self.codes(result))
        self.assertEqual(result.contigs, (("chr1", 100),))

    def test_reference_requires_a_valid_cram_dictionary(self) -> None:
        fai = self.root / "reference.fa.fai"
        fai.write_text("chr1\t100\n", encoding="utf-8")

        result = self.inspect_header(
            "@RG\tID:rg1\tSM:S1\n",
            reference=load_fai(fai),
        )

        self.assertIn("CRAM_CONTIG_DICTIONARY_MISSING", self.codes(result))


if __name__ == "__main__":
    unittest.main()

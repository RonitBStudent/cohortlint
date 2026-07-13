from __future__ import annotations

import gzip
from pathlib import Path
import tempfile
import unittest

from cohortlint.fastq import inspect_fastq_pair
from cohortlint.model import Severity


def _record(name: str, sequence: str = "ACGTN", quality: str = "IIIII") -> str:
    return f"@{name}\n{sequence}\n+\n{quality}\n"


class InspectFastqPairTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def write(self, name: str, content: str) -> Path:
        path = self.root / name
        path.write_text(content, encoding="ascii")
        return path

    def write_gzip(self, name: str, content: str) -> Path:
        path = self.root / name
        with gzip.open(path, "wt", encoding="ascii", newline="") as handle:
            handle.write(content)
        return path

    def codes(self, inspection) -> set[str]:
        return {finding.code for finding in inspection.findings}

    def test_valid_single_end_file_reports_streaming_metrics(self) -> None:
        path = self.write(
            "reads.fastq",
            _record("read-a", "ACGTN", "!!!!!")
            + _record("read-b", "GGCC", "####"),
        )

        result = inspect_fastq_pair(path, None, sample_id="sample-1", full=True)

        self.assertEqual(result.kind, "fastq")
        self.assertEqual(result.sample_names, ("sample-1",))
        self.assertEqual(result.findings, ())
        self.assertEqual(result.metrics["records_1"], 2)
        self.assertEqual(result.metrics["bases_scanned"], 9)
        self.assertEqual(result.metrics["read_length_min"], 4)
        self.assertEqual(result.metrics["read_length_max"], 5)
        self.assertAlmostEqual(result.metrics["gc_fraction"], 6 / 9, places=6)
        self.assertTrue(result.metrics["scan_complete"])

    def test_pairs_slash_and_illumina_style_identifiers(self) -> None:
        first = self.write_gzip(
            "sample_R1.fastq.gz",
            _record("instrument:1:flowcell:2:1101:1:1/1")
            + _record("instrument:1:flowcell:2:1101:1:2 1:N:0:ACGT"),
        )
        second = self.write_gzip(
            "sample_R2.fastq.gz",
            _record("instrument:1:flowcell:2:1101:1:1/2")
            + _record("instrument:1:flowcell:2:1101:1:2 2:N:0:ACGT"),
        )

        result = inspect_fastq_pair(first, second, sample_id="paired", full=True)

        self.assertEqual(result.findings, ())
        self.assertEqual(result.metrics["pair_ids_compared"], 2)
        self.assertEqual(result.metrics["pair_id_mismatches"], 0)
        self.assertTrue(result.metrics["gzip_1"])
        self.assertTrue(result.metrics["gzip_2"])

    def test_mismatched_pair_ids_are_errors(self) -> None:
        first = self.write("R1.fastq", _record("read-a/1") + _record("read-b/1"))
        second = self.write("R2.fastq", _record("read-a/2") + _record("read-X/2"))

        result = inspect_fastq_pair(first, second, sample_id="shifted", full=True)

        mismatch = next(
            finding
            for finding in result.findings
            if finding.code == "FASTQ_PAIR_ID_MISMATCH"
        )
        self.assertEqual(mismatch.severity, Severity.ERROR)
        self.assertIn("read-b", mismatch.detail)
        self.assertIn("read-X", mismatch.detail)
        self.assertEqual(result.metrics["pair_id_mismatches"], 1)

    def test_record_content_and_structure_failures_are_actionable(self) -> None:
        path = self.write(
            "bad.fastq",
            "read-a\nACGT!\nseparator\nIII\n",
        )

        result = inspect_fastq_pair(path, None, sample_id="bad", full=True)

        self.assertTrue(
            {
                "FASTQ_INVALID_HEADER",
                "FASTQ_INVALID_SEPARATOR",
                "FASTQ_INVALID_SEQUENCE_SYMBOL",
                "FASTQ_SEQUENCE_QUALITY_LENGTH_MISMATCH",
            }.issubset(self.codes(result))
        )
        self.assertTrue(all(finding.remediation for finding in result.findings))

    def test_truncated_record_is_detected(self) -> None:
        path = self.write("truncated.fastq", "@read-a\nACGT\n+\n")

        result = inspect_fastq_pair(path, None, sample_id="truncated", full=True)

        self.assertIn("FASTQ_TRUNCATED_RECORD", self.codes(result))
        self.assertEqual(result.metrics["records_1"], 0)

    def test_full_scan_checks_gzip_trailer_integrity(self) -> None:
        path = self.write_gzip(
            "corrupt.fastq.gz",
            "".join(_record(f"read-{index}") for index in range(5000)),
        )
        damaged = path.read_bytes()[:-5]
        path.write_bytes(damaged)

        result = inspect_fastq_pair(path, None, sample_id="corrupt", full=True)

        self.assertIn("FASTQ_GZIP_INTEGRITY_ERROR", self.codes(result))
        self.assertFalse(result.metrics["scan_complete"])

    def test_full_scan_detects_different_pair_counts(self) -> None:
        first = self.write("R1.fastq", _record("a/1") + _record("b/1"))
        second = self.write("R2.fastq", _record("a/2"))

        result = inspect_fastq_pair(first, second, sample_id="orphan", full=True)

        self.assertIn("FASTQ_PAIR_COUNT_MISMATCH", self.codes(result))
        self.assertEqual(result.metrics["records_1"], 2)
        self.assertEqual(result.metrics["records_2"], 1)

    def test_sample_mode_obeys_record_budget(self) -> None:
        path = self.write(
            "long.fastq",
            "".join(_record(f"read-{index}") for index in range(10)),
        )

        result = inspect_fastq_pair(
            path,
            None,
            sample_id="sampled",
            max_records=3,
        )

        self.assertEqual(result.metrics["records_1"], 3)
        self.assertFalse(result.metrics["scan_complete"])

    def test_missing_file_is_reported_without_traceback(self) -> None:
        result = inspect_fastq_pair(
            self.root / "missing.fastq.gz",
            None,
            sample_id="missing",
            full=True,
        )

        self.assertIn("FASTQ_FILE_NOT_FOUND", self.codes(result))
        self.assertEqual(result.metrics["records_1"], 0)

    def test_invalid_record_budget_is_rejected(self) -> None:
        path = self.write("reads.fastq", _record("read-a"))

        with self.assertRaisesRegex(ValueError, "at least 1"):
            inspect_fastq_pair(path, None, sample_id="sample", max_records=0)


if __name__ == "__main__":
    unittest.main()

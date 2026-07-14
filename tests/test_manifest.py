from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from cohortlint.manifest import MANIFEST_FIELDS, discover, load_manifest, write_manifest
from cohortlint.model import ManifestRow, Severity


class ManifestTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def touch(self, relative_path: str) -> Path:
        path = self.root / relative_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
        return path.resolve()

    def write_csv(self, rows: list[list[str]], name: str = "manifest.csv") -> Path:
        path = self.root / name
        with path.open("w", encoding="utf-8", newline="") as handle:
            csv.writer(handle).writerows(rows)
        return path

    @staticmethod
    def codes(findings) -> set[str]:
        return {finding.code for finding in findings}

    def test_load_resolves_relative_paths_against_manifest(self) -> None:
        read_1 = self.touch("reads/subject_R1.fastq.gz")
        read_2 = self.touch("reads/subject_R2.fastq.gz")
        manifest = self.write_csv(
            [
                list(MANIFEST_FIELDS),
                [
                    "subject",
                    "north",
                    "reads/subject_R1.fastq.gz",
                    "reads/subject_R2.fastq.gz",
                    "",
                    "",
                ],
            ]
        )

        rows, findings = load_manifest(manifest)

        self.assertEqual(findings, ())
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].fastq_1, str(read_1))
        self.assertEqual(rows[0].fastq_2, str(read_2))
        self.assertEqual(rows[0].site, "north")

    def test_load_reports_header_id_path_and_ownership_problems(self) -> None:
        shared = self.touch("shared.bam")
        manifest = self.write_csv(
            [
                ["sample_id", "site", "fastq_1", "fastq_2", "alignment", "extra"],
                ["alpha", "a", "", "", str(shared), "ignored"],
                ["ALPHA", "a", "", "", "", "ignored"],
                ["", "b", "", "", "", "ignored"],
                ["beta", "b", "missing.fastq", "", str(shared), "ignored"],
            ]
        )

        rows, findings = load_manifest(manifest)
        codes = self.codes(findings)

        self.assertEqual([row.sample_id for row in rows], ["alpha", "beta"])
        self.assertIn("MANIFEST_HEADERS_MISSING", codes)
        self.assertIn("MANIFEST_HEADERS_EXTRA", codes)
        self.assertIn("MANIFEST_SAMPLE_ID_DUPLICATE", codes)
        self.assertIn("MANIFEST_SAMPLE_ID_BLANK", codes)
        self.assertIn("PATH_NOT_FOUND", codes)
        self.assertIn("DUPLICATE_FILE_OWNERSHIP", codes)
        self.assertTrue(any(f.severity == Severity.ERROR for f in findings))

    def test_load_checks_fastq_mate_names_and_roles(self) -> None:
        read_1 = self.touch("reads/alpha_R1_001.fastq.gz")
        wrong_read_2 = self.touch("reads/beta_R2_002.fastq.gz")
        manifest = self.write_csv(
            [
                list(MANIFEST_FIELDS),
                ["alpha", "", str(read_1), str(wrong_read_2), "", ""],
            ]
        )

        _rows, findings = load_manifest(manifest)

        self.assertIn("FASTQ_MATES_MISMATCH", self.codes(findings))

    def test_load_accepts_an_unmarked_single_end_fastq(self) -> None:
        single = self.touch("reads/alpha.fastq.gz")
        manifest = self.write_csv(
            [list(MANIFEST_FIELDS), ["alpha", "", str(single), "", "", ""]]
        )

        _rows, findings = load_manifest(manifest)

        self.assertNotIn("FASTQ_MATE_2_MISSING", self.codes(findings))

    def test_discover_groups_supported_files_by_conservative_sample_id(self) -> None:
        read_1 = self.touch("input/alpha_R1.fastq.gz")
        read_2 = self.touch("input/alpha_R2.fastq.gz")
        alignment = self.touch("input/alpha.bam")
        variants = self.touch("input/alpha.vcf.gz")
        single = self.touch("input/beta.fq")
        self.touch("input/notes.txt")

        rows, findings = discover(self.root / "input")

        self.assertEqual(findings, ())
        self.assertEqual([row.sample_id for row in rows], ["alpha", "beta"])
        alpha, beta = rows
        self.assertEqual(alpha.fastq_1, str(read_1))
        self.assertEqual(alpha.fastq_2, str(read_2))
        self.assertEqual(alpha.alignment, str(alignment))
        self.assertEqual(alpha.variants, str(variants))
        self.assertEqual(beta.fastq_1, str(single))

    def test_discover_understands_numeric_mates_and_chunk_suffix(self) -> None:
        read_1 = self.touch("input/gamma_1_001.fq.gz")
        read_2 = self.touch("input/gamma_2_001.fq.gz")

        rows, findings = discover(self.root / "input")

        self.assertEqual(findings, ())
        self.assertEqual(rows[0].sample_id, "gamma")
        self.assertEqual(rows[0].fastq_1, str(read_1))
        self.assertEqual(rows[0].fastq_2, str(read_2))

    def test_discover_flags_multilane_fastq_as_unsupported(self) -> None:
        first_r1 = self.touch("input/P1_L001_R1_001.fastq.gz")
        first_r2 = self.touch("input/P1_L001_R2_001.fastq.gz")
        self.touch("input/P1_L002_R1_001.fastq.gz")
        self.touch("input/P1_L002_R2_001.fastq.gz")

        rows, findings = discover(self.root / "input")

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].sample_id, "P1")
        self.assertEqual(rows[0].fastq_1, str(first_r1))
        self.assertEqual(rows[0].fastq_2, str(first_r2))
        self.assertIn("DISCOVERY_MULTIPART_FASTQ_UNSUPPORTED", self.codes(findings))
        self.assertIn("DISCOVERY_ROLE_AMBIGUOUS", self.codes(findings))

    def test_discover_warns_when_one_role_has_multiple_candidates(self) -> None:
        first = self.touch("input/alpha.bam")
        self.touch("input/alpha.cram")

        rows, findings = discover(self.root / "input")

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].alignment, str(first))
        self.assertIn("DISCOVERY_ROLE_AMBIGUOUS", self.codes(findings))

    def test_discover_honors_non_recursive_mode(self) -> None:
        self.touch("input/nested/alpha.bam")

        rows, findings = discover(self.root / "input", recursive=False)

        self.assertEqual(rows, ())
        self.assertIn("DISCOVERY_NO_SUPPORTED_FILES", self.codes(findings))

    def test_write_manifest_uses_stable_schema_and_round_trips(self) -> None:
        alignment = self.touch("alpha.bam")
        destination = self.root / "nested" / "manifest.csv"
        write_manifest(
            [ManifestRow(sample_id="alpha", site="east", alignment=str(alignment))],
            destination,
        )

        with destination.open("r", encoding="utf-8", newline="") as handle:
            records = list(csv.reader(handle))
        self.assertEqual(records[0], list(MANIFEST_FIELDS))

        rows, findings = load_manifest(destination)
        self.assertEqual(findings, ())
        self.assertEqual(rows[0].sample_id, "alpha")
        self.assertEqual(rows[0].alignment, str(alignment))

    def test_joint_called_vcf_can_be_shared_by_manifest_rows(self) -> None:
        shared_vcf = self.root / "cohort.vcf"
        shared_vcf.write_text("##fileformat=VCFv4.3\n", encoding="utf-8")
        manifest = self.root / "cohort.csv"
        manifest.write_text(
            "sample_id,site,fastq_1,fastq_2,alignment,variants\n"
            "A,East,,,,cohort.vcf\n"
            "B,West,,,,cohort.vcf\n",
            encoding="utf-8",
        )

        rows, findings = load_manifest(manifest)

        self.assertEqual(len(rows), 2)
        self.assertNotIn("DUPLICATE_FILE_OWNERSHIP", {finding.code for finding in findings})

    def test_invalid_manifest_path_becomes_a_finding(self) -> None:
        manifest = self.write_csv(
            [list(MANIFEST_FIELDS), ["alpha", "", "bad\0path.fastq", "", "", ""]]
        )

        rows, findings = load_manifest(manifest)

        self.assertEqual(len(rows), 1)
        self.assertIn("MANIFEST_PATH_INVALID", self.codes(findings))


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import contextlib
import csv
import io
import json
from pathlib import Path
import tempfile
import unittest

from cohortlint.cli import main
from cohortlint.demo import write_demo
from cohortlint.engine import check_cohort
from cohortlint.output import render_json


class EngineAndCliTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def invoke(self, *arguments: str) -> tuple[int, str, str]:
        stdout = io.StringIO()
        stderr = io.StringIO()
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            status = main(arguments)
        return status, stdout.getvalue(), stderr.getvalue()

    def test_demo_exposes_cross_lab_failures(self) -> None:
        manifest, reference = write_demo(self.root / "demo")

        report = check_cohort(manifest, reference_fai=reference, full=True)
        codes = {finding.code for finding in report.findings}

        self.assertFalse(report.passed)
        self.assertIn("FASTQ_PAIR_ID_MISMATCH", codes)
        self.assertIn("HTS_VCF_SAMPLE_MISMATCH", codes)
        self.assertIn("REF_LENGTH_MISMATCH", codes)
        self.assertIn("COHORT_CONTIG_SET_MISMATCH", codes)

    def test_compatible_joint_vcf_passes(self) -> None:
        vcf = self.root / "joint.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.3\n"
            "##contig=<ID=chr1,length=1000>\n"
            "##contig=<ID=chr2,length=800>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
            "chr1\t100\t.\tA\tG\t60\tPASS\t.\tGT\t0/1\t0/0\n",
            encoding="utf-8",
        )
        reference = self.root / "reference.fa.fai"
        reference.write_text("chr1\t1000\t0\t60\t61\nchr2\t800\t0\t60\t61\n", encoding="utf-8")
        manifest = self.root / "cohort.csv"
        with manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(("sample_id", "site", "fastq_1", "fastq_2", "alignment", "variants"))
            writer.writerow(("A", "East", "", "", "", "joint.vcf"))
            writer.writerow(("B", "West", "", "", "", "joint.vcf"))

        report = check_cohort(manifest, reference_fai=reference, full=True)

        self.assertTrue(report.passed, [finding.as_dict() for finding in report.findings])
        self.assertEqual(report.errors, 0)
        self.assertEqual(report.sample_count, 2)
        self.assertEqual(report.file_count, 1)

    def test_json_report_has_stable_machine_readable_shape(self) -> None:
        manifest, reference = write_demo(self.root / "json-demo")
        report = check_cohort(manifest, reference_fai=reference, full=True)

        payload = json.loads(render_json(report))

        self.assertEqual(payload["schema_version"], "1.0")
        self.assertEqual(payload["status"], "fail")
        self.assertEqual(payload["summary"]["samples"], 2)
        self.assertTrue(payload["summary"]["reference_dictionary_fingerprint"])
        self.assertTrue(all("code" in finding for finding in payload["findings"]))

    def test_check_command_uses_ci_exit_codes(self) -> None:
        manifest, reference = write_demo(self.root / "cli-demo")

        status, stdout, stderr = self.invoke(
            "check", str(manifest), "--reference", str(reference), "--full",
        )

        self.assertEqual(status, 1, stderr)
        self.assertIn("COHORT NOT READY", stdout)
        self.assertNotIn("Traceback", stderr)

    def test_demo_and_doctor_are_self_contained(self) -> None:
        status, stdout, stderr = self.invoke("doctor")
        self.assertEqual(status, 0, stderr)
        self.assertIn("no runtime dependencies", stdout)

        status, stdout, stderr = self.invoke("demo", "--output", str(self.root / "cli-generated"))
        self.assertEqual(status, 0, stderr)
        self.assertIn("COHORT NOT READY", stdout)
        self.assertTrue((self.root / "cli-generated" / "cohort.csv").exists())

        status, stdout, stderr = self.invoke("demo", "--output", str(self.root / "cli-generated"))
        self.assertEqual(status, 0, stderr)
        self.assertIn("COHORT NOT READY", stdout)
        self.assertIn("reusing demonstration cohort", stderr)

    def test_missing_manifest_is_an_invocation_error(self) -> None:
        status, stdout, stderr = self.invoke("check", str(self.root / "missing.csv"))

        self.assertEqual(status, 2)
        self.assertEqual(stdout, "")
        self.assertIn("cohortlint: error:", stderr)
        self.assertNotIn("Traceback", stderr)


if __name__ == "__main__":
    unittest.main()

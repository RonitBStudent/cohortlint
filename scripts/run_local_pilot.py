#!/usr/bin/env python3
"""Black-box owner-device pilot for an installed CohortLint executable."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import platform
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time


class PilotFailure(RuntimeError):
    """A pilot assertion failed."""


def _run(executable: str, *arguments: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [executable, *arguments],
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )


def _require(
    condition: bool,
    message: str,
    result: subprocess.CompletedProcess[str] | None = None,
) -> None:
    if condition:
        return
    details = ""
    if result is not None:
        details = (
            f"\nexit={result.returncode}"
            f"\nstdout:\n{result.stdout}"
            f"\nstderr:\n{result.stderr}"
        )
    raise PilotFailure(message + details)


def _write_fastq(path: Path, mate: int, count: int) -> None:
    with gzip.open(path, "wt", encoding="ascii", newline="") as handle:
        for index in range(1, count + 1):
            handle.write(
                f"@read-{index}/{mate}\n"
                "ACGTACGT\n"
                "+\n"
                "IIIIIIII\n"
            )


def _write_vcf(path: Path, sample: str, count: int) -> None:
    records = "".join(
        f"chr1\t{index * 10}\t.\tA\tG\t60\tPASS\t.\tGT\t0/1\n"
        for index in range(1, count + 1)
    )
    path.write_text(
        "##fileformat=VCFv4.3\n"
        "##contig=<ID=chr1,length=1000>\n"
        "##contig=<ID=chr2,length=800>\n"
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
        + records,
        encoding="utf-8",
    )


def _write_joint_vcf(path: Path, samples: list[str], count: int) -> None:
    genotype_header = "\t".join(samples)
    genotypes = "\t".join("0/1" for _ in samples)
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write(
            "##fileformat=VCFv4.3\n"
            "##contig=<ID=chr1,length=1000000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            f"{genotype_header}\n"
        )
        for position in range(1, count + 1):
            handle.write(
                f"chr1\t{position}\t.\tA\tG\t60\tPASS\t.\tGT\t{genotypes}\n"
            )


def _write_manifest(
    path: Path,
    *,
    sample: str,
    read_1: Path,
    read_2: Path,
    variants: Path,
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ("sample_id", "site", "fastq_1", "fastq_2", "alignment", "variants")
        )
        writer.writerow(
            (sample, "Pilot-Lab", read_1.name, read_2.name, "", variants.name)
        )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def _doctor_python(stdout: str) -> tuple[str, str]:
    match = re.search(r"^\s*Python\s+(\S+)\s+\((.+)\)\s*$", stdout, re.MULTILINE)
    if match is None:
        raise PilotFailure("doctor output did not identify the target Python")
    return match.group(1), match.group(2)


def _prepare_workspace(requested: Path | None) -> tuple[Path, bool]:
    if requested is None:
        return Path(tempfile.mkdtemp(prefix="cohortlint pilot ")), True
    workspace = requested.expanduser().resolve(strict=False)
    if workspace.exists() and any(workspace.iterdir()):
        raise PilotFailure(f"pilot workspace is not empty: {workspace}")
    workspace.mkdir(parents=True, exist_ok=True)
    return workspace, False


def run_pilot(executable: str, workspace: Path) -> dict[str, object]:
    scenarios: list[str] = []
    started = time.perf_counter()

    version = _run(executable, "--version")
    _require(version.returncode == 0 and version.stdout.startswith("cohortlint "), "version failed", version)
    scenarios.append("version")

    doctor = _run(executable, "doctor")
    _require(doctor.returncode == 0 and "Core       ready" in doctor.stdout, "doctor failed", doctor)
    target_python, target_python_executable = _doctor_python(doctor.stdout)
    scenarios.append("doctor")

    schema = _run(executable, "schema")
    _require(schema.returncode == 0 and "sample_id,site" in schema.stdout, "schema failed", schema)
    scenarios.append("schema")

    demo_directory = workspace / "demo with spaces"
    demo_first = _run(executable, "demo", "--output", str(demo_directory))
    _require(
        demo_first.returncode == 0
        and "CHECK RESULT: BLOCKING FINDINGS" in demo_first.stdout,
        "first demo run failed",
        demo_first,
    )
    demo_second = _run(executable, "demo", "--output", str(demo_directory))
    _require(
        demo_second.returncode == 0 and demo_second.stdout == demo_first.stdout,
        "demo reuse was not deterministic",
        demo_second,
    )
    scenarios.append("repeatable demo and paths with spaces")

    clean = workspace / "clean cohort"
    clean.mkdir()
    read_1 = clean / "S1_R1.fastq.gz"
    read_2 = clean / "S1_R2.fastq.gz"
    variants = clean / "S1.vcf"
    manifest = clean / "cohort.csv"
    reference = clean / "study.fa.fai"
    _write_fastq(read_1, 1, 3)
    _write_fastq(read_2, 2, 3)
    _write_vcf(variants, "S1", 3)
    _write_manifest(
        manifest,
        sample="S1",
        read_1=read_1,
        read_2=read_2,
        variants=variants,
    )
    reference.write_text(
        "chr1\t1000\t0\t60\t61\nchr2\t800\t1017\t60\t61\n",
        encoding="utf-8",
    )
    inputs = (manifest, read_1, read_2, variants, reference)
    hashes_before = {path.name: _sha256(path) for path in inputs}

    full_arguments = (
        "check",
        str(manifest),
        "--reference",
        str(reference),
        "--full",
        "--format",
        "json",
    )
    full_first = _run(executable, *full_arguments)
    _require(full_first.returncode == 0, "clean full check failed", full_first)
    full_payload = json.loads(full_first.stdout)
    _require(
        full_payload["status"] == "pass"
        and full_payload["summary"]["errors"] == 0,
        "clean cohort did not pass",
        full_first,
    )
    full_second = _run(executable, *full_arguments)
    _require(
        full_second.returncode == 0 and full_second.stdout == full_first.stdout,
        "full JSON output was not deterministic",
        full_second,
    )
    scenarios.append("clean full scan and deterministic JSON")

    bounded = _run(
        executable,
        "check",
        str(manifest),
        "--reference",
        str(reference),
        "--max-records",
        "1",
        "--format",
        "json",
    )
    _require(bounded.returncode == 0, "bounded check failed", bounded)
    bounded_payload = json.loads(bounded.stdout)
    _require(
        all(
            inspection["metrics"].get("scan_complete") is False
            for inspection in bounded_payload["inspections"]
        ),
        "bounded scan incorrectly reported completion",
        bounded,
    )
    scenarios.append("bounded scan semantics")

    shared = workspace / "shared VCF"
    shared.mkdir()
    shared_samples = [f"P{index:03d}" for index in range(1, 26)]
    shared_vcf = shared / "joint.vcf"
    shared_manifest = shared / "cohort.csv"
    _write_joint_vcf(shared_vcf, shared_samples, 2_000)
    with shared_manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ("sample_id", "site", "fastq_1", "fastq_2", "alignment", "variants")
        )
        for sample in shared_samples:
            writer.writerow((sample, "Pilot-Lab", "", "", "", shared_vcf.name))
    shared_result = _run(
        executable,
        "check",
        str(shared_manifest),
        "--full",
        "--format",
        "json",
    )
    _require(shared_result.returncode == 0, "shared VCF check failed", shared_result)
    shared_payload = json.loads(shared_result.stdout)
    shared_inspections = [
        inspection
        for inspection in shared_payload["inspections"]
        if inspection["kind"] == "vcf"
    ]
    _require(
        len(shared_inspections) == 1
        and shared_inspections[0]["metrics"]["records_inspected"] == 2_000,
        "shared VCF was not represented as one full inspection",
        shared_result,
    )
    scenarios.append("shared joint VCF single-pass scan")

    bad = workspace / "seeded failures"
    bad.mkdir()
    bad_r1 = bad / "BAD_R1.fastq.gz"
    bad_r2 = bad / "BAD_R2.fastq.gz"
    bad_vcf = bad / "BAD.vcf"
    bad_manifest = bad / "cohort.csv"
    _write_fastq(bad_r1, 1, 3)
    _write_fastq(bad_r2, 2, 2)
    _write_vcf(bad_vcf, "WRONG_SAMPLE", 2)
    _write_manifest(
        bad_manifest,
        sample="BAD",
        read_1=bad_r1,
        read_2=bad_r2,
        variants=bad_vcf,
    )
    seeded = _run(
        executable,
        "check",
        str(bad_manifest),
        "--full",
        "--format",
        "json",
    )
    _require(seeded.returncode == 1, "seeded failures did not fail", seeded)
    seeded_codes = {item["code"] for item in json.loads(seeded.stdout)["findings"]}
    _require(
        {"FASTQ_PAIR_COUNT_MISMATCH", "HTS_VCF_SAMPLE_MISMATCH"}
        <= seeded_codes,
        "seeded failures were not both detected",
        seeded,
    )
    scenarios.append("seeded FASTQ and VCF failures")

    missing = _run(executable, "check", str(workspace / "missing.csv"))
    _require(
        missing.returncode == 2 and "Traceback" not in missing.stderr,
        "missing-input handling failed",
        missing,
    )
    scenarios.append("invocation error contract")

    original_manifest_hash = _sha256(manifest)
    collision = _run(
        executable,
        "check",
        str(manifest),
        "--output",
        str(manifest),
        "--force",
    )
    _require(
        collision.returncode == 2
        and _sha256(manifest) == original_manifest_hash
        and "refusing to overwrite" in collision.stderr,
        "input/output collision was not blocked",
        collision,
    )

    existing_report = workspace / "existing-report.txt"
    existing_report.write_text("sentinel", encoding="utf-8")
    no_clobber = _run(
        executable,
        "check",
        str(manifest),
        "--output",
        str(existing_report),
    )
    _require(
        no_clobber.returncode == 2
        and existing_report.read_text(encoding="utf-8") == "sentinel",
        "existing report was overwritten without --force",
        no_clobber,
    )
    scenarios.append("no-clobber and input protection")

    lanes = workspace / "multi lane"
    lanes.mkdir()
    for lane in ("L001", "L002"):
        for mate in (1, 2):
            (lanes / f"P1_{lane}_R{mate}_001.fastq.gz").touch()
    discovered_manifest = workspace / "multilane.csv"
    discovered = _run(
        executable,
        "discover",
        str(lanes),
        "--output",
        str(discovered_manifest),
    )
    _require(
        discovered.returncode == 1
        and "DISCOVERY_MULTIPART_FASTQ_UNSUPPORTED" in discovered.stderr,
        "multi-lane limitation was not made explicit",
        discovered,
    )
    scenarios.append("multi-lane safety warning")

    hashes_after = {path.name: _sha256(path) for path in inputs}
    _require(hashes_after == hashes_before, "a normal scan modified an input file")
    scenarios.append("input immutability")

    summary: dict[str, object] = {
        "status": "pass",
        "cohortlint": version.stdout.strip(),
        "executable": str(Path(executable).resolve(strict=False)),
        "host": platform.platform(),
        "python": target_python,
        "python_executable": target_python_executable,
        "harness_python": platform.python_version(),
        "scenario_count": len(scenarios),
        "scenarios": scenarios,
        "elapsed_seconds": round(time.perf_counter() - started, 3),
        "synthetic_input_sha256": hashes_after,
    }
    (workspace / "pilot-summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--executable",
        default="cohortlint",
        help="installed CohortLint executable to test (default: cohortlint)",
    )
    parser.add_argument(
        "--workspace",
        type=Path,
        help="empty directory in which to retain synthetic pilot evidence",
    )
    arguments = parser.parse_args()

    executable = shutil.which(arguments.executable) or arguments.executable
    workspace: Path | None = None
    temporary = False
    try:
        workspace, temporary = _prepare_workspace(arguments.workspace)
        summary = run_pilot(executable, workspace)
    except (OSError, subprocess.TimeoutExpired, PilotFailure, ValueError) as error:
        print(f"CohortLint local pilot: FAIL\n{error}", file=sys.stderr)
        if workspace is not None:
            print(f"Workspace: {workspace}", file=sys.stderr)
        return 1

    print(
        f"CohortLint local pilot: PASS\n"
        f"{summary['scenario_count']} scenarios in {summary['elapsed_seconds']}s\n"
        f"Workspace: {workspace}"
    )
    if temporary:
        shutil.rmtree(workspace)
        print("Temporary workspace removed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

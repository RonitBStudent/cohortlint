from __future__ import annotations

import argparse
import platform
import sys
from pathlib import Path
from typing import Sequence

from . import __version__
from .demo import write_demo
from .engine import check_cohort
from .external import tool_version
from .manifest import discover, write_manifest
from .model import CohortLintError, Severity
from .output import render_json, render_text


MANIFEST_SCHEMA = """CohortLint manifest
===================

One row represents one biological sample. Paths may be absolute or relative to
the manifest. Any genomic columns may be blank, but each sample needs at least
one file.

sample_id,site,fastq_1,fastq_2,alignment,variants
P001,Site-East,reads/P001_R1.fastq.gz,reads/P001_R2.fastq.gz,bam/P001.bam,vcf/cohort.vcf.gz

Supported inputs: FASTQ/FASTQ.GZ, BAM, CRAM (with samtools), VCF/VCF.GZ, and
FASTA .fai reference indexes.
"""


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cohortlint",
        description="Validate genomic cohorts before transfer, merge, or analysis.",
        epilog="CohortLint reports interoperability risk; it does not certify scientific validity.",
    )
    parser.add_argument("--version", action="version", version=f"cohortlint {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)

    check = commands.add_parser("check", help="run cohort-level preflight validation")
    check.add_argument("manifest", type=Path, help="cohort CSV manifest")
    check.add_argument("--reference", type=Path, metavar="FASTA.FAI", help="reference FASTA index used by the study")
    check.add_argument("--full", action="store_true", help="scan complete FASTQ and VCF files instead of bounded sampling")
    check.add_argument("--max-records", type=int, default=10_000, metavar="N", help="records sampled per file (default: 10000)")
    check.add_argument("--format", choices=("text", "json"), default="text")
    check.add_argument("--output", type=Path, help="write the report instead of printing it")
    check.add_argument("--verbose", action="store_true", help="include informational findings")
    check.add_argument("--fail-on", choices=("error", "warning"), default="error", help="finding level that produces exit 1")

    discover_parser = commands.add_parser("discover", help="build a draft manifest from a data directory")
    discover_parser.add_argument("directory", type=Path)
    discover_parser.add_argument("--output", type=Path, default=Path("cohort.csv"))
    discover_parser.add_argument("--no-recursive", action="store_true")
    discover_parser.add_argument("--force", action="store_true", help="replace an existing output manifest")

    demo = commands.add_parser("demo", help="create and inspect a deliberately inconsistent two-site cohort")
    demo.add_argument("--output", type=Path, default=Path("cohortlint-demo"), metavar="DIR")
    demo.add_argument("--format", choices=("text", "json"), default="text")

    commands.add_parser("doctor", help="show runtime and optional HTS tool availability")
    commands.add_parser("schema", help="print the cohort manifest contract")
    return parser


def _write_or_print(text: str, path: Path | None) -> None:
    if path is None:
        sys.stdout.write(text)
        return
    try:
        path.write_text(text, encoding="utf-8")
    except OSError as error:
        raise CohortLintError(f"cannot write {path}: {error.strerror or error}") from error
    print(f"cohortlint: wrote {path}", file=sys.stderr)


def _doctor() -> str:
    lines = [
        "CohortLint doctor",
        f"  Python     {platform.python_version()} ({sys.executable})",
        "  Core       ready; no runtime dependencies",
    ]
    for command in ("samtools", "bcftools", "bgzip", "tabix"):
        version = tool_version(command)
        lines.append(f"  {command.ljust(10)} {version or 'not found (optional)'}")
    lines.extend((
        "",
        "samtools is required only for CRAM. Other external tools are reported so",
        "future remediation commands can be checked before use.",
    ))
    return "\n".join(lines) + "\n"


def run(arguments: argparse.Namespace) -> int:
    if arguments.command == "schema":
        sys.stdout.write(MANIFEST_SCHEMA)
        return 0
    if arguments.command == "doctor":
        sys.stdout.write(_doctor())
        return 0
    if arguments.command == "discover":
        if arguments.output.exists() and not arguments.force:
            raise CohortLintError(
                f"output already exists: {arguments.output}; pass --force to replace it"
            )
        rows, findings = discover(arguments.directory, recursive=not arguments.no_recursive)
        fatal = next(
            (finding for finding in findings if finding.severity == Severity.ERROR),
            None,
        )
        if not rows and fatal is not None:
            raise CohortLintError(f"{fatal.message}: {fatal.path or arguments.directory}")
        write_manifest(rows, arguments.output)
        print(f"cohortlint: wrote {arguments.output} with {len(rows)} sample(s)", file=sys.stderr)
        for finding in findings:
            print(f"{finding.severity.label}: {finding.code}: {finding.message}", file=sys.stderr)
        return 1 if any(finding.severity == Severity.ERROR for finding in findings) else 0
    if arguments.command == "demo":
        existing_manifest = arguments.output / "cohort.csv"
        existing_reference = arguments.output / "reference" / "study.fa.fai"
        if existing_manifest.is_file() and existing_reference.is_file():
            manifest, reference = existing_manifest, existing_reference
            print(f"cohortlint: reusing demonstration cohort in {arguments.output}", file=sys.stderr)
        else:
            manifest, reference = write_demo(arguments.output)
            print(f"cohortlint: wrote demonstration cohort to {arguments.output}", file=sys.stderr)
        report = check_cohort(manifest, reference_fai=reference, max_records=10_000, full=True)
        text = render_json(report) if arguments.format == "json" else render_text(report, verbose=True)
        sys.stdout.write(text)
        return 0


    report = check_cohort(
        arguments.manifest,
        reference_fai=arguments.reference,
        max_records=arguments.max_records,
        full=arguments.full,
    )
    text = render_json(report) if arguments.format == "json" else render_text(report, verbose=arguments.verbose)
    _write_or_print(text, arguments.output)
    if report.errors:
        return 1
    if arguments.fail_on == "warning" and report.warnings:
        return 1
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    try:
        return run(arguments)
    except CohortLintError as error:
        print(f"cohortlint: error: {error}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print("cohortlint: interrupted", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())

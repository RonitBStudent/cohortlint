from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from .model import Finding, Inspection, ReferenceIndex, Severity


def tool_path(name: str) -> str | None:
    return shutil.which(name)


def tool_version(name: str) -> str | None:
    executable = tool_path(name)
    if executable is None:
        return None
    try:
        result = subprocess.run(
            [executable, "--version"],
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return "installed (version unavailable)"
    first_line = (result.stdout or result.stderr).splitlines()
    return first_line[0].strip() if first_line else "installed"


def inspect_cram(path: Path, *, sample_id: str, reference: ReferenceIndex | None = None) -> Inspection:
    findings: list[Finding] = []
    executable = tool_path("samtools")
    if executable is None:
        findings.append(Finding(
            code="CRAM_SAMTOOLS_REQUIRED",
            severity=Severity.ERROR,
            message="CRAM validation requires samtools, but it is not available on PATH",
            sample_id=sample_id,
            path=str(path),
            remediation="Install samtools >=1.15 and rerun CohortLint.",
        ))
        return Inspection(kind="cram", path=str(path), sample_names=(sample_id,), findings=tuple(findings))

    try:
        quick = subprocess.run(
            [executable, "quickcheck", "-v", str(path)],
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        findings.append(Finding(
            code="CRAM_CHECK_FAILED",
            severity=Severity.ERROR,
            message=f"samtools could not inspect CRAM: {error}",
            sample_id=sample_id,
            path=str(path),
        ))
        return Inspection(kind="cram", path=str(path), sample_names=(sample_id,), findings=tuple(findings))
    if quick.returncode:
        findings.append(Finding(
            code="CRAM_INTEGRITY_FAILED",
            severity=Severity.ERROR,
            message="samtools quickcheck reported an invalid or truncated CRAM",
            sample_id=sample_id,
            path=str(path),
            detail=(quick.stdout + quick.stderr).strip(),
            remediation="Re-transfer the CRAM and verify its checksum at the source.",
        ))

    header_command = [executable, "view", "-H"]
    if reference is not None:
        fasta = str(reference.path)
        if fasta.endswith(".fai"):
            fasta = fasta[:-4]
        if Path(fasta).exists():
            header_command.extend(("-T", fasta))
    header_command.append(str(path))
    header = subprocess.run(header_command, capture_output=True, text=True, timeout=30, check=False)
    contigs: list[tuple[str, int | None]] = []
    samples: set[str] = set()
    if header.returncode:
        findings.append(Finding(
            code="CRAM_HEADER_UNREADABLE",
            severity=Severity.ERROR,
            message="samtools could not decode the CRAM header",
            sample_id=sample_id,
            path=str(path),
            detail=header.stderr.strip(),
        ))
    else:
        for line in header.stdout.splitlines():
            fields = dict(field.split(":", 1) for field in line.split("\t")[1:] if ":" in field)
            if line.startswith("@SQ") and "SN" in fields:
                length = int(fields["LN"]) if fields.get("LN", "").isdigit() else None
                contigs.append((fields["SN"], length))
            if line.startswith("@RG") and fields.get("SM"):
                samples.add(fields["SM"])
        if samples and sample_id not in samples:
            findings.append(Finding(
                code="CRAM_SAMPLE_MISMATCH",
                severity=Severity.ERROR,
                message=f"manifest sample {sample_id!r} is absent from CRAM read groups",
                sample_id=sample_id,
                path=str(path),
                detail="CRAM sample names: " + ", ".join(sorted(samples)),
                remediation="Correct the manifest or read-group SM field before merging this cohort.",
            ))

    if not any(candidate.exists() for candidate in (Path(str(path) + ".crai"), path.with_suffix(".crai"))):
        findings.append(Finding(
            code="CRAM_INDEX_MISSING",
            severity=Severity.WARNING,
            message="CRAM index (.crai) was not found",
            sample_id=sample_id,
            path=str(path),
            remediation=f"Run: samtools index {path}",
        ))

    return Inspection(
        kind="cram",
        path=str(path),
        sample_names=tuple(sorted(samples)),
        contigs=tuple(contigs),
        metrics={"samtools": tool_version("samtools")},
        findings=tuple(findings),
    )

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


def _cram_reference_findings(
    contigs: tuple[tuple[str, int], ...],
    reference: ReferenceIndex,
    *,
    sample_id: str,
    path: Path,
) -> list[Finding]:
    """Compare a parsed CRAM sequence dictionary with an explicit reference."""

    findings: list[Finding] = []
    observed_names = [name for name, _ in contigs]
    observed_lengths = dict(contigs)
    reference_names = [name for name, _ in reference.contigs]
    reference_lengths = reference.lengths

    unknown = [name for name in observed_names if name not in reference_lengths]
    if unknown:
        findings.append(Finding(
            code="REF_CONTIG_UNKNOWN",
            severity=Severity.ERROR,
            message="CRAM uses contigs that are absent from the selected reference",
            sample_id=sample_id,
            path=str(path),
            detail=", ".join(unknown[:20]),
            remediation="Use a CRAM generated against the selected reference assembly.",
        ))

    missing = [name for name in reference_names if name not in observed_lengths]
    if missing:
        findings.append(Finding(
            code="REF_CONTIG_MISSING",
            severity=Severity.ERROR,
            message="CRAM sequence dictionary omits reference contigs",
            sample_id=sample_id,
            path=str(path),
            detail=", ".join(missing[:20]),
            remediation="Confirm that the CRAM and reference use the same sequence dictionary.",
        ))

    length_mismatches = [
        f"{name}: file={observed_lengths[name]}, reference={reference_lengths[name]}"
        for name in observed_names
        if name in reference_lengths
        and observed_lengths[name] != reference_lengths[name]
    ]
    if length_mismatches:
        findings.append(Finding(
            code="REF_LENGTH_MISMATCH",
            severity=Severity.ERROR,
            message="CRAM contig lengths do not match the selected reference",
            sample_id=sample_id,
            path=str(path),
            detail="; ".join(length_mismatches[:20]),
            remediation="Regenerate the CRAM with the selected reference assembly.",
        ))

    observed_common = [name for name in observed_names if name in reference_lengths]
    expected_common = [name for name in reference_names if name in observed_lengths]
    if observed_common != expected_common and len(observed_common) > 1:
        findings.append(Finding(
            code="REF_CONTIG_ORDER_MISMATCH",
            severity=Severity.ERROR,
            message="CRAM contig order differs from the selected reference",
            sample_id=sample_id,
            path=str(path),
            detail=(
                f"file starts {observed_common[:5]}; "
                f"reference starts {expected_common[:5]}"
            ),
            remediation="Normalize sequence dictionaries before combining cohort files.",
        ))

    return findings


def _header_fields(line: str) -> tuple[dict[str, str], set[str], bool]:
    """Return SAM header fields, duplicate tags, and malformed-field status."""

    fields: dict[str, str] = {}
    duplicates: set[str] = set()
    malformed = False
    for raw_field in line.split("\t")[1:]:
        if ":" not in raw_field:
            malformed = True
            continue
        tag, value = raw_field.split(":", 1)
        if not tag:
            malformed = True
            continue
        if tag in fields:
            duplicates.add(tag)
            continue
        fields[tag] = value
    return fields, duplicates, malformed


def inspect_cram(path: Path, *, sample_id: str, reference: ReferenceIndex | None = None) -> Inspection:
    path = Path(path)
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
            detail=((quick.stdout or "") + (quick.stderr or "")).strip(),
            remediation="Re-transfer the CRAM and verify its checksum at the source.",
        ))

    header_command = [executable, "view", "-H"]
    if reference is not None:
        fasta = str(reference.path)
        if fasta.endswith(".fai"):
            fasta = fasta[:-4]
        if Path(fasta).is_file():
            header_command.extend(("-T", fasta))
    header_command.append(str(path))
    try:
        header = subprocess.run(
            header_command,
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        header = None
        findings.append(Finding(
            code="CRAM_HEADER_UNREADABLE",
            severity=Severity.ERROR,
            message="samtools could not decode the CRAM header",
            sample_id=sample_id,
            path=str(path),
            detail=f"{type(error).__name__}: {error}",
            remediation="Check samtools, the CRAM, and its reference before retrying.",
        ))

    contigs: list[tuple[str, int]] = []
    samples: set[str] = set()
    if header is not None and header.returncode:
        findings.append(Finding(
            code="CRAM_HEADER_UNREADABLE",
            severity=Severity.ERROR,
            message="samtools could not decode the CRAM header",
            sample_id=sample_id,
            path=str(path),
            detail=((header.stdout or "") + (header.stderr or "")).strip(),
            remediation="Verify that samtools can decode this CRAM with its reference.",
        ))
    elif header is not None:
        malformed_contigs: list[str] = []
        duplicate_contigs: set[str] = set()
        seen_contigs: set[str] = set()
        read_groups_without_sample = 0
        malformed_read_groups = 0

        for line in (header.stdout or "").splitlines():
            record_type = line.split("\t", 1)[0]
            if record_type not in {"@SQ", "@RG"}:
                continue
            fields, duplicate_tags, malformed = _header_fields(line)

            if record_type == "@SQ":
                name = fields.get("SN", "")
                raw_length = fields.get("LN", "")
                detail = ""
                if malformed:
                    detail = "contains a field without a TAG:value separator"
                elif duplicate_tags:
                    detail = "duplicate tag(s): " + ", ".join(sorted(duplicate_tags))
                elif not name:
                    detail = "missing or empty SN tag"
                else:
                    try:
                        length = int(raw_length)
                    except ValueError:
                        detail = f"contig {name!r} has invalid LN {raw_length!r}"
                    else:
                        if length <= 0:
                            detail = f"contig {name!r} has non-positive LN {length}"
                        elif name in seen_contigs:
                            duplicate_contigs.add(name)
                        else:
                            seen_contigs.add(name)
                            contigs.append((name, length))
                if detail and len(malformed_contigs) < 10:
                    malformed_contigs.append(detail)

            if record_type == "@RG":
                if malformed or "SM" in duplicate_tags:
                    malformed_read_groups += 1
                    continue
                sample = fields.get("SM", "")
                if sample:
                    samples.add(sample)
                else:
                    read_groups_without_sample += 1

        if malformed_contigs:
            findings.append(Finding(
                code="CRAM_CONTIG_HEADER_INVALID",
                severity=Severity.ERROR,
                message="CRAM header contains malformed @SQ declarations",
                sample_id=sample_id,
                path=str(path),
                detail="; ".join(malformed_contigs),
                remediation="Repair the CRAM sequence dictionary before cohort analysis.",
            ))
        if duplicate_contigs:
            findings.append(Finding(
                code="CRAM_CONTIG_DUPLICATE",
                severity=Severity.ERROR,
                message="CRAM header declares a contig more than once",
                sample_id=sample_id,
                path=str(path),
                detail=", ".join(sorted(duplicate_contigs)),
                remediation="Regenerate the CRAM with one @SQ entry per contig.",
            ))
        if malformed_read_groups:
            findings.append(Finding(
                code="CRAM_READ_GROUP_INVALID",
                severity=Severity.ERROR,
                message="CRAM header contains malformed @RG declarations",
                sample_id=sample_id,
                path=str(path),
                detail=f"malformed read groups: {malformed_read_groups}",
                remediation="Repair duplicate or malformed read-group tags.",
            ))
        if read_groups_without_sample:
            findings.append(Finding(
                code="CRAM_READ_GROUP_SAMPLE_MISSING",
                severity=Severity.WARNING,
                message="CRAM read groups are missing SM tags",
                sample_id=sample_id,
                path=str(path),
                detail=f"read groups without SM: {read_groups_without_sample}",
                remediation="Populate @RG SM tags before cohort analysis.",
            ))

        if sample_id and not samples:
            findings.append(Finding(
                code="CRAM_SAMPLE_MISSING",
                severity=Severity.ERROR,
                message="CRAM header has no usable @RG SM sample name",
                sample_id=sample_id,
                path=str(path),
                remediation="Add read groups with the manifest sample ID in the SM tag.",
            ))
        elif sample_id and sample_id not in samples:
            findings.append(Finding(
                code="CRAM_SAMPLE_MISMATCH",
                severity=Severity.ERROR,
                message=f"manifest sample {sample_id!r} is absent from CRAM read groups",
                sample_id=sample_id,
                path=str(path),
                detail="CRAM sample names: " + ", ".join(sorted(samples)),
                remediation="Correct the manifest or read-group SM field before merging this cohort.",
            ))
        if len(samples) > 1:
            findings.append(Finding(
                code="CRAM_MULTIPLE_SAMPLES",
                severity=Severity.ERROR,
                message="CRAM contains read groups from multiple sample names",
                sample_id=sample_id,
                path=str(path),
                detail=", ".join(sorted(samples)),
                remediation="Split the CRAM by sample before cohort analysis.",
            ))

        if reference is not None:
            if not contigs:
                findings.append(Finding(
                    code="CRAM_CONTIG_DICTIONARY_MISSING",
                    severity=Severity.ERROR,
                    message="CRAM header has no valid @SQ sequence dictionary",
                    sample_id=sample_id,
                    path=str(path),
                    remediation="Regenerate the CRAM with the selected reference dictionary.",
                ))
            else:
                findings.extend(
                    _cram_reference_findings(
                        tuple(contigs),
                        reference,
                        sample_id=sample_id,
                        path=path,
                    )
                )

    if not any(candidate.is_file() for candidate in (Path(str(path) + ".crai"), path.with_suffix(".crai"))):
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

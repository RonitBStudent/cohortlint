from __future__ import annotations

from pathlib import Path

from .external import inspect_cram
from .fastq import inspect_fastq_pair
from .hts import inspect_bam, inspect_vcf
from .manifest import load_manifest
from .model import CohortLintError, Finding, Inspection, Report, Severity
from .reference import load_fai


def _add_cross_file_reference_findings(inspections: list[Inspection], findings: list[Finding]) -> None:
    dictionaries = [inspection for inspection in inspections if inspection.contigs]
    if len(dictionaries) < 2:
        return
    baseline = dictionaries[0]
    baseline_lengths = dict(baseline.contigs)
    baseline_names = set(baseline_lengths)
    for current in dictionaries[1:]:
        current_lengths = dict(current.contigs)
        current_names = set(current_lengths)
        if current_names != baseline_names:
            stripped_baseline = {name.removeprefix("chr") for name in baseline_names}
            stripped_current = {name.removeprefix("chr") for name in current_names}
            convention_only = stripped_baseline == stripped_current
            findings.append(Finding(
                code="COHORT_CONTIG_SET_MISMATCH",
                severity=Severity.ERROR,
                message=(
                    "files use incompatible contig naming conventions"
                    if convention_only else "files do not describe the same contig set"
                ),
                path=current.path,
                detail=f"Compared with {baseline.path}",
                remediation="Normalize all files against one content-identified reference before cohort merging.",
            ))
            continue
        mismatched = sorted(
            name for name in baseline_names
            if baseline_lengths[name] is not None
            and current_lengths[name] is not None
            and baseline_lengths[name] != current_lengths[name]
        )
        if mismatched:
            findings.append(Finding(
                code="COHORT_CONTIG_LENGTH_MISMATCH",
                severity=Severity.ERROR,
                message=f"{len(mismatched)} contig length(s) differ between genomic files",
                path=current.path,
                detail="First mismatches: " + ", ".join(mismatched[:5]),
                remediation="Confirm the exact reference assembly and patch release used at each site.",
            ))


def check_cohort(
    manifest_path: Path,
    *,
    reference_fai: Path | None = None,
    max_records: int = 10_000,
    full: bool = False,
) -> Report:
    if max_records < 1:
        raise CohortLintError("max_records must be at least 1")
    rows, manifest_findings = load_manifest(manifest_path)
    if not rows:
        fatal = next(
            (finding for finding in manifest_findings if finding.severity == Severity.ERROR),
            None,
        )
        if fatal is not None:
            detail = f" ({fatal.detail})" if fatal.detail else ""
            raise CohortLintError(f"{fatal.message}: {fatal.path or manifest_path}{detail}")
        raise CohortLintError(f"manifest contains no usable samples: {manifest_path}")
    findings = list(manifest_findings)
    reference = load_fai(reference_fai) if reference_fai else None
    inspections: list[Inspection] = []

    sites = {row.site for row in rows if row.site}
    if rows and not sites:
        findings.append(Finding(
            code="METADATA_SITE_UNSPECIFIED",
            severity=Severity.INFO,
            message="manifest does not identify contributing sites",
            remediation="Populate the site column for multi-institution provenance.",
        ))

    for row in rows:
        if not row.paths():
            findings.append(Finding(
                code="MANIFEST_SAMPLE_EMPTY",
                severity=Severity.ERROR,
                message="sample has no genomic files",
                sample_id=row.sample_id,
                remediation="Add FASTQ, alignment, or variant paths for this sample.",
            ))
            continue

        if row.fastq_1 and Path(row.fastq_1).exists():
            second = Path(row.fastq_2) if row.fastq_2 and Path(row.fastq_2).exists() else None
            inspection = inspect_fastq_pair(
                Path(row.fastq_1), second,
                sample_id=row.sample_id,
                max_records=max_records,
                full=full,
            )
            inspections.append(inspection)
            findings.extend(inspection.findings)

        if row.alignment and Path(row.alignment).exists():
            alignment = Path(row.alignment)
            lower = alignment.name.lower()
            alignment_inspection: Inspection | None
            if lower.endswith(".bam"):
                alignment_inspection = inspect_bam(alignment, sample_id=row.sample_id, reference=reference)
            elif lower.endswith(".cram"):
                alignment_inspection = inspect_cram(alignment, sample_id=row.sample_id, reference=reference)
            else:
                findings.append(Finding(
                    code="ALIGNMENT_FORMAT_UNSUPPORTED",
                    severity=Severity.ERROR,
                    message="alignment must be BAM or CRAM",
                    sample_id=row.sample_id,
                    path=str(alignment),
                ))
                alignment_inspection = None
            if alignment_inspection is not None:
                inspections.append(alignment_inspection)
                findings.extend(alignment_inspection.findings)

        if row.variants and Path(row.variants).exists():
            variants = Path(row.variants)
            lower = variants.name.lower()
            if lower.endswith((".vcf", ".vcf.gz")):
                inspection = inspect_vcf(
                    variants,
                    sample_id=row.sample_id,
                    reference=reference,
                    max_records=max_records,
                    full=full,
                )
                inspections.append(inspection)
                findings.extend(inspection.findings)
            else:
                findings.append(Finding(
                    code="VARIANT_FORMAT_UNSUPPORTED",
                    severity=Severity.ERROR,
                    message="variant file must be VCF or VCF.GZ",
                    sample_id=row.sample_id,
                    path=str(variants),
                ))

    _add_cross_file_reference_findings(inspections, findings)
    unique = {(item.code, item.severity, item.message, item.sample_id, item.path, item.detail): item for item in findings}
    ordered = sorted(unique.values(), key=lambda item: (-int(item.severity), item.code, item.sample_id, item.path))
    return Report(
        manifest=str(manifest_path.resolve()),
        sample_count=len(rows),
        site_count=len(sites),
        file_count=len({path for row in rows for _, path in row.paths()}),
        reference_fingerprint=reference.fingerprint if reference else "",
        findings=tuple(ordered),
        inspections=tuple(inspections),
    )

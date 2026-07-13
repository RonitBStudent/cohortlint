from __future__ import annotations

import json
from typing import Any

from .model import Finding, Report, Severity


def _context(finding: Finding) -> str:
    values = []
    if finding.sample_id:
        values.append(f"sample={finding.sample_id}")
    if finding.path:
        values.append(f"file={finding.path}")
    return "  " + "  ".join(values) if values else ""


def render_text(report: Report, *, verbose: bool = False) -> str:
    if report.errors:
        state = "NOT READY"
    elif report.warnings:
        state = "READY WITH WARNINGS"
    else:
        state = "READY"
    lines = [
        f"COHORT {state}",
        "═" * (7 + len(state)),
        f"{report.sample_count} samples  {report.site_count} sites  {report.file_count} files",
    ]
    if report.reference_fingerprint:
        lines.append(f"reference dictionary  {report.reference_fingerprint}")
    lines.append(f"{report.errors} errors  {report.warnings} warnings  {report.infos} notes")

    visible = report.findings if verbose else tuple(item for item in report.findings if item.severity != Severity.INFO)
    if visible:
        lines.append("")
    for finding in visible:
        lines.append(f"{finding.severity.label.upper():7} {finding.code}{_context(finding)}")
        lines.append(f"        {finding.message}")
        if finding.detail:
            lines.append(f"        Detail: {finding.detail}")
        if finding.remediation:
            lines.append(f"        Fix: {finding.remediation}")
        lines.append("")
    if not visible:
        lines.extend(("", "No blocking interoperability problems detected.", ""))
    elif not verbose and report.infos:
        lines.append(f"Use --verbose to show {report.infos} informational finding(s).")
    return "\n".join(lines).rstrip() + "\n"


def _inspection_dict(inspection: Any) -> dict[str, Any]:
    return {
        "kind": inspection.kind,
        "path": inspection.path,
        "sample_names": list(inspection.sample_names),
        "contigs": [{"name": name, "length": length} for name, length in inspection.contigs],
        "metrics": inspection.metrics,
    }


def render_json(report: Report) -> str:
    payload = {
        "schema_version": "1.0",
        "status": "fail" if report.errors else "warn" if report.warnings else "pass",
        "summary": {
            "samples": report.sample_count,
            "sites": report.site_count,
            "files": report.file_count,
            "errors": report.errors,
            "warnings": report.warnings,
            "informational": report.infos,
            "reference_dictionary_fingerprint": report.reference_fingerprint or None,
        },
        "manifest": report.manifest,
        "findings": [finding.as_dict() for finding in report.findings],
        "inspections": [_inspection_dict(inspection) for inspection in report.inspections],
    }
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"

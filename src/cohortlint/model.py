from __future__ import annotations

from dataclasses import asdict, dataclass, field
from enum import IntEnum
from pathlib import Path
from typing import Any


class Severity(IntEnum):
    INFO = 10
    WARNING = 20
    ERROR = 30

    @property
    def label(self) -> str:
        return self.name.lower()


@dataclass(frozen=True, slots=True)
class Finding:
    code: str
    severity: Severity
    message: str
    sample_id: str = ""
    path: str = ""
    detail: str = ""
    remediation: str = ""

    def as_dict(self) -> dict[str, Any]:
        result = asdict(self)
        result["severity"] = self.severity.label
        return result


@dataclass(frozen=True, slots=True)
class ManifestRow:
    sample_id: str
    site: str = ""
    fastq_1: str = ""
    fastq_2: str = ""
    alignment: str = ""
    variants: str = ""

    def paths(self) -> tuple[tuple[str, str], ...]:
        return tuple(
            (role, value)
            for role, value in (
                ("fastq_1", self.fastq_1),
                ("fastq_2", self.fastq_2),
                ("alignment", self.alignment),
                ("variants", self.variants),
            )
            if value
        )


@dataclass(frozen=True, slots=True)
class ReferenceIndex:
    path: Path
    contigs: tuple[tuple[str, int], ...]
    fingerprint: str

    @property
    def lengths(self) -> dict[str, int]:
        return dict(self.contigs)


@dataclass(frozen=True, slots=True)
class Inspection:
    kind: str
    path: str
    sample_names: tuple[str, ...] = ()
    contigs: tuple[tuple[str, int | None], ...] = ()
    metrics: dict[str, Any] = field(default_factory=dict)
    findings: tuple[Finding, ...] = ()


@dataclass(frozen=True, slots=True)
class Report:
    manifest: str
    sample_count: int
    site_count: int
    file_count: int
    reference_fingerprint: str
    findings: tuple[Finding, ...]
    inspections: tuple[Inspection, ...] = ()

    @property
    def errors(self) -> int:
        return sum(finding.severity == Severity.ERROR for finding in self.findings)

    @property
    def warnings(self) -> int:
        return sum(finding.severity == Severity.WARNING for finding in self.findings)

    @property
    def infos(self) -> int:
        return sum(finding.severity == Severity.INFO for finding in self.findings)

    @property
    def passed(self) -> bool:
        return self.errors == 0


class CohortLintError(ValueError):
    """Expected user-facing error without a traceback."""

"""Streaming FASTQ validation for CohortLint.

The checker deliberately uses only the Python standard library.  It validates
enough of each record to catch the input failures that otherwise tend to show
up hours into a sequencing workflow, without retaining reads in memory.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import gzip
from pathlib import Path
import re
from typing import BinaryIO, Protocol
import zlib

from .model import Finding, Inspection, Severity


_DNA_SYMBOLS = frozenset(b"ACGTUNRYSWKMBDHVacgtunryswkmbdhv")
_ILLUMINA_MATE_TAG = re.compile(r"^([12]):")
_MAX_EXAMPLES_PER_ISSUE = 5


class _ReadableStream(Protocol):
    def readline(self, size: int = -1) -> bytes: ...
    def close(self) -> None: ...


@dataclass(slots=True)
class _Record:
    number: int
    read_id: str | None
    mate: int | None


class _FindingCollector:
    """Keep representative findings without growing with the input file."""

    def __init__(self, sample_id: str) -> None:
        self.sample_id = sample_id
        self.findings: list[Finding] = []
        self.counts: Counter[str] = Counter()
        self._shown: Counter[tuple[str, str]] = Counter()
        self._suppressed: dict[tuple[str, str], tuple[Severity, str, int]] = {}

    def add(
        self,
        code: str,
        severity: Severity,
        message: str,
        *,
        path: str = "",
        detail: str = "",
        remediation: str = "",
    ) -> None:
        self.counts[code] += 1
        key = (code, path)
        if self._shown[key] < _MAX_EXAMPLES_PER_ISSUE:
            self.findings.append(
                Finding(
                    code=code,
                    severity=severity,
                    message=message,
                    sample_id=self.sample_id,
                    path=path,
                    detail=detail,
                    remediation=remediation,
                )
            )
            self._shown[key] += 1
        else:
            previous = self._suppressed.get(key)
            omitted = 1 if previous is None else previous[2] + 1
            self._suppressed[key] = (severity, remediation, omitted)

    def finish(self) -> tuple[Finding, ...]:
        for (code, path), (severity, remediation, omitted) in self._suppressed.items():
            self.findings.append(
                Finding(
                    code="FASTQ_FINDINGS_SUPPRESSED",
                    severity=severity,
                    message=f"{omitted} additional {code} finding(s) were omitted",
                    sample_id=self.sample_id,
                    path=path,
                    detail="The complete issue count is available in metrics.issue_counts.",
                    remediation=remediation,
                )
            )
        return tuple(self.findings)


class _FastqReader:
    def __init__(
        self,
        path: Path,
        role: str,
        collector: _FindingCollector,
    ) -> None:
        self.path = path
        self.role = role
        self.collector = collector
        self.raw: BinaryIO | None = None
        self.stream: _ReadableStream | None = None
        self.compressed = path.suffix.lower() in {".gz", ".gzip"}
        self.available = False
        self.failed = False
        self.clean_eof = False
        self.records = 0
        self.bases = 0
        self.gc_bases = 0
        self.min_length: int | None = None
        self.max_length: int | None = None
        self.mate_tag_counts: Counter[int] = Counter()

    @property
    def display_path(self) -> str:
        return str(self.path)

    def open(self) -> None:
        try:
            if not self.path.exists():
                self.collector.add(
                    "FASTQ_FILE_NOT_FOUND",
                    Severity.ERROR,
                    f"{self.role} FASTQ file does not exist",
                    path=self.display_path,
                    remediation="Correct the manifest path or restore the missing FASTQ file.",
                )
                return
            if not self.path.is_file():
                self.collector.add(
                    "FASTQ_NOT_A_FILE",
                    Severity.ERROR,
                    f"{self.role} FASTQ path is not a regular file",
                    path=self.display_path,
                    remediation="Point this manifest field at a FASTQ file, not a directory.",
                )
                return

            self.raw = self.path.open("rb")
            magic = self.raw.read(2)
            self.raw.seek(0)
            self.compressed = self.compressed or magic == b"\x1f\x8b"
            if self.compressed:
                self.stream = gzip.GzipFile(fileobj=self.raw, mode="rb")
            else:
                self.stream = self.raw
            self.available = True
        except OSError as exc:
            self._add_read_error(exc)
            self.close()

    def close(self) -> None:
        # GzipFile does not own a fileobj passed by the caller, so both handles
        # need closing.  Avoid closing the same plain-file handle twice.
        if self.stream is not None and self.stream is not self.raw:
            try:
                self.stream.close()
            except OSError:
                pass
        if self.raw is not None:
            try:
                self.raw.close()
            except OSError:
                pass
        self.stream = None
        self.raw = None

    def read_record(self) -> _Record | None:
        if not self.available or self.clean_eof or self.failed:
            return None
        assert self.stream is not None

        lines: list[bytes] = []
        try:
            first = self.stream.readline()
            if first == b"":
                self.clean_eof = True
                return None
            lines.append(first)
            for _ in range(3):
                line = self.stream.readline()
                if line == b"":
                    missing = 4 - len(lines)
                    self.collector.add(
                        "FASTQ_TRUNCATED_RECORD",
                        Severity.ERROR,
                        f"{self.role} ends partway through FASTQ record {self.records + 1}",
                        path=self.display_path,
                        detail=f"Record has {len(lines)} of 4 lines; {missing} line(s) are missing.",
                        remediation="Re-transfer or regenerate this FASTQ file before analysis.",
                    )
                    # The physical end was reached, but this is not a
                    # successful complete scan: the final record is unusable.
                    self.failed = True
                    self.clean_eof = True
                    return None
                lines.append(line)
        except (gzip.BadGzipFile, EOFError, OSError, zlib.error) as exc:
            self._add_read_error(exc)
            return None

        self.records += 1
        number = self.records
        decoded = [self._decode_line(line, number) for line in lines]
        header, sequence, separator, quality = decoded

        read_id: str | None = None
        mate: int | None = None
        if not header.startswith("@"):
            self.collector.add(
                "FASTQ_INVALID_HEADER",
                Severity.ERROR,
                f"{self.role} record {number} header does not start with '@'",
                path=self.display_path,
                detail=_line_preview(header),
                remediation="Regenerate the FASTQ or correct its four-line record framing.",
            )
        else:
            read_id, mate, conflicting_tags = _normalise_read_id(header)
            if not read_id:
                self.collector.add(
                    "FASTQ_EMPTY_READ_ID",
                    Severity.ERROR,
                    f"{self.role} record {number} has an empty read identifier",
                    path=self.display_path,
                    remediation="Ensure every FASTQ header contains a read identifier after '@'.",
                )
            if conflicting_tags:
                self.collector.add(
                    "FASTQ_CONFLICTING_MATE_TAGS",
                    Severity.ERROR,
                    f"{self.role} record {number} has conflicting mate tags",
                    path=self.display_path,
                    detail=_line_preview(header),
                    remediation="Regenerate read headers so slash and Illumina mate tags agree.",
                )
            if mate is not None:
                self.mate_tag_counts[mate] += 1
                expected = 1 if self.role == "R1" else 2
                if mate != expected:
                    self.collector.add(
                        "FASTQ_UNEXPECTED_MATE_TAG",
                        Severity.ERROR,
                        f"{self.role} record {number} is labelled as mate {mate}",
                        path=self.display_path,
                        detail=_line_preview(header),
                        remediation=f"Put mate {expected} reads in the {self.role} file or repair the read headers.",
                    )

        if not separator.startswith("+"):
            self.collector.add(
                "FASTQ_INVALID_SEPARATOR",
                Severity.ERROR,
                f"{self.role} record {number} separator does not start with '+'",
                path=self.display_path,
                detail=_line_preview(separator),
                remediation="Regenerate the FASTQ or correct its four-line record framing.",
            )

        if not sequence:
            self.collector.add(
                "FASTQ_EMPTY_SEQUENCE",
                Severity.ERROR,
                f"{self.role} record {number} has an empty sequence",
                path=self.display_path,
                remediation="Remove the invalid record or regenerate the FASTQ file.",
            )
        invalid_symbols = sorted(set(sequence.encode("ascii", "replace")) - _DNA_SYMBOLS)
        if invalid_symbols:
            rendered = ", ".join(_render_byte(value) for value in invalid_symbols[:8])
            self.collector.add(
                "FASTQ_INVALID_SEQUENCE_SYMBOL",
                Severity.ERROR,
                f"{self.role} record {number} contains non-IUPAC sequence symbols",
                path=self.display_path,
                detail=f"Invalid symbol(s): {rendered}",
                remediation="Confirm the file is FASTQ and regenerate records containing invalid bases.",
            )

        if len(sequence) != len(quality):
            self.collector.add(
                "FASTQ_SEQUENCE_QUALITY_LENGTH_MISMATCH",
                Severity.ERROR,
                f"{self.role} record {number} sequence and quality lengths differ",
                path=self.display_path,
                detail=f"sequence={len(sequence)}, quality={len(quality)}",
                remediation="Re-transfer or regenerate the FASTQ; do not trim sequence and quality lines independently.",
            )
        invalid_quality = sorted(
            {ord(character) for character in quality if not 33 <= ord(character) <= 126}
        )
        if invalid_quality:
            rendered = ", ".join(_render_byte(value) for value in invalid_quality[:8])
            self.collector.add(
                "FASTQ_INVALID_QUALITY_SYMBOL",
                Severity.ERROR,
                f"{self.role} record {number} contains invalid quality characters",
                path=self.display_path,
                detail=f"Invalid ASCII value(s): {rendered}",
                remediation="Regenerate the FASTQ with printable Phred quality characters.",
            )

        length = len(sequence)
        self.bases += length
        self.gc_bases += sum(base in "GCgc" for base in sequence)
        self.min_length = length if self.min_length is None else min(self.min_length, length)
        self.max_length = length if self.max_length is None else max(self.max_length, length)
        return _Record(number=number, read_id=read_id, mate=mate)

    def _decode_line(self, raw_line: bytes, record: int) -> str:
        # Remove only the record terminator.  Other whitespace remains data and
        # is therefore caught by the sequence/quality validators.
        line = raw_line.removesuffix(b"\n").removesuffix(b"\r")
        try:
            return line.decode("ascii")
        except UnicodeDecodeError as exc:
            self.collector.add(
                "FASTQ_NON_ASCII_DATA",
                Severity.ERROR,
                f"{self.role} record {record} contains non-ASCII bytes",
                path=self.display_path,
                detail=f"First invalid byte is at line offset {exc.start}.",
                remediation="Confirm the input is an unmodified FASTQ file and re-transfer it if necessary.",
            )
            return line.decode("ascii", "replace")

    def _add_read_error(self, exc: BaseException) -> None:
        self.failed = True
        if self.compressed:
            code = "FASTQ_GZIP_INTEGRITY_ERROR"
            message = f"{self.role} gzip stream is truncated or corrupt"
            remediation = "Re-transfer the .gz file and verify its checksum before analysis."
        else:
            code = "FASTQ_READ_ERROR"
            message = f"{self.role} FASTQ file could not be read"
            remediation = "Check file permissions and re-transfer the FASTQ if necessary."
        self.collector.add(
            code,
            Severity.ERROR,
            message,
            path=self.display_path,
            detail=f"{type(exc).__name__}: {exc}",
            remediation=remediation,
        )


def inspect_fastq_pair(
    path_1: Path,
    path_2: Path | None,
    *,
    sample_id: str,
    max_records: int = 10_000,
    full: bool = False,
) -> Inspection:
    """Inspect a single-end FASTQ or a synchronized paired-end FASTQ set.

    In the default sampled mode, at most ``max_records`` records from each file
    are consumed.  ``full=True`` streams through both complete files, which also
    forces gzip trailer/CRC validation and produces exact record counts.
    """

    path_1 = Path(path_1)
    path_2 = Path(path_2) if path_2 is not None else None
    if max_records < 1:
        raise ValueError("max_records must be at least 1")

    collector = _FindingCollector(sample_id)
    first = _FastqReader(path_1, "R1", collector)
    second = _FastqReader(path_2, "R2", collector) if path_2 is not None else None
    readers = (first,) if second is None else (first, second)
    for reader in readers:
        reader.open()

    pair_ids_compared = 0
    pair_id_mismatches = 0
    pair_count_mismatch_observed = False
    limit = None if full else max_records
    iterations = 0

    try:
        if second is not None and first.available and second.available:
            while limit is None or iterations < limit:
                record_1 = first.read_record()
                record_2 = second.read_record()
                if record_1 is None and record_2 is None:
                    break
                iterations += 1
                # A clean EOF on one side while a complete record exists on
                # the other proves unequal mate counts immediately.  Record
                # that fact even when this iteration consumes the sampling
                # budget and the longer file cannot be scanned to its EOF.
                if (
                    record_1 is None
                    and record_2 is not None
                    and first.clean_eof
                    and not first.failed
                ) or (
                    record_2 is None
                    and record_1 is not None
                    and second.clean_eof
                    and not second.failed
                ):
                    pair_count_mismatch_observed = True
                if record_1 is not None and record_2 is not None:
                    if record_1.read_id is not None and record_2.read_id is not None:
                        pair_ids_compared += 1
                        if record_1.read_id != record_2.read_id:
                            pair_id_mismatches += 1
                            collector.add(
                                "FASTQ_PAIR_ID_MISMATCH",
                                Severity.ERROR,
                                f"Paired reads disagree at pair {iterations}",
                                path=f"{path_1} | {path_2}",
                                detail=f"R1={record_1.read_id!r}, R2={record_2.read_id!r}",
                                remediation="Re-pair or re-export the FASTQ files; do not analyze positionally shifted mates.",
                            )
        else:
            # If one paired file cannot be opened, still inspect the readable
            # file so the user gets all actionable failures in one invocation.
            for reader in readers:
                if not reader.available:
                    continue
                scanned = 0
                while limit is None or scanned < limit:
                    if reader.read_record() is None:
                        break
                    scanned += 1
    finally:
        for reader in readers:
            reader.close()

    for reader in readers:
        if (
            reader.available
            and reader.records == 0
            and reader.clean_eof
            and not reader.failed
        ):
            collector.add(
                "FASTQ_EMPTY_FILE",
                Severity.ERROR,
                f"{reader.role} FASTQ contains no complete records",
                path=reader.display_path,
                remediation="Provide a non-empty FASTQ file or remove this sample from the manifest.",
            )

    exact_pair_count_mismatch = (
        second is not None
        and first.available
        and second.available
        and first.clean_eof
        and second.clean_eof
        and not first.failed
        and not second.failed
        and first.records != second.records
    )
    if second is not None and (
        pair_count_mismatch_observed or exact_pair_count_mismatch
    ):
        count_1 = f"{first.records:,}" + (
            "" if first.clean_eof and not first.failed else "+"
        )
        count_2 = f"{second.records:,}" + (
            "" if second.clean_eof and not second.failed else "+"
        )
        collector.add(
            "FASTQ_PAIR_COUNT_MISMATCH",
            Severity.ERROR,
            "Paired FASTQ files contain different numbers of records",
            path=f"{path_1} | {path_2}",
            detail=f"R1={count_1}, R2={count_2} ('+' means at least this many records)",
            remediation="Recover the missing mates or re-pair the files before analysis.",
        )

    scan_complete = all(reader.clean_eof and not reader.failed for reader in readers)
    total_bases = sum(reader.bases for reader in readers)
    total_gc = sum(reader.gc_bases for reader in readers)
    read_lengths = [
        value
        for reader in readers
        for value in (reader.min_length, reader.max_length)
        if value is not None
    ]
    metrics = {
        "sample_id": sample_id,
        "paired": second is not None,
        "full_scan": full,
        "scan_complete": scan_complete,
        "max_records": None if full else max_records,
        "records_scanned": min(reader.records for reader in readers),
        "records_1": first.records,
        "records_2": second.records if second is not None else None,
        "bases_scanned": total_bases,
        "bases_1": first.bases,
        "bases_2": second.bases if second is not None else None,
        "read_length_min": min(read_lengths) if read_lengths else None,
        "read_length_max": max(read_lengths) if read_lengths else None,
        "read_length_mean_1": _mean_length(first),
        "read_length_mean_2": _mean_length(second) if second is not None else None,
        "gc_fraction": round(total_gc / total_bases, 6) if total_bases else None,
        "gzip_1": first.compressed,
        "gzip_2": second.compressed if second is not None else None,
        "pair_ids_compared": pair_ids_compared,
        "pair_id_mismatches": pair_id_mismatches,
        "issue_counts": dict(sorted(collector.counts.items())),
    }
    findings = collector.finish()
    # finish() can add only the synthetic suppression summary; issue_counts is
    # intentionally a count of source issues, rather than display rows.
    return Inspection(
        kind="fastq",
        path=str(path_1),
        sample_names=(sample_id,),
        metrics=metrics,
        findings=findings,
    )


def _normalise_read_id(header: str) -> tuple[str | None, int | None, bool]:
    """Return canonical read ID, optional mate number, and tag conflict flag."""

    if not header.startswith("@"):
        return None, None, False
    fields = header[1:].split()
    if not fields:
        return None, None, False

    identifier = fields[0]
    slash_mate: int | None = None
    if identifier.endswith("/1") or identifier.endswith("/2"):
        slash_mate = int(identifier[-1])
        identifier = identifier[:-2]

    illumina_mate: int | None = None
    if len(fields) > 1:
        match = _ILLUMINA_MATE_TAG.match(fields[1])
        if match:
            illumina_mate = int(match.group(1))

    conflict = (
        slash_mate is not None
        and illumina_mate is not None
        and slash_mate != illumina_mate
    )
    return identifier or None, slash_mate or illumina_mate, conflict


def _mean_length(reader: _FastqReader | None) -> float | None:
    if reader is None or reader.records == 0:
        return None
    return round(reader.bases / reader.records, 3)


def _line_preview(line: str, limit: int = 100) -> str:
    rendered = repr(line)
    return rendered if len(rendered) <= limit else rendered[: limit - 3] + "..."


def _render_byte(value: int) -> str:
    if 33 <= value <= 126:
        return repr(chr(value))
    return str(value)

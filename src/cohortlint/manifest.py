"""Manifest I/O and conservative cohort-file discovery.

The manifest module deliberately does not inspect file contents.  It answers the
earlier, project-level questions: can every file be found, is each file owned by
one sample, and do paired FASTQ names agree?  Content inspection lives in the
format-specific checkers.
"""

from __future__ import annotations

import csv
import io
import re
from collections.abc import Iterable, Sequence
from pathlib import Path

from .model import Finding, ManifestRow, Severity


MANIFEST_FIELDS = (
    "sample_id",
    "site",
    "fastq_1",
    "fastq_2",
    "alignment",
    "variants",
)

_PATH_FIELDS = MANIFEST_FIELDS[2:]
_FASTQ_EXTENSION = re.compile(r"(?i)\.(?:fastq|fq)(?:\.gz)?$")
_FASTQ_MATE = re.compile(
    # Examples: sample_R1, sample_R2_001, sample_1, sample_L001_R1_001.
    r"^(?P<sample>.+?)[_.-](?P<read>R?[12])(?:[_.-](?P<chunk>\d+))?$",
    re.IGNORECASE,
)
_FASTQ_LANE = re.compile(
    r"^(?P<sample>.+?)[_.-]L(?P<lane>\d{3})$",
    re.IGNORECASE,
)


def _finding(
    code: str,
    severity: Severity,
    message: str,
    *,
    sample_id: str = "",
    path: str | Path = "",
    detail: str = "",
    remediation: str = "",
) -> Finding:
    return Finding(
        code=code,
        severity=severity,
        message=message,
        sample_id=sample_id,
        path=str(path),
        detail=detail,
        remediation=remediation,
    )


def _resolve_path(value: str, base: Path) -> str:
    candidate = Path(value).expanduser()
    if not candidate.is_absolute():
        candidate = base / candidate
    return str(candidate.resolve(strict=False))


def _without_fastq_extension(name: str) -> str | None:
    match = _FASTQ_EXTENSION.search(name)
    if match is None:
        return None
    return name[: match.start()]


def _fastq_identity(path: str | Path) -> tuple[str, int | None, str]:
    """Return ``(sample key, mate number, optional chunk)`` for a FASTQ name."""

    stem = _without_fastq_extension(Path(path).name)
    if stem is None:
        return "", None, ""
    match = _FASTQ_MATE.fullmatch(stem)
    if match is None:
        return stem, None, ""
    sample = match.group("sample")
    chunk = match.group("chunk") or ""
    lane_match = _FASTQ_LANE.fullmatch(sample)
    if lane_match is not None:
        sample = lane_match.group("sample")
        lane = f"L{lane_match.group('lane')}"
        chunk = f"{lane}_{chunk}" if chunk else lane
    return sample, int(match.group("read")[-1]), chunk


def _validate_fastq_pair(row: ManifestRow) -> list[Finding]:
    findings: list[Finding] = []
    first = _fastq_identity(row.fastq_1) if row.fastq_1 else None
    second = _fastq_identity(row.fastq_2) if row.fastq_2 else None

    if first and first[1] == 2:
        findings.append(
            _finding(
                "FASTQ_MATE_ROLE",
                Severity.ERROR,
                "fastq_1 has a read-2 filename",
                sample_id=row.sample_id,
                path=row.fastq_1,
                remediation="Move this path to fastq_2 or correct the filename.",
            )
        )
    if second and second[1] == 1:
        findings.append(
            _finding(
                "FASTQ_MATE_ROLE",
                Severity.ERROR,
                "fastq_2 has a read-1 filename",
                sample_id=row.sample_id,
                path=row.fastq_2,
                remediation="Move this path to fastq_1 or correct the filename.",
            )
        )

    if row.fastq_2 and not row.fastq_1:
        findings.append(
            _finding(
                "FASTQ_MATE_1_MISSING",
                Severity.ERROR,
                "fastq_2 is present but fastq_1 is missing",
                sample_id=row.sample_id,
                path=row.fastq_2,
                remediation="Add the matching read-1 file or remove fastq_2.",
            )
        )
        return findings

    # An explicitly named read 1 is almost certainly paired data.  An unmarked
    # FASTQ can legitimately be a single-end library and is therefore accepted.
    if row.fastq_1 and not row.fastq_2 and first and first[1] == 1:
        findings.append(
            _finding(
                "FASTQ_MATE_2_MISSING",
                Severity.WARNING,
                "read-1 FASTQ has no matching read-2 path",
                sample_id=row.sample_id,
                path=row.fastq_1,
                remediation=(
                    "Add the read-2 file, or rename the file if this is a "
                    "single-end library."
                ),
            )
        )

    if not (first and second):
        return findings

    if Path(row.fastq_1) == Path(row.fastq_2):
        findings.append(
            _finding(
                "FASTQ_MATES_IDENTICAL",
                Severity.ERROR,
                "fastq_1 and fastq_2 refer to the same file",
                sample_id=row.sample_id,
                path=row.fastq_1,
                remediation="Provide distinct read-1 and read-2 files.",
            )
        )
        return findings

    if first[1] is None or second[1] is None:
        findings.append(
            _finding(
                "FASTQ_MATE_NAMES_AMBIGUOUS",
                Severity.WARNING,
                "FASTQ mate roles cannot be confirmed from both filenames",
                sample_id=row.sample_id,
                detail=f"fastq_1={Path(row.fastq_1).name}; fastq_2={Path(row.fastq_2).name}",
                remediation="Use a conventional _R1/_R2 or _1/_2 filename pair.",
            )
        )
        return findings

    if first[0].casefold() != second[0].casefold() or first[2] != second[2]:
        findings.append(
            _finding(
                "FASTQ_MATES_MISMATCH",
                Severity.ERROR,
                "FASTQ filenames do not describe the same sample or chunk",
                sample_id=row.sample_id,
                detail=(
                    f"read 1 key={first[0]!r}, chunk={first[2] or '-'}; "
                    f"read 2 key={second[0]!r}, chunk={second[2] or '-'}"
                ),
                remediation="Assign a matching read-1/read-2 pair to this row.",
            )
        )
    return findings


def _validate_rows(rows: Sequence[ManifestRow]) -> list[Finding]:
    findings: list[Finding] = []
    ownership: dict[str, tuple[str, str]] = {}

    for row in rows:
        for role, value in row.paths():
            path = Path(value)
            if not path.exists():
                findings.append(
                    _finding(
                        "PATH_NOT_FOUND",
                        Severity.ERROR,
                        f"{role} path does not exist",
                        sample_id=row.sample_id,
                        path=value,
                        remediation="Correct the manifest path or restore the file.",
                    )
                )
            elif not path.is_file():
                findings.append(
                    _finding(
                        "PATH_NOT_FILE",
                        Severity.ERROR,
                        f"{role} path is not a regular file",
                        sample_id=row.sample_id,
                        path=value,
                        remediation="Point the manifest field at a file.",
                    )
                )

            key = str(path.resolve(strict=False))
            previous = ownership.get(key)
            if previous is not None:
                previous_sample, previous_role = previous
                # A joint-called, multi-sample VCF legitimately appears on
                # several sample rows. Raw reads and alignments must retain
                # unambiguous ownership.
                if not (role == "variants" and previous_role == "variants"):
                    findings.append(
                        _finding(
                            "DUPLICATE_FILE_OWNERSHIP",
                            Severity.ERROR,
                            "one input file is assigned more than once",
                            sample_id=row.sample_id,
                            path=value,
                            detail=(
                                f"first assigned to {previous_sample!r} as {previous_role}; "
                                f"also assigned to {row.sample_id!r} as {role}"
                            ),
                            remediation="Assign raw reads and alignments to exactly one biological sample.",
                        )
                    )
            else:
                ownership[key] = (row.sample_id, role)

        findings.extend(_validate_fastq_pair(row))

    return findings


def load_manifest(
    path: Path,
) -> tuple[tuple[ManifestRow, ...], tuple[Finding, ...]]:
    """Load and validate a CohortLint CSV manifest.

    Relative file references are resolved against the manifest's directory, not
    the process working directory.  Invalid records produce findings; a blank-ID
    record and subsequent occurrences of a duplicate ID are not returned because
    downstream checks cannot own them unambiguously.
    """

    source = Path(path).expanduser().resolve(strict=False)
    if not source.exists():
        finding = _finding(
            "MANIFEST_NOT_FOUND",
            Severity.ERROR,
            "manifest does not exist",
            path=source,
            remediation="Provide an existing CSV manifest path.",
        )
        return (), (finding,)
    if not source.is_file():
        finding = _finding(
            "MANIFEST_NOT_FILE",
            Severity.ERROR,
            "manifest path is not a regular file",
            path=source,
            remediation="Provide a CSV file rather than a directory.",
        )
        return (), (finding,)

    findings: list[Finding] = []
    rows: list[ManifestRow] = []

    try:
        handle = source.open("r", encoding="utf-8-sig", newline="")
    except (OSError, UnicodeError) as exc:
        finding = _finding(
            "MANIFEST_UNREADABLE",
            Severity.ERROR,
            "manifest could not be read",
            path=source,
            detail=str(exc),
        )
        return (), (finding,)

    with handle:
        reader = csv.reader(handle, strict=True)
        try:
            raw_header = next(reader)
        except StopIteration:
            finding = _finding(
                "MANIFEST_EMPTY",
                Severity.ERROR,
                "manifest is empty",
                path=source,
                remediation=f"Add the header: {','.join(MANIFEST_FIELDS)}",
            )
            return (), (finding,)
        except (csv.Error, OSError, UnicodeError) as exc:
            finding = _finding(
                "MANIFEST_CSV_INVALID",
                Severity.ERROR,
                "manifest header is not valid CSV",
                path=source,
                detail=str(exc),
            )
            return (), (finding,)

        header = [name.strip() for name in raw_header]
        positions: dict[str, int] = {}
        duplicate_headers: set[str] = set()
        for index, name in enumerate(header):
            if name in positions:
                duplicate_headers.add(name)
            else:
                positions[name] = index

        missing_headers = [field for field in MANIFEST_FIELDS if field not in positions]
        unexpected_headers = [name for name in header if name and name not in MANIFEST_FIELDS]

        if missing_headers:
            findings.append(
                _finding(
                    "MANIFEST_HEADERS_MISSING",
                    Severity.ERROR,
                    "manifest is missing required columns",
                    path=source,
                    detail=", ".join(missing_headers),
                    remediation=f"Use the header: {','.join(MANIFEST_FIELDS)}",
                )
            )
        if duplicate_headers:
            findings.append(
                _finding(
                    "MANIFEST_HEADERS_DUPLICATE",
                    Severity.ERROR,
                    "manifest contains duplicate column names",
                    path=source,
                    detail=", ".join(sorted(duplicate_headers)),
                    remediation="Keep each schema column exactly once.",
                )
            )
        if unexpected_headers:
            findings.append(
                _finding(
                    "MANIFEST_HEADERS_EXTRA",
                    Severity.WARNING,
                    "manifest contains columns CohortLint will ignore",
                    path=source,
                    detail=", ".join(dict.fromkeys(unexpected_headers)),
                )
            )

        if "sample_id" not in positions:
            return (), tuple(findings)

        seen_ids: dict[str, int] = {}
        try:
            for line_number, cells in enumerate(reader, start=2):
                if not cells or all(not cell.strip() for cell in cells):
                    continue

                if len(cells) > len(header):
                    findings.append(
                        _finding(
                            "MANIFEST_ROW_EXTRA_VALUES",
                            Severity.WARNING,
                            "manifest row has values beyond the header",
                            path=source,
                            detail=f"line {line_number}",
                        )
                    )

                def value(field: str) -> str:
                    position = positions.get(field)
                    if position is None or position >= len(cells):
                        return ""
                    return cells[position].strip()

                sample_id = value("sample_id")
                if not sample_id:
                    findings.append(
                        _finding(
                            "MANIFEST_SAMPLE_ID_BLANK",
                            Severity.ERROR,
                            "manifest row has a blank sample_id",
                            path=source,
                            detail=f"line {line_number}",
                            remediation="Give every row a stable, non-blank sample ID.",
                        )
                    )
                    continue

                identity = sample_id.casefold()
                if identity in seen_ids:
                    findings.append(
                        _finding(
                            "MANIFEST_SAMPLE_ID_DUPLICATE",
                            Severity.ERROR,
                            "sample_id appears more than once",
                            sample_id=sample_id,
                            path=source,
                            detail=(
                                f"line {line_number}; first occurrence on "
                                f"line {seen_ids[identity]}"
                            ),
                            remediation="Merge the rows or assign unique sample IDs.",
                        )
                    )
                    continue
                seen_ids[identity] = line_number

                path_values: dict[str, str] = {}
                for field in _PATH_FIELDS:
                    raw_path = value(field)
                    if not raw_path:
                        path_values[field] = ""
                        continue
                    try:
                        path_values[field] = _resolve_path(raw_path, source.parent)
                    except (OSError, ValueError) as exc:
                        findings.append(
                            _finding(
                                "MANIFEST_PATH_INVALID",
                                Severity.ERROR,
                                f"{field} contains an invalid path",
                                sample_id=sample_id,
                                path=source,
                                detail=f"line {line_number}: {exc}",
                                remediation="Replace the path with a valid local filesystem path.",
                            )
                        )
                        path_values[field] = ""
                rows.append(
                    ManifestRow(
                        sample_id=sample_id,
                        site=value("site"),
                        **path_values,
                    )
                )
        except (csv.Error, OSError, UnicodeError) as exc:
            findings.append(
                _finding(
                    "MANIFEST_CSV_INVALID",
                    Severity.ERROR,
                    "manifest contains invalid CSV",
                    path=source,
                    detail=f"near line {reader.line_num}: {exc}",
                )
            )

    findings.extend(_validate_rows(rows))
    return tuple(rows), tuple(findings)


def _classify_discovered(path: Path) -> tuple[str, str, bool] | None:
    """Return ``(sample_id, role, ambiguous)`` for a supported input file."""

    name = path.name
    lower = name.casefold()

    fastq_stem = _without_fastq_extension(name)
    if fastq_stem is not None:
        sample_id, mate, _chunk = _fastq_identity(path)
        if mate == 2:
            return sample_id, "fastq_2", False
        if mate == 1:
            return sample_id, "fastq_1", False
        # An unmarked file can be a valid single-end library.  The role is clear,
        # but its relationship to similarly named files cannot be inferred.
        return sample_id, "fastq_1", False

    for extension in (".bam", ".cram"):
        if lower.endswith(extension):
            return name[: -len(extension)], "alignment", False

    for extension in (".vcf.gz", ".vcf"):
        if lower.endswith(extension):
            return name[: -len(extension)], "variants", False

    return None


def discover(
    directory: Path, *, recursive: bool = True
) -> tuple[tuple[ManifestRow, ...], tuple[Finding, ...]]:
    """Discover supported cohort files and infer a draft manifest.

    The inference intentionally strips only well-known extensions and FASTQ mate
    suffixes.  Pipeline decorations such as ``.sorted`` or ``.filtered`` remain
    part of the sample ID rather than being guessed away.
    """

    root = Path(directory).expanduser().resolve(strict=False)
    if not root.exists():
        finding = _finding(
            "DISCOVERY_DIRECTORY_NOT_FOUND",
            Severity.ERROR,
            "discovery directory does not exist",
            path=root,
            remediation="Provide an existing directory.",
        )
        return (), (finding,)
    if not root.is_dir():
        finding = _finding(
            "DISCOVERY_NOT_DIRECTORY",
            Severity.ERROR,
            "discovery path is not a directory",
            path=root,
            remediation="Provide a directory containing cohort files.",
        )
        return (), (finding,)

    findings: list[Finding] = []
    grouped: dict[str, tuple[str, dict[str, str]]] = {}
    multipart_fastqs: dict[str, set[str]] = {}
    candidates = root.rglob("*") if recursive else root.iterdir()

    try:
        files = sorted(
            (path for path in candidates if path.is_file()),
            key=lambda path: path.as_posix().casefold(),
        )
    except OSError as exc:
        finding = _finding(
            "DISCOVERY_UNREADABLE",
            Severity.ERROR,
            "discovery directory could not be traversed",
            path=root,
            detail=str(exc),
        )
        return (), (finding,)

    recognized = 0
    for path in files:
        classified = _classify_discovered(path)
        if classified is None:
            continue
        recognized += 1
        sample_id, role, ambiguous = classified
        resolved = str(path.resolve(strict=False))
        sample_id = sample_id.strip()

        if role in {"fastq_1", "fastq_2"}:
            _fastq_sample, _mate, part = _fastq_identity(path)
            if part:
                multipart_fastqs.setdefault(sample_id.casefold(), set()).add(part)

        if not sample_id:
            findings.append(
                _finding(
                    "DISCOVERY_SAMPLE_ID_AMBIGUOUS",
                    Severity.WARNING,
                    "a sample ID cannot be inferred from this filename",
                    path=resolved,
                    remediation="Rename the file or add it to a manifest manually.",
                )
            )
            continue
        if ambiguous:
            findings.append(
                _finding(
                    "DISCOVERY_SAMPLE_ID_AMBIGUOUS",
                    Severity.WARNING,
                    "sample ID inference is ambiguous",
                    sample_id=sample_id,
                    path=resolved,
                    remediation="Review the generated manifest before analysis.",
                )
            )

        identity = sample_id.casefold()
        canonical_id, roles = grouped.setdefault(identity, (sample_id, {}))
        previous = roles.get(role)
        if previous is not None:
            findings.append(
                _finding(
                    "DISCOVERY_ROLE_AMBIGUOUS",
                    Severity.WARNING,
                    f"multiple files could fill the {role} field",
                    sample_id=canonical_id,
                    path=resolved,
                    detail=f"already selected {previous}",
                    remediation=(
                        "Choose the intended file manually; only the first sorted "
                        "path is included in the draft manifest."
                    ),
                )
            )
            continue
        roles[role] = resolved

    rows = tuple(
        ManifestRow(sample_id=canonical_id, **roles)
        for canonical_id, roles in sorted(
            grouped.values(), key=lambda item: item[0].casefold()
        )
    )

    for identity, parts in sorted(multipart_fastqs.items()):
        if len(parts) < 2:
            continue
        canonical_id = grouped.get(identity, (identity, {}))[0]
        findings.append(
            _finding(
                "DISCOVERY_MULTIPART_FASTQ_UNSUPPORTED",
                Severity.ERROR,
                "multiple FASTQ lanes or chunks cannot be represented in one manifest row",
                sample_id=canonical_id,
                detail="detected parts: " + ", ".join(sorted(parts)),
                remediation=(
                    "Merge technical FASTQ parts upstream before using the generated manifest; "
                    "do not treat lanes as separate biological samples."
                ),
            )
        )

    if recognized == 0:
        findings.append(
            _finding(
                "DISCOVERY_NO_SUPPORTED_FILES",
                Severity.WARNING,
                "no supported FASTQ, BAM/CRAM, or VCF files were found",
                path=root,
            )
        )

    findings.extend(_validate_rows(rows))
    return rows, tuple(findings)


def write_manifest(rows: Iterable[ManifestRow], path: Path) -> None:
    """Write rows using CohortLint's stable CSV column order."""

    destination = Path(path).expanduser()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(render_manifest(rows), encoding="utf-8", newline="")


def render_manifest(rows: Iterable[ManifestRow]) -> str:
    """Render rows using CohortLint's stable CSV column order."""

    output = io.StringIO(newline="")
    writer = csv.writer(output)
    writer.writerow(MANIFEST_FIELDS)
    for row in rows:
        writer.writerow(
            (
                row.sample_id,
                row.site,
                row.fastq_1,
                row.fastq_2,
                row.alignment,
                row.variants,
            )
        )
    return output.getvalue()

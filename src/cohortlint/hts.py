from __future__ import annotations

import gzip
import re
import struct
from pathlib import Path
from typing import Protocol

from .model import Finding, Inspection, ReferenceIndex, Severity


_VALID_REF = re.compile(r"^[ACGTN]+$", re.IGNORECASE)
_BAM_MAGIC = b"BAM\x01"
_MAX_BAM_HEADER_BYTES = 64 * 1024 * 1024
_MAX_BAM_REFERENCES = 1_000_000
_MAX_BAM_REFERENCE_NAME_BYTES = 1024 * 1024
_BGZF_EOF = bytes.fromhex(
    "1f8b08040000000000ff0600424302001b0003000000000000000000"
)


class _BinaryReader(Protocol):
    def read(self, size: int = -1) -> bytes: ...


def _finding(
    code: str,
    severity: Severity,
    message: str,
    *,
    sample_id: str,
    path: Path,
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


def _has_gzip_magic(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def _index_candidates(path: Path, suffixes: tuple[str, ...]) -> tuple[Path, ...]:
    candidates = [Path(f"{path}{suffix}") for suffix in suffixes]
    if path.suffix:
        candidates.extend(path.with_suffix(suffix) for suffix in suffixes)
    # Preserve order while removing duplicates (e.g. a suffix-less path).
    return tuple(dict.fromkeys(candidates))


def _split_structured_fields(value: str) -> dict[str, str]:
    """Parse the comma-delimited key/value body of a VCF structured header."""

    fields: list[str] = []
    current: list[str] = []
    quoted = False
    escaped = False
    for char in value:
        if escaped:
            current.append(char)
            escaped = False
        elif char == "\\" and quoted:
            current.append(char)
            escaped = True
        elif char == '"':
            current.append(char)
            quoted = not quoted
        elif char == "," and not quoted:
            fields.append("".join(current))
            current = []
        else:
            current.append(char)
    fields.append("".join(current))

    result: dict[str, str] = {}
    for field in fields:
        if "=" not in field:
            continue
        key, raw_value = field.split("=", 1)
        raw_value = raw_value.strip()
        if len(raw_value) >= 2 and raw_value[0] == raw_value[-1] == '"':
            raw_value = raw_value[1:-1].replace(r'\"', '"').replace(r"\\", "\\")
        result[key.strip()] = raw_value
    return result


def _parse_vcf_contig(line: str) -> tuple[str, int | None]:
    prefix = "##contig=<"
    if not line.startswith(prefix) or not line.endswith(">"):
        raise ValueError("malformed ##contig declaration")
    attributes = _split_structured_fields(line[len(prefix) : -1])
    name = attributes.get("ID", "")
    if not name:
        raise ValueError("contig declaration has no ID")

    raw_length = attributes.get("length", attributes.get("Length"))
    if raw_length is None:
        return name, None
    try:
        length = int(raw_length)
    except ValueError as exc:
        raise ValueError(f"contig {name!r} has a non-integer length") from exc
    if length <= 0:
        raise ValueError(f"contig {name!r} has non-positive length {length}")
    return name, length


def _dictionary_findings(
    observed: tuple[tuple[str, int | None], ...],
    reference: ReferenceIndex,
    *,
    sample_id: str,
    path: Path,
    source: str,
    require_all_contigs: bool,
) -> list[Finding]:
    findings: list[Finding] = []
    observed_names = [name for name, _ in observed]
    observed_lengths = dict(observed)
    reference_names = [name for name, _ in reference.contigs]
    reference_lengths = reference.lengths

    unknown = [name for name in observed_names if name not in reference_lengths]
    if unknown:
        findings.append(
            _finding(
                "REF_CONTIG_UNKNOWN",
                Severity.ERROR,
                f"{source} uses contigs that are absent from the selected reference",
                sample_id=sample_id,
                path=path,
                detail=", ".join(unknown[:20]),
                remediation="Use files generated against the same reference assembly.",
            )
        )

    missing = [name for name in reference_names if name not in observed_lengths]
    if missing and observed:
        findings.append(
            _finding(
                "REF_CONTIG_MISSING",
                Severity.ERROR if require_all_contigs else Severity.WARNING,
                f"{source} sequence dictionary omits reference contigs",
                sample_id=sample_id,
                path=path,
                detail=", ".join(missing[:20]),
                remediation="Confirm that the file and reference use the same sequence dictionary.",
            )
        )

    length_mismatches: list[str] = []
    for name in observed_names:
        observed_length = observed_lengths[name]
        if (
            observed_length is not None
            and name in reference_lengths
            and observed_length != reference_lengths[name]
        ):
            length_mismatches.append(
                f"{name}: file={observed_length}, reference={reference_lengths[name]}"
            )
    if length_mismatches:
        findings.append(
            _finding(
                "REF_LENGTH_MISMATCH",
                Severity.ERROR,
                f"{source} contig lengths do not match the selected reference",
                sample_id=sample_id,
                path=path,
                detail="; ".join(length_mismatches[:20]),
                remediation="Regenerate the file with the selected reference assembly.",
            )
        )

    observed_common = [name for name in observed_names if name in reference_lengths]
    expected_common = [name for name in reference_names if name in observed_lengths]
    if observed_common != expected_common and len(observed_common) > 1:
        findings.append(
            _finding(
                "REF_CONTIG_ORDER_MISMATCH",
                Severity.ERROR if require_all_contigs else Severity.WARNING,
                f"{source} contig order differs from the selected reference",
                sample_id=sample_id,
                path=path,
                detail=f"file starts {observed_common[:5]}; reference starts {expected_common[:5]}",
                remediation="Normalize sequence dictionaries before combining cohort files.",
            )
        )

    return findings


def inspect_vcf(
    path: Path,
    *,
    sample_id: str,
    reference: ReferenceIndex | None = None,
    max_records: int = 10_000,
    full: bool = False,
) -> Inspection:
    """Inspect a plain-text or gzip-compressed VCF without external tools.

    By default, record-level checks stop after ``max_records`` variants while
    header and reference-dictionary checks always run.  ``full=True`` scans all
    records.
    """

    path = Path(path)
    findings: list[Finding] = []
    sample_names: tuple[str, ...] = ()
    declared_contigs: list[tuple[str, int | None]] = []
    observed_contigs: list[str] = []
    observed_set: set[str] = set()
    records_inspected = 0
    truncated = False
    compressed = False
    indexed = False

    if not full and (not isinstance(max_records, int) or isinstance(max_records, bool) or max_records <= 0):
        findings.append(
            _finding(
                "HTS_VCF_LIMIT_INVALID",
                Severity.ERROR,
                "VCF record limit must be a positive integer",
                sample_id=sample_id,
                path=path,
                detail=f"max_records={max_records!r}",
            )
        )
        return Inspection(kind="vcf", path=str(path), findings=tuple(findings))

    try:
        compressed = _has_gzip_magic(path)
    except OSError as exc:
        findings.append(
            _finding(
                "HTS_VCF_UNREADABLE",
                Severity.ERROR,
                "VCF could not be opened",
                sample_id=sample_id,
                path=path,
                detail=str(exc),
                remediation="Check the manifest path and file permissions.",
            )
        )
        return Inspection(kind="vcf", path=str(path), findings=tuple(findings))

    if compressed:
        index_paths = _index_candidates(path, (".tbi", ".csi"))
        indexed = any(candidate.is_file() for candidate in index_paths)
        if not indexed:
            findings.append(
                _finding(
                    "HTS_VCF_INDEX_MISSING",
                    Severity.WARNING,
                    "Compressed VCF has no adjacent tabix or CSI index",
                    sample_id=sample_id,
                    path=path,
                    remediation=f"Create {path.name}.tbi or {path.name}.csi before random-access analysis.",
                )
            )

    header_seen = False
    header_columns: tuple[str, ...] = ()
    duplicate_contigs: set[str] = set()
    malformed_contig_headers: list[str] = []
    malformed_records: list[str] = []
    invalid_refs: list[str] = []
    unsorted_records: list[str] = []
    undeclared_record_contigs: set[str] = set()
    unknown_reference_contigs: set[str] = set()
    out_of_range_records: list[str] = []
    sample_column_mismatches: list[str] = []
    last_contig: str | None = None
    last_position: int | None = None
    closed_contigs: set[str] = set()
    declared_rank: dict[str, int] = {}
    reference_rank = (
        {name: index for index, (name, _) in enumerate(reference.contigs)}
        if reference is not None
        else {}
    )

    try:
        if compressed:
            handle = gzip.open(path, "rt", encoding="utf-8", newline="")
        else:
            handle = path.open("rt", encoding="utf-8", newline="")

        with handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.rstrip("\r\n")
                if line_number == 1:
                    line = line.lstrip("\ufeff")

                if line.startswith("##") and not header_seen:
                    if line.startswith("##contig="):
                        try:
                            contig = _parse_vcf_contig(line)
                        except ValueError as exc:
                            if len(malformed_contig_headers) < 10:
                                malformed_contig_headers.append(f"line {line_number}: {exc}")
                        else:
                            name = contig[0]
                            if name in declared_rank:
                                duplicate_contigs.add(name)
                            else:
                                declared_rank[name] = len(declared_contigs)
                                declared_contigs.append(contig)
                    continue

                if line.startswith("#CHROM"):
                    if header_seen:
                        findings.append(
                            _finding(
                                "HTS_VCF_HEADER_DUPLICATE",
                                Severity.ERROR,
                                "VCF contains more than one #CHROM header",
                                sample_id=sample_id,
                                path=path,
                                detail=f"second header at line {line_number}",
                            )
                        )
                        continue
                    header_seen = True
                    header_columns = tuple(line.split("\t"))
                    if len(header_columns) < 8:
                        findings.append(
                            _finding(
                                "HTS_VCF_HEADER_INVALID",
                                Severity.ERROR,
                                "VCF #CHROM header has fewer than eight columns",
                                sample_id=sample_id,
                                path=path,
                                detail=f"line {line_number}",
                            )
                        )
                    sample_names = header_columns[9:] if len(header_columns) > 9 else ()
                    continue

                if line.startswith("#") or not line:
                    continue

                if not full and records_inspected >= max_records:
                    truncated = True
                    break

                records_inspected += 1
                fields = line.split("\t")
                if len(fields) < 8:
                    if len(malformed_records) < 10:
                        malformed_records.append(
                            f"line {line_number}: expected at least 8 columns, found {len(fields)}"
                        )
                    continue

                chrom, raw_position, _record_id, ref = fields[:4]
                try:
                    position = int(raw_position)
                    if position <= 0:
                        raise ValueError
                except ValueError:
                    if len(malformed_records) < 10:
                        malformed_records.append(
                            f"line {line_number}: invalid POS {raw_position!r}"
                        )
                    continue

                if header_seen and len(fields) != len(header_columns):
                    if len(sample_column_mismatches) < 10:
                        sample_column_mismatches.append(
                            f"line {line_number}: header={len(header_columns)}, record={len(fields)}"
                        )

                if chrom not in observed_set:
                    observed_set.add(chrom)
                    observed_contigs.append(chrom)
                if declared_contigs and chrom not in declared_rank:
                    undeclared_record_contigs.add(chrom)

                if reference is not None:
                    if chrom not in reference.lengths:
                        unknown_reference_contigs.add(chrom)
                    elif position > reference.lengths[chrom]:
                        if len(out_of_range_records) < 10:
                            out_of_range_records.append(
                                f"line {line_number}: {chrom}:{position} > {reference.lengths[chrom]}"
                            )

                if not _VALID_REF.fullmatch(ref):
                    if len(invalid_refs) < 10:
                        invalid_refs.append(f"line {line_number}: {chrom}:{position} REF={ref!r}")

                if last_contig is not None:
                    if chrom == last_contig:
                        if last_position is not None and position < last_position:
                            if len(unsorted_records) < 10:
                                unsorted_records.append(
                                    f"line {line_number}: {chrom}:{position} follows {chrom}:{last_position}"
                                )
                    else:
                        order = reference_rank or declared_rank
                        if (
                            chrom in closed_contigs
                            or (
                                last_contig in order
                                and chrom in order
                                and order[chrom] < order[last_contig]
                            )
                        ):
                            if len(unsorted_records) < 10:
                                unsorted_records.append(
                                    f"line {line_number}: contig {chrom!r} appears after {last_contig!r}"
                                )
                        closed_contigs.add(last_contig)
                last_contig = chrom
                last_position = position
    except (OSError, EOFError, UnicodeError, gzip.BadGzipFile) as exc:
        findings.append(
            _finding(
                "HTS_VCF_READ_ERROR",
                Severity.ERROR,
                "VCF is truncated, corrupt, or not valid text",
                sample_id=sample_id,
                path=path,
                detail=str(exc),
                remediation="Re-transfer or regenerate the VCF and its index.",
            )
        )

    if not header_seen:
        findings.append(
            _finding(
                "HTS_VCF_HEADER_MISSING",
                Severity.ERROR,
                "VCF has no #CHROM header",
                sample_id=sample_id,
                path=path,
                remediation="Provide a standards-compliant VCF with a complete header.",
            )
        )
    else:
        if len(set(sample_names)) != len(sample_names):
            findings.append(
                _finding(
                    "HTS_VCF_SAMPLE_DUPLICATE",
                    Severity.ERROR,
                    "VCF contains duplicate sample columns",
                    sample_id=sample_id,
                    path=path,
                    detail=", ".join(sample_names),
                )
            )
        if sample_id and sample_id not in sample_names:
            findings.append(
                _finding(
                    "HTS_VCF_SAMPLE_MISMATCH",
                    Severity.ERROR,
                    "Manifest sample is absent from the VCF genotype columns",
                    sample_id=sample_id,
                    path=path,
                    detail=(
                        f"expected {sample_id!r}; VCF samples: {', '.join(sample_names) or '(none)'}"
                    ),
                    remediation="Correct the manifest sample ID or use the matching VCF.",
                )
            )

    if malformed_contig_headers:
        findings.append(
            _finding(
                "HTS_VCF_CONTIG_HEADER_INVALID",
                Severity.ERROR,
                "VCF contains malformed contig declarations",
                sample_id=sample_id,
                path=path,
                detail="; ".join(malformed_contig_headers),
            )
        )
    if duplicate_contigs:
        findings.append(
            _finding(
                "HTS_VCF_CONTIG_DUPLICATE",
                Severity.ERROR,
                "VCF declares a contig more than once",
                sample_id=sample_id,
                path=path,
                detail=", ".join(sorted(duplicate_contigs)),
            )
        )
    if malformed_records:
        findings.append(
            _finding(
                "HTS_VCF_RECORD_INVALID",
                Severity.ERROR,
                "VCF contains malformed variant records",
                sample_id=sample_id,
                path=path,
                detail="; ".join(malformed_records),
            )
        )
    if sample_column_mismatches:
        findings.append(
            _finding(
                "HTS_VCF_COLUMN_MISMATCH",
                Severity.ERROR,
                "VCF records do not match the header column count",
                sample_id=sample_id,
                path=path,
                detail="; ".join(sample_column_mismatches),
            )
        )
    if invalid_refs:
        findings.append(
            _finding(
                "HTS_VCF_REF_INVALID",
                Severity.ERROR,
                "VCF REF alleles contain invalid bases",
                sample_id=sample_id,
                path=path,
                detail="; ".join(invalid_refs),
                remediation="Normalize or regenerate the affected VCF records.",
            )
        )
    if unsorted_records:
        findings.append(
            _finding(
                "HTS_VCF_UNSORTED",
                Severity.ERROR,
                "VCF records are not coordinate sorted",
                sample_id=sample_id,
                path=path,
                detail="; ".join(unsorted_records),
                remediation="Sort the VCF in reference contig order and rebuild its index.",
            )
        )
    if undeclared_record_contigs:
        findings.append(
            _finding(
                "HTS_VCF_CONTIG_UNDECLARED",
                Severity.WARNING,
                "VCF records use contigs absent from its header dictionary",
                sample_id=sample_id,
                path=path,
                detail=", ".join(sorted(undeclared_record_contigs)),
            )
        )
    if out_of_range_records:
        findings.append(
            _finding(
                "REF_POSITION_OUT_OF_RANGE",
                Severity.ERROR,
                "VCF positions exceed reference contig lengths",
                sample_id=sample_id,
                path=path,
                detail="; ".join(out_of_range_records),
            )
        )

    inspection_contigs: tuple[tuple[str, int | None], ...]
    if declared_contigs:
        inspection_contigs = tuple(declared_contigs)
    else:
        inspection_contigs = tuple((name, None) for name in observed_contigs)

    if reference is not None:
        findings.extend(
            _dictionary_findings(
                tuple(declared_contigs),
                reference,
                sample_id=sample_id,
                path=path,
                source="VCF",
                require_all_contigs=False,
            )
        )
        declared_unknown = {
            name for name, _ in declared_contigs if name not in reference.lengths
        }
        unknown_reference_contigs -= declared_unknown
        if unknown_reference_contigs:
            findings.append(
                _finding(
                    "REF_CONTIG_UNKNOWN",
                    Severity.ERROR,
                    "VCF records use contigs absent from the selected reference",
                    sample_id=sample_id,
                    path=path,
                    detail=", ".join(sorted(unknown_reference_contigs)),
                    remediation="Use a VCF generated against the selected reference assembly.",
                )
            )

    metrics = {
        "records_inspected": records_inspected,
        "scan_complete": not truncated,
        "compressed": compressed,
        "indexed": indexed,
    }
    return Inspection(
        kind="vcf",
        path=str(path),
        sample_names=sample_names,
        contigs=inspection_contigs,
        metrics=metrics,
        findings=tuple(findings),
    )


def _read_exact(handle: _BinaryReader, size: int, description: str) -> bytes:
    data = handle.read(size)
    if len(data) != size:
        raise ValueError(f"truncated while reading {description}")
    return data


def _read_i32(handle: _BinaryReader, description: str) -> int:
    return struct.unpack("<i", _read_exact(handle, 4, description))[0]


def _has_bgzf_extra_field(path: Path) -> bool:
    """Return whether the first gzip member contains the BGZF ``BC`` field."""

    with path.open("rb") as handle:
        fixed = handle.read(12)
        if len(fixed) < 12 or fixed[:3] != b"\x1f\x8b\x08":
            return False
        flags = fixed[3]
        if not flags & 0x04:  # FEXTRA
            return False
        extra_length = struct.unpack("<H", fixed[10:12])[0]
        extra = handle.read(extra_length)
        if len(extra) != extra_length:
            return False

    offset = 0
    while offset + 4 <= len(extra):
        subfield_id = extra[offset : offset + 2]
        subfield_length = struct.unpack("<H", extra[offset + 2 : offset + 4])[0]
        offset += 4
        if offset + subfield_length > len(extra):
            return False
        if subfield_id == b"BC" and subfield_length == 2:
            return True
        offset += subfield_length
    return False


def _has_bgzf_eof(path: Path) -> bool:
    try:
        if path.stat().st_size < len(_BGZF_EOF):
            return False
        with path.open("rb") as handle:
            handle.seek(-len(_BGZF_EOF), 2)
            return handle.read(len(_BGZF_EOF)) == _BGZF_EOF
    except OSError:
        return False


def _bam_sample_names(header_text: str) -> tuple[tuple[str, ...], int]:
    names: list[str] = []
    read_groups_without_sample = 0
    for line in header_text.splitlines():
        if not line.startswith("@RG\t"):
            continue
        fields = line.split("\t")[1:]
        sample_values = [field[3:] for field in fields if field.startswith("SM:")]
        nonempty_values = [value for value in sample_values if value]
        if not nonempty_values:
            read_groups_without_sample += 1
            continue
        for value in nonempty_values:
            if value not in names:
                names.append(value)
    return tuple(names), read_groups_without_sample


def inspect_bam(
    path: Path,
    *,
    sample_id: str,
    reference: ReferenceIndex | None = None,
) -> Inspection:
    """Inspect BAM header identity, reference dictionary, EOF marker, and index."""

    path = Path(path)
    findings: list[Finding] = []
    sample_names: tuple[str, ...] = ()
    contigs: tuple[tuple[str, int | None], ...] = ()
    header_bytes = 0
    read_groups_without_sample = 0

    try:
        is_gzip = _has_gzip_magic(path)
    except OSError as exc:
        findings.append(
            _finding(
                "HTS_BAM_UNREADABLE",
                Severity.ERROR,
                "BAM could not be opened",
                sample_id=sample_id,
                path=path,
                detail=str(exc),
                remediation="Check the manifest path and file permissions.",
            )
        )
        return Inspection(kind="bam", path=str(path), findings=tuple(findings))

    is_bgzf = False
    if is_gzip:
        try:
            is_bgzf = _has_bgzf_extra_field(path)
        except OSError:
            is_bgzf = False
    if not is_gzip:
        findings.append(
            _finding(
                "HTS_BAM_COMPRESSION_INVALID",
                Severity.ERROR,
                "BAM does not have a gzip/BGZF header",
                sample_id=sample_id,
                path=path,
                remediation="Provide a BGZF-compressed BAM file.",
            )
        )
    elif not is_bgzf:
        findings.append(
            _finding(
                "HTS_BAM_NOT_BGZF",
                Severity.ERROR,
                "BAM gzip stream has no BGZF BC extra field",
                sample_id=sample_id,
                path=path,
                remediation="Re-encode the alignment as BGZF-compressed BAM.",
            )
        )

    has_eof = _has_bgzf_eof(path)
    if not has_eof:
        findings.append(
            _finding(
                "HTS_BAM_EOF_MISSING",
                Severity.WARNING,
                "BAM is missing the standard 28-byte BGZF EOF marker",
                sample_id=sample_id,
                path=path,
                remediation="Run an integrity check and regenerate or re-transfer the BAM if truncated.",
            )
        )

    index_paths = _index_candidates(path, (".bai", ".csi"))
    indexed = any(candidate.is_file() for candidate in index_paths)
    if not indexed:
        findings.append(
            _finding(
                "HTS_BAM_INDEX_MISSING",
                Severity.WARNING,
                "BAM has no adjacent BAI or CSI index",
                sample_id=sample_id,
                path=path,
                remediation=f"Create {path.name}.bai or {path.name}.csi before random-access analysis.",
            )
        )

    if is_gzip:
        try:
            with gzip.open(path, "rb") as handle:
                magic = _read_exact(handle, 4, "BAM magic")
                if magic != _BAM_MAGIC:
                    raise ValueError("decompressed stream does not begin with BAM\\x01")

                header_bytes = _read_i32(handle, "BAM header length")
                if not 0 <= header_bytes <= _MAX_BAM_HEADER_BYTES:
                    raise ValueError(f"implausible BAM header length {header_bytes}")
                raw_header = _read_exact(handle, header_bytes, "BAM header text")
                try:
                    header_text = raw_header.decode("utf-8")
                except UnicodeDecodeError as exc:
                    raise ValueError("BAM header text is not valid UTF-8") from exc

                reference_count = _read_i32(handle, "BAM reference count")
                if not 0 <= reference_count <= _MAX_BAM_REFERENCES:
                    raise ValueError(f"implausible BAM reference count {reference_count}")

                parsed_contigs: list[tuple[str, int | None]] = []
                seen_contigs: set[str] = set()
                for index in range(reference_count):
                    name_length = _read_i32(handle, f"reference {index + 1} name length")
                    if not 1 <= name_length <= _MAX_BAM_REFERENCE_NAME_BYTES:
                        raise ValueError(
                            f"invalid BAM reference name length {name_length} at index {index}"
                        )
                    raw_name = _read_exact(
                        handle, name_length, f"reference {index + 1} name"
                    )
                    if not raw_name.endswith(b"\x00"):
                        raise ValueError(f"BAM reference name at index {index} is not NUL terminated")
                    try:
                        name = raw_name[:-1].decode("utf-8")
                    except UnicodeDecodeError as exc:
                        raise ValueError(
                            f"BAM reference name at index {index} is not valid UTF-8"
                        ) from exc
                    if not name or "\x00" in name:
                        raise ValueError(f"invalid BAM reference name at index {index}")
                    if name in seen_contigs:
                        raise ValueError(f"duplicate BAM reference name {name!r}")
                    length = _read_i32(handle, f"reference {name!r} length")
                    if length <= 0:
                        raise ValueError(f"BAM reference {name!r} has invalid length {length}")
                    seen_contigs.add(name)
                    parsed_contigs.append((name, length))

                contigs = tuple(parsed_contigs)
                sample_names, read_groups_without_sample = _bam_sample_names(header_text)
        except (OSError, EOFError, ValueError, struct.error, gzip.BadGzipFile) as exc:
            findings.append(
                _finding(
                    "HTS_BAM_READ_ERROR",
                    Severity.ERROR,
                    "BAM header is truncated or malformed",
                    sample_id=sample_id,
                    path=path,
                    detail=str(exc),
                    remediation="Run samtools quickcheck and regenerate or re-transfer the BAM.",
                )
            )

    if read_groups_without_sample:
        findings.append(
            _finding(
                "HTS_BAM_READ_GROUP_SAMPLE_MISSING",
                Severity.WARNING,
                "BAM read groups are missing SM tags",
                sample_id=sample_id,
                path=path,
                detail=f"read groups without SM: {read_groups_without_sample}",
                remediation="Populate @RG SM tags before cohort analysis.",
            )
        )

    if sample_id:
        if not sample_names:
            findings.append(
                _finding(
                    "HTS_BAM_SAMPLE_MISSING",
                    Severity.ERROR,
                    "BAM header has no usable @RG SM sample name",
                    sample_id=sample_id,
                    path=path,
                    remediation="Add read groups with the manifest sample ID in the SM tag.",
                )
            )
        elif sample_id not in sample_names:
            findings.append(
                _finding(
                    "HTS_BAM_SAMPLE_MISMATCH",
                    Severity.ERROR,
                    "Manifest sample does not match BAM @RG SM tags",
                    sample_id=sample_id,
                    path=path,
                    detail=f"expected {sample_id!r}; BAM samples: {', '.join(sample_names)}",
                    remediation="Correct the manifest sample ID or use the matching BAM.",
                )
            )
        if len(sample_names) > 1:
            findings.append(
                _finding(
                    "HTS_BAM_MULTIPLE_SAMPLES",
                    Severity.ERROR,
                    "BAM contains read groups from multiple sample names",
                    sample_id=sample_id,
                    path=path,
                    detail=", ".join(sample_names),
                    remediation="Split the BAM by sample before cohort analysis.",
                )
            )

    if reference is not None and contigs:
        findings.extend(
            _dictionary_findings(
                contigs,
                reference,
                sample_id=sample_id,
                path=path,
                source="BAM",
                require_all_contigs=True,
            )
        )

    metrics = {
        "header_bytes": header_bytes,
        "reference_count": len(contigs),
        "bgzf": is_bgzf,
        "bgzf_eof": has_eof,
        "indexed": indexed,
    }
    return Inspection(
        kind="bam",
        path=str(path),
        sample_names=sample_names,
        contigs=contigs,
        metrics=metrics,
        findings=tuple(findings),
    )

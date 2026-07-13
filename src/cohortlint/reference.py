from __future__ import annotations

import hashlib
from pathlib import Path

from .model import CohortLintError, ReferenceIndex


def load_fai(path: Path) -> ReferenceIndex:
    """Load a FASTA ``.fai`` file and derive a sequence-dictionary fingerprint.

    Only contig names and lengths participate in the fingerprint.  Byte offsets
    and line widths describe a particular FASTA serialization, not the reference
    assembly itself, and therefore should not make otherwise compatible files
    appear different.
    """

    path = Path(path)
    contigs: list[tuple[str, int]] = []
    seen: set[str] = set()

    try:
        with path.open("rt", encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.rstrip("\r\n")
                if not line:
                    continue

                fields = line.split("\t")
                if len(fields) < 2:
                    raise CohortLintError(
                        f"Invalid FASTA index {path}: line {line_number} has fewer "
                        "than two tab-separated fields"
                    )

                name = fields[0]
                if not name:
                    raise CohortLintError(
                        f"Invalid FASTA index {path}: line {line_number} has an empty "
                        "contig name"
                    )
                if name in seen:
                    raise CohortLintError(
                        f"Invalid FASTA index {path}: duplicate contig {name!r} on "
                        f"line {line_number}"
                    )

                try:
                    length = int(fields[1])
                except ValueError as exc:
                    raise CohortLintError(
                        f"Invalid FASTA index {path}: contig {name!r} has a non-integer "
                        f"length on line {line_number}"
                    ) from exc
                if length <= 0:
                    raise CohortLintError(
                        f"Invalid FASTA index {path}: contig {name!r} has non-positive "
                        f"length {length} on line {line_number}"
                    )

                seen.add(name)
                contigs.append((name, length))
    except CohortLintError:
        raise
    except (OSError, UnicodeError) as exc:
        raise CohortLintError(f"Could not read FASTA index {path}: {exc}") from exc

    if not contigs:
        raise CohortLintError(f"FASTA index {path} contains no contigs")

    canonical_dictionary = "".join(
        f"{name}\t{length}\n" for name, length in contigs
    ).encode("utf-8")
    fingerprint = hashlib.sha256(canonical_dictionary).hexdigest()

    return ReferenceIndex(
        path=path,
        contigs=tuple(contigs),
        fingerprint=fingerprint,
    )

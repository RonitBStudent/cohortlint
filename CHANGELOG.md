# Changelog

All notable changes to CohortLint are documented here.

## 0.1.1 — 2026-07-13

Owner-pilot hardening release.

### Fixed

- Refuse report or discovered-manifest writes that could overwrite inputs.
- Use no-clobber atomic report writes with an explicit `--force` for known reports.
- Convert expected filesystem failures into exit code 2 without Python tracebacks.
- Scan a shared joint VCF once while checking every manifest sample membership.
- Detect unequal FASTQ mate counts proven at a bounded scan boundary.
- Validate VCF fileformat declarations, exact mandatory columns, and compression suffixes.
- Avoid treating observed VCF record contigs as a complete sequence dictionary.
- Mark corrupt VCF scans incomplete in machine-readable metrics.
- Harden CRAM timeouts, sample metadata, sequence dictionaries, and reference comparison.
- Reject invalid manifest paths and explicitly flag unsupported multi-lane/chunk discovery.

### Added

- Repeatable installed-package owner pilot covering 12 black-box scenarios.
- Public HTSlib BAM/VCF fixture validation and an evidence report.
- Focused regression tests, bringing the suite to 53 tests.

## 0.1.0 — 2026-07-13

Initial pilot release.

### Added

- Manifest validation and conservative file discovery.
- Streaming plain/gzip FASTQ structural and paired-read validation.
- BAM header, read-group, reference dictionary, EOF, and index checks.
- Plain/gzip VCF header, sample, contig, ordering, REF, and index checks.
- Optional CRAM validation through samtools.
- Cross-file reference dictionary compatibility checks.
- Human-readable and JSON reports with stable finding codes.
- CI-friendly exit thresholds and a self-contained broken-cohort demo.

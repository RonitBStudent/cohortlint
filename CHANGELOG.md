# Changelog

All notable changes to CohortLint are documented here.

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

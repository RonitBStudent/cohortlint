# CohortLint pilot validation evidence

This document records evidence available for the CohortLint v0.1.1 supervised
pilot as of 2026-07-13. It distinguishes synthetic tests, mocked external-tool
tests, and public-format fixtures because each supports a different level of
confidence.

These results show that the exercised checks behaved as expected in the listed
environments. They do not establish clinical validity, complete file integrity,
performance at production scale, or compatibility with untested platforms and
workflows. See [PILOT.md](PILOT.md) for the lab protocol and research-use
boundary.

## Evidence summary

| Evidence class | Input and environment | Result | What the result supports |
| --- | --- | --- | --- |
| Automated synthetic tests | 53 tests under Python 3.12.2 and Python 3.13.3 | All 53 passed under both interpreters | Regression coverage for the manifest, FASTQ, BAM-header, VCF, reference, reporting, CLI, and error-handling behaviors represented by the synthetic fixtures |
| Mocked external-tool tests | Synthetic CRAM paths and mocked `samtools` subprocess responses within the 53-test suite | Passed | CRAM command construction, header parsing, finding generation, and subprocess error handling; not decoding of a real CRAM file |
| Installed black-box pilot | Fresh wheel installations on macOS 26.5.1 arm64 | All 12 scenarios passed | Packaging, installed-command behavior, deterministic output, bounded/full semantics on synthetic data, exit contracts, input protection, and supported-scope warnings on the owner device |
| Official public BAM fixture | HTSlib `range.bam` with adjacent BAI; expected sample `ERS225193` | Passed | Reading a real BAM header, matching the sample name, recognizing the adjacent index, and accepting the exercised BAM structure |
| Official public single-sample VCF fixture | HTSlib `index.vcf`; expected sample `ERS220911`; 621 records | Full scan passed | Header/sample parsing and record-level checks across all 621 records in this public fixture |
| Official public joint VCF fixture | HTSlib `vcf44_1.vcf`; expected samples `HG00096` and `HG00097`; 28 records | Full scan passed | Shared joint-VCF sample matching and record-level checks across all 28 records |
| Public-fixture negative control | Incompatible VCF dictionaries combined in one cohort check | Correctly returned a blocking dictionary incompatibility finding | Detection of the deliberately established cross-file dictionary mismatch |

## Traceability

The public fixtures were taken from the official HTSlib repository pinned to
commit [`b499e4c7d44eebe85b2f1c13c2d33266f169e8bb`](https://github.com/samtools/htslib/tree/b499e4c7d44eebe85b2f1c13c2d33266f169e8bb/test).
Pinning the source commit prevents a later fixture change from silently altering
the evidence.

The CohortLint owner-device evidence applies to the frozen v0.1.1 pilot build.
Before any external participant starts, the coordinator must record its exact
commit or release tag and the SHA-256 checksum of each distributed wheel. The
validation suite and black-box pilot should be rerun against that exact build.

The locally validated wheel was `cohortlint-0.1.1-py3-none-any.whl` with
SHA-256 `9f164af673cf5abce4e362e7b23e868cb67ba70be77549b7616717dba02d8dc7`.
The official public inputs used in this round had these SHA-256 checksums:

- `range.bam`: `e15d14e3994027d433431c960bf1c5f2d6939f26b5094cd5a86bc6229a5b2661`
- `range.bam.bai`: `f06ef0c00e8ee31d23c16ff78db7e022baec70e430dcb5e0888c6aa94435364b`
- `index.vcf`: `d99c0251010dae47b019b85bb732865fb910cb680e7b43ea3a4b49fcf8216304`
- `vcf44_1.vcf`: `cff5deec04136864c2876a227a8071b9f4811cb84e9c37eb2d10b7b54c3663e6`

## Automated synthetic and mocked evidence

The 53-test suite passed locally with Python 3.12.2 and again with Python
3.13.3. Its fixtures are intentionally small and constructed to exercise known
conditions, including malformed manifests, missing files, FASTQ framing and
pairing problems, malformed VCF records, sample mismatches, BAM header and
reference-dictionary errors, index absence, CLI exit codes, reporting, and
cross-file compatibility.

Synthetic passing controls show that those encoded expectations remain stable;
they are not a substitute for the diversity of files produced by sequencing
instruments and institutional pipelines.

CRAM tests in this suite mock `samtools` output and failures. They support the
Python-side orchestration and parsing logic only. No actual CRAM container was
decoded in this validation round.

## Installed black-box evidence

Fresh wheels were installed and exercised on macOS 26.5.1 arm64. The black-box
pilot invoked the installed `cohortlint` executable rather than importing its
Python internals. All 12 scenarios passed:

1. version output;
2. environment diagnosis;
3. manifest schema output;
4. repeatable demo output and paths containing spaces;
5. a clean synthetic full scan with deterministic JSON;
6. bounded-prefix scan semantics;
7. one-pass inspection of a shared synthetic joint VCF;
8. detection of seeded FASTQ pair-count and VCF sample mismatches;
9. the invalid-invocation exit-code contract without a traceback;
10. report no-clobber behavior and refusal to overwrite an input path;
11. an explicit unsupported warning for multipart/multilane FASTQ discovery;
    and
12. unchanged SHA-256 hashes for the synthetic inputs after scanning.

This is packaging and behavioral evidence on one owner device. It is not a
cross-platform performance study.

## Official HTSlib fixture evidence

The public fixtures add evidence from files maintained outside the CohortLint
project:

- `range.bam` plus its BAI passed with expected sample `ERS225193`. CohortLint
  parsed the real BAM header and recognized the adjacent index. This does not
  mean every alignment record in the BAM was decoded or validated.
- `index.vcf` passed a full 621-record scan with expected sample `ERS220911`.
- `vcf44_1.vcf` passed a full 28-record joint-sample scan for `HG00096` and
  `HG00097`.
- Combining incompatible VCF dictionaries produced the expected blocking
  cross-file error.

These are small compatibility fixtures, not representative cohorts. Their pass
results should not be generalized to all HTSlib-supported encodings, all VCF
versions, or multi-gigabyte production files.

## Explicitly untested or incomplete areas

- **Linux and Python 3.10 CI:** the GitHub Actions jobs are blocked by an account
  billing restriction. No test failure has been observed there, but a blocked
  job is not validation evidence.
- **Real FASTQ data:** FASTQ behavior has been exercised with synthetic plain
  and gzip files, not an instrument- or archive-produced FASTQ corpus.
- **Actual CRAM:** subprocess behavior is mocked in tests; no real CRAM plus
  CRAI has been checked with an installed `samtools` in this round.
- **Multi-gigabyte scale:** runtime, storage behavior, and memory have not been
  characterized on production-size cohorts.
- **Full BAM records:** CohortLint intentionally inspects the BAM header,
  dictionary, index presence, and EOF marker rather than decoding every
  alignment record. A passing BAM result is not a full BAM integrity check.
- **Multipart/multilane FASTQ:** the v0.1.1 manifest cannot represent multiple
  R1/R2 parts for one biological sample. Discovery reports this as unsupported;
  pilot inputs must be safely pre-merged.
- **Cross-institution behavior:** false-positive rates, remediation usefulness,
  and installation friction have not yet been measured independently across
  labs.
- **Other unsupported environments:** Windows, BCF, remote URLs, and cloud
  object stores are outside the v0.1.1 pilot scope.

## Evidence still required before broader beta use

1. Re-run the exact frozen pilot artifact on Linux and Python 3.10.
2. Validate a small, independently sourced real FASTQ pair, including a known
   clean control and a safe negative control.
3. Validate an actual CRAM/CRAI with a recorded `samtools` version and compare
   the result with direct `samtools quickcheck` and header inspection.
4. Measure bounded and full FASTQ/VCF runtime and peak memory on multi-gigabyte
   inputs.
5. Complete supervised pilots at multiple independent institutions and
   adjudicate every error and warning.
6. Record false positives, known misses, remediation outcomes, crashes, and
   willingness to use CohortLint again without overstating results from unknown
   ground truth.

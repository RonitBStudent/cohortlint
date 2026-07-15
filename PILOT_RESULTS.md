# CohortLint pilot validation evidence

This document records evidence available for the CohortLint v0.1.1 supervised
pilot as of 2026-07-14. It distinguishes synthetic tests, mocked external-tool
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
| Automated synthetic tests | 53 tests locally under Python 3.12.2 and 3.13.3, plus GitHub Actions under Python 3.10, 3.12, and 3.13 on Ubuntu and macOS | All 53 passed locally; all six hosted matrix jobs passed | Regression coverage for the manifest, FASTQ, BAM-header, VCF, reference, reporting, CLI, and error-handling behaviors represented by the synthetic fixtures |
| Mocked external-tool tests | Synthetic CRAM paths and mocked `samtools` subprocess responses within the 53-test suite | Passed | CRAM command construction, header parsing, finding generation, and subprocess error handling; not decoding of a real CRAM file |
| Installed black-box pilot | Fresh v0.1.1 wheel installations under Python 3.12.2 and 3.13.3 on macOS 26.5.1 arm64 | All 12 scenarios passed under both interpreters | Packaging, installed-command behavior, deterministic output, bounded/full semantics on synthetic data, exit contracts, input protection, and supported-scope warnings on the owner device |
| Public paired FASTQ | ENA run `DRR000140`; two gzip FASTQs; 2,030,204 36-bp reads per mate | Published MD5s matched; bounded and full scans passed with no findings | Complete streaming and mate-ID validation on an archive-produced Illumina FASTQ pair |
| Actual CRAM/CRAI | HTSlib `range.cram` plus CRAI and matching `ce.fa`/FAI; sample `ERS225193`; samtools 1.22.1 | `samtools quickcheck` and 112-record decode passed; CohortLint passed with no findings | Real CRAM subprocess integration, sample/header parsing, index recognition, and reference-dictionary matching |
| Large streaming benchmark | Deterministic gzip inputs with 3,000,000 150-bp read pairs and 5,000,000 variants; 2.105 GB logical input | Bounded run and two full runs passed; full reports were byte-identical; inputs were unchanged | Low-memory full streaming at the exercised scale, repeatable runtime, deterministic reports, and input immutability |
| Official public BAM fixture | HTSlib `range.bam` with adjacent BAI; expected sample `ERS225193` | Passed | Reading a real BAM header, matching the sample name, recognizing the adjacent index, and accepting the exercised BAM structure |
| Official public single-sample VCF fixture | HTSlib `index.vcf`; expected sample `ERS220911`; 621 records | Full scan passed | Header/sample parsing and record-level checks across all 621 records in this public fixture |
| Official public joint VCF fixture | HTSlib `vcf44_1.vcf`; expected samples `HG00096` and `HG00097`; 28 records | Full scan passed | Shared joint-VCF sample matching and record-level checks across all 28 records |
| Public-fixture negative control | Incompatible VCF dictionaries combined in one cohort check | Correctly returned a blocking dictionary incompatibility finding | Detection of the deliberately established cross-file dictionary mismatch |

## Traceability

The public fixtures were taken from the official HTSlib repository pinned to
commit [`b499e4c7d44eebe85b2f1c13c2d33266f169e8bb`](https://github.com/samtools/htslib/tree/b499e4c7d44eebe85b2f1c13c2d33266f169e8bb/test).
Pinning the source commit prevents a later fixture change from silently altering
the evidence.

The CohortLint owner-device evidence applies to the checksum-pinned v0.1.1
candidate artifact. Before any external participant starts, the coordinator
must record its exact commit or release tag and the SHA-256 checksum of each
distributed wheel. The validation suite and black-box pilot should be rerun
against that exact build.

The validated packaged-product tree was commit
`a5bb61932032150864d65115158b8717cac4edb5`; the evidence and pilot-harness
follow-up does not change packaged product files. The locally validated wheel
was `cohortlint-0.1.1-py3-none-any.whl` with SHA-256
`b8492c63dc697712cbaa088db5814572fb5818393a56eda51097ed7086db5162`.
The official public inputs used in this round had these SHA-256 checksums:

- `range.bam`: `e15d14e3994027d433431c960bf1c5f2d6939f26b5094cd5a86bc6229a5b2661`
- `range.bam.bai`: `f06ef0c00e8ee31d23c16ff78db7e022baec70e430dcb5e0888c6aa94435364b`
- `range.cram`: `ea9217f5a0dd7e57c0f2a94d55d6285d1e8d35cc741de53f12c19eecd0e84326`
- `range.cram.crai`: `fb03d738f4cf48a8cb0469796d6901d3541bfd60a93005607a8fa1abc1b94336`
- `ce.fa`: `5eca163c91918ada9774080ee2274208155f4d1b2d00700ee950cdd7b269508c`
- `ce.fa.fai`: `445a36da04b64dd49b8b964171c6d4b4cafc12e9bff50f8cf4ce33717fec6ff5`
- `index.vcf`: `d99c0251010dae47b019b85bb732865fb910cb680e7b43ea3a4b49fcf8216304`
- `vcf44_1.vcf`: `cff5deec04136864c2876a227a8071b9f4811cb84e9c37eb2d10b7b54c3663e6`

The ENA-published MD5 checksums for `DRR000140_1.fastq.gz` and
`DRR000140_2.fastq.gz` were, respectively,
`3372e8bfe6096adcbdef948d98932ecf` and
`a7b17dbbf0bffcd6f7012c378466f38d`. Both downloaded files matched.
Their local SHA-256 checksums were, respectively,
`92fd844bb652b94f7493f9a34b581b02e282027c6a2812f9131a6ff6670b30ce` and
`31c233c5a9d0c463f3f2460066f38e15afc2ec88a2a65d360b7e3f971913a28d`.

## Automated synthetic and mocked evidence

The 53-test suite passed locally with Python 3.12.2 and again with Python
3.13.3. GitHub Actions run
[`29299922496`](https://github.com/RonitBStudent/cohortlint/actions/runs/29299922496)
then passed the full unit-test, installed-pilot, Ruff, and mypy workflow under
Python 3.10, 3.12, and 3.13 on both Ubuntu and macOS. Its fixtures are
intentionally small and constructed to exercise known conditions, including
malformed manifests, missing files, FASTQ framing and
pairing problems, malformed VCF records, sample mismatches, BAM header and
reference-dictionary errors, index absence, CLI exit codes, reporting, and
cross-file compatibility.

Synthetic passing controls show that those encoded expectations remain stable;
they are not a substitute for the diversity of files produced by sequencing
instruments and institutional pipelines.

CRAM unit tests in this suite mock `samtools` output and failures. They support
the Python-side orchestration and parsing logic; the separate public-fixture
run below supplies actual CRAM decoding evidence.

## Installed black-box evidence

The exact wheel was installed in fresh Python 3.12.2 and 3.13.3 environments
and exercised on macOS 26.5.1 arm64. The black-box pilot invoked the installed
`cohortlint` executable rather than importing its Python internals. All 12
scenarios passed under both interpreters:

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

## Public ENA paired-FASTQ evidence

Run [`DRR000140`](https://www.ebi.ac.uk/ena/browser/view/DRR000140) was selected
through the official ENA Portal API as a public paired Illumina run. ENA listed
compressed sizes of 46,761,080 and 48,764,085 bytes and published the two MD5
checksums recorded above. The downloaded files matched both checksums and
passed independent gzip integrity checks.

The exact installed wheel completed a bounded scan in 0.31 seconds with
26,705,920 bytes maximum resident set size. A full scan completed in 23.27
seconds with 27,049,984 bytes maximum resident set size. It read 2,030,204
records and 73,087,344 bases from each mate, compared every pair identifier,
reported both streams complete, and returned no errors or warnings. The MD5
checksums were unchanged after both scans.

This is one older short-read archive run, not evidence across current
instruments, variable-length reads, quality encodings, or institutional
preprocessing conventions.

## Actual CRAM/CRAI evidence

The official HTSlib `range.cram`, `range.cram.crai`, `ce.fa`, and `ce.fa.fai`
fixtures were downloaded from the same pinned HTSlib commit used for the BAM
and VCF evidence. An isolated Bioconda samtools 1.22.1 environment was used.

`samtools quickcheck` passed, `samtools view -c -T ce.fa range.cram` decoded 112
records, and the decoded header identified sample `ERS225193` and seven
reference sequences. CohortLint then passed the CRAM with no findings, found
the adjacent CRAI, matched the sample, and matched all reference names,
lengths, and order against `ce.fa.fai`.

This is a small compatibility fixture. It does not characterize large CRAM
runtime, reference-cache behavior, every CRAM version, or files emitted by
diverse production pipelines.

## Large streaming benchmark

A deterministic synthetic workload contained 3,000,000 paired 150-bp FASTQ
records and 5,000,000 sorted VCF records. The logical uncompressed input was
2,104,889,003 bytes; the three gzip files occupied 370,539,803 bytes. The VCF
was intentionally left without an index, so `HTS_VCF_INDEX_MISSING` was the
single expected warning.

The exact installed wheel's bounded run completed in 0.44 seconds with
26,836,992 bytes maximum RSS, inspected exactly 10,000 read pairs and 10,000
VCF records, and correctly reported both scans incomplete. Full runs completed
in 88.85 and 89.06 seconds with maximum RSS of 26,886,144 and 27,164,672 bytes.
Each full run inspected all 3,000,000 read pairs and 5,000,000 variants and
reported complete scans.

The two full JSON reports were byte-identical with SHA-256
`8949d7de7ba3fa2ea87d280f4b6ed2b85589dc205494eee7c1b76443e771f571`.
Before/after SHA-256 checksums were identical for the manifest and every input:

- `cohort.csv`: `bc8c2a2b1a32515968e48901d0c0b6a1008eba515ad01282ca8d77145607ae3b`
- `BENCH_R1.fastq.gz`: `3e289f1fdac223b2f4f7c41c072c37142a25b9292453893d04338d6de1f5a403`
- `BENCH_R2.fastq.gz`: `3f967210469443768cd8dd1d1558bc3631af10a6a39f56df3d1260fc9c3848bf`
- `BENCH.vcf.gz`: `326b0ad0f8b4ae77b0936552e63a2f39cbd9eecf9195b91357c46c2e0e483024`

This demonstrates bounded and full streaming behavior at multi-gigabyte
logical scale on one Mac. It is not a production storage, concurrent-cohort,
network-filesystem, or biological-diversity benchmark.

## Official HTSlib fixture evidence

The public fixtures add evidence from files maintained outside the CohortLint
project:

- `range.bam` plus its BAI passed with expected sample `ERS225193`. CohortLint
  parsed the real BAM header and recognized the adjacent index. This does not
  mean every alignment record in the BAM was decoded or validated.
- `range.cram` plus its CRAI passed with expected sample `ERS225193` and the
  matching `ce.fa` dictionary under samtools 1.22.1.
- `index.vcf` passed a full 621-record scan with expected sample `ERS220911`.
- `vcf44_1.vcf` passed a full 28-record joint-sample scan for `HG00096` and
  `HG00097`.
- Combining incompatible VCF dictionaries produced the expected blocking
  cross-file error.

These are small compatibility fixtures, not representative cohorts. Their pass
results should not be generalized to all HTSlib-supported encodings, all VCF
versions, or multi-gigabyte production files.

## Explicitly untested or incomplete areas

- **Real FASTQ diversity:** one archive-produced paired Illumina run passed,
  but current instruments, variable read lengths, alternate header
  conventions, and institutional preprocessing remain uncharacterized.
- **Production CRAM diversity and scale:** one small official CRAM/CRAI fixture
  passed with samtools 1.22.1; large files, other CRAM versions, reference
  caches, and production pipeline outputs remain uncharacterized.
- **Production-scale systems:** the multi-gigabyte logical benchmark was
  synthetic and local. Concurrent cohorts, network filesystems, shared storage,
  and multi-gigabyte real inputs have not been characterized.
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

1. Complete supervised pilots at multiple independent institutions and
   adjudicate every error and warning.
2. Exercise additional independently sourced FASTQ pairs spanning current
   instruments, read lengths, and header conventions, including a safe known
   negative control.
3. Exercise a larger production-origin CRAM/CRAI under the participant's
   approved reference setup and compare it with direct samtools checks.
4. Measure bounded and full runtime on at least one approved multi-gigabyte
   real handoff and a representative shared or network filesystem.
5. Record false positives, known misses, remediation outcomes, crashes, and
   willingness to use CohortLint again without overstating results from unknown
   ground truth.

# CohortLint

**Preflight interoperability checks for multi-site genomic cohorts.**

[![CI](https://github.com/RonitBStudent/cohortlint/actions/workflows/ci.yml/badge.svg)](https://github.com/RonitBStudent/cohortlint/actions/workflows/ci.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-3776ab)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

> **Status:** v0.1.1 research alpha. The checked scope has 53 automated tests and has passed owner-device synthetic and public-fixture pilots. It is not a clinical, production, or release-certification tool.

CohortLint flags selected structural problems before sequencing files are transferred, merged, or submitted to an expensive workflow. It checks relationships between metadata and FASTQ, BAM, CRAM, and VCF files—not only whether each file can be opened.

```text
CHECK RESULT: BLOCKING FINDINGS
════════════════════════════════
96 samples  4 sites  241 files
3 errors  2 warnings  1 note

ERROR   FASTQ_PAIR_ID_MISMATCH  sample=P017
        Read-1 and read-2 identifiers diverge at record 18,291

ERROR   HTS_BAM_SAMPLE_MISMATCH  sample=P042
        Manifest sample 'P042' is absent from BAM read groups

ERROR   COHORT_CONTIG_SET_MISMATCH
        Files use incompatible contig naming conventions
```

## Why it exists

Genomic files are often valid individually but incompatible as a cohort. Common examples include:

- a sample sheet pointing at the wrong mate or subject;
- a BAM read group containing a different sample ID;
- one site using `chr1` and another using `1`;
- files with incompatible reference dictionaries;
- a joint VCF that does not contain the expected samples;
- truncated gzip files or missing BGZF EOF markers discovered only after a workflow starts.

FastQC, samtools, and bcftools remain the authoritative format- and assay-level tools. CohortLint complements them with a project-level question: **do the declared files, identifiers, and reference dictionaries agree in the checks performed?**

## Quick start

CohortLint requires Python 3.10+ and has no core runtime dependencies. For a pilot, install the immutable release tag in a fresh virtual environment:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install "cohortlint @ git+https://github.com/RonitBStudent/cohortlint.git@v0.1.1"
cohortlint doctor
cohortlint demo --output /tmp/cohortlint-demo
```

The source-checkout launcher also works without installation:

```sh
git clone https://github.com/RonitBStudent/cohortlint.git
cd cohortlint
./cohortlint --help
```

## Check a cohort

```sh
./cohortlint check cohort.csv --reference references/GRCh38.fa.fai
```

To stream all FASTQ and VCF records in the checked scope:

```sh
./cohortlint check cohort.csv \
  --reference references/GRCh38.fa.fai \
  --full \
  --format json \
  --output preflight.json
```

Without `--full`, CohortLint checks the **first** 10,000 records per FASTQ/VCF for quick feedback; this is a bounded prefix, not a random statistical sample. `--full` does not decode every BAM record. BAM header structure and BGZF EOF-marker presence are checked, while CRAM validation requires `samtools`.

Report output is no-clobber by default. Use `--force` only to replace a known report; CohortLint refuses to overwrite a manifest, reference index, or discovered genomic input even with `--force`.

## Manifest

One row represents one biological sample. Paths can be absolute or relative to the manifest.

```csv
sample_id,site,fastq_1,fastq_2,alignment,variants
P001,Site-East,reads/P001_R1.fastq.gz,reads/P001_R2.fastq.gz,bam/P001.bam,vcf/cohort.vcf.gz
P002,Site-West,reads/P002_R1.fastq.gz,reads/P002_R2.fastq.gz,bam/P002.bam,vcf/cohort.vcf.gz
```

Each row may provide any combination of one FASTQ pair, one alignment, and one variant file. Multi-sample VCFs may be shared across rows and are scanned once; FASTQ and alignment files must have one unambiguous owner.

The v0.1.1 schema cannot represent multiple FASTQ lanes or chunks for one biological sample. Pre-merge technical FASTQ parts before a pilot. `discover` reports `DISCOVERY_MULTIPART_FASTQ_UNSUPPORTED` rather than silently modeling lanes as separate samples.

Print the contract with:

```sh
./cohortlint schema
```

If no manifest exists, generate a conservative draft from filenames:

```sh
./cohortlint discover ./incoming-data --output cohort.csv
./cohortlint check cohort.csv
```

Discovery records ambiguous roles for a researcher to resolve and refuses to create a header-only manifest when no supported files are found. Always review a generated manifest before running a check.

## Checks in v0.1.1

### Metadata and ownership

- Required schema and unique, non-empty sample identifiers
- Relative-path resolution and missing files
- Duplicate ownership of reads and alignments
- FASTQ mate naming and role consistency
- Presence and count of contributing-site labels

### FASTQ

- Streaming four-line record validation
- Plain or gzip input and complete gzip integrity in `--full` mode
- Sequence/quality length and symbol validation
- `/1`, `/2`, and Illumina mate-tag normalization
- Paired identifier synchronization and pair-count mismatches proven within the scanned prefix or at EOF
- Read count, length, bases, and GC metrics

### BAM and CRAM

- BAM magic/header parsing and BGZF header/EOF-marker presence
- Read-group `SM` agreement with the manifest
- Reference sequence dictionary extraction
- Presence of adjacent BAI/CSI/CRAI indexes; index contents and freshness are not validated
- CRAM `samtools quickcheck` and header inspection

### VCF

- Required `##fileformat`, exact mandatory header columns, and record structure
- Manifest sample presence in single- or multi-sample VCFs
- Contig dictionary and declared length checks
- Position ordering and REF allele syntax
- Missing tabix/CSI indexes for compressed VCFs

### Cohort interoperability

- Cross-file declared contig-name comparison
- Cross-file declared contig-length comparison
- Comparison with an explicit FASTA `.fai`
- Stable reference-dictionary fingerprint in JSON output

## Automation

Text is intended for researchers; JSON is stable machine-readable output for workflows. JSON is not a complete audit artifact and can expose absolute paths, sample/read-group names, VCF sample names, and mismatched read identifiers. Redact it before sharing.

```sh
./cohortlint check cohort.csv --format json --output preflight.json
./cohortlint check cohort.csv --format json --output preflight.json --force
./cohortlint check cohort.csv --fail-on warning
```

Exit codes:

- `0`: no findings at the configured failure threshold
- `1`: cohort findings reached the configured failure threshold
- `2`: CohortLint could not run because an argument or input contract was invalid
- `130`: interrupted

## Test

```sh
PYTHONPATH=src python3 -m unittest discover -s tests -v
```

Fixtures deliberately include truncated gzip files, shifted read pairs, malformed VCF records, sample-ID disagreements, invalid BAM headers, missing indexes, and incompatible reference dictionaries.

Run the installed-package owner pilot and retain its evidence:

```sh
python scripts/run_local_pilot.py \
  --executable cohortlint \
  --workspace ./pilot-evidence
```

See [PILOT_RESULTS.md](PILOT_RESULTS.md) for the current test matrix, public HTSlib fixture results, and what remains unvalidated.

## Pilot testing

We are looking for a small number of sequencing cores, bioinformatics teams, and multi-site studies willing to run a supervised alpha evaluation on public or institution-approved de-identified data. CohortLint contains no telemetry or upload client.

See [PILOT.md](PILOT.md) for the pilot protocol, privacy guidance, and feedback questions. Please do not attach genomic data, protected health information, sample identifiers, or unredacted JSON reports to a public GitHub issue.

## Scientific boundaries

CohortLint checks selected interoperability conditions; it does not validate biological correctness. “No blocking findings” applies only to the performed checks. It does not prove reference-sequence identity, BAM interior-record integrity, index validity, consent, diagnosis, phenotype labels, ancestry inference, variant calling, contamination status, or study design. Run established format and assay QC as well.

## Roadmap toward cross-institution use

- GA4GH refget sequence digests instead of dictionary-only fingerprints
- Multi-lane and multi-library manifest modeling
- Configurable institutional metadata schemas
- BCF and remote object-store support
- SARIF and workflow-engine integrations
- Signed transfer manifests and content checksums
- Real-world validation with public multi-center cohorts

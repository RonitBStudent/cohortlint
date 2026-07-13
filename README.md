# CohortLint

**Preflight interoperability checks for multi-site genomic cohorts.**

[![CI](https://github.com/RonitBStudent/cohortlint/actions/workflows/ci.yml/badge.svg)](https://github.com/RonitBStudent/cohortlint/actions/workflows/ci.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-3776ab)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

> **Status:** v0.1 pilot release. The supported checks are functional and tested, but CohortLint has not yet been validated across enough institutions for production certification.

CohortLint finds structural problems before sequencing files are transferred, merged, or submitted to an expensive workflow. It checks the relationship between metadata and FASTQ, BAM, CRAM, and VCF files—not only whether each file can be opened.

```text
COHORT NOT READY
════════════════
96 samples  4 sites  241 files
3 errors  2 warnings  1 note

ERROR   FASTQ_PAIR_ID_MISMATCH  sample=P017
        Read-1 and read-2 identifiers diverge at record 18,291

ERROR   BAM_SAMPLE_MISMATCH  sample=P042
        Manifest sample 'P042' is absent from BAM read groups

ERROR   COHORT_CONTIG_SET_MISMATCH
        Files use incompatible contig naming conventions
```

## Why it exists

Genomic files are often valid individually but incompatible as a cohort. Common examples include:

- a sample sheet pointing at the wrong mate or subject;
- a BAM read group containing a different sample ID;
- one site using `chr1` and another using `1`;
- files built against different reference patch releases;
- a joint VCF that does not contain the expected samples;
- truncated gzip/BGZF files discovered only after a workflow starts.

FastQC, samtools, and bcftools are excellent format-level tools. CohortLint focuses on the missing project-level question: **do these files, identifiers, and references agree with each other?**

## Quick start

CohortLint requires Python 3.10+ and has no core runtime dependencies.

```sh
cd cohortlint
./cohortlint doctor
./cohortlint demo --output /tmp/cohortlint-demo
```

The source-checkout launcher works without installation. Optionally install the console command:

```sh
python3 -m pip install -e .
cohortlint --help
```

## Check a cohort

```sh
./cohortlint check cohort.csv --reference references/GRCh38.fa.fai
```

For release or archival validation, scan complete FASTQ and VCF files:

```sh
./cohortlint check cohort.csv \
  --reference references/GRCh38.fa.fai \
  --full \
  --format json \
  --output preflight.json
```

Without `--full`, CohortLint samples up to 10,000 records per FASTQ/VCF for quick interactive feedback. BAM headers and BGZF end markers are always checked. CRAM validation uses `samtools` when available.

## Manifest

One row represents one biological sample. Paths can be absolute or relative to the manifest.

```csv
sample_id,site,fastq_1,fastq_2,alignment,variants
P001,Site-East,reads/P001_R1.fastq.gz,reads/P001_R2.fastq.gz,bam/P001.bam,vcf/cohort.vcf.gz
P002,Site-West,reads/P002_R1.fastq.gz,reads/P002_R2.fastq.gz,bam/P002.bam,vcf/cohort.vcf.gz
```

Each row may provide any combination of raw reads, an alignment, and variants. Multi-sample VCFs may be shared across rows; FASTQ and alignment files must have one unambiguous owner.

Print the contract with:

```sh
./cohortlint schema
```

If no manifest exists, generate a conservative draft from filenames:

```sh
./cohortlint discover ./incoming-data --output cohort.csv
./cohortlint check cohort.csv
```

Discovery never silently chooses between ambiguous files. It records the ambiguity for a researcher to resolve.

## Checks in v0.1

### Metadata and ownership

- Required schema and unique, non-empty sample identifiers
- Relative-path resolution and missing files
- Duplicate ownership of reads and alignments
- FASTQ mate naming and role consistency
- Contributing-site provenance

### FASTQ

- Streaming four-line record validation
- Plain or gzip input and complete gzip integrity in `--full` mode
- Sequence/quality length and symbol validation
- `/1`, `/2`, and Illumina mate-tag normalization
- Paired identifier synchronization and complete pair counts
- Read count, length, bases, and GC metrics

### BAM and CRAM

- BAM magic/header and BGZF EOF integrity
- Read-group `SM` agreement with the manifest
- Reference sequence dictionary extraction
- Missing BAI/CSI/CRAI indexes
- CRAM `samtools quickcheck` and header inspection

### VCF

- Required headers and record structure
- Manifest sample presence in single- or multi-sample VCFs
- Contig dictionary and declared length checks
- Position ordering and REF allele syntax
- Missing tabix/CSI indexes for compressed VCFs

### Cohort interoperability

- Cross-file contig-name comparison
- Cross-file contig-length comparison
- Comparison with an explicit FASTA `.fai`
- Stable reference-dictionary fingerprint in JSON output

## Automation

Text is intended for researchers; JSON is intended for workflows and audit systems.

```sh
./cohortlint check cohort.csv --format json --output preflight.json
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

## Pilot testing

We are looking for sequencing cores, bioinformatics teams, and multi-site studies willing to evaluate CohortLint on de-identified project manifests. CohortLint runs locally and never uploads genomic data.

See [PILOT.md](PILOT.md) for the pilot protocol, privacy guidance, and feedback questions. Please do not attach genomic data, protected health information, sample identifiers, or unredacted JSON reports to a public GitHub issue.

## Scientific boundaries

CohortLint validates interoperability, not biological validity. A passing report does not prove correct consent, diagnosis, phenotype labels, ancestry inference, variant calling, contamination status, or study design. It should run before established assay-specific QC—not replace it.

## Roadmap toward cross-institution use

- GA4GH refget sequence digests instead of dictionary-only fingerprints
- Configurable institutional metadata schemas
- BCF and remote object-store support
- SARIF and workflow-engine integrations
- Signed transfer manifests and content checksums
- Real-world validation with public multi-center cohorts

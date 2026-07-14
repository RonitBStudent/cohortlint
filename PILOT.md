# CohortLint supervised pilot protocol

CohortLint v0.1.1 is an alpha research tool for preflight interoperability
checks. This two-week pilot evaluates whether its findings are correct,
actionable, and usable in real genomic file handoffs. It is not a production
certification program.

Current owner-device and public-fixture results are recorded in
[PILOT_RESULTS.md](PILOT_RESULTS.md). Those results establish behavior only for
the environments and inputs listed there.

## Pilot questions

1. Can a new user install the pinned build and produce an interpretable report?
2. Which findings identify confirmed handoff problems, and which are false
   positives in legitimate workflows?
3. Which important, independently known problems are missed?
4. Are the messages and proposed remediations clear enough to act on?
5. Does CohortLint run repeatably, without modifying its inputs, within an
   acceptable resource budget?

## Research-use boundary

CohortLint must not be used for clinical diagnosis, patient management,
regulated release decisions, or certification of a cohort. It does not assess
consent, phenotype accuracy, sample contamination, identity swaps, variant-call
quality, ancestry, study design, or biological validity. A report with no
blocking findings means only that no blocking problem was found by the checks
and scan scope that ran.

The default FASTQ and VCF check examines the first 10,000 records in each
unique file. It is a bounded prefix scan, not a random or statistically
representative sample. `--full` streams every FASTQ and VCF record. BAM checks
cover the header, reference dictionary, adjacent index presence, and BGZF EOF
marker; they do not validate every BAM alignment record. CRAM support depends
on `samtools` and has not yet been exercised with a real CRAM file.

## Who qualifies

A pilot participant should:

- work in a sequencing core, bioinformatics team, multi-site study, or workflow
  team that receives FASTQ, BAM, CRAM, or VCF files;
- have a command-line analyst who can use an isolated Python environment;
- have authority to run local software on the selected system and data;
- use public, synthetic, or locally approved de-identified research data;
- know the expected sample identifiers and reference assembly well enough to
  adjudicate every reported error and warning; and
- commit approximately 60–90 minutes across two weeks for onboarding, one
  representative run, adjudication, and a short closeout interview.

The first cohort of participants should deliberately cover different roles and
formats: at least one sequencing core, one team receiving collaborator data,
one multi-site study, and one workflow maintainer, with FASTQ, BAM, and VCF all
represented. Actual CRAM is optional exploratory work and must not determine
pilot success until it has its own validated baseline.

## Exclusions and unsupported inputs

Do not use the pilot for:

- clinical, diagnostic, or other regulated workflows;
- protected or directly identifying data without explicit institutional
  approval;
- multipart or multilane FASTQ inputs unless each sample's lanes or chunks have
  already been safely merged into one R1/R2 pair;
- BCF, remote URLs, or cloud object-store paths;
- Windows-only workflows; or
- a multi-gigabyte `--full` run until the bounded run has succeeded and the lab
  has approved the compute and time cost.

Multipart/multilane FASTQ discovery is intentionally reported as unsupported;
participants must not accept a partial discovered manifest or select one lane
as if it represented the complete sample.

## Data safety and redaction

CohortLint runs locally and contains no telemetry or upload client. The
`check` command reads genomic inputs and writes only the report path explicitly
requested by the user. Use a new pilot report directory, retain existing source
checksums when available, and follow the lab's normal access, retention, DUA,
IRB, and security requirements. Local execution does not override institutional
policy.

Never send genomic files, manifests, credentials, protected metadata, or an
unredacted report to the maintainer. Treat both text and JSON reports as
potentially sensitive. Before sharing a copy, review and redact all of the
following:

- sample, participant, family, and institution/site identifiers;
- absolute paths, usernames, hostnames, storage locations, and private URLs;
- filenames that encode identifiers;
- BAM/CRAM read-group sample names and VCF sample names;
- FASTQ read identifiers and instrument/run identifiers in finding details;
- the `manifest`, `path`, `sample_id`, `sample_names`, `detail`, `executable`,
  and host/environment fields;
- custom or confidential contig names and reference metadata; and
- tokens, signed URLs, credentials, or other secrets.

Redact a copy, never the original evidence. Use stable pseudonyms such as
`SAMPLE_001` so related findings remain interpretable. Share aggregate counts
and finding codes whenever possible. Public GitHub issues may contain only
de-identified summaries and synthetic reproductions. Use a private GitHub
security advisory only for a suspected software vulnerability, not for routine
pilot feedback.

## Reproducible setup

The pilot coordinator must provide an immutable Git commit or release tag and,
when distributing a wheel, its SHA-256 checksum. Do not pilot from a moving
`main` branch. Record the identifier in the evidence worksheet.

Install the pinned checkout or supplied wheel in a fresh environment:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install /path/to/cohortlint-0.1.1-py3-none-any.whl
cohortlint --version
cohortlint doctor
```

`samtools` is required only for CRAM. Linux and Python 3.10 are useful pilot
targets because they have not yet completed the same owner-device validation as
macOS with Python 3.12 and 3.13.

## Two-week procedure

### Before day 1: qualification

The coordinator confirms the research-use boundary, supported file model,
local approval, operating system, formats, approximate cohort size, and a person
who can adjudicate the findings. Assign an anonymous pilot ID and agree on how
redacted feedback will be returned.

### Days 1–2: installation and orientation

1. Install the pinned build in a fresh environment.
2. Record `cohortlint --version` and a redacted `cohortlint doctor` result.
3. Run the deliberately inconsistent demonstration in a new directory:

   ```sh
   cohortlint demo --output pilot-demo
   ```

4. Confirm that the demonstration displays blocking findings. The `demo`
   command is an orientation exercise and its exit status is not a cohort
   validation result.
5. Select one approved public or de-identified handoff. Start with a small or
   representative subset when practical.

### Days 3–5: bounded real-data run

1. Create the six-column manifest manually, or create a draft and review every
   row before use:

   ```sh
   cohortlint discover ./incoming-data --output cohort.csv
   cohortlint schema
   ```

2. Confirm that each row represents one biological sample and no multipart or
   multilane FASTQ has been omitted.
3. Run the bounded prefix check and record the exit code and elapsed time:

   ```sh
   cohortlint check cohort.csv \
     --reference reference.fa.fai \
     --format json \
     --output pilot-bounded.json
   ```

4. Review every error and warning against source metadata or an established
   tool. Label each finding confirmed, false positive, or unresolved.

### Days 6–9: deeper and negative-control checks

1. If local policy and compute limits permit, scan all FASTQ and VCF records:

   ```sh
   cohortlint check cohort.csv \
     --reference reference.fa.fai \
     --full \
     --format json \
     --output pilot-full.json
   ```

2. Repeat one unchanged command and compare the reports for deterministic
   results.
3. Where a safe synthetic copy can be made, introduce one known structural
   problem, such as a manifest/VCF sample mismatch, and record whether the
   expected finding appears. Never alter the only copy of research data.
4. Compare relevant outcomes with established checks such as gzip integrity,
   `samtools quickcheck`, or bcftools header inspection when those tools are
   already part of the lab workflow.
5. If a confirmed problem can be safely remediated, rerun CohortLint and record
   whether the finding clears.

### Days 10–12: adjudication

Complete the evidence worksheet below. Discuss false positives, missed known
problems, unclear severity, runtime, installation friction, and any workflow
that could not be represented by the manifest. Do not calculate sensitivity or
recall from failures whose ground truth is unknown.

### Days 13–14: closeout

Hold a 20-minute interview and return only the agreed redacted worksheet. The
coordinator aggregates results by finding code and dataset class, records bugs
or unsupported workflows, and decides whether another supervised iteration is
needed before broader recruitment.

## Evidence worksheet

Copy this section into a private local note. Do not paste an unredacted report
into it.

### Participant and build

- Anonymous pilot ID:
- Lab role: sequencing core / receiving team / multi-site study / workflow team
- CohortLint version, commit or tag, and wheel SHA-256:
- OS and architecture:
- Python version:
- Samtools version, if CRAM was used:
- Installation method and minutes to successful `doctor`:

### Dataset description

- Public, synthetic, or locally approved de-identified data:
- Assay context, without participant information:
- FASTQ / BAM / CRAM / VCF files represented:
- Samples, unique files, approximate total bytes, and contributing sites:
- Reference assembly label and source of `.fai`:
- Compression and adjacent-index status:
- Known pre-existing or deliberately seeded failures:

### Run evidence

- Exact command with paths and identifiers redacted:
- Bounded or full scan:
- Exit code:
- Wall time and, if available, peak memory:
- Reported scan coverage or records inspected:
- Repeated output identical: yes / no / not tested
- Input checksums or size/mtime unchanged: yes / no / not tested
- Crash, traceback, hang, or unexpected file write: yes / no; details:

Record one row per finding; add rows as needed:

| Finding code | Severity | Confirmed / false positive / unresolved | Evidence used | Remediation clear? | Cleared after fix? |
| --- | --- | --- | --- | --- | --- |
|  |  |  |  |  |  |

### Missing behavior and usability

- Independently known failures CohortLint missed:
- Legitimate workflow CohortLint could not represent:
- Unclear messages, severities, or remediation:
- Time to first actionable report:
- Maintainer assistance required:
- Would the lab run this on its next non-clinical handoff? Why or why not?
- Highest-priority change before another pilot:

## Pilot success measures

The pilot is informative when every completed run has traceable build and
environment information, every finding is adjudicated, known negative controls
have explicit ground truth, and crashes and false positives are reported rather
than omitted. Advancement to a broader beta should additionally require:

- no input modification or unhandled traceback in completed pilots;
- reproducible installation from an immutable release artifact;
- no unexplained error on clean, independently reviewed controls;
- documented handling of every confirmed false positive;
- successful Linux and Python 3.10 coverage;
- real FASTQ and CRAM evidence, plus multi-gigabyte performance measurements;
  and
- multiple independent labs willing to use the tool again for non-clinical
  preflight checks.

## Recruitment email

**Subject:** Two-week pilot of a local genomic handoff checker

Hi [Name],

I am recruiting a small group of design partners for CohortLint, an open-source
alpha CLI that checks whether manifests, sample identifiers, and reference
dictionaries agree across FASTQ, BAM, CRAM, and VCF handoffs.

The pilot takes about 60–90 minutes over two weeks: install a pinned build, run
a demo and one approved public or de-identified research dataset, adjudicate the
findings, and join a 20-minute closeout. CohortLint runs locally; we will not ask
for genomic files, manifests, or unredacted reports.

This pilot excludes clinical use and multipart/multilane FASTQ unless the lanes
are already safely merged. Participants need Python 3.10+, command-line
familiarity, and someone who knows the expected sample IDs and reference.

If interested, please reply with your OS, file formats, approximate cohort size,
and whether the pilot input would be public or de-identified. Do not include
sample identifiers or data in your reply.

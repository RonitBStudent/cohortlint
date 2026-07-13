# CohortLint pilot protocol

The v0.1 pilot is designed to answer three questions:

1. Does CohortLint catch failures that matter in real genomic handoffs?
2. Which findings are false positives for legitimate institutional workflows?
3. Can a new user reach an actionable report in under five minutes?

## Who should participate

- Sequencing and bioinformatics core facilities
- Teams receiving FASTQ, BAM, CRAM, or VCF files from collaborators
- Multi-site studies harmonizing genomic cohorts
- Workflow maintainers responsible for input validation

## Data safety

CohortLint runs locally. Do not send the project genomic files or protected metadata. Before sharing a report, replace sample IDs, institution names, usernames, and sensitive paths. JSON reports are not automatically de-identified.

## Procedure

1. Install or clone CohortLint.
2. Run `cohortlint doctor` and record the environment.
3. Create or discover a manifest.
4. Run the bounded default check first.
5. Review every error and warning against source metadata.
6. Run `--full` only when local policy and compute limits permit it.
7. Report confirmed findings, false positives, missing checks, runtime, and usability problems.

```sh
cohortlint discover ./incoming-data --output cohort.csv
cohortlint check cohort.csv --reference reference.fa.fai
cohortlint check cohort.csv --reference reference.fa.fai --full --format json --output preflight.json
```

## Feedback template

- Cohort size and file types, without identifying data:
- Operating system and Python version:
- Time to first report:
- Confirmed errors detected:
- False positives:
- Important failures CohortLint missed:
- Unclear messages or remediations:
- Would you run it on the next incoming cohort? Why or why not?

Open a GitHub issue using only de-identified information. For sensitive security concerns, use a private repository security advisory.

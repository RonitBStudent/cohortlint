# Contributing to CohortLint

Contributions from researchers, sequencing facilities, and research software engineers are welcome.

## Before changing code

- Open an issue describing the failure mode or proposed check.
- Include the smallest safe reproduction possible.
- Never upload human genomic data, protected health information, institutional credentials, or unredacted sample identifiers.
- Prefer synthetic fixtures that preserve the structural failure without preserving biological data.

## Development setup

```sh
git clone https://github.com/RonitBStudent/cohortlint.git
cd cohortlint
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -e '.[dev]'
```

Run the complete validation suite:

```sh
python3 -m unittest discover -s tests -v
ruff check src tests cohortlint
mypy src
```

## Adding a check

Every new check should include:

1. A stable, format-prefixed finding code.
2. A severity justified by downstream risk.
3. A concise explanation and actionable remediation.
4. A synthetic failing fixture.
5. A passing control that guards against false positives.
6. Documentation of scientific or format limitations.

Pull requests should remain focused and explain the user-visible behavior, test evidence, and any compatibility implications.

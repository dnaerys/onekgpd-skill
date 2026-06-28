## 1000 Genomes Project Dataset Access

Natural language access to _**1000 Genomes Project dataset**_, hosted online in _[Dnaerys variant store](https://dnaerys.org/)_

Sequenced & aligned by _New York Genome Center_ (_GRCh38_). _3202 samples_: 2504 unrelated samples from phase
three panel + 698 samples from 602 family trios - [dataset details](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38)

### Key Features

- _real-time_ access to _138 044 723_ unique variants and _~442 billion_ individual genotypes

- variant, sample and genotype selection based on coordinates, annotations, zygosity, population

- filtering by VEP (impact, biotype, feature type, variant class, consequences), ClinVar Clinical Significance (202502),
  gnomADe + gnomADg 4.1, AlphaMissense Score & AlphaMissense Class annotations

  - annotated with VEP 115 / GENCODE 49
  - [GENCODE Primary set](https://www.gencodegenes.org/pages/gencode_primary/) transcripts
  - [full annotation composition](https://github.com/dnaerys/onekgpd-mcp/blob/master/docs/annotations.md)

- returned variants annotated with _HGVSp_, _gnomADe + gnomADg_, _AlphaMissense score_ + cohort-wide statistics

  - _**HGVSp**_ annotations are for **_Canonical transcripts_** to reduce LLMs cognitive load

- samples annotated with: _familyId, gender, paternalId, maternalId, relationship, population, superpopulation, phase3 indicator_

## What This Is

OneKGPd is an agentic **Skill**: a `SKILL.md` manifest with agent-facing
instructions, two command-line helper scripts, reference documentation, and a
bundled data asset. An agent reads `SKILL.md`, then invokes the scripts via
`uv run` to answer questions about the 1000 Genomes Project cohort. The scripts
are plain CLIs and can also be run directly.

There is no server to deploy, no API key, no `.env`, and no environment
variables — the only prerequisite is [`uv`](https://docs.astral.sh/uv/), which
provisions each script's dependencies from its inline
[PEP 723](https://peps.python.org/pep-0723/) metadata.

## Repository Layout

```
.
├── README.md
├── LICENSE
├── skills/
│   └── onekgpd/                        # the shipped Skill
│       ├── SKILL.md                    # skill manifest + agent-facing instructions
│       ├── scripts/
│       │   ├── onekgpd_api.py          # online variant/sample/kinship queries (Dnaerys client)
│       │   └── onekgpd_meta.py         # offline sample & population metadata (stdlib only)
│       ├── references/
│       │   ├── onekgpd_commands.md     # full per-command argument tables + output schema
│       │   └── annotation_vocabularies.md  # controlled vocabularies for the CSV filter flags
│       └── assets/
│           └── kgpe.json               # bundled pedigree/population data (3,202 records)
└── tests/                              # black-box test suite (see tests/README.md)
    ├── conftest.py                     # subprocess harness + live-endpoint gating
    ├── test_meta.py                    # offline metadata tests (always run)
    ├── test_cli_contract.py           # offline argparse/validation tests (always run)
    ├── test_e2e_live.py               # live endpoint E2E (auto-skips if unreachable)
    └── test_helpers_unit.py           # opt-in in-process unit tier (--run-unit)
```

## Architecture

The skill is organized in two independent layers, each a self-contained CLI
wrapper invoked through `uv run`. Both write the full result as JSON to a file
(`--output`, default a temp file) and print a concise summary to stdout. Sample
IDs are shared across both layers, so a cohort picked by population can be fed
straight into a variant query.

### Online variant layer — `scripts/onekgpd_api.py`

Wraps the [Dnaerys Python Client](https://dnaerys-python.readthedocs.io/en/latest/)
to query the public Dnaerys variant store at `db.dnaerys.org:443` over gRPC/TLS.
It owns the connection, server-side streaming, pagination, bounded retries on
transient errors, and JSON serialization. Ten subcommands:

- **Variants** — `count-variants`, `select-variants`, and their
  `-in-samples` forms (restricted to a named set of individuals).
- **Samples** — `count-samples`, `select-samples` (individuals carrying a
  matching variant in a region).
- **Homozygous-reference** — `count-samples-hom-ref`, `select-samples-hom-ref`
  at a single position.
- **Kinship** — `kinship` (relatedness degree + KING coefficient between two
  named individuals).
- **Dataset** — `dataset-info` (cohort totals; doubles as a connectivity check).

Its single runtime dependency (`dnaerys`) is declared inline, so no manual
install step is required.

### Offline metadata layer — `scripts/onekgpd_meta.py`

Answers population, sex, pedigree, and superpopulation questions with **no
network, no credentials, and no third-party dependencies** — Python standard
library only. It reads a data file bundled in the skill (`assets/kgpe.json`:
3,202 records, plain JSON, ~1 MB, loaded in-memory) and reproduces the MCP
server's metadata query semantics offline. Six subcommands:

- `sample-metadata` — family, gender, parents, children, population,
  superpopulation, and phase3 status for given sample IDs.
- `list-populations` / `list-superpopulations` — the 26 populations and 5
  superpopulations (`AFR`, `AMR`, `EAS`, `EUR`, `SAS`) with sample counts.
- `population-stats` / `superpopulation-summary` — sex split, phase3, and trio
  membership, per population and per superpopulation.
- `select-samples-by-population` — sample IDs by population and/or
  superpopulation.

Population/superpopulation values match case-insensitively by short code or full
name; sample IDs are case-sensitive.

### Skill design conventions

A few rules are structural to how the agent uses the skill (full detail in
[`skills/onekgpd/SKILL.md`](skills/onekgpd/SKILL.md)):

- **Resolve coordinates first.** The dataset is GRCh38 only; a gene/feature must
  be resolved to verified GRCh38 coordinates against an authoritative source
  before a region query — there is no source-side guard against a misplaced
  region.
- **Count before you select.** Every selection command has a paired counting
  command; size the result set first.
- **Zygosity defaults to both** heterozygous and homozygous carriage; narrow
  with `--het-only` / `--hom-only`.

## Running

```bash
# Online: dataset totals (also a connectivity check)
uv run skills/onekgpd/scripts/onekgpd_api.py dataset-info

# Online: count individuals carrying likely-pathogenic missense variants in a region
uv run skills/onekgpd/scripts/onekgpd_api.py count-samples \
  --chrom chr17 --start 43044292 --end 43170245 \
  --consequence MISSENSE_VARIANT --alpha-missense-class AM_LIKELY_PATHOGENIC

# Offline: per-population sex split, phase3, and trio membership
uv run skills/onekgpd/scripts/onekgpd_meta.py population-stats --populations YRI
```

See [`skills/onekgpd/SKILL.md`](skills/onekgpd/SKILL.md) for the full command guide and
[`skills/onekgpd/references/`](skills/onekgpd/references/) for the per-command argument tables
and output schemas.

## Tests

A black-box test suite lives at the repository root in [`tests/`](tests/) and
ships alongside the skill. Tests drive the skill exactly as an agent would —
invoking `uv run skills/onekgpd/scripts/…` as a subprocess and asserting on exit
code, stdout, stderr, and the JSON output file. They depend only on the standard
library + `pytest` and never import `dnaerys` (the opt-in unit tier is the sole
exception).

- **Offline** (`test_meta.py`, `test_cli_contract.py`) — deterministic, no
  network; always run.
- **Live E2E** (`test_e2e_live.py`) — exercises the commands against the public
  Dnaerys endpoint; auto-skips if it is unreachable and fails on unexpected data.
- **Opt-in unit** (`test_helpers_unit.py`) — in-process helper tests; collected
  only with `--run-unit`.

```bash
# Full suite (offline always; live auto-skips if the endpoint is down)
uv run --with pytest pytest tests -v

# Offline only
uv run --with pytest pytest tests/test_meta.py tests/test_cli_contract.py -v
```

See [`tests/README.md`](tests/README.md) for the detailed tier breakdown and how
to run each one.

## Replica of MCP server

This skill is a replica of [OneKGPd MCP](https://github.com/dnaerys/onekgpd-mcp) server.

## Examples

See [OneKGPd MCP Server examples](https://github.com/dnaerys/onekgpd-mcp/blob/master/examples/README.md).

## Privacy Policy

OneKGPd skill operates as a read-only interface layer for 1000 Genomes Project dataset.
Server does not collect, store, or transmit any user data. No conversation data is recorded.
No personal information is collected. No cookies, tracking mechanisms or authentication are used.

## Support

- Issues and questions: https://github.com/dnaerys/onekgpd-skill/issues
- Email: onekgpd@dnaerys.org

## License

This project is licensed under the MIT license - see the [LICENSE](./LICENSE) file for details.

# AMRScan 🧬

AMRScan v2 is evolving from a BLAST-based resistance gene scanner into a modular
framework for AMR genotyping, phenotype prediction, benchmarking, and
generalisability studies. It implements **AMR determinant calling,
evidence-preserving harmonization, and public AST cohort construction**.
Phenotype prediction and machine learning remain future, separate layers.

A sequence match does not by itself establish phenotypic resistance. A negative
report also does not establish susceptibility. Retain caller methods, database
snapshots, thresholds, and the original report when interpreting results.

## Project links

- [Source repository](https://github.com/taowis/amrscan-dev)
- Historical rendered demonstrations: [R](docs/AMRScan_R.html) and
  [Nextflow](docs/AMRScan_Nextflow.html)
- [Preprint](https://arxiv.org/abs/2507.08062)

The demonstrations and manuscripts describe the legacy version, not the v2
interface. Their historical data and generated reports remain in Git for now.

## Modular workflow

```text
main.nf → workflows/amrscan.nf
          ├─ modules/local/amrfinderplus.nf → raw report + log + provenance
          └─ modules/local/harmonize_amrfinder.nf → harmonized TSV
```

Use assembled nucleotide FASTA (`.fa`, `.fna`, `.fasta`, optionally `.gz`). Raw
FASTQ, protein/GFF combined searches and per-sample organism metadata are not
implemented in this first v2 interface. Sample IDs come from filenames with both
FASTA and gzip suffixes removed. They must be unique and contain only letters,
digits, `_`, `.`, `-`, beginning with a letter or digit. Duplicate IDs fail before
any tasks are scheduled.

Requirements:

- Nextflow **25.04.6** (the CI baseline), Java **17** or a compatible Nextflow JVM.
- Python **3.10+**, standard library only (CI: 3.12).
- AMRFinderPlus **4.0.23** is the documented caller baseline; the adapter targets
  the v4 nucleotide column names. This PR has **not** validated a real caller run.
  Other v4 releases require integration validation before benchmarking.
- A separately installed, indexed AMRFinderPlus database snapshot compatible
  with the caller. Record the release date/version and use an immutable directory.
  No database is downloaded, updated or bundled by this workflow.

Follow the [NCBI installation instructions](https://github.com/ncbi/amr/wiki/Installing-AMRFinder).
Install the caller and its BLAST/HMMER dependencies in your environment. The first
PR uses a local executor; it does not claim to provide a tested container or
Conda environment. Archive environment/package versions with each experiment.

```bash
git clone https://github.com/taowis/amrscan-dev.git
cd amrscan-dev

nextflow run main.nf \
  --input 'assemblies/*.fna' \
  --amrfinder_db /absolute/path/to/versioned-amrfinder-db \
  --outdir results/v2 --threads 4

# For a cohort in the same supported taxonomic group:
nextflow run main.nf \
  --input 'assemblies/*.fna.gz' \
  --amrfinder_db /absolute/path/to/versioned-amrfinder-db \
  --amrfinder_organism Escherichia --outdir results/escherichia
```

The organism option enables taxon-specific behavior, including curated mutation
screening. Without it, the result is not a comprehensive mutation screen. This
mode searches nucleotide assemblies only and does not perform protein HMM
searches. See [NCBI usage](https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus)
and [interpretation](https://github.com/ncbi/amr/wiki/Interpreting-results).

| Parameter | Default | Meaning |
| --- | --- | --- |
| `input` | required | One FASTA or a quoted file glob |
| `outdir` | `results/v2` | Publication directory |
| `threads` | `4` | CPUs per caller task |
| `amrfinder_db` | required | Existing, indexed database directory |
| `amrfinder_bin` | `amrfinder` | Executable name on PATH or absolute executable path |
| `amrfinder_organism` | unset | One supported taxon for all inputs in this run |
| `amrfinder_plus` | `false` | Include caller's plus elements; retain their types |
| `amrfinder_identity` | `-1` | Curated identity thresholds with caller fallback |
| `amrfinder_coverage` | `0.5` | Caller reference coverage cutoff, as a fraction |
| `validate_only` | `false` | Validate paths/parameters, without running processes |

Identity overrides in `[0,1]` replace the caller's curated identity policy; use
only with scientific justification. Provenance records the requested values;
curated per-family settings remain in the identified database snapshot. The
workflow requests all equally scoring references and applies no further top-hit
selection. AMRFinderPlus still applies its own curation and filtering rules.

Outputs per sample:

```text
results/v2/
├── raw/amrfinderplus/<sample>/
│   ├── <sample>.amrfinder.tsv       # Original caller report, unmodified
│   ├── <sample>.amrfinder.log       # Caller stdout/stderr, including runtime versions
│   └── <sample>.provenance.json     # Command, software version, input/db/report hashes
└── harmonized/<sample>.amr.tsv      # Standard columns + complete raw records
```

See [result schema](docs/result-schema.md). Database content is fingerprinted
per sample, including the relative filenames and file SHA-256 hashes; this
costs additional I/O on large cohorts. Keep the database immutable during a run.
A successful no-hit report produces a header-only harmonized TSV, while its raw
report and provenance remain available. Failed calls never produce harmonized
results; task logs and failure provenance remain in the Nextflow work directory.
Use a fresh output directory per experiment to avoid confusing old and new files.

## AST datasets

The separate AST builder archives official NCBI Pathogen Detection tables and
constructs an auditable *E. coli* cohort. It retains submitter phenotypes,
measurements, genome links, metadata, duplicates, conflicts and source hashes.

```bash
python3 bin/ast_dataset.py download --outdir datasets/ecoli/raw
python3 bin/ast_dataset.py build --snapshot datasets/ecoli/raw --outdir datasets/ecoli/processed
python3 bin/ast_dataset.py validate --snapshot datasets/ecoli/raw --processed datasets/ecoli/processed
```

AMRScan does not infer phenotypic resistance from the presence of an AMR gene.
AST phenotype and genomic AMR determinants remain separate evidence layers.
See the [AST builder guide](docs/ast-dataset-builder.md) for offline examples,
source contracts, eligibility rules and validation limits.

## Legacy compatibility

`workflow/AMRScan.nf`, its adjacent config, `scripts/AMRScan.R`, the R Markdown
examples, and historical files are retained. Their algorithms have not been
rewritten in this PR. The legacy Nextflow FASTA/no-hit path has an optional small
smoke test using BLAST+ and R; other historical paths are not newly validated.

```bash
nextflow run workflow/AMRScan.nf \
  --input data/GCF_037966445.1_ASM3796644v1_genomic.fna \
  --card_db db/protein_fasta_protein_homolog_model.fasta \
  --outdir results/legacy --threads 4
```

Legacy tools include BLAST+, R with `dplyr`/`magrittr`, and for FASTQ,
`Biostrings`/`ShortRead`. The old config is not a complete reproducible container
environment. Its highest-bitscore summary is a display ranking, not biological
truth or a phenotype call; inspect the full detailed output. The HTML report
process remains disabled as in the original workflow.

The standalone R script is a **historical demonstration with known defects**:
it changes to a fixed local directory, writes and reads different BLAST output
paths, and does not ensure the final directory exists. Do not treat
`tests/test-run.R` as an automated v2 test. Fixing that interface and expanding
legacy FASTQ validation belong in a separate compatibility PR.

## Tests

```bash
# Offline caller-wrapper and harmonization tests; no third-party Python packages
python3 -m unittest discover -s tests -p test_amrscan.py -v
python3 -m unittest discover -s tests -p test_ast_dataset.py -v

# Nextflow syntax/config, error handling, actual harmonization process, optional legacy smoke
python3 -m unittest discover -s tests -p test_nextflow.py -v

# Entire lightweight suite (missing external tools are explicitly skipped)
python3 -m unittest discover -s tests -p 'test_*.py' -v
```

[Testing details and coverage commands](tests/README.md) describe fixture scope
and distinguish mocks, upstream expected reports, and real integration tests.
GitHub Actions pins Nextflow, checks unit branch coverage, and runs workflow
validation without AMR databases. AMRFinderPlus integration remains pending.

## Directory layout and next steps

```text
bin/                    Python caller/provenance and harmonization CLI
conf/                   Executor/resource defaults
modules/local/          DSL2 caller and harmonization modules
workflows/              Composable v2 orchestration
main.nf                 v2 entry point
nextflow.config         v2 defaults
workflow/               Legacy entry point/config
scripts/                Legacy standalone R demonstration
schemas/                AST raw and normalized schema contracts
training/               AST ingestion and future model-training contracts
benchmark/              Lineage/country/time evaluation contract
tests/                  Small fixtures and automated tests
docs/                   Schema, validation notes and historical reports
paper/                  Historical manuscripts, retained
data/, db/, results/    Historical tracked examples; new generated files ignored
```

Next, validate real AMRFinderPlus runs with a pinned software/database
environment, known positive/negative controls, and cohort accession provenance.
Further work includes adjudicating AST assay/breakpoint metadata, additional
caller adapters (RGI/CARD and ResFinder), and leakage-aware evaluation. See
[training](training/README.md), [benchmark](benchmark/README.md), and
[tracked-data cleanup notes](docs/tracked-data-audit.md).

The accession-verified AMRFinderPlus validation implementation and its current
limits are described in [docs/amrfinder-validation.md](docs/amrfinder-validation.md).
The validation cohort manifest is currently empty because no Assembly-linked
AST isolate could be verified in the network-restricted development environment.

## Citation and license

Lai, K. (2025). *AMRScan: A hybrid R and Nextflow toolkit for rapid antimicrobial
resistance gene detection.* [arXiv:2507.08062](https://arxiv.org/abs/2507.08062).
Also cite each caller and database used in an analysis. Historical source code
is under the repository MIT license; the NCBI test fixture has its own public
domain notice in `tests/fixtures/NCBI-LICENSE.txt`.

# AMRFinderPlus assembly validation

This layer is designed to retrieve only explicitly listed *Escherichia coli*
assemblies, verify their NCBI Assembly/BioSample/organism identity, retain input
hashes and timestamps, and then connect genomic calls to source AST records
without interpreting one as the other. A detected determinant is not a
phenotype; no detected determinant does not establish susceptibility. AST
values remain submitter/source observations. This work does not apply MIC
breakpoints, make predictions, or evaluate genotype and AST jointly.

## Cohort status and selection

`benchmark/cohorts/ecoli_amrfinder_validation.tsv` defines the exact fields for
the intended small cohort, but currently contains only its header. No isolate
is admitted until its versioned PDT accession resolves to exactly one BioSample
and Assembly accession and the NCBI Assembly report independently confirms
that BioSample and the exact organism. Candidate AST observations must be
retained as source records; duplicate/conflicting observations are not
collapsed. Missing, multiple, or incompatible mappings are excluded and
reported, never guessed. Selection is based on availability and identity, not
on AMRFinderPlus calls or agreement with AST.

PR #3 documented a live check of BioSample `SAMN05170351`: its AST rows linked
to SRA `SRR4065706`, but no Assembly accession. During this change, external
NCBI hosts could not be resolved from the development environment. Therefore
no accession-verified AST-linked genomes could be selected, and the manifest
is intentionally unpopulated. This PR does not claim real NCBI assembly
validation.

## Retrieval and identity checks

The standard-library downloader uses NCBI Datasets over HTTPS and only handles
assemblies named in the manifest. It extracts the accession-specific genomic
FASTA and assembly data report from the same package, verifies the reported
Assembly accession, BioSample and organism, validates FASTA structure, and
records the package URL, retrieval time, manifest hash, and FASTA SHA-256.
Downloads use temporary files and refuse to overwrite existing results.

```sh
python3 bin/genome_dataset.py download \
  --manifest benchmark/cohorts/ecoli_amrfinder_validation.tsv \
  --outdir datasets/ecoli_validation/genomes
python3 bin/genome_dataset.py validate \
  --manifest benchmark/cohorts/ecoli_amrfinder_validation.tsv \
  --genomes datasets/ecoli_validation/genomes
```

The current header-only manifest deliberately fails with “no validation
isolates”; this prevents an empty run from looking like a successful retrieval.
Genome data and provenance outputs are ignored by Git. The network transport
and metadata parsing are covered with deterministic offline mocks; NCBI HTTP
retrieval has not run here.

## Caller environment and database

`environment/amrfinderplus.yml` pins the intended caller package at 4.2.7 and
Python 3.12. Create it with `conda env create -f environment/amrfinderplus.yml`
or the equivalent Mamba command. The environment is a software definition,
not a database snapshot. Install an AMRFinderPlus database separately using
the pinned environment's `amrfinder_update -d /path/to/db` command, record the
database release/date reported by the caller and preserve the complete
directory fingerprint emitted by `bin/amrscan.py call`. Do not run database
updates as part of a workflow analysis. Copy the selected database to an
immutable experiment-specific directory: the same software can produce
different calls against different database releases.

AMRFinderPlus is not installed in the development environment. Consequently,
the `Escherichia` organism option, software installation, database setup, and
real caller execution have not been exercised in this change. The existing
`tests/fixtures/amrfinder_v4.0.23.tsv` is an upstream expected-output report;
it is parser input, not a positive-control genome or evidence that the caller
ran. `tests/fixtures/tiny.fna` is synthetic software-test input only; an empty
report from it would not establish a biological negative or susceptibility.
No biological sequence was added or fabricated as a control.

## Implemented and untested boundary

The downloader's schema/accession checks, archive parsing, identity failure
paths, no-overwrite behavior, checksum validation, and FASTA checks have
offline unit coverage with mocked NCBI package bytes. Existing harmonization
tests preserve raw AMRFinderPlus fields independently of AST. This branch has
not performed a real caller run, a real NCBI download, cohort genotyping,
evidence linking, or a real-data summary. Those claims require verified
accessions and an installed caller/database, neither of which is available in
this environment. A small positive caller-control sequence and its expected
call remain to be sourced from an appropriately licensed upstream package.

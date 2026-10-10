# Future model-training interface

No model training or phenotype prediction is implemented in this PR.
Genotype features will consume schema-versioned harmonized calls plus caller
provenance. Observed AST measurements must arrive independently through
`training/ingestion/`; never create labels from gene presence or caller classes.

A future training entry point should accept:

- genotype table and immutable provenance references;
- AST table conforming to `schemas/ast.tsv`, its snapshot manifest and an explicit
  isolate mapping; the [dataset builder](../docs/ast-dataset-builder.md) supplies
  observed phenotypes and conservative endpoint eligibility;
- a frozen split manifest from `benchmark/`;
- a versioned feature/label policy, seed and model configuration.

It should return a model artifact, training-only transformations, prediction
table, calibration/metrics and environment manifest. Missing AST stays missing;
multiple assays, conflicting measurements and breakpoint versions require
explicit adjudication before a label is assigned. No patient data belongs in Git.

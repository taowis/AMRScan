# Future NCBI Pathogen Detection AST ingestion

This directory defines a contract, not an implemented downloader. A future
adapter will accept a versioned source snapshot plus accession-to-sample mapping
and emit the header contract in `schemas/ast.tsv`, a source manifest and rejected
records with reasons. It must record the exact source URL, retrieval time,
snapshot checksum and record identifier so every observation is traceable.

Keep MIC/zone values as strings with comparators, units, assay method and testing
standard/version. Preserve the source's interpretation separately; do not derive
resistance from a BLAST/AMRFinderPlus hit. Use `NA` where metadata is not supplied.
Retain repeated or conflicting assays and unresolved isolate links for review.
The exact NCBI source format and mappings must be verified during implementation.

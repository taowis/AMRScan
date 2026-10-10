# AMRFinderPlus assembly validation

This validation layer keeps source AST observations separate from genomic
caller evidence. It does not infer phenotype from genes, predict R/S, apply
CLSI/EUCAST breakpoints or calculate genotype–phenotype performance. A caller
no-hit result does not establish susceptibility.

## Cohort selection and retrieval

`benchmark/cohorts/ecoli_amrfinder_validation.tsv` defines the intended small
cohort schema. It remains header-only until each exact versioned PDT accession,
BioSample and Assembly relationship is verified and the Assembly report
confirms BioSample and *Escherichia coli*. Candidate AST rows come from the
PR #3 evidence layer and are retained at observation level, including
conflicting and nonbinary categories. Inclusion is based on source availability
and identity, never on AMRFinderPlus calls or agreement with AST. Missing or
ambiguous mappings are excluded with an explicit reason, not guessed.

PR #3's documented BioSample `SAMN05170351` had AST linked to SRA
`SRR4065706` but no Assembly. During this change, two small official Pathogen
Detection `action=retrieve` lookups briefly succeeded on 2026-10-11 UTC and
showed candidate assembly-linked records:

| PDT | BioSample | Assembly | Observation |
| --- | --- | --- | --- |
| `PDT000026734.2` | `SAMN02581401` | `GCA_000633675.2` | First retrieved AST page contained ND observations. |
| `PDT000041778.1` | `SAMN03075588` | `GCA_000770275.1` | First retrieved AST page contained conflicting submitter observations. |

These preliminary lookup responses were not archived as a complete AST
snapshot or independently checked against Assembly reports; neither record is
admitted to the cohort. Subsequent requests failed DNS resolution, so no
complete NCBI AST snapshot or assembly FASTA was published.
The AST request used the official [Pathogen Detection service](https://www.ncbi.nlm.nih.gov/pathogens/pathogens-srv/)
with `collection=ast`, filter `scientific_name:"Escherichia coli"`, and
`limit=20`; the metadata lookup used `collection=isolates` with the two exact
PDT accessions and projected `target_acc,biosample_acc,asm_acc,scientific_name`.
The full table projections and export contract are in
[`schemas/ast_raw_schema.md`](../schemas/ast_raw_schema.md).

The downloader retrieves only explicitly listed accessions from NCBI Datasets,
then checks the returned Assembly accession, BioSample, exact organism, FASTA
structure and SHA-256. Its provenance records the manifest hash, source URL,
retrieval time, package member and FASTA hash. It refuses to overwrite prior
outputs. Genome files and provenance are ignored by Git.

```sh
python3 bin/genome_dataset.py download \
  --manifest benchmark/cohorts/ecoli_amrfinder_validation.tsv \
  --outdir datasets/ecoli_validation/genomes
python3 bin/genome_dataset.py validate \
  --manifest benchmark/cohorts/ecoli_amrfinder_validation.tsv \
  --genomes datasets/ecoli_validation/genomes
```

Once candidate accessions are selected, archive official AST evidence through
the existing builder. The query is the live NCBI browser export interface and
may return a large, changing table; inspect the resulting PDT/BioSample/
Assembly mappings before editing the cohort manifest.

```sh
python3 bin/ast_dataset.py download \
  --outdir datasets/ecoli_validation/ast_raw \
  --query 'scientific_name:"Escherichia coli"' --max-bytes 536870912
python3 bin/ast_dataset.py build \
  --snapshot datasets/ecoli_validation/ast_raw \
  --outdir datasets/ecoli_validation/ast_processed
```

## Caller environment, database and control

`environment/amrfinderplus.yml` pins AMRFinderPlus 4.2.7 and Python 3.12.
The environment could not be created here: Conda could not resolve
`conda.anaconda.org`, and its terms cache is outside the writable workspace.
No installed AMRFinderPlus version or database release/fingerprint is
available. The official NCBI usage documentation lists `Escherichia` as the
organism argument (`-O`); this was checked in the documentation, not against
an installed 4.2.7 binary. Verify it locally with `amrfinder -l` before
calling the cohort.

On a network-enabled host, create the pinned environment and download a
separate database snapshot. Keep the generated dated release directory
immutable; pass that exact directory to `--database` for controls and cohort
runs. `amrfinder --database_version` reports its software and database
versions, while the call provenance records every database file hash.

```sh
conda env create -f environment/amrfinderplus.yml
database_root="$PWD/datasets/amrfinderplus-db"
mkdir -p "$database_root"
conda run -n amrscan-amrfinder-validation amrfinder_update -d "$database_root"
find "$database_root" -mindepth 1 -maxdepth 1 -type d -print
database="<choose-the-dated-release-directory-printed-above>"
conda run -n amrscan-amrfinder-validation amrfinder --database "$database" \
  --database_version
conda run -n amrscan-amrfinder-validation amrfinder -l
```

NCBI's AMRFinderPlus installation test is a reproducible positive **software
control** from the upstream `ncbi/amr` repository at tag `amrfinder_v4.2.7`:
`test_dna.fa`, `test_prot.fa`, `test_prot.gff` and `test_both.expected`. The
expected output includes known AMR elements such as `blaTEM-156`, `blaPDC-114`
and `vanG`. NCBI identifies the software/database repository as public domain.
This is not an E. coli isolate, cohort member or AST observation. Retrieve
The pinned upstream files are [`test_dna.fa`](https://raw.githubusercontent.com/ncbi/amr/amrfinder_v4.2.7/test_dna.fa),
[`test_prot.fa`](https://raw.githubusercontent.com/ncbi/amr/amrfinder_v4.2.7/test_prot.fa),
[`test_prot.gff`](https://raw.githubusercontent.com/ncbi/amr/amrfinder_v4.2.7/test_prot.gff)
and [`test_both.expected`](https://raw.githubusercontent.com/ncbi/amr/amrfinder_v4.2.7/test_both.expected).
The [NCBI usage guide](https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus)
documents this installation test, and the [NCBI AMRFinderPlus repository](https://github.com/ncbi/amr)
contains its public-domain notice. Retrieve those four upstream files into an
ignored scratch directory, preserve each source URL and SHA-256, run the
installation test with the selected immutable database, and compare output
with the expected file. Record
`amrfinder --database_version` and the database file fingerprint. This control
was **not run** here. The tracked v4.0.23 expected-report fixture is parser
input only; the tiny synthetic FASTA is not a biological negative control.
No-hit must not be reported as susceptibility.

On a host with the pinned environment and database, the upstream control can
be retrieved and checked without adding sequence data to Git:

```sh
control_dir=/private/tmp/amrfinderplus-control-4.2.7
mkdir -p "$control_dir"
for name in test_dna.fa test_prot.fa test_prot.gff test_both.expected; do
  curl -fLsS "https://raw.githubusercontent.com/ncbi/amr/amrfinder_v4.2.7/$name" \
    -o "$control_dir/$name"
done
shasum -a 256 "$control_dir"/* > "$control_dir/SHA256SUMS"
amrfinder --database "$database" --plus -n "$control_dir/test_dna.fa" \
  -p "$control_dir/test_prot.fa" -g "$control_dir/test_prot.gff" \
  -O Escherichia --print_node -o "$control_dir/test_both.actual"
diff -u "$control_dir/test_both.expected" "$control_dir/test_both.actual"

# Software no-hit smoke only; this synthetic sequence says nothing about AST.
amrfinder --database "$database" -n tests/fixtures/tiny.fna \
  -O Escherichia -o "$control_dir/no_hit.tsv"
test "$(wc -l < "$control_dir/no_hit.tsv" | tr -d ' ')" = 1
```

For every verified cohort genome, run and harmonize the caller separately;
use the exact file names below, with `pdt`, `fasta` and `database` set from the
manifest, retrieval provenance and immutable database location.

```sh
amrfinder --database_version
python3 bin/amrscan.py call --sample "$pdt" --input "$fasta" \
  --database "$database" --organism Escherichia \
  --log "datasets/ecoli_validation/calls/$pdt.log" \
  --raw "datasets/ecoli_validation/calls/$pdt.amrfinder.tsv" \
  --provenance "datasets/ecoli_validation/calls/$pdt.provenance.json"
python3 bin/amrscan.py harmonize --sample "$pdt" \
  --raw "datasets/ecoli_validation/calls/$pdt.amrfinder.tsv" \
  --provenance "datasets/ecoli_validation/calls/$pdt.provenance.json" \
  --output "datasets/ecoli_validation/calls/$pdt.harmonized.tsv"
```

## Evidence linkage and summary

`bin/amrfinder_validation.py summarize` links exact PDT → BioSample → verified
Assembly → caller result → harmonized determinant evidence → AST observations.
It verifies the downloaded assembly hashes, caller input hash, raw report hash,
sample ID, caller version and database fingerprint. AST rows are joined only
by exact versioned PDT accession, then checked against the cohort BioSample
and organism. Assembly conflicts fail closed. Genotype and AST are stored in
separate evidence tables with shared isolate keys; every AST record is kept.

The deterministic summary reports selected cohort size, valid assemblies,
successfully genotyped isolates, isolates with and without AMR determinant
calls, AST-linked isolates, exclusions and reasons, AMR determinant/gene counts,
caller element types, and AST observations grouped by the normalized submitted
phenotype, including ND/I/SDD/NS/HLAR/NA where present. No sensitivity,
specificity, F1, PPV, NPV, prediction or concordance is calculated. The
summary's provenance JSON fingerprints its manifest, genome provenance, AST
table and caller reports. Results are published to a new directory and never
overwrite an existing summary.

```sh
python3 bin/amrfinder_validation.py summarize \
  --manifest benchmark/cohorts/ecoli_amrfinder_validation.tsv \
  --genomes datasets/ecoli_validation/genomes \
  --ast datasets/ecoli_validation/ast_processed/ecoli_ast.tsv \
  --calls datasets/ecoli_validation/calls \
  --exclusions benchmark/cohorts/ecoli_amrfinder_exclusions.tsv \
  --outdir datasets/ecoli_validation/summary
```

## Validation status

| Stage | Status |
| --- | --- |
| Retrieval, linkage and summary implementation | Implemented. Identity/hash mismatches fail closed. |
| Offline tests | Passed using mocked NCBI packages, synthetic FASTA, upstream report rows and AST fixtures. These validate software behavior only. |
| Official Pathogen Detection lookup | Two small candidate lookups succeeded briefly; candidates are documented above but not admitted. |
| Live AST snapshot and assembly FASTA retrieval | Not validated. Later network requests failed DNS; no source snapshot, assembly report or FASTA was published. |
| AMRFinderPlus 4.2.7 and database | Not installed; no database release/fingerprint. |
| Positive/negative real caller controls | Not run. The NCBI upstream positive software control is identified above; no biological E. coli negative control is claimed. |
| Real cohort, real genotyping and real-data summary | Not validated. The cohort manifest remains header-only; no real cohort counts/calls are claimed. |

Keep PR #4 in Draft until at least one accession-verified real NCBI assembly
has been retrieved and run through the real AMRFinderPlus caller, harmonized,
and included in an auditable provenance chain. Then complete the intended
3–10 isolate cohort and repeat the documented checks. Offline tests do not
meet this scientific completion criterion.

# Future benchmarking and generalisability interface

No benchmark scores or split algorithms are implemented in this PR.
Future evaluation takes frozen genotype/AST snapshots, an isolate mapping and a
split manifest with columns `sample_id`, `split`, `holdout_axis`, `holdout_value`,
`dataset_version`, `split_version`, `seed`. Allowed `split` values are train,
validation and test. Unknown metadata must be accounted for explicitly.

Required evaluation designs:

- **Lineage holdout:** keep all linked isolates/assemblies from a held-out lineage
  outside training and validation; specify the lineage typing method/version.
- **Country holdout:** freeze geographic definitions and handle unknown or
  conflicting provenance separately.
- **Time holdout:** set a cutoff using collection date, preserve date precision,
  and exclude missing/ambiguous dates or evaluate them as a documented cohort.

Prevent identical isolates, near duplicates and linked assemblies from crossing
splits. Fit imputation, scaling and feature selection on training data only.
Keep the test set untouched during tuning, report cohort sizes and exclusions,
per-drug class balance, calibration and uncertainty, and compare to simple
baselines. Publish split manifests and their hashes with caller/database and
AST breakpoint versions so changes in data/callers are distinguishable from
changes in model performance. Random splits alone do not establish transportability.

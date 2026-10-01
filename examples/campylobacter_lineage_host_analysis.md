# Campylobacter accessory methylome: lineage versus host analysis

This analysis addresses the central population-level MOMENTO question for the current 52-isolate Campylobacter panel:

> Do isolates possess accessory methylome repertoires that track genomic lineage, and is there evidence for an additional host/ecology association once lineage is accounted for?

The analysis is deliberately framed as **variance attribution/association**, not causation. Host, lineage and study provenance are partially confounded in this small historical panel.

## Inputs

The script requires only four frozen inputs:

```text
1. IQ-TREE Newick tree
2. accessory_presence_absence.tsv
3. campylobacter_methylotypes_metadata_v0.2.tsv
4. an output directory
```

It does **not** rerun Prokka, Panaroo, IQ-TREE or methylation calling.

Example:

```bash
COHORT=/xdisk/cooperma/bpascoe/CampyMethylome/03_DERIVED_RESULTS/MOMENTO_RESULTS/COHORT_V0.1c
WORK=$COHORT/CORE_PHYLOGENY_V0.1
REPO=/xdisk/cooperma/bpascoe/software/MOMENTO-MethyLome

python "$REPO/scripts/analyse_methylome_lineage_host.py" \
  --tree "$WORK/iqtree/campy52_core.treefile" \
  --methylotype-dir "$COHORT/METHYLOTYPES_V0.1" \
  --metadata "$REPO/metadata/campylobacter_methylotypes_metadata_v0.2.tsv" \
  --output-dir "$COHORT/METHYLOTYPES_V0.1/LINEAGE_HOST_V0.1"
```

## 1. Accessory methylome distance

The response matrix is the complete accessory motif-family presence/absence matrix. For isolates `i` and `j`, the script calculates Jaccard distance:

```text
Jaccard distance = 1 - |features_i ∩ features_j| / |features_i ∪ features_j|
```

This is appropriate for a sparse binary accessory repertoire because shared absences do not make two isolates artificially similar.

Output:

```text
accessory_methylome_jaccard_distance.tsv
```

## 2. Genomic-lineage distance

The script calculates pairwise **patristic distance** from the core-genome IQ-TREE tree: the sum of branch lengths connecting each pair of tips.

This measure is independent of the arbitrary display root of the Newick file.

Output:

```text
core_phylogeny_patristic_distance.tsv
```

## 3. Do genomic and methylome distances track one another?

The script compares the upper triangles of the two distance matrices using both Pearson and Spearman correlations and assesses them with sample-label permutations (999 by default).

Output:

```text
lineage_methylome_distance_correlation.tsv
```

This is a direct answer to the first question: are genetically more distant isolates also more different in their accessory methylome?

A significant positive relationship is evidence of lineage structure in the methylome repertoire. It does **not** imply that every methylation system is vertically inherited; horizontally acquired or phase-variable systems can still create substantial departures from the trend.

## 4. Representing lineage for variance partitioning

A tree-distance matrix cannot be inserted directly into an ordinary design matrix. The script therefore performs principal coordinates analysis (PCoA) on patristic distance and uses the leading positive axes as continuous lineage covariates.

Default rule:

```text
retain enough positive phylogenetic PCoA axes to represent >=80% of
positive-eigenvalue variation, with a maximum of 10 axes
```

The exact number retained is written to:

```text
analysis_settings.tsv
```

This matters for reproducibility: changing the lineage representation can change the apparent residual host effect.

## 5. Which host groups enter the primary test?

The complete metadata remain in descriptive outputs, but the primary host comparison requires at least three isolates per group by default.

Current counts are:

```text
human        31
wild_mammal   8
wild_bird     6
poultry       5
environment   1
unknown       1
```

Therefore the primary host-lineage model uses:

```text
human + wild_mammal + wild_bird + poultry
```

The single environment isolate and unresolved isolate are retained in the distance outputs but excluded from the inferential host partition. With one observation, `environment` has no estimable within-group variation.

Output:

```text
host_group_counts.tsv
analysis_samples.tsv
```

## 6. Distance-based R² models

The script calculates distance-based R² from the Gower-centred accessory Jaccard matrix for the following models:

```text
host
lineage
host + lineage
provenance
lineage + provenance
host + lineage + provenance
```

Both raw and adjusted R² are written. Adjusted R² is particularly useful here because lineage is represented by several continuous axes and host/provenance by multiple dummy variables.

Output:

```text
distance_variance_models.tsv
```

## 7. The headline quantities

The most useful contrasts are written to:

```text
variance_partition_components.tsv
```

They include:

```text
unique_host_after_lineage
unique_lineage_after_host
shared_host_lineage
unique_host_after_lineage_and_provenance
```

The first quantity answers:

> How much additional accessory-methylome variation is associated with host after genomic lineage has already been represented?

The provenance-adjusted quantity is a stricter sensitivity analysis. Because several collections are dominated by particular hosts, host and provenance can be difficult to disentangle completely in this dataset.

### Negative adjusted components

Adjusted-R² differences can occasionally be slightly negative. Do not clamp them to zero silently. A negative estimate is useful information: after penalising model complexity, that component does not add detectable explanatory power in this sample.

## 8. What this analysis does not prove

Even a non-zero host component would not by itself demonstrate host-driven methylome adaptation. Potential explanations include:

- genuine ecological selection;
- residual genomic population structure not captured by the retained axes;
- study/collection effects;
- geographically structured sampling;
- correlated restriction-modification system turnover;
- phase-variable methyltransferases;
- horizontal gene transfer.

The correct next step after a population-level signal is to identify which motif families drive it and map those motifs to candidate MTase/RM loci.

## 9. Sensitivity analyses to add if the signal is interesting

Recommended follow-ups:

1. repeat with a recombination-aware core-genome tree;
2. vary the number of lineage PCoA axes;
3. repeat within well-sampled clonal complexes;
4. test individual motif families with lineage-aware models;
5. project cognate MTase/RM gene presence/absence alongside motif state;
6. examine phase-variable loci separately from stable gene gain/loss;
7. check whether conclusions survive removal of large provenance blocks one at a time.

## 10. Reproducibility rule

Every result directory should retain:

```text
analysis_settings.tsv
analysis_samples.tsv
host_group_counts.tsv
core_phylogeny_patristic_distance.tsv
accessory_methylome_jaccard_distance.tsv
lineage_methylome_distance_correlation.tsv
distance_variance_models.tsv
variance_partition_components.tsv
```

Together with the Git commit SHA and metadata version, those files are enough to reconstruct exactly which samples, variables and distance definitions produced the reported percentages.

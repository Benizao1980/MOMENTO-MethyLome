# Campylobacter methylotype plotting workflow

This example starts from a QC-clean MOMENTO cohort and the outputs of
`momento methylotypes`.

## 1. Generate methylotype analysis outputs

```bash
momento methylotypes \
  --input cohort.momento.tsv \
  --output-dir METHYLOTYPES_V0.1
```

The plotting workflow expects the standard files written by this command,
including `accessory_presence_absence.tsv`, `pcoa.tsv`, `pcoa_variance.tsv`,
`linkage.tsv`, and `motif_family_prevalence.tsv`.

## 2. Curate metadata

Start from:

```text
metadata/campylobacter_methylotypes_metadata_v0.1.tsv
```

The starter table contains the 52 `PASS` Campylobacter isolates used in the
first real-data analysis. Some provenance fields are already filled from known
source information; rows explicitly marked `unknown` or `metadata to curate`
should not be interpreted as complete biological metadata.

For a host-ecology analysis, curate at minimum:

- `host_group`
- `host_detail`
- `provenance_class`
- `species` / `subspecies`
- `st`
- `clonal_complex`

Do not infer an ecological source solely from membership in a historical
reference panel.

## 3. Generate Life Aquatic figures

From the cohort results directory:

```bash
python /path/to/MOMENTO-MethyLome/scripts/plot_methylotypes_life_aquatic.py \
  --methylotype-dir METHYLOTYPES_V0.1 \
  --metadata /path/to/MOMENTO-MethyLome/metadata/campylobacter_methylotypes_metadata_v0.1.tsv \
  --cohort-summary cohort.momento.tsv \
  --output-dir METHYLOTYPES_V0.1/FIGURES_LIFE_AQUATIC
```

Default outputs are written as PNG, PDF and SVG:

```text
figure_accessory_heatmap.*
figure_accessory_heatmap_all_features.*
figure_pcoa.*
figure_accessory_prevalence.*
figure_raatty_backbone.*
heatmap_sample_order.tsv
heatmap_feature_order.tsv
heatmap_all_feature_order.tsv
```

The main heatmap defaults to accessory families present in at least two PASS
isolates. This is only a display filter: the dendrogram still comes from the
full Jaccard/accessory matrix. The all-feature heatmap is written separately so
singleton families remain available for supplementary review.

The default styling is documented in:

```text
docs/figure_style_life_aquatic.md
```

## 4. Useful plotting options

Change the main-heatmap prevalence display threshold:

```bash
--heatmap-min-prevalence 3
```

Skip the supplementary all-feature heatmap while testing:

```bash
--skip-full-heatmap
```

Label every PCoA point for a diagnostic version:

```bash
--label-points
```

Change the PCoA colour / shape metadata:

```bash
--color-by host_group --shape-by provenance_class
```

Change heatmap annotation strips:

```bash
--annotations host_group,provenance_class,st,clonal_complex
```

The publication default is deliberately only:

```text
host_group,provenance_class
```

because the current ST/CC metadata are incomplete.

Disable the shape mapping:

```bash
--shape-by ""
```

Restrict output formats while testing:

```bash
--formats png
```

## Interpretation notes

- The accessory matrix is `PASS`-only and excludes the core RAATTY family.
- Reverse-complement equivalent motif labels are already collapsed by the
  methylotype workflow.
- Nested motifs are retained because the real-data site audit did not identify
  strong same-GFF-record representation redundancy.
- The first 52-isolate analysis contains 50 exact binary accessory profiles;
  therefore the PCoA / clustering should be interpreted as a high-dimensional
  repertoire rather than evidence for a few discrete universal methylotypes.

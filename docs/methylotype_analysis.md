# Accessory methylotype analysis

`momento methylotypes` performs a conservative first-pass population analysis of a canonical MOMENTO cohort table.

The default workflow is designed for sparse bacterial methylome repertoires:

1. keep only samples with `qc_status=PASS`;
2. keep `m6A` and `m4C` as separate modification classes;
3. collapse a motif and its reverse complement to one stable motif-family label;
4. retain motif nesting as an audit flag rather than automatically deleting narrower motifs;
5. exclude the nominated core motif family (`RAATTY` by default) from the accessory matrix;
6. build a binary isolate-by-accessory-family matrix;
7. calculate Jaccard distances;
8. perform classical PCoA and average-linkage hierarchical clustering.

This command intentionally starts with presence/absence rather than methylated fractions. Family-level occupancy needs additional care because reverse-complement summary rows may share genomic opportunity denominators and partial GFFs can otherwise mimic low biological occupancy.

## Run

```bash
momento methylotypes \
  --input cohort.momento.tsv \
  --output-dir methylotypes_v0.1
```

Useful options:

```text
--qc-statuses PASS
--modifications m6A,m4C
--exclude-families RAATTY
--min-called-sites 1
```

`PARTIAL_GFF`, `LOW_DATA` and `FAIL` samples are excluded by the default `PASS` filter.

## Reverse-complement families

Examples such as

```text
CTANNNNNNNGTC <-> GACNNNNNNNTAG
GCANNNNNRTTA  <-> TAAYNNNNNTGC
```

are represented once in the binary matrix. This prevents opposite-strand representations of one recognition specificity from being double weighted in clustering.

## Nested motif audit

IUPAC containment is reported in `nested_motif_pairs.tsv`. For example, `AAAATTY` is sequence-compatible with a narrower subset of the core `RAATTY` recognition language. MOMENTO does **not** automatically delete the child motif: nesting can be a caller representation artefact, but it can also reflect a genuine additional methyltransferase specificity.

The audit table therefore reports parent/child prevalence and co-occurrence, including `child_with_parent_fraction` and `parent_with_child_fraction`. Site-level evidence or MTase annotation should be used before deciding that a nested feature is redundant.

## Outputs

```text
motif_families_long.tsv
motif_family_prevalence.tsv
nested_motif_pairs.tsv
accessory_presence_absence.tsv
jaccard_distance.tsv
pcoa.tsv
pcoa_variance.tsv
cluster_order.tsv
linkage.tsv
```

`accessory_presence_absence.tsv` is the main matrix for first-pass methylotype analysis. Feature names are `modification|motif_family`, so m6A and m4C versions of the same recognition string remain distinct.

`pcoa.tsv` and `cluster_order.tsv` are intended to be joined to isolate metadata (host, species, ST/CC, lineage, geography, platform) rather than interpreted without genomic/population context.

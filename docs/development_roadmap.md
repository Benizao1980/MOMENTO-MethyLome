# Roadmap

## v0.1 — reproduce the current Campylobacter analysis
- stable long-table schema;
- context-aware PacBio GFF parser;
- cognate modified-position inference with auditable site-level filtering;
- FASTA methylatable-target opportunity counter;
- automatic motif summary;
- optional caller-confidence filtering without assuming a universal threshold;
- QC, heatmap, PCA/Jaccard/PCoA;
- validate against published RM1245/RM1477 methylomes and the existing curated spreadsheet.

Acceptance tests:
1. semantically reproduce the published RM1245/RM1477 methylomes from primary GFF + FASTA, investigating version-dependent discrepancies rather than forcing exact equality;
2. rebuild the historical Campylobacter motif matrix without manual Excel editing.

## v0.2 — mechanism
MTase annotation, motif↔MTase linking, gene/promoter mapping, candidate phase-variable tracts.

## v0.3 — population genomics
Lineage overlays, PERMANOVA/variance partitioning, grouped CV, genome vs methylome vs combined models.

## v0.4 — ONT
bedMethyl import and matched-isolate cross-platform validation.

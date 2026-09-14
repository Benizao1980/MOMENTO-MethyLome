# Architecture

## Layers
1. **Import adapters**: PacBio GFF/CSV, ONT bedMethyl/modBAM.
2. **Canonical long table**: platform-independent representation.
3. **Comparative methylomics**: QC, motif opportunities, matrices, ordination, clustering.
4. **Mechanism**: MTase/R-M annotation, phase-variable loci, genomic context.
5. **Population analysis**: metadata, lineage-aware association and grouped CV.

The key design decision is that downstream biology never depends directly on PacBio- or ONT-specific columns.

## Main tables
- `methylation_long.tsv`: sample × motif summaries.
- `samples.tsv`: biological/population metadata.
- `mtases.tsv`: candidate MTase/R-M systems.
- `phase_loci.tsv`: candidate phase states.
- `site_calls.tsv`: optional large position-level table, generally outside Git.

## Cross-platform rule
Analyse PacBio and ONT independently first. Only pool after matched-isolate concordance and platform/caller effects have been assessed.

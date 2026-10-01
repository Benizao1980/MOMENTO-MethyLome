# MOMENTO positioning and related tools

Working notes for software/paper positioning. These are deliberately concise and should be updated as the literature review matures.

## Closest conceptual precedents

### MethylomeMiner / PanMethylomeMiner

MethylomeMiner is an important direct comparator for bacterial comparative methylomics, especially for ONT-derived data. Its population layer is gene/locus and pangenome oriented: methylated sites are associated with coding/intergenic regions and then projected across a pangenome.

MOMENTO should not be framed as the first bacterial pan-methylome tool. The useful distinction is conceptual rather than competitive:

- **MethylomeMiner / PanMethylomeMiner:** which genomic loci or genes are methylated across a bacterial pangenome?
- **MOMENTO:** which methylation systems/states does each isolate possess, how active are they, and how are those states distributed across a population?

The two approaches are complementary and should ultimately be interoperable where practical.

### Population pan-methylome / panepigenome studies

Recent population-scale bacterial methylome studies provide useful precedents for treating methylation as a population-level repertoire rather than only a per-genome annotation problem. Two themes are especially relevant to MOMENTO:

- a conserved/core methylation component can coexist with a large accessory repertoire;
- exact methylation profiles may be highly individualized, so a pan-methylome need not resolve into a few discrete universal methylotypes.

This fits the current Campylobacter result: a conserved RAATTY backbone with a sparse, highly combinatorial accessory methylome.

## MOMENTO's intended niche

MOMENTO is being designed as a **species-agnostic comparative bacterial methylomics framework downstream of methylation callers**.

Key differentiators to preserve:

1. **Caller/platform agnostic data model** rather than a workflow tied to one sequencing platform.
2. **Motif- and methylation-system-centric analysis**, including motif family presence/absence, genomic opportunities, occupancy and score/QC information.
3. **Population-aware analysis**, with explicit separation of methylome structure from genomic lineage.
4. **MTase / restriction-modification integration** as a first-class analysis layer rather than an afterthought.
5. **Phase-variable methylation state** as an explicit target.
6. **Reproducible QC and provenance**, including incomplete/truncated call-file detection.
7. **Cross-platform harmonization** for PacBio and ONT at the standardized-summary level.
8. **Gene/pangenome projection** as a planned downstream layer, complementing rather than replacing motif/system-centric analysis.

A useful hierarchy for the software and eventual manuscript is:

> site -> motif -> methylation system -> isolate methylome profile -> population pan-methylome -> gene/function/ecology

## Terminology

- **Comparative methylomics** is the safest broad description of the field.
- **Pan-methylome / panepigenome** is appropriate for population-scale analysis of shared and accessory methylation features.
- **Methylotype** should be used as a multivariate methylome profile/state, not assumed to mean a small set of discrete classes.
- **Meta-methylomics** is potentially useful for future metagenomic extensions, but is less standardized terminology than metagenomic methylation profiling / methylation-aware metagenomics.

## Roadmap implication

Add a formal gene/pangenome projection layer to the MOMENTO roadmap. This should consume standardized MOMENTO site/motif output plus external genome annotation/pangenome mappings and produce gene-/locus-level methylation summaries without making the core framework dependent on any one pangenome package.

Candidate future command/module names:

- `momento functional`
- `momento project-pangenome`
- `src/momento/functional.py`
- `src/momento/pangenome.py`

The immediate development priority remains lineage-aware comparative methylomics: core-genome phylogeny, aligned methylome tracks and population-structure-aware association testing.

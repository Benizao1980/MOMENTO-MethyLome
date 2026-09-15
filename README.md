# MOMENTO

**MOMENTO** is a species-agnostic toolkit for comparative bacterial methylomics.

It is designed for studies in which methylation has already been called by a sequencing-platform-specific tool (for example PacBio SMRT Analysis or Oxford Nanopore modified-base calling) and the biological question is comparative:

- Which methylation motifs are conserved or variable across isolates?
- Are there reproducible methylation states ("methylotypes")?
- Which methyltransferases (MTases) or restriction-modification systems explain those methylotypes?
- Are phase-variable loci associated with switching methylation states?
- Does methylation associate with host, ecology, disease, phenotype, or population structure?
- Does methylation add predictive information beyond the genome itself?

MOMENTO deliberately does **not** aim to replace raw-signal methylation callers. Instead, it provides a common representation and reproducible downstream analysis layer for bacterial population methylomics.

The name is a subtle nod to *Memento*: bacterial methylation can encode a molecular state or "memory", while the software reconstructs that state across populations.

## Design principles

1. **Species agnostic** — no organism-specific assumptions in the core package.
2. **Platform aware** — PacBio and ONT inputs are imported separately, then harmonised.
3. **Long-table first** — all importers produce one tidy internal schema.
4. **Population-aware** — host/phenotype association should account for genomic lineage.
5. **Interpretable** — motifs, MTases, phase-variable loci, genes and genomic positions remain traceable.
6. **Reproducible** — analyses should be rerunnable from primary caller outputs, not hand-edited spreadsheets.
7. **Teaching friendly** — examples explain *why* each step is performed, not only which command to type.

## Internal data model

The central table has one row per isolate × motif summary, with fields such as `sample_id`, `platform`, `motif`, `modification`, `called_sites`, `genomic_sites`, `methylated_fraction`, `median_score`, `source_file`, and `qc_status`.

## Intended workflow

```text
PacBio motif GFF / modification CSV ─┐
                                    ├─> MOMENTO import
ONT modBAM / bedMethyl ─────────────┘
                                           |
                                           v
                                harmonised long table
                                           |
                    ┌──────────────────────┼──────────────────────┐
                    v                      v                      v
                   QC                 motif matrix          annotation
                    |                      |                 MTase / phase
                    └──────────────┬───────┴───────────────┬─────┘
                                   v                       v
                              methylotypes            associations
                          PCA / PCoA / clustering    host / phenotype
```

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
momento --help
```

## First working importer: PacBio

`momento import-pacbio` converts a PacBio methylation-site GFF plus the matching assembly FASTA into MOMENTO's canonical motif summary. The assembly is used to count genomic motif opportunities, so raw called-site counts can be converted to methylated fractions.

```bash
momento import-pacbio \
  --gff RM1245motifs.gff \
  --fasta CjRM1245.fasta \
  --sample RM1245 \
  --output RM1245.momento.tsv \
  --sites-output RM1245.sites.tsv
```

The importer supports IUPAC motifs such as `RAATTY`, handles reverse complements, preserves site-level GFF attributes, and fails loudly if an explicit motif field cannot be identified instead of guessing from sequence context.

## CLI roadmap

```bash
momento validate
momento import-pacbio   # implemented
momento import-ont      # planned
momento qc
momento matrix
momento cluster
momento mtase           # planned
momento phase           # planned
momento associate       # planned
momento report          # planned
```

The MVP now implements validation, PacBio GFF import, genomic motif-opportunity counting, matrix generation, basic core-motif QC, prevalence summaries, PCA and Jaccard utilities.

## v0.1 acceptance test

The first real-data validation target is to reproduce the published methylome summaries for the *Campylobacter jejuni* HS:19 strains **RM1245** and **RM1477** directly from primary PacBio methylation GFFs plus assembly FASTAs. This provides an external ground truth for motif opportunity counting and methylated fractions before scaling to the full cohort.

## Campylobacter teaching example

See `examples/campylobacter/README.md`. The example uses synthetic/public-safe data and explains:

- why a near-core motif such as `RAATTY` should be handled separately;
- how to distinguish failed methylome runs from genuine biological absence;
- how to define methylotypes;
- how to overlay host and genomic lineage;
- how to connect motifs to candidate MTases and phase-variable loci;
- how to avoid phylogenetic leakage in phenotype prediction.

For the current 57-genome Campylobacter panel, the core-genome phylogeny should be built separately from assembly FASTAs using a bacterial annotation workflow (for example Bakta or Prokka) followed by Panaroo/PIRATE and IQ-TREE. The PacBio `*motifs.gff` files are methylation-call GFFs, **not** gene-annotation GFFs and should not be supplied directly to Panaroo/PIRATE.

## Data policy

Do not commit raw reads, large modification CSVs, unpublished private assemblies, or sensitive metadata. Keep primary data in Box/HPC and commit code, schemas, public-safe derived tables, synthetic examples, documentation, and release-cleared figures.

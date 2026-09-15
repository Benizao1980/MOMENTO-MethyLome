# Worked example: comparative *Campylobacter* methylomics

This tutorial is intentionally detailed. It teaches the logic of population-scale bacterial methylome analysis rather than giving a black-box recipe. The toy data are synthetic and safe to publish.

## Central question

> Do *Campylobacter* isolates possess reproducible methylation states (methylotypes) that associate with host ecology independently of genomic lineage?

Mechanistic extension:

> Are variable methylotypes explained by the presence, absence or phase-variable state of methyltransferase systems?

## Current real-data panel

The working Campylobacter development cohort contains 57 PacBio methylomes assembled from several historical/public strain panels rather than one ecological sampling design. This includes Penner serotype reference strains, published HS:19/GBS strains, wildlife/environment isolates, black-bear isolates and a *C. jejuni* subsp. *doylei* genome. That heterogeneity is useful for developing MOMENTO, but host-association analyses must use carefully defined ecological subsets.

## Why lineage matters

Closely related bacteria share genome content and often share methylation systems. If a cattle-associated lineage carries a motif, a naive model may appear to predict cattle when it is only recognising the lineage. Host association therefore needs within-lineage evidence or grouped/phylogeny-aware validation.

For the 57-genome development panel, build a core-genome phylogeny separately from the assembly FASTAs:

```text
57 assembly FASTAs
      ↓
Bakta or Prokka
      ↓
Panaroo or PIRATE
      ↓
core-gene alignment
      ↓
IQ-TREE
      ↓
phylogeny / lineage metadata for MOMENTO
```

The PacBio `*motifs.gff` files are methylation-call GFFs, **not** the gene-annotation GFFs required by Panaroo/PIRATE.

## RAATTY as a teaching example

`RAATTY` can be near-core and extremely abundant in *C. jejuni*. That means:

1. binary presence/absence may contain little host-discriminating information;
2. raw counts can dominate PCA because the feature is numerically huge;
3. the motif can instead act as a valuable **QC sentinel**;
4. regional/position-level variation may still be biologically informative.

Run methylotype analyses both with all motifs and with a near-core motif removed. This tests whether variable secondary systems drive structure.

## Real project inputs

### Assembly FASTA
Used to count motif opportunities, locate methylated sites, annotate MTases/R-M systems, and inspect candidate repeat tracts.

### PacBio motif GFF
Position-level modification/motif output. The importer must be version-tolerant because SMRT Analysis exports vary.

### PacBio modification CSV
Potentially richer per-site scores. Keep these outside Git.

### Metadata
At minimum: sample, species, host, host group, location/date, ST/CC or another genomic lineage definition, platform.

### MTase table
Candidate enzyme, locus, R-M type, recognition motif, gene status, evidence/confidence.

## Why genomic motif opportunities matter

Raw called-site counts conflate genome composition with methylation state. Prefer, when possible:

```text
methylated_fraction = methylated calls / genomic motif opportunities
```

A strain with 24,000 called RAATTY sites may be ~100% methylated if only 24,100 opportunities exist, or partially methylated if 28,000 opportunities exist.

## Install

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

## Step 0: import a PacBio methylome

The first working MOMENTO adapter converts a PacBio methylation GFF plus the matching assembly FASTA into the canonical motif table:

```bash
momento import-pacbio \
  --gff RM1245motifs.gff \
  --fasta CjRM1245.fasta \
  --sample RM1245 \
  --output RM1245.momento.tsv \
  --sites-output RM1245.sites.tsv
```

The importer counts motif opportunities directly from the assembly, understands degenerate IUPAC motifs such as `RAATTY`, and retains the original GFF attributes in the optional site-level table.

### First real-data acceptance test

Before processing the full cohort, validate MOMENTO against the published RM1245 and RM1477 HS:19 methylomes. MOMENTO should approximately reproduce the published opportunity/call totals for motifs including `RAATTY`, `CATG` and `AGTNNNNNNRTTG` directly from primary FASTA + GFF inputs. Any discrepancy should be investigated before scaling to the remaining isolates.

## Step 1: validate

```bash
momento validate --input examples/campylobacter/data/toy_methylation_long.tsv
```

Before analysis, inspect identifier consistency, missing-value meaning, platform labels and motif representation. Never silently treat unknown as zero.

## Step 2: motif prevalence and core-motif QC

```bash
mkdir -p examples/campylobacter/results
momento qc \
  --input examples/campylobacter/data/toy_methylation_long.tsv \
  --prevalence-output examples/campylobacter/results/motif_prevalence.tsv \
  --core-motif RAATTY \
  --core-output examples/campylobacter/results/raatty_qc.tsv
```

Questions to ask:
- Which motifs are near-core?
- Which occur in intermediate subsets?
- Which are singleton/rare calls?
- Are rare calls real or low-support noise?

If a biologically expected core motif is almost absent in one isolate, investigate sequencing/calling/file matching before declaring a novel methylotype.

## Step 3: variable-motif matrix

```bash
momento matrix \
  --input examples/campylobacter/data/toy_methylation_long.tsv \
  --output examples/campylobacter/results/motif_presence_no_core.tsv \
  --value presence \
  --min-called-sites 100 \
  --exclude RAATTY
```

The `100` threshold is only a transparent teaching threshold. Production QC should use caller confidence, genomic opportunities, methylation fraction and platform validation.

## Step 4: compare representations

Build and compare:
1. presence/absence;
2. raw/log counts;
3. methylated fractions.

If biological conclusions disappear when representation changes, that is a warning that the result may be preprocessing-driven.

## Step 5: ordination

```bash
momento cluster \
  --matrix examples/campylobacter/results/motif_presence_no_core.tsv \
  --output examples/campylobacter/results/pca_presence_no_core.tsv
```

PCA is an accessible first pass. For binary motif states, Jaccard distance plus PCoA is often more natural. Ask whether clusters track host, lineage, species, platform, or failed samples.

## Step 6: add metadata

Join ordination/clustering results to `metadata/toy_samples.tsv`. Compare colouring by host with colouring by lineage. If host and lineage are completely confounded, no model can separate them reliably.

## Step 7: from motif to mechanism

For an interesting variable motif:
1. identify candidate MTase;
2. test MTase gene presence;
3. inspect gene integrity;
4. look for SSR/poly-G/poly-C tracts;
5. estimate provisional ON/OFF state;
6. compare genotype/state with motif activity;
7. map motif sites to genes/promoters.

Start from the methylation phenotype and work backward. Do not scan every repeat tract first and then fish for a story.

## Phase variation caveat

Assemblies may misrepresent homopolymer lengths. Treat assembly-only phase calls as hypotheses. Strong mechanistic claims should be supported with raw reads, especially for ON/OFF tract length.

## Step 8: host/source prediction

Only after QC and population-structure inspection. Use grouped cross-validation by clonal complex/lineage/cgMLST cluster rather than random isolate splits.

Compare:
- Model A: genome-only;
- Model B: methylome-only;
- Model C: genome + methylome.

The scientifically stronger question is whether methylation adds predictive information beyond genomic lineage.

## Position-level extension

Motif summaries are only the first layer. With position-level GFFs, build sample × ortholog/gene/promoter methylation summaries. For a conserved system like RAATTY, *where* methylation occurs may be more informative than simple motif presence.

## Real Campylobacter workflow

### Phase 1 — inventory
Reconcile sample names across FASTA, motif GFF, modification CSV, old spreadsheets, metadata, phylogeny/ST/CC and platform.

### Phase 2 — reproduce published/historical results
First reproduce RM1245/RM1477 published methylome summaries, then programmatically rebuild the old curated motif spreadsheet from primary outputs. These are MOMENTO acceptance tests.

### Phase 3 — QC
Use call counts, confidence, motif opportunities, core-motif behaviour, genome size and failed/empty-file checks.

### Phase 4 — methylotypes
Presence matrix, fraction matrix, heatmap, Jaccard/PCoA, PCA, clustering.

### Phase 5 — mechanism
Host + lineage + MTases + phase state + genomic targets.

### Phase 6 — prediction
Grouped CV and genome-vs-methylome-vs-combined models.

### Phase 7 — ONT validation
Analyse ONT separately first. Use matched isolates to quantify caller/platform effects before pooling.

## Exercises

1. Change the presence threshold; which methylotypes are stable?
2. Cluster with and without RAATTY; explain the change.
3. Compare binary states with methylated fractions.
4. Colour by host, then lineage; which explanation is stronger?
5. For one variable motif, draw: `genotype → MTase state → methylation → genomic targets → phenotype`, and state what evidence supports each arrow.

## Public-repository rule

Commit code, documentation, synthetic/public data and cleared figures. Do not commit private metadata, raw reads, unpublished assemblies or large caller outputs unless release is explicitly authorised.

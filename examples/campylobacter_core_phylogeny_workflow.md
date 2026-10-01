# Campylobacter core-genome phylogeny + methylome workflow

This document is the reproducible worked example for the current Campylobacter MOMENTO population analysis. It explains both **what** is run and **why** each stage is separated.

The current validated cohort contains 52 QC-PASS isolates. The core-genome tree is built independently of the methylome data and is then used as a population-genomic backbone for ordering and analysing accessory methylation patterns.

## Overview

```text
QC-clean assemblies (52)
        |
        v
consistent Prokka annotation
        |
        v
Panaroo strict pangenome + core alignment
        |
        v
IQ-TREE 2 model selection + ML tree
        |
        +------------------------------+
        |                              |
        v                              v
phylogeny-aligned figure        lineage/host analysis
        |                              |
ST + CC + host + provenance     methylome Jaccard distance
RAATTY occupancy                vs patristic distance
accessory motif repertoire      variance partitioning
```

The important design choice is that methylome clustering never determines tree order. Genomic lineage and methylome state remain distinct data layers.

---

## 1. Input cohort and directory layout

Current working cohort on Puma:

```bash
COHORT=/xdisk/cooperma/bpascoe/CampyMethylome/03_DERIVED_RESULTS/MOMENTO_RESULTS/COHORT_V0.1c
WORK=$COHORT/CORE_PHYLOGENY_V0.1
REPO=/xdisk/cooperma/bpascoe/software/MOMENTO-MethyLome
```

The validated cohort has 52 `PASS` isolates. The current exclusions are not silently discarded: they remain in cohort QC outputs and are excluded from the phylogeny because they did not meet the same input criteria.

Stable sample identifiers are critical. The same identifier must appear in:

- the MOMENTO accessory presence/absence matrix;
- the metadata table;
- staged assembly filenames;
- Prokka output prefixes;
- Panaroo alignment sequence labels;
- IQ-TREE tips.

That identity constraint is what allows all tracks to be aligned without manual row ordering.

---

## 2. Stage only the QC-clean assemblies

```bash
python "$REPO/scripts/prepare_phylogeny_inputs.py" \
  --manifest /xdisk/cooperma/bpascoe/CampyMethylome/01_ACTIVE_INPUTS/momento_campylobacter_57_manifest.tsv \
  --qc "$COHORT/cohort.qc.tsv" \
  --output-dir "$WORK"
```

By default the script selects `qc_status=PASS`.

Expected files:

```text
CORE_PHYLOGENY_V0.1/
├── genomes/
│   ├── RM10527.fasta
│   ├── ...
│   └── SKBC94.fasta
├── phylogeny_samples.tsv
└── sample_ids.txt
```

Current expected sample count:

```bash
wc -l "$WORK/sample_ids.txt"
# 52
```

Do not use PacBio methylation-call `*motifs.gff` files as bacterial gene annotations. They describe modified sites, not the coding features Panaroo requires.

---

## 3. Consistent Prokka annotation

All genomes must be annotated with the same tool/version and options. The Puma production script is:

```text
scripts/slurm/01_prokka_array.sh
```

The key command is equivalent to:

```bash
prokka \
  --outdir "prokka/${sample}" \
  --prefix "${sample}" \
  --locustag "${sample}" \
  --cpus 4 \
  "genomes/${sample}.fasta"
```

### Puma batching rule

Puma allows at most **500 tasks in one Slurm array**. `submit_core_phylogeny.sh` therefore calculates `SAMPLES_PER_TASK` automatically. For the current 52-isolate panel this is one genome per task. A larger cohort is grouped into <=500 array tasks rather than submitting an invalid oversized array.

The successful validation run used explicit `standard` partition / `cooperma` account settings. The current Prokka request is deliberately conservative; observed process memory was much lower than the requested allocation, so the request is a known-good setting rather than evidence that Prokka intrinsically requires that amount of RAM.

---

## 4. Panaroo core-genome alignment

The production script is:

```text
scripts/slurm/02_panaroo.sh
```

Equivalent command:

```bash
panaroo \
  -i prokka/*/*.gff \
  -o panaroo \
  --clean-mode strict \
  --alignment core \
  -t 16
```

The current validated run identified 1,387 gene alignments to build before writing the concatenated core alignment.

Expected key output:

```text
panaroo/core_gene_alignment.aln
```

Before tree inference, confirm the number of Prokka GFFs matches the number of intended samples.

---

## 5. IQ-TREE 2

The production script is:

```text
scripts/slurm/03_iqtree.sh
```

Equivalent command:

```bash
iqtree2 \
  -s panaroo/core_gene_alignment.aln \
  -m MFP \
  -B 1000 \
  --alrt 1000 \
  -T 16 \
  --prefix iqtree/campy52_core
```

Key output:

```text
iqtree/campy52_core.treefile
```

The completed Puma validation run finished successfully with 52 tips.

### Scientific interpretation

*Campylobacter* recombines extensively. This Panaroo + IQ-TREE tree is therefore used as a **descriptive population/lineage backbone**, not automatically interpreted as a strictly clonal transmission genealogy.

Patristic distance is still useful for asking whether closely related genomes tend to carry similar methylome repertoires. If conclusions later depend on fine-scale clonal history, add a recombination-aware sensitivity tree and repeat the lineage analysis.

---

## 6. Canonical metadata v0.2

Use:

```text
metadata/campylobacter_methylotypes_metadata_v0.2.tsv
```

This table is derived from the investigator-curated `Methylome Metadata.xlsx` workbook and replaces the placeholder-heavy v0.1 metadata for biological interpretation.

The source workbook used to create v0.2 has SHA256:

```text
a8eab63ca59c8ef3cbaf64ebbe3c89d670c1072a3a5dbdc990baf23716d44b56
```

See `metadata/README.md` for the exact host grouping, ST/CC normalisation and known missing values.

Important distinction:

```text
host_group       = ecology
provenance_class = how/where an isolate entered this dataset
```

Do not use `reference_panel`, `Salinas`, etc. as surrogate host categories.

---

## 7. Revised phylogeny + methylome figure

Run:

```bash
python "$REPO/scripts/plot_phylogeny_methylome_life_aquatic.py" \
  --tree "$WORK/iqtree/campy52_core.treefile" \
  --methylotype-dir "$COHORT/METHYLOTYPES_V0.1" \
  --metadata "$REPO/metadata/campylobacter_methylotypes_metadata_v0.2.tsv" \
  --cohort-summary "$COHORT/cohort.momento.tsv" \
  --output-dir "$COHORT/METHYLOTYPES_V0.1/PHYLOGENY_LIFE_AQUATIC_V0.2"
```

The v0.2 defaults deliberately distinguish lineage from ecology:

- isolate name is drawn at the tree tip;
- `ST` and `CC` are aligned **text** columns;
- host and provenance are categorical colour strips;
- host colours are high-contrast and colour-blind-friendly;
- RAATTY occupancy is continuous;
- accessory motif families are presence/absence;
- the accessory count bar uses the complete accessory matrix even when only the most prevalent families are drawn.

### Host palette

```text
human        #D55E00
wild_bird    #0072B2
wild_mammal  #009E73
poultry      #E69F00
environment  #CC79A7
unknown      #D9D9D9
```

STs are **not** assigned arbitrary colours. There are too many categories, and a separate colour for every ST would be harder to read than aligned text.

### Figure audit outputs

Retain these with the figure:

```text
phylogeny_tip_order.tsv
phylogeny_feature_order.tsv
phylogeny_alignment_audit.tsv
phylogeny_metadata_audit.tsv
```

`phylogeny_metadata_audit.tsv` is particularly important: it freezes the exact ST/CC/host/provenance values used in the plotted result.

---

## 8. One-command Puma workflow

For a new cohort/tree build:

```bash
cd "$REPO"
git pull

bash scripts/slurm/submit_core_phylogeny.sh "$WORK" 12
```

For the current cohort, **do not rerun Prokka/Panaroo/IQ-TREE merely to update metadata or colours**. The tree is already validated. Re-run only the plotting stage when metadata/figure code changes.

To regenerate just the v0.2 plot on Puma:

```bash
cd "$WORK"
sbatch \
  --export=ALL,WORKDIR="$WORK",MOMENTO_REPO="$REPO",COHORT_DIR="$COHORT" \
  "$REPO/scripts/slurm/04_plot_phylogeny.sh"
```

---

## 9. Validated September 2026 run

The first complete Puma dependency chain completed successfully:

```text
Panaroo  23885928  COMPLETED  00:35:29
IQ-TREE  23885929  COMPLETED  01:40:44
plot     23885930  COMPLETED  00:00:16
```

Observed batch MaxRSS was approximately:

```text
Panaroo  ~4.1 GB
IQ-TREE  ~1.6 GB
plot     ~0.65 GB
```

These values document the validation run; they are **not** universal resource requirements for future datasets.

---

## 10. Next analysis: host versus lineage

The next question is quantitative rather than visual:

> How much accessory methylome variation follows genomic lineage, and how much additional variation is associated with host ecology after lineage is represented explicitly?

Use `scripts/analyse_methylome_lineage_host.py` and the separate explainer:

```text
examples/campylobacter_lineage_host_analysis.md
```

That analysis works from the frozen accessory matrix, core-genome tree and metadata v0.2, so it can be rerun without repeating annotation or phylogeny inference.

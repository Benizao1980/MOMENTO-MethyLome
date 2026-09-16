# Campylobacter core-genome phylogeny + methylome workflow

This workflow builds a descriptive core-genome population backbone for the
QC-clean MOMENTO Campylobacter cohort and aligns methylome tracks to that tree.
It is deliberately separate from methylome clustering: the phylogeny determines
isolate order in the final figure.

## 1. Stage the QC-clean assemblies

From the cohort result directory:

```bash
python /path/to/MOMENTO-MethyLome/scripts/prepare_phylogeny_inputs.py \
  --manifest /xdisk/cooperma/bpascoe/CampyMethylome/01_ACTIVE_INPUTS/momento_campylobacter_57_manifest.tsv \
  --qc cohort.qc.tsv \
  --output-dir CORE_PHYLOGENY_V0.1
```

By default only `qc_status=PASS` samples are staged. For the current validated
cohort this should give 52 genomes and exclude SKBC102 (`PARTIAL_GFF`), RM3433
(`LOW_DATA`) and the three failed legacy inputs.

The helper writes:

```text
CORE_PHYLOGENY_V0.1/genomes/<sample_id>.fasta
CORE_PHYLOGENY_V0.1/phylogeny_samples.tsv
CORE_PHYLOGENY_V0.1/sample_ids.txt
```

Stable `sample_id.fasta` names are important because the final Newick tips must
map exactly to MOMENTO sample IDs.

## 2. Annotate all genomes consistently

Use one annotation tool/version for the full set. A Prokka example is shown
below because its GFF output can be used directly by Panaroo.

```bash
cd CORE_PHYLOGENY_V0.1
mkdir -p prokka

while read -r sample; do
  prokka \
    --outdir "prokka/${sample}" \
    --prefix "${sample}" \
    --locustag "${sample}" \
    --cpus 4 \
    "genomes/${sample}.fasta"
done < sample_ids.txt
```

For a production HPC run this loop can be converted to a Slurm array. Keep the
annotation command and tool version identical across isolates.

## 3. Build a core-gene alignment

Example with Panaroo:

```bash
mkdir -p panaroo

panaroo \
  -i prokka/*/*.gff \
  -o panaroo \
  --clean-mode strict \
  --alignment core \
  -t 16
```

The expected alignment is typically:

```text
panaroo/core_gene_alignment.aln
```

Check the exact filename written by the installed Panaroo version before using
it downstream.

## 4. Infer a descriptive core-genome tree

Example with IQ-TREE 2:

```bash
mkdir -p iqtree

# Use iqtree2 or iqtree depending on the installed executable.
iqtree2 \
  -s panaroo/core_gene_alignment.aln \
  -m MFP \
  -B 1000 \
  --alrt 1000 \
  -T AUTO \
  --prefix iqtree/campy52_core
```

The file required by the MOMENTO plotting script is:

```text
iqtree/campy52_core.treefile
```

### Interpretation caveat

Campylobacter undergoes substantial homologous recombination. This initial
Panaroo + IQ-TREE tree is useful as a descriptive lineage backbone for asking
whether methylome patterns track genomic relatedness, but it should not be
presented as a strictly clonal genealogy. If branch-length/clonal-history
interpretation becomes central, add a recombination-aware analysis and retain
both trees for sensitivity analysis.

## 5. Plot the tree with MOMENTO tracks

```bash
python /path/to/MOMENTO-MethyLome/scripts/plot_phylogeny_methylome_life_aquatic.py \
  --tree CORE_PHYLOGENY_V0.1/iqtree/campy52_core.treefile \
  --methylotype-dir METHYLOTYPES_V0.1 \
  --metadata /path/to/MOMENTO-MethyLome/metadata/campylobacter_methylotypes_metadata_v0.1.tsv \
  --cohort-summary cohort.momento.tsv \
  --output-dir METHYLOTYPES_V0.1/PHYLOGENY_LIFE_AQUATIC_V0.1
```

Defaults:

- all 52 MOMENTO PASS samples must occur in the tree;
- extra tree tips are audited and pruned;
- tree tips are ladderized for display;
- host group and provenance are shown as metadata strips;
- RAATTY occupancy is shown as a continuous track;
- accessory features present in at least two isolates are eligible for display;
- the 30 most prevalent eligible accessory features are shown by default;
- the total accessory-feature bar still counts the complete 114-feature matrix.

Useful options:

```bash
# Display more common features
--max-features 40

# Show every feature found in >=2 isolates
--max-features 0

# Restrict to features in >=3 isolates
--min-feature-prevalence 3

# Add better-curated lineage annotations later
--annotations host_group,provenance_class,st,clonal_complex

# Preserve the original Newick child order
--no-ladderize
```

Outputs include PNG/PDF/SVG figures plus:

```text
phylogeny_tip_order.tsv
phylogeny_feature_order.tsv
phylogeny_alignment_audit.tsv
```

These audit files should be retained with manuscript analysis outputs.

## 6. Next population-genetic questions

Once the tree has been checked visually, the main biological question is not
whether isolates form a few universal methylotypes (the current cohort does
not), but how much accessory-methylome similarity can be explained by genomic
lineage and how much remains associated with ecology or other metadata.

Useful follow-up analyses include:

- accessory Jaccard distance versus core-genome genetic distance;
- within-clonal-complex / within-lineage host comparisons;
- PERMANOVA or variance partitioning with lineage represented explicitly;
- motif-family gain/loss mapping on the tree;
- mapping motif families to cognate MTase / restriction-modification genes.

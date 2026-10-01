# MOMENTO core-phylogeny Slurm workflow

These scripts run the validated Puma workflow used to build the first lineage-aware Campylobacter methylome figure and the follow-on host/lineage analysis.

## Environments

The Puma setup uses separate Conda environments:

- `prokka`
- `panaroo`
- `iqtree`
- `momento`

Each batch script activates the required environment explicitly.

## Puma scheduling rules encoded here

The workflow uses the `cooperma` account and `standard` partition explicitly.

Puma permits at most **500 tasks in a single Slurm array**. The submission wrapper handles this automatically. If there are more than 500 isolates, it assigns multiple isolates sequentially to each array task so the array itself remains at or below 500 tasks.

Examples:

```text
52 genomes    -> 52 tasks x 1 genome/task
500 genomes   -> 500 tasks x 1 genome/task
2,217 genomes -> 444 tasks x 5 genomes/task
```

The optional second argument to `submit_core_phylogeny.sh` controls the maximum number of Prokka array tasks running concurrently. It must be between 1 and 500; the default is 12.

## Resource requests and what they mean

| stage | CPUs | memory request | total requested | walltime |
|---|---:|---:|---:|---:|
| Prokka array task | 4 | 4 GB/CPU | 16 GB | 2 h |
| Panaroo | 16 | 2 GB/CPU | 32 GB | 3 h |
| IQ-TREE 2 | 16 | 2 GB/CPU | 32 GB | 6 h |
| MOMENTO phylogeny figure | 2 | 2 GB/CPU | 4 GB | 30 min |
| MOMENTO host/lineage analysis | 2 | 2 GB/CPU | 4 GB | 30 min |

The Prokka allocation is a **known-good Puma configuration**, not a claim that Prokka itself needs 16 GB. In the successful 52-isolate run, Prokka process MaxRSS was generally below 1 GB per task. An earlier run with different Slurm scheduling/account context and an 8 GB total request was OOM-killed during Conda activation, so the current production configuration is kept conservative until there is a reason to optimise it.

The completed validation run used approximately 4.1 GB MaxRSS for Panaroo, 1.6 GB for IQ-TREE and 0.65 GB for plotting. Those observations are useful for future tuning but are not universal requirements.

## Full phylogeny dependency chain

From the MOMENTO repository branch containing these scripts:

```bash
bash scripts/slurm/submit_core_phylogeny.sh \
  /xdisk/cooperma/bpascoe/CampyMethylome/03_DERIVED_RESULTS/MOMENTO_RESULTS/COHORT_V0.1c/CORE_PHYLOGENY_V0.1
```

The wrapper submits:

```text
Prokka array -> Panaroo -> IQ-TREE -> MOMENTO phylogeny figure
```

using Slurm `afterok` dependencies. It writes `phylogeny_slurm_jobs.tsv` in the phylogeny work directory.

The host/lineage analysis is currently a **separate optional fifth job** rather than part of the dependency wrapper. This is deliberate: it lets us validate the new population-statistical stage on the completed 52-isolate tree before making it mandatory for every future phylogeny build.

## Expected inputs

The work directory must already contain the files written by `prepare_phylogeny_inputs.py`:

```text
sample_ids.txt
genomes/<sample_id>.fasta
phylogeny_samples.tsv
```

The plotting and population-analysis jobs also require the canonical metadata table:

```text
metadata/campylobacter_methylotypes_metadata_v0.2.tsv
```

## Key outputs

```text
prokka/<sample>/<sample>.gff
panaroo/core_gene_alignment.aln
iqtree/campy52_core.treefile
../METHYLOTYPES_V0.1/PHYLOGENY_LIFE_AQUATIC_V0.2/
../METHYLOTYPES_V0.1/LINEAGE_HOST_V0.1/
```

## Re-running only what changed

The tree does not need to be rebuilt when metadata or figure aesthetics change.

For a metadata/plot update only:

```bash
cd /xdisk/cooperma/bpascoe/CampyMethylome/03_DERIVED_RESULTS/MOMENTO_RESULTS/COHORT_V0.1c/CORE_PHYLOGENY_V0.1

sbatch \
  --export=ALL,WORKDIR="$PWD",MOMENTO_REPO=/xdisk/cooperma/bpascoe/software/MOMENTO-MethyLome,COHORT_DIR="$(dirname "$PWD")" \
  /xdisk/cooperma/bpascoe/software/MOMENTO-MethyLome/scripts/slurm/04_plot_phylogeny.sh
```

For the lineage/host analysis only:

```bash
sbatch \
  --export=ALL,WORKDIR="$PWD",MOMENTO_REPO=/xdisk/cooperma/bpascoe/software/MOMENTO-MethyLome,COHORT_DIR="$(dirname "$PWD")" \
  /xdisk/cooperma/bpascoe/software/MOMENTO-MethyLome/scripts/slurm/05_lineage_host_analysis.sh
```

## Conservative restart behaviour

- Prokka skips an isolate when its expected GFF already exists.
- Panaroo skips when `core_gene_alignment.aln` already exists.
- IQ-TREE skips when the expected `.treefile` already exists.
- figure/statistics stages write to versioned result directories rather than overwriting the original v0.1 figure.

Move/delete a stage output only when a genuinely clean rerun is intended.

## Interpretation

The IQ-TREE tree is a descriptive core-genome population backbone. *Campylobacter* recombines extensively, so it should not automatically be interpreted as a strictly clonal genealogy. The lineage/host analysis uses patristic distance as a compact representation of genomic relatedness; if conclusions become sensitive to fine-scale ancestry, repeat them with a recombination-aware tree.

For the scientific rationale and every command in context, see:

```text
examples/campylobacter_core_phylogeny_workflow.md
examples/campylobacter_lineage_host_analysis.md
metadata/README.md
```

# MOMENTO core-phylogeny Slurm workflow

These scripts run the first descriptive Campylobacter core-genome phylogeny used by the MOMENTO lineage-aware methylome workflow.

## Environments

The Puma setup used for the validated cohort has separate Conda environments named:

- `prokka`
- `panaroo`
- `iqtree`
- `momento`

Each batch script activates the required environment explicitly.

## Puma scheduling rules encoded here

The workflow uses the `cooperma` account and `standard` partition explicitly.

Puma permits at most **500 tasks in a single Slurm array**. The submission wrapper handles this automatically. If there are more than 500 isolates, it assigns multiple isolates sequentially to each array task so the array itself remains at or below 500 tasks.

Examples:

- 52 genomes -> 52 tasks x 1 genome/task
- 500 genomes -> 500 tasks x 1 genome/task
- 2,217 genomes -> 444 tasks x 5 genomes/task

The optional second argument controls the maximum number of Prokka array tasks running concurrently. It must be between 1 and 500; the default is 12.

Puma Prokka jobs use the resource pattern that has worked reliably on this cluster: 4 CPUs with `--mem-per-cpu=4G`, giving 16 GB total per task. This matters: an earlier 8 GB total request caused Conda activation itself to be OOM-killed before Prokka started.

## Resources

| stage | CPUs | memory request | total memory | walltime |
|---|---:|---:|---:|---:|
| Prokka array task | 4 | 4 GB/CPU | 16 GB | 2 h |
| Panaroo | 16 | 2 GB/CPU | 32 GB | 3 h |
| IQ-TREE 2 | 16 | 2 GB/CPU | 32 GB | 6 h |
| MOMENTO figure | 2 | 2 GB/CPU | 4 GB | 30 min |

## Submit the full dependency chain

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

## Expected inputs

The work directory must already contain the files written by `prepare_phylogeny_inputs.py`:

```text
sample_ids.txt
genomes/<sample_id>.fasta
phylogeny_samples.tsv
```

## Key outputs

```text
prokka/<sample>/<sample>.gff
panaroo/core_gene_alignment.aln
iqtree/campy52_core.treefile
../METHYLOTYPES_V0.1/PHYLOGENY_LIFE_AQUATIC_V0.1/
```

## Re-running

The stage scripts are conservative:

- Prokka skips an isolate when its expected GFF already exists.
- Panaroo skips when `core_gene_alignment.aln` already exists.
- IQ-TREE skips when the expected `.treefile` already exists.

Delete or move a stage output if a clean rerun is required.

## Interpretation

This tree is a descriptive core-genome population backbone. Because Campylobacter recombines extensively, do not treat it automatically as a strictly clonal genealogy. If branch-length or ancestry interpretation becomes central, add a recombination-aware sensitivity analysis.

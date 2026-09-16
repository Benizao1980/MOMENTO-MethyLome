# MOMENTO core-phylogeny Slurm workflow

These scripts run the first descriptive Campylobacter core-genome phylogeny used by the MOMENTO lineage-aware methylome workflow.

## Environments

The Puma setup used for the validated cohort has separate Conda environments named:

- `prokka`
- `panaroo`
- `iqtree`
- `momento`

Each batch script activates the required environment explicitly.

## Resources

| stage | CPUs | memory | walltime |
|---|---:|---:|---:|
| Prokka, per isolate | 4 | 8 GB | 1 h |
| Panaroo | 16 | 32 GB | 3 h |
| IQ-TREE 2 | 16 | 32 GB | 6 h |
| MOMENTO figure | 2 | 4 GB | 30 min |

The submission wrapper defaults to at most 12 simultaneous Prokka tasks. Pass a second argument to change that cap.

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

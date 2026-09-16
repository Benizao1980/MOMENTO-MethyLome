#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 /path/to/CORE_PHYLOGENY_V0.1 [max_concurrent_prokka]" >&2
  exit 2
fi

WORKDIR="$(realpath "$1")"
MAX_CONCURRENT="${2:-12}"
MAX_ARRAY_TASKS=500
MOMENTO_REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
COHORT_DIR="$(dirname "$WORKDIR")"
SLURM_DIR="$MOMENTO_REPO/scripts/slurm"

cd "$WORKDIR"
mkdir -p slurm_logs

if [[ ! -s sample_ids.txt ]]; then
  echo "ERROR: missing $WORKDIR/sample_ids.txt" >&2
  exit 2
fi

N=$(wc -l < sample_ids.txt)
if [[ "$N" -lt 2 ]]; then
  echo "ERROR: expected multiple samples, found $N" >&2
  exit 2
fi

# Puma permits at most 500 tasks in a single Slurm array. Keep this explicit so
# larger cohorts are chunked deliberately rather than failing at submission.
if (( N > MAX_ARRAY_TASKS )); then
  echo "ERROR: Puma array limit is ${MAX_ARRAY_TASKS} tasks; found ${N} samples." >&2
  echo "Split/chunk the Prokka array before submission." >&2
  exit 2
fi

if ! [[ "$MAX_CONCURRENT" =~ ^[0-9]+$ ]] || (( MAX_CONCURRENT < 1 || MAX_CONCURRENT > MAX_ARRAY_TASKS )); then
  echo "ERROR: max_concurrent_prokka must be an integer between 1 and ${MAX_ARRAY_TASKS}" >&2
  exit 2
fi

echo "Submitting MOMENTO core-phylogeny workflow"
echo "  workdir: $WORKDIR"
echo "  samples: $N"
echo "  Puma max array tasks: $MAX_ARRAY_TASKS"
echo "  max concurrent Prokka tasks: $MAX_CONCURRENT"

PROKKA_RAW=$(sbatch --parsable \
  --chdir="$WORKDIR" \
  --array="1-${N}%${MAX_CONCURRENT}" \
  --export=ALL,WORKDIR="$WORKDIR" \
  "$SLURM_DIR/01_prokka_array.sh")
PROKKA_JOB="${PROKKA_RAW%%;*}"

PANAROO_RAW=$(sbatch --parsable \
  --chdir="$WORKDIR" \
  --dependency="afterok:${PROKKA_JOB}" \
  --export=ALL,WORKDIR="$WORKDIR" \
  "$SLURM_DIR/02_panaroo.sh")
PANAROO_JOB="${PANAROO_RAW%%;*}"

IQTREE_RAW=$(sbatch --parsable \
  --chdir="$WORKDIR" \
  --dependency="afterok:${PANAROO_JOB}" \
  --export=ALL,WORKDIR="$WORKDIR" \
  "$SLURM_DIR/03_iqtree.sh")
IQTREE_JOB="${IQTREE_RAW%%;*}"

PLOT_RAW=$(sbatch --parsable \
  --chdir="$WORKDIR" \
  --dependency="afterok:${IQTREE_JOB}" \
  --export=ALL,WORKDIR="$WORKDIR",MOMENTO_REPO="$MOMENTO_REPO",COHORT_DIR="$COHORT_DIR" \
  "$SLURM_DIR/04_plot_phylogeny.sh")
PLOT_JOB="${PLOT_RAW%%;*}"

cat > phylogeny_slurm_jobs.tsv <<EOF
stage\tjob_id\tdependency
prokka\t${PROKKA_JOB}\t-
panaroo\t${PANAROO_JOB}\tafterok:${PROKKA_JOB}
iqtree\t${IQTREE_JOB}\tafterok:${PANAROO_JOB}
plot\t${PLOT_JOB}\tafterok:${IQTREE_JOB}
EOF

echo
echo "Submitted:"
column -t -s $'\t' phylogeny_slurm_jobs.tsv || cat phylogeny_slurm_jobs.tsv

echo
echo "Monitor with:"
echo "  squeue -j ${PROKKA_JOB},${PANAROO_JOB},${IQTREE_JOB},${PLOT_JOB}"
echo "  tail -f slurm_logs/prokka_${PROKKA_JOB}_*.out"
echo
echo "Job table: $WORKDIR/phylogeny_slurm_jobs.tsv"

#!/usr/bin/env bash
#SBATCH --job-name=momento_hostlineage
#SBATCH --account=cooperma
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:30:00
#SBATCH --output=slurm_logs/hostlineage_%j.out
#SBATCH --error=slurm_logs/hostlineage_%j.err

set -euo pipefail
WORKDIR="${WORKDIR:-$PWD}"
MOMENTO_REPO="${MOMENTO_REPO:-/xdisk/cooperma/bpascoe/software/MOMENTO-MethyLome}"
COHORT_DIR="${COHORT_DIR:-$(dirname "$WORKDIR")}"
cd "$WORKDIR"
mkdir -p slurm_logs

TREE="iqtree/campy52_core.treefile"
METADATA="$MOMENTO_REPO/metadata/campylobacter_methylotypes_metadata_v0.2.tsv"
OUTDIR="$COHORT_DIR/METHYLOTYPES_V0.1/LINEAGE_HOST_V0.1"

if [[ ! -s "$TREE" ]]; then
  echo "ERROR: missing tree: $TREE" >&2
  exit 2
fi
if [[ ! -s "$METADATA" ]]; then
  echo "ERROR: missing metadata: $METADATA" >&2
  exit 2
fi

source /home/u12/bpascoe/miniconda3/etc/profile.d/conda.sh
conda activate momento

python "$MOMENTO_REPO/scripts/analyse_methylome_lineage_host.py" \
  --tree "$TREE" \
  --methylotype-dir "$COHORT_DIR/METHYLOTYPES_V0.1" \
  --metadata "$METADATA" \
  --output-dir "$OUTDIR" \
  --permutations 999 \
  --seed 20260917

echo "LINEAGE/HOST ANALYSIS COMPLETE: $OUTDIR"

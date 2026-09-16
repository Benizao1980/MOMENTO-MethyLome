#!/usr/bin/env bash
#SBATCH --job-name=momento_treeplot
#SBATCH --account=cooperma
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:30:00
#SBATCH --output=slurm_logs/treeplot_%j.out
#SBATCH --error=slurm_logs/treeplot_%j.err

set -euo pipefail
WORKDIR="${WORKDIR:-$PWD}"
MOMENTO_REPO="${MOMENTO_REPO:-/xdisk/cooperma/bpascoe/software/MOMENTO-MethyLome}"
COHORT_DIR="${COHORT_DIR:-$(dirname "$WORKDIR")}"
cd "$WORKDIR"
mkdir -p slurm_logs

TREE="iqtree/campy52_core.treefile"
if [[ ! -s "$TREE" ]]; then
  echo "ERROR: missing tree: $TREE" >&2
  exit 2
fi

source /home/u12/bpascoe/miniconda3/etc/profile.d/conda.sh
conda activate momento

python "$MOMENTO_REPO/scripts/plot_phylogeny_methylome_life_aquatic.py" \
  --tree "$TREE" \
  --methylotype-dir "$COHORT_DIR/METHYLOTYPES_V0.1" \
  --metadata "$MOMENTO_REPO/metadata/campylobacter_methylotypes_metadata_v0.1.tsv" \
  --cohort-summary "$COHORT_DIR/cohort.momento.tsv" \
  --output-dir "$COHORT_DIR/METHYLOTYPES_V0.1/PHYLOGENY_LIFE_AQUATIC_V0.1"

echo "PHYLOGENY FIGURE COMPLETE"

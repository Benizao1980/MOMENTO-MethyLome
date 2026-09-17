#!/usr/bin/env bash
#SBATCH --job-name=momento_iqtree
#SBATCH --account=cooperma
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
#SBATCH --time=06:00:00
#SBATCH --output=slurm_logs/iqtree_%j.out
#SBATCH --error=slurm_logs/iqtree_%j.err

set -euo pipefail
WORKDIR="${WORKDIR:-$PWD}"
cd "$WORKDIR"
mkdir -p iqtree slurm_logs

ALIGNMENT="panaroo/core_gene_alignment.aln"
if [[ ! -s "$ALIGNMENT" ]]; then
  echo "ERROR: missing alignment: $ALIGNMENT" >&2
  exit 2
fi

source /home/u12/bpascoe/miniconda3/etc/profile.d/conda.sh
conda activate iqtree

if command -v iqtree2 >/dev/null 2>&1; then
  IQTREE=iqtree2
elif command -v iqtree >/dev/null 2>&1; then
  IQTREE=iqtree
else
  echo "ERROR: iqtree2/iqtree not found in iqtree environment" >&2
  exit 2
fi

echo "IQ-TREE COMMAND: $IQTREE"
"$IQTREE" --version || true

PREFIX="iqtree/campy52_core"
if [[ -s "${PREFIX}.treefile" ]]; then
  echo "IQ-TREE: existing treefile found; skipping"
  exit 0
fi

"$IQTREE" \
  -s "$ALIGNMENT" \
  -m MFP \
  -B 1000 \
  --alrt 1000 \
  -T "${SLURM_CPUS_PER_TASK:-16}" \
  --prefix "$PREFIX"

if [[ ! -s "${PREFIX}.treefile" ]]; then
  echo "ERROR: IQ-TREE completed without ${PREFIX}.treefile" >&2
  exit 3
fi

echo "IQ-TREE COMPLETE: ${PREFIX}.treefile"

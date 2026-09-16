#!/usr/bin/env bash
#SBATCH --job-name=momento_panaroo
#SBATCH --account=cooperma
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
#SBATCH --time=03:00:00
#SBATCH --output=slurm_logs/panaroo_%j.out
#SBATCH --error=slurm_logs/panaroo_%j.err

set -euo pipefail
WORKDIR="${WORKDIR:-$PWD}"
cd "$WORKDIR"
mkdir -p panaroo slurm_logs

source /home/u12/bpascoe/miniconda3/etc/profile.d/conda.sh
conda activate panaroo

echo "PANAROO VERSION: $(panaroo --version 2>&1 | tr '\n' ' ')"

mapfile -t GFFS < <(find prokka -mindepth 2 -maxdepth 2 -type f -name '*.gff' | sort)
EXPECTED=$(wc -l < sample_ids.txt)
FOUND=${#GFFS[@]}
echo "GFF COUNT: ${FOUND}/${EXPECTED}"
if [[ "$FOUND" -ne "$EXPECTED" ]]; then
  echo "ERROR: expected ${EXPECTED} Prokka GFFs but found ${FOUND}" >&2
  exit 2
fi

if [[ -s panaroo/core_gene_alignment.aln ]]; then
  echo "PANAROO: existing core_gene_alignment.aln found; skipping"
  exit 0
fi

panaroo \
  -i "${GFFS[@]}" \
  -o panaroo \
  --clean-mode strict \
  --alignment core \
  -t "${SLURM_CPUS_PER_TASK:-16}"

if [[ ! -s panaroo/core_gene_alignment.aln ]]; then
  echo "ERROR: Panaroo completed but panaroo/core_gene_alignment.aln is missing" >&2
  echo "Contents of panaroo/:"
  ls -lah panaroo
  exit 3
fi

echo "PANAROO COMPLETE: $(du -h panaroo/core_gene_alignment.aln | cut -f1) core alignment"

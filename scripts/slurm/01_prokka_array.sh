#!/usr/bin/env bash
#SBATCH --job-name=momento_prokka
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=slurm_logs/prokka_%A_%a.out
#SBATCH --error=slurm_logs/prokka_%A_%a.err

set -euo pipefail

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID is not set. Submit as an array job." >&2
  exit 2
fi

WORKDIR="${WORKDIR:-$PWD}"
cd "$WORKDIR"
mkdir -p prokka slurm_logs

if [[ ! -s sample_ids.txt ]]; then
  echo "ERROR: missing sample_ids.txt in $WORKDIR" >&2
  exit 2
fi

SAMPLE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" sample_ids.txt)"
if [[ -z "$SAMPLE" ]]; then
  echo "ERROR: no sample for array index ${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

FASTA="genomes/${SAMPLE}.fasta"
OUTDIR="prokka/${SAMPLE}"
if [[ ! -e "$FASTA" ]]; then
  echo "ERROR: missing FASTA: $FASTA" >&2
  exit 2
fi

source /home/u12/bpascoe/miniconda3/etc/profile.d/conda.sh
conda activate prokka

printf 'sample\t%s\n' "$SAMPLE"
printf 'host\t%s\n' "$(hostname)"
printf 'prokka\t%s\n' "$(prokka --version 2>&1 | tr '\n' ' ')"

if [[ -s "${OUTDIR}/${SAMPLE}.gff" ]]; then
  echo "PROKKA: existing GFF found for ${SAMPLE}; skipping"
  exit 0
fi

prokka \
  --outdir "$OUTDIR" \
  --prefix "$SAMPLE" \
  --locustag "$SAMPLE" \
  --cpus "${SLURM_CPUS_PER_TASK:-4}" \
  "$FASTA"

if [[ ! -s "${OUTDIR}/${SAMPLE}.gff" ]]; then
  echo "ERROR: Prokka completed without expected GFF for ${SAMPLE}" >&2
  exit 3
fi

echo "PROKKA COMPLETE: ${SAMPLE}"

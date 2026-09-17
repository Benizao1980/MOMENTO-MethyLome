#!/usr/bin/env bash
#SBATCH --job-name=momento_prokka
#SBATCH --account=cooperma
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --output=slurm_logs/prokka_%A_%a.out
#SBATCH --error=slurm_logs/prokka_%A_%a.err

set -euo pipefail

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID is not set. Submit as an array job." >&2
  exit 2
fi

WORKDIR="${WORKDIR:-$PWD}"
SAMPLES_PER_TASK="${SAMPLES_PER_TASK:-1}"
cd "$WORKDIR"
mkdir -p prokka slurm_logs

if [[ ! -s sample_ids.txt ]]; then
  echo "ERROR: missing sample_ids.txt in $WORKDIR" >&2
  exit 2
fi

if ! [[ "$SAMPLES_PER_TASK" =~ ^[0-9]+$ ]] || (( SAMPLES_PER_TASK < 1 )); then
  echo "ERROR: SAMPLES_PER_TASK must be a positive integer" >&2
  exit 2
fi

TOTAL_SAMPLES=$(wc -l < sample_ids.txt)
START_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) * SAMPLES_PER_TASK + 1 ))
END_INDEX=$(( SLURM_ARRAY_TASK_ID * SAMPLES_PER_TASK ))
if (( END_INDEX > TOTAL_SAMPLES )); then
  END_INDEX=$TOTAL_SAMPLES
fi
if (( START_INDEX > TOTAL_SAMPLES )); then
  echo "ERROR: array task ${SLURM_ARRAY_TASK_ID} maps beyond ${TOTAL_SAMPLES} samples" >&2
  exit 2
fi

source /home/u12/bpascoe/miniconda3/etc/profile.d/conda.sh
conda activate prokka

printf 'host\t%s\n' "$(hostname)"
printf 'prokka\t%s\n' "$(prokka --version 2>&1 | tr '\n' ' ')"
printf 'array_task\t%s\n' "$SLURM_ARRAY_TASK_ID"
printf 'sample_index_range\t%s-%s\n' "$START_INDEX" "$END_INDEX"
printf 'samples_per_task\t%s\n' "$SAMPLES_PER_TASK"

for (( SAMPLE_INDEX=START_INDEX; SAMPLE_INDEX<=END_INDEX; SAMPLE_INDEX++ )); do
  SAMPLE="$(sed -n "${SAMPLE_INDEX}p" sample_ids.txt)"
  if [[ -z "$SAMPLE" ]]; then
    echo "ERROR: no sample for sample index ${SAMPLE_INDEX}" >&2
    exit 2
  fi

  FASTA="genomes/${SAMPLE}.fasta"
  OUTDIR="prokka/${SAMPLE}"
  if [[ ! -e "$FASTA" ]]; then
    echo "ERROR: missing FASTA: $FASTA" >&2
    exit 2
  fi

  printf '\nPROKKA SAMPLE\t%s\tindex=%s\n' "$SAMPLE" "$SAMPLE_INDEX"

  if [[ -s "${OUTDIR}/${SAMPLE}.gff" ]]; then
    echo "PROKKA: existing GFF found for ${SAMPLE}; skipping"
    continue
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
done

echo "PROKKA ARRAY TASK COMPLETE: ${SLURM_ARRAY_TASK_ID}"

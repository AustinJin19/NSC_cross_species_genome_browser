#!/usr/bin/env bash
#SBATCH --job-name=atac-rose
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=24:00:00
set -euo pipefail
pipeline_dir=$1
manifest=$2
project_root=$3
output=$4
shift 4
python "$pipeline_dir/workflow.py" run-sample \
  --manifest "$manifest" --project-root "$project_root" --output "$output" \
  --index "${SLURM_ARRAY_TASK_ID:?Submit as a SLURM array}" "$@"

#!/usr/bin/env bash
set -euo pipefail
if (( $# < 3 )); then
  echo "Usage: $0 samples.tsv PROJECT_ROOT OUTPUT_DIR [run-sample options]" >&2
  exit 2
fi
pipeline_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
manifest=$(realpath "$1")
project_root=$(realpath "$2")
mkdir -p -- "$3/logs"
output=$(realpath "$3")
shift 3
count=$(python - "$manifest" <<'PY'
import csv, sys
with open(sys.argv[1]) as f:
    rows = list(csv.DictReader(f, delimiter='\t'))
if not rows or any(not r.get('bam') for r in rows):
    raise SystemExit('Fill every BAM path, or explicitly subset the manifest before submission.')
print(len(rows))
PY
)
sbatch --array="0-$((count-1))%${MAX_PARALLEL_JOBS:-4}" \
  --output="$output/logs/%A_%a.out" --error="$output/logs/%A_%a.err" \
  "$pipeline_dir/run_array.sbatch" "$pipeline_dir" "$manifest" "$project_root" "$output" "$@"

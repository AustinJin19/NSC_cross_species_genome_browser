#!/usr/bin/env bash
set -euo pipefail
if (( $# < 3 )); then
  echo "Usage: $0 samples.tsv PROJECT_ROOT OUTPUT_DIR [run-sample options]" >&2
  exit 2
fi
pipeline_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
manifest=$1
project_root=$2
output=$3
shift 3
count=$(python - "$manifest" <<'PY'
import csv, sys
with open(sys.argv[1]) as f:
    rows = list(csv.DictReader(f, delimiter='\t'))
if not rows or any(not r.get('bam') for r in rows):
    raise SystemExit('Fill every BAM path or subset the manifest first.')
print(len(rows))
PY
)
for (( index=0; index<count; index++ )); do
  python "$pipeline_dir/workflow.py" run-sample \
    --manifest "$manifest" --project-root "$project_root" --output "$output" \
    --index "$index" "$@"
done
python "$pipeline_dir/workflow.py" summarize --manifest "$manifest" --output "$output"

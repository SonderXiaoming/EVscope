#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$repo_dir"

mkdir -p tests/smoke/work

bash -n EVscope.sh
for sh in bin/*.sh; do bash -n "$sh"; done
python3 -m py_compile bin/*.py

python3 bin/Step_02_calculate_ACC_motif_fraction.py \
  --input_fastq tests/smoke/data/toy_R2.fastq.gz \
  --positions 1-3 \
  --output_tsv tests/smoke/work/acc.tsv

test -s tests/smoke/work/acc.tsv

if python3 - <<'PY'
import importlib.util, sys
sys.exit(0 if importlib.util.find_spec('matplotlib') else 1)
PY
then
  python3 bin/Step_03_plot_fastq_read_length_dist.py \
    --input_fastqs tests/smoke/data/toy_R2.fastq.gz tests/smoke/data/toy_R2.fastq.gz tests/smoke/data/toy_R2.fastq.gz \
    --titles "toy raw" "toy trimmed" "toy clean" \
    --output_pdf tests/smoke/work/read_length.pdf \
    --output_png tests/smoke/work/read_length.png
  test -s tests/smoke/work/read_length.png
  test -s tests/smoke/work/read_length.pdf
else
  echo "WARN: matplotlib unavailable; skipped read length plotting smoke" >&2
fi

while read -r rel; do
  [[ -z "$rel" ]] && continue
  if [[ "$rel" == work/read_length.png && ! -e "tests/smoke/$rel" ]]; then
    continue
  fi
  test -s "tests/smoke/$rel"
done < tests/smoke/expected/expected_files.txt

sha256sum tests/smoke/data/toy_R2.fastq.gz tests/smoke/work/acc.tsv > tests/smoke/work/sha256sums.txt

echo "EVscope smoke tests passed"

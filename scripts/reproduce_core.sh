#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

python scripts/analysis/recalculate_v7_key_results.py --root baseline/v7 --outdir qa/recomputed/v7_python
Rscript scripts/analysis/recalculate_v7_key_results.R baseline/v7 qa/recomputed/v7_r
Rscript scripts/analysis/05_independent_core_statistics_check.R . qa/recomputed/v8_core
python scripts/analysis/04_gse178341_matching_balance_qc.py --analysis-dir analysis/matched_null_gse178341 --outdir qa/recomputed/matched_balance
Rscript scripts/analysis/recalculate_v8_1_corrections.R . qa/recomputed/v8_1_corrections

echo "Core reproducibility checks completed. See qa/recomputed/."


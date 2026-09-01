#!/usr/bin/env python
"""Transparent downstream summaries from the GSE178341 patient-score table."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from scipy import stats


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--patient-scores", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    data = pd.read_csv(args.patient_scores, sep="\t", index_col=0)
    pairs = [
        ("module_all", "T_cell_fraction_all"),
        ("module_all", "stromal_fraction_all"),
        ("module_all", "module_TNKILC"),
        ("module_Epi", "T_cell_fraction_all"),
        ("module_Epi", "stromal_fraction_all"),
        ("module_Epi", "module_TNKILC"),
    ]
    rows = []
    for first, second in pairs:
        frame = data[[first, second]].dropna()
        test = stats.spearmanr(frame[first], frame[second])
        rows.append({"variable_1": first, "variable_2": second, "n": len(frame), "spearman_rho": test.statistic, "p_two_sided": test.pvalue})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
    print(pd.DataFrame(rows).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

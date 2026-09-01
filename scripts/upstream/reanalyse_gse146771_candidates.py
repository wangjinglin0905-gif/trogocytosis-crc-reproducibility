from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd


def wilson(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n == 0:
        return np.nan, np.nan
    p = k / n
    denominator = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denominator
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denominator
    return center - half, center + half


def exact_signflip(differences: np.ndarray) -> float:
    differences = np.asarray(differences, float)
    differences = differences[np.isfinite(differences)]
    if len(differences) == 0:
        return np.nan
    observed = abs(differences.mean())
    if len(differences) <= 20:
        null = [abs(np.mean(differences * np.asarray(signs))) for signs in itertools.product([-1, 1], repeat=len(differences))]
        return float(np.mean(np.asarray(null) >= observed - 1e-15))
    rng = np.random.default_rng(20260831)
    null = np.abs((rng.choice([-1, 1], size=(100000, len(differences))) * differences).mean(axis=1))
    return float((np.sum(null >= observed) + 1) / (len(null) + 1))


def patient_table(cell: pd.DataFrame, candidate: pd.Series, lineage: str | None = None) -> pd.DataFrame:
    d = cell.copy()
    d["candidate"] = candidate.astype(bool)
    if lineage is not None:
        d = d[d["Global_Cluster"].eq(lineage)]
    out = d.groupby(["Sample", "Tissue"]).agg(n=("candidate", "size"), candidates=("candidate", "sum")).reset_index()
    out["rate"] = out["candidates"] / out["n"]
    out["lineage"] = lineage or "All leukocytes"
    return out


def paired_test(summary: pd.DataFrame) -> dict:
    wide = summary[summary["Tissue"].isin(["T", "N"])].pivot(index="Sample", columns="Tissue", values="rate").dropna()
    differences = (wide["T"] - wide["N"]).to_numpy()
    return {
        "paired_patients": len(wide),
        "mean_rate_difference_T_minus_N": float(np.mean(differences)) if len(differences) else np.nan,
        "median_rate_difference_T_minus_N": float(np.median(differences)) if len(differences) else np.nan,
        "exact_signflip_p": exact_signflip(differences),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cell-scores", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    cell = pd.read_csv(args.cell_scores, index_col=0)
    for column in ["epi_score", "CEACAM5_TPM"]:
        cell[column] = pd.to_numeric(cell[column], errors="coerce")
    negative = cell[cell["Tissue"].eq("P")]
    primary_epi_threshold = float(np.nanquantile(negative["epi_score"], 0.995))
    primary = cell["CEACAM5_TPM"].gt(1.0) & cell["epi_score"].le(primary_epi_threshold)
    cell["CEACAM5_candidate_primary"] = primary

    lineages = ["All leukocytes"] + sorted(cell["Global_Cluster"].dropna().unique().tolist())
    summaries = []
    tests = []
    for lineage in lineages:
        summary = patient_table(cell, primary, None if lineage == "All leukocytes" else lineage)
        summaries.append(summary)
        tests.append({"lineage": lineage, **paired_test(summary)})
    sample_summary = pd.concat(summaries, ignore_index=True)
    sample_summary.to_csv(args.outdir / "gse146771_candidate_rates_by_patient_tissue_lineage.csv", index=False)
    paired = pd.DataFrame(tests)
    paired.to_csv(args.outdir / "gse146771_paired_patient_tests.csv", index=False)

    tissue_rows = []
    for tissue, d in cell.groupby("Tissue"):
        k = int(d["CEACAM5_candidate_primary"].sum())
        n = len(d)
        lo, hi = wilson(k, n)
        tissue_rows.append({"tissue": tissue, "n_cells": n, "candidates": k, "rate": k / n, "wilson95_lo": lo, "wilson95_hi": hi})
    pd.DataFrame(tissue_rows).to_csv(args.outdir / "gse146771_candidate_tissue_descriptive.csv", index=False)

    grid_rows = []
    for cea_threshold in [0.0, 1.0, 2.0]:
        for epi_quantile in [0.99, 0.995, 0.999]:
            epi_threshold = float(np.nanquantile(negative["epi_score"], epi_quantile))
            candidate = cell["CEACAM5_TPM"].gt(cea_threshold) & cell["epi_score"].le(epi_threshold)
            summary = patient_table(cell, candidate)
            grid_rows.append({
                "CEACAM5_TPM_threshold": cea_threshold,
                "negative_control_epi_quantile": epi_quantile,
                "epi_threshold": epi_threshold,
                "total_candidates": int(candidate.sum()),
                "candidate_patients": int(cell.loc[candidate, "Sample"].nunique()),
                **paired_test(summary),
            })
    grid = pd.DataFrame(grid_rows)
    grid.to_csv(args.outdir / "gse146771_threshold_sensitivity.csv", index=False)

    primary_all = patient_table(cell, primary)
    paired_patients = primary_all[primary_all["Tissue"].isin(["T", "N"])].pivot(index="Sample", columns="Tissue", values="rate").dropna().index
    loo_rows = []
    for omitted in paired_patients:
        subset = primary_all[primary_all["Sample"].ne(omitted)]
        loo_rows.append({"omitted_patient": omitted, **paired_test(subset)})
    pd.DataFrame(loo_rows).to_csv(args.outdir / "gse146771_leave_one_patient_out.csv", index=False)

    candidate_cells = cell.loc[primary, ["Sample", "Tissue", "Global_Cluster", "Sub_Cluster", "CEACAM5_TPM", "epi_score"]]
    candidate_cells.to_csv(args.outdir / "gse146771_primary_candidate_cells.csv")
    audit = {
        "analysis_unit": "patient for tumor-versus-normal inference; cells are measurement units only",
        "primary_label": "author-annotated leukocyte with CEACAM5 TPM>1 and epithelial-program score at or below the peripheral-blood 99.5th percentile",
        "primary_epi_threshold": primary_epi_threshold,
        "n_cells": len(cell),
        "n_primary_candidates": int(primary.sum()),
        "n_candidate_patients": int(cell.loc[primary, "Sample"].nunique()),
        "terminology": "contamination-filtered CEACAM5-positive leukocyte candidate; not a directly observed trogocytosis event",
        "prohibited_inference": "cell-level Fisher tests are descriptive only and are not used for cohort-level significance",
    }
    (args.outdir / "gse146771_candidate_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(paired.to_string(index=False))
    print(grid.to_string(index=False))
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

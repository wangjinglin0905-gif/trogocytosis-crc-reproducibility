#!/usr/bin/env python
"""Prespecified DepMap 26Q1 replication of the SCD-VPS72 dependency relation."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


SEED = 20260901
TARGETS = ("SCD", "VPS72")


def sha256(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def symbol(column: str) -> str:
    return column.split(" (")[0]


def residual_correlation(scd: np.ndarray, vps72: np.ndarray, covariates: np.ndarray) -> float:
    x_rank = stats.rankdata(scd)
    y_rank = stats.rankdata(vps72)
    cov = np.asarray(covariates, dtype=float)
    cov_mean = np.nanmean(cov, axis=0)
    cov_sd = np.nanstd(cov, axis=0, ddof=1)
    keep = np.isfinite(cov_sd) & (cov_sd > 1e-10)
    cov = (cov[:, keep] - cov_mean[keep]) / cov_sd[keep]
    design = np.column_stack([np.ones(len(scd)), cov])
    x_res = x_rank - design @ np.linalg.lstsq(design, x_rank, rcond=None)[0]
    y_res = y_rank - design @ np.linalg.lstsq(design, y_rank, rcond=None)[0]
    return float(stats.pearsonr(x_res, y_res).statistic)


def adjusted_residuals(scd: np.ndarray, vps72: np.ndarray, covariates: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x_rank = stats.rankdata(scd)
    y_rank = stats.rankdata(vps72)
    cov = np.asarray(covariates, dtype=float)
    cov_mean = np.nanmean(cov, axis=0)
    cov_sd = np.nanstd(cov, axis=0, ddof=1)
    keep = np.isfinite(cov_sd) & (cov_sd > 1e-10)
    cov = (cov[:, keep] - cov_mean[keep]) / cov_sd[keep]
    design = np.column_stack([np.ones(len(scd)), cov])
    return (
        x_rank - design @ np.linalg.lstsq(design, x_rank, rcond=None)[0],
        y_rank - design @ np.linalg.lstsq(design, y_rank, rcond=None)[0],
    )


def bootstrap(
    scd: np.ndarray,
    vps72: np.ndarray,
    covariates: np.ndarray,
    reps: int,
    rng: np.random.Generator,
) -> pd.DataFrame:
    rows = []
    n = len(scd)
    for rep in range(reps):
        idx = rng.integers(0, n, n)
        raw = stats.spearmanr(scd[idx], vps72[idx]).statistic
        adjusted = residual_correlation(scd[idx], vps72[idx], covariates[idx])
        rows.append((rep + 1, raw, adjusted))
    return pd.DataFrame(rows, columns=["replicate", "raw_rho", "adjusted_partial_spearman"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-effect", required=True, type=Path)
    parser.add_argument("--model", required=True, type=Path)
    parser.add_argument("--organoid-qc", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--reps", type=int, default=10000)
    parser.add_argument("--chunk-rows", type=int, default=128)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    model = pd.read_csv(args.model, low_memory=False)
    model_id_column = "ModelID" if "ModelID" in model.columns else "DepMap_ID"
    two_d = model.loc[model["ModelType"].eq("Cell Line")].copy()
    two_d_ids = set(two_d[model_id_column].astype(str))
    crc_ids = set(
        two_d.loc[two_d["OncotreePrimaryDisease"].eq("Colorectal Adenocarcinoma"), model_id_column].astype(str)
    )

    header = pd.read_csv(args.gene_effect, nrows=0).columns.tolist()
    raw_columns = header[1:]
    symbols = np.asarray([symbol(x) for x in raw_columns], dtype=str)
    unique_symbols, symbol_multiplicity = np.unique(symbols, return_counts=True)
    unique_set = set(unique_symbols[symbol_multiplicity == 1])
    keep_columns = [col for col, sym in zip(raw_columns, symbols) if sym in unique_set]
    kept_symbols = np.asarray([symbol(x) for x in keep_columns], dtype=str)
    if any(target not in set(kept_symbols) for target in TARGETS):
        raise RuntimeError("SCD or VPS72 absent/duplicated in DepMap matrix")

    ids: list[str] = []
    matrices: list[np.ndarray] = []
    for chunk in pd.read_csv(args.gene_effect, index_col=0, usecols=[header[0], *keep_columns], chunksize=args.chunk_rows):
        chunk.index = chunk.index.astype(str)
        selected = chunk.loc[chunk.index.isin(two_d_ids)]
        if len(selected):
            ids.extend(selected.index.tolist())
            matrices.append(selected.to_numpy(dtype=np.float32))
    matrix = np.vstack(matrices).astype(np.float64, copy=False)
    ids_array = np.asarray(ids, dtype=str)
    gene_index = {gene: i for i, gene in enumerate(kept_symbols)}

    completeness = np.mean(np.isfinite(matrix), axis=0)
    medians = np.nanmedian(matrix, axis=0)
    prop_dependent = np.nanmean(matrix < -0.5, axis=0)
    common_mask = (completeness >= 0.90) & (medians <= -0.5) & (prop_dependent >= 0.80)
    common_mask[[gene_index[x] for x in TARGETS]] = False
    common_genes = kept_symbols[common_mask]
    common = matrix[:, common_mask].copy()
    common_col_median = np.nanmedian(common, axis=0)
    missing_rows, missing_cols = np.where(~np.isfinite(common))
    common[missing_rows, missing_cols] = common_col_median[missing_cols]
    common_mean = common.mean(axis=0)
    common_sd = common.std(axis=0, ddof=1)
    variable = common_sd > 1e-8
    common = (common[:, variable] - common_mean[variable]) / common_sd[variable]
    u, singular_values, _ = np.linalg.svd(common, full_matrices=False)
    pcs = u[:, :5] * singular_values[:5]

    row_median = np.nanmedian(matrix, axis=1)
    row_mad = np.nanmedian(np.abs(matrix - row_median[:, None]), axis=1)
    covariates_all = np.column_stack([row_median, row_mad, pcs])
    scd_all = matrix[:, gene_index["SCD"]]
    vps_all = matrix[:, gene_index["VPS72"]]
    crc_mask = np.isin(ids_array, list(crc_ids))
    complete = crc_mask & np.isfinite(scd_all) & np.isfinite(vps_all) & np.all(np.isfinite(covariates_all), axis=1)
    scd, vps72, covariates = scd_all[complete], vps_all[complete], covariates_all[complete]
    crc_complete_ids = ids_array[complete]
    if len(scd) < 30:
        raise RuntimeError(f"Too few complete CRC lines: {len(scd)}")

    raw_test = stats.spearmanr(scd, vps72)
    adjusted_rho = residual_correlation(scd, vps72, covariates)
    x_res, y_res = adjusted_residuals(scd, vps72, covariates)
    permuted = np.empty(args.reps, dtype=float)
    for i in range(args.reps):
        permuted[i] = stats.pearsonr(x_res, rng.permutation(y_res)).statistic
    adjusted_p = (1 + np.sum(np.abs(permuted) >= abs(adjusted_rho))) / (args.reps + 1)
    boot = bootstrap(scd, vps72, covariates, args.reps, rng)
    boot.to_csv(args.outdir / "depmap_scd_vps72_bootstrap.tsv.gz", sep="\t", index=False)
    pd.DataFrame({"replicate": np.arange(1, args.reps + 1), "permuted_rho": permuted}).to_csv(
        args.outdir / "depmap_scd_vps72_permutation.tsv.gz", sep="\t", index=False
    )

    organoid = pd.read_csv(args.organoid_qc)
    organoid_row = organoid.loc[organoid["gene"].eq("VPS72")].iloc[0]
    raw_ci = np.quantile(boot["raw_rho"].dropna(), [0.025, 0.975])
    adjusted_ci = np.quantile(boot["adjusted_partial_spearman"].dropna(), [0.025, 0.975])
    replicated = bool(adjusted_rho < 0 and adjusted_p < 0.05 and adjusted_ci[1] < 0)

    result = pd.DataFrame(
        [
            {
                "platform": "CRC organoids",
                "model": "QC-adjusted rank correlation (V7 frozen)",
                "n": 85,
                "rho": float(organoid_row["qc_adjusted_rank_correlation"]),
                "ci95_low": np.nan,
                "ci95_high": np.nan,
                "p": float(organoid_row["qc_adjusted_p"]),
                "interpretation": "discovery",
            },
            {
                "platform": "DepMap 26Q1 CRC 2D",
                "model": "raw Spearman",
                "n": len(scd),
                "rho": float(raw_test.statistic),
                "ci95_low": float(raw_ci[0]),
                "ci95_high": float(raw_ci[1]),
                "p": float(raw_test.pvalue),
                "interpretation": "descriptive replication",
            },
            {
                "platform": "DepMap 26Q1 CRC 2D",
                "model": "partial Spearman: median+MAD+common-essential PC1-5",
                "n": len(scd),
                "rho": adjusted_rho,
                "ci95_low": float(adjusted_ci[0]),
                "ci95_high": float(adjusted_ci[1]),
                "p": adjusted_p,
                "interpretation": "replicated" if replicated else "not replicated",
            },
        ]
    )
    result.to_csv(args.outdir / "depmap_scd_vps72_results.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "ModelID": crc_complete_ids,
            "SCD_gene_effect": scd,
            "VPS72_gene_effect": vps72,
            "global_median_gene_effect": covariates[:, 0],
            "global_mad_gene_effect": covariates[:, 1],
            **{f"common_essential_PC{i + 1}": covariates[:, i + 2] for i in range(5)},
            "SCD_rank_residual": x_res,
            "VPS72_rank_residual": y_res,
        }
    ).to_csv(args.outdir / "depmap_scd_vps72_crc_lines.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "gene": kept_symbols,
            "completeness": completeness,
            "median_gene_effect": medians,
            "proportion_lt_minus_0_5": prop_dependent,
            "common_essential_selected": common_mask,
        }
    ).to_csv(args.outdir / "depmap_common_essential_selection.tsv.gz", sep="\t", index=False)

    non_crc = (~crc_mask) & np.isfinite(scd_all) & np.isfinite(vps_all)
    pan = np.isfinite(scd_all) & np.isfinite(vps_all)
    context = {
        "crc_raw": {"n": int(complete.sum()), "rho": float(raw_test.statistic), "p": float(raw_test.pvalue)},
        "crc_adjusted": {"n": int(complete.sum()), "rho": adjusted_rho, "permutation_p": float(adjusted_p), "bootstrap_ci95": adjusted_ci.tolist()},
        "non_crc_raw": {
            "n": int(non_crc.sum()),
            "rho": float(stats.spearmanr(scd_all[non_crc], vps_all[non_crc]).statistic),
            "p": float(stats.spearmanr(scd_all[non_crc], vps_all[non_crc]).pvalue),
        },
        "all_2d_raw": {
            "n": int(pan.sum()),
            "rho": float(stats.spearmanr(scd_all[pan], vps_all[pan]).statistic),
            "p": float(stats.spearmanr(scd_all[pan], vps_all[pan]).pvalue),
        },
        "organoid_direction": float(organoid_row["qc_adjusted_rank_correlation"]),
        "replication_success": replicated,
    }
    audit = {
        "release": "DepMap Public 26Q1",
        "seed": SEED,
        "two_d_lines_loaded": len(ids_array),
        "crc_lines_complete": len(scd),
        "genes_unique_retained": len(kept_symbols),
        "common_essential_genes_pre_pca": int(common_mask.sum()),
        "common_essential_variable_genes": int(variable.sum()),
        "adjustment_covariates": ["global median", "global MAD", "common-essential PC1-PC5"],
        "results": context,
        "input_sha256": {
            "gene_effect": sha256(args.gene_effect),
            "model": sha256(args.model),
            "organoid_qc": sha256(args.organoid_qc),
        },
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
    }
    (args.outdir / "depmap_scd_vps72_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

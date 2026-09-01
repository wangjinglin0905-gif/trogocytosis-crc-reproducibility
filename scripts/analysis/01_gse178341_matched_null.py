#!/usr/bin/env python
"""Prespecified matched-null analysis for the frozen GSE178341 module.

The independent unit is the patient. Cell-level data are used only to build
gene matching features and patient pseudobulks. No cell-level P value is used.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats


SEED = 20260901
FROZEN = ["CD4", "PTPRC", "CTLA4", "PDCD1", "HAVCR2", "VSIR", "LAG3", "CD38"]
PANEL = {
    "ANO6", "ATF3", "BCAS1", "C3", "CCR7", "CD109", "CD19", "CD22", "CD24",
    "CD274", "CD38", "CD4", "CD47", "CD80", "CD86", "CDH2", "CEACAM5", "CH25H",
    "CLSTN2", "CTLA4", "CTSE", "EGFR", "ERBB2", "FCGR1A", "FCGR2B", "FCGR3A",
    "HAVCR2", "HLA-DRA", "IL6", "KANK4", "LAG3", "MSLN", "PDCD1", "PTPRC",
    "SCD", "SIGLEC10", "SIRPA", "STAT1", "VSIR",
}
T_CELL_MIDWAY = {"TCD4", "TCD8", "Tgd", "TZBTB16"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def robust_z(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    med = np.nanmedian(values, axis=0)
    q25, q75 = np.nanquantile(values, [0.25, 0.75], axis=0)
    scale = q75 - q25
    fallback = np.nanstd(values, axis=0, ddof=1)
    scale = np.where(scale > 0, scale, fallback)
    scale = np.where(scale > 0, scale, 1.0)
    return (values - med) / scale


def genewise_z(log_cpm: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    means = np.nanmean(log_cpm, axis=1, keepdims=True)
    sds = np.nanstd(log_cpm, axis=1, ddof=1, keepdims=True)
    variable = np.isfinite(sds[:, 0]) & (sds[:, 0] > 0)
    z = np.zeros_like(log_cpm, dtype=float)
    z[variable] = (log_cpm[variable] - means[variable]) / sds[variable]
    return z, variable


def spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    result = stats.spearmanr(x, y, nan_policy="omit")
    return float(result.statistic), float(result.pvalue)


def bootstrap_rho(x: np.ndarray, y: np.ndarray, rng: np.random.Generator, reps: int) -> np.ndarray:
    out = np.empty(reps, dtype=float)
    n = len(x)
    for i in range(reps):
        idx = rng.integers(0, n, n)
        out[i] = stats.spearmanr(x[idx], y[idx]).statistic
    return out[np.isfinite(out)]


def aggregate_all_genes(
    h5_path: Path,
    cells: pd.DataFrame,
    patients: list[str],
    chunk_cells: int,
) -> dict[str, object]:
    patient_index = {patient: i for i, patient in enumerate(patients)}
    with h5py.File(h5_path, "r") as handle:
        group = handle["matrix"]
        n_features, n_cells = [int(x) for x in group["shape"][:]]
        if n_cells != len(cells):
            raise RuntimeError(f"H5 cells={n_cells}, annotations={len(cells)}")
        first_barcodes = [x.decode() for x in group["barcodes"][:10]]
        if first_barcodes != cells["sampleID"].iloc[:10].tolist():
            raise RuntimeError("H5 barcode order does not match annotations")

        feature_names = np.asarray([x.decode() for x in group["features/name"][:]], dtype=str)
        symbols, symbol_counts = np.unique(feature_names, return_counts=True)
        unique_symbol_set = set(symbols[symbol_counts == 1])
        unique_feature_idx = np.asarray(
            [i for i, symbol in enumerate(feature_names) if symbol in unique_symbol_set], dtype=int
        )
        unique_names = feature_names[unique_feature_idx]
        patient_counts = np.zeros((len(unique_names), len(patients)), dtype=np.float64)
        total_counts = np.zeros(len(unique_names), dtype=np.float64)
        detected_cells = np.zeros(len(unique_names), dtype=np.int64)
        tnilc_counts = np.zeros(len(unique_names), dtype=np.float64)
        other_counts = np.zeros(len(unique_names), dtype=np.float64)
        patient_cells = np.zeros(len(patients), dtype=np.int64)
        patient_t_cells = np.zeros(len(patients), dtype=np.int64)
        total_library = 0.0
        tnilc_library = 0.0
        other_library = 0.0
        n_tumour_cells = 0
        n_tnilc_cells = 0

        indptr_ds, indices_ds, data_ds = group["indptr"], group["indices"], group["data"]
        for start in range(0, n_cells, chunk_cells):
            end = min(n_cells, start + chunk_cells)
            pointers = indptr_ds[start : end + 1].astype(np.int64)
            lo, hi = int(pointers[0]), int(pointers[-1])
            local = sparse.csc_matrix(
                (data_ds[lo:hi], indices_ds[lo:hi], pointers - lo),
                shape=(n_features, end - start),
            )
            # Ambiguously duplicated symbols are excluded before matching. This
            # makes cell-level detection exact and avoids double-counting a cell.
            local = local[unique_feature_idx, :].tocsc()
            meta = cells.iloc[start:end]
            tumour_mask = meta["SPECIMEN_TYPE"].eq("T") & meta["PID"].isin(patient_index)
            if not tumour_mask.any():
                continue
            tumour_cols = np.flatnonzero(tumour_mask.to_numpy())
            tumour_local = local[:, tumour_cols]
            tumour_meta = meta.iloc[tumour_cols]
            n_tumour_cells += len(tumour_cols)
            total_counts += np.asarray(tumour_local.sum(axis=1)).ravel()
            detected_cells += np.asarray(tumour_local.getnnz(axis=1)).ravel().astype(np.int64)
            total_library += float(tumour_local.sum())

            tnilc_mask = tumour_meta["clTopLevel"].eq("TNKILC").to_numpy()
            if tnilc_mask.any():
                part = tumour_local[:, np.flatnonzero(tnilc_mask)]
                tnilc_counts += np.asarray(part.sum(axis=1)).ravel()
                tnilc_library += float(part.sum())
                n_tnilc_cells += int(tnilc_mask.sum())
            if (~tnilc_mask).any():
                part = tumour_local[:, np.flatnonzero(~tnilc_mask)]
                other_counts += np.asarray(part.sum(axis=1)).ravel()
                other_library += float(part.sum())

            rows = np.arange(len(tumour_cols), dtype=int)
            cols = np.asarray([patient_index[x] for x in tumour_meta["PID"]], dtype=int)
            design = sparse.coo_matrix(
                (np.ones(len(rows)), (rows, cols)),
                shape=(len(rows), len(patients)),
            ).tocsr()
            patient_counts += (tumour_local @ design).toarray()
            for patient, group_meta in tumour_meta.groupby("PID"):
                j = patient_index[patient]
                patient_cells[j] += len(group_meta)
                patient_t_cells[j] += int(group_meta["clMidwayPr"].isin(T_CELL_MIDWAY).sum())
            print(f"aggregated {start:,}-{end:,}/{n_cells:,}", flush=True)

    return {
        "genes": unique_names,
        "patient_counts": patient_counts,
        "total_counts": total_counts,
        "detected_cells": detected_cells,
        "tnilc_counts": tnilc_counts,
        "other_counts": other_counts,
        "patient_cells": patient_cells,
        "patient_t_cells": patient_t_cells,
        "total_library": total_library,
        "tnilc_library": tnilc_library,
        "other_library": other_library,
        "n_tumour_cells": n_tumour_cells,
        "n_tnilc_cells": n_tnilc_cells,
        "raw_features": n_features,
        "duplicate_symbol_features_excluded": int(n_features - len(unique_names)),
    }


def make_modules(
    target_names: list[str],
    target_features: np.ndarray,
    pool_names: np.ndarray,
    pool_features: np.ndarray,
    k: int,
    reps: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, pd.DataFrame]:
    target_scaled = robust_z(np.vstack([pool_features, target_features]))[-len(target_names) :]
    pool_scaled = robust_z(np.vstack([pool_features, target_features]))[: len(pool_names)]
    distances = np.sqrt(((target_scaled[:, None, :] - pool_scaled[None, :, :]) ** 2).sum(axis=2))
    order = np.argsort(distances, axis=1)
    neighbour_rows: list[dict[str, object]] = []
    for i, target in enumerate(target_names):
        for rank, idx in enumerate(order[i, : max(100, k)], start=1):
            neighbour_rows.append(
                {"target_gene": target, "candidate_gene": pool_names[idx], "rank": rank, "distance": distances[i, idx]}
            )

    modules = np.empty((reps, len(target_names)), dtype=int)
    for rep in range(reps):
        used: set[int] = set()
        for i in range(len(target_names)):
            candidates = [int(x) for x in order[i, :k] if int(x) not in used]
            if not candidates:
                candidates = [int(x) for x in order[i] if int(x) not in used]
            chosen = int(rng.choice(candidates))
            modules[rep, i] = chosen
            used.add(chosen)
    return modules, pd.DataFrame(neighbour_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5", required=True, type=Path)
    parser.add_argument("--clusters", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--v7-scores", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--chunk-cells", type=int, default=10000)
    parser.add_argument("--reps", type=int, default=10000)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)

    clusters = pd.read_csv(args.clusters)
    metadata = pd.read_csv(args.metadata)
    if len(clusters) != len(metadata) or not np.array_equal(clusters["sampleID"], metadata["cellID"]):
        raise RuntimeError("Cluster and metadata tables are not row-aligned")
    cells = pd.concat([clusters, metadata.drop(columns=["cellID"])], axis=1)
    tumour = cells.loc[cells["SPECIMEN_TYPE"].eq("T") & cells["MMRStatus"].isin(["MMRd", "MMRp"])]
    patients = sorted(tumour["PID"].unique())
    archived = pd.read_csv(args.v7_scores, sep="\t").set_index("PID").reindex(patients)

    agg = aggregate_all_genes(args.h5, cells, patients, args.chunk_cells)
    genes = np.asarray(agg["genes"], dtype=str)
    patient_counts = np.asarray(agg["patient_counts"], dtype=float)
    total_counts = np.asarray(agg["total_counts"], dtype=float)
    detected_cells = np.asarray(agg["detected_cells"], dtype=float)
    tnilc_counts = np.asarray(agg["tnilc_counts"], dtype=float)
    other_counts = np.asarray(agg["other_counts"], dtype=float)
    patient_cells = np.asarray(agg["patient_cells"], dtype=float)
    patient_t_cells = np.asarray(agg["patient_t_cells"], dtype=float)

    patient_library = patient_counts.sum(axis=0)
    log_cpm = np.log2((patient_counts + 0.5) / (patient_library[None, :] + 1.0) * 1_000_000.0)
    z, variable = genewise_z(log_cpm)
    patient_presence = (patient_counts > 0).sum(axis=1)
    mean_expression = np.log1p(10000.0 * total_counts / float(agg["total_library"]))
    detection_rate = detected_cells / float(agg["n_tumour_cells"])
    t_cpm = tnilc_counts / float(agg["tnilc_library"]) * 1_000_000.0
    other_cpm = other_counts / float(agg["other_library"]) * 1_000_000.0
    specificity = np.log2((t_cpm + 0.1) / (other_cpm + 0.1))

    feature_table = pd.DataFrame(
        {
            "gene": genes,
            "mean_expression_log_cp10k": mean_expression,
            "cell_detection_rate": detection_rate,
            "tnkilc_vs_other_log2_cpm_ratio": specificity,
            "patients_detected": patient_presence,
            "patient_level_variable": variable,
        }
    )
    feature_table.to_csv(args.outdir / "gse178341_gene_matching_features.tsv", sep="\t", index=False)

    lookup = {gene: i for i, gene in enumerate(genes)}
    missing = [gene for gene in FROZEN if gene not in lookup]
    if missing:
        raise RuntimeError(f"Frozen genes missing: {missing}")
    target_idx = np.asarray([lookup[x] for x in FROZEN])
    eligible = (
        (detection_rate >= 0.0005)
        & (patient_presence >= 10)
        & variable
        & np.isfinite(mean_expression)
        & np.isfinite(specificity)
        & ~np.char.startswith(genes, "MT-")
        & ~np.char.startswith(genes, "RPL")
        & ~np.char.startswith(genes, "RPS")
        & ~np.isin(genes, list(PANEL | set(FROZEN)))
    )
    pool_idx = np.flatnonzero(eligible)
    pool_names = genes[pool_idx]
    all_features = np.column_stack([mean_expression, detection_rate, specificity])
    target_features = all_features[target_idx]
    pool_features = all_features[pool_idx]

    t_fraction = patient_t_cells / patient_cells
    target_score = z[target_idx].mean(axis=0)
    observed_rho, observed_p = spearman(target_score, t_fraction)
    archived_rho, archived_p = spearman(archived["module_all"].to_numpy(), archived["T_cell_fraction_all"].to_numpy())
    archived_score_concordance = spearman(target_score, archived["module_all"].to_numpy())
    boot = bootstrap_rho(target_score, t_fraction, rng, args.reps)

    summary_rows: list[dict[str, object]] = []
    null_rows: list[pd.DataFrame] = []
    module_rows: list[pd.DataFrame] = []
    neighbour_frames: list[pd.DataFrame] = []
    for k in (25, 50, 100):
        modules, neighbours = make_modules(
            FROZEN, target_features, pool_names, pool_features, k, args.reps, rng
        )
        neighbours.insert(0, "k", k)
        neighbour_frames.append(neighbours)
        null_rho = np.empty(args.reps, dtype=float)
        for rep in range(args.reps):
            score = z[pool_idx[modules[rep]]].mean(axis=0)
            null_rho[rep] = stats.spearmanr(score, t_fraction).statistic
        valid = null_rho[np.isfinite(null_rho)]
        p_one = (1 + np.sum(valid >= observed_rho)) / (len(valid) + 1)
        p_two = (1 + np.sum(np.abs(valid) >= abs(observed_rho))) / (len(valid) + 1)
        sd = np.std(valid, ddof=1)
        summary_rows.append(
            {
                "k": k,
                "reps_valid": len(valid),
                "observed_rho": observed_rho,
                "observed_asymptotic_p": observed_p,
                "observed_bootstrap_ci95_low": np.quantile(boot, 0.025),
                "observed_bootstrap_ci95_high": np.quantile(boot, 0.975),
                "null_median": np.median(valid),
                "null_ci95_low": np.quantile(valid, 0.025),
                "null_ci95_high": np.quantile(valid, 0.975),
                "observed_percentile": 100.0 * np.mean(valid <= observed_rho),
                "empirical_p_one_sided": p_one,
                "empirical_p_two_sided": p_two,
                "null_z": (observed_rho - np.mean(valid)) / sd if sd > 0 else np.nan,
            }
        )
        null_rows.append(pd.DataFrame({"k": k, "replicate": np.arange(1, len(valid) + 1), "rho": valid}))
        module_frame = pd.DataFrame(pool_names[modules], columns=FROZEN)
        module_frame.insert(0, "replicate", np.arange(1, args.reps + 1))
        module_frame.insert(0, "k", k)
        module_rows.append(module_frame)

    pd.DataFrame(summary_rows).to_csv(args.outdir / "gse178341_matched_null_summary.tsv", sep="\t", index=False)
    pd.concat(null_rows, ignore_index=True).to_csv(args.outdir / "gse178341_matched_null_rhos.tsv.gz", sep="\t", index=False)
    pd.concat(module_rows, ignore_index=True).to_csv(args.outdir / "gse178341_random_module_members.tsv.gz", sep="\t", index=False)
    pd.concat(neighbour_frames, ignore_index=True).to_csv(args.outdir / "gse178341_matching_neighbours.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "PID": patients,
            "t_cell_fraction": t_fraction,
            "frozen_module_score_complete_gene_cpm": target_score,
            "frozen_module_score_v7_archived": archived["module_all"].to_numpy(),
        }
    ).to_csv(args.outdir / "gse178341_patient_scores_matched_null.tsv", sep="\t", index=False)

    audit = {
        "seed": SEED,
        "patients": len(patients),
        "tumour_cells": int(agg["n_tumour_cells"]),
        "tnkilc_cells": int(agg["n_tnilc_cells"]),
        "unique_symbol_features_retained": len(genes),
        "raw_features": int(agg["raw_features"]),
        "duplicate_symbol_features_excluded": int(agg["duplicate_symbol_features_excluded"]),
        "eligible_pool_genes": len(pool_idx),
        "frozen_genes": FROZEN,
        "observed_rho_complete_gene_cpm": observed_rho,
        "observed_p_complete_gene_cpm": observed_p,
        "observed_bootstrap_ci95": [float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))],
        "archived_v7_rho": archived_rho,
        "archived_v7_p": archived_p,
        "new_vs_archived_score_spearman": {"rho": archived_score_concordance[0], "p": archived_score_concordance[1]},
        "matching_features": ["mean_expression", "cell_detection_rate", "TNKILC_vs_other_specificity"],
        "input_sha256": {
            "h5": sha256(args.h5),
            "clusters": sha256(args.clusters),
            "metadata": sha256(args.metadata),
            "v7_scores": sha256(args.v7_scores),
        },
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
    }
    (args.outdir / "gse178341_matched_null_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

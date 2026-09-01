#!/usr/bin/env python
"""Patient-level recalculation for GSE178341.

The script aggregates the official 10x count matrix by patient and tumour-cell
compartment. It never treats cells as biological replicates. CMScaller template
prediction is performed on epithelial pseudobulk; total-tissue and epithelial
scores are reported side by side to expose composition effects.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import rdata
from scipy import sparse, stats


SEED = 20260831
FROZEN_MODULE = ["CD4", "PTPRC", "CTLA4", "PDCD1", "HAVCR2", "VSIR", "LAG3", "CD38"]
EXHAUSTION_PANEL = ["PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "TOX", "CXCL13"]
T_CELL_MIDWAY = {"TCD4", "TCD8", "Tgd", "TZBTB16"}
COMPARTMENTS = ["all", "Epi", "TNKILC", "Strom"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    pvalues = np.asarray(pvalues, dtype=float)
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def quantile_normalize(matrix: np.ndarray) -> np.ndarray:
    """Column-wise quantile normalization with average handling of ties."""
    matrix = np.asarray(matrix, dtype=float)
    order = np.argsort(matrix, axis=0, kind="mergesort")
    sorted_values = np.take_along_axis(matrix, order, axis=0)
    reference = sorted_values.mean(axis=1)
    output = np.empty_like(matrix)
    for col in range(matrix.shape[1]):
        values = sorted_values[:, col]
        assigned = reference.copy()
        starts = np.r_[0, np.flatnonzero(values[1:] != values[:-1]) + 1]
        ends = np.r_[starts[1:], len(values)]
        for start, end in zip(starts, ends):
            if end - start > 1:
                assigned[start:end] = reference[start:end].mean()
        output[order[:, col], col] = assigned
    return output


def collapse_duplicate_symbols(names: np.ndarray, counts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    names = np.asarray(names, dtype=str)
    unique, inverse = np.unique(names, return_inverse=True)
    collapsed = np.zeros((len(unique), counts.shape[1]), dtype=np.float64)
    np.add.at(collapsed, inverse, counts)
    return unique, collapsed


def aggregate_pseudobulk(
    h5_path: Path,
    metadata: pd.DataFrame,
    patients: list[str],
    chunk_cells: int,
    target_symbols: set[str],
    background_genes: int = 10000,
):
    patient_index = {pid: i for i, pid in enumerate(patients)}
    n_patients = len(patients)
    n_groups = n_patients * len(COMPARTMENTS)
    with h5py.File(h5_path, "r") as handle:
        group = handle["matrix"]
        shape = tuple(int(x) for x in group["shape"][:])
        n_genes, n_cells = shape
        if n_cells != len(metadata):
            raise RuntimeError(f"H5 has {n_cells} cells but metadata has {len(metadata)}")
        first_barcodes = [x.decode() for x in group["barcodes"][:10]]
        if first_barcodes != metadata["sampleID"].iloc[:10].tolist():
            raise RuntimeError("H5 barcode order does not match cluster/metadata order")
        feature_names = np.asarray([x.decode() for x in group["features/name"][:]], dtype=str)
        target_idx = np.flatnonzero(np.isin(feature_names, list(target_symbols)))
        remaining_idx = np.setdiff1d(np.arange(n_genes), target_idx, assume_unique=False)
        rng = np.random.default_rng(SEED)
        background_idx = rng.choice(
            remaining_idx, size=min(background_genes, len(remaining_idx)), replace=False
        )
        selected_idx = np.sort(np.unique(np.r_[target_idx, background_idx]))
        selected_names = feature_names[selected_idx]
        counts = np.zeros((len(selected_idx), n_groups), dtype=np.float64)
        cell_counts = np.zeros(n_groups, dtype=np.int64)
        indptr_ds = group["indptr"]
        indices_ds = group["indices"]
        data_ds = group["data"]
        for start in range(0, n_cells, chunk_cells):
            end = min(n_cells, start + chunk_cells)
            pointers = indptr_ds[start : end + 1].astype(np.int64)
            lo, hi = int(pointers[0]), int(pointers[-1])
            local = sparse.csc_matrix(
                (
                    data_ds[lo:hi],
                    indices_ds[lo:hi],
                    pointers - lo,
                ),
                shape=(n_genes, end - start),
            )
            meta = metadata.iloc[start:end]
            rows: list[int] = []
            cols: list[int] = []
            for local_col, row in enumerate(meta.itertuples(index=False)):
                if row.SPECIMEN_TYPE != "T" or row.PID not in patient_index:
                    continue
                pid_col = patient_index[row.PID]
                rows.append(local_col)
                cols.append(pid_col)
                if row.clTopLevel in COMPARTMENTS[1:]:
                    rows.append(local_col)
                    cols.append(COMPARTMENTS.index(row.clTopLevel) * n_patients + pid_col)
            if rows:
                design = sparse.coo_matrix(
                    (np.ones(len(rows), dtype=np.float64), (rows, cols)),
                    shape=(end - start, n_groups),
                ).tocsr()
                counts += (local[selected_idx, :] @ design).toarray()
                cell_counts += np.asarray(design.sum(axis=0)).ravel().astype(np.int64)
            if (start // chunk_cells) % 20 == 0:
                print(f"aggregated cells {start:,}-{end:,}/{n_cells:,}", flush=True)
    columns = pd.MultiIndex.from_product([COMPARTMENTS, patients], names=["compartment", "PID"])
    return selected_names, pd.DataFrame(counts, index=selected_names, columns=columns), pd.Series(cell_counts, index=columns)


def logcpm(counts: pd.DataFrame) -> pd.DataFrame:
    library = counts.sum(axis=0)
    return np.log2((counts + 0.5).divide(library + 1.0, axis=1) * 1_000_000.0)


def gene_zmean(expression: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    matched = [gene for gene in dict.fromkeys(genes) if gene in expression.index]
    if len(matched) < 2:
        raise RuntimeError(f"Too few matched genes: {matched}")
    values = expression.loc[matched].T
    sd = values.std(axis=0, ddof=1)
    keep = sd > 0
    z = (values.loc[:, keep] - values.loc[:, keep].mean(axis=0)) / sd.loc[keep]
    score = z.mean(axis=1)
    score = (score - score.mean()) / score.std(ddof=1)
    return score.rename("score"), list(z.columns)


def cosine_similarity(vector: np.ndarray, templates: np.ndarray) -> np.ndarray:
    denom = np.linalg.norm(vector) * np.linalg.norm(templates, axis=0)
    return np.divide(vector @ templates, denom, out=np.zeros(templates.shape[1]), where=denom > 0)


def cms_ntp(counts: pd.DataFrame, templates: pd.DataFrame, n_perm: int = 1000) -> pd.DataFrame:
    # CMScaller RNA-seq route: log2(count + 0.25), quantile normalize, row scale.
    raw = counts.to_numpy(dtype=float)
    normalized = quantile_normalize(np.log2(raw + 0.25))
    row_mean = normalized.mean(axis=1, keepdims=True)
    row_sd = normalized.std(axis=1, ddof=1, keepdims=True)
    keep = np.isfinite(row_sd[:, 0]) & (row_sd[:, 0] > 0)
    scaled = (normalized[keep] - row_mean[keep]) / row_sd[keep]
    genes = counts.index.to_numpy()[keep]
    lookup = {gene: i for i, gene in enumerate(genes)}
    use = templates.loc[templates["symbol"].isin(lookup)].copy()
    class_names = sorted(use["class"].astype(str).unique())
    template_rows = np.asarray([lookup[x] for x in use["symbol"]], dtype=int)
    class_index = {name: i for i, name in enumerate(class_names)}
    tmat = np.zeros((len(use), len(class_names)), dtype=float)
    for row, cls in enumerate(use["class"].astype(str)):
        tmat[row, class_index[cls]] = 1.0
    rng = np.random.default_rng(SEED)
    rows = []
    for sample_col, sample in enumerate(counts.columns):
        observed = cosine_similarity(scaled[template_rows, sample_col], tmat)
        best = int(np.argmax(observed))
        sampled = scaled[rng.integers(0, scaled.shape[0], size=(len(use), n_perm)), sample_col]
        null_max = np.max(
            (sampled.T @ tmat)
            / (np.linalg.norm(sampled.T, axis=1, keepdims=True) * np.linalg.norm(tmat, axis=0, keepdims=True)),
            axis=1,
        )
        p_value = (1.0 + np.sum(null_max >= observed[best])) / (n_perm + 1.0)
        distance = np.sqrt(0.5 * (1.0 - np.clip(observed, -1, 1)))
        result = {
            "PID": sample,
            "prediction_raw": class_names[best],
            "p_value": p_value,
            "matched_template_genes": len(use),
        }
        result.update({f"distance_{name}": distance[i] for i, name in enumerate(class_names)})
        rows.append(result)
    output = pd.DataFrame(rows)
    output["FDR"] = bh_fdr(output["p_value"].to_numpy())
    output["prediction_FDR05"] = output["prediction_raw"].where(output["FDR"] <= 0.05, "Unclassified")
    return output


def bootstrap_spearman(x: np.ndarray, y: np.ndarray, n_boot: int = 2000) -> tuple[float, float]:
    rng = np.random.default_rng(SEED)
    values = []
    n = len(x)
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        rho = stats.spearmanr(x[idx], y[idx]).statistic
        if np.isfinite(rho):
            values.append(rho)
    return tuple(np.quantile(values, [0.025, 0.975]))


def residual_rank_correlation(x: pd.Series, y: pd.Series, controls: pd.DataFrame):
    frame = pd.concat([x.rename("x"), y.rename("y"), controls], axis=1).dropna()
    design = np.column_stack(
        [np.ones(len(frame)), *[stats.rankdata(frame[col]) for col in controls.columns]]
    )
    xr = stats.rankdata(frame["x"])
    yr = stats.rankdata(frame["y"])
    xres = xr - design @ np.linalg.lstsq(design, xr, rcond=None)[0]
    yres = yr - design @ np.linalg.lstsq(design, yr, rcond=None)[0]
    result = stats.pearsonr(xres, yres)
    return float(result.statistic), float(result.pvalue), len(frame)


def rank_biserial(u: float, n1: int, n2: int) -> float:
    return 2.0 * u / (n1 * n2) - 1.0


def group_test(values: pd.Series, groups: pd.Series, group1: str, group2: str, label: str) -> dict:
    frame = pd.concat([values.rename("value"), groups.rename("group")], axis=1).dropna()
    first = frame.loc[frame["group"] == group1, "value"].to_numpy()
    second = frame.loc[frame["group"] == group2, "value"].to_numpy()
    test = stats.mannwhitneyu(first, second, alternative="two-sided")
    return {
        "analysis": label,
        "group1": group1,
        "group2": group2,
        "n_group1": len(first),
        "n_group2": len(second),
        "median_group1": float(np.median(first)),
        "median_group2": float(np.median(second)),
        "median_difference": float(np.median(first) - np.median(second)),
        "mannwhitney_u": float(test.statistic),
        "rank_biserial_group1_gt_group2": rank_biserial(test.statistic, len(first), len(second)),
        "p_two_sided": float(test.pvalue),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5", required=True, type=Path)
    parser.add_argument("--clusters", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--cmscaller-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--chunk-cells", type=int, default=10000)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    clusters = pd.read_csv(args.clusters)
    metadata = pd.read_csv(args.metadata)
    if len(clusters) != len(metadata) or not np.array_equal(clusters["sampleID"], metadata["cellID"]):
        raise RuntimeError("Cluster and metadata rows are not aligned")
    cells = pd.concat([clusters, metadata.drop(columns=["cellID"])], axis=1)
    tumour = cells.loc[(cells["SPECIMEN_TYPE"] == "T") & cells["MMRStatus"].isin(["MMRd", "MMRp"])].copy()
    patients = sorted(tumour["PID"].unique())
    patient_meta = tumour.groupby("PID", sort=True).agg(
        MMRStatus=("MMRStatus", "first"),
        age=("Age", "first"),
        sex=("Sex", "first"),
        stage=("TumorStage", "first"),
    ).reindex(patients)

    cms_templates = rdata.read_rda(args.cmscaller_dir / "templates.CMS.rda")["templates.CMS"]
    cms_templates = cms_templates[["class", "symbol"]].drop_duplicates().copy()
    gene_sets = rdata.read_rda(args.cmscaller_dir / "geneSets.CMS.rda")["geneSets.CMS"]
    annotation = rdata.read_rda(args.cmscaller_dir / "anno.orgHs.rda")["anno.orgHs"]
    entrez_to_symbol = (
        annotation.dropna(subset=["entrez", "symbol"])
        .drop_duplicates("entrez")
        .set_index("entrez")["symbol"]
        .astype(str)
        .to_dict()
    )
    emt_genes = [entrez_to_symbol[x] for x in map(str, gene_sets["EMT"]) if x in entrez_to_symbol]
    tgfb_genes = [entrez_to_symbol[x] for x in map(str, gene_sets["TGF-Beta"]) if x in entrez_to_symbol]
    target_symbols = set(FROZEN_MODULE + EXHAUSTION_PANEL + emt_genes + tgfb_genes)
    target_symbols.update(cms_templates["symbol"].astype(str))
    features, pseudobulk_raw, cell_counts = aggregate_pseudobulk(
        args.h5,
        cells,
        patients,
        args.chunk_cells,
        target_symbols=target_symbols,
        background_genes=10000,
    )
    unique_features, collapsed = collapse_duplicate_symbols(features, pseudobulk_raw.to_numpy())
    pseudobulk = pd.DataFrame(collapsed, index=unique_features, columns=pseudobulk_raw.columns)
    np.savez_compressed(
        args.output_dir / "gse178341_patient_pseudobulk_counts.npz",
        counts=collapsed,
        genes=unique_features,
        compartments=np.asarray([x[0] for x in pseudobulk.columns]),
        patients=np.asarray([x[1] for x in pseudobulk.columns]),
    )
    cell_counts.rename("n_cells").to_csv(args.output_dir / "gse178341_patient_compartment_cell_counts.tsv", sep="\t")

    pd.DataFrame(
        [("frozen_8_gene_module", x) for x in FROZEN_MODULE]
        + [("CMScaller_CMS4_template", x) for x in cms_templates.loc[cms_templates["class"].astype(str) == "CMS4", "symbol"]]
        + [("CMScaller_EMT", x) for x in emt_genes]
        + [("CMScaller_TGF_Beta", x) for x in tgfb_genes]
        + [("exhaustion_expression_sensitivity", x) for x in EXHAUSTION_PANEL],
        columns=["set", "gene"],
    ).drop_duplicates().to_csv(args.output_dir / "gse178341_frozen_gene_sets.tsv", sep="\t", index=False)

    patient_scores = patient_meta.copy()
    matched_record = []
    for compartment in ["all", "Epi", "TNKILC"]:
        expression = logcpm(pseudobulk[compartment])
        module, matched = gene_zmean(expression, FROZEN_MODULE)
        patient_scores[f"module_{compartment}"] = module
        matched_record.append({"compartment": compartment, "set": "module", "matched": ";".join(matched)})
        if compartment in ["all", "Epi"]:
            emt, matched_emt = gene_zmean(expression, emt_genes)
            tgfb, matched_tgfb = gene_zmean(expression, tgfb_genes)
            patient_scores[f"EMT_{compartment}"] = emt
            patient_scores[f"TGFb_{compartment}"] = tgfb
            matched_record.extend(
                [
                    {"compartment": compartment, "set": "EMT", "matched": ";".join(matched_emt)},
                    {"compartment": compartment, "set": "TGFb", "matched": ";".join(matched_tgfb)},
                ]
            )
        if compartment == "TNKILC":
            exhaustion, matched_exh = gene_zmean(expression, EXHAUSTION_PANEL)
            patient_scores["exhaustion_expression_TNKILC"] = exhaustion
            matched_record.append({"compartment": compartment, "set": "exhaustion", "matched": ";".join(matched_exh)})

    tumour["is_tcell"] = tumour["clMidwayPr"].isin(T_CELL_MIDWAY)
    tumour["is_exhausted_label"] = tumour["cl295v11SubFull"].str.contains("CXCL13|PDCD1", case=False, na=False)
    tumour["is_stromal"] = tumour["clTopLevel"].eq("Strom")
    comp = tumour.groupby("PID").agg(
        total_cells=("sampleID", "size"),
        t_cells=("is_tcell", "sum"),
        exhausted_label_cells=("is_exhausted_label", "sum"),
        stromal_cells=("is_stromal", "sum"),
    ).reindex(patients)
    patient_scores["T_cell_fraction_all"] = comp["t_cells"] / comp["total_cells"]
    patient_scores["exhausted_label_fraction_T"] = comp["exhausted_label_cells"] / comp["t_cells"]
    patient_scores["stromal_fraction_all"] = comp["stromal_cells"] / comp["total_cells"]

    cms_calls = cms_ntp(pseudobulk["Epi"], cms_templates, n_perm=1000).set_index("PID")
    patient_scores = patient_scores.join(cms_calls)
    patient_scores.to_csv(args.output_dir / "gse178341_patient_scores.tsv", sep="\t", index=True)
    cms_calls.to_csv(args.output_dir / "gse178341_cms_calls.tsv", sep="\t", index=True)
    pd.DataFrame(matched_record).to_csv(args.output_dir / "gse178341_gene_matching.tsv", sep="\t", index=False)

    associations = []
    for compartment in ["all", "Epi"]:
        for pathway in ["TGFb", "EMT"]:
            x = patient_scores[f"module_{compartment}"]
            y = patient_scores[f"{pathway}_{compartment}"]
            complete = pd.concat([x, y], axis=1).dropna()
            result = stats.spearmanr(complete.iloc[:, 0], complete.iloc[:, 1])
            ci_low, ci_high = bootstrap_spearman(complete.iloc[:, 0].to_numpy(), complete.iloc[:, 1].to_numpy())
            associations.append(
                {
                    "analysis": f"module_vs_{pathway}",
                    "compartment": compartment,
                    "n": len(complete),
                    "rho": float(result.statistic),
                    "rho_ci95_low_bootstrap": ci_low,
                    "rho_ci95_high_bootstrap": ci_high,
                    "p_two_sided": float(result.pvalue),
                    "adjustment": "none",
                }
            )
            if compartment == "all":
                partial_r, partial_p, partial_n = residual_rank_correlation(
                    x,
                    y,
                    patient_scores[["T_cell_fraction_all", "stromal_fraction_all"]],
                )
                associations.append(
                    {
                        "analysis": f"module_vs_{pathway}",
                        "compartment": compartment,
                        "n": partial_n,
                        "rho": partial_r,
                        "rho_ci95_low_bootstrap": np.nan,
                        "rho_ci95_high_bootstrap": np.nan,
                        "p_two_sided": partial_p,
                        "adjustment": "rank residuals controlling T-cell and stromal fractions",
                    }
                )
    association_df = pd.DataFrame(associations)
    association_df["FDR_BH_within_four_primary_unadjusted"] = np.nan
    primary_mask = association_df["adjustment"].eq("none")
    association_df.loc[primary_mask, "FDR_BH_within_four_primary_unadjusted"] = bh_fdr(
        association_df.loc[primary_mask, "p_two_sided"].to_numpy()
    )
    association_df.to_csv(args.output_dir / "gse178341_pathway_associations.tsv", sep="\t", index=False)

    group_results = [
        group_test(
            patient_scores["T_cell_fraction_all"], patient_scores["MMRStatus"], "MMRd", "MMRp", "MMR_T_cell_fraction"
        ),
        group_test(
            patient_scores["exhausted_label_fraction_T"],
            patient_scores["MMRStatus"],
            "MMRd",
            "MMRp",
            "MMR_exhausted_author_label_fraction",
        ),
        group_test(
            patient_scores["exhaustion_expression_TNKILC"],
            patient_scores["MMRStatus"],
            "MMRd",
            "MMRp",
            "MMR_exhaustion_expression_sensitivity",
        ),
    ]
    for score_name in ["module_all", "module_Epi"]:
        group_results.append(
            group_test(
                patient_scores[score_name],
                patient_scores["prediction_raw"].map(lambda x: "CMS4" if x == "CMS4" else "CMS1-3"),
                "CMS4",
                "CMS1-3",
                f"CMS4_enrichment_{score_name}_raw_call",
            )
        )
        classified = patient_scores["prediction_FDR05"].ne("Unclassified")
        if classified.sum() >= 10 and (patient_scores.loc[classified, "prediction_FDR05"] == "CMS4").sum() >= 3:
            group_results.append(
                group_test(
                    patient_scores.loc[classified, score_name],
                    patient_scores.loc[classified, "prediction_FDR05"].map(lambda x: "CMS4" if x == "CMS4" else "CMS1-3"),
                    "CMS4",
                    "CMS1-3",
                    f"CMS4_enrichment_{score_name}_FDR05_call",
                )
            )
    group_df = pd.DataFrame(group_results)
    group_df["FDR_BH_all_group_tests"] = bh_fdr(group_df["p_two_sided"].to_numpy())
    group_df.to_csv(args.output_dir / "gse178341_group_comparisons.tsv", sep="\t", index=False)

    sample_flow = {
        "all_cells": int(len(cells)),
        "tumour_cells_with_MMR": int(len(tumour)),
        "patients": len(patients),
        "MMRd_patients": int((patient_meta["MMRStatus"] == "MMRd").sum()),
        "MMRp_patients": int((patient_meta["MMRStatus"] == "MMRp").sum()),
        "epithelial_cells": int((tumour["clTopLevel"] == "Epi").sum()),
        "T_cell_definition_cells": int(tumour["is_tcell"].sum()),
        "cms_raw_calls": cms_calls["prediction_raw"].value_counts().to_dict(),
        "cms_FDR05_calls": cms_calls["prediction_FDR05"].value_counts().to_dict(),
    }
    (args.output_dir / "gse178341_sample_flow.json").write_text(
        json.dumps(sample_flow, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    manifest = {
        "analysis_unit": "patient",
        "seed": SEED,
        "frozen_module": FROZEN_MODULE,
        "exhausted_label_definition": "original-author TNKILC cluster label contains CXCL13 or PDCD1",
        "cms_method": "CMScaller official templates; epithelial pseudobulk; log2(count+0.25); quantile normalization; row scaling; cosine NTP; raw call primary; 1000-permutation confidence sensitivity using a seed-locked 10000-gene null background",
        "score_method": "mean gene-wise z scores from patient-level log2 CPM pseudobulk",
        "inputs": {
            str(path): {"bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in [args.h5, args.clusters, args.metadata, args.cmscaller_dir / "templates.CMS.rda", args.cmscaller_dir / "geneSets.CMS.rda", args.cmscaller_dir / "anno.orgHs.rda"]
        },
        "sample_flow": sample_flow,
        "python": sys.version,
    }
    (args.output_dir / "gse178341_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps({"sample_flow": sample_flow, "associations": associations, "groups": group_results}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

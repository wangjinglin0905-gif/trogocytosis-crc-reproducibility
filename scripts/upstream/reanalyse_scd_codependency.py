from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def bh(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, float)
    order = np.argsort(values)
    q = np.minimum.accumulate((values[order] * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1]
    out = np.empty(len(values))
    out[order] = np.clip(q, 0, 1)
    return out


def residualize_rank(values: np.ndarray, covariates: np.ndarray) -> np.ndarray:
    y = stats.rankdata(values)
    x = np.column_stack([np.ones(len(y)), covariates])
    return y - x @ np.linalg.lstsq(x, y, rcond=None)[0]


def enrichment_score(ranked_weights: np.ndarray, hit: np.ndarray) -> float:
    hit = hit.astype(bool)
    weighted_hits = np.abs(ranked_weights) * hit
    hit_norm = weighted_hits.sum()
    miss_norm = (~hit).sum()
    if hit_norm <= 0 or miss_norm <= 0:
        return 0.0
    running = np.cumsum(weighted_hits / hit_norm - (~hit) / miss_norm)
    return float(running[np.argmax(np.abs(running))])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lfc-table", type=Path, required=True)
    parser.add_argument("--selected-libraries", type=Path, required=True)
    parser.add_argument("--gmt", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--permutations", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260831)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    selected = pd.read_csv(args.selected_libraries)
    samples = selected["sample_ID"].astype(str).tolist()
    sample_index = {s: i for i, s in enumerate(samples)}
    pair_set = set(zip(selected["sample_ID"].astype(str), selected["library"].astype(str)))

    values: dict[str, dict[int, float]] = {}
    for chunk in pd.read_csv(args.lfc_table, sep=r"\s+", chunksize=250_000):
        mask = [
            (str(s), str(lib)) in pair_set
            for s, lib in zip(chunk["sample_ID"], chunk["library"])
        ]
        d = chunk.loc[mask, ["sample_ID", "gene", "LFC"]]
        for sample, gene, value in d.itertuples(index=False):
            values.setdefault(str(gene), {})[sample_index[str(sample)]] = float(value)

    genes = sorted(values)
    matrix = np.full((len(genes), len(samples)), np.nan)
    for i, gene in enumerate(genes):
        for j, value in values[gene].items():
            matrix[i, j] = value
    complete = ~np.isnan(matrix).any(axis=1)
    matrix = matrix[complete]
    genes = [g for g, keep in zip(genes, complete) if keep]
    if "SCD" not in genes:
        raise ValueError("SCD absent from complete matrix")

    # Screen-wide shifts can induce spurious codependency. Adjust ranked values
    # for the per-model genome-wide median and the two screen QC summaries.
    global_median = np.median(matrix, axis=0)
    qc = selected[["AUC_ROC", "AUC_PR"]].apply(pd.to_numeric, errors="coerce")
    covariates = np.column_stack([
        stats.rankdata(global_median),
        stats.rankdata(qc["AUC_ROC"].fillna(qc["AUC_ROC"].median())),
        stats.rankdata(qc["AUC_PR"].fillna(qc["AUC_PR"].median())),
    ])
    scd_index = genes.index("SCD")
    scd_resid = residualize_rank(matrix[scd_index], covariates)
    raw_rho = np.empty(len(genes))
    partial_rho = np.empty(len(genes))
    partial_p = np.empty(len(genes))
    for i in range(len(genes)):
        raw_rho[i] = stats.spearmanr(matrix[i], matrix[scd_index]).statistic
        residual = residualize_rank(matrix[i], covariates)
        partial_rho[i], partial_p[i] = stats.pearsonr(residual, scd_resid)
    ranked = pd.DataFrame({
        "gene": genes,
        "raw_spearman_rho": raw_rho,
        "qc_adjusted_rank_correlation": partial_rho,
        "qc_adjusted_p": partial_p,
    })
    ranked = ranked[ranked["gene"] != "SCD"].sort_values("qc_adjusted_rank_correlation", ascending=False)
    ranked["FDR_BH"] = bh(ranked["qc_adjusted_p"].to_numpy())
    ranked.to_csv(args.outdir / "scd_codependency_qc_adjusted.csv", index=False)

    gene_sets: dict[str, set[str]] = {}
    with args.gmt.open(encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                gene_sets[parts[0]] = set(parts[2:])

    ordered = ranked.copy()
    ordered_genes = ordered["gene"].to_numpy()
    ordered_weights = ordered["qc_adjusted_rank_correlation"].to_numpy(float)
    gsea_rows = []
    for name, members in gene_sets.items():
        hit = np.isin(ordered_genes, list(members))
        size = int(hit.sum())
        if not 15 <= size <= 500:
            continue
        observed = enrichment_score(ordered_weights, hit)
        null = np.empty(args.permutations)
        for b in range(args.permutations):
            null[b] = enrichment_score(ordered_weights, rng.permutation(hit))
        same_sign = null[null >= 0] if observed >= 0 else null[null < 0]
        if observed >= 0:
            extreme = np.sum(same_sign >= observed)
        else:
            extreme = np.sum(same_sign <= observed)
        p = (extreme + 1) / (len(same_sign) + 1)
        denominator = np.mean(np.abs(same_sign)) if len(same_sign) else np.nan
        nes = observed / denominator if denominator and np.isfinite(denominator) else np.nan
        gsea_rows.append({"pathway": name, "size": size, "ES": observed, "NES": nes, "p_value": p})
    gsea = pd.DataFrame(gsea_rows)
    gsea["FDR_BH"] = bh(gsea["p_value"].to_numpy())
    gsea = gsea.sort_values("FDR_BH")
    gsea.to_csv(args.outdir / "scd_hallmark_gsea_qc_adjusted.csv", index=False)

    cholesterol = gsea[gsea["pathway"].str.contains("Cholesterol", case=False, na=False)]
    audit = {
        "n_models": len(samples),
        "n_complete_genes": len(genes),
        "correlation": "rank correlation after residualization for genome-wide median LFC, AUC_ROC and AUC_PR",
        "gsea_permutations": args.permutations,
        "gsea_seed": args.seed,
        "n_hallmark_sets_tested": len(gsea),
        "n_hallmark_fdr_lt_0_05": int((gsea["FDR_BH"] < 0.05).sum()),
        "cholesterol_homeostasis": cholesterol.to_dict(orient="records"),
        "interpretation_gate": "Use in the main text only if a prespecified membrane-lipid pathway is FDR<0.05; otherwise report as a supplemental negative/sensitivity analysis.",
    }
    (args.outdir / "scd_codependency_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

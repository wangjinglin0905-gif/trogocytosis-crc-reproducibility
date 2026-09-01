from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


PANEL = [
    "ANO6", "ATF3", "BCAS1", "C3", "CCR7", "CD109", "CD19", "CD22",
    "CD24", "CD274", "CD38", "CD4", "CD47", "CD80", "CD86", "CDH2",
    "CEACAM5", "CH25H", "CLSTN2", "CTLA4", "CTSE", "EGFR", "ERBB2",
    "FCGR1A", "FCGR2B", "FCGR3A", "HAVCR2", "HLA-DRA", "IL6", "KANK4",
    "LAG3", "MSLN", "PDCD1", "PTPRC", "SCD", "SIGLEC10", "SIRPA",
    "STAT1", "VSIR",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def bootstrap_ci(x: np.ndarray, rng: np.random.Generator, reps: int = 10000) -> tuple[float, float]:
    values = np.empty(reps)
    for i in range(reps):
        values[i] = np.median(rng.choice(x, len(x), replace=True))
    return tuple(np.quantile(values, [0.025, 0.975]))


def bootstrap_median_difference(
    x: np.ndarray, y: np.ndarray, rng: np.random.Generator, reps: int = 10000
) -> tuple[float, float]:
    values = np.empty(reps)
    for i in range(reps):
        values[i] = np.median(rng.choice(x, len(x), replace=True)) - np.median(
            rng.choice(y, len(y), replace=True)
        )
    return tuple(np.quantile(values, [0.025, 0.975]))


def bh(p: pd.Series) -> pd.Series:
    x = p.to_numpy(float)
    order = np.argsort(x)
    q = np.minimum.accumulate((x[order] * len(x) / np.arange(1, len(x) + 1))[::-1])[::-1]
    out = np.empty(len(x))
    out[order] = np.clip(q, 0, 1)
    return pd.Series(out, index=p.index)


def symbol(column: str) -> str:
    return column.split(" (")[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-effect", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--organoid-screen", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--release", default="DepMap Public 24Q4")
    parser.add_argument("--seed", type=int, default=20260831)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    models = pd.read_csv(args.model)
    crc_meta = models[
        models["OncotreePrimaryDisease"].eq("Colorectal Adenocarcinoma")
        & models["ModelType"].eq("Cell Line")
    ].copy()
    crc_ids = set(crc_meta["ModelID"])
    all_2d_ids = set(models.loc[models["ModelType"].eq("Cell Line"), "ModelID"])

    header = pd.read_csv(args.gene_effect, nrows=0).columns.tolist()
    sym_to_col = {symbol(c): c for c in header[1:]}
    present_panel = [g for g in PANEL if g in sym_to_col]
    usecols = [header[0]] + [sym_to_col[g] for g in present_panel]
    targets_all = pd.read_csv(args.gene_effect, usecols=usecols, index_col=0)
    targets_all.index = targets_all.index.astype(str)
    targets_all.columns = [symbol(c) for c in targets_all.columns]
    crc = targets_all.loc[targets_all.index.intersection(crc_ids)].copy()
    other = targets_all.loc[
        targets_all.index.isin(all_2d_ids) & ~targets_all.index.isin(crc_ids)
    ].copy()

    target_rows = []
    for gene in ["SCD", "EGFR"]:
        x = crc[gene].dropna().to_numpy()
        y = other[gene].dropna().to_numpy()
        u, p = stats.mannwhitneyu(x, y, alternative="two-sided")
        ci_lo, ci_hi = bootstrap_ci(x, rng)
        diff_lo, diff_hi = bootstrap_median_difference(x, y, rng)
        target_rows.append({
            "gene": gene,
            "n_crc_2d": len(x),
            "median_crc_2d": np.median(x),
            "median_crc_2d_ci95_lo": ci_lo,
            "median_crc_2d_ci95_hi": ci_hi,
            "pct_crc_lt_minus_0_5": 100 * np.mean(x < -0.5),
            "pct_crc_lt_minus_1": 100 * np.mean(x < -1),
            "n_other_2d": len(y),
            "median_other_2d": np.median(y),
            "median_difference_crc_minus_other": np.median(x) - np.median(y),
            "median_difference_ci95_lo": diff_lo,
            "median_difference_ci95_hi": diff_hi,
            "rank_biserial_crc_greater": 2 * u / (len(x) * len(y)) - 1,
            "mannwhitney_p_crc_vs_other": p,
        })
    target_summary = pd.DataFrame(target_rows)
    target_summary["FDR_BH"] = bh(target_summary["mannwhitney_p_crc_vs_other"])
    target_summary.to_csv(args.outdir / "depmap_target_context_summary.csv", index=False)

    panel_rows = []
    for gene in present_panel:
        x = crc[gene].dropna()
        panel_rows.append({
            "gene": gene,
            "n_crc_2d": len(x),
            "median_crc_2d": x.median(),
            "mean_crc_2d": x.mean(),
            "q1_crc_2d": x.quantile(0.25),
            "q3_crc_2d": x.quantile(0.75),
            "pct_crc_lt_minus_0_5": 100 * (x < -0.5).mean(),
            "pct_crc_lt_minus_1": 100 * (x < -1).mean(),
        })
    panel = pd.DataFrame(panel_rows).sort_values("median_crc_2d")
    panel.to_csv(args.outdir / "depmap_crc_panel_summary.csv", index=False)

    # Keep only CRC rows while streaming the full genome-wide matrix. This avoids
    # loading all cell lines and still permits within-release genome-wide ranks.
    crc_full_chunks = []
    for chunk in pd.read_csv(args.gene_effect, index_col=0, chunksize=128):
        keep = chunk.index.astype(str).isin(crc_ids)
        if keep.any():
            crc_full_chunks.append(chunk.loc[keep])
    crc_full = pd.concat(crc_full_chunks, axis=0)
    genome_median = crc_full.median(axis=0, skipna=True)
    genome_rank = genome_median.rank(method="min", ascending=True)
    genome = pd.DataFrame({
        "column": genome_median.index,
        "gene": [symbol(c) for c in genome_median.index],
        "median_crc_2d": genome_median.values,
        "rank_strongest_dependency": genome_rank.values,
    }).sort_values("rank_strongest_dependency")
    genome["percentile_strongest_dependency"] = 100 * (
        1 - (genome["rank_strongest_dependency"] - 1) / max(1, len(genome) - 1)
    )
    genome.to_csv(args.outdir / "depmap_crc_genome_dependency_rank.csv", index=False)

    organoid = pd.read_csv(args.organoid_screen)
    organoid_summary = organoid[organoid["gene"].isin(present_panel)].groupby("gene").agg(
        n_organoid=("LFC", "count"),
        median_organoid=("LFC", "median"),
        q1_organoid=("LFC", lambda s: s.quantile(0.25)),
        q3_organoid=("LFC", lambda s: s.quantile(0.75)),
        pct_officially_depleted=("is_depleted", lambda s: 100 * np.nanmean(s)),
    ).reset_index()
    cross = organoid_summary.merge(panel, on="gene", how="inner")
    rank_rho, rank_p = stats.spearmanr(cross["median_organoid"], cross["median_crc_2d"])
    cross.to_csv(args.outdir / "depmap_organoid_cross_platform_panel.csv", index=False)

    anchors = cross[cross["gene"].isin(["SCD", "EGFR"])].merge(
        genome[["gene", "rank_strongest_dependency", "percentile_strongest_dependency"]],
        on="gene", how="left",
    )
    anchors.to_csv(args.outdir / "depmap_organoid_anchor_comparison.csv", index=False)

    audit = {
        "release": args.release,
        "crc_definition": "OncotreePrimaryDisease=Colorectal Adenocarcinoma and ModelType=Cell Line",
        "metadata_crc_2d_models": int(len(crc_meta)),
        "crc_2d_models_with_gene_effect": int(len(crc)),
        "other_models_with_gene_effect": int(len(other)),
        "panel_genes_present": present_panel,
        "panel_genes_absent": sorted(set(PANEL) - set(present_panel)),
        "panel_cross_platform_spearman_rho": float(rank_rho),
        "panel_cross_platform_spearman_p": float(rank_p),
        "direct_scale_comparison_warning": "Chronos and organoid LFC scales were not tested as exchangeable; cross-platform comparisons are descriptive and rank-based.",
        "gene_effect_sha256": sha256(args.gene_effect),
        "model_sha256": sha256(args.model),
    }
    (args.outdir / "depmap_reanalysis_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

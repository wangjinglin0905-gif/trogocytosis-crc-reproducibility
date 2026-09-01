#!/usr/bin/env python3
"""Post-hoc balance audit for the prespecified GSE178341 matched-null design.

This script does not change matching or inference. It quantifies how closely the
already generated random modules reproduce the frozen module's three matching
features, both target-by-target and at the module-mean level.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd


FEATURES = [
    "mean_expression_log_cp10k",
    "cell_detection_rate",
    "tnkilc_vs_other_log2_cpm_ratio",
]
PRIMARY_K = 50


def load_source_module(source_script: Path):
    spec = importlib.util.spec_from_file_location("matched_null_source", source_script)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Balance audit for an existing GSE178341 matched-null run."
    )
    parser.add_argument("--analysis-dir", required=True, type=Path)
    parser.add_argument("--source-script", type=Path,
                        default=Path(__file__).with_name("01_gse178341_matched_null.py"))
    parser.add_argument("--outdir", type=Path, default=None)
    args = parser.parse_args()
    analysis = args.analysis_dir.resolve()
    out = (args.outdir or analysis).resolve()
    out.mkdir(parents=True, exist_ok=True)

    source = load_source_module(args.source_script.resolve())
    frozen = list(source.FROZEN)
    panel = set(source.PANEL)
    table = pd.read_csv(analysis / "gse178341_gene_matching_features.tsv", sep="\t")
    variable = table["patient_level_variable"].astype(str).str.lower().eq("true")
    eligible = (
        table["cell_detection_rate"].ge(0.0005)
        & table["patients_detected"].ge(10)
        & variable
        & np.isfinite(table["mean_expression_log_cp10k"])
        & np.isfinite(table["tnkilc_vs_other_log2_cpm_ratio"])
        & ~table["gene"].str.startswith("MT-")
        & ~table["gene"].str.startswith("RPL")
        & ~table["gene"].str.startswith("RPS")
        & ~table["gene"].isin(panel | set(frozen))
    )
    pool = table.loc[eligible].copy()
    targets = table.set_index("gene").loc[frozen].reset_index()
    combined = pd.concat([pool, targets], ignore_index=True)
    scaled = source.robust_z(combined[FEATURES].to_numpy(float))
    scaled_table = combined[["gene"]].copy()
    for j, feature in enumerate(FEATURES):
        scaled_table[f"robust_z__{feature}"] = scaled[:, j]
    scaled_lookup = scaled_table.set_index("gene")
    raw_lookup = table.set_index("gene")

    modules = pd.read_csv(analysis / "gse178341_random_module_members.tsv.gz", sep="\t")
    modules = modules.loc[modules["k"].eq(PRIMARY_K)].copy()
    if len(modules) != 10_000:
        raise ValueError(f"Expected 10,000 primary modules, found {len(modules)}")
    long = modules.melt(
        id_vars=["k", "replicate"], value_vars=frozen,
        var_name="target_gene", value_name="selected_gene"
    )

    target_rows = []
    for target in frozen:
        selected = long.loc[long["target_gene"].eq(target), "selected_gene"]
        for feature in FEATURES:
            raw_values = raw_lookup.loc[selected, feature].to_numpy(float)
            z_values = scaled_lookup.loc[selected, f"robust_z__{feature}"].to_numpy(float)
            target_raw = float(raw_lookup.loc[target, feature])
            target_z = float(scaled_lookup.loc[target, f"robust_z__{feature}"])
            target_rows.append(
                {
                    "target_gene": target,
                    "feature": feature,
                    "target_value": target_raw,
                    "selected_median": float(np.median(raw_values)),
                    "selected_ci95_low": float(np.quantile(raw_values, 0.025)),
                    "selected_ci95_high": float(np.quantile(raw_values, 0.975)),
                    "median_signed_robust_z_mismatch": float(np.median(z_values - target_z)),
                    "median_absolute_robust_z_mismatch": float(np.median(np.abs(z_values - target_z))),
                    "p95_absolute_robust_z_mismatch": float(np.quantile(np.abs(z_values - target_z), 0.95)),
                }
            )
    target_balance = pd.DataFrame(target_rows)

    module_rows = []
    module_feature_means = pd.DataFrame({"replicate": modules["replicate"].to_numpy()})
    for feature in FEATURES:
        observed = float(raw_lookup.loc[frozen, feature].mean())
        gene_to_value = raw_lookup[feature].to_dict()
        selected_flat = modules[frozen].to_numpy().ravel()
        mapped = pd.Series(selected_flat, copy=False).map(gene_to_value).to_numpy(float)
        null_means = mapped.reshape(len(modules), len(frozen)).mean(axis=1)
        module_feature_means[feature] = null_means
        sd = float(np.std(null_means, ddof=1))
        module_rows.append(
            {
                "k": PRIMARY_K,
                "feature": feature,
                "observed_module_mean": observed,
                "null_module_mean_median": float(np.median(null_means)),
                "null_module_mean_ci95_low": float(np.quantile(null_means, 0.025)),
                "null_module_mean_ci95_high": float(np.quantile(null_means, 0.975)),
                "observed_percentile": float(100 * np.mean(null_means <= observed)),
                "standardized_difference_vs_null_mean": float(
                    (observed - np.mean(null_means)) / sd if sd > 0 else np.nan
                ),
            }
        )
        pd.DataFrame(
            {
                "replicate": modules["replicate"].to_numpy(),
                "feature": feature,
                "null_module_mean": null_means,
            }
        ).to_csv(
            out / f"gse178341_matching_balance_null_{feature}.tsv.gz",
            sep="\t", index=False, compression="gzip"
        )
    module_balance = pd.DataFrame(module_rows)

    # Outcome-blind post-hoc sensitivity prompted by the balance audit: retain
    # the 10% of primary modules closest to the frozen module in the three
    # module-mean features. This is explicitly secondary and cannot replace the
    # prespecified k=50 result.
    observed_vector = np.array(
        [float(raw_lookup.loc[frozen, feature].mean()) for feature in FEATURES]
    )
    null_matrix = module_feature_means[FEATURES].to_numpy(float)
    null_sd = np.std(null_matrix, axis=0, ddof=1)
    standardized_distance = np.sqrt(
        np.sum(((null_matrix - observed_vector[np.newaxis, :]) / null_sd[np.newaxis, :]) ** 2, axis=1)
    )
    module_feature_means["standardized_distance_to_frozen"] = standardized_distance
    n_closest = max(1, int(round(0.10 * len(module_feature_means))))
    closest = module_feature_means.nsmallest(n_closest, "standardized_distance_to_frozen")
    rhos = pd.read_csv(analysis / "gse178341_matched_null_rhos.tsv.gz", sep="\t")
    rhos = rhos.loc[rhos["k"].eq(PRIMARY_K), ["replicate", "rho"]]
    closest = closest.merge(rhos, on="replicate", how="left", validate="one_to_one")
    observed_rho = float(
        pd.read_csv(analysis / "gse178341_matched_null_summary.tsv", sep="\t")
        .loc[lambda x: x["k"].eq(PRIMARY_K), "observed_rho"].iloc[0]
    )
    valid_rho = closest["rho"].dropna().to_numpy(float)
    conditioned = {
        "selection": "10% smallest standardized Euclidean distance in three module-mean matching features",
        "post_hoc_sensitivity": True,
        "n_modules": int(len(valid_rho)),
        "observed_rho": observed_rho,
        "null_median": float(np.median(valid_rho)),
        "null_ci95_low": float(np.quantile(valid_rho, 0.025)),
        "null_ci95_high": float(np.quantile(valid_rho, 0.975)),
        "observed_percentile": float(100 * np.mean(valid_rho <= observed_rho)),
        "empirical_p_one_sided": float((1 + np.sum(valid_rho >= observed_rho)) / (1 + len(valid_rho))),
        "empirical_p_two_sided": float((1 + np.sum(np.abs(valid_rho) >= abs(observed_rho))) / (1 + len(valid_rho))),
        "distance_cutoff": float(closest["standardized_distance_to_frozen"].max()),
    }

    target_balance.to_csv(out / "gse178341_matching_balance_by_target.tsv", sep="\t", index=False)
    module_balance.to_csv(out / "gse178341_matching_balance_module_means.tsv", sep="\t", index=False)
    module_feature_means.to_csv(out / "gse178341_matching_balance_module_features.tsv.gz", sep="\t", index=False, compression="gzip")
    closest.to_csv(out / "gse178341_balance_conditioned_modules.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.DataFrame([conditioned]).to_csv(out / "gse178341_balance_conditioned_sensitivity.tsv", sep="\t", index=False)
    audit = {
        "primary_k": PRIMARY_K,
        "n_modules": int(len(modules)),
        "n_eligible_pool_genes_reconstructed": int(len(pool)),
        "features": FEATURES,
        "module_balance": module_balance.to_dict(orient="records"),
        "balance_conditioned_sensitivity": conditioned,
        "note": "Post-hoc balance description only; no matching or inferential parameter was changed.",
    }
    (out / "gse178341_matching_balance_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

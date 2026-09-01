from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats


def bh_adjust(p_values: list[float] | np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.clip(adjusted, 0.0, 1.0)
    return output


def cox_components(beta: np.ndarray, x: np.ndarray, time: np.ndarray, event: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    """Efron partial log-likelihood, gradient and Hessian."""
    eta = np.clip(x @ beta, -50.0, 50.0)
    weight = np.exp(eta)
    p = x.shape[1]
    loglik = 0.0
    gradient = np.zeros(p)
    hessian = np.zeros((p, p))
    for event_time in np.unique(time[event == 1]):
        deaths = (time == event_time) & (event == 1)
        risk = time >= event_time
        d = int(deaths.sum())
        x_death = x[deaths]
        w_risk = weight[risk]
        x_risk = x[risk]
        w_death = weight[deaths]
        s0 = float(w_risk.sum())
        s1 = np.einsum("i,ij->j", w_risk, x_risk)
        s2 = np.einsum("i,ij,ik->jk", w_risk, x_risk, x_risk)
        e0 = float(w_death.sum())
        e1 = np.einsum("i,ij->j", w_death, x_death)
        e2 = np.einsum("i,ij,ik->jk", w_death, x_death, x_death)
        # Sum tied-event linear predictors before converting to a scalar.
        # The former ``float(array).sum()`` form fails under current NumPy.
        loglik += float((x_death @ beta).sum())
        gradient += x_death.sum(axis=0)
        for fraction_index in range(d):
            fraction = fraction_index / d
            denom = s0 - fraction * e0
            first = s1 - fraction * e1
            second = s2 - fraction * e2
            loglik -= np.log(denom)
            gradient -= first / denom
            hessian -= second / denom - np.outer(first, first) / (denom * denom)
    return loglik, gradient, hessian


def approximate_schoenfeld_p(beta: np.ndarray, x: np.ndarray, time: np.ndarray, event: np.ndarray) -> np.ndarray:
    eta = np.clip(x @ beta, -50.0, 50.0)
    weight = np.exp(eta)
    residuals: list[np.ndarray] = []
    event_times: list[float] = []
    for event_time in np.unique(time[event == 1]):
        deaths = (time == event_time) & (event == 1)
        risk = time >= event_time
        d = int(deaths.sum())
        x_death = x[deaths]
        w_risk, x_risk = weight[risk], x[risk]
        w_death = weight[deaths]
        s0 = float(w_risk.sum())
        s1 = np.einsum("i,ij->j", w_risk, x_risk)
        e0 = float(w_death.sum())
        e1 = np.einsum("i,ij->j", w_death, x_death)
        expected = np.mean([(s1 - (j / d) * e1) / (s0 - (j / d) * e0) for j in range(d)], axis=0)
        for row in x_death:
            residuals.append(row - expected)
            event_times.append(float(event_time))
    resid = np.vstack(residuals)
    ranked_time = stats.rankdata(np.asarray(event_times))
    p_values = []
    for column in range(resid.shape[1]):
        test = stats.pearsonr(ranked_time, resid[:, column])
        p_values.append(float(test.pvalue))
    return np.asarray(p_values)


def cox_fit(frame: pd.DataFrame, duration: str, event: str, covariates: list[str]) -> dict[str, float | int]:
    data = frame[[duration, event, *covariates]].dropna().copy()
    data = data[data[duration] > 0].copy()
    x = data[covariates].astype(float).to_numpy()
    time = data[duration].astype(float).to_numpy()
    event_values = data[event].astype(int).to_numpy()
    scale = x.std(axis=0, ddof=0)
    scale[scale == 0] = 1.0
    centre = x.mean(axis=0)
    x_scaled = (x - centre) / scale

    def objective(beta: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
        ll, grad, hess = cox_components(beta, x_scaled, time, event_values)
        return -ll, -grad, -hess

    fitted = optimize.minimize(
        lambda beta: objective(beta)[0], np.zeros(x_scaled.shape[1]),
        jac=lambda beta: objective(beta)[1], hess=lambda beta: objective(beta)[2],
        method="trust-exact", options={"gtol": 1e-9, "maxiter": 500},
    )
    if not fitted.success:
        # SciPy 1.18 may reject an otherwise well-behaved trust-region step
        # when the Hessian approximation cannot predict improvement. BFGS
        # uses the same Efron log-likelihood and analytic gradient and is a
        # stable compatibility fallback.
        fitted = optimize.minimize(
            lambda beta: objective(beta)[0], np.zeros(x_scaled.shape[1]),
            jac=lambda beta: objective(beta)[1], method="BFGS",
            options={"gtol": 1e-8, "maxiter": 2000},
        )
    gradient_norm = float(np.linalg.norm(objective(fitted.x)[1], ord=np.inf))
    if not fitted.success and gradient_norm > 1e-6:
        raise RuntimeError(
            f"Cox fit failed: {fitted.message}; gradient_inf={gradient_norm:.3g}"
        )
    beta_scaled = fitted.x
    _, _, hess_scaled = cox_components(beta_scaled, x_scaled, time, event_values)
    covariance_scaled = np.linalg.pinv(-hess_scaled)
    beta = beta_scaled / scale
    covariance = covariance_scaled / np.outer(scale, scale)
    se = np.sqrt(np.diag(covariance))
    module_index = covariates.index("module_z")
    z = beta[module_index] / se[module_index]
    p_value = float(2 * stats.norm.sf(abs(z)))
    ph = approximate_schoenfeld_p(beta_scaled, x_scaled, time, event_values)
    return {
        "n": int(len(data)),
        "events": int(event_values.sum()),
        "hr": float(np.exp(beta[module_index])),
        "ci95_low": float(np.exp(beta[module_index] - 1.96 * se[module_index])),
        "ci95_high": float(np.exp(beta[module_index] + 1.96 * se[module_index])),
        "p": p_value,
        "module_ph_p_rank_approx": float(ph[module_index]),
        "global_ph_min_p_rank_approx": float(ph.min()),
    }


def exact_sign_flip_p(differences: np.ndarray) -> float:
    differences = np.asarray(differences, dtype=float)
    observed = abs(float(differences.mean()))
    values = []
    for signs in itertools.product((-1.0, 1.0), repeat=len(differences)):
        values.append(abs(float(np.mean(differences * np.asarray(signs)))))
    return float(np.mean(np.asarray(values) >= observed - 1e-15))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    root = args.root
    args.outdir.mkdir(parents=True, exist_ok=True)
    results: dict[str, object] = {}

    org = root / "analysis" / "organoid"
    panel = pd.read_csv(org / "organoid_panel_dependency_summary.csv")
    screen = pd.read_csv(org / "organoid_crc_best_library_screen.csv.gz")
    ann = pd.read_csv(org / "organoid_model_annotations.csv")
    screen = screen.merge(
        ann[["sample_ID", "individual_ID", "msStatus", "CMS_prediction", "KRAS_class"]],
        on="sample_ID", how="left", validate="many_to_one"
    )
    anchor_rows: list[dict[str, object]] = []
    for gene in ("SCD", "EGFR"):
        values = screen.loc[screen["gene"] == gene, "LFC"].dropna()
        row = panel.set_index("gene").loc[gene]
        anchor_rows.append({
            "gene": gene,
            "n": int(len(values)),
            "median_lfc_recomputed": float(values.median()),
            "median_lfc_stored": float(row["median_LFC"]),
            "official_pct_depleted": float(row["pct_depleted"]),
        })
    results["organoid_anchors"] = anchor_rows
    scd = screen[screen["gene"] == "SCD"]
    cms_groups = [scd.loc[scd["CMS_prediction"] == cms, "LFC"].dropna() for cms in ("CMS1", "CMS2", "CMS3", "CMS4")]
    results["organoid_scd_cms"] = {
        "counts": [int(len(x)) for x in cms_groups],
        "medians": [float(x.median()) for x in cms_groups],
        "kruskal_p": float(stats.kruskal(*cms_groups).pvalue),
    }
    msi_tests = []
    for gene in sorted(panel["gene"]):
        subset = screen[screen["gene"] == gene]
        msi = subset.loc[subset["msStatus"] == "MSI", "LFC"].dropna()
        mss = subset.loc[subset["msStatus"] == "MSS", "LFC"].dropna()
        test = stats.mannwhitneyu(msi, mss, alternative="two-sided", method="asymptotic")
        msi_tests.append({"gene": gene, "n_msi": len(msi), "n_mss": len(mss), "p": float(test.pvalue)})
    q = bh_adjust([x["p"] for x in msi_tests])
    for row, value in zip(msi_tests, q):
        row["q"] = float(value)
    results["organoid_msi"] = {
        "n_msi": int(ann["msStatus"].eq("MSI").sum()),
        "n_mss": int(ann["msStatus"].eq("MSS").sum()),
        "genes": {x["gene"]: x for x in msi_tests if x["gene"] in {"SCD", "EGFR"}},
        "n_fdr_lt_005": int(sum(x["q"] < 0.05 for x in msi_tests)),
    }

    dep = root / "analysis" / "depmap_26Q1"
    cross = pd.read_csv(dep / "depmap_organoid_cross_platform_panel.csv")
    full = stats.spearmanr(cross["median_organoid"], cross["median_crc_2d"])
    legacy = cross[cross["gene"] != "CD274"]
    legacy_test = stats.spearmanr(legacy["median_organoid"], legacy["median_crc_2d"])
    context = pd.read_csv(dep / "depmap_target_context_summary.csv").set_index("gene")
    results["cross_platform"] = {
        "n_full": int(len(cross)), "rho_full": float(full.statistic), "p_full": float(full.pvalue),
        "n_legacy": int(len(legacy)), "rho_legacy": float(legacy_test.statistic), "p_legacy": float(legacy_test.pvalue),
        "anchors": {
            gene: {k: float(context.loc[gene, k]) for k in (
                "median_crc_2d", "pct_crc_lt_minus_1", "median_difference_crc_minus_other",
                "median_difference_ci95_lo", "median_difference_ci95_hi", "mannwhitney_p_crc_vs_other"
            )} for gene in ("SCD", "EGFR")
        },
    }

    codep = root / "analysis" / "scd_codependency_qc"
    rank = pd.read_csv(codep / "scd_codependency_qc_adjusted.csv")
    gsea = pd.read_csv(codep / "scd_hallmark_raw_vs_qc.csv")
    targets = rank[rank["gene"].isin(["VPS72", "CDX1", "WLS"])].set_index("gene")
    results["scd_codependency"] = {
        "n_genes": int(len(rank)),
        "n_gene_fdr_lt_005": int((rank["FDR_BH"] < 0.05).sum()),
        "target_rows": targets.to_dict(orient="index"),
        "n_hallmark_fdr_lt_005": int((gsea["FDR_BH"] < 0.05).sum()),
        "dna_repair": gsea[gsea["pathway"] == "HALLMARK_DNA_REPAIR"].set_index("ranking").to_dict(orient="index"),
        "cholesterol_homeostasis": gsea[gsea["pathway"] == "HALLMARK_CHOLESTEROL_HOMEOSTASIS"].set_index("ranking").to_dict(orient="index"),
    }

    gse = pd.read_csv(root / "analysis" / "gse39582_recalculation" / "gse39582_module_clinical.tsv", sep="\t")
    gse_univ = cox_fit(gse, "rfs_months", "rfs_event", ["module_z"])
    gse_adj = cox_fit(gse, "rfs_months", "rfs_event", ["module_z", "age", "mmrd", "chemotherapy_yes", "stage_2", "stage_3", "stage_4"])
    results["gse39582"] = {"rfs_univariable": gse_univ, "rfs_adjusted": gse_adj, "os_analysable": False}

    tcga = pd.read_csv(root / "analysis" / "tcga_cbioportal_592" / "tcga_cbioportal_module_clinical.tsv", sep="\t")
    tcga_univ = cox_fit(tcga, "os_months", "os_event", ["module_z"])
    tcga_stage = pd.get_dummies(tcga, columns=["stage_num"], prefix="stage", dtype=float)
    stage_columns = [column for column in ("stage_2.0", "stage_3.0", "stage_4.0") if column in tcga_stage]
    tcga_adj = cox_fit(tcga_stage, "os_months", "os_event", ["module_z", "age", *stage_columns])
    results["tcga"] = {
        "expression_n": int(tcga["module_z"].notna().sum()),
        "os_univariable": tcga_univ,
        "os_age_stage_adjusted": tcga_adj,
    }

    leuko = root / "analysis" / "gse146771"
    rates = pd.read_csv(leuko / "gse146771_candidate_rates_by_patient_tissue_lineage.csv")
    all_rates = rates[rates["lineage"] == "All leukocytes"]
    paired = all_rates.pivot(index="Sample", columns="Tissue", values="rate").dropna(subset=["N", "T"])
    differences = (paired["T"] - paired["N"]).to_numpy()
    tissue = pd.read_csv(leuko / "gse146771_candidate_tissue_descriptive.csv").set_index("tissue")
    sensitivity = pd.read_csv(leuko / "gse146771_threshold_sensitivity.csv")
    results["gse146771"] = {
        "paired_patients": int(len(paired)),
        "mean_difference_T_minus_N": float(differences.mean()),
        "exact_sign_flip_p": exact_sign_flip_p(differences),
        "tissue_counts": tissue[["n_cells", "candidates", "rate"]].to_dict(orient="index"),
        "threshold_candidate_range": {
            str(threshold): [int(group["total_candidates"].min()), int(group["total_candidates"].max())]
            for threshold, group in sensitivity.groupby("CEACAM5_TPM_threshold")
        },
    }

    g178 = pd.read_csv(root / "analysis" / "gse178341_recalculation" / "gse178341_patient_scores.tsv", sep="\t")
    cms4 = g178["prediction_raw"] == "CMS4"
    cms_whole = stats.mannwhitneyu(g178.loc[cms4, "module_all"], g178.loc[~cms4, "module_all"], alternative="two-sided")
    cms_epi = stats.mannwhitneyu(g178.loc[cms4, "module_Epi"], g178.loc[~cms4, "module_Epi"], alternative="two-sided")
    correlations = {}
    for a, b in (("module_all", "TGFb_all"), ("module_all", "EMT_all"), ("module_Epi", "TGFb_Epi"), ("module_Epi", "EMT_Epi"), ("module_all", "T_cell_fraction_all"), ("module_Epi", "T_cell_fraction_all")):
        test = stats.spearmanr(g178[a], g178[b], nan_policy="omit")
        correlations[f"{a}_vs_{b}"] = {"rho": float(test.statistic), "p": float(test.pvalue)}
    mmrd = g178["MMRStatus"] == "MMRd"
    mmr_tests = {}
    for column in ("T_cell_fraction_all", "exhausted_label_fraction_T", "exhaustion_expression_TNKILC"):
        test = stats.mannwhitneyu(g178.loc[mmrd, column], g178.loc[~mmrd, column], alternative="two-sided")
        mmr_tests[column] = {
            "median_mmrd": float(g178.loc[mmrd, column].median()),
            "median_mmrp": float(g178.loc[~mmrd, column].median()),
            "p": float(test.pvalue),
        }
    q = bh_adjust([x["p"] for x in mmr_tests.values()])
    for row, value in zip(mmr_tests.values(), q):
        row["q"] = float(value)
    results["gse178341"] = {
        "patients": int(g178["PID"].nunique()),
        "mmrd": int(mmrd.sum()), "mmrp": int((~mmrd).sum()),
        "cms_raw_counts": g178["prediction_raw"].value_counts().to_dict(),
        "cms4_whole_p": float(cms_whole.pvalue), "cms4_epithelial_p": float(cms_epi.pvalue),
        "correlations": correlations, "mmr_tests": mmr_tests,
    }

    output = args.outdir / "v7_key_results_recalculated.json"
    output.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(results, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

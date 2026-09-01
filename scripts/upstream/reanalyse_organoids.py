from __future__ import annotations

import argparse
import json
import re
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

MSI_FAMILY = [
    "SCD", "EGFR", "ERBB2", "CH25H", "STAT1", "ATF3", "ANO6",
    "CEACAM5", "CD274", "HAVCR2", "CD47",
]

DISPLAY = {"SCD": "SCD"}


def bh(values: pd.Series) -> pd.Series:
    arr = np.asarray(values, dtype=float)
    order = np.argsort(arr)
    ranked = arr[order]
    n = len(arr)
    adjusted = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    out = np.empty(n)
    out[order] = adjusted
    return pd.Series(out, index=values.index)


def bootstrap_median_diff(x: np.ndarray, y: np.ndarray, rng: np.random.Generator,
                          reps: int = 10000) -> tuple[float, float]:
    diffs = np.empty(reps, dtype=float)
    for i in range(reps):
        diffs[i] = np.median(rng.choice(x, len(x), replace=True)) - np.median(
            rng.choice(y, len(y), replace=True)
        )
    return tuple(np.quantile(diffs, [0.025, 0.975]))


def collapse_nonnull(series: pd.Series) -> tuple[object, bool]:
    values = pd.unique(series.dropna().astype(str).str.strip())
    values = [v for v in values if v and v.lower() not in {"nan", "none"}]
    if not values:
        return np.nan, False
    if len(values) == 1:
        return values[0], False
    return "|".join(sorted(values)), True


def classify_kras(variants: pd.Series) -> tuple[str, bool, str]:
    vals = sorted({str(v) for v in variants.dropna() if str(v).lower() != "nan"})
    joined = ";".join(vals)
    classes = []
    for label, pattern in [("G12X", r"G12"), ("G13X", r"G13"), ("Q61X", r"Q61")]:
        if re.search(pattern, joined, flags=re.I):
            classes.append(label)
    other_alteration = bool(joined) and not classes
    if len(classes) == 1:
        return classes[0], False, joined
    if len(classes) > 1:
        return "multiple", True, joined
    return ("other", False, joined) if other_alteration else ("WT", False, joined)


def read_screen_table(path: Path, selected_pairs: pd.MultiIndex, value_column: str) -> pd.DataFrame:
    kept = []
    for chunk in pd.read_csv(path, sep=r"\s+", chunksize=250_000):
        idx = pd.MultiIndex.from_frame(chunk[["sample_ID", "library"]])
        subset = chunk.loc[idx.isin(selected_pairs), ["sample_ID", "library", "gene", value_column]]
        kept.append(subset)
    out = pd.concat(kept, ignore_index=True)
    if out.duplicated(["sample_ID", "library", "gene"]).any():
        raise ValueError(f"duplicate selected records detected in {path.name}")
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260831)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    src = args.source
    rng = np.random.default_rng(args.seed)

    table2 = src / "supplementary_table_2_revision.xlsx"
    table3 = src / "supplementary_table_3_revision.xlsx"
    table4 = src / "supplementary_table_4_revision.xlsx"
    table5 = src / "supplementary_table_5_revision.csv"
    table6 = src / "supplementary_table_6_revision.csv"
    table7 = src / "supplementary_table_7_revision.xlsx"
    table8 = src / "supplementary_table_8_revision.xlsx"

    availability = pd.read_excel(table2, sheet_name="data_availability")
    crc_models = set(
        availability.loc[
            (availability["primary_tumour_type"] == "Colorectal")
            & (availability["CRISPR_available"].astype(str).str.lower() == "yes"),
            "sample_ID",
        ].astype(str)
    )

    qc = pd.read_excel(table4, sheet_name="summary_QC_metrics")
    qc_crc = qc[qc["sample_ID"].astype(str).isin(crc_models)].copy()
    # Source definition: order=1 is poorest, hence maximum order is best.
    best = (
        qc_crc.sort_values(["sample_ID", "order", "AUC_ROC", "AUC_PR"], ascending=[True, False, False, False])
        .drop_duplicates("sample_ID", keep="first")
        .reset_index(drop=True)
    )
    if len(best) != 85 or best["sample_ID"].nunique() != 85:
        raise ValueError(f"expected 85 independent CRC models, got {len(best)}")
    best.to_csv(args.outdir / "organoid_selected_libraries.csv", index=False)
    selected_pairs = pd.MultiIndex.from_frame(best[["sample_ID", "library"]])

    lfc = read_screen_table(table6, selected_pairs, "LFC")
    depleted = read_screen_table(table5, selected_pairs, "is_depleted")
    screen = lfc.merge(depleted, on=["sample_ID", "library", "gene"], how="left", validate="one_to_one")
    screen["is_depleted"] = pd.to_numeric(screen["is_depleted"], errors="coerce")
    screen.to_csv(args.outdir / "organoid_crc_best_library_screen.csv.gz", index=False)

    # Sample-level annotations with conflict flags.
    annotations = availability.loc[availability["sample_ID"].astype(str).isin(best["sample_ID"]),
                                   ["sample_ID", "individual_ID", "primary_tumour_type"]].copy()
    annotations = annotations.drop_duplicates("sample_ID")

    wgs = pd.read_excel(table3, sheet_name="data_obtained_from_WGS")
    wgs = wgs[
        wgs["sample_ID"].astype(str).isin(best["sample_ID"])
        & wgs["model"].astype(str).str.lower().str.startswith("organoid")
    ]
    wgs_rows = []
    for sample_id, group in wgs.groupby("sample_ID"):
        ms, conflict = collapse_nonnull(group["msStatus"])
        wgs_rows.append({"sample_ID": sample_id, "msStatus": ms, "msStatus_conflict": conflict})
    annotations = annotations.merge(pd.DataFrame(wgs_rows), on="sample_ID", how="left")

    cms = pd.read_excel(table3, sheet_name="CMS_CRIS_subtypes")
    cms_rows = []
    for sample_id, group in cms[cms["sample_ID"].astype(str).isin(best["sample_ID"])].groupby("sample_ID"):
        cv, cc = collapse_nonnull(group["CMS_prediction"])
        rv, rc = collapse_nonnull(group["CRIS_prediction"])
        cms_rows.append({"sample_ID": sample_id, "CMS_prediction": cv,
                         "CMS_conflict": cc, "CRIS_prediction": rv, "CRIS_conflict": rc})
    annotations = annotations.merge(pd.DataFrame(cms_rows), on="sample_ID", how="left")

    drivers = pd.read_excel(table3, sheet_name="genomic_driver_events")
    kras = drivers[
        drivers["sample_ID"].astype(str).isin(best["sample_ID"])
        & (drivers["gene"].astype(str).str.upper() == "KRAS")
        & drivers["model"].astype(str).str.lower().str.startswith("organoid")
    ].copy()
    kras_rows = []
    for sample_id in best["sample_ID"]:
        group = kras[kras["sample_ID"] == sample_id]
        source = pd.concat([group.get("variant", pd.Series(dtype=object)),
                            group.get("Genomic_alterations", pd.Series(dtype=object))])
        klass, conflict, raw = classify_kras(source)
        kras_rows.append({"sample_ID": sample_id, "KRAS_class": klass,
                          "KRAS_conflict": conflict, "KRAS_source": raw})
    annotations = annotations.merge(pd.DataFrame(kras_rows), on="sample_ID", how="left")
    annotations.to_csv(args.outdir / "organoid_model_annotations.csv", index=False)

    screen = screen.merge(
        annotations[["sample_ID", "individual_ID", "msStatus", "CMS_prediction", "CRIS_prediction", "KRAS_class"]],
        on="sample_ID", how="left", validate="many_to_one"
    )
    screen["display_gene"] = screen["gene"].replace(DISPLAY)

    full_rank = (
        screen.groupby("gene")
        .agg(n=("LFC", "count"), median_LFC=("LFC", "median"), mean_LFC=("LFC", "mean"),
             q1_LFC=("LFC", lambda x: x.quantile(0.25)), q3_LFC=("LFC", lambda x: x.quantile(0.75)),
             pct_depleted=("is_depleted", lambda x: 100 * np.nanmean(x)))
        .reset_index()
        .sort_values("median_LFC", ascending=True)
        .reset_index(drop=True)
    )
    full_rank["rank_most_dependent"] = np.arange(1, len(full_rank) + 1)
    full_rank["percentile_most_dependent"] = 100 * (1 - (full_rank["rank_most_dependent"] - 1) / max(1, len(full_rank) - 1))
    full_rank.to_csv(args.outdir / "organoid_full_gene_dependency_rank.csv", index=False)

    panel = full_rank[full_rank["gene"].isin(PANEL)].copy()
    panel["panel_gene"] = panel["gene"]
    panel.to_csv(args.outdir / "organoid_panel_dependency_summary.csv", index=False)

    # Official differential-dependency and core-fitness classifications.
    diffdep = pd.read_excel(table7, sheet_name="diff_dep_analysis_colorectal")
    diffdep_panel = diffdep[diffdep["gene"].isin(PANEL)].copy()
    diffdep_panel.to_csv(args.outdir / "organoid_official_diffdep_panel.csv", index=False)
    core = pd.read_excel(table4, sheet_name="organoid_core_fitness_genes")
    core_targets = core[core["gene"].isin(["SCD", "EGFR"])].copy()
    core_targets.to_csv(args.outdir / "organoid_core_fitness_SCD_EGFR.csv", index=False)
    pathways = pd.read_excel(table4, sheet_name="pathway_enrichment_analysis")
    lipid_mask = pathways["pathway"].astype(str).str.contains("sterol|cholesterol|lipid|fatty acid", case=False, regex=True)
    pathways[lipid_mask].to_csv(args.outdir / "organoid_core_fitness_lipid_pathways.csv", index=False)

    # Model-level MSI tests.
    msi_rows = []
    for gene in MSI_FAMILY:
        d = screen[(screen["gene"] == gene) & screen["msStatus"].isin(["MSI", "MSS"])].copy()
        x = d.loc[d["msStatus"] == "MSI", "LFC"].dropna().to_numpy()
        y = d.loc[d["msStatus"] == "MSS", "LFC"].dropna().to_numpy()
        if len(x) == 0 or len(y) == 0:
            continue
        u, p = stats.mannwhitneyu(x, y, alternative="two-sided")
        ci_lo, ci_hi = bootstrap_median_diff(x, y, rng)
        msi_rows.append({
            "gene": gene, "n_MSI": len(x), "n_MSS": len(y),
            "median_MSI": np.median(x), "median_MSS": np.median(y),
            "median_difference_MSI_minus_MSS": np.median(x) - np.median(y),
            "median_diff_CI95_lo": ci_lo, "median_diff_CI95_hi": ci_hi,
            "mannwhitney_U": u, "rank_biserial_MSI_greater": 2 * u / (len(x) * len(y)) - 1,
            "p_value": p,
        })
    msi = pd.DataFrame(msi_rows)
    msi["FDR_BH"] = bh(msi["p_value"])
    msi.to_csv(args.outdir / "organoid_MSI_model_level_tests.csv", index=False)

    # Descriptive subtype summaries plus limited omnibus tests for the two anchors.
    strata_rows = []
    omnibus_rows = []
    for gene in ["SCD", "EGFR"]:
        for variable, allowed in [
            ("CMS_prediction", ["CMS1", "CMS2", "CMS3", "CMS4"]),
            ("KRAS_class", ["WT", "G12X", "G13X", "Q61X", "other"]),
        ]:
            d = screen[(screen["gene"] == gene) & screen[variable].isin(allowed)]
            groups = []
            for level in allowed:
                vals = d.loc[d[variable] == level, "LFC"].dropna().to_numpy()
                if len(vals):
                    groups.append(vals)
                    strata_rows.append({"gene": gene, "stratum_variable": variable,
                                        "stratum": level, "n": len(vals),
                                        "median_LFC": np.median(vals),
                                        "q1_LFC": np.quantile(vals, 0.25),
                                        "q3_LFC": np.quantile(vals, 0.75)})
            if len(groups) >= 2:
                h, p = stats.kruskal(*groups)
                omnibus_rows.append({"gene": gene, "stratum_variable": variable,
                                     "kruskal_H": h, "p_value": p})
    strata = pd.DataFrame(strata_rows)
    strata.to_csv(args.outdir / "organoid_strata_descriptive.csv", index=False)
    omnibus = pd.DataFrame(omnibus_rows)
    omnibus["FDR_BH"] = bh(omnibus["p_value"])
    omnibus.to_csv(args.outdir / "organoid_strata_omnibus_tests.csv", index=False)

    # Donor-level sensitivity: retain the highest-QC model for each donor.
    donor_choice = (
        best.merge(annotations[["sample_ID", "individual_ID"]], on="sample_ID", how="left", validate="one_to_one")
        .sort_values(["individual_ID", "order", "AUC_ROC", "AUC_PR", "sample_ID"], ascending=[True, False, False, False, True])
        .drop_duplicates("individual_ID", keep="first")
    )
    donor_screen = screen[screen["sample_ID"].isin(donor_choice["sample_ID"])].copy()
    donor_msi_rows = []
    for gene in MSI_FAMILY:
        d = donor_screen[(donor_screen["gene"] == gene) & donor_screen["msStatus"].isin(["MSI", "MSS"])]
        x = d.loc[d["msStatus"] == "MSI", "LFC"].dropna().to_numpy()
        y = d.loc[d["msStatus"] == "MSS", "LFC"].dropna().to_numpy()
        if len(x) == 0 or len(y) == 0:
            continue
        u, p = stats.mannwhitneyu(x, y, alternative="two-sided")
        donor_msi_rows.append({
            "gene": gene,
            "n_MSI_donors": len(x),
            "n_MSS_donors": len(y),
            "median_MSI": np.median(x),
            "median_MSS": np.median(y),
            "median_difference_MSI_minus_MSS": np.median(x) - np.median(y),
            "mannwhitney_U": u,
            "rank_biserial_MSI_greater": 2 * u / (len(x) * len(y)) - 1,
            "p_value": p,
        })
    donor_msi = pd.DataFrame(donor_msi_rows)
    donor_msi["FDR_BH"] = bh(donor_msi["p_value"])
    donor_msi.to_csv(args.outdir / "organoid_donor_sensitivity_MSI_tests.csv", index=False)

    donor_omnibus_rows = []
    for gene in ["SCD", "EGFR"]:
        for variable, allowed in [
            ("CMS_prediction", ["CMS1", "CMS2", "CMS3", "CMS4"]),
            ("KRAS_class", ["WT", "G12X", "G13X", "Q61X", "other"]),
        ]:
            d = donor_screen[(donor_screen["gene"] == gene) & donor_screen[variable].isin(allowed)]
            groups = [d.loc[d[variable] == level, "LFC"].dropna().to_numpy() for level in allowed]
            groups = [group for group in groups if len(group)]
            if len(groups) >= 2:
                h, p = stats.kruskal(*groups)
                donor_omnibus_rows.append({
                    "gene": gene,
                    "stratum_variable": variable,
                    "n_donors": int(d["individual_ID"].nunique()),
                    "kruskal_H": h,
                    "p_value": p,
                })
    donor_omnibus = pd.DataFrame(donor_omnibus_rows)
    donor_omnibus["FDR_BH"] = bh(donor_omnibus["p_value"])
    donor_omnibus.to_csv(args.outdir / "organoid_donor_sensitivity_strata_omnibus_tests.csv", index=False)

    # Expression summary for panel genes across independent CRC models with RNA-seq.
    expression = pd.read_excel(table3, sheet_name="RNAseq_TPM")
    expression = expression[expression["gene"].isin(PANEL)].copy()
    crc_columns = [c for c in best["sample_ID"] if c in expression.columns]
    expr_long = expression.set_index("gene")[crc_columns].T.reset_index(names="sample_ID").melt(
        id_vars="sample_ID", var_name="gene", value_name="TPM"
    )
    expr_summary = expr_long.groupby("gene").agg(
        n=("TPM", "count"), median_TPM=("TPM", "median"),
        pct_TPM_gt_0_5=("TPM", lambda x: 100 * np.mean(x > 0.5)),
        pct_TPM_gt_1=("TPM", lambda x: 100 * np.mean(x > 1)),
    ).reset_index()
    expr_summary.to_csv(args.outdir / "organoid_panel_expression_summary.csv", index=False)

    biomarkers = pd.read_excel(table8, sheet_name="significant_biomarker_assoc")
    bm = biomarkers[
        biomarkers["Gene"].isin(["SCD", "EGFR", "CEACAM5"])
        | biomarkers["Feature"].astype(str).str.contains("HSP90AB1|LOX", case=False, regex=True)
    ].copy()
    bm.to_csv(args.outdir / "organoid_claimed_biomarker_records.csv", index=False)

    audit = {
        "qc_order_rule": "maximum order retained; source defines order=1 as poorest",
        "n_crc_models": int(len(best)),
        "n_crc_donors": int(annotations["individual_ID"].nunique()),
        "n_donors_with_multiple_models": int((annotations.groupby("individual_ID")["sample_ID"].nunique() > 1).sum()),
        "donor_sensitivity_models_retained": int(len(donor_choice)),
        "n_crc_models_msi": int((annotations["msStatus"] == "MSI").sum()),
        "n_crc_models_mss": int((annotations["msStatus"] == "MSS").sum()),
        "n_cms_annotated": int(annotations["CMS_prediction"].notna().sum()),
        "n_kras_conflicts": int(annotations["KRAS_conflict"].fillna(False).sum()),
        "n_annotation_conflicts": int(
            annotations[[c for c in annotations if c.endswith("_conflict")]].fillna(False).to_numpy().sum()
        ),
        "screen_rows_selected": int(len(screen)),
        "screen_genes": int(screen["gene"].nunique()),
        "panel_genes_present": sorted(set(panel["gene"])),
        "panel_genes_absent": sorted(set(PANEL) - set(panel["gene"])),
    }
    (args.outdir / "organoid_reanalysis_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

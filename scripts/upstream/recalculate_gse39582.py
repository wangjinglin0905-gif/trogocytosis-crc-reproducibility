#!/usr/bin/env python
"""Recalculate the frozen 8-gene module in the official GSE39582 tumour set.

The Bioconductor curatedCRCData ExpressionSet contains 566 tumour samples and
gene-level fRMA expression.  The GEO author metadata exposes relapse-free
survival (RFS), but not overall survival (OS).  This script therefore models
RFS and audits the unavailable OS endpoint rather than relabelling DFS as OS.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rdata
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import proportional_hazard_test
from scipy import stats


FROZEN_GENES = ["CD4", "PTPRC", "CTLA4", "PDCD1", "HAVCR2", "VSIR", "LAG3", "CD38"]
ALIAS = {"VSIR": "C10orf54"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def extract_field(text: object, field: str) -> str | float:
    match = re.search(rf"(?:^|///)(?:characteristics_ch1(?:\.\d+)?: )?{re.escape(field)}: ([^/]+)", str(text))
    return match.group(1).strip() if match else np.nan


def zscore(series: pd.Series) -> pd.Series:
    return (series - series.mean()) / series.std(ddof=1)


def prepare_model(frame: pd.DataFrame, covariates: list[str]) -> pd.DataFrame:
    required = ["rfs_months", "rfs_event", *covariates]
    data = frame[required].replace([np.inf, -np.inf], np.nan).dropna().copy()
    data = data.loc[data["rfs_months"] > 0].copy()
    return data


def fit_cox(frame: pd.DataFrame, covariates: list[str], label: str) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    data = prepare_model(frame, covariates)
    cph = CoxPHFitter()
    cph.fit(data, duration_col="rfs_months", event_col="rfs_event", show_progress=False)
    summary = cph.summary.reset_index(names="term")
    row = cph.summary.loc["module_z"]
    ph = proportional_hazard_test(cph, data, time_transform="rank").summary.reset_index(names="term")
    global_chi = float(ph["test_statistic"].sum())
    result = {
        "model": label,
        "n": int(len(data)),
        "events": int(data["rfs_event"].sum()),
        "hr_per_sd": float(math.exp(row["coef"])),
        "ci95_low": float(math.exp(row["coef lower 95%"])),
        "ci95_high": float(math.exp(row["coef upper 95%"])),
        "p": float(row["p"]),
        "concordance": float(cph.concordance_index_),
        "module_ph_p_rank": float(ph.loc[ph["term"].eq("module_z"), "p"].iloc[0]),
        "global_ph_component_sum_p": float(stats.chi2.sf(global_chi, df=len(ph))),
    }
    summary.insert(0, "model", label)
    ph.insert(0, "model", label)
    return result, summary, ph


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-rda", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    eset = rdata.read_rda(args.input_rda)["GSE39582_eset"]
    assay = eset.assayData["exprs"]
    genes = pd.Index(assay.coords[assay.dims[0]].values.astype(str), name="feature")
    samples = pd.Index(assay.coords[assay.dims[1]].values.astype(str), name="sample_id")
    matrix = np.asarray(assay.values, dtype=float)

    mapping = []
    selected = {}
    for requested in FROZEN_GENES:
        feature = requested if requested in genes else ALIAS.get(requested)
        if feature is None or feature not in genes:
            raise RuntimeError(f"Frozen module gene unavailable: {requested}")
        positions = np.flatnonzero(genes == feature)
        selected[requested] = matrix[positions].mean(axis=0)
        mapping.append({
            "module_symbol": requested,
            "assay_feature": feature,
            "alias_used": requested != feature,
            "n_collapsed_features": int(len(positions)),
        })
    expr = pd.DataFrame(selected, index=samples)
    gene_z = expr.apply(zscore, axis=0)
    module_raw = gene_z.mean(axis=1, skipna=False)
    module_z = zscore(module_raw)

    pheno = eset.phenoData.data.copy()
    pheno.index = pheno.index.astype(str)
    metadata = pheno["uncurated_author_metadata"].astype(str)
    clinical = pd.DataFrame(index=samples)
    clinical["title"] = metadata.map(lambda x: extract_field(x, "title")).reindex(samples)
    clinical["rfs_event"] = pd.to_numeric(metadata.map(lambda x: extract_field(x, "rfs.event")), errors="coerce").reindex(samples)
    clinical["rfs_months"] = pd.to_numeric(metadata.map(lambda x: extract_field(x, "rfs.delay")), errors="coerce").reindex(samples)
    clinical["age"] = pd.to_numeric(metadata.map(lambda x: extract_field(x, "age.at.diagnosis")), errors="coerce").reindex(samples)
    clinical["stage_num"] = pd.to_numeric(metadata.map(lambda x: extract_field(x, "tnm.stage")), errors="coerce").reindex(samples)
    clinical["mmr"] = metadata.map(lambda x: extract_field(x, "mmr.status")).reindex(samples)
    clinical["chemotherapy"] = metadata.map(lambda x: extract_field(x, "chemotherapy.adjuvant")).reindex(samples)
    clinical["module_z"] = module_z
    clinical["mmrd"] = clinical["mmr"].astype(str).str.lower().eq("dmmr").astype(float)
    clinical["chemotherapy_yes"] = clinical["chemotherapy"].astype(str).str.upper().eq("Y").astype(float)

    # Stage is modelled categorically to avoid an untested linear stage effect.
    stage_dummies = pd.get_dummies(clinical["stage_num"].astype("Int64"), prefix="stage", drop_first=True, dtype=float)
    clinical = clinical.join(stage_dummies)
    stage_terms = list(stage_dummies.columns)

    specifications = [
        ("RFS_univariable", ["module_z"]),
        ("RFS_age_stage_MMR_chemo_adjusted", ["module_z", "age", "mmrd", "chemotherapy_yes", *stage_terms]),
    ]
    results, full, ph_rows = [], [], []
    for label, covariates in specifications:
        result, model_summary, ph = fit_cox(clinical, covariates, label)
        results.append(result)
        full.append(model_summary)
        ph_rows.append(ph)

    # Harmonised curatedCRCData DFS is retained only as an endpoint-audit table.
    # It is not substituted for OS and is not used for the primary inference.
    endpoint_audit = {
        "series_total_samples_reported_by_GEO": 585,
        "tumour_expression_samples_in_curatedCRCData": int(len(samples)),
        "days_to_death_nonmissing": int(pd.Series(pheno["days_to_death"]).notna().sum()),
        "vital_status_nonmissing": int(pd.Series(pheno["vital_status"]).notna().sum()),
        "author_rfs_event_nonmissing": int(clinical["rfs_event"].notna().sum()),
        "author_rfs_delay_nonmissing": int(clinical["rfs_months"].notna().sum()),
        "author_rfs_positive_time": int((clinical["rfs_months"] > 0).sum()),
        "author_rfs_zero_time_excluded": int((clinical["rfs_months"] == 0).sum()),
        "curated_dfs_status_nonmissing": int(pd.Series(pheno["dfs_status"]).notna().sum()),
        "curated_days_to_recurrence_or_death_nonmissing": int(pd.Series(pheno["days_to_recurrence_or_death"]).notna().sum()),
        "os_conclusion": "OS cannot be independently reconstructed from this versioned official object because days_to_death and vital_status are entirely missing.",
    }

    km_data = prepare_model(clinical, ["module_z"])
    km_data["module_group"] = np.where(km_data["module_z"] >= km_data["module_z"].median(), "High", "Low")
    curves = []
    for group, part in km_data.groupby("module_group"):
        kmf = KaplanMeierFitter(label=group).fit(part["rfs_months"], part["rfs_event"])
        curve = kmf.survival_function_.reset_index()
        curve.columns = ["time_months", "survival"]
        ci = kmf.confidence_interval_survival_function_.reset_index(drop=True)
        curve["ci95_low"] = ci.iloc[:, 0].to_numpy()
        curve["ci95_high"] = ci.iloc[:, 1].to_numpy()
        curve["group"] = group
        curve["n_group"] = int(len(part))
        curve["events_group"] = int(part["rfs_event"].sum())
        curves.append(curve)

    expr.to_csv(args.output_dir / "gse39582_frozen_module_expression.tsv", sep="\t")
    gene_z.to_csv(args.output_dir / "gse39582_frozen_module_gene_z.tsv", sep="\t")
    clinical.to_csv(args.output_dir / "gse39582_module_clinical.tsv", sep="\t")
    pd.DataFrame(mapping).to_csv(args.output_dir / "gse39582_gene_mapping.tsv", sep="\t", index=False)
    pd.DataFrame(results).to_csv(args.output_dir / "gse39582_cox_results.tsv", sep="\t", index=False)
    pd.concat(full, ignore_index=True).to_csv(args.output_dir / "gse39582_cox_full.tsv", sep="\t", index=False)
    pd.concat(ph_rows, ignore_index=True).to_csv(args.output_dir / "gse39582_schoenfeld_tests.tsv", sep="\t", index=False)
    pd.concat(curves, ignore_index=True).to_csv(args.output_dir / "gse39582_km_curve.tsv", sep="\t", index=False)
    (args.output_dir / "gse39582_endpoint_audit.json").write_text(json.dumps(endpoint_audit, ensure_ascii=False, indent=2), encoding="utf-8")

    manifest = {
        "source": "Bioconductor curatedCRCData GSE39582_eset",
        "source_url": "https://bioconductor.org/packages/curatedCRCData/",
        "input_file": str(args.input_rda),
        "input_bytes": args.input_rda.stat().st_size,
        "input_sha256": sha256(args.input_rda),
        "expression_preprocessing": "gene-level fRMA values supplied by curatedCRCData",
        "frozen_module": FROZEN_GENES,
        "score_definition": "mean cohort-wise z-score across eight genes, then re-standardised to 1 cohort SD",
        "alias_policy": ALIAS,
        "primary_endpoint": "author-provided RFS event/delay in months; zero follow-up excluded",
        "adjustment": "age, categorical TNM stage, MMR status, adjuvant chemotherapy",
        "os_policy": "not analysed; official days_to_death and vital_status fields are empty",
        "python": sys.version,
        "endpoint_audit": endpoint_audit,
    }
    (args.output_dir / "gse39582_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(pd.DataFrame(results).to_string(index=False))
    print(json.dumps(endpoint_audit, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

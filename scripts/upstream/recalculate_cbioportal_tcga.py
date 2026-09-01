#!/usr/bin/env python
"""Recalculate the frozen 8-gene module in the 592-sample cBioPortal cohort.

The 592 denominator is the RNA-seq sample-list denominator. Survival models use
only one expression sample per patient with complete, positive OS follow-up.
This cohort is analysed independently from the UCSC Xena primary-tumour cohort.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import proportional_hazard_test


BASE_URL = "https://www.cbioportal.org/api"
STUDY_ID = "coadread_tcga_pan_can_atlas_2018"
PROFILE_ID = f"{STUDY_ID}_rna_seq_v2_mrna"
SAMPLE_LIST_ID = PROFILE_ID
FROZEN_GENES = ["CD4", "PTPRC", "CTLA4", "PDCD1", "HAVCR2", "VSIR", "LAG3", "CD38"]


def get_json(session: requests.Session, endpoint: str, params: dict | None = None):
    url = f"{BASE_URL}/{endpoint.lstrip('/')}"
    for attempt in range(5):
        response = session.get(url, params=params, timeout=180)
        if response.ok:
            return response.json(), response.url
        if response.status_code in {429, 500, 502, 503, 504}:
            time.sleep(2**attempt)
            continue
        response.raise_for_status()
    response.raise_for_status()


def save_json(path: Path, payload) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def tidy_stage(value: str) -> float:
    if not isinstance(value, str):
        return np.nan
    cleaned = value.upper().replace("STAGE", "").strip()
    if cleaned.startswith("IV"):
        return 4.0
    if cleaned.startswith("III"):
        return 3.0
    if cleaned.startswith("II"):
        return 2.0
    if cleaned.startswith("I"):
        return 1.0
    return np.nan


def parse_event(value: str) -> float:
    text = str(value).upper()
    if "DECEASED" in text or text.startswith("1"):
        return 1.0
    if "LIVING" in text or text.startswith("0"):
        return 0.0
    return np.nan


def fit_cox(frame: pd.DataFrame, covariates: list[str], label: str) -> tuple[dict, pd.DataFrame]:
    cols = ["os_months", "os_event", *covariates]
    model_df = frame[cols].replace([np.inf, -np.inf], np.nan).dropna().copy()
    model_df = model_df.loc[model_df["os_months"] > 0]
    cph = CoxPHFitter()
    cph.fit(model_df, duration_col="os_months", event_col="os_event", show_progress=False)
    row = cph.summary.loc["module_z"]
    ph = proportional_hazard_test(cph, model_df, time_transform="rank").summary
    result = {
        "model": label,
        "n": int(len(model_df)),
        "events": int(model_df["os_event"].sum()),
        "hr_per_sd": float(math.exp(row["coef"])),
        "ci95_low": float(math.exp(row["coef lower 95%"])),
        "ci95_high": float(math.exp(row["coef upper 95%"])),
        "p": float(row["p"]),
        "concordance": float(cph.concordance_index_),
        "module_ph_p_rank": float(ph.loc["module_z", "p"]),
        "global_ph_min_p_rank": float(ph["p"].min()),
    }
    return result, cph.summary.reset_index(names="term")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--raw-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.raw_dir.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    session.headers.update({"User-Agent": "trogocytosis-crc-recalculation/1.0"})

    sample_list, sample_url = get_json(session, f"sample-lists/{SAMPLE_LIST_ID}")
    if len(sample_list["sampleIds"]) != 592:
        raise RuntimeError(f"Expected 592 samples, received {len(sample_list['sampleIds'])}")
    sample_path = args.raw_dir / "sample_list_592.json"
    save_json(sample_path, sample_list)

    clinical_long, clinical_url = get_json(
        session,
        f"studies/{STUDY_ID}/clinical-data",
        {"clinicalDataType": "PATIENT", "projection": "SUMMARY", "pageSize": 100000},
    )
    clinical_path = args.raw_dir / "patient_clinical_long.json"
    save_json(clinical_path, clinical_long)
    clinical = pd.DataFrame(clinical_long).pivot_table(
        index="patientId", columns="clinicalAttributeId", values="value", aggfunc="first"
    )

    expression = {}
    gene_records = []
    source_urls = {"sample_list": sample_url, "clinical": clinical_url}
    for symbol in FROZEN_GENES:
        gene, gene_url = get_json(session, f"genes/{symbol}")
        entrez = int(gene["entrezGeneId"])
        rows, expression_url = get_json(
            session,
            f"molecular-profiles/{PROFILE_ID}/molecular-data",
            {"sampleListId": SAMPLE_LIST_ID, "entrezGeneId": entrez, "projection": "SUMMARY"},
        )
        save_json(args.raw_dir / f"expression_{symbol}.json", rows)
        source_urls[f"gene_{symbol}"] = gene_url
        source_urls[f"expression_{symbol}"] = expression_url
        gene_records.append({"symbol": symbol, "entrez_gene_id": entrez, "n_records": len(rows)})
        expression[symbol] = pd.Series(
            {row["sampleId"]: pd.to_numeric(row.get("value"), errors="coerce") for row in rows},
            name=symbol,
            dtype=float,
        )

    expr = pd.DataFrame(expression).reindex(sample_list["sampleIds"])
    expr.index.name = "sample_id"
    # cBioPortal values are continuous RSEM abundance. Log2(x+1) limits the
    # leverage of extreme abundance before cohort-wise gene standardisation.
    log_expr = np.log2(expr + 1.0)
    gene_z = (log_expr - log_expr.mean(axis=0)) / log_expr.std(axis=0, ddof=1)
    module = gene_z.mean(axis=1, skipna=False)
    module_z = (module - module.mean()) / module.std(ddof=1)

    cohort = expr.copy()
    cohort.columns = [f"expr_{c}" for c in cohort.columns]
    cohort["module_raw_zmean"] = module
    cohort["module_z"] = module_z
    cohort["patient_id"] = cohort.index.str.slice(0, 12)
    cohort["sample_type_code"] = cohort.index.str.slice(13, 15)
    cohort["os_months"] = pd.to_numeric(cohort["patient_id"].map(clinical.get("OS_MONTHS")), errors="coerce")
    cohort["os_event"] = cohort["patient_id"].map(clinical.get("OS_STATUS")).map(parse_event)
    cohort["age"] = pd.to_numeric(cohort["patient_id"].map(clinical.get("AGE")), errors="coerce")
    cohort["stage_num"] = cohort["patient_id"].map(clinical.get("AJCC_PATHOLOGIC_TUMOR_STAGE")).map(tidy_stage)
    # The profile contains one RNA-seq sample per patient here; retain the first
    # deterministically if a future portal revision introduces duplicates.
    cohort = cohort.sort_index().drop_duplicates("patient_id", keep="first")

    univ, univ_full = fit_cox(cohort, ["module_z"], "OS_univariable")
    adj_age, adj_age_full = fit_cox(cohort, ["module_z", "age"], "OS_age_adjusted")
    adj_stage, adj_stage_full = fit_cox(
        cohort, ["module_z", "age", "stage_num"], "OS_age_stage_adjusted"
    )
    cox_results = pd.DataFrame([univ, adj_age, adj_stage])

    def standardize(score: pd.Series) -> pd.Series:
        return (score - score.mean()) / score.std(ddof=1)

    raw_gene_z = (expr - expr.mean(axis=0)) / expr.std(axis=0, ddof=1)
    score_variants = {
        "primary_log2_then_gene_z_mean": module_z,
        "raw_RSEM_mean": standardize(expr.mean(axis=1)),
        "log2_RSEM_mean": standardize(log_expr.mean(axis=1)),
        "raw_RSEM_gene_z_mean": standardize(raw_gene_z.mean(axis=1)),
    }
    sensitivity_rows = []
    for variant, score in score_variants.items():
        sensitivity_frame = cohort.copy()
        sensitivity_frame["module_z"] = score.reindex(sensitivity_frame.index).to_numpy()
        result, _ = fit_cox(sensitivity_frame, ["module_z"], f"OS_{variant}")
        result["score_variant"] = variant
        result["effect_scale"] = "per 1 cohort SD"
        sensitivity_rows.append(result)
    sensitivity_df = pd.DataFrame(sensitivity_rows)

    km_df = cohort[["os_months", "os_event", "module_z"]].dropna()
    km_df = km_df.loc[km_df["os_months"] > 0].copy()
    km_df["module_group"] = np.where(km_df["module_z"] >= km_df["module_z"].median(), "High", "Low")
    km_rows = []
    for group, part in km_df.groupby("module_group"):
        kmf = KaplanMeierFitter(label=group).fit(part["os_months"], part["os_event"])
        curve = kmf.survival_function_.reset_index()
        curve.columns = ["time_months", "survival"]
        ci = kmf.confidence_interval_survival_function_.reset_index(drop=True)
        curve["ci95_low"] = ci.iloc[:, 0].to_numpy()
        curve["ci95_high"] = ci.iloc[:, 1].to_numpy()
        curve["group"] = group
        km_rows.append(curve)

    expr.to_csv(args.output_dir / "tcga_cbioportal_expression_592.tsv", sep="\t")
    gene_z.to_csv(args.output_dir / "tcga_cbioportal_gene_z_592.tsv", sep="\t")
    cohort.to_csv(args.output_dir / "tcga_cbioportal_module_clinical.tsv", sep="\t", index=True)
    pd.DataFrame(gene_records).to_csv(args.output_dir / "frozen_module_gene_mapping.tsv", sep="\t", index=False)
    cox_results.to_csv(args.output_dir / "tcga_cbioportal_cox_results.tsv", sep="\t", index=False)
    sensitivity_df.to_csv(
        args.output_dir / "tcga_cbioportal_score_sensitivity.tsv", sep="\t", index=False
    )
    pd.concat([univ_full, adj_age_full, adj_stage_full], keys=["univ", "age", "age_stage"]).to_csv(
        args.output_dir / "tcga_cbioportal_cox_full.tsv", sep="\t", index=False
    )
    pd.concat(km_rows, ignore_index=True).to_csv(
        args.output_dir / "tcga_cbioportal_km_curve.tsv", sep="\t", index=False
    )

    sample_flow = {
        "rna_sample_list_n": len(sample_list["sampleIds"]),
        "unique_patients_after_mapping_n": int(cohort["patient_id"].nunique()),
        "complete_8_gene_expression_n": int(cohort["module_z"].notna().sum()),
        "os_complete_positive_time_n": int(
            cohort[["module_z", "os_months", "os_event"]].dropna().query("os_months > 0").shape[0]
        ),
        "os_events_n": int(
            cohort[["module_z", "os_months", "os_event"]].dropna().query("os_months > 0")["os_event"].sum()
        ),
        "primary_tumour_code_01_n": int((cohort["sample_type_code"] == "01").sum()),
    }
    (args.output_dir / "tcga_cbioportal_sample_flow.json").write_text(
        json.dumps(sample_flow, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    manifest = {
        "study_id": STUDY_ID,
        "molecular_profile_id": PROFILE_ID,
        "sample_list_id": SAMPLE_LIST_ID,
        "frozen_genes": FROZEN_GENES,
        "score_definition": "mean cohort-wise gene z-scores after log2(RSEM+1), then module re-standardised",
        "urls": source_urls,
        "raw_files": {
            path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in sorted(args.raw_dir.glob("*.json"))
        },
        "sample_flow": sample_flow,
        "python": sys.version,
    }
    (args.output_dir / "tcga_cbioportal_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "sample_flow": sample_flow,
                "cox": cox_results.to_dict("records"),
                "score_sensitivity": sensitivity_df.to_dict("records"),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

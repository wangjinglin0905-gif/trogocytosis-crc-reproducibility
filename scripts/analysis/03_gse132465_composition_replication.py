#!/usr/bin/env python3
"""Prespecified GSE132465 composition-tracking replication.

The analysis unit is the patient. Raw UMI counts are aggregated into tumour-wide
and author-annotated epithelial pseudobulks. The frozen eight-gene module is
scored as the mean of gene-wise z scores of log2(CPM + 0.25); invariant genes
contribute zero and are explicitly recorded. Cell-level library sizes are reused
from an earlier extraction only after SHA-256 identity of the GEO raw matrix and
annotation has been verified.
"""

from __future__ import annotations

import gzip
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


GENES = ["CD4", "PTPRC", "CTLA4", "PDCD1", "HAVCR2", "VSIR", "LAG3", "CD38"]
GENE_ALIASES = {"C10orf54": "VSIR"}
SEED = 20260901
N_RESAMPLES = 10_000
PSEUDOCOUNT = 0.25


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def extract_target_rows(
    path: Path, genes: list[str]
) -> tuple[list[str], dict[str, np.ndarray], dict[str, str]]:
    wanted = set(genes) | set(GENE_ALIASES)
    found: dict[str, np.ndarray] = {}
    source_symbol: dict[str, str] = {}
    with gzip.open(path, "rb") as handle:
        header = handle.readline().rstrip(b"\r\n").decode("utf-8").split("\t")
        cell_ids = header[1:]
        for line_number, raw_line in enumerate(handle, start=2):
            tab = raw_line.find(b"\t")
            if tab < 0:
                continue
            gene = raw_line[:tab].decode("utf-8")
            if gene not in wanted:
                continue
            target = GENE_ALIASES.get(gene, gene)
            values = np.fromstring(raw_line[tab + 1 :], sep="\t", dtype=np.int64)
            if values.size != len(cell_ids):
                raise ValueError(
                    f"{gene} at line {line_number}: {values.size} values, expected {len(cell_ids)}"
                )
            if target in found:
                # Prefer the current frozen symbol if both current and legacy
                # nomenclature are present; never sum duplicate aliases.
                if gene == target and source_symbol[target] != target:
                    found[target] = values
                    source_symbol[target] = gene
                continue
            found[target] = values
            source_symbol[target] = gene
    missing = sorted(set(genes).difference(found))
    if missing:
        raise ValueError(f"Missing frozen-module genes: {missing}")
    return cell_ids, found, source_symbol


def aggregate(values: np.ndarray, group_codes: np.ndarray, n_groups: int) -> np.ndarray:
    keep = group_codes >= 0
    return np.bincount(
        group_codes[keep], weights=values[keep], minlength=n_groups
    ).astype(float)


def score_module(
    count_matrix: np.ndarray, library_size: np.ndarray, gene_names: list[str]
) -> tuple[np.ndarray, pd.DataFrame, np.ndarray]:
    if np.any(library_size <= 0):
        raise ValueError("Non-positive pseudobulk library size")
    log_cpm = np.log2(count_matrix / library_size[np.newaxis, :] * 1e6 + PSEUDOCOUNT)
    z = np.zeros_like(log_cpm, dtype=float)
    audit_rows = []
    for i, gene in enumerate(gene_names):
        sd = float(np.std(log_cpm[i, :], ddof=0))
        invariant = (not np.isfinite(sd)) or sd == 0.0
        if not invariant:
            z[i, :] = (log_cpm[i, :] - np.mean(log_cpm[i, :])) / sd
        audit_rows.append(
            {
                "gene": gene,
                "mean_log2_cpm": float(np.mean(log_cpm[i, :])),
                "sd_log2_cpm": sd,
                "patients_detected": int(np.sum(count_matrix[i, :] > 0)),
                "invariant_zero_contribution": bool(invariant),
            }
        )
    return np.mean(z, axis=0), pd.DataFrame(audit_rows), log_cpm


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    if np.unique(x).size < 2 or np.unique(y).size < 2:
        return np.nan
    return float(stats.spearmanr(x, y).statistic)


def resampling_inference(
    whole_score: np.ndarray, epithelial_score: np.ndarray, t_fraction: np.ndarray
) -> tuple[dict[str, float], pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(SEED)
    n = len(t_fraction)
    rho_whole = safe_spearman(whole_score, t_fraction)
    rho_epi = safe_spearman(epithelial_score, t_fraction)
    attenuation = abs(rho_whole) - abs(rho_epi)

    boot_whole = np.empty(N_RESAMPLES, dtype=float)
    boot_epi = np.empty(N_RESAMPLES, dtype=float)
    boot_attenuation = np.empty(N_RESAMPLES, dtype=float)
    for i in range(N_RESAMPLES):
        idx = rng.integers(0, n, size=n)
        boot_whole[i] = safe_spearman(whole_score[idx], t_fraction[idx])
        boot_epi[i] = safe_spearman(epithelial_score[idx], t_fraction[idx])
        boot_attenuation[i] = abs(boot_whole[i]) - abs(boot_epi[i])

    perm_whole = np.empty(N_RESAMPLES, dtype=float)
    perm_epi = np.empty(N_RESAMPLES, dtype=float)
    for i in range(N_RESAMPLES):
        shuffled = rng.permutation(t_fraction)
        perm_whole[i] = safe_spearman(whole_score, shuffled)
        perm_epi[i] = safe_spearman(epithelial_score, shuffled)

    finite_bw = boot_whole[np.isfinite(boot_whole)]
    finite_be = boot_epi[np.isfinite(boot_epi)]
    finite_ba = boot_attenuation[np.isfinite(boot_attenuation)]
    finite_pw = perm_whole[np.isfinite(perm_whole)]
    finite_pe = perm_epi[np.isfinite(perm_epi)]

    scipy_whole = stats.spearmanr(whole_score, t_fraction)
    scipy_epi = stats.spearmanr(epithelial_score, t_fraction)
    summary = {
        "n_patients": n,
        "rho_whole": rho_whole,
        "asymptotic_p_whole": float(scipy_whole.pvalue),
        "permutation_p_whole": float((1 + np.sum(np.abs(finite_pw) >= abs(rho_whole))) / (1 + len(finite_pw))),
        "bootstrap_ci_low_whole": float(np.quantile(finite_bw, 0.025)),
        "bootstrap_ci_high_whole": float(np.quantile(finite_bw, 0.975)),
        "rho_epithelial": rho_epi,
        "asymptotic_p_epithelial": float(scipy_epi.pvalue),
        "permutation_p_epithelial": float((1 + np.sum(np.abs(finite_pe) >= abs(rho_epi))) / (1 + len(finite_pe))),
        "bootstrap_ci_low_epithelial": float(np.quantile(finite_be, 0.025)),
        "bootstrap_ci_high_epithelial": float(np.quantile(finite_be, 0.975)),
        "attenuation_abs_rho": attenuation,
        "attenuation_bootstrap_ci_low": float(np.quantile(finite_ba, 0.025)),
        "attenuation_bootstrap_ci_high": float(np.quantile(finite_ba, 0.975)),
        "strong_replication": bool(
            rho_whole > 0
            and ((1 + np.sum(np.abs(finite_pw) >= abs(rho_whole))) / (1 + len(finite_pw))) < 0.05
            and attenuation > 0.20
            and np.quantile(finite_ba, 0.025) > 0
        ),
        "compatible_replication": bool(
            rho_whole > 0
            and ((1 + np.sum(np.abs(finite_pw) >= abs(rho_whole))) / (1 + len(finite_pw))) < 0.05
            and attenuation > 0.20
        ),
    }
    bootstrap = pd.DataFrame(
        {
            "replicate": np.arange(1, N_RESAMPLES + 1),
            "rho_whole": boot_whole,
            "rho_epithelial": boot_epi,
            "attenuation_abs_rho": boot_attenuation,
        }
    )
    permutation = pd.DataFrame(
        {
            "replicate": np.arange(1, N_RESAMPLES + 1),
            "rho_whole": perm_whole,
            "rho_epithelial": perm_epi,
        }
    )
    return summary, bootstrap, permutation


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Patient-level GSE132465 composition-tracking replication."
    )
    parser.add_argument("--umi", required=True, type=Path,
                        help="GEO raw UMI matrix (.txt.gz).")
    parser.add_argument("--annotation", required=True, type=Path,
                        help="GEO author cell annotations (.txt.gz).")
    parser.add_argument("--library-sizes", required=True, type=Path,
                        help="Per-cell total UMI table generated from the same matrix.")
    parser.add_argument("--library-size-manifest", required=True, type=Path,
                        help="Manifest containing source and output SHA-256 values.")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--resamples", type=int, default=10_000)
    args = parser.parse_args()

    global N_RESAMPLES
    N_RESAMPLES = args.resamples
    umi = args.umi.resolve()
    annotation_path = args.annotation.resolve()
    library_path = args.library_sizes.resolve()
    manifest_path = args.library_size_manifest.resolve()
    out = args.outdir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    raw_hash = sha256(umi)
    annotation_hash = sha256(annotation_path)
    prior_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if raw_hash.lower() != prior_manifest["source_umi_sha256"].lower():
        raise ValueError("Library-size file does not derive from the current raw UMI matrix")
    if annotation_hash.lower() != prior_manifest["source_annotation_sha256"].lower():
        raise ValueError("Library-size file does not derive from the current annotation")
    if sha256(library_path).lower() != prior_manifest["library_size_output_sha256"].lower():
        raise ValueError("Library-size file hash mismatch")

    annotation = pd.read_csv(annotation_path, sep="\t")
    cell_ids, rows, source_symbol = extract_target_rows(umi, GENES)
    library = pd.read_csv(library_path, sep="\t")
    if list(annotation["Index"].astype(str)) != cell_ids:
        raise ValueError("Raw UMI header and annotation cell order differ")
    if list(library["cell_id"].astype(str)) != cell_ids:
        raise ValueError("Library-size cell order and raw UMI header differ")

    tumour = annotation["Class"].astype(str).eq("Tumor").to_numpy()
    patients = sorted(annotation.loc[tumour, "Patient"].astype(str).unique())
    if len(patients) < 10:
        raise ValueError("Fewer than ten independent tumour patients")
    patient_map = {patient: i for i, patient in enumerate(patients)}
    patient_codes = np.full(len(annotation), -1, dtype=int)
    patient_codes[tumour] = annotation.loc[tumour, "Patient"].astype(str).map(patient_map).to_numpy()
    epithelial = tumour & annotation["Cell_type"].astype(str).eq("Epithelial cells").to_numpy()
    t_cells = tumour & annotation["Cell_type"].astype(str).eq("T cells").to_numpy()

    whole_codes = np.where(tumour, patient_codes, -1)
    epithelial_codes = np.where(epithelial, patient_codes, -1)
    whole_library = aggregate(library["total_umi"].to_numpy(dtype=float), whole_codes, len(patients))
    epithelial_library = aggregate(library["total_umi"].to_numpy(dtype=float), epithelial_codes, len(patients))
    whole_counts = np.vstack(
        [aggregate(rows[gene], whole_codes, len(patients)) for gene in GENES]
    )
    epithelial_counts = np.vstack(
        [aggregate(rows[gene], epithelial_codes, len(patients)) for gene in GENES]
    )

    n_tumour = aggregate(np.ones(len(annotation)), whole_codes, len(patients))
    n_epithelial = aggregate(np.ones(len(annotation)), epithelial_codes, len(patients))
    t_codes = np.where(t_cells, patient_codes, -1)
    n_t = aggregate(np.ones(len(annotation)), t_codes, len(patients))
    t_fraction = n_t / n_tumour

    whole_score, whole_audit, whole_log_cpm = score_module(whole_counts, whole_library, GENES)
    epithelial_score, epithelial_audit, epithelial_log_cpm = score_module(
        epithelial_counts, epithelial_library, GENES
    )
    whole_audit.insert(0, "compartment", "whole_tumour")
    epithelial_audit.insert(0, "compartment", "epithelial")
    gene_audit = pd.concat([whole_audit, epithelial_audit], ignore_index=True)

    patient_table = pd.DataFrame(
        {
            "patient": patients,
            "n_tumour_cells": n_tumour.astype(int),
            "n_epithelial_cells": n_epithelial.astype(int),
            "n_t_cells": n_t.astype(int),
            "t_cell_fraction": t_fraction,
            "whole_tumour_module_score": whole_score,
            "epithelial_module_score": epithelial_score,
            "whole_tumour_library_size": whole_library.astype(np.int64),
            "epithelial_library_size": epithelial_library.astype(np.int64),
        }
    )
    for i, gene in enumerate(GENES):
        patient_table[f"whole_{gene}_log2CPM"] = whole_log_cpm[i, :]
        patient_table[f"epithelial_{gene}_log2CPM"] = epithelial_log_cpm[i, :]

    summary, bootstrap, permutation = resampling_inference(
        whole_score, epithelial_score, t_fraction
    )

    # Descriptive lineage-stratified scores. Gene-wise standardisation is done
    # jointly across all patient-lineage pseudobulks so cell-type values remain
    # comparable; no cell-level inference is made.
    lineage_blocks = []
    for cell_type in sorted(annotation.loc[tumour, "Cell_type"].astype(str).unique()):
        in_lineage = tumour & annotation["Cell_type"].astype(str).eq(cell_type).to_numpy()
        lineage_codes = np.where(in_lineage, patient_codes, -1)
        counts = np.vstack(
            [aggregate(rows[gene], lineage_codes, len(patients)) for gene in GENES]
        )
        lib = aggregate(library["total_umi"].to_numpy(dtype=float), lineage_codes, len(patients))
        cells = aggregate(np.ones(len(annotation)), lineage_codes, len(patients))
        keep = (cells > 0) & (lib > 0)
        if np.sum(keep) < 3:
            continue
        log_cpm = np.log2(counts[:, keep] / lib[keep][np.newaxis, :] * 1e6 + PSEUDOCOUNT)
        lineage_blocks.append(
            {
                "cell_type": cell_type,
                "patients": np.array(patients)[keep],
                "cells": cells[keep],
                "log_cpm": log_cpm,
            }
        )
    stacked = np.hstack([block["log_cpm"] for block in lineage_blocks])
    stacked_z = np.zeros_like(stacked)
    variable = np.zeros(len(GENES), dtype=bool)
    for i in range(len(GENES)):
        sd = float(np.std(stacked[i, :], ddof=0))
        if np.isfinite(sd) and sd > 0:
            stacked_z[i, :] = (stacked[i, :] - np.mean(stacked[i, :])) / sd
            variable[i] = True
    stacked_score = np.mean(stacked_z, axis=0)
    lineage_rows = []
    offset = 0
    for block in lineage_blocks:
        width = block["log_cpm"].shape[1]
        for patient, n_cells, this_score in zip(
            block["patients"], block["cells"], stacked_score[offset : offset + width]
        ):
            lineage_rows.append(
                {
                    "patient": patient,
                    "cell_type": block["cell_type"],
                    "n_cells": int(n_cells),
                    "module_score_across_lineages": float(this_score),
                    "n_variable_module_genes": int(variable.sum()),
                }
            )
        offset += width
    lineage_table = pd.DataFrame(lineage_rows)

    patient_table.to_csv(out / "gse132465_patient_scores.tsv", sep="\t", index=False)
    gene_audit.to_csv(out / "gse132465_module_gene_audit.tsv", sep="\t", index=False)
    bootstrap.to_csv(out / "gse132465_bootstrap.tsv.gz", sep="\t", index=False, compression="gzip")
    permutation.to_csv(out / "gse132465_permutation.tsv.gz", sep="\t", index=False, compression="gzip")
    lineage_table.to_csv(out / "gse132465_lineage_scores.tsv", sep="\t", index=False)
    pd.DataFrame([summary]).to_csv(out / "gse132465_composition_summary.tsv", sep="\t", index=False)

    extracted = pd.DataFrame({"cell_id": cell_ids})
    for gene in GENES:
        extracted[gene] = rows[gene]
    extracted.to_csv(out / "gse132465_frozen_module_cell_counts.tsv.gz", sep="\t", index=False, compression="gzip")

    audit = {
        "seed": SEED,
        "n_resamples": N_RESAMPLES,
        "raw_umi": umi.name,
        "raw_umi_sha256": raw_hash,
        "annotation": annotation_path.name,
        "annotation_sha256": annotation_hash,
        "cell_library_sizes": library_path.name,
        "cell_library_sizes_sha256": sha256(library_path),
        "n_all_cells": int(len(annotation)),
        "n_tumour_cells": int(tumour.sum()),
        "n_epithelial_tumour_cells": int(epithelial.sum()),
        "n_t_cells": int(t_cells.sum()),
        "n_tumour_patients": int(len(patients)),
        "cell_types_tumour": annotation.loc[tumour, "Cell_type"].value_counts().to_dict(),
        "frozen_genes": GENES,
        "source_symbol_by_frozen_gene": source_symbol,
        "summary": summary,
        "interpretation": (
            "strong_cross_cohort_composition_replication"
            if summary["strong_replication"]
            else "compatible_cross_cohort_composition_replication"
            if summary["compatible_replication"]
            else "prespecified_cross_cohort_replication_not_met"
        ),
    }
    (out / "gse132465_composition_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

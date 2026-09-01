from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


PANEL = [
    "ANO6", "ATF3", "BCAS1", "C3", "CCR7", "CD109", "CD19", "CD22",
    "CD24", "CD274", "CD38", "CD4", "CD47", "CD80", "CD86", "CDH2",
    "CEACAM5", "CH25H", "CLSTN2", "CTLA4", "CTSE", "EGFR", "ERBB2",
    "FCGR1A", "FCGR2B", "FCGR3A", "HAVCR2", "HLA-DRA", "IL6", "KANK4",
    "LAG3", "MSLN", "PDCD1", "PTPRC", "SCD", "SIGLEC10", "SIRPA",
    "STAT1", "VSIR",
]
EPITHELIAL_BACKGROUND = [
    "CDX2", "AGR2", "CLDN4", "CLDN3", "KRT7", "KRT8", "KRT18", "KRT19",
    "EPCAM", "MUC1", "CEACAM6", "KRT20", "MUC13", "TFF3", "CDH1",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--marker-tpm", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.audit.parent.mkdir(parents=True, exist_ok=True)

    tpm = pd.read_csv(args.marker_tpm, sep=r"\s+", index_col=0)
    tpm.index = tpm.index.astype(str).str.strip('"')
    tpm.columns = tpm.columns.astype(str).str.strip('"')
    metadata = pd.read_csv(args.metadata, sep="\t")
    metadata["CellName"] = metadata["CellName"].astype(str).str.strip('"')
    metadata = metadata.set_index("CellName")
    cells = tpm.columns.intersection(metadata.index)
    tpm = tpm[cells]
    cell = metadata.loc[cells].copy()
    module = [gene for gene in PANEL if gene in tpm.index]
    background = [gene for gene in EPITHELIAL_BACKGROUND if gene in tpm.index]
    if "CEACAM5" not in tpm.index:
        raise ValueError("CEACAM5 is absent from the filtered TPM matrix")
    cell["TRG_module"] = np.log2(tpm.loc[module] + 1).mean(axis=0).to_numpy()
    cell["epi_score"] = np.log2(tpm.loc[background].sum(axis=0) + 1).to_numpy()
    cell["CEACAM5_TPM"] = tpm.loc["CEACAM5"].to_numpy()
    keep = ["Sample", "Tissue", "Global_Cluster", "Sub_Cluster", "TRG_module", "epi_score", "CEACAM5_TPM"]
    cell[keep].to_csv(args.output)
    audit = {
        "cells": int(len(cell)),
        "module_genes": module,
        "epithelial_background_genes": background,
        "output_semantics": "inputs for contamination-filtered candidate analysis; no cell is labelled as a trogocytosis event",
    }
    args.audit.write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

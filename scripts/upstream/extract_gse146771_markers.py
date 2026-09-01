from __future__ import annotations

import argparse
import gzip
from pathlib import Path


PANEL = {
    "ANO6", "ATF3", "BCAS1", "C3", "CCR7", "CD109", "CD19", "CD22",
    "CD24", "CD274", "CD38", "CD4", "CD47", "CD80", "CD86", "CDH2",
    "CEACAM5", "CH25H", "CLSTN2", "CTLA4", "CTSE", "EGFR", "ERBB2",
    "FCGR1A", "FCGR2B", "FCGR3A", "HAVCR2", "HLA-DRA", "IL6", "KANK4",
    "LAG3", "MSLN", "PDCD1", "PTPRC", "SCD", "SIGLEC10", "SIRPA",
    "STAT1", "VSIR",
}
EPITHELIAL_BACKGROUND = {
    "CDX2", "AGR2", "CLDN4", "CLDN3", "KRT7", "KRT8", "KRT18", "KRT19",
    "EPCAM", "MUC1", "CEACAM6", "KRT20", "MUC13", "TFF3", "CDH1",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tpm-gzip", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    keep = PANEL | EPITHELIAL_BACKGROUND | {"CEACAM5"}
    retained = []
    with gzip.open(args.tpm_gzip, "rt", encoding="utf-8", errors="replace") as source, args.output.open(
        "w", encoding="utf-8", newline=""
    ) as target:
        target.write(source.readline())
        for line in source:
            gene = line.split(None, 1)[0].strip().strip('"')
            if gene in keep:
                target.write(line)
                retained.append(gene)
    print(f"retained {len(retained)} genes: {','.join(sorted(retained))}")


if __name__ == "__main__":
    main()

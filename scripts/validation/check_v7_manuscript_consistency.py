from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--manuscript", required=True, type=Path)
    args = parser.parse_args()

    root = args.root
    text = args.manuscript.read_text(encoding="utf-8")
    text_normalized = text.replace("−", "-")
    checks: list[dict[str, object]] = []

    def contains(name: str, expected: str, source: str) -> None:
        expected_normalized = expected.replace("−", "-")
        checks.append({"check": name, "source": source, "expected": expected,
                       "status": "PASS" if expected_normalized in text_normalized else "FAIL"})

    org = pd.read_csv(root / "analysis/organoid/organoid_panel_dependency_summary.csv").set_index("gene")
    for gene in ["SCD", "EGFR"]:
        row = org.loc[gene]
        contains(f"{gene} organoid median", f"{row.median_LFC:.3f}",
                 "organoid_panel_dependency_summary.csv")
        contains(f"{gene} organoid depleted", f"{row.pct_depleted:.1f}%",
                 "organoid_panel_dependency_summary.csv")

    dep = json.loads((root / "analysis/depmap_26Q1/depmap_reanalysis_audit.json").read_text(encoding="utf-8"))
    contains("cross-platform rho", f"ρ={dep['panel_cross_platform_spearman_rho']:.3f}",
             "depmap_reanalysis_audit.json")
    contains("cross-platform P", f"P={dep['panel_cross_platform_spearman_p']:.4f}",
             "depmap_reanalysis_audit.json")

    gse = pd.read_csv(root / "analysis/gse39582_recalculation/gse39582_cox_results.tsv", sep="\t")
    tcga = pd.read_csv(root / "analysis/tcga_cbioportal_592/tcga_cbioportal_cox_results.tsv", sep="\t")
    for table, prefix in [(gse, "GSE39582"), (tcga, "TCGA")]:
        for row in table.itertuples():
            if row.model in {"RFS_univariable", "RFS_age_stage_MMR_chemo_adjusted",
                             "OS_univariable", "OS_age_stage_adjusted"}:
                contains(f"{prefix} {row.model} HR", f"{row.hr_per_sd:.3f}",
                         f"{prefix.lower()}_cox_results.tsv")
                contains(f"{prefix} {row.model} P", f"P={row.p:.3f}",
                         f"{prefix.lower()}_cox_results.tsv")

    flow = json.loads((root / "analysis/gse178341_recalculation/gse178341_sample_flow.json").read_text(encoding="utf-8"))
    contains("GSE178341 patient count", "62 patients", "gse178341_sample_flow.json")
    contains("GSE178341 MMR counts", "34 MMRd and 28 MMRp", "gse178341_sample_flow.json")

    required_boundaries = [
        "do not establish it as a regulator of trogocytosis",
        "Candidates are not classified as trogocytosis events",
        "Overall survival was not analysed because the curated endpoint fields were missing",
        "do not support an SCD–sterol co-dependency programme",
        "62 patients (34 MMRd and 28 MMRp)",
        "CD4, PTPRC, CTLA4, PDCD1, HAVCR2, VSIR, LAG3",
        "q=0.890",
    ]
    for item in required_boundaries:
        contains("inference boundary", item, "terminology contract")

    forbidden_current_claims = [
        "SCD/sterol vulnerability",
        "high-confidence trogocytosis event",
        "CMS4 enrichment of the proxy was reproduced",
        "HR 1.21 (95% CI 0.98–1.49)",
        "GSE39582 OS",
        "legacy 33-gene",
    ]
    for phrase in forbidden_current_claims:
        checks.append({"check": "forbidden legacy claim", "source": "V4 audit",
                       "expected": f"absent: {phrase}",
                       "status": "PASS" if phrase not in text else "FAIL"})

    outdir = root / "qa/manuscript_consistency_v7"
    outdir.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(checks)
    frame.to_csv(outdir / "numeric_and_inference_checks.tsv", sep="\t", index=False)
    summary = {"checks": len(frame), "passed": int((frame.status == "PASS").sum()),
               "failed": int((frame.status == "FAIL").sum()),
               "failed_checks": frame.loc[frame.status == "FAIL", "check"].tolist()}
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0 if summary["failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())

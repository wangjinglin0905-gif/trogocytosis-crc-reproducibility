from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import requests


FILES = {
    "GSE146771_metadata": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146771/suppl/GSE146771_CRC.Leukocyte.10x.Metadata.txt.gz",
    "GSE146771_tpm": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146771/suppl/GSE146771_CRC.Leukocyte.10x.TPM.txt.gz",
    "GSE132465_annotation": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132465/suppl/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz",
    "GSE132465_counts": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132465/suppl/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz",
    "GSE144735_annotation": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144735/suppl/GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz",
    "GSE144735_counts": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144735/suppl/GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz",
    "GSE39582_series_matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE39nnn/GSE39582/matrix/GSE39582_series_matrix.txt.gz",
    "GPL570_annotation": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz",
    "GSE178341_cluster": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178341/suppl/GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz",
    "GSE178341_metatables": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178341/suppl/GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz",
    "GSE178341_expression_h5": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178341/suppl/GSE178341_crc10x_full_c295v4_submit.h5",
    # DepMap inputs are intentionally not mirrored here. The analysis uses the
    # official Public 24Q4 release (DOI 10.25452/figshare.plus.27993248.v1),
    # supplied by the caller and verified by SHA-256 in the run manifest.
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, destination: Path, retries: int = 8) -> dict:
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    session = requests.Session()
    session.headers.update({"User-Agent": "Mozilla/5.0 trogocytosis-reanalysis/1.0"})
    for attempt in range(1, retries + 1):
        offset = partial.stat().st_size if partial.exists() else 0
        headers = {"Range": f"bytes={offset}-"} if offset else {}
        try:
            with session.get(url, headers=headers, stream=True, timeout=(30, 60), allow_redirects=True) as response:
                response.raise_for_status()
                if offset and response.status_code != 206:
                    partial.unlink(missing_ok=True)
                    offset = 0
                mode = "ab" if offset else "wb"
                with partial.open(mode) as handle:
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            handle.write(chunk)
                            handle.flush()
            partial.replace(destination)
            return {
                "url": url,
                "path": str(destination),
                "bytes": destination.stat().st_size,
                "sha256": sha256(destination),
            }
        except Exception as exc:
            if attempt == retries:
                raise
            time.sleep(min(30, 2**attempt))
            print(f"retry {attempt}/{retries} for {destination.name}: {type(exc).__name__}", flush=True)
    raise RuntimeError("unreachable")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--only", nargs="*", choices=sorted(FILES))
    args = parser.parse_args()
    selected = args.only or list(FILES)
    records = []
    for key in selected:
        url = FILES[key]
        filename = url.split("/")[-1].split("?")[0].replace("%20", "_")
        destination = args.outdir / filename
        print(f"downloading {key} -> {destination}", flush=True)
        record = {"key": key, **download(url, destination)}
        records.append(record)
        print(json.dumps(record, ensure_ascii=False), flush=True)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(records, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

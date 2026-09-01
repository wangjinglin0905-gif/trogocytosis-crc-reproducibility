#!/usr/bin/env python3
"""Stream a GEO gene-by-cell UMI matrix and derive per-cell library sizes.

The output is deterministic (gzip mtime=0) and accompanied by SHA-256
provenance. It avoids loading the full GSE132465 matrix into memory.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path

import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--umi", required=True, type=Path)
    parser.add_argument("--annotation", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    umi = args.umi.resolve()
    annotation = args.annotation.resolve()
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    output = outdir / "GSE132465_cell_library_sizes.tsv.gz"

    with gzip.open(umi, "rb") as handle:
        header = handle.readline().rstrip(b"\r\n").decode("utf-8").split("\t")
        cell_ids = header[1:]
        totals = np.zeros(len(cell_ids), dtype=np.int64)
        n_genes = 0
        for line_number, raw_line in enumerate(handle, start=2):
            tab = raw_line.find(b"\t")
            if tab < 0:
                continue
            values = np.fromstring(raw_line[tab + 1 :], sep="\t", dtype=np.int64)
            if values.size != totals.size:
                raise ValueError(
                    f"Line {line_number}: {values.size} values, expected {totals.size}"
                )
            totals += values
            n_genes += 1

    with output.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text:
                text.write("cell_id\ttotal_umi\n")
                for cell_id, total in zip(cell_ids, totals):
                    text.write(f"{cell_id}\t{int(total)}\n")

    manifest = {
        "source_umi": umi.name,
        "source_umi_sha256": sha256(umi),
        "source_annotation": annotation.name,
        "source_annotation_sha256": sha256(annotation),
        "library_size_output": output.name,
        "library_size_output_sha256": sha256(output),
        "n_cells": len(cell_ids),
        "n_gene_rows": n_genes,
    }
    (outdir / "extraction_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

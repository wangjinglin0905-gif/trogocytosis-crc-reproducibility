from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path


TEXT_SUFFIXES = {
    "", ".cff", ".csv", ".json", ".md", ".py", ".r", ".sh",
    ".ps1", ".tsv", ".txt", ".yml", ".yaml",
}
SKIP_PARTS = {".git", ".venv", "venv", "__pycache__", "release"}
REQUIRED = {
    "README.md", "CITATION.cff", ".zenodo.json", "LICENSE",
    "LICENSE_DATA.md", "DATA_AVAILABILITY.md", "VERSION",
    "provenance/external_sources.tsv", "provenance/figure_source_map.tsv",
}
FORBIDDEN_RAW_BASENAMES = {
    "CRISPRGeneEffect.csv",
    "Model.csv",
    "GSE178341_crc10x_full_c295v4_submit.h5",
    "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz",
    "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz",
    "GSE39582_eset_curatedCRCData.rda",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit a public reproducibility release.")
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--out", type=Path, default=Path("qa/public_release_audit.json"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    out = args.out if args.out.is_absolute() else root / args.out
    own_path = Path(__file__).resolve()

    errors: list[str] = []
    warnings: list[str] = []
    scanned_files = 0
    total_bytes = 0
    largest: list[tuple[int, str]] = []

    for required in sorted(REQUIRED):
        if not (root / required).is_file():
            errors.append(f"missing required file: {required}")

    local_path_patterns = [
        re.compile(r"(?i)(?<![A-Za-z0-9_])[A-Z]:[\\/](?:Users|KimiData|Documents|Temp)[\\/]"),
        re.compile(r"(?<![A-Za-z0-9_])/(?:Users|home)/[^/\s]+/"),
    ]
    credential_patterns = [
        re.compile(r"gh[pousr]_[A-Za-z0-9]{20,}"),
        re.compile(r"(?i)(?:access[_-]?token|api[_-]?key|client[_-]?secret|password)\s*[:=]\s*['\"][^'\"]{8,}['\"]"),
    ]

    for path in sorted(root.rglob("*")):
        rel = path.relative_to(root)
        if any(part in SKIP_PARTS for part in rel.parts):
            continue
        if path.is_symlink():
            errors.append(f"symlink not permitted in release: {rel.as_posix()}")
            continue
        if not path.is_file():
            continue

        scanned_files += 1
        total_bytes += path.stat().st_size
        largest.append((path.stat().st_size, rel.as_posix()))

        if path.name in FORBIDDEN_RAW_BASENAMES:
            errors.append(f"third-party raw file must not be redistributed: {rel.as_posix()}")
        if path.stat().st_size > 100 * 1024 * 1024:
            errors.append(f"file exceeds GitHub 100 MiB limit: {rel.as_posix()}")

        if path.resolve() == own_path or path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            warnings.append(f"text-like file was not UTF-8 and was skipped: {rel.as_posix()}")
            continue
        for pattern in local_path_patterns:
            if pattern.search(text):
                errors.append(f"local absolute path found in: {rel.as_posix()}")
                break
        for pattern in credential_patterns:
            if pattern.search(text):
                errors.append(f"possible credential found in: {rel.as_posix()}")
                break

    for figure_number in range(1, 8):
        matches = list((root / "figures" / "final").glob(f"Fig{figure_number}_*"))
        suffixes = {item.suffix.lower() for item in matches if item.is_file()}
        missing = {".png", ".tiff", ".pdf", ".svg"} - suffixes
        if missing:
            errors.append(
                f"Figure {figure_number} missing formats: {', '.join(sorted(missing))}"
            )

    report = {
        "status": "PASS" if not errors else "FAIL",
        "root_name": root.name,
        "scanned_files": scanned_files,
        "total_bytes": total_bytes,
        "largest_files": [
            {"path": rel, "bytes": size} for size, rel in sorted(largest, reverse=True)[:10]
        ],
        "required_files_checked": sorted(REQUIRED),
        "errors": sorted(set(errors)),
        "warnings": sorted(set(warnings)),
        "audit_script_sha256": sha256(own_path),
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()


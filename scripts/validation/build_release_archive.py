from __future__ import annotations

import argparse
import hashlib
import os
import zipfile
from pathlib import Path


SKIP_PARTS = {".git", ".venv", "venv", "__pycache__", "release"}
SKIP_QA_PREFIXES = (
    "qa/recomputed/figures_v7/",
    "qa/recomputed/figures_v8_1/",
    "qa/recomputed/figures_seeded_1/",
    "qa/recomputed/figures_seeded_2/",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def include_file(root: Path, path: Path) -> bool:
    rel = path.relative_to(root).as_posix()
    if any(part in SKIP_PARTS for part in Path(rel).parts):
        return False
    if any(rel.startswith(prefix) for prefix in SKIP_QA_PREFIXES):
        return False
    if path.name in {".DS_Store", "Thumbs.db"} or path.suffix.lower() == ".tmp":
        return False
    return path.is_file() and not path.is_symlink()


def collect(root: Path, manifest_path: Path) -> list[Path]:
    return sorted(
        path for path in root.rglob("*")
        if include_file(root, path) and path.resolve() != manifest_path.resolve()
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a deterministic public-release ZIP.")
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--version", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    release_dir = root / "release"
    release_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "provenance" / "file_manifest.tsv"

    files = collect(root, manifest_path)
    rows = ["relative_path\tbytes\tsha256"]
    for path in files:
        rows.append(
            f"{path.relative_to(root).as_posix()}\t{path.stat().st_size}\t{sha256(path)}"
        )
    manifest_path.write_text("\n".join(rows) + "\n", encoding="utf-8", newline="\n")

    files = collect(root, Path("__manifest_is_included__"))
    archive_name = f"trogocytosis-crc-reproducibility-v{args.version}.zip"
    archive_path = release_dir / archive_name
    archive_root = f"trogocytosis-crc-reproducibility-v{args.version}"
    fixed_time = (2026, 9, 2, 0, 0, 0)

    with zipfile.ZipFile(
        archive_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9
    ) as archive:
        for path in files:
            rel = path.relative_to(root).as_posix()
            info = zipfile.ZipInfo(f"{archive_root}/{rel}", date_time=fixed_time)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = (0o755 if path.suffix == ".sh" else 0o644) << 16
            archive.writestr(info, path.read_bytes(), compresslevel=9)

    checksum_path = release_dir / "SHA256SUMS.tsv"
    checksum_path.write_text(
        "file\tbytes\tsha256\n"
        f"{archive_name}\t{archive_path.stat().st_size}\t{sha256(archive_path)}\n",
        encoding="utf-8",
        newline="\n",
    )
    print(f"archive={archive_path}")
    print(f"bytes={archive_path.stat().st_size}")
    print(f"sha256={sha256(archive_path)}")
    print(f"files={len(files)}")


if __name__ == "__main__":
    main()


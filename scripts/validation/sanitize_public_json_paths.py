#!/usr/bin/env python3
"""Replace machine-local absolute paths in JSON provenance with basenames."""

from __future__ import annotations

import argparse
import json
import ntpath
import re
from pathlib import Path


ABSOLUTE_WINDOWS = re.compile(r"^[A-Za-z]:[\\/]")


def clean(value):
    if isinstance(value, dict):
        return {clean_key(key): clean(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean(item) for item in value]
    if isinstance(value, str) and ABSOLUTE_WINDOWS.match(value):
        return f"EXTERNAL_INPUT/{ntpath.basename(value)}"
    return value


def clean_key(value: str) -> str:
    if ABSOLUTE_WINDOWS.match(value):
        return f"EXTERNAL_INPUT/{ntpath.basename(value)}"
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    changed = []
    for path in sorted(args.root.rglob("*.json")):
        if ".git" in path.parts:
            continue
        data = json.loads(path.read_text(encoding="utf-8"))
        cleaned = clean(data)
        if cleaned != data:
            path.write_text(
                json.dumps(cleaned, ensure_ascii=False, indent=2) + "\n",
                encoding="utf-8",
            )
            changed.append(path.relative_to(args.root).as_posix())
    print(json.dumps({"changed": changed, "count": len(changed)}, indent=2))


if __name__ == "__main__":
    main()

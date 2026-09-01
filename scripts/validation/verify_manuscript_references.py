#!/usr/bin/env python
"""Verify DOI and PMID pairs in a Markdown reference list against public APIs."""

from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path
from urllib.parse import quote

import pandas as pd
import requests


def get_retry(session: requests.Session, url: str, *, params: dict | None = None, attempts: int = 4) -> requests.Response:
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            response = session.get(url, params=params, timeout=90)
            response.raise_for_status()
            return response
        except Exception as exc:
            last_error = exc
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"GET failed after {attempts} attempts: {url}") from last_error


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manuscript", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    text = args.manuscript.read_text(encoding="utf-8")
    markers = ["## 参考文献", "## References", "(style=Heading2) References"]
    reference_text = None
    for marker in markers:
        if marker in text:
            reference_text = text.split(marker, 1)[1]
            break
    if reference_text is None:
        raise ValueError("reference heading not found")
    entries = re.findall(
        r"(?ms)^\s*(?:\[P\d+\]\s*)?(\d+)\.\s+(.*?)(?=^\s*(?:\[P\d+\]\s*)?\d+\.\s+|\Z)",
        reference_text,
    )
    session = requests.Session()
    session.headers.update({"User-Agent": "trogocytosis-crc-reference-audit/1.0 (mailto:research@example.org)"})
    rows, raw = [], {}
    for number, entry in entries:
        doi_match = re.search(r"doi:([^\s;]+)", entry, flags=re.I)
        pmid_match = re.search(r"PMID:\s*(\d+)", entry, flags=re.I)
        doi = doi_match.group(1).rstrip(".") if doi_match else ""
        pmid = pmid_match.group(1) if pmid_match else ""
        crossref_ok = False
        datacite_ok = False
        crossref_title = ""
        crossref_doi = ""
        crossref_error = ""
        if doi:
            try:
                response = get_retry(session, f"https://api.crossref.org/works/{quote(doi, safe='')}")
                payload = response.json()["message"]
                raw[f"crossref_{number}"] = payload
                crossref_doi = str(payload.get("DOI", ""))
                crossref_title = " ".join(payload.get("title", [])).strip()
                crossref_ok = crossref_doi.lower() == doi.lower()
            except Exception as exc:  # explicit audit trail
                crossref_error = repr(exc)
                try:
                    response = get_retry(session, f"https://api.datacite.org/dois/{doi}")
                    payload = response.json()["data"]
                    raw[f"datacite_{number}"] = payload
                    datacite_ok = str(payload.get("id", "")).lower() == doi.lower()
                    if datacite_ok:
                        titles = payload.get("attributes", {}).get("titles", [])
                        crossref_title = str(titles[0].get("title", "")) if titles else ""
                except Exception as datacite_exc:
                    crossref_error = f"{crossref_error}; DataCite={datacite_exc!r}"
        pubmed_ok = False
        pubmed_title = ""
        pubmed_doi = ""
        pubmed_error = ""
        if pmid:
            try:
                response = get_retry(
                    session,
                    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                    params={"db": "pubmed", "id": pmid, "retmode": "json"},
                )
                payload = response.json()["result"].get(pmid, {})
                raw[f"pubmed_{number}"] = payload
                pubmed_title = str(payload.get("title", ""))
                pubmed_doi = next(
                    (str(x.get("value")) for x in payload.get("articleids", []) if x.get("idtype") == "doi"),
                    "",
                )
                pubmed_ok = str(payload.get("uid", "")) == pmid
            except Exception as exc:
                pubmed_error = repr(exc)
        rows.append({
            "reference_number": int(number),
            "manuscript_doi": doi,
            "crossref_record_found_and_doi_matches": crossref_ok,
            "datacite_record_found_and_doi_matches": datacite_ok,
            "crossref_title": crossref_title,
            "manuscript_pmid": pmid,
            "pubmed_record_found_and_pmid_matches": pubmed_ok,
            "pubmed_doi": pubmed_doi,
            "doi_matches_pubmed": (not doi or not pubmed_doi or doi.lower() == pubmed_doi.lower()),
            "pubmed_title": pubmed_title,
            "crossref_error": crossref_error,
            "pubmed_error": pubmed_error,
        })
        time.sleep(0.4)

    result = pd.DataFrame(rows).sort_values("reference_number")
    result.to_csv(args.output_dir / "reference_identifier_verification.tsv", sep="\t", index=False)
    (args.output_dir / "reference_api_records.json").write_text(json.dumps(raw, ensure_ascii=False, indent=2), encoding="utf-8")
    summary = {
        "reference_entries": int(len(result)),
        "doi_entries": int(result["manuscript_doi"].ne("").sum()),
        "crossref_doi_matches": int(result["crossref_record_found_and_doi_matches"].sum()),
        "pmid_entries": int(result["manuscript_pmid"].ne("").sum()),
        "pubmed_pmid_matches": int(result["pubmed_record_found_and_pmid_matches"].sum()),
        "doi_pubmed_mismatches": result.loc[~result["doi_matches_pubmed"], "reference_number"].astype(int).tolist(),
        "doi_registry_failures": result.loc[result["manuscript_doi"].ne("") & ~(result["crossref_record_found_and_doi_matches"] | result["datacite_record_found_and_doi_matches"]), "reference_number"].astype(int).tolist(),
        "pubmed_failures": result.loc[result["manuscript_pmid"].ne("") & ~result["pubmed_record_found_and_pmid_matches"], "reference_number"].astype(int).tolist(),
    }
    (args.output_dir / "reference_identifier_verification_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

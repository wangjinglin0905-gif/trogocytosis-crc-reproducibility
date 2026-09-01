# Data availability and provenance

## Public source data

All biological source data used in this project are public and de-identified. The authoritative sources are:

- CRC organoid CRISPR-screen supplementary data: DOI `10.1038/s41586-026-10830-y`.
- DepMap Public 26Q1 CRISPR Chronos gene effects and model metadata.
- MSigDB Hallmark collection, release 2023.2.
- GSE39582 through the Bioconductor `curatedCRCData` resource.
- TCGA COAD/READ PanCancer Atlas data through cBioPortal study `coadread_tcga_pan_can_atlas_2018`.
- GEO series GSE146771, GSE178341 and GSE132465.
- CMScaller templates from the CMScaller package.

Exact landing pages, retrieval dates, expected files and available SHA-256 values are recorded in `provenance/external_sources.tsv`.

## What this repository redistributes

The repository redistributes only author-generated code, parameter contracts, processed numerical derivatives, resampling distributions, source-data tables, quality-control reports and figures. These materials are sufficient to rerun the reported numerical checks and regenerate the final plots without exposing local paths or credentials.

The large external matrices, publisher supplementary workbooks, DepMap downloads, MSigDB gene-set file and curated Bioconductor object are not included. Users should obtain those materials from the named source repositories and follow their original licenses or terms of use.

## Identifiability and ethics

No direct identifiers, protected health information or controlled-access genomic data are included. Patient-level values in the package are derived from public, de-identified research datasets and are limited to the fields needed to reproduce the analyses.

## Versioning

Release `v1.0.0` corresponds to manuscript analysis version V8.1. A Git commit and a version-specific Zenodo record provide immutable version references. The Zenodo DOI will be inserted after final metadata review and publication of the deposition.


# Source data

Large public source files are not duplicated in Git. Place them here, or set
`PCT_DATA_DIR` to another directory.

GLDS-7 reference workflow inputs:

- `GSE56659_series_matrix.txt.gz` and `GPL198.annot.gz` from GEO GSE56659
- `GSE56659_RAW/` containing 36 compressed CEL files
- `ATH1-121501.CDF` for the independent Python RMA workflow
- `arabidopsis_MOESM5.xlsx`, Supplementary Data 4 from Barker et al. (2023)

Context workflow inputs go in `nasa_expansion_data/` and are the processed
GLDS-208 and OSD-251 tables cited in the manuscript. NASA accessions are OSD-7,
OSD-208, and OSD-251. Verify matching downloads with
`source-checksums.sha256`.

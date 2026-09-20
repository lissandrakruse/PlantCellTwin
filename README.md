# PlantCellTwin: NASA OSDR-grounded Arabidopsis digital twin

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22851568.svg)](https://doi.org/10.5281/zenodo.22851568)

This repository accompanies the manuscript “PlantCellTwin reveals a gravity-responsive root defense program in Arabidopsis”. It contains the source code, frozen result tables, software environment records, and checksums needed to inspect and reproduce the computational analyses.

## Project links

- Source code: https://github.com/lissandrakruse/PlantCellTwin
- Archived release and DOI: https://doi.org/10.5281/zenodo.22851568
- Interactive simulator: https://orbital-cell-twin-osdr.fuganti.chatgpt.site

## Authors

- Lissandra Kruse Fuganti (ORCID [0009-0008-8189-112X](https://orcid.org/0009-0008-8189-112X)), Universidade Estadual de Ponta Grossa
- Marcelo Giovanetti Canteri, Universidade Estadual de Ponta Grossa

## Evidence base

- NASA GLDS-7/OSD-7: tissue-resolved discovery analysis
- NASA GLDS-208/OSD-208: root-zone context validation
- NASA OSD-251: fractional-gravity trend validation

## Main components

- `app.js`, `index.html`, `styles.css`, and `plant.css`: interactive simulator
- `analysis/`: R and Python analysis workflows
- `results/`: frozen derived result tables
- Publication figures, manuscript files, supplementary information, and complete environment records are included in the versioned Zenodo archival package

## Quick start

Serve the static interface with `python -m http.server 8000`, then open
`http://localhost:8000`. For analysis, create the environment and obtain the
public inputs listed in `data/README.md`:

```bash
conda env create -f environment.yml
conda activate plantcelltwin
python analysis/analyze_glds7.py
```

Workflows read from `data/` and write to `results/`. Set `PCT_DATA_DIR` and
optionally `PCT_RESULTS_DIR` to use other locations.

## Verify the release

Frozen result tables provide a small open test fixture:

```bash
./run_verification.sh
```

The script checks required assets, expected GLDS-7 R/Python concordance,
OSD-251 panel direction, syntax, and machine-specific paths in the interface.

## Support and contributions

Use the GitHub issue tracker for support and bug reports. Pull requests are
welcome; see `CONTRIBUTING.md`.

## Archived release

Version 1.0.0 is permanently archived on Zenodo:

- DOI: [10.5281/zenodo.22851568](https://doi.org/10.5281/zenodo.22851568)
- GitHub: [lissandrakruse/PlantCellTwin](https://github.com/lissandrakruse/PlantCellTwin)
- Resource type: Software
- Source-code license: MIT
- Manuscript and figure license: CC BY 4.0

## Citation

Fuganti, L. K., & Canteri, M. G. (2026). *PlantCellTwin: NASA OSDR-grounded Arabidopsis digital twin* (Version 1.0.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.22851568

BibTeX:

```bibtex
@software{fuganti_canteri_2026_plantcelltwin,
  author     = {Fuganti, Lissandra Kruse and Canteri, Marcelo Giovanetti},
  title      = {PlantCellTwin: NASA OSDR-grounded Arabidopsis digital twin},
  year       = {2026},
  version    = {1.0.0},
  publisher  = {Zenodo},
  doi        = {10.5281/zenodo.22851568},
  url        = {https://doi.org/10.5281/zenodo.22851568},
  repository = {https://github.com/lissandrakruse/PlantCellTwin}
}
```

## Interpretation

The package supports a reproducible computational reanalysis and hypothesis-generating digital twin. Flight-versus-ground contrasts do not isolate gravity, and cross-dataset validation does not establish causal mechanism. Prospective perturbation experiments remain necessary.

## Keywords

*Arabidopsis thaliana*; NASA GeneLab; OSDR; spaceflight; microgravity; transcriptomics; digital twin; reproducibility.

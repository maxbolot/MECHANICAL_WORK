# Project Architecture Map

This document gives local coding agents a quick map of where things live in this repository.

## Root Overview

- `README.md`: Primary project overview and usage notes.
- `TODO.md`: Work backlog and pending tasks.
- `GLOBALFV3_DATA_LAYOUT_MIGRATION_PLAN.md`: Data layout migration planning notes.
- `pull_data_from_HPC.sh`: Pulls data from HPC storage.
- `send_to_archive.sh`: Sends outputs to archival storage.
- `setup_matlab_proxy.sh`: MATLAB environment/proxy setup helper.

## Directory Map

### `launcher/`
Shell launch scripts for local and async compute runs (mostly no-SLURM wrappers), plus list files used as run inputs.

- Compute launchers:
  - `compute_work_async_*.sh`
  - `compute_prate_threshold*_noslurm.sh`
- Input lists:
  - `launcher/list/*.txt`

Use this area when changing run orchestration, batching, or job argument wiring.

### `local/`
Fortran source and local build entrypoint.

- `local/Makefile`: Local compile/build commands.
- `local/src/*.f90`: Core numerical executables:
  - work and async work computation
  - precipitation-threshold variants
  - latitude-band threshold variants
- `local/README_LAT_BAND_THRESHOLDS.md`: Notes for lat-band threshold logic.

Use this area for algorithmic/compute-kernel changes.

### `matlab/`
Analysis and figure-generation code.

- Legacy scripts:
  - `matlab/figure1_legacy.m`
  - `matlab/figure2_legacy.m`
- Main analysis scripts:
  - `matlab/analysis/*.m` (region masks, thresholded stats, maps, histograms, regression checks, panel generation)
- Supporting material:
  - `matlab/lib/`, `matlab/presets/`, `matlab/entries/`

Use this area for publication figures and post-compute diagnostics in MATLAB.

### `notebooks/`
Jupyter notebooks for exploration, validation, and plotting workflows.

Examples:
- `compare_tropical_work_control_vs_plus4k_prate_thresholded.ipynb`
- `compute_and_plot_work_lift_all_regions.ipynb`
- `plot_work_ratio_lonlat_map.ipynb`

Use this area for exploratory analysis and reproducible notebook-based checks.

### `publication/`
Paper-specific assets and figure source code.

- `publication/paper_dissipation/`
- `publication/paper_work_ratio/`

Within paper folders, expect figure source scripts (including MATLAB figure builders) and publication-ready outputs.

### `python/`
Python utilities and translated logic.

- `python/make_ea_grid.py`: Grid generation helper.
- `python/translated_from_matlab/`: MATLAB-to-Python translated routines.

Use this area for Python tooling or migration away from MATLAB.

### `script/`
Pipeline scripts grouped by stage.

- `script/preprocessing/`
- `script/postprocessing/`

Use this area for data preparation and output post-processing automation.

## Typical Workflow (High Level)

1. Prepare data with `script/preprocessing/` and transfer helpers at repo root.
2. Launch compute runs via `launcher/*.sh`.
3. Execute numerical kernels from `local/src/*.f90` (built via `local/Makefile`).
4. Analyze outputs using `matlab/analysis/`, `python/`, or `notebooks/`.
5. Generate publication figures under `publication/`.

## Quick Navigation for Agents

- Need to change core compute logic: `local/src/`
- Need to change run execution wiring: `launcher/`
- Need to change figure/analysis output: `matlab/analysis/` or `publication/`
- Need Python utility or migration work: `python/`
- Need prep/post pipeline updates: `script/preprocessing/`, `script/postprocessing/`

## Notes

- This repository appears to combine production compute code (Fortran + shell) with analysis/publication pipelines (MATLAB + notebooks + Python).
- Prefer minimal, stage-specific edits: launchers for orchestration, `local/src` for numerical logic, analysis folders for visualization/statistics.

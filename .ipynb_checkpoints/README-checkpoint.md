# MECHANICAL_WORK

End-to-end workflow for computing mechanical work diagnostics from FV3 coarse-grained fields, post-processing outputs into analysis-ready products, and generating regional climatological summaries in MATLAB.

Notebook content under `notebooks/` is provided for dissemination and sharing; MATLAB remains the authoritative analysis stack, and there is no roadmap to move away from MATLAB workflows.

The project is organized as a pipeline:

1. Low-level Fortran kernels compute per-date work products.
2. Launcher scripts run those kernels over many dates (typically via Slurm arrays).
3. Optional preprocessing scripts remap native high-resolution fields onto the analysis grid.
4. Post-processing scripts validate time axes and concatenate outputs.
5. MATLAB entry wrappers call analysis runners backed by the shared lib API to produce diagnostics and plots.
6. Publication submodules manage paper-specific TeX/Bib and final figure assets.

## Project Structure

- `local/`
  - `src/compute_work.f90`: baseline Fortran implementation.
  - `src/compute_work_async.f90`: async/pipelined Fortran implementation (current default binary in launchers).
  - `src/compute_prate_thresholds.f90`: two-pass precipitation-threshold generator.
  - `src/compute_work_async_prate_threshold.f90`: async work/lift fork with precipitation-threshold masking.
  - `Makefile`: build rules for both binaries.
  - `bin/`: compiled executables.
  - `obj/`: object files.
- `launcher/`
  - `compute_work_async_array.sh`: Slurm array launcher using simulation-aware part1/part2 date lists.
  - `compute_work_async_array_inner_parallel.sh`: Slurm array launcher with multiple inner jobs per array task.
  - `compute_work_async_noslurm.sh`: serial fallback runner using simulation-aware part1/part2 date lists.
  - `compute_work_async_test_noslurm.sh`: simple one-date test launcher with one hardcoded source root.
  - `compute_prate_threshold_noslurm.sh`: serial launcher for percentile threshold generation.
  - `compute_work_async_prate_threshold_array.sh`: Slurm array launcher for thresholded async work/lift.
  - `compute_work_async_prate_threshold_noslurm.sh`: serial launcher for thresholded async work/lift.
  - `list/`: simulation-specific date lists for part1/part2 source roots.
- `script/`
  - `preprocessing/`
    - `ea_to_0.25_array.sh`: Slurm array job that remaps native-resolution fields to the 0.25-degree Fortran input grid via conservative equal-area remapping.
  - `postprocessing/`
    - `check_work_time_axis.sh`: checks timestamp drift in per-date work files; optionally removes drifted files.
    - `check_hist_time_axis.sh`: checks timestamp drift in per-date histogram files; optionally removes drifted files.
    - `coarse_grain_and_concat_work.sh`: remaps per-date work files to a target grid, then concatenates into a single date-range file.
    - `coarse_grain_and_concat_precip.sh`: remaps per-date precipitation files, concatenates to a date-range product, and can apply daily aggregation.
    - `concat_histograms.sh`: concatenates per-date histogram files into a single date-range file.
- `python/`
  - `make_ea_grid.py`: generates the CDO equal-area grid description file (`output/grids/ea_25km_grid.txt`) required by `script/preprocessing/ea_to_0.25_array.sh`.
- `matlab/`
  - `entries/`: thin wrapper scripts that define config and call runners.
  - `analysis/`: runner/orchestration scripts.
  - `lib/`: shared MATLAB API helpers.
  - `presets/`: scenario and region presets used by entries/runners.
- `logs/`: Slurm and run logs.
- `publication/`
  - `paper_dissipation/`: paper-specific manuscript and figure workspace tracked as a Git submodule.
  - `.gitmodules` in repo root points to the submodule URL.
- `setup_matlab_proxy.sh`: helper to launch MATLAB proxy environment.

## Low-Level Fortran Logic

### Precipitation Threshold Generator: `compute_prate_thresholds.f90`

This program builds precipitation percentile thresholds from PRATE files in two passes.

- Inputs are read from namelist:
  - `history_root_part1`, `date_list_file_part1`
  - `history_root_part2`, `date_list_file_part2`
  - `output_file`
- Pass 1:
  - Parses all `PRATEsfc_coarse_C3072_1440x720.fre.nc` files listed across both list files.
  - Accumulates area-weighted histogram using `bin_sumarea` from `libhist`.
  - Applies clipped-domain histogram accumulation with hard bounds:
    - latitude: `-30` to `30`
    - longitude: `0` to `360`
- Pass 2:
  - Builds normalized CDF from the accumulated histogram.
  - Computes inverse-CDF precipitation thresholds for requested percentiles using log-log interpolation via `interp1` from `libinterp`.
- Output:
  - ASCII thresholds file.
  - Default launcher paths are simulation-specific:
    - control: `output/thresholds/thresholds_control.txt`
    - warming: `output/thresholds/thresholds_warming.txt`

### Thresholded Async Work/Lift Fork: `compute_work_async_prate_threshold.f90`

This fork preserves async task-pipelined I/O/compute from `compute_work_async.f90` and applies precipitation-threshold masking.

- Reads percentile thresholds from ASCII (`path_thresholds` namelist item).
- Computes work/lift, then applies `PRATE >= threshold` for each percentile.
- Stores zero where threshold is not met.
- Removes histogram outputs/logic from this fork.
- Aggregates native timesteps into daily bins and writes daily means (not daily sums).
- Writes NetCDF `work` and `lift` with leading `percentile` dimension.
- Adds traceability global attributes such as `p99_threshold` and `threshold_file`.
- Adds aggregation metadata:
  - global: `daily_stat=mean`, `daily_aggregation="mean over native timesteps within each day"`
  - variable: `time_stat=daily_mean`, `time_aggregation="mean over native timesteps within each day"` on `work` and `lift`

### Inputs and Outputs

The work kernels read a namelist (`&config`) specifying paths to all required coarse variables (thermodynamic tendencies, hydrometeors, precipitation, etc.) and write:

- A work file with variables:
  - `work` (W m^-2)
  - `lift` (W m^-2)
- A histogram file containing precipitation/work binned diagnostics.

Time aggregation behavior for all `compute_work*.f90` kernels:

- Native model timesteps are grouped into contiguous day bins from the input time axis.
- Per-day outputs are normalized by the number of timesteps in each day.
- NetCDF outputs therefore store daily means.

Each run processes one date directory and writes one `work_YYYYMMDDHH.nc` and one `hist_YYYYMMDDHH.nc`.

### Baseline Kernel: `compute_work.f90`

Core characteristics:

- Uses OpenMP and NetCDF Fortran APIs.
- Processes the horizontal domain in y-chunks (`chunk_size = 144`) with two buffers (`nbuf = 2`).
- Computes daily-mean full-grid work/lift products and selected daily-mean histogram diagnostics.
- Applies clipped-domain bounds for histogram accumulation (currently hardcoded in code):
  - latitude: -30 to 30
  - longitude: 0 to 360
- Writes CF-style coordinate metadata and carries through input time metadata when available.
- Writes explicit daily-aggregation metadata in NetCDF attributes (global and per-variable).

### Async Kernel: `compute_work_async.f90`

This is the same scientific computation, but with overlapped I/O and compute:

- Uses OpenMP tasks and dependency tokens for a two-stage pipeline:
  - read task (serialized NetCDF I/O in a critical section)
  - compute task (operates after read completion)
- Master thread can write completed buffers while worker threads read next chunks.
- Buffer state machine reduces idle time and improves throughput on shared filesystems.

In launchers, `local/bin/compute_work_async` is preferred by default.

## Build System

Build from `local/`:

- Compiler: Intel `ifx`
- OpenMP enabled (`-qopenmp`)
- NetCDF flags from `nf-config`
- External support libraries expected in `$HOME/local`:
  - `-linterp`
  - `-lnc`
  - `-lhist`

Targets:

- `bin/compute_work`
- `bin/compute_work_async`
- `bin/compute_prate_thresholds`
- `bin/compute_work_async_prate_threshold`

## Launch Layer (Date-Oriented Execution)

### `launcher/compute_work_async_array.sh`

Primary Slurm array driver:

- Uses `SIMULATION=${SIMULATION:-control}` to choose between control and warming inputs.
- Reads dates from two simulation-specific list files:
  - `LIST_FILE_PART1`
  - `LIST_FILE_PART2`
- Maps each list to its matching source root:
  - `SOURCE_ROOT_PART1`
  - `SOURCE_ROOT_PART2`
- Submits itself as an array (`--array=1-N`).
- For each task:
  - selects one date from the combined part1+part2 task list
  - resolves the matching source root for that date
  - writes a temporary namelist
  - runs compute binary
  - removes temporary namelist

Key environment knobs:

- `SIMULATION` (`control` or `warming`)
- `NTASKS`
- `CPUS_PER_TASK`
- `MEM_PER_CPU` (default `3G`)
- `LOG_DIR`

### `launcher/compute_work_async_array_inner_parallel.sh`

Array launcher with inner concurrency:

- Uses the same `SIMULATION` selector and combined part1/part2 list logic as `compute_work_async_array.sh`.
- One array task can run multiple dates concurrently (`INNER_PARALLEL`).
- Splits available CPU threads across inner jobs by setting `OMP_NUM_THREADS` per instance.
- Produces per-date logs for easier failure isolation.

### `launcher/compute_work_async_noslurm.sh`

Serial fallback for non-Slurm environments:

- Uses `SIMULATION=${SIMULATION:-control}` to choose control or warming inputs.
- Iterates both part1 and part2 list files in sequence.
- Uses the correct source root for each list.
- Runs one date at a time.
- Writes per-date logs.

### `launcher/compute_work_async_test_noslurm.sh`

Minimal one-date validation launcher:

- Uses one hardcoded source root:
  - `/scratch/cimes/GLOBALFV3/20191020.00Z.C3072.L79x2_pire/history`
- Uses one hardcoded date:
  - `2020010300` (the first entry from `launcher/list/list_control_part1.txt`)
- Writes a single temporary namelist and runs `compute_work_async` once.
- Intended for quick smoke tests before launching full control or warming batches.

### `launcher/compute_prate_threshold_noslurm.sh`

Serial threshold-generation launcher:

- Uses `SIMULATION=${SIMULATION:-control}` to choose control or warming inputs.
- Writes one temporary namelist with two roots and two list files.
- Runs `compute_prate_thresholds` once over all dates.
- Writes threshold ASCII output to a simulation-specific path:
  - control: `output/thresholds/thresholds_control.txt`
  - warming: `output/thresholds/thresholds_warming.txt`

### `launcher/compute_work_async_prate_threshold_array.sh` and `launcher/compute_work_async_prate_threshold_noslurm.sh`

Thresholded async work/lift launchers:

- Use `SIMULATION=${SIMULATION:-control}` to choose control or warming inputs.
- Read dates from paired part1/part2 list files and use the matching source root for each date.
- Pass `path_thresholds` in namelist to `compute_work_async_prate_threshold`.
- Select simulation-specific threshold inputs by default:
  - control: `output/thresholds/thresholds_control.txt`
  - warming: `output/thresholds/thresholds_warming.txt`
- Produce threshold-masked work NetCDF outputs per date.

### Common Launcher Override Pattern

All current launchers default to the control simulation and can be switched manually:

- control: `./launcher/<script>.sh`
- warming: `SIMULATION=warming ./launcher/<script>.sh`

Advanced overrides remain available through environment variables such as:

- `LIST_FILE_PART1`
- `LIST_FILE_PART2`
- `TARGET_DIR_COMPUTE`
- `TARGET_DIR_HISTOGRAMS`
- `THRESHOLD_FILE`

Launcher-generated namelists are now written under `output/namelists/`.
They are deleted after successful runs, but a failed Fortran invocation leaves the namelist behind for inspection.

## Post-Processing Layer

### 1) Time-Axis Quality Checks

`script/postprocessing/check_work_time_axis.sh` scans `work_YYYYMMDDHH.nc` (or thresholded `work_prate_threshold_YYYYMMDDHH.nc`) files and compares:

- expected start timestamp inferred from the filename (`YYYYMMDDHH`)
- actual first timestamp in file (`cdo showtimestamp`)

`script/postprocessing/check_hist_time_axis.sh` performs the same drift check for `hist_YYYYMMDDHH.nc` files.

Both scripts support:

- `MAX_DIFF_DAYS` to control tolerated drift.
- `DRY_RUN=1` to preview actions.
- `REMOVE_SOURCE_FILES=1` to remove files that exceed the drift threshold.

These checks do not modify files unless `REMOVE_SOURCE_FILES=1` is explicitly set.

### 2) Remap Native Fields to 0.25-Degree Inputs (Optional)

`script/preprocessing/ea_to_0.25_array.sh` prepares Fortran-ready coarse inputs from native-resolution fields:

- Reads dates from a list file and processes each date directory.
- Generates conservative remap weights for native -> equal-area and equal-area -> `r1440x720`.
- Produces mean and covariance fields expected by the work kernels.
- Supports split execution by field group (`RUN_GROUP=group1|group2|all`).

### 3) Remap + Concatenate Work Files

`script/postprocessing/coarse_grain_and_concat_work.sh`:

- Collects canonical per-date work files (`work_YYYYMMDDHH.nc`).
- If `*.taxis_repaired.nc` companions are present, uses them instead of original files.
- Remaps each file with CDO (default `remapcon,r360x180`).
- Concatenates remapped files with `ncrcat` into a date-range product:
  - `work_START_END.nc`

### 4) Remap + Concatenate Precipitation Files

`script/postprocessing/coarse_grain_and_concat_precip.sh`:

- Supports six simulation switches:
  - `control`
  - `warming`
  - `control_prate_thresholded`
  - `warming_prate_thresholded`
  - `control_prate_thresholded_by_lat_band`
  - `warming_prate_thresholded_by_lat_band`
- Uses two input modes under one interface:
  - list-driven mode for `control`/`warming` from part1/part2 source roots and date lists.
  - directory mode for thresholded outputs from `work_prate_threshold*.nc` files.
- Extracts precipitation into a unified variable name `precip` across all modes.
- Remaps each file with CDO (default `remapcon,r360x180`).
- Concatenates remapped files with `ncrcat` into:
  - `precip_START_END.nc`

Use this step when analyses require precipitation explicitly (for example Hp/precip diagnostics, Hp delta maps, and decomposition terms involving `P`).
For work/lift-only analyses, this step is optional.

### 5) Concatenate Histograms

`script/postprocessing/concat_histograms.sh`:

- Collects canonical per-date histogram files (`hist_YYYYMMDDHH.nc`).
- Concatenates chronologically to `hist_START_END.nc`.
- If `time` is fixed-size, converts it to a record dimension (`ncks --mk_rec_dmn time`) before `ncrcat`.

## Python Utilities

`python/make_ea_grid.py` generates the CDO grid description file for the equal-area (EA) intermediate grid used during native-to-0.25-degree remapping.

- Grid layout: 1440 × 720, longitude spacing 0.25° uniform, latitude spacing uniform in sin(lat) (equal area).
- Cell centers are computed as midpoints in sine space to ensure exact area conservation.
- Writes `output/grids/ea_25km_grid.txt` in CDO `lonlat` grid format including `xvals`, `xbounds`, `yvals`, `ybounds`.
- Must be run once before generating remap weights with `script/preprocessing/ea_to_0.25_array.sh`.

Run:
```bash
python python/make_ea_grid.py
```

## Notebook Layer (Dissemination)

The `notebooks/` folder contains Python/Jupyter views of selected analyses intended for dissemination, communication, and easier sharing with collaborators.

Important scope note:

1. Notebooks are a dissemination layer, not a replacement for the production MATLAB analysis stack.
2. There is currently no roadmap to move away from MATLAB for analysis workflows.
3. MATLAB (`matlab/entries`, `matlab/analysis`, `matlab/lib`) remains the authoritative analysis implementation.

Translated Python helpers used by notebooks are stored under `python/translated_from_matlab/` to keep notebook orchestration lightweight and reusable.

## MATLAB Analysis Layer

MATLAB consumes concatenated NetCDF products (`work_START_END.nc`, `precip_START_END.nc`, `hist_START_END.nc`) through a three-layer structure:

1. Thin wrapper entries in `matlab/entries`
2. Runner/orchestration scripts in `matlab/analysis`
3. Reusable API helpers in `matlab/lib`

### Thin Wrapper Entries (`matlab/entries`)

Entry scripts are intentionally minimal. They:

1. Set paths (`lib`, `analysis`, `presets` when needed).
2. Build a small `cfg` struct (scenario, mode, plotting switches, dates, thresholds).
3. Call one analysis runner.

Examples:

- `matlab/entries/plot_work_ratio_lonlat_map.m`
- `matlab/entries/plot_hp_precip_delta_lonlat_map.m`
- `matlab/entries/compute_and_plot_work_lift_all_regions.m`
- `matlab/entries/compute_histogram_precip_work_analysis_control.m`
- `matlab/entries/compute_histogram_precip_work_analysis_plus4k.m`

### Runner Scripts (`matlab/analysis`)

Runners contain workflow orchestration and plotting composition. They:

1. Apply defaults to `cfg`.
2. Read inputs using lib API helpers.
3. Execute scenario/region logic.
4. Produce summaries and figures.
5. Return structured outputs when useful for checks/regression.

Representative runners:

- `matlab/analysis/run_work_ratio_lonlat_map_analysis.m`
- `matlab/analysis/run_hp_precip_delta_lonlat_map_analysis.m`
- `matlab/analysis/run_work_lift_all_regions_analysis.m`
- `matlab/analysis/run_histogram_precip_work_analysis.m`

### Lib API (`matlab/lib`)

`matlab/lib` is the shared API surface used by all runners. It centralizes:

- NetCDF input validation and reading helpers
- Dimension reordering/orientation helpers
- Time-axis normalization and datetime conversion helpers
- Time-weighting and missing-step counting helpers
- Shared weighted averaging and histogram metric helpers
- Region/threshold reusable computation helpers

Examples:

- `normalize_time_axis_to_datenum.m`, `netcdf_time_to_datetime.m`, `build_time_axis.m`
- `compute_time_weights_control.m`, `compute_time_weights_plus4k.m`, `compute_time_weights_from_days.m`
- `read_var_as_lon_lat_time.m`, `read_lon_lat_time_from_var.m`, `read_tropical_histogram_diagnostics.m`
- `weighted_nanmean.m`, `weighted_mean_over_time_map.m`, `weighted_time_mean_rows.m`

### Time-Weighting Behavior

Current weighting assumes a uniform 1-day cadence across supported products.
Missing-step counts are computed from the normalized daily time axis.

### Roadmap: Lib Extraction

The lib API is being shaped as a stable shared layer and is expected to be detached into a separate project in the future. The goal is to:

1. Keep wrappers and runners repository-local.
2. Publish versioned reusable MATLAB APIs for NetCDF/time/weighting/stat utilities.
3. Reduce duplication across analysis repositories and publication workflows.

## Typical End-to-End Flow

1. Build Fortran binaries in `local/`.
2. Prepare date lists in `launcher/list/`.
3. Run launcher (Slurm array or serial fallback) to generate per-date `work_*.nc` and `hist_*.nc`.
4. Optional: run `script/preprocessing/ea_to_0.25_array.sh` if native-to-analysis-grid preprocessing is needed.
5. Run `script/postprocessing/check_work_time_axis.sh` and `script/postprocessing/check_hist_time_axis.sh` to detect drift before concatenation.
6. Run `script/postprocessing/coarse_grain_and_concat_work.sh` to produce concatenated `work_START_END.nc`.
7. Optionally run `script/postprocessing/coarse_grain_and_concat_precip.sh` to produce concatenated `precip_START_END.nc` for precipitation-coupled analyses.
8. Optionally run `script/postprocessing/concat_histograms.sh` for histogram products.
9. Run MATLAB entry scripts from `matlab/entries` (regional, map, histogram, or all-regions summaries).
10. Optionally run dissemination notebooks from `notebooks/` for shareable views of selected analyses.

## Publication and Manuscript Workflow

Paper workspaces are managed under `publication/` as separate Git submodules.

Current paper workspace:

- `publication/paper_dissipation/`
  - `main.tex`, `references.bib`, `styles_and_macros.sty`
  - `figures/`, `tables/`
  - `figures_src/`, `logs/`, `manifests/`
  - `Makefile`, `.gitignore`, `.overleafignore`

Overleaf sync is controlled from inside the paper submodule with `.overleafignore`.
By default, non-manuscript assets (`figures_src/`, `logs/`, `manifests/`, `Makefile`, `.gitignore`) are excluded from Overleaf to reduce sync clutter.

## Environment Notes

- Scripts rely on environment modules (Intel, HDF5, NetCDF, CDO, NCO).
- Paths in launchers/scripts are cluster-specific and may need adjustment per system.
- `setup_matlab_proxy.sh` provides one workflow for MATLAB proxy usage.

## Operational Caveats

- Missing-step counts generally indicate missing input dates, but can also reflect malformed metadata if pre-processing is skipped.
- Per-date files with empty time dimensions or invalid timestamps should be repaired or regenerated before concatenation.
- Keep date lists and source directory selection consistent between control and warming pipelines.

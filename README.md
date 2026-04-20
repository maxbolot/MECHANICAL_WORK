# MECHANICAL_WORK

End-to-end workflow for computing mechanical work diagnostics from FV3 coarse-grained fields, post-processing outputs into analysis-ready products, and generating regional climatological summaries in MATLAB.

The project is organized as a pipeline:

1. Low-level Fortran kernels compute per-date work products.
2. Launcher scripts run those kernels over many dates (typically via Slurm arrays).
3. Optional preprocessing scripts remap native high-resolution fields onto the analysis grid.
4. Post-processing scripts validate time axes and concatenate outputs.
5. MATLAB scripts produce regional weighted averages and comparison plots.
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
    - `concat_histograms.sh`: concatenates per-date histogram files into a single date-range file.
- `matlab/`
  - Region-specific analysis scripts for control and warming simulations.
  - `compute_and_plot_work_lift_all_regions.m`: combined multi-region summary + plotting script.
  - `lib/`: shared MATLAB helpers.
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

### 4) Concatenate Histograms

`script/postprocessing/concat_histograms.sh`:

- Collects canonical per-date histogram files (`hist_YYYYMMDDHH.nc`).
- Concatenates chronologically to `hist_START_END.nc`.
- If `time` is fixed-size, converts it to a record dimension (`ncks --mk_rec_dmn time`) before `ncrcat`.

## MATLAB Analysis Layer

MATLAB scripts consume concatenated products (`work_START_END.nc`) and compute regionally averaged diagnostics.

### Shared Helper Library (`matlab/lib`)

- `compute_time_weights_control.m`
- `compute_time_weights_plus4k.m`
- `weighted_nanmean.m`
- `build_time_axis.m`
- `read_var_as_lon_lat_time.m`
- plus utility helpers

### Time-Weighting and Missing-Step Logic

Both weighting helpers:

- Convert time to MATLAB datenum robustly across unit conventions.
- Assign per-sample duration weights according to simulation schedule.
- Detect missing expected outputs.

Current robust missing-step strategy:

1. Prefer parsing canonical input dates from global `history` (ncrcat command paths).
2. Count missing expected start dates against schedule.
3. Fallback to timestamp-axis matching when history is unavailable.

This avoids false missing counts caused by midpoint phase shifts at cadence transitions.

### Region Scripts

Each region has control and warming variants, for example:

- global
- tropics
- maritime continent
- northern midlatitudes
- southern midlatitudes
- southern poles

Each script typically:

1. Reads concatenated work/lift fields.
2. Applies latitude/longitude mask.
3. Computes area-weighted spatial means (`cos(lat)`).
4. Applies schedule-aware weighted temporal mean.
5. Reports work, lift, and lift/work ratio.
6. Plots time series.

### Combined Multi-Region Script

`matlab/compute_and_plot_work_lift_all_regions.m`:

- Processes both control and warming files.
- Computes region-level summary metrics.
- Produces stacked-bar comparison of `lift` and `ke = work - lift`.
- Prints summary table including detected missing timesteps.

## Typical End-to-End Flow

1. Build Fortran binaries in `local/`.
2. Prepare date lists in `launcher/list/`.
3. Run launcher (Slurm array or serial fallback) to generate per-date `work_*.nc` and `hist_*.nc`.
4. Optional: run `script/preprocessing/ea_to_0.25_array.sh` if native-to-analysis-grid preprocessing is needed.
5. Run `script/postprocessing/check_work_time_axis.sh` and `script/postprocessing/check_hist_time_axis.sh` to detect drift before concatenation.
6. Run `script/postprocessing/coarse_grain_and_concat_work.sh` to produce concatenated `work_START_END.nc`.
7. Optionally run `script/postprocessing/concat_histograms.sh` for histogram products.
8. Run MATLAB regional scripts or the all-regions summary script.

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

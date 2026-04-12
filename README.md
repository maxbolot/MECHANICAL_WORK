# MECHANICAL_WORK

End-to-end workflow for computing mechanical work diagnostics from FV3 coarse-grained fields, post-processing outputs into analysis-ready products, and generating regional climatological summaries in MATLAB.

The project is organized as a pipeline:

1. Low-level Fortran kernels compute per-date work products.
2. Launcher scripts run those kernels over many dates (typically via Slurm arrays).
3. Processing scripts repair/normalize time axes and concatenate outputs.
4. MATLAB scripts produce regional weighted averages and comparison plots.

## Project Structure

- `local/`
  - `src/compute_work.f90`: baseline Fortran implementation.
  - `src/compute_work_async.f90`: async/pipelined Fortran implementation (current default binary in launchers).
  - `Makefile`: build rules for both binaries.
  - `bin/`: compiled executables.
  - `obj/`: object files.
- `launcher/`
  - `compute_work_array.sh`: Slurm array launcher (one date per task).
  - `compute_work_array_inner_parallel.sh`: Slurm array launcher with multiple inner jobs per array task.
  - `compute_work_noslurm.sh`: serial fallback runner.
  - `list.txt`: date list used by launchers.
- `script/`
  - `repair_work_time_axis.sh`: detects and repairs inconsistent time coordinates in per-date work files.
  - `coarse_grain_and_concat_work.sh`: remaps per-date files to a target grid, then concatenates.
  - `concat_histograms.sh`: concatenates per-date histogram files.
- `matlab/`
  - Region-specific analysis scripts for control and warming simulations.
  - `compute_and_plot_work_lift_all_regions.m`: combined multi-region summary + plotting script.
  - `lib/`: shared MATLAB helpers.
- `logs/`: Slurm and run logs.
- `setup_matlab_proxy.sh`: helper to launch MATLAB proxy environment.

## Low-Level Fortran Logic

### Inputs and Outputs

Both Fortran programs read a namelist (`&config`) specifying paths to all required coarse variables (thermodynamic tendencies, hydrometeors, precipitation, etc.) and write:

- A work file with variables:
  - `work` (W m^-2)
  - `lift` (W m^-2)
- A histogram file containing precipitation/work binned diagnostics.

Each run processes one date directory and writes one `work_YYYYMMDDHH.nc` and one `hist_YYYYMMDDHH.nc`.

### Baseline Kernel: `compute_work.f90`

Core characteristics:

- Uses OpenMP and NetCDF Fortran APIs.
- Processes the horizontal domain in y-chunks (`chunk_size = 144`) with two buffers (`nbuf = 2`).
- Computes full-grid time-mean work/lift products and selected histogram diagnostics.
- Applies clipped-domain bounds for histogram accumulation (currently hardcoded in code):
  - latitude: -30 to 30
  - longitude: 0 to 360
- Writes CF-style coordinate metadata and carries through input time metadata when available.

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

## Launch Layer (Date-Oriented Execution)

### `launcher/compute_work_array.sh`

Primary Slurm array driver:

- Reads one date per line from `launcher/list.txt` (or `LIST_FILE`).
- Submits itself as an array (`--array=1-N`).
- For each task:
  - selects one date
  - writes a temporary namelist
  - runs compute binary
  - removes temporary namelist

Key environment knobs:

- `NTASKS`
- `CPUS_PER_TASK`
- `MEM_PER_CPU` (default `3G`)
- `LOG_DIR`

### `launcher/compute_work_array_inner_parallel.sh`

Array launcher with inner concurrency:

- One array task can run multiple dates concurrently (`INNER_PARALLEL`).
- Splits available CPU threads across inner jobs by setting `OMP_NUM_THREADS` per instance.
- Produces per-date logs for easier failure isolation.

### `launcher/compute_work_noslurm.sh`

Serial fallback for non-Slurm environments:

- Iterates dates from a list file.
- Runs one date at a time.
- Writes per-date logs.

## Post-Processing Layer

### 1) Time-Axis Repair

`script/repair_work_time_axis.sh` scans `work_YYYYMMDDHH.nc` files and compares:

- expected midpoint timestamp inferred from filename and cadence rules
- actual first timestamp in file (`cdo showtimestamp`)

Default control cadence used for repair script:

- 5-day before `2021-05-27`
- 1-day on/after `2021-05-27`

If absolute drift exceeds `MAX_DIFF_DAYS` (default 5), it writes:

- `work_YYYYMMDDHH.taxis_repaired.nc`

`DRY_RUN=1` prints planned actions without writing files.

### 2) Remap + Concatenate Work Files

`script/coarse_grain_and_concat_work.sh`:

- Collects canonical per-date work files (`work_YYYYMMDDHH.nc`).
- Preferentially substitutes repaired counterparts if present.
- Remaps each file with CDO (default `remapcon,r360x180`).
- Concatenates remapped files with `ncrcat` into a date-range product:
  - `work_START_END.nc`

### 3) Concatenate Histograms

`script/concat_histograms.sh`:

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
2. Prepare date list in `launcher/list.txt`.
3. Run launcher (Slurm array or serial fallback) to generate per-date `work_*.nc` and `hist_*.nc`.
4. Run `script/repair_work_time_axis.sh` if time metadata issues are suspected.
5. Run `script/coarse_grain_and_concat_work.sh` to produce concatenated `work_START_END.nc`.
6. Optionally run `script/concat_histograms.sh` for histogram products.
7. Run MATLAB regional scripts or the all-regions summary script.

## Environment Notes

- Scripts rely on environment modules (Intel, HDF5, NetCDF, CDO, NCO).
- Paths in launchers/scripts are cluster-specific and may need adjustment per system.
- `setup_matlab_proxy.sh` provides one workflow for MATLAB proxy usage.

## Operational Caveats

- Missing-step counts generally indicate missing input dates, but can also reflect malformed metadata if pre-processing is skipped.
- Per-date files with empty time dimensions or invalid timestamps should be repaired or regenerated before concatenation.
- Keep date lists and source directory selection consistent between control and warming pipelines.

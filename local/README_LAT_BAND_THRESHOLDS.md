# Latitude-Band-Specific Precipitation Thresholds

## Overview

This workflow computes precipitation rate thresholds independently for 18 latitude bands (10-degree wide, covering -90 to +90) and applies them per location when computing work/lift.

Current implementation details:
- Threshold masking is latitude-band specific.
- `work` and `lift` are written as daily conditional means over threshold-exceedance events.
- `event_count` is written as daily count of threshold exceedances.
- Output is a single NetCDF file containing `work`, `lift`, and `event_count`.

## Components

### 1. `compute_prate_thresholds_by_lat_band.f90`

**Purpose**: Compute precipitation thresholds independently for each 10-degree latitude band.

**Input**:
- Configuration namelist with:
  - `history_root_part1`, `history_root_part2`: paths to precipitation data directories
  - `date_list_file_part1`, `date_list_file_part2`: date lists for data files
  - `output_file`: output ASCII file path (default: `thresholds_control_by_lat_band.txt`)

**Output**:
- ASCII file with 18 rows (one per latitude band) × 11 columns (one per percentile)
- Format (example):
  ```
  # lat_band p0.50000 p0.75000 p0.90000 ... p0.99999
  1   3.16E-06 1.75E-05 7.07E-05 ... 7.92E-03
  2   4.21E-06 2.10E-05 8.45E-05 ... 8.91E-03
  ...
  18  2.84E-06 1.52E-05 6.12E-05 ... 6.44E-03
  ```

**Latitude bands**:
- Band 1: [-90, -80]
- Band 2: [-80, -70]
- ...
- Band 18: [80, 90]

**Key logic**:
- Area-weighted histogram accumulation per latitude band (matching `compute_prate_thresholds.f90` logic)
- Log-log interpolation in survival-function space (avoiding numerical issues at CDF tails)
- Independent CDF computation per band

### 2. `compute_work_async_prate_threshold_by_lat_band.f90`

**Purpose**: Compute work/lift with latitude-band-specific precipitation masking and write daily conditional means plus daily event counts.

**Input**:
- Configuration namelist with:
  - `path_dz`
  - `path_temp`
  - `path_omega`
  - `path_qv`, `path_qw`, `path_qr`, `path_qi`, `path_qs`, `path_qg`
  - `path_omt`, `path_omqv`, `path_omqw`, `path_omqr`, `path_omqi`, `path_omqs`, `path_omqg`
  - `path_pr`
  - `path_work_out` (single output NetCDF path)
  - `path_thresholds`: ASCII threshold file from `compute_prate_thresholds_by_lat_band`

**Output**:
- One NetCDF file with dimensions `(percentile, lon, lat, time)`.
- Variables:
  - `work`: daily conditional mean over exceedance events; `_FillValue=-9999` when `event_count=0`.
  - `lift`: daily conditional mean over exceedance events; `_FillValue=-9999` when `event_count=0`.
  - `event_count`: daily count of threshold exceedance events.
- Global attributes for threshold traceability:
  - `bandNN_p..._threshold` for all 18x11 (for example `band01_p99_threshold`).

**Key logic**:
- For each grid point:
  1. Determine latitude band based on grid-point latitude
  2. Compute work and lift using all pressure levels
  3. Apply latitude-band-specific thresholds
  4. Accumulate masked sums and event counts per day
  5. Compute conditional means as `sum(masked_values)/event_count` when `event_count > 0`
- Parallel I/O with buffering (same async structure as `compute_work_async_prate_threshold.f90`)

**Latitude band determination**:
```fortran
ilat = int((lat_value + 90.0d0) / 10.0d0) + 1
ilat = max(1, min(nlat_bands, ilat))
```

## Usage Example

### 1. Compute thresholds:
```bash
./bin/compute_prate_thresholds_by_lat_band config_thresholds.nml
# Produces: output/thresholds/thresholds_control_by_lat_band.txt
```

Configuration file (`config_thresholds.nml`):
```namelist
&config
  history_root_part1 = '/scratch/cimes/GLOBALFV3/20191020.00Z.C3072.L79x2_pire/history'
  history_root_part2 = '/scratch/cimes/GLOBALFV3/stellar_run/processed_new/20191020.00Z.C3072.L79x2_pire/pp'
  date_list_file_part1 = 'launcher/list/list_control_part1.txt'
  date_list_file_part2 = 'launcher/list/list_control_part2.txt'
  output_file = 'output/thresholds/thresholds_control_by_lat_band.txt'
/
```

### 2. Compute work/lift:
```bash
./bin/compute_work_async_prate_threshold_by_lat_band config_work.nml
# Produces: one NetCDF with work, lift, and event_count
```

### 3. Coarse-grain and concatenate precipitation (lat-band thresholded):
```bash
# Control simulation
SIMULATION=control_prate_thresholded_by_lat_band \
  ./script/postprocessing/coarse_grain_and_concat_precip.sh

# Warming simulation
SIMULATION=warming_prate_thresholded_by_lat_band \
  ./script/postprocessing/coarse_grain_and_concat_precip.sh
```

These modes use the lat-band threshold files automatically:
- `output/thresholds/thresholds_control_by_lat_band.txt`
- `output/thresholds/thresholds_warming_by_lat_band.txt`

Outputs are written to:
- `/scratch/gpfs/mbolot/results/GLOBALFV3/precip_coarse_C3072_360x180_prate_thresholded_by_lat_band`
- `/scratch/gpfs/mbolot/results/GLOBALFV3/precip_coarse_C3072_360x180_PLUS_4K_CO2_1270ppmv_prate_thresholded_by_lat_band`

Configuration file (`config_work.nml`):
```namelist
&config
  path_dz = '/path/to/DZ.nc'
  path_temp = '/path/to/temp.nc'
  path_omega = '/path/to/omega.nc'
  path_qv = '/path/to/qv.nc'
  path_qw = '/path/to/qw.nc'
  path_qr = '/path/to/qr.nc'
  path_qi = '/path/to/qi.nc'
  path_qs = '/path/to/qs.nc'
  path_qg = '/path/to/qg.nc'
  path_omt = '/path/to/omt.nc'
  path_omqv = '/path/to/omqv.nc'
  path_omqw = '/path/to/omqw.nc'
  path_omqr = '/path/to/omqr.nc'
  path_omqi = '/path/to/omqi.nc'
  path_omqs = '/path/to/omqs.nc'
  path_omqg = '/path/to/omqg.nc'
  path_pr = '/path/to/prate.nc'
  path_work_out = 'work_lift_prate_band.nc'
  path_thresholds = 'output/thresholds/thresholds_control_by_lat_band.txt'
/
```

## Makefile Targets

```bash
make clean
make all  # Builds all targets including new ones

# Individual targets:
make bin/compute_prate_thresholds_by_lat_band
make bin/compute_work_async_prate_threshold_by_lat_band
```

## Performance Considerations

1. **Threshold computation**: ~18× more independent CDF calculations (one per band), but each histogram is smaller (fewer grid points per band). Overall runtime similar to single-domain version.

2. **Work/lift computation**: Same async I/O structure; per-grid-point latitude band lookup is O(1), minimal overhead.

3. **Output size**: Same as non-banded version (11 percentiles × domain grid). Thresholds stored as global attributes (~18KB per file for 18×11 values).

## Verification

To verify correctness:

1. **Threshold file structure**:
   ```bash
   wc -l thresholds_control_by_lat_band.txt
   # Should show 19 lines (1 header + 18 data)
   head -3 thresholds_control_by_lat_band.txt
   # Should show header and first 2 bands
   ```

2. **Output attributes** (using `ncdump`):
   ```bash
  ncdump -h work_lift_prate_band.nc | grep band..
   # Should show 18×11=198 global attributes
   ```

3. **Output variables and stats metadata**:
  ```bash
  ncdump -h work_lift_prate_band.nc | egrep 'double (work|lift|event_count)|time_stat|time_aggregation'
  ```
  Expected:
  - `work` and `lift` have conditional-mean metadata.
  - `event_count` has daily-count metadata.

## Debugging

- **Threshold file parsing errors**: Check ASCII file format (18 rows, first column is band index, remaining 11 columns are thresholds)
- **Grid latitude mismatch**: Verify that `grid_yt_coarse` is present and consistent with the input files
- **Band index out of range**: Ensure latitude values are in [-90, 90]; clamping is applied but may indicate data issue
- **No-event regions**: `work/lift` are set to `_FillValue` where `event_count=0`; this is expected behavior for conditional means

## Related Programs

- `compute_prate_thresholds.f90`: Single-domain threshold computation (original)
- `compute_work_async_prate_threshold.f90`: Single-domain work/lift computation (original)

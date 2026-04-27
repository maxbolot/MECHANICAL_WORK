# GLOBALFV3 Data Layout Migration Plan

Date: 2026-04-25
Owner: mbolot
Target window: about 2 weeks

Warming dataset note: in this project, "warming" corresponds to the simulation with +4K surface temperature and CO2 set to 1270 ppmv.

## Objectives

1. Migrate data root from:
   - `/scratch/gpfs/mbolot/results/GLOBALFV3`
   to:
   - `/scratch/gpfs/mbolot/processed_data/GLOBALFV3`
2. Reorganize tree under `GLOBALFV3` into:
   - `C3072_1440x720/`
   - `C3072_360x180/`
   - `histogram/`
3. Ensure level-2 organization under each level-1 directory:
   - `work/`
   - `precip/`
4. Add region-aware organization for histogram products so multi-region diagnostics (including latitude-band partitions) are first-class in the layout.
4. Safely migrate all hardcoded paths across the project with rollback options.

## Target Layout

```text
/scratch/gpfs/mbolot/processed_data/GLOBALFV3/
  C3072_1440x720/
    work/
      control/
      warming/
      control_prate_thresholded/
      warming_prate_thresholded/
      control_prate_thresholded_by_lat_band/
      warming_prate_thresholded_by_lat_band/
      control_ea_to_025/               (optional legacy)
      test/                            (optional legacy)
    precip/
      (if produced at this resolution in future)

  C3072_360x180/
    work/
      control/
      warming/
      control_prate_thresholded/
      warming_prate_thresholded/
      control_prate_thresholded_by_lat_band/
      warming_prate_thresholded_by_lat_band/
    precip/
      control/
      warming/
      control_prate_thresholded/
      warming_prate_thresholded/
      control_prate_thresholded_by_lat_band/
      warming_prate_thresholded_by_lat_band/

  histogram/
    work/
      control/
        global/
          all/
        lat_band/
          tropics/
          nh_midlat/
          sh_midlat/
      warming/
        global/
          all/
        lat_band/
          tropics/
          nh_midlat/
          sh_midlat/
    precip/
      (empty placeholder for schema consistency)
```

Notes for histogram region structure:

1. `region_family` examples: `global`, `lat_band`, `custom_box`, `mask`.
2. `region_id` identifies one concrete region in that family (for example `all`, `tropics`, `band_30S_0`).
3. Keep one histogram file per scenario + region family + region id + date range.

## Old to New Mapping Rules

### Resolution and data type mapping

1. `work_coarse_C3072_1440x720*` -> `C3072_1440x720/work/<scenario>/`
2. `work_coarse_C3072_360x180*` -> `C3072_360x180/work/<scenario>/`
3. `precip_coarse_C3072_360x180*` -> `C3072_360x180/precip/<scenario>/`
4. `work_histograms*` -> `histogram/work/<scenario>/`
  - region-aware target: `histogram/work/<scenario>/<region_family>/<region_id>/`

### Scenario suffix mapping

1. no special suffix -> `control`
2. `_PLUS_4K_CO2_1270ppmv` -> `warming`
3. `_prate_thresholded` -> `control_prate_thresholded`
4. `_PLUS_4K_CO2_1270ppmv_prate_thresholded` -> `warming_prate_thresholded`
5. `_prate_thresholded_by_lat_band` -> `control_prate_thresholded_by_lat_band`
6. `_PLUS_4K_CO2_1270ppmv_prate_thresholded_by_lat_band` -> `warming_prate_thresholded_by_lat_band`

### Region mapping for histogram products

1. Existing tropical histogram products map to:
  - `<region_family>=lat_band`
  - `<region_id>=tropics`
2. Existing global histogram products (if/when present) map to:
  - `<region_family>=global`
  - `<region_id>=all`
3. Future regional products should not create new scenario names; they should use new region-family/region-id directories.

## Migration Phases

## Phase 0: Freeze + inventory (Day 1)

1. Stop launchers/postprocessing jobs writing to old root.
2. Capture inventory for each old folder:
   - file count
   - total bytes
   - representative netCDF metadata snapshots (`ncdump -h`, `cdo sinfov`)
3. Save inventory logs under project (suggested path: `output/migration_logs/`).

## Phase 1: Create destination skeleton (Day 1)

1. Create full destination directories under `/scratch/gpfs/mbolot/processed_data/GLOBALFV3`.
2. Set permissions/ownership to match old tree.
3. Pre-create region-family and region-id folders under `histogram/work/<scenario>/` for currently supported regions (at minimum `lat_band/tropics`).

## Phase 2: Data copy + verify (Days 2-7)

1. Copy in batches using `rsync -aH --info=progress2`.
2. Copy largest/high-risk folders first:
   - thresholded/by-lat-band work at 1440x720
   - thresholded work/precip at 360x180
3. After each batch, verify:
   - file count parity
   - byte parity
   - metadata spot checks on representative files
4. Re-run rsync for each batch until no differences.

Suggested copy pattern:

```bash
rsync -aH --info=progress2 \
  /scratch/gpfs/mbolot/results/GLOBALFV3/work_coarse_C3072_360x180_prate_thresholded/ \
  /scratch/gpfs/mbolot/processed_data/GLOBALFV3/C3072_360x180/work/control_prate_thresholded/
```

## Phase 3: Compatibility bridge + cutover (Days 8-9)

1. Keep old tree intact (read-only if possible).
2. Optionally create temporary compatibility symlink:

```bash
# Only after moving old root out of the way or replacing it intentionally
ln -s /scratch/gpfs/mbolot/processed_data/GLOBALFV3 /scratch/gpfs/mbolot/results/GLOBALFV3
```

3. Switch default paths in code to new root.
4. Run smoke tests before any deletion.

## Phase 4: Code path migration (Days 8-12)

### Priority 1: launchers and postprocessing scripts

- `launcher/compute_work_async_array.sh`
- `launcher/compute_work_async_array_inner_parallel.sh`
- `launcher/compute_work_async_noslurm.sh`
- `launcher/compute_work_async_test_noslurm.sh`
- `launcher/compute_work_async_prate_threshold_array.sh`
- `launcher/compute_work_async_prate_threshold_noslurm.sh`
- `launcher/compute_work_async_prate_threshold_by_lat_band_array.sh`
- `launcher/compute_work_async_prate_threshold_by_lat_band_noslurm.sh`
- `script/postprocessing/coarse_grain_and_concat_work.sh`
- `script/postprocessing/coarse_grain_and_concat_precip.sh`
- `script/postprocessing/check_work_time_axis.sh`
- `script/postprocessing/concat_histograms.sh`
- `script/postprocessing/check_hist_time_axis.sh`

Histogram-specific roadmap note:

- Update histogram concat/check scripts to resolve scenario + region family + region id paths, not scenario-only paths.

### Priority 2: MATLAB/Python defaults and presets

- `matlab/presets/scenario_control.m`
- `matlab/presets/scenario_warming.m` (or keep current filename and rename internals to warming)
- `matlab/analysis/run_work_lift_all_regions_analysis.m`
- `matlab/analysis/run_histogram_precip_work_analysis.m`
- `matlab/entries/compute_histogram_precip_work_analysis_control.m`
- `matlab/entries/compute_histogram_precip_work_analysis_warming.m` (or keep current filename and rename internals to warming)
- `python/translated_from_matlab/presets/scenarios.py`
- `python/translated_from_matlab/analysis/work_lift_all_regions.py`

### Priority 3: publication and notebook references

- `publication/paper_dissipation/figures_src/*.m`
- `publication/paper_dissipation/manifests/*.yaml`
- notebooks that contain hardcoded printed paths

### Generated artifacts note

`output/namelists/*.nml` contains many old-path references and should be regenerated after changing source defaults. Do not hand-edit these generated files.

## Path Management Refactor (Recommended)

To avoid future churn, centralize the data root path:

1. Bash: use `GLOBALFV3_DATA_ROOT` env var with default:
   - `/scratch/gpfs/mbolot/processed_data/GLOBALFV3`
2. MATLAB: helper function reading `getenv('GLOBALFV3_DATA_ROOT')`.
3. Python: helper constant/function reading `os.getenv('GLOBALFV3_DATA_ROOT', default)`.

Then construct all scenario paths from this root.

## Verification Gates

## Gate A: Data integrity

1. Source vs destination file count parity per mapped folder.
2. Source vs destination byte size parity.
3. NetCDF spot checks (`ncdump -h`, `cdo sinfov`) for representative files.
4. Histogram metadata checks include region identity attributes (for example `region_family`, `region_id`, and region bounds or equivalent region descriptors).

## Gate B: Code integrity

1. `grep` scan for old root in active source directories returns zero intentional references.
2. Scripts pass syntax checks (`bash -n`).
3. MATLAB/Python analyses run with new defaults.

## Gate C: Functional smoke tests

1. Control standard workflow (work + precip + histogram).
2. Warming standard workflow.
3. Control thresholded workflow.
4. Warming thresholded-by-lat-band workflow.
5. Histogram retrieval and analysis for at least two distinct regions (for example `lat_band/tropics` and `global/all` when available).

## Rollback Plan

1. Keep old tree available and unchanged during validation window.
2. If failures occur, switch defaults back to old root (or use compatibility symlink).
3. Re-run failed step, fix, and retry cutover.
4. Do not delete old root until all gates pass for several days of normal usage.

## Two-Week Checklist

### Week 1

1. Freeze writers and collect inventories.
2. Create destination tree.
3. Copy data in batches with verification after each batch.
4. Resolve any copy mismatches.

### Week 2

1. Introduce centralized root variable in scripts/MATLAB/Python.
2. Migrate hardcoded paths by priority.
3. Regenerate namelists.
4. Run smoke tests and publication/notebook checks.
5. Enable temporary compatibility symlink if needed.
6. Final acceptance and deprecate old references.

## Acceptance Criteria

1. All expected data available under `/scratch/gpfs/mbolot/processed_data/GLOBALFV3` in the new layout.
2. Project scripts/analyses run without relying on `/scratch/gpfs/mbolot/results/GLOBALFV3`.
3. Validation logs and smoke tests are clean.
4. Rollback procedure remains available until final sign-off.

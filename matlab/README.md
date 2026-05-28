# MATLAB Analyses

This folder contains MATLAB entry scripts, analysis orchestrators, reusable utilities, and scenario/region presets.

## Directory Layout

- `entries/`: runnable entry-point scripts
- `analysis/`: orchestration and plotting routines
- `lib/`: shared helpers and utilities
- `presets/`: scenario and region presets

## GLOBALFV3 Data Root Auto-Detection

MATLAB scripts resolve the GLOBALFV3 root through:

- `lib/get_globalfv3_data_root.m`

Resolution order:

1. Use `GLOBALFV3_DATA_ROOT` if set.
2. Otherwise use OS defaults:
   - Windows (`ispc`): `C:\climate_processed_data\GLOBALFV3`
   - Unix/Linux (`isunix`): `/scratch/gpfs/mbolot/results/GLOBALFV3`

Scenario presets and analysis wrappers build file paths with `fullfile(...)`, so only the root changes by environment while filenames remain unchanged.

## Recommended Overrides

HPC shell:

```bash
export GLOBALFV3_DATA_ROOT=/scratch/gpfs/mbolot/results/GLOBALFV3
```

Windows MATLAB session:

```matlab
setenv('GLOBALFV3_DATA_ROOT', 'C:\climate_processed_data\GLOBALFV3');
```

WSL/Linux MATLAB session reading Windows files:

```bash
export GLOBALFV3_DATA_ROOT=/mnt/c/climate_processed_data/GLOBALFV3
```

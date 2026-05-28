# Python Translated Analyses

This package contains Python translations of selected MATLAB analysis workflows.

## Package Layout

- `analysis/`: orchestration and plotting routines
- `presets/`: scenario and region presets
- `lib/`: shared helpers and utilities

## GLOBALFV3 Data Root Auto-Detection

Translated workflows resolve the GLOBALFV3 root with:

- `translated_from_matlab.lib.paths.get_globalfv3_data_root()`

Resolution order:

1. Use `GLOBALFV3_DATA_ROOT` if set.
2. Otherwise use platform defaults:
   - Windows: `C:\climate_processed_data\GLOBALFV3`
   - WSL: `/mnt/c/climate_processed_data/GLOBALFV3`
   - Linux/HPC: `/scratch/gpfs/mbolot/results/GLOBALFV3`

Scenario presets and analysis defaults build filenames from this root, so only the root changes by environment.

## Recommended Overrides

HPC shell:

```bash
export GLOBALFV3_DATA_ROOT=/scratch/gpfs/mbolot/results/GLOBALFV3
```

WSL shell:

```bash
export GLOBALFV3_DATA_ROOT=/mnt/c/climate_processed_data/GLOBALFV3
```

Windows PowerShell:

```powershell
$env:GLOBALFV3_DATA_ROOT = 'C:\climate_processed_data\GLOBALFV3'
```

#!/usr/bin/env bash
# pull_data_from_HPC.sh - Pull processed NetCDF files from the HPC

HPC_HOST="mbolot@stellar-vis1.princeton.edu"
HPC_DIR="/scratch/gpfs/mbolot/results/GLOBALFV3/"
LOCAL_DIR="/mnt/c/climate_processed_data/GLOBALFV3/"

echo "Syncing processed NetCDF files from HPC to local environment..."

# Create the local directory if it doesn't exist yet
mkdir -p "${LOCAL_DIR}"

rsync -avP \
  --include="*.nc" \
  --exclude="*1440x720*" \
  "${HPC_HOST}:${HPC_DIR}" "${LOCAL_DIR}"

echo -e "\nSync complete! Files are natively accessible on Windows at C:\\climate_processed_data\\GLOBALFV3\\"
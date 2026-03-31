#!/bin/bash

set -euo pipefail

LOG_FILE="/home/mbolot/Projects/MECHANICAL_WORK/.vscode/ifx-lint.log"

{
  echo "=== $(date '+%Y-%m-%d %H:%M:%S') ==="
  echo "PWD=$PWD"
  echo "ARGS=$*"
} >> "$LOG_FILE"

source /usr/share/Modules/init/bash >/dev/null 2>&1 || source /etc/profile.d/modules.sh >/dev/null 2>&1
module load intel-oneapi/2024.2 hdf5/oneapi-2024.2/1.14.4 netcdf/oneapi-2024.2/hdf5-1.14.4/4.9.2 >/dev/null 2>&1

exec /opt/intel/oneapi/compiler/2024.2/bin/ifx "$@"
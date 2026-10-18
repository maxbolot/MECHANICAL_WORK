#!/usr/bin/env bash
# env.sh - Sourced environment configuration for MECHANICAL_WORK
# Works seamlessly on both Princeton HPC (Della/Stellar/Spirit) and Thinkpad P51 workstation.

# Locate project root
_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PROJECT_ROOT="${PROJECT_ROOT:-$_SCRIPT_DIR}"

# 1. Load local .env or .env.local if present
if [[ -f "$PROJECT_ROOT/.env.local" ]]; then
    # shellcheck disable=SC1091
    source "$PROJECT_ROOT/.env.local"
elif [[ -f "$PROJECT_ROOT/.env" ]]; then
    # shellcheck disable=SC1091
    source "$PROJECT_ROOT/.env"
fi

# 2. Host detection
_HOSTNAME=$(hostname 2>/dev/null || echo "unknown")

if [[ "$_HOSTNAME" =~ (della|stellar|spirit|oondemand) ]]; then
    # On Princeton HPC
    export HPC_ENV="true"
    export GLOBALFV3_DATA_ROOT="${GLOBALFV3_DATA_ROOT:-/scratch/gpfs/mbolot/results/GLOBALFV3}"
    export GLOBALFV3_RAW_ROOT="${GLOBALFV3_RAW_ROOT:-/scratch/cimes/GLOBALFV3}"
else
    # On Workstation (Thinkpad P51 or other Linux workstation)
    export HPC_ENV="false"
    export GLOBALFV3_DATA_ROOT="${GLOBALFV3_DATA_ROOT:-/media/max/data/climate_work/processed}"
    export GLOBALFV3_RAW_ROOT="${GLOBALFV3_RAW_ROOT:-/media/max/data/climate_work/raw}"
    export HPC_HOST="${HPC_HOST:-stellar-vis1.princeton.edu}"
    export HPC_USER="${HPC_USER:-mbolot}"

    # Auto-initialize Intel oneAPI if compiler not in PATH and setvars exists
    if ! command -v ifx >/dev/null 2>&1; then
        _ONEAPI_SETVARS="${ONEAPI_SETVARS:-/opt/intel/oneapi/setvars.sh}"
        if [[ -f "$_ONEAPI_SETVARS" ]]; then
            _RESTORE_E=0
            [[ $- == *e* ]] && _RESTORE_E=1
            _RESTORE_U=0
            [[ $- == *u* ]] && _RESTORE_U=1
            set +eu
            # shellcheck disable=SC1091
            source "$_ONEAPI_SETVARS" --force >/dev/null 2>&1 || true
            [[ $_RESTORE_E -eq 1 ]] && set -e
            [[ $_RESTORE_U -eq 1 ]] && set -u
        fi
    fi
    export FC="${FC:-ifx}"

    # Support external libraries path (default ~/local)
    export EXT_INC_DIR="${EXT_INC_DIR:-$HOME/local/include}"
    export EXT_LIB_DIR="${EXT_LIB_DIR:-$HOME/local/lib}"
    if [[ -d "$EXT_LIB_DIR" ]]; then
        export LD_LIBRARY_PATH="$EXT_LIB_DIR:${LD_LIBRARY_PATH:-}"
    fi
fi

# Add project python directory to PYTHONPATH
if [[ -d "$PROJECT_ROOT/python" ]]; then
    export PYTHONPATH="$PROJECT_ROOT/python:${PYTHONPATH:-}"
fi


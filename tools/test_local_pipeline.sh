#!/usr/bin/env bash
# test_local_pipeline.sh - Smoke test runner for local sandbox & AI iteration
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Source environment
# shellcheck disable=SC1091
source "$PROJECT_ROOT/env.sh"

# Re-enforce strict error handling
set -euo pipefail

echo "========================================================================"
echo " MECHANICAL_WORK: Local Sandbox & AI Pipeline Validation"
echo "========================================================================"
echo "Project Root:         $PROJECT_ROOT"
echo "Processed Data Root:  $GLOBALFV3_DATA_ROOT"
echo "Raw Data Root:        $GLOBALFV3_RAW_ROOT"
echo "Fortran Compiler:     $(command -v ifx || command -v gfortran || echo 'NOT FOUND')"
echo "CDO:                  $(command -v cdo || echo 'NOT FOUND (activate conda geo)')"
echo "NCO (ncks):           $(command -v ncks || echo 'NOT FOUND (activate conda geo)')"
echo "========================================================================"

SAMPLE_DATE="2020080500"
SAMPLE_RAW_DIR="$GLOBALFV3_RAW_ROOT/20191020.00Z.C3072.L79x2_pire/history/$SAMPLE_DATE"
COMPUTE_BIN="$PROJECT_ROOT/local/bin/compute_work_async"

# 1. Check Toolchain
if ! command -v ifx >/dev/null 2>&1 && ! command -v gfortran >/dev/null 2>&1; then
    echo "Warning: Fortran compiler not found in PATH."
    if [[ -f /opt/intel/oneapi/setvars.sh ]]; then
        echo "Found Intel oneAPI at /opt/intel/oneapi/setvars.sh. Sourcing..."
        _RESTORE_E=0
        [[ $- == *e* ]] && _RESTORE_E=1
        set +eu
        # shellcheck disable=SC1091
        source /opt/intel/oneapi/setvars.sh --force >/dev/null 2>&1 || true
        [[ $_RESTORE_E -eq 1 ]] && set -e
    fi
fi

# 2. Check Support Libraries
LIB_INTERP="$HOME/local/lib/libinterp.a"
if [[ ! -f "$LIB_INTERP" ]]; then
    echo -e "\n[Notice] Support library $LIB_INTERP not found locally."
    echo "To pull compiled libraries and headers from HPC, run:"
    echo "  ./tools/sync_hpc_data.sh --profile vendor-libs"
fi

# 3. Check Raw Sample Data
if [[ ! -d "$SAMPLE_RAW_DIR" ]]; then
    echo -e "\n[Notice] Benchmark raw data for $SAMPLE_DATE not found at:"
    echo "  $SAMPLE_RAW_DIR"
    echo "To pull benchmark sample data from HPC, run:"
    echo "  ./tools/sync_hpc_data.sh --profile sample-raw --sim control --date $SAMPLE_DATE"
    echo ""
    echo "To run with simulated dry-run or verify other parts of the pipeline, ensure data is synced first."
    exit 0
fi

# 4. Compile binary if needed
if [[ ! -x "$COMPUTE_BIN" ]]; then
    echo -e "\nCompiling $COMPUTE_BIN..."
    LIBRARY_PATH="" make -C "$PROJECT_ROOT/local" bin/compute_work_async
fi

# 5. Run test execution
echo -e "\nRunning test computation on sample date $SAMPLE_DATE..."
TEST_RUN_ID="local_test_$(date +%Y%m%d_%H%M%S)"
DATE="$SAMPLE_DATE" \
RUN_ID_ROOT="$TEST_RUN_ID" \
"$PROJECT_ROOT/launcher/compute_work_async_test_noslurm.sh"

echo -e "\n========================================================================"
echo " Pipeline smoke test completed successfully!"
echo " Output files written under: $GLOBALFV3_DATA_ROOT/work_coarse_C3072_1440x720_test/$TEST_RUN_ID"
echo "========================================================================"


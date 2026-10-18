#!/usr/bin/env bash
# pull_data_from_HPC.sh - Legacy wrapper; forwards to tools/sync_hpc_data.sh
echo "Notice: pull_data_from_HPC.sh is now replaced by tools/sync_hpc_data.sh."
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "$SCRIPT_DIR/tools/sync_hpc_data.sh" --profile analysis "$@"
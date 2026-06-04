#!/bin/bash
# =============================================================================
#  transRNA pipeline v1.2 — login-node launcher
#
#  Usage:
#    bash run.sh           # fresh run from scratch
#    bash run.sh --resume  # resume a previous run (reuse cached steps)
# =============================================================================

set -euo pipefail

# ── Edit these two if needed ──────────────────────────────────────────────────
PARTITION="icelake"
PROJECT="MAORI-SL2-CPU"
# ─────────────────────────────────────────────────────────────────────────────

# Parse optional --resume flag
RESUME_FLAG=""
for arg in "$@"; do
    if [[ "$arg" == "--resume" ]]; then
        RESUME_FLAG="-resume"
    fi
done

# Suppress tput warnings from Nextflow in non-interactive sessions
export TERM=xterm

PIPELINE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
mkdir -p "${PIPELINE_DIR}/logs" "${PIPELINE_DIR}/pipeline_reports"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG="${PIPELINE_DIR}/logs/nextflow_${TIMESTAMP}.log"

echo "============================================================"
echo "  transRNA pipeline v1.2"
echo "  $(date)"
echo "  Partition : ${PARTITION}"
echo "  Project   : ${PROJECT}"
echo "  Mode      : ${RESUME_FLAG:-(fresh run)}"
echo "  Nextflow  : $(nextflow -version 2>&1 | head -1)"
echo "  Log       : ${LOG}"
echo "============================================================"

nextflow run "${PIPELINE_DIR}/main.nf" \
    -profile cambridge \
    -params-file "${PIPELINE_DIR}/params.yml" \
    --partition "${PARTITION}" \
    --project   "${PROJECT}" \
    ${RESUME_FLAG} \
    2>&1 | tee "${LOG}"

echo "============================================================"
echo "  Finished: $(date)"
echo "  Log: ${LOG}"
echo "============================================================"

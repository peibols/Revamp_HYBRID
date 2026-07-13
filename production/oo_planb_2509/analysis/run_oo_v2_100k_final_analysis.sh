#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${ROOT:-/raid5/data/yjlee/hybrid_dev}"
FIRST_WORK="${FIRST_WORK:-${ROOT}/test/oo5360_v2_500hydro_50k_20260711}"
CONTINUATION_WORK="${CONTINUATION_WORK:-${ROOT}/test/oo5360_v2_500hydro_cont50k_to100k_20260712}"

export ROOT
export WORK="${WORK:-${FIRST_WORK}}"
export LOCAL_EOS="${LOCAL_EOS:-${FIRST_WORK}/eos_snapshot}"
export ADDITIONAL_LOCAL_EOS="${ADDITIONAL_LOCAL_EOS:-${CONTINUATION_WORK}/eos_snapshot}"
export STRICT_MARKER="${STRICT_MARKER:-${FIRST_WORK}/v2_50k_strict_complete.txt}"
export ADDITIONAL_STRICT_MARKER="${ADDITIONAL_STRICT_MARKER:-${CONTINUATION_WORK}/v2_50k_strict_complete.txt}"
export TASK_MANIFEST="${TASK_MANIFEST:-${CONTINUATION_WORK}/hydro_prepared/aa_task_manifest_combined_100k.tsv}"
export OUT="${OUT:-${FIRST_WORK}/final_analysis_100k}"
export BUILD_DIR="${BUILD_DIR:-${ROOT}/test/oo5360_v2_100k_final_root_build}"
export TARGET=100000
export TAG=100k
export PREFIX=oo5360_v2_100k
export FINAL_MARKER="${OUT}/v2_100k_final_analysis_complete.txt"

exec "${SCRIPT_DIR}/run_oo_v2_final_analysis.sh"

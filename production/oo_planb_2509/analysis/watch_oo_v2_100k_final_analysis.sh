#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${ROOT:-/raid5/data/yjlee/hybrid_dev}"
FIRST_WORK="${FIRST_WORK:-${ROOT}/test/oo5360_v2_500hydro_50k_20260711}"
CONTINUATION_WORK="${CONTINUATION_WORK:-${ROOT}/test/oo5360_v2_500hydro_cont50k_to100k_20260712}"
POLL_SECONDS="${POLL_SECONDS:-300}"

FIRST_MARKER="${FIRST_WORK}/v2_50k_strict_complete.txt"
CONTINUATION_MARKER="${CONTINUATION_WORK}/v2_50k_strict_complete.txt"
FINAL_MARKER="${FIRST_WORK}/final_analysis_100k/v2_100k_final_analysis_complete.txt"

marker_is_complete() {
    local marker="$1"
    [[ -s "${marker}" ]] && grep -qx 'accepted_pairs=50000' "${marker}"
}

if [[ -s "${FINAL_MARKER}" ]] && grep -qx 'status=PASS' "${FINAL_MARKER}"; then
    printf '%s final 100k analysis already passed: %s\n' \
        "$(date --iso-8601=seconds)" "${FINAL_MARKER}"
    exit 0
fi

printf '%s waiting for both audited 50k strict markers\n' \
    "$(date --iso-8601=seconds)"
while ! marker_is_complete "${FIRST_MARKER}" || \
      ! marker_is_complete "${CONTINUATION_MARKER}"; do
    sleep "${POLL_SECONDS}"
done

printf '%s both strict markers passed; starting combined 100k analysis\n' \
    "$(date --iso-8601=seconds)"
exec "${SCRIPT_DIR}/run_oo_v2_100k_final_analysis.sh"

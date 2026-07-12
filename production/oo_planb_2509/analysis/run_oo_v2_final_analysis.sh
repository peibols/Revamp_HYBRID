#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
SOURCE="${SOURCE:-$(cd -- "${SCRIPT_DIR}/../../.." && pwd)}"
ROOT="${ROOT:-/raid5/data/yjlee/hybrid_dev}"
WORK="${WORK:-${ROOT}/test/oo5360_v2_500hydro_50k_20260711}"
LOCAL_EOS="${LOCAL_EOS:-${WORK}/eos_snapshot}"
PP_LOCAL_EOS="${PP_LOCAL_EOS:-${ROOT}/test/tmp_oo_10k_prehydro_raa_20260709/local_eos}"
TASK_MANIFEST="${TASK_MANIFEST:-${WORK}/hydro_prepared/aa_task_manifest.tsv}"
STRICT_MARKER="${STRICT_MARKER:-${WORK}/v2_50k_strict_complete.txt}"
OUT="${OUT:-${WORK}/final_analysis}"
BUILD_DIR="${BUILD_DIR:-${ROOT}/test/tmp_oo_v2_final_root_build}"
OVERWRITE="${OVERWRITE:-false}"

TARGET=50000
RAA_OUT="${OUT}/hadron_raa"
ROOT_OUT="${OUT}/root/oo5360_c0_5_planb_v2_50k.root"
JET20_OUT="${OUT}/jets/inclusive_pt20"
JET30_OUT="${OUT}/jets/inclusive_pt30"
JET_SLICE_OUT="${OUT}/jets/pt_slices"
FINAL_MARKER="${OUT}/v2_final_analysis_complete.txt"

test -s "${STRICT_MARKER}"
test -d "${LOCAL_EOS}/status/aa"
test -d "${LOCAL_EOS}/outputs/aa"
test -d "${PP_LOCAL_EOS}/status/pp"
test -d "${PP_LOCAL_EOS}/outputs/pp"
test -s "${TASK_MANIFEST}"

mkdir -p "${RAA_OUT}" "$(dirname -- "${ROOT_OUT}")" \
  "${JET20_OUT}" "${JET30_OUT}" "${JET_SLICE_OUT}" "${BUILD_DIR}"

python3 "${SCRIPT_DIR}/analyze_oo_prehydro_pair.py" \
  --local-eos "${LOCAL_EOS}" \
  --pp-local-eos "${PP_LOCAL_EOS}" \
  --out-dir "${RAA_OUT}" \
  --require-paired-aa \
  --aa-events-per-chunk 1 \
  --aa-task-manifest "${TASK_MANIFEST}" \
  --aa-task-limit "${TARGET}" \
  --require-complete-aa-prefix

converter_args=(
  --source "v2=${LOCAL_EOS}"
  --aa-task-manifest "${TASK_MANIFEST}"
  --output "${ROOT_OUT}"
  --build-dir "${BUILD_DIR}"
  --progress-every 500
)
case "${OVERWRITE}" in
  1|true|TRUE) converter_args+=(--overwrite) ;;
  0|false|FALSE) ;;
  *) echo "OVERWRITE must be true or false, got: ${OVERWRITE}" >&2; exit 1 ;;
esac
python3 "${SCRIPT_DIR}/convert_oo_paired_to_root.py" "${converter_args[@]}"

python3 "${SCRIPT_DIR}/plot_oo_jet_variables.py" \
  --input-root "${ROOT_OUT}" --out-dir "${JET20_OUT}" \
  --pt-min 20 --prefix oo5360_v2_50k_jet_variables_pt20
python3 "${SCRIPT_DIR}/plot_oo_jet_variables.py" \
  --input-root "${ROOT_OUT}" --out-dir "${JET30_OUT}" \
  --pt-min 30 --prefix oo5360_v2_50k_jet_variables_pt30
python3 "${SCRIPT_DIR}/plot_oo_jet_variables.py" \
  --input-root "${ROOT_OUT}" --out-dir "${JET_SLICE_OUT}" \
  --pt-min 20 --pt-max 30 --prefix oo5360_v2_50k_jet_variables_pt20to30
python3 "${SCRIPT_DIR}/plot_oo_jet_variables.py" \
  --input-root "${ROOT_OUT}" --out-dir "${JET_SLICE_OUT}" \
  --pt-min 30 --pt-max 50 --prefix oo5360_v2_50k_jet_variables_pt30to50
python3 "${SCRIPT_DIR}/plot_oo_jet_variables.py" \
  --input-root "${ROOT_OUT}" --out-dir "${JET_SLICE_OUT}" \
  --pt-min 50 --pt-max 80 --prefix oo5360_v2_50k_jet_variables_pt50to80
python3 "${SCRIPT_DIR}/plot_oo_jet_variables.py" \
  --input-root "${ROOT_OUT}" --out-dir "${JET_SLICE_OUT}" \
  --pt-min 80 --prefix oo5360_v2_50k_jet_variables_pt80plus

python3 "${SCRIPT_DIR}/summarize_oo_jet_pt_slices.py" \
  --slice-metadata \
    "${JET_SLICE_OUT}/oo5360_v2_50k_jet_variables_pt20to30_metadata.json" \
    "${JET_SLICE_OUT}/oo5360_v2_50k_jet_variables_pt30to50_metadata.json" \
    "${JET_SLICE_OUT}/oo5360_v2_50k_jet_variables_pt50to80_metadata.json" \
    "${JET_SLICE_OUT}/oo5360_v2_50k_jet_variables_pt80plus_metadata.json" \
  --inclusive-metadata \
    "${JET20_OUT}/oo5360_v2_50k_jet_variables_pt20_metadata.json" \
  --out-dir "${JET_SLICE_OUT}" \
  --prefix oo5360_v2_50k_jet_pt_slice_summary

python3 - \
  "${ROOT_OUT%.root}.summary.json" \
  "${JET_SLICE_OUT}/oo5360_v2_50k_jet_pt_slice_summary_validation.json" \
  "${TASK_MANIFEST}" \
  "${TARGET}" <<'PY'
import hashlib
import json
import math
from pathlib import Path
import sys

conversion_path, closure_path, manifest_path = map(Path, sys.argv[1:4])
target = int(sys.argv[4])
conversion = json.loads(conversion_path.read_text())
closure = json.loads(closure_path.read_text())
expected_manifest_sha = hashlib.sha256(manifest_path.read_bytes()).hexdigest()

if conversion["acceptedPairs"] != target:
    raise SystemExit(
        f"ROOT conversion accepted {conversion['acceptedPairs']} pairs, expected {target}"
    )
if conversion["rejectedArchives"] != 0 or conversion["skippedChunkRecords"] != 0:
    raise SystemExit(
        "ROOT conversion has rejected or skipped records: "
        f"{conversion['rejectedArchives']} rejected, "
        f"{conversion['skippedChunkRecords']} skipped"
    )
manifest = conversion.get("aaTaskManifest")
if manifest is None or manifest["rows"] != target:
    raise SystemExit("ROOT conversion did not record the complete task manifest")
if manifest["sha256"] != expected_manifest_sha:
    raise SystemExit("ROOT conversion task-manifest SHA256 does not match")
if closure.get("status") != "PASS":
    raise SystemExit("jet-pT slice closure did not pass")
if any(
    not math.isfinite(difference) or abs(difference) > 1.0e-12
    for radius in closure["crossSectionClosureDifferenceMb"].values()
    for difference in radius.values()
):
    raise SystemExit("jet-pT slice closure exceeds 1e-12 mb")
PY

source_commit="$(git -C "${SOURCE}" rev-parse HEAD)"
root_sha256="$(sha256sum "${ROOT_OUT}" | awk '{print $1}')"
manifest_sha256="$(sha256sum "${TASK_MANIFEST}" | awk '{print $1}')"
cat > "${FINAL_MARKER}" <<EOF
date=$(date -Is)
accepted_pairs=${TARGET}
source_commit=${source_commit}
task_manifest_sha256=${manifest_sha256}
root_file=${ROOT_OUT}
root_sha256=${root_sha256}
raa_dir=${RAA_OUT}
jet_pt20_dir=${JET20_OUT}
jet_pt30_dir=${JET30_OUT}
jet_slice_dir=${JET_SLICE_OUT}
status=PASS
EOF
echo "Wrote ${FINAL_MARKER}"

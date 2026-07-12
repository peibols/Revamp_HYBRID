#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
SOURCE="${SOURCE:-$(cd -- "${SCRIPT_DIR}/../../.." && pwd)}"
ROOT="${ROOT:-/raid5/data/yjlee/hybrid_dev}"
WORK="${WORK:-${ROOT}/test/oo5360_v2_500hydro_50k_20260711}"
LOCAL_EOS="${LOCAL_EOS:-${WORK}/eos_snapshot}"
PP_LOCAL_EOS="${PP_LOCAL_EOS:-${ROOT}/test/tmp_oo_10k_prehydro_raa_20260709/local_eos}"
PP_JET_CACHE="${PP_JET_CACHE:-${PP_LOCAL_EOS}/analysis_cache/oo5360_pp1m_jet_spectrum_v1.tsv}"
TASK_MANIFEST="${TASK_MANIFEST:-${WORK}/hydro_prepared/aa_task_manifest.tsv}"
STRICT_MARKER="${STRICT_MARKER:-${WORK}/v2_50k_strict_complete.txt}"
OUT="${OUT:-${WORK}/final_analysis}"
BUILD_DIR="${BUILD_DIR:-${ROOT}/test/tmp_oo_v2_final_root_build}"
OVERWRITE="${OVERWRITE:-false}"

TARGET=50000
RAA_OUT="${OUT}/hadron_raa"
ROOT_OUT="${OUT}/root/oo5360_c0_5_planb_v2_50k.root"
JET20_OUT="${OUT}/jets/inclusive_pt20"
JET_RAA_OUT="${OUT}/jets/raa"
JET30_OUT="${OUT}/jets/inclusive_pt30"
JET_SLICE_OUT="${OUT}/jets/pt_slices"
JET_CHARGE_OUT="${OUT}/jets/charge_response"
JET_PAIRED_OUT="${OUT}/jets/paired_substructure"
FINAL_MARKER="${OUT}/v2_final_analysis_complete.txt"

test -s "${STRICT_MARKER}"
test -d "${LOCAL_EOS}/status/aa"
test -d "${LOCAL_EOS}/outputs/aa"
test -d "${PP_LOCAL_EOS}/status/pp"
test -d "${PP_LOCAL_EOS}/outputs/pp"
test -s "${TASK_MANIFEST}"

mkdir -p "${RAA_OUT}" "$(dirname -- "${ROOT_OUT}")" \
  "${JET20_OUT}" "${JET30_OUT}" "${JET_RAA_OUT}" "${JET_SLICE_OUT}" "${JET_CHARGE_OUT}" \
  "${JET_PAIRED_OUT}" "${BUILD_DIR}"

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

python3 "${SCRIPT_DIR}/analyze_oo_jet_raa.py" \
  --input-root "${ROOT_OUT}" \
  --pp-local-eos "${PP_LOCAL_EOS}" \
  --out-dir "${JET_RAA_OUT}" \
  --pp-cache "${PP_JET_CACHE}" \
  --build-dir "${BUILD_DIR}" \
  --expected-aa-events "${TARGET}"

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

python3 "${SCRIPT_DIR}/plot_oo_jet_charge_response.py" \
  --input-root "${ROOT_OUT}" --out-dir "${JET_CHARGE_OUT}" \
  --pt-min 20 --prefix oo5360_v2_50k_jet_charge_response_pt20
python3 "${SCRIPT_DIR}/plot_oo_jet_charge_response.py" \
  --input-root "${ROOT_OUT}" --out-dir "${JET_CHARGE_OUT}" \
  --pt-min 30 --prefix oo5360_v2_50k_jet_charge_response_pt30
python3 "${SCRIPT_DIR}/plot_oo_jet_charge_response.py" \
  --input-root "${ROOT_OUT}" --out-dir "${JET_CHARGE_OUT}" \
  --pt-min 20 --pt-max 30 --prefix oo5360_v2_50k_jet_charge_response_pt20to30
python3 "${SCRIPT_DIR}/plot_oo_jet_charge_response.py" \
  --input-root "${ROOT_OUT}" --out-dir "${JET_CHARGE_OUT}" \
  --pt-min 30 --pt-max 50 --prefix oo5360_v2_50k_jet_charge_response_pt30to50
python3 "${SCRIPT_DIR}/plot_oo_jet_charge_response.py" \
  --input-root "${ROOT_OUT}" --out-dir "${JET_CHARGE_OUT}" \
  --pt-min 50 --pt-max 80 --prefix oo5360_v2_50k_jet_charge_response_pt50to80
python3 "${SCRIPT_DIR}/plot_oo_jet_charge_response.py" \
  --input-root "${ROOT_OUT}" --out-dir "${JET_CHARGE_OUT}" \
  --pt-min 80 --prefix oo5360_v2_50k_jet_charge_response_pt80plus

python3 "${SCRIPT_DIR}/plot_oo_jet_paired_substructure.py" \
  --input-root "${ROOT_OUT}" --out-dir "${JET_PAIRED_OUT}" \
  --pt-min 20 --prefix oo5360_v2_50k_jet_paired_substructure_pt20
python3 "${SCRIPT_DIR}/plot_oo_jet_paired_substructure.py" \
  --input-root "${ROOT_OUT}" --out-dir "${JET_PAIRED_OUT}" \
  --pt-min 30 --prefix oo5360_v2_50k_jet_paired_substructure_pt30
python3 "${SCRIPT_DIR}/plot_oo_jet_paired_substructure.py" \
  --input-root "${ROOT_OUT}" --out-dir "${JET_PAIRED_OUT}" \
  --pt-min 20 --pt-max 30 \
  --prefix oo5360_v2_50k_jet_paired_substructure_pt20to30
python3 "${SCRIPT_DIR}/plot_oo_jet_paired_substructure.py" \
  --input-root "${ROOT_OUT}" --out-dir "${JET_PAIRED_OUT}" \
  --pt-min 30 --pt-max 50 \
  --prefix oo5360_v2_50k_jet_paired_substructure_pt30to50
python3 "${SCRIPT_DIR}/plot_oo_jet_paired_substructure.py" \
  --input-root "${ROOT_OUT}" --out-dir "${JET_PAIRED_OUT}" \
  --pt-min 50 --pt-max 80 \
  --prefix oo5360_v2_50k_jet_paired_substructure_pt50to80
python3 "${SCRIPT_DIR}/plot_oo_jet_paired_substructure.py" \
  --input-root "${ROOT_OUT}" --out-dir "${JET_PAIRED_OUT}" \
  --pt-min 80 --prefix oo5360_v2_50k_jet_paired_substructure_pt80plus

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
  "${TARGET}" \
  "${JET_PAIRED_OUT}" \
  "${JET_RAA_OUT}/oo5360_c0_5_jet_raa_R020408_metadata.json" \
  "${JET_RAA_OUT}/oo5360_c0_5_jet_raa_R020408.tsv" <<'PY'
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

conversion_path, closure_path, manifest_path = map(Path, sys.argv[1:4])
target = int(sys.argv[4])
paired_dir = Path(sys.argv[5])
jet_raa_metadata_path = Path(sys.argv[6])
jet_raa_table_path = Path(sys.argv[7])
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
paired_metadata = sorted(paired_dir.glob("*_metadata.json"))
if len(paired_metadata) != 6:
    raise SystemExit(
        f"expected six paired-substructure metadata files, found "
        f"{len(paired_metadata)}"
    )
for metadata_path in paired_metadata:
    paired = json.loads(metadata_path.read_text())
    if paired.get("pairCount") != target:
        raise SystemExit(
            f"{metadata_path.name} used {paired.get('pairCount')} pairs, "
            f"expected {target}"
        )
    if any(
        radius.get("pairPtAudit") != "PASS"
        or radius.get("pairHardTagAudit") != "PASS"
        or radius.get("oneToOneMatchAudit") != "PASS"
        for radius in paired.get("radii", {}).values()
    ):
        raise SystemExit(f"paired match audit failed in {metadata_path.name}")
jet_raa = json.loads(jet_raa_metadata_path.read_text())
with jet_raa_table_path.open() as handle:
    jet_raa_rows = list(csv.DictReader(handle, delimiter="\t"))
if (
    jet_raa.get("status") != "PASS"
    or jet_raa.get("aaEvents") != target
    or jet_raa.get("ppEvents") != 1_000_000
    or len(jet_raa_rows) != 54
    or {float(row["radius"]) for row in jet_raa_rows} != {0.2, 0.4, 0.8}
    or any(
        not math.isfinite(float(row[key]))
        for row in jet_raa_rows
        for key in ("raa", "stat_error", "aa_spectrum_mb_per_gev", "pp_spectrum_mb_per_gev")
    )
):
    raise SystemExit("jet RAA output failed event-count, row-count, radius, or finite-value audit")
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
jet_raa_dir=${JET_RAA_OUT}
jet_slice_dir=${JET_SLICE_OUT}
jet_charge_response_dir=${JET_CHARGE_OUT}
jet_paired_substructure_dir=${JET_PAIRED_OUT}
status=PASS
EOF
echo "Wrote ${FINAL_MARKER}"

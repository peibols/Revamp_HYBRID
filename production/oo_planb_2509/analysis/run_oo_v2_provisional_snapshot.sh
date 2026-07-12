#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
SOURCE="${SOURCE:-$(cd -- "${SCRIPT_DIR}/../../.." && pwd)}"
ROOT="${ROOT:-/raid5/data/yjlee/hybrid_dev}"
SNAPSHOT="${1:?usage: run_oo_v2_provisional_snapshot.sh SNAPSHOT_DIR}"
PP_LOCAL_EOS="${PP_LOCAL_EOS:-${ROOT}/test/tmp_oo_10k_prehydro_raa_20260709/local_eos}"
PP_JET_CACHE="${PP_JET_CACHE:-${PP_LOCAL_EOS}/analysis_cache/oo5360_pp1m_jet_spectrum_r01020408_v2.tsv}"
BUILD_DIR="${BUILD_DIR:-${SNAPSHOT}/build}"
OVERWRITE="${OVERWRITE:-false}"

LOCAL_EOS="${SNAPSHOT}/local_eos"
TASK_MANIFEST="${SNAPSHOT}/aa_task_manifest.tsv"
ACCEPTED_IDS="${SNAPSHOT}/accepted_task_ids.txt"
COUNT="$(wc -l < "${ACCEPTED_IDS}")"
PREFIX="oo5360_v2_provisional${COUNT}"
RAA_OUT="${SNAPSHOT}/hadron_raa"
ROOT_OUT="${SNAPSHOT}/root/oo5360_c0_5_planb_v2_provisional${COUNT}.root"
JET_OUT="${SNAPSHOT}/jet_variables"
JET_RAA_OUT="${SNAPSHOT}/jet_raa"
CHARGE_OUT="${SNAPSHOT}/charge_response"
PAIRED_OUT="${SNAPSHOT}/paired_substructure"
VALIDATION_OUT="${SNAPSHOT}/validation"
FINAL_MARKER="${SNAPSHOT}/provisional_analysis_complete.txt"

test -s "${SNAPSHOT}/snapshot_metadata.tsv"
test -s "${TASK_MANIFEST}"
test -d "${LOCAL_EOS}/status/aa"
test -d "${LOCAL_EOS}/outputs/aa"
test "${COUNT}" -gt 0
mkdir -p "${RAA_OUT}" "$(dirname -- "${ROOT_OUT}")" "${JET_OUT}" \
  "${JET_RAA_OUT}" "${CHARGE_OUT}" "${PAIRED_OUT}" "${VALIDATION_OUT}" \
  "${BUILD_DIR}"

python3 "${SCRIPT_DIR}/analyze_oo_prehydro_pair.py" \
  --local-eos "${LOCAL_EOS}" \
  --pp-local-eos "${PP_LOCAL_EOS}" \
  --out-dir "${RAA_OUT}" \
  --require-paired-aa \
  --aa-events-per-chunk 1 \
  --aa-task-manifest "${TASK_MANIFEST}"

converter_args=(
  --source "v2_fullstats=${LOCAL_EOS}"
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
  --expected-aa-events "${COUNT}"

run_jet_variables() {
  local suffix="$1"
  shift
  python3 "${SCRIPT_DIR}/plot_oo_jet_variables.py" \
    --input-root "${ROOT_OUT}" --out-dir "${JET_OUT}" \
    "$@" --prefix "${PREFIX}_jet_variables_${suffix}"
}
run_jet_variables pt20 --pt-min 20
run_jet_variables pt30 --pt-min 30
run_jet_variables pt20to30 --pt-min 20 --pt-max 30
run_jet_variables pt30to50 --pt-min 30 --pt-max 50
run_jet_variables pt50to80 --pt-min 50 --pt-max 80
run_jet_variables pt80plus --pt-min 80

python3 "${SCRIPT_DIR}/summarize_oo_jet_pt_slices.py" \
  --slice-metadata \
    "${JET_OUT}/${PREFIX}_jet_variables_pt20to30_metadata.json" \
    "${JET_OUT}/${PREFIX}_jet_variables_pt30to50_metadata.json" \
    "${JET_OUT}/${PREFIX}_jet_variables_pt50to80_metadata.json" \
    "${JET_OUT}/${PREFIX}_jet_variables_pt80plus_metadata.json" \
  --inclusive-metadata "${JET_OUT}/${PREFIX}_jet_variables_pt20_metadata.json" \
  --out-dir "${JET_OUT}" --prefix "${PREFIX}_jet_pt_slice_summary"

run_charge_response() {
  local suffix="$1"
  shift
  python3 "${SCRIPT_DIR}/plot_oo_jet_charge_response.py" \
    --input-root "${ROOT_OUT}" --out-dir "${CHARGE_OUT}" \
    "$@" --prefix "${PREFIX}_jet_charge_response_${suffix}"
}
run_charge_response pt20 --pt-min 20
run_charge_response pt30 --pt-min 30
run_charge_response pt20to30 --pt-min 20 --pt-max 30
run_charge_response pt30to50 --pt-min 30 --pt-max 50
run_charge_response pt50to80 --pt-min 50 --pt-max 80
run_charge_response pt80plus --pt-min 80
run_charge_response pt30_strictRquarter \
  --pt-min 30 --max-pair-match-dr-fraction 0.25
run_charge_response pt30to50_strictRquarter \
  --pt-min 30 --pt-max 50 --max-pair-match-dr-fraction 0.25

run_paired_substructure() {
  local suffix="$1"
  shift
  python3 "${SCRIPT_DIR}/plot_oo_jet_paired_substructure.py" \
    --input-root "${ROOT_OUT}" --out-dir "${PAIRED_OUT}" \
    "$@" --prefix "${PREFIX}_paired_${suffix}"
}
run_paired_substructure pt20 --pt-min 20
run_paired_substructure pt30 --pt-min 30
run_paired_substructure pt20to30 --pt-min 20 --pt-max 30
run_paired_substructure pt30to50 --pt-min 30 --pt-max 50
run_paired_substructure pt50to80 --pt-min 50 --pt-max 80
run_paired_substructure pt80plus --pt-min 80

python3 - \
  "${SNAPSHOT}" "${ROOT_OUT%.root}.summary.json" \
  "${JET_OUT}/${PREFIX}_jet_pt_slice_summary_validation.json" \
  "${JET_RAA_OUT}/oo5360_c0_5_jet_raa_R01020408_metadata.json" \
  "${JET_RAA_OUT}/oo5360_c0_5_jet_raa_R01020408.tsv" \
  "${COUNT}" "${PREFIX}" <<'PY'
import csv
import json
import math
from pathlib import Path
import sys

snapshot, conversion_path, closure_path, jet_raa_metadata_path, jet_raa_table_path = map(
    Path, sys.argv[1:6]
)
expected = int(sys.argv[6])
prefix = sys.argv[7]
conversion = json.loads(conversion_path.read_text())
closure = json.loads(closure_path.read_text())
jet_raa_metadata = json.loads(jet_raa_metadata_path.read_text())
with (snapshot / "hadron_raa/oo5360_c0_5_prehydro_overlay_raa.tsv").open() as handle:
    raa_rows = list(csv.DictReader(handle, delimiter="\t"))
with jet_raa_table_path.open() as handle:
    jet_raa_rows = list(csv.DictReader(handle, delimiter="\t"))

if conversion["acceptedPairs"] != expected:
    raise SystemExit("ROOT accepted-pair count does not match the frozen snapshot")
if conversion["rejectedArchives"] or conversion["skippedChunkRecords"]:
    raise SystemExit("ROOT conversion rejected or skipped a frozen strict archive")
if len(raa_rows) != 18 or {int(row["aa_events"]) for row in raa_rows} != {expected}:
    raise SystemExit("RAA table does not contain the expected frozen event count")
if (
    jet_raa_metadata.get("status") != "PASS"
    or jet_raa_metadata.get("aaEvents") != expected
    or jet_raa_metadata.get("ppEvents") != 1_000_000
):
    raise SystemExit("jet RAA metadata does not contain the expected AA/pp event counts")
if (
    len(jet_raa_rows) != 72
    or {row["variant"] for row in jet_raa_rows} != {"noPrehydro", "withPrehydro"}
    or {float(row["radius"]) for row in jet_raa_rows} != {0.1, 0.2, 0.4, 0.8}
    or any(
        not math.isfinite(float(row[key]))
        for row in jet_raa_rows
        for key in ("raa", "stat_error", "aa_spectrum_mb_per_gev", "pp_spectrum_mb_per_gev")
    )
):
    raise SystemExit("jet RAA table failed row-count, radius, variant, or finite-value audit")
if closure.get("status") != "PASS":
    raise SystemExit("jet pT-slice closure failed")
residuals = [
    abs(value)
    for radius in closure["crossSectionClosureDifferenceMb"].values()
    for value in radius.values()
]
if any(not math.isfinite(value) or value > 1.0e-12 for value in residuals):
    raise SystemExit("jet pT-slice closure exceeds 1e-12 mb")

metadata_files = [
    *sorted((snapshot / "jet_variables").glob(f"{prefix}_jet_variables_*_metadata.json")),
    *sorted((snapshot / "charge_response").glob(f"{prefix}_jet_charge_response_*_metadata.json")),
    *sorted((snapshot / "paired_substructure").glob(f"{prefix}_paired_*_metadata.json")),
]
for path in metadata_files:
    metadata = json.loads(path.read_text())
    if metadata.get("pairCount") not in (None, expected):
        raise SystemExit(f"{path.name}: pairCount does not match snapshot")

report = {
    "status": "PASS",
    "publicationStatus": "PROVISIONAL_DIAGNOSTIC_NOT_UNBIASED",
    "acceptedPairs": expected,
    "rootAcceptedPairs": conversion["acceptedPairs"],
    "rootRejectedArchives": conversion["rejectedArchives"],
    "rootSkippedChunkRecords": conversion["skippedChunkRecords"],
    "ppEvents": int(raa_rows[0]["pp_events"]),
    "raaRows": len(raa_rows),
    "jetRaaRows": len(jet_raa_rows),
    "jetMetadataFiles": len(metadata_files),
    "pdfPlots": len(list(snapshot.rglob("*.pdf"))),
    "pngPlots": len(list(snapshot.rglob("*.png"))),
    "maxJetSliceClosureResidualMb": max(residuals, default=0.0),
    "jetSliceClosureToleranceMb": 1.0e-12,
}
(snapshot / "validation/validation_report.json").write_text(
    json.dumps(report, indent=2, sort_keys=True) + "\n"
)
print(json.dumps(report, indent=2, sort_keys=True))
PY

source_commit="$(git -C "${SOURCE}" rev-parse HEAD)"
root_sha256="$(sha256sum "${ROOT_OUT}" | awk '{print $1}')"
cat > "${FINAL_MARKER}" <<EOF
date=$(date -Is)
accepted_pairs=${COUNT}
selection=strict-complete-at-frozen-boundary_completion-order-selected
publication_status=PROVISIONAL_DIAGNOSTIC_NOT_UNBIASED
source_commit=${source_commit}
root_file=${ROOT_OUT}
root_sha256=${root_sha256}
jet_raa_dir=${JET_RAA_OUT}
status=PASS
EOF
echo "Wrote ${FINAL_MARKER}"

#!/usr/bin/env bash
set -euo pipefail
ROOT="${ROOT:-/raid5/data/yjlee/hybrid_dev}"
WORK="${WORK:-${ROOT}/test/oo5360_no_moliere_raa_20260606}"
SUPPORT="${WORK}/cern_support"
CERNCTL="${CERNCTL:-/data/yjlee/cernLxplus/cernctl}"
CERN_REMOTE="${CERN_REMOTE:-lxplus}"
CAMPAIGN="${CAMPAIGN:-hybrid_oo5360_c0_5_no_moliere_paired_prehydro_pthat4_unbounded_20260709}"
EOS_BASE="${EOS_BASE:-/eos/user/y/yjlee/${CAMPAIGN}}"
JOB_PAYLOAD_EOS_BASE="${JOB_PAYLOAD_EOS_BASE:-${EOS_BASE}}"
JOB_RUNTIME_PAYLOAD_EOS_BASE="${JOB_RUNTIME_PAYLOAD_EOS_BASE:-${EOS_BASE}}"
AFS_WORK="${AFS_WORK:-/afs/cern.ch/user/y/yjlee/cernLxplus_jobs/${CAMPAIGN}}"
TMP_REMOTE="${TMP_REMOTE:-/tmp/yjlee_${CAMPAIGN}_payloads}"
DRY_RUN="${DRY_RUN:-false}"
SUBMIT_PP="${SUBMIT_PP:-true}"
PP_REFERENCE_CAMPAIGN="${PP_REFERENCE_CAMPAIGN:-}"
PYTHIA_SOURCE="${PYTHIA_SOURCE:-/eos/user/y/yjlee/hybrid_mmli_smoke_20260505_120203/payloads/pythia8315_alma9_install.tar.gz}"
MMLI_SOURCE_ROOT="${MMLI_SOURCE_ROOT:-${ROOT}/wt_main_moliere_lres_integration_clean}"
AA_CHUNKS="${AA_CHUNKS:-220}"
PP_CHUNKS="${PP_CHUNKS:-100}"
EVENTS="${EVENTS:-1}"
AA_EVENTS="${AA_EVENTS:-1}"
PP_EVENTS="${PP_EVENTS:-${EVENTS}}"
AA_SEED_OFFSET="${AA_SEED_OFFSET:-610000}"
PP_SEED_OFFSET="${PP_SEED_OFFSET:-710000}"
BUILD_JOB_FLAVOUR="${BUILD_JOB_FLAVOUR:-workday}"
JOB_FLAVOUR="${JOB_FLAVOUR:-workday}"
TIMEOUT_S="${TIMEOUT_S:-7200}"
PTHAT_MIN="${PTHAT_MIN:-4}"
PTHAT_MAX="${PTHAT_MAX:--1}"
PDF_MODE="${PDF_MODE:-lhapdf}"
LHAPDF_SET="${LHAPDF_SET:-EPPS21nlo_CT18Anlo_O16/0}"
AA_PDF_MODE="${AA_PDF_MODE:-${PDF_MODE}}"
AA_LHAPDF_SET="${AA_LHAPDF_SET:-${LHAPDF_SET}}"
PP_PDF_MODE="${PP_PDF_MODE:-off}"
PP_LHAPDF_SET="${PP_LHAPDF_SET:-}"
LHAPDF_DATA_ROOT="${LHAPDF_DATA_ROOT:-${ROOT}/test/oo5360_no_moliere_local_clear_20260607/lhapdf_data}"
LHAPDF_CVMFS_VIEW="${LHAPDF_CVMFS_VIEW:-/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/lhapdf/6.5.3-3fa11/x86_64-centos7-gcc11-opt}"
AA_CENTRALITY_INDEX="${AA_CENTRALITY_INDEX:-0}"
AA_CENTRALITY_LABEL="${AA_CENTRALITY_LABEL:-C0-5}"
RUN_PREHYDRO_PAIR="${RUN_PREHYDRO_PAIR:-true}"
NO_PREHYDRO_ALPHA="${NO_PREHYDRO_ALPHA:-0.37}"
PREHYDRO_ALPHA="${PREHYDRO_ALPHA:-0.37}"
BROADENING_K="${BROADENING_K:-15.0}"
PREHYDRO_TAU_MIN="${PREHYDRO_TAU_MIN:-0.24}"
PREHYDRO_TAU_GRID="${PREHYDRO_TAU_GRID:-0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.399}"
PREHYDRO_ETA_OVER_S="${PREHYDRO_ETA_OVER_S:-0.12}"
PREHYDRO_EOS_FACTOR="${PREHYDRO_EOS_FACTOR:-15.626873635058152}"
TOLERATE_CHUNK_FAILURES="${TOLERATE_CHUNK_FAILURES:-true}"
DEFAULT_PREHYDRO_ATTRACTOR_TABLE="${WORK}/reference_data/qcd_kinetic_attractor_lambda10_Cinf0p87.tsv"
PREHYDRO_ATTRACTOR_TABLE="${PREHYDRO_ATTRACTOR_TABLE:-${DEFAULT_PREHYDRO_ATTRACTOR_TABLE}}"
PREHYDRO_VISCOUS_ANCHOR="${PREHYDRO_VISCOUS_ANCHOR:-true}"
AA_TASK_MANIFEST="${AA_TASK_MANIFEST:-}"
HYDRO_MANIFEST="${HYDRO_MANIFEST:-}"
HYDRO_PAYLOAD_DIR="${HYDRO_PAYLOAD_DIR:-}"
HYDRO_PREPARATION_SUMMARY="${HYDRO_PREPARATION_SUMMARY:-}"
STORE_PREHYDRO_TABLE="${STORE_PREHYDRO_TABLE:-true}"
V2_MULTI_HYDRO=false
if [[ -n "${AA_TASK_MANIFEST}${HYDRO_MANIFEST}${HYDRO_PAYLOAD_DIR}" ]]; then
  V2_MULTI_HYDRO=true
  for REQUIRED in AA_TASK_MANIFEST HYDRO_MANIFEST HYDRO_PAYLOAD_DIR; do
    if [[ -z "${!REQUIRED}" ]]; then
      echo "${REQUIRED} is required for v2 multi-hydro submission" >&2
      exit 1
    fi
  done
fi

for VALUE_NAME in AA_CHUNKS PP_CHUNKS AA_EVENTS PP_EVENTS; do
  VALUE="${!VALUE_NAME}"
  if [[ ! "${VALUE}" =~ ^[1-9][0-9]*$ ]]; then
    echo "${VALUE_NAME} must be a positive integer, got: ${VALUE}" >&2
    exit 1
  fi
done
python3 - "${NO_PREHYDRO_ALPHA}" "${PREHYDRO_ALPHA}" "${BROADENING_K}" <<'PY'
import math
import sys

names = ("NO_PREHYDRO_ALPHA", "PREHYDRO_ALPHA", "BROADENING_K")
for name, raw in zip(names, sys.argv[1:]):
    try:
        value = float(raw)
    except ValueError as exc:
        raise SystemExit(f"{name} must be numeric, got: {raw}") from exc
    if not math.isfinite(value) or value < 0.0:
        raise SystemExit(f"{name} must be finite and nonnegative, got: {raw}")
PY
case "${SUBMIT_PP}" in
  1|true|TRUE) SUBMIT_PP=true ;;
  0|false|FALSE) SUBMIT_PP=false ;;
  *) echo "SUBMIT_PP must be true or false, got: ${SUBMIT_PP}" >&2; exit 1 ;;
esac
case "${TOLERATE_CHUNK_FAILURES}" in
  1|true|TRUE) TOLERATE_CHUNK_FAILURES=true ;;
  0|false|FALSE) TOLERATE_CHUNK_FAILURES=false ;;
  *) echo "TOLERATE_CHUNK_FAILURES must be true or false, got: ${TOLERATE_CHUNK_FAILURES}" >&2; exit 1 ;;
esac
if [[ "${SUBMIT_PP}" == "false" && -z "${PP_REFERENCE_CAMPAIGN}" ]]; then
  echo "PP_REFERENCE_CAMPAIGN is required when SUBMIT_PP=false" >&2
  exit 1
fi
if [[ "${AA_EVENTS}" != "1" ]]; then
  echo "AA_EVENTS must be 1: each AA job runs one hard event in both paired variants" >&2
  exit 1
fi
if [[ "${V2_MULTI_HYDRO}" == "true" ]]; then
  python3 - "${AA_TASK_MANIFEST}" "${HYDRO_MANIFEST}" "${HYDRO_PAYLOAD_DIR}" "${AA_CHUNKS}" "${AA_SEED_OFFSET}" <<'PY'
import csv
import hashlib
from pathlib import Path
import sys

task_path, hydro_path, payload_dir = map(Path, sys.argv[1:4])
expected_tasks, seed_offset = map(int, sys.argv[4:6])
for path in (task_path, hydro_path):
    if not path.is_file():
        raise SystemExit(f"v2 manifest not found: {path}")
if not payload_dir.is_dir():
    raise SystemExit(f"v2 hydro payload directory not found: {payload_dir}")
with hydro_path.open(newline="") as handle:
    hydros = list(csv.DictReader(handle, delimiter="\t"))
with task_path.open(newline="") as handle:
    tasks = list(csv.DictReader(handle, delimiter="\t"))
if len(hydros) != 500:
    raise SystemExit(f"expected 500 hydro rows, found {len(hydros)}")
if len(tasks) != expected_tasks:
    raise SystemExit(f"AA_CHUNKS={expected_tasks}, task manifest has {len(tasks)} rows")
if [int(row["task_id"]) for row in tasks] != list(range(expected_tasks)):
    raise SystemExit("task IDs must be exactly 0..AA_CHUNKS-1")
if [int(row["hard_seed"]) for row in tasks] != list(range(seed_offset, seed_offset + expected_tasks)):
    raise SystemExit("task hard seeds do not match AA_SEED_OFFSET")
by_slot = {int(row["hydro_slot"]): row for row in hydros}
if set(by_slot) != set(range(500)):
    raise SystemExit("hydro slots must be exactly 0..499")
for row in tasks:
    hydro = by_slot[int(row["hydro_slot"])]
    expected = {
        "hydro_event_id": hydro["event_id"],
        "hydro_ncoll": hydro["ncoll"],
        "hydro_dir": hydro["hydro_dir"],
        "hydro_payload_key": hydro["payload_key"],
        "hydro_payload_sha256": hydro["payload_sha256"],
    }
    for key, value in expected.items():
        if row[key] != value:
            raise SystemExit(f"task {row['task_id']} {key} does not match hydro manifest")
for hydro in hydros:
    payload = payload_dir / Path(hydro["payload_key"]).name
    if not payload.is_file() or payload.stat().st_size != int(hydro["payload_size"]):
        raise SystemExit(f"missing or wrong-size payload: {payload}")
    digest = hashlib.sha256(payload.read_bytes()).hexdigest()
    if digest != hydro["payload_sha256"]:
        raise SystemExit(f"payload checksum mismatch: {payload}")
print(f"validated v2 manifests: {len(tasks)} tasks, {len(hydros)} hydro payloads")
PY
fi
if [[ "${PREHYDRO_ETA_OVER_S}" != "0.12" ]]; then
  echo "PREHYDRO_ETA_OVER_S must be 0.12 for arXiv:2509.19430v2 Plan B" >&2
  exit 1
fi
case "${PREHYDRO_VISCOUS_ANCHOR}" in
  1|true|TRUE) ;;
  *)
    echo "PREHYDRO_VISCOUS_ANCHOR must be true for arXiv:2509.19430v2 Plan B" >&2
    exit 1
    ;;
esac
if [[ ! -f "${PREHYDRO_ATTRACTOR_TABLE}" ]]; then
  echo "Published prehydro attractor table not found: ${PREHYDRO_ATTRACTOR_TABLE}" >&2
  exit 1
fi
PREHYDRO_ATTRACTOR_SHA256="$(sha256sum "${PREHYDRO_ATTRACTOR_TABLE}" | awk '{print $1}')"
EXPECTED_PREHYDRO_ATTRACTOR_SHA256="1bea7289d3dc8ed95819eaa86cf4c489442a054c14aae47eff010cf45155eba0"
if [[ "${PREHYDRO_ATTRACTOR_SHA256}" != "${EXPECTED_PREHYDRO_ATTRACTOR_SHA256}" ]]; then
  echo "Unexpected prehydro attractor SHA256: ${PREHYDRO_ATTRACTOR_SHA256}" >&2
  echo "Expected published table: ${EXPECTED_PREHYDRO_ATTRACTOR_SHA256}" >&2
  exit 1
fi

PAIR_TAG="single"
if [[ "${RUN_PREHYDRO_PAIR}" == "1" || "${RUN_PREHYDRO_PAIR}" == "true" || "${RUN_PREHYDRO_PAIR}" == "TRUE" ]]; then
  PREHYDRO_TAU_TAG="${PREHYDRO_TAU_MIN//./p}"
  PREHYDRO_ALPHA_TAG="${PREHYDRO_ALPHA//./p}"
  PAIR_TAG="pairedPrehydroTau${PREHYDRO_TAU_TAG}Alpha${PREHYDRO_ALPHA_TAG}"
fi
CENTRALITY_TAG="${AA_CENTRALITY_LABEL//-/_}"
PDF_TAG="aa${AA_PDF_MODE}_pp${PP_PDF_MODE}"
RUN_NAME="${RUN_NAME:-runs/oo5360_no_moliere_${CENTRALITY_TAG}_${PDF_TAG}_pthat${PTHAT_MIN}_${PTHAT_MAX}_${PAIR_TAG}_aa${AA_CHUNKS}x${AA_EVENTS}_pp${PP_CHUNKS}x${PP_EVENTS}}"
PAYLOAD_DIR="${WORK}/payloads"
mkdir -p "${PAYLOAD_DIR}" "${WORK}/logs"
AA_TOTAL_EVENTS=$((AA_CHUNKS * AA_EVENTS))
PP_TOTAL_EVENTS=$((PP_CHUNKS * PP_EVENTS))
PP_SUBMITTED_EVENTS="${PP_TOTAL_EVENTS}"
if [[ "${SUBMIT_PP}" == "false" ]]; then
  PP_SUBMITTED_EVENTS=0
fi
HYDRO_PAYLOAD_DESCRIPTION="staged_hydro one event per Zenodo centrality bin C0-5...C90-100"
ANALYSIS_WEIGHT_DESCRIPTION="centrality_width*ncoll"
if [[ "${V2_MULTI_HYDRO}" == "true" ]]; then
  HYDRO_PAYLOAD_DESCRIPTION="500 event-by-event Zenodo C0-5 MUSIC hydros; one checksum-pinned assigned payload fetched per AA job"
  ANALYSIS_WEIGHT_DESCRIPTION="no extra Ncoll factor: task sampling is proportional to hydro-event Ncoll within every 5000-task block"
fi

for MODE in "${AA_PDF_MODE}" "${PP_PDF_MODE}"; do
  case "${MODE}" in
    native_npdf|off|lhapdf) ;;
    *) echo "Unsupported PDF mode: ${MODE}" >&2; exit 1 ;;
  esac
done
if [[ "${AA_PDF_MODE}" == "lhapdf" && -z "${AA_LHAPDF_SET}" ]]; then
  echo "AA_LHAPDF_SET is required when AA_PDF_MODE=lhapdf" >&2
  exit 1
fi
if [[ "${PP_PDF_MODE}" == "lhapdf" && -z "${PP_LHAPDF_SET}" ]]; then
  echo "PP_LHAPDF_SET is required when PP_PDF_MODE=lhapdf" >&2
  exit 1
fi

echo "campaign=${CAMPAIGN}"
echo "eos_base=${EOS_BASE}"
echo "job_payload_eos_base=${JOB_PAYLOAD_EOS_BASE}"
echo "job_runtime_payload_eos_base=${JOB_RUNTIME_PAYLOAD_EOS_BASE}"
echo "afs_work=${AFS_WORK}"
echo "mmli_source=${MMLI_SOURCE_ROOT}"
echo "dry_run=${DRY_RUN}"
echo "aa_centrality=${AA_CENTRALITY_LABEL} aa_centrality_index=${AA_CENTRALITY_INDEX}"
echo "aa_chunks=${AA_CHUNKS} aa_events_per_chunk=${AA_EVENTS} aa_total_events=${AA_TOTAL_EVENTS}"
echo "submit_pp=${SUBMIT_PP} pp_chunks=${PP_CHUNKS} pp_events_per_chunk=${PP_EVENTS} pp_total_events=${PP_TOTAL_EVENTS} pp_reference_campaign=${PP_REFERENCE_CAMPAIGN:-none}"
echo "job_flavour=${JOB_FLAVOUR} timeout_s=${TIMEOUT_S} pthat_min=${PTHAT_MIN} pthat_max=${PTHAT_MAX}"
echo "tolerate_chunk_failures=${TOLERATE_CHUNK_FAILURES}"
echo "no_prehydro_alpha=${NO_PREHYDRO_ALPHA} prehydro_alpha=${PREHYDRO_ALPHA} broadening_k=${BROADENING_K}"
echo "aa_pdf_mode=${AA_PDF_MODE} aa_lhapdf_set=${AA_LHAPDF_SET}"
echo "pp_pdf_mode=${PP_PDF_MODE} pp_lhapdf_set=${PP_LHAPDF_SET:-none}"
echo "run_prehydro_pair=${RUN_PREHYDRO_PAIR} prehydro_tau_min=${PREHYDRO_TAU_MIN} prehydro_tau_grid=${PREHYDRO_TAU_GRID} prehydro_attractor_table=${PREHYDRO_ATTRACTOR_TABLE}"
echo "prehydro_attractor_sha256=${PREHYDRO_ATTRACTOR_SHA256}"
echo "v2_multi_hydro=${V2_MULTI_HYDRO} aa_task_manifest=${AA_TASK_MANIFEST:-none} hydro_manifest=${HYDRO_MANIFEST:-none}"

if ! git -C "${MMLI_SOURCE_ROOT}" diff --quiet ||
   ! git -C "${MMLI_SOURCE_ROOT}" diff --cached --quiet; then
  echo "Tracked HYBRID source changes must be committed before packaging" >&2
  exit 1
fi
SOURCE_GIT_BRANCH="$(git -C "${MMLI_SOURCE_ROOT}" branch --show-current)"
SOURCE_GIT_COMMIT="$(git -C "${MMLI_SOURCE_ROOT}" rev-parse HEAD)"
git -C "${MMLI_SOURCE_ROOT}" archive --format=tar "${SOURCE_GIT_COMMIT}" | gzip -n > "${PAYLOAD_DIR}/mmli_source.tar.gz"
SOURCE_PAYLOAD_SHA256="$(sha256sum "${PAYLOAD_DIR}/mmli_source.tar.gz" | awk '{print $1}')"
rm -rf "${PAYLOAD_DIR}/runtime_payload"
mkdir -p "${PAYLOAD_DIR}/runtime_payload/runtime"
cp "${MMLI_SOURCE_ROOT}/test/setup_pythia.cmnd" "${PAYLOAD_DIR}/runtime_payload/runtime/setup_pythia.cmnd"
cp "${SUPPORT}/run_oo_validation_chunk.py" "${PAYLOAD_DIR}/runtime_payload/runtime/"
PREHYDRO_ATTRACTOR_TABLE_RUNTIME=""
if [[ -n "${PREHYDRO_ATTRACTOR_TABLE}" ]]; then
  test -f "${PREHYDRO_ATTRACTOR_TABLE}"
  cp "${PREHYDRO_ATTRACTOR_TABLE}" "${PAYLOAD_DIR}/runtime_payload/runtime/prehydro_attractor_table.dat"
  PREHYDRO_ATTRACTOR_TABLE_RUNTIME="runtime/prehydro_attractor_table.dat"
fi
copy_lhapdf_set() {
  local SET_NAME="$1"
  local SET_DIR="${SET_NAME%%/*}"
  test -d "${LHAPDF_DATA_ROOT}/${SET_DIR}"
  mkdir -p "${PAYLOAD_DIR}/runtime_payload/runtime/lhapdf_data"
  if [[ ! -d "${PAYLOAD_DIR}/runtime_payload/runtime/lhapdf_data/${SET_DIR}" ]]; then
    cp -a "${LHAPDF_DATA_ROOT}/${SET_DIR}" "${PAYLOAD_DIR}/runtime_payload/runtime/lhapdf_data/"
  fi
}
if [[ "${AA_PDF_MODE}" == "lhapdf" ]]; then
  copy_lhapdf_set "${AA_LHAPDF_SET}"
fi
if [[ "${PP_PDF_MODE}" == "lhapdf" ]]; then
  copy_lhapdf_set "${PP_LHAPDF_SET}"
fi
if [[ "${AA_PDF_MODE}" == "lhapdf" || "${PP_PDF_MODE}" == "lhapdf" ]]; then
  mkdir -p "${PAYLOAD_DIR}/runtime_payload/runtime/lhapdf_data"
  cp -a "${LHAPDF_DATA_ROOT}/pdfsets.index" "${PAYLOAD_DIR}/runtime_payload/runtime/lhapdf_data/"
  echo "${LHAPDF_CVMFS_VIEW}" > "${PAYLOAD_DIR}/runtime_payload/runtime/lhapdf_cvmfs_view.txt"
fi
AA_TASK_MANIFEST_RUNTIME=""
HYDRO_MANIFEST_SHA256="none"
AA_TASK_MANIFEST_SHA256="none"
if [[ "${V2_MULTI_HYDRO}" == "true" ]]; then
  cp "${AA_TASK_MANIFEST}" "${PAYLOAD_DIR}/runtime_payload/runtime/aa_task_manifest.tsv"
  cp "${HYDRO_MANIFEST}" "${PAYLOAD_DIR}/runtime_payload/runtime/hydro_manifest.tsv"
  if [[ -n "${HYDRO_PREPARATION_SUMMARY}" ]]; then
    cp "${HYDRO_PREPARATION_SUMMARY}" "${PAYLOAD_DIR}/runtime_payload/runtime/hydro_preparation_summary.json"
  fi
  mkdir -p "${PAYLOAD_DIR}/runtime_payload/runtime/staged_hydro"
  AA_TASK_MANIFEST_RUNTIME="runtime/aa_task_manifest.tsv"
  HYDRO_MANIFEST_SHA256="$(sha256sum "${HYDRO_MANIFEST}" | awk '{print $1}')"
  AA_TASK_MANIFEST_SHA256="$(sha256sum "${AA_TASK_MANIFEST}" | awk '{print $1}')"
else
  tar -czf "${PAYLOAD_DIR}/runtime_payload/runtime/staged_hydro.tar.gz" -C "${WORK}" staged_hydro hydro_manifest.tsv
  (
    cd "${PAYLOAD_DIR}/runtime_payload/runtime"
    tar -xzf staged_hydro.tar.gz
    rm staged_hydro.tar.gz
  )
fi
tar -czf "${PAYLOAD_DIR}/runtime_payload.tar.gz" -C "${PAYLOAD_DIR}/runtime_payload" .
RUNTIME_PAYLOAD_SHA256="$(sha256sum "${PAYLOAD_DIR}/runtime_payload.tar.gz" | awk '{print $1}')"
seq 0 $((AA_CHUNKS - 1)) > "${WORK}/aa_chunk_ids.txt"
if [[ "${SUBMIT_PP}" == "true" ]]; then
  seq 0 $((PP_CHUNKS - 1)) > "${WORK}/pp_chunk_ids.txt"
else
  : > "${WORK}/pp_chunk_ids.txt"
fi

if [[ "${DRY_RUN}" == "1" || "${DRY_RUN}" == "true" || "${DRY_RUN}" == "TRUE" ]]; then
  echo "dry_run_note=skipping remote mkdir/scp/EOS payload staging"
else
  "${CERNCTL}" run mkdir -p "${AFS_WORK}/log" "${EOS_BASE}/payloads" "${EOS_BASE}/outputs/aa" "${EOS_BASE}/outputs/pp" "${EOS_BASE}/status/aa" "${EOS_BASE}/status/pp" "${TMP_REMOTE}"
  scp -q -o BatchMode=yes "${PAYLOAD_DIR}/mmli_source.tar.gz" "${PAYLOAD_DIR}/runtime_payload.tar.gz" "${CERN_REMOTE}:${TMP_REMOTE}/"
  if [[ "${V2_MULTI_HYDRO}" == "true" && "${JOB_PAYLOAD_EOS_BASE}" == "${EOS_BASE}" ]]; then
    HYDRO_BUNDLE="${PAYLOAD_DIR}/hydro_payloads_v2.tar"
    HYDRO_CHECKSUMS="${PAYLOAD_DIR}/hydro_payloads_v2.sha256"
    rm -f "${HYDRO_BUNDLE}" "${HYDRO_CHECKSUMS}"
    tar -cf "${HYDRO_BUNDLE}" --transform='s,^,hydro/C0-5/,' -C "${HYDRO_PAYLOAD_DIR}" .
    awk -F '\t' 'NR > 1 {print $8 "  " $6}' "${HYDRO_MANIFEST}" > "${HYDRO_CHECKSUMS}"
    scp -q -o BatchMode=yes "${HYDRO_BUNDLE}" "${HYDRO_CHECKSUMS}" "${CERN_REMOTE}:${TMP_REMOTE}/"
  fi
  scp -q -o BatchMode=yes "${SUPPORT}/build_runtime_job.sh" "${SUPPORT}/run_chunk_job.sh" "${CERN_REMOTE}:${AFS_WORK}/"
  "${CERNCTL}" run bash -lc "cp ${TMP_REMOTE}/mmli_source.tar.gz ${TMP_REMOTE}/runtime_payload.tar.gz ${PYTHIA_SOURCE} ${EOS_BASE}/payloads/ && chmod +x ${AFS_WORK}/build_runtime_job.sh ${AFS_WORK}/run_chunk_job.sh"
  if [[ "${V2_MULTI_HYDRO}" == "true" && "${JOB_PAYLOAD_EOS_BASE}" == "${EOS_BASE}" ]]; then
    "${CERNCTL}" run bash -lc "cd ${EOS_BASE}/payloads && tar -xf ${TMP_REMOTE}/hydro_payloads_v2.tar && sha256sum -c ${TMP_REMOTE}/hydro_payloads_v2.sha256"
  fi
fi

cat > "${WORK}/build_oo_no_moliere.sub" <<EOF_SUB
executable = build_runtime_job.sh
arguments =
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
output = build_oo_no_moliere.\$(ClusterId).out
error = build_oo_no_moliere.\$(ClusterId).err
log = log/build_oo_no_moliere.\$(ClusterId).log
environment = "EOS_BASE=${EOS_BASE} LHAPDF_CVMFS_VIEW=${LHAPDF_CVMFS_VIEW}"
+JobFlavour = "${BUILD_JOB_FLAVOUR}"
request_cpus = 1
request_memory = 4000
request_disk = 20000000
queue 1
EOF_SUB

cat > "${WORK}/oo_no_moliere_aa.sub" <<EOF_SUB
executable = run_chunk_job.sh
arguments = \$(chunk_id)
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_output_files = ""
# Physics archives and status records are uploaded to EOS by the wrapper.
output = /dev/null
error = /dev/null
log = log/oo_no_moliere_aa.\$(ClusterId).log
environment = "KIND=aa EOS_BASE=${EOS_BASE} PAYLOAD_EOS_BASE=${JOB_PAYLOAD_EOS_BASE} RUNTIME_PAYLOAD_EOS_BASE=${JOB_RUNTIME_PAYLOAD_EOS_BASE} SEED_OFFSET=${AA_SEED_OFFSET} EVENTS=${AA_EVENTS} RUN_NAME=${RUN_NAME} PTHAT_MIN=${PTHAT_MIN} PTHAT_MAX=${PTHAT_MAX} PDF_MODE=${AA_PDF_MODE} LHAPDF_SET=${AA_LHAPDF_SET} LHAPDF_CVMFS_VIEW=${LHAPDF_CVMFS_VIEW} AA_CENTRALITY_INDEX=${AA_CENTRALITY_INDEX} AA_TASK_MANIFEST=${AA_TASK_MANIFEST_RUNTIME} RUN_PREHYDRO_PAIR=${RUN_PREHYDRO_PAIR} NO_PREHYDRO_ALPHA=${NO_PREHYDRO_ALPHA} PREHYDRO_ALPHA=${PREHYDRO_ALPHA} BROADENING_K=${BROADENING_K} PREHYDRO_TAU_MIN=${PREHYDRO_TAU_MIN} PREHYDRO_TAU_GRID=${PREHYDRO_TAU_GRID} PREHYDRO_ETA_OVER_S=${PREHYDRO_ETA_OVER_S} PREHYDRO_EOS_FACTOR=${PREHYDRO_EOS_FACTOR} PREHYDRO_ATTRACTOR_TABLE=${PREHYDRO_ATTRACTOR_TABLE_RUNTIME} PREHYDRO_VISCOUS_ANCHOR=${PREHYDRO_VISCOUS_ANCHOR} STORE_PREHYDRO_TABLE=${STORE_PREHYDRO_TABLE} TOLERATE_CHUNK_FAILURE=${TOLERATE_CHUNK_FAILURES} TIMEOUT_S=${TIMEOUT_S}"
+JobFlavour = "${JOB_FLAVOUR}"
request_cpus = 1
request_memory = 4000
request_disk = 4000000
queue chunk_id from aa_chunk_ids.txt
EOF_SUB

cat > "${WORK}/oo_no_moliere_pp.sub" <<EOF_SUB
executable = run_chunk_job.sh
arguments = \$(chunk_id)
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_output_files = ""
# Physics archives and status records are uploaded to EOS by the wrapper.
output = /dev/null
error = /dev/null
log = log/oo_no_moliere_pp.\$(ClusterId).log
environment = "KIND=pp EOS_BASE=${EOS_BASE} SEED_OFFSET=${PP_SEED_OFFSET} EVENTS=${PP_EVENTS} RUN_NAME=${RUN_NAME} PTHAT_MIN=${PTHAT_MIN} PTHAT_MAX=${PTHAT_MAX} PDF_MODE=${PP_PDF_MODE} LHAPDF_SET=${PP_LHAPDF_SET} LHAPDF_CVMFS_VIEW=${LHAPDF_CVMFS_VIEW} RUN_PREHYDRO_PAIR=false TOLERATE_CHUNK_FAILURE=${TOLERATE_CHUNK_FAILURES} TIMEOUT_S=${TIMEOUT_S}"
+JobFlavour = "${JOB_FLAVOUR}"
request_cpus = 1
request_memory = 4000
request_disk = 4000000
queue chunk_id from pp_chunk_ids.txt
EOF_SUB

if [[ "${SUBMIT_PP}" == "true" ]]; then
cat > "${WORK}/campaign.dag" <<'EOF_DAG'
JOB BUILD build_oo_no_moliere.sub
JOB AA oo_no_moliere_aa.sub
JOB PP oo_no_moliere_pp.sub
PARENT BUILD CHILD AA PP
EOF_DAG
else
cat > "${WORK}/campaign.dag" <<'EOF_DAG'
JOB BUILD build_oo_no_moliere.sub
JOB AA oo_no_moliere_aa.sub
PARENT BUILD CHILD AA
EOF_DAG
fi

if [[ "${DRY_RUN}" == "1" || "${DRY_RUN}" == "true" || "${DRY_RUN}" == "TRUE" ]]; then
  echo "dry_run_note=skipping remote submit-file copy and condor_submit_dag"
else
  scp -q -o BatchMode=yes "${WORK}/build_oo_no_moliere.sub" "${WORK}/oo_no_moliere_aa.sub" "${WORK}/oo_no_moliere_pp.sub" "${WORK}/aa_chunk_ids.txt" "${WORK}/pp_chunk_ids.txt" "${WORK}/campaign.dag" "${CERN_REMOTE}:${AFS_WORK}/"
  "${CERNCTL}" run bash -lc "source /etc/profile.d/modules.sh 2>/dev/null || true; module load lxbatch/eossubmit >/dev/null 2>&1; cd ${AFS_WORK} && condor_submit_dag -force campaign.dag" | tee "${WORK}/logs/submit_cern_oo_no_moliere_validation.log"
fi
cat > "${WORK}/campaign_manifest.txt" <<EOF_MANIFEST
date=$(date -Is)
campaign=${CAMPAIGN}
eos_base=${EOS_BASE}
job_payload_eos_base=${JOB_PAYLOAD_EOS_BASE}
job_runtime_payload_eos_base=${JOB_RUNTIME_PAYLOAD_EOS_BASE}
afs_work=${AFS_WORK}
mmli_source=${MMLI_SOURCE_ROOT}
source_git_branch=${SOURCE_GIT_BRANCH}
source_git_commit=${SOURCE_GIT_COMMIT}
source_payload_sha256=${SOURCE_PAYLOAD_SHA256}
runtime_payload_sha256=${RUNTIME_PAYLOAD_SHA256}
aa_chunks=${AA_CHUNKS}
pp_chunks=${PP_CHUNKS}
aa_events_per_chunk=${AA_EVENTS}
pp_events_per_chunk=${PP_EVENTS}
aa_total_events_per_variant=${AA_TOTAL_EVENTS}
pp_total_events=${PP_TOTAL_EVENTS}
pp_submitted_events=${PP_SUBMITTED_EVENTS}
submit_pp=${SUBMIT_PP}
pp_reference_campaign=${PP_REFERENCE_CAMPAIGN:-none}
aa_seed_offset=${AA_SEED_OFFSET}
pp_seed_offset=${PP_SEED_OFFSET}
build_job_flavour=${BUILD_JOB_FLAVOUR}
job_flavour=${JOB_FLAVOUR}
timeout_s=${TIMEOUT_S}
tolerate_chunk_failures=${TOLERATE_CHUNK_FAILURES}
aa_centrality=${AA_CENTRALITY_LABEL}
aa_centrality_index=${AA_CENTRALITY_INDEX}
pthat_min=${PTHAT_MIN}
pthat_max=${PTHAT_MAX}
dry_run=${DRY_RUN}
pythia=8.315
beams_eCM=5360
aa_pdf_mode=${AA_PDF_MODE}
aa_lhapdf_set=${AA_LHAPDF_SET}
pp_pdf_mode=${PP_PDF_MODE}
pp_lhapdf_set=${PP_LHAPDF_SET:-none}
pythia_lhapdf_plugin=built_on_cern_from_Pythia8Plugins_LHAPDF6_h
lhapdf_cvmfs_view=${LHAPDF_CVMFS_VIEW}
aa_nPDF=EPPS21 O16 full nuclear PDF through LHAPDF6 when AA_PDF_MODE=lhapdf, or PYTHIA native O16 nPDF when AA_PDF_MODE=native_npdf
pp_reference_pdf=proton PDF:pSet=13 with no nPDF when PP_PDF_MODE=off
no_prehydro_alpha_kappa_sc=${NO_PREHYDRO_ALPHA}
prehydro_alpha_kappa_sc=${PREHYDRO_ALPHA}
broadening_K=${BROADENING_K}
do_elastic=false
do_lres=false
hydro_payload=${HYDRO_PAYLOAD_DESCRIPTION}
analysis_weight=${ANALYSIS_WEIGHT_DESCRIPTION}
run_prehydro_pair=${RUN_PREHYDRO_PAIR}
prehydro_tau_min=${PREHYDRO_TAU_MIN}
prehydro_tau_grid=${PREHYDRO_TAU_GRID}
prehydro_eta_over_s=${PREHYDRO_ETA_OVER_S}
prehydro_eos_factor=${PREHYDRO_EOS_FACTOR}
prehydro_attractor_table=${PREHYDRO_ATTRACTOR_TABLE}
prehydro_attractor_runtime=runtime/prehydro_attractor_table.dat
prehydro_attractor_sha256=${PREHYDRO_ATTRACTOR_SHA256}
prehydro_attractor_source=DOI:10.4119/unibi/2939684 QCD lambda=10 C_inf=0.87 curve
prehydro_viscous_anchor=${PREHYDRO_VISCOUS_ANCHOR}
prehydro_reference=arXiv:2509.19430v2 Eqs. (2)-(4), published QCD kinetic attractor table, three-flavor conformal EOS, natural-unit tau conversion, linear transverse pre-flow, and Bjorken longitudinal flow
paired_job_layout=AA task_NNNNN is no-prehydro baseline; task_NNNNN_prehydro is the same seed and hydro with 2509.19430 prehydro
v2_multi_hydro=${V2_MULTI_HYDRO}
aa_task_manifest=${AA_TASK_MANIFEST:-none}
aa_task_manifest_sha256=${AA_TASK_MANIFEST_SHA256}
hydro_manifest=${HYDRO_MANIFEST:-none}
hydro_manifest_sha256=${HYDRO_MANIFEST_SHA256}
hydro_payload_dir=${HYDRO_PAYLOAD_DIR:-none}
hydro_preparation_summary=${HYDRO_PREPARATION_SUMMARY:-none}
hydro_assignment=each 5000-task milestone block is allocated by largest remainder proportional to event Ncoll and deterministically shuffled
hydro_events=500 when v2_multi_hydro=true
store_prehydro_table=${STORE_PREHYDRO_TABLE}
EOF_MANIFEST

echo "Wrote ${WORK}/campaign_manifest.txt"

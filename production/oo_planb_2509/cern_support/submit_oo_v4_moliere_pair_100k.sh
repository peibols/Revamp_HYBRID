#!/usr/bin/env bash
set -Eeuo pipefail

ROOT="${ROOT:-/raid5/data/yjlee/hybrid_dev}"
SOURCE="${SOURCE:-${ROOT}/wt_main_moliere_lres_integration_clean}"
SUPPORT="${SOURCE}/production/oo_planb_2509/cern_support"
WORK="${WORK:-${ROOT}/test/oo5360_v4_moliere_pair_alpha0335_100k_20260715}"
CAMPAIGN="${CAMPAIGN:-hybrid_oo5360_c0_5_500hydro_moliere_paired_alpha0335_v4_100kAA_20260715}"
EOS_BASE="${EOS_BASE:-/eos/user/y/yjlee/${CAMPAIGN}}"
AFS_WORK="${AFS_WORK:-/afs/cern.ch/user/y/yjlee/oo_v4_moliere_pair_alpha0335_100k_20260715}"
CERN_REMOTE="${CERN_REMOTE:-lxplus}"
SCHEDD="${SCHEDD:-bigbird101.cern.ch}"
TMP_REMOTE="${TMP_REMOTE:-/tmp/yjlee_${CAMPAIGN}}"

SHARED_PAYLOAD_EOS_BASE="${SHARED_PAYLOAD_EOS_BASE:-/eos/user/y/yjlee/hybrid_oo5360_c0_5_500hydro_no_moliere_paired_public2509_planB_v2_50kAA_20260711}"
BASE_RUNTIME="${BASE_RUNTIME:-${ROOT}/test/oo5360_v2_500hydro_cont50k_to100k_20260712/payloads/mmli_runtime_alma9.tar.gz}"
TASK_MANIFEST="${TASK_MANIFEST:-${ROOT}/test/oo5360_v2_500hydro_cont50k_to100k_20260712/hydro_prepared/aa_task_manifest_combined_100k.tsv}"
HYDRO_MANIFEST="${HYDRO_MANIFEST:-${ROOT}/test/oo5360_v2_500hydro_50k_20260711/hydro_prepared/hydro_manifest.tsv}"
ATTRACTOR_TABLE="${ATTRACTOR_TABLE:-${SOURCE}/production/oo_planb_2509/reference_data/qcd_kinetic_attractor_lambda10_Cinf0p87.tsv}"
MOLIERE_TABLE_ARCHIVE="${MOLIERE_TABLE_ARCHIVE:-${ROOT}/test/tmp_redownload_tables_redo/a10_tables.zip}"
MOLIERE_TABLES_EOS_BASE="${MOLIERE_TABLES_EOS_BASE:-/eos/user/y/yjlee/hybrid_shared_moliere_a10_tables_20260715}"
MOLIERE_TABLES_KEY="${MOLIERE_TABLES_KEY:-payloads/a10_tables.zip}"

FULL_TASK_COUNT=100000
SUBMIT_TARGET="${SUBMIT_TARGET:-100000}"
TASK_ID_START="${TASK_ID_START:-0}"
MAX_JOBS_PER_SUBMIT="${MAX_JOBS_PER_SUBMIT:-10000}"
MAX_MATERIALIZE_PER_CLUSTER="${MAX_MATERIALIZE_PER_CLUSTER:-250}"
MAX_IDLE_PER_CLUSTER="${MAX_IDLE_PER_CLUSTER:-100}"
SEED_OFFSET="${SEED_OFFSET:-900000}"
ALPHA="${ALPHA:-0.335}"
BROADENING_K="${BROADENING_K:-15.0}"
TIMEOUT_S="${TIMEOUT_S:-72000}"
JOB_FLAVOUR="${JOB_FLAVOUR:-tomorrow}"
JOB_PRIORITY="${JOB_PRIORITY:-50}"
REQUEST_MEMORY_MB="${REQUEST_MEMORY_MB:-8000}"
REQUEST_DISK_KB="${REQUEST_DISK_KB:-8000000}"
STORE_PREHYDRO_TABLE="${STORE_PREHYDRO_TABLE:-false}"
DRY_RUN="${DRY_RUN:-false}"

EXPECTED_BASE_RUNTIME_SHA256="a208df01f48867ca21560f074ebd1cc2768c06ebd37ed7cbb30bd9dfca71c3e2"
EXPECTED_MAIN_SHA256="0d4c68b6ded87379e598e41670dd0bf8f8a39c234f0a9176b1b73ce6e2d88649"
EXPECTED_TASK_MANIFEST_SHA256="2317bd840ee8fffcc71f4b216a34ec4c22d5d7fa91eddb78dc81c1cd7f0a9a5e"
EXPECTED_HYDRO_MANIFEST_SHA256="293a1fb65a1c9641d14b413aa7e48f233235faecc9dfb5c0a4f014fa8dcca85b"
EXPECTED_ATTRACTOR_SHA256="1bea7289d3dc8ed95819eaa86cf4c489442a054c14aae47eff010cf45155eba0"
EXPECTED_MOLIERE_TABLES_SHA256="33974c43c35f2adac1634afc7c27c8d11985f46465f6b60e3ad4ccfa268f8b9a"

normalize_bool() {
  case "$2" in
    1|true|TRUE) printf -v "$1" '%s' true ;;
    0|false|FALSE) printf -v "$1" '%s' false ;;
    *) echo "$1 must be true or false, got: $2" >&2; exit 1 ;;
  esac
}
normalize_bool DRY_RUN "${DRY_RUN}"
normalize_bool STORE_PREHYDRO_TABLE "${STORE_PREHYDRO_TABLE}"

for name in SUBMIT_TARGET TASK_ID_START MAX_JOBS_PER_SUBMIT MAX_MATERIALIZE_PER_CLUSTER MAX_IDLE_PER_CLUSTER REQUEST_MEMORY_MB REQUEST_DISK_KB; do
  value="${!name}"
  if [[ ! "${value}" =~ ^[0-9]+$ ]]; then
    echo "${name} must be a nonnegative integer, got: ${value}" >&2
    exit 1
  fi
done
if (( SUBMIT_TARGET <= 0 || MAX_JOBS_PER_SUBMIT <= 0 || MAX_MATERIALIZE_PER_CLUSTER <= 0 || MAX_IDLE_PER_CLUSTER <= 0 || REQUEST_MEMORY_MB <= 0 || REQUEST_DISK_KB <= 0 )); then
  echo "submission target, queue limits, memory, and disk must be positive" >&2
  exit 1
fi
if (( MAX_IDLE_PER_CLUSTER > MAX_MATERIALIZE_PER_CLUSTER )); then
  echo "MAX_IDLE_PER_CLUSTER cannot exceed MAX_MATERIALIZE_PER_CLUSTER" >&2
  exit 1
fi
if (( TASK_ID_START + SUBMIT_TARGET > FULL_TASK_COUNT )); then
  echo "requested task interval exceeds ${FULL_TASK_COUNT} tasks" >&2
  exit 1
fi
if [[ ! "${JOB_PRIORITY}" =~ ^-?[0-9]+$ ]]; then
  echo "JOB_PRIORITY must be an integer" >&2
  exit 1
fi
python3 - "${ALPHA}" "${BROADENING_K}" <<'PY'
import math
import sys

alpha, broadening = map(float, sys.argv[1:])
if not math.isclose(alpha, 0.335, rel_tol=0.0, abs_tol=1e-12):
    raise SystemExit(f"ALPHA must be exactly 0.335, got {alpha}")
if not math.isclose(broadening, 15.0, rel_tol=0.0, abs_tol=1e-12):
    raise SystemExit(f"BROADENING_K must be exactly 15.0, got {broadening}")
PY

check_sha256() {
  local expected="$1"
  local path="$2"
  local actual
  actual="$(sha256sum "${path}" | awk '{print $1}')"
  if [[ "${actual}" != "${expected}" ]]; then
    echo "SHA256 mismatch for ${path}: ${actual}, expected ${expected}" >&2
    exit 1
  fi
}
for path in "${BASE_RUNTIME}" "${TASK_MANIFEST}" "${HYDRO_MANIFEST}" "${ATTRACTOR_TABLE}" "${MOLIERE_TABLE_ARCHIVE}"; do
  test -f "${path}"
done
if ! git -C "${SOURCE}" diff --quiet || ! git -C "${SOURCE}" diff --cached --quiet; then
  echo "Tracked source changes must be committed before campaign packaging" >&2
  exit 1
fi
check_sha256 "${EXPECTED_BASE_RUNTIME_SHA256}" "${BASE_RUNTIME}"
check_sha256 "${EXPECTED_TASK_MANIFEST_SHA256}" "${TASK_MANIFEST}"
check_sha256 "${EXPECTED_HYDRO_MANIFEST_SHA256}" "${HYDRO_MANIFEST}"
check_sha256 "${EXPECTED_ATTRACTOR_SHA256}" "${ATTRACTOR_TABLE}"
check_sha256 "${EXPECTED_MOLIERE_TABLES_SHA256}" "${MOLIERE_TABLE_ARCHIVE}"

python3 - "${TASK_MANIFEST}" "${FULL_TASK_COUNT}" "${SEED_OFFSET}" <<'PY'
import csv
from pathlib import Path
import sys

path = Path(sys.argv[1])
target, seed_offset = map(int, sys.argv[2:])
with path.open(newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
if len(rows) != target:
    raise SystemExit(f"expected {target} task rows, found {len(rows)}")
if [int(row["task_id"]) for row in rows] != list(range(target)):
    raise SystemExit("task IDs are not exactly 0..99999")
if [int(row["hard_seed"]) for row in rows] != list(range(seed_offset, seed_offset + target)):
    raise SystemExit("hard seeds do not match the V2 seed range")
if any(len(row["hydro_payload_sha256"]) != 64 for row in rows):
    raise SystemExit("task manifest contains an invalid hydro checksum")
print(f"validated {len(rows)} task identities")
PY

PAYLOAD_DIR="${WORK}/payloads"
RUNTIME_ROOT="${PAYLOAD_DIR}/runtime_moliere_pair"
RUNTIME_ARCHIVE="${PAYLOAD_DIR}/mmli_runtime_alma9.tar.gz"
mkdir -p "${PAYLOAD_DIR}" "${WORK}/logs" "${WORK}/cern_support"
rm -rf "${RUNTIME_ROOT}"
mkdir -p "${RUNTIME_ROOT}"
tar -xzf "${BASE_RUNTIME}" -C "${RUNTIME_ROOT}"
check_sha256 "${EXPECTED_MAIN_SHA256}" "${RUNTIME_ROOT}/bin/main"
cp "${SUPPORT}/run_oo_validation_chunk.py" "${RUNTIME_ROOT}/runtime/run_oo_validation_chunk.py"
cp "${TASK_MANIFEST}" "${RUNTIME_ROOT}/runtime/aa_task_manifest.tsv"
cp "${ATTRACTOR_TABLE}" "${RUNTIME_ROOT}/runtime/prehydro_attractor_table.dat"
cat > "${RUNTIME_ROOT}/runtime/moliere_pair_campaign.tsv" <<EOF_PROVENANCE
key	value
campaign	${CAMPAIGN}
no_prehydro_alpha	${ALPHA}
prehydro_alpha	${ALPHA}
moliere_mode	legacy_resolved
moliere_tables_sha256	${EXPECTED_MOLIERE_TABLES_SHA256}
task_manifest_sha256	${EXPECTED_TASK_MANIFEST_SHA256}
hydro_manifest_sha256	${EXPECTED_HYDRO_MANIFEST_SHA256}
EOF_PROVENANCE
tar -czf "${RUNTIME_ARCHIVE}" -C "${RUNTIME_ROOT}" .
RUNTIME_SHA256="$(sha256sum "${RUNTIME_ARCHIVE}" | awk '{print $1}')"
PACKAGED_MAIN_SHA256="$(tar -xOzf "${RUNTIME_ARCHIVE}" ./bin/main | sha256sum | awk '{print $1}')"
if [[ "${PACKAGED_MAIN_SHA256}" != "${EXPECTED_MAIN_SHA256}" ]]; then
  echo "Packaged executable changed: ${PACKAGED_MAIN_SHA256}" >&2
  exit 1
fi

cp "${SUPPORT}/run_chunk_job.sh" "${WORK}/cern_support/run_chunk_job.sh"
cp "${SUPPORT}/supervise_oo_moliere_pair.py" "${WORK}/cern_support/supervise_oo_moliere_pair.py"

python3 - "${WORK}" "${TASK_ID_START}" "${SUBMIT_TARGET}" "${MAX_JOBS_PER_SUBMIT}" <<'PY'
from pathlib import Path
import sys

work = Path(sys.argv[1])
start, target, batch = map(int, sys.argv[2:])
for old in work.glob("aa_chunk_ids_part_*.txt"):
    old.unlink()
for part, first in enumerate(range(start, start + target, batch)):
    stop = min(first + batch, start + target)
    (work / f"aa_chunk_ids_part_{part:02d}.txt").write_text(
        "".join(f"{task_id}\n" for task_id in range(first, stop))
    )
PY

PART_COUNT=$(( (SUBMIT_TARGET + MAX_JOBS_PER_SUBMIT - 1) / MAX_JOBS_PER_SUBMIT ))
RUN_NAME="runs/oo5360_c0_5_moliere_paired_alpha0335_v4"
for ((part = 0; part < PART_COUNT; ++part)); do
  printf -v part_tag '%02d' "${part}"
  cat > "${WORK}/oo_moliere_aa_part_${part_tag}.sub" <<EOF_SUB
executable = cern_support/run_chunk_job.sh
arguments = \$(chunk_id)
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_output_files = ""
output = /dev/null
error = /dev/null
log = /dev/null
environment = "KIND=aa EOS_BASE=${EOS_BASE} PAYLOAD_EOS_BASE=${SHARED_PAYLOAD_EOS_BASE} RUNTIME_PAYLOAD_EOS_BASE=${EOS_BASE} SEED_OFFSET=${SEED_OFFSET} EVENTS=1 RUN_NAME=${RUN_NAME} PTHAT_MIN=4 PTHAT_MAX=-1 PDF_MODE=lhapdf LHAPDF_SET=EPPS21nlo_CT18Anlo_O16/0 LHAPDF_CVMFS_VIEW=/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/lhapdf/6.5.3-3fa11/x86_64-centos7-gcc11-opt AA_CENTRALITY_INDEX=-1 AA_TASK_MANIFEST=runtime/aa_task_manifest.tsv RUN_PREHYDRO_PAIR=true RUN_PREHYDRO_ONLY=false NO_PREHYDRO_ALPHA=${ALPHA} PREHYDRO_ALPHA=${ALPHA} BROADENING_K=${BROADENING_K} PREHYDRO_TAU_MIN=0.24 PREHYDRO_ETA_OVER_S=0.12 PREHYDRO_EOS_FACTOR=15.626873635058152 PREHYDRO_ATTRACTOR_TABLE=runtime/prehydro_attractor_table.dat PREHYDRO_VISCOUS_ANCHOR=true STORE_PREHYDRO_TABLE=${STORE_PREHYDRO_TABLE} DO_MOLIERE=true MOLIERE_TABLES_EOS_BASE=${MOLIERE_TABLES_EOS_BASE} MOLIERE_TABLES_KEY=${MOLIERE_TABLES_KEY} MOLIERE_TABLES_SHA256=${EXPECTED_MOLIERE_TABLES_SHA256} TOLERATE_CHUNK_FAILURE=true TIMEOUT_S=${TIMEOUT_S}"
+JobFlavour = "${JOB_FLAVOUR}"
+JobBatchName = "OO_moliere_pair_a0335_${part_tag}"
+OOMoliereCampaign = "${CAMPAIGN}"
priority = ${JOB_PRIORITY}
request_cpus = 1
request_memory = ${REQUEST_MEMORY_MB}
request_disk = ${REQUEST_DISK_KB}
max_materialize = ${MAX_MATERIALIZE_PER_CLUSTER}
max_idle = ${MAX_IDLE_PER_CLUSTER}
queue chunk_id from aa_chunk_ids_part_${part_tag}.txt
EOF_SUB
done
cp "${WORK}/oo_moliere_aa_part_00.sub" "${WORK}/oo_moliere_aa.sub"
cp "${WORK}/oo_moliere_aa_part_00.sub" "${WORK}/oo_no_moliere_aa.sub"

SOURCE_COMMIT="$(git -C "${SOURCE}" rev-parse HEAD)"
RUNNER_SHA256="$(sha256sum "${SUPPORT}/run_oo_validation_chunk.py" | awk '{print $1}')"
WRAPPER_SHA256="$(sha256sum "${SUPPORT}/run_chunk_job.sh" | awk '{print $1}')"
cat > "${WORK}/campaign_manifest.txt" <<EOF_MANIFEST
date=$(date -Is)
campaign=${CAMPAIGN}
eos_base=${EOS_BASE}
afs_work=${AFS_WORK}
schedd=${SCHEDD}
source_commit=${SOURCE_COMMIT}
physics_main_sha256=${PACKAGED_MAIN_SHA256}
runtime_payload_sha256=${RUNTIME_SHA256}
runner_sha256=${RUNNER_SHA256}
wrapper_sha256=${WRAPPER_SHA256}
task_id_start=${TASK_ID_START}
submitted_tasks=${SUBMIT_TARGET}
full_task_count=${FULL_TASK_COUNT}
hard_seed_offset=${SEED_OFFSET}
events_per_job=1
variants_per_job=2
paired_job_layout=same hard seed and assigned hydro for no-prehydro and prehydro
no_prehydro_alpha=${ALPHA}
prehydro_alpha=${ALPHA}
broadening_K=${BROADENING_K}
prehydro_tau_min=0.24
prehydro_tau_hyd=0.4
do_elastic=true
do_lres=false
do_moliere=true
moliere_mode=legacy_resolved
moliere_unresolved_modes=false
compat_moliere_legacy_hydro=true
hadro_type=1
medium_response=true
moliere_tables_archive=${MOLIERE_TABLE_ARCHIVE}
moliere_tables_eos_base=${MOLIERE_TABLES_EOS_BASE}
moliere_tables_key=${MOLIERE_TABLES_KEY}
moliere_tables_sha256=${EXPECTED_MOLIERE_TABLES_SHA256}
task_manifest=${TASK_MANIFEST}
task_manifest_sha256=${EXPECTED_TASK_MANIFEST_SHA256}
hydro_manifest=${HYDRO_MANIFEST}
hydro_manifest_sha256=${EXPECTED_HYDRO_MANIFEST_SHA256}
hydro_payload_eos_base=${SHARED_PAYLOAD_EOS_BASE}
pythia=8.315
pthat_min=4
pthat_max=-1
bias2Selection=on
bias2SelectionPow=4
bias2SelectionRef=10
aa_pdf=EPPS21nlo_CT18Anlo_O16/0
job_flavour=${JOB_FLAVOUR}
job_priority=${JOB_PRIORITY}
request_memory_mb=${REQUEST_MEMORY_MB}
request_disk_kb=${REQUEST_DISK_KB}
timeout_s=${TIMEOUT_S}
submit_parts=${PART_COUNT}
jobs_per_submit_at_most=${MAX_JOBS_PER_SUBMIT}
max_materialize_per_cluster=${MAX_MATERIALIZE_PER_CLUSTER}
max_idle_per_cluster=${MAX_IDLE_PER_CLUSTER}
pp_generated=false
pp_reference=reuse the existing 1M pp reference because Moliere is a medium interaction
EOF_MANIFEST

echo "Prepared ${WORK}"
echo "runtime_sha256=${RUNTIME_SHA256}"
echo "main_sha256=${PACKAGED_MAIN_SHA256}"
if [[ "${DRY_RUN}" == "true" ]]; then
  echo "DRY_RUN=true: remote staging and submission skipped"
  exit 0
fi

ACTIVE="$({ ssh -o BatchMode=yes "${CERN_REMOTE}" \
  "condor_q -name ${SCHEDD} -constraint 'OOMoliereCampaign == \"${CAMPAIGN}\"' -af ClusterId" || true; } | sed '/^[[:space:]]*$/d')"
if [[ -n "${ACTIVE}" ]]; then
  echo "Campaign already has active jobs on ${SCHEDD}; refusing duplicate submission" >&2
  exit 1
fi

ssh -o BatchMode=yes "${CERN_REMOTE}" \
  "mkdir -p ${AFS_WORK}/cern_support ${EOS_BASE}/payloads ${EOS_BASE}/outputs/aa ${EOS_BASE}/status/aa ${TMP_REMOTE} ${MOLIERE_TABLES_EOS_BASE}/payloads"

REMOTE_TABLE_SHA="$({ ssh -o BatchMode=yes "${CERN_REMOTE}" \
  "test -f ${MOLIERE_TABLES_EOS_BASE}/${MOLIERE_TABLES_KEY} && sha256sum ${MOLIERE_TABLES_EOS_BASE}/${MOLIERE_TABLES_KEY}" || true; } | awk '{print $1}')"
if [[ "${REMOTE_TABLE_SHA}" != "${EXPECTED_MOLIERE_TABLES_SHA256}" ]]; then
  scp -q -o BatchMode=yes "${MOLIERE_TABLE_ARCHIVE}" \
    "${CERN_REMOTE}:${TMP_REMOTE}/a10_tables.zip"
  ssh -o BatchMode=yes "${CERN_REMOTE}" \
    "cp ${TMP_REMOTE}/a10_tables.zip ${MOLIERE_TABLES_EOS_BASE}/${MOLIERE_TABLES_KEY} && echo '${EXPECTED_MOLIERE_TABLES_SHA256}  ${MOLIERE_TABLES_EOS_BASE}/${MOLIERE_TABLES_KEY}' | sha256sum -c -"
fi

scp -q -o BatchMode=yes "${RUNTIME_ARCHIVE}" \
  "${CERN_REMOTE}:${TMP_REMOTE}/mmli_runtime_alma9.tar.gz"
scp -q -o BatchMode=yes "${WORK}/cern_support/run_chunk_job.sh" \
  "${CERN_REMOTE}:${AFS_WORK}/cern_support/run_chunk_job.sh"
scp -q -o BatchMode=yes "${WORK}"/aa_chunk_ids_part_*.txt \
  "${WORK}"/oo_moliere_aa_part_*.sub "${WORK}/oo_moliere_aa.sub" \
  "${WORK}/oo_no_moliere_aa.sub" "${CERN_REMOTE}:${AFS_WORK}/"
ssh -o BatchMode=yes "${CERN_REMOTE}" \
  "cp ${TMP_REMOTE}/mmli_runtime_alma9.tar.gz ${EOS_BASE}/payloads/mmli_runtime_alma9.tar.gz && chmod +x ${AFS_WORK}/cern_support/run_chunk_job.sh"

: > "${WORK}/submitted_clusters.tsv"
for ((part = 0; part < PART_COUNT; ++part)); do
  printf -v part_tag '%02d' "${part}"
  result="$(ssh -o BatchMode=yes "${CERN_REMOTE}" \
    "cd ${AFS_WORK} && condor_submit -name ${SCHEDD} oo_moliere_aa_part_${part_tag}.sub")"
  printf '%s\n' "${result}" | tee -a "${WORK}/logs/submit.log"
  cluster="$(printf '%s\n' "${result}" | sed -n 's/.*submitted to cluster \([0-9][0-9]*\).*/\1/p' | tail -n 1)"
  if [[ -z "${cluster}" ]]; then
    echo "Could not parse cluster ID for part ${part_tag}" >&2
    exit 1
  fi
  first=$((TASK_ID_START + part * MAX_JOBS_PER_SUBMIT))
  count=$((SUBMIT_TARGET - part * MAX_JOBS_PER_SUBMIT))
  if (( count > MAX_JOBS_PER_SUBMIT )); then count="${MAX_JOBS_PER_SUBMIT}"; fi
  printf '%s\t%s\t%s\t%s\n' "${part_tag}" "${cluster}" "${first}" "${count}" >> "${WORK}/submitted_clusters.tsv"
done

echo "Submitted ${SUBMIT_TARGET} paired Moliere jobs in ${PART_COUNT} cluster(s)"

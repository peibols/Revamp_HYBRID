#!/usr/bin/env bash
set -Eeuo pipefail

ROOT="${ROOT:-/raid5/data/yjlee/hybrid_dev}"
SOURCE="${SOURCE:-${ROOT}/wt_main_moliere_lres_integration_clean}"
WORK="${WORK:-${ROOT}/test/oo5360_v2_prehydro_only_alpha0335_100k_20260713}"
SUPPORT="${SOURCE}/production/oo_planb_2509/cern_support"
CAMPAIGN="${CAMPAIGN:-hybrid_oo5360_c0_5_500hydro_planB_alpha0335_only_v2_100kAA_20260713}"
EOS_BASE="${EOS_BASE:-/eos/user/y/yjlee/${CAMPAIGN}}"
AFS_WORK="${AFS_WORK:-/afs/cern.ch/user/y/yjlee/oo_v2_alpha0335_only_100k_20260713}"
CERN_REMOTE="${CERN_REMOTE:-lxplus}"
SCHEDD="${SCHEDD:-bigbird101.cern.ch}"
TMP_REMOTE="${TMP_REMOTE:-/tmp/yjlee_${CAMPAIGN}}"
SHARED_PAYLOAD_EOS_BASE="${SHARED_PAYLOAD_EOS_BASE:-/eos/user/y/yjlee/hybrid_oo5360_c0_5_500hydro_no_moliere_paired_public2509_planB_v2_50kAA_20260711}"
BASE_RUNTIME="${BASE_RUNTIME:-${ROOT}/test/oo5360_v2_500hydro_cont50k_to100k_20260712/payloads/mmli_runtime_alma9.tar.gz}"
TASK_MANIFEST="${TASK_MANIFEST:-${ROOT}/test/oo5360_v2_500hydro_cont50k_to100k_20260712/hydro_prepared/aa_task_manifest_combined_100k.tsv}"
HYDRO_MANIFEST="${HYDRO_MANIFEST:-${ROOT}/test/oo5360_v2_500hydro_50k_20260711/hydro_prepared/hydro_manifest.tsv}"
ATTRACTOR_TABLE="${ATTRACTOR_TABLE:-${SOURCE}/production/oo_planb_2509/reference_data/qcd_kinetic_attractor_lambda10_Cinf0p87.tsv}"
ALPHA="${ALPHA:-0.335}"
BROADENING_K="${BROADENING_K:-15.0}"
TARGET="${TARGET:-100000}"
MAX_JOBS_PER_SUBMIT="${MAX_JOBS_PER_SUBMIT:-10000}"
SEED_OFFSET="${SEED_OFFSET:-900000}"
TIMEOUT_S="${TIMEOUT_S:-72000}"
JOB_FLAVOUR="${JOB_FLAVOUR:-tomorrow}"
JOB_PRIORITY="${JOB_PRIORITY:-100}"
DRY_RUN="${DRY_RUN:-false}"
STORE_PREHYDRO_TABLE="${STORE_PREHYDRO_TABLE:-false}"

EXPECTED_BASE_RUNTIME_SHA256="a208df01f48867ca21560f074ebd1cc2768c06ebd37ed7cbb30bd9dfca71c3e2"
EXPECTED_MAIN_SHA256="0d4c68b6ded87379e598e41670dd0bf8f8a39c234f0a9176b1b73ce6e2d88649"
EXPECTED_TASK_MANIFEST_SHA256="2317bd840ee8fffcc71f4b216a34ec4c22d5d7fa91eddb78dc81c1cd7f0a9a5e"
EXPECTED_HYDRO_MANIFEST_SHA256="293a1fb65a1c9641d14b413aa7e48f233235faecc9dfb5c0a4f014fa8dcca85b"
EXPECTED_ATTRACTOR_SHA256="1bea7289d3dc8ed95819eaa86cf4c489442a054c14aae47eff010cf45155eba0"

case "${DRY_RUN}" in
  1|true|TRUE) DRY_RUN=true ;;
  0|false|FALSE) DRY_RUN=false ;;
  *) echo "DRY_RUN must be true or false" >&2; exit 1 ;;
esac
case "${STORE_PREHYDRO_TABLE}" in
  1|true|TRUE) STORE_PREHYDRO_TABLE=true ;;
  0|false|FALSE) STORE_PREHYDRO_TABLE=false ;;
  *) echo "STORE_PREHYDRO_TABLE must be true or false" >&2; exit 1 ;;
esac
if [[ "${TARGET}" != "100000" || "${MAX_JOBS_PER_SUBMIT}" != "10000" ]]; then
  echo "This campaign is pinned to TARGET=100000 and MAX_JOBS_PER_SUBMIT=10000" >&2
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

for path in "${BASE_RUNTIME}" "${TASK_MANIFEST}" "${HYDRO_MANIFEST}" "${ATTRACTOR_TABLE}"; do
  test -f "${path}"
done
if ! git -C "${SOURCE}" diff --quiet || ! git -C "${SOURCE}" diff --cached --quiet; then
  echo "Tracked source changes must be committed before campaign packaging" >&2
  exit 1
fi

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
check_sha256 "${EXPECTED_BASE_RUNTIME_SHA256}" "${BASE_RUNTIME}"
check_sha256 "${EXPECTED_TASK_MANIFEST_SHA256}" "${TASK_MANIFEST}"
check_sha256 "${EXPECTED_HYDRO_MANIFEST_SHA256}" "${HYDRO_MANIFEST}"
check_sha256 "${EXPECTED_ATTRACTOR_SHA256}" "${ATTRACTOR_TABLE}"

python3 - "${TASK_MANIFEST}" "${TARGET}" "${SEED_OFFSET}" <<'PY'
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
    raise SystemExit("hard seeds do not match the existing V2 seed range")
for row in rows:
    if len(row["hydro_payload_sha256"]) != 64:
        raise SystemExit(f"task {row['task_id']} has an invalid hydro checksum")
print(f"validated {len(rows)} task identities")
PY

PAYLOAD_DIR="${WORK}/payloads"
RUNTIME_ROOT="${PAYLOAD_DIR}/runtime_alpha0335"
RUNTIME_ARCHIVE="${PAYLOAD_DIR}/mmli_runtime_alma9.tar.gz"
mkdir -p "${PAYLOAD_DIR}" "${WORK}/logs" "${WORK}/cern_support"
rm -rf "${RUNTIME_ROOT}"
mkdir -p "${RUNTIME_ROOT}"
tar -xzf "${BASE_RUNTIME}" -C "${RUNTIME_ROOT}"
check_sha256 "${EXPECTED_MAIN_SHA256}" "${RUNTIME_ROOT}/bin/main"
cp "${SUPPORT}/run_oo_validation_chunk.py" "${RUNTIME_ROOT}/runtime/run_oo_validation_chunk.py"
cp "${TASK_MANIFEST}" "${RUNTIME_ROOT}/runtime/aa_task_manifest.tsv"
cp "${ATTRACTOR_TABLE}" "${RUNTIME_ROOT}/runtime/prehydro_attractor_table.dat"
cat > "${RUNTIME_ROOT}/runtime/prehydro_only_campaign.tsv" <<EOF_PROVENANCE
key	value
campaign	${CAMPAIGN}
prehydro_alpha	${ALPHA}
reference_main_sha256	${EXPECTED_MAIN_SHA256}
reference_runtime_sha256	${EXPECTED_BASE_RUNTIME_SHA256}
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
cp "${SUPPORT}/supervise_oo_v2.py" "${WORK}/cern_support/supervise_oo_v2.py"
cp "${SUPPORT}/supervise_oo_prehydro_only.py" "${WORK}/cern_support/supervise_oo_prehydro_only.py"

python3 - "${WORK}" "${TARGET}" "${MAX_JOBS_PER_SUBMIT}" <<'PY'
from pathlib import Path
import sys

work = Path(sys.argv[1])
target, batch = map(int, sys.argv[2:])
for part, start in enumerate(range(0, target, batch)):
    stop = min(start + batch, target)
    path = work / f"aa_chunk_ids_part_{part:02d}.txt"
    path.write_text("".join(f"{task_id}\n" for task_id in range(start, stop)))
PY

RUN_NAME="runs/oo5360_c0_5_planB_alpha0335_only_v2_100k"
for part in $(seq -f '%02g' 0 9); do
  cat > "${WORK}/oo_alpha0335_aa_part_${part}.sub" <<EOF_SUB
executable = cern_support/run_chunk_job.sh
arguments = \$(chunk_id)
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_output_files = ""
output = /dev/null
error = /dev/null
log = /dev/null
environment = "KIND=aa EOS_BASE=${EOS_BASE} PAYLOAD_EOS_BASE=${SHARED_PAYLOAD_EOS_BASE} RUNTIME_PAYLOAD_EOS_BASE=${EOS_BASE} SEED_OFFSET=${SEED_OFFSET} EVENTS=1 RUN_NAME=${RUN_NAME} PTHAT_MIN=4 PTHAT_MAX=-1 PDF_MODE=lhapdf LHAPDF_SET=EPPS21nlo_CT18Anlo_O16/0 LHAPDF_CVMFS_VIEW=/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/lhapdf/6.5.3-3fa11/x86_64-centos7-gcc11-opt AA_CENTRALITY_INDEX=-1 AA_TASK_MANIFEST=runtime/aa_task_manifest.tsv RUN_PREHYDRO_PAIR=false RUN_PREHYDRO_ONLY=true NO_PREHYDRO_ALPHA=0.37 PREHYDRO_ALPHA=${ALPHA} BROADENING_K=${BROADENING_K} PREHYDRO_TAU_MIN=0.24 PREHYDRO_ETA_OVER_S=0.12 PREHYDRO_EOS_FACTOR=15.626873635058152 PREHYDRO_ATTRACTOR_TABLE=runtime/prehydro_attractor_table.dat PREHYDRO_VISCOUS_ANCHOR=true STORE_PREHYDRO_TABLE=${STORE_PREHYDRO_TABLE} TOLERATE_CHUNK_FAILURE=true TIMEOUT_S=${TIMEOUT_S}"
+JobFlavour = "${JOB_FLAVOUR}"
+JobBatchName = "OO_alpha0335_only_100k_part_${part}"
+OOAlphaCampaign = "${CAMPAIGN}"
priority = ${JOB_PRIORITY}
request_cpus = 1
request_memory = 4000
request_disk = 4000000
queue chunk_id from aa_chunk_ids_part_${part}.txt
EOF_SUB
done

# The strict supervisor uses this canonical template for missing-task retries.
cp "${WORK}/oo_alpha0335_aa_part_00.sub" "${WORK}/oo_no_moliere_aa.sub"

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
run_mode=prehydro_only
prehydro_alpha=0.335
no_prehydro_generated=false
reference_no_prehydro_alpha=0.37
broadening_K=15.0
prehydro_tau_min=0.24
prehydro_tau_hyd=0.4
prehydro_eta_over_s=0.12
prehydro_viscous_anchor=true
do_elastic=false
do_lres=false
do_moliere=false
medium_response=true
pythia=8.315
pthat_min=4
pthat_max=-1
bias2Selection=on
bias2SelectionPow=4
bias2SelectionRef=10
aa_pdf=EPPS21nlo_CT18Anlo_O16/0
task_ids=0-99999
hard_seeds=900000-999999
events_per_job=1
target_events=100000
task_manifest=${TASK_MANIFEST}
task_manifest_sha256=${EXPECTED_TASK_MANIFEST_SHA256}
hydro_manifest=${HYDRO_MANIFEST}
hydro_manifest_sha256=${EXPECTED_HYDRO_MANIFEST_SHA256}
hydro_payload_eos_base=${SHARED_PAYLOAD_EOS_BASE}
reference_runtime=${BASE_RUNTIME}
reference_runtime_sha256=${EXPECTED_BASE_RUNTIME_SHA256}
reference_physics_source_commit=aa02143a27dcbabb7ee183fabac635feebb8861f
reference_main_sha256=${EXPECTED_MAIN_SHA256}
packaged_main_sha256=${PACKAGED_MAIN_SHA256}
runtime_payload_sha256=${RUNTIME_SHA256}
runner_sha256=${RUNNER_SHA256}
wrapper_sha256=${WRAPPER_SHA256}
attractor_sha256=${EXPECTED_ATTRACTOR_SHA256}
submit_parts=10
jobs_per_submit=10000
job_priority=${JOB_PRIORITY}
condor_logs=/dev/null
pp_generated=false
reference_v2_first50k_eos=${SHARED_PAYLOAD_EOS_BASE}
reference_v2_second50k_eos=/eos/user/y/yjlee/hybrid_oo5360_c0_5_500hydro_no_moliere_paired_public2509_planB_v2_cont50k_to100k_20260712
analysis_join=join by task_id against the existing V2 no-prehydro and alpha=0.37 archives; require exact seed, hydro metadata, payload SHA, event weight, and hard marker before RAA, jet spectra, or substructure comparison
EOF_MANIFEST

echo "Prepared ${WORK}"
echo "runtime_sha256=${RUNTIME_SHA256}"
echo "main_sha256=${PACKAGED_MAIN_SHA256}"
if [[ "${DRY_RUN}" == "true" ]]; then
  echo "DRY_RUN=true: remote staging and submission skipped"
  exit 0
fi

ACTIVE="$({ ssh -o BatchMode=yes "${CERN_REMOTE}" \
  "condor_q -name ${SCHEDD} -constraint 'OOAlphaCampaign == \"${CAMPAIGN}\"' -af ClusterId" || true; } | sed '/^[[:space:]]*$/d')"
if [[ -n "${ACTIVE}" ]]; then
  echo "Campaign already has active jobs on ${SCHEDD}; refusing duplicate submission" >&2
  exit 1
fi

ssh -o BatchMode=yes "${CERN_REMOTE}" \
  "mkdir -p ${AFS_WORK}/cern_support ${EOS_BASE}/payloads ${EOS_BASE}/outputs/aa ${EOS_BASE}/status/aa ${TMP_REMOTE}"
scp -q -o BatchMode=yes "${RUNTIME_ARCHIVE}" \
  "${CERN_REMOTE}:${TMP_REMOTE}/mmli_runtime_alma9.tar.gz"
scp -q -o BatchMode=yes \
  "${WORK}/cern_support/run_chunk_job.sh" \
  "${CERN_REMOTE}:${AFS_WORK}/cern_support/run_chunk_job.sh"
scp -q -o BatchMode=yes "${WORK}"/aa_chunk_ids_part_*.txt \
  "${WORK}"/oo_alpha0335_aa_part_*.sub \
  "${WORK}/oo_no_moliere_aa.sub" \
  "${CERN_REMOTE}:${AFS_WORK}/"
ssh -o BatchMode=yes "${CERN_REMOTE}" \
  "cp ${TMP_REMOTE}/mmli_runtime_alma9.tar.gz ${EOS_BASE}/payloads/mmli_runtime_alma9.tar.gz && chmod +x ${AFS_WORK}/cern_support/run_chunk_job.sh"

: > "${WORK}/submitted_clusters.tsv"
for part in $(seq -f '%02g' 0 9); do
  result="$(ssh -o BatchMode=yes "${CERN_REMOTE}" \
    "cd ${AFS_WORK} && condor_submit -name ${SCHEDD} oo_alpha0335_aa_part_${part}.sub")"
  printf '%s\n' "${result}" | tee -a "${WORK}/logs/submit_100k.log"
  cluster="$(printf '%s\n' "${result}" | sed -n 's/.*submitted to cluster \([0-9][0-9]*\).*/\1/p' | tail -n 1)"
  if [[ -z "${cluster}" ]]; then
    echo "Could not parse cluster ID for part ${part}" >&2
    exit 1
  fi
  printf '%s\t%s\t10000\n' "${part}" "${cluster}" >> "${WORK}/submitted_clusters.tsv"
done

echo "Submitted ${TARGET} prehydro-only jobs in 10 clusters"

#!/usr/bin/env bash
set -euo pipefail
ROOT="${ROOT:-/raid5/data/yjlee/hybrid_dev}"
WORK="${WORK:-${ROOT}/test/oo5360_no_moliere_raa_20260606}"
SUPPORT="${WORK}/cern_support"
CERNCTL="${CERNCTL:-/data/yjlee/cernLxplus/cernctl}"
CERN_REMOTE="${CERN_REMOTE:-lxplus}"
CAMPAIGN="${CAMPAIGN:-hybrid_oo5360_c0_5_no_moliere_paired_prehydro_pthat4_unbounded_20260709}"
EOS_BASE="${EOS_BASE:-/eos/user/y/yjlee/${CAMPAIGN}}"
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
PREHYDRO_TAU_MIN="${PREHYDRO_TAU_MIN:-0.24}"
PREHYDRO_TAU_GRID="${PREHYDRO_TAU_GRID:-0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.399}"
PREHYDRO_ETA_OVER_S="${PREHYDRO_ETA_OVER_S:-0.12}"
PREHYDRO_EOS_FACTOR="${PREHYDRO_EOS_FACTOR:-15.626873635058152}"
DEFAULT_PREHYDRO_ATTRACTOR_TABLE="${WORK}/reference_data/qcd_kinetic_attractor_lambda10_Cinf0p87.tsv"
PREHYDRO_ATTRACTOR_TABLE="${PREHYDRO_ATTRACTOR_TABLE:-${DEFAULT_PREHYDRO_ATTRACTOR_TABLE}}"
PREHYDRO_VISCOUS_ANCHOR="${PREHYDRO_VISCOUS_ANCHOR:-true}"

for VALUE_NAME in AA_CHUNKS PP_CHUNKS AA_EVENTS PP_EVENTS; do
  VALUE="${!VALUE_NAME}"
  if [[ ! "${VALUE}" =~ ^[1-9][0-9]*$ ]]; then
    echo "${VALUE_NAME} must be a positive integer, got: ${VALUE}" >&2
    exit 1
  fi
done
case "${SUBMIT_PP}" in
  1|true|TRUE) SUBMIT_PP=true ;;
  0|false|FALSE) SUBMIT_PP=false ;;
  *) echo "SUBMIT_PP must be true or false, got: ${SUBMIT_PP}" >&2; exit 1 ;;
esac
if [[ "${SUBMIT_PP}" == "false" && -z "${PP_REFERENCE_CAMPAIGN}" ]]; then
  echo "PP_REFERENCE_CAMPAIGN is required when SUBMIT_PP=false" >&2
  exit 1
fi
if [[ "${AA_EVENTS}" != "1" ]]; then
  echo "AA_EVENTS must be 1: each AA job runs one hard event in both paired variants" >&2
  exit 1
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
  PAIR_TAG="pairedPrehydroTau${PREHYDRO_TAU_TAG}"
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
echo "afs_work=${AFS_WORK}"
echo "mmli_source=${MMLI_SOURCE_ROOT}"
echo "dry_run=${DRY_RUN}"
echo "aa_centrality=${AA_CENTRALITY_LABEL} aa_centrality_index=${AA_CENTRALITY_INDEX}"
echo "aa_chunks=${AA_CHUNKS} aa_events_per_chunk=${AA_EVENTS} aa_total_events=${AA_TOTAL_EVENTS}"
echo "submit_pp=${SUBMIT_PP} pp_chunks=${PP_CHUNKS} pp_events_per_chunk=${PP_EVENTS} pp_total_events=${PP_TOTAL_EVENTS} pp_reference_campaign=${PP_REFERENCE_CAMPAIGN:-none}"
echo "job_flavour=${JOB_FLAVOUR} timeout_s=${TIMEOUT_S} pthat_min=${PTHAT_MIN} pthat_max=${PTHAT_MAX}"
echo "aa_pdf_mode=${AA_PDF_MODE} aa_lhapdf_set=${AA_LHAPDF_SET}"
echo "pp_pdf_mode=${PP_PDF_MODE} pp_lhapdf_set=${PP_LHAPDF_SET:-none}"
echo "run_prehydro_pair=${RUN_PREHYDRO_PAIR} prehydro_tau_min=${PREHYDRO_TAU_MIN} prehydro_tau_grid=${PREHYDRO_TAU_GRID} prehydro_attractor_table=${PREHYDRO_ATTRACTOR_TABLE}"
echo "prehydro_attractor_sha256=${PREHYDRO_ATTRACTOR_SHA256}"

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
tar -czf "${PAYLOAD_DIR}/runtime_payload/runtime/staged_hydro.tar.gz" -C "${WORK}" staged_hydro hydro_manifest.tsv
(
  cd "${PAYLOAD_DIR}/runtime_payload/runtime"
  tar -xzf staged_hydro.tar.gz
  rm staged_hydro.tar.gz
)
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
  scp -q -o BatchMode=yes "${SUPPORT}/build_runtime_job.sh" "${SUPPORT}/run_chunk_job.sh" "${CERN_REMOTE}:${AFS_WORK}/"
  "${CERNCTL}" run bash -lc "cp ${TMP_REMOTE}/mmli_source.tar.gz ${TMP_REMOTE}/runtime_payload.tar.gz ${PYTHIA_SOURCE} ${EOS_BASE}/payloads/ && chmod +x ${AFS_WORK}/build_runtime_job.sh ${AFS_WORK}/run_chunk_job.sh"
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
output = log/oo_no_moliere_aa.\$(ClusterId).\$(ProcId).out
error = log/oo_no_moliere_aa.\$(ClusterId).\$(ProcId).err
log = log/oo_no_moliere_aa.\$(ClusterId).log
environment = "KIND=aa EOS_BASE=${EOS_BASE} SEED_OFFSET=${AA_SEED_OFFSET} EVENTS=${AA_EVENTS} RUN_NAME=${RUN_NAME} PTHAT_MIN=${PTHAT_MIN} PTHAT_MAX=${PTHAT_MAX} PDF_MODE=${AA_PDF_MODE} LHAPDF_SET=${AA_LHAPDF_SET} LHAPDF_CVMFS_VIEW=${LHAPDF_CVMFS_VIEW} AA_CENTRALITY_INDEX=${AA_CENTRALITY_INDEX} RUN_PREHYDRO_PAIR=${RUN_PREHYDRO_PAIR} PREHYDRO_TAU_MIN=${PREHYDRO_TAU_MIN} PREHYDRO_TAU_GRID=${PREHYDRO_TAU_GRID} PREHYDRO_ETA_OVER_S=${PREHYDRO_ETA_OVER_S} PREHYDRO_EOS_FACTOR=${PREHYDRO_EOS_FACTOR} PREHYDRO_ATTRACTOR_TABLE=${PREHYDRO_ATTRACTOR_TABLE_RUNTIME} PREHYDRO_VISCOUS_ANCHOR=${PREHYDRO_VISCOUS_ANCHOR} TIMEOUT_S=${TIMEOUT_S}"
+JobFlavour = "${JOB_FLAVOUR}"
request_cpus = 1
request_memory = 4000
request_disk = 20000000
queue chunk_id from aa_chunk_ids.txt
EOF_SUB

cat > "${WORK}/oo_no_moliere_pp.sub" <<EOF_SUB
executable = run_chunk_job.sh
arguments = \$(chunk_id)
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
output = log/oo_no_moliere_pp.\$(ClusterId).\$(ProcId).out
error = log/oo_no_moliere_pp.\$(ClusterId).\$(ProcId).err
log = log/oo_no_moliere_pp.\$(ClusterId).log
environment = "KIND=pp EOS_BASE=${EOS_BASE} SEED_OFFSET=${PP_SEED_OFFSET} EVENTS=${PP_EVENTS} RUN_NAME=${RUN_NAME} PTHAT_MIN=${PTHAT_MIN} PTHAT_MAX=${PTHAT_MAX} PDF_MODE=${PP_PDF_MODE} LHAPDF_SET=${PP_LHAPDF_SET} LHAPDF_CVMFS_VIEW=${LHAPDF_CVMFS_VIEW} RUN_PREHYDRO_PAIR=false TIMEOUT_S=${TIMEOUT_S}"
+JobFlavour = "${JOB_FLAVOUR}"
request_cpus = 1
request_memory = 4000
request_disk = 20000000
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
alpha_kappa_sc=0.37
broadening_K=15.0
do_elastic=false
do_lres=false
hydro_payload=staged_hydro one event per Zenodo centrality bin C0-5...C90-100
analysis_weight=centrality_width*ncoll
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
EOF_MANIFEST

echo "Wrote ${WORK}/campaign_manifest.txt"

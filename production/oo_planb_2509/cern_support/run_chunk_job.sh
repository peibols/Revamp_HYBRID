#!/usr/bin/env bash
set -Eeuo pipefail
EOS_BASE="${EOS_BASE:?EOS_BASE is required}"
KIND="${KIND:?KIND is required: aa or pp}"
SEED_OFFSET="${SEED_OFFSET:?SEED_OFFSET is required}"
EVENTS="${EVENTS:?EVENTS is required}"
RUN_NAME="${RUN_NAME:?RUN_NAME is required}"
TIMEOUT_S="${TIMEOUT_S:-7200}"
PTHAT_MIN="${PTHAT_MIN:-4}"
PTHAT_MAX="${PTHAT_MAX:--1}"
if [[ -z "${PDF_MODE:-}" ]]; then
  if [[ "${KIND}" == "aa" ]]; then
    PDF_MODE="native_npdf"
  else
    PDF_MODE="off"
  fi
fi
LHAPDF_SET="${LHAPDF_SET:-EPPS21nlo_CT18Anlo_O16/0}"
AA_CENTRALITY_INDEX="${AA_CENTRALITY_INDEX:--1}"
RUN_PREHYDRO_PAIR="${RUN_PREHYDRO_PAIR:-false}"
PREHYDRO_TAU_MIN="${PREHYDRO_TAU_MIN:-0.24}"
PREHYDRO_TAU_GRID="${PREHYDRO_TAU_GRID:-0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.399}"
PREHYDRO_ETA_OVER_S="${PREHYDRO_ETA_OVER_S:-0.12}"
PREHYDRO_EOS_FACTOR="${PREHYDRO_EOS_FACTOR:-15.626873635058152}"
PREHYDRO_ATTRACTOR_TABLE="${PREHYDRO_ATTRACTOR_TABLE:-runtime/prehydro_attractor_table.dat}"
PREHYDRO_VISCOUS_ANCHOR="${PREHYDRO_VISCOUS_ANCHOR:-true}"
TASK_ID="${1:-${_CONDOR_PROCNO:-0}}"
if [[ "${RUN_PREHYDRO_PAIR}" == "1" || "${RUN_PREHYDRO_PAIR}" == "true" || "${RUN_PREHYDRO_PAIR}" == "TRUE" ]]; then
  if [[ -z "${PREHYDRO_ATTRACTOR_TABLE}" ]]; then
    echo "PREHYDRO_ATTRACTOR_TABLE is required for paired prehydro production" >&2
    exit 1
  fi
fi
WORK="$PWD/work"
INITIAL_DIR="$PWD"
fetch() { xrdcp -f "root://eosuser.cern.ch/${EOS_BASE}/$1" "$2"; }
put() { xrdcp -f "$1" "root://eosuser.cern.ch/${EOS_BASE}/$2"; }
finalize() {
  rc=$?; set +e
  echo "kind=$KIND" > "$INITIAL_DIR/chunk_status.txt"
  echo "task_id=$TASK_ID" >> "$INITIAL_DIR/chunk_status.txt"
  echo "status=$([[ $rc -eq 0 ]] && echo success || echo failed)" >> "$INITIAL_DIR/chunk_status.txt"
  echo "exit_code=$rc" >> "$INITIAL_DIR/chunk_status.txt"
  echo "date=$(date -Is)" >> "$INITIAL_DIR/chunk_status.txt"
  if [[ -d "$WORK/$RUN_NAME" ]]; then
    (
      cd "$WORK"
      find "$RUN_NAME" -type f \( \
        -name HYBRID_Hadrons.out -o \
        -name summary.tsv -o \
        -name stdout.log -o \
        -name hybrid_input.dat -o \
        -name setup_pythia.cmnd -o \
        -name prehydro_table.tsv -o \
        -name 'task_*_pair_summary.tsv' \
      \) -print0 | tar --null --ignore-failed-read -czf "$INITIAL_DIR/chunk_output.tar.gz" --files-from -
    ) 2>/dev/null || tar -czf "$INITIAL_DIR/chunk_output.tar.gz" -C "$INITIAL_DIR" chunk_status.txt
  else
    tar -czf "$INITIAL_DIR/chunk_output.tar.gz" -C "$INITIAL_DIR" chunk_status.txt
  fi
  put "$INITIAL_DIR/chunk_status.txt" "status/${KIND}/chunk_${TASK_ID}.txt" >/dev/null 2>&1 || true
  put "$INITIAL_DIR/chunk_output.tar.gz" "outputs/${KIND}/chunk_${TASK_ID}.tar.gz" >/dev/null 2>&1 || true
  exit $rc
}
trap finalize EXIT
mkdir -p "$WORK"
cd "$WORK"
fetch payloads/pythia8315_alma9_install.tar.gz pythia8315_alma9_install.tar.gz
fetch payloads/mmli_runtime_alma9.tar.gz mmli_runtime_alma9.tar.gz
tar -xzf pythia8315_alma9_install.tar.gz
tar -xzf mmli_runtime_alma9.tar.gz
PYTHIA_ROOT="$WORK/pythia8315_alma9_install"
export PYTHIA8DATA="$PYTHIA_ROOT/share/Pythia8/xmldoc"
export LD_LIBRARY_PATH="$PYTHIA_ROOT/lib:${LD_LIBRARY_PATH:-}"
if [[ -d "$WORK/runtime/root_lib" ]]; then
  export LD_LIBRARY_PATH="$WORK/runtime/root_lib:$LD_LIBRARY_PATH"
fi
ROOT_LIBDIR="${ROOT_LIBDIR:-}"
if [[ -z "$ROOT_LIBDIR" && -f "$WORK/runtime/root_libdir.txt" ]]; then
  ROOT_LIBDIR="$(cat "$WORK/runtime/root_libdir.txt")"
fi
if [[ -z "$ROOT_LIBDIR" ]] && command -v root-config >/dev/null 2>&1; then
  ROOT_LIBDIR="$(root-config --libdir)"
fi
if [[ -z "$ROOT_LIBDIR" && -d /usr/lib64/root ]]; then
  ROOT_LIBDIR="/usr/lib64/root"
fi
if [[ -n "$ROOT_LIBDIR" && -d "$ROOT_LIBDIR" ]]; then
  export LD_LIBRARY_PATH="$ROOT_LIBDIR:$LD_LIBRARY_PATH"
fi
LHAPDF_VIEW="${LHAPDF_CVMFS_VIEW:-}"
if [[ -z "$LHAPDF_VIEW" && -f "$WORK/runtime/lhapdf_cvmfs_view.txt" ]]; then
  LHAPDF_VIEW="$(cat "$WORK/runtime/lhapdf_cvmfs_view.txt")"
fi
if [[ -n "$LHAPDF_VIEW" ]]; then
  export LD_LIBRARY_PATH="$LHAPDF_VIEW/lib:$LD_LIBRARY_PATH"
fi
if [[ -d "$WORK/runtime/lhapdf_plugin" ]]; then
  export LD_LIBRARY_PATH="$WORK/runtime/lhapdf_plugin:$LD_LIBRARY_PATH"
fi
if [[ -d "$WORK/runtime/lhapdf_data" ]]; then
  export LHAPDF_DATA_PATH="$WORK/runtime/lhapdf_data${LHAPDF_DATA_PATH:+:$LHAPDF_DATA_PATH}"
fi
if [[ -n "$LHAPDF_VIEW" && -d "$LHAPDF_VIEW/share/LHAPDF" ]]; then
  export LHAPDF_DATA_PATH="${LHAPDF_DATA_PATH:+$LHAPDF_DATA_PATH:}$LHAPDF_VIEW/share/LHAPDF"
fi
export MMLI_BIN="$WORK/bin/main"
export MMLI_PYTHIA_TEMPLATE="$WORK/runtime/setup_pythia.cmnd"
export OO_HYDRO_ROOT="$WORK/runtime/staged_hydro"
args=(
  runtime/run_oo_validation_chunk.py
  --kind "$KIND" \
  --task-id "$TASK_ID" \
  --seed-offset "$SEED_OFFSET" \
  --events "$EVENTS" \
  --pthat-min "$PTHAT_MIN" \
  --pthat-max "$PTHAT_MAX" \
  --pdf-mode "$PDF_MODE" \
  --lhapdf-set "$LHAPDF_SET" \
  --timeout-s "$TIMEOUT_S" \
  --run-name "$RUN_NAME" \
  --aa-centrality-index "$AA_CENTRALITY_INDEX"
)
if [[ "$RUN_PREHYDRO_PAIR" == "1" || "$RUN_PREHYDRO_PAIR" == "true" || "$RUN_PREHYDRO_PAIR" == "TRUE" ]]; then
  args+=(--run-prehydro-pair)
fi
args+=(
  --prehydro-tau-min "$PREHYDRO_TAU_MIN"
  --prehydro-tau-grid "$PREHYDRO_TAU_GRID"
  --prehydro-eta-over-s "$PREHYDRO_ETA_OVER_S"
  --prehydro-eos-factor "$PREHYDRO_EOS_FACTOR"
)
if [[ -n "$PREHYDRO_ATTRACTOR_TABLE" ]]; then
  args+=(--prehydro-attractor-table "$PREHYDRO_ATTRACTOR_TABLE")
fi
if [[ "$PREHYDRO_VISCOUS_ANCHOR" == "0" || "$PREHYDRO_VISCOUS_ANCHOR" == "false" || "$PREHYDRO_VISCOUS_ANCHOR" == "FALSE" ]]; then
  args+=(--no-prehydro-viscous-anchor)
elif [[ "$PREHYDRO_VISCOUS_ANCHOR" == "1" || "$PREHYDRO_VISCOUS_ANCHOR" == "true" || "$PREHYDRO_VISCOUS_ANCHOR" == "TRUE" ]]; then
  args+=(--prehydro-viscous-anchor)
fi
python3 "${args[@]}"

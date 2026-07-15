#!/usr/bin/env bash
set -Eeuo pipefail
EOS_BASE="${EOS_BASE:?EOS_BASE is required}"
PAYLOAD_EOS_BASE="${PAYLOAD_EOS_BASE:-$EOS_BASE}"
RUNTIME_PAYLOAD_EOS_BASE="${RUNTIME_PAYLOAD_EOS_BASE:-$PAYLOAD_EOS_BASE}"
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
RUN_PREHYDRO_ONLY="${RUN_PREHYDRO_ONLY:-false}"
NO_PREHYDRO_ALPHA="${NO_PREHYDRO_ALPHA:-0.37}"
PREHYDRO_ALPHA="${PREHYDRO_ALPHA:-0.37}"
BROADENING_K="${BROADENING_K:-15.0}"
PREHYDRO_TAU_MIN="${PREHYDRO_TAU_MIN:-0.24}"
PREHYDRO_TAU_GRID="${PREHYDRO_TAU_GRID:-0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.399}"
PREHYDRO_ETA_OVER_S="${PREHYDRO_ETA_OVER_S:-0.12}"
PREHYDRO_EOS_FACTOR="${PREHYDRO_EOS_FACTOR:-15.626873635058152}"
PREHYDRO_ATTRACTOR_TABLE="${PREHYDRO_ATTRACTOR_TABLE:-runtime/prehydro_attractor_table.dat}"
PREHYDRO_VISCOUS_ANCHOR="${PREHYDRO_VISCOUS_ANCHOR:-true}"
TOLERATE_CHUNK_FAILURE="${TOLERATE_CHUNK_FAILURE:-false}"
AA_TASK_MANIFEST="${AA_TASK_MANIFEST:-}"
STORE_PREHYDRO_TABLE="${STORE_PREHYDRO_TABLE:-true}"
DO_MOLIERE="${DO_MOLIERE:-false}"
MOLIERE_TABLES_EOS_BASE="${MOLIERE_TABLES_EOS_BASE:-}"
MOLIERE_TABLES_KEY="${MOLIERE_TABLES_KEY:-}"
MOLIERE_TABLES_SHA256="${MOLIERE_TABLES_SHA256:-}"
TASK_ID="${1:-${_CONDOR_PROCNO:-0}}"
HYDRO_SLOT=""
HYDRO_EVENT_ID=""
HYDRO_NCOLL=""
HYDRO_DIR=""
HYDRO_PAYLOAD_KEY=""
HYDRO_PAYLOAD_SHA256=""
case "${TOLERATE_CHUNK_FAILURE}" in
  1|true|TRUE) TOLERATE_CHUNK_FAILURE=true ;;
  0|false|FALSE) TOLERATE_CHUNK_FAILURE=false ;;
  *) echo "TOLERATE_CHUNK_FAILURE must be true or false" >&2; exit 1 ;;
esac
case "${RUN_PREHYDRO_PAIR}" in
  1|true|TRUE) RUN_PREHYDRO_PAIR=true ;;
  0|false|FALSE) RUN_PREHYDRO_PAIR=false ;;
  *) echo "RUN_PREHYDRO_PAIR must be true or false" >&2; exit 1 ;;
esac
case "${RUN_PREHYDRO_ONLY}" in
  1|true|TRUE) RUN_PREHYDRO_ONLY=true ;;
  0|false|FALSE) RUN_PREHYDRO_ONLY=false ;;
  *) echo "RUN_PREHYDRO_ONLY must be true or false" >&2; exit 1 ;;
esac
case "${DO_MOLIERE}" in
  1|true|TRUE) DO_MOLIERE=true ;;
  0|false|FALSE) DO_MOLIERE=false ;;
  *) echo "DO_MOLIERE must be true or false" >&2; exit 1 ;;
esac
if [[ "${RUN_PREHYDRO_PAIR}" == "true" && "${RUN_PREHYDRO_ONLY}" == "true" ]]; then
  echo "RUN_PREHYDRO_PAIR and RUN_PREHYDRO_ONLY are mutually exclusive" >&2
  exit 1
fi
if [[ "${DO_MOLIERE}" == "true" ]]; then
  if [[ "${KIND}" != "aa" ]]; then
    echo "Moliere scattering is only supported for AA jobs" >&2
    exit 1
  fi
  if [[ -z "${MOLIERE_TABLES_EOS_BASE}" || -z "${MOLIERE_TABLES_KEY}" || ! "${MOLIERE_TABLES_SHA256}" =~ ^[0-9a-f]{64}$ ]]; then
    echo "MOLIERE_TABLES_EOS_BASE, MOLIERE_TABLES_KEY, and a lowercase SHA256 are required" >&2
    exit 1
  fi
fi
if [[ "${RUN_PREHYDRO_PAIR}" == "true" || "${RUN_PREHYDRO_ONLY}" == "true" ]]; then
  if [[ -z "${PREHYDRO_ATTRACTOR_TABLE}" ]]; then
    echo "PREHYDRO_ATTRACTOR_TABLE is required for prehydro production" >&2
    exit 1
  fi
fi
WORK="$PWD/work"
INITIAL_DIR="$PWD"
fetch_from() {
  local base="$1"
  local key="$2"
  local destination="$3"
  xrdcp -f "root://eosuser.cern.ch/${base}/$key" "$destination"
}
put() { xrdcp -f "$1" "root://eosuser.cern.ch/${EOS_BASE}/$2"; }
finalize() {
  rc=$?; set +e
  echo "kind=$KIND" > "$INITIAL_DIR/chunk_status.txt"
  echo "task_id=$TASK_ID" >> "$INITIAL_DIR/chunk_status.txt"
  echo "status=$([[ $rc -eq 0 ]] && echo success || echo failed)" >> "$INITIAL_DIR/chunk_status.txt"
  echo "exit_code=$rc" >> "$INITIAL_DIR/chunk_status.txt"
  echo "date=$(date -Is)" >> "$INITIAL_DIR/chunk_status.txt"
  echo "run_prehydro_pair=${RUN_PREHYDRO_PAIR}" >> "$INITIAL_DIR/chunk_status.txt"
  echo "run_prehydro_only=${RUN_PREHYDRO_ONLY}" >> "$INITIAL_DIR/chunk_status.txt"
  echo "no_prehydro_alpha=${NO_PREHYDRO_ALPHA}" >> "$INITIAL_DIR/chunk_status.txt"
  echo "prehydro_alpha=${PREHYDRO_ALPHA}" >> "$INITIAL_DIR/chunk_status.txt"
  echo "broadening_k=${BROADENING_K}" >> "$INITIAL_DIR/chunk_status.txt"
  echo "do_moliere=${DO_MOLIERE}" >> "$INITIAL_DIR/chunk_status.txt"
  echo "moliere_mode=$([[ "${DO_MOLIERE}" == "true" ]] && echo legacy_resolved || echo disabled)" >> "$INITIAL_DIR/chunk_status.txt"
  echo "moliere_tables_sha256=${MOLIERE_TABLES_SHA256}" >> "$INITIAL_DIR/chunk_status.txt"
  if [[ -n "${HYDRO_EVENT_ID}" ]]; then
    echo "hydro_slot=${HYDRO_SLOT}" >> "$INITIAL_DIR/chunk_status.txt"
    echo "hydro_event_id=${HYDRO_EVENT_ID}" >> "$INITIAL_DIR/chunk_status.txt"
    echo "hydro_ncoll=${HYDRO_NCOLL}" >> "$INITIAL_DIR/chunk_status.txt"
    echo "hydro_dir=${HYDRO_DIR}" >> "$INITIAL_DIR/chunk_status.txt"
    echo "hydro_payload_key=${HYDRO_PAYLOAD_KEY}" >> "$INITIAL_DIR/chunk_status.txt"
    echo "hydro_payload_sha256=${HYDRO_PAYLOAD_SHA256}" >> "$INITIAL_DIR/chunk_status.txt"
  fi
  if [[ -d "$WORK/$RUN_NAME" ]]; then
    if [[ "${STORE_PREHYDRO_TABLE}" != "1" && "${STORE_PREHYDRO_TABLE}" != "true" && "${STORE_PREHYDRO_TABLE}" != "TRUE" ]]; then
      find "$WORK/$RUN_NAME" -type f -name prehydro_table.tsv -delete
    fi
    (
      cd "$WORK"
      find "$RUN_NAME" -type f \( \
        -name HYBRID_Hadrons.out -o \
        -name summary.tsv -o \
        -name stdout.log -o \
        -name hybrid_input.dat -o \
        -name setup_pythia.cmnd -o \
        -name prehydro_table.tsv -o \
        -name 'task_*_pair_summary.tsv' -o \
        -name 'task_*_prehydro_only_summary.tsv' \
      \) -print0 | tar --null --ignore-failed-read -czf "$INITIAL_DIR/chunk_output.tar.gz" --files-from -
    ) 2>/dev/null || tar -czf "$INITIAL_DIR/chunk_output.tar.gz" -C "$INITIAL_DIR" chunk_status.txt
  else
    tar -czf "$INITIAL_DIR/chunk_output.tar.gz" -C "$INITIAL_DIR" chunk_status.txt
  fi
  put "$INITIAL_DIR/chunk_status.txt" "status/${KIND}/chunk_${TASK_ID}.txt" >/dev/null 2>&1 || true
  put "$INITIAL_DIR/chunk_output.tar.gz" "outputs/${KIND}/chunk_${TASK_ID}.tar.gz" >/dev/null 2>&1 || true
  if [[ $rc -ne 0 && "${TOLERATE_CHUNK_FAILURE}" == "true" ]]; then
    echo "chunk ${TASK_ID} failed with exit code ${rc}; failure retained in EOS status" >&2
    exit 0
  fi
  exit $rc
}
trap finalize EXIT
mkdir -p "$WORK"
cd "$WORK"
fetch_from "$PAYLOAD_EOS_BASE" \
  payloads/pythia8315_alma9_install.tar.gz \
  pythia8315_alma9_install.tar.gz
fetch_from "$RUNTIME_PAYLOAD_EOS_BASE" \
  payloads/mmli_runtime_alma9.tar.gz \
  mmli_runtime_alma9.tar.gz
tar -xzf pythia8315_alma9_install.tar.gz
tar -xzf mmli_runtime_alma9.tar.gz
MOLIERE_TABLES_PATH=""
if [[ "${DO_MOLIERE}" == "true" ]]; then
  fetch_from "${MOLIERE_TABLES_EOS_BASE}" "${MOLIERE_TABLES_KEY}" moliere_a10_tables.zip
  echo "${MOLIERE_TABLES_SHA256}  moliere_a10_tables.zip" | sha256sum -c -
  mkdir -p runtime/moliere_tables
  python3 - moliere_a10_tables.zip runtime/moliere_tables <<'PY'
from pathlib import Path, PurePosixPath
import stat
import sys
import zipfile

archive_path, destination = map(Path, sys.argv[1:])
with zipfile.ZipFile(archive_path) as archive:
    members = archive.infolist()
    if not members:
        raise SystemExit("Moliere table archive is empty")
    for member in members:
        path = PurePosixPath(member.filename)
        mode = member.external_attr >> 16
        if (
            path.is_absolute()
            or ".." in path.parts
            or not path.parts
            or path.parts[0] != "a10_tables"
            or stat.S_ISLNK(mode)
        ):
            raise SystemExit(f"unsafe Moliere archive member: {member.filename}")
    archive.extractall(destination)

root = destination / "a10_tables"
for species in ("quark_tables", "gluon_tables"):
    tables = list((root / species).glob("m*x_g_*_d_*_n_*.dat"))
    if len(tables) != 476 or any(path.stat().st_size <= 0 for path in tables):
        raise SystemExit(
            f"{species}: expected 476 nonempty Moliere tables, found {len(tables)}"
        )
PY
  MOLIERE_TABLES_PATH="$WORK/runtime/moliere_tables/a10_tables"
fi
if [[ "${KIND}" == "aa" && -n "${AA_TASK_MANIFEST}" ]]; then
  if [[ ! -f "${AA_TASK_MANIFEST}" ]]; then
    echo "AA task manifest not found in runtime: ${AA_TASK_MANIFEST}" >&2
    exit 1
  fi
  IFS=$'\t' read -r HYDRO_SLOT HYDRO_EVENT_ID HYDRO_NCOLL HYDRO_DIR HYDRO_PAYLOAD_KEY HYDRO_PAYLOAD_SHA256 < <(
    python3 - "${AA_TASK_MANIFEST}" "${TASK_ID}" <<'PY'
import csv
import sys

path, requested = sys.argv[1], int(sys.argv[2])
with open(path, newline="") as handle:
    rows = [row for row in csv.DictReader(handle, delimiter="\t") if int(row["task_id"]) == requested]
if len(rows) != 1:
    raise SystemExit(f"{path}: expected one row for task {requested}, found {len(rows)}")
row = rows[0]
fields = (
    "hydro_slot",
    "hydro_event_id",
    "hydro_ncoll",
    "hydro_dir",
    "hydro_payload_key",
    "hydro_payload_sha256",
)
print("\t".join(row[field] for field in fields))
PY
  )
  if [[ -z "${HYDRO_PAYLOAD_KEY}" || ! "${HYDRO_PAYLOAD_SHA256}" =~ ^[0-9a-f]{64}$ ]]; then
    echo "Invalid hydro assignment for task ${TASK_ID}" >&2
    exit 1
  fi
  fetch_from "$PAYLOAD_EOS_BASE" \
    "payloads/${HYDRO_PAYLOAD_KEY}" \
    assigned_hydro.tar.gz
  echo "${HYDRO_PAYLOAD_SHA256}  assigned_hydro.tar.gz" | sha256sum -c -
  mkdir -p runtime/staged_hydro
  python3 - assigned_hydro.tar.gz runtime/staged_hydro "${HYDRO_DIR}" <<'PY'
from pathlib import PurePosixPath
import sys
import tarfile

archive_path, destination, expected_dir = sys.argv[1:]
with tarfile.open(archive_path, "r:gz") as archive:
    members = archive.getmembers()
    if not members:
        raise SystemExit("assigned hydro archive is empty")
    for member in members:
        path = PurePosixPath(member.name)
        if path.is_absolute() or ".." in path.parts or not path.parts or path.parts[0] != expected_dir:
            raise SystemExit(f"unsafe or unexpected hydro archive member: {member.name}")
        if member.issym() or member.islnk() or member.isdev():
            raise SystemExit(f"unsupported hydro archive member type: {member.name}")
    archive.extractall(destination)
PY
  test -s "runtime/staged_hydro/${HYDRO_DIR}/evolution_all_xyeta.dat"
  test -s "runtime/staged_hydro/${HYDRO_DIR}/NcollList.dat"
fi
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
  --aa-centrality-index "$AA_CENTRALITY_INDEX" \
  --no-prehydro-alpha "$NO_PREHYDRO_ALPHA" \
  --prehydro-alpha "$PREHYDRO_ALPHA" \
  --broadening-k "$BROADENING_K"
)
if [[ "$RUN_PREHYDRO_PAIR" == "true" ]]; then
  args+=(--run-prehydro-pair)
fi
if [[ "$RUN_PREHYDRO_ONLY" == "true" ]]; then
  args+=(--run-prehydro-only)
fi
if [[ "$DO_MOLIERE" == "true" ]]; then
  args+=(
    --do-moliere
    --moliere-tables-path "$MOLIERE_TABLES_PATH"
    --moliere-tables-sha256 "$MOLIERE_TABLES_SHA256"
  )
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
if [[ -n "$AA_TASK_MANIFEST" ]]; then
  args+=(--aa-task-manifest "$AA_TASK_MANIFEST")
fi
if [[ "$PREHYDRO_VISCOUS_ANCHOR" == "0" || "$PREHYDRO_VISCOUS_ANCHOR" == "false" || "$PREHYDRO_VISCOUS_ANCHOR" == "FALSE" ]]; then
  args+=(--no-prehydro-viscous-anchor)
elif [[ "$PREHYDRO_VISCOUS_ANCHOR" == "1" || "$PREHYDRO_VISCOUS_ANCHOR" == "true" || "$PREHYDRO_VISCOUS_ANCHOR" == "TRUE" ]]; then
  args+=(--prehydro-viscous-anchor)
fi
python3 "${args[@]}"

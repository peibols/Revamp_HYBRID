#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "usage: $0 ALPHA037_ROOT ALPHA0335_ROOT OUT_DIR EXPECTED_EVENTS" >&2
  exit 2
fi

alpha037_root=$(realpath "$1")
alpha0335_root=$(realpath "$2")
out_dir=$(realpath -m "$3")
expected_events=$4
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
python=${PYTHON:-python3}
SUBSTRUCTURE_REBIN_FACTOR=${SUBSTRUCTURE_REBIN_FACTOR:-2}

if [[ ! "${SUBSTRUCTURE_REBIN_FACTOR}" =~ ^[1-9][0-9]*$ ]]; then
  echo "SUBSTRUCTURE_REBIN_FACTOR must be a positive integer" >&2
  exit 1
fi

mkdir -p "$out_dir"

run_pair_analysis() {
  local root=$1
  local pair_name=$2
  local tag=$3
  local pt_min=$4
  local pt_max=$5
  local destination="$out_dir/$pair_name/$tag"
  local prefix="${pair_name}_${tag}"
  local args=(
    "$python" "$script_dir/plot_oo_jet_variables.py"
    --input-root "$root"
    --out-dir "$destination"
    --pt-min "$pt_min"
    --substructure-rebin-factor "$SUBSTRUCTURE_REBIN_FACTOR"
    --prefix "$prefix"
  )
  if [[ -n "$pt_max" ]]; then
    args+=(--pt-max "$pt_max")
  fi
  "${args[@]}"
}

merge_selection() {
  local tag=$1
  local pt_min=$2
  local pt_max=$3
  local prefix="oo5360_v3_matched_jet_variables_${tag}"
  local destination="$out_dir/merged/$tag"
  local args=(
    "$python" "$script_dir/plot_oo_v3_jet_variables.py"
    --hist-alpha037 "$out_dir/pair_alpha037/$tag/pair_alpha037_${tag}_histograms.tsv"
    --hist-alpha0335 "$out_dir/pair_alpha0335/$tag/pair_alpha0335_${tag}_histograms.tsv"
    --summary-alpha037 "$out_dir/pair_alpha037/$tag/pair_alpha037_${tag}_summary.tsv"
    --summary-alpha0335 "$out_dir/pair_alpha0335/$tag/pair_alpha0335_${tag}_summary.tsv"
    --pt-min "$pt_min"
    --substructure-rebin-factor "$SUBSTRUCTURE_REBIN_FACTOR"
    --expected-events "$expected_events"
    --out-dir "$destination"
    --prefix "$prefix"
  )
  if [[ -n "$pt_max" ]]; then
    args+=(--pt-max "$pt_max")
  fi
  "${args[@]}"
}

selections=(
  "pt20:20:"
  "pt30:30:"
  "pt20to30:20:30"
  "pt30to50:30:50"
  "pt50to80:50:80"
  "pt80plus:80:"
)

for selection in "${selections[@]}"; do
  IFS=: read -r tag pt_min pt_max <<<"$selection"
  run_pair_analysis "$alpha037_root" pair_alpha037 "$tag" "$pt_min" "$pt_max"
  run_pair_analysis "$alpha0335_root" pair_alpha0335 "$tag" "$pt_min" "$pt_max"
  merge_selection "$tag" "$pt_min" "$pt_max"
done

summary_dir="$out_dir/merged/summary"
mkdir -p "$summary_dir"
"$python" "$script_dir/summarize_oo_v3_jet_pt_slices.py" \
  --slice-metadata \
  "$out_dir/merged/pt20to30/oo5360_v3_matched_jet_variables_pt20to30_metadata.json" \
  "$out_dir/merged/pt30to50/oo5360_v3_matched_jet_variables_pt30to50_metadata.json" \
  "$out_dir/merged/pt50to80/oo5360_v3_matched_jet_variables_pt50to80_metadata.json" \
  "$out_dir/merged/pt80plus/oo5360_v3_matched_jet_variables_pt80plus_metadata.json" \
  --inclusive-metadata \
  "$out_dir/merged/pt20/oo5360_v3_matched_jet_variables_pt20_metadata.json" \
  --out-dir "$summary_dir" \
  --prefix oo5360_v3_matched_jet_pt_slice_summary

echo "V3 jet spectra complete: $out_dir"

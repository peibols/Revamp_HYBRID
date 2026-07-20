#!/bin/bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
validation_root="${repo_root}/validation/mmlhi"
workspace="$(cd "${repo_root}/.." && pwd)"

parent_repo="${PARENT_REPO:-${workspace}/wt_main_moliere_lres_integration_clean}"
pythia_home="${PYTHIA_HOME:-/data/yjlee/pythia/pythia8/pythia8315}"
tables_path="${MOLIERE_TABLES:-${workspace}/moliere_table_bundle/a10_tables}"
reference_run="${MMLHI_REFERENCE_RUN:-${workspace}/test/events_per_seed_hybrid_timing_20260505/mmli_pbpb_noelastic_events1_n1/no_moliere_no_wake/seed_860001}"
output_root="${MMLHI_VALIDATION_OUT:-${workspace}/test/mmhli_validation_current}"

hydro_file="${reference_run}/hydroinfoPlaintxtHuichaoFormat.dat"
tab_file="${reference_run}/TAb2LL.dat"

for path in "${parent_repo}" "${pythia_home}/include" "${pythia_home}/lib" \
            "${tables_path}" "${hydro_file}" "${tab_file}"; do
    if [[ ! -e "${path}" ]]; then
        echo "Required validation input is missing: ${path}" >&2
        exit 2
    fi
done

mkdir -p "${output_root}/config" "${output_root}/runs" "${output_root}/logs" "${output_root}/bin"

build_repo() {
    local repo="$1"
    (
        cd "${repo}"
        PYTHIA_INCLUDE="${pythia_home}/include" \
        PYTHIA_LIB="${pythia_home}/lib" \
        ./compiler.sh main
    )
}

render_config() {
    local output="$1"
    local seed="$2"
    local nev="$3"
    local heavy_mode="$4"
    local do_elastic="$5"
    local do_lres="$6"
    local lres_mode="$7"
    local c_res="$8"
    local dump_history="$9"
    local rpower="${10:-2.0}"
    local mode_b=false mode_c=false mode_d=false mode_e=false hadro_type=0
    if [[ "${do_elastic}" == true ]]; then hadro_type=1; fi
    case "${lres_mode}" in
        A) ;;
        B) mode_b=true ;;
        C) mode_c=true ;;
        D) mode_d=true ;;
        E) mode_e=true ;;
        *) echo "Unknown LRES/Moliere mode: ${lres_mode}" >&2; exit 2 ;;
    esac
    sed \
        -e "s|@SEED@|${seed}|g" \
        -e "s|@NEV@|${nev}|g" \
        -e "s|@HEAVY_MODE@|${heavy_mode}|g" \
        -e "s|@DO_ELASTIC@|${do_elastic}|g" \
        -e "s|@DO_LRES@|${do_lres}|g" \
        -e "s|@MODE_B@|${mode_b}|g" \
        -e "s|@MODE_C@|${mode_c}|g" \
        -e "s|@MODE_D@|${mode_d}|g" \
        -e "s|@MODE_E@|${mode_e}|g" \
        -e "s|@C_RES@|${c_res}|g" \
        -e "s|@RPOWER@|${rpower}|g" \
        -e "s|@HADRO_TYPE@|${hadro_type}|g" \
        -e "s|@TABLES_PATH@|${tables_path}/|g" \
        -e "s|@DUMP_HISTORY@|${dump_history}|g" \
        "${validation_root}/hybrid.input.in" > "${output}"
}

run_one() {
    local name="$1"
    local executable="$2"
    local config="$3"
    local pythia_card="$4"
    local run_dir="${output_root}/runs/${name}"
    local expected_events
    expected_events="$(awk -F= '$1 ~ /^Nev / {gsub(/[[:space:]]/, "", $2); print $2}' "${config}")"
    mkdir -p "${run_dir}"
    if [[ "${MMLHI_RESUME:-0}" == 1 && -s "${run_dir}/wall.txt" &&
          -s "${run_dir}/result_Hadrons.out" && -s "${run_dir}/result_Partons.out" ]]; then
        local completed_events
        completed_events="$(grep -c '^# event ' "${run_dir}/result_Hadrons.out")"
        if [[ "${completed_events}" == "${expected_events}" ]]; then
            return
        fi
    fi
    rm -f "${run_dir}/result_Hadrons.out" "${run_dir}/result_Partons.out" \
          "${run_dir}/result_history.tsv" "${run_dir}/run.log" "${run_dir}/wall.txt"
    ln -sfn "${config}" "${run_dir}/hybrid_input.dat"
    ln -sfn "${pythia_card}" "${run_dir}/setup_pythia.cmnd"
    ln -sfn "${hydro_file}" "${run_dir}/hydroinfoPlaintxtHuichaoFormat.dat"
    ln -sfn "${tab_file}" "${run_dir}/TAb2LL.dat"
    (
        cd "${run_dir}"
        /usr/bin/time -f "wall_seconds=%e" -o wall.txt \
            "${executable}" hybrid_input.dat > run.log 2>&1
    )
}

assert_same_output() {
    local first="$1"
    local second="$2"
    cmp "${output_root}/runs/${first}/result_Hadrons.out" \
        "${output_root}/runs/${second}/result_Hadrons.out"
    cmp "${output_root}/runs/${first}/result_Partons.out" \
        "${output_root}/runs/${second}/result_Partons.out"
}

run_parent_child() {
    local label="$1"
    local config="$2"
    local card="$3"
    run_one "parent_${label}" "${parent_repo}/main" "${config}" "${card}" &
    local parent_pid=$!
    run_one "child_${label}" "${repo_root}/main" "${config}" "${card}" &
    local child_pid=$!
    wait "${parent_pid}"
    wait "${child_pid}"
    assert_same_output "parent_${label}" "child_${label}"
}

echo "Building parent and child with pinned PYTHIA 8.315"
build_repo "${parent_repo}" > "${output_root}/logs/build_parent.log" 2>&1
build_repo "${repo_root}" > "${output_root}/logs/build_child.log" 2>&1

echo "Running heavy-kernel unit and sanitizer tests"
"${repo_root}/test/run_heavy_quark_unit.sh" > "${output_root}/logs/unit.log" 2>&1
g++ -std=c++17 -O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined \
    -I"${repo_root}" \
    "${repo_root}/test/test_heavy_quark_energy_loss.cc" \
    "${repo_root}/HeavyQuarkEnergyLoss.cc" "${repo_root}/Random.cc" \
    -o "${output_root}/bin/test_heavy_quark_energy_loss_sanitized"
ASAN_OPTIONS=detect_leaks=1 \
    "${output_root}/bin/test_heavy_quark_energy_loss_sanitized" \
    > "${output_root}/logs/unit_sanitized.log" 2>&1

generic_card="${validation_root}/setup_generic.cmnd"
charm_card="${validation_root}/setup_charm.cmnd"
bottom_card="${validation_root}/setup_bottom.cmnd"

echo "Running exact feature-off closure against the MMLI parent"
render_config "${output_root}/config/closure_standard.input" 860001 1 0 false false A 1.0 false
render_config "${output_root}/config/closure_lres.input" 860001 1 0 false true A 1.0 false
render_config "${output_root}/config/closure_moliere.input" 860001 1 0 true false A 1.0 false
run_parent_child closure_standard "${output_root}/config/closure_standard.input" "${generic_card}"
run_parent_child closure_lres "${output_root}/config/closure_lres.input" "${generic_card}"
run_parent_child closure_moliere "${output_root}/config/closure_moliere.input" "${generic_card}"

for lres_mode in A B C D E; do
    mode_lower="$(printf '%s' "${lres_mode}" | tr '[:upper:]' '[:lower:]')"
    config="${output_root}/config/closure_mode_${mode_lower}.input"
    render_config "${config}" 0 7 0 true true "${lres_mode}" 0.0 true
    run_parent_child "closure_mode_${mode_lower}" "${config}" "${generic_card}"
done

# The compact matrix reaches real candidates in C-E. These two targeted cases
# additionally require the resolving branch itself, including Mode E's live
# recursive tree update, to execute identically in parent and child.
render_config "${output_root}/config/closure_mode_c_resolving.input" 0 14 0 true true C 0.0 true
run_parent_child closure_mode_c_resolving \
    "${output_root}/config/closure_mode_c_resolving.input" "${generic_card}"
render_config "${output_root}/config/closure_mode_e_resolving.input" 0 14 0 true true E 0.0 true 0.2
run_parent_child closure_mode_e_resolving \
    "${output_root}/config/closure_mode_e_resolving.input" "${generic_card}"

for resolving_run in child_closure_mode_d child_closure_mode_c_resolving \
                     child_closure_mode_e_resolving; do
    resolving_count="$(sed -n 's/.*n_unresolved_resolving_scatters= \([0-9][0-9]*\).*/\1/p' \
        "${output_root}/runs/${resolving_run}/run.log" | tail -1)"
    if [[ -z "${resolving_count}" || "${resolving_count}" -le 0 ]]; then
        echo "${resolving_run} did not exercise a resolving scattering" >&2
        exit 3
    fi
done

python3 "${validation_root}/validate_history.py" C \
    "${output_root}/runs/child_closure_mode_c_resolving/result_history.tsv" \
    > "${output_root}/logs/history_mode_c.log"
python3 "${validation_root}/validate_history.py" D \
    "${output_root}/runs/child_closure_mode_d/result_history.tsv" \
    > "${output_root}/logs/history_mode_d.log"
python3 "${validation_root}/validate_history.py" E \
    "${output_root}/runs/child_closure_mode_e_resolving/result_history.tsv" \
    > "${output_root}/logs/history_mode_e.log"

echo "Running feature-on charm and bottom coverage"
for heavy_mode in 1 2 3; do
    config="${output_root}/config/charm_mode_${heavy_mode}.input"
    render_config "${config}" 870001 1 "${heavy_mode}" false false A 1.0 false
    run_one "charm_mode_${heavy_mode}" "${repo_root}/main" "${config}" "${charm_card}"
done

render_config "${output_root}/config/charm_lres_mode_2.input" 870001 1 2 false true A 1.0 false
run_one charm_lres_mode_2 "${repo_root}/main" \
    "${output_root}/config/charm_lres_mode_2.input" "${charm_card}"

render_config "${output_root}/config/charm_moliere_mode_2.input" 870001 100 2 true false A 1.0 true
run_one charm_moliere_mode_2 "${repo_root}/main" \
    "${output_root}/config/charm_moliere_mode_2.input" "${charm_card}"

for lres_mode in A B C D E; do
    mode_lower="$(printf '%s' "${lres_mode}" | tr '[:upper:]' '[:lower:]')"
    config="${output_root}/config/charm_lres_moliere_mode_${mode_lower}.input"
    render_config "${config}" 870001 20 2 true true "${lres_mode}" 0.0 true
    run_one "charm_lres_moliere_mode_${mode_lower}" "${repo_root}/main" "${config}" "${charm_card}"
done

render_config "${output_root}/config/charm_lres_moliere_mode_e_hard100.input" \
    870001 100 2 true true E 0.0 true
run_one charm_lres_moliere_mode_e_hard100 "${repo_root}/main" \
    "${output_root}/config/charm_lres_moliere_mode_e_hard100.input" "${charm_card}"

render_config "${output_root}/config/bottom_mode_2.input" 880001 1 2 false false A 1.0 false
render_config "${output_root}/config/bottom_moliere_mode_2.input" 880001 1 2 true false A 1.0 true
run_one bottom_mode_2 "${repo_root}/main" \
    "${output_root}/config/bottom_mode_2.input" "${bottom_card}"
run_one bottom_moliere_mode_2 "${repo_root}/main" \
    "${output_root}/config/bottom_moliere_mode_2.input" "${bottom_card}"

hard_charm_count="$(awk 'NF == 6 && ($5 == 4 || $5 == -4) && $6 == 1 {n++} END {print n+0}' \
    "${output_root}/runs/charm_moliere_mode_2/result_Partons.out")"
if [[ "${hard_charm_count}" -le 0 ]]; then
    echo "The charm Moliere coverage run did not retain a hard-scattered charm" >&2
    exit 3
fi

echo "Checking deterministic reruns"
run_one charm_mode_2_repeat "${repo_root}/main" \
    "${output_root}/config/charm_mode_2.input" "${charm_card}"
assert_same_output charm_mode_2 charm_mode_2_repeat
run_one charm_lres_moliere_mode_e_repeat "${repo_root}/main" \
    "${output_root}/config/charm_lres_moliere_mode_e.input" "${charm_card}"
assert_same_output charm_lres_moliere_mode_e charm_lres_moliere_mode_e_repeat

echo "Checking required-parameter failure"
invalid_config="${output_root}/config/invalid_missing_lambda.input"
render_config "${invalid_config}" 1 1 2 false false A 1.0 false
sed -i '/^heavy_quark_lambda =/d' "${invalid_config}"
invalid_dir="${output_root}/runs/invalid_missing_lambda"
mkdir -p "${invalid_dir}"
set +e
(cd "${invalid_dir}" && "${repo_root}/main" "${invalid_config}" > run.log 2>&1)
invalid_status=$?
set -e
if [[ ${invalid_status} -ne 3 ]] || ! grep -q "heavy_quark_lambda must be set explicitly" "${invalid_dir}/run.log"; then
    echo "Missing-lambda validation did not fail with the expected diagnostic" >&2
    exit 3
fi

python3 "${validation_root}/summarize.py" "${output_root}/runs" \
    > "${output_root}/summary.tsv"

{
    echo "parent_commit=$(git -C "${parent_repo}" rev-parse HEAD)"
    echo "child_commit=$(git -C "${repo_root}" rev-parse HEAD)"
    echo "pythia_home=${pythia_home}"
    echo "tables_path=${tables_path}"
    echo "reference_run=${reference_run}"
} > "${output_root}/provenance.txt"

echo "MMLHI validation passed. Summary: ${output_root}/summary.tsv"

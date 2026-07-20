#!/bin/bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="$(mktemp -d "${repo_root}/test/tmp_modee_policy_unit.XXXXXX")"
trap 'rm -rf "${build_dir}"' EXIT

"${CXX:-g++}" -std=c++17 -Wall -Wextra -Werror \
    -I"${repo_root}" "${repo_root}/test/test_modee_policy.cc" \
    -o "${build_dir}/test_modee_policy"
"${build_dir}/test_modee_policy"

#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
binary="${TMPDIR:-/tmp}/mmli_config_strict_unit"
input="${TMPDIR:-/tmp}/mmli_config_strict_unit.input"

g++ -std=c++17 -Wall -Wextra -pedantic -I"${repo_dir}" \
  "${repo_dir}/test/test_config_strict.cc" "${repo_dir}/Config.cc" -o "${binary}"
"${binary}" "${input}"

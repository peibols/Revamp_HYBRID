#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
binary="${TMPDIR:-/tmp}/mmli_response_ledger_unit"

g++ -std=c++17 -Wall -Wextra -pedantic -I"${repo_dir}" \
  "${repo_dir}/test/test_response_ledger.cc" \
  "${repo_dir}/Parton.cc" "${repo_dir}/Quench.cc" -o "${binary}"
"${binary}"
